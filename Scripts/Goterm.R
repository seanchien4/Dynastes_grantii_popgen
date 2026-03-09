# Go enrichment analysis
library(dplyr)
library(readr)
library(tidyr)
library(topGO)
library(Rgraphviz)
library(Gviz)

annotated_candidates <- read.csv("../Data/GEA/annotated_high_confidence_candidates.csv")
sca <- sort(unique(annotated_candidates$chromosome))[1:11]

PD.snps <- annotated_candidates[grepl("PD", annotated_candidates$lfmm_association) & annotated_candidates$chromosome %in% sca,][,c(3,4)]
write.table(PD.snps,
            file = "../Data/GEA/PD_SNP.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
PW.snps <- annotated_candidates[grepl("PW", annotated_candidates$lfmm_association) & annotated_candidates$chromosome %in% sca,][,c(3,4)]
write.table(PW.snps,
            file = "../Data/GEA/PW_SNP.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
PC.snps <- annotated_candidates[grepl("PC", annotated_candidates$lfmm_association) & annotated_candidates$chromosome %in% sca,][,c(3,4)]
write.table(PC.snps,
            file = "../Data/GEA/PC_SNP.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)



interproscan_file <- read_tsv('../Data/GEA/braker_isoform_removed_nostop.aa.tsv', col_names = FALSE, show_col_types = FALSE)
go_col_index <- 14 
universe_genes <- interproscan_file %>%
  pull(X1) %>%                     
  gsub("\\.t\\d+$", "", .) %>%
  unique()
write.table(universe_genes, 
            file = '../Data/GEA/universe_output_file', 
            row.names = FALSE, 
            col.names = FALSE, 
            quote = FALSE)

gene_to_go_map <- interproscan_file %>%
  dplyr::select(gene_id_raw = X1, go_terms = all_of(go_col_index)) %>%
  
  filter(!is.na(go_terms)) %>%
  
  mutate(gene_id = gsub("\\.t\\d+$", "", gene_id_raw)) %>%
  
  separate_rows(go_terms, sep = "\\|") %>%
  
  # Remove invalid hyphen rows ***
  filter(go_terms != "-") %>%
  
  # Remove the (InterPro) suffix ***
  mutate(go_terms = gsub("\\(.*?\\)", "", go_terms)) %>%
  
  dplyr::select(gene_id, go_terms) %>%
  
  distinct()
write.table(gene_to_go_map, 
            file = '../Data/GEA/map_output_file', 
            sep = "\t", 
            row.names = FALSE, 
            col.names = FALSE, 
            quote = FALSE)


# read files
candidate_gene_ids <- read.table('../Data/GEA/PW_gene.list') 
candidate_gene_ids <- candidate_gene_ids$V1

# B) The gene-to-GO mapping file you created from InterProScan
gene_to_go_map <- read_tsv("../Data/GEA/map_output_file",col_names = FALSE, comment = "#")
colnames(gene_to_go_map) <- c("gene_id", "go_term")

# C) The background/universe gene list
universe_gene_ids <- readLines("../Data/GEA/universe_output_file")

# --- topGO requires a specific format for the mapping ---
# It needs a named list where names are gene IDs and values are GO terms.
# This code converts your data frame into that format.
geneID2GO <- tapply(gene_to_go_map$go_term, gene_to_go_map$gene_id, function(x) x)

# --- topGO also needs a list of all genes with a 1/0 for candidate status ---
# This is a named factor where names are all genes and values are 1 if the
# gene is a candidate, and 0 otherwise.
gene_list_for_topgo <- factor(as.integer(universe_gene_ids %in% candidate_gene_ids))
names(gene_list_for_topgo) <- universe_gene_ids

# 2. CREATE THE topGO DATA OBJECT
# --------------------------------
# We will focus on "BP" - Biological Process, which is usually the most informative.
# You can also run this for "MF" (Molecular Function) and "CC" (Cellular Component).

GOdata <- new("topGOdata",
              ontology = "BP",  # Can be "BP", "MF", or "CC"
              allGenes = gene_list_for_topgo,
              annot = annFUN.gene2GO,
              gene2GO = geneID2GO)

# The GOdata object now contains all your information in the correct structure.


# 3. RUN THE ENRICHMENT TEST
# --------------------------
# topGO offers several algorithms. 'weight01' is a good modern choice that
# accounts for the GO hierarchy. 'classic' is a simpler gene-by-gene count.
# We will use Fisher's exact test.

resultFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")

# The resultFisher object now contains the p-values for every GO term.

# 4. VIEW THE RESULTS
# --------------------
# GenTable creates a data frame of the results.
# We ask for the top 50 significant terms.
# The p-values here are NOT corrected for multiple testing yet.
results_table <- GenTable(GOdata, 
                          Fisher = resultFisher,
                          orderBy = "Fisher", 
                          ranksOf = "Fisher", 
                          topNodes = 50)

# --- Add Multiple-Testing Correction (FDR) ---
# This is a crucial step!
results_table$FDR <- p.adjust(results_table$Fisher, method = "fdr")

# Now, filter for your truly significant results
significant_results <- results_table %>%
  filter(FDR < 0.05)

# View your final, high-confidence table
print(significant_results)
significant_go_ids <- significant_results$GO.ID
genes_for_significant_terms <- list()

for (go_id in significant_go_ids) {
  # Get all genes in the universe associated with this GO term
  all_genes_in_term <- genesInTerm(GOdata, go_id)[[1]]
  
  # Now, find which of those genes are also in your candidate list
  # This finds the intersection of the two lists.
  candidate_genes_in_term <- intersect(all_genes_in_term, candidate_gene_ids)
  
  # Store the list of candidate genes for this significant term
  genes_for_significant_terms[[go_id]] <- candidate_genes_in_term
}

# 4. PRINT AND INSPECT THE RESULTS
# ---------------------------------
# `genes_for_significant_terms` is now a named list containing your answer.

print(genes_for_significant_terms)

########
# Plot # 
########
plot_data <- significant_results[order(significant_results$Significant, decreasing = FALSE), ]

# Create a -log10(FDR) column for the color scale
plot_data$log10FDR <- -log10(plot_data$FDR)

# 2. CREATE A CUSTOM COLOR PALETTE
# --------------------------------
# We will create a color gradient from light blue to dark blue.
# `colorRampPalette` creates a function that generates colors along a gradient.
num_colors <- 100
color_palette_func <- colorRampPalette(c("#DEEBF7", "#08519C"))
palette <- color_palette_func(num_colors)

# Assign a color to each bar based on its -log10(FDR) value
# We scale the log10FDR values to be an index from 1 to num_colors
fdr_range <- range(plot_data$log10FDR)
# Handle case where there's only one term to avoid division by zero
if (diff(fdr_range) == 0) {
  color_indices <- rep(num_colors, nrow(plot_data))
} else {
  color_indices <- 1 + round((num_colors - 1) * (plot_data$log10FDR - fdr_range[1]) / (fdr_range[2] - fdr_range[1]))
}
bar_colors <- palette[color_indices]

bar_centers <- barplot(
  height = plot_data$Significant,
  names.arg = plot_data$Term,
  cex.names = 0.6,
  horiz = TRUE,
  las = 1,
  col = bar_colors,
  border = NA,
  main = "Enriched GO Terms for Climate Adaptation Candidates",
  xlab = "Number of Candidate Genes",
)

# 
num_sig_terms <- nrow(significant_results)
showSigOfNodes(GOdata, 
               score(resultFisher), 
               firstSigNodes = num_sig_terms, 
               useInfo = 'all')

