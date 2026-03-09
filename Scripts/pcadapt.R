library(pcadapt)
library(vcfR)
library(qvalue)
library(gtools)


genotypes <- read.pcadapt('../Data/dynastes_ped.ped', type = "ped")
pca_results <- pcadapt(input = genotypes, K = 5)
plot(pca_results, option = "screeplot")
populations <- c(rep('Lemmon',5),rep('Portal',5),rep('Utah',5),rep('Mogollon',12)) 

# Plot PC1 vs PC2
plot(pca_results, option = "scores", pop = populations, i = 1, j = 2)
optimal_K <- 4 
final_results <- pcadapt(input = genotypes, K = optimal_K)
plot(final_results, option = "manhattan")
p_values <- final_results$pvalues
# Convert p-values to q-values
q_values <- qvalue(p = p_values)$qvalues

alpha <- 0.001
outlier_indices <- which(q_values < alpha)
cat("Number of candidate SNPs found:", length(outlier_indices))
outlier_indices <- which(q_values < 0.05) # Or your chosen alpha

snp_metadata <- read.table("../Data/dynastes_ped.map", 
                           col.names = c("CHR", "SNP_ID", "Genetic_Dist", "Position"))
candidate_snps_info <- snp_metadata[outlier_indices, ]
cat("Found", nrow(candidate_snps_info), "candidate SNPs. Top hits:\n")
unique(candidate_snps_info$CHR)

# plot 
plot_data <- data.frame(
  CHR = snp_metadata$CHR,
  Position = snp_metadata$Position,
  p_value = final_results$pvalues
)

plot_data$log_p <- -log10(plot_data$p_value)
plot_data$log_p[is.infinite(plot_data$log_p)] <- max(plot_data$log_p[is.finite(plot_data$log_p)], na.rm = TRUE) * 1.1
plot_data <- plot_data[order(mixedsort(plot_data$CHR), plot_data$Position), ]

chr_lengths <- tapply(plot_data$Position, plot_data$CHR, max)
chr_offsets <- c(0, cumsum(as.numeric(chr_lengths))[-length(chr_lengths)])
names(chr_offsets) <- names(chr_lengths)

# Add the cumulative offset to each SNP's position
plot_data$genome_pos <- plot_data$Position + chr_offsets[plot_data$CHR]
# save the plot data
# write.table(plot_data, file = '../Data/pcadapt_plot.txt')
# 4. Prepare chromosome colors (adapting your code)
chr_levels <- unique(plot_data$CHR)
cols <- c("#3030D0", "#9090F0") # Your requested colors (Dark Blue, Light Blue)
chr_levels <- chr_levels[order(chr_levels)]
chr_col_map <- setNames(rep(cols, length.out = length(chr_levels)), chr_levels)
plot_data$color <- chr_col_map[plot_data$CHR]

# 5. Calculate the FDR threshold line position
alpha <- 0.05
# Find the maximum p-value that is still significant (q < alpha)
# This is more accurate than using a Bonferroni-corrected p-value
significant_p_values <- plot_data$p_value[q_values < alpha]
if (length(significant_p_values) > 0) {
  threshold_p <- max(significant_p_values, na.rm = TRUE)
  threshold_logp <- -log10(threshold_p)
} else {
  # If no SNPs are significant, place the line off-screen or at a default high value
  threshold_logp <- -log10(0.05 / nrow(plot_data)) # Fallback to simple Bonferroni
  cat("No SNPs found below FDR threshold. Drawing a default Bonferroni line.\n")
}

# Set plotting parameters for a single, clean plot
par(mar = c(5, 5, 4, 2) + 0.1) # Good margins for labels (bottom, left, top, right)

# Create the plot
plot(plot_data$genome_pos, plot_data$log_p,
     pch = 16, cex = 0.5,                                   
     col = plot_data$color,                                  
     ylim = c(0, max(plot_data$log_p, na.rm = TRUE) + 2),     
     xaxt = "n",                                            
     xaxs = "i", yaxs = "i",                                 
     xlab = "Genomic Position",                             
     ylab = expression(-log[10](italic(p))),
     main = "",
     cex.lab = 1.2, cex.axis = 1.1) 
chr_centers <- chr_offsets + chr_lengths / 2
# axis_labels <- gsub("scaffold_|_", "", names(chr_centers))
axis_labels <- c(1:11)
axis(1, at = chr_centers, labels = axis_labels, las = 1, tick = FALSE, cex.axis = 0.8)

# Add the significance threshold line
abline(h = -log10(alpha), col = "red", lty = 2, lwd = 1.5)

# -------------------------------------------------------------------------
# NEW CODE TO IDENTIFY AND COLOR THE PEAK
# -------------------------------------------------------------------------

# 1. Define the criteria for the peak region
peak_chromosome <- "JAROYF020000008.1"
peak_threshold <- 100    # A -log10(p) value threshold to define the peak

# 2. Add a new column to 'plot_data' for custom coloring
# First, use the existing chromosome color as the default
plot_data$peak_color <- plot_data$color # 'color' was created in the previous script

# 3. Identify the SNPs that meet our peak criteria
# These are the SNPs on the correct chromosome AND above our significance threshold
is_in_peak <- (plot_data$CHR == peak_chromosome) & (plot_data$log_p > peak_threshold)

# 4. Assign a new, highlight color to just those peak SNPs
# Let's choose a vibrant red color to make it stand out
plot_data$peak_color[is_in_peak] <- "red"

# Optional: Check how many SNPs you've highlighted
cat("Number of SNPs highlighted in the peak on chromosome", peak_chromosome, ":", sum(is_in_peak), "\n")