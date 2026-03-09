library(terra)
library(readr)
library(usdm)
library(dplyr)
library(vegan)
library(lfmm)
library(viridis)
gps_data <- read.csv("../Data/GEA/coords.txt")
coords <- gps_data[, c("lon", "lat")]
climate_files <- list.files(path = "../Data//wc2.1_30s_bio/",
                            pattern = "\\.tif$", 
                            full.names = TRUE)
climate_stack <- rast(climate_files)
names(climate_stack) <- paste0("BIO", 1:19)
print(climate_stack)

# Extract climate data
extracted_values <- terra::extract(climate_stack, coords)
extracted_values <- extracted_values[, -1] 
env_data <- cbind(gps_data, extracted_values)
head(env_data)

write.csv(env_data, file = "../Data/GEA/environmental_data_raw.csv")

env_data_raw <- read.csv("../Data/GEA/environmental_data_raw.csv")
all_bio_vars <- env_data_raw[, paste0("BIO", 1:19)]

# selected_bio_vars <- env_data_raw[, c("BIO1", "BIO4", "BIO5", "BIO6",
#                                      "BIO12", "BIO14", "BIO15", "BIO18")]
selected_bio_vars <- env_data_raw[,-c(1:4)]

# calculate variance inflation factors for each variable
# only keep the variable with VIF < 10

vif_results <- vifstep(selected_bio_vars, th = 10)
print(vif_results)
final_vars_names <- vif_results@results$Variables


# Create a new data frame with only the non-correlated variables,
# plus your essential sample identifiers.
final_env_data <- env_data_raw[, c("SampleID", "lat", "lon", final_vars_names)]
colnames(final_env_data) <- c("SampleID", "lat", "lon", "PD", "PW", "PC")
# Inspect the final dataset
head(final_env_data)
str(final_env_data)
# write.csv(final_env_data, "../Data/GEA/environmental_data_final_VIF.csv")

############
# Analysis #
############
raw_file_path <- "../Data/GEA/genotype_data.raw"
geno_raw <- read.delim(raw_file_path, sep = " ")

# process the genotype matrix
# The .raw file has 6 columns of sample info (FID, IID, PAT, MAT, SEX, PHENOTYPE)
# followed by columns for each SNP.
# We need to separate the sample info from the genotype data.

# Extract the sample IDs (the IID column is what we usually need)
sample_ids_from_plink <- geno_raw$IID

# Extract the genotype data (columns 7 to the end) and convert to a matrix
genotype_matrix <- as.matrix(geno_raw[, 7:ncol(geno_raw)])

# Set the row names of the matrix to be the sample IDs for easy matching
rownames(genotype_matrix) <- sample_ids_from_plink

# The column names are the SNP IDs from PLINK
cat("Loaded genotype matrix with", nrow(genotype_matrix), "samples and", ncol(genotype_matrix), "SNPs.\n")

# Handle missing data
# PLINK command coded missing genotypes as '9'. We need to convert these
# to NA and then impute them, which is a common requirement for RDA.

genotype_matrix[genotype_matrix == 9] <- NA

# Impute missing data with the mean allele frequency for that SNP (column-wise).
# This is a more accurate approach than using a single global mean.
for(i in 1:ncol(genotype_matrix)) {
  mean_val <- mean(genotype_matrix[, i], na.rm = TRUE)
  genotype_matrix[is.na(genotype_matrix[, i]), i] <- mean_val
}

# Verify that there are no NAs left
if(any(is.na(genotype_matrix))) {
  warning("Warning: NAs still present after imputation!")
} else {
  cat("Missing data successfully imputed.\n")
}


# Synchronize with env data
# The order of samples in your genotype matrix MUST match the order in your
# environmental data frame for the GEA models to work correctly.

# Load final, VIF-filtered environmental data
env_data_final <- read.csv("../Data/GEA/environmental_data_final_VIF.csv")

# Get the required sample order from the environmental data
sample_order <- env_data_final$SampleID

# Check if all samples are present in both datasets
if (!all(sample_order %in% rownames(genotype_matrix))) {
  stop("ERROR: Not all samples from the environmental file are in the genotype file!")
}

# Reorder the rows of the genotype matrix to match the environmental data
genotype_matrix <- genotype_matrix[sample_order, ]

# Prepare the final obj for GEA 
#   1. genotype_matrix
#   2. env_predictors (created from env_data_final)

env_predictors <- env_data_final %>% select(PD,PW,PC)

# Let's inspect our final, ready-to-use objects
print("Dimensions of final genotype matrix:")
print(dim(genotype_matrix))
print("Dimensions of final environmental predictors:")
print(dim(env_predictors))

# Run the RDA 
# The formula specifies that the genotype matrix is the response variable,
# and all variables in env_predictors are explanatory.
# scale = TRUE is important as it standardizes the environmental variables.
rda_model <- rda(genotype_matrix ~ ., data = env_predictors, scale = TRUE)

# Model significance
# We test the significance of the overall model using a permutation test.
# This can take a few minutes.
anova_result <- anova.cca(rda_model, parallel = getOption("mc.cores"))


# Also, check the proportion of variance explained.
r_squared <- RsquareAdj(rda_model)
cat("\nAdjusted R-squared for RDA model:", r_squared$adj.r.squared, "\n")
summary(rda_model)$concont

# Identify outlier loci
# Outliers are SNPs with strong loadings on the RDA axes. These are our candidates.
# We'll define outliers as those falling > 3 standard deviations from the mean loading.

get_rda_outliers <- function(rda_object, axis_number, num_sds = 3) {
  loadings <- scores(rda_object, choices = axis_number, display = "species")
  mean_loading <- mean(loadings)
  sd_loading <- sd(loadings)
  threshold_upper <- mean_loading + num_sds * sd_loading
  threshold_lower <- mean_loading - num_sds * sd_loading
  outlier_indices <- which(loadings > threshold_upper | loadings < threshold_lower)
  outlier_names <- rownames(loadings)[outlier_indices]
  return(outlier_names)
}

# Get outliers from the first two RDA axes (these usually explain the most variance)
rda_candidates_ax1 <- get_rda_outliers(rda_model, 1)
rda_candidates_ax2 <- get_rda_outliers(rda_model, 2)

# Combine into a single, unique list of candidate SNPs from RDA
rda_candidates <- unique(c(rda_candidates_ax1, rda_candidates_ax2))
cat("Found", length(rda_candidates), "candidate SNPs using RDA (>3 SD).\n")

# -------------------------------------------------------------------
# LATENT FACTOR MIXED MODELS (LFMM)
# -------------------------------------------------------------------

# CHOOSE THE NUMBER OF LATENT FACTORS (K)

# This is a critical step. A PCA scree plot of your genetic data is a good way
# to estimate the number of major axes of population structure.

K <- 4

# Run LFMM 
# We must run the model separately for each environmental predictor.
lfmm_results <- list()

# This loop can be time-consuming.
for (env_var_name in colnames(env_predictors)) {
  
  cat("\nRunning LFMM for:", env_var_name, "...\n")
  
  # The LFMM model requires the environmental variable as a matrix
  current_env_matrix <- as.matrix(env_predictors[[env_var_name]])
  
  # Run the LFMM model
  lfmm_model <- lfmm_ridge(Y = genotype_matrix,
                           X = current_env_matrix,
                           K = K)
  
  # Perform association testing to get p-values.
  # The GIF calibration helps control inflation of test statistics.
  pv <- lfmm_test(Y = genotype_matrix,
                  X = current_env_matrix,
                  lfmm = lfmm_model,
                  calibrate = "gif")
  lfmm_results[[env_var_name]] <- pv$pvalue[,1]
}

# Adjust p-values 

# We must correct for testing thousands of SNPs. We use the False Discovery Rate (FDR).
fdr_threshold <- 0.05
lfmm_candidates <- list()

for (env_var_name in names(lfmm_results)) {
  
  # Adjust p-values using the Benjamini-Hochberg method (FDR)
  adjusted_p_values <- p.adjust(lfmm_results[[env_var_name]], method = "fdr")
  
  # Identify SNPs that pass the FDR threshold
  significant_snp_indices <- which(adjusted_p_values < fdr_threshold)
  
  if (length(significant_snp_indices) > 0) {
    # Get the names of the significant SNPs
    lfmm_candidates[[env_var_name]] <- colnames(genotype_matrix)[significant_snp_indices]
  }
}

# Print a summary of LFMM findings
cat("\n--- LFMM Results Summary ---\n")
for (env_var_name in names(lfmm_candidates)) {
  cat("Found", length(lfmm_candidates[[env_var_name]]), 
      "candidate SNPs for", env_var_name, "at FDR <", fdr_threshold, "\n")
}

######################
# Synthesize results # 
######################

# concordant candidates
# all unique SNPs found by LFMM across all variables.
all_lfmm_candidates <- unique(unlist(lfmm_candidates))

# find the intersection between the RDA list and the LFMM list.
concordant_candidates <- intersect(rda_candidates, all_lfmm_candidates)

cat("\n--- Analysis Synthesis ---\n")
cat("Total unique candidates from RDA:", length(rda_candidates), "\n")
cat("Total unique candidates from LFMM:", length(all_lfmm_candidates), "\n")
cat("Number of high-confidence CONCORDANT candidates:", length(concordant_candidates), "\n")

# This list is your final, most reliable set of candidate SNPs for adaptation.
print("Concordant SNP IDs:")
print(concordant_candidates)


#############
# RDA plots # 
#############
screeplot(rda_model)
# This plot shows samples, environmental variables, and candidate SNPs.
ind <- final_env_data$SampleID
cols <- viridis(4)
bg <- c(rep(cols[2],5), rep(cols[3], 5), rep(cols[4],5), rep(cols[1],12))

summary(rda_model)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
# axes  1 & 2 
plot(rda_model, scaling = 3, xlim=c(-3,12), ylim=c(-10,7), 
     main="", xlab = "RAD1 (46.8%)", ylab = "RAD2 (22.9%)")
points(rda_model, display="species", pch=21, cex=1, col="gray32", scaling=3)
points(rda_model, display="sites", pch=21, cex=1.3, col='black', bg=bg, scaling=3)
text(rda_model, scaling=3, display="bp", col="#0868ac", cex=1)
x = 9.5
y = 7
s = 0.8
cex=1.3
cols <- viridis(4)
points(x,y, bg = cols[1], col='black', pch =21, cex = cex )
text(x,y, pos = 4, labels = 'Mogollon Rim')
points(x,y-s, bg = cols[2], col='black', pch =21, cex = cex)
text(x,y-s, pos = 4, labels = 'Mt. Lemmon')
points(x,y-2*s, bg = cols[3], col='black', pch =21, cex = cex)
text(x,y-2*s, pos = 4, labels = 'Portal')
points(x,y-3*s, bg = cols[4], col='black', pch =21, cex = cex)
text(x,y-3*s, pos = 4, labels = 'Utah')

# axes  1 & 3 
plot(rda_model, scaling = 3, xlim=c(-4,13), ylim=c(-8,9), main="",choices=c(1,3),
     xlab = "RAD1 (46.8%)", ylab = "RAD3 (17.2%)")
points(rda_model, display="species", pch=21, cex=1, col="gray32", scaling=3,choices=c(1,3))
points(rda_model, display="sites", pch=21, cex=1.3, col='black', bg=bg, scaling=3,choices=c(1,3))
text(rda_model, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))
x = 9.2
y = 9
s = 0.7
cex=1.3
cols <- viridis(4)
points(x,y, bg = cols[1], col='black', pch =21, cex = cex )
text(x,y, pos = 4, labels = 'Mogollon Rim')
points(x,y-s, bg = cols[2], col='black', pch =21, cex = cex)
text(x,y-s, pos = 4, labels = 'Mt. Lemmon')
points(x,y-2*s, bg = cols[3], col='black', pch =21, cex = cex)
text(x,y-2*s, pos = 4, labels = 'Portal')
points(x,y-3*s, bg = cols[4], col='black', pch =21, cex = cex)
text(x,y-3*s, pos = 4, labels = 'Utah')

# plot SNPs
env_var_names <- colnames(env_predictors)


# Final list of candidate snps
final_candidates_df <- data.frame(
  snp_id = concordant_candidates
)
color_palette <- setNames(
  c("#E69F00", "#56B4E9", "#009E73"),
  env_var_names
)
color_palette["Multiple"] <- "#000000"
candidates_for_plot <- final_candidates_df %>%
  mutate(main_predictor = sapply(strsplit(lfmm_association, ";"), function(x) {
    if (length(x) == 1) return(x) else return("Multiple")
  })) %>%
  # Add a 'color' column based on the main predictor
  mutate(color = color_palette[main_predictor])

# 2. CREATE A COLOR VECTOR FOR ALL SNPS
all_rda_snps <- rownames(scores(rda_model, display = "species"))
snp_colors <- rep("gray80", length(all_rda_snps))
snp_sizes <- rep(0.5, length(all_rda_snps))

for (i in 1:nrow(candidates_for_plot)) {
  # Get the candidate SNP name and its assigned color
  snp_name <- candidates_for_plot$snp_id[i]
  snp_col <- candidates_for_plot$color[i]
  
  # Find the position (index) of this candidate in the full list of RDA SNPs
  match_index <- match(snp_name, all_rda_snps)
  
  # If a match is found, update the color and size at that position
  if (!is.na(match_index)) {
    snp_colors[match_index] <- snp_col
    snp_sizes[match_index] <- 0.5 # Make candidates larger
  }
}
r_squared <- RsquareAdj(rda_model)
plot(rda_model, type = "n", scaling = 3,
     main = "",
     xlab = paste0("RDA1 (46.8%)"),
     ylab = paste0("RDA2 (22.9%)"),
     xlim = c(-0.25, 0.25), ylim = c(-0.25,0.25))
points(rda_model, display = "species", pch = 21,
       bg = snp_colors, col = 'gray32',cex = 1, scaling = 3)
points(rda_model, display = "sites", pch = 4, col = "gray30", cex = 0.8, scaling = 3)
text(rda_model, display = "bp", col = "#0868ac", cex = 1, scaling = 3, font=1)
x = -0.29
y = 0.25
s = 0.02
cex=1
cols <- c("#E69F00", "#56B4E9", "#009E73")
points(x,y, bg = cols[1], col='black', pch =21, cex = cex )
text(x,y, pos = 4, labels = 'PD')
points(x,y-s, bg = cols[2], col='black', pch =21, cex = cex)
text(x,y-s, pos = 4, labels = 'PW')
points(x,y-2*s, bg = cols[3], col='black', pch =21, cex = cex)
text(x,y-2*s, pos = 4, labels = 'PC')
points(x,y-3*s, bg = '#000000', col='black', pch =21, cex = cex)
text(x,y-3*s, pos = 4, labels = 'Mixed')

# 1 & 3
plot(rda_model, type = "n", scaling = 3,
     main = "",
     xlab = paste0("RDA1 (46.8%)"),
     ylab = paste0("RDA3 (17.2%)"),
     xlim = c(-0.35, 0.35), ylim = c(-0.35,0.35),choices=c(1,3))
points(rda_model, display = "species", pch = 21,
       bg = snp_colors, col = 'gray32',cex = 1, scaling = 3,choices=c(1,3))
text(rda_model, display = "bp", col = "#0868ac", cex = 1, scaling = 3, font=1,choices=c(1,3))
# legend 
x = -0.4
y = 0.34
s = 0.03
cex=1
cols <- c("#E69F00", "#56B4E9", "#009E73")
points(x,y, bg = cols[1], col='black', pch =21, cex = cex )
text(x,y, pos = 4, labels = 'PD')
points(x,y-s, bg = cols[2], col='black', pch =21, cex = cex)
text(x,y-s, pos = 4, labels = 'PW')
points(x,y-2*s, bg = cols[3], col='black', pch =21, cex = cex)
text(x,y-2*s, pos = 4, labels = 'PC')
points(x,y-3*s, bg = '#000000', col='black', pch =21, cex = cex)
text(x,y-3*s, pos = 4, labels = 'Mixed')

##############################
#### RDA plot for manuscript #
##############################
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
# axes  1 & 2 
plot(rda_model, scaling = 3, xlim=c(-7,14), ylim=c(-10,7), 
     main="", xlab = "RAD1 (46.8%)", ylab = "RAD2 (22.9%)")
points(rda_model, display="species", pch=21, cex=1, col="gray32", scaling=3)
points(rda_model, display="sites", pch=21, cex=1, col='black', bg=bg, scaling=3)
text(rda_model, scaling=3, display="bp", col="#0868ac", cex = 0.7)
title(main = "(A)", adj = 0, line = 0.5)
x = 6.8
y = 8
s = 1
cex.point = 0.8
cex.text = 0.8
pch = 16 
cols <- viridis(4)
points(x,y, col = cols[1], pch = pch, cex = cex.point )
text(x,y, pos = 4, labels = 'Mogollon Rim', cex = cex.text)
points(x,y-s, col = cols[2], pch = pch, cex = cex.point)
text(x,y-s, pos = 4, labels = 'Mt. Lemmon', cex = cex.text)
points(x,y-2*s, col = cols[3], pch = pch, cex = cex.point)
text(x,y-2*s, pos = 4, labels = 'Portal', cex = cex.text)
points(x,y-3*s, col = cols[4], pch = pch, cex = cex.point)
text(x,y-3*s, pos = 4, labels = 'Utah', cex = cex.text)

# SNPs
plot(rda_model, type = "n", scaling = 3,
     main = "",
     xlab = paste0("RDA1 (46.8%)"),
     ylab = paste0("RDA2 (22.9%)"),
     xlim = c(-0.25, 0.25), ylim = c(-0.25,0.25))
points(rda_model, display = "species", pch = 21,
       bg = snp_colors, col = 'gray32',cex = 1, scaling = 3)
points(rda_model, display = "sites", pch = 4, col = "gray30", cex = 0.8, scaling = 3)
text(rda_model, display = "bp", col = "#0868ac", cex = 0.7, scaling = 3, font=1)
title(main = "(B)", adj = 0, line = 0.5)
x = -0.26
y = 0.24
s = 0.03
cols <- c("#E69F00", "#56B4E9", "#009E73")
points(x,y, bg = cols[1], pch =21, cex = cex.point )
text(x,y, pos = 4, labels = 'PD',cex = cex.text)
points(x,y-s, bg = cols[2], pch =21, cex = cex.point)
text(x,y-s, pos = 4, labels = 'PW', cex = cex.text)
points(x,y-2*s, bg = cols[3], pch =21, cex = cex.point)
text(x,y-2*s, pos = 4, labels = 'PC', cex = cex.text)
points(x,y-3*s, bg = '#000000', pch =21, cex = cex.point)
text(x,y-3*s, pos = 4, labels = 'Mixed', cex = cex.text)

# axes  1 & 3 
plot(rda_model, scaling = 3, xlim=c(-6,14), ylim=c(-8,9), main="",choices=c(1,3),
     xlab = "RAD1 (46.8%)", ylab = "RAD3 (17.2%)")
points(rda_model, display="species", pch=21, cex = 1, col="gray32", scaling=3,choices=c(1,3))
points(rda_model, display="sites", pch=21, cex = 1, col='black', bg=bg, scaling=3,choices=c(1,3))
text(rda_model, scaling=3, display="bp", col="#0868ac", cex = 0.7, choices=c(1,3))
title(main = "(C)", adj = 0, line = 0.5)
x = 7
y = 9.5
s = 0.9
cols <- viridis(4)
points(x,y, col = cols[1], pch = pch, cex = cex.point )
text(x,y, pos = 4, labels = 'Mogollon Rim', cex = cex.text)
points(x,y-s, col = cols[2], pch = pch, cex = cex.point)
text(x,y-s, pos = 4, labels = 'Mt. Lemmon', cex = cex.text)
points(x,y-2*s, col = cols[3], pch = pch, cex = cex.point)
text(x,y-2*s, pos = 4, labels = 'Portal', cex = cex.text)
points(x,y-3*s, col = cols[4], pch = pch, cex = cex.point)
text(x,y-3*s, pos = 4, labels = 'Utah', cex = cex.text)
# SNPs
plot(rda_model, type = "n", scaling = 3,
     main = "",
     xlab = paste0("RDA1 (46.8%)"),
     ylab = paste0("RDA3 (17.2%)"),
     xlim = c(-0.35, 0.35), ylim = c(-0.35,0.35),choices=c(1,3))
points(rda_model, display = "species", pch = 21,
       bg = snp_colors, col = 'gray32',cex = 1, scaling = 3,choices=c(1,3))
text(rda_model, display = "bp", col = "#0868ac", cex = 0.7 , scaling = 3, font=1,choices=c(1,3))
title(main = "(D)", adj = 0, line = 0.5)
x = -0.38
y = 0.34
s = 0.04
cex=1
cols <- c("#E69F00", "#56B4E9", "#009E73")
points(x,y, bg = cols[1], pch =21, cex = cex.point )
text(x,y, pos = 4, labels = 'PD', cex = cex.text)
points(x,y-s, bg = cols[2], pch =21, cex = cex.point)
text(x,y-s, pos = 4, labels = 'PW', cex = cex.text)
points(x,y-2*s, bg = cols[3], pch =21, cex = cex.point)
text(x,y-2*s, pos = 4, labels = 'PC', cex = cex.text)
points(x,y-3*s, bg = '#000000', pch =21, cex = cex.point)
text(x,y-3*s, pos = 4, labels = 'Mixed', cex = cex.text)

# Optional: Add which variable(s) each SNP was associated with in LFMM
lfmm_associations <- sapply(concordant_candidates, function(snp) {
  associated_vars <- c()
  for (var_name in names(lfmm_candidates)) {
    if (snp %in% lfmm_candidates[[var_name]]) {
      associated_vars <- c(associated_vars, var_name)
    }
  }
  return(paste(associated_vars, collapse = ";"))
})
final_candidates_df$lfmm_association <- lfmm_associations

# Save to CSV
write.csv(final_candidates_df, "../Data/GEA/high_confidence_candidate_snps.csv", row.names = FALSE)
# annotation
candidate_snps <- read_csv("../Data/GEA/high_confidence_candidate_snps.csv")
map_file <- read_tsv("../Data/GEA/genotype_data.map", col_names = FALSE)
colnames(map_file) <- c("chromosome", "snp_id_plink", "genetic_dist", "position")
raw_header <- geno_raw
raw_snp_names <- colnames(raw_header)[7:length(colnames(raw_header))]
if (length(raw_snp_names) != nrow(map_file)) {
  stop("FATAL ERROR: The number of SNPs in the .raw file header (", length(raw_snp_names), 
       ") does not match the number of rows in the .map file (", nrow(map_file), 
       ").\nCheck that you used the exact same filters for both PLINK commands.")
}
snp_lookup_table <- data.frame(
  snp_id_r = raw_snp_names,  
  chromosome = map_file$chromosome,
  position = map_file$position
)

# Merge candidate list with genomic coordinates
candidate_snps_to_join <- candidate_snps %>%
  rename(snp_id_r = snp_id)
annotated_candidates <- left_join(candidate_snps_to_join, snp_lookup_table, by = "snp_id_r")
write_csv(annotated_candidates, "../Data/GEA/annotated_high_confidence_candidates.csv", na = "NA")

# 4.2. PLOT 2: LFMM MANHATTAN PLOTS

########################
# LFMM Manhattan Plots # 
########################
fdr_threshold <- 0.05
fdr_threshold <- 0.02
fdr_threshold <- 0.01
chr_colors <- c("grey39", "grey83")
sig_colors <- c("#3030D0", "#9090F0")
par(mfrow = c(3, 1),
    mar   = c(1.5, 4, 1.5, 3),
    oma   = c(0, 0, 0, 0),
    mgp   = c(3, 0.3, 0))
variables_to_plot <- names(lfmm_results)
cex = 0.5
# Combine all p-values into one large data frame
all_pvals_df <- do.call(rbind, lapply(names(lfmm_results), function(var_name) {
  data.frame(snp_id_r = names(lfmm_results[[var_name]]), 
             p_value = lfmm_results[[var_name]], 
             variable = var_name)
}))
# Merge with SNP coordinate information
plot_data_full <- inner_join(all_pvals_df, snp_lookup_table, by = "snp_id_r")
plot_data_full <- plot_data_full %>% filter(!is.na(chromosome) & !is.na(position))
# Identify and filter for the 11 largest chromosomes
top_chromosomes <- plot_data_full %>% count(chromosome, sort = TRUE) %>% slice(1:11) %>% pull(chromosome)
plot_data_full <- plot_data_full %>% filter(chromosome %in% top_chromosomes)
plot_data_full$chromosome <- factor(plot_data_full$chromosome, levels = top_chromosomes)
plot_data_prepared <- plot_data_full %>%
  group_by(chromosome) %>%
  summarise(chr_len = max(position)) %>%
  mutate(total_len = cumsum(as.numeric(chr_len)) - chr_len) %>%
  select(-chr_len) %>%
  left_join(plot_data_full, ., by = "chromosome") %>%
  arrange(chromosome, position) %>%
  mutate(bp_cumulative = position + total_len)

# Calculate axis positions (chromosome midpoints and boundaries) for the shared x-axis
axis_df <- plot_data_prepared %>% 
  group_by(chromosome) %>% 
  summarize(center = mean(bp_cumulative), start = min(bp_cumulative), end = max(bp_cumulative))

for (i in 1:length(variables_to_plot)) {
  
  env_var_name <- variables_to_plot[i]
  current_plot_data <- filter(plot_data_prepared, variable == env_var_name)
  
  # Calculate -log10(p) and identify significant SNPs for THIS variable
  current_plot_data$log10_p <- -log10(current_plot_data$p_value)
  significant_snp_names <- lfmm_candidates[[env_var_name]] 
  current_plot_data$is_significant <- current_plot_data$snp_id_r %in% significant_snp_names
  if (any(current_plot_data$is_significant)) {
    threshold_log10 <- -log10(max(current_plot_data$p_value[current_plot_data$is_significant], na.rm = TRUE))
  } else { threshold_log10 <- NA }
  plot(x = current_plot_data$bp_cumulative, y = current_plot_data$log10_p,
       type = "n", xaxt = "n", yaxt = "s",
       ylab = "-log10(p-value)", xlab = "",
       main = paste(env_var_name),
       cex.lab = 1, cex.main = 1)
  axis(1, at = axis_df$center, labels = c(1:11), tick = F, cex.axis = 0.8)
  for (j in 1:length(top_chromosomes)) {
    bg_subset <- filter(current_plot_data, chromosome == top_chromosomes[j])
    points(x = bg_subset$bp_cumulative, y = bg_subset$log10_p, col = chr_colors[j %% 2 + 1], pch = 16, cex = cex)
    sig_subset <- filter(current_plot_data, chromosome == top_chromosomes[j] & is_significant == TRUE)
    if (nrow(sig_subset) > 0) {
      points(x = sig_subset$bp_cumulative, y = sig_subset$log10_p, col = sig_colors[j %% 2 + 1], pch = 16, cex = cex)
    }
  }
  if (!is.na(threshold_log10)) {
    abline(h = threshold_log10, col = "red", lty = 2, lwd = 1.5)
  }
}



##############################
# RDA & LFMM Manhattan Plots #
##############################
# Define the color scheme: Muted grays for background, one bright color for candidates
background_colors <- c("#BDBDBD", "#E0E0E0") # Two shades of gray
candidate_color <- rgb(1, 0, 0, 0.3)                # A strong red for high-confidence candidates
cex = 0.8
# Define which environmental variables to plot
variables_to_plot <- names(lfmm_results) # Or specify a subset, e.g., c("PD", "PW", "PC")

# Re-create the full candidate lists if they are not in the environment
all_lfmm_candidates <- unique(unlist(lfmm_candidates))
concordant_candidates <- intersect(rda_candidates, all_lfmm_candidates)
par(mfrow = c(3, 1),
    mar   = c(2, 5, 2, 2),
    oma   = c(3, 0, 2, 0),
    mgp   = c(3, 0.5, 0))

all_pvals_df <- do.call(rbind, lapply(names(lfmm_results), function(var_name) {
  data.frame(snp_id_r = names(lfmm_results[[var_name]]), p_value = lfmm_results[[var_name]], variable = var_name)
}))
plot_data_full <- inner_join(all_pvals_df, snp_lookup_table, by = "snp_id_r")
plot_data_full <- plot_data_full %>% filter(!is.na(chromosome) & !is.na(position))
top_chromosomes <- plot_data_full %>% count(chromosome, sort = TRUE) %>% slice(1:11) %>% pull(chromosome)
plot_data_full <- plot_data_full %>% filter(chromosome %in% top_chromosomes)
plot_data_full$chromosome <- factor(plot_data_full$chromosome, levels = top_chromosomes)
plot_data_prepared <- plot_data_full %>%
  group_by(chromosome) %>% summarise(chr_len = max(position)) %>%
  mutate(total_len = cumsum(as.numeric(chr_len)) - chr_len) %>%
  select(-chr_len) %>%
  left_join(plot_data_full, ., by = "chromosome") %>%
  arrange(chromosome, position) %>%
  mutate(bp_cumulative = position + total_len)
axis_df <- plot_data_prepared %>% group_by(chromosome) %>% 
  summarize(center = mean(bp_cumulative), start = min(bp_cumulative), end = max(bp_cumulative))

for (i in 1:length(variables_to_plot)) {
  
  env_var_name <- variables_to_plot[i]
  current_plot_data <- filter(plot_data_prepared, variable == env_var_name)
  current_plot_data$log10_p <- -log10(current_plot_data$p_value)
  
  # *** START OF REVISED LOGIC ***
  
  # Condition 1: Which SNPs are LFMM-significant for THIS specific variable?
  lfmm_sig_for_this_var <- lfmm_candidates[[env_var_name]]
  
  # Condition 2: From that list, which are ALSO on the global concordant list?
  # This is the list of SNPs we will highlight in THIS panel.
  highlight_candidates_for_this_panel <- intersect(concordant_candidates, lfmm_sig_for_this_var)
  
  # Create a flag in our data frame based on this strict list.
  current_plot_data$is_highlighted <- current_plot_data$snp_id_r %in% highlight_candidates_for_this_panel
  
  # *** END OF REVISED LOGIC ***
  
  # The FDR threshold line is still based on all LFMM hits for this variable.
  if (!is.null(lfmm_sig_for_this_var)) {
    threshold_log10 <- -log10(max(current_plot_data$p_value[current_plot_data$snp_id_r %in% lfmm_sig_for_this_var], na.rm = TRUE))
  } else { 
    threshold_log10 <- NA 
  }
  
  # Initialize an empty plot
  plot(x = current_plot_data$bp_cumulative, y = current_plot_data$log10_p,
       type = "n", xaxt = "n", yaxt = "s", ylab = "-log10(p-value)", xlab = "",
       main = paste(env_var_name))
  axis(1, at = axis_df$center, labels = c(1:11), tick = F, cex.axis = 0.8)
  # Layer 1: Plot ALL points as the gray background
  for (j in 1:length(top_chromosomes)) {
    chr_subset <- filter(current_plot_data, chromosome == top_chromosomes[j])
    points(x = chr_subset$bp_cumulative, y = chr_subset$log10_p, 
           col = background_colors[j %% 2 + 1], pch = 16, cex = cex)
  }
  # Layer 2: Plot the strictly filtered highlighted points
  highlight_points <- filter(current_plot_data, is_highlighted == TRUE)
  if (nrow(highlight_points) > 0) {
    points(
      x = highlight_points$bp_cumulative,
      y = highlight_points$log10_p,
      col = 'white',
      pch = 16, 
      cex = cex
    )
    points(
      x = highlight_points$bp_cumulative,
      y = highlight_points$log10_p,
      col = candidate_color,
      pch = 16, 
      cex = cex
    )
  }
  
  # Layer 3: Add the FDR threshold line
  if (!is.na(threshold_log10)) {
    abline(h = threshold_log10, col = "black", lty = 2, lwd = 1)
  }
}

# chromosome 8 
# Initialize an empty plot
current_plot_data <- current_plot_data[current_plot_data$chromosome == "JAROYF020000008.1",]
plot(x = current_plot_data$bp_cumulative, y = current_plot_data$log10_p,
     type = "n", xaxt = "n", yaxt = "s", ylab = "-log10(p-value)", xlab = "",
     main = paste(env_var_name),
     xlim = c(66847,398579)+447883963)
axis(1, at = axis_df$center, labels = c(1:11), tick = F, cex.axis = 0.8)
# Layer 1: Plot ALL points as the gray background
for (j in 1:length(top_chromosomes)) {
  chr_subset <- filter(current_plot_data, chromosome == top_chromosomes[j])
  points(x = chr_subset$bp_cumulative, y = chr_subset$log10_p, 
         col = background_colors[j %% 2 + 1], pch = 16, cex = cex)
}
# Layer 2: Plot the strictly filtered highlighted points
highlight_points <- filter(current_plot_data, is_highlighted == TRUE)
if (nrow(highlight_points) > 0) {
  points(
    x = highlight_points$bp_cumulative,
    y = highlight_points$log10_p,
    col = 'white',
    pch = 16, 
    cex = cex
  )
  points(
    x = highlight_points$bp_cumulative,
    y = highlight_points$log10_p,
    col = candidate_color,
    pch = 16, 
    cex = cex
  )
}

# Layer 3: Add the FDR threshold line
if (!is.na(threshold_log10)) {
  abline(h = threshold_log10, col = "black", lty = 2, lwd = 1)
}

abline(v = c(66847,69198,82562,85084,120819,123772)+447883963)
abline(v = c(152532,159598,346642,398579)+447883963)

 # plot