# ============================================================
# GENOTYPE-ENVIRONMENT ASSOCIATION ANALYSIS (GEA)
# RDA + LFMM Pipeline
# ============================================================
library(terra)
library(readr)
library(usdm)
library(dplyr)
library(vegan)
library(lfmm)
library(viridis)

# ============================================================
# PART 1: ENVIRONMENTAL DATA PREPARATION
# ============================================================

# Load GPS coordinates
gps_data <- read.csv("~/Desktop/Projects/Dynastes_grantii_pop/Data/GEA/coords.txt")
coords   <- gps_data[, c("lon", "lat")]

# Load WorldClim bioclimatic rasters
climate_files        <- list.files("~/Desktop/Projects/Dynastes_grantii_pop/Data/wc2.1_30s_bio/", pattern = "\\.tif$", full.names = TRUE)
climate_stack        <- rast(climate_files)
names(climate_stack) <- paste0("BIO", 1:19)

# Extract climate values at sample locations
extracted_values <- terra::extract(climate_stack, coords)[, -1]
env_data         <- cbind(gps_data, extracted_values)

# Save raw environmental data (run once, then comment out)
# write.csv(env_data, "../Data/GEA/environmental_data_raw.csv", row.names = FALSE)

# ── VIF filtering ──────────────────────────────────────────
env_data_raw      <- read.csv("~/Desktop/Projects/Dynastes_grantii_pop/Data/GEA/environmental_data_raw.csv")
selected_bio_vars <- env_data_raw[, -c(1:4)]          # drop non-BIO columns

vif_results      <- vifstep(selected_bio_vars, th = 10)
print(vif_results)
final_vars_names <- vif_results@results$Variables

# Build final environmental data frame (rename VIF-selected vars to PD, PW, PC)
final_env_data           <- env_data_raw[, c("SampleID", "lat", "lon", final_vars_names)]
colnames(final_env_data) <- c("SampleID", "lat", "lon", "PD", "PW", "PC")

# Save VIF-filtered data (run once, then comment out)
# write.csv(final_env_data, "../Data/GEA/environmental_data_final_VIF.csv", row.names = FALSE)

# ============================================================
# PART 2: GENOTYPE DATA PREPARATION
# ============================================================

# Load PLINK .raw file (6 info columns + SNP columns)
geno_raw              <- read.delim("~/Desktop/Projects/Dynastes_grantii_pop/Data/GEA/genotype_data.raw", sep = " ")
sample_ids_from_plink <- geno_raw$IID

# Build genotype matrix
genotype_matrix           <- as.matrix(geno_raw[, 7:ncol(geno_raw)])
rownames(genotype_matrix) <- sample_ids_from_plink
cat("Loaded:", nrow(genotype_matrix), "samples,", ncol(genotype_matrix), "SNPs.\n")

# Recode PLINK missing value (9) → NA, then impute with per-SNP mean
genotype_matrix[genotype_matrix == 9] <- NA

for (i in seq_len(ncol(genotype_matrix))) {
  col_mean <- mean(genotype_matrix[, i], na.rm = TRUE)
  genotype_matrix[is.na(genotype_matrix[, i]), i] <- col_mean
}

if (any(is.na(genotype_matrix))) {
  warning("NAs still present after imputation!")
} else {
  cat("Missing data successfully imputed.\n")
}

# ============================================================
# PART 3: SYNCHRONISE GENOTYPE & ENVIRONMENTAL DATA
# ============================================================

env_data_final <- read.csv("~/Desktop/Projects/Dynastes_grantii_pop/Data/GEA/environmental_data_final_VIF.csv")
sample_order   <- env_data_final$SampleID

if (!all(sample_order %in% rownames(genotype_matrix))) {
  stop("ERROR: Not all samples from the environmental file are in the genotype file!")
}

genotype_matrix <- genotype_matrix[sample_order, ]
env_predictors  <- env_data_final %>% select(PD, PW, PC)

cat("Genotype matrix:  ", nrow(genotype_matrix), "samples x", ncol(genotype_matrix), "SNPs\n")
cat("Env. predictors:  ", nrow(env_predictors),  "samples x", ncol(env_predictors),  "vars\n")

# ============================================================
# PART 4: REDUNDANCY ANALYSIS (RDA)
# ============================================================

rda_model    <- rda(genotype_matrix ~ ., data = env_predictors, scale = TRUE)
anova_result <- anova.cca(rda_model, parallel = getOption("mc.cores"))
print(anova_result)

r_squared <- RsquareAdj(rda_model)
cat("\nAdjusted R² for RDA model:", r_squared$adj.r.squared, "\n")
print(summary(rda_model)$concont)

# ── Identify outlier SNPs (> 3 SD from mean loading on each axis) ──
get_rda_outliers <- function(rda_object, axis_number, num_sds = 3) {
  loadings   <- scores(rda_object, choices = axis_number, display = "species")
  upper      <- mean(loadings) + num_sds * sd(loadings)
  lower      <- mean(loadings) - num_sds * sd(loadings)
  idx        <- which(loadings > upper | loadings < lower)
  rownames(loadings)[idx]
}

rda_candidates <- unique(c(
  get_rda_outliers(rda_model, 1),
  get_rda_outliers(rda_model, 2)
))
cat("RDA candidate SNPs (> 3 SD):", length(rda_candidates), "\n")
# 1411
# ============================================================
# PART 5: LATENT FACTOR MIXED MODELS (LFMM)
# ============================================================

K             <- 4      # number of latent factors (estimated from PCA scree plot)
lfmm_results  <- list()

for (env_var_name in colnames(env_predictors)) {
  cat("Running LFMM for:", env_var_name, "...\n")
  
  X          <- as.matrix(env_predictors[[env_var_name]])
  lfmm_fit   <- lfmm_ridge(Y = genotype_matrix, X = X, K = K)
  pv         <- lfmm_test(Y = genotype_matrix, X = X, lfmm = lfmm_fit, calibrate = "gif")
  lfmm_results[[env_var_name]] <- pv$pvalue[, 1]
}

# ── FDR correction ────────────────────────────────────────
fdr_threshold   <- 0.01
lfmm_candidates <- list()

for (env_var_name in names(lfmm_results)) {
  adj_p   <- p.adjust(lfmm_results[[env_var_name]], method = "fdr")
  sig_idx <- which(adj_p < fdr_threshold)
  if (length(sig_idx) > 0) {
    lfmm_candidates[[env_var_name]] <- colnames(genotype_matrix)[sig_idx]
  }
}

cat("\n--- LFMM Results Summary ---\n")
for (env_var_name in names(lfmm_candidates)) {
  cat("  ", env_var_name, ":", length(lfmm_candidates[[env_var_name]]),
      "candidates at FDR <", fdr_threshold, "\n")
}
# PD : 950 candidates at FDR < 0.01 
# PW : 8462 candidates at FDR < 0.01 
# PC : 169 candidates at FDR < 0.01 

# ============================================================
# PART 6: SYNTHESISE RDA + LFMM RESULTS
# ============================================================

all_lfmm_candidates   <- unique(unlist(lfmm_candidates))
concordant_candidates <- intersect(rda_candidates, all_lfmm_candidates)

cat("\n--- Analysis Synthesis ---\n")
cat("  RDA candidates:         ", length(rda_candidates),          "\n")
cat("  LFMM candidates:        ", length(all_lfmm_candidates),     "\n")
cat("  Concordant candidates:  ", length(concordant_candidates),   "\n")
# --- Analysis Synthesis ---
# RDA candidates:          1411  
# LFMM candidates:         9319 
# Concordant candidates:   1275 

# ── Build annotated candidate table ───────────────────────
lfmm_association_vec <- sapply(concordant_candidates, function(snp) {
  vars <- names(lfmm_candidates)[sapply(lfmm_candidates, function(x) snp %in% x)]
  paste(vars, collapse = ";")
})

final_candidates_df <- data.frame(
  snp_id           = concordant_candidates,
  lfmm_association = lfmm_association_vec,
  stringsAsFactors = FALSE
)

# ── Build colour lookup for plotting ──────────────────────
env_var_names <- colnames(env_predictors)

color_palette <- setNames(c("#E69F00", "#56B4E9", "#009E73"), env_var_names)
color_palette["Multiple"] <- "#000000"

candidates_for_plot <- final_candidates_df %>%
  mutate(
    main_predictor = sapply(strsplit(lfmm_association, ";"), function(x) {
      if (length(x) == 1) x else "Multiple"
    }),
    color = color_palette[main_predictor]
  )

# ── SNP colour vector for RDA biplots ─────────────────────
all_rda_snps <- rownames(scores(rda_model, display = "species"))
snp_colors   <- rep("gray80", length(all_rda_snps))

for (i in seq_len(nrow(candidates_for_plot))) {
  idx <- match(candidates_for_plot$snp_id[i], all_rda_snps)
  if (!is.na(idx)) snp_colors[idx] <- candidates_for_plot$color[i]
}

# ── Add genomic coordinates and save ──────────────────────
map_file        <- read_tsv("~/Desktop/Projects/Dynastes_grantii_pop/Data/GEA/genotype_data.map", col_names = FALSE)
colnames(map_file) <- c("chromosome", "snp_id_plink", "genetic_dist", "position")

raw_snp_names <- colnames(geno_raw)[7:ncol(geno_raw)]
if (length(raw_snp_names) != nrow(map_file)) {
  stop("SNP count mismatch between .raw header and .map file!")
}

snp_lookup_table <- data.frame(
  snp_id_r   = raw_snp_names,
  chromosome = map_file$chromosome,
  position   = map_file$position
)

annotated_candidates <- final_candidates_df %>%
  rename(snp_id_r = snp_id) %>%
  left_join(snp_lookup_table, by = "snp_id_r")

write_csv(annotated_candidates, "~/Desktop/Projects/Dynastes_grantii_pop/Data/GEA/annotated_high_confidence_candidates.csv", na = "NA")

saveRDS(genotype_matrix,        "~/Desktop/Projects/Dynastes_grantii_pop/Data/GEA/genotype_matrix.rds")
saveRDS(concordant_candidates,  "~/Desktop/Projects/Dynastes_grantii_pop/Data/GEA/concordant_candidates.rds")
saveRDS(rda_candidates,         "~/Desktop/Projects/Dynastes_grantii_pop/Data/GEA/rda_candidates.rds")
saveRDS(rda_model,              "~/Desktop/Projects/Dynastes_grantii_pop/Data/GEA/rda_model.rds")
saveRDS(lfmm_candidates,        "~/Desktop/Projects/Dynastes_grantii_pop/Data/GEA/lfmm_candidates.rds")
saveRDS(env_predictors,         "~/Desktop/Projects/Dynastes_grantii_pop/Data/GEA/env_predictors.rds")
saveRDS(final_env_data,         "~/Desktop/Projects/Dynastes_grantii_pop/Data/GEA/final_env_data.rds")

cat("GEA objects saved for genetic offset analysis.\n")

# ============================================================
# PART 7: HELPER FUNCTIONS FOR PLOTS
# ============================================================

add_pop_legend <- function(x, y, s, cex_pt = 0.8, cex_txt = 0.8, pch_val = 16) {
  cols <- viridis(5)
  labels <- c("Mogollon Rim (Payson)", "Mogollon Rim (Reserve)", "Mt. Lemmon", "Chiricahua Mt.", "Utah")
  for (k in 1:5) {
    points(x, y - (k - 1) * s, col = cols[k], pch = pch_val, cex = cex_pt)
    text(x,   y - (k - 1) * s, pos = 4, labels = labels[k], cex = cex_txt)
  }
}

add_predictor_legend <- function(x, y, s, cex_pt = 0.8, cex_txt = 0.8) {
  cols   <- c("#E69F00", "#56B4E9", "#009E73", "#000000")
  labels <- c("PD", "PW", "PC", "Mixed")
  for (k in 1:4) {
    points(x, y - (k - 1) * s, bg = cols[k], col = "black", pch = 21, cex = cex_pt)
    text(x,   y - (k - 1) * s, pos = 4, labels = labels[k], cex = cex_txt)
  }
}

# ============================================================
# PART 8: RDA PLOTS
# ============================================================

screeplot(rda_model)

# Population colour assignments
cols_pop <- viridis(5)
bg_pop   <- c(rep(cols_pop[3], 5),   # Mt. Lemmon
              rep(cols_pop[4], 5),   # Portal (Chiricahua)
              rep(cols_pop[5], 5),   # Utah
              rep(cols_pop[1], 6),
              rep(cols_pop[2], 6)
)  # Mogollon Rim

# ── 2 × 2 panel figure ────────────────────────────────────
par(mfrow = c(2, 2), 
    mar = c(3, 3, 1.5, 1),
    mgp   = c(2, 0.5, 0))
# Panel A – sites, axes 1 & 2
plot(rda_model, scaling = 3, xlim = c(-7, 14), ylim = c(-10, 7),
     main = "", xlab = "RDA1 (46.8%)", ylab = "RDA2 (30.4%)")
points(rda_model, display = "species", pch = 21, cex = 1, col = "gray32", scaling = 3)
points(rda_model, display = "sites",   pch = 21, cex = 1, col = "black",  bg = bg_pop, scaling = 3)
text(rda_model,   scaling = 3, display = "bp", col = "#0868ac", cex = 0.7)
title(main = "(A)", adj = 0, line = 0.5)
add_pop_legend(x = 2, y = -6.8, s = 1.1, cex_pt = 0.8, cex_txt = 0.8, pch_val = 16)

# Panel B – SNPs, axes 1 & 2
plot(rda_model, type = "n", scaling = 3,
     main = "", xlab = "RDA1 (46.8%)", ylab = "RDA2 (30.4%)",
     xlim = c(-0.25, 0.25), ylim = c(-0.25, 0.25))
points(rda_model, display = "species", pch = 21, bg = snp_colors, col = "gray32", cex = 1, scaling = 3)
points(rda_model, display = "sites",   pch = 4,  col = "gray30",  cex = 0.8,      scaling = 3)
text(rda_model,   display = "bp", col = "#0868ac", cex = 0.7, scaling = 3, font = 1)
title(main = "(B)", adj = 0, line = 0.5)
add_predictor_legend(x = -0.26, y = 0.24, s = 0.03)

# Panel C – sites, axes 1 & 3
plot(rda_model, scaling = 3, xlim = c(-6.2, 14), ylim = c(-8, 9),
     main = "", xlab = "RDA1 (46.8%)", ylab = "RDA3 (22.8%)", choices = c(1, 3))
points(rda_model, display = "species", pch = 21, cex = 1, col = "gray32", scaling = 3, choices = c(1, 3))
points(rda_model, display = "sites",   pch = 21, cex = 1, col = "black",  bg = bg_pop, scaling = 3, choices = c(1, 3))
text(rda_model,   scaling = 3, display = "bp", col = "#0868ac", cex = 0.7, choices = c(1, 3))
title(main = "(C)", adj = 0, line = 0.5)
add_pop_legend(x = 2.5, y = -4.8, s = 1, cex_pt = 0.8, cex_txt = 0.8, pch_val = 16)

# Panel D – SNPs, axes 1 & 3
plot(rda_model, type = "n", scaling = 3,
     main = "", xlab = "RDA1 (46.8%)", ylab = "RDA3 (22.8%)",
     xlim = c(-0.35, 0.35), ylim = c(-0.35, 0.35), choices = c(1, 3))
points(rda_model, display = "species", pch = 21, bg = snp_colors, col = "gray32", cex = 1, scaling = 3, choices = c(1, 3))
text(rda_model,   display = "bp", col = "#0868ac", cex = 0.7, scaling = 3, font = 1, choices = c(1, 3))
title(main = "(D)", adj = 0, line = 0.5)
add_predictor_legend(x = -0.35, y = 0.34, s = 0.04)

# ============================================================
# PART 9: LFMM MANHATTAN PLOTS
# ============================================================
# Color scheme
background_colors <- c("#BDBDBD", "#E0E0E0")  # alternating grey for chromosomes
candidate_color   <- rgb(1, 0, 0, 0.3)         # translucent red for concordant candidates
cex               <- 0.8

# Variables to plot
variables_to_plot     <- names(lfmm_results)
all_lfmm_candidates   <- unique(unlist(lfmm_candidates))
concordant_candidates <- intersect(rda_candidates, all_lfmm_candidates)

# ── Build cumulative-position data frame ──────────────────
all_pvals_df <- do.call(rbind, lapply(names(lfmm_results), function(var_name) {
  data.frame(
    snp_id_r = names(lfmm_results[[var_name]]),
    p_value  = lfmm_results[[var_name]],
    variable = var_name
  )
}))

plot_data_full <- inner_join(all_pvals_df, snp_lookup_table, by = "snp_id_r") %>%
  filter(!is.na(chromosome), !is.na(position))

top_chromosomes        <- plot_data_full %>% count(chromosome, sort = TRUE) %>% slice(1:11) %>% pull(chromosome)
plot_data_full         <- plot_data_full %>% filter(chromosome %in% top_chromosomes)
plot_data_full$chromosome <- factor(plot_data_full$chromosome, levels = top_chromosomes)

plot_data_prepared <- plot_data_full %>%
  group_by(chromosome) %>%
  summarise(chr_len = max(position)) %>%
  mutate(total_len = cumsum(as.numeric(chr_len)) - chr_len) %>%
  select(-chr_len) %>%
  left_join(plot_data_full, ., by = "chromosome") %>%
  arrange(chromosome, position) %>%
  mutate(bp_cumulative = position + total_len)

axis_df <- plot_data_prepared %>%
  group_by(chromosome) %>%
  summarize(center = mean(bp_cumulative),
            start  = min(bp_cumulative),
            end    = max(bp_cumulative))

# ── Three-panel Manhattan plot ─────────────────────────────
par(mfrow = c(3, 1),
    mar   = c(1, 3, 1.5, 1),
    oma   = c(3, 0, 2, 0),
    mgp   = c(2, 0.5, 0))
titles <- c('E','F','G')
for (i in 1:length(variables_to_plot)) {
  
  env_var_name      <- variables_to_plot[i]
  current_plot_data <- filter(plot_data_prepared, variable == env_var_name)
  current_plot_data$log10_p <- -log10(current_plot_data$p_value)
  
  # SNPs significant in LFMM for this variable AND concordant with RDA
  lfmm_sig_for_this_var           <- lfmm_candidates[[env_var_name]]
  highlight_candidates_for_this_panel <- intersect(concordant_candidates, lfmm_sig_for_this_var)
  current_plot_data$is_highlighted <- current_plot_data$snp_id_r %in% highlight_candidates_for_this_panel
  
  # FDR threshold line (based on all LFMM hits for this variable)
  if (!is.null(lfmm_sig_for_this_var)) {
    threshold_log10 <- -log10(max(current_plot_data$p_value[current_plot_data$snp_id_r %in% lfmm_sig_for_this_var], na.rm = TRUE))
  } else {
    threshold_log10 <- NA
  }
  
  # Empty canvas
  plot(x = current_plot_data$bp_cumulative, y = current_plot_data$log10_p,
       type = "n", xaxt = "n", yaxt = "s",
       ylab = "-log10(p-value)", xlab = "",
       main = '')
  title(main = paste0('(',titles[i],')'), adj = 0, line = 0.5)
  axis(1, at = axis_df$center, labels = 1:11, tick = FALSE, cex.axis = 1)
  
  # Layer 1: all SNPs as grey background
  for (j in 1:length(top_chromosomes)) {
    chr_subset <- filter(current_plot_data, chromosome == top_chromosomes[j])
    points(x = chr_subset$bp_cumulative, y = chr_subset$log10_p,
           col = background_colors[j %% 2 + 1], pch = 16, cex = cex)
  }
  
  # Layer 2: concordant candidates (white mask then red)
  highlight_points <- filter(current_plot_data, is_highlighted == TRUE)
  if (nrow(highlight_points) > 0) {
    points(x = highlight_points$bp_cumulative, y = highlight_points$log10_p,
           col = "white", pch = 16, cex = cex)
    points(x = highlight_points$bp_cumulative, y = highlight_points$log10_p,
           col = candidate_color, pch = 16, cex = cex)
  }
  
  # Layer 3: FDR threshold line
  if (!is.na(threshold_log10)) {
    abline(h = threshold_log10, col = "black", lty = 2, lwd = 1)
  }
}

mtext("Chromosome", side = 1, outer = TRUE, line = 1.5, cex = 1)
# order is variables_to_plot
# "PD" "PW" "PC"

# ── Chromosome 8 zoom (last variable in loop) ─────────────
chr8_id     <- "JAROYF020000008.1"
chr8_offset <- 447883963
xlim_chr8   <- c(66847, 398579) + chr8_offset

current_plot_data <- current_plot_data[current_plot_data$chromosome == chr8_id, ]

# Empty canvas
plot(x = current_plot_data$bp_cumulative, y = current_plot_data$log10_p,
     type = "n", xaxt = "n", yaxt = "s",
     ylab = "-log10(p-value)", xlab = "",
     main = env_var_name,
     xlim = xlim_chr8)
axis(1, at = axis_df$center, labels = 1:11, tick = FALSE, cex.axis = 0.8)

# Layer 1: all SNPs as grey background
for (j in 1:length(top_chromosomes)) {
  chr_subset <- filter(current_plot_data, chromosome == top_chromosomes[j])
  points(x = chr_subset$bp_cumulative, y = chr_subset$log10_p,
         col = background_colors[j %% 2 + 1], pch = 16, cex = cex)
}

# Layer 2: concordant candidates (white mask then red)
highlight_points <- filter(current_plot_data, is_highlighted == TRUE)
if (nrow(highlight_points) > 0) {
  points(x = highlight_points$bp_cumulative, y = highlight_points$log10_p,
         col = "white", pch = 16, cex = cex)
  points(x = highlight_points$bp_cumulative, y = highlight_points$log10_p,
         col = candidate_color, pch = 16, cex = cex)
}

# Layer 3: FDR threshold line
if (!is.na(threshold_log10)) {
  abline(h = threshold_log10, col = "black", lty = 2, lwd = 1)
}

# Gene structure annotations
abline(v = c(66847, 69198,  82562,  85084, 120819, 123772) + chr8_offset)  # exon boundaries
abline(v = c(152532, 159598, 346642, 398579) + chr8_offset)                # intron/region boundaries
