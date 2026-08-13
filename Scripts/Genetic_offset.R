# ============================================================
# GENETIC OFFSET ANALYSIS
# Dynastes Sky Island Project
# Builds on: RDA + LFMM GEA pipeline
# Method: RDA-based genetic offset (Capblancq & Forester 2021)
# ============================================================

library(terra)
library(dplyr)
library(vegan)
library(ggplot2)
library(viridis)
library(tidyr)
library(readr)
library(usdm)
library(maps)

# ============================================================
# PART 1: LOAD GEA OUTPUTS
# ============================================================

genotype_matrix       <- readRDS("~/Desktop/Projects/Dynastes_grantii_pop/Data/GEA/genotype_matrix.rds")
concordant_candidates <- readRDS("~/Desktop/Projects/Dynastes_grantii_pop/Data/GEA/concordant_candidates.rds")
rda_candidates        <- readRDS("~/Desktop/Projects/Dynastes_grantii_pop/Data/GEA/rda_candidates.rds")
lfmm_candidates       <- readRDS("~/Desktop/Projects/Dynastes_grantii_pop/Data/GEA/lfmm_candidates.rds")
env_predictors        <- readRDS("~/Desktop/Projects/Dynastes_grantii_pop/Data/GEA/env_predictors.rds")
final_env_data        <- readRDS("~/Desktop/Projects/Dynastes_grantii_pop/Data/GEA/final_env_data.rds")
rda_model             <- readRDS("~/Desktop/Projects/Dynastes_grantii_pop/Data/GEA/rda_model.rds")
gps_data              <- read.csv("~/Desktop/Projects/Dynastes_grantii_pop/Data/GEA/coords.txt")

cat("Concordant candidate SNPs:", length(concordant_candidates), "\n")
cat("Samples:                  ", nrow(genotype_matrix), "\n")
pop_assignments <- c(
  rep("Mt. Lemmon",             5),
  rep("Chiricahua Mt.",         5),
  rep("Utah",                   5),
  rep("Mogollon Rim ",  12)
)
stopifnot(length(pop_assignments) == nrow(final_env_data))

# ============================================================
# PART 2: IDENTIFY BIO VARIABLES USED IN GEA (PD / PW / PC)
# ============================================================
# Re-run VIF to confirm which BIO variables became PD, PW, PC
# Position 1 = PD, 2 = PW, 3 = PC

env_data_raw   <- read.csv("~/Desktop/Projects/Dynastes_grantii_pop/Data/GEA/environmental_data_raw.csv")
selected_bio   <- env_data_raw[, -c(1:4)]          # drop non-BIO columns
vif_results    <- vifstep(selected_bio, th = 10)
final_vars_bio <- vif_results@results$Variables     # e.g. c("BIO4","BIO12","BIO15")

cat("BIO variables (= PD, PW, PC):", final_vars_bio, "\n")

# ============================================================
# PART 3: LOAD AND PREPARE CLIMATE RASTERS
# ============================================================

# ── Current climate ───────────────────────────────────────
climate_files <- list.files("~/Desktop/Projects/Dynastes_grantii_pop/Data/wc2.1_30s_bio/",
                            pattern = "\\.tif$", full.names = TRUE)
climate_stack <- rast(climate_files)
names(climate_stack) <- paste0("BIO", 1:19)

# ── Future climate (WorldClim 2.1 / BCC-CSM2-MR) ─────────
future_2041_ssp245 <- rast("~/Downloads/wc2.1_30s_bioc_BCC-CSM2-MR_ssp245_2041-2060.tif")
future_2041_ssp585 <- rast("~/Downloads/wc2.1_30s_bioc_BCC-CSM2-MR_ssp585_2041-2060.tif")
future_2061_ssp245 <- rast("~/Downloads/wc2.1_30s_bioc_BCC-CSM2-MR_ssp245_2061-2080.tif")
future_2061_ssp585 <- rast("~/Downloads/wc2.1_30s_bioc_BCC-CSM2-MR_ssp585_2061-2080.tif")
future_2081_ssp245 <- rast("~/Downloads/wc2.1_30s_bioc_BCC-CSM2-MR_ssp245_2081-2100.tif")
future_2081_ssp585 <- rast("~/Downloads/wc2.1_30s_bioc_BCC-CSM2-MR_ssp585_2081-2100.tif")

names(future_2041_ssp245) <- paste0("BIO", 1:19)
names(future_2041_ssp585) <- paste0("BIO", 1:19)
names(future_2061_ssp245) <- paste0("BIO", 1:19)
names(future_2061_ssp585) <- paste0("BIO", 1:19)
names(future_2081_ssp245) <- paste0("BIO", 1:19)
names(future_2081_ssp585) <- paste0("BIO", 1:19)

# ── Helper: extract the 3 target BIO layers → rename PD/PW/PC ──
prep_climate_stack <- function(stack, bio_vars) {
  s <- stack[[bio_vars]]
  names(s) <- c("PD", "PW", "PC")
  return(s)
}

current_3v      <- prep_climate_stack(climate_stack,      final_vars_bio)
fut_245_2041_3v <- prep_climate_stack(future_2041_ssp245, final_vars_bio)
fut_585_2041_3v <- prep_climate_stack(future_2041_ssp585, final_vars_bio)
fut_245_2061_3v <- prep_climate_stack(future_2061_ssp245, final_vars_bio)
fut_585_2061_3v <- prep_climate_stack(future_2061_ssp585, final_vars_bio)
fut_245_2081_3v <- prep_climate_stack(future_2081_ssp245, final_vars_bio)
fut_585_2081_3v <- prep_climate_stack(future_2081_ssp585, final_vars_bio)

# ── Crop all rasters to sky islands study extent ──────────
coords_vect  <- vect(gps_data, geom = c("lon", "lat"), crs = "EPSG:4326")
study_extent <- ext(coords_vect) + 2.5
current_3v      <- crop(current_3v,      study_extent)
fut_245_2041_3v <- crop(fut_245_2041_3v, study_extent)
fut_585_2041_3v <- crop(fut_585_2041_3v, study_extent)
fut_245_2061_3v <- crop(fut_245_2061_3v, study_extent)
fut_585_2061_3v <- crop(fut_585_2061_3v, study_extent)
fut_245_2081_3v <- crop(fut_245_2081_3v, study_extent)
fut_585_2081_3v <- crop(fut_585_2081_3v, study_extent)

# Sanity check — sample points should sit on land
plot(current_3v[["PD"]], main = "PD (Current) — verify sample locations")
map('state', fill = FALSE, add = TRUE, col = "black")
map.scale(x=-114.5, y=30.5, ratio=FALSE, relwidth=0.2, cex = 0.5)
points(coords_vect, pch = 16, col = "red", cex = 1)

# ============================================================
# PART 4: RDA-BASED GENETIC OFFSET FUNCTION
# ============================================================
# Logic:
#   predict(rda, newdata = current_climate) → current genetic position
#   predict(rda, newdata = future_climate)  → future  genetic position
#   offset = Euclidean distance between the two in RDA space

calc_rda_offset_raster <- function(rda_object,
                                   current_rast,
                                   future_rast,
                                   n_axes = 3) {
  
  # Rasters → data frames (keep xy for reconstruction)
  curr_df  <- as.data.frame(current_rast, xy = TRUE, na.rm = FALSE)
  fut_df   <- as.data.frame(future_rast,  xy = TRUE, na.rm = FALSE)
  
  xy       <- curr_df[, c("x", "y")]
  curr_env <- curr_df[, c("PD", "PW", "PC")]
  fut_env  <- fut_df[,  c("PD", "PW", "PC")]
  
  # Process only cells that are complete in BOTH layers
  complete <- complete.cases(curr_env) & complete.cases(fut_env)
  cat("  Processing", sum(complete), "non-NA cells...\n")
  
  # Project climate onto RDA axes (linear combination scores)
  scores_current <- predict(rda_object,
                            newdata = curr_env[complete, ],
                            type    = "lc",
                            scaling = 1)[, 1:n_axes]
  
  scores_future  <- predict(rda_object,
                            newdata = fut_env[complete, ],
                            type    = "lc",
                            scaling = 1)[, 1:n_axes]
  
  # Euclidean distance in RDA space = genetic offset
  offset_vals <- sqrt(rowSums((scores_future - scores_current)^2))
  
  # Reconstruct raster
  offset_rast         <- rast(current_rast[[1]])
  values(offset_rast) <- NA
  cell_ids            <- cellFromXY(offset_rast, xy[complete, ])
  offset_rast[cell_ids] <- offset_vals
  
  return(offset_rast)
}

# ============================================================
# PART 5: CALCULATE LANDSCAPE-LEVEL GENETIC OFFSET
# ============================================================

cat("Calculating landscape genetic offsets...\n")

cat("SSP2-4.5  2041-2060...\n")
offset_rast_245_2041 <- calc_rda_offset_raster(rda_model, current_3v, fut_245_2041_3v)

cat("SSP5-8.5  2041-2060...\n")
offset_rast_585_2041 <- calc_rda_offset_raster(rda_model, current_3v, fut_585_2041_3v)

cat("SSP2-4.5  2061-2080...\n")
offset_rast_245_2061 <- calc_rda_offset_raster(rda_model, current_3v, fut_245_2061_3v)

cat("SSP5-8.5  2061-2080...\n")
offset_rast_585_2061 <- calc_rda_offset_raster(rda_model, current_3v, fut_585_2061_3v)

cat("SSP2-4.5  2081-2100...\n")
offset_rast_245_2081 <- calc_rda_offset_raster(rda_model, current_3v, fut_245_2081_3v)

cat("SSP5-8.5  2081-2100...\n")
offset_rast_585_2081 <- calc_rda_offset_raster(rda_model, current_3v, fut_585_2081_3v)


# ── Save rasters ──────────────────────────────────────────
writeRaster(offset_rast_245_2041,
            "~/Desktop/Projects/Dynastes_grantii_pop/Data/GEA/GeneticOffset/rda_offset_ssp245_2041-2060.tif", overwrite = TRUE)
writeRaster(offset_rast_585_2041,
            "~/Desktop/Projects/Dynastes_grantii_pop/Data/GEA/GeneticOffset/rda_offset_ssp585_2041-2060.tif", overwrite = TRUE)
writeRaster(offset_rast_245_2061,
            "~/Desktop/Projects/Dynastes_grantii_pop/Data/GEA/GeneticOffset/rda_offset_ssp245_2061-2080.tif", overwrite = TRUE)
writeRaster(offset_rast_585_2061,
            "~/Desktop/Projects/Dynastes_grantii_pop/Data/GEA/GeneticOffset/rda_offset_ssp585_2061-2080.tif", overwrite = TRUE)
writeRaster(offset_rast_245_2081,
            "~/Desktop/Projects/Dynastes_grantii_pop/Data/GEA/GeneticOffset/rda_offset_ssp245_2081-2100.tif", overwrite = TRUE)
writeRaster(offset_rast_585_2081,
            "~/Desktop/Projects/Dynastes_grantii_pop/Data/GEA/GeneticOffset/rda_offset_ssp585_2081-2100.tif", overwrite = TRUE)

# ============================================================
# PART 6: POPULATION-LEVEL GENETIC OFFSET
# ============================================================

coords_mat <- as.matrix(gps_data[, c("lon", "lat")])

# Predict RDA scores at sample locations under current climate
scores_current_samples <- predict(rda_model,
                                  newdata = env_predictors,
                                  type    = "lc",
                                  scaling = 1)[, 1:3]

# Helper: extract future climate at samples → predict RDA scores → offset
calc_sample_offset <- function(future_rast, scores_current, rda_object,
                               coords_mat) {
  fut_env <- terra::extract(future_rast, coords_mat)[, c("PD", "PW", "PC")]
  scores_fut <- predict(rda_object,
                        newdata = fut_env,
                        type    = "lc",
                        scaling = 1)[, 1:3]
  sqrt(rowSums((scores_fut - scores_current)^2))
}

ind_offset <- data.frame(
  SampleID        = final_env_data$SampleID,
  Population      = pop_assignments,
  lon             = gps_data$lon,
  lat             = gps_data$lat,
  PD_current      = env_predictors$PD,
  PW_current      = env_predictors$PW,
  PC_current      = env_predictors$PC,
  offset_245_2041 = calc_sample_offset(fut_245_2041_3v, scores_current_samples,
                                       rda_model, coords_mat),
  offset_585_2041 = calc_sample_offset(fut_585_2041_3v, scores_current_samples,
                                       rda_model, coords_mat),
  offset_245_2061 = calc_sample_offset(fut_245_2061_3v, scores_current_samples,
                                       rda_model, coords_mat),
  offset_585_2061 = calc_sample_offset(fut_585_2061_3v, scores_current_samples,
                                       rda_model, coords_mat),
  offset_245_2081 = calc_sample_offset(fut_245_2081_3v, scores_current_samples,
                                       rda_model, coords_mat),
  offset_585_2081 = calc_sample_offset(fut_585_2081_3v, scores_current_samples,
                                       rda_model, coords_mat)

)

# ── Population summary ────────────────────────────────────
pop_offset_summary <- ind_offset %>%
  group_by(Population) %>%
  summarise(
    n                  = n(),
    mean_lon           = mean(lon),
    mean_lat           = mean(lat),
    off_245_2041_mean  = mean(offset_245_2041),
    off_245_2041_sd    = sd(offset_245_2041),
    off_585_2041_mean  = mean(offset_585_2041),
    off_585_2041_sd    = sd(offset_585_2041),
    off_245_2061_mean  = mean(offset_245_2061),
    off_245_2061_sd    = sd(offset_245_2061),
    off_585_2061_mean  = mean(offset_585_2061),
    off_585_2061_sd    = sd(offset_585_2061),
    off_245_2081_mean  = mean(offset_245_2081),
    off_245_2081_sd    = sd(offset_245_2081),
    off_585_2081_mean  = mean(offset_585_2081),
    off_585_2081_sd    = sd(offset_585_2081),
    .groups = "drop"
  ) %>%
  arrange(desc(off_585_2061_mean))   # rank: highest offset first

cat("\n--- Population Genetic Offset Ranking (worst-case top) ---\n")
print(pop_offset_summary)

write_csv(ind_offset,         "~/Desktop/Projects/Dynastes_grantii_pop/Data/GEA/GeneticOffset/individual_offsets.csv")
write_csv(pop_offset_summary, "~/Desktop/Projects/Dynastes_grantii_pop/Data/GEA/GeneticOffset/population_offset_summary.csv")

# ============================================================
# PART 7: VISUALISATION — LANDSCAPE MAPS
# ============================================================

# Population colours — must match your sample order
cols_pop <- viridis(4)
pop_colors_pts <- c(
  rep(cols_pop[3], 5),    # Mt. Lemmon
  rep(cols_pop[2], 5),    # Chiricahua
  rep(cols_pop[4], 5),    # Utah
  rep(cols_pop[1], 12)    # Mogollon Rim (Payson + Reserve)
)

# Shared z-scale across ALL 6 panels for fair comparison
all_vals <- c(
  values(offset_rast_245_2041), values(offset_rast_585_2041),
  values(offset_rast_245_2061), values(offset_rast_585_2061),
  values(offset_rast_245_2081), values(offset_rast_585_2081)
)
zlim_share <- c(0, quantile(all_vals, 0.99, na.rm = TRUE))

scenarios <- list(
  list(rast = offset_rast_245_2041, title = "(A) SSP2-4.5 | 2041-2060"),
  list(rast = offset_rast_585_2041, title = "(B) SSP5-8.5 | 2041-2060"),
  list(rast = offset_rast_245_2061, title = "(C) SSP2-4.5 | 2061-2080"),
  list(rast = offset_rast_585_2061, title = "(D) SSP5-8.5 | 2061-2080"),
  list(rast = offset_rast_245_2081, title = "(E) SSP2-4.5 | 2081-2100"),
  list(rast = offset_rast_585_2081, title = "(F) SSP5-8.5 | 2081-2100")
)
par(mfrow = c(3, 2),
    mar   = c(0, 0, 0, 0),
    oma   = c(0, 0, 0, 0))

for (sc in scenarios) {
  plot(sc$rast,
       col    = plasma(100),
       zlim   = zlim_share,
       main   = sc$title,
       axes   = FALSE,
       legend = TRUE)
  map('state', fill = FALSE, add = TRUE, col = "black")
  # map.scale(x=-109.8, y=30.6, ratio=FALSE, relwidth=0.2, cex = 1)
  map.scale(x=-109.8, y=31.2, ratio=FALSE, relwidth=0.1, cex = 0.35)
  points(coords_mat,
         pch = 21,
         bg  = pop_colors_pts,
         col = "black",
         cex = 1.2,
         lwd = 0.5)
}

# Main title
mtext("Genetic Offset — Dynastes Sky Islands",
      outer = TRUE, side = 3,
      line  = 0.5, cex = 1.2, font = 2)

# Shared x-axis label
mtext("Higher values = more evolutionary change required",
      outer = TRUE, side = 1,
      line  = 0.5, cex = 0.8, col = "grey40")
# ============================================================
# PART 8: VISUALISATION — POPULATION BAR CHART
# ============================================================

# ── Clean population names (remove trailing spaces) ────────
pop_offset_summary$Population <- trimws(pop_offset_summary$Population)

# ── Build matrix: rows = scenarios, cols = populations ─────
# Population order: highest offset first (worst case SSP5-8.5 2061)
pop_order <- pop_offset_summary$Population[
  order(pop_offset_summary$off_585_2061_mean, decreasing = TRUE)
]

# Extract offset values in correct population order
get_offsets <- function(scenario_col) {
  vals <- pop_offset_summary[[scenario_col]]
  names(vals) <- pop_offset_summary$Population
  vals[pop_order]   # reorder to match pop_order
}

pop_matrix <- rbind(
  get_offsets("off_245_2041_mean"),
  get_offsets("off_585_2041_mean"),
  get_offsets("off_245_2061_mean"),
  get_offsets("off_585_2061_mean"),
  get_offsets("off_245_2081_mean"),
  get_offsets("off_585_2081_mean")
)

colnames(pop_matrix) <- pop_order
rownames(pop_matrix) <- c(
  "SSP2-4.5  2041-2060",   # ← fixed typo SSP5-4.5 → SSP2-4.5
  "SSP5-8.5  2041-2060",
  "SSP2-4.5  2061-2080",
  "SSP5-8.5  2061-2080",
  "SSP2-4.5  2081-2100",
  "SSP5-8.5  2081-2100"
)
print(round(pop_matrix, 1))

bar_colors <- c(
  "#cb181d",   # SSP2-4.5 2041  dark blue
  "#fb6a4a",   # SSP5-8.5 2041  light blue
  "#238b45",   # SSP2-4.5 2061  dark red
  "#a1d99b",   # SSP5-8.5 2061  light red
  "#2171b5",   # SSP2-4.5 2081  dark green
  "#9ecae1"    # SSP5-8.5 2081  light green
)


par(mfrow = c(1, 1),
    mar   = c(3, 3, 1, 1))
bp <- barplot(
  pop_matrix,
  beside    = TRUE,
  horiz     = F,
  col       = bar_colors,
  border    = NA,
  names.arg = colnames(pop_matrix),
  las       = 1,
  cex.names = 1,
  cex.axis  = 0.75,
  xlab      = "",
  #xlim      = c(0, max(pop_matrix, na.rm = TRUE) * 1.05)
)

legend(
  x      = 20,
  y      = 1100,
  legend = rownames(pop_matrix),
  fill   = bar_colors,
  border = NA,
  bty    = "n",
  cex    = 0.72,
  title  = expression(bold("")),
  xpd    = TRUE,
  y.intersp = 1.1
)

# ============================================================
# CUSTOM FUNCTIONS REPLACING RDAforest/QuantGenResources
# Only needs: vegan, terra
# ============================================================

# ── Function 1: adapt_scale ───────────────────────────────
# Computes 90th percentile of pairwise distances in gPC space
# Used as scaling factor for env_mismatch

adapt_scale <- function(predictions) {
  pred_mat <- as.matrix(predictions)
  
  # Sample max 500 rows for speed (pairwise = n^2)
  if (nrow(pred_mat) > 500) {
    set.seed(42)
    pred_mat <- pred_mat[sample(nrow(pred_mat), 500), ]
  }
  
  # All pairwise Euclidean distances
  dists <- as.vector(dist(pred_mat, method = "euclidean"))
  
  return(c(
    mean_dist   = mean(dists,             na.rm = TRUE),
    scale_90pct = quantile(dists, 0.90,   na.rm = TRUE)
  ))
}

# ── Function 2: env_mismatch ──────────────────────────────
# Distance from a TARGET (site/individual) to every landscape pixel
# Low value = good genomic match = good rescue donor source

env_mismatch <- function(X,     # target: future gPCs needed (vector or 1-row df)
                         Y,     # oj model (list with $predictions.direct, $goodrows)
                         sy,    # landscape coordinates (lon, lat)
                         sc) {  # scaling factor from adapt_scale()[2]
  
  # Ensure X is numeric vector
  if (is.data.frame(X) | is.matrix(X)) X <- as.numeric(X[1, ])
  X <- as.numeric(X)
  
  # Get landscape predictions
  preds <- as.matrix(Y$predictions.direct)
  
  # Filter to stable rows
  if (!is.null(Y$goodrows)) {
    preds  <- preds[Y$goodrows, ]
    coords <- sy[Y$goodrows, ]
  } else {
    coords <- sy
  }
  
  # Euclidean distance from target to each pixel, scaled
  mismatch <- apply(preds, 1, function(row) {
    sqrt(sum((row - X)^2)) / sc
  })
  
  return(data.frame(
    x            = as.numeric(coords[, 1]),
    y            = as.numeric(coords[, 2]),
    env.mismatch = mismatch
  ))
}

# ── Function 3: ordinationJackknife ──────────────────────
# Core model: maps genotype → environment across landscape
# Jackknife resampling = stable, robust predictions
# Equivalent to the tutorial's ordinationJackknife()

ordinationJackknife <- function(Y,
                                X,
                                newX,
                                covariates,
                                nreps   = 20,
                                top.pcs = 4,
                                extra   = 0.1,
                                verbose = TRUE) {
  library(vegan)
  X          <- as.data.frame(X)
  newX       <- as.data.frame(newX)
  covariates <- as.matrix(covariates)
  n          <- nrow(X)
  n_pixels   <- nrow(newX)
  
  cat("=== Input dimensions ===\n")
  cat("Y (dist matrix):  ", nrow(as.matrix(Y)), "x", ncol(as.matrix(Y)), "\n")
  cat("X (env samples):  ", nrow(X), "x", ncol(X), "\n")
  cat("newX (landscape): ", nrow(newX), "x", ncol(newX), "\n")
  cat("covariates:       ", nrow(covariates), "x", ncol(covariates), "\n")
  cat("top.pcs:          ", top.pcs, "\n")
  cat("========================\n")
  
  if (!all(colnames(X) == colnames(newX))) {
    colnames(newX) <- colnames(X)
    cat("Column names forced to match.\n")
  }
  
  # Storage
  all_preds <- vector("list", nreps)
  good_reps <- logical(nreps)
  ref_pred  <- NULL    # ← initialise HERE, outside the loop
  
  if (verbose) cat("Running", nreps, "jackknife replicates...\n")
  
  for (rep in seq_len(nreps)) {
    if (verbose) cat("  Rep", rep, "/", nreps, "\n")
    
    tryCatch({
      # ── Step A: Jackknife resample ─────────────────────────
      leave_out <- sample(seq_len(n), max(1, floor(n / nreps)))
      keep_idx  <- setdiff(seq_len(n), leave_out)
      
      Y_rep   <- as.matrix(Y)[keep_idx, keep_idx]
      X_rep   <- X[keep_idx, , drop = FALSE]
      cov_rep <- covariates[keep_idx, , drop = FALSE]
      
      # ── Step B: db-RDA (Euclidean → no Cailliez correction) ──
      ord_rep <- capscale(
        as.dist(Y_rep) ~ . + Condition(cov_rep),
        data = X_rep,
        add  = FALSE    # ← FALSE now that we use Euclidean distance
      )
      
      n_constrained <- length(ord_rep$CCA$eig)
      cat("    Constrained axes available:", n_constrained, "\n")
      
      n_axes <- min(top.pcs, n_constrained)
      if (n_axes < top.pcs) cat("    NOTE: using", n_axes, "axes instead of", top.pcs, "\n")
      if (n_axes == 0) stop("No constrained axes produced!")
      
      sc_sites <- scores(ord_rep,
                         scaling = 1,
                         display = "sites",
                         choices = 1:n_axes)
      if (is.vector(sc_sites)) sc_sites <- matrix(sc_sites, ncol = 1)
      
      # ── Step C: RDA of gPC scores ~ environment ───────────
      rda_rep     <- rda(sc_sites ~ ., data = X_rep, scale = FALSE)
      n_rda_axes  <- length(rda_rep$CCA$eig)
      n_pred_axes <- min(n_axes, n_rda_axes)
      
      # ── Step D: Predict at landscape pixels ───────────────
      pixel_preds <- predict(rda_rep,
                             newdata = newX,
                             type    = "lc",
                             scaling = 1)[, 1:n_pred_axes, drop = FALSE]
      
      # ── Step E: Procrustes alignment to first good rep ────
      if (is.null(ref_pred)) {
        # First successful rep becomes the reference — store as-is
        ref_pred         <- pixel_preds
        all_preds[[rep]] <- pixel_preds
        cat("    Rep", rep, "set as Procrustes reference\n")
        
      } else {
        # Align this rep to the reference
        # Subsample for speed (Procrustes only needs rotation, not all pixels)
        n_sub  <- min(500, nrow(pixel_preds))
        set.seed(rep)                              # reproducible subsample
        sub_idx <- sample(nrow(pixel_preds), n_sub)
        
        pro <- procrustes(
          ref_pred[sub_idx, , drop = FALSE],       # reference (fixed)
          pixel_preds[sub_idx, , drop = FALSE],    # this rep (rotated)
          scale     = FALSE,
          symmetric = FALSE
        )
        
        # Apply the rotation matrix to ALL pixels
        pixel_preds_aligned  <- as.matrix(pixel_preds) %*% pro$rotation
        all_preds[[rep]]     <- pixel_preds_aligned
        cat("    Rep", rep, "aligned via Procrustes\n")
      }
      
      good_reps[rep] <- TRUE
      cat("    Rep", rep, "SUCCESS — pred dim:", dim(all_preds[[rep]]), "\n")
      
    }, error = function(e) {
      cat("  Rep", rep, "failed:", conditionMessage(e), "\n")
      good_reps[rep] <<- FALSE
    })
  }
  
  # ── Collect results ──────────────────────────────────────
  good_idx  <- which(good_reps)
  all_preds <- all_preds[good_idx]
  
  cat("\nSuccessful reps:", length(good_idx), "/", nreps, "\n")
  if (length(good_idx) == 0) stop("All replicates failed!")
  
  # Standardise axes across reps
  min_axes  <- min(sapply(all_preds, ncol))
  all_preds <- lapply(all_preds, function(p) p[, 1:min_axes, drop = FALSE])
  
  # ── Average predictions ──────────────────────────────────
  pred_array <- array(
    data = unlist(lapply(all_preds, as.vector)),
    dim  = c(n_pixels, min_axes, length(good_idx))
  )
  mean_preds <- apply(pred_array, c(1, 2), mean, na.rm = TRUE)
  
  # ── Stability filter ─────────────────────────────────────
  pred_var   <- apply(pred_array, c(1, 2), var,  na.rm = TRUE)
  mean_var   <- rowMeans(pred_var, na.rm = TRUE)
  var_thresh <- quantile(mean_var, 1 - extra, na.rm = TRUE)
  good_rows  <- which(mean_var <= var_thresh)
  
  cat("Stable landscape pixels:", length(good_rows), "/", n_pixels, "\n")
  
  colnames(mean_preds) <- paste0("gPC", 1:min_axes)
  
  return(list(
    predictions.direct = as.data.frame(mean_preds),
    goodrows           = good_rows,
    nreps_total        = nreps,
    nreps_good         = length(good_idx),
    top.pcs            = min_axes
  ))
}

stopifnot(length(pop_assignments) == nrow(final_env_data))

# Use concordant adaptive SNPs only
geno <- genotype_matrix[, concordant_candidates]
env  <- as.matrix(env_predictors)   # PD, PW, PC

cat("Samples:        ", nrow(geno), "\n")
cat("Candidate SNPs: ", ncol(geno), "\n")

# ============================================================
# PART 2: BUILD LANDSCAPE PREDICTION GRIDS
# ============================================================

# Current climate grid
envc <- as.data.frame(current_3v, xy = TRUE, na.rm = TRUE)
colnames(envc) <- c("lon", "lat", "PD", "PW", "PC")

# Future climate grid (worst case SSP5-8.5 2061-2080)
envf <- as.data.frame(fut_585_2061_3v, xy = TRUE, na.rm = TRUE)
colnames(envf) <- c("lon", "lat", "PD", "PW", "PC")

cat("Landscape pixels: ", nrow(envc), "\n")

# ============================================================
# PART 3: BUILD THE LANDSCAPE GENOMIC MODEL
# ============================================================

tokeep      <- 4    # matches your K=4 from LFMM

# ── Check how many zero-variance SNPs you have ────────────
snp_variance  <- apply(genotype_matrix, 2, var)
zero_var_snps <- sum(snp_variance == 0)
cat("Zero-variance SNPs found:", zero_var_snps, "\n")
cat("Total SNPs before filter:", ncol(genotype_matrix), "\n")
# ── Remove them ───────────────────────────────────────────
geno_matrix_filtered <- genotype_matrix[, snp_variance > 0]
cat("SNPs after filter:       ", ncol(geno_matrix_filtered), "\n")

# ── Now run PCA safely ────────────────────────────────────
pcs_neutral <- prcomp(geno_matrix_filtered, scale = TRUE)$x[, 1:tokeep]

# Genetic distance matrix
cat("Computing genetic distance matrix...\n")
cordist <- dist(geno, method = "euclidean")

# ── Run ordinationJackknife ───────────────────────────────
cat("Running ordinationJackknife...\n")
cat("Go get a coffee — this takes several minutes!\n")

set.seed(4)
oj <- ordinationJackknife(
  Y          = cordist,
  X          = env,
  newX       = as.matrix(envc[, c("PD", "PW", "PC")]),
  covariates = pcs_neutral,
  nreps      = 20,
  top.pcs    = tokeep,
  extra      = 0.1
)

# Save immediately — expensive to rerun!
saveRDS(oj, "~/Desktop/Projects/Dynastes_grantii_pop/Data/GEA/GeneticOffset/oj_model.rds")
# oj <- readRDS("../Data/GeneticOffset/oj_model.rds")   # reload next session

cat("Model complete!\n")
cat("Stable pixels:", length(oj$goodrows), "/", nrow(envc), "\n")

# ── Scaling factor ────────────────────────────────────────
sc <- adapt_scale(oj$predictions.direct)[2]
cat("Scale factor (90th pct):", round(sc, 4), "\n")

# ── Future gPC raster (predicted genomics under future climate)
rfut <- terra::rast(
  data.frame(cbind(
    envf[oj$goodrows, 1:2],                          
    oj$predictions.direct[oj$goodrows, ]            
  ))
)

# ============================================================
# PART 4: SITE SUITABILITY — BEST DONOR FOR EACH POPULATION
# ============================================================
# Question: "For each sky island population facing high offset,
#            which landscape pixels have the best donor genetics?"

# Population centroid coordinates
pop_coords <- data.frame(
  Population = c("Mt. Lemmon", "Chiricahua Mt.", "Utah",
                 "Mogollon Rim"),
  lon = as.numeric(tapply(gps_data$lon, pop_assignments, mean)),
  lat = as.numeric(tapply(gps_data$lat, pop_assignments, mean))
)

coords_mat <- as.matrix(gps_data[, c("lon", "lat")])

# ── Calculate suitability for each population ─────────────
site_suit_list <- list()

for (i in 1:nrow(pop_coords)) {
  
  pop_name <- pop_coords$Population[i]
  site_xy  <- c(pop_coords$lon[i], pop_coords$lat[i])
  cat("Processing site suitability for:", pop_name, "\n")
  
  # Future gPCs needed at this site
  site_future_gpcs <- terra::extract(
    x      = rfut,
    y      = data.frame(rbind(site_xy)),
    method = "bilinear"
  )[-1]
  
  # Mismatch between future needs and present landscape genetics
  agf <- env_mismatch(
    X  = site_future_gpcs,
    Y  = oj,
    sy = envc[, 1:2],
    sc = sc
  )
  
  agf$Population   <- pop_name
  site_suit_list[[pop_name]] <- agf
}

# ── Plot — 5 population maps ──────────────────────────────

par(mfrow = c(2, 2),
    mar   = c(0, 0, 0, 0),
    oma   = c(0, 0, 0, 0))
points_star <- function(x0, y0, size = 0.1, col = "gold", border = "black", npoints = 5) {
  angles <- seq(0, 2*pi, length.out = 2*npoints + 1)
  radii  <- rep(c(size, size/2), length.out = length(angles))
  x <- x0 + radii * cos(angles)
  y <- y0 + radii * sin(angles)
  polygon(x, y, col = col, border = border)
}

for (pop_name in names(site_suit_list)) {
  
  agf      <- site_suit_list[[pop_name]]
  agf_rast <- terra::rast(agf[, c("x", "y", "env.mismatch")])
  
  # Invert: high = genetically suitable donor
  agf_inv  <- max(values(agf_rast), na.rm = TRUE) - agf_rast
  
  plot(agf_inv,
       col    = viridis(100, option = "inferno"),
       main   = pop_name,
       axes   = FALSE,
       legend = TRUE)
  map('state', fill = FALSE, add = TRUE, col = "black")
  # map.scale(x=-109.8, y=30.6, ratio=FALSE, relwidth=0.2, cex = 1)
  map.scale(x=-109.8, y=30.9, ratio=FALSE, relwidth=0.2, cex = 0.6)
  # Target population (cyan square)
  target <- pop_coords[pop_coords$Population == pop_name, ]
  points(target$lon, target$lat,
         pch = 0, col = "chartreuse", cex = 3.5, lwd = 2)
  
  # All sample locations (white dots)
  # points(coords_mat,
  #        pch = 16, col = "white", cex = 0.9)
  for (j in 1:27){
    points_star(coords_mat[,1][j], coords_mat[,2][j], col = "skyblue", size = 0.3)
  }
  
}

# High suitability pixels should make biological sense
# e.g. similar elevation, forest cover, not desert lowlands
plot(agf_inv, main = "",xlim = c(-114.5, -106.5), ylim = c(30, 38))
map('state', fill = FALSE, add = TRUE, col = "black")
map.scale(x=-109.8, y=30.9, ratio=FALSE, relwidth=0.2, cex = 0.6)


# ============================================================
# PART 6: SUMMARY TABLE — CONSERVATION RECOMMENDATIONS
# ============================================================

pop_offset_summary <- read.csv("~/Desktop/Projects/Dynastes_grantii_pop/Data/GEA/GeneticOffset/population_offset_summary.csv")
pop_offset_summary$Population <- c('Mt. Lemmon',"Chiricahua Mt.", "Mogollon Rim","Utah"  )
# For each high-offset population, find the best donor pixel
recommendations <- data.frame()

for (pop_name in names(site_suit_list)) {
  
  agf      <- site_suit_list[[pop_name]]
  
  # Best donor = pixel with LOWEST mismatch
  best_idx <- which.min(agf$env.mismatch)
  
  # Offset rank from your previous analysis
  offset_rank <- which(pop_offset_summary$Population == pop_name)
  
  rec <- data.frame(
    Population        = pop_name,
    Offset_rank       = offset_rank,
    Offset_585_2061   = pop_offset_summary$off_585_2061_mean[offset_rank],
    Best_donor_lon    = round(agf$x[best_idx], 4),
    Best_donor_lat    = round(agf$y[best_idx], 4),
    Min_mismatch      = round(agf$env.mismatch[best_idx], 4)
  )
  recommendations <- rbind(recommendations, rec)
}

recommendations <- recommendations %>%
  arrange(Offset_rank)

cat("\n=== CONSERVATION RECOMMENDATIONS ===\n")
print(recommendations)
for (pop_name in names(site_suit_list)) {
  agf <- site_suit_list[[pop_name]]
  cat(pop_name, "\n")
  cat("  Top 5% suitability pixels — mean lon/lat:\n")
  top5 <- agf[agf$env.mismatch <= quantile(agf$env.mismatch, 0.05), ]
  cat("  lon range:", round(range(top5$x), 2), "\n")
  cat("  lat range:", round(range(top5$y), 2), "\n\n")
}
write.csv(recommendations,
          "~/Desktop/Projects/Dynastes_grantii_pop/Data/GEA/GeneticOffset/rescue_recommendations.csv",
          row.names = FALSE)

