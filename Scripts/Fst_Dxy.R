library(dplyr)
#######
# Fst #
#######
fst_data_raw <- read.table("../Data/Dxy/pixy_fst.txt", header = TRUE, stringsAsFactors = FALSE)
fst_data <- fst_data_raw[is.finite(fst_data_raw$avg_wc_fst), ]
chrom <- unique(fst_data$chromosome)[order(unique(fst_data$chromosome))][1:11]
fst_data <- fst_data[fst_data$chromosome %in% chrom,]
chr_levels <- sort(unique(fst_data$chromosome))
fst_data <- fst_data[order(factor(fst_data$chromosome, levels = chr_levels), fst_data$window_pos_1), ]
fst_data$mid <- (fst_data$window_pos_1 + fst_data$window_pos_2) / 2
chr_max <- tapply(fst_data$window_pos_2, fst_data$chromosome, max)
chr_off <- c(0, cumsum(as.numeric(chr_max[chr_levels])))
chr_off <- chr_off[-length(chr_off)]
names(chr_off) <- chr_levels
fst_data$genome_pos <- fst_data$mid + chr_off[fst_data$chromosome]

list <- unique(paste0(fst_data$pop1, '_', fst_data$pop2))
com.list <- c()
for (i in 1:length(list)){
  pop1 <- strsplit(list[i], split = '_')[[1]][1]
  pop2 <- strsplit(list[i], split = '_')[[1]][2]
  com.name <- paste0(strsplit(pop1, split = '')[[1]][1], strsplit(pop2, split = '')[[1]][1])
  com.list <- c(com.list, com.name)
  temp.dat <- fst_data[fst_data$pop1 == pop1 & fst_data$pop2 == pop2, ]
  assign(com.name, temp.dat, envir = .GlobalEnv)
}

cols <- c("#3030D0", "#9090F0")
par(mfrow = c(6, 1),            
    mar   = c(1.5, 4, 0.5, 3),
    oma   = c(0, 0, 0, 0),
    mgp   = c(2.2, 0.6, 0),
    xaxs  = "i", yaxs = "i")
chr_col <- setNames(cols[(seq_along(chr_levels)-1) %% 2 + 1], chr_levels)
for (i in 1:length(com.list)){
  df <- get(com.list[i])
  plot(df$genome_pos, df$avg_wc_fst,
       pch = 16, cex = 0.45,
       col = chr_col[df$chromosome],
       ylim = c(0, 1.1),
       xlab = "Genome position", ylab = "Fst", xaxt = "n")
  mtext(paste0(strsplit(list[i], split = '_')[[1]][1],' vs. ',strsplit(list[i], split = '_')[[1]][2]), side = 4, line = 0.5, las = 0, cex = 0.6)
  # centers <- chr_off + chr_max[chr_levels] / 2
  # axis(1, at = centers, labels = chr_levels, las = 2)
}
# plot highlight data points
highlight_keys <- paste(
  fully_overlapping_windows$chromosome,
  fully_overlapping_windows$window_pos_1
)
for (i in 1:length(com.list)){
  df <- get(com.list[i])
  
  # --- CHANGE 1: PLOT ALL POINTS NORMALLY FIRST ---
  # This is your original plotting call, which draws the background.
  plot(df$genome_pos, df$avg_wc_fst,
       pch = 16, cex = 0.45,
       col = chr_col[df$chromosome], # Use the original alternating blue colors
       ylim = c(0, 1.1),
       xlab = "Genome position", ylab = "Fst", xaxt = "n")
  
  # --- CHANGE 2: SUBSET AND PLOT THE HIGHLIGHT POINTS ON TOP ---
  # 1. Create a logical index (TRUE/FALSE) for the highlight points.
  is_highlight <- paste(df$chromosome, df$window_pos_1) %in% highlight_keys
  
  # 2. Subset the data frame to get only the highlight points.
  highlight_df <- df[is_highlight, ]
  
  # 3. Use points() to draw them on top of the existing plot.
  # We can also make them slightly larger (e.g., cex = 0.6) to stand out.
  if (nrow(highlight_df) > 0) {
    points(highlight_df$genome_pos, highlight_df$avg_wc_fst,
           pch = 16,
           cex = 0.6, # Make them a bit bigger
           col = "red")
  }
  mtext(paste0(strsplit(list[i], split = '_')[[1]][1],' vs. ',strsplit(list[i], split = '_')[[1]][2]), side = 4, line = 0.5, las = 0, cex = 0.6)
}


# Focus on Utah to others
Utah <- com.list[c(1,3,6)]
Utah.list <- list[c(1,3,6)]
par(mfrow = c(3, 1),            
    mar   = c(1.5, 4, 0.5, 3),
    oma   = c(0, 0, 0, 0),
    mgp   = c(2.2, 0.6, 0),
    xaxs  = "i", yaxs = "i")
for (i in 1:length(Utah)){
  df <- get(Utah[i])
  plot(df$genome_pos, df$avg_wc_fst,
       pch = 16, cex = 0.45,
       col = chr_col[df$chromosome],
       ylim = c(0, 1.1),
       xlab = "Genome position", ylab = "Fst", xaxt = "n")
  mtext(Utah.list[i], side = 4, line = 0.5, las = 0, cex = 0.6)
  # centers <- chr_off + chr_max[chr_levels] / 2
  # axis(1, at = centers, labels = chr_levels, las = 2)
  chr14 <- df[df$chromosome == "JAROYF020000014.1",]
  hit  <- chr14[chr14$avg_wc_fst >= 0.9,]
  print(hit)
  points(hit$genome_pos, hit$avg_wc_fst,
         col = "red", pch = 16, cex = 0.7)
}
cutoff <- length(df[df$chromosome %in% unique(df$chromosome)[1:11],]$avg_wc_fst)/100*0.1
df[df$chromosome %in% unique(df$chromosome)[1:11],]$avg_wc_fst[order(df[df$chromosome %in% unique(df$chromosome)[1:11],]$avg_wc_fst, decreasing = T)][1:cutoff]
#######
# Dxy #
#######

Dxy_data_raw <- read.table("../Data/Dxy/pixy_dxy.txt", header = TRUE, stringsAsFactors = FALSE)
Dxy_data <- Dxy_data_raw[is.finite(Dxy_data_raw$avg_dxy), ]
chrom <- unique(Dxy_data$chromosome)[order(unique(Dxy_data$chromosome))][1:11]
Dxy_data <- Dxy_data[Dxy_data$chromosome %in% chrom,]
chr_levels <- sort(unique(Dxy_data$chromosome))
Dxy_data <- Dxy_data[order(factor(Dxy_data$chromosome, levels = chr_levels), Dxy_data$window_pos_1), ]
Dxy_data$mid <- (Dxy_data$window_pos_1 + Dxy_data$window_pos_2) / 2
chr_max <- tapply(Dxy_data$window_pos_2, Dxy_data$chromosome, max)
chr_off <- c(0, cumsum(as.numeric(chr_max[chr_levels])))
chr_off <- chr_off[-length(chr_off)]
names(chr_off) <- chr_levels
Dxy_data$genome_pos <- Dxy_data$mid + chr_off[Dxy_data$chromosome]
list <- unique(paste0(Dxy_data$pop1,'_',Dxy_data$pop2))
com.list <- c()

for (i in 1:length(list)){
  pop1 <- strsplit(list[i], split = '_')[[1]][1]
  pop2 <- strsplit(list[i], split = '_')[[1]][2]
  com.name <- paste0(strsplit(pop1, split = '')[[1]][1], strsplit(pop2, split = '')[[1]][1])
  com.list <- c(com.list, com.name)
  temp.dat <- Dxy_data[Dxy_data$pop1 == pop1 & Dxy_data$pop2 == pop2, ]
  assign(com.name, temp.dat, envir = .GlobalEnv)
}

cols <- c("#3030D0", "#9090F0")
par(mfrow = c(6, 1),            
    mar   = c(1.5, 4, 0.5, 3),
    oma   = c(0, 0, 0, 0),
    mgp   = c(2.2, 0.6, 0),
    xaxs  = "i", yaxs = "i")
chr_col <- setNames(cols[(seq_along(chr_levels)-1) %% 2 + 1], chr_levels)
for (i in 1:length(com.list)){
  df <- get(com.list[i])
  
  plot(df$genome_pos, df$avg_dxy,
       pch = 16, cex = 0.45,
       col = chr_col[df$chromosome],
       ylim = c(0, 1.1),
       xlab = "Genome position", ylab = "Dxy", xaxt = "n")
  mtext(list[i], side = 4, line = 0.5, las = 0, cex = 0.6)
  # centers <- chr_off + chr_max[chr_levels] / 2
  # axis(1, at = centers, labels = chr_levels, las = 2)
}



#############
# Fst - Dxy #
#############

combined_data <- inner_join(
  fst_data %>% select(pop1, pop2, chromosome, window_pos_1, window_pos_2, avg_wc_fst),
  Dxy_data %>% select(pop1, pop2, chromosome, window_pos_1, window_pos_2, avg_dxy),
  by = c("pop1", "pop2", "chromosome", "window_pos_1", "window_pos_2")
)
# Remove any rows with non-finite values (NA, NaN, Inf) in Fst or Dxy
combined_data <- combined_data %>%
  filter(is.finite(avg_wc_fst) & is.finite(avg_dxy))

# Get a list of unique population comparisons
pop_pairs <- unique(paste0(combined_data$pop1, "_", combined_data$pop2))
par(mfrow = c(6, 1),            
    mar   = c(1.5, 4, 0.5, 3),
    oma   = c(0, 0, 0, 0),
    mgp   = c(2.2, 0.6, 0),
    xaxs  = "i", yaxs = "i")
for (pair in pop_pairs) {
  
  # Subset the data for the current population pair
  pop1_name <- strsplit(pair, "_")[[1]][1]
  pop2_name <- strsplit(pair, "_")[[1]][2]
  
  subset_data <- combined_data %>%
    filter(pop1 == pop1_name & pop2 == pop2_name)
  
  # --- Step 3: Identify Outlier Windows (Top 1% Fst) ---
  
  # Calculate the Fst threshold for the top 1%
  # We use na.rm = TRUE just in case, though we already filtered.
  fst_threshold <- quantile(subset_data$avg_wc_fst, probs = 0.99, na.rm = TRUE)
  
  # Identify the outlier windows
  outlier_windows <- subset_data %>%
    filter(avg_wc_fst >= fst_threshold) %>%
    arrange(desc(avg_wc_fst)) # Order by Fst for clarity
  
  # (Optional) If you want to save these outlier regions to a file:
  # write.csv(outlier_windows, file = paste0("outliers_", pair, ".csv"), row.names = FALSE)
  
  
  # --- Step 4: Create Scatterplot of Fst vs. Dxy ---
  
  # Add a new column to the subset_data to mark outliers for plotting
  subset_data <- subset_data %>%
    mutate(is_outlier = avg_wc_fst >= fst_threshold)
  
  # Create the plot title
  plot_title <- paste("Fst vs. Dxy for", pair)
  plot_colors <- ifelse(subset_data$is_outlier, "black", "gray")
  
  plot(
    x = subset_data$avg_dxy,
    y = subset_data$avg_wc_fst,
    xlim = c(0,1), ylim = c(0,1),
    main = '',
    xlab = "Dxy",
    ylab = "Fst",
    pch = 16,          # Solid circle points
    cex = 0.8,         # Adjust point size
    col = plot_colors  # Apply the custom colors
  )
  mtext(paste0(strsplit(pair, split = '_')[[1]][1],' v.s. ',strsplit(pair, split = '_')[[1]][2]), side = 4, line = 0.5, las = 0, cex = 0.6)
  abline(h = fst_threshold, col = "blue", lty = 2)
  
}

# Loop through each pair to find candidates
meta.dat <- c()
for (pair in pop_pairs) {
  
  # --- Setup for the current population pair ---
  cat("==============================================================\n")
  cat("CANDIDATE ANALYSIS FOR POPULATION PAIR:", pair, "\n")
  cat("==============================================================\n\n")
  
  pop1_name <- strsplit(pair, "_")[[1]][1]
  pop2_name <- strsplit(pair, "_")[[1]][2]
  
  subset_data <- combined_data %>%
    filter(pop1 == pop1_name & pop2 == pop2_name)
  
  # --- Step 1: Define High-Fst Candidates ---
  
  # Define the Fst percentile threshold (e.g., 99th percentile for top 1%)
  fst_threshold_high <- quantile(subset_data$avg_wc_fst, probs = 0.99, na.rm = TRUE)
  
  # Get the initial pool of high-Fst outlier windows
  high_fst_candidates <- subset_data %>%
    filter(avg_wc_fst >= fst_threshold_high)
  
  cat("Total windows with Fst in the top 1% ( >=", round(fst_threshold_high, 4), "): ", nrow(high_fst_candidates), "\n\n")
  
  
  # --- METHOD 1: Two-Step Percentile ---

  # From the high-Fst candidates, find the Dxy threshold for the bottom 10%
  dxy_threshold_low <- quantile(high_fst_candidates$avg_dxy, probs = 0.10, na.rm = TRUE)

  # Filter for the final list of "High Fst / Low Dxy" candidates
  candidates_percentile <- high_fst_candidates %>%
    filter(avg_dxy <= dxy_threshold_low) %>%
    arrange(desc(avg_wc_fst)) # Order by Fst for clarity

  cat("--- Method 1: Two-Step Percentile Candidates ---\n")
  cat("(Top 1% Fst AND Bottom 10% Dxy among those outliers)\n")
  cat("Dxy threshold for this group (<= ", round(dxy_threshold_low, 4), ")\n")
  print(candidates_percentile)
  meta.dat <- rbind(meta.dat,candidates_percentile)
  cat("\n")
  
  
  # # --- METHOD 2: Z-Score ---
  # 
  # # Calculate genome-wide Dxy stats for this pair
  # mean_dxy <- mean(subset_data$avg_dxy, na.rm = TRUE)
  # sd_dxy <- sd(subset_data$avg_dxy, na.rm = TRUE)
  # 
  # # Calculate Z-score for all windows and then filter
  # candidates_zscore <- subset_data %>%
  #   mutate(z_dxy = (avg_dxy - mean_dxy) / sd_dxy) %>%
  #   filter(avg_wc_fst >= fst_threshold_high & z_dxy < -1.5) %>% # Apply both filters at once
  #   arrange(desc(avg_wc_fst)) # Order by Fst
  # 
  # cat("--- Method 2: Z-Score Candidates ---\n")
  # cat("(Top 1% Fst AND Dxy Z-score < -1.5)\n")
  # cat("Genome-wide Mean Dxy:", round(mean_dxy, 4), "\n")
  # print(candidates_zscore)
  # cat("\n\n")
  # meta.dat <- rbind(meta.dat,candidates_zscore)
}

utah_candidates <- meta.dat %>%
  filter(pop1 == 'Utah' | pop2 == 'Utah')

# Step 2: Determine how many total comparisons 'Utah' is a part of.
# This gives us the number we need for a complete overlap.
# For example, if Utah was compared to Portal, Lemmon, and Mogollon, this number will be 3.
num_utah_comparisons <- utah_candidates %>%
  # We create a temporary column to uniquely identify each pair
  mutate(pair_id = paste(pmin(pop1, pop2), pmax(pop1, pop2), sep = "_")) %>%
  summarise(n_distinct(pair_id)) %>%
  pull() # pull() extracts the single value (the count)

cat("Number of comparisons involving 'Utah':", num_utah_comparisons, "\n\n")


# Step 3 & 4: Group by genomic window and count the occurrences.
window_counts <- utah_candidates %>%
  group_by(chromosome, window_pos_1, window_pos_2) %>%
  summarise(overlap_count = n(), .groups = 'drop') %>%
  arrange(desc(overlap_count)) # Sort to see the most frequent overlaps first

# You can inspect this intermediate result
cat("Count of how many times each candidate window appears:\n")
print(head(window_counts))

# Step 5: Filter for windows that appear in ALL comparisons.
# These are your highest-confidence candidates for selection on the Utah lineage.
fully_overlapping_windows <- window_counts %>%
  filter(overlap_count == num_utah_comparisons)

cat("--- Highest-Confidence Candidates: Windows significant in ALL", num_utah_comparisons, "Utah comparisons ---\n")
print(fully_overlapping_windows)


# Step 6 (Optional but Recommended): Get the full data for these top windows.
# This shows you the Fst, Dxy, etc. values for the overlapping windows in each comparison.
if (nrow(fully_overlapping_windows) > 0) {
  top_candidate_details <- semi_join(utah_candidates, fully_overlapping_windows, 
                                     by = c("chromosome", "window_pos_1", "window_pos_2")) %>%
    arrange(chromosome, window_pos_1, pop1) # Arrange for easy reading
  
  cat("--- Full details for the highest-confidence candidate windows ---\n")
  print(top_candidate_details)
} else {
  cat("No windows were found to be significant across all", num_utah_comparisons, "comparisons.\n")
}

# plot with those highlight data points

if (!exists("fully_overlapping_windows") || nrow(fully_overlapping_windows) == 0) {
  stop("The 'fully_overlapping_windows' object is empty or does not exist. Please run the overlap analysis first.")
}

top_candidate_keys <- paste(
  fully_overlapping_windows$chromosome,
  fully_overlapping_windows$window_pos_1,
  sep = "_"
)

# --- Step 2: Loop Through Each Population Pair to Plot ---

# Get a list of unique population comparisons from your merged data
pop_pairs <- unique(paste0(combined_data$pop1, "_", combined_data$pop2))

# Loop through each pair
for (pair in pop_pairs) {
  
  # Subset the data for the current population pair
  pop1_name <- strsplit(pair, "_")[[1]][1]
  pop2_name <- strsplit(pair, "_")[[1]][2]
  
  subset_data <- combined_data %>%
    filter(pop1 == pop1_name & pop2 == pop2_name)
  
  # --- Step 3: Create a new 'highlight_status' column for coloring ---
  
  # Calculate the Fst threshold for this specific pair
  fst_threshold_high <- quantile(subset_data$avg_wc_fst, probs = 0.99, na.rm = TRUE)
  
  # Create a unique key for each window in the current subset
  subset_data$key <- paste(subset_data$chromosome, subset_data$window_pos_1, sep = "_")
  
  # Use case_when() to assign each window to a category for plotting
  subset_data <- subset_data %>%
    mutate(
      highlight_status = case_when(
        key %in% top_candidate_keys ~ "Top Utah Candidate",
        avg_wc_fst >= fst_threshold_high ~ "Other Top 1% Outlier",
        TRUE ~ "Background"
      )
    )
  
  # Convert to a factor to control the plotting order (important!)
  # This ensures 'Top Utah Candidate' points are plotted last (on top).
  subset_data$highlight_status <- factor(
    subset_data$highlight_status,
    levels = c("Background", "Other Top 1% Outlier", "Top Utah Candidate")
  )
  
  # Re-order the data frame based on this factor
  subset_data <- subset_data %>%
    arrange(highlight_status)
  
  # --- Step 4: Create the Scatterplot with new colors ---
  
  # Define the colors that correspond to the factor levels
  plot_colors <- c(
    "Background" = "grey", 
    "Other Top 1% Outlier" = "black", 
    "Top Utah Candidate" = "red"
  )
  
  # Create the plot
  plot(
    x = subset_data$avg_dxy,
    y = subset_data$avg_wc_fst,
    main = '',
    xlab = "Dxy",
    ylab = "Fst",
    xlim = c(0,1), ylim = c(0,1),
    pch = 16, # Use solid circles for better visibility
    cex = 0.8,
    col = plot_colors[subset_data$highlight_status] # Assign colors based on status
  )
  mtext(paste0(strsplit(pair, split = '_')[[1]][1],' v.s. ',strsplit(pair, split = '_')[[1]][2]), side = 4, line = 0.5, las = 0, cex = 0.6)
  # Add the threshold line
  abline(h = fst_threshold_high, col = "blue", lty = 2)
}

#######################
# Mann-Whitney U test #
#######################

# making the dataframe
processed_data <- combined_data %>%
  # Explicitly tell R to use the 'rename' function FROM the 'dplyr' package
  dplyr::rename(
    Fst_value = avg_wc_fst,
    Dxy_value = avg_dxy
  ) %>%
  # Be explicit with mutate as well, just to be safe
  dplyr::mutate(
    ComparisonType = case_when(
      pop1 == "Utah" | pop2 == "Utah" ~ "Utah-involved",
      TRUE ~ "South-only"
    )
  )

#------------------------------------------------------------------
# STEP 1: PERFORM THE STATISTICAL TESTS
#------------------------------------------------------------------

# --- Test 1: Compare Fst Distributions ---
# Hypothesis: Fst in 'Utah-involved' pairs is greater than in 'South-only' pairs.
# This tests for the effect of drift/isolation.

# We use the formula interface for wilcox.test, which is clean and easy to read.
# The formula y ~ x means "test variable y grouped by variable x".
fst_test_result <- wilcox.test(Fst_value ~ ComparisonType, 
                               data = processed_data, 
                               alternative = "less")

# Print a clear summary of the Fst test result
cat("--- Fst 'Drift' Test ---\n")
print(fst_test_result)
cat("\nThis p-value should be extremely small, confirming the 'drift' story.\n")


# --- Test 2: Compare Dxy Distributions ---
# Hypothesis: Dxy distributions are NOT different between the two groups.
# This supports the 'congruent history' story.

# We use a 'two.sided' test because we're looking for any difference.
dxy_test_result <- wilcox.test(Dxy_value ~ ComparisonType, 
                               data = processed_data, 
                               alternative = "two.sided")

# Print a clear summary of the Dxy test result
cat("\n--- Dxy 'Congruent History' Test ---\n")
print(dxy_test_result)
cat("\nThis p-value should be non-significant (e.g., > 0.05), supporting the 'congruent history' story.\n")

#------------------------------------------------------------------
# STEP 1: SET UP THE PLOTTING AREA
#------------------------------------------------------------------

# Set up a 1x2 plotting grid
par(mfrow = c(1, 2), mar = c(5, 5, 4, 1))

#------------------------------------------------------------------
# STEP 2: CREATE THE FST PLOT WITH POINTS (LEFT PANEL)
#------------------------------------------------------------------

# First, draw the boxplot but WITHOUT the outlier points (outline = FALSE)
boxplot(Fst_value ~ factor(ComparisonType, levels = c("South-only", "Utah-involved")),
        data = processed_data,
        main = "",
        xlab = "",
        ylab = "Fst",
        col = "lightgray",
        outline = FALSE,
        cex.axis = 0.8 )
boxplot(Dxy_value ~ factor(ComparisonType, levels = c("South-only", "Utah-involved")),
        data = processed_data,
        main = "",
        xlab = "",
        ylab = "Dxy",
        col = "lightgray",
        outline = FALSE,
        cex.axis = 0.8)


######
# Pi #
######
Pi_data_raw <- read.table("../Data/Dxy/pixy_pi.txt", header = TRUE, stringsAsFactors = FALSE)
Pi_data <- Pi_data_raw[is.finite(Pi_data_raw$avg_pi), ]
keep <- sort(unique(Pi_data$chromosome))[1:11]
Pi_data <- Pi_data[Pi_data$chromosome %in% keep,]
pop <- unique(Pi_data$pop)
for (i in 1:length(pop)){
  temp.dat <- Pi_data[Pi_data$pop == pop[i],]
  temp.dat <- temp.dat[order(temp.dat$chromosome, temp.dat$window_pos_1), ]
  temp.dat$mid <- (temp.dat$window_pos_1 + temp.dat$window_pos_2) / 2
  temp.dat <- temp.dat[order(temp.dat$chromosome, temp.dat$window_pos_1), ]
  temp.dat$mid <- (temp.dat$window_pos_1 + temp.dat$window_pos_2) / 2
  chr_levels <- unique(temp.dat$chromosome)
  chr_max <- tapply(temp.dat$window_pos_2, temp.dat$chromosome, max)
  chr_off <- c(0, cumsum(chr_max[chr_levels]))
  chr_off <- chr_off[seq_along(chr_levels)]
  names(chr_off) <- chr_levels
  temp.dat$genome_pos <- temp.dat$mid + chr_off[temp.dat$chromosome]
  assign(paste0(pop[i],'.pi'), temp.dat, envir = .GlobalEnv)
}
cols <- c("#3030D0", "#9090F0")
par(mfrow = c(4, 1),            
    mar   = c(1.5, 4, 0.5, 3),
    oma   = c(0, 0, 0, 0),
    mgp   = c(2.2, 0.6, 0),
    xaxs  = "i", yaxs = "i")
chr_col <- setNames(cols[(seq_along(chr_levels)-1) %% 2 + 1], chr_levels)
for (i in 1:length(pop)){
  df <- get(paste0(pop[i],'.pi'))
  
  plot(df$genome_pos, df$avg_pi,
       pch = 16, cex = 0.45,
       col = chr_col[df$chromosome],
       ylim = c(0, 0.65),
       xlab = "", ylab = "Pi", xaxt = "n")
  mtext(pop[i], side = 4, line = 0.5, las = 0, cex = 0.7)
 
  # centers <- chr_off + chr_max[chr_levels] / 2
  # axis(1, at = centers, labels = chr_levels, las = 2)
}

df <- Utah.pi
df <- df[df$chromosome == "JAROYF020000003.1",]
plot(df$window_pos_1, df$avg_pi, type = 'l',
     pch = 16, cex = 0.45,
     col = chr_col[df$chromosome],
     ylim = c(0, 0.65),
     xlab = "", ylab = "Pi", xaxt = "n",
     xlim = c(58500000,59000000))
abline(v=c(58790000,58795000,58910001,58915000))
df[df$window_pos_1 == 58790001,]
df[df$window_pos_1 == 58910001,]
