# Nucleotide diversity & Tajuma's D 
library(viridis)
library(vioplot)
library(dplyr)
library(coda)

results_dir <- "../Data/tajima" 
population_names <- c("Utah", "Lemmon", "Portal", "Mogollon")
pi_data_list <- list()
tajimaD_data_list <- list()

for (i in 1:length(population_names)){
  pi_file_path <- file.path(results_dir, paste0(population_names[i], "_pi.windowed.pi"))
  pi_table <- read.table(pi_file_path, header = TRUE, stringsAsFactors = FALSE)
  pi_values <- pi_table$PI[is.finite(pi_table$PI)]
  pi_data_list[[population_names[i]]] <- pi_values
  tajimaD_file_path <- file.path(results_dir, paste0(population_names[i], "_tajimaD.Tajima.D"))
  tajimaD_table <- read.table(tajimaD_file_path, header = TRUE, stringsAsFactors = FALSE)
  tajimaD_values <- tajimaD_table$TajimaD[is.finite(tajimaD_table$TajimaD)]
  tajimaD_data_list[[population_names[i]]] <- tajimaD_values
}

# Plot 
results_dir <- "../Data/tajima"
population_names <- c("Utah", "Lemmon", "Portal", "Mogollon")

list_of_data_frames <- list()

for (pop_name in population_names) {
  # Construct the file path
  file_path <- file.path(results_dir, paste0(pop_name, "_tajimaD.Tajima.D"))
  
  # Read the table
  # Use tryCatch to handle cases where a file might be missing or empty
  tryCatch({
    pop_data <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE)
    
    # *** The Key Step: Add a column for the population name ***
    pop_data$pop <- pop_name
    
    # Add the complete data frame to our list
    list_of_data_frames[[pop_name]] <- pop_data
  }, error = function(e) {
    message(paste("Could not read or process file for population:", pop_name))
    message("Error:", e$message)
  })
}

combined_df <- bind_rows(list_of_data_frames)

combined_df <- combined_df %>%
  rename(
    chromosome = CHROM,
    window_start = BIN_START,
    n_snps = N_SNPS,
    tajima_d = TajimaD
  )

window_size <- 5000 
clean_df <- combined_df %>%
  filter(is.finite(tajima_d)) %>%
  mutate(window_mid = window_start + (window_size / 2))

head(clean_df)

pop_summaries <- clean_df %>%
  group_by(pop) %>%
  summarise(
    mean_d = mean(tajima_d),
    sd_d = sd(tajima_d),
    .groups = 'drop' # Good practice to drop grouping after summarising
  )

print(pop_summaries)

plot_list <- list()

top_11_chroms <- clean_df %>%
  count(chromosome, sort = TRUE) %>%
  slice(1:11) %>%
  pull(chromosome)

population_names <- c("Mogollon", "Lemmon", "Portal", "Utah")
plot_colors <- viridis(4,alpha = 0.6)
names(plot_colors) <- population_names
par(mfrow = c(4, 3), mar = c(4, 4, 2, 1))

for (chrom_name in top_11_chroms) {
  
  cat(paste("Generating base R plot for chromosome:", chrom_name, "\n"))
  
  # Get all data for the current chromosome
  chrom_data <- subset(clean_df, chromosome == chrom_name)
  
  # --- CRITICAL: Determine plot limits BEFORE plotting ---
  # In base R, you must define the full range of the plot first
  # so that all overlaid lines will fit.
  xlim_range <- range(chrom_data$window_mid / 1e6, na.rm = TRUE)
  ylim_range <- range(chrom_data$tajima_d, na.rm = TRUE)
  
  # Add a little buffer to the y-axis
  ylim_range[1] <- ylim_range[1] - 0.2
  ylim_range[2] <- ylim_range[2] + 0.2
  
  # --- 1. Initialize an empty plot ---
  # 'type = "n"' creates the axes and labels but plots no data.
  # This sets up our "canvas" for this chromosome.
  plot(1, 
       type = "n", 
       xlim = xlim_range, 
       ylim = ylim_range,
       main = '',
       xlab = "Position (Mb)",
       ylab = "Tajima's D",
       cex.main = 1, # Adjust title size
       cex.lab = 0.7, # Adjust axis label size
       cex.axis = 0.7 # Adjust axis tick mark size
  )
  for (pop_name in population_names) {
    pop_chrom_data <- subset(chrom_data, pop == pop_name)
    points(pop_chrom_data$window_mid / 1e6, 
           pop_chrom_data$tajima_d, 
           col = plot_colors[pop_name], 
           pch =16, cex = 0.1)
    abline(h = 0, lty = "dotted", col = "black")
  }
  
}

par(mfrow = c(1, 1))

results_with_z <- clean_df %>%
  left_join(pop_summaries, by = "pop") %>%
  mutate(z_score = (tajima_d - mean_d) / sd_d)

manhattan_data <- results_with_z %>%
  
  group_by(chromosome) %>%
  summarise(chr_len = max(window_mid), .groups = 'drop') %>%
  arrange(chromosome) %>% 
  mutate(
    total_len = cumsum(as.numeric(chr_len)) - chr_len
  ) %>%
  
  right_join(results_with_z, by = "chromosome") %>%
  mutate(bp_cumulative = window_mid + total_len)


axis_df <- manhattan_data %>% 
  group_by(chromosome) %>% 
  summarize(center = mean(bp_cumulative), .groups = 'drop')


par(mar = c(6, 4.5, 3, 2))

# --- Initialize an empty plot ---
plot(
  x = manhattan_data$bp_cumulative, 
  y = abs(manhattan_data$z_score),
  type = "n",
  xaxt = "n", # Suppress the default x-axis
  main = "",
  xlab = "",  # X-axis label will be added manually
  ylab = "Absolute Z-score",
  cex.main = 1.2,
  cex.lab = 1
)

abline(h = 3, lty = "dashed", col = "red", lwd = 1)
population_names <- c("AdmixedZone", "Lemmon", "Portal", "Outgroup")
plot_colors <- viridis(4,alpha = 0.3)
names(plot_colors) <- population_names

for (pop_name in population_names) {
  pop_data_subset <- subset(manhattan_data, pop == pop_name)
  points(
    x = pop_data_subset$bp_cumulative,
    y = abs(pop_data_subset$z_score),
    col = plot_colors[pop_name],
    pch = 16, # Solid circle
    cex = 0.7
  )
}

# --- Manually draw the custom x-axis ---
axis(side = 1, at = axis_df$center, labels = FALSE)
text(
  x = axis_df$center[1:11],
  y = par("usr")[3] - (par("usr")[4] * 0.05), # Position labels below the plot area
  labels = axis_df$chromosome[1:11],
  srt = 45, # Rotate labels 45 degrees
  adj = 1,
  xpd = TRUE, # Allow drawing outside the plot region
  cex = 0.8
)
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
