# Plot 
# Pi 
# Dxy
# Fst
library(viridis)
library(vioplot)
pi_file <- "../Data/Dxy/pixy_pi.txt"
dxy_file <- "../Data/Dxy/pixy_dxy.txt"
fst_file <- "../Data/Dxy/pixy_fst.txt"

# Load the data into R data frames
pi_data_raw <- read.table(pi_file, header = TRUE, stringsAsFactors = FALSE)
dxy_data_raw <- read.table(dxy_file, header = TRUE, stringsAsFactors = FALSE)
fst_data_raw <- read.table(fst_file, header = TRUE, stringsAsFactors = FALSE)

# --- NEW ROBUST DATA CLEANING ---
# We will filter each data frame to keep only the rows where the
# statistic of interest (e.g., avg_pi, avg_fst) is a finite number.
# is.finite() correctly handles NA, NaN, and Inf all at once.

pi_data <- pi_data_raw[is.finite(pi_data_raw$avg_pi), ]
dxy_data <- dxy_data_raw[is.finite(dxy_data_raw$avg_dxy), ]
fst_data <- fst_data_raw[is.finite(fst_data_raw$avg_wc_fst), ]
pi_data$pop <- factor(pi_data$pop, levels = c('Mogollon','Lemmon','Portal', 'Utah'))
# --- Create the pop_pair column AFTER cleaning ---
dxy_data$pop_pair <- paste(dxy_data$pop1, dxy_data$pop2, sep = "-")
fst_data$pop_pair <- paste(fst_data$pop1, fst_data$pop2, sep = "-")

# For Pi
avg_pi <- tapply(pi_data$avg_pi * pi_data$no_sites, pi_data$pop, sum) / tapply(pi_data$no_sites, pi_data$pop, sum)
cat("Genome-wide average Nucleotide Diversity (pi):\n")
print(avg_pi)

# For dxy
avg_dxy <- tapply(dxy_data$avg_dxy * dxy_data$no_sites, dxy_data$pop_pair, sum) / tapply(dxy_data$no_sites, dxy_data$pop_pair, sum)
cat("\nGenome-wide average Absolute Divergence (Dxy):\n")
print(avg_dxy)

# For fst (This command should now succeed)
avg_fst <- tapply(fst_data$avg_wc_fst * fst_data$no_snps, fst_data$pop_pair, sum) / tapply(fst_data$no_snps, fst_data$pop_pair, sum)
cat("\nGenome-wide average Relative Divergence (Fst):\n")
print(avg_fst)


cols <- viridis(4, alpha = 0.6)
# violoin plot
vioplot(avg_pi ~ pop, data = pi_data,
        main = "",
        xlab = "",
        ylab = "Nucleotide Diversity (π)",
        col = cols,          # Apply the custom colors
        cex.axis = 1,
        cex.lab = 1,
        cex.main = 1,
        las = 1)
boxplot(avg_pi ~ pop, data = pi_data,
        add = TRUE,
        boxwex = 0.2,  
        col = "white", 
        border = "black",
        outline = FALSE, 
        names = NA,
        axes = FALSE)
unique(pi_data$chromosome)[1]
lemmon_1 <- pi_data[pi_data$chromosome == "JAROYF020000001.1",][pi_data[pi_data$chromosome == "JAROYF020000001.1",]$pop == 'Lemmon',]
mogollon_1 <- pi_data[pi_data$chromosome == "JAROYF020000001.1",][pi_data[pi_data$chromosome == "JAROYF020000001.1",]$pop == 'Mogollon',]
utah_1 <- pi_data[pi_data$chromosome == "JAROYF020000001.1",][pi_data[pi_data$chromosome == "JAROYF020000001.1",]$pop == 'Utah',]
portal_1 <- pi_data[pi_data$chromosome == "JAROYF020000001.1",][pi_data[pi_data$chromosome == "JAROYF020000001.1",]$pop == 'Portal',]

cols <- viridis(4, alpha = 0.3)
cex = 0.3
plot(lemmon_1$window_pos_2 + 5000 - 1, lemmon_1$avg_pi, cex = cex, col = cols[1], type = 'l')
points(mogollon_1$window_pos_2 + 5000 - 1, mogollon_1$avg_pi, cex = cex, col = cols[2], type = 'l')
points(utah_1$window_pos_2 + 5000 - 1, utah_1$avg_pi, cex = cex, col = cols[3], type = 'l')
points(portal_1$window_pos_2 + 5000 - 1, portal_1$avg_pi, cex = cex, col = cols[4], type = 'l')


#############################
# Absolute Divergence (Dxy) #
#############################
pair_order <- c("Utah-Mogollon", "Lemmon-Utah", "Portal-Utah", 
                "Lemmon-Mogollon", "Portal-Mogollon", "Lemmon-Portal")

# Convert the 'pop_pair' column into an ordered factor
dxy_data$pop_pair <- factor(dxy_data$pop_pair, levels = pair_order)

# Create the main violin plot
vioplot(avg_dxy ~ pop_pair, data = dxy_data,
        main = "Distribution of Absolute Divergence (Dxy)",
        xlab = "Population Pair",
        ylab = "Average Dxy",
        col = "lightgray", # Use a neutral color for all pairs
        las = 2) # Rotate x-axis labels for readability

# Overlay a boxplot on top
boxplot(avg_dxy ~ pop_pair, data = dxy_data,
        add = TRUE,
        boxwex = 0.2,       # Make the boxes narrow
        col = "white",
        border = "black",
        outline = FALSE,    # Hide outliers
        axes = FALSE,       # Don't re-draw axes
        names = NA)
#------------------------#
# load libraries
library(tidyverse)
library(ggridges)
library(patchwork)

# load data
dxy_data <- read.table("../Data/Dxy/pixy_dxy.txt", header=T)
fst_data <- read.table("../Data/Dxy/pixy_fst.txt", header=T)

# add some columns
dxy_data$data_type <- "Dxy"
dxy_data <- mutate(dxy_data,
                   comparison = paste(pop1, pop2, sep = '_v_'))

fst_data$data_type <- "Fst"
fst_data <- mutate(fst_data,
                   comparison = paste(pop1, pop2, sep = '_v_'))

# subset and merge dataframes
dxy_data_sub <- dxy_data %>% select(comparison,data_type,chromosome,window_pos_1,window_pos_2,avg_dxy)
colnames(dxy_data_sub) <- c("comparison","data_type","chromosome","window_pos_1","window_pos_2","value")

fst_data_sub <- fst_data %>% select(comparison,data_type,chromosome,window_pos_1,window_pos_2,avg_wc_fst)
colnames(fst_data_sub) <- c("comparison","data_type","chromosome","window_pos_1","window_pos_2","value")

# cheeky fix to get matching rows from both datasets
tmp_dxy <- semi_join(dxy_data_sub, fst_data_sub, by=c("comparison", "chromosome", "window_pos_1", "window_pos_2"))
tmp_fst <- semi_join(fst_data_sub, dxy_data_sub,  by=c("comparison", "chromosome", "window_pos_1", "window_pos_2"))

# join the datasets together to create a single dataframe
data <- full_join(tmp_dxy, tmp_fst)

# add numbering to help with plotting
data <- data %>%
  filter(!str_detect(chromosome, "^TJAROYF02")) %>%
  group_by(comparison, data_type) %>%
  mutate(position = 1:n())

# add sex chromosome information
data <- data %>%
  mutate(chr_type = ifelse(str_detect(chromosome, "^TJAROYF02"), "sexchr", "autosome"))

# get position for vertical lines used in plot to delineate the linkage groups
data %>%
  group_by(chromosome) %>%
  summarise(max = max(position, na.rm = TRUE))

# summarise median data for Fst and Dxy
data %>%
  group_by(comparison,data_type) %>%
  summarise(median = median(value, na.rm = TRUE))

# # subset data to get only fst values
data_fst <- data %>% filter(data_type=="Fst" & 
                              (comparison=="Lemmon_v_Portal" | comparison=="Lemmon_v_Utah" | comparison=="Portal_v_Mogollon" 
                               | comparison=="Lemmon_v_Mogollon" | comparison=="Portal_v_Utah" | comparison=="Utah_v_Mogollon"))

# median values for density plots
data_fst_median <- data_fst %>%
  group_by(comparison) %>%
  summarise(median = median(value, na.rm = TRUE))

# genomewide plot of fst values for each pairwise comparison
plot_fst_gw <- ggplot(data_fst, aes(position*20000, value, col=chromosome)) +
  geom_point(size=0.5) +
  scale_colour_cyclical(values = c("#3030D0", "#9090F0")) +
  facet_grid(comparison~.) +
  ylim(0,1) +
  theme_bw() +
  labs(x="Genomic Position", y="Genetic differentiation")

plot_fst_density <- ggplot(data_fst, aes(value, col=comparison, fill=comparison)) +
  geom_density(show.legend = FALSE) +
  facet_grid(comparison~.) +
  theme_bw() + xlim(0,1) +
  geom_vline(data=data_fst_median,aes(xintercept=median),linetype="dashed")+
  labs(x="Fst", y="Density") +
  theme()

plot_fst_gw + plot_fst_density +  plot_layout(widths = c(5, 1))
data_fst[data_fst$comparison == 'Utah_v_Mogollon',]$value
max(data_fst[data_fst$comparison == 'Utah_v_Mogollon',]$value, na.rm = T)

data_fst[data_fst$comparison == 'Utah_v_Mogollon',][which(data_fst[data_fst$comparison == 'Utah_v_Mogollon',]$value == 1),]
data_fst[data_fst$comparison == 'Lemmon_v_Utah',][which(data_fst[data_fst$comparison == 'Lemmon_v_Utah',]$value == 1),]
data_fst[data_fst$comparison == 'Portal_v_Utah',][which(data_fst[data_fst$comparison == 'Portal_v_Utah',]$value == 1),]

data_fst[data_fst$comparison == 'Portal_v_Utah',][which(data_fst[data_fst$window_pos_1 == 7105001,]),]
data_fst[data_fst$comparison == 'Portal_v_Utah',][data_fst[data_fst$comparison == 'Portal_v_Utah',]$window_pos_1 == 7105001,]
data_fst[data_fst$comparison == 'Portal_v_Utah',]$value[order(data_fst[data_fst$comparison == 'Portal_v_Utah',]$value,decreasing = T)]
length(data_fst[data_fst$comparison == 'Portal_v_Utah',]$value)/100

# # dxy data
# # subset data to get only dxy values
data_dxy <- data %>% filter(data_type=="Dxy" & (comparison=="Lemmon_v_Portal" | comparison=="Lemmon_v_Utah" | comparison=="Portal_v_Mogollon" 
                                                | comparison=="Lemmon_v_Mogollon" | comparison=="Portal_v_Utah" | comparison=="Utah_v_Mogollon"))

# median values for density plots
data_dxy_median <- data_dxy %>%
  group_by(comparison) %>%
  summarise(median = median(value, na.rm = TRUE))

# genomewide plot of dxy values for each pairwise comparison
plot_dxy_gw <- ggplot(data_dxy, aes(position*20000, value, col=chromosome)) +
  geom_point(size=0.5) +
  scale_colour_cyclical(values = c("#3030D0", "#9090F0")) +
  facet_grid(comparison~.) +
  ylim(0,1) +
  theme_bw() +
  labs(x="Genomic Position", y="Genetic differentiation")

plot_dxy_density <- ggplot(data_dxy, aes(value, col=comparison, fill=comparison)) +
  geom_density(show.legend = FALSE) +
  facet_grid(comparison~.) +
  theme_bw() + xlim(0,1) +
  geom_vline(data=data_dxy_median,aes(xintercept=median),linetype="dashed")+
  labs(x="Dxy", y="Density") +
  theme()

plot_dxy_gw + plot_dxy_density +  plot_layout(widths = c(5, 1))
