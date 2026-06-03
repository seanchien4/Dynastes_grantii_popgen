library(tidyverse)
library(viridis)
# ============================================================
# LOAD ALL POPULATIONS
# ============================================================

load_roh_data <- function(file_path, pop_name) {
  data <- read.table(file_path)
  data <- data[, c(2:6)]
  colnames(data) <- c("Sample", "Chrom", "Start", "End", "Length")
  data$Length_Mb <- data$Length / 1e6
  data$Population <- pop_name
  return(data)
}

# Load all four populations
mogollon <- load_roh_data("~/Desktop/Projects/Dynastes_grantii_pop/Data/ROH/Mogollon_roh_regions.txt", "Mogollon")
lemmon <- load_roh_data("~/Desktop/Projects/Dynastes_grantii_pop/Data/ROH/Lemmon_roh_regions.txt", "Mt. Lemmon")
lemmon$Population <- 'Lemmon'
portal <- load_roh_data("~/Desktop/Projects/Dynastes_grantii_pop/Data/ROH/Portal_roh_regions.txt", "Portal")
portal$Population <- 'Chiricahua'
utah <- load_roh_data("~/Desktop/Projects/Dynastes_grantii_pop/Data/ROH/Utah_roh_regions.txt", "Utah")

# Combine
roh_all <- rbind(mogollon, lemmon, portal, utah)

# ============================================================
# DEFINE POPULATION ORDER AND COLORS
# ============================================================
cols <- viridis(4)
pop_order  <- c("Utah", "Chiricahua", "Lemmon", "Mogollon")  
pop_colors <- c("Utah"     = cols[4],
                "Chiricahua"   = cols[2], 
                "Lemmon"   = cols[3], 
                "Mogollon" = cols[1]) 
cols_light <- viridis(4,alpha = 0.2)
pop_colors_light <- c("Utah"     = cols_light[4],   
                      "Chiricahua"   = cols_light[2], 
                      "Lemmon"   = cols_light[3], 
                      "Mogollon" = cols_light[1])  

roh_all$Population <- factor(roh_all$Population, levels = pop_order)



# ============================================================
# ROH LENGTH CATEGORIES
# ============================================================

roh_all <- roh_all %>%
  mutate(ROH_class = case_when(
    Length_Mb < 0.1               ~ "< 100 kb",
    Length_Mb >= 0.1 & Length_Mb < 1 ~ "100 kb – 1 Mb",
    Length_Mb >= 1   & Length_Mb < 5 ~ "1 – 5 Mb",
    Length_Mb >= 5                ~ "> 5 Mb"
  ),
  ROH_class = factor(ROH_class, 
                     levels = c("< 100 kb", 
                                "100 kb – 1 Mb",
                                "1 – 5 Mb", 
                                "> 5 Mb")))
# ============================================================
# PER-INDIVIDUAL STATISTICS
# ============================================================

# You need genome size for FROH
# Replace with your actual autosomal genome size
genome_size <- 6.787e8
ind_summary <- roh_all %>%
  group_by(Population, Sample) %>%
  summarise(
    NROH    = n(),                          # number of ROH
    SROH    = sum(Length),                  # total ROH length (bp)
    SROH_Mb = sum(Length_Mb),              # total ROH in Mb
    FROH    = sum(Length) / genome_size,   # inbreeding coefficient
    Mean_ROH_length = mean(Length_Mb),
    Max_ROH_length  = max(Length_Mb),
    .groups = "drop"
  )

# Quick check
ind_summary %>%
  group_by(Population) %>%
  summarise(
    n_ind       = n(),
    mean_FROH   = mean(FROH),
    sd_FROH     = sd(FROH),
    mean_NROH   = mean(NROH),
    mean_SROH   = mean(SROH_Mb)
  )
# ==============================================
# ROH Length Distribution (histogram/density)
# ==============================================
pops <- unique(roh_all$Population)
par(
  mfrow = c(2, 2),   
  mar = c(4, 4, 2, 1)
)

# Loop through populations
labels <- c("(A)", "(B)", "(C)", "(D)")
i <- 1
for (pop in pops) {
  
  # Subset data
  dat <- roh_all$Length_Mb[roh_all$Population == pop]
  
  # Log-transform (equivalent to scale_x_log10)
  log_dat <- log10(dat)
  
  # Draw histogram
  hist(
    log_dat,
    breaks = 80,
    col = pop_colors[pop],
    border = "white",
    main = '',
    xlab = "",
    ylab = "Count",
    xaxt = "n",
    xlim = c(-4,0.1)
  )
  # Custom x-axis to match your breaks
  ticks <- c(0.01, 0.1, 1, 10)   # reduced number of ticks
  axis(
    1,
    at = log10(ticks),
    labels = ticks,
    cex.axis = 0.8
  )
  mtext("ROH Length (Mb, log scale)", side = 1, line = 2.5, cex = 0.8)
  title(
    main = labels[i],
    adj = 0,      # <-- left
    line = 0.5,
    cex.main = 1
  )
  
  i <- i + 1
}


# =======================================
# ROH BURDEN (SROH) per individual
# =======================================
cex.axis = 0.7
par(
  mfrow = c(2, 2),   
  mar = c(4, 4, 2, 1)
)
# par(mfrow = c(1, 1), mar = c(4, 4.5, 1, 1))
ind_summary$Population <- as.factor(ind_summary$Population)
# Basic boxplot (no outliers plotted)
boxplot(
  SROH_Mb ~ Population,
  data = ind_summary,
  col = pop_colors_light[levels(ind_summary$Population)],
  border = "black",
  outline = FALSE,   # removes outliers like ggplot version
  ylab = "",
  xlab = "",
  main = "",
  xaxt = "n"
)
title(ylab = expression("S"["ROH"]*" (Mb)"), line = 2.5)
title(xlab = "Population", line = 2.5)
title(main = '(A)', adj = 0, line = 0.5,
      cex.main = 1
)
# Add jittered points (beeswarm-like)
for (i in seq_along(levels(ind_summary$Population))) {
  
  pop <- levels(ind_summary$Population)[i]
  
  y <- ind_summary$SROH_Mb[ind_summary$Population == pop]
  
  # jitter x positions around box center
  x <- jitter(rep(i, length(y)), amount = 0.15)
  
  points(
    x, y,
    bg = pop_colors[pop],
    pch = 21,
    cex = 1
  )
}
axis(1, at = 1:4, labels = c("Utah", "Chiricahua", "Lemmon", "Mogollon"),cex.axis = cex.axis)
# =====================
# FROH per individual
# =====================
ind_summary$Population <- as.factor(ind_summary$Population)

# Boxplot (no default x-axis so you can customize if needed)
boxplot(
  FROH ~ Population,
  data = ind_summary,
  col = pop_colors_light[levels(ind_summary$Population)],
  border = "black",
  outline = FALSE,
  xlab = "",
  ylab = "",
  main = '',
  xaxt = "n"
)
# Add jittered points (beeswarm-like)
for (i in seq_along(levels(ind_summary$Population))) {
  
  pop <- levels(ind_summary$Population)[i]
  
  y <- ind_summary$FROH[ind_summary$Population == pop]
  
  x <- jitter(rep(i, length(y)), amount = 0.15)
  
  points(
    x, y,
    bg = pop_colors[pop],
    pch = 21,
    cex = 1
  )
}

axis(1, at = 1:4, labels = c("Utah", "Chiricahua", "Lemmon", "Mogollon"),cex.axis = cex.axis)
title(xlab = "Population", line = 2.5)
title(ylab = expression(F[ROH]), line = 2.5)
title(main = '(B)', adj = 0, line = 0.5,
      cex.main = 1)
# ============================================================
# PLOT 5: NROH vs SROH (scatter) — classic ROH plot
# ============================================================
ind_summary$Population <- as.factor(ind_summary$Population)

# Get levels and colors
pops <- levels(ind_summary$Population)

# Base plot
plot(
  ind_summary$NROH,
  ind_summary$SROH_Mb,
  type = "n",   # empty plot first
  xlab = "",
  ylab = "",
  main = '',
  xaxt = 'n'
  
)
axis(1, at = seq(from = 6000, to = 14000, by = 2000), labels = seq(from = 6000, to = 14000, by = 2000),cex.axis = cex.axis)
title(xlab = expression("N"["ROH"]*""), line = 2.5)
title(ylab = expression("S"["ROH"]*"(Mb)"), line = 2.5)
title(main = '(C)', adj = 0, line = 0.5,
      cex.main = 1)
# Add points by population
for (pop in pops) {
  
  idx <- ind_summary$Population == pop
  
  points(
    ind_summary$NROH[idx],
    ind_summary$SROH_Mb[idx],
    bg = pop_colors[pop],
    pch = 21,   # same shape for all (like ggplot default dot)
    cex = 1
  )
}
s = 6.5
cex = 0.8
points(4500, 125, pch = 21, bg = pop_colors[4], cex = cex )
text(4500, 125, labels = 'Mogollon Rim', pos = 4)
points(4500, 125-s, pch = 21, bg = pop_colors[2], cex = cex )
text(4500, 125-s, labels = 'Mt. Lemmon', pos = 4)
points(4500, 125-2*s, pch = 21, bg = pop_colors[3], cex = cex )
text(4500, 125-2*s, labels = 'Chiricahua Mt.', pos = 4)
points(4500, 125-3*s, pch = 21, bg = pop_colors[1], cex = cex )
text(4500, 125-3*s, labels = 'Utah', pos = 4)
# ============================================================
# PLOT 6: SROH broken down by ROH class per individual
# This shows ancient vs recent inbreeding contribution
# ============================================================
sroh_by_class <- roh_all %>% group_by(Population, Sample, ROH_class) %>%
  summarise(SROH_class = sum(Length_Mb), .groups = "drop")
sroh_by_class$Population <- as.factor(sroh_by_class$Population)
sroh_by_class$Sample <- as.factor(sroh_by_class$Sample)
sroh_by_class$ROH_class <- as.factor(sroh_by_class$ROH_class)

pops <- unique(sroh_by_class$Population)
classes <- levels(sroh_by_class$ROH_class)

cols <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdYlBu")))(length(classes))

# GLOBAL Y LIMIT
ymax <- max(sroh_by_class$SROH_class, na.rm = TRUE)

par(
  mfrow = c(2, 2),
  mar = c(6, 4, 3, 1)
)

for (pop in pops) {
  
  dat_pop <- sroh_by_class[sroh_by_class$Population == pop, ]
  samples <- unique(dat_pop$Sample)
  
  mat <- matrix(0, nrow = length(classes), ncol = length(samples))
  rownames(mat) <- classes
  colnames(mat) <- samples
  
  for (i in seq_along(samples)) {
    
    tmp <- dat_pop[dat_pop$Sample == samples[i], ]
    mat[match(tmp$ROH_class, classes), i] <- tmp$SROH_class
  }
  
  barplot(
    mat,
    col = cols,
    border = NA,
    space = 0,
    las = 2,
    ylim = c(0, ymax),   # <<< KEY FIX
    main = pop,
    ylab = "Total ROH Length (Mb)",
    cex.names = 0.6
  )
}
# ============================================================
# STATISTICAL COMPARISONS
# ============================================================
pop_summary <- ind_summary %>%
  group_by(Population) %>%
  summarise(
    n           = n(),
    # SROH
    mean_SROH   = round(mean(SROH_Mb), 2),
    median_SROH = round(median(SROH_Mb), 2),
    sd_SROH     = round(sd(SROH_Mb), 2),
    min_SROH    = round(min(SROH_Mb), 2),
    max_SROH    = round(max(SROH_Mb), 2),
    # FROH
    mean_FROH   = round(mean(FROH), 4),
    median_FROH = round(median(FROH), 4),
    sd_FROH     = round(sd(FROH), 4),
    # NROH
    mean_NROH   = round(mean(NROH), 1),
    median_NROH = round(median(NROH), 1),
    sd_NROH     = round(sd(NROH), 1),
    .groups = "drop"
  )

print(pop_summary)

# Also get per individual so you can report outliers
print(ind_summary %>% 
        arrange(Population, desc(SROH_Mb)) %>%
        select(Population, Sample, SROH_Mb, FROH, NROH))
# Kruskal-Wallis (non-parametric, does not assume normality)
kruskal.test(SROH_Mb ~ Population, data = ind_summary)
kruskal.test(FROH    ~ Population, data = ind_summary)
kruskal.test(NROH    ~ Population, data = ind_summary)

# Pairwise Wilcoxon with correction
pairwise.wilcox.test(ind_summary$SROH_Mb, 
                     ind_summary$Population,
                     p.adjust.method = "bonferroni")

pairwise.wilcox.test(ind_summary$FROH, 
                     ind_summary$Population,
                     p.adjust.method = "bonferroni")

