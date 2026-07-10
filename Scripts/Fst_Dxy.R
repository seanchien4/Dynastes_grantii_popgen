library(dplyr)
library(viridis)
#######
# Dxy #
#######
Dxy_data_raw <- read.table("~/Desktop/Projects/Dynastes_grantii_pop/Data/Dxy/pixy_dxy.txt", header = TRUE, stringsAsFactors = FALSE)
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

dxy_list <- lapply(mget(com.list), function(df) df$avg_dxy)
names(dxy_list) <- com.list
ord_names <- c("PM", "LM", "LP", "UM", "PU", "LU")
dxy_list <- dxy_list[ord_names]
names(dxy_list) <- c('Mogollon-Chiricahua', 'Mogollon-Lemmon','Chiricahua-Lemmon',
                     'Utah-Mogollon','Utah-Chiricahua','Utah-Lemmon')

temp_cols <- viridis(2, begin = 0.5)
cols <- ifelse(grepl("U", ord_names), temp_cols[2], temp_cols[1])
# 2. Setup densities and scaling
densities <- lapply(dxy_list, density, na.rm = TRUE)
max_dens <- max(sapply(densities, function(d) max(d$y)))
n <- length(dxy_list)
# 3. Create the empty plot frame
# par(mar = c(5, 4, 1, 1))
par(mfrow = c(1, 2), mar = c(6, 4, 1.5, 0.5))
plot(1, type = "n", 
     xlim = c(0.5, n + 0.5), ylim = c(0,1.1),
     xaxt = "n", 
     xlab = "", 
     ylab = expression(italic(D)[XY]),
     main = "")

text(x = 1:n, 
     y = par("usr")[3] - 0.03,
     labels = names(dxy_list), 
     srt = 45,       
     adj = 1,        
     cex = 0.8,
     xpd = TRUE)


# 5. Loop to draw each violin and median
for(i in 1:n) {
  # Get the density coordinates for this specific pair
  d <- densities[[i]]
  
  # Calculate the width of the 'violin' belly
  # 0.4 ensures violins don't overlap into next lane
  x_offset <- (d$y / max_dens) * 0.4
  
  # Draw the violin using polygon()
  polygon(c(i - x_offset, rev(i + x_offset)), 
          c(d$x, rev(d$x)), 
          border = "black",col = cols[i])
  
  # Draw the white median point
  points(i, median(dxy_list[[i]], na.rm = TRUE), 
         pch = 21, 
         bg = "white", 
         cex = 1)
}

points(0.5,1.1, pch = 16, col = cols[1])
text(0.5, 1.1, pos = 4, label = 'Non-Utah involved')
points(0.5,1.05, pch = 16, col = cols[4])
text(0.5, 1.05, pos = 4, label = 'Utah involved')
title(main = "(A)", adj = 0, line = 0.5)
# statistic 
# Does any pair differ from the others?
# 1. Get all unique combinations of your 6 pairs (15 combinations)
pair_names <- names(dxy_list)
combos <- combn(pair_names, 2)

# 2. Create a function to run the tests
run_comparison <- function(p1, p2) {
  test <- wilcox.test(dxy_list[[p1]], dxy_list[[p2]])
  med_diff <- median(dxy_list[[p1]], na.rm=TRUE) - median(dxy_list[[p2]], na.rm=TRUE)
  return(c(p_val = test$p.value, diff = med_diff))
}

# 3. Apply the function across all combinations
results_raw <- apply(combos, 2, function(x) run_comparison(x[1], x[2]))

# 4. Format into a clean table
full_comp <- data.frame(
  Comparison = apply(combos, 2, paste, collapse = " vs "),
  P_Value = results_raw["p_val", ],
  Median_Diff = results_raw["diff", ]
)

# 5. Apply FDR (Benjamini-Hochberg)
full_comp$FDR <- p.adjust(full_comp$P_Value, method = "BH")

# 6. Apply Biological Threshold (0.01)
threshold <- 0.01
full_comp$Significant <- full_comp$FDR < 0.05 & abs(full_comp$Median_Diff) > threshold

# 7. Print sorted by difference
full_comp <- full_comp[order(abs(full_comp$Median_Diff), decreasing = TRUE), ]
print(full_comp)

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


#######
# Fst #
#######
fst_data_raw <- read.table("~/Desktop/Projects/Dynastes_grantii_pop/Data/Dxy/pixy_fst.txt", header = TRUE, stringsAsFactors = FALSE)
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

fst_list <- lapply(mget(com.list), function(df) df$avg_wc_fst)
names(fst_list) <- com.list



# 1. Reorder the list so Utah pairs are grouped together
# This makes the "step up" in Fst visible
ord_names <- c("PM", "LM", "LP", "UM", "PU", "LU")
plot_data_fst <- fst_list[ord_names]
no_utah_med <- median(unlist(fst_list[c("PM", "LM", "LP")]), na.rm=TRUE)
names(fst_list) <- c('Mogollon-Chiricahua', 'Mogollon-Lemmon','Chiricahua-Lemmon',
                     'Utah-Mogollon','Utah-Chiricahua','Utah-Lemmon')

temp_cols <- viridis(2, begin = 0.5)
cols <- ifelse(grepl("U", ord_names), temp_cols[2], temp_cols[1])
# 2. Setup Plot
densities <- lapply(plot_data_fst, density, na.rm = TRUE)
max_dens <- max(sapply(densities, function(d) max(d$y)))
n <- length(plot_data_fst)
plot(1, type = "n", xlim = c(0.5, n + 0.5), ylim = c(-0.2, 1.1),
     xaxt = "n", xlab = "", ylab = expression(italic(F)[ST]))
text(x = 1:n, 
     y = par("usr")[3] - 0.03,
     labels = names(fst_list), 
     srt = 45,       
     adj = 1,        
     cex = 0.8,
     xpd = TRUE)

# 3. Add a dashed baseline at the median of non-Utah pairs 

# abline(h = no_utah_med, col = "gray", lty = 2, lwd = 1.5)

# 4. Loop to draw
for(i in 1:n) {
  d <- densities[[i]]
  x_offset <- (d$y / max_dens) * 0.4
  polygon(c(i - x_offset, rev(i + x_offset)), c(d$x, rev(d$x)), 
          col = cols[i], border = "black")
  points(i, median(plot_data_fst[[i]], na.rm = TRUE), pch = 21, bg = "white")
  
  # 5. ADD TEXT: Label the Median Difference for Utah pairs
  if(grepl("U", ord_names[i])) {
    text(i, 1.1, labels = "***", cex = 1)
  }
}
title(main = "(B)", adj = 0, line = 0.5)
# points(0.5,1.1, pch = 16, col = cols[1])
# text(0.5, 1.1, pos = 4, label = 'Non-Utah involved')
# points(0.5,1.04, pch = 16, col = cols[4])
# text(0.5, 1.04, pos = 4, label = 'Utah involved')
# statistic
# 1. Get all unique combinations of your 6 Fst pairs (15 total)
pair_names_fst <- names(fst_list)
combos_fst <- combn(pair_names_fst, 2)

# 2. Function to run Wilcoxon test and calculate median difference
run_fst_comparison <- function(p1, p2) {
  test <- wilcox.test(fst_list[[p1]], fst_list[[p2]])
  med_diff <- median(fst_list[[p1]], na.rm=TRUE) - median(fst_list[[p2]], na.rm=TRUE)
  return(c(p_val = test$p.value, diff = med_diff))
}

# 3. Apply the function
results_fst_raw <- apply(combos_fst, 2, function(x) run_fst_comparison(x[1], x[2]))

# 4. Create the summary table
fst_comp_results <- data.frame(
  Comparison = apply(combos_fst, 2, paste, collapse = " vs "),
  P_Value = results_fst_raw["p_val", ],
  Median_Diff = results_fst_raw["diff", ]
)

# 5. Apply FDR (Benjamini-Hochberg)
fst_comp_results$FDR <- p.adjust(fst_comp_results$P_Value, method = "BH")

# 6. Apply Biological Threshold (0.01)
# You may even consider a higher threshold for Fst (like 0.05) if drift is strong
threshold <- 0.01 
fst_comp_results$Significant <- fst_comp_results$FDR < 0.05 & abs(fst_comp_results$Median_Diff) > threshold

# 7. Print sorted by the magnitude of the difference
fst_comp_results <- fst_comp_results[order(abs(fst_comp_results$Median_Diff), decreasing = TRUE), ]
print(fst_comp_results)


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
