# tajima's D 
library(dplyr)
library(viridis)
results_dir <- "../Data/tajima" 
population_names <- c("Utah", "Lemmon", "Portal", "Mogollon")
for (i in 1:length(population_names)){
  tajimaD_file_path <- file.path(results_dir, paste0(population_names[i], "_tajimaD.Tajima.D"))
  temp.name <- paste0(population_names[i],'.tajima.dat') 
  temp.dat <- read.table(tajimaD_file_path, header = TRUE, stringsAsFactors = FALSE)
  genome_wide_mean_d <- mean(temp.dat$TajimaD, na.rm = TRUE)
  genome_wide_sd_d <- sd(temp.dat$TajimaD, na.rm = TRUE)
  temp.dat <- temp.dat %>%
    mutate(z_score = (TajimaD - genome_wide_mean_d) / genome_wide_sd_d)
  temp.dat <- temp.dat[order(temp.dat$CHROM, temp.dat$BIN_START), ]
  assign(temp.name,temp.dat)
}
pseudo <- sort(unique(temp.dat$CHROM))[1:11]
unplace <- sort(unique(temp.dat$CHROM))[12:length(unique(temp.dat$CHROM))]

########
cex = 0.3
pch = 19
cols <- c("#3030D0", "#9090F0")
cols <- c(viridis(4)[1],viridis(4,alpha = .5)[1])
dat <- Mogollon.tajima.dat
bin_size <- Mogollon.tajima.dat$BIN_START[2]
dat$CHROM <- factor(dat$CHROM, levels = unique(dat$CHROM))
chr_len <- tapply(dat$BIN_START, dat$CHROM, max) + bin_size
chr_offset <- c(0, cumsum(as.numeric(chr_len))[-length(chr_len)])
names(chr_offset) <- names(chr_len)
dat$x <- dat$BIN_START + (bin_size/2) + chr_offset[as.character(dat$CHROM)]
chr_mid <- chr_offset + as.numeric(chr_len)/2
# 
# cols <- c("#3030D0", "#9090F0")
cols <- c(viridis(4)[1],viridis(4,alpha = .2)[1])
n_colored <- 11
chr_levels <- levels(dat$CHROM)
chr_colmap <- setNames(rep("grey70", length(chr_levels)), chr_levels)
chr_colmap[chr_levels[1:min(n_colored, length(chr_levels))]] <-
  rep(cols, length.out = min(n_colored, length(chr_levels)))

dat$col <- chr_colmap[as.character(dat$CHROM)]

plot(dat$x, dat$z_score, pch=pch, cex=cex,col=dat$col,
     xlab="Position", ylab="Z score", xaxt="n")
# Mt. Lemmone
dat <- Lemmon.tajima.dat
bin_size <- dat$BIN_START[2]
dat$CHROM <- factor(dat$CHROM, levels = unique(dat$CHROM))
chr_len <- tapply(dat$BIN_START, dat$CHROM, max) + bin_size
chr_offset <- c(0, cumsum(as.numeric(chr_len))[-length(chr_len)])
names(chr_offset) <- names(chr_len)
dat$x <- dat$BIN_START + (bin_size/2) + chr_offset[as.character(dat$CHROM)]
chr_mid <- chr_offset + as.numeric(chr_len)/2
# 
cols <- c(viridis(4)[2],viridis(4,alpha = .5)[2])
n_colored <- 11
chr_levels <- levels(dat$CHROM)
chr_colmap <- setNames(rep("grey70", length(chr_levels)), chr_levels)
chr_colmap[chr_levels[1:min(n_colored, length(chr_levels))]] <-
  rep(cols, length.out = min(n_colored, length(chr_levels)))
dat$col <- chr_colmap[as.character(dat$CHROM)]
points(dat$x, dat$z_score, pch=pch, cex=cex,col=dat$col)
#
# Portal
dat <- Portal.tajima.dat
bin_size <- dat$BIN_START[2]
dat$CHROM <- factor(dat$CHROM, levels = unique(dat$CHROM))
chr_len <- tapply(dat$BIN_START, dat$CHROM, max) + bin_size
chr_offset <- c(0, cumsum(as.numeric(chr_len))[-length(chr_len)])
names(chr_offset) <- names(chr_len)
dat$x <- dat$BIN_START + (bin_size/2) + chr_offset[as.character(dat$CHROM)]
chr_mid <- chr_offset + as.numeric(chr_len)/2
# 
cols <- c(viridis(4)[3],viridis(4,alpha = .5)[3])
n_colored <- 11
chr_levels <- levels(dat$CHROM)
chr_colmap <- setNames(rep("grey70", length(chr_levels)), chr_levels)
chr_colmap[chr_levels[1:min(n_colored, length(chr_levels))]] <-
  rep(cols, length.out = min(n_colored, length(chr_levels)))
dat$col <- chr_colmap[as.character(dat$CHROM)]
points(dat$x, dat$z_score, pch=pch, cex=cex,col=dat$col)

axis(1, at=chr_mid[1:11], labels=names(chr_len)[1:11], las=2, cex.axis=0.7)
abline(h=3, col = 'red')

