# male and female read depth distribution
male_density <- read.table('../Data/sexchrome/male_rd_density.txt')
plot(male_density$V1, male_density$V2, type = 'l',
     xlim = c(0,50),
     xlab = 'Read depth', ylab = 'Density')
abline(v=male_density$V1[which(male_density$V2 == max(male_density$V2[10:30]))], col = 'red')
text(42,5.2e+07, labels = 'Mode read depth = 16')

female_density <- read.table('../Data/sexchrome/female_rd_density.txt')
plot(female_density$V1, female_density$V2, type = 'l',
    xlim = c(0,50),
    xlab = 'Read depth', ylab = 'Density')
abline(v=female_density$V1[which(female_density$V2 == max(female_density$V2[10:30]))], col = 'red')
text(42,5.3e+07, labels = 'Mode read depth = 17')

files <- list.files("../Data/sexchrome/scaffold_rd_dist_5/", full.names = TRUE)

par(
  mfrow = c(4, 4),
  mar = c(2, 2, 1.5, 0.5)
)

for (f in files) {
  dat <- read.table(f)
  
  # select y values where x is between 0.5 and 1.5
  y_max <- max(dat$V2[dat$V1 >= 0.5 & dat$V1 <= 1.5], na.rm = TRUE)
  
  plot(dat$V1, dat$V2,
       type = "l",
       xlim = c(0, 2.5),
       ylim = c(0, y_max),
       main =strsplit(basename(f), split = '_')[[1]][1])
}

