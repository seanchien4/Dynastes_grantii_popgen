library(viridis)
library(scales)
# par(mfrow = c(2, 2))
par(mar = c(4,4,1,1))
sim <- 100
plot(NA,
     xlim = c(0, 14),
     ylim = c(9, 13),
     xlab = "Yeasrs ago",
     ylab = expression(italic(N[e])),
     xaxt = "n",
     yaxt = "n",
)
x_ticks_orig <- c(1, 10, 50, 500, 5000, 50000, 500000)
x_ticks_log  <- log(x_ticks_orig)
axis(1, at = x_ticks_log, labels = c('1','10','50','500','5K','50K','500K'))
y_ticks_orig <- c(1e4, 1e5, 1e6)
y_ticks_log  <- log(y_ticks_orig)
axis(2, at = y_ticks_log, labels = c('10K','100K','1M'))
# 
abline(v=c(log(11700),log(115000)), col=alpha("black", 0.4))
y = 9.1
cex = 0.8
text(13,y,"LGP",cex=cex)
text(10.5,y,"Holocene",cex=cex)


lwd = 0.5
cols_BS <- viridis(4, alpha = 0.05)
cols <- viridis(4)
lwd_main <- 2.5
# Mogollon
file0 <- read.csv('~/Desktop/Projects/Dynastes_grantii_pop/Data/Ne/csv/Mo_1.csv')
x <- log(file0$x)
y_mat <- matrix(NA, nrow = length(x), ncol = sim)
for (i in 1:sim) {
  file <- paste('~/Desktop/Projects/Dynastes_grantii_pop/Data/Ne/csv/Mo_', i, '.csv', sep = '')
  dat <- read.csv(file)
  y <- log(dat$y)
  lines(x, y, col = cols_BS[1], lwd = lwd)
  y_mat[, i] <- y
}
mean_y <- apply(y_mat, 1, mean, na.rm = TRUE)
lines(x, mean_y, col = cols[1], lwd = lwd_main)

# Mt. Lemmon
file0 <- read.csv('~/Desktop/Projects/Dynastes_grantii_pop/Data/Ne/csv/ML_1.csv')
x <- log(file0$x)
y_mat <- matrix(NA, nrow = length(x), ncol = sim)
for (i in 1:sim) {
  file <- paste('~/Desktop/Projects/Dynastes_grantii_pop/Data/Ne/csv/ML_', i, '.csv', sep = '')
  dat <- read.csv(file)
  y <- log(dat$y)
  # draw individual simulation
  lines(x, y, col = cols_BS[2], lwd = lwd)
  y_mat[, i] <- y
}
mean_y <- apply(y_mat, 1, mean, na.rm = TRUE)
lines(x, mean_y, col = cols[2], lwd = lwd_main)

# Portal
file0 <- read.csv('~/Desktop/Projects/Dynastes_grantii_pop/Data/Ne/csv/PT_1.csv')
x <- log(file0$x)
y_mat <- matrix(NA, nrow = length(x), ncol = sim)
for (i in 1:sim) {
  file <- paste('~/Desktop/Projects/Dynastes_grantii_pop/Data/Ne/csv/PT_', i, '.csv', sep = '')
  dat <- read.csv(file)
  y <- log(dat$y)
  lines(x, y, col = cols_BS[3], lwd = lwd)
  y_mat[, i] <- y
}
mean_y <- apply(y_mat, 1, mean, na.rm = TRUE)
lines(x, mean_y, col = cols[3], lwd = lwd_main)

# Utah
file0 <- read.csv('~/Desktop/Projects/Dynastes_grantii_pop/Data/Ne/csv/UT_1.csv')
x <- log(file0$x)
y_mat <- matrix(NA, nrow = length(x), ncol = sim)
for (i in 1:sim) {
  file <- paste('~/Desktop/Projects/Dynastes_grantii_pop/Data/Ne/csv/UT_', i, '.csv', sep = '')
  dat <- read.csv(file)
  y <- log(dat$y)
  lines(x, y, col = cols_BS[4], lwd = lwd)
  y_mat[, i] <- y
}
mean_y <- apply(y_mat, 1, mean, na.rm = TRUE)
lines(x, mean_y, col = cols[4], lwd = lwd_main)

# Unknown
# file0 <- read.csv('../Data/Ne/csv/UN_1.csv')
# x <- log(file0$x)
# y_mat <- matrix(NA, nrow = length(x), ncol = sim)
# for (i in 1:sim) {
#   file <- paste('../Data/Ne/csv/UN_', i, '.csv', sep = '')
#   dat <- read.csv(file)
#   y <- log(dat$y)
#   lines(x, y, col = cols_BS[5], lwd = lwd)
#   y_mat[, i] <- y
# }
# mean_y <- apply(y_mat, 1, mean, na.rm = TRUE)
# lines(x, mean_y, col = cols[5], lwd = lwd_main)


x = 0.05
y = 9.6
s = 0.15
points(x,y, col = cols[1], pch =16)
text(x,y, pos = 4, labels = 'Mongollon Rim')
points(x,y-s, col = cols[2], pch =16)
text(x,y-s, pos = 4, labels = 'Mt. Lemmon')
points(x,y-2*s, col = cols[3], pch =16)
text(x,y-2*s, pos = 4, labels = 'Chiricahua Mt.')
points(x,y-3*s, col = cols[4], pch =16)
text(x,y-3*s, pos = 4, labels = 'Utah')
# points(x,y-4*s, col = cols[5], pch =16)
# text(x,y-4*s, pos = 4, labels = 'Unknown')
