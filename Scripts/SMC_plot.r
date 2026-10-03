library(viridis)
library(scales)
# par(mar = c(4,4,2,1))
# par(mfrow = c(2, 2), 
#     mar = c(3, 3, 1.5, 1),
#     mgp   = c(2, 0.5, 0))
sim <- 100
plot(NA,
     xlim = c(1, 14),
     ylim = c(9, 14),
     xlab = "",
     ylab = '',
     xaxt = "n",
     yaxt = "n"
)
# title(main = "(A)", adj = 0, line = 0.5)
# 5.8
# plot(NA,
#      xlim = c(1, 14),
#      ylim = c(9, 13),
#      xlab = "",
#      ylab = '',
#      xaxt = "n",
#      yaxt = "n",
# )
line = 2.5
title(ylab = expression(italic(N[e])), line = line)
title(xlab = "Yeasrs ago", line = line)
# for small plot 
cex.axis = 1
x_ticks_orig <- c(1, 10, 50, 500, 5000, 50000, 500000)
x_ticks_log  <- log(x_ticks_orig)
axis(1, at = x_ticks_log, labels = c('1','10','50','500','5K','50K','500K'), cex.axis = cex.axis)
y_ticks_orig <- c(1e4, 1e5, 1e6)
y_ticks_log  <- log(y_ticks_orig)
axis(2, at = y_ticks_log, labels = c('10K','100K','1M'), cex.axis = cex.axis)
abline(v=c(log(11700),log(115000)), col=alpha("black", 0.4))
y = 9
# cex = 0.8
# text(13,y,"LGP",cex=cex)
# text(10.5,y,"Holocene",cex=cex)
cex = 1
# cex = 0.6
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

# for larger plot
x = 1
y = 9.5
s = 0.18
# for 4 plots in onc
# y = 9.9
# s = 0.29

pch = 21
points(x,y, bg = cols[1], pch = pch)
text(x,y, pos = 4, labels = 'Mogollon Rim')
points(x,y-s, bg = cols[2], pch = pch)
text(x,y-s, pos = 4, labels = 'Mt. Lemmon')
points(x,y-2*s, bg = cols[3], pch = pch)
text(x,y-2*s, pos = 4, labels = 'Chiricahua Mt.')
points(x,y-3*s, bg = cols[4], pch = pch)
text(x,y-3*s, pos = 4, labels = 'Utah')
