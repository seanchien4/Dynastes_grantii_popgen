library(viridis)
library(scales)
# par(mfrow = c(2, 2))
sim <- 100
plot(NA,
     xlim = c(0, 14),
     ylim = c(10, 14),
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
y = 10
cex = 0.8
text(13,y,"LGP",cex=cex)
text(10.5,y,"Holocene",cex=cex)


lwd = 0.5
cols <- viridis(4, alpha = 0.05)
lwd_main <- 2.5
# Mogollon
file0 <- read.csv('../Data/csv/Mo_1.csv')
x <- log(file0$x)
y_mat <- matrix(NA, nrow = length(x), ncol = sim)
for (i in 1:sim) {
  file <- paste('../Data/csv/Mo_', i, '.csv', sep = '')
  dat <- read.csv(file)
  y <- log(dat$y)
  lines(x, y, col = cols[1], lwd = lwd)
  y_mat[, i] <- y
}
mean_y <- apply(y_mat, 1, mean, na.rm = TRUE)
lines(x, mean_y, col = viridis(4)[1], lwd = lwd_main)

# Mt. Lemmon
file0 <- read.csv('../Data/csv/ML_1.csv')
x <- log(file0$x)
y_mat <- matrix(NA, nrow = length(x), ncol = sim)
for (i in 1:sim) {
  file <- paste('../Data/csv/ML_', i, '.csv', sep = '')
  dat <- read.csv(file)
  y <- log(dat$y)
  # draw individual simulation
  lines(x, y, col = cols[2], lwd = lwd)
  y_mat[, i] <- y
}
mean_y <- apply(y_mat, 1, mean, na.rm = TRUE)
lines(x, mean_y, col = viridis(4)[2], lwd = lwd_main)

# Portal
file0 <- read.csv('../Data/csv/PT_1.csv')
x <- log(file0$x)
y_mat <- matrix(NA, nrow = length(x), ncol = sim)
for (i in 1:sim) {
  file <- paste('../Data/csv/PT_', i, '.csv', sep = '')
  dat <- read.csv(file)
  y <- log(dat$y)
  lines(x, y, col = cols[3], lwd = lwd)
  y_mat[, i] <- y
}
mean_y <- apply(y_mat, 1, mean, na.rm = TRUE)
lines(x, mean_y, col = viridis(4)[3], lwd = lwd_main)

# Utah
file0 <- read.csv('../Data/csv/UT_1.csv')
x <- log(file0$x)
y_mat <- matrix(NA, nrow = length(x), ncol = sim)
for (i in 1:sim) {
  file <- paste('../Data/csv/UT_', i, '.csv', sep = '')
  dat <- read.csv(file)
  y <- log(dat$y)
  lines(x, y, col = cols[4], lwd = lwd)
  y_mat[, i] <- y
}
mean_y <- apply(y_mat, 1, mean, na.rm = TRUE)
lines(x, mean_y, col = viridis(4)[4], lwd = lwd_main)



x = 0.05
y = 10.6
s = 0.2
cols <- viridis(4)
points(x,y, col = cols[1], pch =16)
text(x,y, pos = 4, labels = 'Mongollon Rim')
points(x,y-s, col = cols[2], pch =16)
text(x,y-s, pos = 4, labels = 'Mt. Lemmon')
points(x,y-2*s, col = cols[3], pch =16)
text(x,y-2*s, pos = 4, labels = 'Portal')
points(x,y-3*s, col = cols[4], pch =16)
text(x,y-3*s, pos = 4, labels = 'Utah')

