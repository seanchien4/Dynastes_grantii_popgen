# GONE2 Ne
library(RColorBrewer)
library(scico)
library(viridis)
cols <- scico(4,palette = 'roma')
cols <- viridis(4)
dat <- read.table('../Data/GONE_Ne/GONE2_Mogollon.ped_GONE2_Ne')
colnames(dat) <- c(dat[1,1], dat[1,2])
dat <- dat[-1,]
plot(dat$Generation, as.numeric(dat$Ne_diploids), type = 'l', col = cols[1], lwd = 2.5,
     xaxt = "n", yaxt = "n",
     ylab = expression(italic(N[e])),
     xlab = "Years ago",
     ylim = c(-100,40000)
)
g <- seq(from = 0, to = 150, by = 20)
axis(side = 1, at = g, labels = g*2)
ticks_y  <- c(1000,1e+04, 2e+04, 3e+04, 4e+04)
labels_y <- c("1k","10K", "20K", "30K", "40K")
axis(2, at = ticks_y, labels = labels_y, las = 1)

dat <- read.table('../Data/GONE_Ne/GONE2_Lemmon.ped_GONE2_Ne', header = TRUE)
lines(dat$Generation,as.numeric(dat$Ne_diploids),type = 'l', col = cols[2], lwd = 2.5)
dat <- read.table('../Data/GONE_Ne/GONE2_Portal.ped_GONE2_Ne', header = TRUE)
lines(dat$Generation,as.numeric(dat$Ne_diploids),type = 'l', col = cols[3], lwd = 2.5)
dat <- read.table('../Data/GONE_Ne/GONE2_Utah.ped_GONE2_Ne', header = TRUE)
lines(dat$Generation,as.numeric(dat$Ne_diploids),type = 'l', col = cols[4], lwd = 2.5)

s=2000
labels <- c('Mogollon','Mt. Lemmon','Portal','Utah')
for (i in 1:4){
  x = 0
  y = 42000-i*s
  points(x,y, col = cols[i], pch = 16 )
  text(x,y, labels[i], pos = 4)
}