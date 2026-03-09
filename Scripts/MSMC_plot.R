# msmc2 
# Sean Chien
mu <- 3e-9
gen <- 2

# read data 
DgMo <- read.table("../Data/DgMo.msmc2.final.txt", header=TRUE)
DgMo.loop <- read.table('../Data/DgMo.msmc2.loop')
colnames(DgMo.loop) <- c('time','pop')

DgUT<-read.table("../Data/DgUT.msmc2.final.txt", header=TRUE)
DgUT.loop <- read.table('../Data/DgUT.msmc2.loop')
colnames(DgUT.loop) <- c('time','pop')

DgSky <- read.table("../Data/Dgmiddle.msmc2.final.txt", header=TRUE)
DgSky.loop <- read.table('../Data/Dgmiddle.msmc2.loop')
colnames(DgSky.loop) <- c('time','pop')

DgSR <- read.table("../Data/DgSR.msmc2.final.txt", header=TRUE)
DgSR.loop <- read.table('../Data/DgSR.msmc2.loop')
colnames(DgSR.loop) <- c('time','pop')

# plot
plot(DgMo$left_time_boundary/mu*gen, (1/DgMo$lambda)/(2*mu), 
     log="x",
     xlim = c(2000,400000),
     ylim=c(0,200000),
     type="n", xlab="Years ago", ylab=expression(italic(N[e])),
     xaxt = "n") 

axis(1, at = c(1e3,1e4, 1e5, 1e6, 1e7), labels = c("1,000","10,000", "100,000", "1,000,000", "10,000,000"))
# axis(2, at = c(1e3,1e4,1e5), labels = expression(10^3,10^4,10^5))
axis(1, at = c(2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000), labels = FALSE, tcl = -0.3)
axis(1, at = c(20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000), labels = FALSE, tcl = -0.3)
axis(1, at = c(200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000), labels = FALSE, tcl = -0.3)
# axis(2, las = 2)  # vertical labels for y-axis
# Pleistocene time frame
rect(xleft = 11700, xright = 2600000,
     ybottom = par("usr")[3],  # bottom of y-axis
     ytop = par("usr")[4],     # top of y-axis
     col = rgb(0.9, 0.9, 0.9, 0.3), border = NA)
abline(v=11700, col = 'gray')
abline(v=2600000, col = 'gray')

# iteration 
lwd = 0.8
for (i in 0:19) {
  start <- i * 32 + 1
  end <- start + 31
  chunk <- DgMo.loop[start:end, ]
  lines(chunk[[1]]/mu*gen, (1/chunk[[2]])/(2*mu), type="s", col="#fee0d2",lwd = lwd )
}

for (i in 0:19) {
  start <- i * 32 + 1
  end <- start + 31
  chunk <- DgUT.loop[start:end, ]
  lines(chunk[[1]]/mu*gen, (1/chunk[[2]])/(2*mu), type="s", col="#deebf7",lwd = lwd )
}

for (i in 0:19) {
  start <- i * 32 + 1
  end <- start + 31
  chunk <- DgSky.loop[start:end, ]
  lines(chunk[[1]]/mu*gen, (1/chunk[[2]])/(2*mu), type="s", col="#e5f5e0",lwd = lwd )
}

for (i in 0:19) {
  start <- i * 32 + 1
  end <- start + 31
  chunk <- DgSR.loop[start:end, ]
  lines(chunk[[1]]/mu*gen, (1/chunk[[2]])/(2*mu), type="s", col="#efedf5",lwd = lwd )
}

lwd = 1.5
lines(DgMo$left_time_boundary/mu*gen, (1/DgMo$lambda)/(2*mu), type="s", col="#ef3b2c", lwd = lwd)
lines(DgUT$left_time_boundary/mu*gen, (1/DgUT$lambda)/(2*mu), type="s", col="#4292c6", lwd = lwd)
lines(DgSky$left_time_boundary/mu*gen, (1/DgSky$lambda)/(2*mu), type="s", col="#41ab5d", lwd = lwd)
lines(DgSR$left_time_boundary/mu*gen, (1/DgSR$lambda)/(2*mu), type="s", col="#807dba", lwd = lwd)

y = 200000
x = 2000
s = 10000
points(x,y, col="#ef3b2c", pch =16)
text(x=x, y, labels = "Mogollon Rim", pos = 4)
points(x,y-1*s, col="#4292c6", pch =16)
text(x=x, y-1*s, labels = "Utah", pos = 4)
points(x,y-2*s, col="#41ab5d", pch =16)
text(x=x, y-2*s, labels = "Sky Islands", pos = 4)
points(x,y-3*s, col="#807dba", pch =16)
text(x=x, y-3*s, labels = "Mexico + Unknown", pos = 4)


text(x=4500, 0, labels = "Holocene")
text(x=70000, 0, labels = "Pleistocene")

