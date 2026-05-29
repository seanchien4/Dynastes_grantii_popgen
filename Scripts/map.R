library(maps)
library(raster)
library(geodata)
library(terra)
library(paletteer)
library(viridis)
library(RColorBrewer)
# inset map
map('state', fill = F, xlim = c(-115, -106), ylim = c(29.5, 38), xlab = "lon", ylab = "lat")
map.axes(cex.axis=0.8)
map.scale(x=-114.5, y=30.5, ratio=FALSE, relwidth=0.2, cex = 0.5)
# Inmap
# par(usr=c(-125, 9.8, 25, 150))
par(usr=c(-216, -63, 22, 144))
rect(xleft =-126.2,ybottom = 23.8,xright = -65.5,ytop = 50.6,col = "white")
map("usa", xlim=c(-126.2,-65.5), ylim=c(23.8,50.6),add=T)
map("state", add = T)
map("state", region = "arizona", fill = T, add = T)
### elevation map
# download map data
usa.dem1 <- elevation_3s(lon = -112, lat = 31,path=tempdir())
usa.dem2 <- elevation_3s(lon = -110, lat = 31,path=tempdir())
usa.dem3 <- elevation_3s(lon = -112, lat = 38,path=tempdir())
usa.dem4 <- elevation_3s(lon = -110, lat = 38,path=tempdir())
usa.dem5 <- elevation_3s(lon = -114, lat = 31,path=tempdir())
usa.dem6 <- elevation_3s(lon = -114, lat = 38,path=tempdir())
merge <- mosaic(usa.dem1, usa.dem2, usa.dem3, usa.dem4, usa.dem5, usa.dem6)

#############################
### sky islands elevation ###
#############################
# pdf("../plots/skyislands_elevationL.pdf", width = 6, height = 6)
par(mar = c(0, 0, 0, 0))
plot(merge, xlim = c(-115, -106.5), ylim = c(30, 38), col = terrain.colors(4000, rev = T))
map('state', fill = FALSE, xlim = c(-115, -106), ylim = c(29.5, 38.5), xlab = "lon", ylab = "lat", add= T)
map.scale(x=-114.5, y=30.4, ratio=FALSE, relwidth=0.2, cex = 0.5)
# Inmap
par(usr=c(-280, -21, 10.5, 200))
rect(xleft =-126.2,ybottom = 23.8,xright = -65.5,ytop = 50.6,col = "white")
map("usa", xlim=c(-126.2,-65.5), ylim=c(23.8,50.6),add=T)
map("state", add = T)
map("state", region = "arizona", fill = T, add = T)
# dev.off()

##################
### tree cover ###
##################
tree.cover <- landcover('trees', path = tempdir())
# pdf("../plots/studysite_dist.pdf", width = 6, height = 6)
par(mar = c(0, 0, 0, 0))
aoi <- ext(-115, -106, 29.5, 38)
tree.crop <- crop(tree.cover, aoi)
plot(tree.crop, xlim = c(-115, -106), ylim = c(29.5, 38), col = terrain.colors(100, rev = T))
# study site
map('state', fill = FALSE, xlim = c(-115, -106), ylim = c(29.5, 38.5), xlab = "lon", ylab = "lat", add= T, col = '#636363')
map.scale(x=-114.5, y=30, ratio=FALSE, relwidth=0.2, cex = 0.5)
### species distribution
grantii <- read.csv('~/Desktop/Projects/Dynastes_grantii_pop/data/grantii_gbif.csv')
## all species from GBIF
points(x=grantii$lon, y = grantii$lat, pch = 16, cex = 1, col = rgb(37,37,37, maxColorValue = 255, alpha = 100))
# location 
cols <- brewer.pal(8,name = 'Set1')
points_star <- function(x0, y0, size = 0.1, col = "gold", border = "black", npoints = 5) {
  angles <- seq(0, 2*pi, length.out = 2*npoints + 1)
  radii <- rep(c(size, size/2), npoints)
  x <- x0 + radii * cos(angles)
  y <- y0 + radii * sin(angles)
  polygon(x, y, col = col, border = border)  # <- filled
}
size = 0.15
# Utah
points_star(-112.68175,37.222528, col = cols[3], size = size )
# pine, payson n = 4
points_star(-111.25847,34.25504,col = cols[6], size = size )
points_star(-111.45514,34.38447, col = cols[6], size = size)
# researve, NM n = 10
points_star(-108.756496974, 33.705163846, col = cols[6], size = size)
# Portal
points_star(-109.21473, 31.89373, col = cols[2], size = size)
# Santa Rita mountains madera canyon 
points_star(-110.87374, 31.71387, col = cols[2], size = size)
# mt lemmon n = 3 
points_star(-110.69176, 32.37330, col = cols[2], size = size)
# mt Graham
points_star(-109.817017, 32.644867, col = cols[2], size = size)

# legend
# points(-114.5, 31.7, col = rgb(37,37,37, maxColorValue = 255, alpha = 160), pch =16)
# text(-114.5, 31.7, pos=4, label = 'GBIF database')
x = -114.5
y = 31.6
s = 0.4
size = 0.15
points_star(x, y, col = cols[6], size = size )
text(x, y, pos=4, label = 'Mogollon Rim')
points_star(x, y-s, col = cols[3],  size = size )
text(x, y-s, pos=4, label = 'Utah')
points_star(x, y-2*s, col = cols[2], size = size )
text(x, y-2*s, pos=4, label = 'Sky islands')

points(x, y-3*s, col = rgb(37,37,37, maxColorValue = 255, alpha = 255), pch =16, cex = 1.5)
text(x, y-3*s, pos=4, label = 'GBIF database')


# Inmap
par(usr=c(-280, -21, 10.5, 200))
rect(xleft =-126.2,ybottom = 23.8,xright = -65.5,ytop = 50.6,col = "white")
map("usa", xlim=c(-126.2,-65.5), ylim=c(23.8,50.6),add=T)
map("state", add = T)
map("state", region = "arizona", fill = T, add = T)
# dev.off()

