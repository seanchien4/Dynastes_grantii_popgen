library(phytools)
library(mapdata)
map('state', fill = F, xlim = c(-115, -106), ylim = c(29.5, 38), xlab = "lon", ylab = "lat")
tree <- rcoal(n = 9)
# coorinates matrix
coor <- read.csv('../Data/sample_location.csv', header = T)
rownames(coor) <- coor$X
coor <- coor[,-1]
tree$tip.label <- rownames(coor)
obj <- phylo.to.map(tree, coor, plot = F,direction="rightwards",
                    database="state",
                    regions = c("new mexico", "arizona", "utah"))
# Plot the map, zooming into the specific states
plot(obj, direction = "rightwards" ,ftype = 'i', xlim = c(-115, -106), ylim = c(29.5, 38),
     fsize = 0.0003,lty="dashed",colors="chartreuse")

data(tortoise.tree)
data(tortoise.geog)
tortoise.phymap<-phylo.to.map(tortoise.tree,
                              tortoise.geog,plot=FALSE,direction="rightwards",
                              regions="Ecuador")
plot(tortoise.phymap,direction="rightwards",pts=FALSE,
     xlim=c(-92.25,-89.25),ylim=c(-1.8,0.75),ftype="i",
     fsize=0.8,lty="dashed",map.bg="lightgreen",
     colors="slategrey")
## reset margins
par(mar=c(5.1,4.1,4.1,2.1))