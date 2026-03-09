library(phytools)
library(mapdata)
tree <- read.tree('Newick Export.nwk')
plot(tree)
ultrametric_tree <- chronos(tree)
tips <- ultrametric_tree$tip.label
tips <- gsub('DgSR', 'SR', tips)
tips <- gsub('DgrP', 'P', tips)
tips <- gsub('DgrR', 'R', tips)
tips <- gsub('DgUT', 'UT', tips)
tips <- gsub('DgrCO0', 'C', tips)
tips <- gsub('DGS0', 'M', tips)
ultrametric_tree$tip.label <- tips


cols <- viridis(4, end = 1)
pop.cols <- tips 
pop.cols[20:24] <- cols[1]
pop.cols[c(1:4,10)] <- cols[2]
pop.cols[c(5:6, 8:9,11:19)] <- cols[4]
pop.cols[7] <- cols[3]
plot.phylo(ultrametric_tree, tip.color = pop.cols )

# create coordination matrix
loc.matrix <- matrix(NA, nrow = length(ultrametric_tree$tip.label), ncol = 2)
rownames(loc.matrix) <- ultrametric_tree$tip.label
colnames(loc.matrix) <- c('lat', 'long')
for (i in 1:nrow(loc.matrix)){
  if (grepl('R',rownames(loc.matrix)[i])){
    loc.matrix[i,] <- c(33.70516385,-108.756497)
  }
  if (grepl('UT',rownames(loc.matrix)[i])){
    loc.matrix[i,] <- c(37.222528,-112.68175)
  }
  if (grepl('P',rownames(loc.matrix)[i])){
    loc.matrix[i,] <- c(34.38447,-111.45514)
  }
  if (grepl('SR',rownames(loc.matrix)[i])){
    loc.matrix[i,] <- c(31.8259,-110.7748)
  }
  if (grepl('C',rownames(loc.matrix)[i])){
    loc.matrix[i,] <- c(31.8494,-109.6392)
  }
  if (grepl('M',rownames(loc.matrix)[i])){
    loc.matrix[i,] <- c(29.59,-109.87)
  }
}
obj <- phylo.to.map(ultrametric_tree, loc.matrix, plot = F,
                    database="state",
                    regions = c("new mexico", "arizona", "utah"))
# Plot the map, zooming into the specific states
plot(obj, direction = "rightwards" ,ftype = 'i',
     fsize = 1,lty="dashed",colors="chartreuse", delimit_map = T)


