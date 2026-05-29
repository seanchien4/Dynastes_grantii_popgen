library(tidyverse)
library(viridis)
library(RColorBrewer)
library(car)
library(dplyr)
library(plotly)
library(scatterplot3d)
options(rgl.useNULL = TRUE)
library(rgl)
library(magick)
pca <- read.table('~/Desktop/Projects/Dynastes_grantii_pop/Data/PCA/final.eigenvec')
eigenval <- scan("~/Desktop/Projects/Dynastes_grantii_pop/Data/PCA/final.eigenval")
pca <- pca[,-1]
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
# location data
loc <- rep(NA, length(pca$ind))
loc[grep("DgUT", pca$ind)] <- 'Utah'
loc[grep("DgSR", pca$ind)] <- 'Unknown'
loc[grep("DGS01", pca$ind)] <- 'Mexico'
loc[grep("Dgr", pca$ind)] <- 'Mogollon_Rim'
loc[grep("DgrCO", pca$ind)] <- 'Cochise_County'

loc[grep("DgPT", pca$ind)] <- 'Portal'
loc[grep("DgML", pca$ind)] <- 'Mt_Lemmon'
loc[grep("DgMC", pca$ind)] <- 'Madera_canyon'
loc[grep("DgMG", pca$ind)] <- 'Mt_Graham'

pca <- as.tibble(data.frame(pca,loc))
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
par(mfrow = c(2, 2))
# barplot to show the % of variance each principal component explains
bp <- barplot(pve$pve, ylim = c(0, 10), xaxt = "n")
axis(1, at = bp, labels = paste0("PC", 1:20), las = 2, cex.axis = 0.6)

# plot
# cols = viridis(8, alpha = 0.8)
cols <- brewer.pal(9,name = 'Set1')
loc_col <- rep(NA, length(pca$ind))
loc_col[grep("Mogollon_Rim", pca$loc)] <- cols[6]
loc_col[grep("Madera_canyon", pca$loc)] <- cols[1]
loc_col[grep("Portal", pca$loc)] <- cols[2]
loc_col[grep("Mt_Lemmon", pca$loc)] <- cols[7]
loc_col[grep("Mt_Graham", pca$loc)] <- cols[5]
loc_col[grep("Unknown", pca$loc)] <- cols[4]
loc_col[grep("Mexico", pca$loc)] <- cols[8]
loc_col[grep("Utah", pca$loc)] <- cols[3]
loc_col[grep("Cochise_County", pca$loc)] <- cols[9]

lwd = 0.4
cex = 0.8
plot(x = pca$PC1, y = pca$PC2, pch = 21, bg = loc_col,lwd = lwd, cex = cex, 
     xlab = 'PC1(9.9%)', ylab = 'PC2(6.5%)',)
     #xlim = rev(range(pca$PC2)))
x = 0.28
y = 0.25
s = -0.02
pch = 21
lwd = 0.2
cex = 0.5

points(x,y, pch = pch, bg = cols[6],lwd = lwd, cex = cex)
text(x, y, pos =4 , labels = c("Mogollon Rim"), cex = cex)
points(x,y+s*1, pch = pch, bg = cols[1],lwd = lwd, cex = cex)
text(x, y+s*1, pos =4 , labels = c("Madera Canyon"), cex = cex)
points(x,y+s*2, pch = pch, bg = cols[2],lwd = lwd, cex = cex)
text(x, y+s*2, pos =4 , labels = c("Portal"), cex = cex)
points(x,y+s*3, pch = pch, bg = cols[7],lwd = lwd, cex = cex)
text(x, y+s*3, pos =4 , labels = c("Mt. Lemmon"), cex = cex)
points(x,y+s*4, pch = pch, bg = cols[5],lwd = lwd, cex = cex)
text(x, y+s*4, pos =4 , labels = c("Mt. Graham"), cex = cex)
points(x,y+s*5, pch = pch, bg = cols[9],lwd = lwd, cex = cex)
text(x, y+s*5, pos =4 , labels = c("Cochise County"), cex = cex)
points(x,y+s*6, pch = pch, bg = cols[4],lwd = lwd, cex = cex)
text(x, y+s*6, pos =4 , labels = c("Unknown"), cex = cex)
points(x,y+s*7, pch = pch, bg = cols[8],lwd = lwd, cex = cex)
text(x, y+s*7, pos =4 , labels = c("Mexico"), cex = cex)
points(x,y+s*8, pch = pch, bg = cols[3],lwd = lwd, cex = cex)
text(x, y+s*8, pos =4 , labels = c("Utah"), cex = cex)

############
# Other PC #
############
lwd = 0.4
cex = 0.8
plot(x = pca$PC1, y = pca$PC3, pch = 21, bg = loc_col,lwd = lwd, cex = cex, 
     xlab = 'PC1(9.9%)', ylab = 'PC3(5.7%)',)
x = 0.28
y = 0.34
s = -0.028
pch = 21
lwd = 0.2
cex = 0.5
points(x,y, pch = pch, bg = cols[6],lwd = lwd, cex = cex)
text(x, y, pos =4 , labels = c("Mogollon Rim"), cex = cex)
points(x,y+s*1, pch = pch, bg = cols[1],lwd = lwd, cex = cex)
text(x, y+s*1, pos =4 , labels = c("Madera Canyon"), cex = cex)
points(x,y+s*2, pch = pch, bg = cols[2],lwd = lwd, cex = cex)
text(x, y+s*2, pos =4 , labels = c("Portal"), cex = cex)
points(x,y+s*3, pch = pch, bg = cols[7],lwd = lwd, cex = cex)
text(x, y+s*3, pos =4 , labels = c("Mt. Lemmon"), cex = cex)
points(x,y+s*4, pch = pch, bg = cols[5],lwd = lwd, cex = cex)
text(x, y+s*4, pos =4 , labels = c("Mt. Graham"), cex = cex)
points(x,y+s*5, pch = pch, bg = cols[9],lwd = lwd, cex = cex)
text(x, y+s*5, pos =4 , labels = c("Cochise County"), cex = cex)
points(x,y+s*6, pch = pch, bg = cols[4],lwd = lwd, cex = cex)
text(x, y+s*6, pos =4 , labels = c("Unknown"), cex = cex)
points(x,y+s*7, pch = pch, bg = cols[8],lwd = lwd, cex = cex)
text(x, y+s*7, pos =4 , labels = c("Mexico"), cex = cex)
points(x,y+s*8, pch = pch, bg = cols[3],lwd = lwd, cex = cex)
text(x, y+s*8, pos =4 , labels = c("Utah"), cex = cex)

lwd = 0.4
cex = 0.8
plot(x = pca$PC2, y = pca$PC3, pch = 21, bg = loc_col,lwd = lwd, cex = cex, 
     xlab = 'PC2(6.5%)', ylab = 'PC3(5.7%)',)
x = 0.12
y = 0.34
s = -0.028
pch = 21
lwd = 0.2
cex = 0.5
points(x,y, pch = pch, bg = cols[6],lwd = lwd, cex = cex)
text(x, y, pos =4 , labels = c("Mogollon Rim"), cex = cex)
points(x,y+s*1, pch = pch, bg = cols[1],lwd = lwd, cex = cex)
text(x, y+s*1, pos =4 , labels = c("Madera Canyon"), cex = cex)
points(x,y+s*2, pch = pch, bg = cols[2],lwd = lwd, cex = cex)
text(x, y+s*2, pos =4 , labels = c("Portal"), cex = cex)
points(x,y+s*3, pch = pch, bg = cols[7],lwd = lwd, cex = cex)
text(x, y+s*3, pos =4 , labels = c("Mt. Lemmon"), cex = cex)
points(x,y+s*4, pch = pch, bg = cols[5],lwd = lwd, cex = cex)
text(x, y+s*4, pos =4 , labels = c("Mt. Graham"), cex = cex)
points(x,y+s*5, pch = pch, bg = cols[9],lwd = lwd, cex = cex)
text(x, y+s*5, pos =4 , labels = c("Cochise County"), cex = cex)
points(x,y+s*6, pch = pch, bg = cols[4],lwd = lwd, cex = cex)
text(x, y+s*6, pos =4 , labels = c("Unknown"), cex = cex)
points(x,y+s*7, pch = pch, bg = cols[8],lwd = lwd, cex = cex)
text(x, y+s*7, pos =4 , labels = c("Mexico"), cex = cex)
points(x,y+s*8, pch = pch, bg = cols[3],lwd = lwd, cex = cex)
text(x, y+s*8, pos =4 , labels = c("Utah"), cex = cex)

# for (group in unique(pca$loc)) {
#   data_subset <- subset(pca, loc == group)
#   if (nrow(data_subset) > 2) {
#     dataEllipse(data_subset$PC1, data_subset$PC2, levels = 0.95, plot.points = FALSE,
#                 center.pch = FALSE, col = loc_col[which(pca$loc == group)][1], lwd = 2)
#   }
# }


# plot 3D 
# Create the 3D plot
# Use 'bg' for the fill color and 'color' for the border color
s3d <- scatterplot3d(
  x = pca$PC1,
  y = pca$PC2,
  z = pca$PC3,
  xlab = paste0("PC1 (", round(pve$pve[1], 1), "%)"),
  ylab = paste0("PC2 (", round(pve$pve[2], 1), "%)"),
  zlab = paste0("PC3 (", round(pve$pve[3], 1), "%)"),
  pch = 21,        # Use a filled circle shape
  bg = loc_col,    # Use your loc_col vector for the fill (background) color
  color = "black", # Use black for the point borders
  lwd = 0.5,
  cex.symbols = 1.5,
  main = "",
  angle = 40,
  color.grid = "lightgray",
  col.axis = "gray50" # Makes the main box edges gray instead of black
)
# locator(1)
# Add a legend
x_legend = 1.5
y_legend = 6.8
s_legend = -0.32  
pch_legend = 22
lwd_legend = 1
cex_legend = 1.2

# Add legend points and text manually
points(x_legend, y_legend, pch = pch_legend, bg = cols[6], lwd = lwd_legend, cex = cex_legend)
text(x_legend, y_legend, pos = 4, labels = "Mogollon Rim", cex = 0.9)

points(x_legend, y_legend + s_legend * 1, pch = pch_legend, bg = cols[1], lwd = lwd_legend, cex = cex_legend)
text(x_legend, y_legend + s_legend * 1, pos = 4, labels = "Madera Canyon", cex = 0.9)

points(x_legend, y_legend + s_legend * 2, pch = pch_legend, bg = cols[2], lwd = lwd_legend, cex = cex_legend)
text(x_legend, y_legend + s_legend * 2, pos = 4, labels = "Portal", cex = 0.9)

points(x_legend, y_legend + s_legend * 3, pch = pch_legend, bg = cols[7], lwd = lwd_legend, cex = cex_legend)
text(x_legend, y_legend + s_legend * 3, pos = 4, labels = "Mt. Lemmon", cex = 0.9)

points(x_legend, y_legend + s_legend * 4, pch = pch_legend, bg = cols[5], lwd = lwd_legend, cex = cex_legend)
text(x_legend, y_legend + s_legend * 4, pos = 4, labels = "Mt. Graham", cex = 0.9)

points(x_legend, y_legend + s_legend * 5, pch = pch_legend, bg = cols[9], lwd = lwd_legend, cex = cex_legend)
text(x_legend, y_legend + s_legend * 5, pos = 4, labels = "Cochise County", cex = 0.9)

points(x_legend, y_legend + s_legend * 6, pch = pch_legend, bg = cols[4], lwd = lwd_legend, cex = cex_legend)
text(x_legend, y_legend + s_legend * 6, pos = 4, labels = "Unknown", cex = 0.9)

points(x_legend, y_legend + s_legend * 7, pch = pch_legend, bg = cols[8], lwd = lwd_legend, cex = cex_legend)
text(x_legend, y_legend + s_legend * 7, pos = 4, labels = "Mexico", cex = 0.9)

points(x_legend, y_legend + s_legend * 8, pch = pch_legend, bg = cols[3], lwd = lwd_legend, cex = cex_legend)
text(x_legend, y_legend + s_legend * 8, pos = 4, labels = "Utah", cex = 0.9)

## animation 
options(rgl.printRglwidget = TRUE)
# Static chart
par3d(windowRect = c(0, 0, 1000, 1000))
par3d(cex.lab = 1, cex.axis = 0.7)
plot3d( pca$PC1, pca$PC2, pca$PC3, col = loc_col, type = "s", 
        radius = .015, size = 1,
        xlab = "PC1 (9.9%)",  
        ylab = "PC2 (6.5%)",  
        zlab = "PC3 (5.7%)" )
# play3d( spin3d( axis = c(0, 0, 1), rpm = 10,dev = cur3d()),startTime = 0, duration = 10 )

frame_files <- movie3d(
  movie = "~/Desktop/Projects/Dynastes_grantii_pop/plots/PCA_animation/3dAnimatedScatterplot",
  spin3d(axis = c(0, 0, 1), rpm = 2),
  duration = 30,
  dir = ".",
  type = "png", 
  clean = FALSE,
  convert = FALSE
)

png_files <- list.files(path = "../plots/PCA_animation/", 
                        pattern = "3dAnimatedScatterplot.*\\.png$", 
                        full.names = TRUE)
if (length(png_files) > 0) {
  animation <- image_read(png_files)
  animation_joined <- image_join(animation)
  
  # Assuming 100 frames over 10 seconds = 10 fps
  animation_final <- image_animate(animation_joined, fps = 10)
  image_write(animation_final, "../plots/PCA_animation/3dAnimatedScatterplot_final.gif")
} 
file.remove(png_files)



##########
# Reduce #
##########
pca <- read.table('../Data/PCA/reduce.eigenvec')
eigenval <- scan("../Data/PCA/reduce.eigenval")
pca <- pca[,-1]
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
# location data
loc <- rep(NA, length(pca$ind))
loc[grep("DgUT", pca$ind)] <- 'Utah'
loc[grep("Dgr", pca$ind)] <- 'Mogollon_Rim'
loc[grep("DgPT", pca$ind)] <- 'Portal'
loc[grep("DgML", pca$ind)] <- 'Mt_Lemmon'

pca <- as.tibble(data.frame(pca,loc))
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
par(mfrow = c(2, 2))
# barplot to show the % of variance each principal component explains
bp <- barplot(pve$pve, ylim = c(0, 10), xaxt = "n")
axis(1, at = bp, labels = paste0("PC", 1:20), las = 2, cex.axis = 0.6)

# plot
# cols = viridis(8, alpha = 0.8)
cols <- brewer.pal(9,name = 'Set1')
loc_col <- rep(NA, length(pca$ind))
loc_col[grep("Mogollon_Rim", pca$loc)] <- cols[6]
loc_col[grep("Portal", pca$loc)] <- cols[2]
loc_col[grep("Mt_Lemmon", pca$loc)] <- cols[7]
loc_col[grep("Utah", pca$loc)] <- cols[3]

lwd = 0.4
cex = 0.8
plot(x = pca$PC1, y = pca$PC2, pch = 21, bg = loc_col,lwd = lwd, cex = cex, 
     xlab = 'PC1(9.9%)', ylab = 'PC2(6.5%)',)
#xlim = rev(range(pca$PC2)))
x = 0.28
y = 0.25
s = -0.02
pch = 21
lwd = 0.2
cex = 0.5

points(x,y, pch = pch, bg = cols[6],lwd = lwd, cex = cex)
text(x, y, pos =4 , labels = c("Mogollon Rim"), cex = cex)
points(x,y+s*1, pch = pch, bg = cols[1],lwd = lwd, cex = cex)
text(x, y+s*1, pos =4 , labels = c("Madera Canyon"), cex = cex)
points(x,y+s*2, pch = pch, bg = cols[2],lwd = lwd, cex = cex)
text(x, y+s*2, pos =4 , labels = c("Portal"), cex = cex)
points(x,y+s*3, pch = pch, bg = cols[7],lwd = lwd, cex = cex)
text(x, y+s*3, pos =4 , labels = c("Mt. Lemmon"), cex = cex)
points(x,y+s*4, pch = pch, bg = cols[5],lwd = lwd, cex = cex)
text(x, y+s*4, pos =4 , labels = c("Mt. Graham"), cex = cex)
points(x,y+s*5, pch = pch, bg = cols[9],lwd = lwd, cex = cex)
text(x, y+s*5, pos =4 , labels = c("Cochise County"), cex = cex)
points(x,y+s*6, pch = pch, bg = cols[4],lwd = lwd, cex = cex)
text(x, y+s*6, pos =4 , labels = c("Unknown"), cex = cex)
points(x,y+s*7, pch = pch, bg = cols[8],lwd = lwd, cex = cex)
text(x, y+s*7, pos =4 , labels = c("Mexico"), cex = cex)
points(x,y+s*8, pch = pch, bg = cols[3],lwd = lwd, cex = cex)
text(x, y+s*8, pos =4 , labels = c("Utah"), cex = cex)

############
# Other PC #
############
lwd = 0.4
cex = 0.8
plot(x = pca$PC1, y = pca$PC3, pch = 21, bg = loc_col,lwd = lwd, cex = cex, 
     xlab = 'PC1(9.9%)', ylab = 'PC3(5.7%)',)
x = 0.28
y = 0.34
s = -0.028
pch = 21
lwd = 0.2
cex = 0.5
points(x,y, pch = pch, bg = cols[6],lwd = lwd, cex = cex)
text(x, y, pos =4 , labels = c("Mogollon Rim"), cex = cex)
points(x,y+s*1, pch = pch, bg = cols[1],lwd = lwd, cex = cex)
text(x, y+s*1, pos =4 , labels = c("Madera Canyon"), cex = cex)
points(x,y+s*2, pch = pch, bg = cols[2],lwd = lwd, cex = cex)
text(x, y+s*2, pos =4 , labels = c("Portal"), cex = cex)
points(x,y+s*3, pch = pch, bg = cols[7],lwd = lwd, cex = cex)
text(x, y+s*3, pos =4 , labels = c("Mt. Lemmon"), cex = cex)
points(x,y+s*4, pch = pch, bg = cols[5],lwd = lwd, cex = cex)
text(x, y+s*4, pos =4 , labels = c("Mt. Graham"), cex = cex)
points(x,y+s*5, pch = pch, bg = cols[9],lwd = lwd, cex = cex)
text(x, y+s*5, pos =4 , labels = c("Cochise County"), cex = cex)
points(x,y+s*6, pch = pch, bg = cols[4],lwd = lwd, cex = cex)
text(x, y+s*6, pos =4 , labels = c("Unknown"), cex = cex)
points(x,y+s*7, pch = pch, bg = cols[8],lwd = lwd, cex = cex)
text(x, y+s*7, pos =4 , labels = c("Mexico"), cex = cex)
points(x,y+s*8, pch = pch, bg = cols[3],lwd = lwd, cex = cex)
text(x, y+s*8, pos =4 , labels = c("Utah"), cex = cex)

lwd = 0.4
cex = 0.8
plot(x = pca$PC2, y = pca$PC3, pch = 21, bg = loc_col,lwd = lwd, cex = cex, 
     xlab = 'PC2(6.5%)', ylab = 'PC3(5.7%)',)
x = 0.12
y = 0.34
s = -0.028
pch = 21
lwd = 0.2
cex = 0.5
points(x,y, pch = pch, bg = cols[6],lwd = lwd, cex = cex)
text(x, y, pos =4 , labels = c("Mogollon Rim"), cex = cex)
points(x,y+s*1, pch = pch, bg = cols[1],lwd = lwd, cex = cex)
text(x, y+s*1, pos =4 , labels = c("Madera Canyon"), cex = cex)
points(x,y+s*2, pch = pch, bg = cols[2],lwd = lwd, cex = cex)
text(x, y+s*2, pos =4 , labels = c("Portal"), cex = cex)
points(x,y+s*3, pch = pch, bg = cols[7],lwd = lwd, cex = cex)
text(x, y+s*3, pos =4 , labels = c("Mt. Lemmon"), cex = cex)
points(x,y+s*4, pch = pch, bg = cols[5],lwd = lwd, cex = cex)
text(x, y+s*4, pos =4 , labels = c("Mt. Graham"), cex = cex)
points(x,y+s*5, pch = pch, bg = cols[9],lwd = lwd, cex = cex)
text(x, y+s*5, pos =4 , labels = c("Cochise County"), cex = cex)
points(x,y+s*6, pch = pch, bg = cols[4],lwd = lwd, cex = cex)
text(x, y+s*6, pos =4 , labels = c("Unknown"), cex = cex)
points(x,y+s*7, pch = pch, bg = cols[8],lwd = lwd, cex = cex)
text(x, y+s*7, pos =4 , labels = c("Mexico"), cex = cex)
points(x,y+s*8, pch = pch, bg = cols[3],lwd = lwd, cex = cex)
text(x, y+s*8, pos =4 , labels = c("Utah"), cex = cex)

