library(geosphere)
library(vegan)
library(viridis)
library(terra)
library(geosphere)
# IBD test
#######
# ALL #
#######
set.seed(4)
coords <- read.table("~/Desktop/Projects/Dynastes_grantii_pop/Data/IBD/coords.txt", header = FALSE, col.names = c("ID", "Lat", "Lon"))
geo_dist_matrix <- distm(coords[, c("Lon", "Lat")], fun = distHaversine)
geo_dist_matrix_km <- geo_dist_matrix / 1000
rownames(geo_dist_matrix_km) <- coords$ID
colnames(geo_dist_matrix_km) <- coords$ID

gen_dist_full <- as.matrix(read.table("~/Desktop/Projects/Dynastes_grantii_pop/Data/IBD/all_dist.dist", header=FALSE))
rownames(gen_dist_full) <- coords$ID
colnames(gen_dist_full) <- coords$ID
# Ensure the order of individuals is the same in both matrices
common_ids <- intersect(rownames(gen_dist_full), rownames(geo_dist_matrix_km))
gen_dist_ordered <- gen_dist_full[common_ids, common_ids]
geo_dist_ordered <- geo_dist_matrix_km[common_ids, common_ids]
# Convert to 'dist' objects for the mantel test
gen_dist_obj <- as.dist(gen_dist_ordered)
geo_dist_obj <- as.dist(geo_dist_ordered)

mantel_result <- mantel(gen_dist_obj, geo_dist_obj, method = "pearson", permutations = 9999)
print(mantel_result)


# Calculate R² explicitly
r_value <- mantel_result$statistic
r_squared <- r_value^2
print(paste("R² =", round(r_squared, 3)))
# Report: "r=0.786 (R²=0.62)

# plot
col = viridis(1, begin = 0.35, alpha = 0.6)
par(mar = c(4,4,1,1))
plot(
  x = geo_dist_obj,
  y = gen_dist_obj,
  main = "",
  xlab = "Geographic Distance (km)",
  ylab = "Genetic Distance",
  pch = 20, # Use solid circles for points
  col = col # Use a semi-transparent blue for points
)
abline(
  lm(gen_dist_obj ~ geo_dist_obj), 
  col = "red", 
  lwd = 1
)
text(250, 120000,labels = expression("Mantel statistic "*italic(r)*" = 0.79, " * R^2 == 0.62),pos = 4)
text(250, 118000, label = expression(italic(p-value)* '< 0.0001'),pos = 4)
