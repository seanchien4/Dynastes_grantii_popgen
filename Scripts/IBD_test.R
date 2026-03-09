library(geosphere)
library(vegan)
library(viridis)
# IBD test
#######
# ALL #
#######
set.seed(4)
coords <- read.table("../Data/IBD/coords.txt", header = FALSE, col.names = c("ID", "Lat", "Lon"))
geo_dist_matrix <- distm(coords[, c("Lon", "Lat")], fun = distHaversine)
geo_dist_matrix_km <- geo_dist_matrix / 1000
rownames(geo_dist_matrix_km) <- coords$ID
colnames(geo_dist_matrix_km) <- coords$ID

gen_dist_full <- as.matrix(read.table("../Data/IBD/pop1_all_dist.dist", header=FALSE))
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

# plot
col = viridis(1, begin = 0.35, alpha = 0.6)
plot(
  x = geo_dist_obj,
  y = gen_dist_obj,
  main = "",
  xlab = "Geographic Distance (km)",
  ylab = "Genetic Distance (IBS)",
  pch = 20, # Use solid circles for points
  col = col # Use a semi-transparent blue for points
)
abline(
  lm(gen_dist_obj ~ geo_dist_obj), # Calculate the linear model on the fly
  col = "red", # Set the line color to red
  lwd = 2 # Set the line width to be thicker
)
text(250, 125000, label = expression("Mantel statistic "*italic(r)*" = 0.41"), pos = 4)
text(250, 124000, label = expression(italic(p-value)* '= 0.0012'),pos = 4)


###############
# Only Payson # 
###############
set.seed(4)
coords <- read.table("../Data/IBD/coords_payson.txt", header = FALSE, col.names = c("ID", "Lat", "Lon"))
geo_dist_matrix <- distm(coords[, c("Lon", "Lat")], fun = distHaversine)
geo_dist_matrix_km <- geo_dist_matrix / 1000
rownames(geo_dist_matrix_km) <- coords$ID
colnames(geo_dist_matrix_km) <- coords$ID

gen_dist_full <- as.matrix(read.table("../Data/IBD/pop1_payson_dist.dist", header=FALSE))
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
# plot
plot(
  x = geo_dist_obj,
  y = gen_dist_obj,
  main = "",
  xlab = "Geographic Distance (km)",
  ylab = "Genetic Distance (IBS)",
  pch = 20,
  col = col 
)
abline(
  lm(gen_dist_obj ~ geo_dist_obj), 
  col = "red",
  lwd = 2
)

################
# Only Reserve #
################
set.seed(4)
coords <- read.table("../Data/IBD/coords_reserve.txt", header = FALSE, col.names = c("ID", "Lat", "Lon"))
geo_dist_matrix <- distm(coords[, c("Lon", "Lat")], fun = distHaversine)
geo_dist_matrix_km <- geo_dist_matrix / 1000
rownames(geo_dist_matrix_km) <- coords$ID
colnames(geo_dist_matrix_km) <- coords$ID
gen_dist_full <- as.matrix(read.table("../Data/IBD/pop1_reserve_dist.dist", header=FALSE))
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
# plot
plot(
  x = geo_dist_obj,
  y = gen_dist_obj,
  main = "",
  xlab = "Geographic Distance (km)",
  ylab = "Genetic Distance (IBS)",
  pch = 20,
  col = col
)
abline(
  lm(gen_dist_obj ~ geo_dist_obj), 
  col = "red",
  lwd = 2
)



########
# UTAH #
########
set.seed(4)
coords <- read.table("../Data/IBD/utah_coords.txt", header = FALSE, col.names = c("ID", "Lat", "Lon"))
geo_dist_matrix <- distm(coords[, c("Lon", "Lat")], fun = distHaversine)
geo_dist_matrix_km <- geo_dist_matrix / 1000
rownames(geo_dist_matrix_km) <- coords$ID
colnames(geo_dist_matrix_km) <- coords$ID

gen_dist_full <- as.matrix(read.table("../Data/IBD/utah_dist.dist", header=FALSE))
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

# partail IBD test 


library(terra)

# 1. Load your GPS coordinates (IDs, Lat, Lon)
# Ensure this matches the order of your genetic distance matrix
coords <- read.table("../Data/IBD/utah_coords.txt", header = FALSE, col.names = c("ID", "Lat", "Lon"))

# 2. Load the WorldClim 2 .tif files (BIO12, BIO14, etc., or your PD, PW, PC layers)
# Replace the paths with your actual .tif file locations
pd_raster <- rast("../Data/wc2.1_30s_bio/wc2.1_30s_bio_17.tif") # PD (Precip of Driest Month/Quarter)
pw_raster <- rast("../Data/wc2.1_30s_bio/wc2.1_30s_bio_16.tif") # PW (Precip of Wettest Month/Quarter)
pc_raster <- rast("../Data/wc2.1_30s_bio/wc2.1_30s_bio_19.tif") # PC (Precip of Coldest Quarter)

# 3. Combine into a multi-layer raster object
climate_stack <- c(pd_raster, pw_raster, pc_raster)
names(climate_stack) <- c("PD", "PW", "PC")

# 4. Extract values based on your coordinates
# We create a simple data frame of Lon/Lat for extraction
points <- coords[, c("Lon", "Lat")]
extracted_values <- terra::extract(climate_stack, points)

# 5. Combine with IDs and save
utah_climate <- data.frame(ID = coords$ID, extracted_values[, -1]) # -1 removes the ID column terra adds
write.table(utah_climate, "../Data/IBD/utah_climate.txt", sep="\t", row.names=FALSE, quote=FALSE)

library(geosphere)

# 1. Load your climate data (PD, PW, PC values for each individual)
# Assuming a file with columns: ID, PD, PW, PC
env_data <- read.table("../Data/IBD/utah_climate.txt", header = TRUE)
rownames(env_data) <- env_data$ID

# 2. Create the Climate Distance Matrix
# We use Euclidean distance for environmental variables
env_dist_matrix <- dist(env_data[common_ids, c("PD", "PW", "PC")], method = "euclidean")

# 3. Perform the Partial Mantel Tests
# Question A: Is it Geography (IBD) after controlling for Climate?
partial_mantel_geo <- mantel.partial(gen_dist_obj, geo_dist_obj, env_dist_matrix, 
                                     method = "pearson", permutations = 9999)

# Question B: Is it Climate (IBE) after controlling for Geography?
partial_mantel_env <- mantel.partial(gen_dist_obj, env_dist_matrix, geo_dist_obj, 
                                     method = "pearson", permutations = 9999)

print("--- Partial Mantel: Geography (controlling for Climate) ---")
print(partial_mantel_geo)

print("--- Partial Mantel: Climate (controlling for Geography) ---")
print(partial_mantel_env)
