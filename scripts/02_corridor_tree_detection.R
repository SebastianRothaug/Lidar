library(lidR)  # LiDAR reading / processing
library(terra) # Modern raster work-horse (SpatRaster)
library(tmap)  # Nice quick maps (optional)
library(RCSF)
library(sf)

# ---------------------------------------
### Reguired Data:
# AOI as LAS file
# Corridor in type Line as clipping reference to minimize storage
# ---------------------------------------


####  Settings Input

# ---- INPUT LAS/LAZ of the AOI -------------------------------------------------
las_file <- "C:/Users/User/Desktop/Lidar_Course/lidar_holzkirchen.las"   # <-- change if needed
# READ LAS --------------------------------------------------
las <- readLAS(las_file)  
st_crs(las) <- 25832
# Zeigt an, wie viele Punkte pro Klasse im Original existieren
table(las$Classification)

# ---- INPUT Line = Base Line for the Corridor -------------------------------------------------
corridor <- st_read("C:/Users/User/Desktop/Lidar_Course/powerline_project/data/trassenlinie.gpkg") # <-- change if needed
#### ------------------



####  Settings Output

# ---- RESOLUTION OF THE OUTPUT RASTERS -----------------------------
res_chm    <- 0.5   # cell size in metres (smaller -> more detail, larger RAM)
res_dtm    <- 1.0   # DTM can be coarser - you decide

# ---- output folder ---------------------------------------
out_dir    <- "C:/Users/User/Desktop/Lidar_Course"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
#### ------------------



#### Processing 

###1 Corridor: -------------------------------
# Clip the extend of the to improve processing speed and minimize data size
#(you could skip this part and just continue las = full AOI but speed will decrease)

# buffer to define AOI
buffered_line <- st_buffer(corridor, dist = 100) # 100m Corridor from Base Line

# LAS clip
# Transform buffer to the CRS of the LAS
buffered_line_proj <- st_transform(buffered_line, st_crs(las))
# if line has different feature segments it is gonna be united to one
buffered_line_proj <- st_union(buffered_line_proj)
# Clipping of the LAS-Pointcloud with reference to the buffer
las_final <- clip_roi(las, buffered_line_proj)
# (Optional) check if it is a LAS"-object:
class(las_final)

# plot the result
plot(las_final)

# Save the clipped area of the AOI as LAS
aoi_clip_file <- file.path(out_dir, "lidar_holzkirchen_clipped.las")
writeLAS(las_final, aoi_clip_file)

# redefine the las to the now smaller AOI in our environment
las <- las_final
## (you could have skip this part and continued with previous las = full AOI but speed will decrease)
#### ------------------



###2 Pointcloud Processing different Ground, Height Models

###2.1 GROUND CLASSIFICATION = DTM -------------------------------

# Classifying ground points using Cloth-Simulation-Filter
las <- classify_ground(las, csf())

# Building Digital Terrain Model (DTM)
dtm <- rasterize_terrain(las,
                         algorithm = knnidw(k = 10L, p = 2),   # simple IDW
                         res = res_dtm)                        # cell size


###2.2 Normalisation = CANOPY HEIGHT MODEL CHM -------------------------------
# NORMALISE HEIGHTS (Z -> height above ground) -------------
# Difference from just created DTM to higher Points of the LAS pointcloud
las_norm <- normalize_height(las, dtm)   # Z now = "height above ground"

# CREATE CANOPY HEIGHT MODEL (CHM)
#Rasterising normalised points -> CHM 
chm <- rasterize_canopy(las_norm,
                        res = res_chm,
                        algorithm = p2r(0.15))   # point-to-raster (bilinear),`p2r()` is fast and gives a smooth surface.

# SAVE RASTERS TO DISK

dtm_file <- file.path(out_dir, "DTM.tif")
chm_file <- file.path(out_dir, "CHM_TreeHeight.tif")

writeRaster(dtm, filename = dtm_file,
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))

writeRaster(chm, filename = chm_file,
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))

# Ploting Result ----------------------------------------
tmap_mode("plot")   # static map (use "view" for interactive leaflet)

tm_shape(dtm) +
  tm_raster(palette = "-viridis", title = "DTM (m)", alpha = 0.6) +
  tm_shape(chm) +
  tm_raster(palette = "YlOrRd", title = "Tree height (CHM, m)",
            style = "cont") +
  tm_layout(main.title = "LiDAR-derived Tree Height Raster",
            legend.outside = TRUE)

# Generate kernel and smooth chm
kernel <- matrix(1, 3, 3)
schm <- terra::focal(x = chm, w = kernel, fun = median, na.rm = TRUE)
plot(schm)
#### ------------------




####3 Tree Detection

###3.1 Detection of Trees, based on the maxima of the canopy height model
# Moving window to Identify trees: dynamic size aproach depending on the canopy height
f <- function(x) {abs(10*sin(0.03*x)+2)}
heights <- seq(0,60,5)
ws <- f(heights)
plot(heights, ws, type = "l", ylim = c(0,20), xlim = c(0,60))

# Detection of maxima as tree tops
ttops <- locate_trees(las = schm, algorithm = lmf(ws=f,hmin = 3.5)) # windowsize = f dynamic function and minal trees 3.5m
ttops

# plot the detection result
plot(chm)
plot(ttops, col = "red", add = TRUE, cex = 0.5)

# Segment trees using dalponte
las <- segment_trees(las = las, algorithm = dalponte2016(chm = schm, treetops = ttops))

# Count number of trees detected and segmented
length(unique(las$treeID) |> na.omit())

# Result is a a point layer


###3.2 Extraction of height information for the detected point

# check if SpatRaster for extraction
if (!inherits(dtm, "SpatRaster")) {
  dtm <- terra::rast(dtm)
}
# check if SpatVector for extraction
if (inherits(ttops, "SpatVector")) {
  ttops_vect <- ttops
} else if (inherits(ttops, "sf")) {
  ttops_vect <- terra::vect(ttops)
} else {
  ttops_vect <- terra::vect(sf::st_as_sf(ttops))
}


# DTM values extracted at the tree top
dgm_extract <- terra::extract(dtm, ttops_vect)
dgm_height  <- dgm_extract[[2]]  # 2. Spalte = Rasterwert des DTM

# crate tree top output point layer
ttops_out <- ttops
if (!inherits(ttops_out, "sf")) {
  ttops_out <- sf::st_as_sf(ttops_out)
}

# Column Tree Height: height above ground from Normilazation derrived
ttops_out$tree_height   <- ttops$Z
# column Terrain Height: height of the Terrain at the tree location
ttops_out$dgm_height    <- dgm_height
# column Canopy height: Terrain height + Tree Height
ttops_out$canopy_height <- ttops_out$dgm_height + ttops_out$tree_height 

# check CRS
if (is.na(sf::st_crs(ttops_out))) {
  sf::st_crs(ttops_out) <- 25832
}

# Trees outside of the Raster: NA values
n_na <- sum(is.na(ttops_out$dgm_height))
if (n_na > 0) {
  warning(sprintf("%d von %d Baeumen liegen ausserhalb des DTM-Rasters (dgm_height = NA).",
                  n_na, nrow(ttops_out)))
}

# only keep this relevant columns in the new tree top output file
keep_cols <- intersect(c("treeID", "dgm_height", "tree_height", "canopy_height"),
                       names(ttops_out))
ttops_out <- ttops_out[, keep_cols]

## Save as GeoPackage Point Layer resampling tree with all 3 height information ------
----
gpkg_path <- file.path(out_dir, "tree_tops_EPSG25832.gpkg")

sf::st_write(
  obj          = ttops_out,
  dsn          = gpkg_path,
  layer        = "tree_tops",
  delete_layer = TRUE,
  quiet        = FALSE
)
#-------
cat("\n[OK] GeoPackage saved:\n", gpkg_path, "\n")
