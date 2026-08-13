# If you don't have them yet, uncomment the install lines:
# install.packages(c("lidR", "terra", "tmap", "RCSF", "sf"))

library(lidR)  # LiDAR reading / processing
library(terra) # Modern raster work-horse (SpatRaster)
library(tmap)  # Nice quick maps (optional)
library(RCSF)
library(sf)

### Reguired Data:
# AOI as LAS file
# Corridor in type Line as clipping reference to minimize storage
####


# ---- INPUT LAS/LAZ -------------------------------------------------
las_file <- "C:/Users/User/Desktop/Lidar_Course/lidar_holzkirchen.las"   # <-- change if needed

# ---- RESOLUTION OF THE OUTPUT RASTERS -----------------------------
res_chm    <- 0.5   # cell size in metres (smaller -> more detail, larger RAM)
res_dtm    <- 1.0   # DTM can be coarser - you decide

# ---- OPTIONAL: output folder ---------------------------------------
out_dir    <- "C:/Users/User/Desktop/Lidar_Course"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# READ LAS --------------------------------------------------
las <- readLAS(las_file)   # keep only X,Y,Z (faster)
st_crs(las) <- 25832
if (is.empty(las)) {
  stop("[FEHLER] Could not read the LAS/LAZ file - check the path and format.")
}

# Zeigt an, wie viele Punkte pro Klasse im Original existieren
table(las$Classification)

# Clip to extend to improve processing speed and minimize data size

corridor <- st_read("C:/Users/User/Desktop/Lidar_Course/powerline_project/data/trassenlinie.gpkg")

# buffer to define AOI
buffered_line <- st_buffer(corridor, dist = 100)

# LAS clip
# Den sf-Puffer in das Koordinatensystem der LAS-Datei transformieren
buffered_line_proj <- st_transform(buffered_line, st_crs(las))
# Verschmilzt alle Segmente zu einem zusammenhaengenden Puffer
buffered_line_proj <- st_union(buffered_line_proj)

# Punktwolke auf den Puffer zuschneiden
las_final <- clip_roi(las, buffered_line_proj)

# (Optional) Zur Sicherheit pruefen, ob es jetzt ein "LAS"-Objekt ist:
class(las_final)

# 2. Die Punktwolke in 3D plotten
plot(las_final)

# 3. Die zugeschnittene Datei speichern
writeLAS(las_final, "C:/Users/User/Desktop/Lidar_Course/lidar_holzkirchen_clipped.las")

# rename
las <- las_final


# GROUND CLASSIFICATION & DTM -------------------------------

cat("Classifying ground points ...\n")
# A fast, robust classifier (CSF works well on most terrains)
las <- classify_ground(las, csf())

cat("Building Digital Terrain Model (DTM) ...\n")
dtm <- rasterize_terrain(las,
                         algorithm = knnidw(k = 10L, p = 2),   # simple IDW
                         res = res_dtm)                        # cell size


# NORMALISE HEIGHTS (Z -> height above ground) -------------

cat("Normalising point heights to the ground ...\n")
las_norm <- normalize_height(las, dtm)   # Z now = "height above ground"

# CREATE CANOPY HEIGHT MODEL (CHM) -------------------------

cat("Rasterising normalised points -> CHM ...\n")

## Several rasterisation algorithms exist; `p2r()` is fast and gives a smooth surface.
chm <- rasterize_canopy(las_norm,
                        res = res_chm,
                        algorithm = p2r(0.15))   # point-to-raster (bilinear)


# At this stage *chm* already contains "tree height above ground".
# If you prefer to compute it as DSM - DTM, you could also do:
# dsm <- rasterize_canopy(las, res = res_chm, algorithm = p2r())
# chm_alt <- dsm - dtm   # should be identical (modulo edge effects)


# QUICK VISUAL CHECK ----------------------------------------
cat("Plotting the result ...\n")
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


### Detect trees -----
#f <- function(x) {x * 0.25 + 5}
f <- function(x) {abs(10*sin(0.03*x)+2)}
heights <- seq(0,60,5)
ws <- f(heights)
plot(heights, ws, type = "l", ylim = c(0,20), xlim = c(0,60))


ttops <- locate_trees(las = schm, algorithm = lmf(ws=f,hmin = 3.5))
ttops
plot(chm)
plot(ttops, col = "red", add = TRUE, cex = 0.5)

# Segment trees using dalponte
las <- segment_trees(las = las, algorithm = dalponte2016(chm = schm, treetops = ttops))

# Count number of trees detected and segmented
length(unique(las$treeID) |> na.omit())


# SAVE RASTERS TO DISK -----------------------------------------------

dtm_file <- file.path(out_dir, "DTM.tif")
chm_file <- file.path(out_dir, "CHM_TreeHeight.tif")

cat("Writing DTM ->", dtm_file, "\n")
writeRaster(dtm, filename = dtm_file,
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))

cat("Writing CHM (tree-height) ->", chm_file, "\n")
writeRaster(chm, filename = chm_file,
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))


## BAUMHOEHEN EXTRAHIEREN: dgm_height, tree_height, canopy_height -----------
# ttops ist bereits ein sf-Objekt (aus locate_trees()) mit Spalten
# treeID, Z, geometry. Z = Hoehe ueber Grund, da schm hoehennormalisiert ist.

cat("Extrahiere DGM-Hoehe an den Baumspitzen-Positionen ...\n")

# Absichern: dtm muss ein SpatRaster sein (aeltere lidR-Funktionen wie
# grid_terrain() liefern noch ein RasterLayer aus dem raster-Paket)
if (!inherits(dtm, "SpatRaster")) {
  dtm <- terra::rast(dtm)
}

# Absichern: ttops muss ein SpatVector sein (kann je nach lidR-Version
# als sf, sp/SpatialPointsDataFrame oder bereits SpatVector vorliegen)
if (inherits(ttops, "SpatVector")) {
  ttops_vect <- ttops
} else if (inherits(ttops, "sf")) {
  ttops_vect <- terra::vect(ttops)
} else {
  ttops_vect <- terra::vect(sf::st_as_sf(ttops))
}

# DGM-Wert (Gelaendehoehe) an jeder Baumspitzen-Position extrahieren
dgm_extract <- terra::extract(dtm, ttops_vect)
dgm_height  <- dgm_extract[[2]]  # 2. Spalte = Rasterwert des DTM

ttops_out <- ttops
if (!inherits(ttops_out, "sf")) {
  ttops_out <- sf::st_as_sf(ttops_out)
}
ttops_out$tree_height   <- ttops$Z                                       # Hoehe ueber Grund (aus CHM)
ttops_out$dgm_height    <- dgm_height                                    # Gelaendehoehe (DTM) am Baumstandort
ttops_out$canopy_height <- ttops_out$dgm_height + ttops_out$tree_height  # DGM + Baumhoehe

# CRS absichern, falls ttops (noch) keine CRS-Info traegt
if (is.na(sf::st_crs(ttops_out))) {
  sf::st_crs(ttops_out) <- 25832
}

# Baeume ausserhalb des DTM-Rasters (dgm_height = NA) melden
n_na <- sum(is.na(ttops_out$dgm_height))
if (n_na > 0) {
  warning(sprintf("%d von %d Baeumen liegen ausserhalb des DTM-Rasters (dgm_height = NA).",
                  n_na, nrow(ttops_out)))
}

# Nur relevante Spalten behalten (treeID nur, falls vorhanden)
keep_cols <- intersect(c("treeID", "dgm_height", "tree_height", "canopy_height"),
                       names(ttops_out))
ttops_out <- ttops_out[, keep_cols]

## Als GeoPackage speichern --------------------------------------------------

gpkg_path <- file.path(out_dir, "tree_tops_EPSG25832.gpkg")

sf::st_write(
  obj          = ttops_out,
  dsn          = gpkg_path,
  layer        = "tree_tops",
  delete_layer = TRUE,
  quiet        = FALSE
)

cat("\n[OK] GeoPackage erfolgreich erstellt:\n", gpkg_path, "\n")
cat("   Felder: treeID, dgm_height, tree_height, canopy_height (= dgm_height + tree_height)\n")
