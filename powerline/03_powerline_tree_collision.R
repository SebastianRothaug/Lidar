# collision detection
library(terra)
library(sf)
library(ggplot2)
library(tidyterra)   # fuer geom_spatraster() -> SpatRaster in ggplot darstellen

# ---- Eingabe-Pfade (bitte anpassen) ---------------------------------------
tiff_path <- "C:/Users/User/Desktop/Python_VS/output/leitungsdurchhang_hoehenmodell.tif"
gpkg_path <- "C:/Users/User/Desktop/Lidar_Course/tree_tops_EPSG25832.gpkg"
gpkg_layer <- NULL   # z.B. "tree_tops", falls mehrere Layer im GPKG vorhanden sind

out_gpkg_path <- file.path(
  dirname(gpkg_path),
  paste0(tools::file_path_sans_ext(basename(gpkg_path)), "_collision.gpkg")
)

# ---- 1) Baumpunkte einlesen -------------------------------------------------
trees <- if (!is.null(gpkg_layer)) {
  st_read(gpkg_path, layer = gpkg_layer, quiet = FALSE)
} else {
  st_read(gpkg_path, quiet = FALSE)
}

# --- Karte: eingelesene Baumpunkte -----------------------------------------
p_trees <- ggplot() +
  geom_sf(data = trees, color = "darkgreen", size = 1, alpha = 0.7) +
  labs(title = "Eingelesene Baumstandorte", x = NULL, y = NULL) +
  theme_minimal()
print(p_trees)

# ---- 2) Leitungshoehenmodell (TIFF) oeffnen --------------------------------
dtm_leitung <- terra::rast(tiff_path)

# --- Karte: Leitungshoehenmodell --------------------------------------------
p_dtm <- ggplot() +
  geom_spatraster(data = dtm_leitung) +
  scale_fill_viridis_c(name = "Hoehe (m)", na.value = NA) +
  labs(title = "Leitungsdurchhangshoehenmodell", x = NULL, y = NULL) +
  theme_minimal()
print(p_dtm)

# CRS angleichen, falls noetig
raster_crs <- terra::crs(dtm_leitung, describe = TRUE)$code
trees_epsg <- sf::st_crs(trees)$epsg


# ---- 3) Leitungshoehe an den Baumpositionen auslesen -----------------------
trees_vect <- terra::vect(trees)
sampled <- terra::extract(dtm_leitung, trees_vect)
leitung_height <- sampled[[2]]  # 2. Spalte = Rasterwert

trees$leitung_height <- leitung_height

# ---- 4) Kollisionslogik ----------------------------------------------------
canopy <- trees$canopy_height
tree_h <- trees$tree_height
line_h <- trees$leitung_height

no_data_mask <- is.na(line_h) | is.na(canopy)
collision_bool <- (canopy >= line_h) & !no_data_mask

diff <- canopy - line_h  # positiv = Baum ragt in/ueber die Leitung hinein
collision_max_tree_height <- tree_h - diff

# Nur dort befuellen, wo tatsaechlich eine Kollision vorliegt
diff[!collision_bool] <- NA
collision_max_tree_height[!collision_bool] <- NA

collision_str <- ifelse(
  no_data_mask, "no data",
  ifelse(canopy >= line_h, "yes", "no")
)

trees$collision <- collision_str
trees$collision_diff_m <- diff
trees$collision_max_tree_height <- collision_max_tree_height

# ---- 5) Ergebnis schreiben --------------------------------------------------
sf::st_write(
  obj          = trees,
  dsn          = out_gpkg_path,
  layer        = "tree_tops_collision",
  delete_layer = TRUE,
  quiet        = FALSE
)

n_total    <- nrow(trees)
n_collision <- sum(collision_str == "yes")
n_ok        <- sum(collision_str == "no")
n_nodata    <- sum(collision_str == "no data")

cat(sprintf("[OK] %d Baeume geprueft.\n", n_total))
cat(sprintf("     Kollisionen (ja): %d\n", n_collision))
cat(sprintf("     Keine Kollision (nein): %d\n", n_ok))
cat(sprintf("     Ausserhalb des Leitungsmodells (kein Wert): %d\n", n_nodata))
cat(sprintf("[OK] Ergebnis geschrieben nach: %s\n", out_gpkg_path))

# --------------------------------------------------------------------------
# Kartenausgabe: Kollisionsstatus je Baum (yes = rot, no = gruen, no data = blau)
# --------------------------------------------------------------------------
collision_colors <- c("yes" = "red", "no" = "green", "no data" = "blue")

p_collision <- ggplot() +
  geom_sf(data = trees, aes(color = collision), size = 1.5, alpha = 0.8) +
  scale_color_manual(values = collision_colors, name = "Kollision") +
  labs(title = "Kollisionsstatus der Baeume", x = NULL, y = NULL) +
  theme_minimal()
print(p_collision)

# --------------------------------------------------------------------------
# Kurzes Balkendiagramm: Anzahl Baeume je Kollisionsstatus
# --------------------------------------------------------------------------
p_bar <- ggplot(trees, aes(x = collision, fill = collision)) +
  geom_bar() +
  scale_fill_manual(values = collision_colors, name = "Kollision") +
  labs(title = "Anzahl Baeume je Kollisionsstatus", x = NULL, y = "Anzahl") +
  theme_minimal()
print(p_bar)
