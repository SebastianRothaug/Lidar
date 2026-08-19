import numpy as np
import pandas as pd
import fiona
from fiona.crs import from_epsg
from shapely.geometry import LineString, mapping
from pathlib import Path

import rasterio
from rasterio.transform import from_origin, rowcol, xy as transform_xy
from rasterio.features import geometry_mask
from scipy.interpolate import griddata

import laspy

####  Erforderliche Input-Daten: Mast- bzw Spannfeldliste in bestimmter Struktur:
    ### Mast1: x0, y0, z0 // Mast2: x2, y2, z2 // minimaler Durchhangpunkt: z1 (= Höhe minimaler Durchang: Dieser Wert wird von Bauingenieuren/Statikern errechnet, er lässt sich durch die Höhe der beiden Masten + angelegte Kabelspannung bestimmen)
    ## Durch diese zwei 3D-Mastkoordinaten und den Durchhangpunkt (entspicht Nullpunkt der Parabel) kann die Parabel in 3D eindeutig bestimmt und durch Punkte modelliert werden.

# --------------------------------------------------------------------------
# Input Definition
# --------------------------------------------------------------------------
CSV_PATH = Path.cwd() / "i01_powerline_location.csv"  # Anpassen
if not CSV_PATH.is_file():
    raise FileNotFoundError(f"CSV-Datei nicht gefunden: {CSV_PATH}")

df = pd.read_csv(CSV_PATH, sep=';', dtype={'group_id': str})
required_cols = {'group_id', 'x0', 'y0', 'z0', 'x2', 'y2', 'z2', 'z1'}
if not required_cols.issubset(df.columns):
    raise ValueError(f"CSV muss die Spalten {required_cols} enthalten.")


# --------------------------------------------------------------------------
# Output Definition
# --------------------------------------------------------------------------
OUT_DIR = Path.cwd() / "output"
OUT_DIR.mkdir(parents=True, exist_ok=True)

GPKG_PATH = OUT_DIR / "01_leitungsdurchhang_points.gpkg" # Leitungsparabel (repräsentiert durch punkte) einfach zum einladen im GIS Geopackage mit Z-Wert
LAS_PATH = OUT_DIR / "01_leitungsdurchhang_korridor_points.las" # Leitungsparabel auf Korridor(/Schutzstreifen) Breite direkt als .las file
TIFF_PATH = OUT_DIR / "01_leitungsdurchhang_hoehenmodell.tif" # Leitungsparabel als Höhenmodell = Rasterisierung der Korridor Points

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------
CRS_EPSG = 25832

# --- Breite des Korridors ---------------------------------------------------
CORRIDOR_HALFWIDTH_M = 40.0   # in der Regel Schutzstreifenbreite 40m
RASTER_RESOLUTION_M = 1.0     # Pixelgroesse des Ausgabe-TIFFs

# --- Intervall der erzeugten Punktwolke innerhalb des Korridors ---------------------------------------------------
INTERVAL_M = 1.0              # Abstand der Stuetzpunkte entlang der Linie (in Linien Richtung = von Mast zu Mast)
CROSS_INTERVAL_M = 1.0        # Abstand der Stuetzpunkte quer zur Linie (Breite also rechtwinkling zur Linie; repräsentiert Schutzstreifen)




# --------------------------------------------------------------------------
# Grundfunktionen (unveraendert aus dem Ursprungsskript)
# --------------------------------------------------------------------------
def interpolate_mid_point(p0, p2):
    x_mid = (p0[0] + p2[0]) / 2.0
    y_mid = (p0[1] + p2[1]) / 2.0
    return x_mid, y_mid


def fit_quadratic(t_vals, coord_vals):
    return np.polyfit(t_vals, coord_vals, 2)


def planar_arc_length(coeffs_x, coeffs_y, n_samples=100_000):
    t_dense = np.linspace(0.0, 2.0, n_samples)

    der_x = np.polyder(coeffs_x)
    der_y = np.polyder(coeffs_y)

    dx = np.polyval(der_x, t_dense)
    dy = np.polyval(der_y, t_dense)

    ds_dt = np.sqrt(dx**2 + dy**2)
    s2d_dense = np.cumsum(ds_dt) * (t_dense[1] - t_dense[0])
    s2d_dense -= s2d_dense[0]
    return t_dense, s2d_dense


# --------------------------------------------------------------------------
# Punkte entlang der Parabel in 3D erzeugen (inkl. Tangenten- und Normalenrichtung, wichtig fuer die Korridor-Breite)
def points_and_tangents_at_interval(t_dense, s2d_dense, coeffs_xyz, interval=1.0):
    """Wie points_at_interval im Ursprungsskript, gibt zusaetzlich die
    (normierte) Tangenten- und Normalenrichtung je Punkt zurueck."""
    total_len = s2d_dense[-1]
    target_dists = np.arange(0.0, np.floor(total_len) + interval, interval)

    t_at_dist = np.interp(target_dists, s2d_dense, t_dense)

    x = np.polyval(coeffs_xyz['x'], t_at_dist)
    y = np.polyval(coeffs_xyz['y'], t_at_dist)
    z = np.polyval(coeffs_xyz['z'], t_at_dist)

    der_x = np.polyder(coeffs_xyz['x'])
    der_y = np.polyder(coeffs_xyz['y'])
    dx = np.polyval(der_x, t_at_dist)
    dy = np.polyval(der_y, t_at_dist)
    norm = np.sqrt(dx**2 + dy**2)
    norm[norm == 0] = 1e-9  # Division durch 0 vermeiden

    tangent_x, tangent_y = dx / norm, dy / norm
    # Normale = Tangente um 90 Grad gedreht
    normal_x, normal_y = -tangent_y, tangent_x

    return list(zip(x, y, z, target_dists, normal_x, normal_y))


# --------------------------------------------------------------------------
# Korridorpunkte erzeugen (seitlich aufgeweitete Punktwolke)
def build_corridor_points(centerline_pts, halfwidth, cross_step):
    """Erzeugt aus den Mittellinienpunkten (mit Normalenrichtung) die
    seitlich aufgeweitete Punktwolke. Der Z-Wert wird ueber die Breite
    konstant gehalten (Extrusion der Durchhangslinie)."""
    n_cross = int(round((2 * halfwidth) / cross_step)) + 1
    offsets = np.linspace(-halfwidth, halfwidth, n_cross)

    xs, ys, zs = [], [], []
    for (x, y, z, dist, nx, ny) in centerline_pts:
        px = x + offsets * nx
        py = y + offsets * ny
        pz = np.full_like(offsets, z)
        xs.append(px)
        ys.append(py)
        zs.append(pz)

    return np.concatenate(xs), np.concatenate(ys), np.concatenate(zs)


def write_las(xs, ys, zs, group_index, group_id_lookup, path, crs_epsg):
    """Schreibt die gesamte Korridor-Punktwolke (alle Punkte, aus denen das
    Geländemodell interpoliert wurde) als georeferenzierte LAS-Datei.
    group_index: int-Array je Punkt -> Index in group_id_lookup (fuer die
    urspruengliche GROUP_ID als String, da LAS keine String-Extra-Dims kann)."""
    header = laspy.LasHeader(point_format=3, version="1.2")
    header.offsets = [float(np.min(xs)), float(np.min(ys)), float(np.min(zs))]
    header.scales = [0.001, 0.001, 0.001]

    try:
        import pyproj
        header.add_crs(pyproj.CRS.from_epsg(crs_epsg))
    except Exception as e:
        print(f"[LAS] Warnung: CRS konnte nicht in den LAS-Header geschrieben werden ({e})")

    header.add_extra_dim(laspy.ExtraBytesParams(name="group_idx", type=np.uint16))

    las = laspy.LasData(header)
    las.x = xs
    las.y = ys
    las.z = zs
    las.group_idx = group_index.astype(np.uint16)
    las.write(str(path))

    print(f"[LAS] {len(xs)} Punkte -> {path}")
    print(f"[LAS] group_idx -> GROUP_ID Zuordnung: {group_id_lookup}")


def process_group(group_row):
    gid = group_row['group_id']

    p0 = np.array([group_row['x0'], group_row['y0'], group_row['z0']])
    p2 = np.array([group_row['x2'], group_row['y2'], group_row['z2']])

    x_mid, y_mid = interpolate_mid_point(p0, p2)
    p1 = np.array([x_mid, y_mid, group_row['z1']])

    P = np.vstack([p0, p1, p2])
    t_vals = np.arange(3)  # 0,1,2
    coeffs_xyz = {
        'x': fit_quadratic(t_vals, P[:, 0]),
        'y': fit_quadratic(t_vals, P[:, 1]),
        'z': fit_quadratic(t_vals, P[:, 2]),
    }

    t_dense, s2d_dense = planar_arc_length(coeffs_xyz['x'], coeffs_xyz['y'])
    centerline_pts = points_and_tangents_at_interval(
        t_dense, s2d_dense, coeffs_xyz, interval=INTERVAL_M
    )

    # Fiona-Features fuer das GeoPackage (nur Mittellinie, wie bisher)
    features = []
    for i, (xp, yp, zp, dist, nx, ny) in enumerate(centerline_pts):
        geom = {"type": "Point", "coordinates": (float(xp), float(yp), float(zp))}
        feat = {
            "geometry": geom,
            "properties": {
                "ID": int(i),
                "GROUP_ID": str(gid),
                "DIST_2D": float(dist),
                "X": float(xp),
                "Y": float(yp),
                "Z": float(zp),
            },
        }
        features.append(feat)

    return gid, centerline_pts, features


# --------------------------------------------------------------------------
# GeoPackage-Schema (Mittellinie, unveraendert)
# --------------------------------------------------------------------------
schema = {
    "geometry": "Point",
    "properties": {
        "ID": "int",
        "GROUP_ID": "str",
        "DIST_2D": "float",
        "X": "float",
        "Y": "float",
        "Z": "float",
    },
}


# --------------------------------------------------------------------------
# Hauptlauf
# --------------------------------------------------------------------------
def main():
    all_groups = []  # (gid, centerline_pts) fuer die Korridor-Berechnung

    # ---- 1) Mittellinien-Punkte je Gruppe berechnen & ins GPKG schreiben ----
    with fiona.open(
        GPKG_PATH,
        mode='w',
        driver='GPKG',
        crs=from_epsg(CRS_EPSG),
        schema=schema,
        layer='points_1m',
    ) as dst:
        for idx, row in df.iterrows():
            gid, centerline_pts, features = process_group(row)
            all_groups.append((gid, centerline_pts))
            for f in features:
                dst.write(f)

    print(f"[GPKG] {len(df)} Gruppen -> {GPKG_PATH}")

    # ---- 2) Korridor-Punktwolke je Gruppe erzeugen --------------------------
    corridor_data = []  # (gid, xs, ys, zs, buffer_polygon)
    all_x, all_y = [], []

    # Fuer den LAS-Export: alle Korridorpunkte aller Gruppen inkl. Gruppen-Index
    las_xs, las_ys, las_zs, las_group_idx = [], [], [], []
    group_id_lookup = {}

    for group_index, (gid, centerline_pts) in enumerate(all_groups):
        xs, ys, zs = build_corridor_points(
            centerline_pts, CORRIDOR_HALFWIDTH_M, CROSS_INTERVAL_M
        )
        line = LineString([(p[0], p[1]) for p in centerline_pts])
        poly = line.buffer(CORRIDOR_HALFWIDTH_M, cap_style="flat")

        corridor_data.append((gid, xs, ys, zs, poly))
        all_x.append(xs)
        all_y.append(ys)

        las_xs.append(xs)
        las_ys.append(ys)
        las_zs.append(zs)
        las_group_idx.append(np.full(xs.shape, group_index))
        group_id_lookup[group_index] = str(gid)

    all_x = np.concatenate(all_x)
    all_y = np.concatenate(all_y)

    las_xs = np.concatenate(las_xs)
    las_ys = np.concatenate(las_ys)
    las_zs = np.concatenate(las_zs)
    las_group_idx = np.concatenate(las_group_idx)

    # ---- 3) Gemeinsames Raster (Ausgabe-Grid) festlegen ---------------------
    margin = RASTER_RESOLUTION_M * 2
    xmin, xmax = all_x.min() - margin, all_x.max() + margin
    ymin, ymax = all_y.min() - margin, all_y.max() + margin

    width_px = int(np.ceil((xmax - xmin) / RASTER_RESOLUTION_M))
    height_px = int(np.ceil((ymax - ymin) / RASTER_RESOLUTION_M))

    transform = from_origin(xmin, ymax, RASTER_RESOLUTION_M, RASTER_RESOLUTION_M)
    raster = np.full((height_px, width_px), np.nan, dtype="float32")

    # ---- 4) Je Gruppe interpolieren und ins Gesamtraster einsetzen ----------
    for gid, xs, ys, zs, poly in corridor_data:
        gx_min, gx_max = xs.min(), xs.max()
        gy_min, gy_max = ys.min(), ys.max()

        row_max, col_min = rowcol(transform, gx_min, gy_min)
        row_min, col_max = rowcol(transform, gx_max, gy_max)

        row_min = max(int(row_min) - 1, 0)
        row_max = min(int(row_max) + 1, height_px)
        col_min = max(int(col_min) - 1, 0)
        col_max = min(int(col_max) + 1, width_px)

        if row_min >= row_max or col_min >= col_max:
            continue

        sub_rows = np.arange(row_min, row_max)
        sub_cols = np.arange(col_min, col_max)
        cc, rr = np.meshgrid(sub_cols, sub_rows)
        mesh_x, mesh_y = transform_xy(transform, rr, cc)
        # rasterio.transform.xy liefert bei Array-Eingabe flache Listen zurueck
        # -> auf die urspruengliche 2D-Grid-Form (rows, cols) zurueckformen
        mesh_x = np.array(mesh_x).reshape(rr.shape)
        mesh_y = np.array(mesh_y).reshape(rr.shape)

        interp_vals = griddata(
            points=np.column_stack([xs, ys]),
            values=zs,
            xi=(mesh_x, mesh_y),
            method="linear",
        )

        # Auf das tatsaechliche Korridor-Polygon maskieren (keine
        # rechteckigen Kunstartefakte an den Trassen-Enden)
        window_transform = rasterio.transform.from_origin(
            *transform_xy(transform, row_min, col_min, offset="ul"),
            RASTER_RESOLUTION_M,
            RASTER_RESOLUTION_M,
        )
        mask_outside = geometry_mask(
            [mapping(poly)],
            out_shape=(row_max - row_min, col_max - col_min),
            transform=window_transform,
            invert=False,  # True = innerhalb Polygon markieren -> hier: False liefert True fuer Aussen
        )

        valid = (~mask_outside) & (~np.isnan(interp_vals))
        sub_raster = raster[row_min:row_max, col_min:col_max]
        sub_raster[valid] = interp_vals[valid]
        raster[row_min:row_max, col_min:col_max] = sub_raster

    # ---- 5) GeoTIFF schreiben ------------------------------------------------
    with rasterio.open(
        TIFF_PATH,
        "w",
        driver="GTiff",
        height=height_px,
        width=width_px,
        count=1,
        dtype="float32",
        crs=f"EPSG:{CRS_EPSG}",
        transform=transform,
        nodata=np.nan,
        compress="deflate",
    ) as dst:
        dst.write(raster, 1)

    print(f"[TIFF] Korridorbreite gesamt: {2 * CORRIDOR_HALFWIDTH_M:.0f} m")
    print(f"[TIFF] Aufloesung: {RASTER_RESOLUTION_M} m/Pixel, Groesse: {width_px} x {height_px} px")
    print(f"[TIFF] Geschrieben nach: {TIFF_PATH}")

    # ---- 6) Korridor-Punktwolke als LAS schreiben ----------------------------
    write_las(las_xs, las_ys, las_zs, las_group_idx, group_id_lookup, LAS_PATH, CRS_EPSG)

    print("\n✅  Fertig!")


if __name__ == "__main__":
    main()
