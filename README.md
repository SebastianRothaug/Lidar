# Lidar

This repository contains scripts for processing LiDAR point cloud data, mainly focusing on application involing vegetation analysis regarding height, structure, coverage.

The repository is structured du to different projects

## Project: powerline
**Clearance Analysis for Powerline Planning - Forest Impact and Vegetation Cut-Off Analysis in Powerline Corridor**
Description Text

**Script structure:**
* **01_powerline_modeling.py**
  * Modeling a powerline mathematically based on a symetrical parable calculation (2x pole location xyz + the zero point height; provided in a initial CSV).
    Resulting in an output of a powerline height model representing the electricity cable.
* **02_corridor_tree_detection.R**
  * Preprocessing of the poingt cloud - clipping the powerline corridor to minimize data size and increase computing performance.
  * Point cloud processing - classification and rasterisation
    * Ground = Digital Terreain Model (DTM) derived via Cloth-Simulation-Filter (if no previous classification information)
    * Canopy = Canopy Height Model (CHM) derived from hormalized height calculation: original_las_height - DTM
  * Tree detection based on identifying the maxima of the Canopy Height Model - using a moving window approach with a dynamic window size for hetereogenous forest density and structure.
    Resulting in an output of a tree layer involving tree height informations within the corridor.
* **03_powerline_tree_collision**
  * Identification of critical collision, cut of length and max. tree height - by comparing the tree detection layer heights information (from 02) to the powerline height model (from 01). 

**Data folder structure:**
* **Initial (external) data**, to initiate the scripts - marked with: **i** + the script number needs to be implemented.
  * **i01_powerline_location.csv** - CSV table of the pole locations (dt. Mast- bzw Spannfeldliste). *Each row representing a span between to poles with 2x pole location in xyz + the zero point height. Data and structure/nomenclature provided by the structural engineers (they calculate the zero point heiht based on materialistic varibels of the cable just as taking into account the physical factors like tractive power and span power)
  * **i02_lidar_full_aoi.las** - LiDaR point cloud of the area of interest.  **!! due to data size issues that specific example dataset couldn't be uloaded to the repository (contact me or alternatively continue with the already clipped AOI file: 02_lidar_aoi_clipped.las !!**
  * **i02_corridor_line.gpkg** - Line layer representing the base/middle line of the corridor for the tree detection.

* Processed (internal) data, created by the scripts - named with the script number it is exported from.
  * Some of these internal processed data is not a final product, but required to initiate the next scripts.

**My example data and project:**
* The data LiDaR data I used in my example project is provided by "Bayerische Vermessungsverwaltung – www.geodaten.bayern.de" and openly accessibly.
* **The project is fictional created!** neither the AOI location is in a real planing involved nor adminatrative and legal factors (just as the Federal Land Utilization Ordinance) are taken into account.

## Requirements

These scripts are written in Python. Depending on the specific script, you may need standard spatial and point cloud libraries such as `laspy`, `pdal`, `rasterio`, and `scipy`.
