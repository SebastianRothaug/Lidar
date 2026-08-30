# Lidar

This repository contains scripts for processing LiDAR point cloud data, mainly focusing on applications involving vegetation analysis.

The repository is structured according to different projects.

## Project: powerline
**Clearance Analysis for Powerline Planning - Forest Impact and Vegetation Cut-Off Analysis in Powerline Corridor** <br>
An automated workflow combining mathematical powerline modeling and LiDAR-based tree detection to assess vegetation risks and calculate clearance requirements in utility corridors. By comparing individual tree heights against simulated cable sags, the project identifies critical collision zones and calculates the exact cut-off length or max. height needed for each tree to ensure infrastructure safety.

### My example data and project:
* The LiDAR data I used in my example project is provided by "Bayerische Vermessungsverwaltung – www.geodaten.bayern.de" and is openly accessible.
* **The project is fictionally created!** Neither is the AOI location involved in a real planning process, nor are administrative and legal factors (such as the Federal Land Utilization Ordinance) taken into account in this infrastructure planning scenario.

### Script structure:
* **01_powerline_modeling.py**
  * Modeling a powerline mathematically based on a symmetrical parabola calculation (2x pole location xyz + the zero point height; provided in an initial CSV).
    <br>This results in an output of a powerline height model representing the electricity cable.
* **02_corridor_tree_detection.R**
  * Preprocessing of the point cloud - clipping the powerline corridor to minimize data size and increase computing performance.
  * Point cloud processing - classification and rasterization
    * Ground = Digital Terrain Model (DTM) derived via Cloth-Simulation-Filter (if no previous classification information is available).
    * Canopy = Canopy Height Model (CHM) derived from normalized height calculation: original_las_height - DTM.
  * Tree detection based on identifying the maxima of the Canopy Height Model - using a moving window approach with a dynamic window size for heterogeneous forest density and structure.
    <br>This results in an output of a tree layer involving tree height information within the corridor.
* **03_powerline_tree_collision**
  * Identification of critical collisions, cut-off length, and max. tree height - by comparing the tree detection layer height information (from 02) to the powerline height model (from 01). 

### Data folder structure:
* **Initial (external) data**, to initiate the scripts - marked with: **i** + the script number it needs to be implemented in.
  * **i01_powerline_location.csv** - CSV table of the pole locations (dt. Mast- bzw. Spannfeldliste). *Each row represents a span between two poles with 2x pole location in xyz + the zero point height. Data and structure/nomenclature are provided by the structural engineers (they calculate the zero point height based on material variables of the cable, as well as taking into account physical factors like tractive power and span power).*
  * **i02_lidar_full_aoi.las** - LiDAR point cloud of the area of interest.
    <br>**!** *Due to data size issues, this specific example dataset couldn't be uploaded to the repository (contact me or alternatively continue with the already clipped AOI file: 02_lidar_aoi_clipped.las).* 
  * **i02_corridor_line.gpkg** - Line layer representing the base/middle line of the corridor for tree detection.

* Processed (internal) data, created by the scripts - named with the script number it is exported from.
  * Some of this internally processed data is not a final product, but is required to initiate the next scripts.
