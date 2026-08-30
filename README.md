# Lidar

This repository contains scripts for processing LiDAR point cloud data, mainly focusing on application involing vegetation analysis regarding height, structure, coverage.

The repository is structured du to different projects

## Project: powerline
**Clearance Analysis for Powerline Planning - Forest Impact and Vegetation Cut-Off Analysis in Powerline Corridor**

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


## Requirements

These scripts are written in Python. Depending on the specific script, you may need standard spatial and point cloud libraries such as `laspy`, `pdal`, `rasterio`, and `scipy`.
