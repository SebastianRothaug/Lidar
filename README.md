# Lidar

This repository contains scripts for processing LiDAR point cloud data, mainly focusing on application involing vegetation analysis regarding height, structure, coverage.

The repository is structured du to different projects

## Project: powerline
**Clearance Analysis for Powerline Planning - Forest Impact and Vegetation Cut-Off Analysis in Powerline Corridor**

**Script structure:**
* **01_powerline_modeling.py** Modeling a powerline mathematically based on a symetrical parable calculation (2x pole location xyz + the zero point height; provided in CSV). 
* **02_corridor_tree_detection.R** 
* **03_powerline_tree_collision** 


## Requirements

These scripts are written in Python. Depending on the specific script, you may need standard spatial and point cloud libraries such as `laspy`, `pdal`, `rasterio`, and `scipy`.
