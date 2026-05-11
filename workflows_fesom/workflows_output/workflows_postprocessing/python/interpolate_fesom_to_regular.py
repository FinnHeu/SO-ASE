"""
Interpolate FESOM2 data on nodes to a regular lat/lon mesh

This script interpolates FESOM2 data on nodes to a regular latitude-longitude grid using a modified
`pyfesom2.fesom2regular()`. It processes all time steps and vertical levels in the input data
for a given year and saves the interpolated output as NetCDF files.

Author: Finn Heukamp
Updated: June 2025

Usage:
    python interpolate_fesom_to_regular.py <YEAR> <BASE_PATH_MODEL> <BASE_PATH_MESH> \
        <MIN_LAT> <MAX_LAT> <MIN_LON> <MAX_LON> <LAT_INCREMENT> <LON_INCREMENT> <VARIABLES...>

Arguments:
    <YEAR>            Year to process (e.g., 1995)
    <BASE_PATH_MODEL> Base path where the model output files are located
    <BASE_PATH_MESH>  Path to the FESOM2 mesh directory
    <MIN_LAT>         Minimum latitude of the regular grid
    <MAX_LAT>         Maximum latitude of the regular grid
    <MIN_LON>         Minimum longitude of the regular grid
    <MAX_LON>         Maximum longitude of the regular grid
    <LAT_INCREMENT>   Latitude increment (e.g., 1.0)
    <LON_INCREMENT>   Longitude increment (e.g., 1.0)
    <VARIABLES...>    One or more variable names to interpolate (e.g., unod vnod temp salt)

Dependencies:
    - os
    - sys
    - numpy
    - xarray
    - pyfesom2

Input:
    - FESOM2 node-based NetCDF files for the specified year and variables
    - FESOM2 mesh and meshdiagnostics file (typically: 'fesom.mesh.diag.nc')

Output:
    - Interpolated NetCDF files named:
        <variable>.interp.<lon-range>E.<lat-range>N.fesom.<year>.nc

Notes:
    - Designed for parallel execution (e.g., over multiple years)
    - Assumes data follows the naming pattern: <variable>.fesom.<year>.nc
"""

import sys
import xarray as xr
import pyfesom2 as pf
import numpy as np
import os
import joblib
import scipy
from os.path import isdir, isfile
from os import makedirs
import so_ase as so

# Get inputs from command line arguments
# Parse fixed arguments
year = int(sys.argv[1])
base_path_model = sys.argv[2]
base_path_mesh = sys.argv[3]

min_lat = float(sys.argv[4])
max_lat = float(sys.argv[5])
min_lon = float(sys.argv[6])
max_lon = float(sys.argv[7])
lat_increment = float(sys.argv[8])
lon_increment = float(sys.argv[9])
radius_of_influence = float(sys.argv[10])

# Destination path for interpolated data
dest_path = sys.argv[11] + "/gridded/"

raw_variables = sys.argv[12:]  # Remaining arguments are variable names
variables = [
    v.strip()
    for part in raw_variables
    for v in part.split(",")
    if v.strip()
]

# Print a nice header
header = f"""
{'='*78}
{' '*20}FESOM2 to Regular Grid Interpolation{' '*20}
{'='*78}
Year: {year}
Model path: {base_path_model}
Mesh path: {base_path_mesh}
Output path: {dest_path}
Grid bounds: [{min_lat}, {max_lat}] latitude, [{min_lon}, {max_lon}] longitude
Grid resolution: {lat_increment}° x {lon_increment}°
Radius of influence: {radius_of_influence/1000:.1f} km
Variables: {', '.join(variables)}
{'='*78}
"""
print(header)

# Ensure the output directory exists
if not isdir(dest_path):
    print(f"Creating output directory: {dest_path}")
    makedirs(dest_path, exist_ok=True)

so.fesom2regular_nd(
    base_path_model, 
    base_path_mesh, 
    reg_lat=(min_lat, max_lat), 
    reg_lon=(min_lon, max_lon), 
    lat_increment=lat_increment, 
    lon_increment=lon_increment,
    variables=variables, 
    years=(year, year), 
    radius_of_influence=100000, 
    dest_path=dest_path,
    dumpfile=True,
    log=True
)

print("Interpolation completed for all variables.")
