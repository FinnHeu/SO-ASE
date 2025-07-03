"""
Interpolate FESOM2 data on nodes to a regular lat/lon mesh

This script interpolates FESOM2 data on nodes to a regular latitude-longitude grid using
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
from os.path import isdir
from os import makedirs

# Get inputs from command line arguments
# Parse fixed arguments
year = sys.argv[1]
base_path_model = sys.argv[2]
base_path_mesh = sys.argv[3]

min_lat = float(sys.argv[4])
max_lat = float(sys.argv[5])
min_lon = float(sys.argv[6])
max_lon = float(sys.argv[7])
lat_increment = float(sys.argv[8])
lon_increment = float(sys.argv[9])

variables = sys.argv[10:]  # Remaining arguments are variable names

# Destination path for interpolated data
dest_path = base_path_model + 'interpolated/'

# Logging for confirmation
print(f"Processing year: {year}")
print(f"Base path model: {base_path_model}")
print(f"Base path mesh: {base_path_mesh}")
print(f"Base path output: {dest_path}")
print(f"Latitude bounds: {min_lat} to {max_lat} (step {lat_increment})")
print(f"Longitude bounds: {min_lon} to {max_lon} (step {lon_increment})")
print(f"Variables to interpolate: {variables}")

# Ensure the output directory exists
if not isdir(dest_path):
    print(f"Creating output directory: {dest_path}")
    makedirs(dest_path)

# Generate the lat/lon grid
lat = np.arange(min_lat, max_lat + lat_increment, lat_increment)
lon = np.arange(min_lon, max_lon + lon_increment, lon_increment)
lon_grid, lat_grid = np.meshgrid(lon, lat)

print(f"Generated grid: {len(lat)} latitudes x {len(lon)} longitudes")

# Load Mesh
mesh = pf.load_mesh(base_path_mesh, usepickle=False)
mesh_diag = xr.open_dataset(f"{base_path_model}fesom.mesh.diag.nc")

# Get number of vertical levels
nlevs = len(mesh_diag.nz1)

# Interpolate
for variable in variables:
    print(f"Interpolating variable: {variable}")

    # Load the variable data
    data_file = f"{base_path_model}{variable}.fesom.{year}.nc"
    ds = xr.open_dataset(data_file)

    # allocate array
    interp = np.zeros((len(lat), len(lon), nlevs,
                      ds.dims['time']), dtype=np.float32) * np.nan

    for i in range(len(ds.time)):
        print(
            f"Processing time step {i+1}/{len(ds.time)} for variable {variable}")
        for j in range(nlevs):
            print(
                f"Processing vertical level {j+1}/{nlevs} for variable {variable} at time step {i+1}")

            # Extract the variable data for the current time step and vertical level
            data = ds[variable].isel(time=i, nz1=j).values

            # Interpolate to regular grid
            interp_2D = pf.fesom2regular(
                data,
                mesh,
                lon_grid,
                lat_grid,
                how='nn',
                dumpfile=True
            )

            # Add vertical level and time dimension
            interp[:, :, j, i] = interp_2D[:, :]

    # Create a new xarray Dataset for the interpolated data)
    interp_ds = xr.Dataset(
        {
            variable: (['lat', 'lon', 'nz1', 'time'], interp)
        },
        coords={
            'lon': lon,
            'lat': lat,
            'nz1': mesh_diag.nz1,
            'time': ds.time
        }
    )

    # Save the interpolated data
    output_file = f"{dest_path}{variable}.interp.{min_lon}E.{max_lon}E.{lon_increment}E.{min_lat}N.{min_lat}N.{lat_increment}N.fesom.{year}.nc"
    interp_ds.to_netcdf(output_file)

    print(f"Saved interpolated data to: {output_file}")

print("Interpolation completed for all variables.")
