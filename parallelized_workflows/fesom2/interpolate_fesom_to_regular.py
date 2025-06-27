#!/usr/bin/env python3
"""
Interpolate FESOM2 data on nodes to regular mesh

This script interpolates FESOM2 data on nodes to a regular lat/lon mesh using `pyfesom2.fesom2regular()`.
The interpolation is performed for all time steps and vertical levels in the input data.

Author: Finn Heukamp
Date: June 2024

Usage:
    python interpolate_fesom_to_regular.py <YEAR> <BASE_PATH_MODEL>

Arguments:
    <YEAR>  Year to process (e.g., 1995). Must match available input files:
            - u.fesom.<YEAR>.nc
            - v.fesom.<YEAR>.nc
    
    <BASE_PATH_MODEL>  Base path where the model output files are located.

Dependencies:
    - xarray
    - numpy
    - pyfesom2
    - netCDF4 (for xarray backend)
    - FESOM2 mesh diagnostics file: 'fesom.mesh.diag.nc'

Input:
    - FESOM files on nodes for the specified year
    - Mesh diagnostics file with longitude and latitude information

Output:
    - Interpolated files:
        - <variable>.interp.fesom.<YEAR>.nc

Notes:
    - This script is designed to be run in parallel for multiple years,
      ideally on a single node with sufficient CPUs.
"""
import sys
import xarray as xr
import pyfesom2 as pf
import numpy as np

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


# Logging for confirmation
print(f"Processing year: {year}")
print(f"Base path model: {base_path_model}")
print(f"Base path mesh: {base_path_mesh}")
print(f"Latitude bounds: {min_lat} to {max_lat} (step {lat_increment})")
print(f"Longitude bounds: {min_lon} to {max_lon} (step {lon_increment})")
print(f"Variables to interpolate: {variables}")

# Example: Generate the lat/lon grid
lat_grid = np.arange(min_lat, max_lat + lat_increment, lat_increment)
lon_grid = np.arange(min_lon, max_lon + lon_increment, lon_increment)

print(f"Generated grid: {len(lat_grid)} latitudes x {len(lon_grid)} longitudes")

# Interpolate each variable
mesh = pf.load_mesh(base_path_mesh)
mesh_diag = xr.open_dataset(f"{base_path_model}fesom.mesh.diag.nc")

# Get number of vertical levels
nlevs = len(mesh_diag.nz1)

for variable in variables:
    print(f"Interpolating variable: {variable}")
    
    # Load the variable data
    data_file = f"{base_path_model}{variable}.fesom.{year}.nc"
    ds = xr.open_dataset(data_file)

    # allocate array
    interp = np.zeros((len(lon_grid), len(lat_grid), nlevs, ds.dims['time']), dtype=np.float32) * np.nan
    
    for i in range(len(ds.time)):
        print(f"Processing time step {i+1}/{len(ds.time)} for variable {variable}")
        for j in range(nlevs):
            print(f"Processing vertical level {j+1}/{nlevs} for variable {variable} at time step {i+1}")
            
            # Extract the variable data for the current time step and vertical level
            data = ds.isel(time=i, nz1=j).values.squeeze()
            
            # Interpolate to regular grid
            interp_2D = pf.fesom2regular(
                data,
                mesh,
                lon_grid,
                lat_grid,
                how='nn'
            )

            interp[:, :, j, i] = interp_2D[:,:, np.newaxis, np.newaxis]  # Add vertical level and time dimension
        
    # Create a new xarray Dataset for the interpolated data)
    interp_ds = xr.Dataset(
        {
            variable: (['lon', 'lat', 'nz1', 'time'], interp)
        },
        coords={
            'lon': lon_grid,
            'lat': lat_grid,
            'nz1': mesh_diag.nz1,
            'time': ds.time
        }
    )
    
    # Save the interpolated data
    output_file = f"{base_path_model}{variable}.interp.{min_lon}E.{max_lon}E.{lon_increment}E.{min_lat}N.{min_lat}N.{lat_increment}N.fesom.{year}.nc"
    interp_ds.to_netcdf(output_file)
    
    print(f"Saved interpolated data to: {output_file}")

print("Interpolation completed for all variables.")


    
                           