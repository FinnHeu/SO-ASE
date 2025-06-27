#!/usr/bin/env python3
"""
Unrotate FESOM2 Element Velocities

This script loads rotated velocity components (`u`, `v`) on elements from FESOM2 output,
unrotates them using `pyfesom2.vec_rotate_r2g`, and writes the unrotated
velocities (`u_unrot`, `v_unrot`) to NetCDF files.

The unrotation is performed for all time steps and vertical levels in the input data,
based on the average longitude and latitude of the element centers.

Author: Finn Heukamp
Date: June 2024

Usage:
    python unrotate_elem.py <YEAR>

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
    - Rotated velocity files for `u` and `v` for the specified year
    - Mesh diagnostics file with longitude and latitude information

Output:
    - Unrotated velocity files:
        - u_unrot.fesom.<YEAR>.nc
        - v_unrot.fesom.<YEAR>.nc

Notes:
    - This script is designed to be run in parallel for multiple years,
      ideally on a single node with sufficient CPUs.
"""
import sys
import xarray as xr
import pyfesom2 as pf
import numpy as np

# Get inputs from command line arguments
year = sys.argv[1] 
year = int(year)
print('Unrotating velocities for year: ' + str(year), flush=True)

base_path_model = sys.argv[2]
base_path_model = str(base_path_model)
print('Base path model: ' + base_path_model, flush=True)
print('Output is written to the same path as the input files.', flush=True)


# load mesh_diag
mesh_diag = xr.open_dataset(base_path_model + 'fesom.mesh.diag.nc')

# Files
src_file_u = base_path_model + 'u.fesom.' + str(year) + '.nc'
src_file_v = base_path_model + 'v.fesom.' + str(year) + '.nc'

# Load rotated u and v
ds_urot = xr.open_dataset(src_file_u)
ds_vrot = xr.open_dataset(src_file_v)

# Extract UV as numpy arrays
u_rot = ds_urot.u.values
v_rot = ds_vrot.v.values

# Extract LON, LAT (average 3 nods to center of element)
lons = mesh_diag.lon.values[mesh_diag.face_nodes.values-1].mean(axis=0)
lats = mesh_diag.lat.values[mesh_diag.face_nodes.values-1].mean(axis=0)

# Allocate
u_unrot_array = np.ones_like(u_rot) * np.nan
v_unrot_array = np.ones_like(u_rot) * np.nan

# Unrotate in Loop since pyfesom2.vec_rotate_r2g does only support 1D arrays
# Loop over time and depth (nz1) and unrotate u and v
for m in range(len(ds_urot.time)):
    print('Unrotating time step: ' + str(m+1) + '/' + str(len(ds_urot.time)), flush=True)
    for d in range(len(ds_urot.nz1)):

        u_unrot, v_unrot = pf.vec_rotate_r2g(50, 15, -90, lons, lats, u_rot[m,d,:], v_rot[m,d,:], flag=1)

        u_unrot_array[m, d, :] = u_unrot[np.newaxis, np.newaxis,:]
        v_unrot_array[m, d, :] = v_unrot[np.newaxis, np.newaxis,:]

# write to dataset
ds_urot['u_unrot'] = (('time','nz1','elem'), u_unrot_array)
ds_vrot['v_unrot'] = (('time','nz1','elem'), v_unrot_array)

# drop rotated velocities
ds_u_unrot = ds_urot.drop_vars('u')
ds_v_unrot = ds_vrot.drop_vars('v')

# to netcdf
dst_file_u = base_path_model + 'u_unrot.fesom.' + str(year) + '.nc'
dst_file_v = base_path_model + 'v_unrot.fesom.' + str(year) + '.nc'

ds_u_unrot.to_netcdf(dst_file_u)
ds_v_unrot.to_netcdf(dst_file_v)

print('Unrotated velocities saved to:', flush=True)
print(dst_file_u, flush=True)
print(dst_file_v, flush=True)

    
                           