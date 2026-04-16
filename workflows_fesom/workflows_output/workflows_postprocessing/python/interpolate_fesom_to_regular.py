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
from scipy.spatial import cKDTree

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
radius_of_influence = float(sys.argv[10])

raw_variables = sys.argv[11:]  # Remaining arguments are variable names
variables = [
    v.strip()
    for part in raw_variables
    for v in part.split(",")
    if v.strip()
]

# Destination path for interpolated data
dest_path = base_path_model + "/gridded/"

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

# Generate the lat/lon grid
lat = np.arange(min_lat, max_lat + lat_increment, lat_increment)
lon = np.arange(min_lon, max_lon + lon_increment, lon_increment)
lon_grid, lat_grid = np.meshgrid(lon, lat)

print(f"Generated grid: {len(lat)} latitudes x {len(lon)} longitudes")

# Load Mesh
mesh = pf.load_mesh(base_path_mesh, usepickle=False)
mesh_diag = xr.open_dataset(f"{base_path_mesh}fesom.mesh.diag.nc")

# Define time encoding
time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)

# Get number of vertical levels
nlevs = len(mesh_diag.nz1)

# Define helper functions
def lon_lat_to_cartesian(lon, lat, R=6371000):
    """
    calculates lon, lat coordinates of a point on a sphere with
    radius R. Taken from http://earthpy.org/interpolation_between_grids_with_ckdtree.html
    """
    lon_r = np.radians(lon)
    lat_r = np.radians(lat)

    x = R * np.cos(lat_r) * np.cos(lon_r)
    y = R * np.cos(lat_r) * np.sin(lon_r)
    z = R * np.sin(lat_r)
    return x, y, z

def create_indexes_and_distances(mesh, lons, lats, k=1):
    """
    Creates KDTree object and query it for indexes of points in FESOM mesh that are close to the
    points of the target grid. Also return distances of the original points to target points.

    Parameters
    ----------
    mesh : fesom_mesh object
        pyfesom mesh representation
    lons/lats : array
        2d arrays with target grid values.
    k : int
        k-th nearest neighbors to return.
    n_jobs : int, optional
        Number of jobs to schedule for parallel processing. If -1 is given
        all processors are used. Default: 1.

    Returns
    -------
    distances : array of floats
        The distances to the nearest neighbors.
    inds : ndarray of ints
        The locations of the neighbors in data.

    """
    xs, ys, zs = lon_lat_to_cartesian(mesh.x2, mesh.y2)
    xt, yt, zt = lon_lat_to_cartesian(lons.flatten(), lats.flatten())

    tree = cKDTree(list(zip(xs, ys, zs)))

    border_version = "1.6.0"
    current_version = scipy.__version__
    v1_parts = list(map(int, border_version.split(".")))
    v2_parts = list(map(int, current_version.split(".")))

    if v2_parts > v1_parts:
        distances, inds = tree.query(list(zip(xt, yt, zt)), k=k)
    else:
        distances, inds = tree.query(list(zip(xt, yt, zt)), k=k)

    return distances, inds

def fesom2regular_nn(
    data,
    mesh,
    lons,
    lats,
    distances_path=None,
    inds_path=None,
    radius_of_influence=100000,
    dumpfile=False,
):
    """
    Nearest-neighbor interpolation from FESOM mesh to target grid.
    """

    left, right = np.min(lons), np.max(lons)
    down, up = np.min(lats), np.max(lats)
    lonNumber, latNumber = lons.shape[1], lats.shape[0]

    kk = 1  # nearest neighbor

    MESH_BASE = os.path.basename(mesh.path)
    CACHE_DIR = os.environ.get(
        "PYFESOM_CACHE",
        os.path.join(os.getcwd(), "MESH_cache"),
    )
    CACHE_DIR = os.path.join(CACHE_DIR, MESH_BASE)
    os.makedirs(CACHE_DIR, exist_ok=True)

    distances_file = (
        f"distances_{mesh.n2d}_{left}_{right}_{down}_{up}_"
        f"{lonNumber}_{latNumber}_{kk}"
    )
    inds_file = (
        f"inds_{mesh.n2d}_{left}_{right}_{down}_{up}_"
        f"{lonNumber}_{latNumber}_{kk}"
    )

    distances_paths = [
        distances_path,
        os.path.join(mesh.path, distances_file),
        os.path.join(CACHE_DIR, distances_file),
    ]
    inds_paths = [
        inds_path,
        os.path.join(mesh.path, inds_file),
        os.path.join(CACHE_DIR, inds_file),
    ]

    distances = inds = None

    # Try loading cached distances
    for path in filter(None, distances_paths):
        if os.path.isfile(path):
            try:
                distances = joblib.load(path)
                break
            except PermissionError:
                pass

    # Try loading cached indices
    for path in filter(None, inds_paths):
        if os.path.isfile(path):
            try:
                inds = joblib.load(path)
                break
            except PermissionError:
                pass

    # Compute if not available
    if distances is None or inds is None:
        distances, inds = create_indexes_and_distances(
            mesh, lons, lats, k=kk
        )

        if dumpfile:
            for path in distances_paths:
                if path is None:
                    continue
                try:
                    joblib.dump(distances, path)
                    break
                except PermissionError:
                    pass

            for path in inds_paths:
                if path is None:
                    continue
                try:
                    joblib.dump(inds, path)
                    break
                except PermissionError:
                    pass

    # Nearest-neighbor interpolation
    data_interpolated = data[inds]
    data_interpolated[distances >= radius_of_influence] = np.nan
    data_interpolated = data_interpolated.reshape(lons.shape)

    return np.ma.masked_invalid(data_interpolated)

# Loop over variables and interpolate
for variable in variables:
    print(f"Interpolating variable: {variable}")

    output_file = f"{dest_path}{variable}.interp.{min_lon}E.{max_lon}E.{lon_increment}E.{min_lat}N.{max_lat}N.{lat_increment}N.fesom.{year}.nc"
    if isfile(output_file):
        print(f"Output file already exists, skipping: {output_file}")
        continue

    # Load the variable data
    data_file = f"{base_path_model}{variable}.fesom.{year}.nc"
    ds = xr.open_dataset(data_file, decode_times=time_coder)

    # allocate array for 2D and 3D variables
    ndims = len(list(ds.dims))
    if ndims == 2:
        interp = (
            np.zeros((len(lat), len(lon), ds.dims["time"]), dtype=np.float32)
            * np.nan
        )

    elif ndims == 3:
        interp = (
            np.zeros((len(lat), len(lon), nlevs, ds.dims["time"]), dtype=np.float32)
            * np.nan
        )

    if ndims == 2:
        for i in range(len(ds.time)):
            print(f"Processing time step {i+1}/{len(ds.time)} for 2D variable {variable}")
        
            # Extract the variable data for the current time step and vertical level
            data = ds[variable].isel(time=i).values

            # Interpolate to regular grid
            interp_2D = fesom2regular_nn(
                data, mesh, lon_grid, lat_grid, radius_of_influence=radius_of_influence, dumpfile=False
            )

            # Add vertical level and time dimension
            interp[:, :, i] = interp_2D[:, :]
        
            # Create a new xarray Dataset for the interpolated data)
        interp_ds = xr.Dataset(
        {variable: (["lat", "lon", "time"], interp)},
        coords={"lon": lon, "lat": lat, "time": ds.time},
        )


    elif ndims == 3:
        for i in range(len(ds.time)):
            print(f"Processing time step {i+1}/{len(ds.time)} for 3D variable {variable}")
            for j in range(nlevs):
                print(
                    f"Processing vertical level {j+1}/{nlevs} for variable {variable} at time step {i+1}"
                )

                # Extract the variable data for the current time step and vertical level
                data = ds[variable].isel(time=i, nz1=j).values

                # Interpolate to regular grid
                interp_2D = fesom2regular_nn(
                    data, mesh, lon_grid, lat_grid, radius_of_influence=radius_of_influence, dumpfile=False
                )

                # Add vertical level and time dimension
                interp[:, :, j, i] = interp_2D[:, :]

        # Create a new xarray Dataset for the interpolated data)
        interp_ds = xr.Dataset(
            {variable: (["lat", "lon", "nz1", "time"], interp)},
            coords={"lon": lon, "lat": lat, "nz1": mesh_diag.nz1, "time": ds.time},
        )

    # Save the interpolated data
    interp_ds.to_netcdf(output_file)

    print(f"Saved interpolated data to: {output_file}")

print("Interpolation completed for all variables.")
