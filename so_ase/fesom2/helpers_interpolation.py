# so_ase.helpers_interpolation.py

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
from ..miscellaneous.helpers_misc import lon_lat_to_cartesian

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

def fesom2regular_1d(
    data,
    mesh,
    lons,
    lats,
    distances_path=None,
    inds_path=None,
    radius_of_influence=100000,
    dumpfile=True,
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

def fesom2regular_nd(
    datapath, 
    meshpath, 
    reg_lat=(-90, 90), 
    reg_lon=(-180, 180), 
    lat_increment=1,
    lon_increment=1,
    variables=['temp'], 
    years=(1850, 1860), 
    radius_of_influence=100000, 
    dest_path='./',
    dumpfile=False,
    log=False
):

    mesh = pf.load_mesh(meshpath, usepickle=False)
    mesh_diag = xr.open_dataset(f"{meshpath}fesom.mesh.diag.nc")

    reg_lats = np.arange(reg_lat[0], reg_lat[1], lat_increment)
    reg_lons = np.arange(reg_lon[0], reg_lon[1], lon_increment)
    lon_grid, lat_grid = np.meshgrid(reg_lons, reg_lats)
    if log:
        print(f"Generated grid: {len(reg_lats)} latitudes x {len(reg_lons)} longitudes")
        
    # Get number of vertical levels
    nlevs = len(mesh_diag.nz1)

    # Define time encoding
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)

    years = np.arange(years[0], years[1]+1)

    if log:
        print(f"Interpolating for years - {years}")

    for year in years:
        # Loop over variables and interpolate
        for variable in variables:
            if log:
                print(f"Interpolating variable: {variable}")
        
            output_file = f"{dest_path}{variable}.interp.fesom.{year}.nc"
            if isfile(output_file):
                if log:
                    print(f"Output file already exists, skipping: {output_file}")
                continue
        
            # Load the variable data
            data_file = f"{datapath}{variable}.fesom.{year}.nc"
            ds = xr.open_dataset(data_file, decode_times=time_coder)
        
            # allocate array for 2D and 3D variables
            ndims = len(list(ds.dims))
            if ndims == 2:
                interp = (
                    np.zeros((len(reg_lats), len(reg_lons), ds.sizes["time"]), dtype=np.float32)
                    * np.nan
                )
        
            elif ndims == 3:
                interp = (
                    np.zeros((len(reg_lats), len(reg_lons), nlevs, ds.sizes["time"]), dtype=np.float32)
                    * np.nan
                )
        
            if ndims == 2:
                for i in range(len(ds.time)):
                    if log:
                        print(f"Processing time step {i+1}/{len(ds.time)} for 2D variable {variable}")
                
                    # Extract the variable data for the current time step and vertical level
                    data = ds[variable].isel(time=i).values
        
                    # Interpolate to regular grid
                    interp_2D = fesom2regular_1d(
                        data, mesh, lon_grid, lat_grid, radius_of_influence=radius_of_influence, dumpfile=dumpfile
                    )
        
                    # Add vertical level and time dimension
                    interp[:, :, i] = interp_2D[:, :]
                
                    # Create a new xarray Dataset for the interpolated data)
                interp_ds = xr.Dataset(
                {variable: (["lat", "lon", "time"], interp)},
                coords={"lon": reg_lons, "lat": reg_lats, "time": ds.time},
                )
        
        
            elif ndims == 3:
                for i in range(len(ds.time)):
                    if log:
                        print(f"Processing time step {i+1}/{len(ds.time)} for 3D variable {variable}")
                    for j in range(nlevs):
                        # if log:
                        #     print(f"Processing vertical level {j+1}/{nlevs} for variable {variable} at time step {i+1}")
        
                        # Extract the variable data for the current time step and vertical level
                        data = ds[variable].isel(time=i, nz1=j).values
        
                        # Interpolate to regular grid
                        interp_2D = fesom2regular_1d(
                            data, mesh, lon_grid, lat_grid, radius_of_influence=radius_of_influence, dumpfile=dumpfile
                        )
        
                        # Add vertical level and time dimension
                        interp[:, :, j, i] = interp_2D[:, :]
        
                # Create a new xarray Dataset for the interpolated data)
                interp_ds = xr.Dataset(
                    {variable: (["lat", "lon", "nz1", "time"], interp)},
                    coords={"lon": reg_lons, "lat": reg_lats, "nz1": mesh_diag.nz1, "time": ds.time},
                )
        
            # Save the interpolated data
            interp_ds.to_netcdf(output_file)
            if log:
                print(f"Saved to {output_file}")

            