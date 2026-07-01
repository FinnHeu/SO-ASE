# so_ase/fesom2/helpers_interpolation.py

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
    Interpolate 1D FESOM data to a regular lat/lon grid using nearest-neighbor interpolation.

    This function maps unstructured FESOM mesh data onto a regular grid by finding the
    nearest mesh node for each target grid point. It supports caching of KDTree query
    results (distances and indices) to speed up repeated interpolations on the same grid.

    Parameters
    ----------
    data : array-like
        1D array of values on FESOM mesh nodes to be interpolated.
    mesh : fesom_mesh object
        pyfesom mesh representation containing mesh geometry (x2, y2 coordinates).
    lons : ndarray
        2D array of target grid longitudes (e.g., from np.meshgrid).
    lats : ndarray
        2D array of target grid latitudes (e.g., from np.meshgrid).
    distances_path : str, optional
        Path to a cached distances file. If provided and exists, distances are loaded
        from this file instead of being computed.
    inds_path : str, optional
        Path to a cached indices file. If provided and exists, indices are loaded
        from this file instead of being computed.
    radius_of_influence : float, optional
        Maximum distance (in meters) for valid interpolation. Grid points farther
        than this from any mesh node are set to NaN. Default is 100000 (100 km).
    dumpfile : bool, optional
        If True, cache computed distances and indices to disk for future use.
        Default is True.

    Returns
    -------
    np.ma.MaskedArray
        2D masked array with shape matching `lons`/`lats`, containing interpolated
        values. Points outside the radius of influence are masked.

    Notes
    -----
    Cache files are searched in the following order:
    1. User-specified path (distances_path/inds_path)
    2. Mesh directory
    3. PYFESOM_CACHE environment variable or ./MESH_cache
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
    """
    Batch interpolate FESOM output files to a regular lat/lon grid.

    This function processes multiple years and variables of FESOM output, interpolating
    each from the unstructured mesh to a user-defined regular grid. It handles both
    2D (time, nodes) and 3D (time, depth, nodes) variables automatically.

    Parameters
    ----------
    datapath : str
        Path to directory containing FESOM output files. Files are expected to follow
        the naming convention: `{variable}.fesom.{year}.nc`.
    meshpath : str
        Path to FESOM mesh directory containing mesh files and `fesom.mesh.diag.nc`.
    reg_lat : tuple of float, optional
        Latitude range (min, max) for the target grid. Default is (-90, 90).
    reg_lon : tuple of float, optional
        Longitude range (min, max) for the target grid. Default is (-180, 180).
    lat_increment : float, optional
        Latitude spacing of the target grid in degrees. Default is 1.
    lon_increment : float, optional
        Longitude spacing of the target grid in degrees. Default is 1.
    variables : list of str, optional
        List of variable names to interpolate. Default is ['temp'].
    years : tuple of int, optional
        Year range (start, end) inclusive to process. Default is (1850, 1860).
    radius_of_influence : float, optional
        Maximum distance (in meters) for valid interpolation. Default is 100000.
    dest_path : str, optional
        Output directory for interpolated NetCDF files. Default is './'.
    dumpfile : bool, optional
        If True, cache KDTree distances/indices for faster subsequent runs.
        Default is False.
    log : bool, optional
        If True, print progress messages. Default is False.

    Returns
    -------
    None
        Interpolated data is written to NetCDF files at `dest_path` with naming
        convention: `{variable}.interp.fesom.{year}.nc`.

    Notes
    -----
    - Existing output files are skipped (not overwritten).
    - 2D variables produce output with dimensions (lat, lon, time).
    - 3D variables produce output with dimensions (lat, lon, nz1, time).
    - Vertical levels are taken from the mesh diagnostic file.
    """

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


def fesom2regular_binned(
    data,
    meshpath,
    which='node',
    reg_lon=(-180, 180),
    reg_lat=(-90, 90),
    lon_increment=1.0,
    lat_increment=1.0,
):
    """
    Regrid FESOM data to a regular lon/lat grid using area-weighted binning.

    For each grid cell, finds all FESOM nodes/elements whose coordinates fall
    within the cell, computes the area-weighted mean of their values, and
    returns the result on the regular grid.

    Parameters
    ----------
    data : array-like
        1D array of values on FESOM mesh nodes or elements.
    meshpath : str
        Path to the directory containing `fesom.mesh.diag.nc`.
    which : str, optional
        'node' for node-based data (uses nod_area), 'element' for element-based
        data (uses elem_area). Default is 'node'.
    reg_lon : tuple of float, optional
        Longitude range (min, max) for the target grid. Default is (-180, 180).
    reg_lat : tuple of float, optional
        Latitude range (min, max) for the target grid. Default is (-90, 90).
    lon_increment : float, optional
        Longitude spacing of the target grid in degrees. Default is 1.0.
    lat_increment : float, optional
        Latitude spacing of the target grid in degrees. Default is 1.0.

    Returns
    -------
    xarray.Dataset
        Dataset with the binned variable 'data' on dimensions (lat, lon).
        Grid cells with no data are set to NaN.

    Example
    -------
    >>> binned = fesom2regular_binned(
    ...     ds['temp'].isel(time=0, nz1=0).values,
    ...     '/path/to/mesh/',
    ...     which='node',
    ...     reg_lon=(-180, 180),
    ...     reg_lat=(-90, -60),
    ...     lon_increment=0.5,
    ...     lat_increment=0.5
    ... )
    """
    data = np.asarray(data)

    # Load mesh diagnostics
    mesh_diag = xr.open_dataset(f"{meshpath}fesom.mesh.diag.nc")

    if which == 'node':
        lon = mesh_diag.lon.values
        lat = mesh_diag.lat.values
        area = mesh_diag.nod_area.max(dim='nz').values
    elif which == 'element':
        lon = mesh_diag.elem_lon.values
        lat = mesh_diag.elem_lat.values
        area = mesh_diag.elem_area.values
    else:
        mesh_diag.close()
        raise ValueError(f"which must be 'node' or 'element', got {which}")

    mesh_diag.close()

    # Define grid edges
    lon_edges = np.arange(reg_lon[0], reg_lon[1] + lon_increment, lon_increment)
    lat_edges = np.arange(reg_lat[0], reg_lat[1] + lat_increment, lat_increment)
    
    # Grid cell centers
    lon_centers = (lon_edges[:-1] + lon_edges[1:]) / 2
    lat_centers = (lat_edges[:-1] + lat_edges[1:]) / 2
    
    nlat = len(lat_centers)
    nlon = len(lon_centers)

    # Compute bin indices for each point
    lon_idx = np.digitize(lon, lon_edges) - 1
    lat_idx = np.digitize(lat, lat_edges) - 1

    # Mask points outside the grid
    valid = (
        (lon_idx >= 0) & (lon_idx < nlon) &
        (lat_idx >= 0) & (lat_idx < nlat) &
        np.isfinite(data)
    )

    # Compute linear bin index for valid points
    bin_idx = lat_idx[valid] * nlon + lon_idx[valid]
    data_valid = data[valid]
    area_valid = area[valid]

    # Weighted sum and total area per bin
    weighted_sum = np.bincount(bin_idx, weights=data_valid * area_valid, minlength=nlat * nlon)
    total_area = np.bincount(bin_idx, weights=area_valid, minlength=nlat * nlon)

    # Compute area-weighted mean
    with np.errstate(divide='ignore', invalid='ignore'):
        binned = weighted_sum / total_area
    binned[total_area == 0] = np.nan

    # Reshape to 2D grid
    binned = binned.reshape((nlat, nlon))

    # Create xarray Dataset
    ds_out = xr.Dataset(
        data_vars={
            'data': (['lat', 'lon'], binned.astype(np.float32)),
            'cell_area': (['lat', 'lon'], total_area.reshape((nlat, nlon)).astype(np.float32)),
        },
        coords={
            'lon': lon_centers,
            'lat': lat_centers,
        },
        attrs={
            'description': 'Area-weighted binned FESOM data on regular grid',
            'lon_increment': lon_increment,
            'lat_increment': lat_increment,
        }
    )

    return ds_out