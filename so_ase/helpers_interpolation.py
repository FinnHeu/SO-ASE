# so_ase.helpers_interpolation.py

import numpy as np
import xarray as xr
from scipy.interpolate import griddata
from so_ase.helpers_mesh import read_nodes, build_land_sea_mask

def fesom_to_gridded(meshpath, data, varname, lon_grid, lat_grid, method="nearest", mask_land=True):
    """
    Interpolate FESOM2 data to a regular grid.
    
    Parameters:
    -----------
    meshpath : str
        Path to the FESOM2 mesh folder.
    data : xarray.DataArray
        FESOM2 data to interpolate. Only single time step supported.
    varname : str
        Name of the variable.
    lon_grid : numpy.ndarray
        1D array of longitudes for the regular grid.
    lat_grid : numpy.ndarray
        1D array of latitudes for the regular grid.
    method : str, optional
        Interpolation method. Default is "nearest". See scipy.interpolate.griddata for options.
    mask_land : bool, optional
        If True, mask land cells with NaN. Default is True.
    
    Returns:
    --------
    xarray.Dataset
        Interpolated data on the regular grid.
    """

    # extract lon and lat from fesom mesh
    lon, lat, _, _, = read_nodes(meshpath)
    lon, lat = lon, lat
        
    # define regular grid
    lon2d, lat2d = np.meshgrid(lon_grid, lat_grid)

    # interpolate
    data_gridded = griddata(
        (lon, lat),
        data.values,
        (lon2d, lat2d),
        method=method
    )

    # mask land cells
    if mask_land:
        ds_mask = build_land_sea_mask(meshpath, nlon=len(lon_grid), nlat=len(lat_grid), has_cavity=False, cavity_is_land=False)
        data_gridded = np.where(ds_mask.mask.values == 0, np.nan, data_gridded)
    
    
    ds = xr.Dataset(
        data_vars={
            varname: (("lat", "lon"), data_gridded)
        },
        coords={
            "lon": lon_grid,
            "lat": lat_grid
        },
        attrs={
            "description": "Interpolated FESOM2 " + varname,
        }
    )

    return ds