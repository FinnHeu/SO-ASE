# so_ase.helpers_interpolation.py

import numpy as np
import xarray as xr
from scipy.interpolate import griddata

def fesom2gridded(meshpath, data, varname, lon_grid, lat_grid, method="nearest"):
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
    
    Returns:
    --------
    xarray.Dataset
        Interpolated data on the regular grid.
    """

    # extract lon and lat from fesom mesh
    lon, lat, _, _, = so.read_nodes(meshpath)
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