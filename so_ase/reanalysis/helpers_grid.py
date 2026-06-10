# so_ase/reanalysis/helpers_grid.py

import xarray as xr
import numpy as np
import math

from pyproj import Proj, Transformer

def gridcell_area_hadley(ds, R=6371.0):
    """
    Compute and add grid cell area (km²) to an xarray.Dataset on a regular 1° lat-lon grid.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with 'lat' and 'lon' coordinates.
    R : float, optional
        Radius of the Earth in kilometers (default: 6371 km).

    Returns
    -------
    ds_out : xarray.Dataset
        New dataset with an additional variable 'cell_area'.
    """
    # Broadcast 2D lat-lon grids
    lat2d, lon2d = xr.broadcast(ds.latitude, ds.longitude)

    # Convert to radians
    dlat = np.deg2rad(1.0)
    dlon = np.deg2rad(1.0)
    lat_rad = np.deg2rad(lat2d)

    # Compute latitude bounds
    lat1 = lat_rad - dlat / 2
    lat2 = lat_rad + dlat / 2

    # Area formula for a spherical Earth
    area = (R ** 2) * dlon * (np.sin(lat2) - np.sin(lat1))
    area = np.abs(area) * 1e6

    # Create DataArray for area
    area_da = xr.DataArray(
        area,
        coords=lat2d.coords,
        dims=lat2d.dims,
        name="cell_area",
        attrs={"units": "m^2", "description": "Grid cell area assuming spherical Earth"}
    )

    # Add to dataset
    ds_out = ds.copy()
    ds_out["cell_area"] = area_da

    return ds_out

def reproject_to_latlon(ds, input_proj="+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84"):
    """
    Transforms x/y coordinates of a polar-stereographic dataset (NSIDC datasets) into latitude/longitude and adds them as 2D coordinate variables.

    Parameters:
        ds (xr.Dataset): Dataset in projected coordinates (e.g., EPSG:3412).
        proj_str (str): PROJ string for the dataset's current projection.

    Returns:
        xr.Dataset: Dataset with 2D 'lat' and 'lon' variables added.
    """
    x = ds['x'].values
    y = ds['y'].values
    x2d, y2d = np.meshgrid(x, y)

    # Define projections explicitly
    proj_stereo = Proj(input_proj)
    proj_geo = Proj(proj='latlong', datum='WGS84')
    
    transformer = Transformer.from_proj(proj_stereo, proj_geo, always_xy=True)

    lon2d, lat2d = transformer.transform(x2d, y2d)

    ds['lat'] = (('y','x'), lat2d)
    ds['lon'] = (('y','x'), lon2d)
    ds['lat'].attrs['units'] = 'degrees_north'
    ds['lon'].attrs['units'] = 'degrees_east'
    
    return ds
    
