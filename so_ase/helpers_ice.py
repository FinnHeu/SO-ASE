# so_ase.helpers_ice.py

import xarray as xr
import numpy as np
from pyproj import Proj, Transformer

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