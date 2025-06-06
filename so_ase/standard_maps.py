# so_ase/maps.py

import matplotlib.pyplot as plt
import cartopy.crs as ccrs

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.path as mpath
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import numpy as np

def create_map(ax, extent='global', land=True, coastline=True, lon_inc=30, lat_inc=5, tick_labels=True, circular=True):
    """
    Creates a map with optional land and coastline features, 
    gridlines.

    Parameters:
        ax (matplotlib.axes._subplots.AxesSubplot): 
            The axes object to plot on.
        extent (list or str, optional): 
            The geographic extent of the map in the form [lon_min, lon_max, lat_min, lat_max] or str 'global', 'southern_ocean' or 'arctic_ocean'. 
            Defaults to 'global'.
        land (bool, optional): 
            Whether to add land features to the map. Defaults to True.
        coastline (bool, optional): 
            Whether to add coastlines to the map. Defaults to True.
        lon_inc (int, optional): 
            Longitude grid interval. Defaults to 30.
        lat_inc (int, optional): 
            Latitude grid interval. Defaults to 5.
        tick_labels (bool, optional): 
            Whether to display longitude and latitude labels on the gridlines. Defaults to True.
        circular (bool, optional):
            Whether to apply a circular shape to the map. Defaults to True.

    Returns:
        matplotlib.axes._subplots.AxesSubplot: 
            The modified axes object with the Southern Ocean map.
    """

    # handle extent and circular input
    if isinstance(extent, str):
        if extent == 'global':
            extent = [-180, 180, -90, 90]
            circular = False  # Global maps are not typically circular
        elif extent == 'southern_ocean':
            extent = [-180, 180, -90, -40]
        elif extent == 'arctic_ocean':
            extent = [-180, 180, 65, 90]
        else:
            raise ValueError("Invalid extent string. Use 'global', 'southern_ocean', or 'arctic_ocean'.")
    
    # validate extent input
    elif not isinstance(extent, (list, tuple)) or len(extent) != 4:
        raise ValueError("Extent must be a list or tuple of four values: [lon_min, lon_max, lat_min, lat_max].")
    
    if not all(isinstance(x, (int, float)) for x in extent):
        raise ValueError("Extent values must be numeric (int or float).")
    
    if extent[0] >= extent[1] or extent[2] >= extent[3]:
        raise ValueError("Invalid extent: lon_min must be less than lon_max and lat_min must be less than lat_max.")
    
    # set map extent
    ax.set_extent(extent, crs=ccrs.PlateCarree())

    # add land/coastlines
    if land:
        ax.add_feature(cfeature.LAND, color='lightgrey')
    if coastline:
        ax.add_feature(cfeature.COASTLINE, color='black', linewidth=0.5)

    # add grid and grid labels
    gl = ax.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=tick_labels,
        linewidth=1,
        color="gray",
        alpha=0.5,
        linestyle="--",
        x_inline=False,
        y_inline=True,
    )

    gl.xlocator = mticker.FixedLocator(range(-180, 180, lon_inc))
    gl.ylocator = mticker.FixedLocator(range(-90, 90, lat_inc))
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {"size": 8, "rotation": 0}
    gl.ylabel_style = {"size": 8, "rotation": 30}

    # make plot circular
    if circular:
        ax = circular_shape(ax)
    
    return ax

def circular_shape(ax):
    """
    Clips the plotting area to a circular shape, typically used to create 
    polar plots or maps with circular boundaries.

    Parameters:
        ax (matplotlib.axes._subplots.AxesSubplot): 
            The axes object to apply the circular boundary to.

    Returns:
        matplotlib.axes._subplots.AxesSubplot: 
            The modified axes object with a circular clipping boundary.
    """
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)

    ax.set_boundary(circle, transform=ax.transAxes)
    return ax
