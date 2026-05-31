# so_ase/plotting/plotting_restarts.py

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import cartopy.crs as ccrs

from .plotting_maps import create_map


def plot_mapper(mapper, lon_src, lat_src, lon_tgt, lat_tgt, horiz, path_dst_plots, n=10000):
    """
    Plot the mapper.
    
    Parameters
    ----------
    mapper : np.ndarray
        Mapper array.
    lon_src : np.ndarray
        Source longitudes.
    lat_src : np.ndarray
        Source latitudes.
    lon_tgt : np.ndarray
        Target longitudes.
    lat_tgt : np.ndarray
        Target latitudes.
    horiz : str
        Horizontal grid type.
    path_dst_plots : str
        Path to the plots directory.
    n : int, optional
        Number of points to plot, by default 10000
    """
    fig, ax = plt.subplots(1,1, figsize=(40,20), subplot_kw=dict(projection=ccrs.PlateCarree()))

    box = [-180, 180, -90, 90]
    
    for i in range(len(mapper))[::n]:
        lon_ori = lon_src[mapper[i]]
        lat_ori = lat_src[mapper[i]]

        lon_q = lon_tgt[i]
        lat_q = lat_tgt[i] 
        
        color = np.random.rand(3,)

        ax.plot([lon_ori, lon_q], [lat_ori, lat_q], color=color, marker='d', linestyle='-', transform=ccrs.PlateCarree())
        ax.plot(lon_q, lat_q, color=color, marker='o', markeredgecolor='k', markersize=5, linestyle='-', transform=ccrs.PlateCarree())


    create_map(ax, extent=box, land=False, coastline=True, circular=False)

    if horiz == 'elem':
        filename = 'mapper_elem_plot.png'
    elif horiz == 'node':
        filename = 'mapper_node_plot.png'
    
    plt.savefig(f"{path_dst_plots}{filename}", bbox_inches='tight', dpi=300)
    plt.close()


def plot_interpolated_extrapolated_field(path_src, path_int, path_restart_tgt, varname, lon_src, lat_src, lon_int, lat_int, path_dst_plots, n=20, level=15):
    """
    Plot the comparison of source and interpolated/extrapolated field.
    
    Parameters
    ----------
    path_src : str
        Path to the source directory.
    path_int : str
        Path to the interpolated directory.
    path_restart_tgt : str
        Path to the target restart directory.
    varname : str
        Name of the variable to plot.
    lon_src : np.ndarray
        Source longitudes.
    lat_src : np.ndarray
        Source latitudes.
    lon_int : np.ndarray
        Interpolated longitudes.
    lat_int : np.ndarray
        Interpolated latitudes.
    path_dst_plots : str
        Path to the plots directory.
    n : int, optional
        Number of points to plot, by default 20
    level : int, optional
        Level to plot, by default 15
    """
    print(f"Plotting comparison of source and interpolated/extrapolated field for variable: {varname}")
    print(' ')

    data_src = xr.open_dataset(f"{path_src}{varname}.nc").isel(time=-1)[varname].values
    data_int = xr.open_dataset(f"{path_int}{varname}.nc").isel(time=-1)[varname].values
    data_tgt = xr.open_dataset(f"{path_restart_tgt}{varname}.nc").isel(time=-1)[varname].values

    if len(data_src.shape) == 2:
        data_src = data_src[level,:]
        data_int = data_int[level,:]
        data_tgt = data_tgt[level,:]

    box = [-180, 180, -90, 90]
    
    fig, ax = plt.subplots(1,3, figsize=(40,15), subplot_kw=dict(projection=ccrs.PlateCarree()))
    
    for axis in ax:
        create_map(axis, extent=box, land=True, coastline=True, circular=False)

    min_src = np.min(data_src)
    max_src = np.max(data_src)
    min_int = np.min(data_int)
    max_int = np.max(data_int)

    if min_src < min_int:
        vmin = min_src
    else:
        vmin = min_int

    if max_src > max_int:
        vmax = max_src
    else:
        vmax = max_int

    lev = np.linspace(vmin, vmax, 20)

    try:
        cb = ax[0].tricontourf(lon_src[::n], lat_src[::n], data_src[::n], levels=lev, cmap='viridis', transform=ccrs.PlateCarree())
        ax[1].tricontourf(lon_int[::n], lat_int[::n], data_int[::n], levels=lev, cmap='viridis', transform=ccrs.PlateCarree())
        ax[2].tricontourf(lon_int[::n], lat_int[::n], data_tgt[::n], levels=lev, cmap='viridis', transform=ccrs.PlateCarree())

        plt.colorbar(cb, ax=ax, shrink=.7)

        ax[0].set_title('Source')
        ax[1].set_title('Interpolated/Extrapolated')
        ax[2].set_title('Target')

    except:
        print('Plotting failed!')

    filename = f"compare_{varname}.png"
    plt.savefig(f"{path_dst_plots}{filename}", bbox_inches='tight', dpi=300)
    plt.close()


def plot_refill_comparison(varname, path_restart_src, path_restart_dst, path_restart_cavity_fill, 
                           lon, lat, path_dst_plots, n=50):
    """
    Plot comparison of the dataset to be filled, the filling dataset, and the filled result.
    
    For 2D data: plots all three directly.
    For 3D data: uses the first active layer.
    Third subplot shows binary colors indicating whether filled data matches 
    the to-be-filled dataset or the filling dataset.
    
    Parameters
    ----------
    varname : str
        Name of the variable to plot.
    path_restart_src : str
        Path to the source restart directory (dataset to be filled).
    path_restart_dst : str
        Path to the destination restart directory (filled dataset).
    path_restart_cavity_fill : str
        Path to the cavity fill restart directory (filling dataset).
    lon : np.ndarray
        Longitudes of grid points.
    lat : np.ndarray
        Latitudes of grid points.
    path_dst_plots : str
        Path to save the plots.
    n : int, optional
        Subsampling factor for plotting, by default 50.
    """
    print(f"Plotting refill comparison for variable: {varname}")
    
    # Load datasets
    ds_to_fill = xr.open_dataset(f"{path_restart_src}{varname}.nc").isel(time=-1)
    ds_fill = xr.open_dataset(f"{path_restart_cavity_fill}{varname}.nc").isel(time=-1)
    ds_filled = xr.open_dataset(f"{path_restart_dst}{varname}.nc").isel(time=-1)
    
    # Extract data arrays
    data_to_fill = ds_to_fill[varname].values
    data_fill = ds_fill[varname].values
    data_filled = ds_filled[varname].values
    
    # Determine if data has depth dimension and extract appropriate layer
    has_depth = len(data_to_fill.shape) == 2
    
    if has_depth:
        # Get first active layer by finding first non-zero value along depth dimension
        n_points = data_to_fill.shape[1]
        first_layer_idx = np.argmax(data_to_fill != 0, axis=0)
        
        # Extract data at first active layer for each point
        data_to_fill_2d = data_to_fill[first_layer_idx, np.arange(n_points)]
        data_fill_2d = data_fill[first_layer_idx, np.arange(n_points)]
        data_filled_2d = data_filled[first_layer_idx, np.arange(n_points)]
    else:
        # 2D data - use directly
        data_to_fill_2d = data_to_fill
        data_fill_2d = data_fill
        data_filled_2d = data_filled
    
    # Create figure with 3 subplots
    box = [-180, 180, -90, 90]
    fig, ax = plt.subplots(1, 3, figsize=(40, 15), subplot_kw=dict(projection=ccrs.PlateCarree()))
    
    for axis in ax:
        create_map(axis, extent=box, land=True, coastline=True, circular=False)
    
    # Compute common color limits for first two plots
    vmin = min(np.nanmin(data_to_fill_2d), np.nanmin(data_fill_2d))
    vmax = max(np.nanmax(data_to_fill_2d), np.nanmax(data_fill_2d))
    
    try:
        # Plot 1: Dataset to be filled
        sc1 = ax[0].scatter(lon[::n], lat[::n], c=data_to_fill_2d[::n], 
                            s=1, cmap='viridis', vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())
        ax[0].set_title('To Be Filled (Source)', fontsize=14)
        
        # Plot 2: Filling dataset
        ax[1].scatter(lon[::n], lat[::n], c=data_fill_2d[::n], 
                      s=1, cmap='viridis', vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())
        ax[1].set_title('Filling Dataset (Cavity Source)', fontsize=14)
        
        plt.colorbar(sc1, ax=ax[:2], shrink=0.7, orientation='horizontal', pad=0.05)
        
        # Plot 3: Binary comparison - which source does the filled data match?
        # 1 = matches to_fill, 2 = matches fill_source
        tol = 1e-10
        matches_to_fill = np.abs(data_filled_2d - data_to_fill_2d) < tol
        matches_fill = np.abs(data_filled_2d - data_fill_2d) < tol
        
        # Create binary array: 0 = matches to_fill, 1 = matches fill_source
        binary_match = np.zeros_like(data_filled_2d)
        binary_match[matches_fill & ~matches_to_fill] = 1
        
        # Custom colormap for binary plot
        binary_cmap = ListedColormap(['#1f77b4', '#ff7f0e'])  # Blue for to_fill, Orange for fill
        
        scatter = ax[2].scatter(lon[::n], lat[::n], c=binary_match[::n], 
                                cmap=binary_cmap, s=1, transform=ccrs.PlateCarree(),
                                vmin=0, vmax=1)
        ax[2].set_title('Match Source (Blue=Original, Orange=Cavity Fill)', fontsize=14)
        
        # Add colorbar for binary plot
        cbar = plt.colorbar(scatter, ax=ax[2], shrink=0.7, orientation='horizontal', pad=0.05)
        cbar.set_ticks([0.25, 0.75])
        cbar.set_ticklabels(['Original', 'Cavity Fill'])
        
    except Exception as e:
        print(f'Plotting failed for {varname}: {e}')
    
    filename = f"refill_compare_{varname}.png"
    plt.savefig(f"{path_dst_plots}{filename}", bbox_inches='tight', dpi=300)
    plt.close()
    print(f"Saved plot: {path_dst_plots}{filename}")
