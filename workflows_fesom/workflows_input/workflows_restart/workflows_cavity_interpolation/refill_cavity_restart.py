
# =============================================================================
# ============================== READ MODULES =================================
# =============================================================================
import warnings
warnings.filterwarnings("ignore")

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import cartopy.crs as ccrs
from scipy.spatial import cKDTree
import os
import so_ase as so


# =============================================================================
# ============================== USER SETTINGS ================================
# =============================================================================

# Mesh directories
path_mesh = "/work/ab0995/a270186/model_inputs/fesom2/mesh/DARS2cav/"

# Restart Source (already generated restarts on DARS2cav grid for the branchoff initialization run)
restart_year = 1599
path_restart_src_oce = f"/work/ba1550/a270186/simulations/awiesm3-v3.4.1/restarts/DARS2_to_DARS2cav/CAVini/fesom.{restart_year}.oce.restart/"
path_restart_src_ice = f"/work/ba1550/a270186/simulations/awiesm3-v3.4.1/restarts/DARS2_to_DARS2cav/CAVini/fesom.{restart_year}.ice.restart/"

# Cavity Restart (restarts from which to take the cavity values)
path_restart_cavity_fill_oce = "/work/ba1550/a270301/runtime/awiesm3-v3.4.1/branchoff_DARS2cav/restart/fesom/fesom.1611.oce.restart/"
path_restart_cavity_fill_ice = "/work/ba1550/a270301/runtime/awiesm3-v3.4.1/branchoff_DARS2cav/restart/fesom/fesom.1611.ice.restart/"

# Restart Destination (location where restart files are going to be stored)
path_restart_dst_oce = f"/work/ba1550/a270186/simulations/awiesm3-v3.4.1/restarts/DARS2_to_DARS2cav/CAV/fesom.{restart_year}.oce.restart/"
path_restart_dst_ice = f"/work/ba1550/a270186/simulations/awiesm3-v3.4.1/restarts/DARS2_to_DARS2cav/CAV/fesom.{restart_year}.ice.restart/"

# Plots
plot = True
path_dst_plots = "./plots/refill/"

# =============================================================================
# ================================= FUNCTIONS =================================
# =============================================================================

def fill_cavities_from_existing_restart(varname, path_restart_src, path_restart_dst, path_restart_cavity_fill, verbose=True):
    """ 
    Replace the cavity values in the restart files with the values from another restart files.

    Parameters
    ----------
    varname : str
        Name of the variable to fill.
    path_restart_src_oce : str
        Path to the source restart file.
    path_restart_cavity_fill_oce : str
        Path to the cavity fill source restart file.
    path_restart_cavity_fill_ice : str
        Path to the cavity fill source ice restart file.
    """

    if verbose:
        print('\n============ fill_cavities_from_existing_restart.py ==============')
        print(f'Filling cavity values from existing restart file:')
        print(f'Source: {path_restart_dst_oce}{varname}.nc')
        print(f'Cavity fill source: {path_restart_cavity_fill_oce}{varname}.nc')
        print(f'Target: {path_restart_dst_oce}{varname}.nc')
        print(f'Variable name: {varname}.nc')
    
    
    # Load the restart to be filled with new cavity values
    print(f'Loading {path_restart_src}{varname}.nc')
    ds_to_fill = xr.open_dataset(f"{path_restart_src}{varname}.nc")

    # Load the dataset to fill the cavities from 
    ds_fill = xr.open_dataset(f"{path_restart_cavity_fill}{varname}.nc").isel(time=-1)

    # Build a horizontal cavity mask by checking if the first layer is zero
    mesh_diag = xr.open_dataset(f"{path_mesh}fesom.mesh.diag.nc")
    if 'node' in ds_to_fill.dims:
        print("Using node dimension")
        cavity_mask = mesh_diag.zbar_n_surface != 0
        cavity_mask = cavity_mask.rename({"nod2": "node"})
    elif 'elem' in ds_to_fill.dims:
        print("Using elem dimension")
        cavity_mask = mesh_diag.zbar_e_surface != 0
        # Rename mesh_diag dimension to match restart file dimension
        if 'ne' in cavity_mask.dims:
            cavity_mask = cavity_mask.rename({"ne": "elem"})
        elif 'elem2' in cavity_mask.dims:
            cavity_mask = cavity_mask.rename({"elem2": "elem"})
    else:
        raise ValueError("No node or elem dimension found in the dataset")
    
    
    print("Cavity mask shape:", cavity_mask.shape)
    print("Dataset to fill shape:", ds_to_fill[varname].shape)
    print("Fill dataset shape:", ds_fill[varname].shape)

    # Replace the cavity values with the values from the fill dataset
    ds_to_fill[varname] = ds_to_fill[varname].where(~cavity_mask, ds_fill[varname])

    # Save the filled dataset
    ds_to_fill.to_netcdf(f"{path_restart_dst}{varname}.nc")

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
        Subsampling factor for plotting, by default 20.
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
        so.create_map(axis, extent=box, land=True, coastline=True, circular=False)
    
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

if __name__ == "__main__":
    # Create destination directories if they do not exist
    if not os.path.isdir(path_restart_dst_oce):
        os.makedirs(path_restart_dst_oce)
    if not os.path.isdir(path_restart_dst_ice):
        os.makedirs(path_restart_dst_ice)
    if plot and not os.path.isdir(path_dst_plots):
        os.makedirs(path_dst_plots)

    # Read mesh coordinates for plotting
    node_lon, node_lat, node_idx, _ = so.read_nodes(path_mesh)
    elem, elem_lon, elem_lat = so.read_elements(path_mesh, return_coordinates=True)

    # Variables to process/copy
    vars_oce = ['salt', 'temp', 'temp_AB', 'salt_AB', 'temp_M1', 'salt_M1']
    vars_oce_copy = ['ssh', 'ssh_rhs_old', 'hbar', 'hnode', 'u', 'v', 'vrhs_AB', 'urhs_AB', 'urhs_AB3', 'vrhs_AB3', 'w', 'w_impl', 'w_expl'] 
    vars_ice_copy = ['area', 'hsnow', 'hice', 'uice', 'vice', 'ice_temp', 'ice_albedo']
    

    # Process required ocean variables which need cavity refill
    for var in vars_oce:
        fill_cavities_from_existing_restart(var, path_restart_src_oce, path_restart_dst_oce, path_restart_cavity_fill_oce)
        if plot:
            plot_refill_comparison(var, path_restart_src_oce, path_restart_dst_oce, path_restart_cavity_fill_oce,
                                   node_lon, node_lat, path_dst_plots)

    # Copy variables that don't need cavity filling
    for var in vars_oce_copy:
        print(f'Copying {var}...')
        os.system(f'cp {path_restart_src_oce}/{var}.nc {path_restart_dst_oce}/{var}.nc')
    
    for var in vars_ice_copy:
        print(f'Copying {var}...')
        os.system(f'cp {path_restart_src_ice}/{var}.nc {path_restart_dst_ice}/{var}.nc')

    print('=============================================================================')
    print('============================ REFILL COMPLETE ================================')
    print('=============================================================================')  