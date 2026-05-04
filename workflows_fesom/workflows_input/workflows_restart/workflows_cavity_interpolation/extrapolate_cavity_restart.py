# =============================================================================
# DARS2/CORE2 to DARS2cav/CORE2ice: Generate Restarts
# =============================================================================
#
# Here, FESOM2 restarts are taken from a DARS2/CORE2 simulation and extrapolated
# onto an existing, however empty, DARS2cav/CORE2ice restart.
#
# The aim is to avoid a full cold start with the DARS2cav mesh and instead
# branch off a DARS2 simulation at the end of the restart.
#
# Extrapolation strategy:
#
# Cavity:
#   - Temperature and Salinity: nearest neighbor
#   - Velocity: set to 0 (cold start) in all cavity elements
#
# Open Ocean:
#   - Temperature, Salinity, Velocity: nearest neighbor
#
# List of all restart variables
# ----------------------------------------------------------------------------- 
#
# Ocean:
#
# 2D node: (time, node)
#   - ssh.nc          (time, node) [x]
#   - ssh_rhs_old.nc  (time, node)
#   - hbar.nc         (time, node)
#
# 3D node nz1: (time, nz_1, node)
#   - hnode.nc    (time, nz_1, node)
#   - salt.nc     (time, nz_1, node)
#   - temp.nc     (time, nz_1, node)
#   - temp_AB.nc  (time, nz_1, node)
#   - salt_AB.nc  (time, nz_1, node)
#   - temp_M1.nc  (time, nz_1, node)
#   - salt_M1.nc  (time, nz_1, node)
#
# 3D node nz: (time, nz, node)
#   - w_impl.nc  (time, nz, node)
#   - w_expl.nc  (time, nz, node)
#   - w.nc       (time, nz_1, node)
#
# 3D element: (time, nz_1, elem)
#   - u.nc        (time, nz_1, elem)
#   - v.nc        (time, nz_1, elem)
#   - vrhs_AB.nc  (time, nz_1, elem)
#   - urhs_AB.nc  (time, nz_1, elem)
#
# Ice:
#   - area.nc
#   - hice.nc
#   - hsnow.nc
#   - uice.nc
#   - vice.nc
#
# =============================================================================

# =============================================================================
# ============================== READ MODULES =================================
# =============================================================================
import warnings
warnings.filterwarnings("ignore")

import xarray as xr
import numpy as np
import so_ase as so
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import pyfesom2 as pf
import cmocean.cm as cmo
from scipy.spatial import cKDTree
import os
import sys


# =============================================================================
# ============================== SET PATHS =================================
# =============================================================================

# Mesh directories
path_mesh_src = "/work/ab0995/a270186/model_inputs/fesom2/mesh/DARS2/"
path_mesh_tgt = "/work/ab0995/a270186/model_inputs/fesom2/mesh/DARS2cav/"

# Restart files on DARS2 mesh which are supposed to be used on DARS2cav mesh
path_restart_src_oce = f"/work/bb1469/a270089/runtime/awiesm3-v3.4.1/AWI-ESM3-VEG-HR-CMIP7-Spinup_cont1/restart/fesom/fesom.1475.oce.restart/"
path_restart_src_ice = f"/work/bb1469/a270089/runtime/awiesm3-v3.4.1/AWI-ESM3-VEG-HR-CMIP7-Spinup_cont1/restart/fesom/fesom.1475.ice.restart/"

# Restart files on DARS2cav mesh (template files of which the file structure is taken)
path_restart_tgt_oce = f"/work/ab0995/a270186/model_inputs/awicm3/pool/restarts/templates/DARS2cav/v2.7.1/fesom.1850.oce.restart/"
path_restart_tgt_ice = f"/work/ab0995/a270186/model_inputs/awicm3/pool/restarts/templates/DARS2cav/v2.7.1/fesom.1850.ice.restart/"

# Restart Destination (destination of the generated restart files)
restart_year = 1475
path_dst_restarts_oce = f"/work/ba1550/a270186/simulations/awicm3-develop/restarts_DARS2_to_DARS2_mod_blacksea/fesom.{restart_year}.oce.restart/"
path_dst_restarts_ice = f"/work/ba1550/a270186/simulations/awicm3-develop/restarts_DARS2_to_DARS2_mod_blacksea/fesom.{restart_year}.ice.restart/"

# Plots Destination
plot = True
path_dst_plots = "./plots/"

# Coupled Model (AWI-CM3 with FESOM2.7) also requires ice_temp.nc and ice_albedo.nc
is_coupled = True

# Path to restart files on target mesh ( ---> fesom v2.7 <---  ) for masking cavities
path_restart_tgt_oce_v27 = path_restart_tgt_oce

# Fill cavities from existing restart files
fill_cavities = False
path_restart_cavity_fill_oce = ""

# =============================================================================
# ============================ SET LOG FILES ==================================
# =============================================================================

# Create destination directories if they do not exist
if not os.path.isdir(path_dst_restarts_oce):
    os.makedirs(path_dst_restarts_oce)
if not os.path.isdir(path_dst_restarts_ice):
    os.makedirs(path_dst_restarts_ice)

# Determine the parent folder of path_dst_restarts_oce
log_file = os.path.join(os.path.dirname(path_dst_restarts_oce.rstrip('/')), 'extrapolate_cavity_restart.log')

# Tee class to write to both terminal and log file
class Tee:
    def __init__(self, *files):
        self.files = files
    def write(self, obj):
        for f in self.files:
            f.write(obj)
            f.flush()
    def flush(self):
        for f in self.files:
            f.flush()

# Redirect stdout and stderr to both terminal and log file
log_file_handle = open(log_file, 'w')
sys.stdout = Tee(sys.__stdout__, log_file_handle)
sys.stderr = Tee(sys.__stderr__, log_file_handle)

# =============================================================================
# ================================= FUNCTIONS =================================
# =============================================================================

def build_spherical_nn_mapper(lon_source, lat_source,
                              lon_target, lat_target):
    """
    Build nearest-neighbor mapping between two sets of points
    on a sphere (lon/lat in degrees). The KDTree is built on the source grid.

    Parameters
    ----------
    lon_source : array-like
        Longitudes of source grid points (degrees).
    lat_source : array-like
        Latitudes of source grid points (degrees).
    lon_target : array-like
        Longitudes of target grid points (degrees).
    lat_target : array-like
        Latitudes of target grid points (degrees).

    Returns
    -------
    mapping : ndarray
        For each target point, the index of the nearest source point.
    distances : ndarray
        The spherical Euclidean distance in 3D Cartesian space.
    """
    # --- Convert degrees to radians ---
    lon1 = np.radians(lon_source)
    lat1 = np.radians(lat_source)
    lon2 = np.radians(lon_target)
    lat2 = np.radians(lat_target)

    # --- Convert spherical coordinates to Cartesian ---
    def sph2cart(lon, lat):
        x = np.cos(lat) * np.cos(lon)
        y = np.cos(lat) * np.sin(lon)
        z = np.sin(lat)
        return np.column_stack((x, y, z))

    coords1 = sph2cart(lon1, lat1)
    coords2 = sph2cart(lon2, lat2)

    # --- KDTree on source grid ---
    tree = cKDTree(coords1)

    # --- Query nearest neighbor for every target point ---
    distances, indices = tree.query(coords2, k=1)

    return indices, distances

def interpolate_extrapolate_2D(varname, path_restart_src, path_restart_tgt, mapper_nodes, path_out, mask_file, mask_file_varname, t_step=-1, verbose=True):
    """
    Interpolate/extrapolate 2D variables from source to target grid.
    
    Parameters
    ----------
    varname : str
        Name of the variable to interpolate/extrapolate.
    path_restart_src : str
        Path to the source restart file.
    path_restart_tgt : str
        Path to the target restart file.
    mapper_nodes : callable
        Function to map source grid to target grid.
    path_out : str
        Path to the output file.
    mask_file : str
        Path to the mask file.
    mask_file_varname : str
        Name of the variable in the mask file.
    t_step : int, optional
        Time step to use, by default -1 (last time step).
    verbose : bool, optional
        Whether to print verbose output, by default True
    """
    
    if verbose:
        print('============ interpolate_extrapolate_node2D.py ==============')
        print(f'Inter/Extrapolating restart file:')
        print(f'Source: {path_restart_src}{varname}.nc')
        print(f'Target: {path_restart_tgt}{varname}.nc')
        print(f'Variable name: {varname}.nc')
        print('')

    # Open files
    ds_src = xr.open_dataset(f"{path_restart_src}{varname}.nc").isel(time=-1)
    ds_tgt = xr.open_dataset(f"{path_restart_tgt}{varname}.nc").isel(time=-1)
    ds_mask = xr.open_dataset(mask_file).isel(time=-1)

    # Extract array of source data and destination data
    data_src = ds_src[varname].values
    data_mask = ds_mask[mask_file_varname].values
    
    # Map source data (without cavity) to destination grid (with cavity)
    data_int = data_src[mapper_nodes]
    
    # Force bathymetry and cavity nodes to 0
    data_int[data_mask == 0] = 0 # cavity and bathymetry are all 0
    
    # Make a deep copy of destination restart
    ds_int = ds_tgt.copy()
    
    # Add a time dimension
    time = xr.open_dataset(f"{path_restart_src}{varname}.nc").time.values[t_step]
    ds_int = ds_int.assign_coords(time=time).expand_dims('time')
    
    # Extract iter
    iteration = xr.open_dataset(f"{path_restart_src}{varname}.nc").iter.values[t_step]
    
    # Write iter to deep copy
    ds_int['iter'] = (('time'), np.atleast_1d(iteration))
    
    # Write intermediate data to deep copy
    ds_int[varname] = (('time', 'node'), data_int[np.newaxis,:])
    
    # Add metadata
    for var in ds_src.data_vars:
        if var in ds_int:
            ds_int[var].attrs = ds_src[var].attrs.copy()
    
    ds_int.attrs['source restart (data)'] = f'{path_restart_src}{varname}.nc'
    ds_int.attrs['target restart (grid)'] = f'{path_restart_tgt}{varname}.nc'

    # Save
    ds_int.to_netcdf(f'{path_out}{varname}.nc')

    return

def interpolate_extrapolate_3D(varname, path_restart_src, path_restart_tgt, mapper_nodes, path_out, t_step=-1, verbose=True):
    """
    Interpolate/extrapolate 3D variables from source to target grid.
    
    Parameters
    ----------
    varname : str
        Name of the variable to interpolate/extrapolate.
    path_restart_src : str
        Path to the source restart file.
    path_restart_tgt : str
        Path to the target restart file.
    mapper_nodes : callable
        Function to map source grid to target grid.
    path_out : str
        Path to the output file.
    t_step : int, optional
        Time step to use, by default -1 (last time step).
    verbose : bool, optional
        Whether to print verbose output, by default True
    """
    if verbose:
        print('============ interpolate_extrapolate_node3D.py ==============')
        print(f'Inter/Extrapolating restart file:')
        print(f'Source: {path_restart_src}{varname}.nc')
        print(f'Target: {path_restart_tgt}{varname}.nc')
        print(f'Variable name: {varname}.nc')

    # Open files
    ds_src = xr.open_dataset(f"{path_restart_src}{varname}.nc").isel(time=t_step)
    ds_tgt = xr.open_dataset(f"{path_restart_tgt}{varname}.nc").isel(time=t_step)

    if varname == 'hnode':
        print('Special case: hnode.nc')
        print('Only copy the file from the target grid and modify time, iter, metadata: contains active layer thickness.')
        
        # Make a deep copy of destination restart
        ds_int = ds_tgt.copy()
        
        # Add a time dimension
        time = xr.open_dataset(f"{path_restart_src}{varname}.nc").time.values[t_step]
        ds_int = ds_int.assign_coords(time=time).expand_dims('time')
        
        # Extract iter
        iteration = xr.open_dataset(f"{path_restart_src}{varname}.nc").iter.values[t_step]
        
        # Write iter to deep copy
        ds_int['iter'] = (('time'), np.atleast_1d(iteration))
        
        ds_int.attrs['information'] = 'The hnode content of this restart file was copied from <source restart>' 
        ds_int.attrs['source restart (data)'] = f'{path_restart_src}{varname}.nc'
        ds_int.attrs['target restart (grid)'] = f'{path_restart_tgt}{varname}.nc'
    
        # Save
        ds_int.to_netcdf(f'{path_out}{varname}.nc')
        

    else:
        ###---> Extract array of source data and destination data
        data_src = ds_src[varname].values
        data_dst = ds_tgt[varname].values

        ###---> Make sure the are no NaN values in the source data
        data_src[~np.isfinite(data_src)] = 0
        
        ###---> Apply mapping to source data
        data_int = data_src[:, mapper_nodes]
    
        ###---> Fill seafloor with deepest active layer
        # Create a mask of nonzero entries
        mask = (data_int != 0)
        
        # Indices along axis 0
        idx = np.arange(data_int.shape[0])[:, None]
        
        # For each column, find the last index where data != 0
        # np.where(mask, idx, -1) gives -1 where data == 0, so max gives last nonzero index
        last_nonzero_idx = np.max(np.where(mask, idx, -1), axis=0)
        
        # Get the last nonzero values for each column
        last_nonzero_val = data_int[last_nonzero_idx, np.arange(data_int.shape[1])]
        
        # Fill zeros with the last nonzero value
        data_int_filled = np.where(mask, data_int, last_nonzero_val)
        
        ###---> Force cavity and bathymetry nodes back to 0
        data_int_filled[data_dst == 0] = 0

        ###---> For velocity variables, set in-cavity values to 0 to force cold start
        if varname in ['u', 'v', 'urhs_AB', 'vrhs_AB', 'urhs_AB3', 'vrhs_AB3', 'w', 'w_impl', 'w_expl']:
            print('Setting cavity velocities to 0 to force cavity cold start.')

            np.where(data_dst[0,:]) == 0
            data_int_filled[:, data_dst[0,:] == 0] = 0

        ###---> Make a deep copy of destination restart
        ds_int = ds_tgt.copy()
        
        ###---> Add a time dimension
        time = xr.open_dataset(f"{path_restart_src}{varname}.nc").time.values[t_step]
        ds_int = ds_int.assign_coords(time=time).expand_dims('time')
        
        ###---> Extract iter
        iteration = xr.open_dataset(f"{path_restart_src}{varname}.nc").iter.values[t_step]
        
        ###---> Write iter to deep copy
        ds_int['iter'] = (('time'), np.atleast_1d(iteration))
        
        ###---> Write intermediate data to deep copy
        if 'node' in ds_tgt.dims:
            horiz = 'node'
        elif 'elem' in ds_tgt.dims:
            horiz = 'elem'
        else:
            raise ValueError('Neither node nor elem are found in the dataset dimensions!')
        
        if "nz_1" in ds_tgt.dims:
            vert = 'nz_1'
        elif "nz" in ds_tgt.dims:
            vert = 'nz'
        else:
            raise ValueError('Neither nz1 nor nz are found in the dataset dimensions!')
        
        ds_int[varname] = (('time', vert, horiz), data_int_filled[np.newaxis,:])
        
        ###---> Add metadata
        for var in ds_src.data_vars:
            if var in ds_int:
                ds_int[var].attrs = ds_src[var].attrs.copy()
        
        ds_int.attrs['source restart (data)'] = f'{path_restart_src}{varname}.nc'
        ds_int.attrs['target restart (grid)'] = f'{path_restart_src}{varname}.nc'

        ###---> Final checks
        if varname not in ['u', 'v', 'urhs_AB', 'vrhs_AB', 'urhs_AB3', 'vrhs_AB3', 'w', 'w_impl', 'w_expl']:
            finite_sum = np.isfinite(ds_int[varname].values).sum()
            nansum = np.isnan(ds_int[varname].values).sum()
            zerosum = (ds_int[varname].values == 0).sum()
            
            if finite_sum != np.prod(np.array(ds_tgt[varname].values.shape)):
                raise ValueError('There are non-finite values in the dataset!')
            elif nansum != 0:
                raise ValueError('There are non-finite values in the dataset!')
            elif zerosum != (ds_tgt[varname].values == 0).sum():
                raise ValueError('There are too few cavity/topography values in the dataset!')
            
        ###---> Save
        ds_int.to_netcdf(f'{path_out}{varname}.nc')

    print(' ')
    return

def fill_cavities_from_existing_restart(varname, path_restart_tgt_oce, path_restart_cavity_fill_oce):
    """ 
    Replace the cavity values in the destination restart files with the values from the source restart files.

    Parameters
    ----------
    varname : str
        Name of the variable to fill.
    path_restart_tgt_oce : str
        Path to the target restart file.
    path_restart_cavity_fill_oce : str
        Path to the source restart file.
    """

    if verbose:
        print('============ fill_cavities_from_existing_restart.py ==============')
        print(f'Filling cavity values from existing restart file:')
        print(f'Source: {path_restart_cavity_fill_oce}{varname}.nc')
        print(f'Target: {path_restart_tgt_oce}{varname}.nc')
        print(f'Variable name: {varname}.nc')
    
    # Rename the originally computed restart file
    os.rename(f"{path_restart_tgt_oce}{varname}.nc", f"{path_restart_tgt_oce}{varname}_original.nc")
    
    # Load the dataset to be filled (just created)
    ds_to_fill = xr.open_dataset(f"{path_restart_tgt_oce}{varname}_original.nc").isel(time=-1).squeeze()

    # Load the dataset to fill from (another simulation)
    ds_fill = xr.open_dataset(f"{path_restart_cavity_fill_oce}{varname}.nc").isel(time=-1).squeeze()

    # Build a horizontal cavity mask by checking if the first layer is finite
    cavity_mask = np.isfinite(ds_to_fill[varname].values[0,:])
    print(f'Number of cavity cells: {cavity_mask.sum()}')

    # Replace the cavity values with the values from the fill dataset
    ds_to_fill[varname].values[:, cavity_mask] = ds_fill[varname].values[:, cavity_mask]

    # Save the filled dataset
    ds_to_fill.to_netcdf(f"{path_restart_tgt_oce}{varname}.nc")
    
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


    so.create_map(ax, extent=box, land=False, coastline=True, circular=False)

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
        so.create_map(axis, extent=box, land=True, coastline=True, circular=False)

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
    

if __name__ == "__main__":
    print('=============================================================================')
    print('======================== EXTRAPOLATE CAVITY RESTARTS ========================')
    print('=============================================================================')
    print(' ')

    # =============================================================================
    # ============================ READ MESHES =================================
    # =============================================================================

    # DARS2/CORE2
    node_lon_src, node_lat_src, node_id_src, node_coastmask_src = so.read_nodes(path_mesh_src)
    elem_src = so.read_elements(path_mesh_src)
    elem_lon_src = node_lon_src[elem_src].mean(axis=1)
    elem_lat_src = node_lat_src[elem_src].mean(axis=1)

    # DARS2cav/CORE2ice
    node_lon_tgt, node_lat_tgt, node_id_tgt, node_coastmask_tgt = so.read_nodes(path_mesh_tgt)
    elem_tgt = so.read_elements(path_mesh_tgt)
    elem_lon_tgt = node_lon_tgt[elem_tgt].mean(axis=1)
    elem_lat_tgt = node_lat_tgt[elem_tgt].mean(axis=1)


    # =============================================================================
    # ============================ BUILD MAPPERS ==================================
    # =============================================================================
    print("###---> Building Nearest Neighbor Mappers...")
    print(" ")
    mapper_elements, distance_elements = build_spherical_nn_mapper(elem_lon_src, elem_lat_src, elem_lon_tgt, elem_lat_tgt)
    mapper_nodes, distance_nodes = build_spherical_nn_mapper(node_lon_src, node_lat_src, node_lon_tgt, node_lat_tgt)

    if plot:
        print("###---> Plotting Nearest Neighbor Mappers...")
        print(" ")
        plot_mapper(mapper_elements, elem_lon_src, elem_lat_src, elem_lon_tgt, elem_lat_tgt, 'elem', path_dst_plots)
        plot_mapper(mapper_nodes, node_lon_src, node_lat_src, node_lon_tgt, node_lat_tgt, 'node', path_dst_plots)

    # =============================================================================
    # ======================== 2D nodal fiels (time, node) ========================
    # =============================================================================
    varnames_2D_node = ['ssh', 'ssh_rhs_old', 'hbar']
    mask_file = f"{path_restart_tgt_oce_v27}ssh.nc"

    for varname in varnames_2D_node:
        interpolate_extrapolate_2D(varname, path_restart_src_oce, path_restart_tgt_oce, mapper_nodes, path_dst_restarts_oce, mask_file, mask_file_varname='ssh', t_step=-1, verbose=True)
        if plot:
            plot_interpolated_extrapolated_field(path_restart_src_oce, path_dst_restarts_oce, path_restart_tgt_oce, varname, node_lon_src, node_lat_src, node_lon_tgt, node_lat_tgt, path_dst_plots)

    varnames_2D_node = ['area', 'hsnow', 'hice', 'uice', 'vice']
    if is_coupled:
        varnames_2D_node.extend(['ice_temp', 'ice_albedo'])

    for varname in varnames_2D_node:
        interpolate_extrapolate_2D(varname, path_restart_src_ice, path_restart_tgt_ice, mapper_nodes, path_dst_restarts_ice, mask_file, mask_file_varname='ssh', t_step=-1, verbose=True)
        if plot:
            plot_interpolated_extrapolated_field(path_restart_src_ice, path_dst_restarts_ice, path_restart_tgt_ice, varname, node_lon_src, node_lat_src, node_lon_tgt, node_lat_tgt, path_dst_plots)

    # =============================================================================
    # ======================== 3D nodal/elem fiels (time, nz_1/nz, node/elem) =====
    # =============================================================================
    varnames_3D_node = ['hnode', 'salt', 'temp', 'temp_AB', 'salt_AB', 'temp_M1', 'salt_M1', 'w', 'w_impl', 'w_expl']

    for varname in varnames_3D_node:
        interpolate_extrapolate_3D(varname, path_restart_src_oce, path_restart_tgt_oce, mapper_nodes, path_dst_restarts_oce, t_step=-1, verbose=True)
        if plot:
            plot_interpolated_extrapolated_field(path_restart_src_oce, path_dst_restarts_oce, path_restart_tgt_oce, varname, node_lon_src, node_lat_src, node_lon_tgt, node_lat_tgt, path_dst_plots)


    varnames_3D_element = ['u', 'v', 'vrhs_AB', 'urhs_AB', 'urhs_AB3', 'vrhs_AB3']
    for varname in varnames_3D_element:
        interpolate_extrapolate_3D(varname, path_restart_src_oce, path_restart_tgt_oce, mapper_elements, path_dst_restarts_oce, t_step=-1, verbose=True)
        if plot:
            plot_interpolated_extrapolated_field(path_restart_src_oce, path_dst_restarts_oce, path_restart_tgt_oce, varname, elem_lon_src, elem_lat_src, elem_lon_tgt, elem_lat_tgt, path_dst_plots)

    # =============================================================================
    # ============== Fill cavities from existing restart files  ===================
    # =============================================================================
    if fill_cavities:
        varnames_3D = ['salt', 'temp', 'temp_AB', 'salt_AB', 'temp_M1', 'salt_M1', 'w', 'w_impl', 'w_expl', 'u', 'v', 'vrhs_AB', 'urhs_AB', 'urhs_AB3', 'vrhs_AB3']
        for varname in varnames_3D:
            fill_cavities_from_existing_restart(varname, path_restart_tgt_oce, path_restart_cavity_fill_oce)

    print('=============================================================================')
    print('======================== EXTRAPOLATION COMPLETE ============================')
    print('=============================================================================')