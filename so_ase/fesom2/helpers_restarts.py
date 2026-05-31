# so_ase/fesom2/helpers_restarts.py

import xarray as xr
import numpy as np
from scipy.spatial import cKDTree


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


def fill_cavities_from_existing_restart(varname, path_restart_src, path_restart_dst, path_restart_cavity_fill, path_mesh, verbose=True):
    """ 
    Replace the cavity values in the restart files with the values from another restart files.

    Parameters
    ----------
    varname : str
        Name of the variable to fill.
    path_restart_src : str
        Path to the source restart file.
    path_restart_dst : str
        Path to the destination restart file.
    path_restart_cavity_fill : str
        Path to the cavity fill source restart file.
    path_mesh : str
        Path to the mesh directory containing fesom.mesh.diag.nc.
    verbose : bool, optional
        Whether to print verbose output, by default True
    """

    if verbose:
        print('\n============ fill_cavities_from_existing_restart.py ==============')
        print(f'Filling cavity values from existing restart file:')
        print(f'Source: {path_restart_src}{varname}.nc')
        print(f'Cavity fill source: {path_restart_cavity_fill}{varname}.nc')
        print(f'Target: {path_restart_dst}{varname}.nc')
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
