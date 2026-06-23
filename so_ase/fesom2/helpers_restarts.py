# so_ase/fesom2/helpers_restarts.py

import xarray as xr
import numpy as np
from scipy.spatial import cKDTree
import os
from .helpers_mesh import read_nodes, read_elements, read_element_levels, build_cavity_mask, find_element_for_points    
from .eval_icebergs import read_iceberg_restart_file

###--------> 1. FESOM ocean/ice restart files

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


###--------> 2. FESOM iceberg restart files

def icbdat_from_ismrestart(srcpath, dstpath, meshpath=None, verbose=True):
    """
    Convert an iceberg.restart.ISM file into a set of corresponding initial icb_*.dat files.
    
    Parameters
    ----------
    srcpath : str
        Path to the directory containing iceberg.restart.ISM file.
    dstpath : str
        Path to the output directory for icb_*.dat files.
    meshpath : str, optional
        Path to the target mesh directory. If provided, icebergs are filtered to only
        include those in valid open ocean elements (not cavity, not coastal nodes).
    verbose : bool, optional
        Whether to print progress information, by default True.
    
    Notes
    -----
    Creates the following files in dstpath:
        - icb_longitude.dat : Longitude in geographical coordinates
        - icb_latitude.dat : Latitude in geographical coordinates
        - icb_height.dat : Iceberg height
        - icb_length.dat : Iceberg length
        - icb_scaling.dat : Scaling factor
        - icb_calving_day.dat : Calving day
    """
    
    if verbose:
        print('============ icbdat_from_restart ==============')
        print(f'Source: {srcpath}/iceberg.restart.ISM')
        print(f'Destination: {dstpath}')
        if meshpath:
            print(f'Target mesh: {meshpath}')
    
    # Read the restart file (unrotate coordinates for output)
    ds_icb = read_iceberg_restart_file(f"{srcpath}/iceberg.restart.ISM", unrotate=True)
    
    n_icebergs = len(ds_icb.ib)
    if verbose:
        print(f'Number of icebergs in restart: {n_icebergs}')
    
    # Optional mesh element check
    if meshpath is not None:
        if verbose:
            print('Checking iceberg locations on target mesh...')
        
        # Read mesh data
        node_lon, node_lat, node_idx, node_coast = read_nodes(meshpath)
        elements = read_elements(meshpath)
        cavity_elvls = read_element_levels(meshpath, which='cavity', raw=False, python_indexing=False)
        
        # Build node cavity mask (node is cavity if all its elements are cavity)
        node_cavity_mask = build_cavity_mask(meshpath, which='node')
        
        # Load forbidden elements from cache file if available
        # For restart conversion, we want icebergs only in pure open ocean:
        # NOT in cavity elements, NOT in calving front elements, NOT in seeding elements
        cache_file = os.path.join(meshpath, "elem_icb_seeding.npz")
        forbidden_elems = set()
        if os.path.exists(cache_file):
            if verbose:
                print(f'  Loading seeding constraints from: {cache_file}')
            cache_data = np.load(cache_file, allow_pickle=True)
            # Calving front elements are forbidden
            calving_front_elems = set(cache_data['all_calving_front'])
            # Seeding elements (near calving front) are also forbidden for restart conversion
            seeding_elems = set(cache_data['all_seeding_elems'])
            # Combine both into forbidden set
            forbidden_elems = calving_front_elems | seeding_elems
            if verbose:
                print(f'  Calving front elements: {len(calving_front_elems)}')
                print(f'  Seeding elements (near calving front): {len(seeding_elems)}')
                print(f'  Total forbidden elements: {len(forbidden_elems)}')
        else:
            if verbose:
                print(f'  Warning: Cache file not found: {cache_file}')
                print(f'  Skipping forbidden element check.')
        
        # Find element for each iceberg
        elem_indices = find_element_for_points(
            ds_icb['lon_deg'].values,
            ds_icb['lat_deg'].values,
            node_lon, node_lat, elements
        )
        
        # Check validity of each iceberg location
        # Only pure open ocean elements are allowed (not cavity, not calving front, not seeding)
        valid_mask = np.ones(n_icebergs, dtype=bool)
        n_not_in_elem = 0
        n_in_cavity = 0
        n_near_coast_cavity = 0
        n_in_forbidden = 0
        
        for i in range(n_icebergs):
            eidx = elem_indices[i]
            
            # Check if iceberg is in any element
            if eidx == -1:
                valid_mask[i] = False
                n_not_in_elem += 1
                continue
            
            # Check if element is not a cavity element (cavity_elvls == 1 means open ocean)
            if cavity_elvls[eidx] > 1:
                valid_mask[i] = False
                n_in_cavity += 1
                continue
            
            # Check if all three nodes are valid (not coastal and not cavity)
            n1, n2, n3 = elements[eidx]
            node_invalid = False
            for node in [n1, n2, n3]:
                # Check if node is coastal
                if node_coast[node] == 1:
                    node_invalid = True
                    break
                # Check if node is in cavity
                if node_cavity_mask[node]:
                    node_invalid = True
                    break
            if node_invalid:
                valid_mask[i] = False
                n_near_coast_cavity += 1
                continue
            
            # Check if element is in forbidden elements (calving front or seeding elements)
            if eidx in forbidden_elems:
                valid_mask[i] = False
                n_in_forbidden += 1
                continue
        
        # Filter dataset to valid icebergs
        n_valid = valid_mask.sum()
        
        if verbose:
            print(f'  Icebergs not in any element: {n_not_in_elem}')
            print(f'  Icebergs in cavity elements: {n_in_cavity}')
            print(f'  Icebergs near coast/cavity nodes: {n_near_coast_cavity}')
            print(f'  Icebergs in forbidden elements (calving front + seeding): {n_in_forbidden}')
            print(f'  Valid icebergs (pure open ocean): {n_valid} / {n_icebergs}')
        
        # Apply filter
        ds_icb = ds_icb.isel(ib=valid_mask)
        n_icebergs = n_valid
    
    # Create output directory if it doesn't exist
    os.makedirs(dstpath, exist_ok=True)
    
    # Write icb_longitude.dat (geographical coordinates)
    with open(f'{dstpath}/icb_longitude.dat', 'w') as f:
        for val in ds_icb['lon_deg'].values:
            f.write(f'{val}\n')
    if verbose:
        print(f'  Written: icb_longitude.dat')
    
    # Write icb_latitude.dat (geographical coordinates)
    with open(f'{dstpath}/icb_latitude.dat', 'w') as f:
        for val in ds_icb['lat_deg'].values:
            f.write(f'{val}\n')
    if verbose:
        print(f'  Written: icb_latitude.dat')
    
    # Write icb_height.dat
    with open(f'{dstpath}/icb_height.dat', 'w') as f:
        for val in ds_icb['height_ib'].values:
            f.write(f'{val}\n')
    if verbose:
        print(f'  Written: icb_height.dat')
    
    # Write icb_length.dat
    with open(f'{dstpath}/icb_length.dat', 'w') as f:
        for val in ds_icb['length_ib'].values:
            f.write(f'{val}\n')
    if verbose:
        print(f'  Written: icb_length.dat')
    
    # Write icb_scaling.dat
    with open(f'{dstpath}/icb_scaling.dat', 'w') as f:
        for val in ds_icb['scaling'].values:
            f.write(f'{int(val)}\n')
    if verbose:
        print(f'  Written: icb_scaling.dat')
    
    # Write icb_calving_day.dat
    with open(f'{dstpath}/icb_calving_day.dat', 'w') as f:
        for val in ds_icb['calving_day'].values:
            f.write(f'{int(val)}\n')
    if verbose:
        print(f'  Written: icb_calving_day.dat')
    
    if verbose:
        print('Done.')
    
    return