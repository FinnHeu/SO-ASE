# so_ase/fesom2/helpers_extraction.py

import xarray as xr
import numpy as np
from os.path import isfile

def fesom_extract_bottom_values(
    src_path,
    dest_path,
    varname='temp',
    years=(1979, 2015),
    log=True
):
    """
    Extract bottom (deepest valid) values from 3D FESOM output and save as 2D fields.

    For each grid node and time step, finds the deepest non-NaN value in the water
    column. This handles both open ocean profiles (valid values followed by NaN 
    bathymetry) and cavity profiles (NaN ice shelf, valid water, NaN bathymetry).

    Parameters
    ----------
    src_path : str
        Path to the directory containing FESOM output files.
        Files are expected as `{varname}.fesom.{year}.nc`.
    dest_path : str
        Path to the directory where output files will be saved.
        Output files are named `bottom_{varname}.fesom.{year}.nc`.
    varname : str, optional
        Variable name to process. Supported: 'temp', 'salt', 'sigma0'.
        Default is 'temp'.
    years : tuple of int, optional
        Year range (start, end) to process. Default is (1979, 2015).
    log : bool, optional
        If True, print progress messages. Default is True.

    Returns
    -------
    None
        Output is written to NetCDF files at dest_path.

    Notes
    -----
    - Existing output files are skipped (not overwritten).
    - The function handles profiles with ice shelf cavities where the top
      of the water column may also be NaN.
    - If a profile has no valid values (all NaN), the bottom value is NaN.

    Example
    -------
    >>> fesom_extract_bottom_values(
    ...     src_path='/path/to/fesom/output/',
    ...     dest_path='/path/to/output/',
    ...     varname='temp',
    ...     years=(2000, 2010)
    ... )
    """
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)

    for year in range(years[0], years[-1] + 1):
        input_file = f"{src_path}{varname}.fesom.{year}.nc"
        output_file = f"{dest_path}bottom_{varname}.fesom.{year}.nc"

        if isfile(output_file):
            if log:
                print(f"Skipping (exists): {output_file}", flush=True)
            continue

        if not isfile(input_file):
            if log:
                print(f"Input file not found: {input_file}", flush=True)
            continue

        if log:
            print(f"Processing: {input_file}", flush=True)

        ds = xr.open_dataset(input_file, decode_times=time_coder)
        
        # Get dimension names and find their positions
        var_dims = ds[varname].dims  # e.g., ('time', 'nod2', 'nz')
        depth_dim = [d for d in var_dims if d.startswith('nz')][0]
        node_dim = [d for d in var_dims if d.startswith('nod')][0]
        depth_axis = var_dims.index(depth_dim)
        
        ntime = ds.dims['time']
        nz = ds.dims[depth_dim]
        nnodes = ds.dims[node_dim]
        
        data = ds[varname].values  # shape depends on file, e.g., (time, nod2, nz)

        # Find bottom values: last non-NaN along depth axis
        # Create mask of valid (non-NaN) values
        valid_mask = np.isfinite(data)

        # Find the index of the last valid value along depth axis
        # We reverse along depth, find first valid, then convert back to original index
        reversed_mask = np.flip(valid_mask, axis=depth_axis)
        first_valid_reversed = np.argmax(reversed_mask, axis=depth_axis)
        has_valid = np.any(valid_mask, axis=depth_axis)

        # Convert reversed index to original index
        bottom_idx = (nz - 1) - first_valid_reversed  # shape: (time, nod2)

        # Extract bottom values using take_along_axis
        bottom_idx_expanded = np.expand_dims(bottom_idx, axis=depth_axis)
        bottom_values = np.take_along_axis(data, bottom_idx_expanded, axis=depth_axis)
        bottom_values = np.squeeze(bottom_values, axis=depth_axis)

        # Set to NaN where no valid values exist
        bottom_values[~has_valid] = np.nan

        # Get nod2 coordinate from dataset
        if node_dim in ds.coords:
            nod2_coord = ds[node_dim]
        elif node_dim in ds:
            nod2_coord = ds[node_dim]
        else:
            nod2_coord = np.arange(nnodes)

        # Create output dataset
        ds_out = xr.Dataset(
            {
                f"bottom_{varname}": (["time", node_dim], bottom_values.astype(np.float32)),
            },
            coords={
                "time": ds.time,
                node_dim: nod2_coord,
            },
            attrs={
                "description": f"Bottom (deepest valid) {varname} values from FESOM output",
                "source_file": input_file,
            }
        )

        ds_out.to_netcdf(output_file)
        if log:
            print(f"Saved: {output_file}", flush=True)

        ds.close()

    if log:
        print("Done.", flush=True)


def fesom_depth_average(
    src_path,
    dest_path,
    mesh_diag_path,
    varname='temp',
    depth_range=(0, 50),
    years=(1979, 2015),
    log=True
):
    """
    Compute layer-thickness-weighted depth averages from 3D FESOM output.

    For each grid node and time step, computes the weighted mean of the variable
    over the specified depth range, using layer thickness as weights.

    Parameters
    ----------
    src_path : str
        Path to the directory containing FESOM output files.
        Files are expected as `{varname}.fesom.{year}.nc`.
    dest_path : str
        Path to the directory where output files will be saved.
        Output files are named `{varname}_{depth_min}m_{depth_max}m.fesom.{year}.nc`.
    mesh_diag_path : str
        Path to the directory containing `fesom.mesh.diag.nc`.
    varname : str, optional
        Variable name to process. Default is 'temp'.
    depth_range : tuple of float, optional
        Depth range (min, max) in meters for averaging. Default is (0, 50).
    years : tuple of int, optional
        Year range (start, end) to process. Default is (1979, 2015).
    log : bool, optional
        If True, print progress messages. Default is True.

    Returns
    -------
    None
        Output is written to NetCDF files at dest_path.

    Notes
    -----
    - Existing output files are skipped (not overwritten).
    - Layer thickness is computed from mesh_diag.nz (depth levels).
    - The depth selection uses the nz coordinate values in meters.

    Example
    -------
    >>> fesom_depth_average(
    ...     src_path='/path/to/fesom/output/',
    ...     dest_path='/path/to/output/',
    ...     mesh_diag_path='/path/to/mesh/',
    ...     varname='temp',
    ...     depth_range=(0, 50),
    ...     years=(2000, 2010)
    ... )
    """
    import os
    os.makedirs(dest_path, exist_ok=True)
    
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
    
    # Load mesh diagnostics for depth levels
    mesh_diag = xr.open_dataset(f"{mesh_diag_path}fesom.mesh.diag.nc")
    
    depth_min, depth_max = depth_range
    
    for year in range(years[0], years[-1] + 1):
        input_file = f"{src_path}{varname}.fesom.{year}.nc"
        output_file = f"{dest_path}{varname}_{int(depth_min)}m_{int(depth_max)}m.fesom.{year}.nc"

        if isfile(output_file):
            if log:
                print(f"Skipping (exists): {output_file}", flush=True)
            continue

        if not isfile(input_file):
            if log:
                print(f"Input file not found: {input_file}", flush=True)
            continue

        if log:
            print(f"Processing: {input_file}", flush=True)

        ds = xr.open_dataset(input_file, decode_times=time_coder)
        
        # Get depth coordinate name (nz or nz1)
        depth_dim = None
        for dim in ds[varname].dims:
            if dim.startswith('nz'):
                depth_dim = dim
                break
        if depth_dim is None:
            raise ValueError(f"Could not find depth dimension (nz*) in {ds[varname].dims}")
        
        # Compute layer thickness from mesh_diag.nz (always use this)
        # mesh_diag.nz has 57 levels, np.diff gives 56 layer thicknesses matching the data
        layer_thickness = np.abs(np.diff(mesh_diag.nz.values))
        
        # Create DataArray for layer thickness
        ds['layer_thickness'] = (depth_dim, layer_thickness)
        
        # Get depth values from the dataset for masking (56 levels)
        depths = ds[depth_dim].values
        
        # Select depth range and compute weighted mean
        depth_mask = (depths >= depth_min) & (depths <= depth_max)
        ds_subset = ds.isel({depth_dim: depth_mask})
        
        if ds_subset.dims[depth_dim] == 0:
            if log:
                print(f"Warning: No depth levels in range {depth_min}-{depth_max}m", flush=True)
            ds.close()
            continue
        
        # Compute weighted mean
        depth_avg = ds_subset[varname].weighted(ds_subset['layer_thickness']).mean(dim=depth_dim)

        # Create output dataset
        out_varname = f"{varname}_{int(depth_min)}m_{int(depth_max)}m"
        ds_out = xr.Dataset(
            {
                out_varname: depth_avg.astype(np.float32),
            },
            attrs={
                "description": f"Layer-thickness-weighted depth average of {varname} from {depth_min}m to {depth_max}m",
                "source_file": input_file,
                "depth_range_m": f"{depth_min}-{depth_max}",
            }
        )

        ds_out.to_netcdf(output_file)
        if log:
            print(f"Saved: {output_file}", flush=True)

        ds.close()

    mesh_diag.close()
    
    if log:
        print("Done.", flush=True)


def fesom_extract_at_depth(
    src_path,
    dest_path,
    mesh_diag_path,
    varname='temp',
    depth=50,
    years=(1979, 2015),
    log=True
):
    """
    Extract FESOM2 data at a specific depth level and save to NetCDF files.
    
    Uses nearest-neighbor selection to find the closest available depth level
    to the requested depth. The actual extracted depth is included in the
    output filename.
    
    Parameters
    ----------
    src_path : str
        Path to directory containing FESOM2 output files named
        `{varname}.fesom.{year}.nc`.
    dest_path : str
        Path to directory where output files will be saved.
    mesh_diag_path : str
        Path to directory containing `fesom.mesh.diag.nc`.
    varname : str, optional
        Variable name to extract (e.g., 'temp', 'salt'). Default is 'temp'.
    depth : float, optional
        Target depth in meters. The nearest available depth level will be
        selected. Default is 50.
    years : tuple of int, optional
        Start and end year (inclusive) to process. Default is (1979, 2015).
    log : bool, optional
        If True, print progress messages. Default is True.
    
    Returns
    -------
    None
        Results are written to NetCDF files named
        `{varname}_{actual_depth}m.fesom.{year}.nc`.
    
    Example
    -------
    >>> fesom_extract_at_depth(
    ...     src_path='/path/to/fesom/output/',
    ...     dest_path='/path/to/output/',
    ...     mesh_diag_path='/path/to/mesh/',
    ...     varname='temp',
    ...     depth=100,
    ...     years=(2000, 2010)
    ... )
    """
    import os
    os.makedirs(dest_path, exist_ok=True)
    
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
    
    for year in range(years[0], years[-1] + 1):
        input_file = f"{src_path}{varname}.fesom.{year}.nc"
        
        if not isfile(input_file):
            if log:
                print(f"Input file not found: {input_file}", flush=True)
            continue
        
        ds = xr.open_dataset(input_file, decode_times=time_coder)
        
        # Get depth coordinate name (nz or nz1)
        depth_dim = None
        for dim in ds[varname].dims:
            if dim.startswith('nz'):
                depth_dim = dim
                break
        if depth_dim is None:
            raise ValueError(f"Could not find depth dimension (nz*) in {ds[varname].dims}")
        
        # Select nearest depth level
        ds_at_depth = ds.sel({depth_dim: depth}, method='nearest')
        
        # Get the actual depth value that was selected
        actual_depth = float(ds_at_depth[depth_dim].values)
        actual_depth_int = int(round(actual_depth))
        
        # Build output filename with actual depth
        output_file = f"{dest_path}{varname}_{actual_depth_int}m.fesom.{year}.nc"
        
        if isfile(output_file):
            if log:
                print(f"Skipping (exists): {output_file}", flush=True)
            ds.close()
            continue
        
        if log:
            print(f"Processing: {input_file} -> depth {actual_depth_int}m", flush=True)
        
        # Extract the variable at the selected depth
        data_at_depth = ds_at_depth[varname]
        
        # Create output dataset
        out_varname = f"{varname}_{actual_depth_int}m"
        ds_out = xr.Dataset(
            {
                out_varname: data_at_depth.astype(np.float32),
            },
            attrs={
                "description": f"{varname} extracted at {actual_depth_int}m depth",
                "source_file": input_file,
                "requested_depth_m": depth,
                "actual_depth_m": actual_depth,
            }
        )
        
        ds_out.to_netcdf(output_file)
        if log:
            print(f"Saved: {output_file}", flush=True)
        
        ds.close()
    
    if log:
        print("Done.", flush=True)


def fesom_region_mean_vertical_profile(
    src_path,
    dest_path,
    mesh_diag_path,
    mask_file,
    varname='temp',
    years=(1979, 2015),
    log=True
):
    """
    Compute area-weighted mean vertical profiles for a masked region from FESOM2 output.
    
    For each year and time step, computes the area-weighted horizontal mean of 
    vertical profiles within the masked region. Outputs time x depth datasets.
    
    Parameters
    ----------
    src_path : str
        Path to directory containing FESOM2 output files named
        `{varname}.fesom.{year}.nc`.
    dest_path : str
        Path to directory where output files will be saved.
    mesh_diag_path : str
        Path to directory containing `fesom.mesh.diag.nc`.
    mask_file : str
        Path to NetCDF file containing a 1D binary mask variable.
        The mask should have dimension matching nod2 with values 0/1 or True/False.
    varname : str, optional
        Variable name to extract (e.g., 'temp', 'salt'). Default is 'temp'.
    years : tuple of int, optional
        Start and end year (inclusive) to process. Default is (1979, 2015).
    log : bool, optional
        If True, print progress messages. Default is True.
    
    Returns
    -------
    None
        Results are written to NetCDF files named
        `{varname}_region_profile.{mask_name}.{year}.nc`.
    
    Example
    -------
    >>> fesom_region_mean_vertical_profile(
    ...     src_path='/path/to/fesom/output/',
    ...     dest_path='/path/to/output/',
    ...     mesh_diag_path='/path/to/mesh/',
    ...     mask_file='/path/to/mesh/mask.AAshelf.max1000m.nc',
    ...     varname='temp',
    ...     years=(2000, 2010)
    ... )
    """
    import os
    os.makedirs(dest_path, exist_ok=True)
    
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
    
    # Load mesh diagnostics for nodal areas
    mesh_diag = xr.open_dataset(f"{mesh_diag_path}fesom.mesh.diag.nc")
    
    # Load the mask file and extract mask variable
    mask_ds = xr.open_dataset(mask_file)
    mask_varname = 'zbar_n_bottom'
    mask = mask_ds[mask_varname].values
    
    # Get mask name from filename for output naming
    mask_basename = os.path.basename(mask_file)
    mask_name = mask_basename.replace('.nc', '').replace('mask.', '')
    
    # Find indices where mask is True/1
    inds = np.where(mask > 0)[0]
    if log:
        print(f"Mask '{mask_name}': {len(inds)} nodes selected", flush=True)
    
    # Get nodal area for weighting (max over nz dimension, then select nodes)
    nodal_area = mesh_diag.nod_area.max(dim='nz').isel(nod2=inds)
    
    for year in range(years[0], years[-1] + 1):
        input_file = f"{src_path}{varname}.fesom.{year}.nc"
        output_file = f"{dest_path}{varname}_region_profile.{mask_name}.{year}.nc"
        
        if isfile(output_file):
            if log:
                print(f"Skipping (exists): {output_file}", flush=True)
            continue
        
        if not isfile(input_file):
            if log:
                print(f"Input file not found: {input_file}", flush=True)
            continue
        
        if log:
            print(f"Processing: {input_file}", flush=True)
        
        ds = xr.open_mfdataset(input_file, decode_times=time_coder)
        
        # Select nodes within the masked region
        ds_region = ds.isel(nod2=inds)
        
        # Compute area-weighted horizontal mean (preserves time and nz dimensions)
        mean_profile = ds_region[varname].weighted(nodal_area).mean(dim='nod2')
        
        # Create output dataset
        ds_out = xr.Dataset(
            {
                varname: mean_profile.astype(np.float32),
            },
            attrs={
                "description": f"Area-weighted mean vertical profile of {varname} in region {mask_name}",
                "source_file": input_file,
                "mask_file": mask_file,
                "mask_name": mask_name,
            }
        )
        
        ds_out.to_netcdf(output_file)
        if log:
            print(f"Saved: {output_file}", flush=True)
        
        ds.close()
    
    mesh_diag.close()
    mask_ds.close()
    
    if log:
        print("Done.", flush=True)