# so_ase/fesom2/eval_cavity.py

import os
import xarray as xr
import numpy as np
import glob
from os.path import isfile
from .helpers_mesh import build_cavity_mask, build_cavity_regional_mask, get_first_active_layer_area
from ..miscellaneous.helpers_misc import seconds_per_month

def fesom_subshelf_freshwaterflux(src_path, mesh_diag_path, mesh_path, mask, year, log=False, savepath='./'):
    """
    Compute and save integrated subshelf basal melt time series from FESOM 
    freshwater flux output [m3/s].

    This function loads monthly freshwater flux fields (`fw`) from the FESOM
    output file of a single year, multiplies them by the mean horizontal node
    area from a mesh diagnostics file, applies a cavity mask, and integrates the
    melt signal over all masked nodes to obtain a subshelf melt time series for
    that year. The result is saved as a standalone NetCDF file.

    Parameters
    ----------
    src_path : str
        Directory containing the annual FESOM freshwater flux files 
        named ``fw.fesom.<year>.nc``.
    mesh_diag_path : str
        Directory containing ``fesom.mesh.diag.nc`` with `nod_area` and 
        related mesh diagnostics.
    mesh_path : str
        Path to the FESOM mesh directory (used for building the cavity mask).
    mask : dict or array-like
        Identifier for which ocean cavity mask to apply.  
        If ``dict['name'] == 'all'`` is given, the mask is generated using 
        ``build_cavity_mask(mesh_path, which='node')``.  
        If ``dict['name'] == 'Amery'`` is given a .kml file called ``Amery.kml`` 
        will be loaded from ``dict['kml_path']`` from which the boolean mask will be generated.
        If an array-like mask is supplied, it is used directly for ``isel``.
    year : int
        Year to process.
    log : bool, default False
        If True, print a message whenever a file is saved or skipped.
    savepath : str, default './'
        Directory in which to save the output files. The files are created as
        ``subshelf_melt_<mask>.<year>.nc``.
    """
    
    print("\n--> Compute fesom subshelf freshwaterflux...")

    out_file = f"{savepath}subshelf_melt_{mask['name']}.{year}.nc"
    if os.path.exists(out_file):
        if log:
            print(f"Skipping existing file: {out_file}")
        return

    os.makedirs(savepath, exist_ok=True)

    # Build mask
    mesh_diag = xr.open_dataset(f"{mesh_diag_path}fesom.mesh.diag.nc")
    if isinstance(mask, dict):
        if mask['name'] == 'all':
            node_mask = build_cavity_mask(mesh_path, which='node')
        else:
            node_mask = build_cavity_regional_mask(mesh_path, mask['kml_path'], name=mask['name'], which='node')
    elif isinstance(mask, list) or isinstance(mask, np.array):
        node_mask = mask
    else:
        raise ValueError('Mask type is not supported!')

    # Open file with cftime decoder
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
    ds_fw = xr.open_mfdataset([f"{src_path}fw.fesom.{year}.nc"], decode_times=time_coder, chunks={"time": 12})

    # Compute subshelf melt
    SSM = (ds_fw.fw * get_first_active_layer_area(mesh_diag)).isel(nod2=node_mask).sum(dim='nod2')

    # Create output dataset
    ds_out = xr.Dataset(
        {
            "subshelf_melt": SSM
        },
        coords={
            "time": SSM.time
        }
    )

    # Add units attribute
    ds_out["subshelf_melt"].attrs = {
        "units": "m3/s",
        "description": "sum(subshelf_melt(i) [m/s] * area(i) [m2]) over cavity region"
    }

    # Select the monthly data of this year and save
    ds_year = ds_out.sel(time=ds_out.time.dt.year == year)
    ds_year.to_netcdf(out_file)
    if log:
        print(f"Saved: {out_file}")

    if log:
        print("All done!")

    return

def fesom_subshelf_heatflux(src_path, mesh_diag_path, mesh_path, mask, year, log=False, savepath='./'):

    """
    Compute and save integrated subshelf heat flux time series from FESOM 
    heat flux output.

    This function loads monthly subshelf heat flux fields (`fh`) from the
    annual FESOM output file of a single year, multiplies them by the mean
    horizontal node area from a mesh diagnostics file, applies a cavity mask,
    and integrates the heat flux over all masked nodes to obtain an
    area-integrated time series for that year. The result is saved as a
    standalone NetCDF file.

    Parameters
    ----------
    src_path : str
        Directory containing the annual FESOM heat flux files 
        named ``fh.fesom.<year>.nc``.
    mesh_diag_path : str
        Directory containing ``fesom.mesh.diag.nc``, which includes 
        `nod_area` and other mesh diagnostics.
    mesh_path : str
        Path to the FESOM mesh directory. Used to construct the subshelf 
        cavity mask.
    mask : dict or array-like
        Identifier for which ocean cavity mask to apply.  
        If ``dict['name'] == 'all'`` is given, the mask is generated using 
        ``build_cavity_mask(mesh_path, which='node')``.  
        If ``dict['name'] == 'Amery'`` is given a .kml file called ``Amery.kml`` 
        will be loaded from ``dict['kml_path']`` from which the boolean mask will be generated.
        If an array-like mask is supplied, it is used directly for ``isel``.
    year : int
        Year to process.
    log : bool, default False
        If True, print status messages when saving or skipping files.
    savepath : str, default './'
        Directory into which output files will be written.  
        Files are saved as ``subshelf_heatflux_<mask>.<year>.nc``.

    Output Variables
    ----------------
    subshelf_heatflux : DataArray (time)
        Integrated subshelf heat flux time series for the given year.

    Returns
    -------
    None
        The function writes NetCDF files but does not return a value.
    """
    print("\n--> Compute fesom subshelf heatflux...")

    out_file = f"{savepath}subshelf_heatflux_{mask['name']}.{year}.nc"
    if os.path.exists(out_file):
        if log:
            print(f"Skipping existing file: {out_file}")
        return

    os.makedirs(savepath, exist_ok=True)

    # load mesh diagnostics
    mesh_diag = xr.open_dataset(f"{mesh_diag_path}fesom.mesh.diag.nc")

    # Build mask
    if isinstance(mask, dict):
        if mask['name'] == 'all':
            node_mask = build_cavity_mask(mesh_path, which='node')
        else:
            node_mask = build_cavity_regional_mask(mesh_path, mask['kml_path'], name=mask['name'], which='node')
    elif isinstance(mask, list) or isinstance(mask, np.array):
        node_mask = mask
    else:
        raise ValueError('Mask type is not supported!')

    # Open file with cftime decoder
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
    ds_fh = xr.open_mfdataset([f"{src_path}fh.fesom.{year}.nc"], decode_times=time_coder, chunks={"time": 12})

    # Compute subshelf heatflux
    SSHF = (ds_fh.fh * get_first_active_layer_area(mesh_diag)).isel(nod2=node_mask).sum(dim='nod2')

    # Create output dataset
    ds_out = xr.Dataset(
        {
            "subshelf_heatflux": SSHF
        },
        coords={
            "time": SSHF.time
        }
    )

    # Select the monthly data of this year and save
    ds_year = ds_out.sel(time=ds_out.time.dt.year == year)
    ds_year.to_netcdf(out_file)
    if log:
        print(f"Saved: {out_file}")

    if log:
        print("All done!")

    return
    
def freshwaterflux_to_massflux_Gty(src_path, dst_path, rho_fw=1000, year=None, log=True):
    """
    Convert monthly mean subshelf melt time series (m³/s) into annual integrated
    mass fluxes (Gt/yr).

    This function reads the output files generated by
    ``fesom_subshelf_freshwaterflux()``, which contain a variable
    ``subshelf_melt`` representing monthly mean freshwater melt rates in
    m³/s. For each file, the function converts the freshwater flux into an
    ice-mass-equivalent flux (kg/s), multiplies each monthly mean by the exact number
    of seconds in that month (including leap years), and sums over each year
    to obtain a total melt rate in gigatonnes per year (Gt/yr). A new NetCDF
    file with the suffix ``_GTY`` is written for each processed input file.

    Parameters
    ----------
    src_path : str
        Directory containing the input files produced by
        ``fesom_subshelf_melt()``.  
    dst_path : str
        Directory for writing the output files
    rho_fw : float, default 1000
        Density of freshwater in kg/m³.  
        Used to convert freshwater volume flux (m³/s) into mass flux (kg/s).
    year : int, optional
        If specified, only process files for this year. If None, process all files.
    log : bool, default True
        If True, print progress messages when opening, skipping, or saving files.
    """

    print("\n--> Convert freshwaterflux to massflux Gty...")

    if year is not None:
        files2process = np.sort(glob.glob(f"{src_path}subshelf_melt*.{year}.nc"))
    else:
        files2process = np.sort(glob.glob(f"{src_path}subshelf_melt*.nc"))
    
    # Open files with cftime decoder
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
    
    for file in files2process:
        outfile = dst_path + file.split('/')[-1].split('.')[0] + '_GTY.' + file.split('/')[-1].split('.')[1] + '.nc'
        if isfile(outfile):
            if log:
                print(f"Skipping: {file}")
        else:
            if log:
                print(f"Opening: {file}")
            
            ds_subshelf_melt = xr.open_dataset(file, decode_times=time_coder)
            massflux_water = ds_subshelf_melt.subshelf_melt * rho_fw # m3/s * kg/m3 = kg/s 
            massflux_ice = massflux_water.copy() # mass is conserved: ice mass melting == freshwater mass
        
            year = int(massflux_ice.groupby('time.year').mean().year.values)
            seconds = seconds_per_month(year) # compute seconds of each month
            
            Gty = ((massflux_ice * seconds).groupby('time.year').sum() * 1e-12) # convert monthly mean kg/s to Gt/year
        
            ds_out = xr.Dataset(
                    {
                        "subshelf_melt_GTY": Gty
                    },
                    coords={
                        "year": Gty.year
                    }
                    )

            ds_out["subshelf_melt_GTY"].attrs["units"] = "Gt/y"
            ds_out["subshelf_melt_GTY"].attrs["long_name"] = (
                "Freshwater flux inside cavities in Gigatons per year"
            )
    
            ds_out.to_netcdf(outfile)
            if log:
                print(f"Saved: {outfile}")

    return 

def freshwaterflux_to_massflux_Gtm(src_path, dst_path, rho_fw=1000, year=None, log=True):
    """
    Convert monthly mean subshelf melt time series (m³/s) into monthly integrated
    mass fluxes (Gt/month).

    This function reads the output files generated by
    ``fesom_subshelf_freshwaterflux()``, which contain a variable
    ``subshelf_melt`` representing monthly mean freshwater melt rates in
    m³/s. For each file, the function converts the freshwater flux into an
    ice-mass-equivalent flux (kg/s), multiplies each monthly mean by the exact number
    of seconds in that month (including leap years), 
    to obtain a total melt rate in gigatonnes per month (Gt/m). A new NetCDF
    file with the suffix ``_GTM`` is written for each processed input file.

    Parameters
    ----------
    src_path : str
        Directory containing the input files produced by
        ``fesom_subshelf_melt()``.  
    dst_path : str
        Directory for writing the output files
    rho_fw : float, default 1000
        Density of freshwater in kg/m³.  
        Used to convert freshwater volume flux (m³/s) into mass flux (kg/s).
    year : int, optional
        If specified, only process files for this year. If None, process all files.
    log : bool, default True
        If True, print progress messages when opening, skipping, or saving files.
    """

    print("\n--> Convert freshwaterflux to massflux Gtm...")

    if year is not None:
        files2process = np.sort(glob.glob(f"{src_path}subshelf_melt*.{year}.nc"))
    else:
        files2process = np.sort(glob.glob(f"{src_path}subshelf_melt*.nc"))
    
    # Open files with cftime decoder
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
    
    for file in files2process:
        outfile = dst_path + file.split('/')[-1].split('.')[0] + '_GTM.' + file.split('/')[-1].split('.')[1] + '.nc'
        if isfile(outfile):
            if log:
                print(f"Skipping: {file}")
        else:
            if log:
                print(f"Opening: {file}")
            
            ds_subshelf_melt = xr.open_dataset(file, decode_times=time_coder)
            massflux_water = ds_subshelf_melt.subshelf_melt * rho_fw # m3/s * kg/m3 = kg/s 
            massflux_ice = massflux_water.copy() # mass is conserved: ice mass melting == freshwater mass
        
            year = int(massflux_ice.groupby('time.year').mean().year.values)
            seconds = seconds_per_month(year) # compute seconds of each month
            
            Gtm = (massflux_ice * seconds * 1e-12) # convert monthly mean kg/s to Gt/year
        
            ds_out = xr.Dataset(
                    {
                        "subshelf_melt_GTM": Gtm
                    },
                    coords={
                        "year": Gtm.time
                    }
                    )
            ds_out["subshelf_melt_GTM"].attrs["units"] = "Gt/m"
            ds_out["subshelf_melt_GTM"].attrs["long_name"] = (
                "Freshwater flux inside cavities in Gigatons per month"
            )
    
    
            ds_out.to_netcdf(outfile)
            if log:
                print(f"Saved: {outfile}")

    return 

def fesom_subshelf_hydrography(src_path, mesh_diag_path, mesh_path, mask, year, variable='temp', mean=True, log=False, savepath='./'):
    """
    Extract subshelf hydrography fields from FESOM output for a given node mask
    and save the subsetted fields to new NetCDF files.

    Parameters
    ----------
    src_path : str
        Directory containing yearly FESOM variable files, e.g.
        ``/data/fesom/output/``.  
        Each file is expected to follow the format:
        ``{variable}.fesom.{year}.nc``.
    mesh_diag_path : str
        Directory containing the FESOM mesh diagnostic file
        ``fesom.mesh.diag.nc``. This file must include the field
        ``nod_area`` used to add node area information to the output.
    mesh_path : str
        Directory containing the FESOM mesh files. Required by the
        cavity mask builders.
    mask : dict, list, or numpy.ndarray
        Node mask used to subset the FESOM mesh:
        - If a **dict**, it must contain at least a `'name'` key:
            * If ``mask['name'] == 'all'``: use the full cavity mask.
            * Otherwise: use a regional mask defined by a KML file
              provided in ``mask['kml_path']``.
        - If a **list** or **np.array**, it is interpreted as a boolean or
          integer mask directly indexable along the ``nod2`` dimension.
    year : int
        Year to process. Only ``{variable}.fesom.<year>.nc`` is loaded.
    variable : str or list of str, optional
        One or more variable names to extract from the source files.
    mean : bool, optional
        If True, compute the horizontal mean of the variable over the mask.
    log : bool, optional
        If True, print progress messages.
    savepath : str, optional
        Directory where the subsetted NetCDF files will be written.
    """
    
    print("\n--> Compute fesom subshelf hydrography...")

    if isinstance(variable, str):
        variable = [variable]

    mask_name = mask['name'] if isinstance(mask, dict) else "custom"
    todo = []
    for var in variable:
        out_file = f"{savepath}{var}_{mask_name}.{year}.nc"
        if os.path.exists(out_file):
            if log:
                print(f"Skipping existing file: {out_file}", flush=True)
        else:
            todo.append(var)

    if not todo:
        return

    os.makedirs(savepath, exist_ok=True)

    # --- Build mask ---
    if isinstance(mask, dict):
        if mask['name'] == 'all':
            node_mask = build_cavity_mask(mesh_path, which='node')
        else:
            node_mask = build_cavity_regional_mask(
                mesh_path, mask['kml_path'],
                name=mask['name'], which='node'
            )
    elif isinstance(mask, (list, np.ndarray)):
        node_mask = mask
    else:
        raise ValueError('Mask type is not supported!')

    # --- Load mesh diag once ---
    mesh_diag = xr.open_dataset(f"{mesh_diag_path}fesom.mesh.diag.nc")

    # Precompute nodal area (lazy, small)
    nod_area = mesh_diag.nod_area.max(dim='nz').isel(nod2=node_mask)

    # Use cftime decoder
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)

    for var in todo:

        if log:
            print(f"Opening file for {var}, year {year}...")

        # OPEN ONLY THE REQUESTED YEAR
        ds = xr.open_mfdataset(
            [f"{src_path}{var}.fesom.{year}.nc"],
            combine='by_coords',
            parallel=True,
            decode_times=time_coder,
            chunks={'time': 12, 'nod2': 50000}
        )

        # --- subset nodes early ---
        ds = ds.isel(nod2=node_mask)

        # --- add weights ---
        ds['nod_area'] = nod_area

        # --- compute ---
        if mean:
            result = ds[var].weighted(ds['nod_area']).mean(dim='nod2')
        else:
            result = ds[var]

        # IMPORTANT: force graph creation ONCE (optional but often faster downstream)
        result = result.chunk({'time': 12})

        # --- write ---
        outfile = f"{savepath}{var}_{mask_name}.{year}.nc"

        if log:
            print(f"Writing {outfile}...")

        result.sel(time=str(year)).to_netcdf(
            outfile,
            #engine='h5netcdf',
            #encoding={var: {'zlib': True, 'complevel': 4}}
        )

        if log:
            print(f"Saved: {outfile}")

    return
