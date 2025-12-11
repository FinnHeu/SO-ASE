# so_ase.eval_cavity.py

import xarray as xr
import numpy as np
import glob
from os.path import isfile
from .helpers_mesh import build_cavity_mask, build_cavity_regional_mask
from .helpers_misc import seconds_per_month

def fesom_subshelf_freshwaterflux(src_path, mesh_diag_path, mesh_path, mask, years=(1979, 2015), log=False, savepath='./'):
    """
    Compute and save integrated subshelf basal melt time series from FESOM 
    freshwater flux output [m3/s].

    This function loads monthly freshwater flux fields (`fw`) from a sequence 
    of FESOM output files, multiplies them by the mean horizontal node area 
    from a mesh diagnostics file, applies a cavity mask, and integrates the 
    melt signal over all masked nodes to obtain a single subshelf melt 
    time series for each year. The result for each year is saved as a 
    standalone NetCDF file.

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
    years : tuple of int, default (1979, 2015)
        Start and end years. The function processes files for all years in 
        ``range(years[0], years[-1])``. (Note: the end year is excluded.)
    log : bool, default False
        If True, print a message whenever a file is saved or skipped.
    savepath : str, default './'
        Directory in which to save the output files. The files are created as
        ``subshelf_melt_<mask>.<year>.nc``.

    # load mesh diagnostics
    mesh_diag = xr.open_dataset(f"{mesh_diag_path}fesom.mesh.diag.nc")
    """
    # load mesh diag
    mesh_diag = xr.open_dataset(f"{mesh_diag_path}fesom.mesh.diag.nc")

    # Build mask
    if isinstance(mask, dict):
        if mask['name'] == 'all':
            node_mask = build_cavity_mask(mesh_path, which='node')
        else:
            node_mask = build_cavity_regional_mask(mesh_path, mask['kml_path'], which=mask['name'])
    elif isinstance(mask, list) or isinstance(mask, np.array):
        node_mask = mask
    else:
        raise ValueError('Mask type is not supported!')
        
    # Build list of all input files
    years_list = list(range(years[0], years[-1]))
    files = [f"{src_path}fw.fesom.{y}.nc" for y in years_list]

    for i, file in enumerate(files):

        file2save = f"{savepath}subshelf_melt_{mask['name']}.{years_list[i]}.nc"
        if not isfile(file2save):
            
            ds_fw = xr.open_dataset(file)
            SSM = (ds_fw.fw * mesh_diag.nod_area.mean(dim='nz')).isel(nod2=node_mask).sum(dim='nod2')

            ds_out = xr.Dataset(
            {
                "subshelf_melt": SSM
            },
            coords={
                "time": SSM.time
            }
            )
    
            ds_out.to_netcdf(file2save)
            if log:
                print(f"Saved: {file2save}")
        else:
            if log:
                print(f"Skipped: {file2save}")

    return

def fesom_subshelf_heatflux(src_path, mesh_diag_path, mesh_path, mask, years=(1979, 2015), log=False, savepath='./'):

    """
    Compute and save integrated subshelf heat flux time series from FESOM 
    heat flux output.

    This function loads monthly subshelf heat flux fields (`fh`) from a 
    sequence of annual FESOM output files, multiplies them by the mean 
    horizontal node area from a mesh diagnostics file, applies a cavity mask, 
    and integrates the heat flux over all masked nodes to obtain a single 
    area-integrated time series for each year. The result for each year is 
    saved as a standalone NetCDF file.

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
    years : tuple of int, default (1979, 2015)
        Start and end years. The function processes files for all years in 
        ``range(years[0], years[-1])``. (Note: the end year is excluded.)
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
    # load mesh diagnostics
    mesh_diag = xr.open_dataset(f"{mesh_diag_path}fesom.mesh.diag.nc")

    # Build mask
    if isinstance(mask, dict):
        if mask['name'] == 'all':
            node_mask = build_cavity_mask(mesh_path, which='node')
        else:
            node_mask = build_cavity_regional_mask(mesh_path, mask['kml_path'], which=mask['name'])
    elif isinstance(mask, list) or isinstance(mask, np.array):
        node_mask = mask
    else:
        raise ValueError('Mask type is not supported!')
        
    # Build list of all input files
    years_list = list(range(years[0], years[-1]))
    files = [f"{src_path}fh.fesom.{y}.nc" for y in years_list]

    for i, file in enumerate(files):

        file2save = f"{savepath}subshelf_melt_{mask['name']}.{years_list[i]}.nc"
        if not isfile(file2save):
            
            ds_fh = xr.open_dataset(file)
            SSHF = (ds_fh.fh * mesh_diag.nod_area.mean(dim='nz')).isel(nod2=node_mask).sum(dim='nod2')

            ds_out = xr.Dataset(
            {
                "subshelf_heatflux": SSHF
            },
            coords={
                "time": SSHF.time
            }
            )
    
            ds_out.to_netcdf(file2save)
            if log:
                print(f"Saved: {file2save}")
        else:
            if log:
                print(f"Skipped: {file2save}")

    return
    
def freshwaterflux_to_massflux_Gty(src_path, dst_path, rho_fw=1000, log=True):
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
    log : bool, default True
        If True, print progress messages when opening, skipping, or saving files.
    """

    files2process = np.sort(glob.glob(f"{src_path}*.nc"))
    
    for file in files2process:
        outfile = dst_path + file.split('/')[-1].split('.')[0] + '_GTY.' + file.split('/')[-1].split('.')[1] + '.nc'
        if isfile(outfile):
            if log:
                print(f"Skipping: {file}")
        else:
            if log:
                print(f"Opening: {file}")
            
            ds_subshelf_melt = xr.open_dataset(file)
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
    
    
            ds_out.to_netcdf(outfile)
            if log:
                print(f"Saved: {outfile}")

    return 

def fesom_subshelf_hydrography(src_path, mesh_diag_path, mesh_path, mask, years=(1979,1980), variable='temp', log=False, savepath='./'):
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
    years : tuple of int, optional
        Two-element tuple defining the year range to process.
        The function processes all years in ``range(years[0], years[1])``.
        For example, ``years=(1979, 1980)`` processes only 1979.
    variable : str or list of str, optional
        One or more variable names to extract from the source files.
    log : bool, optional
        If True, print progress messages.
    savepath : str, optional
        Directory where the subsetted NetCDF files will be written.
    """
    # Build mask
    if log:
        print("Building mask...")
    if isinstance(mask, dict):
        if mask['name'] == 'all':
            node_mask = build_cavity_mask(mesh_path, which='node')
        else:
            node_mask = build_cavity_regional_mask(mesh_path, mask['kml_path'], which=mask['name'])
    elif isinstance(mask, list) or isinstance(mask, np.array):
        node_mask = mask
    else:
        raise ValueError('Mask type is not supported!')

    if log:
        print("Loading mesh diag...")
    mesh_diag = xr.open_dataset(f"{mesh_diag_path}fesom.mesh.diag.nc")
    years_list = list(range(years[0], years[-1]))

    if isinstance(variable, str):
        variable = [variable]    
    
    for var in variable:
        for y in years_list:
            infile = f"{src_path}{var}.fesom.{y}.nc"
            outfile = f"{savepath}{var}_{mask['name']}_{y}.nc"
            if not isfile(outfile):
                ds = xr.open_dataset(infile).load().isel(nod2=node_mask)
                ds['nod_area'] = mesh_diag.nod_area.isel(nod2=node_mask, nz=0).squeeze()
                ds.to_netcdf(outfile)
                if log:
                    print(f"Saved: {outfile}")
            else:
                if log:
                    print(f"Skipped: {outfile}")

    return