# so_ase.eval_cavity.py

import xarray as xr
import numpy as np
import glob
from os.path import isfile
from .helpers_mesh import build_cavity_mask
from .helpers_misc import seconds_per_month

def fesom_subshelf_freshwaterflux(src_path, mesh_diag_path, mesh_path, years=(1979, 2015), mask='all', log=False, savepath='./'):
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
    years : tuple of int, default (1979, 2015)
        Start and end years. The function processes files for all years in 
        ``range(years[0], years[-1])``. (Note: the end year is excluded.)
    mask : str or array-like, default 'all'
        Identifier for which ocean cavity mask to apply.  
        If ``'all'`` is given, the mask is generated using 
        ``build_cavity_mask(mesh_path, which='node')``.  
        If an array-like mask is supplied, it is used directly for ``isel``.
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
    if mask == 'all':
        node_mask = build_cavity_mask(mesh_path, which='node')
        
    # Build list of all input files
    years_list = list(range(years[0], years[-1]))
    files = [f"{src_path}fw.fesom.{y}.nc" for y in years_list]

    for i, file in enumerate(files):

        file2save = f"{savepath}subshelf_melt_{mask}.{years_list[i]}.nc"
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

def fesom_subshelf_heatflux(src_path, mesh_diag_path, mesh_path, years=(1979, 2015), mask='all', log=False, savepath='./'):

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
    years : tuple of int, default (1979, 2015)
        Start and end years. The function processes files for all years in
        ``range(years[0], years[-1])`` (end year excluded).
    mask : str or array-like, default 'all'
        Mask specifying which nodes to include.  
        If ``'all'`` is provided, a node-based cavity mask is generated using  
        ``build_cavity_mask(mesh_path, which='node')``.  
        If an array-like mask is supplied, it is applied directly.
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
    if mask == 'all':
        node_mask = build_cavity_mask(mesh_path, which='node')
        
    # Build list of all input files
    years_list = list(range(years[0], years[-1]))
    files = [f"{src_path}fh.fesom.{y}.nc" for y in years_list]

    for i, file in enumerate(files):

        file2save = f"{savepath}subshelf_heatflux_{mask}.{years_list[i]}.nc"
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
    
def freshwaterflux_to_massflux_Gty(src_path, rho_fw=1000, filename='subshelf_melt_all.*.nc', log=True):
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
        Expected filenames match the pattern given in ``filename``.
    rho_fw : float, default 1000
        Density of freshwater in kg/m³.  
        Used to convert freshwater volume flux (m³/s) into mass flux (kg/s).
    filename : str, default 'subshelf_melt_all.*.nc'
        Glob pattern indicating which NetCDF files to process.
        Each file must contain a DataArray named ``subshelf_melt``.
    log : bool, default True
        If True, print progress messages when opening, skipping, or saving files.
    """

    files2process = np.sort(glob.glob(f"{src_path}{filename}"))
    
    for file in files2process:
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

        file2save = src_path + filename.split('.')[0] +'_GTY' + '.' + str(year) + '.nc'

        if isfile(file2save):
            if log:
                print(f"Skipped: {file2save}")
        else:
            ds_out.to_netcdf(file2save)
            if log:
                print(f"Saved: {file2save}")
    
    return 
