# so_ase.eval_calving.py

import xarray as xr
from so_ase.helpers_misc import seconds_per_month
from .helpers_mesh import build_runoff_basin_mask

def fesom_calving_flux(src_path, meshpath, mesh_diag_path, runoff_maps, basin=66, which='runoff_solid', years=(1979, 2015), log=False, savepath='./'):
    """
    Calculate calving flux from FESOM runoff data.
    
    Parameters:
    -----------
    src_path : str
        Path to the directory containing the runoff files.
    meshpath : str
        Path to the mesh directory.
    mesh_diag_path : str
        Path to the mesh diagnostic file.
    runoff_maps : str
        Path to the runoff maps file.
    basin : int, optional
        Basin ID to filter the runoff data. Default is 66.
    which : str, optional
        Which runoff type to use. Default is 'runoff_solid'.
    years : tuple, optional
        Tuple of start and end years to process. Default is (1979, 2015).
    log : bool, optional
        Whether to log the processing steps. Default is False.
    savepath : str, optional
        Path to save the output files. Default is './'.
    
    Returns:
    --------
    None
    """

    mesh_diag = xr.open_dataset(f"{mesh_diag_path}/fesom.mesh.diag.nc")
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)

    if which == 'runoff_solid':
        # load files
        files = [f"{src_path}/runoff_solid.fesom.{y}.nc" for y in range(years[0], years[1])]
        ds_calving = xr.open_mfdataset(files, combine='by_coords', decode_times=time_coder)

        # build mask
        ds_mask = build_runoff_basin_mask(meshpath, runoff_maps, which='solid')
        calving = (ds_calving["runoff_solid"].where(ds_mask.basin_id == basin, 0) * mesh_diag["nod_area"].isel(nz=0)).sum(dim='nod2')# m3/s montly mean

    # convert to Gt/year
    sec = seconds_per_month([y for y in range(years[0], years[1])])
    calving = calving * sec * 1e-9 # m3/s monthly mean to Gt per month

    # Split the monthly data into chunks for each individual year and save as xarray dataset
    for y in range(years[0], years[1]):
        calving_year = calving.sel(time=calving.time.dt.year == y)
        calving_year = calving_year.to_dataset(name='calving_GTM')
        calving_year['calving_GTM'].attrs['units'] = 'Gt per month'
        calving_year.to_netcdf(f"{savepath}/calving_basin_{basin}_GTM_{y}.nc")
    
    return

        
    
    
    
