# so_ase.eval_ocean.py

import xarray as xr
import numpy as np
import pyfesom2 as pf
from scipy.stats import linregress
from scipy.interpolate import griddata
from os.path import isfile
from os import remove
from .helpers_mesh import find_nodes_in_box, add_element_volumes, build_cavity_mask

# To Do:
# ocean temperature 50S-65S horizontal mean


def fesom_ocean_heat_transport_as_residual(
    src_path,
    mesh_diag_path,
    ref_date="2000-01-31",
    eval_date="2002-01-31",
    box=[-180, 180, -90, -60],
    rho=1028,
    cp=4190,
    log=True,
):
    """
    Compute the ocean heat transport (OHT) into a region of interest during a particular period as the residual of the total surface heat flux during the period and ocean heat
    content change between the first and last timestep of the period, using FESOM2 model output.

    Parameters
    ----------
    src_path : str
        Directory path to the FESOM2 NetCDF output files (e.g., `fh.fesom.YYYY.nc`, `temp.fesom.YYYY.nc`).
    mesh_diag_path : str
        Directory path to the mesh diagnostic NetCDF file `fesom.mesh.diag.nc`.
    ref_date : str, optional
        Start date of the analysis period in 'YYYY-MM-DD' format. Default is '2000-01-31'.
    eval_date : str, optional
        End date of the analysis period in 'YYYY-MM-DD' format. Default is '2002-01-31'.
    box : list of float, optional
        Geographic bounding box `[lon_min, lon_max, lat_min, lat_max]` to subset the spatial region of interest.
        Default is `[-180, 180, -90, -60]` (Southern Ocean).
    rho : float, optional
        Seawater density in kg/m³. Default is 1028.
    cp : float, optional
        Seawater specific heat capacity in J/(kg·K). Default is 4190.
    log : bool, optional
        If True, prints logging information during processing. Default is True.

    Returns
    -------
    float
        The estimated ocean heat transport (OHT) in joules (J) over the specified period and region.

    Notes
    -----
    - The OHT is estimated as the residual of:

      where:
        - ΔOHC is the change in ocean heat content between `ref_date` and `eval_date`
        - SHF is the surface heat flux (in W/m²)
        - The surface heat flux integral includes time (in seconds) and nodal area (m²)

    - This function assumes monthly-averaged data and uses fixed calendar days per month (non-leap year).
    """
    # Preprocess inputs
    years = (int(ref_date.split("-")[0]), int(eval_date.split("-")[0]))

    # Load mesh dignostic file
    mesh_diag = xr.open_dataset(f"{mesh_diag_path}fesom.mesh.diag.nc")
    if log:
        print("Mesh diagnostics loaded:", flush=True)
        print(f"{mesh_diag_path}fesom.mesh.diag.nc", flush=True)

    # Load files for surface heat flux from src_path
    files2load = [f"{src_path}fh.fesom.{y}.nc" for y in range(years[0], years[1])]
    # Open files with cftime decoder
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
    ds_shf = xr.open_mfdataset(files2load, decode_times=time_coder).load()
    if log:
        print("Files loaded:", flush=True)
        for f in files2load:
            print(f"{f}", flush=True)

    # Crop to examination period
    ds_shf = ds_shf.sel(time=slice(ref_date, eval_date))
    if log:
        print(f"Cropped to: {ref_date} to {eval_date}", flush=True)

    # Load files for ocean temperature from src_path (only for the ref and eval year)
    files2load = [f"{src_path}temp.fesom.{y}.nc" for y in [years[0], years[1]]]
    ds_temp = xr.open_mfdataset(files2load, decode_times=time_coder).load()
    if log:
        print("Files loaded:", flush=True)
        for f in files2load:
            print(f"{f}", flush=True)

    # Crop to examination period
    ds_temp = ds_temp.sel(time=[ref_date, eval_date], method="nearest")
    if log:
        print(f"Cropped to: {ref_date} / {eval_date}", flush=True)
        print(f"Timestamps: {ds_temp.time[0].values}, {ds_temp.time[1].values}")

    # Add layerthickness and volumes to mesh_diag
    mesh_diag["layer_thickness"] = (("nz1"), np.diff(mesh_diag.nz))
    mesh_diag["volumes"] = mesh_diag.nod_area.isel(nz=0) * mesh_diag.layer_thickness

    # Find indices of nodes within the specified box
    inds = find_nodes_in_box(mesh_diag_path, box=box, log=log)

    # Crop datasets to specified box
    ds_temp_cropped = ds_temp.isel(nod2=inds)
    ds_shf_cropped = ds_shf.isel(nod2=inds)
    mesh_diag_cropped = mesh_diag.isel(nod2=inds)

    # Compute ocean heat content change (OHC(eval_date) - OHC(ref_date))
    OHC = (ds_temp_cropped.temp * rho * cp * mesh_diag_cropped.volumes).sum(
        dim=("nod2", "nz1")
    )
    deltaOHC = OHC.isel(time=-1) - OHC.isel(time=0)

    # Compute total number of seconds per month for total SHF computation
    days_in_month = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
    month_indices = (
        ds_shf_cropped.time.dt.month.values - 1
    )  # values already numpy array
    seconds_per_month = days_in_month[month_indices] * 86400
    ds_shf_cropped["seconds_per_month"] = (("time"), seconds_per_month)

    # Compute total sum of the surface heat flux over the full examination period (SHF [W/m2] * seconds_per_month [s] * nodal area [m2])
    SHF = (
        ds_shf_cropped.fh
        * ds_shf_cropped.seconds_per_month
        * mesh_diag_cropped.nod_area.isel(nz=0)
    ).sum(dim="nod2")
    sumSHF = SHF.sum(dim="time")

    # Compute the ocean heat transport as a residual (deltaOHC - sumSHF)
    OHT = deltaOHC - sumSHF

    return OHT.values

def fesom_timeseries_of_mean_vertical_profile_in_region(
    src_path,
    mesh_diag_path,
    years=(1979, 2015),
    box=[-180, 180, -90, -60],
    varname='temp',
    log=True
):
    """
    Compute a time series of area-weighted mean vertical profiles for a specified region 
    from FESOM2 model output.

    Parameters:
    -----------
    src_path : str
        Path to the directory containing annual FESOM2 output files (e.g., temp.fesom.YYYY.nc).
    
    mesh_diag_path : str
        Path to the directory containing the mesh diagnostic file (fesom.mesh.diag.nc).
    
    years : tuple of int, optional
        Start and end year for the time series. Default is (1979, 2015).
    
    box : list of float, optional
        Geographic bounds of the region of interest in the format [lon_min, lon_max, lat_min, lat_max].
        Default is global Southern Ocean: [-180, 180, -90, -60].
    
    varname : str, optional
        Name of the variable to extract and average (must match variable name in NetCDF files).
        Default is 'temp'.
    
    log : bool, optional
        Whether to print progress information. Default is True.

    Returns:
    --------
    xarray.DataSet
        A time series of area-weighted mean vertical profiles in the specified region.
        The vertical levels are preserved along the 'nz' dimension.
    """
    
    # Find indices of nodes within the specified box
    if log:
        print(f"Box {box[0]}E {box[1]}E {box[2]}N {box[3]}N")
    inds = find_nodes_in_box(mesh_diag_path, box=box, log=log)

    mesh_diag = xr.open_dataset(f"{mesh_diag_path}fesom.mesh.diag.nc")
    nodal_area = mesh_diag.nod_area.isel(nod2=inds, nz=0)
    
    mean_profiles = []

    # Open files with cftime decoder
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)

    for year in range(years[0], years[-1]):
        # Load file for each single year
        file2load = f"{src_path}{varname}.fesom.{year}.nc"
        ds = xr.open_mfdataset(file2load, decode_times=time_coder).isel(nod2=inds).load()
        if log:
            print(f"File loaded: {file2load}", flush=True)
        
        # Compute the area-weighted horizontal mean of all vertical profiles 
        mean_profiles.append(ds.weighted(nodal_area).mean(dim='nod2'))
        
    # Concatenate to one Dataset
    ds_out = xr.concat(mean_profiles, dim='time')

    return ds_out
    
def fesom_total_kinetic_energy(src_path, mesh_diag_path, meshpath, years=(1979, 2015), mask='cavity', log=False, savepath='./'):  

    # load mesh diagnostics
    mesh_diag = xr.open_dataset(f"{mesh_diag_path}fesom.mesh.diag.nc")
    mesh_diag = add_element_volumes(mesh_diag)

    # Build mask
    if mask == 'all':
        element_mask = np.ones_like(len(mesh_diag.elem_area))
    elif mask == 'cavity':
        element_mask = build_cavity_mask(meshpath, which='element')
    elif mask == 'open_ocean':
        mask = build_cavity_mask(meshpath, which='element')
        element_mask = ~mask
    else:
        pass
        
    # Build list of all input files
    years_list = list(range(years[0], years[-1]))
    files_u = [f"{src_path}u.fesom.{y}.nc" for y in years_list]
    files_v = [f"{src_path}v.fesom.{y}.nc" for y in years_list]

    for i, (file_u, file_v) in enumerate(zip(files_u, files_v)):

        file2save = f"{savepath}kinetic_energy_{mask}.{years_list[i]}.nc"
        if not isfile(file2save):
            
            # Open files with cftime decoder
            time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
            
            ds_u = xr.open_dataset(file_u, decode_times=time_coder)
            ds_v = xr.open_dataset(file_v, decode_times=time_coder)
        
            if log:
                print(f"Opened:{file_u}, {file_v}", flush=True)
        
            # ---- Compute |u|^2 from u, v ----
            da_uv2 = ds_u.u**2 + ds_v.v**2
        
            # ---- Multiply by constant element volumes ----
            KE = 0.5 * mesh_diag.elem_volume * 1030 * da_uv2 # 1/2 * m * |u|^2
        
            # ---- Sum along vertical ----
            KE = KE.sum(dim="nz1")
        
            # ---- Apply element mask ----
            if mask is not None:
                KE = KE.isel(elem=element_mask).sum(dim="elem")
            else:
                KE = KE.sum(dim="elem")
        
            ds_out = xr.Dataset(
            {
                "KE": KE
            },
            coords={
                "time": KE.time
            }
            )
    
            ds_out.to_netcdf(file2save)
            if log:
                print(f"Saved: {file2save}")
        else:
            if log:
                print(f"Skipped: {file2save}")

    return
        

def fesom_total_runoff(src_path, meshpath, basin_mask_file, basins=[66], years=(1979, 2015), which='solid', ori=False, log=False, savepath='./', replace=False):
    """
    Calculate total runoff for specified basins from FESOM output files.
    
    This function processes FESOM runoff data files, applies basin masks to isolate
    runoff from specific drainage basins, and saves the total runoff time series
    for each year as separate NetCDF files.
    
    Parameters:
        src_path (str): Path to directory containing FESOM runoff files.
        meshpath (str): Path to directory containing FESOM mesh files.
        basin_mask_file (str): Path to NetCDF file containing basin mask data.
        basins (list or int): List of basin IDs to include in calculation. Defaults to [66].
        years (tuple): Tuple of (start_year, end_year) for processing period. Defaults to (1979, 2015).
        which (str): Type of runoff to process ('solid' or 'liquid'). Defaults to 'solid'.
        ori (bool): If True, process original runoff files with '_ori' suffix. Defaults to False.
        log (bool): If True, print progress messages. Defaults to False.
        savepath (str): Directory path for output files. Defaults to './'.
        replace (bool): If True, overwrite existing output files. Defaults to False.
    
    Returns:
        None
    
    Notes:
        - Input files follow pattern: runoff_{which}.fesom.{year}.nc or runoff_{which}_ori.fesom.{year}.nc
        - Output files follow pattern: {which}_runoff_basins_{basin_ids}.{year}.nc
        - Runoff is calculated as sum over nodes in specified basins, weighted by nodal area
        - Uses mesh diagnostics from fesom.mesh.diag.nc for nodal area information
    """
    # Load mesh diagnostics
    mesh_diag = xr.open_dataset(f"{meshpath}fesom.mesh.diag.nc")

    # Build basin mask
    if log:
        print(f"Read nodal basin mask from: {basin_mask_file}")
    
    if not isinstance(basins, list):
        basins = [basins]
        
    basin_mask = xr.open_dataset(basin_mask_file).basin_id
    node_mask = np.isin(basin_mask, basins)

    # Build list of all input files
    years_list = list(range(years[0], years[-1]))
    if ori:
        files = [f"{src_path}runoff_{which}_ori.fesom.{y}.nc" for y in years_list]
    else:
        files = [f"{src_path}runoff_{which}.fesom.{y}.nc" for y in years_list]

    for i, file in enumerate(files):

        if ori:
            file2save = f"{savepath}{which}_runoff_ori_basins_{"_".join(map(str, basins))}.{years_list[i]}.nc"
        else:
            file2save = f"{savepath}{which}_runoff_basins_{"_".join(map(str, basins))}.{years_list[i]}.nc"
            
        if not isfile(file2save):
            time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
            ds_rnf = xr.open_dataset(file, decode_times=time_coder)
            RNF = (ds_rnf[f"runoff_{which}"] * mesh_diag.nod_area.isel(nz=0)).isel(nod2=node_mask).sum(dim='nod2')

            ds_out = xr.Dataset(
            {
                f"runoff_{which}": RNF
            },
            coords={
                "time": RNF.time
            }
            )
    
            ds_out.to_netcdf(file2save)
            if log:
                print(f"Saved: {file2save}")
        else:
            if replace:
                remove(file2save)
                if log:
                    print(f"Removed existing: {file2save}")
                    
                time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
                ds_rnf = xr.open_dataset(file, decode_times=time_coder)
                RNF = (ds_rnf[f"runoff_{which}"] * mesh_diag.nod_area.isel(nz=0)).isel(nod2=node_mask).sum(dim='nod2')
    
                ds_out = xr.Dataset(
                {
                    f"runoff_{which}": RNF
                },
                coords={
                    "time": RNF.time
                }
                )
        
                ds_out.to_netcdf(file2save)
                if log:
                    print(f"Saved: {file2save}")
            else:
                if log:
                    print(f"Skipped: {file2save}")
    return
    
        

        

