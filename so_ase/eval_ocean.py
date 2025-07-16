# so_ase.eval_ocean.py

import xarray as xr
import numpy as np
import pyfesom2 as pf
from .helpers_mesh import find_nodes_in_box

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
    ds_shf = xr.open_mfdataset(files2load).load()
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
    ds_temp = xr.open_mfdataset(files2load).load()
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
    seconds_per_month = days_in_month[month_indices] * 86500
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

    for year in years:
        # Load file for each single year
        file2load = f"{src_path}{varname}.fesom.{year}.nc"
        ds = xr.open_dataset(file2load).isel(nod2=inds).load()
        if log:
            print(f"File loaded for year {year}:", flush=True)
        
        # Compute the area-weighted horizontal mean of all vertical profiles 
        mean_profiles.append(ds.weighted(nodal_area).mean(dim='nod2'))
    
    n_samples_per_depth = (ds[varname].isel(time=0) != 0).sum(dim='nod2')
    # Concatenate to one Dataset
    ds_out = xr.concat(mean_profiles, dim='time')

    return ds_out, n_samples_per_depth
        

        
        
        
        
        

        

        

