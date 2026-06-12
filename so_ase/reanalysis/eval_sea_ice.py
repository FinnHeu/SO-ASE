# so_ase/reanalysis/eval_sea_ice.py

import xarray as xr
import numpy as np
import os
import glob

from .helpers_grid import *

def nsidc_ice_diag(src_path,
    years=(1979, 2015),
    box=[-180, 180, -90, -50],
    siconc_threshold=0.15,
    grouping='annual.mean',
    version=5,
    diag='area',               
    log=True):

    """
    Compute total sea ice area from NSIDC CDR data within a geographic region.

    This function loads NSIDC sea ice concentration data, reprojects it from 
    polar stereographic coordinates to regular lat/lon, and calculates the 
    total sea ice area where concentration exceeds a given threshold.

    Parameters
    ----------
    src_path : str
        Path to directory containing SIC NetCDF files named as 'sic.<year>.nc'.
    years : tuple of int, optional
        Start and end year (exclusive) for processing, e.g., (1979, 2015).
    box : list of float, optional
        Geographic bounding box [lon_min, lon_max, lat_min, lat_max] for area calculation.
    siconc_threshold : float, optional
        Sea ice concentration threshold above which a grid cell is counted as ice-covered (default is 0.15).
    grouping : {'annual.mean', 'annual.max', 'annual.min', 'monthly.mean'}, optional
        Temporal aggregation applied to the computed sea-ice diagnostic. The default is 'annual.mean'.
    version : {5, 6}, optional
        NSIDC product version. Determines the variable name used during calculation.
    diag : {'area', 'extent'}, optional
        Which diagnostic to compute after filtering all cells exceeding the siconc_threshold: 'area' (concentration × cell area) or 'extent'(area of ALL cells exceeding siconc_threshold)
    log : bool, optional
        If True, print progress messages during processing.

    Returns
    -------
    xarray.DataArray
        Time series of total sea ice area (in km²) per time step within the specified region.

    Notes
    -----
    - Assumes a constant 25 km x 25 km grid cell size.

    """

    files2load = [f"{src_path}siconc.{y}.nc" for y in range(years[0], years[1])]

    if version==6:
        var = 'cdr_seaice_conc_monthly'
    elif version==5:
        var= 'siconc'

    # Open files with cftime decoder
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)

    result = []
    for file in files2load:    
    
        ds = xr.open_dataset(file, decode_times=time_coder).load()
        if log:
            print(f"File loaded: {file}", flush=True)

        # transform projection
        ds = reproject_to_latlon(ds)

        # Compute sea ice area
        mask_geo = (
            (ds.lat >= box[2]) & (ds.lat <= box[3]) &
            (ds.lon >= box[0]) & (ds.lon <= box[1])
        )

        mask_isice = ds[var] > siconc_threshold
    
        gridcell_area = 25000 * 25000 # m^2

        if diag=='area':
            ice_area = ((ds[var] * mask_geo * mask_isice) * gridcell_area).sum(dim=('x','y'))
            name = 'sea_ice_area'
        elif diag=='extent':
            ice_area = ((mask_geo * mask_isice) * gridcell_area).sum(dim=('x','y'))
            name = 'sea_ice_extent'
        
        ice_area = ice_area.where(ice_area>0,np.nan)

        result.append(ice_area)

    result = xr.concat(result, dim='time').rename(name)

    result.attrs['units'] = 'm^2'
    result.attrs['long_name'] = f'total {name}'
    result.attrs['bounding box'] = (
        f"Longitude: {box[0]}E to {box[1]}E, Latitude: {box[2]}N to {box[3]}N"
    )

    if grouping == 'annual.mean':
        result = result.groupby('time.year').mean('time')
    elif grouping == 'annual.max':
        result = result.groupby('time.year').max('time')
    elif grouping == 'annual.min':
        result = result.groupby('time.year').min('time')
    elif grouping == 'monthly.mean':
        pass
    
    print('Done!')
    return result
    
def hadlsst_ice_area(src_path,
    years=(1979, 2015),
    box=[-180, 180, -90, -50],
    siconc_threshold=0.15,
    grouping='annual.mean',
    log=True):

    """
    Compute total sea ice area from HadlSST_ice data within a geographic region.

    This function loads HadlSST sea ice concentration data, adds gridd cell areea, and calculates the 
    total sea ice area where concentration exceeds a given threshold.

    Parameters
    ----------
    src_path : str
        Path to directory containing SIC NetCDF files named as 'sic.<year>.nc'.
    years : tuple of int, optional
        Start and end year (exclusive) for processing, e.g., (1979, 2015).
    box : list of float, optional
        Geographic bounding box [lon_min, lon_max, lat_min, lat_max] for area calculation.
    siconc_threshold : float, optional
        Sea ice concentration threshold above which a grid cell is counted as ice-covered (default is 0.15).
    grouping : {'annual.mean', 'annual.max', 'annual.min', 'monthly.mean'}, optional
        Temporal aggregation applied to the computed sea-ice diagnostic. The default is 'annual.mean'.
    log : bool, optional
        If True, print progress messages during processing.

    Returns
    -------
    xarray.DataArray
        Time series of total sea ice area (in km²) per time step within the specified region.

    Notes
    -----
    - Compute grid cell size from lat/lon
    """

    files2load = [f"{src_path}siconc.{y}.nc" for y in range(years[0], years[1])]

    # Open files with cftime decoder
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)

    result = []
    for file in files2load:    
    
        ds = xr.open_dataset(file, decode_times=time_coder)
        if log:
            print(f"File loaded: {file}", flush=True)

        # add gridcell area
        ds = gridcell_area_hadley(ds, R=6371.0)

        # Crop to box
        ds = ds.sel(latitude=slice(box[2], box[3]), longitude=slice(box[0], box[1]))

        # commented line is ice extent, we want ice area
        #ice_area = ds.cell_area.where(ds.siconc > siconc_threshold, 0).sum(dim=('longitude','latitude'))

        mask_isice = ds.siconc > siconc_threshold
        ice_area = (ds.siconc * mask_isice * ds.cell_area).sum(dim=('longitude','latitude'))

        result.append(ice_area)

    result = xr.concat(result, dim='time').rename('sea_ice_area')

    result.attrs['units'] = 'm^2'
    result.attrs['long_name'] = 'total sea ice area'
    result.attrs['bounding box'] = (
        f"Longitude: {box[0]}E to {box[1]}E, Latitude: {box[2]}N to {box[3]}N"
    )

    if grouping == 'annual.mean':
        result = result.groupby('time.year').mean('time')
    elif grouping == 'annual.max':
        result = result.groupby('time.year').max('time')
    elif grouping == 'annual.min':
        result = result.groupby('time.year').min('time')
    elif grouping == 'monthly.mean':
        pass
    
    print('Done!')
    return result

def osisaf_ice_diag(src_path,
    years=(1979, 2015),
    box=[-180, 180, -90, -50],
    siconc_threshold=0.15,
    grouping='annual.mean',
    diag='area',               
    log=True):

    """
    Compute total sea ice area from OSISAF_v3p0 data within a geographic region.

    This function loads OSISAF_v3p0 sea ice concentration data, reprojects it from 
    polar stereographic coordinates to regular lat/lon, and calculates the 
    total sea ice area where concentration exceeds a given threshold.

    Parameters
    ----------
    src_path : str
        Path to directory containing SIC NetCDF files named, e.g. 'ice_conc_sh_ease2-250_cdr-v3p0_199504.nc'.
    years : tuple of int, optional
        Start and end year (exclusive) for processing, e.g., (1979, 2015).
    box : list of float, optional
        Geographic bounding box [lon_min, lon_max, lat_min, lat_max] for area calculation.
    siconc_threshold : float, optional
        Sea ice concentration threshold above which a grid cell is counted as ice-covered (default is 0.15).
    grouping : {'annual.mean', 'annual.max', 'annual.min', 'monthly.mean'}, optional
        Temporal aggregation applied to the computed sea-ice diagnostic. The default is 'annual.mean'.
    diag : {'area', 'extent'}, optional
        Which diagnostic to compute after filtering all cells exceeding the siconc_threshold: 'area' (concentration × cell area) or 'extent'(area of ALL cells exceeding siconc_threshold)
    log : bool, optional
        If True, print progress messages during processing.

    Returns
    -------
    xarray.DataArray
        Time series of total sea ice area (in km²) per time step within the specified region.

    Notes
    -----
    - Assumes a constant 25 km x 25 km grid cell size.
    - OSISAF SIC is stored as percent values (0–100) and is internally converted to fractional concentration (0–1) for area calculations.

    """
    
    
    files2load = []

    for y in range(years[0], years[1]):
        files2load.extend(glob.glob(f"{src_path}/ice_conc_sh*{y}*.nc"))

    files2load = sorted(files2load) 
    
    var = 'ice_conc'

    # Open files with cftime decoder
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)

    result = []
    for file in files2load:    
    
        ds = xr.open_dataset(file, decode_times=time_coder).load()
        if log:
            print(f"File loaded: {file}", flush=True)

        # Compute sea ice area
        mask_geo = (
            (ds.lat >= box[2]) & (ds.lat <= box[3]) &
            (ds.lon >= box[0]) & (ds.lon <= box[1])
        )

        mask_isice = ds[var] > (siconc_threshold*100) # ice_conc is in %
    
        gridcell_area = 25000 * 25000 # m^2

        if diag=='area':
            ice_area = (((ds[var]/100) * mask_geo * mask_isice) * gridcell_area).sum(dim=('xc','yc'))
            name = 'sea_ice_area'
        elif diag=='extent':
            ice_area = ((mask_geo * mask_isice) * gridcell_area).sum(dim=('xc','yc'))
            name = 'sea_ice_extent'
        
        ice_area = ice_area.where(ice_area>0,np.nan)

        result.append(ice_area)

    result = xr.concat(result, dim='time').rename(name)

    result.attrs['units'] = 'm^2'
    result.attrs['long_name'] = f'total {name}'
    result.attrs['bounding box'] = (
        f"Longitude: {box[0]}E to {box[1]}E, Latitude: {box[2]}N to {box[3]}N"
    )

    if grouping == 'annual.mean':
        result = result.groupby('time.year').mean('time')
    elif grouping == 'annual.max':
        result = result.groupby('time.year').max('time')
    elif grouping == 'annual.min':
        result = result.groupby('time.year').min('time')
    elif grouping == 'monthly.mean':
        pass
    
    print('Done!')
    return result