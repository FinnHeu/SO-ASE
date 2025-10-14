import numpy as np
import xarray as xr
import pandas as pd
from eofs.standard import Eof
from datetime import datetime
import matplotlib.pyplot as plt

def compute_ZW3_index(
    src_paths,
    level=500,
    lat_range=(-70, -40),
    time_range=('1979', '2020'),
    neofs=6,
    save_to=None,
):
    """
        Compute the Zonal Wave 3 Index after Goyal et al. 2022 (https://doi.org/10.1175/JCLI-D-21-0927.1)
        -> Monthly anomalies of meridional winds at 500 hPa and from 40S - 70S

        Parameters
        ----------
        src_paths : list or str
            List of file paths (for xr.open_mfdataset) OR a single netcdf file with meriodional winds at 
            different pressure levels.
        level : int, default=500
            Pressure level in hPa to extract from the meridional wind field.
        lat_range : tuple, default=(-70, -40)
            Latitude bounds (south to north).
        time_range : tuple of str, default=('1979','2020')
            Time slice (format 'YYYY').
        neofs : int, default=6
            Number of EOFs/PCs to compute.
        save_to : str or None
            If provided, save output .csv-File to this path.

        Returns
        -------
        xr.Dataset
            Dataset containing ZW3 magnitude, phase, PC1, PC2, and optionally EOFs & variance.
    
        Notes
         -------
        
    """

    # --- Read & preprocess data ---
    print(">>> Loading dataset...")
    ds = xr.open_mfdataset(src_paths, combine='by_coords')
    print(ds)
    
    # --- Normalize coordinate names ---

    # Time coordinate
    if 'valid_time' in ds.coords:
        ds = ds.rename({'valid_time': 'time'})
    # Latitude
    if 'latitude' in ds.coords:
        ds = ds.rename({'latitude': 'lat'})
    # Longitude
    if 'longitude' in ds.coords:
        ds = ds.rename({'longitude': 'lon'})
    # Pressure level
    lev_names = ['pressure_level', 'plev', 'isobaricInhPa', 'lev']
    for lname in lev_names:
         if lname in ds.coords:
            ds = ds.rename({lname: 'level'})
    # Final check
    for var in ['time', 'lat', 'lon', 'level']:
        if var not in ds.coords:
            raise ValueError(f"Required coordinate '{var}' not found in dataset.")

    # --- Print diagnostics ---
    print("\n>>> Coordinate checks:")
    for var in ['time', 'lat', 'lon', 'level']:
        print(f" - {var} present:", var in ds.coords)

    print(">>> Latitude sample:", ds.lat.values[:5], "...", ds.lat.values[-5:])
    print(">>> Available pressure levels:", ds.level.values if "level" in ds else "None")
    print(">>> Variables:", list(ds.data_vars))

    # --- Subsample data according the Definition ---
    if 'v' not in ds.data_vars:
        raise ValueError("Dataset must contain 'v' variable for meridional wind.")

    print(f"\n>>> Selecting level={level}, time={time_range}, lat={lat_range}")
    # Auto-detect slicing direction
    if ds.lat.values[0] > ds.lat.values[-1]:
        lat_slice = slice(lat_range[1], lat_range[0])
    else:
        lat_slice = slice(lat_range[0], lat_range[1])
    v = ds.v.sel(level=level, time=slice(*time_range)).sel(lat=lat_slice)

    # --- Compute Monthly anomalies ---
    v_anom = v.groupby('time.month') - v.groupby('time.month').mean(dim='time')

    # --- EOF computation ---
    lat         = v_anom.lat
    coslat      = np.cos(np.deg2rad(lat.values)).clip(0.,1.)
    wgts        = np.sqrt(coslat)[..., np.newaxis]
    solver      = Eof(v_anom.values, weights=wgts)
    pcs         = solver.pcs(npcs=neofs, pcscaling=1)

    # --- Build PC into xarray ---
    pc1 = xr.DataArray(pcs[:, 0], coords=[v_anom.time], name='pc1')
    pc2 = xr.DataArray(pcs[:, 1], coords=[v_anom.time], name='pc2')

    # --- Compute indices (vectorized) ---
    # Strength on ZW3:          magnitude = sqrt((PC1)²+(PC2)²)
    magnitude = np.sqrt(pc1**2 + pc2**2)
    # Longitudinal Orientation: phase     = arctan(PC2/PC1)
    phase = np.arctan2(pc2, pc1) * 180 / np.pi

    # --- Normalize PC1 & PC2 --- 
    pc1_mean = pc1.mean(dim='time')
    pc1_std  = pc1.std(dim='time')
    pc1_normalized = (pc1 - pc1_mean) / pc1_std

    pc2_mean = pc2.mean(dim='time')
    pc2_std  = pc2.std(dim='time')
    pc2_normalized = (pc2 - pc2_mean) / pc2_std
    

   # --- Save outputs ---
    out = xr.Dataset({
        'zw3_magnitude': magnitude,
        'zw3_phase': phase, 
        'pc1': pc1,
        'pc1_normalized': pc1_normalized,
        'pc2': pc2,
        'pc2_normalized': pc2_normalized,
        'variance_fraction': xr.DataArray(solver.varianceFraction(), dims=['mode']),
    })

    # --- Save to .csv ---
    if save_to:
        df = out[['zw3_magnitude', 'zw3_phase', 'pc1', 'pc1_normalized', 'pc2', 'pc2_normalized']].to_dataframe()
        df = df.reset_index()
        df.to_csv(save_to, index=False)
        print(f">>> Saved CSV to {save_to}")

    print(">>> Done.")

    return out

def compute_SAM_index(
        src_paths,
        ref_period=(datetime(1982, 1, 1), datetime(2015, 12, 31)),
        save_to=None
):
    """
        Compute the Southern Annular Mode (SAM) Index from SLP following Marshall, 2003:

             SAM Index = (SLP40_anom / std40) - (SLP65_anom / std65)

        Parameters
        ----------
        src_paths : list or str
            List of NetCDF file paths or a single file containing monthly SLP ('msl').
        ref_period : tuple of datetime, default=(1982-01-01, 2015-12-31)
            Reference period for climatology.
        save_to : str or None
            If provided, save output to a CSV file.

        Returns
        -------
        xr.Dataset
            Dataset with slp anomalies and monthly SAM index.
        
        Notes
        -------
     """
    # --- Read & preprocess data ---
    print(">>> Loading dataset...")
    ds = xr.open_mfdataset(src_paths, combine='by_coords')
    print(ds)

    if 'msl' not in ds.data_vars:
        raise ValueError("Dataset must contain a 'msl' variable for mean sea level pressure.")

    # --- Normalize coordinate names ---

    # Time coordinate
    if 'valid_time' in ds.coords:
        ds = ds.rename({'valid_time': 'time'})
    # Latitude
    if 'latitude' in ds.coords:
        ds = ds.rename({'latitude': 'lat'})
    # Longitude
    if 'longitude' in ds.coords:
        ds = ds.rename({'longitude': 'lon'})

    # --- Extract zonal mean SLP at 40S and 65S ---
    slp40 = ds['msl'].sel(lat=-40, method='nearest').mean(dim='lon').drop_vars('lat')
    slp65 = ds['msl'].sel(lat=-65, method='nearest').mean(dim='lon').drop_vars('lat')
    print(f">>> Extracted zonal mean SLP: slp40 shape {slp40.shape}, slp65 shape {slp65.shape}")

    # --- Select reference period ---
    slp40_ref = slp40.sel(time=slice(ref_period[0], ref_period[1]))
    slp65_ref = slp65.sel(time=slice(ref_period[0], ref_period[1]))
    print(f">>> Reference period selected: {ref_period[0]} to {ref_period[1]}")
    print(f"    slp40_ref shape: {slp40_ref.shape}, slp65_ref shape: {slp65_ref.shape}")

    # --- Compute monthly climatology ---
    print(">>> Computing Climatology...")
    slp40_clim = slp40_ref.groupby('time.month').mean('time')
    slp65_clim = slp65_ref.groupby('time.month').mean('time')

    # --- Compute anomalies ---
    print(">>> Computing Anomalies..")
    slp40_anom = slp40.groupby('time.month') - slp40_clim
    slp65_anom = slp65.groupby('time.month') - slp65_clim

    # --- Compute reference std for normalization ---
    slp40_std = slp40_anom.sel(time=slice(ref_period[0], ref_period[1])).std(dim='time')
    slp65_std = slp65_anom.sel(time=slice(ref_period[0], ref_period[1])).std(dim='time')

    # --- SAM Index ---
    print(">>> Computing SAM index...")
    SAM_index = (slp40_anom / slp40_std) - (slp65_anom / slp65_std)
    SAM_index.name = 'SAM_index'
    

    out = xr.Dataset({
    'slp40_anom': slp40_anom,
    'slp65_anom': slp65_anom,
    'SAM_index': SAM_index
    })

    # --- Save to CSV if requested ---
    if save_to:
        df = SAM_index.to_dataframe().reset_index()
        df.to_csv(save_to, index=False)
        print(f">>> Saved CSV to {save_to}")

    print(">>> Done.")

    return out

def compute_nino34_index(
        src_paths,
        ref_period=(datetime(1982, 1, 1), datetime(2015, 12, 31)),
        save_to=None
):
    """
        Compute the Southern Annular Mode (SAM) Index from SLP following Marshall, 2003:

             SAM Index = (SLP40_anom / std40) - (SLP65_anom / std65)

        Parameters
        ----------
        src_paths : list or str
            List of NetCDF file paths or a single file containing monthly SLP ('msl').
        ref_period : tuple of datetime, default=(1982-01-01, 2015-12-31)
            Reference period for climatology.
        save_to : str or None
            If provided, save output to a CSV file.

        Returns
        -------
        xr.Dataset
            Dataset with slp anomalies and monthly SAM index.
        
        Notes
        -------
     """
    # --- Read & preprocess data ---
    ds = xr.open_mfdataset(src_paths, combine='by_coords')

    if 'msl' not in ds.data_vars:
        raise ValueError("Dataset must contain a 'msl' variable for mean sea level pressure.")

    # --- Normalize coordinate names ---

    # Time coordinate
    if 'valid_time' in ds.coords:
        ds = ds.rename({'valid_time': 'time'})
    # Latitude
    if 'latitude' in ds.coords:
        ds = ds.rename({'latitude': 'lat'})
    # Longitude
    if 'longitude' in ds.coords:
        ds = ds.rename({'longitude': 'lon'})

    # --- Extract zonal mean SLP at 40S and 65S ---
    slp40 = ds['msl'].sel(lat=-40, method='nearest').mean(dim='lon')
    slp65 = ds['msl'].sel(lat=-65, method='nearest').mean(dim='lon')

    # --- Select reference period ---
    slp40_ref = slp40.sel(time=slice(ref_period[0], ref_period[1]))
    slp65_ref = slp65.sel(time=slice(ref_period[0], ref_period[1]))

    # --- Compute monthly climatology ---
    slp40_clim = slp40_ref.groupby('time.month').mean('time')
    slp65_clim = slp65_ref.groupby('time.month').mean('time')

    # --- Compute anomalies ---
    slp40_anom = slp40.groupby('time.month') - slp40_clim
    slp65_anom = slp65.groupby('time.month') - slp65_clim

    # --- Compute reference std for normalization ---
    slp40_std = slp40_anom.sel(time=slice(ref_period[0], ref_period[1])).std(dim='time')
    slp65_std = slp65_anom.sel(time=slice(ref_period[0], ref_period[1])).std(dim='time')

    # --- SAM Index ---
    SAM_index = (slp40_anom / slp40_std) - (slp65_anom / slp65_std)
    SAM_index.name = 'SAM_index'

    out = xr.Dataset({
    'slp40_anom': slp40_anom,
    'slp65_anom': slp65_anom,
    'SAM_index': SAM_index
    })

    # --- Save to CSV if requested ---
    if save_to:
        df = SAM_index.to_dataframe().reset_index()
        df.to_csv(save_to, index=False)

    return out
