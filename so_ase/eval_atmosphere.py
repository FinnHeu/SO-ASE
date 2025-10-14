import numpy as np
import xarray as xr
import pandas as pd
from eofs.standard import Eof
from datetime import datetime
import matplotlib.pyplot as plt
# __________________________________________________________________________________________________________________________________________

# Functions to compute (and plot) different Atmospheric Modes from reanalysis data sets:

# 1.1 Compute Zonal Wave 3 Index
# 1.2 Plot Zonal Wave 3 Index
# 2. SAM Index
# 3.1 Compute Nino3.4 Index
# 3.2 Plot Nino3.4 Index
# To-Do:
# Plotting SAM, Compute ASL

# 1.1 Compute Zonal Wave 3 Index
#__________________________________________________________________________________________________________________________________________
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

# 1.2 Plot Zonal Wave 3 Index
# __________________________________________________________________________________________________________________________________________
def plot_ZW3_index(out_ds):
    """
    Plot Zonal Wave 3 (ZW3) index components with 3-month rolling means.
    
    Parameters
    ----------
    out_ds : xr.Dataset
        Output from compute_ZW3_index(), must contain
        'pc1_normalized', 'pc2_normalized', 'zw3_magnitude', 'zw3_phase'.
    """

    # Extract normalized PCs and indices
    pc1_normalized = out_ds['pc1_normalized']
    pc2_normalized = out_ds['pc2_normalized']
    magnitude_index = out_ds['zw3_magnitude']
    phase_index = out_ds['zw3_phase']

    # Compute 3-month rolling means
    pc1_roll = pc1_normalized.rolling(time=3, center=True).mean()
    pc2_roll = pc2_normalized.rolling(time=3, center=True).mean()
    mag_roll = magnitude_index.rolling(time=3, center=True).mean()
    phase_roll = phase_index.rolling(time=3, center=True).mean()

    # Time vector
    time = pc1_normalized.time
    time_min, time_max = time_min, time_max = time[0], time[-1]

    # Function to add positive/negative shading
    def shade_pos_neg(ax, data, color_pos='red', color_neg='blue', alpha=0.8):
        ax.fill_between(time, 0, data, where=(data>0), facecolor=color_pos, alpha=alpha, interpolate=True)
        ax.fill_between(time, 0, data, where=(data<0), facecolor=color_neg, alpha=alpha, interpolate=True)

    # --- Plot setup ---
    fig, axes = plt.subplots(4, 1, figsize=(15, 10), sharex=True)

    # PC1
    axes[0].plot(time, pc1_normalized, color='grey', label='PC1', linewidth=0.5)
    axes[0].plot(time, pc1_roll, color='k', label='PC1 3-month rolling mean')
    shade_pos_neg(axes[0], pc1_roll)
    axes[0].set_ylabel('PC1 (normalized)')
    axes[0].set_ylim(-2.5, 2.5)
    axes[0].set_xlim(time_min, time_max)
    axes[0].grid(True, linestyle='--', linewidth=0.5, alpha=0.7)

    # PC2
    axes[1].plot(time, pc2_normalized, color='grey', label='PC2', linewidth=0.5)
    axes[1].plot(time, pc2_roll, color='k', label='PC2 3-month rolling mean')
    shade_pos_neg(axes[1], pc2_roll)
    axes[1].set_ylabel('PC2 (normalized)')
    axes[1].set_ylim(-2.5, 2.5)
    axes[1].set_xlim(time_min, time_max)
    axes[1].grid(True, linestyle='--', linewidth=0.5, alpha=0.7)

    # ZW3 magnitude
    axes[2].plot(time, magnitude_index, color='grey', label='ZW3 Magnitude', linewidth=0.5)
    axes[2].plot(time, mag_roll, color='k', label='ZW3 Magnitude 3-month rolling mean')
    axes[2].fill_between(time, mag_roll, 1.5, color='grey', alpha=0.95, where=(mag_roll > 1.5))
    axes[2].set_ylabel('ZW3 index (magnitude)')
    axes[2].set_xlim(time_min, time_max)
    axes[2].grid(True, linestyle='--', linewidth=0.5, alpha=0.7)

    # ZW3 phase
    axes[3].plot(time, phase_index, color='grey', label='ZW3 Phase', linewidth=0.5)
    axes[3].plot(time, phase_roll, color='k', label='ZW3 Phase 3-month rolling mean')
    axes[3].set_ylabel('ZW3 phase (deg)')
    axes[3].set_xlabel('Time')
    axes[3].set_ylim(-175, 175)
    axes[0].set_xlim(time_min, time_max)
    axes[3].grid(True, linestyle='--', linewidth=0.5, alpha=0.7)

    plt.tight_layout()
    plt.show()

# 2. SAM Index
# __________________________________________________________________________________________________________________________________________
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

# 3.1 Compute Nino3.4 Index
# __________________________________________________________________________________________________________________________________________
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

# 3.2 Plot Nino3.4 Index
# __________________________________________________________________________________________________________________________________________

def plot_nino34_index(nino34_index, title="Niño 3.4 SST Anomaly Index", highlight=True, threshold=0.5):
    """
    Plot the Niño 3.4 index with ENSO event shading and labels.

    Parameters
    ----------
    nino34_index : xr.DataArray
        Output from compute_nino34_index()
    title : str
        Plot title
    highlight : bool
        Whether to shade El Niño / La Niña periods
    threshold : float
        Threshold for moderate ENSO events (default 0.5°C)
    """
    import matplotlib.pyplot as plt
    import pandas as pd

    # Convert to pandas for easier plotting
    ts = nino34_index.to_pandas()

    plt.figure(figsize=(14,4))
    plt.plot(ts.index, ts.values, label='Niño 3.4 SST Anomaly', color='black', linewidth=1.5)

    # Plot thresholds
    plt.axhline(threshold, linestyle='--', color='red', linewidth=1)
    plt.axhline(-threshold, linestyle='--', color='blue', linewidth=1)

    if highlight:
        # Moderate events shading
        plt.fill_between(ts.index, ts.values, threshold, where=(ts.values >= threshold), color='red', alpha=0.5)
        plt.fill_between(ts.index, ts.values, -threshold, where=(ts.values <= -threshold), color='blue', alpha=0.5)

    plt.title(title)
    plt.ylabel("SST Anomaly (°C)")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.xlim(ts.index[0], ts.index[-1])
    plt.tight_layout()
    plt.show()