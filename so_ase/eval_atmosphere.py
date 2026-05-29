# so_ase.eval_atmosphere.py

import xarray as xr
import numpy as np
from eofs.xarray import Eof as Eof_xarray
from eofs.standard import Eof as Eof_standard
from datetime import datetime
import matplotlib.pyplot as plt 

def norm_180_lon(ds):
    """
    Roll longitudes to the –180 to 180° convention if the grid uses 0–360°.
    Works for any coordinate named 'lon' or 'longitude'.
    """
    if 'longitude' in ds.coords:
        ds = ds.rename({'longitude': 'lon'})
        ds = ds.rename({'latitude': 'lat'})
    lon = ds["lon"].values
    if lon.max() > 180:
        ds = ds.assign_coords(lon=((lon + 180) % 360) - 180)
        ds = ds.sortby("lon")
    return ds

def SAM_idx_EOF(src_paths, ref_period=(datetime(1979, 1, 1), datetime(2000, 12, 31)), savepath='./'):
    '''
    Compute monthly area-weighted EOF1 and PC1 of mean sea level pressure south of 20°S according
    to the NOAA SAM definition.
 
    See: https://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/history/method.shtml

    Parameters
        ----------
        src_paths : list or str
            List of file paths (for xr.open_mfdataset) OR a single netcdf with monthly means of mean sea level pressure. 
        ref_period : tuple of datetime, optional
            Start and end of the climatology reference period.
            Default is (1979-01-01, 2000-12-31) matching the NOAA reference dataset.
        savepath : str, optional
        If provided, saves the output dataset (EOF1 + PC1) to this path. Default is './'.
    
    Returns
    -------
    eof1 : xr.DataArray
        First EOF of MSL (as correlation).
    pc1 : xr.DataArray
        First principal component (SAM index).
    
    '''
    # --- Load data set & normalise variable names ---
    print(">>> Loading data from:", src_paths) 
    ds = xr.open_mfdataset(src_paths, combine='by_coords')
    
    if 'valid_time' in ds.coords:
        ds = ds.rename({'valid_time': 'time'})
    if 'psl' in ds.data_vars:
        ds = ds.rename({'psl': 'msl'})
    if 'latitude' in ds.coords:
        ds = ds.rename({'latitude': 'lat'})
    if 'longitude' in ds.coords:
        ds = ds.rename({'longitude': 'lon'})
    
    ds = ds.resample(time="MS").mean()
    ds = norm_180_lon(ds)
    print(f'Longitude range after normalisation: {ds.lon.values.min():.1f} -> {ds.lon.values.max():.1f}')
    print('Time Range:', ds.time.min().values, ds.time.max().values)
    print('Reference period: {} to {}'.format(ref_period[0].date(), ref_period[1].date()))

    # ---  Extract subset (90°S – 20°S) --- 
    ds = ds.sortby('lat').sel(lat=slice(-90, -20))
    msl = ds['msl']

    # --- Climatology & anomalies ---
    msl_ref = msl.sel(time=slice(ref_period[0], ref_period[1]))
    clim = msl_ref.mean(dim='time')
    ds_anom = (msl - clim).transpose('time', 'lat', 'lon') # # Eof requires time as the first axis
    
    # --- Area weights ---
    coslat = np.cos(np.deg2rad(ds_anom.lat.values)) 
    wgts = np.sqrt(coslat)[np.newaxis, :, np.newaxis]

    # --- EOF analysis ---
    print(">>> Performing EOF analysis...")
    ds_anom = ds_anom.chunk({'time': -1}) # necessary for EOF analysis
    solver = Eof_xarray(ds_anom, weights=wgts)
    eof1 = solver.eofsAsCorrelation(neofs=1)
    pc1 = solver.pcs(npcs=1, pcscaling=1)

    # --- Sign convention ---
    # if high latitude (65S) value of 1st EOF in South Pac sector is positive, multiply EOFs by -1 to show negative values over pole
    lat_idx = np.argmin(np.abs(ds_anom.lat.values - (-65)))
    lon_idx = np.argmin(np.abs(ds_anom.lon.values - 125))  # 125°E ~ S. Pacific or 235º for 0-360º
    # Flip sign if S. Pacific EOF > 0
    if eof1[0, lat_idx, lon_idx] > 0:
        eof1 = -eof1
        pc1 = -pc1
    
    # --- Metadata ---
    eof1.name = "EOF1"
    pc1.name = "PC1"
    eof1.attrs["long_name"] = "First EOF of MSL (correlation)"
    eof1.attrs["reference_period"] = f"{ref_period[0]}–{ref_period[1]}"
    pc1.attrs["long_name"] = "First principal component (SAM index)"
    pc1.attrs["pcscaling"] = 1

    # --- Output dataset ---
    ds_out = xr.Dataset({"EOF1": eof1, "PC1": pc1})

    # Plotting
    # --- EOF map ---
    print(">>> Plotting EOF1 and PC1...")
    ax1 = plt.subplot(2, 1, 1)
    eof1.isel(mode=0).plot(ax=ax1, cmap="RdBu_r", vmin=-1, vmax=1, cbar_kwargs={"label": "Correlation"},)
    ax1.set_title("EOF1 of MSL (SAM pattern)")

    # --- PC time series ---
    ax2 = plt.subplot(2, 1, 2)
    pc1.isel(mode=0).plot(ax=ax2, color="k", lw=1)
    ax2.axhline(0, color="gray", lw=0.8)
    ax2.set_ylabel("SAM Index (PC1)")
    ax2.set_title("Principal Component 1")
    ax2.grid(True)

    plt.tight_layout()
    plt.show()

    # --- Save ---
    filename = f'{savepath}/SAM_Index.nc'
    ds_out.to_netcdf(filename)
    print(f">>> Saved NetCDF to {filename}")

    print(">>> EOF analysis complete.")
    
    return eof1, pc1

def ZW3_index(src_paths, level=500, lat_range=(-70, -40), ref_period=(datetime(1979, 1, 1), datetime(2020, 12, 31)), savepath='./'):
    """
        Compute the Zonal Wave 3 Index after Goyal et al. 2022 (https://doi.org/10.1175/JCLI-D-21-0927.1)
        -> Monthly anomalies of meridional winds at 500 hPa and from 40S - 70S

        Parameters
        ----------
        src_paths : list or str
            List of file paths (for xr.open_mfdataset) OR a single netcdf file with meriodional winds (v) at 
            different pressure levels.
        level : int, default=500
            Pressure level in hPa to extract from the meridional wind field.
        lat_range : tuple, default=(-70, -40)
            Latitude bounds (south to north).
        time_range : tuple of str, default=('1979','2020')
            Time slice (format 'YYYY').
        savepath : str, optional
            If provided, saves the output dataset (EOF1 + PC1) to this path. Default is './'.

        Returns
        -------
        xr.Dataset
            Dataset containing ZW3 magnitude, phase, PC1, PC2 & variance.
        
    """

    # --- Read & preprocess data ---
    print(">>> Loading dataset...")
    ds = xr.open_mfdataset(src_paths, combine='by_coords')
    
    # --- Normalize coordinate names ---
    if 'valid_time' in ds.coords:
        ds = ds.rename({'valid_time': 'time'})
    if 'latitude' in ds.coords:
        ds = ds.rename({'latitude': 'lat'})
    if 'longitude' in ds.coords:
        ds = ds.rename({'longitude': 'lon'})
    lev_names = ['pressure_level', 'plev', 'isobaricInhPa', 'lev']
    for lname in lev_names:
         if lname in ds.coords:
            ds = ds.rename({lname: 'level'})
    for var in ['time', 'lat', 'lon', 'level']:
        if var not in ds.coords:
            raise ValueError(f"Required coordinate '{var}' not found in dataset.")
    
    ds = norm_180_lon(ds)

    print(">>> Variables:", list(ds.data_vars))
    print(f">>> Time:{str(ds.time.values[0])[:10]} to {str(ds.time.values[-1])[:10]}")
    print(f">>> Lat range : {float(ds.lat.min()):.1f} to {float(ds.lat.max()):.1f}")
    print(">>> Available pressure levels:", ds.level.values if "level" in ds else "None")

    # --- Subsample data according the Definition ---
    if 'v' not in ds.data_vars:
        raise ValueError("Dataset must contain 'v' variable for meridional wind.")

    print(f"\n>>> Selecting level={level}, ref_period={ref_period[0]} to {ref_period[1]}, lat={lat_range}")
    # Auto-detect slicing direction00
    if ds.lat.values[0] > ds.lat.values[-1]:
        lat_slice = slice(lat_range[1], lat_range[0])
    else:
        lat_slice = slice(lat_range[0], lat_range[1])
    v = ds.v.sel(level=level).sel(lat=lat_slice)

    # --- Compute Monthly anomalies ---
    v_ref = v.sel(time=slice(ref_period[0], ref_period[1]))
    v_clim = v_ref.groupby('time.month').mean(dim='time')
    v_anom = v.groupby('time.month') - v_clim

    # --- EOF computation ---
    lat         = v_anom.lat
    coslat      = np.cos(np.deg2rad(lat.values)).clip(0.,1.)
    wgts        = np.sqrt(coslat)[..., np.newaxis]
    print(">>> Running EOF analysis...")
    solver      = Eof_standard(v_anom.values, weights=wgts)
    pcs         = solver.pcs(npcs=6, pcscaling=1)
    var         = solver.varianceFraction() 
    print(f">>> Variance explained -- EOF1: {var[0]*100:.1f}%  EOF2: {var[1]*100:.1f}%")

    # --- Magnitude = sqrt((PC1)²+(PC2)²)---
    pc1 = xr.DataArray(pcs[:,0], coords=[v.time], name='pc1')
    pc2 = xr.DataArray(pcs[:,1], coords=[v.time], name='pc2')
    magnitude_index = np.sqrt((pc1**2 + pc2**2))
    magnitude_index.name = 'zw3_magnitude'
    
    # --- phase = arctan(PC2/PC1) ---
    ntime = len(v_anom.time)
    phase = np.zeros(len(v[:,0,0])) * np.nan # Prepare array
    
    for i in range(ntime):
        p1, p2 = pcs[i, 0], pcs[i, 1]
        if   p1 > 0 and p2 > 0:
            phase[i] =  np.arctan(p2 / p1) * 180 / np.pi
        elif p1 < 0 and p2 > 0:
            phase[i] = (np.arctan(p2 / p1) * 180 / np.pi) + 180
        elif p1 > 0 and p2 < 0:
            phase[i] =  np.arctan(p2 / p1) * 180 / np.pi
        elif p1 < 0 and p2 < 0:
            phase[i] = (np.arctan(p2 / p1) * 180 / np.pi) - 180
        # if p1 == 0 or p2 == 0: stays NaN (edge case, rare for continuous data)
    phase_index = xr.DataArray(phase, coords=[v_anom.time], name='zw3_phase')
    
   # --- Save ---
    ds_out = xr.Dataset({
        'zw3_magnitude': magnitude_index,
        'zw3_phase': phase_index, 
        'pc1': pc1,
        'pc2': pc2,
        'variance_fraction': xr.DataArray(var, dims=['mode'])})
    ds_out['zw3_magnitude'].attrs = {'long_name': 'ZW3 magnitude index (Goyal et al. 2022)'}
    ds_out['zw3_phase'].attrs = {'long_name': 'ZW3 phase index (Goyal et al. 2022)', 'units': 'degrees'}
    filename = f'{savepath}/ZW3_Index.nc'
    ds_out.to_netcdf(filename)
    print(f">>> Saved NetCDF to {filename}")

    print(">>> Done.")

    return ds_out

def nino34_index( src_paths, ref_period=(datetime(1979, 1, 1), datetime(2010, 12, 31)), nino_box=[-170, -120, -5, 5], savepath='./'):
    """
    Compute the Niño 3.4 SST Index. Average SST anomalies within the region 5°S–5°N, 170°W–120°W.

    Parameters
    ----------
    src_paths : list or str
        List of NetCDF files or a single file containing SST ('sst' or 'tos').
    ref_period : tuple of datetime
        Reference period for climatology (mean SST).
    nino_box : list
        [lon_min, lon_max, lat_min, lat_max] defining the region. Default is Niño3.4 region.
    savepath : str, optional
        If provided, saves the output dataset to this path. Default is './'.

    Returns
    -------
    xr.DataArray
        Time series of Niño 3.4 SST anomalies.
    """

    # --- Load data ---
    print(">>> Loading dataset...")
    ds = xr.open_mfdataset(src_paths, combine='by_coords')

    # --- Normalize coordinate names ---
    if 'valid_time' in ds.coords:
        ds = ds.rename({'valid_time': 'time'})
    if 'latitude' in ds.coords:
        ds = ds.rename({'latitude': 'lat'})
    if 'longitude' in ds.coords:
        ds = ds.rename({'longitude': 'lon'})
    ds = norm_180_lon(ds)
    if 'sst' in ds.data_vars:
        sst = ds['sst']
    elif 'tos' in ds.data_vars:
        sst = ds['tos'].rename('sst')
    else:
        raise ValueError("Dataset must contain 'sst' or 'tos' variable.")
    print(f">>> Variable  : {sst.name}")
    print(f">>> Time:{str(sst.time.values[0])[:10]} to {str(sst.time.values[-1])[:10]}")

    # --- Extract Niño region ---
    lon_min, lon_max, lat_min, lat_max = nino_box
    if sst.lat.values[0] > sst.lat.values[-1]:
        lat_slice = slice(lat_max, lat_min)
    else:
        lat_slice = slice(lat_min, lat_max)
    sst_region = sst.sel(
        lon=slice(lon_min, lon_max),
        lat=lat_slice)
    
    print(f"\n>>> Niño 3.4 box extracted:")
    print(f">>> Lat range : {float(sst_region.lat.min()):.1f} to {float(sst_region.lat.max()):.1f}")
    print(f">>> Lon range : {float(sst_region.lon.min()):.1f} to {float(sst_region.lon.max()):.1f}")
    if sst_region.lon.size == 0 or sst_region.lat.size == 0:
        raise ValueError(f"Niño 3.4 box {nino_box} returned empty selection. ")

    # --- Compute area-weighted mean ---
    wgts = np.cos(np.deg2rad(sst_region.lat))
    nino34 = sst_region.weighted(wgts).mean(dim=['lat', 'lon'])
    
    # --- Reference period climatology ---
    ref = nino34.sel(time=slice(ref_period[0], ref_period[1]))
    clim = ref.groupby('time.month').mean('time')
    nino34_index = nino34.groupby('time.month') - clim # Anomaly
    nino34_index.name = 'nino34_index'
    nino34_index.attrs = {
        'long_name'  : 'Niño 3.4 SST anomaly index',
        'units'      : sst.attrs.get('units', 'K or °C'),
        'ref_period' : f"{str(ref_period[0])[:10]} to {str(ref_period[1])[:10]}",
        'nino_box'   : str(nino_box)}
    coords_to_drop = [c for c in nino34_index.coords 
                  if c not in ('time',)]
    nino34_index = nino34_index.drop_vars(coords_to_drop)
    print(f">>> Ref_period : {ref.time.min()} to {ref.time.max()}")

    # --- Save ---
    filename = f'{savepath}/Nino34_Index.nc'
    nino34_index.to_netcdf(filename)
    print(f">>> Saved NetCDF to {filename}")

    return nino34_index 