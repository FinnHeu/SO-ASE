# so_ase/fesom2/eval_icebergs.py

import numpy as np
import xarray as xr
import pandas as pd
import re
import io
import os
from os.path import isfile
import glob
    
from .helpers_mesh import unrotate_coordinates
from ..miscellaneous.helpers_misc import seconds_per_month

def read_iceberg_initial_files(icebergpath):
    """Reads iceberg initial condition files from the specified directory.
    
    Parameters:
        icebergpath (str): Path to the directory containing iceberg initial condition files (icb*.dat).
        
    Returns:
        xarray.Dataset: Dataset containing iceberg initial conditions with dimension 'ib'.
    """

    def read_dat(filepath):
        data = []
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if line:
                    data.append(float(line))
        return np.array(data, dtype=float)
    
    icb_lon = read_dat(icebergpath + "/icb_longitude.dat")
    icb_lat = read_dat(icebergpath + "/icb_latitude.dat")
    icb_height = read_dat(icebergpath + "/icb_height.dat")
    icb_length = read_dat(icebergpath + "/icb_length.dat")
    icb_scaling = read_dat(icebergpath + "/icb_scaling.dat").astype(int)
    icb_calving_day = read_dat(icebergpath + "/icb_calving_day.dat").astype(int)

    # Create xarray dataset
    ds = xr.Dataset(
        data_vars={
            'lon_deg': (['ib'], icb_lon),
            'lat_deg': (['ib'], icb_lat),
            'height_ib': (['ib'], icb_height),
            'length_ib': (['ib'], icb_length),
            'scaling': (['ib'], icb_scaling),
            'calving_day': (['ib'], icb_calving_day)
        },
        coords={'ib': np.arange(1, len(icb_lon) + 1)},
        attrs={
            'description': 'Iceberg initial condition data',
            'source_path': icebergpath,
            'format': 'FESOM iceberg initial files'
        }
    )
    
    # Add variable attributes
    ds['lon_deg'].attrs = {'long_name': 'geographical longitude (unrotated)', 'units': 'degrees'}
    ds['lat_deg'].attrs = {'long_name': 'geographical latitude (unrotated)', 'units': 'degrees'}
    ds['height_ib'].attrs = {'long_name': 'iceberg height', 'units': 'm'}
    ds['length_ib'].attrs = {'long_name': 'iceberg length', 'units': 'm'}
    ds['scaling'].attrs = {'long_name': 'scaling factor', 'units': '1'}
    ds['calving_day'].attrs = {'long_name': 'calving day', 'units': 'days'}
    
    return ds


def iceberg_occurrence_heatmap(trackfile, lon_res=1.0, lat_res=1.0, lon_range=(-180, 180), lat_range=(-90, 90)):
    """Generate a heatmap of iceberg occurrence from trajectory data.
    
    Counts the number of iceberg occurrences per grid cell. If multiple icebergs
    are in the same cell at the same timestep, each is counted separately.
    
    Parameters:
        trackfile (str): Path to iceberg track NetCDF file (icb_track.nc) containing
            pos_lon_deg and pos_lat_deg variables with dims (time, number_tracer).
        lon_res (float): Longitude resolution in degrees. Default is 1.0.
        lat_res (float): Latitude resolution in degrees. Default is 1.0.
        lon_range (tuple): (min, max) longitude range. Default is (-180, 180).
        lat_range (tuple): (min, max) latitude range. Default is (-90, 90).
        
    Returns:
        xarray.Dataset: Dataset with 'iceberg_count' variable on regular lon-lat grid.
    """
    ds = xr.open_dataset(trackfile, decode_times=False)
    
    lon = ds['pos_lon_deg'].values.flatten()
    lat = ds['pos_lat_deg'].values.flatten()
    
    valid = np.isfinite(lon) & np.isfinite(lat) & ~((lon == 0) & (lat == 0))
    lon = lon[valid]
    lat = lat[valid]
    
    lon_bins = np.arange(lon_range[0], lon_range[1] + lon_res, lon_res)
    lat_bins = np.arange(lat_range[0], lat_range[1] + lat_res, lat_res)
    
    counts, _, _ = np.histogram2d(lon, lat, bins=[lon_bins, lat_bins])
    counts = counts.T.astype(np.int32)
    
    lon_centers = lon_bins[:-1] + lon_res / 2
    lat_centers = lat_bins[:-1] + lat_res / 2
    
    ds_out = xr.Dataset(
        data_vars={
            'iceberg_count': (['lat', 'lon'], counts)
        },
        coords={
            'lon': lon_centers,
            'lat': lat_centers
        },
        attrs={
            'description': 'Iceberg occurrence heatmap',
            'source_file': trackfile,
            'lon_resolution': lon_res,
            'lat_resolution': lat_res
        }
    )
    
    ds_out['iceberg_count'].attrs = {
        'long_name': 'iceberg occurrence count',
        'units': '1',
        'description': 'Number of iceberg occurrences per grid cell (summed over all timesteps and icebergs)'
    }
    
    ds.close()
    return ds_out


def read_iceberg_restart_file(icebergpath, unrotate=True):
    """Reads iceberg restart file and returns content as xarray dataset.
    
    Parameters:
        icebergpath (str): Path to the iceberg.restart.ISM file.
        unrotate (bool): Whether to unrotate coordinates to regular lat/lon (default: True)
        
    Returns:
        xarray.Dataset: Dataset containing all iceberg restart variables with dimension 'ib'.
    """
    # Define column names and their data types
    column_names = [
        'height_ib', 'length_ib', 'width_ib', 'lon_deg', 'lat_deg', 'Co', 'Ca', 'Ci',
        'Cdo_skin', 'Cda_skin', 'rho_icb', 'conc_sill', 'P_sill', 'rho_h2o', 'rho_air',
        'rho_ice', 'u_ib', 'v_ib', 'iceberg_elem', 'find_iceberg_elem', 'f_u_ib_old',
        'f_v_ib_old', 'calving_day', 'grounded', 'scaling', 'melted'
    ]
    
    # Read the file using numpy's loadtxt with structured dtype
    dtype = [
        ('height_ib', 'f8'), ('length_ib', 'f8'), ('width_ib', 'f8'),
        ('lon_deg', 'f8'), ('lat_deg', 'f8'), ('Co', 'f8'), ('Ca', 'f8'), ('Ci', 'f8'),
        ('Cdo_skin', 'f8'), ('Cda_skin', 'f8'), ('rho_icb', 'f8'), ('conc_sill', 'f8'),
        ('P_sill', 'f8'), ('rho_h2o', 'f8'), ('rho_air', 'f8'), ('rho_ice', 'f8'),
        ('u_ib', 'f8'), ('v_ib', 'f8'), ('iceberg_elem', 'i8'), ('find_iceberg_elem', 'U1'),
        ('f_u_ib_old', 'f8'), ('f_v_ib_old', 'f8'), ('calving_day', 'i8'),
        ('grounded', 'U1'), ('scaling', 'i8'), ('melted', 'U1')
    ]
    
    with open(icebergpath, 'r') as f:
        content = f.read()
    # Fix patterns like "0.1234567-310" -> "0.1234567E-310" and "0.1234567+100" -> "0.1234567E+100"
    content = re.sub(r'(\d)([+-])(\d{2,3})(\s|$)', r'\1E\2\3\4', content)
    
    # Read the data from the fixed content
    data = np.loadtxt(io.StringIO(content), dtype=dtype)
    
    # Convert boolean flags from characters to bool
    find_iceberg_elem = data['find_iceberg_elem'] == b'F'
    grounded = data['grounded'] == b'F'
    melted = data['melted'] == b'F'
    
    # Create data dictionary for xarray
    data_vars = {}
    for i, name in enumerate(column_names):
        if name in ['find_iceberg_elem', 'grounded', 'melted']:
            # Handle boolean flags
            if name == 'find_iceberg_elem':
                data_vars[name] = (['ib'], find_iceberg_elem)
            elif name == 'grounded':
                data_vars[name] = (['ib'], grounded)
            elif name == 'melted':
                data_vars[name] = (['ib'], melted)
        else:
            # Handle numeric values
            data_vars[name] = (['ib'], data[name])
    
    # Create xarray dataset
    ds = xr.Dataset(
        data_vars=data_vars,
        coords={'ib': np.arange(len(data))},
        attrs={
            'description': 'Iceberg restart file data',
            'source_file': icebergpath,
            'format': 'FESOM iceberg.restart.ISM'
        }
    )
    
    # Add variable attributes
    ds['height_ib'].attrs = {'long_name': 'iceberg height', 'units': 'm'}
    ds['length_ib'].attrs = {'long_name': 'iceberg length', 'units': 'm'}
    ds['width_ib'].attrs = {'long_name': 'iceberg width', 'units': 'm'}
    ds['lon_deg'].attrs = {'long_name': 'longitude', 'units': 'degrees'}
    ds['lat_deg'].attrs = {'long_name': 'latitude', 'units': 'degrees'}
    ds['Co'].attrs = {'long_name': 'drag coefficient', 'units': '1'}
    ds['Ca'].attrs = {'long_name': 'added mass coefficient', 'units': '1'}
    ds['Ci'].attrs = {'long_name': 'inertia coefficient', 'units': '1'}
    ds['Cdo_skin'].attrs = {'long_name': 'skin drag coefficient (ocean)', 'units': '1'}
    ds['Cda_skin'].attrs = {'long_name': 'skin drag coefficient (air)', 'units': '1'}
    ds['rho_icb'].attrs = {'long_name': 'iceberg density', 'units': 'kg/m^3'}
    ds['conc_sill'].attrs = {'long_name': 'concentration sill', 'units': '1'}
    ds['P_sill'].attrs = {'long_name': 'pressure sill', 'units': 'Pa'}
    ds['rho_h2o'].attrs = {'long_name': 'water density', 'units': 'kg/m^3'}
    ds['rho_air'].attrs = {'long_name': 'air density', 'units': 'kg/m^3'}
    ds['rho_ice'].attrs = {'long_name': 'ice density', 'units': 'kg/m^3'}
    ds['u_ib'].attrs = {'long_name': 'iceberg u-velocity', 'units': 'm/s'}
    ds['v_ib'].attrs = {'long_name': 'iceberg v-velocity', 'units': 'm/s'}
    ds['iceberg_elem'].attrs = {'long_name': 'iceberg element index', 'units': '1'}
    ds['find_iceberg_elem'].attrs = {'long_name': 'find iceberg element flag', 'units': '1'}
    ds['f_u_ib_old'].attrs = {'long_name': 'old u-force', 'units': 'N'}
    ds['f_v_ib_old'].attrs = {'long_name': 'old v-force', 'units': 'N'}
    ds['calving_day'].attrs = {'long_name': 'calving day', 'units': 'days'}
    ds['grounded'].attrs = {'long_name': 'grounded flag', 'units': '1'}
    ds['scaling'].attrs = {'long_name': 'scaling factor', 'units': '1'}
    ds['melted'].attrs = {'long_name': 'melted flag', 'units': '1'}
    
    # Unrotate coordinates if requested
    if unrotate:
        lon_unrot, lat_unrot = unrotate_coordinates(50.0, 15.0, -90.0, ds['lon_deg'].values, ds['lat_deg'].values)
        
        # Replace original coordinates with unrotated ones
        ds['lon_deg'] = (['ib'], lon_unrot)
        ds['lat_deg'] = (['ib'], lat_unrot)
        
        # Update attributes to indicate these are now unrotated coordinates
        ds['lon_deg'].attrs = {'long_name': 'longitude (unrotated)', 'units': 'degrees'}
        ds['lat_deg'].attrs = {'long_name': 'latitude (unrotated)', 'units': 'degrees'}
    
    return ds


def fesom_iceberg_heatflux_vertical_integral(
    src_path,
    dest_path,
    years=(1979, 2015),
    log=True
):
    """
    Vertically integrate iceberg heat flux from 3D FESOM output.

    For each grid node and time step, computes the vertical integral of the
    iceberg heat flux (ibhf) using layer thickness as weights.

    Parameters
    ----------
    src_path : str
        Path to the directory containing FESOM output files.
        Files are expected as `ibhf.fesom.{year}.nc`.
    dest_path : str
        Path to the directory where output files will be saved.
        Output files are named `ibhf_vertint.fesom.{year}.nc`.
    mesh_diag_path : str
        Path to the directory containing `fesom.mesh.diag.nc`.
    years : tuple of int, optional
        Year range (start, end) to process. Default is (1979, 2015).
    log : bool, optional
        If True, print progress messages. Default is True.

    Returns
    -------
    None
        Output is written to NetCDF files at dest_path.

    Notes
    -----
    - Existing output files are skipped (not overwritten).
    - Layer thickness is computed from mesh_diag.nz (depth levels).
    - The vertical integral sums ibhf * layer_thickness over all depth levels.

    Example
    -------
    >>> fesom_iceberg_heatflux_vertical_integral(
    ...     src_path='/path/to/fesom/output/',
    ...     dest_path='/path/to/output/',
    ...     mesh_diag_path='/path/to/mesh/',
    ...     years=(2000, 2010)
    ... )
    """
    os.makedirs(dest_path, exist_ok=True)
    
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
    
    for year in range(years[0], years[-1] + 1):
        input_file = f"{src_path}ibhf.fesom.{year}.nc"
        output_file = f"{dest_path}ibhf_vertint.fesom.{year}.nc"

        if isfile(output_file):
            if log:
                print(f"Skipping (exists): {output_file}", flush=True)
            continue

        if not isfile(input_file):
            if log:
                print(f"Input file not found: {input_file}", flush=True)
            continue

        if log:
            print(f"Processing: {input_file}", flush=True)

        ds = xr.open_dataset(input_file, decode_times=time_coder)
        
        # Compute vertical integral over depth
        vertint = ds['ibhf'].sum(dim='nz1')

        # Create output dataset
        ds_out = xr.Dataset(
            {
                'ibhf': vertint.astype(np.float32),
            },
            attrs={
                "description": f"Vertically integrated iceberg heat flux",
                "source_file": input_file,
            }
        )
        
        ds_out['ibhf'].attrs = {
            'long_name': 'vertically integrated iceberg heat flux',
            'units': 'W/m^2',
        }

        ds_out.to_netcdf(output_file)
        if log:
            print(f"Saved: {output_file}", flush=True)

        ds.close()
    
    if log:
        print("Done.", flush=True)


def fesom_total_iceberg_volume(
    src_path,
    dest_path,
    years=(1979, 2015),
    log=True
):
    """
    Compute total iceberg volume from buoys_track files.

    Loads buoys_track*.nc files, adds proper time dimension (12-hourly data),
    resamples to monthly means, and computes total iceberg volume as
    sum(height * length^2) over all icebergs.

    Parameters
    ----------
    src_path : str
        Path to the directory containing FESOM output files.
        Files are expected as `buoys_track.nc_{year}0101-{year}1231`.
    dest_path : str
        Path to the directory where output files will be saved.
        Output files are named `icb_vol.fesom.{year}.nc`.
    years : tuple of int, optional
        Year range (start, end) to process. Default is (1979, 2015).
    log : bool, optional
        If True, print progress messages. Default is True.

    Returns
    -------
    None
        Output is written to NetCDF files at dest_path.

    Notes
    -----
    - Existing output files are skipped (not overwritten).
    - Data is assumed to be 12-hourly and is resampled to monthly means.
    - Volume is computed as height * length^2 (assuming square cross-section).

    Example
    -------
    >>> fesom_total_iceberg_volume(
    ...     src_path='/path/to/fesom/output/',
    ...     dest_path='/path/to/output/',
    ...     years=(2000, 2010)
    ... )
    """
    os.makedirs(dest_path, exist_ok=True)
    
    for year in range(years[0], years[-1] + 1):
        input_pattern = f"{src_path}buoys_track.nc_{year}0101-{year}1231"
        output_file = f"{dest_path}icb_vol.fesom.{year}.nc"

        if isfile(output_file):
            if log:
                print(f"Skipping (exists): {output_file}", flush=True)
            continue

        # Find matching files
        input_files = glob.glob(input_pattern)
        if not input_files:
            if log:
                print(f"Input file not found: {input_pattern}", flush=True)
            continue

        if log:
            print(f"Processing: {input_files[0]}", flush=True)

        ds = xr.open_mfdataset(input_pattern, decode_times=False)
        
        # Add proper time dimension (12-hourly data)
        ds['time'] = pd.date_range(f'{year}-01-01', f'{year}-12-31T23:59:00', freq='12h')
        
        # Resample to monthly means
        ds_res = ds.resample(time='ME').mean()
        
        # Compute total iceberg volume: sum(height * length^2) over all icebergs
        icb_vol = (ds_res.height * ds_res.length**2).sum(dim='number_tracer')

        # Create output dataset
        ds_out = xr.Dataset(
            {
                'icb_vol': icb_vol.astype(np.float32),
            },
            attrs={
                "description": "Total iceberg volume (monthly means)",
                "source_file": input_files[0],
            }
        )
        
        ds_out['icb_vol'].attrs = {
            'long_name': 'total iceberg volume',
            'units': 'm^3',
        }

        ds_out.to_netcdf(output_file)
        if log:
            print(f"Saved: {output_file}", flush=True)

        ds.close()
    
    if log:
        print("Done.", flush=True)


def fesom_iceberg_combine_components(src_path, dest_path, years=(1979, 2015), log=True):
    """
    Combine iceberg freshwater flux components (ibfwb, ibfwbv, ibfwe, ibfwl) 
    into a single ibfw file per year.

    Parameters
    ----------
    src_path : str
        Directory containing the annual FESOM iceberg component files 
        (e.g., ibfwb.fesom.<year>.nc).
    dest_path : str
        Directory where the combined ibfw files will be saved.
    years : tuple of int, default (1979, 2015)
        Start and end years (inclusive) to process.
    log : bool, default True
        If True, print progress messages.

    Returns
    -------
    None
        Results are saved as NetCDF files named `ibfw.fesom.<year>.nc`.
    """
    os.makedirs(dest_path, exist_ok=True)
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
    
    for year in range(years[0], years[-1] + 1):
        outfile = f"{dest_path}ibfw.fesom.{year}.nc"
        if isfile(outfile):
            if log:
                print(f"Skipping (exists): {outfile}", flush=True)
            continue
            
        components = ['ibfwb', 'ibfwbv', 'ibfwe', 'ibfwl']
        ds_sum = None
        
        found_any = False
        for comp in components:
            infile = f"{src_path}{comp}.fesom.{year}.nc"
            if not isfile(infile):
                continue
            
            found_any = True
            ds = xr.open_dataset(infile, decode_times=time_coder)
            if ds_sum is None:
                ds_sum = ds.copy()
                ds_sum = ds_sum.rename({comp: 'ibfw'})
            else:
                ds_sum['ibfw'] = ds_sum['ibfw'] + ds[comp]
            ds.close()
            
        if found_any and ds_sum is not None:
            ds_sum['ibfw'].attrs = {
                'long_name': 'total iceberg freshwater flux',
                'units': 'm/s',
                'description': 'Sum of ibfwb, ibfwbv, ibfwe, and ibfwl'
            }
            ds_sum.to_netcdf(outfile)
            if log:
                print(f"Saved: {outfile}", flush=True)
            ds_sum.close()


def fesom_iceberg_integrated_flux(
    src_path, 
    mesh_diag_path, 
    varname='ibfw', 
    mask=None, 
    mask_name='global',
    years=(1979, 2015), 
    log=True, 
    savepath='./'
):
    """
    Compute and save integrated iceberg flux (freshwater or heat) time series.

    Integrates the specified flux variable over all nodes or a subset of nodes 
    defined by a mask. The flux is multiplied by the nodal area before summing.

    Parameters
    ----------
    src_path : str
        Directory containing the annual FESOM output files.
    mesh_diag_path : str
        Directory containing `fesom.mesh.diag.nc` for nodal areas.
    varname : str, default 'ibfw'
        Variable name to integrate (e.g., 'ibfw', 'ibhf').
    mask : array-like, optional
        Boolean mask or array of indices along the `nod2` dimension for 
        regional integration. If None, computes global integration.
    mask_name : str, default 'global'
        Descriptive name of the mask used in the output filename.
    years : tuple of int, default (1979, 2015)
        Start and end years (inclusive) to process.
    log : bool, default True
        If True, print progress messages.
    savepath : str, default './'
        Directory where output NetCDF files will be saved.

    Returns
    -------
    None
        Results are saved as NetCDF files named `iceberg_<varname>_<mask_name>.<year>.nc`.
    """
    os.makedirs(savepath, exist_ok=True)
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
    
    mesh_diag = xr.open_dataset(f"{mesh_diag_path}fesom.mesh.diag.nc")
    nod_area = mesh_diag.nod_area.max(dim='nz')
    
    for year in range(years[0], years[-1] + 1):
        infile = f"{src_path}{varname}.fesom.{year}.nc"
        outfile = f"{savepath}iceberg_{varname}_{mask_name}.{year}.nc"
        
        if isfile(outfile):
            if log:
                print(f"Skipping (exists): {outfile}", flush=True)
            continue
            
        if not isfile(infile):
            continue
            
        if log:
            print(f"Integrating {varname} ({mask_name}) for {year}...", flush=True)
            
        ds = xr.open_dataset(infile, decode_times=time_coder)
        data = ds[varname]
        
        # Handle vertical dimension if present
        if 'nz1' in data.dims:
            data = data.sum(dim='nz1')
        elif 'nz' in data.dims:
            data = data.sum(dim='nz')
            
        if mask is not None:
            integrated = (data * nod_area).isel(nod2=mask).sum(dim='nod2')
        else:
            integrated = (data * nod_area).sum(dim='nod2')
            
        ds_out = xr.Dataset(
            {
                f"{varname}_integrated": integrated.astype(np.float32),
            },
            attrs={
                "description": f"Integrated iceberg {varname} ({mask_name})",
                "source_file": infile,
                "mask_name": mask_name,
            }
        )
        
        ds_out.to_netcdf(outfile)
        if log:
            print(f"Saved: {outfile}", flush=True)
        ds.close()
    
    mesh_diag.close()


def icebergflux_to_massflux_Gty(src_path, dst_path, varname='ibfw', rho_fw=1000, year=None, log=True):
    """
    Convert spatially integrated iceberg flux time series into annual totals.

    For freshwater flux (`ibfw`), converts volume flux (m³/s) to annual 
    integrated mass flux (Gt/yr). For heat flux, computes the annual mean rate (W).

    Parameters
    ----------
    src_path : str
        Directory containing the integrated iceberg flux files.
    dst_path : str
        Directory where the annual total files will be saved.
    varname : str, default 'ibfw'
        Variable name processed ('ibfw' or 'ibhf').
    rho_fw : float, default 1000
        Density of freshwater in kg/m³, used to convert volume to mass.
    year : int, optional
        If specified, only process files for this year.
    log : bool, default True
        If True, print progress messages.

    Returns
    -------
    None
        Results are saved as NetCDF files with the suffix `_GTY`.
    """
    os.makedirs(dst_path, exist_ok=True)
    pattern = f"iceberg_{varname}_*.{year}.nc" if year else f"iceberg_{varname}_*.nc"
    files = np.sort(glob.glob(os.path.join(src_path, pattern)))
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
    
    for file in files:
        basename = os.path.basename(file)
        if "_GTY" in basename or "_GTM" in basename: continue
        
        outfile = os.path.join(dst_path, basename.replace(".nc", "_GTY.nc"))
        if isfile(outfile):
            continue
            
        ds = xr.open_dataset(file, decode_times=time_coder)
        data = ds[f"{varname}_integrated"]
        
        if varname == 'ibfw':
            massflux_water = data * rho_fw  # kg/s
            year_val = int(data.groupby('time.year').mean().year.values)
            seconds = seconds_per_month(year_val)
            Gty = ((massflux_water * seconds).groupby('time.year').sum() * 1e-12)
            ds_out = xr.Dataset({"iceberg_melt_GTY": Gty}, coords={"year": Gty.year})
            ds_out["iceberg_melt_GTY"].attrs = {"units": "Gt/y", "long_name": "Iceberg freshwater flux in Gt/yr"}
        else:
            annual_mean = data.groupby('time.year').mean()
            ds_out = xr.Dataset({f"{varname}_GTY": annual_mean}, coords={"year": annual_mean.year})
            ds_out[f"{varname}_GTY"].attrs = {"units": "W", "long_name": f"Annual mean iceberg {varname}"}

        ds_out.to_netcdf(outfile)
        if log: print(f"Saved: {outfile}", flush=True)
        ds.close()


def icebergflux_to_massflux_Gtm(src_path, dst_path, varname='ibfw', rho_fw=1000, year=None, log=True):
    """
    Convert spatially integrated iceberg flux time series into monthly totals.

    For freshwater flux (`ibfw`), converts volume flux (m³/s) to monthly 
    integrated mass flux (Gt/month). For heat flux, keeps the monthly mean rate (W).

    Parameters
    ----------
    src_path : str
        Directory containing the integrated iceberg flux files.
    dst_path : str
        Directory where the monthly total files will be saved.
    varname : str, default 'ibfw'
        Variable name processed ('ibfw' or 'ibhf').
    rho_fw : float, default 1000
        Density of freshwater in kg/m³, used to convert volume to mass.
    year : int, optional
        If specified, only process files for this year.
    log : bool, default True
        If True, print progress messages.

    Returns
    -------
    None
        Results are saved as NetCDF files with the suffix `_GTM`.
    """
    os.makedirs(dst_path, exist_ok=True)
    pattern = f"iceberg_{varname}_*.{year}.nc" if year else f"iceberg_{varname}_*.nc"
    files = np.sort(glob.glob(os.path.join(src_path, pattern)))
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
    
    for file in files:
        basename = os.path.basename(file)
        if "_GTY" in basename or "_GTM" in basename: continue
        
        outfile = os.path.join(dst_path, basename.replace(".nc", "_GTM.nc"))
        if isfile(outfile):
            continue
            
        ds = xr.open_dataset(file, decode_times=time_coder)
        data = ds[f"{varname}_integrated"]
        
        if varname == 'ibfw':
            massflux_water = data * rho_fw  # kg/s
            year_val = int(data.groupby('time.year').mean().year.values)
            seconds = seconds_per_month(year_val)
            Gtm = (massflux_water * seconds * 1e-12)
            ds_out = xr.Dataset({"iceberg_melt_GTM": Gtm}, coords={"time": Gtm.time})
            ds_out["iceberg_melt_GTM"].attrs = {"units": "Gt/m", "long_name": "Iceberg freshwater flux in Gt/month"}
        else:
            ds_out = xr.Dataset({f"{varname}_GTM": data}, coords={"time": data.time})
            ds_out[f"{varname}_GTM"].attrs = {"units": "W", "long_name": f"Monthly mean iceberg {varname}"}

        ds_out.to_netcdf(outfile)
        if log: print(f"Saved: {outfile}", flush=True)
        ds.close()
