# so_ase.eval_icebergs.py

import numpy as np
import xarray as xr
from .helpers_mesh import unrotate_coordinates

def read_iceberg_initial_files(icebergpath):
    """Reads iceberg initial condition files from the specified directory.
    
    Parameters:
        icebergpath (str): Path to the directory containing iceberg initial condition files.
        
    Returns:
        array: Arrays containing iceberg longitude, latitude, height, face element indices,
               length, and scaling factors.
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
    icb_felem = read_dat(icebergpath + "/icb_felem.dat").astype(int)
    icb_length = read_dat(icebergpath + "/icb_length.dat")
    icb_scaling = read_dat(icebergpath + "/icb_scaling.dat")

    return icb_lon, icb_lat, icb_height, icb_felem, icb_length, icb_scaling


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
    
    # Read the data
    data = np.loadtxt(icebergpath, dtype=dtype)
    
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

