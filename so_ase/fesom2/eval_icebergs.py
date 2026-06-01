# so_ase.eval_icebergs.py

import numpy as np
import xarray as xr

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