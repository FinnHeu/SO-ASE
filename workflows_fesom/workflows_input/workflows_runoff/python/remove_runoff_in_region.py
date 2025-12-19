import xarray as xr
import so_ase as so
import numpy as np
import glob

def remove_runoff_in_region(src_path_runoff, dst_path_runoff, src_path_mesh, box=[-180, 180, -90, -60]):
    """Remove runoff in a specified rectangular region defined by box = [lon_min, lon_max, lat_min, lat_max].
    Parameters
    ----------
    src_path_runoff : str
        Path to the source runoff files.
    dst_path_runoff : str
        Path to the destination runoff files.
    src_path_mesh : str
        Path to the mesh files.
    box : list
        List defining the rectangular region [lon_min, lon_max, lat_min, lat_max].
    """

    files = glob.glob(f"{src_path_runoff}*.nc")

    node_lon, mode_lat, node_idx, node_coast = so.read_nodes(src_path_mesh)

    idx_in_box = np.where((node_lon >= box[0]) & (node_lon <= box[1]) & (node_lat >= box[2]) & (node_lat <= box[3]))
    
    for file in files:
        ds = xr.open_dataset(file)
        print(f'Reading file: {file}')
        runoff_data = ds.runoff.values
        runoff_data[:,idx_in_box] = 0 
        filename = dst_path_runoff + file.split('/')[-1]
        print(f'Writing file: {filename}')
        print(' ')
        ds.to_netcdf(filename)

    return

# Example usage
src_path_mesh = '/albedo/work/user/fheukamp/PostDoc2/CORE2/inputs/mesh/'
src_path_runoff = '/albedo/work/user/fheukamp/PostDoc2/CORE2/inputs/runoff/JRA55/'
dst_path_runoff = '/albedo/work/user/fheukamp/PostDoc2/CORE2/inputs/runoff/JRA55_noSO/'

remove_runoff_in_region(src_path_runoff, dst_path_runoff, src_path_mesh, box=[-180, 180, -90, -60])
    