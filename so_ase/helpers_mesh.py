# so_ase.helpers_mesh.py

import xarray as xr
import numpy as np
from pyproj import Proj, Transformer

def read_nodes(meshpath):
    """
    Reads 2D node coordinates from a `nod2d.out` file in a given mesh path.

    Parameters:
        meshpath (str): Path to the directory containing `nod2d.out`.

    Returns:
        array: A list of (longitude, latitude) tuples for each node.
    """
    with open(f'{meshpath}nod2d.out', 'r') as f:
        num_nodes = int(f.readline())
        nodes_lon = []
        nodes_lat = []
        nodes_idx = []
        nodes_coast = []
        
        for _ in range(num_nodes):
            parts = f.readline().split()
            idx = int(parts[0]) - 1
            lon = float(parts[1])
            lat = float(parts[2])
            coast = bool(parts[3])
            
            nodes_lon.append(lon)
            nodes_lat.append(lat)
            nodes_idx.append(idx)
            nodes_coast.append(coast)
            
    return np.array(nodes_lon), np.array(nodes_lat), np.array(nodes_idx), np.array(nodes_coast)

def read_elements(meshpath):
    """
    Reads element connectivity information from a `elem2d.out` file.

    Parameters:
        meshpath (str): Path to the directory containing `elem2d.out`.

    Returns:
        array: A list of elements, where each element is a tuple
                       of 0-based node indices (n1, n2, n3).
    """
    with open(f'{meshpath}elem2d.out', 'r') as f:
        num_elems = int(f.readline())
        elements = []
        for _ in range(num_elems):
            parts = f.readline().split()
            # Convert to 0-based indexing
            n1 = int(parts[0]) - 1
            n2 = int(parts[1]) - 1
            n3 = int(parts[2]) - 1
            elements.append((n1, n2, n3))
    return np.array(elements)

def read_aux3d(meshpath):
    """
    Reads vertical level information (e.g., depths) from an `aux3d.out` file,
    skipping the first `num_levels` lines that typically contain header or
    level-wise data, and returns one depth value per node.

    Parameters:
        meshpath (str): Path to the directory containing `aux3d.out` and `nod2d.out`.

    Returns:
        array of float: A list of depth values, one for each node.
    """
    with open(f'{meshpath}nod2d.out', 'r') as f:
        num_nodes = int(f.readline())
    with open(f'{meshpath}aux3d.out', 'r') as f:
        num_levels = int(f.readline())
        depths = []
        for _ in range(num_nodes + num_levels):
            parts = f.readline()
            d = float(parts)
            depths.append(d)
    depths = depths[num_levels:]  # Skip level header values
    return np.array(depths)

def read_cavity_depth_at_node(meshpath):
    """
    Reads ice base depth information from cavity_depth@node.out file,

    Parameters:
        meshpath (str): Path to the directory containing `aux3d.out` and `nod2d.out`.

    Returns:
        list of float: A list of depth values, one for each node.
    """
    with open(f'{meshpath}nod2d.out', 'r') as f:
        num_nodes = int(f.readline())
    with open(f'{meshpath}cavity_depth@node.out', 'r') as f:
        depths = []
        for _ in range(num_nodes):
            parts = f.readline()
            d = float(parts)
            depths.append(d)
    return np.array(depths)

def read_element_levels(meshpath, which='seafloor', raw=False):
    """
    Reads vertical level information (first active/last active layer index) from elvls.out/cavity_elvls.out file,

    Parameters:
        meshpath (str): Path to the directory containing elvls/cavity_elvls.out.
        which (str): either <seafloor> or <cavity> for last active layer or first active layer. 

    Returns:
        array: A list of level indices values, one for each element.
    """
    with open(f'{meshpath}elem2d.out', 'r') as f:
        num_elem = int(f.readline())
    if which == 'seafloor':
        filename = 'elvls.out'
        if raw:
            filename='elvls_raw.out'
    elif which == 'cavity':
        filename = 'cavity_elvls.out'
        if raw:
            filename='cavity_elvls_raw.out'
    with open(f'{meshpath}{filename}', 'r') as f:
        elvls = []
        for _ in range(num_elem):
            parts = f.readline()
            elvls.append(int(parts))
    return np.array(elvls)

def find_nodes_in_box(
        mesh_diag_path,
        box=[-180, 180, -90, -60],
        log=True
):
    """ """
    mesh_diag = xr.open_dataset(f"{mesh_diag_path}fesom.mesh.diag.nc")
    
    # Find indices of nodes within the specified box
    inds = np.where(
        (mesh_diag.lon > box[0])
        & (mesh_diag.lon < box[1])
        & (mesh_diag.lat > box[2])
        & (mesh_diag.lat < box[3])
    )[0]
    
    if log:
        print(f"Found {len(inds)} nodes in the specified box.", flush=True)
    
    return inds

def gridcell_area_hadley(ds, R=6371.0):
    """
    Compute and add grid cell area (km²) to an xarray.Dataset on a regular 1° lat-lon grid.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with 'lat' and 'lon' coordinates.
    R : float, optional
        Radius of the Earth in kilometers (default: 6371 km).

    Returns
    -------
    ds_out : xarray.Dataset
        New dataset with an additional variable 'cell_area'.
    """
    # Broadcast 2D lat-lon grids
    lat2d, lon2d = xr.broadcast(ds.latitude, ds.longitude)

    # Convert to radians
    dlat = np.deg2rad(1.0)
    dlon = np.deg2rad(1.0)
    lat_rad = np.deg2rad(lat2d)

    # Compute latitude bounds
    lat1 = lat_rad - dlat / 2
    lat2 = lat_rad + dlat / 2

    # Area formula for a spherical Earth
    area = (R ** 2) * dlon * (np.sin(lat2) - np.sin(lat1))
    area = np.abs(area) * 1e6

    # Create DataArray for area
    area_da = xr.DataArray(
        area,
        coords=lat2d.coords,
        dims=lat2d.dims,
        name="cell_area",
        attrs={"units": "m^2", "description": "Grid cell area assuming spherical Earth"}
    )

    # Add to dataset
    ds_out = ds.copy()
    ds_out["cell_area"] = area_da

    return ds_out

def reproject_to_latlon(ds, input_proj="+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84"):
    """
    Transforms x/y coordinates of a polar-stereographic dataset (NSIDC datasets) into latitude/longitude and adds them as 2D coordinate variables.

    Parameters:
        ds (xr.Dataset): Dataset in projected coordinates (e.g., EPSG:3412).
        proj_str (str): PROJ string for the dataset's current projection.

    Returns:
        xr.Dataset: Dataset with 2D 'lat' and 'lon' variables added.
    """
    x = ds['x'].values
    y = ds['y'].values
    x2d, y2d = np.meshgrid(x, y)

    # Define projections explicitly
    proj_stereo = Proj(input_proj)
    proj_geo = Proj(proj='latlong', datum='WGS84')
    
    transformer = Transformer.from_proj(proj_stereo, proj_geo, always_xy=True)

    lon2d, lat2d = transformer.transform(x2d, y2d)

    ds['lat'] = (('y','x'), lat2d)
    ds['lon'] = (('y','x'), lon2d)
    ds['lat'].attrs['units'] = 'degrees_north'
    ds['lon'].attrs['units'] = 'degrees_east'
    
    return ds
