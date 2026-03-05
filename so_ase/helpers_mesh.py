# so_ase.helpers_mesh.py

import xarray as xr
import numpy as np
import shapely.geometry as sh
import math
from pyproj import Proj, Transformer
from collections import defaultdict, deque
from .helpers_misc import lon_to_360, read_kml_coords

###--------> 1. Read raw mesh files

def read_nodes(meshpath):
    """
    Reads 2D node coordinates from a `nod2d.out` file in a given mesh path.

    Parameters:
        meshpath (str): Path to the directory containing `nod2d.out`.

    Returns:
        node_lon (array): Array of node longitudes.
        node_lat (array): Array of node latitudes.
        node_idx (array): Array of node indices (0-based).
        node_coast (array): Array indicating if a node is coastal (1) or not (0).
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
            coast = int(parts[3])
            
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

def read_depth_zlev(meshpath):
    """
    Reads vertical layer depth information from depth_zlev.out file,
    
    Parameters:
        meshpath (str): Path to the directory containing `depth_zlev.out`.
    
    Returns:
        array of float: A list of depth values for each vertical layer.
        int: Number of vertical layers.
    """
    depth = []
    with open(f"{meshpath}depth_zlev.out", 'r') as f:
        num_layers = int(f.readline())
        for i in range(num_layers):
            depth.append(float(f.readline()))

    return np.array(depth), num_layers

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

def level_idx_to_depth(meshpath, which='seafloor', raw=False):
    """
    Converts level indices to depth values for either seafloor or cavity levels.
    
    Parameters:
        meshpath (str): Path to the directory containing mesh files.
        which (str): either <seafloor> or <cavity> for last active layer or first active layer. 
        raw (bool): If True, reads raw level indices. If False, reads processed level indices.
    
    Returns:
        array of float: A list of depth values, one for each element.
    """
    if which == 'seafloor':
        idx_layer = read_element_levels(meshpath, which='seafloor', raw=raw, python_indexing=True)
    elif which == 'cavity':
        idx_layer = read_element_levels(meshpath, which='cavity', raw=raw, python_indexing=True)
    
    depth_levels, num_layer = read_depth_zlev(meshpath)

    depth = depth_levels[idx_layer]

    return depth

def read_element_levels(meshpath, which='seafloor', raw=False, python_indexing=False):
    """
    Reads vertical level information (first active/last active layer index) from elvls.out/cavity_elvls.out file,

    Parameters:
        meshpath (str): Path to the directory containing elvls/cavity_elvls.out.
        which (str): either <seafloor> or <cavity> for last active layer or first active layer. 
        raw (bool): If True, reads raw level indices. If False, reads processed level indices.
        python_indexing (bool): If True, converts indices to 0-based indexing.

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
            if python_indexing:
                elvls.append(int(parts)-1)
            else:
                elvls.append(int(parts))
    return np.array(elvls)


###--------> 2. Build neighbors of nodes and elements

def build_element_neighbors(elements):
    """
    Builds a list of neighboring elements for each triangle in the mesh.
    Parameters:
        elements (array): An array of shape (ntri, 3) containing the node indices for each triangle.
    Returns:
        list: A list of lists, where each sublist contains the indices of neighboring triangles for the corresponding triangle.
    """
    ntri = elements.shape[0]
    edge_to_tri = defaultdict(list)
    
    for tidx, tri in enumerate(elements):
        a, b, c = tri
        edges = [(a, b), (b, c), (c, a)]
        for e1, e2 in edges:
            edge = tuple(sorted((e1, e2)))
            edge_to_tri[edge].append(tidx)
    
    ###---> Build Neighbors
    neighbors = [[] for _ in range(ntri)]
    for tidx, tri in enumerate(elements):
        a, b, c = tri
        edges = [(a, b), (b, c), (c, a)]
    
        for (e1, e2) in edges:
            edge = tuple(sorted((e1, e2)))
            attached = edge_to_tri[edge]
    
            # If edge is shared by two triangles
            if len(attached) == 2:
                neigh = attached[1] if attached[0] == tidx else attached[0]
                neighbors[tidx].append(neigh)

    return neighbors

def build_node_neighbors(elements, node_idx):
    """
    Builds neighboring nodes for each node in the mesh.
    
    Parameters:
        elements (array): An array of shape (ntri, 3) containing the node indices for each triangle.
        node_idx (array): An array of node indices.
    Returns:
        list: A list of lists, where each sublist contains the neighboring node indices for the corresponding node.
    """    
    N = node_idx.size
    neighbors = [set() for _ in range(N)]
    
    for a, b, c in elements:
        neighbors[a].update((b, c))
        neighbors[b].update((a, c))
        neighbors[c].update((a, b))

    neighbors = [list(s) for s in neighbors]

    return neighbors

def build_node_k_ring_neighbors(elements, node_idx, k):
    """
    Builds k-ring neighboring nodes for each node in the mesh.
    Parameters:
        elements (array): An array of shape (ntri, 3) containing the node indices for each triangle.
        node_idx (array): An array of node indices.
        k (int): The ring number to compute (e.g., 1 for 1-ring neighbors).
    Returns:
        list: A list of sets, where each set contains the k-ring neighboring node indices for the corresponding node.
    """

    N = node_idx.size
    neighbors_1 = [set() for _ in range(N)]
    
    for a, b, c in elements:
        neighbors_1[a].update((b, c))
        neighbors_1[b].update((a, c))
        neighbors_1[c].update((a, b))
    
    N = len(neighbors_1)
    k_ring = [set() for _ in range(N)]

    for node in range(N):
        visited = {node}            # avoid including the center node
        queue = deque([(node, 0)])

        while queue:
            current, dist = queue.popleft()
            if dist == k:     # stop expanding
                continue

            for nb in neighbors_1[current]:
                if nb not in visited:
                    visited.add(nb)
                    k_ring[node].add(nb)   # add to result
                    queue.append((nb, dist + 1))
    
    return k_ring

def build_elements_of_nodes(elements, node_idx):
    """
    For each node in the mesh, return the list of element indices
    (triangles) that contain that node.

    Parameters:
        elements (array): (ntri, 3) array of triangle vertex indices.
        node_idx (array): Array of node indices (e.g. np.arange(N)).

    Returns:
        list: A list of lists. The i-th list contains the triangle indices
              of all elements that include node i.
    """
    N = node_idx.size
    ntri = elements.shape[0]

    # Initialize a list of empty lists, one per node.
    elems_of_node = [[] for _ in range(N)]

    # Loop over triangles
    for tidx in range(ntri):
        a, b, c = elements[tidx]
        elems_of_node[a].append(tidx)
        elems_of_node[b].append(tidx)
        elems_of_node[c].append(tidx)

    return elems_of_node


###--------> 3. Select particular nodes/elements or build masks

def find_nodes_in_box(
        mesh_diag_path,
        box=[-180, 180, -90, -60],
        log=True
):
    """ 
    Finds node indices within a specified geographical box.
    Parameters:
        mesh_diag_path (str): Path to the directory containing `fesom.mesh.diag.nc`.
        box (list): List of [lon_min, lon_max, lat_min, lat_max].
        log (bool): If True, prints the number of found nodes.
    Returns:
        array: Array of node indices within the specified box.
    """
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

def build_cavity_mask(meshpath, which='element'):
    """
    Builds a cavity mask for either elements or nodes based on cavity levels.
    A node is considered in the cavity if all its connected elements are cavity elements.
    Inverting the mask gives the open ocean mask.

    Parameters:
        meshpath (str): Path to the directory containing mesh files.
        which (str): 'element' to build mask for elements, 'node' for nodes.

    Returns:
        array of bool: Cavity mask array.
    """
    lev_cav = read_element_levels(meshpath, which='cavity', raw=False, python_indexing=False)
    
    if which == 'element':
        cavity_mask = lev_cav > 1
    elif which == 'node':
        elements = read_elements(meshpath)
        node_stats = read_nodes(meshpath)
        elems_of_node = build_elements_of_nodes(elements, node_stats[2])

        cavity_mask = np.zeros(len(elems_of_node), dtype='bool')
        for i, e in enumerate(elems_of_node):
            if all(lev_cav[e] > 1):
                cavity_mask[i] = True
    
    return cavity_mask

def build_cavity_regional_mask(meshpath, kml_path, name='Filchner-Ronne', which='node'):
    """
    Builds a regional cavity mask for nodes based on a KML-defined polygon.
    
    Parameters:
        meshpath (str): Path to the directory containing mesh files.
        kml_path (str): Path to the directory containing KML files.
        name (str): Name of the region (used to select the KML file).
        which (str): 'element' to build mask for elements, 'node' for nodes.
    
    Returns:
        array of bool: Regional cavity mask.
    """

    
    if which == 'node':
        lon, lat, node_idx, node_coast = read_nodes(meshpath)
        mask_cavity = build_cavity_mask(meshpath, which='node')
    elif which == 'element':
        elem = read_elements(meshpath)
        node_lon, node_lat, node_idx, node_coast = read_nodes(meshpath)
        lon, lat  = node_lon[elem].mean(axis=1), node_lat[elem].mean(axis=1)
        lev_cav = read_element_levels(meshpath, which='cavity', raw=False, python_indexing=False)
        mask_cavity = lev_cav > 1
    else:
        raise ValueError(f"which must be 'node' or 'element', got {which}")

    kml_file = f"{kml_path}{name}.kml"        
    coords = read_kml_coords(kml_file, close_ring=True)

    if name == 'Ross':
        lon = lon_to_360(lon)
    #    coords = [(i+360, j) for i, j in coords if i < 0]
        
    polygon = sh.Polygon(coords)
    mask_region = mask_cavity.copy()
    for i, b in enumerate(mask_cavity):
        if b:
            point = sh.Point((lon[i], lat[i]))
            if not polygon.contains(point):
                mask_region[i] = False

    return mask_region

def build_runoff_basin_mask(meshpath, runoff_file, which='liquid'):
    """
    Assigns FESOM mesh nodes to runoff basins defined in a runoff_maps.nc file.
    
    Parameters:
        meshpath (str): Path to the directory containing mesh files (nod2d.out).
        runoff_file (str): Path to the runoff_maps.nc file containing basin definitions.
        which (str): Either 'liquid' (for arrival_point_id) or 'solid' (for calving_point_id).
    
    Returns:
        xarray.Dataset: Dataset containing basin IDs with node coordinates as dimensions.
                       Variables include 'basin_id' with basin ID for each mesh node.
                       Nodes not in any basin will have ID 0.
    """
    
    node_lon, node_lat, node_idx, node_coast = read_nodes(meshpath)
    
    ds_runoff = xr.open_dataset(runoff_file)
    
    if which == 'liquid':
        basin_var = 'arrival_point_id'
    elif which == 'solid':
        basin_var = 'calving_point_id'
    else:
        raise ValueError("which must be either 'liquid' or 'solid'")
    
    basin_data = ds_runoff[basin_var]
    basin_lon = ds_runoff['lon'].values
    basin_lat = ds_runoff['lat'].values

    node_lon = lon_to_360(node_lon) # basin_lon is defined on 0-360E, node_lon is on -180-180E
    
    basin_ids = np.zeros(len(node_lon), dtype=int)
    
    for i in range(len(node_lon)):
        lon_node = node_lon[i]
        lat_node = node_lat[i]
        
        lon_diff = np.abs(basin_lon - lon_node)
        lat_diff = np.abs(basin_lat - lat_node)
        
        lon_idx = np.argmin(lon_diff)
        lat_idx = np.argmin(lat_diff)
        
        basin_id = int(basin_data.isel(lon=lon_idx, lat=lat_idx).values)
        basin_ids[i] = basin_id
    
    ds_runoff.close()
    
    # Create xarray dataset with basin IDs
    ds_basin = xr.Dataset(
        {
            'basin_id': (['nod2'], basin_ids)
        },
        coords={
            'nod2': node_idx,
            'lon': (['nod2'], node_lon),
            'lat': (['nod2'], node_lat),
            'coast': (['nod2'], node_coast)
        },
        attrs={
            'description': f'Runoff basin mask for {which} runoff',
            'basin_type': which,
            'mesh_path': meshpath,
            'runoff_file': runoff_file
        }
    )
    
    # Add attributes to basin_id variable
    ds_basin['basin_id'].attrs = {
        'long_name': 'Runoff basin identifier',
        'units': '1',
        'description': f'Basin ID for {which} runoff. Nodes not in any basin have ID 0.'
    }
    
    return ds_basin


###--------> 4. Add mesh diagnostics to fesom.mesh.diag.nc

def add_nodal_volumes(mesh_diag):
    """
    Adds layer thickness and nodal volume to the mesh_diag xarray.Dataset.
    Parameters:
        mesh_diag (xarray.Dataset): Dataset containing 'nz' and 'nod_area'.
    Returns:
        xarray.Dataset: Updated dataset with 'layer_thickness' and 'nod_volume'.
    """
    mesh_diag['layer_thickness'] = (('nz1'), np.diff(mesh_diag.nz))
    mesh_diag['nod_volume'] = mesh_diag.nod_area.isel(nz=0) * mesh_diag.layer_thickness
    return mesh_diag

def add_element_volumes(mesh_diag):
    """
    Adds layer thickness and element volume to the mesh_diag xarray.Dataset.
    Parameters:
        mesh_diag (xarray.Dataset): Dataset containing 'nz' and 'nod_area'.
    Returns:
        xarray.Dataset: Updated dataset with 'layer_thickness' and 'elem_volume'.
    """
    mesh_diag['layer_thickness'] = (('nz1'), np.diff(mesh_diag.nz))
    mesh_diag['elem_volume'] = mesh_diag.elem_area * mesh_diag.layer_thickness
    return mesh_diag


###--------> 5. Miscellaneous mesh helper functions

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

def mesh2vtk(meshpath, which='seafloor'):
    """
    Converts a FESOM2 mesh defined by `nod2d.out` and `elem2d.out` files into VTK format.
    
    Parameters:
        meshpath (str): Path to the directory containing `nod2d.out` and `elem2d.out`.
        which (str): 'seafloor' to include last active layer, 'cavity' for first active layer, 'none' for none of both.
    
    Returns:
        None: Writes a VTK file named 'mesh_output.vtk' in the current directory.
    """


    EARTH_RADIUS_KM = 6371.0

    def read_nodes_for_vtk(filename):
        with open(filename, 'r') as f:
            num_nodes = int(f.readline())
            nodes = []
            for _ in range(num_nodes):
                parts = f.readline().split()
                node_id = int(parts[0])
                lon = float(parts[1])
                lat = float(parts[2])
                # Convert to radians
                lon_rad = math.radians(lon)
                lat_rad = math.radians(lat)
                # Convert to Cartesian
                x = EARTH_RADIUS_KM * math.cos(lat_rad) * math.cos(lon_rad)
                y = EARTH_RADIUS_KM * math.cos(lat_rad) * math.sin(lon_rad)
                z = EARTH_RADIUS_KM * math.sin(lat_rad)
                nodes.append((x, y, z))
        return nodes

    def read_elements_for_vtk(filename):
        with open(filename, 'r') as f:
            num_elems = int(f.readline())
            elements = []
            for _ in range(num_elems):
                parts = f.readline().split()
                # Convert to 0-based indexing
                n1 = int(parts[0]) - 1
                n2 = int(parts[1]) - 1
                n3 = int(parts[2]) - 1
                elements.append((n1, n2, n3))
        return elements

    def write_vtk(filename, nodes, elements, cell_data=None):
        num_cells = len(elements)

        # --- Basic mesh output ---
        with open(filename, 'w') as f:
            f.write("# vtk DataFile Version 3.0\n")
            f.write("Mesh converted from nod2d and elem2d\n")
            f.write("ASCII\n")
            f.write("DATASET UNSTRUCTURED_GRID\n")

            # Points
            f.write(f"POINTS {len(nodes)} float\n")
            for x, y, z in nodes:
                f.write(f"{x} {y} {z}\n")

            # Triangles
            size = num_cells * 4  # 3 vertices + size indicator
            f.write(f"CELLS {num_cells} {size}\n")
            for n1, n2, n3 in elements:
                f.write(f"3 {n1} {n2} {n3}\n")

            f.write(f"CELL_TYPES {num_cells}\n")
            for _ in range(num_cells):
                f.write("5\n")  # VTK_TRIANGLE = 5

            # --- CELL_DATA section: optional multiple arrays ---
            if cell_data is not None and len(cell_data) > 0:
                f.write(f"CELL_DATA {num_cells}\n")

                for name, values in cell_data.items():
                    if len(values) != num_cells:
                        raise ValueError(
                            f"Cell data array '{name}' has length {len(values)}, "
                            f"but mesh has {num_cells} cells."
                        )

                    f.write(f"SCALARS {name} float\n")
                    f.write("LOOKUP_TABLE default\n")

                    for v in values:
                        f.write(f"{float(v)}\n")

    # Main execution
    nodes = read_nodes_for_vtk("f{meshpath}nod2d.out")
    elements = read_elements_for_vtk("f{meshpath}elem2d.out")
    if which == 'seafloor':
        depth = read_element_levels(meshpath, which='seafloor', raw=False)
        cell_data = {
            "last active layer": depth
            }
    elif which == 'cavity':
        depth = read_element_levels(meshpath, which='cavity', raw=False)
        cell_data = {
            "first active layer": depth
            }
    elif which == 'none':
        cell_data = None

    write_vtk("mesh_output_DARS2ice_seafloor.vtk", nodes, elements, cell_data)
    return