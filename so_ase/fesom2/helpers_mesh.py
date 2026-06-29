# so_ase/fesom2/helpers_mesh.py

import xarray as xr
import numpy as np
import shapely.geometry as sh
import shapely
import math

from scipy.spatial import cKDTree
from pyproj import Proj, Transformer
from collections import defaultdict, deque
from scipy.interpolate import griddata
from ..miscellaneous.helpers_misc import lon_to_360, read_kml_coords

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

def read_elements(meshpath, return_coordinates=False):
    """
    Reads element connectivity information from a `elem2d.out` file.

    Parameters:
        meshpath (str): Path to the directory containing `elem2d.out`.
        return_coordinates (bool): If True, also return centroid coordinates of elements.

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
    
    elements = np.array(elements)
    
    if return_coordinates:
        # Get node coordinates
        node_lon, node_lat, _, _ = read_nodes(meshpath)
        # Calculate centroid coordinates for each element
        elem_lon = node_lon[elements].mean(axis=1)
        elem_lat = node_lat[elements].mean(axis=1)
        return elements, elem_lon, elem_lat
    else:
        return elements

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

def level_idx_to_depth(meshpath, which='seafloor', raw=False, write=False):
    """
    Converts level indices to depth values for either seafloor or cavity levels.
    
    Parameters:
        meshpath (str): Path to the directory containing mesh files.
        which (str): either <seafloor> or <cavity> for last active layer or first active layer. 
        raw (bool): If True, reads raw level indices. If False, reads processed level indices.
        write (bool): If True, writes depths line-by-line to depth@elem.out (seafloor) or cavity_depth@elem.out (cavity).
    
    Returns:
        array of float: A list of depth values, one for each element.
    """
    if which == 'seafloor':
        idx_layer = read_element_levels(meshpath, which='seafloor', raw=raw, python_indexing=True)
    elif which == 'cavity':
        idx_layer = read_element_levels(meshpath, which='cavity', raw=raw, python_indexing=True)
    
    depth_levels, num_layer = read_depth_zlev(meshpath)

    depth = depth_levels[idx_layer]

    if write:
        if which == 'seafloor':
            filename = f'{meshpath}depth@elem.out'
        elif which == 'cavity':
            filename = f'{meshpath}cavity_depth@elem.out'
        with open(filename, 'w') as f:
            for d in depth:
                f.write(f'{d}\n')
        print(f'Wrote: {filename}')

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

def build_elem_k_ring_neighbors(elements, k):
    """
    Builds k-ring neighboring elements for each element in the mesh.
    Parameters:
        elements (array): An array of shape (ntri, 3) containing the node indices for each triangle.
        k (int): The ring number to compute (e.g., 1 for 1-ring neighbors).
    Returns:
        list: A list of sets, where each set contains the k-ring neighboring element indices for the corresponding element.
    """
    neighbors_1 = build_element_neighbors(elements)
    
    N = len(neighbors_1)
    k_ring = [set() for _ in range(N)]

    for elem in range(N):
        visited = {elem}            # avoid including the center element
        queue = deque([(elem, 0)])

        while queue:
            current, dist = queue.popleft()
            if dist == k:     # stop expanding
                continue

            for nb in neighbors_1[current]:
                if nb not in visited:
                    visited.add(nb)
                    k_ring[elem].add(nb)   # add to result
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

def find_element_for_points(lon_points, lat_points, node_lon, node_lat, elements):
    """
    Find the mesh element index for each point (lon, lat).
    
    Parameters
    ----------
    lon_points : array-like
        Longitudes of points in degrees.
    lat_points : array-like
        Latitudes of points in degrees.
    node_lon : array
        Longitudes of mesh nodes.
    node_lat : array
        Latitudes of mesh nodes.
    elements : array
        Element connectivity array (n_elem, 3).
    
    Returns
    -------
    elem_indices : ndarray
        Element index (0-based) for each point. Returns -1 if point is not found in any element.
    """
    from scipy.spatial import cKDTree
    
    lon_points = np.asarray(lon_points)
    lat_points = np.asarray(lat_points)
    n_points = len(lon_points)
    
    elem_indices = np.full(n_points, -1, dtype=np.int64)
    
    def point_in_triangle(px, py, x1, y1, x2, y2, x3, y3):
        denom = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)
        if abs(denom) < 1e-12:
            return False
        a = ((y2 - y3) * (px - x3) + (x3 - x2) * (py - y3)) / denom
        b = ((y3 - y1) * (px - x3) + (x1 - x3) * (py - y3)) / denom
        c = 1 - a - b
        return (a >= -1e-10) and (b >= -1e-10) and (c >= -1e-10)
    
    # Build element centroids for quick filtering
    elem_lon = node_lon[elements].mean(axis=1)
    elem_lat = node_lat[elements].mean(axis=1)
    
    # Build KDTree on element centroids for fast nearest-neighbor search
    tree = cKDTree(np.column_stack([elem_lon, elem_lat]))
    
    for i in range(n_points):
        px, py = lon_points[i], lat_points[i]
        
        # Find nearest elements to check
        _, candidate_elems = tree.query([px, py], k=min(100, len(elements)))
        
        for eidx in candidate_elems:
            n1, n2, n3 = elements[eidx]
            x1, y1 = node_lon[n1], node_lat[n1]
            x2, y2 = node_lon[n2], node_lat[n2]
            x3, y3 = node_lon[n3], node_lat[n3]
            
            if point_in_triangle(px, py, x1, y1, x2, y2, x3, y3):
                elem_indices[i] = eidx
                break
    
    return elem_indices

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

def build_regional_mask(meshpath, kml_file, which='node', close_ring=True, polar='south'):
    """
    Builds a binary mask for FESOM nodes or elements based on a KML-defined polygon.
    
    Points/elements inside the polygon are marked True, outside are marked False.
    The polygon can be non-closed in the KML file; it will be closed automatically
    if close_ring=True.
    
    For polar regions (especially polygons wrapping around the pole or crossing the
    antimeridian), the function projects coordinates to a stereographic projection
    before performing the point-in-polygon test.
    
    Parameters:
        meshpath (str): Path to the directory containing mesh files (nod2d.out, elem2d.out).
        kml_file (str): Full path to the KML file containing the polygon definition.
        which (str): 'node' to build mask for nodes, 'element' for element centroids.
        close_ring (bool): If True, closes the polygon by appending the first coordinate
                           at the end if not already closed. Default is True.
        polar (str or None): 'south' for Antarctic stereographic projection,
                             'north' for Arctic stereographic projection,
                             None to use lon/lat directly (for non-polar regions).
    
    Returns:
        array of bool: Boolean mask array. True if inside the polygon, False otherwise.
                       Length equals number of nodes (if which='node') or elements (if which='element').
    
    Example:
        >>> mask = build_regional_mask('/path/to/mesh/', '/path/to/AA_shelf.kml', which='node', polar='south')
        >>> # mask[i] is True if node i is inside the polygon
    """
    if which == 'node':
        lon, lat, node_idx, _ = read_nodes(meshpath)
    elif which == 'element':
        elements, lon, lat = read_elements(meshpath, return_coordinates=True)
    else:
        raise ValueError(f"which must be 'node' or 'element', got {which}")

    coords = read_kml_coords(kml_file, close_ring=close_ring)
    poly_lon = np.array([c[0] for c in coords])
    poly_lat = np.array([c[1] for c in coords])

    if polar is not None:
        # Use stereographic projection to handle polar regions correctly
        if polar == 'south':
            proj = Proj(proj='stere', lat_0=-90, lon_0=0, ellps='WGS84')
        elif polar == 'north':
            proj = Proj(proj='stere', lat_0=90, lon_0=0, ellps='WGS84')
        else:
            raise ValueError(f"polar must be 'south', 'north', or None, got {polar}")
        
        # Transform polygon coordinates
        poly_x, poly_y = proj(poly_lon, poly_lat)
        polygon = sh.Polygon(zip(poly_x, poly_y))
        
        # Transform mesh coordinates
        mesh_x, mesh_y = proj(lon, lat)
        points = shapely.points(mesh_x, mesh_y)
    else:
        # Use lon/lat directly for non-polar regions
        polygon = sh.Polygon(coords)
        points = shapely.points(lon, lat)

    mask = shapely.contains(polygon, points)

    return mask

def build_runoff_basin_mask(meshpath, runoff_maps_file, which='solid'):
    
    node_lon, node_lat, node_idx, node_coast = read_nodes(meshpath)

    ds_runoff = xr.open_dataset(runoff_maps_file)

    if which == 'liquid':
        basin_var = 'arrival_point_id'
    elif which == 'solid':
        basin_var = 'calving_point_id'
    else:
        raise ValueError("which must be either 'liquid' or 'solid'")

    basin_data = ds_runoff[basin_var].values  # <-- load as numpy array
    basin_lon = ds_runoff['lon'].values
    basin_lat = ds_runoff['lat'].values

    node_lon = lon_to_360(node_lon)

    # --- Compute indices vectorized ---
    # assumes regular grid!
    dlon = basin_lon[1] - basin_lon[0]
    dlat = basin_lat[1] - basin_lat[0]

    lon_idx = np.clip(((node_lon - basin_lon[0]) / dlon).round().astype(int), 0, len(basin_lon) - 1)
    lat_idx = np.clip(((node_lat - basin_lat[0]) / dlat).round().astype(int), 0, len(basin_lat) - 1)

    # --- Vectorized lookup ---
    basin_ids = basin_data[lat_idx, lon_idx]

    ds_runoff.close()

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
            'runoff_file': runoff_maps_file
        }
    )

    return ds_basin

def find_calving_front_elements(meshpath):
    """
    Finds all ocean elements that have at least one neighboring cavity element.
    
    Ocean elements are defined as elements with cavity_elvls == 1.
    Cavity elements are defined as elements with cavity_elvls > 1.
    
    Parameters:
        meshpath (str): Path to the directory containing mesh files.
    
    Returns:
        array of bool: Boolean mask of length num_elements, True for ocean elements
                       that have at least one neighboring cavity element.
    """
    elements = read_elements(meshpath)
    node_lon, node_lat, node_idx, _ = read_nodes(meshpath)
    lev_cav = read_element_levels(meshpath, which='cavity', raw=False, python_indexing=False)
    
    is_ocean = lev_cav == 1
    is_cavity = lev_cav > 1
    
    elems_of_node = build_elements_of_nodes(elements, node_idx)
    
    ocean_near_cavity = np.zeros(len(elements), dtype=bool)
    for eidx in range(len(elements)):
        if is_ocean[eidx]:
            for node in elements[eidx]:
                if any(is_cavity[e] for e in elems_of_node[node]):
                    ocean_near_cavity[eidx] = True
                    break
    
    return ocean_near_cavity


###--------> 4. Compute mesh resolution

def compute_mesh_resolution(meshpath, R=6371000.0, interpolate2regular=False):
    """
    Computes the local mesh resolution for a FESOM2 unstructured triangular mesh.
    
    Resolution is defined as the great-circle distance between neighboring nodes.
    For each unique edge in the mesh, the function returns the edge length (resolution)
    and the midpoint coordinates.
    
    Parameters:
        meshpath (str): Path to the directory containing `nod2d.out` and `elem2d.out`.
        R (float): Earth's radius in meters (default: 6371000 m).
        interpolate2regular (bool): If True, interpolate the resolution to a regular grid.
    
    Returns:
        resolution (np.ndarray): 1D array of distances between neighboring nodes (meters).
        lon_mid (np.ndarray): 1D array of edge midpoint longitudes (degrees).
        lat_mid (np.ndarray): 1D array of edge midpoint latitudes (degrees).
    """
    # Read mesh data using existing infrastructure
    node_lon, node_lat, node_idx, _ = read_nodes(meshpath)
    elements = read_elements(meshpath)
    
    # Build edges from triangles: (n1,n2), (n2,n3), (n3,n1)
    n1 = elements[:, 0]
    n2 = elements[:, 1]
    n3 = elements[:, 2]
    
    edges_a = np.column_stack([n1, n2])
    edges_b = np.column_stack([n2, n3])
    edges_c = np.column_stack([n3, n1])
    
    all_edges = np.vstack([edges_a, edges_b, edges_c])
    
    # Sort each edge so (i,j) and (j,i) become identical
    sorted_edges = np.sort(all_edges, axis=1)
    
    # Remove duplicates
    unique_edges = np.unique(sorted_edges, axis=0)
    
    # Extract node indices for each edge
    idx_i = unique_edges[:, 0]
    idx_j = unique_edges[:, 1]
    
    # Get coordinates for each node pair
    lon_i = node_lon[idx_i]
    lat_i = node_lat[idx_i]
    lon_j = node_lon[idx_j]
    lat_j = node_lat[idx_j]
    
    # Convert to radians for great-circle calculation
    lon_i_rad = np.deg2rad(lon_i)
    lat_i_rad = np.deg2rad(lat_i)
    lon_j_rad = np.deg2rad(lon_j)
    lat_j_rad = np.deg2rad(lat_j)
    
    # Haversine formula for great-circle distance
    dlat = lat_j_rad - lat_i_rad
    dlon = lon_j_rad - lon_i_rad
    
    a = np.sin(dlat / 2.0)**2 + np.cos(lat_i_rad) * np.cos(lat_j_rad) * np.sin(dlon / 2.0)**2
    c = 2.0 * np.arcsin(np.sqrt(a))
    resolution = R * c
    
    # Compute midpoint coordinates (simple arithmetic mean in degrees)
    lon_mid = (lon_i + lon_j) / 2.0
    lat_mid = (lat_i + lat_j) / 2.0

    if interpolate2regular:
        # define regular grid
        lon_grid = np.linspace(-180, 180, 1440)
        lat_grid = np.linspace(-90, 90, 720)

        lon2d, lat2d = np.meshgrid(lon_grid, lat_grid)

        # interpolate
        res_grid = griddata(
            (lon_mid, lat_mid),
            resolution,
            (lon2d, lat2d),
            method="nearest"
        )

        # create xarray Dataset
        ds = xr.Dataset(
            data_vars={
                "resolution": (("lat", "lon"), res_grid)
            },
            coords={
                "lon": lon_grid,
                "lat": lat_grid
            },
            attrs={
                "description": "Interpolated FESOM2 mesh resolution",
                "units": "meters"
            }
        )
        return ds
    else:
        # create xarray Dataset with edge-based data
        ds = xr.Dataset(
            data_vars={
                "resolution": (("edge",), resolution),
                "lon_mid": (("edge",), lon_mid),
                "lat_mid": (("edge",), lat_mid)
            },
            coords={
                "edge": np.arange(len(resolution))
            },
            attrs={
                "description": "FESOM2 mesh resolution at edge midpoints",
                "units": "meters"
            }
        )
        return ds


###--------> 5. Add mesh diagnostics to fesom.mesh.diag.nc

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

###--------> 6. Interpolation/Extrapolation
def build_spherical_nn_mapper(lon_source, lat_source,
                              lon_target, lat_target):
    """
    Build nearest-neighbor mapping between two sets of points
    on a sphere (lon/lat in degrees). The KDTree is built on the source grid.

    Parameters
    ----------
    lon_source : array-like
        Longitudes of source grid points (degrees).
    lat_source : array-like
        Latitudes of source grid points (degrees).
    lon_target : array-like
        Longitudes of target grid points (degrees).
    lat_target : array-like
        Latitudes of target grid points (degrees).

    Returns
    -------
    mapping : ndarray
        For each target point, the index of the nearest source point.
    distances : ndarray
        The spherical Euclidean distance in 3D Cartesian space.
    """
    # --- Convert degrees to radians ---
    lon1 = np.radians(lon_source)
    lat1 = np.radians(lat_source)
    lon2 = np.radians(lon_target)
    lat2 = np.radians(lat_target)

    # --- Convert spherical coordinates to Cartesian ---
    def sph2cart(lon, lat):
        x = np.cos(lat) * np.cos(lon)
        y = np.cos(lat) * np.sin(lon)
        z = np.sin(lat)
        return np.column_stack((x, y, z))

    coords1 = sph2cart(lon1, lat1)
    coords2 = sph2cart(lon2, lat2)

    # --- KDTree on source grid ---
    tree = cKDTree(coords1)

    # --- Query nearest neighbor for every target point ---
    distances, indices = tree.query(coords2, k=1)

    return indices, distances

###--------> 7. Coordinate Transformations

def unrotate_coordinates(al, be, ga, rlon, rlat):
    """
    Converts rotated coordinates to geographical coordinates.

    Parameters
    ----------
    al : float
        alpha Euler angle
    be : float
        beta Euler angle
    ga : float
        gamma Euler angle
    rlon : array
        1d array of longitudes in rotated coordinates
    rlat : array
        1d araay of latitudes in rotated coordinates

    Returns
    -------
    lon : array
        1d array of longitudes in geographical coordinates
    lat : array
        1d array of latitudes in geographical coordinates

    """
    import math as mt
    
    rad = mt.pi / 180
    al = al * rad
    be = be * rad
    ga = ga * rad
    rotate_matrix = np.zeros(shape=(3, 3))
    rotate_matrix[0, 0] = np.cos(ga) * np.cos(al) - np.sin(ga) * np.cos(be) * np.sin(al)
    rotate_matrix[0, 1] = np.cos(ga) * np.sin(al) + np.sin(ga) * np.cos(be) * np.cos(al)
    rotate_matrix[0, 2] = np.sin(ga) * np.sin(be)
    rotate_matrix[1, 0] = -np.sin(ga) * np.cos(al) - np.cos(ga) * np.cos(be) * np.sin(
        al
    )
    rotate_matrix[1, 1] = -np.sin(ga) * np.sin(al) + np.cos(ga) * np.cos(be) * np.cos(
        al
    )
    rotate_matrix[1, 2] = np.cos(ga) * np.sin(be)
    rotate_matrix[2, 0] = np.sin(be) * np.sin(al)
    rotate_matrix[2, 1] = -np.sin(be) * np.cos(al)
    rotate_matrix[2, 2] = np.cos(be)

    rotate_matrix = np.linalg.pinv(rotate_matrix)

    rlat = rlat * rad
    rlon = rlon * rad

    # Rotated Cartesian coordinates:
    xr = np.cos(rlat) * np.cos(rlon)
    yr = np.cos(rlat) * np.sin(rlon)
    zr = np.sin(rlat)

    # Geographical Cartesian coordinates:
    xg = rotate_matrix[0, 0] * xr + rotate_matrix[0, 1] * yr + rotate_matrix[0, 2] * zr
    yg = rotate_matrix[1, 0] * xr + rotate_matrix[1, 1] * yr + rotate_matrix[1, 2] * zr
    zg = (
        rotate_matrix[2, 0] * xr + rotate_matrix[2, 1] * yr + rotate_matrix[2, 2] * zr
    )

    # Geographical coordinates:
    lat = np.arcsin(zg)
    lon = np.arctan2(yg, xg)

    a = np.where((np.abs(xg) + np.abs(yg)) == 0)
    if a:
        lon[a] = 0

    lat = lat / rad
    lon = lon / rad

    return lon, lat

###--------> 8. Build Land Sea Mask

def build_land_sea_mask(meshpath, nlon=1440, nlat=720, has_cavity=True, cavity_is_land=True):
    """
    Build a regular lon-lat land-sea mask from a FESOM2 unstructured mesh.
    
    The mesh only contains ocean cells (land/continents are holes). Grid points
    falling inside mesh triangles are marked as ocean, others as land.
    
    Parameters
    ----------
    meshpath : str
        Path to mesh directory containing nod2d.out, elem2d.out, and cavity_elvls.out.
    nlon : int
        Number of longitude grid points (default: 1440 for 0.25° resolution).
    nlat : int
        Number of latitude grid points (default: 720 for 0.25° resolution).
    has_cavity : bool
        If True, the mesh contains cavity information (cavity_elvls.out) and 
        cavity_is_land is used to decide how to treat cavity cells.
        If False, all FESOM grid cells are treated as ocean.
    cavity_is_land : bool
        Only used when has_cavity=True. If True, sub-ice cavity cells 
        (cavity_elvls > 1) are treated as land.
    
    Returns
    -------
    xarray.Dataset
        Dataset with 'mask' variable (1 = ocean, 0 = land) on regular lon-lat grid.
    """
    # Read mesh data
    node_lon, node_lat, _, _ = read_nodes(meshpath)
    elements = read_elements(meshpath)
    
    # Handle cavity information
    if has_cavity:
        # Read cavity levels (values > 1 indicate cavity elements)
        cavity_elvls = read_element_levels(meshpath, which='cavity', raw=False, python_indexing=False)
        is_cavity = cavity_elvls > 1
    else:
        # No cavity info: treat all elements as non-cavity (ocean)
        is_cavity = np.zeros(len(elements), dtype=bool)
    
    # Create regular grid
    lon_grid = np.linspace(-180, 180, nlon, endpoint=False)
    lat_grid = np.linspace(-90, 90, nlat, endpoint=False)
    dlon = lon_grid[1] - lon_grid[0]
    dlat = lat_grid[1] - lat_grid[0]
    # Center grid cells
    lon_grid = lon_grid + dlon / 2
    lat_grid = lat_grid + dlat / 2
    
    # Initialize mask as land (0)
    mask = np.zeros((nlat, nlon), dtype=np.int8)
    
    # Helper function: check if point is inside triangle using barycentric coordinates
    def point_in_triangle(px, py, x1, y1, x2, y2, x3, y3):
        denom = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)
        if abs(denom) < 1e-12:
            return False
        a = ((y2 - y3) * (px - x3) + (x3 - x2) * (py - y3)) / denom
        b = ((y3 - y1) * (px - x3) + (x1 - x3) * (py - y3)) / denom
        c = 1 - a - b
        return (a >= 0) and (b >= 0) and (c >= 0)
    
    # Process each element
    for eidx, (n1, n2, n3) in enumerate(elements):
        # Skip cavity elements if they should be treated as land
        if cavity_is_land and is_cavity[eidx]:
            continue
        
        # Get triangle vertices
        lons = np.array([node_lon[n1], node_lon[n2], node_lon[n3]])
        lats = np.array([node_lat[n1], node_lat[n2], node_lat[n3]])
        
        # Handle triangles crossing the antimeridian (lon jump > 180°)
        lon_range = lons.max() - lons.min()
        if lon_range > 180:
            # Shift longitudes to 0-360 range for this triangle
            lons = np.where(lons < 0, lons + 360, lons)
            lon_min, lon_max = lons.min(), lons.max()
            lat_min, lat_max = lats.min(), lats.max()
            
            # Find grid cells that could intersect this triangle
            i_min = max(0, int((lat_min - lat_grid[0]) / dlat) - 1)
            i_max = min(nlat - 1, int((lat_max - lat_grid[0]) / dlat) + 1)
            
            for i in range(i_min, i_max + 1):
                lat_pt = lat_grid[i]
                for j in range(nlon):
                    lon_pt = lon_grid[j]
                    # Shift grid lon to 0-360 for comparison
                    lon_pt_shifted = lon_pt if lon_pt >= 0 else lon_pt + 360
                    if lon_pt_shifted < lon_min - dlon or lon_pt_shifted > lon_max + dlon:
                        continue
                    if point_in_triangle(lon_pt_shifted, lat_pt, 
                                         lons[0], lats[0], lons[1], lats[1], lons[2], lats[2]):
                        mask[i, j] = 1
        else:
            # Normal triangle
            lon_min, lon_max = lons.min(), lons.max()
            lat_min, lat_max = lats.min(), lats.max()
            
            # Find grid cell index range
            j_min = max(0, int((lon_min - lon_grid[0]) / dlon) - 1)
            j_max = min(nlon - 1, int((lon_max - lon_grid[0]) / dlon) + 1)
            i_min = max(0, int((lat_min - lat_grid[0]) / dlat) - 1)
            i_max = min(nlat - 1, int((lat_max - lat_grid[0]) / dlat) + 1)
            
            # Check grid cells within bounding box
            for i in range(i_min, i_max + 1):
                lat_pt = lat_grid[i]
                for j in range(j_min, j_max + 1):
                    lon_pt = lon_grid[j]
                    if point_in_triangle(lon_pt, lat_pt, 
                                         lons[0], lats[0], lons[1], lats[1], lons[2], lats[2]):
                        mask[i, j] = 1
    
    # Create xarray Dataset
    ds = xr.Dataset(
        data_vars={
            'mask': (['lat', 'lon'], mask)
        },
        coords={
            'lon': lon_grid,
            'lat': lat_grid
        },
        attrs={
            'description': 'Land-sea mask from FESOM2 mesh',
            'mask_convention': '1 = ocean, 0 = land',
            'has_cavity': str(has_cavity),
            'cavity_is_land': str(cavity_is_land) if has_cavity else 'N/A',
            'mesh_path': meshpath
        }
    )
    
    ds['mask'].attrs = {
        'long_name': 'land-sea mask',
        'units': '1',
        'flag_values': [0, 1],
        'flag_meanings': 'land ocean'
    }
    
    return ds



def mesh2vtk(meshpath, dstfile, which='seafloor'):
    """
    Converts a FESOM2 mesh defined by `nod2d.out` and `elem2d.out` files into VTK format.
    
    Parameters:
        meshpath (str): Path to the directory containing `nod2d.out` and `elem2d.out`.
        dstfile (str): File name of the output vtk file (without .vtk)
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

    write_vtk(f"{dstfile}.vtk", nodes, elements, cell_data)
    return
