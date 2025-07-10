# so_ase.helpers_mesh.py

def read_nodes(meshpath):
    """
    Reads 2D node coordinates from a `nod2d.out` file in a given mesh path.

    Parameters:
        meshpath (str): Path to the directory containing `nod2d.out`.

    Returns:
        list of tuple: A list of (longitude, latitude) tuples for each node.
    """
    with open(f'{meshpath}nod2d.out', 'r') as f:
        num_nodes = int(f.readline())
        nodes = []
        for _ in range(num_nodes):
            parts = f.readline().split()
            node_id = int(parts[0])  # unused
            lon = float(parts[1])
            lat = float(parts[2])
            nodes.append((lon, lat))
    return nodes


def read_elements(meshpath):
    """
    Reads element connectivity information from a `elem2d.out` file.

    Parameters:
        meshpath (str): Path to the directory containing `elem2d.out`.

    Returns:
        list of tuple: A list of elements, where each element is a tuple
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
    return elements


def read_aux3d(meshpath):
    """
    Reads vertical level information (e.g., depths) from an `aux3d.out` file,
    skipping the first `num_levels` lines that typically contain header or
    level-wise data, and returns one depth value per node.

    Parameters:
        meshpath (str): Path to the directory containing `aux3d.out` and `nod2d.out`.

    Returns:
        list of float: A list of depth values, one for each node.
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
    return depths