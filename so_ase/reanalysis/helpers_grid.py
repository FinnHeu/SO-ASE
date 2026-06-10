# so_ase/reanalysis/helpers_grid.py

import xarray as xr
import numpy as np
import math

from pyproj import Proj, Transformer

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