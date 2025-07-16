# so_ase.eval_sea_ice.py

# To Do:
# Sea Ice Production Maps
#
#
#
#
#

import xarray as xr
import numpy as np
from .helpers_mesh import find_nodes_in_box


def fesom_ice_area(
    src_path,
    mesh_diag_path,
    years=(1979, 2015),
    box=[-180, 180, -90, -60],
    aice_threshhold=0.15,
    log=True,
):
    """
    Compute total sea ice area within a specified geographic bounding box using FESOM2 output.

    This function loads sea ice concentration (`a_ice`) from FESOM2 output files over a range
    of years and calculates the sea ice area by summing the nodal areas where sea ice
    concentration exceeds a given threshold, within a user-defined geographic region.

    Parameters
    ----------
    src_path : str
        Path to the directory containing FESOM2 output NetCDF files named
        `a_ice.fesom.{year}.nc`.
    mesh_diag_path : str
        Path to the directory containing the mesh diagnostic file `fesom.mesh.diag.nc`,
        which provides nodal coordinates and areas.
    years : tuple of int, optional
        Start and end year (exclusive) for the time series. Defaults to (1979, 2015).
    box : list of float, optional
        Geographic bounding box as [lon_min, lon_max, lat_min, lat_max] to define the region
        of interest. Defaults to `[-180, 180, -90, -60]` (entire Southern Hemisphere).
    aice_threshhold : float, optional
        Minimum sea ice concentration (0–1) above which a node is considered ice-covered.
        Defaults to 0.15.
    log : bool, optional
        If True, prints progress and file loading information to standard output.

    Returns
    -------
    xarray.Dataset
        Dataset containing the time series variable `sea_ice_area` (in square meters)
        representing total sea ice area within the defined region. The variables `a_ice`
        and the intermediate `ice_mask` are removed from the returned dataset.

    Notes
    -----
    A binary sea ice mask is created where `a_ice > aice_threshhold`. The total sea ice area
    is computed by summing the nodal areas (`nod_area`) over the masked nodes.
    """

    # Load mesh dignostic file
    mesh_diag = xr.open_dataset(f"{mesh_diag_path}fesom.mesh.diag.nc")
    if log:
        print("Mesh diagnostics loaded:", flush=True)
        print(f"{mesh_diag_path}fesom.mesh.diag.nc", flush=True)

    # Find indices of nodes within the specified box
    inds = find_nodes_in_box(mesh_diag_path, box=box, log=log)
      
    # Load files for sea ice concentration from src_path
    files2load = [f"{src_path}a_ice.fesom.{y}.nc" for y in range(years[0], years[1])]

    result = []
    
    for file in files2load:
        ds = xr.open_dataset(file).load()
        if log:
            print(f"File loaded: {file}", flush=True)
        
        # Crop datasets
        ds_cropped = ds.isel(nod2=inds)
        mesh_diag_cropped = mesh_diag.isel(nod2=inds)
    
        # Create sea ice mask
        ds_cropped["ice_mask"] = (
            ("time", "nod2"),
            np.where(ds_cropped.a_ice > aice_threshhold, True, False),
        )
    
        # Sum over non-masked nodal areas
        ds_cropped["sea_ice_area"] = (
            ("time"),
            (ds_cropped.ice_mask * mesh_diag_cropped.nod_area.isel(nz=0))
            .sum(dim="nod2")
            .values,
        )
    
        # Prepare output
        ds_cropped = ds_cropped.drop_vars(["a_ice", "ice_mask"])
        ds_cropped.sea_ice_area.attrs["units"] = "$m^2$"
        ds_cropped.sea_ice_area.attrs["long_name"] = "total sea ice area"
        ds_cropped.sea_ice_area.attrs["bounding box"] = (
            f"Longitude: {box[0]}E to {box[1]}E, Latitude: {box[2]}N to {box[3]}N"
        )

        result.append(ds_cropped)

    result = xr.concat(result, dim='time')

    return result


def fesom_ice_volume(
    src_path, 
    mesh_diag_path, 
    years=(1979, 2015), 
    box=[-180, 180, -90, -60], 
    log=True
):
    """
    Compute total sea ice volume within a specified geographic bounding box using FESOM2 output.

    This function loads sea ice concentration (`a_ice`) and thickness (`m_ice`) variables
    from FESOM2 output files over a range of years, and calculates the sea ice volume by
    summing the product of sea ice thickness, concentration, and nodal area over the nodes
    that fall within a user-defined geographic region.

    Parameters
    ----------
    src_path : str
        Path to the directory containing FESOM2 output NetCDF files named
        `a_ice.fesom.{year}.nc` and `m_ice.fesom.{year}.nc`.
    mesh_diag_path : str
        Path to the directory containing the mesh diagnostic file `fesom.mesh.diag.nc`,
        which provides nodal coordinates and areas.
    years : tuple of int, optional
        Start and end year (exclusive) for the time series. Defaults to (1979, 2015).
    box : list of float, optional
        Geographic bounding box as [lon_min, lon_max, lat_min, lat_max] to define the region
        of interest. Defaults to `[-180, 180, -90, -60]` (entire Southern Hemisphere).
    log : bool, optional
        If True, prints progress and file loading information to standard output.

    Returns
    -------
    xarray.Dataset
        Dataset containing the time series variable `sea_ice_volume` (in cubic meters)
        representing total sea ice volume within the defined region. The original `m_ice`
        variable is dropped from the returned dataset.

    Notes
    -----
    Sea ice volume is calculated as:

        sea_ice_volume = sum_over_nodes(area * m_ice * a_ice)
    """

    # Load mesh dignostic file
    mesh_diag = xr.open_dataset(f"{mesh_diag_path}fesom.mesh.diag.nc")
    if log:
        print("Mesh diagnostics loaded:", flush=True)
        print(f"{mesh_diag_path}fesom.mesh.diag.nc", flush=True)

    # Find indices of nodes within the specified box
    inds = find_nodes_in_box(mesh_diag_path, box=box, log=log)

    files2load_aice = [f"{src_path}a_ice.fesom.{y}.nc" for y in range(years[0], years[1])]
    files2load_mice = [f"{src_path}m_ice.fesom.{y}.nc" for y in range(years[0], years[1])]

    result = []
    
    for file1, file2 in zip(files2load_aice, files2load_mice):
        # Load files for sea ice concentration from src_path
        ds = xr.open_mfdataset([file1, file2]).load()
        if log:
            print(f"File loaded: {file1}", flush=True)
            print(f"File loaded: {file2}", flush=True)
            
        # Crop datasets
        ds_cropped = ds.isel(nod2=inds)
        mesh_diag_cropped = mesh_diag.isel(nod2=inds, nz1=0)
    
        # Compute total ice volume (nodal area * nodal ice height * nodal ice concentration)
        sea_ice_volume = (
            mesh_diag_cropped.nod_area
            * ds_cropped.m_ice
            * ds_cropped.a_ice
        ).sum(dim="nod2")
    
        # Sum over nodal areas
        ds_cropped["sea_ice_volume"] = sea_ice_volume
    
        # Prepare output
        ds_cropped = ds_cropped.drop_vars(["m_ice", "a_ice"])
        ds_cropped.sea_ice_volume.attrs["units"] = "$m^3$"
        ds_cropped.sea_ice_volume.attrs["long_name"] = "total sea ice volume"
        ds_cropped.sea_ice_volume.attrs["bounding box"] = (
            f"Longitude: {box[0]}E to {box[1]}E, Latitude: {box[2]}N to {box[3]}N"
        )

        result.append(ds_cropped)

    result = xr.concat(result, dim='time')

    return result
