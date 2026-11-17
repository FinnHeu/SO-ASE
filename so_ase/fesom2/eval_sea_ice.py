# so_ase/fesom2/eval_sea_ice.py

import xarray as xr
import numpy as np
import os
from .helpers_mesh import find_nodes_in_box


def _split_box(box):
    if isinstance(box, dict):
        return box["box"], box.get("name")
    return box, None


def fesom_sea_ice_area(
    src_path,
    mesh_diag_path,
    year,
    box=[-180, 180, -90, -40],
    siconc_threshold=0.0,
    savepath='./',
    log=True
):
    """
    Compute and save total sea ice area time series from FESOM2 output within a specified
    geographic bounding box.

    This function processes the yearly FESOM2 sea ice concentration file
    (`a_ice.fesom.<year>.nc`), computes the total sea ice area by summing the
    nodal areas where sea ice concentration exceeds a given threshold, and
    writes the result to a NetCDF file for that year. If the output file already
    exists, the function returns immediately without loading any data.

    The sea ice area is calculated only for mesh nodes that fall within the
    user-defined geographic bounding box. Spatial masking is based on the
    FESOM mesh diagnostic file (`fesom.mesh.diag.nc`).

    Output filenames are self-describing and include:
      - the processed year,
      - the geographic bounding box formatted using N/S/E/W notation.

    Parameters
    ----------
    src_path : str
        Path to the directory containing FESOM2 sea ice concentration files
        named `a_ice.fesom.<year>.nc`.
    mesh_diag_path : str
        Path to the directory containing the FESOM mesh diagnostic file
        `fesom.mesh.diag.nc`.
    year : int
        Year to process.
    box : list of float or dict, optional
        Geographic bounding box specified as
        [lon_min, lon_max, lat_min, lat_max] in degrees.
        Longitudes are expected in degrees east, latitudes in degrees north.
        Alternatively, a dict with key "box" (the bounding box above) and optional key "name",
        which is used as the box part of the output filenames instead of the N/S/E/W string.
        Default is [-180, 180, -90, -60].
    siconc_threshold : float, optional
        Minimum sea ice concentration (range 0–1) required for a node to be
        considered ice-covered. Default is 0.15.
    savepath : str
        Directory where the output NetCDF files will be written. The directory
        is created if it does not already exist.
    log : bool, optional
        If True, print progress messages, including file loading, skipping,
        and saving information. Default is True.

    Returns
    -------
    Sea ice area dataset for `year`, or None if the output file already exists
        
    """

    def format_lat(lat):
        return f"{abs(lat)}{'S' if lat < 0 else 'N'}"

    def format_lon(lon):
        return f"{abs(lon)}{'W' if lon < 0 else 'E'}"

    # Prepare filename-safe strings
    box, name = _split_box(box)
    box_str = name or (
        f"{format_lon(box[0])}_{format_lon(box[1])}_"
        f"{format_lat(box[2])}_{format_lat(box[3])}"
    )

    out_file = f"{savepath}/sea_ice_area_{year}_{box_str}.nc"
    if os.path.exists(out_file):
        if log:
            print(f"Skipping existing file: {out_file}", flush=True)
        return

    os.makedirs(savepath, exist_ok=True)

    # Load mesh diagnostic file once
    mesh_diag = xr.open_dataset(f"{mesh_diag_path}fesom.mesh.diag.nc")
    if log:
        print("Mesh diagnostics loaded:", flush=True)
        print(f"{mesh_diag_path}fesom.mesh.diag.nc", flush=True)

    # Find indices of nodes within the specified box
    inds = find_nodes_in_box(mesh_diag_path, box=box, log=log)

    # Open file with cftime decoder
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
    ds = xr.open_mfdataset([f"{src_path}a_ice.fesom.{year}.nc"], decode_times=time_coder, chunks={"time": 12})
    if 'nz' in ds.dims:
        ds = ds.isel(nz=0)
        ds = ds.drop_vars('nz')

    # Crop datasets
    ds_cropped = ds.isel(nod2=inds)
    mesh_diag_cropped = mesh_diag.isel(nod2=inds)

    # Create sea ice mask
    ice_mask = ds_cropped.a_ice > siconc_threshold

    # Compute sea ice area time series
    sea_ice_area = (
        ds_cropped.a_ice * ice_mask * mesh_diag_cropped.nod_area.isel(nz=0)
    ).sum(dim="nod2")

    # Create output dataset
    ds_out = sea_ice_area.to_dataset(name="sea_ice_area")

    # Variable metadata
    ds_out.sea_ice_area.attrs["units"] = "m^2"
    ds_out.sea_ice_area.attrs["long_name"] = "total sea ice area"
    ds_out.sea_ice_area.attrs["bounding_box"] = (
        f"Longitude: {box[0]}E to {box[1]}E, "
        f"Latitude: {box[2]}N to {box[3]}N"
    )

    # Global metadata
    ds_out.attrs["source"] = "FESOM2"
    ds_out.attrs["siconc_threshold"] = siconc_threshold

    # Select the monthly data of this year and save
    sia_year = ds_out.sel(time=ds_out.time.dt.year == year)
    sia_year.to_netcdf(out_file)

    if log:
        print(f"Saved: {out_file}", flush=True)
        print("All done!", flush=True)
    return ds_out

def fesom_sea_ice_volume(
    src_path, 
    mesh_diag_path,
    year,
    box=[-180, 180, -90, -40], 
    savepath='./',
    log=True
):
    """
    Compute and save total sea ice volume time series from FESOM2 output within a specified
    geographic bounding box.

    This function loads sea ice concentration (`a_ice`) and thickness (`m_ice`) variables
    from the FESOM2 output files of a single year, and calculates the sea ice volume by
    summing the product of sea ice thickness, concentration, and nodal area over the nodes
    that fall within a user-defined geographic region. One output file is written per
    year. If the output file already exists, the function returns immediately without loading any data.

    Parameters
    ----------
    src_path : str
        Path to the directory containing FESOM2 output NetCDF files named
        `a_ice.fesom.{year}.nc` and `m_ice.fesom.{year}.nc`.
    mesh_diag_path : str
        Path to the directory containing the mesh diagnostic file `fesom.mesh.diag.nc`,
        which provides nodal coordinates and areas.
    year : int
        Year to process.
    box : list of float or dict, optional
        Geographic bounding box as [lon_min, lon_max, lat_min, lat_max] to define the region
        of interest. Alternatively, a dict with key "box" (the bounding box above) and optional key "name",
        which is used as the box part of the output filenames instead of the N/S/E/W string.
        Defaults to `[-180, 180, -90, -60]` (entire Southern Hemisphere).
    savepath : str
        Directory where the output NetCDF files will be written. The directory
        is created if it does not already exist.
    log : bool, optional
        If True, prints progress and file loading information to standard output.

    Returns
    -------
    Sea ice volume dataset for `year`, or None if the output file already exists

    Notes
    -----
    Sea ice volume is calculated as:

        sea_ice_volume = sum_over_nodes(area * m_ice * a_ice)
    """

    def format_lat(lat):
        return f"{abs(lat)}{'S' if lat < 0 else 'N'}"

    def format_lon(lon):
        return f"{abs(lon)}{'W' if lon < 0 else 'E'}"

    # Prepare filename-safe strings
    box, name = _split_box(box)
    box_str = name or (
        f"{format_lon(box[0])}_{format_lon(box[1])}_"
        f"{format_lat(box[2])}_{format_lat(box[3])}"
    )

    out_file = f"{savepath}/sea_ice_volume_{year}_{box_str}.nc"
    if os.path.exists(out_file):
        if log:
            print(f"Skipping existing file: {out_file}", flush=True)
        return

    os.makedirs(savepath, exist_ok=True)

    # Load mesh diagnostic file
    mesh_diag = xr.open_dataset(f"{mesh_diag_path}fesom.mesh.diag.nc")
    if log:
        print("Mesh diagnostics loaded:", flush=True)
        print(f"{mesh_diag_path}fesom.mesh.diag.nc", flush=True)

    # Find indices of nodes within the specified box
    inds = find_nodes_in_box(mesh_diag_path, box=box, log=log)

    # Open files with cftime decoder
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
    files2load = [f"{src_path}a_ice.fesom.{year}.nc", f"{src_path}m_ice.fesom.{year}.nc"]
    ds = xr.open_mfdataset(files2load, decode_times=time_coder, chunks={"time": 12})
    if 'nz' in ds.dims:
        ds = ds.isel(nz=0)
        ds = ds.drop_vars('nz')

    # Crop datasets
    ds_cropped = ds.isel(nod2=inds)
    mesh_diag_cropped = mesh_diag.isel(nod2=inds, nz=0)

    # Compute total ice volume time series (nodal area * nodal ice height * nodal ice concentration)
    sea_ice_volume = (
        mesh_diag_cropped.nod_area
        * ds_cropped.m_ice
        * ds_cropped.a_ice
    ).sum(dim="nod2")

    # Create output dataset
    ds_out = sea_ice_volume.to_dataset(name="sea_ice_volume")

    # Variable metadata
    ds_out.sea_ice_volume.attrs["units"] = "m^3"
    ds_out.sea_ice_volume.attrs["long_name"] = "total sea ice volume"
    ds_out.sea_ice_volume.attrs["bounding_box"] = (
        f"Longitude: {box[0]}E to {box[1]}E, "
        f"Latitude: {box[2]}N to {box[3]}N"
    )

    # Global metadata
    ds_out.attrs["source"] = "FESOM2"

    # Select the monthly data of this year and save
    siv_year = ds_out.sel(time=ds_out.time.dt.year == year)
    siv_year.to_netcdf(out_file)

    if log:
        print(f"Saved: {out_file}", flush=True)
        print("All done!", flush=True)

    return ds_out

def fesom_sea_ice_extent(
    src_path,
    mesh_diag_path,
    year,
    box=[-180, 180, -90, -40],
    siconc_threshold=0.15,
    savepath='./',
    log=True
):
    """
    Compute and save total sea ice extent time series from FESOM2 output within a specified
    geographic bounding box.

    This function processes the yearly FESOM2 sea ice concentration file
    (`a_ice.fesom.<year>.nc`), computes the total sea ice extent by summing the
    nodal areas where sea ice concentration exceeds a given threshold, and
    writes the result to a NetCDF file for that year. If the output file already
    exists, the function returns immediately without loading any data.

    The sea ice extent is calculated only for mesh nodes that fall within the
    user-defined geographic bounding box. Spatial masking is based on the
    FESOM mesh diagnostic file (`fesom.mesh.diag.nc`).

    Output filenames are self-describing and include:
      - the processed year,
      - the geographic bounding box formatted using N/S/E/W notation.

    Parameters
    ----------
    src_path : str
        Path to the directory containing FESOM2 sea ice concentration files
        named `a_ice.fesom.<year>.nc`.
    mesh_diag_path : str
        Path to the directory containing the FESOM mesh diagnostic file
        `fesom.mesh.diag.nc`.
    year : int
        Year to process.
    box : list of float or dict, optional
        Geographic bounding box specified as
        [lon_min, lon_max, lat_min, lat_max] in degrees.
        Longitudes are expected in degrees east, latitudes in degrees north.
        Alternatively, a dict with key "box" (the bounding box above) and optional key "name",
        which is used as the box part of the output filenames instead of the N/S/E/W string.
        Default is [-180, 180, -90, -60].
    siconc_threshold : float, optional
        Minimum sea ice concentration (range 0–1) required for a node to be
        considered ice-covered. Default is 0.15.
    savepath : str
        Directory where the output NetCDF files will be written. The directory
        is created if it does not already exist.
    log : bool, optional
        If True, print progress messages, including file loading, skipping,
        and saving information. Default is True.

    Returns
    -------
    Sea ice extent dataset for `year`, or None if the output file already exists

    """

    def format_lat(lat):
        return f"{abs(lat)}{'S' if lat < 0 else 'N'}"

    def format_lon(lon):
        return f"{abs(lon)}{'W' if lon < 0 else 'E'}"

    # Prepare filename-safe strings
    box, name = _split_box(box)
    box_str = name or (
        f"{format_lon(box[0])}_{format_lon(box[1])}_"
        f"{format_lat(box[2])}_{format_lat(box[3])}"
    )

    out_file = f"{savepath}/sea_ice_extent_{year}_{box_str}.nc"
    if os.path.exists(out_file):
        if log:
            print(f"Skipping existing file: {out_file}", flush=True)
        return

    os.makedirs(savepath, exist_ok=True)

    # Load mesh diagnostic file once
    mesh_diag = xr.open_dataset(f"{mesh_diag_path}fesom.mesh.diag.nc")
    if log:
        print("Mesh diagnostics loaded:", flush=True)
        print(f"{mesh_diag_path}fesom.mesh.diag.nc", flush=True)

    # Find indices of nodes within the specified box
    inds = find_nodes_in_box(mesh_diag_path, box=box, log=log)

    # Open file with cftime decoder
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
    ds = xr.open_mfdataset([f"{src_path}a_ice.fesom.{year}.nc"], decode_times=time_coder, chunks={"time": 12})
    if 'nz' in ds.dims:
        ds = ds.isel(nz=0)
        ds = ds.drop_vars('nz')

    # Crop datasets
    ds_cropped = ds.isel(nod2=inds)
    mesh_diag_cropped = mesh_diag.isel(nod2=inds)

    # Create sea ice mask
    ice_mask = ds_cropped.a_ice > siconc_threshold

    # Compute sea ice extent time series
    sea_ice_extent = (
        ice_mask * mesh_diag_cropped.nod_area.isel(nz=0)
    ).sum(dim="nod2")

    # Create output dataset
    ds_out = sea_ice_extent.to_dataset(name="sea_ice_extent")

    # Variable metadata
    ds_out.sea_ice_extent.attrs["units"] = "m^2"
    ds_out.sea_ice_extent.attrs["long_name"] = "total sea ice extent"
    ds_out.sea_ice_extent.attrs["bounding_box"] = (
        f"Longitude: {box[0]}E to {box[1]}E, "
        f"Latitude: {box[2]}N to {box[3]}N"
    )

    # Global metadata
    ds_out.attrs["source"] = "FESOM2"
    ds_out.attrs["siconc_threshold"] = siconc_threshold

    # Select the monthly data of this year and save
    sie_year = ds_out.sel(time=ds_out.time.dt.year == year)
    sie_year.to_netcdf(out_file)

    if log:
        print(f"Saved: {out_file}", flush=True)
        print("All done!", flush=True)

    return ds_out
