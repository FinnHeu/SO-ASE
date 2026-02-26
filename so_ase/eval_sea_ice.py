# so_ase.eval_sea_ice.py

import xarray as xr
import numpy as np
from .helpers_mesh import find_nodes_in_box, gridcell_area_hadley, reproject_to_latlon


import os
import numpy as np
import xarray as xr

def fesom_sea_ice_area(
    src_path,
    mesh_diag_path,
    years=(1979, 2025),
    box=[-180, 180, -90, -60],
    siconc_threshold=0.15,
    savepath='./',
    log=True
):
    """
    Compute and save total sea ice area time series from FESOM2 output within a specified
    geographic bounding box.

    This function processes yearly FESOM2 sea ice concentration files
    (`a_ice.fesom.<year>.nc`), computes the total sea ice area by summing the
    nodal areas where sea ice concentration exceeds a given threshold, and
    writes the result to disk as NetCDF files. One output file is written per
    year containing the full time series for that year. If an output file already exists, that year is skipped.

    The sea ice area is calculated only for mesh nodes that fall within the
    user-defined geographic bounding box. Spatial masking is based on the
    FESOM mesh diagnostic file (`fesom.mesh.diag.nc`), which is loaded once
    and reused for all years.

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
    years : tuple of int, optional
        Start and end year (end year exclusive) defining the range of years
        to process. Default is (1979, 2025).
    box : list of float, optional
        Geographic bounding box specified as
        [lon_min, lon_max, lat_min, lat_max] in degrees.
        Longitudes are expected in degrees east, latitudes in degrees north.
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
    None
        The function does not return any objects. Results are written directly
        to disk as NetCDF files.
    """

    def format_lat(lat):
        return f"{abs(lat)}{'S' if lat < 0 else 'N'}"

    def format_lon(lon):
        return f"{abs(lon)}{'W' if lon < 0 else 'E'}"
        
    os.makedirs(savepath, exist_ok=True)

    # Load mesh diagnostic file once
    mesh_diag = xr.open_dataset(f"{mesh_diag_path}fesom.mesh.diag.nc")
    if log:
        print("Mesh diagnostics loaded:", flush=True)
        print(f"{mesh_diag_path}fesom.mesh.diag.nc", flush=True)

    # Find indices of nodes within the specified box
    inds = find_nodes_in_box(mesh_diag_path, box=box, log=log)

    # Prepare filename-safe strings
    box_str = (
        f"{format_lon(box[0])}_{format_lon(box[1])}_"
        f"{format_lat(box[2])}_{format_lat(box[3])}"
    )

    for year in range(years[0], years[1]):

        out_file = f"{savepath}/sea_ice_area_{year}_{box_str}.nc"

        if os.path.exists(out_file):
            if log:
                print(f"Skipping existing file: {out_file}", flush=True)
            continue

        in_file = f"{src_path}a_ice.fesom.{year}.nc"
        
        # Open files with cftime decoder
        time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
        ds = xr.open_dataset(in_file, decode_times=time_coder).load()

        # Crop datasets
        ds_cropped = ds.isel(nod2=inds)
        mesh_diag_cropped = mesh_diag.isel(nod2=inds)

        # Create sea ice mask
        ice_mask = ds_cropped.a_ice > siconc_threshold

        # Compute sea ice area time series
        sea_ice_area = (
            ice_mask * mesh_diag_cropped.nod_area.isel(nz=0)
        ).sum(dim="nod2")

        # Create output dataset with full time series
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

        # Save to disk
        ds_out.to_netcdf(out_file)

        if log:
            print(f"Saved: {out_file}", flush=True)

    if log:
        print("All done!", flush=True)

def fesom_ice_volume(
    src_path, 
    mesh_diag_path,
    years=(1979, 2015), 
    box=[-180, 180, -90, -60], 
    savepath='./',
    log=True
):
    """
    Compute and save total sea ice volume time series from FESOM2 output within a specified
    geographic bounding box.

    This function loads sea ice concentration (`a_ice`) and thickness (`m_ice`) variables
    from FESOM2 output files over a range of years, and calculates the sea ice volume by
    summing the product of sea ice thickness, concentration, and nodal area over the nodes
    that fall within a user-defined geographic region. One output file is written per
    year containing the full time series for that year. If an output file already exists, that year is skipped.

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
    savepath : str
        Directory where the output NetCDF files will be written. The directory
        is created if it does not already exist.
    log : bool, optional
        If True, prints progress and file loading information to standard output.

    Returns
    -------
    None
        The function does not return any objects. Results are written directly
        to disk as NetCDF files.

    Notes
    -----
    Sea ice volume is calculated as:

        sea_ice_volume = sum_over_nodes(area * m_ice * a_ice)
    """

    def format_lat(lat):
        return f"{abs(lat)}{'S' if lat < 0 else 'N'}"

    def format_lon(lon):
        return f"{abs(lon)}{'W' if lon < 0 else 'E'}"
        
    os.makedirs(savepath, exist_ok=True)

    # Load mesh diagnostic file
    mesh_diag = xr.open_dataset(f"{mesh_diag_path}fesom.mesh.diag.nc")
    if log:
        print("Mesh diagnostics loaded:", flush=True)
        print(f"{mesh_diag_path}fesom.mesh.diag.nc", flush=True)

    # Find indices of nodes within the specified box
    inds = find_nodes_in_box(mesh_diag_path, box=box, log=log)

    # Prepare filename-safe strings
    box_str = (
        f"{format_lon(box[0])}_{format_lon(box[1])}_"
        f"{format_lat(box[2])}_{format_lat(box[3])}"
    )

    for year in range(years[0], years[1]):

        out_file = f"{savepath}/sea_ice_volume_{year}_{box_str}.nc"

        if os.path.exists(out_file):
            if log:
                print(f"Skipping existing file: {out_file}", flush=True)
            continue

        file_aice = f"{src_path}a_ice.fesom.{year}.nc"
        file_mice = f"{src_path}m_ice.fesom.{year}.nc"
            
        # Load files for sea ice concentration and thickness
        time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
        ds = xr.open_mfdataset([file_aice, file_mice], decode_times=time_coder).load()
        
        # Crop datasets
        ds_cropped = ds.isel(nod2=inds)
        mesh_diag_cropped = mesh_diag.isel(nod2=inds, nz1=0)
    
        # Compute total ice volume time series (nodal area * nodal ice height * nodal ice concentration)
        sea_ice_volume = (
            mesh_diag_cropped.nod_area
            * ds_cropped.m_ice
            * ds_cropped.a_ice
        ).sum(dim="nod2")
    
        # Create output dataset with full time series
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

        # Save to disk
        ds_out.to_netcdf(out_file)

        if log:
            print(f"Saved: {out_file}", flush=True)

    if log:
        print("All done!", flush=True)

def nsidc_ice_area(src_path,
    years=(1979, 2015),
    box=[-180, 180, -90, -50],
    siconc_threshold=0.15,
    grouping='annual.mean',
    log=True):

    """
    Compute total sea ice area from NSIDC CDR data within a geographic region.

    This function loads NSIDC sea ice concentration data, reprojects it from 
    polar stereographic coordinates to regular lat/lon, and calculates the 
    total sea ice area where concentration exceeds a given threshold.

    Parameters
    ----------
    src_path : str
        Path to directory containing SIC NetCDF files named as 'sic.<year>.nc'.
    years : tuple of int, optional
        Start and end year (exclusive) for processing, e.g., (1979, 2015).
    box : list of float, optional
        Geographic bounding box [lon_min, lon_max, lat_min, lat_max] for area calculation.
    siconc_threshold : float, optional
        Sea ice concentration threshold above which a grid cell is counted as ice-covered (default is 0.15).
    log : bool, optional
        If True, print progress messages during processing.

    Returns
    -------
    xarray.DataArray
        Time series of total sea ice area (in km²) per time step within the specified region.

    Notes
    -----
    - Assumes a constant 25 km x 25 km grid cell size.

    """

    files2load = [f"{src_path}siconc.{y}.nc" for y in range(years[0], years[1])]

    # Open files with cftime decoder
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)

    result = []
    for file in files2load:    
    
        ds = xr.open_dataset(file, decode_times=time_coder).load()
        if log:
            print(f"File loaded: {file}", flush=True)

        # transform projection
        ds = reproject_to_latlon(ds)

        # Compute sea ice area
        mask_geo = (
            (ds.lat >= box[2]) & (ds.lat <= box[3]) &
            (ds.lon >= box[0]) & (ds.lon <= box[1])
        )

        mask_isice = ds.siconc > siconc_threshold
    
        gridcell_area = 25000 * 25000 # m^2
        ice_area = ((ds.siconc * mask_geo * mask_isice) * gridcell_area).sum(dim=('x','y'))
        ice_area = ice_area.where(ice_area>0,np.nan)

        result.append(ice_area)

    result = xr.concat(result, dim='time').rename('sea_ice_area')

    result.attrs['units'] = 'm^2'
    result.attrs['long_name'] = 'total sea ice area'
    result.attrs['bounding box'] = (
        f"Longitude: {box[0]}E to {box[1]}E, Latitude: {box[2]}N to {box[3]}N"
    )

    if grouping == 'annual.mean':
        result = result.groupby('time.year').mean('time')
    elif grouping == 'annual.max':
        result = result.groupby('time.year').max('time')
    elif grouping == 'annual.min':
        result = result.groupby('time.year').min('time')
    elif grouping == 'monthly.mean':
        pass
    
    print('Done!')
    return result
    
def hadlsst_ice_area(src_path,
    years=(1979, 2015),
    box=[-180, 180, -90, -50],
    siconc_threshold=0.15,
    grouping='annual.mean',
    log=True):

    """
    Compute total sea ice area from HadlSST_ice data within a geographic region.

    This function loads HadlSST sea ice concentration data, adds gridd cell areea, and calculates the 
    total sea ice area where concentration exceeds a given threshold.

    Parameters
    ----------
    src_path : str
        Path to directory containing SIC NetCDF files named as 'sic.<year>.nc'.
    years : tuple of int, optional
        Start and end year (exclusive) for processing, e.g., (1979, 2015).
    box : list of float, optional
        Geographic bounding box [lon_min, lon_max, lat_min, lat_max] for area calculation.
    siconc_threshold : float, optional
        Sea ice concentration threshold above which a grid cell is counted as ice-covered (default is 0.15).
    log : bool, optional
        If True, print progress messages during processing.

    Returns
    -------
    xarray.DataArray
        Time series of total sea ice area (in km²) per time step within the specified region.

    Notes
    -----
    - Compute grid cell size from lat/lon
    """

    files2load = [f"{src_path}siconc.{y}.nc" for y in range(years[0], years[1])]

    # Open files with cftime decoder
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)

    result = []
    for file in files2load:    
    
        ds = xr.open_dataset(file, decode_times=time_coder).load()
        if log:
            print(f"File loaded: {file}", flush=True)

        # add gridcell area
        ds = gridcell_area_hadley(ds, R=6371.0)

        # Crop to box
        ds = ds.sel(latitude=slice(box[2], box[3]), longitude=slice(box[0], box[1]))

        # commented line is ice extent, we want ice area
        #ice_area = ds.cell_area.where(ds.siconc > siconc_threshold, 0).sum(dim=('longitude','latitude'))

        mask_isice = ds.siconc > siconc_threshold
        ice_area = (ds.siconc * mask_isice * ds.cell_area).sum(dim=('longitude','latitude'))

        result.append(ice_area)

    result = xr.concat(result, dim='time').rename('sea_ice_area')

    result.attrs['units'] = 'm^2'
    result.attrs['long_name'] = 'total sea ice area'
    result.attrs['bounding box'] = (
        f"Longitude: {box[0]}E to {box[1]}E, Latitude: {box[2]}N to {box[3]}N"
    )

    if grouping == 'annual.mean':
        result = result.groupby('time.year').mean('time')
    elif grouping == 'annual.max':
        result = result.groupby('time.year').max('time')
    elif grouping == 'annual.min':
        result = result.groupby('time.year').min('time')
    elif grouping == 'monthly.mean':
        pass
    
    print('Done!')
    return result
