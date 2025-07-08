# so_ase.eval_sea_ice.py

import xarray as xr
import numpy as np


def fesom_ice_area(ds, mesh_diag, box=[-180, 180, -90, -60], aice_threshhold=0.15):
    """
    Calculate total sea ice area within a specified geographic bounding box using FESOM2 output.

    Parameters
    ----------
    ds : xarray.Dataset
        FESOM2 model output dataset containing the variable `a_ice` (sea ice concentration)
        on unstructured nodes over time with dimension `nod2`.
    mesh_diag : xarray.Dataset
        Corresponding mesh diagnostic dataset containing nodal coordinates (`lon`, `lat`)
        and nodal area (`nod_area`) with dimensions `nod2` and `nz`.
    box : list of float, optional
        Geographic bounding box `[lon_min, lon_max, lat_min, lat_max]` used to subset
        the region of interest. Defaults to `[-180, 180, -90, -60]` (entire Southern Hemisphere).
    aice_threshhold : float, optional
        Sea ice concentration threshold above which nodes are considered ice-covered.
        Default is `0.15`.

    Returns
    -------
    xarray.Dataset
        Dataset containing a single time series variable `sea_ice_area`, representing the
        total sea ice area (in square meters) over time within the specified region.
        The original `a_ice` and intermediate `ice_mask` variables are removed from the output.

    Notes
    -----
    The sea ice mask is created where `a_ice > aice_threshhold`.
    Total area is computed as the sum of nodal areas (from `mesh_diag.nod_area`) where ice is present.
    """

    # Find indices of nodes within the specified box
    inds = np.where(
        (mesh_diag.lon > box[0])
        & (mesh_diag.lon < box[1])
        & (mesh_diag.lat > box[2])
        & (mesh_diag.lat < box[3])
    )[0]

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

    return ds_cropped


def fesom_ice_volume(ds_aice, ds_mice, mesh_diag, box=[-180, 180, -90, -60]):
    """
    Calculate total sea ice volume within a specified geographic bounding box using FESOM2 output.

    Parameters
    ----------
    ds_aice : xarray.Dataset
        FESOM2 model output dataset containing sea ice concentration variable `a_ice`
        on unstructured nodes with dimension `nod2`.
    ds_mice : xarray.Dataset
        FESOM2 model output dataset containing sea ice thickness variable `m_ice`
        on unstructured nodes with dimension `nod2`.
    mesh_diag : xarray.Dataset
        Mesh diagnostic dataset containing nodal coordinates (`lon`, `lat`) and nodal area (`nod_area`)
        with dimensions `nod2` and `nz`.
    box : list of float, optional
        Geographic bounding box `[lon_min, lon_max, lat_min, lat_max]` used to subset
        the region of interest. Defaults to `[-180, 180, -90, -60]` (entire Southern Hemisphere).

    Returns
    -------
    xarray.Dataset
        Dataset containing a single time series variable `sea_ice_volume`, representing the
        total sea ice volume (in cubic meters) over time within the specified region.
        The original variable `m_ice` is dropped from the output dataset.

    Notes
    -----
    The total sea ice volume is calculated as the sum over all nodes of the product:
    nodal area × sea ice thickness (`m_ice`) × sea ice concentration (`a_ice`).
    """

    # Find indices of nmodes within the specified box
    inds = np.where(
        (mesh_diag.lon > box[0])
        & (mesh_diag.lon < box[1])
        & (mesh_diag.lat > box[2])
        & (mesh_diag.lat < box[3])
    )[0]

    # Crop datasets
    ds_aice_cropped = ds_aice.isel(nod2=inds)
    ds_mice_cropped = ds_mice.isel(nod2=inds)
    mesh_diag_cropped = mesh_diag.isel(nod2=inds)

    # Compute total ice volume (nodal area * nodal ice height * nodal ice concentration)
    sea_ice_volume = (
        mesh_diag_cropped.nod_area.isel(nz=0)
        * ds_mice_cropped.m_ice
        * ds_aice_cropped.a_ice
    ).sum(dim="nod2")

    # Sum over nodal areas
    ds_mice_cropped["sea_ice_volume"] = sea_ice_volume

    # Prepare output
    ds_mice_cropped = ds_mice_cropped.drop_vars(["m_ice"])
    ds_mice_cropped.sea_ice_volume.attrs["units"] = "$m^3$"
    ds_mice_cropped.sea_ice_volume.attrs["long_name"] = "total sea ice volume"
    ds_mice_cropped.sea_ice_volume.attrs["bounding box"] = (
        f"Longitude: {box[0]}E to {box[1]}E, Latitude: {box[2]}N to {box[3]}N"
    )

    return ds_mice_cropped
