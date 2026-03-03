# so_ase.eval_ocean.fesom_timeseries_of_mean_vertical_profile_in_region

### so_ase.eval_ocean.fesom_timeseries_of_mean_vertical_profile_in_region(src_path, mesh_diag_path, years=(1979, 2015), box=[-180, 180, -90, -60], varname='temp', log=True)

Compute a time series of area-weighted mean vertical profiles for a specified region
from FESOM2 model output.

## Parameters:

src_path
: Path to the directory containing annual FESOM2 output files (e.g., temp.fesom.YYYY.nc).

mesh_diag_path
: Path to the directory containing the mesh diagnostic file (fesom.mesh.diag.nc).

years
: Start and end year for the time series. Default is (1979, 2015).

box
: Geographic bounds of the region of interest in the format [lon_min, lon_max, lat_min, lat_max].
  Default is global Southern Ocean: [-180, 180, -90, -60].

varname
: Name of the variable to extract and average (must match variable name in NetCDF files).
  Default is ‘temp’.

log
: Whether to print progress information. Default is True.

## Returns:

xarray.DataSet
: A time series of area-weighted mean vertical profiles in the specified region.
  The vertical levels are preserved along the ‘nz’ dimension.
