# so_ase package

## Submodules

## so_ase.eval_atmosphere module

## so_ase.eval_ocean module

### so_ase.eval_ocean.fesom_ocean_heat_transport_as_residual(src_path, mesh_diag_path, ref_date='2000-01-31', eval_date='2002-01-31', box=[-180, 180, -90, -60], rho=1028, cp=4190, log=True)

Compute the ocean heat transport (OHT) into a region of interest during a particular period as the residual of the total surface heat flux during the period and ocean heat
content change between the first and last timestep of the period, using FESOM2 model output.

* **Parameters:**
  * **src_path** (*str*) – Directory path to the FESOM2 NetCDF output files (e.g., fh.fesom.YYYY.nc, temp.fesom.YYYY.nc).
  * **mesh_diag_path** (*str*) – Directory path to the mesh diagnostic NetCDF file fesom.mesh.diag.nc.
  * **ref_date** (*str* *,* *optional*) – Start date of the analysis period in ‘YYYY-MM-DD’ format. Default is ‘2000-01-31’.
  * **eval_date** (*str* *,* *optional*) – End date of the analysis period in ‘YYYY-MM-DD’ format. Default is ‘2002-01-31’.
  * **box** (*list* *of* *float* *,* *optional*) – Geographic bounding box [lon_min, lon_max, lat_min, lat_max] to subset the spatial region of interest.
    Default is [-180, 180, -90, -60] (Southern Ocean).
  * **rho** (*float* *,* *optional*) – Seawater density in kg/m³. Default is 1028.
  * **cp** (*float* *,* *optional*) – Seawater specific heat capacity in J/(kg·K). Default is 4190.
  * **log** (*bool* *,* *optional*) – If True, prints logging information during processing. Default is True.
* **Returns:**
  The estimated ocean heat transport (OHT) in joules (J) over the specified period and region.
* **Return type:**
  float

### Notes

- The OHT is estimated as the residual of:

  where:
  : - ΔOHC is the change in ocean heat content between ref_date and eval_date
    - SHF is the surface heat flux (in W/m²)
    - The surface heat flux integral includes time (in seconds) and nodal area (m²)
- This function assumes monthly-averaged data and uses fixed calendar days per month (non-leap year).

### so_ase.eval_ocean.fesom_timeseries_of_mean_vertical_profile_in_region(src_path, mesh_diag_path, years=(1979, 2015), box=[-180, 180, -90, -60], varname='temp', log=True)

Compute a time series of area-weighted mean vertical profiles for a specified region
from FESOM2 model output.

### Parameters:

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

### Returns:

xarray.DataSet
: A time series of area-weighted mean vertical profiles in the specified region.
  The vertical levels are preserved along the ‘nz’ dimension.

## so_ase.eval_sea_ice module

### so_ase.eval_sea_ice.fesom_ice_area(src_path, mesh_diag_path, years=(1979, 2015), box=[-180, 180, -90, -60], siconc_threshold=0.15, grouping='annual.mean', log=True)

Compute total sea ice area within a specified geographic bounding box using FESOM2 output.

This function loads sea ice concentration (a_ice) from FESOM2 output files over a range
of years and calculates the sea ice area by summing the nodal areas where sea ice
concentration exceeds a given threshold, within a user-defined geographic region.

* **Parameters:**
  * **src_path** (*str*) – Path to the directory containing FESOM2 output NetCDF files named
    a_ice.fesom.{year}.nc.
  * **mesh_diag_path** (*str*) – Path to the directory containing the mesh diagnostic file fesom.mesh.diag.nc,
    which provides nodal coordinates and areas.
  * **years** (*tuple* *of* *int* *,* *optional*) – Start and end year (exclusive) for the time series. Defaults to (1979, 2015).
  * **box** (*list* *of* *float* *,* *optional*) – Geographic bounding box as [lon_min, lon_max, lat_min, lat_max] to define the region
    of interest. Defaults to [-180, 180, -90, -60] (entire Southern Hemisphere).
  * **aice_threshhold** (*float* *,* *optional*) – Minimum sea ice concentration (0–1) above which a node is considered ice-covered.
    Defaults to 0.15.
  * **grouping** (*str* *,* *optional*) – Type of aggregation for the time series. Options are ‘annual.mean’, ‘annual.max’,
    ‘annual.min’, or ‘monthly.mean’. Defaults to ‘annual.mean’.
  * **log** (*bool* *,* *optional*) – If True, prints progress and file loading information to standard output.
* **Returns:**
  Dataset containing the time series variable sea_ice_area (in square meters)
  representing total sea ice area within the defined region. The variables a_ice
  and the intermediate ice_mask are removed from the returned dataset.
* **Return type:**
  xarray.Dataset

### Notes

A binary sea ice mask is created where a_ice > aice_threshhold. The total sea ice area
is computed by summing the nodal areas (nod_area) over the masked nodes.

### so_ase.eval_sea_ice.fesom_ice_volume(src_path, mesh_diag_path, years=(1979, 2015), box=[-180, 180, -90, -60], grouping='annual.mean', log=True)

Compute total sea ice volume within a specified geographic bounding box using FESOM2 output.

This function loads sea ice concentration (a_ice) and thickness (m_ice) variables
from FESOM2 output files over a range of years, and calculates the sea ice volume by
summing the product of sea ice thickness, concentration, and nodal area over the nodes
that fall within a user-defined geographic region.

* **Parameters:**
  * **src_path** (*str*) – Path to the directory containing FESOM2 output NetCDF files named
    a_ice.fesom.{year}.nc and m_ice.fesom.{year}.nc.
  * **mesh_diag_path** (*str*) – Path to the directory containing the mesh diagnostic file fesom.mesh.diag.nc,
    which provides nodal coordinates and areas.
  * **years** (*tuple* *of* *int* *,* *optional*) – Start and end year (exclusive) for the time series. Defaults to (1979, 2015).
  * **box** (*list* *of* *float* *,* *optional*) – Geographic bounding box as [lon_min, lon_max, lat_min, lat_max] to define the region
    of interest. Defaults to [-180, 180, -90, -60] (entire Southern Hemisphere).
  * **log** (*bool* *,* *optional*) – If True, prints progress and file loading information to standard output.
* **Returns:**
  Dataset containing the time series variable sea_ice_volume (in cubic meters)
  representing total sea ice volume within the defined region. The original m_ice
  variable is dropped from the returned dataset.
* **Return type:**
  xarray.Dataset

### Notes

Sea ice volume is calculated as:

> sea_ice_volume = sum_over_nodes(area \* m_ice \* a_ice)

### so_ase.eval_sea_ice.hadlsst_ice_area(src_path, years=(1979, 2015), box=[-180, 180, -90, -50], siconc_threshold=0.15, grouping='annual.mean', log=True)

Compute total sea ice area from HadlSST_ice data within a geographic region.

This function loads HadlSST sea ice concentration data, adds gridd cell areea, and calculates the
total sea ice area where concentration exceeds a given threshold.

* **Parameters:**
  * **src_path** (*str*) – Path to directory containing SIC NetCDF files named as ‘sic.<year>.nc’.
  * **years** (*tuple* *of* *int* *,* *optional*) – Start and end year (exclusive) for processing, e.g., (1979, 2015).
  * **box** (*list* *of* *float* *,* *optional*) – Geographic bounding box [lon_min, lon_max, lat_min, lat_max] for area calculation.
  * **siconc_threshold** (*float* *,* *optional*) – Sea ice concentration threshold above which a grid cell is counted as ice-covered (default is 0.15).
  * **log** (*bool* *,* *optional*) – If True, print progress messages during processing.
* **Returns:**
  Time series of total sea ice area (in km²) per time step within the specified region.
* **Return type:**
  xarray.DataArray

### Notes

- Compute grid cell size from lat/lon

### so_ase.eval_sea_ice.nsidc_ice_area(src_path, years=(1979, 2015), box=[-180, 180, -90, -50], siconc_threshold=0.15, grouping='annual.mean', log=True)

Compute total sea ice area from NSIDC CDR data within a geographic region.

This function loads NSIDC sea ice concentration data, reprojects it from
polar stereographic coordinates to regular lat/lon, and calculates the
total sea ice area where concentration exceeds a given threshold.

* **Parameters:**
  * **src_path** (*str*) – Path to directory containing SIC NetCDF files named as ‘sic.<year>.nc’.
  * **years** (*tuple* *of* *int* *,* *optional*) – Start and end year (exclusive) for processing, e.g., (1979, 2015).
  * **box** (*list* *of* *float* *,* *optional*) – Geographic bounding box [lon_min, lon_max, lat_min, lat_max] for area calculation.
  * **siconc_threshold** (*float* *,* *optional*) – Sea ice concentration threshold above which a grid cell is counted as ice-covered (default is 0.15).
  * **log** (*bool* *,* *optional*) – If True, print progress messages during processing.
* **Returns:**
  Time series of total sea ice area (in km²) per time step within the specified region.
* **Return type:**
  xarray.DataArray

### Notes

- Assumes a constant 25 km x 25 km grid cell size.

## so_ase.helpers_mesh module

### so_ase.helpers_mesh.build_element_neighbors(elements)

Builds a list of neighboring elements for each triangle in the mesh.
Parameters:

> elements (array): An array of shape (ntri, 3) containing the node indices for each triangle.

Returns:
: list: A list of lists, where each sublist contains the indices of neighboring triangles for the corresponding triangle.

### so_ase.helpers_mesh.build_node_k_ring_neighbors(elements, node_idx, k)

Builds k-ring neighboring nodes for each node in the mesh.
Parameters:

> elements (array): An array of shape (ntri, 3) containing the node indices for each triangle.
> node_idx (array): An array of node indices.
> k (int): The ring number to compute (e.g., 1 for 1-ring neighbors).

Returns:
: list: A list of sets, where each set contains the k-ring neighboring node indices for the corresponding node.

### so_ase.helpers_mesh.build_node_neighbors(node_idx, elements)

### so_ase.helpers_mesh.find_nodes_in_box(mesh_diag_path, box=[-180, 180, -90, -60], log=True)

### so_ase.helpers_mesh.gridcell_area_hadley(ds, R=6371.0)

Compute and add grid cell area (km²) to an xarray.Dataset on a regular 1° lat-lon grid.

* **Parameters:**
  * **ds** (*xarray.Dataset*) – Dataset with ‘lat’ and ‘lon’ coordinates.
  * **R** (*float* *,* *optional*) – Radius of the Earth in kilometers (default: 6371 km).
* **Returns:**
  **ds_out** – New dataset with an additional variable ‘cell_area’.
* **Return type:**
  xarray.Dataset

### so_ase.helpers_mesh.read_aux3d(meshpath)

Reads vertical level information (e.g., depths) from an aux3d.out file,
skipping the first num_levels lines that typically contain header or
level-wise data, and returns one depth value per node.

Parameters:
: meshpath (str): Path to the directory containing aux3d.out and nod2d.out.

Returns:
: array of float: A list of depth values, one for each node.

### so_ase.helpers_mesh.read_cavity_depth_at_node(meshpath)

Reads ice base depth information from [cavity_depth@node.out](mailto:cavity_depth@node.out) file,

Parameters:
: meshpath (str): Path to the directory containing aux3d.out and nod2d.out.

Returns:
: list of float: A list of depth values, one for each node.

### so_ase.helpers_mesh.read_element_levels(meshpath, which='seafloor', raw=False)

Reads vertical level information (first active/last active layer index) from elvls.out/cavity_elvls.out file,

Parameters:
: meshpath (str): Path to the directory containing elvls/cavity_elvls.out.
  which (str): either <seafloor> or <cavity> for last active layer or first active layer.

Returns:
: array: A list of level indices values, one for each element.

### so_ase.helpers_mesh.read_elements(meshpath)

Reads element connectivity information from a elem2d.out file.

Parameters:
: meshpath (str): Path to the directory containing elem2d.out.

Returns:
: array: A list of elements, where each element is a tuple
  : of 0-based node indices (n1, n2, n3).

### so_ase.helpers_mesh.read_nodes(meshpath)

Reads 2D node coordinates from a nod2d.out file in a given mesh path.

Parameters:
: meshpath (str): Path to the directory containing nod2d.out.

Returns:
: node_lon (array): Array of node longitudes.
  node_lat (array): Array of node latitudes.
  node_idx (array): Array of node indices (0-based).
  node_coast (array): Array indicating if a node is coastal (1) or not (0).

### so_ase.helpers_mesh.reproject_to_latlon(ds, input_proj='+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84')

Transforms x/y coordinates of a polar-stereographic dataset (NSIDC datasets) into latitude/longitude and adds them as 2D coordinate variables.

Parameters:
: ds (xr.Dataset): Dataset in projected coordinates (e.g., EPSG:3412).
  proj_str (str): PROJ string for the dataset’s current projection.

Returns:
: xr.Dataset: Dataset with 2D ‘lat’ and ‘lon’ variables added.

## so_ase.helpers_plots module

### so_ase.helpers_plots.add_notebook_path_to_fig(fig, y_position=-0.05, fontsize=8)

Adds the absolute path of the Jupyter Notebook to an existing figure.

Parameters:
fig (matplotlib.figure.Figure): The existing figure to modify.
y_position (float): The vertical position of the text (default: 0.01, near the bottom).
fontsize (int): Font size of the path text (default: 8).

## so_ase.helpers_regression module

### so_ase.helpers_regression.anomaly2D_fesom(src_path, mesh_diag_path, ref_period=(2011, 2024), box=[-180, 180, -65, -55], depth=None, varname='a_ice', grouping='annual.mean', grid_data=True, grid_inc=0.25, log=True)

Compute 2D anomalies of FESOM data over a specified reference period, with optional horizontal interpolation to a regular lat-lon grid.

* **Parameters:**
  * **src_path** (*str*) – Path to the folder containing FESOM netCDF output files.
  * **mesh_diag_path** (*str*) – Path to the mesh diagnostic file (fesom.mesh.diag.nc).
  * **ref_period** (*tuple* *of* *int* *,* *optional*) – Start and end year (inclusive start, exclusive end) used to define the anomaly reference period. Default is (2011, 2024).
  * **box** (*list* *of* *float* *,* *optional*) – Geographic bounding box [lon_min, lon_max, lat_min, lat_max] used to subset the spatial domain. Default is global Southern Ocean sector.
  * **depth** (*float* *or* *None* *,* *optional*) – Vertical level (in meters) to extract if variable has a vertical extent (temp, salt, u, v, …). If None, it is assumed that surface data (sst, sic, hf, …) is used. Default is None.
  * **varname** (*str* *,* *optional*) – Variable name to process (e.g., ‘sst’, ‘a_ice’). Default is ‘a_ice’.
  * **grouping** (*str* *,* *optional*) – Temporal aggregation strategy, either ‘annual.mean’, ‘annual.max’, ‘annual.min’, or ‘monthly.mean’. Only ‘mean’ is supported for monthly grouping (no further grouping of monthly mean output). Default is ‘annual.mean’.
  * **grid_data** (*bool* *,* *optional*) – If True, interpolate unstructured data to a regular lat-lon grid. If False, return unstructured anomalies. Default is True.
  * **grid_inc** (*float* *,* *optional*) – Grid spacing in degrees for regular interpolation. Default is 0.25°.
  * **log** (*bool* *,* *optional*) – If True, print progress messages to standard output. Default is True.
* **Returns:**
  **ds_out** – Dataset containing the computed anomalies. If grid_data=True, the anomalies are on a regular (lat, lon) grid; otherwise they are on the unstructured mesh (nod2).
* **Return type:**
  xarray.Dataset

### so_ase.helpers_regression.dataset_regression_on_time_1D(ds, varname, split=None, log=True)

Performs linear regression of a variable in a Dataset over its time axis. The variable must only have a time dimension.
If split is provided, the time series is split and two regressions are performed.

* **Parameters:**
  * **ds** (*xarray.Dataset*) – The input dataset.
  * **varname** (*str*) – Name of the variable to regress.
  * **split** (*int* *,* *optional*) – Year to split the regression. If provided, runs two regressions: one for
    time <= split, one for time > split.
* **Returns:**
  One or two dictionaries, each containing:
  - ‘year’: the time values
  - ‘fit’: fitted linear trend (mx + c)
  - ‘slope’: slope (m)
  - ‘intercept’: intercept (c)
  - ‘pvalue’: p-value of slope
  - ‘rvalue’: correlation coefficient
* **Return type:**
  List[dict]

### so_ase.helpers_regression.grid_fesom(src_path, fname, mesh_diag_path, ref_period=(2011, 2024), box=[-180, 180, -65, -55], depth=None, varname='a_ice', grouping='annual.mean', grid_inc=0.25, log=True)

Read in 2D fesom data and do horizontal interpolation to a regular lat-lon grid.

* **Parameters:**
  * **src_path** (*str*) – Path to the folder containing FESOM netCDF output files.
  * **mesh_diag_path** (*str*) – Path to the mesh diagnostic file (fesom.mesh.diag.nc).
  * **ref_period** (*tuple* *of* *int* *,* *optional*) – Start and end year (inclusive start, exclusive end) used to define the anomaly reference period. Default is (2011, 2024).
  * **box** (*list* *of* *float* *,* *optional*) – Geographic bounding box [lon_min, lon_max, lat_min, lat_max] used to subset the spatial domain. Default is global Southern Ocean sector.
  * **depth** (*float* *or* *None* *,* *optional*) – Vertical level (in meters) to extract if variable has a vertical extent (temp, salt, u, v, …). If None, it is assumed that surface data (sst, sic, hf, …) is used. Default is None.
  * **varname** (*str* *,* *optional*) – Variable name to process (e.g., ‘sst’, ‘a_ice’). Default is ‘a_ice’.
  * **grouping** (*str* *,* *optional*) – Temporal aggregation strategy, either ‘annual.mean’, ‘annual.max’, ‘annual.min’, or ‘monthly.mean’. Only ‘mean’ is supported for monthly grouping (no further grouping of monthly mean output). Default is ‘annual.mean’.
  * **grid_inc** (*float* *,* *optional*) – Grid spacing in degrees for regular interpolation. Default is 0.25°.
  * **log** (*bool* *,* *optional*) – If True, print progress messages to standard output. Default is True.
* **Returns:**
  **ds_out** – Dataset containing the computed anomalies on a regular (lat, lon) grid
* **Return type:**
  xarray.Dataset

### so_ase.helpers_regression.regression2D_fesom(src_path, mesh_diag_path, years=(2011, 2024), box=[-180, 180, -65, -55], depth=None, varname='a_ice', grouping='annual.mean', grid_data=True, log=True)

Perform linear regression on a FESOM2 variable over time at each unstructured grid node,
and optionally interpolate the regression results to a regular lon-lat grid.

* **Parameters:**
  * **src_path** (*str*) – Path to the directory containing the FESOM2 NetCDF files.
  * **mesh_diag_path** (*str*) – Path to the fesom.mesh.diag.nc file containing mesh diagnostics (node coordinates).
  * **years** (*tuple* *of* *int* *,* *optional*) – Start and end year for the analysis (inclusive start, exclusive end). Default is (2011, 2024).
  * **box** (*list* *of* *float* *,* *optional*) – Geographic bounding box [lon_min, lon_max, lat_min, lat_max] used to subset the mesh. Default is global Southern Ocean sector.
  * **depth** (*float* *or* *None* *,* *optional*) – If specified, selects a vertical level (in meters) using nearest match for 3D variables.
  * **varname** (*str* *,* *optional*) – Variable name in the NetCDF files to be analyzed. Default is ‘sic’.
  * **grouping** (*str* *,* *optional*) – Temporal aggregation of the data before regression. Options are:
    - ‘annual.mean’
    - ‘annual.max’
    - ‘annual.min’
    - ‘monthly.mean’ (automatically deseasonalized)
  * **grid_data** (*bool* *,* *optional*) – If True, interpolates regression results to a regular lon-lat grid. Default is True.
  * **log** (*bool* *,* *optional*) – If True, prints progress messages. Default is True.
* **Returns:**
  **result** – Dataset containing regression parameters:
  - slope
  - intercept
  - rvalue
  - pvalue
  - stderr
  - intercept_stderr
  - lon and lat coordinates (if grid_data is False).

  If grid_data is True, the dataset is on regular lat-lon grid;
  otherwise, it’s indexed by the unstructured mesh node coordinate nod2.
* **Return type:**
  xarray.Dataset

### Notes

- This function uses scipy.stats.linregress to perform regression at each node.
- Data is deseasonalized when using monthly.mean.
- When grid_data=True, scipy.interpolate.griddata is used for nearest neighbor spatial interpolation.

### so_ase.helpers_regression.regression2D_hadlsst(src_path, years=(2011, 2024), box=[-180, 180, -65, -55], depth=None, varname='siconc', grouping='annual.mean', log=True)

Perform 2D linear regression on Hadley Centre SST or SIC data over a specified region and time range.

This function reads HadISST NetCDF data for sea surface temperature (SST) or sea ice concentration (SIC),
subsets it spatially and temporally, optionally aggregates or deseasonalizes it, and performs pixel-wise
linear regression to estimate temporal trends and statistical significance.

* **Parameters:**
  * **src_path** (*str*) – Path to the directory containing HadISST NetCDF files.
  * **years** (*tuple* *of* *int* *,* *optional*) – Start and end years (exclusive) for the analysis period. Default is (2011, 2024).
  * **box** (*list* *of* *float* *,* *optional*) – Geographic bounding box [lon_min, lon_max, lat_min, lat_max]. Default is [-180, 180, -65, -55].
  * **depth** (*int* *or* *None* *,* *optional*) – Reserved for compatibility; not used for HadISST data. Default is None.
  * **varname** (*str* *,* *optional*) – Variable name to analyze. Should be either ‘siconc’ for sea ice concentration or ‘sst’ for temperature.
    Default is ‘siconc’.
  * **grouping** (*str* *,* *optional*) – Temporal grouping method. Options:
    - ‘annual.mean’
    - ‘annual.max’
    - ‘annual.min’
    - ‘monthly.mean’ (deseasonalized)
  * **log** (*bool* *,* *optional*) – If True, print progress messages to stdout. Default is True.
* **Returns:**
  A dataset containing the following regression statistics for each grid cell:
  - slope : Trend over time (per year)
  - intercept : Y-intercept of regression line
  - rvalue : Pearson correlation coefficient
  - pvalue : Two-sided p-value for a hypothesis test whose null hypothesis is that the slope is zero
  - stderr : Standard error of the estimated slope
  - intercept_stderr : Standard error of the estimated intercept
* **Return type:**
  xarray.Dataset

### Notes

- Input NetCDF files must follow the naming convention {varname}.{year}.nc.
- This function assumes the HadISST grid format with dimensions: time, latitude, longitude, nv (ignored).
- Grid is assumed to be regular 1-degree.
- Monthly data is deseasonalized before regression using monthly climatology.

### so_ase.helpers_regression.regression2D_nsidc(src_path, years=(2011, 2024), box=[-180, 180, -65, -55], depth=None, varname='siconc', grouping='annual.mean', log=True)

Perform 2D linear regression on NSIDC SIC data over a specified region and time range.

This function reads HadISST NetCDF data for sea surface temperature (SST) or sea ice concentration (SIC),
subsets it spatially and temporally, optionally aggregates or deseasonalizes it, and performs pixel-wise
linear regression to estimate temporal trends and statistical significance.

* **Parameters:**
  * **src_path** (*str*) – Path to the directory containing HadISST NetCDF files.
  * **years** (*tuple* *of* *int* *,* *optional*) – Start and end years (exclusive) for the analysis period. Default is (2011, 2024).
  * **box** (*list* *of* *float* *,* *optional*) – Geographic bounding box [lon_min, lon_max, lat_min, lat_max]. Default is [-180, 180, -65, -55].
  * **depth** (*int* *or* *None* *,* *optional*) – Reserved for compatibility; not used for HadISST data. Default is None.
  * **varname** (*str* *,* *optional*) – Variable name to analyze. Should be either ‘siconc’ for sea ice concentration or ‘sst’ for temperature.
    Default is ‘siconc’.
  * **grouping** (*str* *,* *optional*) – Temporal grouping method. Options:
    - ‘annual.mean’
    - ‘annual.max’
    - ‘annual.min’
    - ‘monthly.mean’ (deseasonalized)
  * **log** (*bool* *,* *optional*) – If True, print progress messages to stdout. Default is True.
* **Returns:**
  A dataset containing the following regression statistics for each grid cell:
  - slope : Trend over time (per year)
  - intercept : Y-intercept of regression line
  - rvalue : Pearson correlation coefficient
  - pvalue : Two-sided p-value for a hypothesis test whose null hypothesis is that the slope is zero
  - stderr : Standard error of the estimated slope
  - intercept_stderr : Standard error of the estimated intercept
* **Return type:**
  xarray.Dataset

### Notes

- Input NetCDF files must follow the naming convention {varname}.{year}.nc.
- This function assumes the HadISST grid format with dimensions: time, latitude, longitude, nv (ignored).
- Grid is assumed to be regular 1-degree.
- Monthly data is deseasonalized before regression using monthly climatology.

### so_ase.helpers_regression.xrlinregress(first_samples, second_samples, dim)

## so_ase.plotting_maps module

### so_ase.plotting_maps.circular_shape(ax)

Clips the plotting area to a circular shape, typically used to create
polar plots or maps with circular boundaries.

Parameters:
: ax (matplotlib.axes._subplots.AxesSubplot):
  : The axes object to apply the circular boundary to.

Returns:
: matplotlib.axes._subplots.AxesSubplot:
  : The modified axes object with a circular clipping boundary.

### so_ase.plotting_maps.create_map(ax, extent='global', land=True, coastline=True, lon_inc=30, lat_inc=5, tick_labels=True, circular=True, zorder=1000)

Creates a map with optional land and coastline features,
gridlines.

Parameters:
: ax (matplotlib.axes._subplots.AxesSubplot):
  : The axes object to plot on.
  <br/>
  extent (list or str, optional):
  : The geographic extent of the map in the form [lon_min, lon_max, lat_min, lat_max] or str ‘global’, ‘southern_ocean’ or ‘arctic_ocean’.
    Defaults to ‘global’.
  <br/>
  land (bool, optional):
  : Whether to add land features to the map. Defaults to True.
  <br/>
  coastline (bool, optional):
  : Whether to add coastlines to the map. Defaults to True.
  <br/>
  lon_inc (int, optional):
  : Longitude grid interval. Defaults to 30.
  <br/>
  lat_inc (int, optional):
  : Latitude grid interval. Defaults to 5.
  <br/>
  tick_labels (bool, optional):
  : Whether to display longitude and latitude labels on the gridlines. Defaults to True.
  <br/>
  circular (bool, optional):
  : Whether to apply a circular shape to the map. Defaults to True.
  <br/>
  zorder (int, optional):
  : The z-order for the map features. Defaults to 1000.

Returns:
: matplotlib.axes._subplots.AxesSubplot:
  : The modified axes object with the Southern Ocean map.

## Module contents
