# SO-ASE

## Ocean

<!-- !! processed by numpydoc !! -->

### fesom_ocean_heat_transport_as_residual(src_path, mesh_diag_path, ref_date='2000-01-31', eval_date='2002-01-31', box=[-180, 180, -90, -60], rho=1028, cp=4190, log=True)

Compute the ocean heat transport (OHT) into a region of interest during a particular period as the residual of the total surface heat flux during the period and ocean heat
content change between the first and last timestep of the period, using FESOM2 model output.

* **Parameters:**
  **src_path**
  : Directory path to the FESOM2 NetCDF output files (e.g., fh.fesom.YYYY.nc, temp.fesom.YYYY.nc).

  **mesh_diag_path**
  : Directory path to the mesh diagnostic NetCDF file fesom.mesh.diag.nc.

  **ref_date**
  : Start date of the analysis period in ‘YYYY-MM-DD’ format. Default is ‘2000-01-31’.

  **eval_date**
  : End date of the analysis period in ‘YYYY-MM-DD’ format. Default is ‘2002-01-31’.

  **box**
  : Geographic bounding box [lon_min, lon_max, lat_min, lat_max] to subset the spatial region of interest.
    Default is [-180, 180, -90, -60] (Southern Ocean).

  **rho**
  : Seawater density in kg/m³. Default is 1028.

  **cp**
  : Seawater specific heat capacity in J/(kg·K). Default is 4190.

  **log**
  : If True, prints logging information during processing. Default is True.
* **Returns:**
  float
  : The estimated ocean heat transport (OHT) in joules (J) over the specified period and region.

### Notes

- The OHT is estimated as the residual of:

  where:
  : - ΔOHC is the change in ocean heat content between ref_date and eval_date
    - SHF is the surface heat flux (in W/m²)
    - The surface heat flux integral includes time (in seconds) and nodal area (m²)
- This function assumes monthly-averaged data and uses fixed calendar days per month (non-leap year).

<!-- !! processed by numpydoc !! -->

### fesom_timeseries_of_mean_vertical_profile_in_region(src_path, mesh_diag_path, years=(1979, 2015), box=[-180, 180, -90, -60], varname='temp', log=True)

Compute a time series of area-weighted mean vertical profiles for a specified region 
from FESOM2 model output.

* **Parameters:**
  **src_path**
  : Path to the directory containing annual FESOM2 output files (e.g., temp.fesom.YYYY.nc).

  **mesh_diag_path**
  : Path to the directory containing the mesh diagnostic file (fesom.mesh.diag.nc).

  **years**
  : Start and end year for the time series. Default is (1979, 2015).

  **box**
  : Geographic bounds of the region of interest in the format [lon_min, lon_max, lat_min, lat_max].
    Default is global Southern Ocean: [-180, 180, -90, -60].

  **varname**
  : Name of the variable to extract and average (must match variable name in NetCDF files).
    Default is ‘temp’.

  **log**
  : Whether to print progress information. Default is True.
* **Returns:**
  xarray.DataSet
  : A time series of area-weighted mean vertical profiles in the specified region.
    The vertical levels are preserved along the ‘nz’ dimension.

<!-- !! processed by numpydoc !! -->

### fesom_total_kinetic_energy(src_path, mesh_diag_path, meshpath, years=(1979, 2015), mask='cavity', log=False, savepath='./')

<!-- !! processed by numpydoc !! -->

### fesom_total_runoff(src_path, meshpath, basin_mask_file, basins=[66], years=(1979, 2015), which='solid', ori=False, log=False, savepath='./', replace=False)

Calculate total runoff for specified basins from FESOM output files.

This function processes FESOM runoff data files, applies basin masks to isolate
runoff from specific drainage basins, and saves the total runoff time series
for each year as separate NetCDF files.

Parameters:
: src_path (str): Path to directory containing FESOM runoff files.
  meshpath (str): Path to directory containing FESOM mesh files.
  basin_mask_file (str): Path to NetCDF file containing basin mask data.
  basins (list or int): List of basin IDs to include in calculation. Defaults to [66].
  years (tuple): Tuple of (start_year, end_year) for processing period. Defaults to (1979, 2015).
  which (str): Type of runoff to process (‘solid’ or ‘liquid’). Defaults to ‘solid’.
  ori (bool): If True, process original runoff files with ‘_ori’ suffix. Defaults to False.
  log (bool): If True, print progress messages. Defaults to False.
  savepath (str): Directory path for output files. Defaults to ‘./’.
  replace (bool): If True, overwrite existing output files. Defaults to False.

Returns:
: None

Notes:
: - Input files follow pattern: runoff_{which}.fesom.{year}.nc or runoff_{which}_ori.fesom.{year}.nc
  - Output files follow pattern: {which}_runoff_basins_{basin_ids}.{year}.nc
  - Runoff is calculated as sum over nodes in specified basins, weighted by nodal area
  - Uses mesh diagnostics from fesom.mesh.diag.nc for nodal area information

<!-- !! processed by numpydoc !! -->

## Sea Ice

<!-- !! processed by numpydoc !! -->

### fesom_ice_volume(src_path, mesh_diag_path, years=(1979, 2015), box=[-180, 180, -90, -60], savepath='./', log=True)

Compute and save total sea ice volume time series from FESOM2 output within a specified
geographic bounding box.

This function loads sea ice concentration (a_ice) and thickness (m_ice) variables
from FESOM2 output files over a range of years, and calculates the sea ice volume by
summing the product of sea ice thickness, concentration, and nodal area over the nodes
that fall within a user-defined geographic region. One output file is written per
year containing the full time series for that year. If an output file already exists, that year is skipped.

* **Parameters:**
  **src_path**
  : Path to the directory containing FESOM2 output NetCDF files named
    a_ice.fesom.{year}.nc and m_ice.fesom.{year}.nc.

  **mesh_diag_path**
  : Path to the directory containing the mesh diagnostic file fesom.mesh.diag.nc,
    which provides nodal coordinates and areas.

  **years**
  : Start and end year (exclusive) for the time series. Defaults to (1979, 2015).

  **box**
  : Geographic bounding box as [lon_min, lon_max, lat_min, lat_max] to define the region
    of interest. Defaults to [-180, 180, -90, -60] (entire Southern Hemisphere).

  **savepath**
  : Directory where the output NetCDF files will be written. The directory
    is created if it does not already exist.

  **log**
  : If True, prints progress and file loading information to standard output.
* **Returns:**
  None
  : The function does not return any objects. Results are written directly
    to disk as NetCDF files.

### Notes

Sea ice volume is calculated as:

> sea_ice_volume = sum_over_nodes(area \* m_ice \* a_ice)
<!-- !! processed by numpydoc !! -->

### fesom_sea_ice_area(src_path, mesh_diag_path, years=(1979, 2025), box=[-180, 180, -90, -60], siconc_threshold=0.15, savepath='./', log=True)

Compute and save total sea ice area time series from FESOM2 output within a specified
geographic bounding box.

This function processes yearly FESOM2 sea ice concentration files
(a_ice.fesom.<year>.nc), computes the total sea ice area by summing the
nodal areas where sea ice concentration exceeds a given threshold, and
writes the result to disk as NetCDF files. One output file is written per
year containing the full time series for that year. If an output file already exists, that year is skipped.

The sea ice area is calculated only for mesh nodes that fall within the
user-defined geographic bounding box. Spatial masking is based on the
FESOM mesh diagnostic file (fesom.mesh.diag.nc), which is loaded once
and reused for all years.

Output filenames are self-describing and include:
: - the processed year,
  - the geographic bounding box formatted using N/S/E/W notation.

* **Parameters:**
  **src_path**
  : Path to the directory containing FESOM2 sea ice concentration files
    named a_ice.fesom.<year>.nc.

  **mesh_diag_path**
  : Path to the directory containing the FESOM mesh diagnostic file
    fesom.mesh.diag.nc.

  **years**
  : Start and end year (end year exclusive) defining the range of years
    to process. Default is (1979, 2025).

  **box**
  : Geographic bounding box specified as
    [lon_min, lon_max, lat_min, lat_max] in degrees.
    Longitudes are expected in degrees east, latitudes in degrees north.
    Default is [-180, 180, -90, -60].

  **siconc_threshold**
  : Minimum sea ice concentration (range 0–1) required for a node to be
    considered ice-covered. Default is 0.15.

  **savepath**
  : Directory where the output NetCDF files will be written. The directory
    is created if it does not already exist.

  **log**
  : If True, print progress messages, including file loading, skipping,
    and saving information. Default is True.
* **Returns:**
  None
  : The function does not return any objects. Results are written directly
    to disk as NetCDF files.

<!-- !! processed by numpydoc !! -->

### hadlsst_ice_area(src_path, years=(1979, 2015), box=[-180, 180, -90, -50], siconc_threshold=0.15, grouping='annual.mean', log=True)

Compute total sea ice area from HadlSST_ice data within a geographic region.

This function loads HadlSST sea ice concentration data, adds gridd cell areea, and calculates the 
total sea ice area where concentration exceeds a given threshold.

* **Parameters:**
  **src_path**
  : Path to directory containing SIC NetCDF files named as ‘sic.<year>.nc’.

  **years**
  : Start and end year (exclusive) for processing, e.g., (1979, 2015).

  **box**
  : Geographic bounding box [lon_min, lon_max, lat_min, lat_max] for area calculation.

  **siconc_threshold**
  : Sea ice concentration threshold above which a grid cell is counted as ice-covered (default is 0.15).

  **log**
  : If True, print progress messages during processing.
* **Returns:**
  xarray.DataArray
  : Time series of total sea ice area (in km²) per time step within the specified region.

### Notes

- Compute grid cell size from lat/lon

<!-- !! processed by numpydoc !! -->

### nsidc_ice_area(src_path, years=(1979, 2015), box=[-180, 180, -90, -50], siconc_threshold=0.15, grouping='annual.mean', log=True)

Compute total sea ice area from NSIDC CDR data within a geographic region.

This function loads NSIDC sea ice concentration data, reprojects it from 
polar stereographic coordinates to regular lat/lon, and calculates the 
total sea ice area where concentration exceeds a given threshold.

* **Parameters:**
  **src_path**
  : Path to directory containing SIC NetCDF files named as ‘sic.<year>.nc’.

  **years**
  : Start and end year (exclusive) for processing, e.g., (1979, 2015).

  **box**
  : Geographic bounding box [lon_min, lon_max, lat_min, lat_max] for area calculation.

  **siconc_threshold**
  : Sea ice concentration threshold above which a grid cell is counted as ice-covered (default is 0.15).

  **log**
  : If True, print progress messages during processing.
* **Returns:**
  xarray.DataArray
  : Time series of total sea ice area (in km²) per time step within the specified region.

### Notes

- Assumes a constant 25 km x 25 km grid cell size.

<!-- !! processed by numpydoc !! -->

## Icebergs

<!-- !! processed by numpydoc !! -->

### read_iceberg_initial_files(icebergpath)

Reads iceberg initial condition files from the specified directory.

Parameters:
: icebergpath (str): Path to the directory containing iceberg initial condition files.

Returns:
: array: Arrays containing iceberg longitude, latitude, height, face element indices,
  : length, and scaling factors.

<!-- !! processed by numpydoc !! -->

## Ice cavities

<!-- !! processed by numpydoc !! -->

### fesom_subshelf_freshwaterflux(src_path, mesh_diag_path, mesh_path, mask, years=(1979, 2015), log=False, savepath='./')

Compute and save integrated subshelf basal melt time series from FESOM 
freshwater flux output [m3/s].

This function loads monthly freshwater flux fields (fw) from a sequence 
of FESOM output files, multiplies them by the mean horizontal node area 
from a mesh diagnostics file, applies a cavity mask, and integrates the 
melt signal over all masked nodes to obtain a single subshelf melt 
time series for each year. The result for each year is saved as a 
standalone NetCDF file.

* **Parameters:**
  **src_path**
  : Directory containing the annual FESOM freshwater flux files 
    named `fw.fesom.<year>.nc`.

  **mesh_diag_path**
  : Directory containing `fesom.mesh.diag.nc` with nod_area and 
    related mesh diagnostics.

  **mesh_path**
  : Path to the FESOM mesh directory (used for building the cavity mask).

  **mask**
  : Identifier for which ocean cavity mask to apply.  
    If `dict['name'] == 'all'` is given, the mask is generated using 
    `build_cavity_mask(mesh_path, which='node')`.  
    If `dict['name'] == 'Amery'` is given a .kml file called `Amery.kml` 
    will be loaded from `dict['kml_path']` from which the boolean mask will be generated.
    If an array-like mask is supplied, it is used directly for `isel`.

  **years**
  : Start and end years. The function processes files for all years in 
    `range(years[0], years[-1])`. (Note: the end year is excluded.)

  **log**
  : If True, print a message whenever a file is saved or skipped.

  **savepath**
  : Directory in which to save the output files. The files are created as
    `subshelf_melt_<mask>.<year>.nc`.

<!-- !! processed by numpydoc !! -->

### fesom_subshelf_heatflux(src_path, mesh_diag_path, mesh_path, mask, years=(1979, 2015), log=False, savepath='./')

Compute and save integrated subshelf heat flux time series from FESOM 
heat flux output.

This function loads monthly subshelf heat flux fields (fh) from a 
sequence of annual FESOM output files, multiplies them by the mean 
horizontal node area from a mesh diagnostics file, applies a cavity mask, 
and integrates the heat flux over all masked nodes to obtain a single 
area-integrated time series for each year. The result for each year is 
saved as a standalone NetCDF file.

* **Parameters:**
  **src_path**
  : Directory containing the annual FESOM heat flux files 
    named `fh.fesom.<year>.nc`.

  **mesh_diag_path**
  : Directory containing `fesom.mesh.diag.nc`, which includes 
    nod_area and other mesh diagnostics.

  **mesh_path**
  : Path to the FESOM mesh directory. Used to construct the subshelf 
    cavity mask.

  **mask**
  : Identifier for which ocean cavity mask to apply.  
    If `dict['name'] == 'all'` is given, the mask is generated using 
    `build_cavity_mask(mesh_path, which='node')`.  
    If `dict['name'] == 'Amery'` is given a .kml file called `Amery.kml` 
    will be loaded from `dict['kml_path']` from which the boolean mask will be generated.
    If an array-like mask is supplied, it is used directly for `isel`.

  **years**
  : Start and end years. The function processes files for all years in 
    `range(years[0], years[-1])`. (Note: the end year is excluded.)

  **log**
  : If True, print status messages when saving or skipping files.

  **savepath**
  : Directory into which output files will be written.  
    Files are saved as `subshelf_heatflux_<mask>.<year>.nc`.
* **Returns:**
  None
  : The function writes NetCDF files but does not return a value.

<!-- !! processed by numpydoc !! -->

### fesom_subshelf_hydrography(src_path, mesh_diag_path, mesh_path, mask, years=(1979, 1980), variable='temp', log=False, savepath='./')

Extract subshelf hydrography fields from FESOM output for a given node mask
and save the subsetted fields to new NetCDF files.

* **Parameters:**
  **src_path**
  : Directory containing yearly FESOM variable files, e.g.
    `/data/fesom/output/`.  
    Each file is expected to follow the format:
    `{variable}.fesom.{year}.nc`.

  **mesh_diag_path**
  : Directory containing the FESOM mesh diagnostic file
    `fesom.mesh.diag.nc`. This file must include the field
    `nod_area` used to add node area information to the output.

  **mesh_path**
  : Directory containing the FESOM mesh files. Required by the
    cavity mask builders.

  **mask**
  : Node mask used to subset the FESOM mesh:
    - If a **dict**, it must contain at least a ‘name’ key:
    > * If `mask['name'] == 'all'`: use the full cavity mask.
    > * Otherwise: use a regional mask defined by a KML file
    >   provided in `mask['kml_path']`.
    - If a **list** or **np.array**, it is interpreted as a boolean or
      integer mask directly indexable along the `nod2` dimension.

  **years**
  : Two-element tuple defining the year range to process.
    The function processes all years in `range(years[0], years[1])`.
    For example, `years=(1979, 1980)` processes only 1979.

  **variable**
  : One or more variable names to extract from the source files.

  **log**
  : If True, print progress messages.

  **savepath**
  : Directory where the subsetted NetCDF files will be written.

<!-- !! processed by numpydoc !! -->

### freshwaterflux_to_massflux_Gtm(src_path, dst_path, rho_fw=1000, log=True)

Convert monthly mean subshelf melt time series (m³/s) into monthly integrated
mass fluxes (Gt/month).

This function reads the output files generated by
`fesom_subshelf_freshwaterflux()`, which contain a variable
`subshelf_melt` representing monthly mean freshwater melt rates in
m³/s. For each file, the function converts the freshwater flux into an
ice-mass-equivalent flux (kg/s), multiplies each monthly mean by the exact number
of seconds in that month (including leap years), 
to obtain a total melt rate in gigatonnes per month (Gt/m). A new NetCDF
file with the suffix `_GTM` is written for each processed input file.

* **Parameters:**
  **src_path**
  : Directory containing the input files produced by
    `fesom_subshelf_melt()`.

  **dst_path**
  : Directory for writing the output files

  **rho_fw**
  : Density of freshwater in kg/m³.  
    Used to convert freshwater volume flux (m³/s) into mass flux (kg/s).

  **log**
  : If True, print progress messages when opening, skipping, or saving files.

<!-- !! processed by numpydoc !! -->

### freshwaterflux_to_massflux_Gty(src_path, dst_path, rho_fw=1000, log=True)

Convert monthly mean subshelf melt time series (m³/s) into annual integrated
mass fluxes (Gt/yr).

This function reads the output files generated by
`fesom_subshelf_freshwaterflux()`, which contain a variable
`subshelf_melt` representing monthly mean freshwater melt rates in
m³/s. For each file, the function converts the freshwater flux into an
ice-mass-equivalent flux (kg/s), multiplies each monthly mean by the exact number
of seconds in that month (including leap years), and sums over each year
to obtain a total melt rate in gigatonnes per year (Gt/yr). A new NetCDF
file with the suffix `_GTY` is written for each processed input file.

* **Parameters:**
  **src_path**
  : Directory containing the input files produced by
    `fesom_subshelf_melt()`.

  **dst_path**
  : Directory for writing the output files

  **rho_fw**
  : Density of freshwater in kg/m³.  
    Used to convert freshwater volume flux (m³/s) into mass flux (kg/s).

  **log**
  : If True, print progress messages when opening, skipping, or saving files.

<!-- !! processed by numpydoc !! -->

## Mesh

<!-- !! processed by numpydoc !! -->

### add_element_volumes(mesh_diag)

Adds layer thickness and element volume to the mesh_diag xarray.Dataset.
Parameters:

> mesh_diag (xarray.Dataset): Dataset containing ‘nz’ and ‘nod_area’.

Returns:
: xarray.Dataset: Updated dataset with ‘layer_thickness’ and ‘elem_volume’.

<!-- !! processed by numpydoc !! -->

### add_nodal_volumes(mesh_diag)

Adds layer thickness and nodal volume to the mesh_diag xarray.Dataset.
Parameters:

> mesh_diag (xarray.Dataset): Dataset containing ‘nz’ and ‘nod_area’.

Returns:
: xarray.Dataset: Updated dataset with ‘layer_thickness’ and ‘nod_volume’.

<!-- !! processed by numpydoc !! -->

### build_cavity_mask(meshpath, which='element')

Builds a cavity mask for either elements or nodes based on cavity levels.
A node is considered in the cavity if all its connected elements are cavity elements.
Inverting the mask gives the open ocean mask.

Parameters:
: meshpath (str): Path to the directory containing mesh files.
  which (str): ‘element’ to build mask for elements, ‘node’ for nodes.

Returns:
: array of bool: Cavity mask array.

<!-- !! processed by numpydoc !! -->

### build_cavity_regional_mask(meshpath, kml_path, which='Filchner-Ronne')

Builds a regional cavity mask for nodes based on a KML-defined polygon.

Parameters:
: meshpath (str): Path to the directory containing mesh files.
  kml_path (str): Path to the directory containing KML files.
  which (str): Name of the region (used to select the KML file).

Returns:
: array of bool: Regional cavity mask for nodes.

<!-- !! processed by numpydoc !! -->

### build_element_neighbors(elements)

Builds a list of neighboring elements for each triangle in the mesh.
Parameters:

> elements (array): An array of shape (ntri, 3) containing the node indices for each triangle.

Returns:
: list: A list of lists, where each sublist contains the indices of neighboring triangles for the corresponding triangle.

<!-- !! processed by numpydoc !! -->

### build_elements_of_nodes(elements, node_idx)

For each node in the mesh, return the list of element indices
(triangles) that contain that node.

Parameters:
: elements (array): (ntri, 3) array of triangle vertex indices.
  node_idx (array): Array of node indices (e.g. np.arange(N)).

Returns:
: list: A list of lists. The i-th list contains the triangle indices
  : of all elements that include node i.

<!-- !! processed by numpydoc !! -->

### build_node_k_ring_neighbors(elements, node_idx, k)

Builds k-ring neighboring nodes for each node in the mesh.
Parameters:

> elements (array): An array of shape (ntri, 3) containing the node indices for each triangle.
> node_idx (array): An array of node indices.
> k (int): The ring number to compute (e.g., 1 for 1-ring neighbors).

Returns:
: list: A list of sets, where each set contains the k-ring neighboring node indices for the corresponding node.

<!-- !! processed by numpydoc !! -->

### build_node_neighbors(elements, node_idx)

Builds neighboring nodes for each node in the mesh.

Parameters:
: elements (array): An array of shape (ntri, 3) containing the node indices for each triangle.
  node_idx (array): An array of node indices.

Returns:
: list: A list of lists, where each sublist contains the neighboring node indices for the corresponding node.

<!-- !! processed by numpydoc !! -->

### build_runoff_basin_mask(meshpath, runoff_file, which='liquid')

Assigns FESOM mesh nodes to runoff basins defined in a runoff_maps.nc file.

Parameters:
: meshpath (str): Path to the directory containing mesh files (nod2d.out).
  runoff_file (str): Path to the runoff_maps.nc file containing basin definitions.
  which (str): Either ‘liquid’ (for arrival_point_id) or ‘solid’ (for calving_point_id).

Returns:
: xarray.Dataset: Dataset containing basin IDs with node coordinates as dimensions.
  : Variables include ‘basin_id’ with basin ID for each mesh node.
    Nodes not in any basin will have ID 0.

<!-- !! processed by numpydoc !! -->

### find_nodes_in_box(mesh_diag_path, box=[-180, 180, -90, -60], log=True)

Finds node indices within a specified geographical box.
Parameters:

> mesh_diag_path (str): Path to the directory containing fesom.mesh.diag.nc.
> box (list): List of [lon_min, lon_max, lat_min, lat_max].
> log (bool): If True, prints the number of found nodes.

Returns:
: array: Array of node indices within the specified box.

<!-- !! processed by numpydoc !! -->

### gridcell_area_hadley(ds, R=6371.0)

Compute and add grid cell area (km²) to an xarray.Dataset on a regular 1° lat-lon grid.

* **Parameters:**
  **ds**
  : Dataset with ‘lat’ and ‘lon’ coordinates.

  **R**
  : Radius of the Earth in kilometers (default: 6371 km).
* **Returns:**
  **ds_out**
  : New dataset with an additional variable ‘cell_area’.

<!-- !! processed by numpydoc !! -->

### level_idx_to_depth(meshpath, which='seafloor', raw=False)

Converts level indices to depth values for either seafloor or cavity levels.

Parameters:
: meshpath (str): Path to the directory containing mesh files.
  which (str): either <seafloor> or <cavity> for last active layer or first active layer. 
  raw (bool): If True, reads raw level indices. If False, reads processed level indices.

Returns:
: array of float: A list of depth values, one for each element.

<!-- !! processed by numpydoc !! -->

### mesh2vtk(meshpath, which='seafloor')

Converts a FESOM2 mesh defined by nod2d.out and elem2d.out files into VTK format.

Parameters:
: meshpath (str): Path to the directory containing nod2d.out and elem2d.out.
  which (str): ‘seafloor’ to include last active layer, ‘cavity’ for first active layer, ‘none’ for none of both.

Returns:
: None: Writes a VTK file named ‘mesh_output.vtk’ in the current directory.

<!-- !! processed by numpydoc !! -->

### read_aux3d(meshpath)

Reads vertical level information (e.g., depths) from an aux3d.out file,
skipping the first num_levels lines that typically contain header or
level-wise data, and returns one depth value per node.

Parameters:
: meshpath (str): Path to the directory containing aux3d.out and nod2d.out.

Returns:
: array of float: A list of depth values, one for each node.

<!-- !! processed by numpydoc !! -->

### read_cavity_depth_at_node(meshpath)

Reads ice base depth information from [cavity_depth@node.out](mailto:cavity_depth@node.out) file,

Parameters:
: meshpath (str): Path to the directory containing aux3d.out and nod2d.out.

Returns:
: list of float: A list of depth values, one for each node.

<!-- !! processed by numpydoc !! -->

### read_depth_zlev(meshpath)

Reads vertical layer depth information from depth_zlev.out file,

Parameters:
: meshpath (str): Path to the directory containing depth_zlev.out.

Returns:
: array of float: A list of depth values for each vertical layer.
  int: Number of vertical layers.

<!-- !! processed by numpydoc !! -->

### read_element_levels(meshpath, which='seafloor', raw=False, python_indexing=False)

Reads vertical level information (first active/last active layer index) from elvls.out/cavity_elvls.out file,

Parameters:
: meshpath (str): Path to the directory containing elvls/cavity_elvls.out.
  which (str): either <seafloor> or <cavity> for last active layer or first active layer. 
  raw (bool): If True, reads raw level indices. If False, reads processed level indices.
  python_indexing (bool): If True, converts indices to 0-based indexing.

Returns:
: array: A list of level indices values, one for each element.

<!-- !! processed by numpydoc !! -->

### read_elements(meshpath)

Reads element connectivity information from a elem2d.out file.

Parameters:
: meshpath (str): Path to the directory containing elem2d.out.

Returns:
: array: A list of elements, where each element is a tuple
  : of 0-based node indices (n1, n2, n3).

<!-- !! processed by numpydoc !! -->

### read_nodes(meshpath)

Reads 2D node coordinates from a nod2d.out file in a given mesh path.

Parameters:
: meshpath (str): Path to the directory containing nod2d.out.

Returns:
: node_lon (array): Array of node longitudes.
  node_lat (array): Array of node latitudes.
  node_idx (array): Array of node indices (0-based).
  node_coast (array): Array indicating if a node is coastal (1) or not (0).

<!-- !! processed by numpydoc !! -->

### reproject_to_latlon(ds, input_proj='+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84')

Transforms x/y coordinates of a polar-stereographic dataset (NSIDC datasets) into latitude/longitude and adds them as 2D coordinate variables.

Parameters:
: ds (xr.Dataset): Dataset in projected coordinates (e.g., EPSG:3412).
  proj_str (str): PROJ string for the dataset’s current projection.

Returns:
: xr.Dataset: Dataset with 2D ‘lat’ and ‘lon’ variables added.

<!-- !! processed by numpydoc !! -->

## Plots

<!-- !! processed by numpydoc !! -->

### add_notebook_path_to_fig(fig, y_position=-0.05, fontsize=8)

Adds the absolute path of the Jupyter Notebook to an existing figure.

Parameters:
fig (matplotlib.figure.Figure): The existing figure to modify.
y_position (float): The vertical position of the text (default: 0.01, near the bottom).
fontsize (int): Font size of the path text (default: 8).

<!-- !! processed by numpydoc !! -->

### remove_axes_frame(ax, left=True, right=True, top=True, bottom=True, ticks=True, labels=True)

Remove or hide selected sides of a matplotlib Axes frame.

* **Parameters:**
  **ax**
  : The axes to modify.

  **left, right, top, bottom**
  : If True, remove the corresponding spine.

  **ticks**
  : If True, remove ticks on removed spines.

  **labels**
  : If True, remove tick labels on removed spines.

<!-- !! processed by numpydoc !! -->
<!-- !! processed by numpydoc !! -->

### circular_shape(ax)

Clips the plotting area to a circular shape, typically used to create
polar plots or maps with circular boundaries.

Parameters:
: ax (matplotlib.axes._subplots.AxesSubplot):
  : The axes object to apply the circular boundary to.

Returns:
: matplotlib.axes._subplots.AxesSubplot:
  : The modified axes object with a circular clipping boundary.

<!-- !! processed by numpydoc !! -->

### create_map(ax, extent='global', land=True, coastline=True, lon_inc=30, lat_inc=5, tick_labels=True, circular=True, zorder=1000)

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

<!-- !! processed by numpydoc !! -->

### plot_on_elements(ax, lon, lat, elements, data, mask, vmin='None', vmax='None', cmap='RdBu', zorder=1)

Plots data on a triangular mesh defined by elements.

Parameters:
: ax (matplotlib.axes._subplots.AxesSubplot):
  : The axes object to plot on.
  <br/>
  lon (np.ndarray):
  : Array of longitudes of the mesh nodes.
  <br/>
  lat (np.ndarray):
  : Array of latitudes of the mesh nodes.
  <br/>
  elements (np.ndarray):
  : Array defining the triangular elements by indices of the nodes.
  <br/>
  data (np.ndarray):
  : Array of data values at each node.
  <br/>
  mask (np.ndarray):
  : Boolean array indicating which elements to plot.
  <br/>
  vmin (float or str, optional):
  : Minimum data value for colormap scaling. If ‘None’, uses min of data. Defaults to ‘None’.
  <br/>
  vmax (float or str, optional):
  : Maximum data value for colormap scaling. If ‘None’, uses max of data. Defaults to ‘None’.
  <br/>
  cmap (str, optional):
  : Colormap to use for plotting. Defaults to “RdBu”.
  <br/>
  zorder (int, optional):
  : The z-order for the plot. Defaults to 1.

Returns: 
: matplotlib.collections.Tripcolor:
  : The tripcolor object created by the plot.

<!-- !! processed by numpydoc !! -->

## Regression

<!-- !! processed by numpydoc !! -->

### anomaly2D_fesom(src_path, mesh_diag_path, ref_period=(2011, 2024), box=[-180, 180, -65, -55], depth=None, varname='a_ice', grouping='annual.mean', grid_data=True, grid_inc=0.25, log=True)

Compute 2D anomalies of FESOM data over a specified reference period, with optional horizontal interpolation to a regular lat-lon grid.

* **Parameters:**
  **src_path**
  : Path to the folder containing FESOM netCDF output files.

  **mesh_diag_path**
  : Path to the mesh diagnostic file (fesom.mesh.diag.nc).

  **ref_period**
  : Start and end year (inclusive start, exclusive end) used to define the anomaly reference period. Default is (2011, 2024).

  **box**
  : Geographic bounding box [lon_min, lon_max, lat_min, lat_max] used to subset the spatial domain. Default is global Southern Ocean sector.

  **depth**
  : Vertical level (in meters) to extract if variable has a vertical extent (temp, salt, u, v, …). If None, it is assumed that surface data (sst, sic, hf, …) is used. Default is None.

  **varname**
  : Variable name to process (e.g., ‘sst’, ‘a_ice’). Default is ‘a_ice’.

  **grouping**
  : Temporal aggregation strategy, either ‘annual.mean’, ‘annual.max’, ‘annual.min’, or ‘monthly.mean’. Only ‘mean’ is supported for monthly grouping (no further grouping of monthly mean output). Default is ‘annual.mean’.

  **grid_data**
  : If True, interpolate unstructured data to a regular lat-lon grid. If False, return unstructured anomalies. Default is True.

  **grid_inc**
  : Grid spacing in degrees for regular interpolation. Default is 0.25°.

  **log**
  : If True, print progress messages to standard output. Default is True.
* **Returns:**
  **ds_out**
  : Dataset containing the computed anomalies. If grid_data=True, the anomalies are on a regular (lat, lon) grid; otherwise they are on the unstructured mesh (nod2).

<!-- !! processed by numpydoc !! -->

### dataset_regression_on_time_1D(ds, varname, split=None, log=True)

Performs linear regression of a variable in a Dataset over its time axis. The variable must only have a time dimension.
If split is provided, the time series is split and two regressions are performed.

* **Parameters:**
  **ds**
  : The input dataset.

  **varname**
  : Name of the variable to regress.

  **split**
  : Year to split the regression. If provided, runs two regressions: one for
    time <= split, one for time > split.
* **Returns:**
  List[dict]
  : One or two dictionaries, each containing:
    - ‘year’: the time values
    - ‘fit’: fitted linear trend (mx + c)
    - ‘slope’: slope (m)
    - ‘intercept’: intercept (c)
    - ‘pvalue’: p-value of slope
    - ‘rvalue’: correlation coefficient

<!-- !! processed by numpydoc !! -->

### regression2D_fesom(src_path, mesh_diag_path, years=(2011, 2024), box=[-180, 180, -65, -55], depth=None, varname='a_ice', grouping='annual.mean', grid_data=True, log=True)

Perform linear regression on a FESOM2 variable over time at each unstructured grid node,
and optionally interpolate the regression results to a regular lon-lat grid.

* **Parameters:**
  **src_path**
  : Path to the directory containing the FESOM2 NetCDF files.

  **mesh_diag_path**
  : Path to the fesom.mesh.diag.nc file containing mesh diagnostics (node coordinates).

  **years**
  : Start and end year for the analysis (inclusive start, exclusive end). Default is (2011, 2024).

  **box**
  : Geographic bounding box [lon_min, lon_max, lat_min, lat_max] used to subset the mesh. Default is global Southern Ocean sector.

  **depth**
  : If specified, selects a vertical level (in meters) using nearest match for 3D variables.

  **varname**
  : Variable name in the NetCDF files to be analyzed. Default is ‘sic’.

  **grouping**
  : Temporal aggregation of the data before regression. Options are:
    - ‘annual.mean’
    - ‘annual.max’
    - ‘annual.min’
    - ‘monthly.mean’ (automatically deseasonalized)

  **grid_data**
  : If True, interpolates regression results to a regular lon-lat grid. Default is True.

  **log**
  : If True, prints progress messages. Default is True.
* **Returns:**
  **result**
  : Dataset containing regression parameters:
    - slope
    - intercept
    - rvalue
    - pvalue
    - stderr
    - intercept_stderr
    - lon and lat coordinates (if grid_data is False).
    <br/>
    If grid_data is True, the dataset is on regular lat-lon grid;
    otherwise, it’s indexed by the unstructured mesh node coordinate nod2.

### Notes

- This function uses scipy.stats.linregress to perform regression at each node.
- Data is deseasonalized when using monthly.mean.
- When grid_data=True, scipy.interpolate.griddata is used for nearest neighbor spatial interpolation.

<!-- !! processed by numpydoc !! -->

### regression2D_hadlsst(src_path, years=(2011, 2024), box=[-180, 180, -65, -55], depth=None, varname='siconc', grouping='annual.mean', log=True)

Perform 2D linear regression on Hadley Centre SST or SIC data over a specified region and time range.

This function reads HadISST NetCDF data for sea surface temperature (SST) or sea ice concentration (SIC),
subsets it spatially and temporally, optionally aggregates or deseasonalizes it, and performs pixel-wise
linear regression to estimate temporal trends and statistical significance.

* **Parameters:**
  **src_path**
  : Path to the directory containing HadISST NetCDF files.

  **years**
  : Start and end years (exclusive) for the analysis period. Default is (2011, 2024).

  **box**
  : Geographic bounding box [lon_min, lon_max, lat_min, lat_max]. Default is [-180, 180, -65, -55].

  **depth**
  : Reserved for compatibility; not used for HadISST data. Default is None.

  **varname**
  : Variable name to analyze. Should be either ‘siconc’ for sea ice concentration or ‘sst’ for temperature.
    Default is ‘siconc’.

  **grouping**
  : Temporal grouping method. Options:
    - ‘annual.mean’
    - ‘annual.max’
    - ‘annual.min’
    - ‘monthly.mean’ (deseasonalized)

  **log**
  : If True, print progress messages to stdout. Default is True.
* **Returns:**
  xarray.Dataset
  : A dataset containing the following regression statistics for each grid cell:
    - slope : Trend over time (per year)
    - intercept : Y-intercept of regression line
    - rvalue : Pearson correlation coefficient
    - pvalue : Two-sided p-value for a hypothesis test whose null hypothesis is that the slope is zero
    - stderr : Standard error of the estimated slope
    - intercept_stderr : Standard error of the estimated intercept

### Notes

- Input NetCDF files must follow the naming convention {varname}.{year}.nc.
- This function assumes the HadISST grid format with dimensions: time, latitude, longitude, nv (ignored).
- Grid is assumed to be regular 1-degree.
- Monthly data is deseasonalized before regression using monthly climatology.

<!-- !! processed by numpydoc !! -->

### regression2D_nsidc(src_path, years=(2011, 2024), box=[-180, 180, -65, -55], depth=None, varname='siconc', grouping='annual.mean', log=True)

Perform 2D linear regression on NSIDC SIC data over a specified region and time range.

This function reads HadISST NetCDF data for sea surface temperature (SST) or sea ice concentration (SIC),
subsets it spatially and temporally, optionally aggregates or deseasonalizes it, and performs pixel-wise
linear regression to estimate temporal trends and statistical significance.

* **Parameters:**
  **src_path**
  : Path to the directory containing HadISST NetCDF files.

  **years**
  : Start and end years (exclusive) for the analysis period. Default is (2011, 2024).

  **box**
  : Geographic bounding box [lon_min, lon_max, lat_min, lat_max]. Default is [-180, 180, -65, -55].

  **depth**
  : Reserved for compatibility; not used for HadISST data. Default is None.

  **varname**
  : Variable name to analyze. Should be either ‘siconc’ for sea ice concentration or ‘sst’ for temperature.
    Default is ‘siconc’.

  **grouping**
  : Temporal grouping method. Options:
    - ‘annual.mean’
    - ‘annual.max’
    - ‘annual.min’
    - ‘monthly.mean’ (deseasonalized)

  **log**
  : If True, print progress messages to stdout. Default is True.
* **Returns:**
  xarray.Dataset
  : A dataset containing the following regression statistics for each grid cell:
    - slope : Trend over time (per year)
    - intercept : Y-intercept of regression line
    - rvalue : Pearson correlation coefficient
    - pvalue : Two-sided p-value for a hypothesis test whose null hypothesis is that the slope is zero
    - stderr : Standard error of the estimated slope
    - intercept_stderr : Standard error of the estimated intercept

### Notes

- Input NetCDF files must follow the naming convention {varname}.{year}.nc.
- This function assumes the HadISST grid format with dimensions: time, latitude, longitude, nv (ignored).
- Grid is assumed to be regular 1-degree.
- Monthly data is deseasonalized before regression using monthly climatology.

<!-- !! processed by numpydoc !! -->

### xrlinregress(first_samples, second_samples, dim)

<!-- !! processed by numpydoc !! -->

## Miscellaneous

<!-- !! processed by numpydoc !! -->

### lon_to_360(lon)

Converts longitude in -180/180E convention to 0/360E convention.

* **Parameters:**
  **lon**
  : longitude array in -180/180E convention
* **Returns:**
  **lon_360**

<!-- !! processed by numpydoc !! -->

### read_kml_coords(path, close_ring=True)

Read coordinates from a KML file and return a list of (lon, lat) tuples.
If the coordinate sequence is not closed, the first coordinate is appended at the end.

* **Parameters:**
  **path**
  : Path to the .kml file.
* **Returns:**
  **coords**
  : List of (lon, lat) tuples; first coordinate appended at end to close the ring.

<!-- !! processed by numpydoc !! -->

### seconds_per_month(years)

Return the number of seconds in each month for given year(s).

* **Parameters:**
  **years**
  : Year or list of years (e.g., 1980 or [1980, 1981, 1982])
* **Returns:**
  dict
  : {year: [seconds_in_Jan, seconds_in_Feb, …, seconds_in_Dec]}

<!-- !! processed by numpydoc !! -->
