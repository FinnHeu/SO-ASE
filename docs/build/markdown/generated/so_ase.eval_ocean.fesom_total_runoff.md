# so_ase.eval_ocean.fesom_total_runoff

### fesom_total_runoff(src_path, meshpath, basin_mask_file, basins=[66], years=(1979, 2015), which='solid', ori=False, log=False, savepath='./', replace=False)

Calculate total runoff for specified basins from FESOM output files.

This function processes FESOM runoff data files, applies basin masks to isolate
runoff from specific drainage basins, and saves the total runoff time series
for each year as separate NetCDF files.

* **Parameters:**
  * **src_path** (*str*) – Path to directory containing FESOM runoff files.
  * **meshpath** (*str*) – Path to directory containing FESOM mesh files.
  * **basin_mask_file** (*str*) – Path to NetCDF file containing basin mask data.
  * **basins** (*list* *or* *int*) – List of basin IDs to include in calculation. Defaults to [66].
  * **years** (*tuple*) – Tuple of (start_year, end_year) for processing period. Defaults to (1979, 2015).
  * **which** (*str*) – Type of runoff to process (‘solid’ or ‘liquid’). Defaults to ‘solid’.
  * **ori** (*bool*) – If True, process original runoff files with ‘_ori’ suffix. Defaults to False.
  * **log** (*bool*) – If True, print progress messages. Defaults to False.
  * **savepath** (*str*) – Directory path for output files. Defaults to ‘./’.
  * **replace** (*bool*) – If True, overwrite existing output files. Defaults to False.
* **Returns:**
  None

### Notes

- Input files follow pattern: runoff_{which}.fesom.{year}.nc or runoff_{which}_ori.fesom.{year}.nc
- Output files follow pattern: {which}_runoff_basins_{basin_ids}.{year}.nc
- Runoff is calculated as sum over nodes in specified basins, weighted by nodal area
- Uses mesh diagnostics from fesom.mesh.diag.nc for nodal area information
