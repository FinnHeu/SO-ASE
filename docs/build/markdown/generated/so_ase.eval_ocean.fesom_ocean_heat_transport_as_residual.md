# so_ase.eval_ocean.fesom_ocean_heat_transport_as_residual

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
