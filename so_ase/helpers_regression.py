# helpers_regression.py

import xarray as xr
import numpy as np
from scipy.stats import linregress
from scipy.interpolate import griddata
from .helpers_mesh import find_nodes_in_box, reproject_to_latlon

def regression2D_fesom(src_path, mesh_diag_path, years=(2011, 2024), box=[-180, 180, -65, -55], depth=None, varname='a_ice', grouping='annual.mean', grid_data=True, log=True):
    """
    Perform linear regression on a FESOM2 variable over time at each unstructured grid node,
    and optionally interpolate the regression results to a regular lon-lat grid.

    Parameters
    ----------
    src_path : str
        Path to the directory containing the FESOM2 NetCDF files.
    mesh_diag_path : str
        Path to the `fesom.mesh.diag.nc` file containing mesh diagnostics (node coordinates).
    years : tuple of int, optional
        Start and end year for the analysis (inclusive start, exclusive end). Default is (2011, 2024).
    box : list of float, optional
        Geographic bounding box [lon_min, lon_max, lat_min, lat_max] used to subset the mesh. Default is global Southern Ocean sector.
    depth : float or None, optional
        If specified, selects a vertical level (in meters) using nearest match for 3D variables.
    varname : str, optional
        Variable name in the NetCDF files to be analyzed. Default is 'sic'.
    grouping : str, optional
        Temporal aggregation of the data before regression. Options are:
        - 'annual.mean'
        - 'annual.max'
        - 'annual.min'
        - 'monthly.mean' (automatically deseasonalized)
    grid_data : bool, optional
        If True, interpolates regression results to a regular lon-lat grid. Default is True.
    log : bool, optional
        If True, prints progress messages. Default is True.

    Returns
    -------
    result : xarray.Dataset
        Dataset containing regression parameters:
        - slope
        - intercept
        - rvalue
        - pvalue
        - stderr
        - intercept_stderr
        - lon and lat coordinates (if `grid_data` is False).

        If `grid_data` is True, the dataset is on regular lat-lon grid;
        otherwise, it's indexed by the unstructured mesh node coordinate `nod2`.

    Notes
    -----
    - This function uses `scipy.stats.linregress` to perform regression at each node.
    - Data is deseasonalized when using `monthly.mean`.
    - When `grid_data=True`, `scipy.interpolate.griddata` is used for nearest neighbor spatial interpolation.

    """
    files2open = [f'{src_path}{varname}.fesom.{y}.nc' for y in range(years[0], years[-1])]
    
    inds = find_nodes_in_box(mesh_diag_path, box)

    if log:
        print('Loading files...')
        if depth is not None:
            print(f'Chosen depth level: {depth}m')
            ds = xr.open_mfdataset(files2open).isel(nod2=inds).sel(nz1=depth, method='nearest').squeeze().load()
        else:        
            ds = xr.open_mfdataset(files2open).isel(nod2=inds).load()
    if log:
        for f in files2open:
            print(f'Files loaded: {f}')
    
    ds = ds.transpose('time','nod2')
    
    freq, how = grouping.split('.')
    if freq == 'annual':
        if how == 'mean':
            ds = ds.groupby('time.year').mean()
        elif how == 'max':
            ds = ds.groupby('time.year').max()
        elif how == 'min':
            ds = ds.groupby('time.year').min()
        # Get time as float
        x = ds.year.values
    
    elif freq == 'monthly':
        if how != 'mean':
            raise Warning('For monthly data no further grouping is allowed.')
        # Deaseason
        ds = ds.groupby('time.month') - ds.groupby('time.month').mean()
        # Get time as float
        x = ds.time.dt.year.values + np.tile(np.arange(0,1,1/12), len(np.unique(ds.time.dt.year.values)))
    
    data = ds[varname].values
    
    # Perform regression
    slope = []
    intercept = []
    rvalue = []
    pvalue =[]
    stderr = []
    intercept_stderr = []
    
    for j in range(data.shape[-1]):
        reg = linregress(x, data[:,j])
        slope.append(reg.slope)
        intercept.append(reg.intercept)
        rvalue.append(reg.rvalue)
        pvalue.append(reg.pvalue)
        stderr.append(reg.stderr)
        intercept_stderr.append(reg.intercept_stderr)

    if grid_data:
        mesh_diag = xr.open_dataset(f'{src_path}fesom.mesh.diag.nc')
        points = np.column_stack((mesh_diag.lon.values[inds], mesh_diag.lat.values[inds]))
        
        # Define target regular grid
        lon_grid = np.arange(box[0], box[1], .5)
        lat_grid = np.arange(box[2], box[3], .5)
        lon_mesh, lat_mesh = np.meshgrid(lon_grid, lat_grid)
        
        #result = {'lon': lon_grid, 'lat': lat_grid}
        data_vars = {}
        for data, name in zip([slope, intercept, rvalue, pvalue, stderr, intercept_stderr], 
                ['slope', 'intercept', 'rvalue', 'pvalue', 'stderr', 'intercept_stderr']):
        
            # Interpolate to grid
            data_grid = griddata(
                points,
                data,
                (lon_mesh, lat_mesh),
                method='nearest')
            
            #result.update({name: data_grid})
            # Add to dataset
            data_vars[name] = (("lat", "lon"), data_grid)

            # Create the xarray dataset
            result = xr.Dataset(
            data_vars=data_vars,
            coords={
                "lon": lon_grid,
                "lat": lat_grid,
            }
        )
            
    else:
        # Create the dataset
        result = xr.Dataset(
            data_vars={
                'lon': ('nod2', mesh_diag.lon.isel(nod2=inds).values),
                'lat': ('nod2', mesh_diag.lat.isel(nod2=inds).values),
                
                'slope': ('nod2', slope),
                'intercept': ('nod2', intercept),
                'rvalue': ('nod2', rvalue),
                'pvalue': ('nod2', pvalue),
                'stderr': ('nod2', stderr),
                'intercept_stderr': ('nod2', intercept_stderr),
            },
            coords={
                'nod2': range(len(lon_grid))
            }
    )
    if log:
        print('Done!')
    return result

def regression2D_hadlsst(src_path, years=(2011, 2024), box=[-180, 180, -65, -55], depth=None, varname='siconc', grouping='annual.mean', log=True):
    """
    Perform 2D linear regression on Hadley Centre SST or SIC data over a specified region and time range.

    This function reads HadISST NetCDF data for sea surface temperature (SST) or sea ice concentration (SIC),
    subsets it spatially and temporally, optionally aggregates or deseasonalizes it, and performs pixel-wise
    linear regression to estimate temporal trends and statistical significance.

    Parameters
    ----------
    src_path : str
        Path to the directory containing HadISST NetCDF files.
    years : tuple of int, optional
        Start and end years (exclusive) for the analysis period. Default is (2011, 2024).
    box : list of float, optional
        Geographic bounding box [lon_min, lon_max, lat_min, lat_max]. Default is [-180, 180, -65, -55].
    depth : int or None, optional
        Reserved for compatibility; not used for HadISST data. Default is None.
    varname : str, optional
        Variable name to analyze. Should be either 'siconc' for sea ice concentration or 'sst' for temperature.
        Default is 'siconc'.
    grouping : str, optional
        Temporal grouping method. Options:
        - 'annual.mean'
        - 'annual.max'
        - 'annual.min'
        - 'monthly.mean' (deseasonalized)
    log : bool, optional
        If True, print progress messages to stdout. Default is True.

    Returns
    -------
    xarray.Dataset
        A dataset containing the following regression statistics for each grid cell:
        - slope : Trend over time (per year)
        - intercept : Y-intercept of regression line
        - rvalue : Pearson correlation coefficient
        - pvalue : Two-sided p-value for a hypothesis test whose null hypothesis is that the slope is zero
        - stderr : Standard error of the estimated slope
        - intercept_stderr : Standard error of the estimated intercept

    Notes
    -----
    - Input NetCDF files must follow the naming convention `{varname}.{year}.nc`.
    - This function assumes the HadISST grid format with dimensions: time, latitude, longitude, nv (ignored).
    - Grid is assumed to be regular 1-degree.
    - Monthly data is deseasonalized before regression using monthly climatology.
    """
    
    files2open = [f'{src_path}{varname}.{y}.nc' for y in range(years[0], years[-1])]
    

    if log:
        print('Loading files...')
        ds = xr.open_mfdataset(files2open).sel(longitude=slice(box[0], box[1]), latitude=slice(box[2], box[3])).load()
    if log:
        for f in files2open:
            print(f'Files loaded: {f}')
    
    ds = ds.transpose('time','longitude', 'latitude', 'nv')
    
    freq, how = grouping.split('.')
    if freq == 'annual':
        if how == 'mean':
            ds = ds.groupby('time.year').mean()
        elif how == 'max':
            ds = ds.groupby('time.year').max()
        elif how == 'min':
            ds = ds.groupby('time.year').min()
        # Get time as float
        x = ds.year.values
    
    elif freq == 'monthly':
        if how != 'mean':
            raise Warning('For monthly data no further grouping is allowed.')
        # Deaseason
        ds = ds.groupby('time.month') - ds.groupby('time.month').mean()
        # Get time as float
        x = ds.time.dt.year.values + np.tile(np.arange(0,1,1/12), len(np.unique(ds.time.dt.year.values)))
    
    data = ds[varname].values
    
    # Perform regression
    slope = np.ones((data.shape[1], data.shape[2]))
    intercept = np.ones_like(slope)
    rvalue = np.ones_like(slope)
    pvalue =np.ones_like(slope)
    stderr = np.ones_like(slope)
    intercept_stderr = np.ones_like(slope)
    
    for i in range(data.shape[1]):
        for j in range(data.shape[2]):
            reg = linregress(x, data[:,i,j])
            slope[i,j] = reg.slope
            intercept[i,j] = reg.intercept
            rvalue[i,j] = reg.rvalue
            pvalue[i,j] = reg.pvalue
            stderr[i,j] = reg.stderr
            intercept_stderr[i,j] = reg.intercept_stderr

    # Create the dataset
    result = xr.Dataset(
    data_vars={
        'slope': (('longitude','latitude'), slope),
        'intercept': (('longitude','latitude'), intercept),
        'rvalue': (('longitude','latitude'), rvalue),
        'pvalue': (('longitude','latitude'), pvalue),
        'stderr': (('longitude','latitude'), stderr),
        'intercept_stderr': (('longitude','latitude'), intercept_stderr),
    },
    coords={
        'longitude': ds.longitude,
        'latitude': ds.latitude
        }
)
    if log:
        print('Done!')
    return result

def regression2D_nsidc(src_path, years=(2011, 2024), box=[-180, 180, -65, -55], depth=None, varname='siconc', grouping='annual.mean', log=True):
    """
    Perform 2D linear regression on NSIDC SIC data over a specified region and time range.

    This function reads HadISST NetCDF data for sea surface temperature (SST) or sea ice concentration (SIC),
    subsets it spatially and temporally, optionally aggregates or deseasonalizes it, and performs pixel-wise
    linear regression to estimate temporal trends and statistical significance.

    Parameters
    ----------
    src_path : str
        Path to the directory containing HadISST NetCDF files.
    years : tuple of int, optional
        Start and end years (exclusive) for the analysis period. Default is (2011, 2024).
    box : list of float, optional
        Geographic bounding box [lon_min, lon_max, lat_min, lat_max]. Default is [-180, 180, -65, -55].
    depth : int or None, optional
        Reserved for compatibility; not used for HadISST data. Default is None.
    varname : str, optional
        Variable name to analyze. Should be either 'siconc' for sea ice concentration or 'sst' for temperature.
        Default is 'siconc'.
    grouping : str, optional
        Temporal grouping method. Options:
        - 'annual.mean'
        - 'annual.max'
        - 'annual.min'
        - 'monthly.mean' (deseasonalized)
    log : bool, optional
        If True, print progress messages to stdout. Default is True.

    Returns
    -------
    xarray.Dataset
        A dataset containing the following regression statistics for each grid cell:
        - slope : Trend over time (per year)
        - intercept : Y-intercept of regression line
        - rvalue : Pearson correlation coefficient
        - pvalue : Two-sided p-value for a hypothesis test whose null hypothesis is that the slope is zero
        - stderr : Standard error of the estimated slope
        - intercept_stderr : Standard error of the estimated intercept

    Notes
    -----
    - Input NetCDF files must follow the naming convention `{varname}.{year}.nc`.
    - This function assumes the HadISST grid format with dimensions: time, latitude, longitude, nv (ignored).
    - Grid is assumed to be regular 1-degree.
    - Monthly data is deseasonalized before regression using monthly climatology.
    """
    
    files2open = [f'{src_path}{varname}.{y}.nc' for y in range(years[0], years[-1])]
    

    if log:
        print('Loading files...')
        ds = xr.open_mfdataset(files2open)#.sel(longitude=slice(box[0], box[1]), latitude=slice(box[2], box[3])).load()
    if log:
        for f in files2open:
            print(f'Files loaded: {f}')

    ds = reproject_to_latlon(ds)
    #ds = ds.transpose('time','longitude', 'latitude', 'nv')
    
    freq, how = grouping.split('.')
    if freq == 'annual':
        if how == 'mean':
            ds = ds.groupby('time.year').mean()
        elif how == 'max':
            ds = ds.groupby('time.year').max()
        elif how == 'min':
            ds = ds.groupby('time.year').min()
        # Get time as float
        x = ds.year.values
    
    elif freq == 'monthly':
        if how != 'mean':
            raise Warning('For monthly data no further grouping is allowed.')
        # Deaseason
        ds = ds.groupby('time.month') - ds.groupby('time.month').mean()
        # Get time as float
        x = ds.time.dt.year.values + np.tile(np.arange(0,1,1/12), len(np.unique(ds.time.dt.year.values)))
    
    data = ds[varname].values
    
    # Perform regression
    slope = np.ones((data.shape[1], data.shape[2]))
    intercept = np.ones_like(slope)
    rvalue = np.ones_like(slope)
    pvalue =np.ones_like(slope)
    stderr = np.ones_like(slope)
    intercept_stderr = np.ones_like(slope)
    
    for i in range(data.shape[1]):
        for j in range(data.shape[2]):
            reg = linregress(x, data[:,i,j])
            slope[i,j] = reg.slope
            intercept[i,j] = reg.intercept
            rvalue[i,j] = reg.rvalue
            pvalue[i,j] = reg.pvalue
            stderr[i,j] = reg.stderr
            intercept_stderr[i,j] = reg.intercept_stderr

    # Create the dataset
    result = xr.Dataset(
    data_vars={
        'slope': (('longitude','latitude'), slope),
        'intercept': (('longitude','latitude'), intercept),
        'rvalue': (('longitude','latitude'), rvalue),
        'pvalue': (('longitude','latitude'), pvalue),
        'stderr': (('longitude','latitude'), stderr),
        'intercept_stderr': (('longitude','latitude'), intercept_stderr),
    },
    coords={
        'longitude': ds.x,
        'latitude': ds.y
        }
)
    if log:
        print('Done!')
    return result

def anomaly2D_fesom(src_path, mesh_diag_path, ref_period=(2011, 2024), box=[-180, 180, -65, -55], depth=None, varname='a_ice', grouping='annual.mean', grid_data=True, grid_inc=.25, log=True):

    """
    Compute 2D anomalies of FESOM data over a specified reference period, with optional horizontal interpolation to a regular lat-lon grid.

    Parameters
    ----------
    src_path : str
        Path to the folder containing FESOM netCDF output files.
    mesh_diag_path : str
        Path to the mesh diagnostic file (`fesom.mesh.diag.nc`).
    ref_period : tuple of int, optional
        Start and end year (inclusive start, exclusive end) used to define the anomaly reference period. Default is (2011, 2024).
    box : list of float, optional
        Geographic bounding box [lon_min, lon_max, lat_min, lat_max] used to subset the spatial domain. Default is global Southern Ocean sector.
    depth : float or None, optional
        Vertical level (in meters) to extract if variable has a vertical extent (temp, salt, u, v, ...). If None, it is assumed that surface data (sst, sic, hf, ...) is used. Default is None.
    varname : str, optional
        Variable name to process (e.g., 'sst', 'a_ice'). Default is 'a_ice'.
    grouping : str, optional
        Temporal aggregation strategy, either 'annual.mean', 'annual.max', 'annual.min', or 'monthly.mean'. Only 'mean' is supported for monthly grouping (no further grouping of monthly mean output). Default is 'annual.mean'.
    grid_data : bool, optional
        If True, interpolate unstructured data to a regular lat-lon grid. If False, return unstructured anomalies. Default is True.
    grid_inc : float, optional
        Grid spacing in degrees for regular interpolation. Default is 0.25°.
    log : bool, optional
        If True, print progress messages to standard output. Default is True.

    Returns
    -------
    ds_out : xarray.Dataset
        Dataset containing the computed anomalies. If `grid_data=True`, the anomalies are on a regular (lat, lon) grid; otherwise they are on the unstructured mesh (nod2).
    """
    
    files2open = [f'{src_path}{varname}.fesom.{y}.nc' for y in range(ref_period[0], ref_period[-1])]
    
    inds = find_nodes_in_box(mesh_diag_path, box)

    if log:
        print('Loading files...')
        if depth is not None:
            print(f'Chosen depth level: {depth}m')
            ds = xr.open_mfdataset(files2open).isel(nod2=inds).sel(nz1=depth, method='nearest').squeeze().load()
        else:        
            ds = xr.open_mfdataset(files2open).isel(nod2=inds).load()
    if log:
        for f in files2open:
            print(f'Files loaded: {f}')
    
    ds = ds.transpose('time','nod2')
    
    freq, how = grouping.split('.')
    if freq == 'annual':
        if how == 'mean':
            ds = ds.groupby('time.year').mean()
        elif how == 'max':
            ds = ds.groupby('time.year').max()
        elif how == 'min':
            ds = ds.groupby('time.year').min()

        # Remove the long-term mean from the annual means
        ds = ds - ds.mean(dim='year') 
        
    elif freq == 'monthly':
        if how != 'mean':
            raise Warning('For monthly data no further grouping is allowed.')
        
        # Remove the climatology from the monthly mean data
        ds = ds.groupby('time.month') - ds.groupby('time.month').mean()

    if grid_data:
        if log:
            print('Gridding data...')
            
        mesh_diag = xr.open_dataset(f'{src_path}fesom.mesh.diag.nc')
        points = np.column_stack((mesh_diag.lon.values[inds], mesh_diag.lat.values[inds]))
        
        # Define target regular grid
        lon_grid = np.arange(box[0], box[1] + grid_inc, grid_inc)
        lat_grid = np.arange(box[2], box[3] + grid_inc, grid_inc)
        lon_mesh, lat_mesh = np.meshgrid(lon_grid, lat_grid)
        
        data_vars = {}
        data_gridded = np.zeros((ds[varname].shape[0], len(lat_grid), len(lon_grid))) * np.nan
        
        #for data, name in zip(ds[varname].value, varname):
        data = ds[varname].values
        for i in range(ds[varname].shape[0]):
            # Interpolate to grid
            data_gridded[i,:,:] = griddata(
                points,
                data[i,:],
                (lon_mesh, lat_mesh),
                method='nearest')
            
            # Add to dataset
            time_dim_name = list(ds[varname].dims)[0]
            data_vars[varname] = ((time_dim_name, "lat", "lon"), data_gridded)

            # Create the xarray dataset
            ds_out = xr.Dataset(
            data_vars=data_vars,
            coords={
                time_dim_name: ds[time_dim_name],
                "lon": lon_grid,
                "lat": lat_grid,
            }
        )

    else:
        ds_out = ds

    if log:
        print('Done!')
        
    return ds_out

def dataset_regression_on_time_1D(ds, varname, split=None, log=True):
    """
    Performs linear regression of a variable in a Dataset over its time axis. The variable must only have a time dimension.
    If `split` is provided, the time series is split and two regressions are performed.

    Parameters
    ----------
    ds : xarray.Dataset
        The input dataset.

    varname : str
        Name of the variable to regress.

    split : int, optional
        Year to split the regression. If provided, runs two regressions: one for
        time <= split, one for time > split.

    Returns
    -------
    List[dict]
        One or two dictionaries, each containing:
        - 'year': the time values
        - 'fit': fitted linear trend (mx + c)
        - 'slope': slope (m)
        - 'intercept': intercept (c)
        - 'pvalue': p-value of slope
        - 'rvalue': correlation coefficient
    """
    da = ds[varname]

    # Identify time dimension
    if "time" in da.dims:
        time_dim = "time"
    elif "year" in da.dims:
        time_dim = "year"
    else:
        raise ValueError("Time dimension must be 'time' or 'year'.")

    time = ds[time_dim].values
    other_dims = [dim for dim in da.dims if dim != time_dim]

    if other_dims:
        raise ValueError("Only 1D time series (no extra dimensions) supported in this version.")

    y = da.values
    output = []

    if split is None:
        reg = linregress(time, y)
        fit = reg.slope * time + reg.intercept
        output.append({
            "year": time,
            "fit": fit,
            "slope": reg.slope,
            "intercept": reg.intercept,
            "pvalue": reg.pvalue,
            "rvalue": reg.rvalue
        })
    else:
        mask1 = time <= split
        mask2 = time > split

        if mask1.any():
            reg1 = linregress(time[mask1], y[mask1])
            fit1 = reg1.slope * time[mask1] + reg1.intercept
            output.append({
                "year": time[mask1],
                "fit": fit1,
                "slope": reg1.slope,
                "intercept": reg1.intercept,
                "pvalue": reg1.pvalue,
                "rvalue": reg1.rvalue
            })

        if mask2.any():
            reg2 = linregress(time[mask2], y[mask2])
            fit2 = reg2.slope * time[mask2] + reg2.intercept
            output.append({
                "year": time[mask2],
                "fit": fit2,
                "slope": reg2.slope,
                "intercept": reg2.intercept,
                "pvalue": reg2.pvalue,
                "rvalue": reg2.rvalue
            })

    if log:
        print('Done!')
    return output