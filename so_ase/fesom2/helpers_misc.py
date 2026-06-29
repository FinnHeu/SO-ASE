# so_ase/fesom2/helpers_misc.py

import xarray as xr
import numpy as np
from ..miscellaneous import seconds_per_month

def total_annual_from_monthly_mean(ds, var='fw'):
    """
    Convert monthly mean flux rates to total annual flux.

    Takes a dataset with monthly mean flux values (in units per second) and
    computes the total annual flux by multiplying each month's value by the
    number of seconds in that month, then summing over the time dimension.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset containing the flux variable with a 'time' dimension.
        Time coordinates should be monthly.
    var : str, optional
        Name of the flux variable to convert. Default is 'fw'.

    Returns
    -------
    xarray.DataArray
        Total annual flux summed over all months. NaN values in the input
        are treated as zero before summation.

    Notes
    -----
    Uses `so.seconds_per_month()` to compute the number of seconds in each
    month, accounting for varying month lengths and leap years.
    """
    ds['spm'] = (('time'), seconds_per_month(ds.groupby('time.year').mean().year.values))
    return (ds[var].fillna(0) * ds['spm']).sum(dim='time')