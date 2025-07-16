# helpers_regression.py

import xarray as xr
import numpy as np
from scipy.stats import linregress

def dataset_regression_on_time_1D(ds, varname, split=None):
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

    return output