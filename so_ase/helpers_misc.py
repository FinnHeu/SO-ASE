import calendar
import numpy as np
from datetime import datetime, timedelta

def seconds_per_month(years):
    """
    Return the number of seconds in each month for given year(s).

    Parameters
    ----------
    years : int or iterable of int
        Year or list of years (e.g., 1980 or [1980, 1981, 1982])

    Returns
    -------
    dict
        {year: [seconds_in_Jan, seconds_in_Feb, ..., seconds_in_Dec]}
    """
    
    # Accept single year or iterable
    if isinstance(years, int):
        years = [years]
    
    secs = []
    
    for y in years:
        for m in range(1, 13):
            # start of month
            start = datetime(y, m, 1)
            # start of next month (handle December→January rollover)
            if m == 12:
                end = datetime(y + 1, 1, 1)
            else:
                end = datetime(y, m + 1, 1)
            secs.append((end - start).total_seconds())
    
    return np.array(secs)