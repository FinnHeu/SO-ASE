import calendar
import numpy as np
import xml.etree.ElementTree as ET
from pathlib import Path
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

def read_kml_coords(path, close_ring=True):
    """
    Read coordinates from a KML file and return a list of (lon, lat) tuples.
    If the coordinate sequence is not closed, the first coordinate is appended at the end.

    Parameters
    ----------
    path : str or Path
        Path to the .kml file.

    Returns
    -------
    coords : list of (float, float)
        List of (lon, lat) tuples; first coordinate appended at end to close the ring.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"KML file not found: {path}")

    tree = ET.parse(str(path))
    root = tree.getroot()

    # Find all <coordinates> elements regardless of namespace by checking tag suffix
    coord_texts = []
    for elem in root.iter():
        if elem.tag is None:
            continue
        if elem.tag.split('}')[-1].lower() == 'coordinates':
            if elem.text and elem.text.strip():
                coord_texts.append(elem.text.strip())

    if not coord_texts:
        # No coordinates found
        return []

    # Parse all coordinate blocks and concatenate them
    coords = []
    for txt in coord_texts:
        # coordinates are whitespace-separated tuples like "lon,lat[,alt]"
        parts = [p for p in txt.replace('\\n', ' ').split() if p.strip()]
        for p in parts:
            pieces = p.split(',')
            if len(pieces) < 2:
                continue
            try:
                lon = float(pieces[0])
                lat = float(pieces[1])
            except ValueError:
                # Fallback: try replacing comma decimal separators (very rare in KML)
                try:
                    lon = float(pieces[0].replace(',', '.'))
                    lat = float(pieces[1].replace(',', '.'))
                except Exception:
                    continue
            coords.append((lon, lat))

    if not coords:
        return []

    # Close the ring by appending first coordinate if necessary
    if close_ring:
        if coords[0] != coords[-1]:
            coords.append(coords[0])

    return coords

def lon_to_360(lon):
    """
    Converts longitude in -180/180E convention to 0/360E convention.
    
    Parameters
    ----------
    lon : array
        longitude array in -180/180E convention

    Returns
    -------
    lon_360 : array of floats in 0/360E convention

    """
    return np.where(lon < 0, lon + 360, lon)