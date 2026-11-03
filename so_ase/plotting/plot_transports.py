# so_ase/plotting/plot_transports.py

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
from matplotlib.collections import PolyCollection, LineCollection
from .plotting_maps import create_map


def _unwrap_lon(lon):
    """Shifts each row of longitudes by multiples of 360 so that it is continuous (antimeridian safe)."""
    ref = lon[..., :1]
    return ref + ((lon - ref + 180) % 360 - 180)


def plot_section_profile(ds, var='normal_velocity', time=None, ax=None, vmax=None, cmap='RdBu_r'):
    """
    Plots a variable of a section transport dataset as a vertical section:
    x = distance along the section (km), y = depth (m, increasing downward), colour = `var`.

    Parameters
    ----------
    ds : xarray.Dataset
        Output of `fesom_section_transport`.
    var : str, optional
        'normal_velocity' (positive to the right of travel), 'u' (eastward) or 'v' (northward).
        Default is 'normal_velocity'.
    time : optional
        Time to plot (value accepted by `.sel(time=..., method='nearest')`). Default is the first record.
    ax : matplotlib.axes.Axes, optional
        Axes to draw into. A new figure is created if None.
    vmax : float, optional
        Symmetric colour limit. Default is the 99th percentile of |var|.
    cmap : str, optional
        Colormap name. Default is 'RdBu_r'.

    Returns
    -------
    fig, ax : matplotlib Figure and Axes
    """
    da = ds[var]
    da = da.isel(time=0) if time is None else da.sel(time=time, method='nearest')
    da = da.sortby('dist')
    values = da.values
    if vmax is None:
        vmax = np.nanpercentile(np.abs(values), 99)

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 4))
    else:
        fig = ax.figure
    pm = ax.pcolormesh(da.dist.values, da.nz.values, np.ma.masked_invalid(values),
                       cmap=cmap, vmin=-vmax, vmax=vmax, shading='auto')
    ax.invert_yaxis()
    ax.set_xlabel('distance along section (km)')
    ax.set_ylabel('depth (m)')
    cb = fig.colorbar(pm, ax=ax)
    cb.set_label(f"{var} ({da.attrs.get('units', '')})")
    ax.set_title(f"{var}, {str(da.time.values)[:10]}")
    return fig, ax


def plot_section_map(ds, time=None, ax=None, projection=None, extent=None, vmax=None, cmap='RdBu_r'):
    """
    Plots the section, the mesh triangles used for the transport calculation and the crossed
    edges on a map. Triangles are coloured by their depth-integrated transport (Sv).

    Parameters
    ----------
    ds : xarray.Dataset
        Output of `fesom_section_transport`.
    time : optional
        Time record to use for the triangle transports. Default None uses the time mean.
    ax : cartopy GeoAxes, optional
        Axes to draw into. Created with `projection` if None.
    projection : cartopy.crs projection, optional
        Default is SouthPolarStereo for sections south of 40S, PlateCarree otherwise.
    extent : list of float, optional
        Map extent [lon_min, lon_max, lat_min, lat_max] to zoom into (e.g. to resolve single triangles).
        Default is the section bounding box, or the Southern Ocean for ring-like sections.
    vmax : float, optional
        Symmetric colour limit in Sv. Default is the 99th percentile of |transport|.
    cmap : str, optional
        Colormap name. Default is 'RdBu_r'.

    Returns
    -------
    fig, ax : matplotlib Figure and GeoAxes
    """
    if projection is None:
        projection = ccrs.SouthPolarStereo() if float(ds.path_lat.max()) < -40 else ccrs.PlateCarree()
    polar = isinstance(projection, ccrs.SouthPolarStereo)
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 7), subplot_kw=dict(projection=projection))
    else:
        fig = ax.figure

    lon_u = _unwrap_lon(ds.path_lon.values[None, :])[0]
    if extent is not None:
        create_map(ax, extent=list(extent))
    elif polar and lon_u.max() - lon_u.min() > 300:
        create_map(ax, extent='southern_ocean', circular=True)
    else:
        pad = 3.0
        bbox = [float(lon_u.min() - pad), float(lon_u.max() + pad),
                float(max(ds.path_lat.min() - pad, -90)), float(min(ds.path_lat.max() + pad, 90))]
        create_map(ax, extent=bbox)

    tr = ds.transport.sel(time=time, method='nearest') if time is not None else ds.transport.mean('time')
    val = tr.sum('nz').values / 1e6
    if vmax is None:
        vmax = np.nanpercentile(np.abs(val), 99)

    verts = np.stack([_unwrap_lon(ds.vertex_lon.values), ds.vertex_lat.values], axis=-1)
    pc = PolyCollection(verts, array=val, cmap=cmap, norm=mcolors.Normalize(-vmax, vmax),
                        edgecolors='0.25', linewidths=0.4, zorder=5)
    pc.set_transform(ccrs.PlateCarree()._as_mpl_transform(ax))
    ax.add_collection(pc)

    edges = np.stack([_unwrap_lon(ds.edge_lon.values), ds.edge_lat.values], axis=-1)
    lc = LineCollection(edges, colors='0.1', linewidths=0.8, zorder=6)
    lc.set_transform(ccrs.PlateCarree()._as_mpl_transform(ax))
    ax.add_collection(lc)

    path_lon = np.degrees(np.unwrap(np.radians(ds.path_lon.values)))
    ax.plot(path_lon, ds.path_lat.values, transform=ccrs.PlateCarree(), color='k', lw=1.5, zorder=7)

    cb = fig.colorbar(pc, ax=ax, shrink=0.7)
    cb.set_label('triangle transport (Sv)')
    label = 'time mean' if time is None else str(ds.time.sel(time=time, method='nearest').values)[:10]
    ax.set_title(f"section triangles and crossed edges ({label})")
    return fig, ax


def plot_transport_timeseries(ds, ax=None):
    """
    Plots the net and gross section transport (Sv) for every time record.

    Parameters
    ----------
    ds : xarray.Dataset
        Output of `fesom_section_transport`.
    ax : matplotlib.axes.Axes, optional
        Axes to draw into. A new figure is created if None.

    Returns
    -------
    fig, ax : matplotlib Figure and Axes
    """
    net = ds.transport.sum(('nz', 'elem')).values / 1e6
    gross = ds.transport_gross.sum('nz').values / 1e6
    x = np.arange(len(net))

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 3.5))
    else:
        fig = ax.figure
    ax.plot(x, net, 'o-', label='net')
    ax.plot(x, gross, 's--', color='0.5', label='gross (sum of |contributions|)')
    ax.axhline(0, color='k', lw=0.5)
    ax.set_xticks(x, [t.strftime('%b') for t in ds.time.values])
    ax.set_ylabel('transport (Sv)')
    ax.legend(frameon=False)
    return fig, ax
