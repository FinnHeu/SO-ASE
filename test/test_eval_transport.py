# test/test_eval_transport.py
#
# Tests for so_ase.fesom2.eval_transport (section volume transports on the FESOM2 mesh)
# and so_ase.plotting.plot_transports. Run with:  pytest test/test_eval_transport.py -v
# The tests read the FESOM2 example output and mesh diagnostics from the paths below.

import functools

import numpy as np
import pytest
import xarray as xr
from pyproj import Transformer, Geod

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import so_ase as so
from so_ase.fesom2.eval_transport import (
    load_section_mesh, find_section_crossings, _section_elements, EULER_ANGLES, R_EARTH,
)

SRC = '/work/ba1550/a270186/simulations/awiesm3-v3.4.2-CAV-ICB/production/CAV-INI/outdata/fesom/'
MESH = '/work/ab0995/a270186/model_inputs/fesom2/mesh/DARS2cav/'
YEARS = (1853, 1853)

DRAKE = [(-57.5, -61.8), (-66.5, -57.0)]

to3031 = Transformer.from_crs('EPSG:4326', 'EPSG:3031', always_xy=True)
to4326 = Transformer.from_crs('EPSG:3031', 'EPSG:4326', always_xy=True)
geod = Geod(ellps='sphere')


def antarctic_ring(lat0=-60.0, amp_km=0.0, n_waves=0, n_pts=360):
    """Closed ring in Antarctic stereographic coordinates with r(theta) = r0 + amp * sin(n theta)."""
    x0, y0 = to3031.transform(0.0, lat0)
    r0 = np.hypot(x0, y0)
    th = np.linspace(0, 2 * np.pi, n_pts, endpoint=False)
    r = r0 + 1e3 * amp_km * np.sin(n_waves * th)
    lon, lat = to4326.transform(r * np.cos(th), r * np.sin(th))
    return np.vstack([np.column_stack([lon, lat]), [lon[0], lat[0]]])


def geodesic_circle(lon_c, lat_c, radius_km, n_pts=72):
    """Closed small circle of geodesic radius radius_km around (lon_c, lat_c)."""
    brg = np.linspace(0, 360, n_pts, endpoint=False)
    lon, lat, _ = geod.fwd(np.full(n_pts, lon_c), np.full(n_pts, lat_c), brg, np.full(n_pts, radius_km * 1e3))
    return np.vstack([np.column_stack([lon, lat]), [lon[0], lat[0]]])


@functools.lru_cache(maxsize=1)
def mesh():
    return load_section_mesh(MESH)


def transport(waypoints):
    return so.fesom_section_transport(SRC, MESH, waypoints, years=YEARS, log=False)


def annual_net_Sv(ds):
    return float(ds.transport.mean('time').sum(('nz', 'elem')) / 1e6)


def annual_gross_Sv(ds):
    return float(ds.transport_gross.mean('time').sum('nz') / 1e6)


def max_monthly_net_Sv(ds):
    return float(np.abs(ds.transport.sum(('nz', 'elem')).values).max() / 1e6)


@pytest.fixture(scope='module')
def drake_ds():
    return transport(DRAKE)


def test_streamfunction_flow_has_zero_net_flux():
    """A non-divergent flow (streamfunction psi = sin(3 lon) cos^2(lat)) must give ~zero net flux through a closed ring."""
    m = mesh()
    wp = antarctic_ring()
    sec = find_section_crossings(m, wp)
    el = _section_elements(m, sec, EULER_ANGLES)
    lam, phi = np.radians(el['lon']), np.radians(el['lat'])
    dpsi_dlam = 3 * np.cos(3 * lam) * np.cos(phi) ** 2 * 1e6
    dpsi_dphi = np.sin(3 * lam) * (-2 * np.cos(phi) * np.sin(phi)) * 1e6
    u, v = -dpsi_dphi / R_EARTH, dpsi_dlam / (R_EARTH * np.cos(phi))
    net = np.sum(u * el['S'][:, 0] + v * el['S'][:, 1])
    gross = np.abs(u[el['half_out']] * el['half_S'][:, 0] + v[el['half_out']] * el['half_S'][:, 1]).sum()
    assert abs(net) / gross < 1e-6


def test_circle_ring_closes():
    ds = transport(antarctic_ring())
    assert bool(ds.attrs['closed'])
    assert max_monthly_net_Sv(ds) < 1.0
    assert abs(annual_net_Sv(ds)) / annual_gross_Sv(ds) < 1e-4


@pytest.mark.parametrize('amp_km, n_waves', [(200, 6), (300, 3)])
def test_sine_wave_ring_closes(amp_km, n_waves):
    ds = transport(antarctic_ring(amp_km=amp_km, n_waves=n_waves))
    assert bool(ds.attrs['closed'])
    assert max_monthly_net_Sv(ds) < 1.0
    assert abs(annual_net_Sv(ds)) / annual_gross_Sv(ds) < 1e-4


def test_unrotate_velocities_flag_changes_closure():
    """Treating the u/v output as already geographical (flag False) breaks the closure of the 60S ring."""
    ds = so.fesom_section_transport(SRC, MESH, antarctic_ring(), years=YEARS, unrotate_velocities=False, log=False)
    assert abs(annual_net_Sv(ds)) > 10.0


def test_small_closed_loop_closes():
    ds = transport(geodesic_circle(-120.0, -60.0, 300))
    assert bool(ds.attrs['closed'])
    assert max_monthly_net_Sv(ds) < 0.1
    assert abs(annual_net_Sv(ds)) / annual_gross_Sv(ds) < 1e-4


def test_reversed_section_flips_sign(drake_ds):
    rev = transport(DRAKE[::-1])
    assert annual_net_Sv(rev) == pytest.approx(-annual_net_Sv(drake_ds), abs=1e-6)


def test_antimeridian_longitude_representation_invariant():
    a = transport([(179.0, -58.0), (-179.0, -58.0)])
    b = transport([(179.0, -58.0), (181.0, -58.0)])
    assert annual_net_Sv(a) == pytest.approx(annual_net_Sv(b), abs=1e-9)


def test_drake_passage_magnitude_and_direction(drake_ds):
    """Drake Passage transport is eastward and of order 100-200 Sv (published ~134 +- 6 Sv)."""
    net = annual_net_Sv(drake_ds)
    assert 100.0 < net < 200.0


def test_drake_passage_u_and_v_stored(drake_ds):
    assert drake_ds.u.dims == ('time', 'nz', 'elem')
    assert drake_ds.v.dims == ('time', 'nz', 'elem')
    assert bool((drake_ds.u.isnull() == drake_ds.normal_velocity.isnull()).all())


def test_plot_functions_run(drake_ds):
    fig, ax = so.plot_section_profile(drake_ds, var='u', time=drake_ds.time.values[0])
    assert ax is not None
    plt.close(fig)
    fig, ax = so.plot_section_profile(drake_ds, var='normal_velocity')
    assert ax is not None
    plt.close(fig)
    fig, ax = so.plot_section_map(drake_ds, extent=[-68, -56, -64, -54])
    assert ax is not None
    plt.close(fig)
    fig, ax = so.plot_transport_timeseries(drake_ds)
    assert ax is not None
    plt.close(fig)


def test_drake_temp_salt_stored(drake_ds):
    for name in ('temp', 'salt'):
        assert drake_ds[name].dims == ('time', 'nz', 'elem')
        assert bool((drake_ds[name].isnull() == drake_ds.u.isnull()).all())


def test_drake_temp_salt_plausible(drake_ds):
    t = drake_ds.temp.values
    s = drake_ds.salt.values
    assert -2.5 < np.nanmin(t) and np.nanmax(t) < 10.0
    # Lower bound 31 psu allows for coastal surface meltwater in austral summer (~31.5 psu near the South Shetlands).
    assert 31.0 < np.nanmin(s) and np.nanmax(s) < 35.0


def test_temp_salt_equal_corner_node_mean(drake_ds):
    """Element temperature and salinity equal the mean of the 3 corner nodes read directly from the files."""
    m = mesh()
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
    ds_t = xr.open_dataset(SRC + f'temp.fesom.{YEARS[0]}.nc', decode_times=time_coder)
    ds_s = xr.open_dataset(SRC + f'salt.fesom.{YEARS[0]}.nc', decode_times=time_coder)
    k = len(drake_ds.elem) // 2
    elem_idx = int(drake_ds.elem_idx[k])
    corners = m['face_nodes'][elem_idx]
    t_nodes = ds_t.temp.isel(time=0, nod2=corners).load().values.astype(float)
    s_nodes = ds_s.salt.isel(time=0, nod2=corners).load().values.astype(float)
    t_ref = np.nanmean(t_nodes, axis=0)
    s_ref = np.nanmean(s_nodes, axis=0)
    j = np.arange(t_ref.size)
    active = (j >= m['ulevels'][elem_idx] - 1) & (j <= m['nlevels'][elem_idx] - 2)
    t_ds = drake_ds.temp.isel(time=0, elem=k).values
    s_ds = drake_ds.salt.isel(time=0, elem=k).values
    assert np.allclose(t_ds[active], t_ref[active], atol=1e-5)
    assert np.allclose(s_ds[active], s_ref[active], atol=1e-5)
    ds_t.close()
    ds_s.close()


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
