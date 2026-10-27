# so_ase/fesom2/eval_transport.py

import os
import numpy as np
import xarray as xr
import pyfesom2 as pf
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree
from ..miscellaneous.helpers_misc import lon_lat_to_cartesian

R_EARTH = 6371000.0
EULER_ANGLES = (50, 15, -90)


###--------> 1. Mesh geometry

def load_section_mesh(mesh_diag_path):
    """
    Loads the FESOM2 mesh geometry needed for section transports from `fesom.mesh.diag.nc`.

    All indices are converted to 0-based. Edge-to-element links that are NaN (boundary edges)
    are set to -1. The vectors in `edge_cross_dxdy` are stored in the rotated frame of the model
    and are rotated to geographical components in `fesom_section_transport`.

    Parameters
    ----------
    mesh_diag_path : str
        Directory containing `fesom.mesh.diag.nc`.

    Returns
    -------
    dict
        Mesh arrays, unit node vectors and a KD-tree over the nodes.
    """
    ds = xr.open_dataset(os.path.join(mesh_diag_path, 'fesom.mesh.diag.nc'))

    edge_elems = ds.edge_face_links.values.T
    edge_elems = np.where(np.isnan(edge_elems), -1, np.nan_to_num(edge_elems) - 1).astype(np.int64)
    face_nodes = (ds.face_nodes.values.T - 1).astype(np.int64)
    edge_nodes = (ds.edge_nodes.values.T - 1).astype(np.int64)

    e_idx = np.repeat(np.arange(len(edge_elems)), 2)
    t_idx = edge_elems.ravel()
    m = t_idx >= 0
    elem_edges = e_idx[m][np.argsort(t_idx[m], kind='stable')].reshape(len(face_nodes), 3)

    X = np.column_stack(lon_lat_to_cartesian(ds.lon.values, ds.lat.values, R=1))
    mesh = dict(
        X=X,
        lon=ds.lon.values,
        lat=ds.lat.values,
        face_nodes=face_nodes,
        edge_nodes=edge_nodes,
        edge_elems=edge_elems,
        elem_edges=elem_edges,
        edge_cross=ds.edge_cross_dxdy.values.T,
        nod_in_elem=ds.nod_in_elem2D.values.T - 1,
        nod_in_elem_num=ds.nod_in_elem2D_num.values,
        ulevels=ds.ulevels.values,
        nlevels=ds.nlevels.values,
        nz=ds.nz.values,
        tree=cKDTree(X),
    )
    return mesh


###--------> 2. Section geometry and element crossings

def _projector(c):
    """Oblique stereographic projection (unit sphere) from the antipode of c onto the tangent plane at c."""
    a = np.eye(3)[np.argmin(np.abs(c))]
    e1 = a - (a @ c) * c
    e1 /= np.linalg.norm(e1)
    E = np.column_stack([e1, np.cross(c, e1)])
    return lambda X: 2 * (X @ E) / (1 + X @ c)[:, None]


def _orient_waypoints(lon, lat, closed=None):
    """
    Returns waypoints and the centre vector c of the section. Closed sections are oriented
    counterclockwise seen from c, so that the region around c lies on the left and positive
    transport is outward from it.
    """
    W = np.column_stack(lon_lat_to_cartesian(lon, lat, R=1))
    if closed is None:
        closed = np.linalg.norm(W[0] - W[-1]) < 1e-9
    if closed:
        if np.linalg.norm(W[0] - W[-1]) > 1e-9:
            lon, lat = np.append(lon, lon[0]), np.append(lat, lat[0])
            W = np.column_stack(lon_lat_to_cartesian(lon, lat, R=1))
        c = W[:-1].sum(0)
        c /= np.linalg.norm(c)
        if np.sum(np.cross(W[:-1], W[1:]) @ c) < 0:
            lon, lat, W = lon[::-1], lat[::-1], W[::-1]
    else:
        c = W.sum(0)
        c /= np.linalg.norm(c)
    return lon, lat, closed, c


def _densify(lon, lat, max_step_km):
    """Great-circle densification of waypoints. Returns unit vectors and cumulative arc length (km)."""
    W = np.column_stack(lon_lat_to_cartesian(lon, lat, R=1))
    pts = []
    for A, B in zip(W[:-1], W[1:]):
        th = np.arccos(np.clip(A @ B, -1, 1))
        n = max(int(np.ceil(th * R_EARTH / 1e3 / max_step_km)), 1)
        t = np.arange(n) / n
        pts.append((np.sin((1 - t) * th)[:, None] * A + np.sin(t * th)[:, None] * B) / np.sin(th))
    Xd = np.vstack(pts + [W[-1:]])
    seg = 2 * np.arcsin(np.linalg.norm(np.diff(Xd, axis=0), axis=1) / 2) * R_EARTH / 1e3
    return Xd, np.concatenate([[0.0], np.cumsum(seg)])


def _in_triangle(tri, p, tol=1e-9):
    (ax, ay), (bx, by), (cx, cy) = tri
    d = (by - cy) * (ax - cx) + (cx - bx) * (ay - cy)
    l1 = ((by - cy) * (p[0] - cx) + (cx - bx) * (p[1] - cy)) / d
    l2 = ((cy - ay) * (p[0] - cx) + (ax - cx) * (p[1] - cy)) / d
    return min(l1, l2, 1 - l1 - l2) >= -tol


def _locate(mesh, x, proj):
    """Element containing unit vector x, searched among the elements of the nearest nodes."""
    _, nodes = mesh['tree'].query(x, k=12)
    p = proj(x[None])[0]
    for n in nodes:
        for T in mesh['nod_in_elem'][n, :mesh['nod_in_elem_num'][n]]:
            if _in_triangle(proj(mesh['X'][mesh['face_nodes'][T]]), p):
                return T
    raise ValueError('Section start point is not inside any mesh element (land or outside the mesh).')


def _seg_hit(p, q, A, B):
    """Parameters (s along p->q, u along A->B) of the intersection of two segments, or None if parallel."""
    r, w, ap = q - p, B - A, A - p
    den = r[0] * w[1] - r[1] * w[0]
    if abs(den) < 1e-15:
        return None
    return (ap[0] * w[1] - ap[1] * w[0]) / den, (ap[0] * r[1] - ap[1] * r[0]) / den


def _walk(mesh, proj, Xd):
    """Walks element to element along the densified section; returns element sequence and crossed edges."""
    P = proj(Xd)
    T = _locate(mesh, Xd[0], proj)
    seq, edges, segs, fracs = [T], [], [], []
    for i in range(len(P) - 1):
        p, q = P[i], P[i + 1]
        entry = -1
        for _ in range(100):
            ed = mesh['elem_edges'][T]
            ends = proj(mesh['X'][mesh['edge_nodes'][ed]].reshape(-1, 3)).reshape(3, 2, 2)
            best = None
            for k, e in enumerate(ed):
                if e == entry:
                    continue
                h = _seg_hit(p, q, ends[k, 0], ends[k, 1])
                if h and 1e-9 < h[0] <= 1 + 1e-9 and -1e-9 <= h[1] <= 1 + 1e-9 and (best is None or h[0] < best[0]):
                    best = (h[0], e)
            if best is None:
                break
            s, e = best
            ee = mesh['edge_elems'][e]
            nxt = ee[1] if ee[0] == T else ee[0]
            if nxt < 0:
                raise ValueError('Section leaves the mesh through a coastline (adjust the waypoints).')
            seq.append(nxt)
            edges.append(e)
            segs.append(i)
            fracs.append(s)
            p = p + s * (q - p)
            T, entry = nxt, e
        else:
            raise RuntimeError('Element walk did not terminate.')
    return np.array(seq), np.array(edges), np.array(segs), np.array(fracs)


def find_section_crossings(mesh, waypoints, max_step_km=5.0, closed=None):
    """
    Finds the ordered sequence of mesh elements and edges crossed by a section.

    Waypoints are connected by great-circle arcs, densified to `max_step_km`, and walked
    element-to-element in a stereographic projection centred on the section.

    Parameters
    ----------
    mesh : dict
        Mesh returned by `load_section_mesh`.
    waypoints : array-like
        List of (lon, lat) pairs in degrees.
    max_step_km : float, optional
        Maximum spacing of densified section points. Default is 5 km.
    closed : bool or None, optional
        If None, the section is closed when the last waypoint equals the first one.

    Returns
    -------
    dict
        'seq': element indices visited (len n+1), 'edge': crossed edge per crossing (len n),
        'dist_km': distance along the section at each crossing (len n), 'sigma': +1 if the
        edge's A->B direction points to the right of the travel direction, else -1,
        'total_km': section length, 'closed': bool.
    """
    wp = np.asarray(waypoints, dtype=float)
    lon, lat, closed, c = _orient_waypoints(wp[:, 0], wp[:, 1], closed)
    Xd, cum = _densify(lon, lat, max_step_km)
    seq, edges, segs, fracs = _walk(mesh, _projector(c), Xd)

    dist = cum[segs] + fracs * (cum[segs + 1] - cum[segs])
    d = Xd[segs + 1] - Xd[segs]
    xc = Xd[segs] + fracs[:, None] * d
    xc /= np.linalg.norm(xc, axis=1)[:, None]
    en = mesh['edge_nodes'][edges]
    sigma = np.sign(np.einsum('ij,ij->i', mesh['X'][en[:, 1]] - mesh['X'][en[:, 0]], np.cross(d, xc)))
    path_lon = np.degrees(np.arctan2(Xd[:, 1], Xd[:, 0]))
    path_lat = np.degrees(np.arcsin(np.clip(Xd[:, 2], -1, 1)))
    return dict(seq=seq, edge=edges, dist_km=dist, sigma=sigma, total_km=cum[-1], closed=closed,
                path_lon=path_lon, path_lat=path_lat, path_dist=cum)


def _section_elements(mesh, sec, euler):
    """
    Builds, for every element visited by the section, the summed signed dual-face area vector
    (east, north; m) and the total half-segment length (m), plus the distance along the section.
    """
    seq, edges, sigma = sec['seq'], sec['edge'], sec['sigma']
    n = len(edges)

    en = mesh['edge_nodes'][edges]
    A, B = en[:, 0], en[:, 1]
    M = mesh['X'][A] + mesh['X'][B]
    M /= np.linalg.norm(M, axis=1)[:, None]
    mlon = np.degrees(np.arctan2(M[:, 1], M[:, 0]))
    mlat = np.degrees(np.arcsin(M[:, 2]))
    dlon = (mesh['lon'][B] - mesh['lon'][A] + 180) % 360 - 180
    eE = R_EARTH * np.cos(np.radians(mlat)) * np.radians(dlon)
    eN = R_EARTH * np.radians(mesh['lat'][B] - mesh['lat'][A])

    S = []
    for h in (0, 1):
        D = mesh['edge_cross'][edges, 2 * h:2 * h + 2]
        De, Dn = pf.vec_rotate_r2g(*euler, mlon, mlat, D[:, 0], D[:, 1], 1)
        s = np.sign(Dn * eE - De * eN)
        S.append(np.column_stack([s * Dn, -s * De]))

    ee = mesh['edge_elems'][edges]
    hf = (ee[:, 0] != seq[:-1]).astype(int)
    ht = (ee[:, 0] != seq[1:]).astype(int)
    Sf = np.where(hf[:, None] == 0, S[0], S[1]) * sigma[:, None]
    St = np.where(ht[:, None] == 0, S[0], S[1]) * sigma[:, None]
    half_elem = np.concatenate([seq[:-1], seq[1:]])
    half_S = np.vstack([Sf, St])
    half_L = np.linalg.norm(half_S, axis=1)

    uniq, first = np.unique(seq, return_index=True)
    order = np.argsort(first)
    elem_out = uniq[order]
    rank = np.empty(len(uniq), dtype=int)
    rank[order] = np.arange(len(uniq))
    half_out = rank[np.searchsorted(uniq, half_elem)]
    agg = csr_matrix((np.ones(len(half_out)), (half_out, np.arange(len(half_out)))), shape=(len(uniq), len(half_out)))
    S_out = agg @ half_S
    L_out = agg @ half_L

    entry, exit_ = {seq[0]: 0.0}, {}
    for k in range(n):
        exit_.setdefault(seq[k], sec['dist_km'][k])
        entry.setdefault(seq[k + 1], sec['dist_km'][k])
    exit_.setdefault(seq[-1], sec['total_km'])
    dist_out = np.array([0.5 * (entry[T] + exit_[T]) for T in elem_out])

    C = mesh['X'][mesh['face_nodes'][elem_out]].mean(1)
    C /= np.linalg.norm(C, axis=1)[:, None]
    lon_c = np.degrees(np.arctan2(C[:, 1], C[:, 0]))
    lat_c = np.degrees(np.arcsin(C[:, 2]))

    return dict(elem_idx=elem_out, S=S_out, L=L_out, dist_km=dist_out, lon=lon_c, lat=lat_c,
                half_out=half_out, half_S=half_S)


###--------> 3. Transport

def fesom_section_transport(src_path, mesh_diag_path, waypoints, years=(1979, 2015), euler=(50, 15, -90),
                            unrotate_velocities=True, max_step_km=5.0, closed=None, log=True):
    """
    Computes the volume transport across an arbitrary section of the FESOM2 mesh from monthly
    u and v output, using the dual (median-dual) faces of the crossed edges.

    For every crossed edge, the flux through its dual face is evaluated from the two half-segments
    (midpoint to element centroid) using the velocity of the element that contains each half-segment.
    Contributions are summed per element and vertical level.

    Sign convention: positive transport is flow to the right of the travel direction along the
    section. For closed sections the waypoints are oriented counterclockwise seen from the region
    they enclose, so positive is outward from that region.

    Parameters
    ----------
    src_path : str
        Directory containing `u.fesom.YYYY.nc` and `v.fesom.YYYY.nc`.
    mesh_diag_path : str
        Directory containing `fesom.mesh.diag.nc`.
    waypoints : array-like
        List of (lon, lat) pairs in degrees defining the section.
    years : tuple of int, optional
        First and last year to process. Default is (1979, 2015).
    euler : tuple of float, optional
        Euler angles (alpha, beta, gamma) of the rotated frame of `edge_cross_dxdy` and of the u/v output.
        Default is (50, 15, -90).
    unrotate_velocities : bool, optional
        If True, the u/v output is in the rotated frame and is rotated to geographical components at the
        element centroids using `euler`. Set False if u/v are already geographical. Default True.
    max_step_km : float, optional
        Maximum spacing of densified section points. Default is 5 km.
    closed : bool or None, optional
        Force open/closed interpretation. Default None (closed if first and last waypoint match).
    log : bool, optional
        Print progress information. Default is True.

    Returns
    -------
    xarray.Dataset
        'transport' (time, nz, elem) in m3/s, 'normal_velocity' (time, nz, elem) in m/s,
        'transport_gross' (time, nz) in m3/s (sum of absolute half-segment contributions).
        Coordinates: elem_idx, lon, lat, dist (km at the element's crossing points), nz (mid-level depth),
        depth_top, depth_bottom, layer_thickness.
    """
    mesh = load_section_mesh(mesh_diag_path)
    if log:
        print(f"Mesh loaded: {mesh_diag_path}", flush=True)

    wp = np.asarray(waypoints, dtype=float)
    sec = find_section_crossings(mesh, wp, max_step_km=max_step_km, closed=closed)
    if log:
        print(f"Section: {len(sec['edge'])} crossed edges, {len(sec['seq']) - 1} element steps, "
              f"length {sec['total_km']:.1f} km, closed={sec['closed']}", flush=True)
    el = _section_elements(mesh, sec, euler)
    elem = el['elem_idx']

    files_u = [os.path.join(src_path, f"u.fesom.{y}.nc") for y in range(years[0], years[1] + 1)]
    files_v = [os.path.join(src_path, f"v.fesom.{y}.nc") for y in range(years[0], years[1] + 1)]
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
    ds_u = xr.open_mfdataset(files_u, combine='by_coords', decode_times=time_coder)
    ds_v = xr.open_mfdataset(files_v, combine='by_coords', decode_times=time_coder)
    u = ds_u.u.isel(elem=elem).load().values.astype(float)
    v = ds_v.v.isel(elem=elem).load().values.astype(float)
    if log:
        print(f"Loaded u/v for {len(elem)} elements, {u.shape[0]} time steps", flush=True)

    nlay = u.shape[2]
    j = np.arange(nlay)
    active = (j[None, :] >= mesh['ulevels'][elem][:, None] - 1) & (j[None, :] <= mesh['nlevels'][elem][:, None] - 2)
    u = np.where(active[None], np.nan_to_num(u), 0.0)
    v = np.where(active[None], np.nan_to_num(v), 0.0)
    if unrotate_velocities:
        shape = u.shape
        lon2 = np.broadcast_to(el['lon'][None, :, None], shape).ravel()
        lat2 = np.broadcast_to(el['lat'][None, :, None], shape).ravel()
        u, v = pf.vec_rotate_r2g(*euler, lon2, lat2, u.ravel(), v.ravel(), 1)
        u, v = np.reshape(u, shape), np.reshape(v, shape)
    u_out, v_out = np.where(active[None], u, np.nan), np.where(active[None], v, np.nan)

    corner = mesh['face_nodes'][elem]
    nodes = np.unique(corner)
    pos = np.searchsorted(nodes, corner)

    def elem_mean(da):
        nv = da.isel(nod2=nodes).load().values.astype(float)[:, pos, :]
        s, c = np.nansum(nv, axis=2), np.sum(~np.isnan(nv), axis=2)
        return np.where(c > 0, s / np.maximum(c, 1), np.nan)

    files_t = [os.path.join(src_path, f"temp.fesom.{y}.nc") for y in range(years[0], years[1] + 1)]
    files_s = [os.path.join(src_path, f"salt.fesom.{y}.nc") for y in range(years[0], years[1] + 1)]
    ds_t = xr.open_mfdataset(files_t, combine='by_coords', decode_times=time_coder)
    ds_s = xr.open_mfdataset(files_s, combine='by_coords', decode_times=time_coder)
    temp_out = np.where(active[None], elem_mean(ds_t.temp), np.nan)
    salt_out = np.where(active[None], elem_mean(ds_s.salt), np.nan)
    if log:
        print(f"Loaded temp/salt at {len(nodes)} nodes", flush=True)

    dz = np.diff(mesh['nz'])

    flux = u * el['S'][None, :, 0, None] + v * el['S'][None, :, 1, None]
    normal_velocity = np.where(active[None], flux / el['L'][None, :, None], np.nan)
    transport = flux * dz[None, None, :]

    uh, vh = u[:, el['half_out'], :], v[:, el['half_out'], :]
    gross = np.abs(uh * el['half_S'][None, :, 0, None] + vh * el['half_S'][None, :, 1, None]).sum(1) * dz[None, :]

    ds = xr.Dataset(
        data_vars=dict(
            transport=(('time', 'nz', 'elem'), np.transpose(transport, (0, 2, 1)),
                       dict(long_name='volume transport across section', units='m3 s-1')),
            normal_velocity=(('time', 'nz', 'elem'), np.transpose(normal_velocity, (0, 2, 1)),
                             dict(long_name='normal velocity (positive to the right of travel)', units='m s-1')),
            transport_gross=(('time', 'nz'), gross,
                             dict(long_name='sum of absolute half-segment transports', units='m3 s-1')),
            u=(('time', 'nz', 'elem'), np.transpose(u_out, (0, 2, 1)),
               dict(long_name='eastward velocity at element (geographical)', units='m s-1')),
            v=(('time', 'nz', 'elem'), np.transpose(v_out, (0, 2, 1)),
               dict(long_name='northward velocity at element (geographical)', units='m s-1')),
            temp=(('time', 'nz', 'elem'), np.transpose(temp_out, (0, 2, 1)),
                  dict(long_name='temperature at element (mean of corner nodes)', units='degC')),
            salt=(('time', 'nz', 'elem'), np.transpose(salt_out, (0, 2, 1)),
                  dict(long_name='practical salinity at element (mean of corner nodes)', units='psu')),
        ),
        coords=dict(
            time=ds_u.time.values,
            nz=ds_u.nz.values,
            elem=np.arange(len(elem)),
            elem_idx=('elem', elem),
            lon=('elem', el['lon'], dict(units='degrees_east')),
            lat=('elem', el['lat'], dict(units='degrees_north')),
            dist=('elem', el['dist_km'], dict(long_name='distance along section', units='km')),
            path_lon=('path', sec['path_lon'], dict(units='degrees_east')),
            path_lat=('path', sec['path_lat'], dict(units='degrees_north')),
            path_dist=('path', sec['path_dist'], dict(long_name='distance along densified section', units='km')),
            vertex_lon=(('elem', 'vertex'), mesh['lon'][mesh['face_nodes'][elem]], dict(units='degrees_east')),
            vertex_lat=(('elem', 'vertex'), mesh['lat'][mesh['face_nodes'][elem]], dict(units='degrees_north')),
            crossing_dist=('crossing', sec['dist_km'], dict(long_name='distance along section at crossing', units='km')),
            edge_lon=(('crossing', 'endpoint'), mesh['lon'][mesh['edge_nodes'][sec['edge']]], dict(units='degrees_east')),
            edge_lat=(('crossing', 'endpoint'), mesh['lat'][mesh['edge_nodes'][sec['edge']]], dict(units='degrees_north')),
            depth_top=('nz', mesh['nz'][:-1], dict(units='m')),
            depth_bottom=('nz', mesh['nz'][1:], dict(units='m')),
            layer_thickness=('nz', dz, dict(units='m')),
        ),
        attrs=dict(
            sign_convention='positive = flow to the right of travel; closed sections oriented CCW around enclosed region (positive = outward)',
            closed=int(sec['closed']),
            section_length_km=float(sec['total_km']),
            n_crossed_edges=len(sec['edge']),
        ),
    )
    ds_u.close()
    ds_v.close()
    ds_t.close()
    ds_s.close()
    return ds

