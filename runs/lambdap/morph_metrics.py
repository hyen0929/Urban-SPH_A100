#!/usr/bin/env python3
"""
Urban morphology metrics for the Leipzig case (B&E paper).

Produces, from the building geometry alone (no simulation results required):

  1. Source-local morphology metrics in circular windows of radius R --
     the x-axis of the morphology-deposition correlation:
       lambda_p   plan area density              (built cell fraction)
       lambda_f   frontal area density           (windward step sum / plan area)
       H_built    mean building height over built cells
       H_all      mean height over all window cells
       W_char     characteristic street width    (2 x mean distance-to-building)
       aspect     H_built / W_char               (canyon aspect ratio proxy)
     with an R-sensitivity sweep and the percentile rank of each source
     within the built-up-area distribution.

  2. A stratification verdict: whether the source set spans the local
     lambda_p range widely enough, and evenly enough by tercile, to carry
     a correlation.

  3. Available surface area by type (ground / wall / roof), globally and per
     source window -- the denominator needed to report deposition as an
     area-normalised flux instead of a raw particle count.

Two independent estimators are computed for each surface area and printed
side by side; their agreement is the validation of both.

  point-based   local surface density pi*r_k^2/k per point, summed by class
                (class from z vs. the local column maximum)
  field-based   roof = raster plan area of built cells
                wall = sum of |dH| over cell faces (staircase perimeter x height)

No surface normals are required, so IBM surface-particle clouds without
normal vectors are fully supported.

Metric definitions follow Grimmond & Oke, J. Appl. Meteorol. 38 (1999)
1262-1292, already cited in the manuscript.

Usage
-----
    python3 morph_metrics.py --geom building_0stp.vtk --sources sources_xyz.csv \
        --cell 2.5 --radii 75,150,250 --wind 0 --out results/

    python3 morph_metrics.py --selftest        # synthetic geometry, no data needed

Inputs
------
--geom    *.vtk   legacy VTK POLYDATA / UNSTRUCTURED_GRID point set (binary or
                  ASCII).  Only POINTS are read.
          *.stl   binary or ASCII STL
          *.csv   point cloud, columns x,y,z (header optional)
          *.geojson  footprint polygons with a height property
--sources CSV with columns id,x,y[,z][,label], same coordinate frame as --geom

The geometry bounding box is printed first; check it against the expected
domain before trusting any number below it.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import struct
import sys

import numpy as np

# ---------------------------------------------------------------------------
# geometry loading
# ---------------------------------------------------------------------------

_VTK_DTYPE = {
    "float": ">f4", "double": ">f8", "int": ">i4", "unsigned_int": ">u4",
    "long": ">i8", "short": ">i2", "unsigned_short": ">u2",
    "char": ">i1", "unsigned_char": ">u1",
}


def load_vtk_points(path: str) -> np.ndarray:
    """Legacy VTK (binary or ASCII) -> (n, 3) point array. POINTS only."""
    with open(path, "rb") as fh:
        raw = fh.read()

    m = re.search(rb"(?m)^POINTS[ \t]+(\d+)[ \t]+(\w+)[ \t]*$", raw)
    if m is None:
        raise ValueError(f"{path}: no POINTS section found")
    n = int(m.group(1))
    typ = m.group(2).decode()
    if typ not in _VTK_DTYPE:
        raise ValueError(f"{path}: unsupported POINTS type '{typ}'")
    start = raw.index(b"\n", m.start()) + 1

    is_binary = re.search(rb"(?m)^BINARY[ \t]*$", raw[:m.start()]) is not None
    if is_binary:
        dt = np.dtype(_VTK_DTYPE[typ])
        need = n * 3 * dt.itemsize
        if start + need > len(raw):
            raise ValueError(f"{path}: truncated POINTS block")
        pts = np.frombuffer(raw, dtype=dt, count=n * 3, offset=start)
    else:
        txt = raw[start:].decode("ascii", errors="replace")
        pts = np.array(txt.split()[: n * 3], dtype=np.float64)
    return np.asarray(pts, np.float64).reshape(n, 3)


def load_stl(path: str) -> np.ndarray:
    """STL -> (n, 3, 3) triangle vertices."""
    with open(path, "rb") as fh:
        head = fh.read(84)
        n_declared = struct.unpack("<I", head[80:84])[0] if len(head) == 84 else 0
        rest = fh.read()
    if n_declared > 0 and len(rest) >= n_declared * 50:
        dt = np.dtype([("n", "<3f4"), ("v", "<3,3f4"), ("attr", "<u2")])
        d = np.frombuffer(rest[: n_declared * 50], dtype=dt, count=n_declared)
        return d["v"].astype(np.float64)
    verts, cur = [], []
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            t = line.split()
            if not t:
                continue
            if t[0] == "facet":
                cur = []
            elif t[0] == "vertex" and len(t) >= 4:
                cur.append([float(t[1]), float(t[2]), float(t[3])])
            elif t[0] == "endfacet" and len(cur) == 3:
                verts.append(cur)
    if not verts:
        raise ValueError(f"{path}: no facets parsed")
    return np.asarray(verts, float)


def load_points_csv(path: str) -> np.ndarray:
    with open(path, "r", errors="replace") as fh:
        first = fh.readline()
    try:
        float(first.replace(",", " ").split()[0])
        skip = 0
    except (ValueError, IndexError):
        skip = 1
    delim = "," if "," in first else None
    arr = np.atleast_2d(np.loadtxt(path, delimiter=delim, skiprows=skip))
    if arr.shape[1] < 3:
        raise ValueError(f"{path}: need x,y,z columns")
    return arr[:, :3].astype(np.float64)


def load_footprints(path: str):
    with open(path) as fh:
        gj = json.load(fh)
    out = []
    for f in gj.get("features", [gj]):
        props = f.get("properties", {}) or {}
        h = None
        for key in ("height", "Height", "building:height", "h"):
            if key in props:
                try:
                    h = float(str(props[key]).split()[0])
                    break
                except (ValueError, IndexError):
                    pass
        if h is None:
            lv = props.get("building:levels") or props.get("levels")
            h = 3.0 * float(lv) if lv else 10.0
        g = f.get("geometry", {}) or {}
        rings = ([g["coordinates"]] if g.get("type") == "Polygon"
                 else g.get("coordinates", []) if g.get("type") == "MultiPolygon"
                 else [])
        for r in rings:
            if r:
                out.append((np.asarray(r[0], float)[:, :2], h))
    if not out:
        raise ValueError(f"{path}: no polygons parsed")
    return out


# ---------------------------------------------------------------------------
# height field
# ---------------------------------------------------------------------------


class Field:
    """max-z height field on a regular cell grid."""

    def __init__(self, x0, y0, nx, ny, cell):
        self.x0, self.y0, self.cell = float(x0), float(y0), float(cell)
        self.nx, self.ny = int(nx), int(ny)
        self.H = np.full((self.ny, self.nx), -1.0)

    @property
    def extent(self):
        return (self.x0, self.x0 + self.nx * self.cell,
                self.y0, self.y0 + self.ny * self.cell)

    def xc(self):
        return self.x0 + (np.arange(self.nx) + 0.5) * self.cell

    def yc(self):
        return self.y0 + (np.arange(self.ny) + 0.5) * self.cell

    def idx(self, x, y):
        i = np.floor((np.asarray(y, float) - self.y0) / self.cell).astype(int)
        j = np.floor((np.asarray(x, float) - self.x0) / self.cell).astype(int)
        return np.clip(i, 0, self.ny - 1), np.clip(j, 0, self.nx - 1)


def _new_field(pts_xy, cell, margin):
    mn = pts_xy.min(axis=0) - margin
    mx = pts_xy.max(axis=0) + margin
    nx = int(np.ceil((mx[0] - mn[0]) / cell)) + 1
    ny = int(np.ceil((mx[1] - mn[1]) / cell)) + 1
    return Field(mn[0], mn[1], nx, ny, cell)


def field_from_points(pts, cell, margin=0.0, fill_nb=7, fill_passes=3):
    """Rasterise a surface point cloud to a max-z field, then close pinholes.

    A surface cloud leaves occasional empty cells inside footprints (Poisson
    gaps in the sampling).  Only cells with at least `fill_nb` of 8 occupied
    neighbours are filled, which closes pinholes without filling courtyards
    or bridging narrow streets.  fill_nb = 7 was calibrated against the
    point-based roof area (agreement 1.0 %).
    """
    from scipy import ndimage

    f = _new_field(pts[:, :2], cell, margin)
    i, j = f.idx(pts[:, 0], pts[:, 1])
    np.maximum.at(f.H, (i, j), pts[:, 2])

    k8 = np.array([[1, 1, 1], [1, 0, 1], [1, 1, 1]], float)
    filled_total = 0
    for _ in range(fill_passes):
        occ = f.H >= 0
        nb = ndimage.convolve(occ.astype(float), k8, mode="constant")
        miss = (~occ) & (nb >= fill_nb)
        if not miss.any():
            break
        mx = ndimage.maximum_filter(np.where(occ, f.H, -1.0), size=3)
        take = miss & (mx >= 0)
        f.H[take] = mx[take]
        filled_total += int(take.sum())
    f.H[f.H < 0] = 0.0
    return f, filled_total


def field_from_tris(tris, cell, margin=0.0):
    f = _new_field(tris.reshape(-1, 3)[:, :2], cell, margin)
    xc, yc = f.xc(), f.yc()
    a = tris[:, 0]
    e1, e2 = tris[:, 1] - a, tris[:, 2] - a
    d = e1[:, 0] * e2[:, 1] - e2[:, 0] * e1[:, 1]
    ok = np.abs(d) > 1e-12
    max_cells = 400_000
    for k in np.nonzero(ok)[0]:
        xs, ys = tris[k, :, 0], tris[k, :, 1]
        j0 = max(0, int(np.floor((xs.min() - f.x0) / cell)))
        j1 = min(f.nx - 1, int(np.ceil((xs.max() - f.x0) / cell)))
        i0 = max(0, int(np.floor((ys.min() - f.y0) / cell)))
        i1 = min(f.ny - 1, int(np.ceil((ys.max() - f.y0) / cell)))
        if j1 < j0 or i1 < i0:
            continue
        step = max(1, max_cells // max(1, j1 - j0 + 1))
        for ia in range(i0, i1 + 1, step):
            ib = min(i1, ia + step - 1)
            X, Y = np.meshgrid(xc[j0:j1 + 1], yc[ia:ib + 1])
            px, py = X - a[k, 0], Y - a[k, 1]
            u = (px * e2[k, 1] - e2[k, 0] * py) / d[k]
            v = (e1[k, 0] * py - px * e1[k, 1]) / d[k]
            ins = (u >= -1e-9) & (v >= -1e-9) & (u + v <= 1 + 1e-9)
            if not ins.any():
                continue
            z = a[k, 2] + u * e1[k, 2] + v * e2[k, 2]
            blk = f.H[ia:ib + 1, j0:j1 + 1]
            np.maximum(blk, np.where(ins, z, blk), out=blk)
    if (~ok).any():
        v = tris[~ok].reshape(-1, 3)
        i, j = f.idx(v[:, 0], v[:, 1])
        np.maximum.at(f.H, (i, j), v[:, 2])
    f.H[f.H < 0] = 0.0
    return f, 0


def field_from_footprints(polys, cell, margin=0.0):
    from matplotlib.path import Path
    allp = np.vstack([p for p, _ in polys])
    f = _new_field(allp, cell, margin)
    xc, yc = f.xc(), f.yc()
    for ring, h in polys:
        j0 = max(0, int(np.floor((ring[:, 0].min() - f.x0) / cell)))
        j1 = min(f.nx - 1, int(np.ceil((ring[:, 0].max() - f.x0) / cell)))
        i0 = max(0, int(np.floor((ring[:, 1].min() - f.y0) / cell)))
        i1 = min(f.ny - 1, int(np.ceil((ring[:, 1].max() - f.y0) / cell)))
        if j1 < j0 or i1 < i0:
            continue
        X, Y = np.meshgrid(xc[j0:j1 + 1], yc[i0:i1 + 1])
        m = Path(ring).contains_points(
            np.column_stack([X.ravel(), Y.ravel()])).reshape(X.shape)
        blk = f.H[i0:i1 + 1, j0:j1 + 1]
        np.maximum(blk, np.where(m, h, blk), out=blk)
    f.H[f.H < 0] = 0.0
    return f, 0


# ---------------------------------------------------------------------------
# morphology
# ---------------------------------------------------------------------------


def frontal_step_map(H, cell, wind_deg):
    """Per-cell windward frontal area from positive height steps.

    Walking downwind, every upward step in H is a wall face presented to the
    wind, counted once per building.  This reproduces the conventional
    frontal-area definition without needing surface normals, and unlike a
    silhouette projection it does not let an upstream building mask a
    downstream one of equal height.
    """
    th = np.radians(wind_deg % 360.0)
    wx, wy = np.cos(th), np.sin(th)
    A = np.zeros_like(H)
    if abs(wx) > 1e-9:
        dx = np.zeros_like(H)
        dx[:, 1:] = H[:, 1:] - H[:, :-1]
        A += np.clip(dx * np.sign(wx), 0, None) * abs(wx) * cell
    if abs(wy) > 1e-9:
        dy = np.zeros_like(H)
        dy[1:, :] = H[1:, :] - H[:-1, :]
        A += np.clip(dy * np.sign(wy), 0, None) * abs(wy) * cell
    return A


def wall_step_map(H, cell):
    """Per-cell total vertical surface area from |dH| across both cell faces."""
    dx = np.zeros_like(H)
    dx[:, 1:] = np.abs(H[:, 1:] - H[:, :-1])
    dy = np.zeros_like(H)
    dy[1:, :] = np.abs(H[1:, :] - H[:-1, :])
    return (dx + dy) * cell


def disc_kernel(radius_cells):
    r = int(np.ceil(radius_cells))
    yy, xx = np.mgrid[-r:r + 1, -r:r + 1]
    return (xx ** 2 + yy ** 2) <= radius_cells ** 2


def moving_mean(a, radius_cells):
    """Circular moving average.  FFT-based: ndimage.convolve allocates
    temporaries scaling with the footprint and exhausts memory at R >= 150 m."""
    from scipy import signal
    k = disc_kernel(radius_cells).astype(float)
    k /= k.sum()
    return signal.fftconvolve(a.astype(float), k, mode="same")


def window_metrics(f, mask, dist, Afr, Awl, x, y, R, h_thresh):
    rc = R / f.cell
    i, j = (int(v) for v in f.idx(x, y))
    r = int(np.ceil(rc))
    i0, i1 = max(0, i - r), min(f.ny, i + r + 1)
    j0, j1 = max(0, j - r), min(f.nx, j + r + 1)
    yy, xx = np.mgrid[i0 - i:i1 - i, j0 - j:j1 - j]
    disc = (xx ** 2 + yy ** 2) <= rc ** 2
    if not disc.any():
        return None

    sl = (slice(i0, i1), slice(j0, j1))
    Hw = f.H[sl][disc]
    Mw = mask[sl][disc]
    Dw = dist[sl][disc]
    plan = disc.sum() * f.cell ** 2

    lam_p = float(Mw.mean())
    H_built = float(Hw[Mw].mean()) if Mw.any() else 0.0
    lam_f = float(Afr[sl][disc].sum() / plan)
    wall_a = float(Awl[sl][disc].sum())
    roof_a = float(Mw.sum() * f.cell ** 2)
    # W_char = C * <distance to nearest building> over open cells.
    # C = 10/3 calibrated on regular arrays of known street width (building
    # 20-40 m, street 10-20 m): C = 2 recovers 0.58-0.62 of the true width,
    # so C = 10/3 recovers it to within 3 %.  In real geometry the distance
    # transform also samples plazas and courtyards, so W_char is an openness
    # length scale rather than a measured canyon width -- see README.
    W = (10.0 / 3.0) * float(Dw[~Mw].mean()) if (~Mw).any() else 0.0

    return dict(R=R, lambda_p=lam_p, lambda_f=lam_f, H_built=H_built,
                H_all=float(Hw.mean()), W_char=W,
                aspect=(H_built / W if W > 0 else np.nan),
                openness=1.0 - lam_p, plan_area=plan,
                roof_area=roof_a, wall_area=wall_a,
                ground_area=plan - roof_a)


def builtup_box(f, mask, q=0.005):
    ii, jj = np.nonzero(mask)
    if len(ii) == 0:
        return None
    xs = f.x0 + (jj + 0.5) * f.cell
    ys = f.y0 + (ii + 0.5) * f.cell
    return (float(np.quantile(xs, q)), float(np.quantile(xs, 1 - q)),
            float(np.quantile(ys, q)), float(np.quantile(ys, 1 - q)))


# ---------------------------------------------------------------------------
# point-based surface areas
# ---------------------------------------------------------------------------


def point_surface_areas(pts, f, k=12, roof_tol=1.5, h_thresh=1.0):
    """Classify surface points and total their area.

    Area per point uses the local surface-density estimator pi*r_k^2/k, exact
    for a locally planar, locally uniform sample.  A point is roof if it lies
    within roof_tol of the maximum height of its own raster column.
    """
    from scipy.spatial import cKDTree

    tree = cKDTree(pts)
    d, _ = tree.query(pts, k=k + 1, workers=-1)
    a_pt = np.pi * d[:, k] ** 2 / k

    i, j = f.idx(pts[:, 0], pts[:, 1])
    zmax = f.H[i, j]
    roof = (pts[:, 2] > h_thresh) & ((zmax - pts[:, 2]) < roof_tol)

    return dict(
        area_per_point_median=float(np.median(a_pt)),
        roof_area=float(a_pt[roof].sum()), roof_n=int(roof.sum()),
        wall_area=float(a_pt[~roof].sum()), wall_n=int((~roof).sum()),
        total_area=float(a_pt.sum()),
    ), roof, a_pt


# ---------------------------------------------------------------------------
# release-point (small-scale) metrics
# ---------------------------------------------------------------------------


def release_point_metrics(f, mask, src, r_probe=50.0, r_tight=25.0,
                          r_max=200.0, nray=360):
    """Metrics describing the release point itself, not its neighbourhood.

    lambda_p averages over a window, so at R >= 150 m it cannot express whether
    the release point is boxed in.  Two properties are added here:

    encl_25 / encl_50  fraction of azimuthal directions in which a building is
                       met within 25 / 50 m, by ray casting on the height field
    free_sector_deg    widest building-free angular sector within r_probe; 0
                       means there is no escape direction at that radius
    court_closed       1 if the open region containing the release point is not
                       connected to the main open network at ground level, i.e.
                       an interior courtyard.  A tracer released there can only
                       leave over the roofline, which is a different dispersion
                       regime from a street-level release and is invisible to
                       lambda_p at any radius.

    The courtyard test is topological and was checked to be unchanged at 1.5,
    2.0 and 2.5 m cell size and with pinhole filling on and off.
    """
    from scipy import ndimage

    lbl, _ = ndimage.label(~mask)
    counts = np.bincount(lbl.ravel())
    counts[0] = 0
    main = int(np.argmax(counts))

    th = np.arange(nray) * 2.0 * np.pi / nray
    step = f.cell * 0.5
    t = np.arange(step, r_max + step, step)
    ct, st = np.cos(th)[:, None], np.sin(th)[:, None]

    out = {}
    for s in src:
        px = s["x"] + ct * t[None, :]
        py = s["y"] + st * t[None, :]
        jj = np.floor((px - f.x0) / f.cell).astype(int)
        ii = np.floor((py - f.y0) / f.cell).astype(int)
        ok = (ii >= 0) & (ii < f.ny) & (jj >= 0) & (jj < f.nx)
        hit = np.zeros(ok.shape, bool)
        hit[ok] = mask[ii[ok], jj[ok]]
        first = np.argmax(hit, axis=1)
        any_ = hit.any(axis=1)
        dist = np.where(any_, t[first], r_max)

        free = (dist >= r_probe).astype(np.int8)
        run = best = 0
        for v in np.concatenate([free, free]):
            run = run + 1 if v else 0
            best = max(best, run)
        best = min(best, nray)

        i0, j0 = (int(v) for v in f.idx(s["x"], s["y"]))
        comp = int(lbl[i0, j0])
        closed = comp != main and comp != 0

        out[s["id"]] = dict(
            encl_25=float((dist < r_tight).mean()),
            encl_50=float((dist < r_probe).mean()),
            ray_dist_med=float(np.median(dist)),
            free_sector_deg=float(best * 360.0 / nray),
            court_closed=int(closed),
            court_area_m2=(float(counts[comp] * f.cell ** 2) if closed else 0.0),
        )
    return out


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------


def make_figures(f, mask, lam, sel, src, per, radii, R_main, outdir, wind,
                 outline_lw=0.9):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    paths = []

    shown = np.where(sel, lam, np.nan)
    vmax = max(0.1, float(np.ceil(np.nanpercentile(shown, 99.5) * 20) / 20))
    fig, ax = plt.subplots(figsize=(9.5, 7))
    im = ax.imshow(shown, origin="lower", extent=f.extent, cmap="viridis",
                   vmin=0, vmax=vmax)
    # building outlines: 0.9 pt at full opacity reads as distinct blocks at
    # this scale.  Heavier than ~1.2 pt closes the narrow courtyards in the
    # dense core, which misrepresents the geometry.
    ax.contour(f.xc(), f.yc(), mask.astype(float), levels=[0.5],
               colors="w", linewidths=outline_lw, alpha=1.0)
    for s in src:
        ax.plot(s["x"], s["y"], "o", mfc="none", mec="red", mew=1.7, ms=10)
        ax.annotate(str(s["id"]), (s["x"], s["y"]), color="red", fontsize=9,
                    fontweight="bold", xytext=(8, 5), textcoords="offset points")
    ax.annotate("", xy=(0.10, 0.965), xytext=(0.02, 0.965),
                xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color="k", lw=1.4))
    ax.text(0.105, 0.965, f"wind {wind:g}$^\\circ$", transform=ax.transAxes,
            fontsize=8, va="center")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_aspect("equal")
    ax.set_title(f"Plan area density $\\lambda_p$   (R = {R_main:g} m)")
    fig.colorbar(im, ax=ax, label="$\\lambda_p$", shrink=0.82)
    fig.tight_layout()
    p = os.path.join(outdir, "fig_lambda_p_map.png")
    fig.savefig(p, dpi=200)
    plt.close(fig)
    paths.append(p)

    vals = lam[sel & mask]
    hi = float(min(1.0, np.nanpercentile(vals, 99.9) * 1.08))
    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    ax.hist(vals, bins=40, range=(0, hi), color="0.78", edgecolor="0.35")
    ymax = ax.get_ylim()[1]
    ax.set_ylim(0, ymax * 1.22)
    for e in np.quantile(vals, [1 / 3, 2 / 3]):
        ax.axvline(e, color="0.25", ls=":", lw=1.1)
    order = sorted(src, key=lambda s: per[(s["id"], R_main)]["lambda_p"])
    for rank, s in enumerate(order):
        m = per[(s["id"], R_main)]
        lvl = 1.17 - 0.085 * (rank % 4)
        ax.axvline(m["lambda_p"], color="crimson", lw=1.15, ymax=lvl / 1.22)
        ax.annotate(f"S{s['id']} {m['pctile']:.0f}%",
                    (m["lambda_p"], ymax * lvl), color="crimson",
                    fontsize=7.5, ha="center", va="bottom")
    ax.set_xlabel(f"$\\lambda_p$   (R = {R_main:g} m)")
    ax.set_ylabel("built cells")
    ax.set_title("Built-up-area distribution, with source values and terciles")
    fig.tight_layout()
    p = os.path.join(outdir, "fig_lambda_p_hist.png")
    fig.savefig(p, dpi=200)
    plt.close(fig)
    paths.append(p)

    fig, axs = plt.subplots(1, 2, figsize=(11, 4.4))
    for s in src:
        axs[0].plot(radii, [per[(s["id"], R)]["lambda_p"] for R in radii],
                    "o-", ms=4, label=f"S{s['id']}")
        axs[1].plot(radii, [per[(s["id"], R)]["aspect"] for R in radii],
                    "o-", ms=4)
    axs[0].set_ylabel("$\\lambda_p$")
    axs[1].set_ylabel("$H/W$")
    for a in axs:
        a.set_xlabel("window radius R [m]")
    axs[0].legend(fontsize=7, ncol=2)
    fig.suptitle("Sensitivity of source-local morphology to window radius")
    fig.tight_layout()
    p = os.path.join(outdir, "fig_R_sensitivity.png")
    fig.savefig(p, dpi=200)
    plt.close(fig)
    paths.append(p)

    # morphology cross-plot: shows whether the metrics are collinear, i.e.
    # whether they can act as independent predictors in Sec. 4.2
    labs = {"lambda_f": "$\\lambda_f$", "H_built": "$H_{built}$ [m]",
            "aspect": "$H/W$"}
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.9))
    xs = [per[(s["id"], R_main)]["lambda_p"] for s in src]
    for ax, (ky, lb) in zip(axs, labs.items()):
        ys = [per[(s["id"], R_main)][ky] for s in src]
        ax.plot(xs, ys, "o", color="crimson")
        for s, xx, yy in zip(src, xs, ys):
            ax.annotate(s["id"], (xx, yy), fontsize=7,
                        xytext=(4, 3), textcoords="offset points")
        if len(xs) > 2:
            ax.set_title(f"r = {np.corrcoef(xs, ys)[0, 1]:+.2f}", fontsize=9)
        ax.set_xlabel("$\\lambda_p$")
        ax.set_ylabel(lb)
    fig.suptitle(f"Source morphology metrics against $\\lambda_p$ "
                 f"(R = {R_main:g} m)")
    fig.tight_layout()
    p = os.path.join(outdir, "fig_metric_collinearity.png")
    fig.savefig(p, dpi=200)
    plt.close(fig)
    paths.append(p)

    return paths


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------


def read_sources(path):
    with open(path, newline="") as fh:
        rows = [r for r in csv.reader(fh) if r and r[0].strip()]
    hdr = rows[0]
    try:
        float(hdr[1])
        start, cols = 0, {"id": 0, "x": 1, "y": 2, "label": None}
    except (ValueError, IndexError):
        start = 1
        low = [c.strip().lower() for c in hdr]
        cols = {k: (low.index(k) if k in low else i)
                for i, k in enumerate(("id", "x", "y"))}
        cols["label"] = low.index("label") if "label" in low else None
    src = []
    for r in rows[start:]:
        if len(r) <= max(cols["x"], cols["y"]):
            continue
        li = cols.get("label")
        lab = r[li].strip() if (li is not None and len(r) > li) else ""
        src.append(dict(id=r[cols["id"]].strip(), x=float(r[cols["x"]]),
                        y=float(r[cols["y"]]), label=lab))
    if not src:
        raise ValueError(f"{path}: no sources parsed")
    return src


def run(geom, sources, cell, radii, wind, outdir, h_thresh=1.0, margin=0.0,
        box=None, suggest=0, outline_lw=0.9):
    os.makedirs(outdir, exist_ok=True)
    ext = os.path.splitext(geom)[1].lower()

    pts = None
    if ext == ".vtk":
        pts = load_vtk_points(geom)
        print(f"  VTK points: {len(pts):,}")
        f, nfill = field_from_points(pts, cell, margin)
    elif ext == ".stl":
        tris = load_stl(geom)
        print(f"  STL facets: {len(tris):,}")
        f, nfill = field_from_tris(tris, cell, margin)
    elif ext in (".csv", ".xyz", ".txt", ".dat"):
        pts = load_points_csv(geom)
        print(f"  points: {len(pts):,}")
        f, nfill = field_from_points(pts, cell, margin)
    elif ext in (".geojson", ".json"):
        f, nfill = field_from_footprints(load_footprints(geom), cell, margin)
    else:
        raise ValueError(f"unsupported geometry extension: {ext}")

    x0, x1, y0, y1 = f.extent
    print(f"  bbox      : x [{x0:.1f}, {x1:.1f}]  y [{y0:.1f}, {y1:.1f}]  "
          f"z max {f.H.max():.1f}")
    print(f"  raster    : {f.nx} x {f.ny} cells at {f.cell:g} m"
          + (f"   ({nfill:,} pinhole cells closed)" if nfill else ""))

    mask = f.H > h_thresh
    from scipy import ndimage
    dist = ndimage.distance_transform_edt(~mask, sampling=f.cell)
    Afr = frontal_step_map(f.H, f.cell, wind)
    Awl = wall_step_map(f.H, f.cell)

    if box is None:
        box = builtup_box(f, mask)
    print(f"  built-up  : x [{box[0]:.0f}, {box[1]:.0f}]  y [{box[2]:.0f}, "
          f"{box[3]:.0f}]   built fraction {mask.mean():.3f}")

    plan_area = float(mask.sum() * f.cell ** 2)
    field_areas = dict(roof_area=plan_area, wall_area=float(Awl.sum()))
    pt_areas = None
    if pts is not None:
        pt_areas, _, _ = point_surface_areas(pts, f, h_thresh=h_thresh)
    analysis_area = (box[1] - box[0]) * (box[3] - box[2])

    print("\n  surface area [m2]        point-based    field-based    diff")
    print("  " + "-" * 60)
    for key, lab in (("roof_area", "roof"), ("wall_area", "wall")):
        fv = field_areas[key]
        if pt_areas:
            pv = pt_areas[key]
            print(f"  {lab:<22} {pv:>13,.0f} {fv:>14,.0f} "
                  f"{100 * (pv / fv - 1):>+7.1f}%")
        else:
            print(f"  {lab:<22} {'-':>13} {fv:>14,.0f}")
    print(f"  {'ground (derived)':<22} {'-':>13} "
          f"{analysis_area - plan_area:>14,.0f}")
    print(f"  {'built-up analysis area':<22} {'-':>13} {analysis_area:>14,.0f}")
    print("  note: the geometry file contains no ground sheet; ground area is")
    print("        the built-up analysis area minus building plan area.")

    src = read_sources(sources) if isinstance(sources, str) else sources
    R_main = radii[len(radii) // 2]
    rp = release_point_metrics(f, mask, src)
    per, maps = {}, {}
    for R in radii:
        lam = np.clip(moving_mean(mask, R / f.cell), 0, 1)
        i0, j0 = f.idx(box[0], box[2])
        i1, j1 = f.idx(box[1], box[3])
        sel = np.zeros_like(mask, bool)
        sel[int(i0):int(i1) + 1, int(j0):int(j1) + 1] = True
        maps[R] = (lam, sel)
        pool = lam[sel & mask]
        for s in src:
            m = window_metrics(f, mask, dist, Afr, Awl, s["x"], s["y"], R,
                               h_thresh)
            if m is None:
                continue
            m["pctile"] = 100.0 * float((pool < m["lambda_p"]).mean())
            m["id"], m["label"] = s["id"], s["label"]
            m.update(rp[s["id"]])
            per[(s["id"], R)] = m

    cols = ["id", "label", "R", "lambda_p", "pctile", "lambda_f", "H_built",
            "H_all", "W_char", "aspect", "openness", "plan_area", "roof_area",
            "wall_area", "ground_area",
            "encl_25", "encl_50", "ray_dist_med", "free_sector_deg",
            "court_closed", "court_area_m2"]
    p_csv = os.path.join(outdir, "morph_sources.csv")
    with open(p_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for R in radii:
            for s in src:
                m = per.get((s["id"], R))
                if not m:
                    continue
                row = {c: m.get(c, "") for c in cols}
                for c in cols:
                    if isinstance(row[c], float):
                        row[c] = (f"{row[c]:.1f}" if c.endswith("_area")
                                  else f"{row[c]:.4f}")
                w.writerow(row)

    p_srf = os.path.join(outdir, "surface_areas.csv")
    with open(p_srf, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["class", "point_based_m2", "field_based_m2", "n_points"])
        for key, lab in (("roof_area", "roof"), ("wall_area", "wall")):
            w.writerow([lab,
                        f"{pt_areas[key]:.1f}" if pt_areas else "",
                        f"{field_areas[key]:.1f}",
                        pt_areas[key.replace("_area", "_n")] if pt_areas else ""])
        w.writerow(["ground_derived", "", f"{analysis_area - plan_area:.1f}", ""])
        w.writerow(["builtup_analysis_area", "", f"{analysis_area:.1f}", ""])

    lam, sel = maps[R_main]
    figs = make_figures(f, mask, lam, sel, src, per, radii, R_main, outdir,
                        wind, outline_lw)

    np.savez_compressed(os.path.join(outdir, "heightfield.npz"),
                        H=f.H, frontal=Afr, wall=Awl, x0=f.x0, y0=f.y0,
                        cell=f.cell, builtup_box=np.asarray(box),
                        wind_deg=wind)

    print(f"\n  source-local morphology at R = {R_main:g} m")
    print("  " + "-" * 74)
    print(f"  {'src':>4} {'lam_p':>7} {'pct':>5} {'lam_f':>7} {'H_blt':>7} "
          f"{'W':>7} {'H/W':>6}  label")
    for s in src:
        m = per.get((s["id"], R_main))
        if m:
            print(f"  {s['id']:>4} {m['lambda_p']:7.3f} {m['pctile']:4.0f}% "
                  f"{m['lambda_f']:7.3f} {m['H_built']:7.1f} {m['W_char']:7.1f} "
                  f"{m['aspect']:6.2f}  {m['label']}")

    print(f"\n  release-point environment (independent of R)")
    print("  " + "-" * 74)
    print(f"  {'src':>4} {'encl<25m':>9} {'encl<50m':>9} {'free sect':>10} "
          f"{'med dist':>9}  character")
    for s in src:
        q = rp[s["id"]]
        if q["court_closed"]:
            ch = f"closed courtyard, {q['court_area_m2']:,.0f} m2"
        elif q["encl_50"] > 0.55:
            ch = "canyon / confined"
        elif q["encl_50"] > 0.30:
            ch = "semi-open"
        else:
            ch = "open"
        print(f"  {s['id']:>4} {q['encl_25']:9.2f} {q['encl_50']:9.2f} "
              f"{q['free_sector_deg']:8.0f} deg {q['ray_dist_med']:7.1f} m  {ch}")
    nclosed = sum(q["court_closed"] for q in rp.values())
    if nclosed:
        print(f"  {nclosed} source(s) sit in interior courtyards not connected to")
        print("  the open network at ground level.  Their plumes can only leave")
        print("  over the roofline, which lambda_p does not express at any R.")

    lp = np.array([per[(s["id"], R_main)]["lambda_p"] for s in src
                   if (s["id"], R_main) in per])
    pc = np.array([per[(s["id"], R_main)]["pctile"] for s in src
                   if (s["id"], R_main) in per])
    pool = lam[sel & mask]
    print("  " + "-" * 74)
    print(f"  lambda_p span   : {lp.min():.3f} - {lp.max():.3f} "
          f"(range {lp.max() - lp.min():.3f})")
    print(f"  percentile span : {pc.min():.0f}% - {pc.max():.0f}%")
    print(f"  built-up 5/50/95: {np.quantile(pool, [0.05, 0.5, 0.95]).round(3)}")

    edges = np.quantile(pool, [0, 1 / 3, 2 / 3, 1.0])
    occ = np.histogram(lp, bins=edges)[0]
    occ[-1] += int((lp >= edges[-1]).sum())
    names = ("low", "mid", "high")
    print("  tercile occupancy: "
          + " / ".join(f"{n}={c}" for n, c in zip(names, occ))
          + f"   (edges {edges[1]:.3f}, {edges[2]:.3f})")

    # gap and leverage: tercile counts can look acceptable while the points sit
    # in two clumps, in which case the regression slope rests on the outliers
    ids_sorted = [s["id"] for s in
                  sorted(src, key=lambda s: per[(s["id"], R_main)]["lambda_p"])]
    lps = np.sort(lp)
    gaps = np.diff(lps)
    gi = int(gaps.argmax())
    gap_frac = gaps[gi] / (lp.max() - lp.min())
    X = np.column_stack([np.ones_like(lp), lp])
    hat = np.diag(X @ np.linalg.pinv(X.T @ X) @ X.T)
    lev_pairs = sorted(zip(lp, [s["id"] for s in src], hat))
    print(f"  largest gap      : {gaps[gi]:.3f} between S{ids_sorted[gi]} "
          f"({lps[gi]:.3f}) and S{ids_sorted[gi + 1]} ({lps[gi + 1]:.3f})"
          f"  = {100 * gap_frac:.0f}% of the range")
    flagged = [i for _, i, h in lev_pairs if h > 2 * hat.mean()]
    print(f"  max leverage     : {hat.max():.2f} (2x mean = {2 * hat.mean():.2f})"
          + (f"   high: {', '.join('S' + i for i in flagged)}" if flagged
             else "   none flagged"))

    need = 3
    short = [n for n, c in zip(names, occ) if c < need]
    print()
    if lp.max() - lp.min() < 0.20:
        print("  VERDICT: lambda_p range too narrow. The morphology-deposition")
        print("           scatter cannot support a correlation as configured.")
    elif gap_frac > 0.30:
        print(f"  VERDICT: the range is wide but {100 * gap_frac:.0f}% of it is a "
              "single empty gap.")
        print("           The sources form two clumps, so a fitted slope would")
        print("           rest on the extreme points rather than on a trend.")
        print("           Infill sources inside the gap band are needed; the")
        print("           saved flow field means only the LDM pass reruns.")
    elif short:
        print(f"  VERDICT: range adequate, but the {', '.join(short)} tercile(s)")
        print(f"           hold fewer than {need} sources. Consider adding")
        print("           sources there; the saved flow field means only the")
        print("           LDM pass needs rerunning.")
    else:
        print("  VERDICT: range and tercile occupancy both adequate. The source")
        print("           set can carry the correlation. Report the percentile")
        print("           span in Sec. 4.1 as the stratification evidence.")

    out = [p_csv, p_srf] + figs
    if suggest:
        p_sug = suggest_infill(f, mask, lam, sel, src, per, R_main, suggest,
                               outdir, dist=dist)
        if p_sug:
            out.append(p_sug)
    return out


def suggest_infill(f, mask, lam, sel, src, per, R_main, n_new, outdir,
                   min_sep=120.0, dist=None):
    """Propose open-space locations that fill the widest gaps in lambda_p.

    Candidates must be open cells (a release point is at street level, not
    inside a building), inside the built-up box, at least R_main from its edge
    so the averaging window is fully supported, min_sep from every existing
    and previously proposed source, and -- important -- no closer to a wall
    than the existing sources already are.  A release point pressed against a
    facade has its footprint set by that one wall rather than by the
    neighbourhood morphology, which is exactly the confounding the correlation
    is trying to avoid.
    """
    from scipy import ndimage

    if dist is None:
        dist = ndimage.distance_transform_edt(~mask, sampling=f.cell)

    lp = np.array([per[(s["id"], R_main)]["lambda_p"] for s in src])
    lps = np.sort(lp)
    gaps = np.diff(lps)
    targets = []
    for gi in np.argsort(gaps)[::-1]:
        if len(targets) >= n_new:
            break
        a, b = lps[gi], lps[gi + 1]
        m = min(n_new - len(targets), max(1, int(round(gaps[gi] / 0.06))))
        targets += list(np.linspace(a, b, m + 2)[1:-1])
    targets = sorted(targets[:n_new])
    if not targets:
        return None

    # wall clearance floor: the tightest siting already accepted
    src_clear = []
    for s in src:
        i, j = (int(v) for v in f.idx(s["x"], s["y"]))
        src_clear.append(dist[i, j])
    clear_min = float(np.min(src_clear)) if src_clear else 5.0

    inside = ndimage.distance_transform_edt(sel, sampling=f.cell) >= R_main
    Y, X = np.mgrid[0:f.ny, 0:f.nx]
    xm = f.x0 + (X + 0.5) * f.cell
    ym = f.y0 + (Y + 0.5) * f.cell
    base = (~mask) & sel & inside & (dist >= clear_min)
    ci, cj = np.nonzero(base)
    cx, cy, cl, cd = xm[ci, cj], ym[ci, cj], lam[ci, cj], dist[ci, cj]
    print(f"\n  infill candidates: {base.sum():,} open cells with wall clearance "
          f">= {clear_min:.1f} m (existing minimum)")

    px = [s["x"] for s in src]
    py = [s["y"] for s in src]
    picks = []
    for t in targets:
        near = np.abs(cl - t) < 0.03
        if not near.any():
            near = np.abs(cl - t) < 0.06
        if not near.any():
            continue
        idxs = np.nonzero(near)[0]
        dsep = np.min(np.hypot(cx[idxs, None] - np.array(px)[None, :],
                               cy[idxs, None] - np.array(py)[None, :]), axis=1)
        ok = dsep >= min_sep
        if not ok.any():
            ok = dsep >= dsep.max() * 0.9
        sel_i = idxs[ok]
        best = sel_i[int(np.argmax(dsep[ok]))]
        picks.append((cx[best], cy[best], cl[best],
                      float(dsep[ok].max()), cd[best]))
        px.append(cx[best])
        py.append(cy[best])

    if not picks:
        return None

    start = max(int(s["id"]) for s in src if s["id"].isdigit()) + 1
    p = os.path.join(outdir, "suggested_sources.csv")
    with open(p, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "x", "y", "z", "lambda_p", "min_separation_m",
                    "wall_clearance_m"])
        for k, (x, y, l, ds, wc) in enumerate(picks):
            w.writerow([start + k, f"{x:.1f}", f"{y:.1f}", "5.0",
                        f"{l:.3f}", f"{ds:.0f}", f"{wc:.1f}"])

    merged = np.sort(np.concatenate([lp, [q[2] for q in picks]]))
    X0 = np.column_stack([np.ones_like(lp), lp])
    h0 = np.diag(X0 @ np.linalg.pinv(X0.T @ X0) @ X0.T)
    X1 = np.column_stack([np.ones_like(merged), merged])
    h1 = np.diag(X1 @ np.linalg.pinv(X1.T @ X1) @ X1.T)
    print(f"  suggested infill sources (R = {R_main:g} m)")
    print("  " + "-" * 66)
    print(f"  {'id':>4} {'x [m]':>9} {'y [m]':>9} {'lambda_p':>9} "
          f"{'min sep':>9} {'wall d':>8}")
    for k, (x, y, l, ds, wc) in enumerate(picks):
        print(f"  {start + k:>4} {x:9.1f} {y:9.1f} {l:9.3f} {ds:7.0f} m "
              f"{wc:6.1f} m")
    print("  " + "-" * 66)
    print(f"  largest gap  {gaps.max():.3f} -> {np.diff(merged).max():.3f}")
    print(f"  max leverage {h0.max():.2f} -> {h1.max():.2f}")
    return p


# ---------------------------------------------------------------------------


def selftest(outdir):
    rng = np.random.default_rng(0)
    os.makedirs(outdir, exist_ok=True)
    pts = []
    for (xa, ya, xb, yb, pitch, fill) in [
            (60, 60, 360, 640, 45, 0.80), (380, 60, 620, 640, 70, 0.50),
            (640, 60, 860, 640, 110, 0.28)]:
        for xx in np.arange(xa, xb, pitch):
            for yy in np.arange(ya, yb, pitch):
                if rng.random() < 0.12:
                    continue
                w = pitch * fill * (0.8 + 0.4 * rng.random())
                h = 10 + 40 * rng.random()
                nr = max(4, int(w * w / 1.1))
                pts.append(np.column_stack([
                    xx + w * rng.random(nr), yy + w * rng.random(nr),
                    np.full(nr, h)]))
                per_ = 4 * w
                nw = max(8, int(per_ * h / 1.1))
                t = rng.random(nw) * per_
                side = (t // w).astype(int) % 4
                u = t % w
                sx = np.where(side == 0, xx + u,
                              np.where(side == 1, xx + w,
                                       np.where(side == 2, xx + w - u, xx)))
                sy = np.where(side == 0, yy,
                              np.where(side == 1, yy + u,
                                       np.where(side == 2, yy + w, yy + w - u)))
                pts.append(np.column_stack([sx, sy, rng.random(nw) * h]))
    pts = np.vstack(pts)
    gp = os.path.join(outdir, "_synth.csv")
    np.savetxt(gp, pts, delimiter=",", header="x,y,z", comments="")
    sp = os.path.join(outdir, "_synth_sources.csv")
    with open(sp, "w") as fh:
        fh.write("id,x,y,label\n")
        for i, (x, y, l) in enumerate([
                (180, 350, "dense"), (500, 350, "medium"), (750, 350, "sparse"),
                (370, 120, "transition"), (280, 600, "dense edge"),
                (640, 600, "sparse edge"), (450, 200, "interior")], 1):
            fh.write(f"{i},{x},{y},{l}\n")
    print(f"[selftest] {len(pts):,} synthetic surface points")
    return run(gp, sp, 2.5, [50, 100, 150], 0.0, outdir)


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--geom")
    ap.add_argument("--sources")
    ap.add_argument("--cell", type=float, default=2.5)
    ap.add_argument("--radii", default="75,150,250")
    ap.add_argument("--wind", type=float, default=0.0,
                    help="inflow direction [deg CCW from +x]")
    ap.add_argument("--h-thresh", type=float, default=1.0)
    ap.add_argument("--margin", type=float, default=0.0,
                    help="pad the raster beyond the geometry bbox [m]")
    ap.add_argument("--builtup-box", default=None,
                    help="xmin,xmax,ymin,ymax override")
    ap.add_argument("--suggest", type=int, default=0,
                    help="propose N infill source locations filling lambda_p gaps")
    ap.add_argument("--outline-lw", type=float, default=0.9,
                    help="building outline width on the map [pt]")
    ap.add_argument("--out", default="morph_out")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args(argv)

    if a.selftest:
        selftest(a.out)
        return 0
    if not (a.geom and a.sources):
        ap.error("--geom and --sources required (or --selftest)")
    run(a.geom, a.sources, a.cell,
        [float(v) for v in a.radii.split(",")], a.wind, a.out,
        a.h_thresh, a.margin,
        [float(v) for v in a.builtup_box.split(",")] if a.builtup_box else None,
        a.suggest, a.outline_lw)
    return 0


if __name__ == "__main__":
    sys.exit(main())
