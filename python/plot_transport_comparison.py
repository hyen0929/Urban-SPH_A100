#!/usr/bin/env python3
"""
Scalar transport verification: Analytical vs SPH vs LDM
6-panel figure (2x3):
  1) mean_x   2) mean_y, mean_z   3) var_x
  4) x-dist t=10s  5) x-dist t=25s  6) x-dist t=50s
"""
import numpy as np
import matplotlib.pyplot as plt
import os, re, sys

# ============================================================
# Parameters
# ============================================================
U  = 0.3       # advection velocity (m/s)
D  = 0.01      # diffusion coefficient (m^2/s)
x0 = 5.0       # source x
y0 = 0.0       # source y
z0 = 3.0       # source z
dt_base = 0.01  # timestep

vtk_dir = "Transport_advection"
csv_path = os.path.join(vtk_dir, "compare_analytical_LDM_SPH.csv")

# ============================================================
# VTK parsers (from verify_compare.py)
# ============================================================
DTYPE_MAP = {
    "float":  (">f4", 4),
    "double": (">f8", 8),
    "int":    (">i4", 4),
}

def parse_vtk_LDM(filepath):
    with open(filepath, "rb") as f:
        content = f.read()
    m = re.search(rb"POINTS\s+(\d+)\s+(\w+)\n", content)
    if not m:
        return None
    n_pts = int(m.group(1))
    dt_str = m.group(2).decode("ascii")
    pts_offset = m.end()
    dt, sz = DTYPE_MAP.get(dt_str, (">f4", 4))
    pts = np.frombuffer(content, dtype=dt, count=n_pts*3, offset=pts_offset).reshape(n_pts, 3)
    return pts.astype(np.float64)

def parse_vtk_fluid(filepath):
    with open(filepath, "rb") as f:
        content = f.read()
    m = re.search(rb"POINTS\s+(\d+)\s+(\w+)\n", content)
    if not m:
        return None, None
    n_pts = int(m.group(1))
    dt_str = m.group(2).decode("ascii")
    pts_offset = m.end()
    dt, sz = DTYPE_MAP.get(dt_str, (">f4", 4))
    pts = np.frombuffer(content, dtype=dt, count=n_pts*3, offset=pts_offset).reshape(n_pts, 3).astype(np.float64)

    m_field = re.search(rb"FIELD FieldData\s+(\d+)\n", content)
    if not m_field:
        return pts, None
    n_fields = int(m_field.group(1))
    offset = m_field.end()
    for _ in range(n_fields):
        newline = content.index(b"\n", offset)
        header = content[offset:newline].decode("ascii").split()
        fname, ncomp, ntuple, dtype_s = header[0], int(header[1]), int(header[2]), header[3]
        offset = newline + 1
        dt2, itemsize = DTYPE_MAP.get(dtype_s, (">f4", 4))
        nbytes = ncomp * ntuple * itemsize
        if fname == "concn":
            concn = np.frombuffer(content, dtype=dt2, count=ncomp*ntuple, offset=offset).astype(np.float64)
            return pts, concn
        offset += nbytes
    return pts, None

# ============================================================
# 1. Load CSV for time series (panels 1-3)
# ============================================================
data = np.genfromtxt(csv_path, delimiter=",", names=True)

t       = data["time"]
ana_mx  = data["mean_x_analytical"]
ldm_mx  = data["mean_x_LDM"]
ldm_my  = data["mean_y_LDM"]
ldm_mz  = data["mean_z_LDM"]
sph_mx  = data["mean_x_SPH"]
sph_my  = data["mean_y_SPH"]
sph_mz  = data["mean_z_SPH"]
ana_var  = data["var_analytical"]
ldm_vx  = data["var_x_LDM"]
sph_vx  = data["var_x_SPH"]

# Analytical mean y, z are constants
ana_my = np.full_like(t, y0)
ana_mz = np.full_like(t, z0)

# ============================================================
# 2. Load VTK for x-distributions (panels 4-6)
# ============================================================
snap_times = [10, 25, 50]
snap_steps = [int(ts / dt_base) for ts in snap_times]  # 1000, 2500, 5000

# Initial variance from CSV (t=0)
var0 = ana_var[0]

def analytical_pdf_x(x_arr, t_val):
    """1D Gaussian PDF in x-direction"""
    mu = x0 + U * t_val
    sig2 = var0 + 2.0 * D * t_val
    return (1.0 / np.sqrt(2 * np.pi * sig2)) * np.exp(-(x_arr - mu)**2 / (2 * sig2))

dist_data = {}
for ts, step in zip(snap_times, snap_steps):
    ldm_file = os.path.join(vtk_dir, f"LDM_{step}stp.vtk")
    fluid_file = os.path.join(vtk_dir, f"fluid_{step}stp.vtk")

    # LDM
    ldm_pts = parse_vtk_LDM(ldm_file)
    ldm_x = ldm_pts[:, 0] if ldm_pts is not None else np.array([])

    # SPH
    sph_pts, sph_concn = parse_vtk_fluid(fluid_file)

    dist_data[ts] = (ldm_x, sph_pts, sph_concn)
    print(f"  t={ts}s (step={step}): LDM={len(ldm_x)} particles, SPH={len(sph_pts)} particles")

# ============================================================
# 3. Plot
# ============================================================
fig, axes = plt.subplots(2, 3, figsize=(15, 9))

for row in axes:
    for a in row:
        a.grid(True, linestyle=':', color='gray', alpha=0.6)
        a.set_axisbelow(True)

ms = 4  # marker size

# ----- Panel 1: mean x -----
ax = axes[0, 0]
ax.plot(t, ana_mx, "k-", linewidth=1.5, label="Analytical")
ax.plot(t, sph_mx, "bs", markersize=ms, linestyle="none", label="SPH")
ax.plot(t, ldm_mx, "r^", markersize=ms, linestyle="none", label="LDM")
ax.set_xlabel("Time (s)")
ax.set_ylabel("Mean x (m)")
ax.set_title("(a) Mean x")
ax.legend(fontsize=9)

# ----- Panel 2: mean y, mean z -----
ax = axes[0, 1]
ax.plot(t, ana_my, "k-", linewidth=1.5, label="Analytical y")
ax.plot(t, ana_mz, "k--", linewidth=1.5, label="Analytical z")
ax.plot(t, sph_my, "bs", markersize=ms, linestyle="none", label="SPH y")
ax.plot(t, sph_mz, "bD", markersize=ms, linestyle="none", label="SPH z")
ax.plot(t, ldm_my, "r^", markersize=ms, linestyle="none", label="LDM y")
ax.plot(t, ldm_mz, "rv", markersize=ms, linestyle="none", label="LDM z")
ax.set_xlabel("Time (s)")
ax.set_ylabel("Mean y, z (m)")
ax.set_title("(b) Mean y, z")
ax.legend(fontsize=8, ncol=2)

# ----- Panel 3: variance x -----
ax = axes[0, 2]
ax.plot(t, ana_var, "k-", linewidth=1.5, label="Analytical")
ax.plot(t, sph_vx, "bs", markersize=ms, linestyle="none", label="SPH")
ax.plot(t, ldm_vx, "r^", markersize=ms, linestyle="none", label="LDM")
ax.set_xlabel("Time (s)")
ax.set_ylabel(r"Variance $\sigma_x^2$ (m$^2$)")
ax.set_title(r"(c) Variance $\sigma_x^2$")
ax.legend(fontsize=9)

# ----- Panels 4-6: x-distributions at t=10,25,50 -----
panel_labels = ["(d)", "(e)", "(f)"]
dx_sph = 0.1  # SPH particle spacing

for idx, ts in enumerate(snap_times):
    ax = axes[1, idx]
    ldm_x, sph_pts, sph_concn = dist_data[ts]

    # x-range: center on analytical mean, +/- 4 sigma
    mu = x0 + U * ts
    sig = np.sqrt(var0 + 2.0 * D * ts)
    x_lo = mu - 4 * sig
    x_hi = mu + 4 * sig
    x_theory = np.linspace(x_lo, x_hi, 300)

    # Analytical
    pdf_theory = analytical_pdf_x(x_theory, ts)
    ax.plot(x_theory, pdf_theory, "k-", linewidth=1.5, label="Analytical")

    # Use bin width >= 2*dx_sph to avoid aliasing
    n_bins = max(20, int((x_hi - x_lo) / (2.0 * dx_sph)))
    bin_edges = np.linspace(x_lo, x_hi, n_bins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_width = bin_edges[1] - bin_edges[0]

    # LDM histogram as PDF
    ldm_counts, _ = np.histogram(ldm_x, bins=bin_edges)
    ldm_pdf = ldm_counts / (len(ldm_x) * bin_width)
    ax.plot(bin_centers, ldm_pdf, "r^", markersize=ms, linestyle="none", label="LDM")

    # SPH: group by unique x-coordinate, sum concentration per x-slice, then bin
    sph_pdf_interp = None
    if sph_concn is not None:
        sph_x_all = sph_pts[:, 0]
        total_C = sph_concn.sum()

        # Round to grid spacing to group particles on the same x-plane
        sph_x_round = np.round(sph_x_all / dx_sph) * dx_sph
        unique_x = np.unique(sph_x_round)

        # Sum concentration per x-plane
        sph_xvals = []
        sph_cvals = []
        for ux in unique_x:
            if ux < x_lo or ux > x_hi:
                continue
            mask = (sph_x_round == ux)
            c_sum = sph_concn[mask].sum()
            sph_xvals.append(ux)
            sph_cvals.append(c_sum)

        sph_xvals = np.array(sph_xvals)
        sph_cvals = np.array(sph_cvals)

        # Normalize to PDF (area = 1)
        area = np.trapezoid(sph_cvals, sph_xvals) if len(sph_xvals) > 1 else 1.0
        sph_pdf = sph_cvals / area
        ax.plot(sph_xvals, sph_pdf, "bs", markersize=ms, linestyle="none", label="SPH")

        # Interpolate SPH PDF onto bin_centers for error calculation
        sph_pdf_interp = np.interp(bin_centers, sph_xvals, sph_pdf, left=0, right=0)

    # NRMSE: normalized by analytical peak
    ana_at_bins = analytical_pdf_x(bin_centers, ts)
    peak_ana = ana_at_bins.max()

    ldm_nrmse = np.sqrt(np.mean((ldm_pdf - ana_at_bins)**2)) / peak_ana * 100
    err_txt = f"LDM: {ldm_nrmse:.1f}%"

    if sph_pdf_interp is not None:
        sph_nrmse = np.sqrt(np.mean((sph_pdf_interp - ana_at_bins)**2)) / peak_ana * 100
        err_txt += f"\nSPH: {sph_nrmse:.1f}%"

    ax.text(0.97, 0.95, f"NRMSE\n{err_txt}", transform=ax.transAxes,
            fontsize=8, verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

    ax.set_xlabel("x (m)")
    if idx == 0:
        ax.set_ylabel("PDF")
    ax.set_title(f"{panel_labels[idx]} x-distribution, t={ts}s")
    ax.legend(fontsize=9, loc='upper left')

plt.tight_layout()
out_path = os.path.join(vtk_dir, "transport_comparison_6panel.png")
plt.savefig(out_path, dpi=200, bbox_inches="tight")
print(f"\nSaved: {out_path}")
plt.close()
