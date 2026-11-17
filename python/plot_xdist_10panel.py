#!/usr/bin/env python3
"""
x-distribution comparison at 5s intervals (t=5..50s): 2x5 panel figure
Analytical (black line) vs SPH (blue symbols) vs LDM (red symbols)
"""
import numpy as np
import matplotlib.pyplot as plt
import os, re

# ============================================================
# Parameters
# ============================================================
U  = 0.3
D  = 0.01
x0 = 5.0
dt_base = 0.01
dx_sph = 0.1

vtk_dir = "Transport_advection"
csv_path = os.path.join(vtk_dir, "compare_analytical_LDM_SPH.csv")

# Initial variance from CSV
data = np.genfromtxt(csv_path, delimiter=",", names=True)
var0 = data["var_analytical"][0]

# ============================================================
# VTK parsers
# ============================================================
DTYPE_MAP = {"float": (">f4", 4), "double": (">f8", 8), "int": (">i4", 4)}

def parse_vtk_LDM(filepath):
    with open(filepath, "rb") as f:
        content = f.read()
    m = re.search(rb"POINTS\s+(\d+)\s+(\w+)\n", content)
    if not m: return None
    n = int(m.group(1))
    dt, sz = DTYPE_MAP.get(m.group(2).decode(), (">f4", 4))
    return np.frombuffer(content, dtype=dt, count=n*3, offset=m.end()).reshape(n, 3).astype(np.float64)

def parse_vtk_fluid(filepath):
    with open(filepath, "rb") as f:
        content = f.read()
    m = re.search(rb"POINTS\s+(\d+)\s+(\w+)\n", content)
    if not m: return None, None
    n = int(m.group(1))
    dt, sz = DTYPE_MAP.get(m.group(2).decode(), (">f4", 4))
    pts = np.frombuffer(content, dtype=dt, count=n*3, offset=m.end()).reshape(n, 3).astype(np.float64)
    m2 = re.search(rb"FIELD FieldData\s+(\d+)\n", content)
    if not m2: return pts, None
    offset = m2.end()
    for _ in range(int(m2.group(1))):
        nl = content.index(b"\n", offset)
        hdr = content[offset:nl].decode().split()
        fname, ncomp, ntuple, dts = hdr[0], int(hdr[1]), int(hdr[2]), hdr[3]
        offset = nl + 1
        dt2, isz = DTYPE_MAP.get(dts, (">f4", 4))
        if fname == "concn":
            return pts, np.frombuffer(content, dtype=dt2, count=ncomp*ntuple, offset=offset).astype(np.float64)
        offset += ncomp * ntuple * isz
    return pts, None

def analytical_pdf_x(x_arr, t_val):
    mu = x0 + U * t_val
    sig2 = var0 + 2.0 * D * t_val
    return (1.0 / np.sqrt(2 * np.pi * sig2)) * np.exp(-(x_arr - mu)**2 / (2 * sig2))

# ============================================================
# Load VTK data
# ============================================================
snap_times = list(range(5, 55, 5))  # 5, 10, ..., 50
snap_steps = [int(ts / dt_base) for ts in snap_times]

dist_data = {}
for ts, step in zip(snap_times, snap_steps):
    ldm_file = os.path.join(vtk_dir, f"LDM_{step}stp.vtk")
    fluid_file = os.path.join(vtk_dir, f"fluid_{step}stp.vtk")

    ldm_pts = parse_vtk_LDM(ldm_file) if os.path.exists(ldm_file) else None
    ldm_x = ldm_pts[:, 0] if ldm_pts is not None else np.array([])

    sph_pts, sph_concn = parse_vtk_fluid(fluid_file) if os.path.exists(fluid_file) else (None, None)

    dist_data[ts] = (ldm_x, sph_pts, sph_concn)
    n_ldm = len(ldm_x)
    n_sph = len(sph_pts) if sph_pts is not None else 0
    print(f"  t={ts}s (step={step}): LDM={n_ldm}, SPH={n_sph}")

# ============================================================
# Plot: 2x5
# ============================================================
fig, axes = plt.subplots(2, 5, figsize=(22, 8))
ms = 4

for idx, ts in enumerate(snap_times):
    row = idx // 5
    col = idx % 5
    ax = axes[row, col]

    ax.grid(True, linestyle=':', color='gray', alpha=0.6)
    ax.set_axisbelow(True)

    ldm_x, sph_pts, sph_concn = dist_data[ts]

    mu = x0 + U * ts
    sig = np.sqrt(var0 + 2.0 * D * ts)
    x_lo = mu - 4 * sig
    x_hi = mu + 4 * sig
    x_theory = np.linspace(x_lo, x_hi, 300)

    # Analytical
    ax.plot(x_theory, analytical_pdf_x(x_theory, ts), "k-", linewidth=1.5, label="Analytical")

    # LDM
    n_bins = max(20, int((x_hi - x_lo) / (2.0 * dx_sph)))
    bin_edges = np.linspace(x_lo, x_hi, n_bins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_width = bin_edges[1] - bin_edges[0]

    if len(ldm_x) > 0:
        ldm_counts, _ = np.histogram(ldm_x, bins=bin_edges)
        ldm_pdf = ldm_counts / (len(ldm_x) * bin_width)
        ax.plot(bin_centers, ldm_pdf, "r^", markersize=ms, linestyle="none", label="LDM")

    # SPH
    if sph_pts is not None and sph_concn is not None:
        sph_x_round = np.round(sph_pts[:, 0] / dx_sph) * dx_sph
        unique_x = np.unique(sph_x_round)

        xv, cv = [], []
        for ux in unique_x:
            if ux < x_lo or ux > x_hi:
                continue
            mask = (sph_x_round == ux)
            xv.append(ux)
            cv.append(sph_concn[mask].sum())
        xv = np.array(xv)
        cv = np.array(cv)

        if len(xv) > 1:
            area = np.trapezoid(cv, xv)
            if area > 1e-10:
                ax.plot(xv, cv / area, "bs", markersize=ms, linestyle="none", label="SPH")

    ax.set_xlabel("x (m)")
    if col == 0:
        ax.set_ylabel("PDF")
    ax.set_title(f"t = {ts} s")
    if idx == 0:
        ax.legend(fontsize=8)

plt.tight_layout()
out_path = os.path.join(vtk_dir, "xdist_10panel.png")
plt.savefig(out_path, dpi=200, bbox_inches="tight")
print(f"\nSaved: {out_path}")
plt.close()
