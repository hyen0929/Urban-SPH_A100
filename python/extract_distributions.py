#!/usr/bin/env python3
"""
Extract x, y, z distributions at 5s intervals for LDM, SPH, and Analytical.
Output: CSV with columns [direction, time, bin_center, pdf_analytical, pdf_LDM, pdf_SPH]
"""
import numpy as np
import os, re, glob

# ============================================================
# Parameters
# ============================================================
U  = 0.3
D  = 0.01
x0, y0, z0 = 5.0, 0.0, 3.0
dt_base = 0.01
dx_sph = 0.1

vtk_dir = "Transport_advection"
out_csv = os.path.join(vtk_dir, "distributions_5s_xyz.csv")

# Initial variance from compare CSV
compare = np.genfromtxt(os.path.join(vtk_dir, "compare_analytical_LDM_SPH.csv"),
                        delimiter=",", names=True)
var0 = compare["var_analytical"][0]

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

# ============================================================
# Extract distributions
# ============================================================
times = list(range(0, 55, 5))  # 0, 5, 10, ..., 50
steps = [int(t / dt_base) for t in times]
directions = ["x", "y", "z"]
dir_idx = {"x": 0, "y": 1, "z": 2}
means = {"x": lambda t: x0 + U*t, "y": lambda t: y0, "z": lambda t: z0}

n_bins = 50

rows = []

for t_val, step in zip(times, steps):
    print(f"Processing t={t_val}s (step={step})...")

    ldm_file = os.path.join(vtk_dir, f"LDM_{step}stp.vtk")
    fluid_file = os.path.join(vtk_dir, f"fluid_{step}stp.vtk")

    ldm_pts = parse_vtk_LDM(ldm_file) if os.path.exists(ldm_file) else None
    sph_pts, sph_concn = parse_vtk_fluid(fluid_file) if os.path.exists(fluid_file) else (None, None)

    sig2 = var0 + 2.0 * D * t_val

    for d in directions:
        di = dir_idx[d]
        mu = means[d](t_val)
        sig = np.sqrt(sig2)

        # Bin range: mu +/- 4*sigma (at least 1m wide)
        half = max(4 * sig, 0.5)
        lo = mu - half
        hi = mu + half

        bin_edges = np.linspace(lo, hi, n_bins + 1)
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        bw = bin_edges[1] - bin_edges[0]

        # Analytical PDF
        pdf_ana = (1.0 / np.sqrt(2 * np.pi * sig2)) * np.exp(-(bin_centers - mu)**2 / (2 * sig2))

        # LDM PDF
        pdf_ldm = np.full(n_bins, np.nan)
        if ldm_pts is not None:
            counts, _ = np.histogram(ldm_pts[:, di], bins=bin_edges)
            pdf_ldm = counts / (len(ldm_pts) * bw)

        # SPH PDF
        pdf_sph = np.full(n_bins, np.nan)
        if sph_pts is not None and sph_concn is not None:
            coords = sph_pts[:, di]
            coords_round = np.round(coords / dx_sph) * dx_sph
            unique_vals = np.unique(coords_round)

            xv, cv = [], []
            for uv in unique_vals:
                if uv < lo or uv > hi:
                    continue
                mask = (coords_round == uv)
                xv.append(uv)
                cv.append(sph_concn[mask].sum())
            xv = np.array(xv)
            cv = np.array(cv)

            if len(xv) > 1:
                area = np.trapezoid(cv, xv)
                if area > 1e-10:
                    # Interpolate onto bin_centers
                    sph_interp = np.interp(bin_centers, xv, cv / area, left=0, right=0)
                    pdf_sph = sph_interp

        for i in range(n_bins):
            rows.append((d, t_val, bin_centers[i], pdf_ana[i], pdf_ldm[i], pdf_sph[i]))

# Write CSV
with open(out_csv, "w") as f:
    f.write("direction,time,bin_center,pdf_analytical,pdf_LDM,pdf_SPH\n")
    for r in rows:
        f.write(f"{r[0]},{r[1]:.1f},{r[2]:.6f},{r[3]:.8e},{r[4]:.8e},{r[5]:.8e}\n")

print(f"\nSaved: {out_csv}")
print(f"Total rows: {len(rows)}")
