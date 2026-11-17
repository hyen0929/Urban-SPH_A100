#!/usr/bin/env python3
"""
해석해 vs LDM vs SPH scalar transport 비교 CSV 생성
사용법: python3 verify_compare.py <vtk_dir> [U] [D] [x0] [y0] [z0]
"""
import numpy as np
import os, sys, re, glob

DTYPE_MAP = {
    "float":  (">f4", 4),
    "double": (">f8", 8),
    "int":    (">i4", 4),
    "unsigned_int": (">u4", 4),
    "short":  (">i2", 2),
    "char":   (">i1", 1),
}

# ============================================================
# VTK parsers
# ============================================================
def parse_vtk_LDM(filepath):
    with open(filepath, "rb") as f:
        content = f.read()
    m = re.search(rb"POINTS\s+(\d+)\s+(\w+)\n", content)
    if not m:
        return None
    n_pts = int(m.group(1))
    dtype_str = m.group(2).decode("ascii")
    pts_offset = m.end()
    dt, sz = DTYPE_MAP.get(dtype_str, (">f4", 4))
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


def get_step(fpath, pattern):
    m = re.search(pattern, fpath)
    return int(m.group(1)) if m else -1


# ============================================================
# Main
# ============================================================
def analyze(vtk_dir, U, D, x0, y0, z0, dt_base):
    # --- LDM files ---
    ldm_files = sorted(glob.glob(os.path.join(vtk_dir, "LDM_[0-9]*stp.vtk")),
                       key=lambda f: get_step(f, r"LDM_(\d+)stp"))
    # --- Fluid files ---
    fluid_files = sorted(glob.glob(os.path.join(vtk_dir, "fluid_[0-9]*stp.vtk")),
                         key=lambda f: get_step(f, r"fluid_(\d+)stp"))

    print(f"LDM files: {len(ldm_files)}, Fluid files: {len(fluid_files)}")
    print(f"U={U}, D={D}, source=({x0},{y0},{z0}), dt={dt_base}")

    # --- Build step-keyed dictionaries ---
    ldm_data = {}
    for fpath in ldm_files:
        step = get_step(fpath, r"LDM_(\d+)stp")
        pts = parse_vtk_LDM(fpath)
        if pts is None or len(pts) == 0:
            continue
        N = len(pts)
        mx, my, mz = pts[:,0].mean(), pts[:,1].mean(), pts[:,2].mean()
        vx, vy, vz = pts[:,0].var(), pts[:,1].var(), pts[:,2].var()
        ldm_data[step] = (N, mx, my, mz, vx, vy, vz)

    sph_data = {}
    for fpath in fluid_files:
        step = get_step(fpath, r"fluid_(\d+)stp")
        pts, concn = parse_vtk_fluid(fpath)
        if pts is None or concn is None:
            continue
        total_C = concn.sum()
        if total_C < 1e-10:
            continue
        w = concn / total_C
        mx = np.sum(w * pts[:,0])
        my = np.sum(w * pts[:,1])
        mz = np.sum(w * pts[:,2])
        vx = np.sum(w * (pts[:,0] - mx)**2)
        vy = np.sum(w * (pts[:,1] - my)**2)
        vz = np.sum(w * (pts[:,2] - mz)**2)
        idx_peak = np.argmax(concn)
        peak_c = concn[idx_peak]
        sph_data[step] = (mx, my, mz, vx, vy, vz, total_C, peak_c)

    # --- Merge on common steps ---
    all_steps = sorted(set(ldm_data.keys()) & set(sph_data.keys()))
    if not all_steps:
        print("No common steps found!")
        return

    print(f"Common steps: {len(all_steps)} ({all_steps[0]} ~ {all_steps[-1]})")

    # Initial variances (해석해 기준: 각각의 t=0 분산에서 출발)
    s0 = all_steps[0]
    ldm0 = ldm_data[s0]
    sph0 = sph_data[s0]
    ldm_var0 = (ldm0[4] + ldm0[5] + ldm0[6]) / 3.0
    sph_var0 = (sph0[3] + sph0[4] + sph0[5]) / 3.0

    print(f"\nInitial variance (step={s0}):")
    print(f"  LDM:  var0 = {ldm_var0:.6f} (σ₀={np.sqrt(ldm_var0):.4f})")
    print(f"  SPH:  var0 = {sph_var0:.6f} (σ₀={np.sqrt(sph_var0):.4f})")

    # --- CSV output ---
    csv_path = os.path.join(vtk_dir, "compare_analytical_LDM_SPH.csv")
    with open(csv_path, "w") as f:
        f.write("time,"
                "mean_x_analytical,"
                "mean_x_LDM,mean_y_LDM,mean_z_LDM,"
                "mean_x_SPH,mean_y_SPH,mean_z_SPH,"
                "var_analytical,var_x_LDM,var_y_LDM,var_z_LDM,"
                "var_x_SPH,var_y_SPH,var_z_SPH,"
                "N_LDM,total_C_SPH,peak_C_SPH\n")

        for step in all_steps:
            t = step * dt_base
            ldm = ldm_data[step]
            sph = sph_data[step]

            # Analytical (공통: σ₀² + 2Dt, σ₀는 LDM 기준)
            ana_mean_x = x0 + U * t
            ana_var = ldm_var0 + 2.0 * D * t

            # LDM: N, mx, my, mz, vx, vy, vz
            # SPH: mx, my, mz, vx, vy, vz, total_C, peak_c
            f.write(f"{t:.4f},"
                    f"{ana_mean_x:.6f},"
                    f"{ldm[1]:.6f},{ldm[2]:.6f},{ldm[3]:.6f},"
                    f"{sph[0]:.6f},{sph[1]:.6f},{sph[2]:.6f},"
                    f"{ana_var:.6f},{ldm[4]:.6f},{ldm[5]:.6f},{ldm[6]:.6f},"
                    f"{sph[3]:.6f},{sph[4]:.6f},{sph[5]:.6f},"
                    f"{ldm[0]},{sph[6]:.2f},{sph[7]:.6f}\n")

    print(f"\nCSV saved: {csv_path}")

    # --- Summary table ---
    print(f"\n{'step':>6} | {'time':>5} | {'<x> ana':>8} | {'<x> LDM':>8} | {'<x> SPH':>8} | {'Var ana':>9} | {'Var_x LDM':>9} | {'Var_x SPH':>9}")
    print("-" * 90)
    for step in all_steps:
        t = step * dt_base
        ldm = ldm_data[step]
        sph = sph_data[step]
        ana_mx = x0 + U * t
        ana_var = ldm_var0 + 2.0 * D * t
        print(f"{step:>6} | {t:>5.1f} | {ana_mx:>8.3f} | {ldm[1]:>8.3f} | {sph[0]:>8.3f} | {ana_var:>9.5f} | {ldm[4]:>9.5f} | {sph[3]:>9.5f}")

    # Summary
    last = all_steps[-1]
    t_last = last * dt_base
    ldm_last = ldm_data[last]
    sph_last = sph_data[last]
    ana_var_last = ldm_var0 + 2.0 * D * t_last

    print(f"\n{'='*60}")
    print(f"SUMMARY (t={t_last:.1f}s)")
    print(f"{'='*60}")
    print(f"  Analytical mean_x = {x0 + U*t_last:.3f}")
    print(f"  LDM mean_x        = {ldm_last[1]:.3f}  (err={ldm_last[1]-(x0+U*t_last):+.4f})")
    print(f"  SPH mean_x        = {sph_last[0]:.3f}  (err={sph_last[0]-(x0+U*t_last):+.4f})")
    print()
    print(f"  Analytical Var     = {ana_var_last:.5f}")
    print(f"  LDM Var_x          = {ldm_last[4]:.5f}  (ratio={ldm_last[4]/ana_var_last:.3f})")
    print(f"  SPH Var_x          = {sph_last[3]:.5f}  (ratio={sph_last[3]/ana_var_last:.3f})")

    D_num_ldm = (ldm_last[4] - ldm_var0) / (2.0 * t_last)
    D_num_sph = (sph_last[3] - sph_var0) / (2.0 * t_last)
    print(f"\n  Physical D         = {D:.6e}")
    print(f"  LDM effective D_x  = {D_num_ldm:.6e}  (ratio={D_num_ldm/D:.2f}x)")
    print(f"  SPH effective D_x  = {D_num_sph:.6e}  (ratio={D_num_sph/D:.2f}x)")


if __name__ == "__main__":
    vtk_dir = sys.argv[1] if len(sys.argv) > 1 else "Transport_advection"
    U  = float(sys.argv[2]) if len(sys.argv) > 2 else 0.3
    D  = float(sys.argv[3]) if len(sys.argv) > 3 else 0.01
    x0 = float(sys.argv[4]) if len(sys.argv) > 4 else 5.0
    y0 = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0
    z0 = float(sys.argv[6]) if len(sys.argv) > 6 else 3.0
    dt_base = 0.01

    analyze(vtk_dir, U, D, x0, y0, z0, dt_base)
