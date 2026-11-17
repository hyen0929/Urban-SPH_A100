#!/usr/bin/env python3
"""
Eulerian scalar transport 검증 스크립트: 시뮬레이션 결과를 해석해와 비교
사용법: python3 verify_scalar.py <vtk_dir> [U] [D] [x0] [y0] [z0]

Parameters:
  U  = 0.3 m/s (uniform advection velocity in x)
  D  = 1e-7 m²/s (molecular diffusion, essentially zero)
  Source at (5.0, 0.0, 3.0)
"""
import numpy as np
import os, sys, re, glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DTYPE_MAP = {
    "float":  (">f4", 4),
    "double": (">f8", 8),
    "int":    (">i4", 4),
    "unsigned_int": (">u4", 4),
    "short":  (">i2", 2),
    "char":   (">i1", 1),
}

# ============================================================
# VTK parser (POLYDATA binary FIELD)
# ============================================================
def parse_vtk_fluid(filepath, fields_to_read=("concn",)):
    """
    Fluid VTK binary (POLYDATA + FIELD) 형식에서
    points와 요청된 field 배열을 반환
    """
    with open(filepath, "rb") as f:
        content = f.read()

    # --- POINTS ---
    m = re.search(rb"POINTS\s+(\d+)\s+(\w+)\n", content)
    if not m:
        return None, {}
    n_pts = int(m.group(1))
    dtype_str = m.group(2).decode("ascii")
    pts_offset = m.end()

    dt_pts, sz_pts = DTYPE_MAP.get(dtype_str, (">f4", 4))
    pts = np.frombuffer(content, dtype=dt_pts, count=n_pts*3, offset=pts_offset).reshape(n_pts, 3)

    # --- FIELD FieldData ---
    m_field = re.search(rb"FIELD FieldData\s+(\d+)\n", content)
    if not m_field:
        return pts.astype(np.float64), {}

    n_fields = int(m_field.group(1))
    offset = m_field.end()
    result = {}

    for _ in range(n_fields):
        newline = content.index(b"\n", offset)
        header = content[offset:newline].decode("ascii").split()
        fname = header[0]
        ncomp = int(header[1])
        ntuple = int(header[2])
        dtype_s = header[3]
        offset = newline + 1

        dt, itemsize = DTYPE_MAP.get(dtype_s, (">f4", 4))
        nbytes = ncomp * ntuple * itemsize

        if fname in fields_to_read:
            arr = np.frombuffer(content, dtype=dt, count=ncomp*ntuple, offset=offset)
            if ncomp > 1:
                arr = arr.reshape(ntuple, ncomp)
            result[fname] = arr.astype(np.float64)

        offset += nbytes

    return pts.astype(np.float64), result


def get_step(fpath):
    m = re.search(r"fluid_(\d+)stp", fpath)
    return int(m.group(1)) if m else -1


# ============================================================
# Analytical solutions
# ============================================================
def analytical_mean(x0, U, t):
    return x0 + U * t

def analytical_variance(D, t):
    return 2.0 * D * t


# ============================================================
# Main analysis
# ============================================================
def analyze(vtk_dir, U, D, x0, y0, z0, dt_base, freq_output):
    files = sorted(glob.glob(os.path.join(vtk_dir, "fluid_[0-9]*stp.vtk")),
                   key=get_step)
    if not files:
        print(f"No fluid VTK files found in {vtk_dir}")
        return

    print(f"Found {len(files)} fluid VTK files")
    print(f"Parameters: U={U} m/s, D={D} m²/s, source=({x0},{y0},{z0})")
    print(f"dt={dt_base}s, freq_output={freq_output}")
    print()

    steps = []
    times = []
    mean_x_list, mean_y_list, mean_z_list = [], [], []
    var_x_list, var_y_list, var_z_list = [], [], []
    total_mass_list = []
    peak_concn_list = []
    peak_x_list = []

    for fpath in files:
        step = get_step(fpath)
        t = step * dt_base

        pts, fields = parse_vtk_fluid(fpath, fields_to_read=("concn",))
        if pts is None or "concn" not in fields:
            continue

        concn = fields["concn"]
        total_C = concn.sum()

        if total_C < 1e-10:
            continue

        # concentration-weighted statistics
        w = concn / total_C  # normalized weights

        mx = np.sum(w * pts[:, 0])
        my = np.sum(w * pts[:, 1])
        mz = np.sum(w * pts[:, 2])

        vx = np.sum(w * (pts[:, 0] - mx)**2)
        vy = np.sum(w * (pts[:, 1] - my)**2)
        vz = np.sum(w * (pts[:, 2] - mz)**2)

        # peak tracking
        idx_peak = np.argmax(concn)
        peak_c = concn[idx_peak]
        peak_x = pts[idx_peak, 0]

        steps.append(step)
        times.append(t)
        mean_x_list.append(mx)
        mean_y_list.append(my)
        mean_z_list.append(mz)
        var_x_list.append(vx)
        var_y_list.append(vy)
        var_z_list.append(vz)
        total_mass_list.append(total_C)
        peak_concn_list.append(peak_c)
        peak_x_list.append(peak_x)

    if not steps:
        print("No valid data.")
        return

    times = np.array(times)
    mean_x = np.array(mean_x_list)
    mean_y = np.array(mean_y_list)
    mean_z = np.array(mean_z_list)
    var_x = np.array(var_x_list)
    var_y = np.array(var_y_list)
    var_z = np.array(var_z_list)
    total_mass = np.array(total_mass_list)
    peak_concn = np.array(peak_concn_list)
    peak_x = np.array(peak_x_list)

    # Initial variance (from first valid frame)
    var0_x = var_x[0]
    var0_y = var_y[0]
    var0_z = var_z[0]
    print(f"Initial state (t={times[0]:.2f}s, step={steps[0]}):")
    print(f"  mean = ({mean_x[0]:.4f}, {mean_y[0]:.4f}, {mean_z[0]:.4f})")
    print(f"  var  = ({var0_x:.6f}, {var0_y:.6f}, {var0_z:.6f})")
    print(f"  σ₀   = ({np.sqrt(var0_x):.4f}, {np.sqrt(var0_y):.4f}, {np.sqrt(var0_z):.4f})")
    print(f"  total_C = {total_mass[0]:.2f}, peak = {peak_concn[0]:.4f}")
    print()

    # Analytical: Var(t) = σ₀² + 2Dt
    theo_mean_x = analytical_mean(x0, U, times)
    theo_var_x = var0_x + analytical_variance(D, times)
    theo_var_y = var0_y + analytical_variance(D, times)
    theo_var_z = var0_z + analytical_variance(D, times)

    # ============================================================
    # Table output
    # ============================================================
    print("=== Eulerian Scalar Transport Verification ===\n")

    # Mean position
    print(f"{'step':>6} | {'time':>6} | {'<x>_sim':>10} | {'<x>_thy':>10} | {'err_x':>10} | "
          f"{'peak_x':>8} | {'peak_C':>10} | {'total_C':>10}")
    print("-" * 95)
    for i in range(len(steps)):
        err_x = mean_x[i] - theo_mean_x[i]
        print(f"{steps[i]:>6} | {times[i]:>6.2f} | {mean_x[i]:>10.4f} | {theo_mean_x[i]:>10.4f} | "
              f"{err_x:>+10.4f} | {peak_x[i]:>8.4f} | {peak_concn[i]:>10.4f} | {total_mass[i]:>10.2f}")

    # Variance
    print()
    print(f"{'step':>6} | {'time':>6} | {'Var_x':>10} | {'Var_x thy':>10} | "
          f"{'Var_y':>10} | {'Var_y thy':>10} | {'Var_z':>10} | {'Var_z thy':>10}")
    print("-" * 100)
    for i in range(len(steps)):
        print(f"{steps[i]:>6} | {times[i]:>6.2f} | {var_x[i]:>10.6f} | {theo_var_x[i]:>10.6f} | "
              f"{var_y[i]:>10.6f} | {theo_var_y[i]:>10.6f} | {var_z[i]:>10.6f} | {theo_var_z[i]:>10.6f}")

    # Numerical diffusion estimate
    print("\n=== Numerical Diffusion Estimate ===")
    print("D_numerical = (Var(t) - Var(0)) / (2t)\n")
    print(f"{'step':>6} | {'time':>6} | {'D_num_x':>12} | {'D_num_y':>12} | {'D_num_z':>12} | {'D_physical':>12}")
    print("-" * 75)
    for i in range(len(steps)):
        t = times[i]
        if t > 1e-10:
            Dnx = (var_x[i] - var0_x) / (2.0 * t)
            Dny = (var_y[i] - var0_y) / (2.0 * t)
            Dnz = (var_z[i] - var0_z) / (2.0 * t)
        else:
            Dnx = Dny = Dnz = 0.0
        print(f"{steps[i]:>6} | {times[i]:>6.2f} | {Dnx:>12.6e} | {Dny:>12.6e} | {Dnz:>12.6e} | {D:>12.6e}")

    # ============================================================
    # CSV output
    # ============================================================
    csv_path = os.path.join(vtk_dir, "scalar_transport_stats.csv")
    with open(csv_path, "w") as f:
        f.write("time,mean_x_sim,mean_x_theory,mean_y_sim,mean_z_sim,"
                "var_x_sim,var_y_sim,var_z_sim,var_x_theory,var_y_theory,var_z_theory,"
                "peak_x,peak_concn,total_concn\n")
        for i in range(len(times)):
            f.write(f"{times[i]:.4f},{mean_x[i]:.6f},{theo_mean_x[i]:.6f},"
                    f"{mean_y[i]:.6f},{mean_z[i]:.6f},"
                    f"{var_x[i]:.6f},{var_y[i]:.6f},{var_z[i]:.6f},"
                    f"{theo_var_x[i]:.6f},{theo_var_y[i]:.6f},{theo_var_z[i]:.6f},"
                    f"{peak_x[i]:.6f},{peak_concn[i]:.6f},{total_mass[i]:.2f}\n")
    print(f"\nCSV saved: {csv_path}")

    # ============================================================
    # Plots
    # ============================================================
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle(f"Eulerian Scalar Transport Verification\n(U={U}, D={D})", fontsize=14)

    # (0,0) Mean x-position
    ax = axes[0, 0]
    ax.plot(times, mean_x, 'bo-', label='Simulation <x>', markersize=3)
    ax.plot(times, theo_mean_x, 'r--', label='Theory: x₀ + U·t', linewidth=2)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Mean x [m]')
    ax.set_title('Concentration-weighted mean x')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # (0,1) Peak x-position
    ax = axes[0, 1]
    ax.plot(times, peak_x, 'go-', label='Peak x (sim)', markersize=3)
    ax.plot(times, theo_mean_x, 'r--', label='Theory: x₀ + U·t', linewidth=2)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Peak x [m]')
    ax.set_title('Peak concentration x-position')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # (0,2) Mean y, z
    ax = axes[0, 2]
    ax.plot(times, mean_y, 'g^-', label='<y>', markersize=3)
    ax.plot(times, mean_z, 'ms-', label='<z>', markersize=3)
    ax.axhline(y=y0, color='g', linestyle='--', alpha=0.5, label=f'y₀={y0}')
    ax.axhline(y=z0, color='m', linestyle='--', alpha=0.5, label=f'z₀={z0}')
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Mean position [m]')
    ax.set_title('Mean y, z (should stay constant)')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # (1,0) Variance x
    ax = axes[1, 0]
    ax.plot(times, var_x, 'bo-', label='Var(x) sim', markersize=3)
    ax.plot(times, theo_var_x, 'r--', label=f'Theory: σ₀²+2Dt', linewidth=2)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Variance [m²]')
    ax.set_title('Variance in x (numerical diffusion)')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # (1,1) Variance y, z
    ax = axes[1, 1]
    ax.plot(times, var_y, 'g^-', label='Var(y) sim', markersize=3)
    ax.plot(times, var_z, 'ms-', label='Var(z) sim', markersize=3)
    ax.plot(times, theo_var_y, 'r--', label='Theory', linewidth=2)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Variance [m²]')
    ax.set_title('Variance in y, z')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # (1,2) Peak concentration & total mass
    ax = axes[1, 2]
    ax2 = ax.twinx()
    l1 = ax.plot(times, peak_concn, 'bo-', label='Peak C', markersize=3)
    l2 = ax2.plot(times, total_mass, 'r^-', label='Total C', markersize=3)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Peak concentration', color='b')
    ax2.set_ylabel('Total concentration (mass)', color='r')
    ax.set_title('Peak & Total concentration')
    lines = l1 + l2
    labels = [l.get_label() for l in lines]
    ax.legend(lines, labels)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plot_name = os.path.join(vtk_dir, "scalar_transport_verify.png")
    plt.savefig(plot_name, dpi=150)
    print(f"Plot saved: {plot_name}")

    # ============================================================
    # Summary statistics
    # ============================================================
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)

    # Mean position error
    if len(times) > 1:
        err_mean_x = mean_x - theo_mean_x
        mask = times > 0
        if mask.any():
            rms_err = np.sqrt(np.mean(err_mean_x[mask]**2))
            max_err = np.max(np.abs(err_mean_x[mask]))
            rel_travel = theo_mean_x[-1] - x0
            print(f"\nMean x-position:")
            print(f"  RMS error     = {rms_err:.6f} m")
            print(f"  Max |error|   = {max_err:.6f} m")
            print(f"  Total travel  = {rel_travel:.4f} m")
            if rel_travel > 0:
                print(f"  Relative error= {max_err/rel_travel*100:.2f}%")

    # Numerical diffusion
    if len(times) > 2:
        t_last = times[-1]
        if t_last > 0:
            Dnx_final = (var_x[-1] - var0_x) / (2.0 * t_last)
            Dny_final = (var_y[-1] - var0_y) / (2.0 * t_last)
            Dnz_final = (var_z[-1] - var0_z) / (2.0 * t_last)
            print(f"\nNumerical diffusion (from variance growth):")
            print(f"  D_num_x = {Dnx_final:.6e} m²/s")
            print(f"  D_num_y = {Dny_final:.6e} m²/s")
            print(f"  D_num_z = {Dnz_final:.6e} m²/s")
            print(f"  D_physical = {D:.6e} m²/s")
            if D > 0:
                print(f"  Ratio D_num_x/D_phys = {Dnx_final/D:.1f}x")
            print(f"\n  Var_x growth: {var0_x:.6f} → {var_x[-1]:.6f} m² (+{var_x[-1]-var0_x:.6f})")
            print(f"  Var_y growth: {var0_y:.6f} → {var_y[-1]:.6f} m² (+{var_y[-1]-var0_y:.6f})")
            print(f"  Var_z growth: {var0_z:.6f} → {var_z[-1]:.6f} m² (+{var_z[-1]-var0_z:.6f})")

    # Peak decay
    if len(peak_concn) > 1 and peak_concn[0] > 0:
        print(f"\nPeak concentration decay:")
        print(f"  Initial peak = {peak_concn[0]:.4f}")
        print(f"  Final peak   = {peak_concn[-1]:.4f}")
        print(f"  Ratio        = {peak_concn[-1]/peak_concn[0]:.4f} ({peak_concn[-1]/peak_concn[0]*100:.1f}%)")

    # Mass conservation
    if len(total_mass) > 1 and total_mass[0] > 0:
        print(f"\nMass conservation:")
        print(f"  Initial total = {total_mass[0]:.2f}")
        print(f"  Final total   = {total_mass[-1]:.2f}")
        print(f"  Change        = {(total_mass[-1]-total_mass[0])/total_mass[0]*100:+.2f}%")


# ============================================================
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 verify_scalar.py <vtk_dir> [U] [D] [x0] [y0] [z0]")
        sys.exit(1)

    vtk_dir = sys.argv[1]
    U  = float(sys.argv[2]) if len(sys.argv) > 2 else 0.3
    D  = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-7
    x0 = float(sys.argv[4]) if len(sys.argv) > 4 else 5.0
    y0 = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0
    z0 = float(sys.argv[6]) if len(sys.argv) > 6 else 3.0

    dt_base = 0.01
    freq_output = 100

    analyze(vtk_dir, U, D, x0, y0, z0, dt_base, freq_output)
