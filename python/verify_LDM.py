#!/usr/bin/env python3
"""
LDM 검증 스크립트: 시뮬레이션 결과를 해석해와 비교
사용법: python3 verify_LDM.py <vtk_dir> <case>

<case>:
  1  순수 Diffusion  (U=0,   D=0.01)
  2  순수 Advection  (U=0.3, D=0)
  3  Advection+Diffusion (U=0.3, D=0.01)
"""
import numpy as np
import os, sys, re, glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ============================================================
# VTK 파서 (POLYDATA binary)
# ============================================================
def parse_vtk_LDM(filepath):
    """LDM VTK binary (POLYDATA) 형식에서 입자 좌표를 반환"""
    with open(filepath, "rb") as f:
        content = f.read()

    # POINTS 위치 찾기
    m = re.search(rb"POINTS\s+(\d+)\s+(\w+)\n", content)
    if not m:
        return None
    n_pts = int(m.group(1))
    dtype_str = m.group(2).decode("ascii")
    pts_offset = m.end()

    if dtype_str == "float":
        pts = np.frombuffer(content, dtype=">f4", count=n_pts * 3, offset=pts_offset).reshape(n_pts, 3)
    elif dtype_str == "double":
        pts = np.frombuffer(content, dtype=">f8", count=n_pts * 3, offset=pts_offset).reshape(n_pts, 3)
    else:
        pts = np.frombuffer(content, dtype=">f4", count=n_pts * 3, offset=pts_offset).reshape(n_pts, 3)

    return pts.astype(np.float64)


def get_step_from_filename(fpath):
    m = re.search(r"LDM_(\d+)stp", fpath)
    return int(m.group(1)) if m else -1


# ============================================================
# 해석해 (Analytical solutions)
# ============================================================
def analytical_mean(x0, U, t):
    """평균 위치: <x> = x0 + U*t"""
    return x0 + U * t

def analytical_variance(D, t):
    """위치 분산: Var(x) = 2*D*t (각 축 독립)"""
    return 2.0 * D * t

def analytical_concentration_1D(x, x0, U, D, t, N_total=1.0):
    """
    1D 농도 프로필 (x 방향):
    C(x,t) = N / sqrt(4*pi*D*t) * exp(-(x - x0 - U*t)^2 / (4*D*t))
    """
    if D <= 0 or t <= 0:
        return np.zeros_like(x)
    return N_total / np.sqrt(4 * np.pi * D * t) * np.exp(-(x - x0 - U * t)**2 / (4 * D * t))


# ============================================================
# 메인 분석
# ============================================================
def analyze(vtk_dir, case_id, U, D, x0, y0, z0, dt_base, freq_LDM, freq_output):
    """
    case_id: 1=diffusion, 2=advection, 3=both
    """
    files = sorted(glob.glob(os.path.join(vtk_dir, "LDM_[0-9]*stp.vtk")),
                   key=get_step_from_filename)

    if not files:
        print(f"No LDM VTK files found in {vtk_dir}")
        return

    print(f"Found {len(files)} LDM VTK files")
    print(f"Parameters: U={U} m/s, D={D} m²/s, source=({x0},{y0},{z0})")
    print(f"dt={dt_base}s, freq_LDM={freq_LDM}, freq_output={freq_output}")
    print()

    steps = []
    times = []
    mean_x_list, mean_y_list, mean_z_list = [], [], []
    var_x_list, var_y_list, var_z_list = [], [], []
    n_particles = []

    for fpath in files:
        step = get_step_from_filename(fpath)
        t = step * dt_base

        pts = parse_vtk_LDM(fpath)
        if pts is None or len(pts) == 0:
            continue

        N = len(pts)
        mx, my, mz = pts[:, 0].mean(), pts[:, 1].mean(), pts[:, 2].mean()
        vx = pts[:, 0].var()
        vy = pts[:, 1].var()
        vz = pts[:, 2].var()

        steps.append(step)
        times.append(t)
        mean_x_list.append(mx)
        mean_y_list.append(my)
        mean_z_list.append(mz)
        var_x_list.append(vx)
        var_y_list.append(vy)
        var_z_list.append(vz)
        n_particles.append(N)

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

    # 초기 분산 (t=0 데이터에서 측정, 점 소스가 아닌 경우 σ₀² > 0)
    var0_x = var_x[0]
    var0_y = var_y[0]
    var0_z = var_z[0]
    var0_avg = (var0_x + var0_y + var0_z) / 3.0
    print(f"Initial variance (t={times[0]:.2f}s): var0_x={var0_x:.6f}, var0_y={var0_y:.6f}, var0_z={var0_z:.6f}")
    print(f"  -> σ₀ ≈ {np.sqrt(var0_avg):.4f} m")
    print()

    # 해석해: Var(t) = σ₀² + 2Dt
    theo_mean_x = analytical_mean(x0, U, times)
    theo_var = var0_avg + analytical_variance(D, times)

    # ============================================================
    # 결과 출력 (표)
    # ============================================================
    case_names = {1: "Pure Diffusion", 2: "Pure Advection", 3: "Advection+Diffusion"}
    print(f"=== Case {case_id}: {case_names[case_id]} ===\n")

    if case_id in [2, 3]:
        print(f"{'step':>6} | {'time':>6} | {'<x>_sim':>10} | {'<x>_theory':>10} | {'err_x':>10} | {'N':>6}")
        print("-" * 65)
        for i in range(len(steps)):
            err_x = mean_x[i] - theo_mean_x[i]
            print(f"{steps[i]:>6} | {times[i]:>6.2f} | {mean_x[i]:>10.4f} | {theo_mean_x[i]:>10.4f} | {err_x:>+10.4f} | {n_particles[i]:>6}")

    if case_id in [1, 3]:
        print()
        print(f"{'step':>6} | {'time':>6} | {'Var_x sim':>10} | {'Var_x thy':>10} | {'ratio':>7} | {'Var_y sim':>10} | {'Var_z sim':>10}")
        print("-" * 85)
        for i in range(len(steps)):
            ratio = var_x[i] / theo_var[i] if theo_var[i] > 1e-15 else float('nan')
            print(f"{steps[i]:>6} | {times[i]:>6.2f} | {var_x[i]:>10.6f} | {theo_var[i]:>10.6f} | {ratio:>7.3f} | {var_y[i]:>10.6f} | {var_z[i]:>10.6f}")

    # ============================================================
    # CSV 출력 (Origin 용)
    # ============================================================
    csv_path = os.path.join(vtk_dir, f"LDM_case{case_id}_stats.csv")
    with open(csv_path, "w") as f:
        f.write("time,mean_x_sim,mean_x_theory,mean_y_sim,mean_z_sim,var_x_sim,var_y_sim,var_z_sim,var_theory,N\n")
        for i in range(len(times)):
            f.write(f"{times[i]:.4f},{mean_x[i]:.6f},{theo_mean_x[i]:.6f},{mean_y[i]:.6f},{mean_z[i]:.6f},"
                    f"{var_x[i]:.6f},{var_y[i]:.6f},{var_z[i]:.6f},{theo_var[i]:.6f},{n_particles[i]}\n")
    print(f"\nCSV saved: {csv_path}")

    # ============================================================
    # 플롯
    # ============================================================
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f"LDM Verification — Case {case_id}: {case_names[case_id]}\n(U={U}, D={D})", fontsize=14)

    # (0,0) Mean x-position
    ax = axes[0, 0]
    ax.plot(times, mean_x, 'bo-', label='Simulation <x>', markersize=3)
    ax.plot(times, theo_mean_x, 'r--', label=f'Theory: x₀ + U·t', linewidth=2)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Mean x [m]')
    ax.set_title('Mean x-position')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # (0,1) Mean y, z position
    ax = axes[0, 1]
    ax.plot(times, mean_y, 'g^-', label='<y>', markersize=3)
    ax.plot(times, mean_z, 'ms-', label='<z>', markersize=3)
    ax.axhline(y=y0, color='g', linestyle='--', alpha=0.5, label=f'y₀={y0}')
    ax.axhline(y=z0, color='m', linestyle='--', alpha=0.5, label=f'z₀={z0}')
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Mean position [m]')
    ax.set_title('Mean y, z position (should stay constant)')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # (1,0) Variance (x)
    ax = axes[1, 0]
    ax.plot(times, var_x, 'bo-', label='Var(x) sim', markersize=3)
    ax.plot(times, var_y, 'g^-', label='Var(y) sim', markersize=3)
    ax.plot(times, var_z, 'ms-', label='Var(z) sim', markersize=3)
    if D > 0:
        ax.plot(times, theo_var, 'r--', label=f'Theory: 2Dt', linewidth=2)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Variance [m²]')
    ax.set_title('Position variance')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # (1,1) Histogram of x at last timestep + analytical
    ax = axes[1, 1]
    last_file = files[-1]
    last_pts = parse_vtk_LDM(last_file)
    last_step = get_step_from_filename(last_file)
    last_t = last_step * dt_base

    if last_pts is not None and len(last_pts) > 0:
        x_data = last_pts[:, 0]
        ax.hist(x_data, bins=50, density=True, alpha=0.6, color='steelblue', label=f'Simulation (t={last_t:.1f}s)')

        if D > 0:
            x_range = np.linspace(x_data.min() - 0.5, x_data.max() + 0.5, 300)
            C_theory = analytical_concentration_1D(x_range, x0, U, D, last_t)
            ax.plot(x_range, C_theory, 'r-', linewidth=2, label='Analytical')

        ax.set_xlabel('x [m]')
        ax.set_ylabel('Probability density')
        ax.set_title(f'x-distribution at t={last_t:.1f}s')
        ax.legend()
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plot_name = os.path.join(vtk_dir, f"LDM_verify_case{case_id}.png")
    plt.savefig(plot_name, dpi=150)
    print(f"\nPlot saved: {plot_name}")


# ============================================================
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 verify_LDM.py <vtk_dir> <case_id> [U] [D] [x0] [y0] [z0]")
        print("  case_id: 1=diffusion, 2=advection, 3=both")
        sys.exit(1)

    vtk_dir = sys.argv[1]
    case_id = int(sys.argv[2])

    # 기본값 (ISPH_Calc.cuh의 LDM_U_mean, LDM_D와 일치시킬 것)
    defaults = {
        1: (0.0, 0.01),   # pure diffusion
        2: (0.3, 0.0),    # pure advection
        3: (0.3, 0.01),   # advection + diffusion
    }
    U_def, D_def = defaults.get(case_id, (0.3, 0.01))

    U  = float(sys.argv[3]) if len(sys.argv) > 3 else U_def
    D  = float(sys.argv[4]) if len(sys.argv) > 4 else D_def
    x0 = float(sys.argv[5]) if len(sys.argv) > 5 else 5.0
    y0 = float(sys.argv[6]) if len(sys.argv) > 6 else 0.0
    z0 = float(sys.argv[7]) if len(sys.argv) > 7 else 3.0

    dt_base = 0.01       # time step
    freq_LDM = 10        # LDM frequency (from solv file)
    freq_output = 100    # output frequency (from solv file)

    analyze(vtk_dir, case_id, U, D, x0, y0, z0, dt_base, freq_LDM, freq_output)
