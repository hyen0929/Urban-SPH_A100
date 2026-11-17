"""
MOST Analytical Velocity & TKE Profiles
Compare simulation results with Monin-Obukhov similarity theory predictions.

Usage:
  python MOST_analytical_profile.py [vtk_file]
  - Without argument: plots analytical profiles only
  - With VTK file: overlays simulation data
"""

import numpy as np
import matplotlib.pyplot as plt
import sys, os

# ============================================================
# Parameters (match Parameters.cuh)
# ============================================================
u_star  = 0.5       # Friction velocity [m/s]
kappa   = 0.4       # von Karman constant
z0      = 0.1       # Roughness length [m]
z0h     = 0.01      # Roughness length for heat [m]
ABL_h   = 125.0     # ABL height [m]
T_ini   = 300.0     # Initial temperature [K]
g       = 9.81      # Gravitational acceleration

# ============================================================
# Stability correction functions (Businger-Dyer)
# ============================================================
def psi_m(zeta):
    """Momentum stability correction"""
    zeta = np.asarray(zeta, dtype=float)
    result = np.zeros_like(zeta)
    # Unstable (zeta < 0)
    mask = zeta < 0
    x = (1.0 - 16.0 * zeta[mask])**0.25
    result[mask] = (np.pi/2 - 2*np.arctan(x)
                    + np.log((1+x)**2 * (1+x**2) / 8))
    # Stable (zeta >= 0)
    result[~mask] = -5.0 * zeta[~mask]
    return result

def psi_h(zeta):
    """Heat stability correction"""
    zeta = np.asarray(zeta, dtype=float)
    result = np.zeros_like(zeta)
    # Unstable
    mask = zeta < 0
    y = (1.0 - 16.0 * zeta[mask])**0.5
    result[mask] = 2.0 * np.log(0.5 * (1 + y))
    # Stable
    result[~mask] = -5.0 * zeta[~mask]
    return result

# ============================================================
# Analytical profiles
# ============================================================
def u_profile(z, L=None):
    """
    MOST velocity profile: u(z) = (u*/κ) [ln(z/z0) - ψ_m(z/L) + ψ_m(z0/L)]
    L=None or L=inf → neutral
    """
    z = np.maximum(z, z0)
    if L is None or np.isinf(L):
        return (u_star / kappa) * np.log(z / z0)
    else:
        return (u_star / kappa) * (np.log(z / z0) - psi_m(z / L) + psi_m(z0 / L))

def theta_profile(z, L, theta_s, theta_star):
    """
    MOST temperature profile:
    θ(z) = θ_s + (θ*/κ) [ln(z/z0h) - ψ_h(z/L) + ψ_h(z0h/L)]
    """
    z = np.maximum(z, z0h)
    return theta_s + (theta_star / kappa) * (np.log(z / z0h) - psi_h(z / L) + psi_h(z0h / L))

def k_profile_neutral(z):
    """
    Neutral TKE profile: k(z) = u*²/√Cμ × (1 - z/h)²
    Cμ ≈ 0.09 → √Cμ ≈ 0.3
    """
    ratio = np.clip(1.0 - z / ABL_h, 0, 1)
    k_surf = u_star**2 / 0.3
    return k_surf * ratio**2

def nut_profile_neutral(z, delta=2.5):
    """
    Neutral eddy viscosity: ν_t = κ u* z (1 - z/h)
    Or from DDF model: ν_t = Ck Δ √k
    """
    z_clip = np.clip(z, 0, ABL_h)
    return kappa * u_star * z_clip * np.clip(1.0 - z_clip / ABL_h, 0, 1)

def L_from_heat_flux(Q_H):
    """
    Obukhov length from sensible heat flux:
    L = -ρ cp u*³ T / (κ g Q_H)
    Q_H > 0 (upward) → L < 0 (unstable)
    Q_H < 0 (downward) → L > 0 (stable)
    """
    rho = 1.225
    cp = 1005.0
    if abs(Q_H) < 1e-10:
        return np.inf
    return -rho * cp * u_star**3 * T_ini / (kappa * g * Q_H)


# ============================================================
# Read VTK binary (POLYDATA) if provided
# ============================================================
def read_vtk_binary(filepath):
    """Read SOPHIA VTK binary output (FIELD format), return dict of arrays"""
    import struct
    data = {}
    with open(filepath, 'rb') as f:
        # Header (5 lines: version, title, BINARY, DATASET POLYDATA, POINTS)
        for _ in range(4):
            f.readline()

        # POINTS N float
        line = f.readline().decode().strip()
        parts = line.split()
        npts = int(parts[1])
        raw = f.read(npts * 3 * 4)
        coords = np.array(struct.unpack('>' + 'f' * (npts * 3), raw)).reshape(npts, 3)
        data['x'] = coords[:, 0]
        data['y'] = coords[:, 1]
        data['z'] = coords[:, 2]

        # Scan for POINT_DATA
        while True:
            line = f.readline()
            if not line:
                return data
            try:
                line_str = line.decode('ascii', errors='ignore').strip()
            except:
                continue
            if line_str.startswith('POINT_DATA'):
                break

        # Read FIELD header: "FIELD FieldData N"
        line = f.readline().decode().strip()
        parts = line.split()
        nfields = int(parts[2])

        # Read each field: "name ncomp ntuples dtype\n<binary>"
        for _ in range(nfields):
            line = f.readline().decode().strip()
            parts = line.split()
            name = parts[0]
            ncomp = int(parts[1])
            ntuples = int(parts[2])
            dtype = parts[3]  # "float" or "int"

            nbytes = ntuples * ncomp * 4
            raw = f.read(nbytes)

            if dtype == 'float':
                arr = np.array(struct.unpack('>' + 'f' * (ntuples * ncomp), raw))
            else:  # int
                arr = np.array(struct.unpack('>' + 'i' * (ntuples * ncomp), raw))

            if ncomp == 1:
                data[name] = arr
            elif ncomp == 3:
                arr = arr.reshape(ntuples, 3)
                data[name + '_x'] = arr[:, 0]
                data[name + '_y'] = arr[:, 1]
                data[name + '_z'] = arr[:, 2]
            else:
                data[name] = arr.reshape(ntuples, ncomp)

    return data


def extract_profile(data, x_target, y_target=0.0, tol_x=2.0, tol_y=2.0):
    """Extract vertical profile at (x_target, y_target) from VTK data"""
    mask = ((np.abs(data['x'] - x_target) < tol_x) &
            (np.abs(data['y'] - y_target) < tol_y))
    z = data['z'][mask]
    ux = data.get('velocity_x', np.zeros_like(data['x']))[mask]
    speed = np.sqrt(ux**2) if 'velocity_x' in data else None

    # Sort by z
    idx = np.argsort(z)
    result = {'z': z[idx]}
    if 'velocity_x' in data:
        result['ux'] = data['velocity_x'][mask][idx]
    if 'velocity_y' in data:
        result['uy'] = data['velocity_y'][mask][idx]
    if 'k_turb' in data:
        result['k_turb'] = data['k_turb'][mask][idx]
    if 'e_turb' in data:
        result['e_turb'] = data['e_turb'][mask][idx]
    return result


# ============================================================
# Plotting
# ============================================================
def plot_profiles(vtk_file=None, x_locations=None):
    z = np.linspace(0.5, ABL_h, 500)

    # --- Stability cases ---
    cases = {
        'Neutral (L=∞)':         {'L': np.inf,  'color': 'black',  'ls': '-'},
        'Weakly unstable (L=-200m)':  {'L': -200.0, 'color': 'red',    'ls': '--'},
        'Unstable (L=-50m)':     {'L': -50.0,  'color': 'darkred', 'ls': '-.'},
        'Weakly stable (L=200m)':     {'L': 200.0,  'color': 'blue',   'ls': '--'},
        'Stable (L=50m)':        {'L': 50.0,   'color': 'navy',   'ls': '-.'},
    }

    fig, axes = plt.subplots(1, 3, figsize=(15, 7))

    # --- (a) Velocity profiles ---
    ax = axes[0]
    for label, c in cases.items():
        u = u_profile(z, L=c['L'])
        ax.plot(u, z, color=c['color'], ls=c['ls'], lw=1.5, label=label)

    # Overlay simulation data
    if vtk_file and os.path.exists(vtk_file):
        data = read_vtk_binary(vtk_file)
        if x_locations is None:
            x_locations = [50, 150, 250, 400]
        colors_sim = ['green', 'orange', 'purple', 'brown']
        for x_loc, col in zip(x_locations, colors_sim):
            prof = extract_profile(data, x_loc)
            if len(prof['z']) > 0 and 'ux' in prof:
                ax.plot(prof['ux'], prof['z'], 'o', color=col, ms=3,
                        alpha=0.6, label=f'Sim x={x_loc}m')

    ax.set_xlabel('u [m/s]', fontsize=12)
    ax.set_ylabel('z [m]', fontsize=12)
    ax.set_title('(a) Velocity Profile', fontsize=13)
    ax.set_xlim(0, 15)
    ax.set_ylim(0, ABL_h)
    ax.legend(fontsize=8, loc='upper left')
    ax.grid(True, alpha=0.3)

    # --- (b) Velocity in semi-log (z in log scale) ---
    ax = axes[1]
    z_log = np.logspace(np.log10(0.5), np.log10(ABL_h), 500)
    for label, c in cases.items():
        u = u_profile(z_log, L=c['L'])
        ax.semilogx(z_log, u, color=c['color'], ls=c['ls'], lw=1.5, label=label)

    if vtk_file and os.path.exists(vtk_file):
        for x_loc, col in zip(x_locations, colors_sim):
            prof = extract_profile(data, x_loc)
            if len(prof['z']) > 0 and 'ux' in prof:
                ax.semilogx(prof['z'], prof['ux'], 'o', color=col, ms=3,
                            alpha=0.6, label=f'Sim x={x_loc}m')

    ax.set_xlabel('z [m] (log scale)', fontsize=12)
    ax.set_ylabel('u [m/s]', fontsize=12)
    ax.set_title('(b) Velocity (semi-log)', fontsize=13)
    ax.set_ylim(0, 15)
    ax.legend(fontsize=8, loc='upper left')
    ax.grid(True, alpha=0.3, which='both')

    # --- (c) TKE profile (neutral only) ---
    ax = axes[2]
    k_analytical = k_profile_neutral(z)
    ax.plot(k_analytical, z, 'k-', lw=2, label='Analytical (neutral)')

    if vtk_file and os.path.exists(vtk_file):
        for x_loc, col in zip(x_locations, colors_sim):
            prof = extract_profile(data, x_loc)
            if len(prof['z']) > 0 and 'k_turb' in prof:
                ax.plot(prof['k_turb'], prof['z'], 'o', color=col, ms=3,
                        alpha=0.6, label=f'Sim x={x_loc}m')

    ax.set_xlabel('k [m²/s²]', fontsize=12)
    ax.set_ylabel('z [m]', fontsize=12)
    ax.set_title('(c) TKE Profile (neutral)', fontsize=13)
    ax.set_xlim(0, 1.5)
    ax.set_ylim(0, ABL_h)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    plt.suptitle(f'MOST Analytical Profiles  (u*={u_star}, z₀={z0}, h={ABL_h}m)',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('MOST_profiles.png', dpi=200, bbox_inches='tight')
    plt.show()
    print("Saved: MOST_profiles.png")


def print_table():
    """Print reference values at selected heights"""
    heights = [1.25, 2.5, 5, 10, 20, 50, 75, 100, 125]
    Ls = [np.inf, -200, -50, 200, 50]
    L_names = ['Neutral', 'L=-200', 'L=-50', 'L=200', 'L=50']

    print("\n" + "="*85)
    print(f"  MOST Velocity Profile [m/s]  (u*={u_star}, κ={kappa}, z0={z0}m)")
    print("="*85)
    header = f"{'z [m]':>8s}"
    for name in L_names:
        header += f"  {name:>12s}"
    header += f"  {'k_neutral':>10s}"
    print(header)
    print("-"*85)

    for z_val in heights:
        row = f"{z_val:8.2f}"
        for L in Ls:
            u = u_profile(np.array([z_val]), L=L)[0]
            row += f"  {u:12.4f}"
        k = k_profile_neutral(np.array([z_val]))[0]
        row += f"  {k:10.4f}"
        print(row)
    print("="*85)

    print(f"\n  Key reference values:")
    print(f"  u(z=2.5m, neutral)  = {u_profile(np.array([2.5]))[0]:.4f} m/s")
    print(f"  u(z=125m, neutral)  = {u_profile(np.array([125.0]))[0]:.4f} m/s  (geostrophic)")
    print(f"  k(z=2.5m, neutral)  = {k_profile_neutral(np.array([2.5]))[0]:.4f} m²/s²")
    print(f"  ν_t(z=2.5m)         = {nut_profile_neutral(np.array([2.5]))[0]:.4f} m²/s")
    print(f"  f_geo = u*²/h       = {u_star**2/ABL_h:.6f} m/s²")
    print()


# ============================================================
if __name__ == '__main__':
    print_table()

    vtk_file = sys.argv[1] if len(sys.argv) > 1 else None
    if vtk_file:
        print(f"Loading VTK: {vtk_file}")
    plot_profiles(vtk_file)
