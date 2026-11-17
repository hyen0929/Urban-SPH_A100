"""Analyze VTK binary files to find where velocity blowup occurs."""
import struct
import numpy as np
import sys
import glob
import os

def read_vtk_binary(filepath):
    """Read VTK binary POLYDATA file, extract positions and scalar/vector fields."""
    with open(filepath, 'rb') as f:
        data = f.read()

    # Find POINTS section
    idx = data.find(b'POINTS')
    if idx < 0:
        return None
    header_end = data.index(b'\n', idx) + 1
    header_line = data[idx:header_end].decode().strip()
    parts = header_line.split()
    n_points = int(parts[1])
    dtype = parts[2]  # float or double

    if n_points == 0:
        return None

    # Read positions (3 floats per point, big-endian)
    pos_start = header_end
    pos_bytes = n_points * 3 * 4  # float32
    pos_data = np.frombuffer(data[pos_start:pos_start+pos_bytes], dtype='>f4').reshape(n_points, 3)

    result = {'n': n_points, 'x': pos_data[:,0], 'y': pos_data[:,1], 'z': pos_data[:,2]}

    # Find POINT_DATA and read fields
    pd_idx = data.find(b'POINT_DATA', pos_start + pos_bytes)
    if pd_idx < 0:
        return result

    # Parse fields after POINT_DATA
    cursor = data.index(b'\n', pd_idx) + 1
    while cursor < len(data):
        line_end = data.find(b'\n', cursor)
        if line_end < 0:
            break
        line = data[cursor:line_end].decode(errors='ignore').strip()

        if line.startswith('SCALARS'):
            parts = line.split()
            name = parts[1]
            cursor = line_end + 1
            lt_end = data.find(b'\n', cursor)
            cursor = lt_end + 1
            nbytes = n_points * 4
            if cursor + nbytes <= len(data):
                arr = np.frombuffer(data[cursor:cursor+nbytes], dtype='>f4')
                result[name] = arr
            cursor += nbytes

        elif line.startswith('VECTORS'):
            parts = line.split()
            name = parts[1]
            cursor = line_end + 1
            nbytes = n_points * 3 * 4
            if cursor + nbytes <= len(data):
                arr = np.frombuffer(data[cursor:cursor+nbytes], dtype='>f4').reshape(n_points, 3)
                result[name + '_x'] = arr[:,0]
                result[name + '_y'] = arr[:,1]
                result[name + '_z'] = arr[:,2]
            cursor += nbytes

        elif line.startswith('FIELD') or line.startswith('POINT_DATA'):
            cursor = line_end + 1

        else:
            # Try FIELD FieldData format: "name\tn_comp\tn_tuples\tdtype"
            parts = line.split('\t')
            if len(parts) >= 4:
                name = parts[0].strip()
                try:
                    n_comp = int(parts[1])
                    n_tuples = int(parts[2])
                    dtype_str = parts[3].strip()
                except (ValueError, IndexError):
                    cursor = line_end + 1
                    continue
                cursor = line_end + 1
                if dtype_str == 'float':
                    nbytes = n_tuples * n_comp * 4
                    if cursor + nbytes <= len(data):
                        arr = np.frombuffer(data[cursor:cursor+nbytes], dtype='>f4')
                        if n_comp == 1:
                            result[name] = arr
                        elif n_comp == 3:
                            arr = arr.reshape(n_tuples, 3)
                            result[name+'_x'] = arr[:,0]
                            result[name+'_y'] = arr[:,1]
                            result[name+'_z'] = arr[:,2]
                    cursor += nbytes
                elif dtype_str == 'int':
                    nbytes = n_tuples * n_comp * 4
                    if cursor + nbytes <= len(data):
                        arr = np.frombuffer(data[cursor:cursor+nbytes], dtype='>i4')
                        if n_comp == 1:
                            result[name] = arr
                    cursor += nbytes
                else:
                    cursor = line_end + 1
            else:
                cursor = line_end + 1

    return result


def analyze_file(filepath):
    """Print statistics for a VTK file."""
    basename = os.path.basename(filepath)
    size = os.path.getsize(filepath)

    if size < 500:
        print(f"{basename:40s}  EMPTY ({size} bytes)")
        return

    d = read_vtk_binary(filepath)
    if d is None or d['n'] == 0:
        print(f"{basename:40s}  No particles")
        return

    n = d['n']
    vmag = np.sqrt(d['x']*0 + 0)  # placeholder

    # Check for velocity fields
    if 'velocity_x' in d:
        ux, uy, uz = d['velocity_x'], d['velocity_y'], d['velocity_z']
        vmag = np.sqrt(ux**2 + uy**2 + uz**2)
    elif 'vorticity_x' in d:
        wx, wy, wz = d['vorticity_x'], d['vorticity_y'], d['vorticity_z']
        vmag = np.sqrt(wx**2 + wy**2 + wz**2)

    # Find max velocity particle
    imax = np.argmax(vmag)

    print(f"{basename:40s}  N={n:>8d}  |v|_max={vmag[imax]:>10.3f}"
          f"  @ ({d['x'][imax]:.1f}, {d['y'][imax]:.1f}, {d['z'][imax]:.1f})")

    # Count particles with extreme values
    n_over5 = np.sum(vmag > 5)
    n_over10 = np.sum(vmag > 10)
    n_over20 = np.sum(vmag > 20)
    if n_over5 > 0:
        print(f"{'':40s}  |v|>5: {n_over5}  |v|>10: {n_over10}  |v|>20: {n_over20}")

        mask = vmag > 5
        print(f"{'':40s}  High-v region: "
              f"x=[{d['x'][mask].min():.1f},{d['x'][mask].max():.1f}] "
              f"y=[{d['y'][mask].min():.1f},{d['y'][mask].max():.1f}] "
              f"z=[{d['z'][mask].min():.1f},{d['z'][mask].max():.1f}]")

        # Z distribution of high-velocity particles
        z_vals = d['z'][mask]
        z_bins = [0, 5, 10, 20, 50, 100, 500]
        for j in range(len(z_bins)-1):
            cnt = np.sum((z_vals >= z_bins[j]) & (z_vals < z_bins[j+1]))
            if cnt > 0:
                print(f"{'':40s}    z=[{z_bins[j]},{z_bins[j+1]}): {cnt}")

    # Ground-level check
    ground = (d['z'] >= 0) & (d['z'] < 5)
    if np.sum(ground) > 0:
        vmag_g = vmag[ground]
        ig = np.argmax(vmag_g)
        gi = np.where(ground)[0][ig]
        if vmag_g[ig] > 1.0:
            print(f"{'':40s}  Ground(z<5): max|v|={vmag_g[ig]:.3f} @ ({d['x'][gi]:.1f},{d['y'][gi]:.1f},{d['z'][gi]:.1f})")

    # Pressure
    if 'pressure' in d:
        p = d['pressure']
        print(f"{'':40s}  pres: [{p.min():.1f}, {p.max():.1f}] std={p.std():.2f}")

    # buffer_type check
    if 'buffer_type' in d:
        bt = d['buffer_type']
        if 99 in bt:
            mask99 = bt == 99
            v99 = vmag[mask99] if len(vmag) == len(bt) else np.zeros(1)
            print(f"{'':40s}  interior(bt=99): N={np.sum(mask99)} max|v|={v99.max():.4f}")


if __name__ == "__main__":
    vtk_dir = sys.argv[1] if len(sys.argv) > 1 else "/home/kdh/Urban-SPH/Leipzig_point"
    pattern = sys.argv[2] if len(sys.argv) > 2 else "vorticity_*stp.vtk"

    files = sorted(glob.glob(os.path.join(vtk_dir, pattern)),
                   key=lambda f: int(f.split('_')[-1].replace('stp.vtk','')))

    print(f"Analyzing {len(files)} files in {vtk_dir}")
    print("=" * 100)
    for f in files:
        analyze_file(f)
