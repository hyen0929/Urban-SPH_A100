"""
Flag interior fluid particles in SPH input file using STL ray casting.

Reads the SPH input file and STL building mesh. For each fluid particle
(p_type==1), checks whether it lies inside a building by casting a ray
upward and counting roof hits. Interior particles get p_type changed to 99.

The CUDA code then converts p_type=99 → buffer_type=99, p_type=1 at init.

Usage:
    python flag_interior_fluid.py <input.txt> <building.stl> [output.txt] [dx]

    If output.txt is omitted, overwrites input.txt.
    dx (default=2.5) is the particle spacing used to build a spatial hash
    for fast lookup of interior positions.
"""
import numpy as np
import trimesh
import sys
import time
import os


def load_interior_positions(stl_path, dx, min_z_layers=3):
    """
    Generate set of interior grid positions using z-direction ray casting.
    Returns a dict: (ix, iy) -> z_top  for fast interior checking.
    """
    print(f"[*] Loading STL: {stl_path}")
    mesh = trimesh.load(stl_path)
    print(f"[*] Bounds: {mesh.bounds[0]} to {mesh.bounds[1]}")

    xmin, ymin, zmin = mesh.bounds[0]
    xmax, ymax, zmax = mesh.bounds[1]

    # Build 2D grid
    xs = np.arange(xmin + dx / 2, xmax, dx)
    ys = np.arange(ymin + dx / 2, ymax, dx)
    print(f"[*] XY grid: {len(xs)} x {len(ys)} = {len(xs)*len(ys)} columns")

    xx, yy = np.meshgrid(xs, ys)
    xx = xx.ravel()
    yy = yy.ravel()
    n_cols = len(xx)

    # Ray cast downward from above
    ray_z = zmax + 10.0
    origins = np.column_stack([xx, yy, np.full(n_cols, ray_z)])
    directions = np.tile([0.0, 0.0, -1.0], (n_cols, 1))

    print(f"[*] Casting {n_cols} rays ...")
    t0 = time.time()
    hit_points, ray_indices, _ = mesh.ray.intersects_location(
        ray_origins=origins, ray_directions=directions
    )
    print(f"[*] Ray casting: {time.time()-t0:.2f}s, {len(hit_points)} hits")

    if len(hit_points) == 0:
        print("[!] No hits — check STL.")
        return {}

    # For each column, find highest hit z (= roof)
    roof_z = {}
    for pt, ri in zip(hit_points, ray_indices):
        z_hit = pt[2]
        if ri not in roof_z or z_hit > roof_z[ri]:
            roof_z[ri] = z_hit

    # Build set of interior (x, y, z_top) per column
    # Use snapped grid coordinates as keys
    interior = {}  # (round_x, round_y) -> z_top
    for ri, z_top in roof_z.items():
        x_col = float(xx[ri])
        y_col = float(yy[ri])
        interior[(round(x_col, 4), round(y_col, 4))] = z_top

    print(f"[*] {len(interior)} interior columns found")
    return interior


def flag_interior_in_sph(input_path, stl_path, output_path, dx=2.5):
    """
    Read SPH input, flag interior fluid particles as p_type=99.
    Uses a KD-tree of interior column centers for fast spatial lookup.
    """
    # Build interior column info from STL
    interior_cols = load_interior_positions(stl_path, dx)
    if not interior_cols:
        print("[!] No interior columns. Aborting.")
        return

    # Build arrays for vectorized lookup
    col_xy = np.array(list(interior_cols.keys()))  # (N_cols, 2)
    col_ztop = np.array(list(interior_cols.values()))  # (N_cols,)

    from scipy.spatial import cKDTree
    tree = cKDTree(col_xy)
    search_radius = dx * 0.6  # slightly more than dx/2 to catch snapping errors

    # Parse header
    print(f"[*] Reading SPH input: {input_path}")
    t0 = time.time()

    with open(input_path, 'r') as f:
        header = f.readline().strip()
        lines = f.readlines()

    col_indices = list(map(int, header.split()))
    # Find which data column holds p_type (inp_ptype = 8)
    try:
        ptype_col = col_indices.index(8)
    except ValueError:
        print("[!] Column 8 (p_type) not found in header!")
        return

    # Find x, y, z columns (inp_x=1, inp_y=2, inp_z=3)
    x_col_idx = col_indices.index(1)
    y_col_idx = col_indices.index(2)
    z_col_idx = col_indices.index(3)

    print(f"[*] {len(lines)} particles, header columns: {col_indices}")
    print(f"[*] p_type at data col {ptype_col}, xyz at cols {x_col_idx},{y_col_idx},{z_col_idx}")

    # Parse all particle data
    n_flagged = 0
    n_fluid = 0
    output_lines = [header + '\n']

    # Process in chunks for speed
    # First pass: extract fluid positions
    fluid_indices = []
    fluid_xy = []
    fluid_z = []

    for i, line in enumerate(lines):
        fields = line.split('\t')
        ptype = int(float(fields[ptype_col]))
        if ptype == 1:
            x = float(fields[x_col_idx])
            y = float(fields[y_col_idx])
            z = float(fields[z_col_idx])
            fluid_indices.append(i)
            fluid_xy.append([x, y])
            fluid_z.append(z)

    n_fluid = len(fluid_indices)
    print(f"[*] Parsing done in {time.time()-t0:.1f}s. {n_fluid} fluid particles.")

    if n_fluid == 0:
        print("[!] No fluid particles found.")
        return

    fluid_xy = np.array(fluid_xy)
    fluid_z = np.array(fluid_z)

    # Query KD-tree: for each fluid particle, find nearest interior column
    print("[*] Querying KD-tree for interior matches ...")
    t0 = time.time()
    dists, idxs = tree.query(fluid_xy, k=1, workers=-1)

    # Flag: within search_radius AND z < z_top of that column
    is_interior = np.zeros(n_fluid, dtype=bool)
    for fi in range(n_fluid):
        if dists[fi] < search_radius:
            ci = idxs[fi]
            if fluid_z[fi] < col_ztop[ci] - dx * 0.3:
                is_interior[fi] = True

    n_flagged = np.sum(is_interior)
    print(f"[*] KD-tree query: {time.time()-t0:.1f}s. Flagged {n_flagged} interior fluid particles.")

    # Build flag set for the fluid indices
    flag_set = set()
    for fi in range(n_fluid):
        if is_interior[fi]:
            flag_set.add(fluid_indices[fi])

    # Write output
    print(f"[*] Writing output: {output_path}")
    t0 = time.time()
    with open(output_path, 'w') as f:
        f.write(header + '\n')
        for i, line in enumerate(lines):
            if i in flag_set:
                fields = line.split('\t')
                fields[ptype_col] = '99'
                f.write('\t'.join(fields))
            else:
                f.write(line)

    print(f"[*] Writing done in {time.time()-t0:.1f}s")
    print(f"[*] Summary: {n_flagged}/{n_fluid} fluid particles flagged as interior (p_type=99)")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python flag_interior_fluid.py <input.txt> <building.stl> [output.txt] [dx]")
        sys.exit(1)

    input_file = sys.argv[1]
    stl_file = sys.argv[2]
    output_file = sys.argv[3] if len(sys.argv) >= 4 else input_file
    dx = float(sys.argv[4]) if len(sys.argv) >= 5 else 2.5

    flag_interior_in_sph(input_file, stl_file, output_file, dx)
