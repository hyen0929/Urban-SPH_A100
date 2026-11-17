import numpy as np
import trimesh
import time
import sys


def generate_interior_particles(stl_path, dx, min_z_layers=3,
                                output_filename="interior_particles.txt"):
    """
    Generate interior building particles using z-direction ray casting.
    Works with open-bottom (non-watertight) building STLs.

    For each (x,y) grid point, cast a ray downward from above.
    If the ray hits the mesh, fill from z=0 up to the hit point.
    Columns with fewer than min_z_layers particles are discarded.
    """
    print(f"[*] Loading STL file: {stl_path} ...")
    mesh = trimesh.load(stl_path)

    print(f"[*] Watertight: {mesh.is_watertight}")
    print(f"[*] Bounds (min): {mesh.bounds[0]}")
    print(f"[*] Bounds (max): {mesh.bounds[1]}")

    xmin, ymin, zmin = mesh.bounds[0]
    xmax, ymax, zmax = mesh.bounds[1]

    # Build 2D grid of (x, y) candidate points
    xs = np.arange(xmin + dx / 2, xmax, dx)
    ys = np.arange(ymin + dx / 2, ymax, dx)
    print(f"[*] XY grid: {len(xs)} x {len(ys)} = {len(xs)*len(ys)} columns")

    xx, yy = np.meshgrid(xs, ys)
    xx = xx.ravel()
    yy = yy.ravel()
    n_cols = len(xx)

    # Ray origins: start above the mesh, shoot downward (-z)
    ray_z = zmax + 10.0
    origins = np.column_stack([xx, yy, np.full(n_cols, ray_z)])
    directions = np.tile([0.0, 0.0, -1.0], (n_cols, 1))

    print(f"[*] Casting {n_cols} rays downward ...")
    start_time = time.time()

    # ray.intersects_location returns (hit_points, ray_indices, face_indices)
    hit_points, ray_indices, _ = mesh.ray.intersects_location(
        ray_origins=origins,
        ray_directions=directions
    )

    elapsed = time.time() - start_time
    print(f"[*] Ray casting complete in {elapsed:.2f} sec. Hits: {len(hit_points)}")

    if len(hit_points) == 0:
        print("[!] No ray hits. Check STL orientation/position.")
        return None

    # For each column (ray), find the highest hit z (= roof/top of building)
    # A column is "inside building" if a ray hits the mesh
    roof_z = {}
    for pt, ri in zip(hit_points, ray_indices):
        z_hit = pt[2]
        if ri not in roof_z:
            roof_z[ri] = z_hit
        else:
            # Take the highest hit (roof), since ray goes downward
            roof_z[ri] = max(roof_z[ri], z_hit)

    # Generate interior particles: z from dx/2 up to (roof_z - dx/2)
    zs_base = np.arange(dx / 2, zmax, dx)
    particles = []
    skipped_cols = 0

    for ri, z_top in roof_z.items():
        x_col = xx[ri]
        y_col = yy[ri]
        # Fill from ground to just below roof
        valid_z = zs_base[zs_base < z_top - dx / 2]
        # Skip thin columns (pointy rooftops etc.)
        if len(valid_z) < min_z_layers:
            skipped_cols += 1
            continue
        for z in valid_z:
            particles.append([x_col, y_col, z])

    print(f"[*] Skipped {skipped_cols} thin columns (< {min_z_layers} layers)")

    if len(particles) == 0:
        print("[!] No interior particles generated.")
        return None

    particles = np.array(particles)
    n = particles.shape[0]
    print(f"[*] Generated {n} interior particles.")

    # Save
    print(f"[*] Saving to {output_filename} ...")
    with open(output_filename, 'w') as f:
        f.write(f"{n}\n")
        for p in particles:
            f.write(f"{p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")

    print("[*] Done.")
    return particles


if __name__ == "__main__":
    STL_FILE = "Leipzig_buildings_shifted.stl"
    DX = 2.5

    if len(sys.argv) >= 2:
        STL_FILE = sys.argv[1]
    if len(sys.argv) >= 3:
        DX = float(sys.argv[2])

    pts = generate_interior_particles(STL_FILE, DX)

    if pts is not None:
        print(f"\n=== SUMMARY ===")
        print(f"STL file   : {STL_FILE}")
        print(f"dx (pitch) : {DX} m")
        print(f"Total pts  : {len(pts)}")
        print(f"X range    : [{pts[:,0].min():.3f}, {pts[:,0].max():.3f}]")
        print(f"Y range    : [{pts[:,1].min():.3f}, {pts[:,1].max():.3f}]")
        print(f"Z range    : [{pts[:,2].min():.3f}, {pts[:,2].max():.3f}]")
