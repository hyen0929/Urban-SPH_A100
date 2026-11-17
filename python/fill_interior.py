#!/usr/bin/env python3
"""
Leipzig_surface.txt 의 건물 표면 내부를
(0,0,0) 기준 2.5m 정규 격자점으로 채우는 스크립트.

Approach: voxelize surface → dilate → BFS flood-fill from outside → interior extraction
(scipy 없이 순수 numpy + collections.deque 사용)
"""
import numpy as np
from collections import deque
import time

# ============================================================
# Parameters
# ============================================================
dx = 2.5           # particle spacing (m)
surface_file = "Leipzig_surface.txt"
output_file = "Leipzig_interior_2p5m.txt"

# ============================================================
# 1. Read surface points
# ============================================================
t0 = time.time()
print("Reading surface points...")
with open(surface_file, "r") as f:
    n_pts = int(f.readline().strip())
    pts = np.loadtxt(f, max_rows=n_pts)
print(f"  {n_pts} points loaded ({time.time()-t0:.1f}s)")
print(f"  X: [{pts[:,0].min():.1f}, {pts[:,0].max():.1f}]")
print(f"  Y: [{pts[:,1].min():.1f}, {pts[:,1].max():.1f}]")
print(f"  Z: [{pts[:,2].min():.1f}, {pts[:,2].max():.1f}]")

# ============================================================
# 2. Build voxel grid aligned to (0,0,0) with spacing dx
# ============================================================
pad = 2  # padding in voxel units

i_min = int(np.floor(pts[:,0].min() / dx)) - pad
i_max = int(np.ceil(pts[:,0].max() / dx)) + pad
j_min = int(np.floor(pts[:,1].min() / dx)) - pad
j_max = int(np.ceil(pts[:,1].max() / dx)) + pad
k_min = int(np.floor(pts[:,2].min() / dx)) - pad
k_max = int(np.ceil(pts[:,2].max() / dx)) + pad

# Ensure k_min covers below ground for sealing
if k_min > -2:
    k_min = -2

nx = i_max - i_min + 1
ny = j_max - j_min + 1
nz = k_max - k_min + 1
print(f"\nVoxel grid: {nx} x {ny} x {nz} = {nx*ny*nz:,} voxels")

# ============================================================
# 3. Mark surface voxels
# ============================================================
print("Marking surface voxels...")
# 0 = unknown, 1 = blocked (surface/shell), 2 = outside
grid = np.zeros((nx, ny, nz), dtype=np.uint8)

si = np.round(pts[:,0] / dx).astype(int) - i_min
sj = np.round(pts[:,1] / dx).astype(int) - j_min
sk = np.round(pts[:,2] / dx).astype(int) - k_min

valid = (si >= 0) & (si < nx) & (sj >= 0) & (sj < ny) & (sk >= 0) & (sk < nz)
si, sj, sk = si[valid], sj[valid], sk[valid]
grid[si, sj, sk] = 1

n_surf_voxels = (grid == 1).sum()
print(f"  Surface voxels: {n_surf_voxels:,}")

# ============================================================
# 4. Seal building bottoms at z=0
# ============================================================
print("Sealing building bottoms at z=0...")
k_z0 = round(0.0 / dx) - k_min  # voxel index for z=0

surface_mask = (grid == 1)
has_surface = surface_mask.any(axis=2)  # (nx, ny)

for ii in range(nx):
    for jj in range(ny):
        if has_surface[ii, jj]:
            col = surface_mask[ii, jj, :]
            k_indices = np.where(col)[0]
            k_lowest = k_indices.min()
            if k_lowest > 0 and k_z0 < k_lowest:
                grid[ii, jj, k_z0:k_lowest] = 1

n_after_seal = (grid == 1).sum()
print(f"  Voxels after sealing: {n_after_seal:,} (+{n_after_seal - n_surf_voxels:,})")

# ============================================================
# 5. Dilate surface shell (6-connectivity, 1 iteration)
#    Using numpy roll for each of the 6 face-neighbors
# ============================================================
print("Dilating surface shell (6-connectivity)...")
blocked = (grid == 1)
dilated = blocked.copy()

# 6 face-neighbor shifts
dilated |= np.roll(blocked, 1, axis=0)
dilated |= np.roll(blocked, -1, axis=0)
dilated |= np.roll(blocked, 1, axis=1)
dilated |= np.roll(blocked, -1, axis=1)
dilated |= np.roll(blocked, 1, axis=2)
dilated |= np.roll(blocked, -1, axis=2)

# Fix roll wrap-around edges
dilated[0, :, :] |= blocked[0, :, :]
dilated[-1, :, :] |= blocked[-1, :, :]
dilated[:, 0, :] |= blocked[:, 0, :]
dilated[:, -1, :] |= blocked[:, -1, :]
dilated[:, :, 0] |= blocked[:, :, 0]
dilated[:, :, -1] |= blocked[:, :, -1]

grid[dilated] = 1
n_shell = (grid == 1).sum()
print(f"  Shell voxels: {n_shell:,}")

# ============================================================
# 6. BFS flood fill from outside
#    Start from all boundary voxels that are NOT blocked
# ============================================================
print("BFS flood fill from outside...")
t1 = time.time()

# Seed: all boundary faces of the grid that are not blocked
queue = deque()

def seed_face(idx_array):
    """Add non-blocked boundary voxels to queue."""
    for ii, jj, kk in idx_array:
        if grid[ii, jj, kk] == 0:
            grid[ii, jj, kk] = 2
            queue.append((ii, jj, kk))

# Seed all 6 boundary faces
faces = [
    (0, slice(None), slice(None)),       # i=0
    (nx-1, slice(None), slice(None)),    # i=nx-1
    (slice(None), 0, slice(None)),       # j=0
    (slice(None), ny-1, slice(None)),    # j=ny-1
    (slice(None), slice(None), 0),       # k=0
    (slice(None), slice(None), nz-1),    # k=nz-1
]
for face in faces:
    face_grid = grid[face]
    idxs = np.argwhere(face_grid == 0)
    for idx in idxs:
        # Map 2D face index back to 3D
        coord = [0, 0, 0]
        dim = 0
        for d in range(3):
            if isinstance(face[d], int):
                coord[d] = face[d]
            else:
                coord[d] = idx[dim]
                dim += 1
        ci, cj, ck = coord
        if grid[ci, cj, ck] == 0:
            grid[ci, cj, ck] = 2
            queue.append((ci, cj, ck))

print(f"  Seed voxels: {len(queue):,}")

# BFS with 6-connectivity
n_visited = len(queue)
while queue:
    ci, cj, ck = queue.popleft()
    for di, dj, dk in [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]:
        ni, nj, nk = ci+di, cj+dj, ck+dk
        if 0 <= ni < nx and 0 <= nj < ny and 0 <= nk < nz:
            if grid[ni, nj, nk] == 0:
                grid[ni, nj, nk] = 2
                queue.append((ni, nj, nk))
                n_visited += 1

print(f"  Outside voxels: {n_visited:,} ({time.time()-t1:.1f}s)")

# ============================================================
# 7. Interior = voxels still marked as 0 (not blocked, not outside)
# ============================================================
interior = (grid == 0)
n_interior = interior.sum()
print(f"  Interior voxels: {n_interior:,}")

# ============================================================
# 8. Extract coordinates, filter z >= 0
# ============================================================
print("Extracting coordinates...")
ki, kj, kk = np.where(interior)
x_coords = (ki + i_min) * dx
y_coords = (kj + j_min) * dx
z_coords = (kk + k_min) * dx

above = z_coords >= 0.0
x_coords = x_coords[above]
y_coords = y_coords[above]
z_coords = z_coords[above]

n_final = len(x_coords)
print(f"\n  Interior particles (z>=0): {n_final:,}")

if n_final > 0:
    print(f"  X: [{x_coords.min():.1f}, {x_coords.max():.1f}]")
    print(f"  Y: [{y_coords.min():.1f}, {y_coords.max():.1f}]")
    print(f"  Z: [{z_coords.min():.1f}, {z_coords.max():.1f}]")

# ============================================================
# 9. Save output
# ============================================================
print(f"\nSaving {output_file}...")
with open(output_file, "w") as f:
    f.write(f"{n_final}\n")
    for xi, yi, zi in zip(x_coords, y_coords, z_coords):
        f.write(f"{xi:.6f} {yi:.6f} {zi:.6f}\n")

print(f"Done! ({time.time()-t0:.1f}s total)")
