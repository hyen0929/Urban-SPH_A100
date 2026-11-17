#!/usr/bin/env python3
"""
BFS flood-fill interior detection using STL-based voxelization.
Works on non-watertight meshes from OSM.

Usage:
    python flag_interior_bfs.py <particles.txt> <buildings.stl> [output.txt] [dx]
"""
import sys, time
import numpy as np
import trimesh
from collections import deque


def build_interior_mask(stl_path, dx=2.5, verbose=True):
    """
    STL 삼각형 기반 복셀화 → BFS flood fill → interior 마스크 반환.
    Returns: (grid, origin, dx)
      grid: uint8 array, 0=interior, 1=wall, 2=outside
      origin: world coord of voxel (0,0,0)
    """
    t0 = time.time()
    mesh = trimesh.load(stl_path)
    if verbose:
        print(f"[*] STL loaded: {len(mesh.faces)} faces, watertight={mesh.is_watertight}")
        print(f"[*] Bounds: {mesh.bounds[0]} ~ {mesh.bounds[1]}")

    # STL 삼각형 기반 복셀화
    vox = trimesh.voxel.creation.voxelize(mesh, pitch=dx)
    origin = vox.transform[:3, 3].copy()
    grid = vox.matrix.copy().astype(np.uint8)  # 1=wall, 0=unknown
    nx, ny, nz = grid.shape

    # z 방향으로 1층 확장 (바닥 봉합용 가상층)
    floor_layer = np.ones((nx, ny, 1), dtype=np.uint8)  # 전체 막음
    grid = np.concatenate([floor_layer, grid], axis=2)  # k=0이 z=-dx 가상층
    nz += 1
    origin[2] -= dx  # origin을 한 층 아래로

    if verbose:
        print(f"[*] Voxel grid: {nx}x{ny}x{nz} (z extended by 1), wall={grid.sum()}")

    # BFS seeds: 도메인 6면 경계
    OUTSIDE = 2
    queue = deque()

    def seed(i, j, k):
        if 0 <= i < nx and 0 <= j < ny and 0 <= k < nz and grid[i, j, k] == 0:
            grid[i, j, k] = OUTSIDE
            queue.append((i, j, k))

    for j in range(ny):
        for k in range(nz): seed(0, j, k);   seed(nx-1, j, k)
    for i in range(nx):
        for k in range(nz): seed(i, 0, k);   seed(i, ny-1, k)
    for i in range(nx):
        for j in range(ny): seed(i, j, 0);   seed(i, j, nz-1)

    if verbose:
        print(f"[*] BFS seeds: {len(queue)}, starting flood fill...")

    t1 = time.time()
    nb = [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]
    while queue:
        ci, cj, ck = queue.popleft()
        for di, dj, dk in nb:
            ni, nj, nk = ci+di, cj+dj, ck+dk
            if 0 <= ni < nx and 0 <= nj < ny and 0 <= nk < nz and grid[ni, nj, nk] == 0:
                grid[ni, nj, nk] = OUTSIDE
                queue.append((ni, nj, nk))

    n_interior = (grid == 0).sum()
    if verbose:
        print(f"[*] BFS done in {time.time()-t1:.1f}s | interior voxels: {n_interior}")
        print(f"[*] Total: {time.time()-t0:.1f}s")

    return grid, origin, dx


def is_interior_voxel(grid, origin, dx, points):
    """
    점들이 interior 복셀에 속하는지 반환 (bool array).
    Uses floor division (not round) for z to avoid floor-layer collision.
    """
    nx, ny, nz = grid.shape
    # XY: nearest neighbor (round), Z: floor division
    vi = np.round((points[:, 0] - origin[0]) / dx).astype(int)
    vj = np.round((points[:, 1] - origin[1]) / dx).astype(int)
    vk = np.floor((points[:, 2] - origin[2]) / dx).astype(int)

    valid = (vi >= 0) & (vi < nx) & (vj >= 0) & (vj < ny) & (vk >= 0) & (vk < nz)
    result = np.zeros(len(points), dtype=bool)
    idx = np.where(valid)[0]
    result[idx] = (grid[vi[idx], vj[idx], vk[idx]] == 0)
    return result


def flag_interior(input_path, stl_path, output_path, dx=2.5):
    grid, origin, dx = build_interior_mask(stl_path, dx)

    print(f"[*] Reading: {input_path}")
    with open(input_path, 'r') as f:
        header = f.readline().strip()
        lines = f.readlines()

    col_indices = list(map(int, header.split()))
    ptype_col = col_indices.index(8)
    x_col = col_indices.index(1)
    y_col = col_indices.index(2)
    z_col = col_indices.index(3)

    # 유체 입자만 추출
    fluid_idx, fluid_pts = [], []
    for i, line in enumerate(lines):
        f = line.split('\t')
        if int(float(f[ptype_col])) == 1:
            fluid_idx.append(i)
            fluid_pts.append([float(f[x_col]), float(f[y_col]), float(f[z_col])])

    fluid_pts = np.array(fluid_pts)
    print(f"[*] Fluid particles: {len(fluid_pts)}")

    inside = is_interior_voxel(grid, origin, dx, fluid_pts)
    n_flagged = inside.sum()
    print(f"[*] Interior particles: {n_flagged} / {len(fluid_pts)}")

    flag_set = {fluid_idx[i] for i in np.where(inside)[0]}

    with open(output_path, 'w') as f:
        f.write(header + '\n')
        for i, line in enumerate(lines):
            if i in flag_set:
                parts = line.split('\t')
                parts[ptype_col] = '99'
                f.write('\t'.join(parts))
            else:
                f.write(line)

    print(f"[*] Written: {output_path}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python flag_interior_bfs.py <particles.txt> <buildings.stl> [output.txt] [dx]")
        sys.exit(1)

    inp   = sys.argv[1]
    stl   = sys.argv[2]
    out   = sys.argv[3] if len(sys.argv) >= 4 else inp
    dx    = float(sys.argv[4]) if len(sys.argv) >= 5 else 2.5
    flag_interior(inp, stl, out, dx)
