#!/usr/bin/env python3
"""
STL 없이 IBM 입자 위치 기반으로 내부 유체 입자를 마킹.

단순 박스 형상: IBM bounding box 기반 직접 판별
복잡 형상: IBM 입자 voxel BFS 방식 (fill_interior.py 로직 재사용)

Usage:
    python flag_interior_noSTL.py <particles.txt> [output.txt] [method: box|bfs] [margin]
"""
import sys, time
import numpy as np
from collections import deque


def flag_by_bounding_box(data, col, margin=0.6):
    """IBM bounding box 기반 내부 판별. 박스 건물에 적합."""
    ptypes = data[:, col['p_type']].astype(int)
    ibm = data[ptypes >= 1000]
    fluid_mask = ptypes == 1

    xlo, xhi = ibm[:, col['x']].min() + margin, ibm[:, col['x']].max() - margin
    ylo, yhi = ibm[:, col['y']].min() + margin, ibm[:, col['y']].max() - margin
    zlo, zhi = 0.0,                              ibm[:, col['z']].max() - margin

    print(f"[*] IBM bounds: x=[{xlo:.2f},{xhi:.2f}] y=[{ylo:.2f},{yhi:.2f}] z=[{zlo:.2f},{zhi:.2f}]")

    x, y, z = data[:, col['x']], data[:, col['y']], data[:, col['z']]
    inside = (fluid_mask
              & (x > xlo) & (x < xhi)
              & (y > ylo) & (y < yhi)
              & (z > zlo) & (z < zhi))
    return inside


def flag_by_bfs(data, col, dx=None, margin=0.3):
    """IBM 입자를 벽 복셀로 사용해 BFS flood fill."""
    ptypes = data[:, col['p_type']].astype(int)
    ibm = data[ptypes >= 1000]

    if dx is None:
        # IBM 입자 간격으로 dx 추정
        xs = np.sort(np.unique(np.round(ibm[:, col['x']], 3)))
        dx = float(np.median(np.diff(xs))) if len(xs) > 1 else 1.0
    print(f"[*] Estimated dx={dx:.3f}m")

    surf_xyz = ibm[:, [col['x'], col['y'], col['z']]]

    pad = 2
    origin = surf_xyz.min(axis=0) - pad * dx
    sz = np.ceil((surf_xyz.max(axis=0) - origin) / dx).astype(int) + pad + 1
    nx, ny, nz = sz

    # 벽 복셀 마킹
    si = np.round((surf_xyz[:, 0] - origin[0]) / dx).astype(int)
    sj = np.round((surf_xyz[:, 1] - origin[1]) / dx).astype(int)
    sk = np.round((surf_xyz[:, 2] - origin[2]) / dx).astype(int)
    valid = (si>=0)&(si<nx)&(sj>=0)&(sj<ny)&(sk>=0)&(sk<nz)
    grid = np.zeros((nx, ny, nz), dtype=np.uint8)
    grid[si[valid], sj[valid], sk[valid]] = 1

    # z=0 바닥 봉합
    k_floor = max(0, int(round((0.0 - origin[2]) / dx)) - 1)
    grid[:, :, k_floor] = 1
    print(f"[*] Grid: {nx}x{ny}x{nz}, wall voxels: {grid.sum()}, floor at k={k_floor}")

    # BFS from boundary
    OUTSIDE = 2
    queue = deque()
    def seed(i, j, k):
        if 0<=i<nx and 0<=j<ny and 0<=k<nz and grid[i,j,k]==0:
            grid[i,j,k]=OUTSIDE; queue.append((i,j,k))

    for j in range(ny):
        for k in range(nz): seed(0,j,k); seed(nx-1,j,k)
    for i in range(nx):
        for k in range(nz): seed(i,0,k); seed(i,ny-1,k)
    for i in range(nx):
        for j in range(ny): seed(i,j,0); seed(i,j,nz-1)

    nb = [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]
    while queue:
        ci,cj,ck = queue.popleft()
        for di,dj,dk in nb:
            ni,nj,nk = ci+di,cj+dj,ck+dk
            if 0<=ni<nx and 0<=nj<ny and 0<=nk<nz and grid[ni,nj,nk]==0:
                grid[ni,nj,nk]=OUTSIDE; queue.append((ni,nj,nk))

    print(f"[*] Interior voxels: {(grid==0).sum()}")

    # 유체 입자 판별
    fluid_mask = ptypes == 1
    fx = data[fluid_mask, col['x']]
    fy = data[fluid_mask, col['y']]
    fz = data[fluid_mask, col['z']]
    vi = np.round((fx - origin[0]) / dx).astype(int)
    vj = np.round((fy - origin[1]) / dx).astype(int)
    vk = np.floor((fz - origin[2]) / dx).astype(int)
    valid_f = (vi>=0)&(vi<nx)&(vj>=0)&(vj<ny)&(vk>=0)&(vk<nz)
    inside_f = np.zeros(fluid_mask.sum(), dtype=bool)
    idx = np.where(valid_f)[0]
    inside_f[idx] = (grid[vi[idx], vj[idx], vk[idx]] == 0)

    inside_all = np.zeros(len(data), dtype=bool)
    inside_all[np.where(fluid_mask)[0]] = inside_f
    return inside_all


def run(input_path, output_path, method='bfs', margin=0.6):
    t0 = time.time()
    print(f"[*] Reading: {input_path}")
    with open(input_path) as f:
        header = f.readline().strip()
        lines = f.readlines()

    col_indices = list(map(int, header.split()))
    col = {
        'x':      col_indices.index(1),
        'y':      col_indices.index(2),
        'z':      col_indices.index(3),
        'p_type': col_indices.index(8),
    }
    data = np.array([[float(v) for v in l.split('\t')] for l in lines])
    print(f"[*] {len(data)} particles loaded ({time.time()-t0:.1f}s)")

    if method == 'box':
        inside = flag_by_bounding_box(data, col, margin=margin)
    else:
        inside = flag_by_bfs(data, col, margin=margin)

    n = inside.sum()
    print(f"[*] Interior fluid particles: {n} / {(data[:,col['p_type']]==1).sum()}")

    with open(output_path, 'w') as f:
        f.write(header + '\n')
        for i, line in enumerate(lines):
            if inside[i]:
                parts = line.split('\t')
                parts[col['p_type']] = '99'
                f.write('\t'.join(parts))
            else:
                f.write(line)

    print(f"[*] Written: {output_path}  ({time.time()-t0:.1f}s)")


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python flag_interior_noSTL.py <particles.txt> [output.txt] [box|bfs] [margin]")
        sys.exit(1)
    inp    = sys.argv[1]
    out    = sys.argv[2] if len(sys.argv) >= 3 else inp
    method = sys.argv[3] if len(sys.argv) >= 4 else 'bfs'
    margin = float(sys.argv[4]) if len(sys.argv) >= 5 else 0.6
    run(inp, out, method, margin)
