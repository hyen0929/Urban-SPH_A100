#!/usr/bin/env python3
"""Detailed analysis: building corner blowup + IBM filt quality."""
import numpy as np
import re

DTYPE_MAP = {
    "float": (">f4", 4), "double": (">f8", 8),
    "int": (">i4", 4), "unsigned_int": (">u4", 4),
}

def parse_vtk(filepath):
    with open(filepath, "rb") as f:
        content = f.read()
    m = re.search(rb"POINTS\s+(\d+)\s+(\w+)\n", content)
    n_pts = int(m.group(1))
    dt, sz = DTYPE_MAP[m.group(2).decode("ascii")]
    pts = np.frombuffer(content, dtype=dt, count=n_pts*3, offset=m.end()).reshape(n_pts, 3).astype(np.float64)
    fields = {}
    m_field = re.search(rb"FIELD FieldData\s+(\d+)\n", content)
    if not m_field:
        return pts, fields
    offset = m_field.end()
    for _ in range(int(m_field.group(1))):
        newline = content.index(b"\n", offset)
        header = content[offset:newline].decode("ascii").split()
        fname, ncomp, ntuple, dtype_s = header[0], int(header[1]), int(header[2]), header[3]
        offset = newline + 1
        dt2, itemsize = DTYPE_MAP.get(dtype_s, (">f4", 4))
        arr = np.frombuffer(content, dtype=dt2, count=ncomp*ntuple, offset=offset).astype(np.float64)
        if ncomp > 1: arr = arr.reshape(ntuple, ncomp)
        fields[fname] = arr
        offset += ncomp * ntuple * itemsize
    return pts, fields

fpath = "Leipzig_small_point/fluid_condition_3000stp.vtk"
pts, fields = parse_vtk(fpath)
vel = fields["velocity"]
ptype = fields["p_type"].astype(int)
speed = np.sqrt(vel[:,0]**2 + vel[:,1]**2 + vel[:,2]**2)

# IBM geometry
ibm = ptype >= 1000
ibm_pts = pts[ibm]
IBM_X = (ibm_pts[:,0].min(), ibm_pts[:,0].max())  # ~116, ~347
IBM_Y = (ibm_pts[:,1].min(), ibm_pts[:,1].max())  # ~0, ~133
IBM_Z = (ibm_pts[:,2].min(), ibm_pts[:,2].max())  # ~0, ~29

print(f"Building: x=[{IBM_X[0]:.1f},{IBM_X[1]:.1f}], y=[{IBM_Y[0]:.1f},{IBM_Y[1]:.1f}], z=[{IBM_Z[0]:.1f},{IBM_Z[1]:.1f}]")

fluid = ptype == 1
fp = pts[fluid]
fv = vel[fluid]
fs = speed[fluid]

# ========================================
# Issue 1: z=1.2 layer 분석 (비물리적 고속도)
# ========================================
print(f"\n{'='*80}")
print("Issue 1: z=1.2 층 (최하단 유체) 비물리적 고속도")
print(f"{'='*80}")

z_low = fp[:, 2] < 2.0
fp_low = fp[z_low]
fv_low = fv[z_low]
fs_low = fs[z_low]

# 건물 footprint 안 vs 밖
inside_footprint = (fp_low[:,0] > IBM_X[0]) & (fp_low[:,0] < IBM_X[1]) & \
                   (fp_low[:,1] > IBM_Y[0]) & (fp_low[:,1] < IBM_Y[1])
outside_footprint = ~inside_footprint

print(f"\nz<2.0 유체 입자: N={z_low.sum()}")
print(f"  건물 footprint 안: N={inside_footprint.sum()}, |u|_mean={fs_low[inside_footprint].mean():.2f}, max={fs_low[inside_footprint].max():.2f}")
print(f"  건물 footprint 밖: N={outside_footprint.sum()}, |u|_mean={fs_low[outside_footprint].mean():.2f}, max={fs_low[outside_footprint].max():.2f}")

# 건물 edge 근처 세분화 (footprint 안, edge 근처 vs 중심)
if inside_footprint.any():
    edge_dist = 10.0  # edge에서의 거리
    fp_in = fp_low[inside_footprint]
    fv_in = fv_low[inside_footprint]
    fs_in = fs_low[inside_footprint]

    near_x_edge = (fp_in[:,0] < IBM_X[0]+edge_dist) | (fp_in[:,0] > IBM_X[1]-edge_dist)
    near_y_edge = (fp_in[:,1] < IBM_Y[0]+edge_dist) | (fp_in[:,1] > IBM_Y[1]-edge_dist)
    near_edge = near_x_edge | near_y_edge
    core = ~near_edge

    print(f"\n  건물 안쪽 z<2.0:")
    print(f"    Edge 근처 (<{edge_dist}m): N={near_edge.sum()}, |u|_mean={fs_in[near_edge].mean():.2f}, max={fs_in[near_edge].max():.2f}")
    if core.any():
        print(f"    Core (중앙):    N={core.sum()}, |u|_mean={fs_in[core].mean():.2f}, max={fs_in[core].max():.2f}")

# 건물 밖 edge 근처 (flow separation zone)
if outside_footprint.any():
    fp_out = fp_low[outside_footprint]
    fv_out = fv_low[outside_footprint]
    fs_out = fs_low[outside_footprint]

    # y 방향 (측면)
    near_y_out = (fp_out[:,1] > IBM_Y[1]) & (fp_out[:,1] < IBM_Y[1]+30) & \
                 (fp_out[:,0] > IBM_X[0]) & (fp_out[:,0] < IBM_X[1])
    # x 방향 downstream
    downstream = (fp_out[:,0] > IBM_X[1]) & (fp_out[:,0] < IBM_X[1]+100) & \
                 (fp_out[:,1] > IBM_Y[0]) & (fp_out[:,1] < IBM_Y[1])
    # corners
    corner = ((fp_out[:,0] > IBM_X[1]) | (fp_out[:,1] > IBM_Y[1])) & \
             (fp_out[:,0] > IBM_X[0]-10) & (fp_out[:,0] < IBM_X[1]+30) & \
             (fp_out[:,1] > IBM_Y[0]-10) & (fp_out[:,1] < IBM_Y[1]+30)

    print(f"\n  건물 밖 z<2.0 (flow separation):")
    if near_y_out.any():
        print(f"    측면(y>{IBM_Y[1]:.0f}): N={near_y_out.sum()}, |u|_mean={fs_out[near_y_out].mean():.2f}, max={fs_out[near_y_out].max():.2f}")
    if downstream.any():
        print(f"    후류(x>{IBM_X[1]:.0f}): N={downstream.sum()}, |u|_mean={fs_out[downstream].mean():.2f}, max={fs_out[downstream].max():.2f}")

# 고속도 입자의 시간 발전 위치 추적 (y, x edge 인접성)
print(f"\n  |u| > 30 m/s 입자 위치 분석:")
extreme = fs_low > 30
if extreme.any():
    ep = fp_low[extreme]
    ev = fv_low[extreme]
    es = fs_low[extreme]
    print(f"    N={extreme.sum()}")
    print(f"    x: [{ep[:,0].min():.1f}, {ep[:,0].max():.1f}] (IBM: {IBM_X[0]:.1f}~{IBM_X[1]:.1f})")
    print(f"    y: [{ep[:,1].min():.1f}, {ep[:,1].max():.1f}] (IBM: {IBM_Y[0]:.1f}~{IBM_Y[1]:.1f})")
    print(f"    z: [{ep[:,2].min():.2f}, {ep[:,2].max():.2f}]")

    # IBM footprint 기준 위치
    in_x = (ep[:,0] >= IBM_X[0]) & (ep[:,0] <= IBM_X[1])
    in_y = (ep[:,1] >= IBM_Y[0]) & (ep[:,1] <= IBM_Y[1])
    print(f"    x in IBM: {in_x.sum()}/{len(ep)}")
    print(f"    y in IBM: {in_y.sum()}/{len(ep)}")
    print(f"    both in IBM footprint: {(in_x & in_y).sum()}/{len(ep)}")
    print(f"    y > IBM_y_max ({IBM_Y[1]:.0f}): {(ep[:,1]>IBM_Y[1]).sum()}/{len(ep)}")

# ========================================
# Issue 2: IBM 경계 아래쪽 (유체와 먼 쪽) 분석
# ========================================
print(f"\n{'='*80}")
print("Issue 2: IBM 입자 위치별 주변 유체 속도")
print(f"{'='*80}")

# IBM 입자를 z층별, surface/interior별로 나누고 가장 가까운 유체 속도 확인
ibm_vel = vel[ibm]
ibm_sp = speed[ibm]

# IBM 입자의 surface 판별: x, y, z edge에서 가까운지
ibm_edge_dist = 5.0  # 입자간격 기준
ibm_near_top = ibm_pts[:,2] > IBM_Z[1] - ibm_edge_dist
ibm_near_bot = ibm_pts[:,2] < IBM_Z[0] + ibm_edge_dist
ibm_near_xmin = ibm_pts[:,0] < IBM_X[0] + ibm_edge_dist
ibm_near_xmax = ibm_pts[:,0] > IBM_X[1] - ibm_edge_dist
ibm_near_ymin = ibm_pts[:,1] < IBM_Y[0] + ibm_edge_dist
ibm_near_ymax = ibm_pts[:,1] > IBM_Y[1] - ibm_edge_dist
ibm_surface = ibm_near_top | ibm_near_bot | ibm_near_xmin | ibm_near_xmax | ibm_near_ymin | ibm_near_ymax
ibm_interior = ~ibm_surface

print(f"\nIBM 입자 분류:")
print(f"  Surface: N={ibm_surface.sum()}")
print(f"  Interior: N={ibm_interior.sum()}")

# IBM 입자 z별 프로필 - 가장 가까운 유체 입자 속도
# (full nearest neighbor search는 비싸므로, z-layer별 대표 분석)
from scipy.spatial import cKDTree

# z layer별 IBM 입자 주변 유체 속도 (건물 안 유체 + 건물 밖 유체)
print(f"\nIBM z-layer별 가장 가까운 유체 입자 속도 (10개 샘플씩):")
z_layers = np.arange(IBM_Z[0], IBM_Z[1]+1, 5.0)

# fluid kdtree
fluid_tree = cKDTree(fp)

for z_lo in z_layers:
    z_hi = z_lo + 5.0
    layer_mask = (ibm_pts[:,2] >= z_lo) & (ibm_pts[:,2] < z_hi)
    if layer_mask.sum() == 0:
        continue

    # surface vs interior in this layer
    layer_surf = layer_mask & ibm_surface
    layer_int = layer_mask & ibm_interior

    # Sample 10 interior IBM particles
    for label, mask in [("surface", layer_surf), ("interior", layer_int)]:
        if mask.sum() == 0:
            continue
        sample_idx = np.where(mask)[0][:min(20, mask.sum())]
        dists_all = []
        speeds_all = []
        for idx in sample_idx:
            d, j = fluid_tree.query(ibm_pts[idx], k=5)
            dists_all.append(d[0])
            speeds_all.append(fs[j].mean())

        avg_dist = np.mean(dists_all)
        avg_speed = np.mean(speeds_all)
        max_speed = np.max(speeds_all)
        print(f"  z=[{z_lo:.0f},{z_hi:.0f}] {label:>8}: N={mask.sum():>5}, "
              f"nearest_fluid_dist={avg_dist:.1f}m, "
              f"fluid_|u|_mean={avg_speed:.2f}, max={max_speed:.2f}")

# ========================================
# p_type=0 (buffer) 입자 분석
# ========================================
print(f"\n{'='*80}")
print("p_type=0 (buffer) 입자 분석")
print(f"{'='*80}")
buf = ptype == 0
bp = pts[buf]
bv = vel[buf]
bs = speed[buf]
print(f"N={buf.sum()}, |u|_max={bs.max():.1f}")
print(f"x: [{bp[:,0].min():.1f}, {bp[:,0].max():.1f}]")
print(f"y: [{bp[:,1].min():.1f}, {bp[:,1].max():.1f}]")
print(f"z: [{bp[:,2].min():.1f}, {bp[:,2].max():.1f}]")

# buffer 입자 중 IBM 근처
buf_near_ibm = (bp[:,0] > IBM_X[0]-20) & (bp[:,0] < IBM_X[1]+20) & \
               (bp[:,1] > IBM_Y[0]-20) & (bp[:,1] < IBM_Y[1]+20)
if buf_near_ibm.any():
    print(f"IBM 근처 buffer: N={buf_near_ibm.sum()}, |u|_max={bs[buf_near_ibm].max():.1f}")
