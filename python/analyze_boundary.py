#!/usr/bin/env python3
"""Analyze boundary particles (p_type=-1, 0) and z<0 region."""
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
    if not m_field: return pts, fields
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

pts, fields = parse_vtk("Leipzig_small_point/fluid_condition_3000stp.vtk")
vel = fields["velocity"]
ptype = fields["p_type"].astype(int)
speed = np.sqrt(vel[:,0]**2 + vel[:,1]**2 + vel[:,2]**2)

# ========================================
# 전체 p_type별 z-분포
# ========================================
print("="*80)
print("전체 p_type별 z-분포")
print("="*80)
for pt in sorted(np.unique(ptype)):
    m = ptype == pt
    z = pts[m, 2]
    sp = speed[m]
    print(f"  p_type={pt:>5}: N={m.sum():>8}, z=[{z.min():>8.2f}, {z.max():>8.2f}], "
          f"|u|: mean={sp.mean():.3f}, max={sp.max():.3f}")

# ========================================
# z<0 영역 전체 분석
# ========================================
print(f"\n{'='*80}")
print("z<0 영역 상세 분석")
print(f"{'='*80}")

below_ground = pts[:, 2] < 0
bp = pts[below_ground]
bv = vel[below_ground]
bs = speed[below_ground]
bpt = ptype[below_ground]

print(f"z<0 전체: N={below_ground.sum()}")
for pt in sorted(np.unique(bpt)):
    m = bpt == pt
    print(f"  p_type={pt}: N={m.sum()}, z=[{bp[m,2].min():.2f},{bp[m,2].max():.2f}], "
          f"|u|_mean={bs[m].mean():.3f}, max={bs[m].max():.3f}")

# ========================================
# Boundary(p_type=-1) z-layer별 속도 프로필
# ========================================
print(f"\n{'='*80}")
print("p_type=-1 (boundary) z-layer별 속도")
print(f"{'='*80}")

bnd = ptype == -1
bnd_pts = pts[bnd]
bnd_vel = vel[bnd]
bnd_sp = speed[bnd]

z_unique = np.unique(np.round(bnd_pts[:,2], 2))
print(f"  z 고유값: {z_unique}")

print(f"\n  {'z':>8} | {'N':>6} | {'|u|_mean':>10} | {'|u|_max':>10} | {'ux_mean':>10} | {'uy_mean':>10} | {'uz_mean':>10}")
for z_val in sorted(z_unique):
    zm = np.abs(bnd_pts[:,2] - z_val) < 0.1
    if zm.sum() > 0:
        sp_z = bnd_sp[zm]
        v_z = bnd_vel[zm]
        print(f"  {z_val:>8.2f} | {zm.sum():>6} | {sp_z.mean():>10.3f} | {sp_z.max():>10.3f} | "
              f"{v_z[:,0].mean():>10.3f} | {v_z[:,1].mean():>10.3f} | {v_z[:,2].mean():>10.3f}")

# ========================================
# z=-8.75 상세 분석
# ========================================
print(f"\n{'='*80}")
print("z=-8.75 층 상세 분석")
print(f"{'='*80}")

z_target = -8.75
z_tol = 0.5
layer = np.abs(pts[:,2] - z_target) < z_tol

for pt in sorted(np.unique(ptype[layer])):
    m = layer & (ptype == pt)
    if m.sum() == 0: continue
    sp_m = speed[m]
    v_m = vel[m]
    p_m = pts[m]
    print(f"\n  p_type={pt}: N={m.sum()}")
    print(f"    |u|: mean={sp_m.mean():.3f}, max={sp_m.max():.3f}")
    print(f"    ux: mean={v_m[:,0].mean():.3f}, max={v_m[:,0].max():.3f}, min={v_m[:,0].min():.3f}")
    print(f"    x: [{p_m[:,0].min():.1f}, {p_m[:,0].max():.1f}]")
    print(f"    y: [{p_m[:,1].min():.1f}, {p_m[:,1].max():.1f}]")

    # IBM 건물 영역과 겹치는지
    ibm = ptype >= 1000
    IBM_X = (pts[ibm,0].min(), pts[ibm,0].max())
    IBM_Y = (pts[ibm,1].min(), pts[ibm,1].max())
    in_ibm_xy = (p_m[:,0] > IBM_X[0]) & (p_m[:,0] < IBM_X[1]) & \
                (p_m[:,1] > IBM_Y[0]) & (p_m[:,1] < IBM_Y[1])
    out_ibm_xy = ~in_ibm_xy
    if in_ibm_xy.any():
        print(f"    건물 footprint 안: N={in_ibm_xy.sum()}, |u|_mean={sp_m[in_ibm_xy].mean():.3f}, max={sp_m[in_ibm_xy].max():.3f}")
    if out_ibm_xy.any():
        print(f"    건물 footprint 밖: N={out_ibm_xy.sum()}, |u|_mean={sp_m[out_ibm_xy].mean():.3f}, max={sp_m[out_ibm_xy].max():.3f}")

    # 고속도 입자 위치
    if sp_m.max() > 15:
        high = sp_m > sp_m.max() * 0.5
        hp = p_m[high]
        hv = v_m[high]
        hs = sp_m[high]
        print(f"\n    |u| > {sp_m.max()*0.5:.1f} 입자 (N={high.sum()}):")
        print(f"      x: [{hp[:,0].min():.1f}, {hp[:,0].max():.1f}]")
        print(f"      y: [{hp[:,1].min():.1f}, {hp[:,1].max():.1f}]")
        # Top 5
        top5 = np.argsort(hs)[-min(5,len(hs)):][::-1]
        for idx in top5:
            print(f"      |u|={hs[idx]:.2f}, ux={hv[idx,0]:.2f}, uy={hv[idx,1]:.2f}, uz={hv[idx,2]:.2f}, "
                  f"pos=({hp[idx,0]:.1f},{hp[idx,1]:.1f},{hp[idx,2]:.2f})")

# ========================================
# boundary 입자 velocity 설정 확인: BC 커널에서 어떻게 처리되는지 추적
# ========================================
print(f"\n{'='*80}")
print("Boundary vs fluid 속도 비교 (z 방향)")
print(f"{'='*80}")

# z-layer별 p_type=1과 p_type=-1 비교
z_vals_all = np.arange(-10, 15, 2.5)
print(f"  {'z_center':>8} | {'fluid N':>8} | {'fluid |u|':>10} | {'bnd N':>8} | {'bnd |u|':>10} | {'diff':>8}")
for i in range(len(z_vals_all)-1):
    zlo, zhi = z_vals_all[i], z_vals_all[i+1]
    zc = (zlo+zhi)/2

    fm = (ptype == 1) & (pts[:,2] >= zlo) & (pts[:,2] < zhi)
    bm = (ptype == -1) & (pts[:,2] >= zlo) & (pts[:,2] < zhi)

    if fm.sum() > 0 or bm.sum() > 0:
        f_sp = speed[fm].mean() if fm.sum() > 0 else 0
        b_sp = speed[bm].mean() if bm.sum() > 0 else 0
        f_n = fm.sum()
        b_n = bm.sum()
        diff = b_sp - f_sp if fm.sum() > 0 and bm.sum() > 0 else 0
        print(f"  {zc:>8.2f} | {f_n:>8} | {f_sp:>10.3f} | {b_n:>8} | {b_sp:>10.3f} | {diff:>+8.3f}")

# ========================================
# BC kernel 확인: boundary 속도가 어떻게 결정되는지
# ========================================
print(f"\n{'='*80}")
print("z=-8.75 boundary 입자 - z=1.25 유체 입자 속도 관계")
print(f"{'='*80}")

# z=-8.75 boundary 입자의 x,y 위치와 대응하는 z=1.25 유체 입자 속도 비교
bnd_layer = (ptype == -1) & (np.abs(pts[:,2] - (-8.75)) < 0.5)
fluid_layer = (ptype == 1) & (np.abs(pts[:,2] - 1.25) < 0.5)

bp_l = pts[bnd_layer]
bv_l = vel[bnd_layer]
bs_l = speed[bnd_layer]

fp_l = pts[fluid_layer]
fv_l = vel[fluid_layer]
fs_l = speed[fluid_layer]

# x-bin별 비교
x_bins = np.arange(-250, 500, 50)
print(f"  {'x_range':>15} | {'bnd_N':>6} | {'bnd_ux':>10} | {'fluid_N':>6} | {'fluid_ux':>10}")
for i in range(len(x_bins)-1):
    xlo, xhi = x_bins[i], x_bins[i+1]
    bm = (bp_l[:,0] >= xlo) & (bp_l[:,0] < xhi)
    fm = (fp_l[:,0] >= xlo) & (fp_l[:,0] < xhi)
    if bm.sum() > 0:
        print(f"  [{xlo:>6.0f},{xhi:>6.0f}) | {bm.sum():>6} | {bv_l[bm,0].mean():>10.3f} | "
              f"{fm.sum():>6} | {fv_l[fm,0].mean():>10.3f}" if fm.sum() > 0 else
              f"  [{xlo:>6.0f},{xhi:>6.0f}) | {bm.sum():>6} | {bv_l[bm,0].mean():>10.3f} | {'---':>6} | {'---':>10}")
