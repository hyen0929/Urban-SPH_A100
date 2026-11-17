#!/usr/bin/env python3
"""IBM interior vs surface: nearest fluid velocity (no scipy)."""
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

ibm = ptype >= 1000
fluid = ptype == 1
ibm_pts = pts[ibm]
fp = pts[fluid]
fv = vel[fluid]
fs = speed[fluid]

IBM_X = (ibm_pts[:,0].min(), ibm_pts[:,0].max())
IBM_Y = (ibm_pts[:,1].min(), ibm_pts[:,1].max())
IBM_Z = (ibm_pts[:,2].min(), ibm_pts[:,2].max())

# ========================================
# IBM 입자 z-layer별, surface/interior별 분석
# IBM 입자에서 가까운 유체 속도가 어떤지 확인
# ========================================
print("IBM 입자 위치별 가장 가까운 유체 분석")
print(f"Building: x=[{IBM_X[0]:.1f},{IBM_X[1]:.1f}], y=[{IBM_Y[0]:.1f},{IBM_Y[1]:.1f}], z=[{IBM_Z[0]:.1f},{IBM_Z[1]:.1f}]")

edge_dist = 5.0
ibm_surface = (ibm_pts[:,2] > IBM_Z[1]-edge_dist) | (ibm_pts[:,2] < IBM_Z[0]+edge_dist) | \
              (ibm_pts[:,0] < IBM_X[0]+edge_dist) | (ibm_pts[:,0] > IBM_X[1]-edge_dist) | \
              (ibm_pts[:,1] < IBM_Y[0]+edge_dist) | (ibm_pts[:,1] > IBM_Y[1]-edge_dist)

# z-layer별 분석
z_layers = np.arange(0, 30, 5.0)
print(f"\n{'z_range':>12} | {'type':>8} | {'N_ibm':>6} | {'near_fluid_dist':>16} | {'fluid_|u|_mean':>16} | {'fluid_|u|_max':>16}")
print("-"*100)

for z_lo in z_layers:
    z_hi = z_lo + 5.0
    layer = (ibm_pts[:,2] >= z_lo) & (ibm_pts[:,2] < z_hi)

    for label, type_mask in [("surface", ibm_surface & layer), ("interior", ~ibm_surface & layer)]:
        if type_mask.sum() == 0:
            continue
        # Sample 10 IBM particles
        idx_all = np.where(type_mask)[0]
        np.random.seed(42)
        sample = idx_all[np.random.choice(len(idx_all), min(15, len(idx_all)), replace=False)]

        dists = []
        near_speeds = []
        for si in sample:
            # Brute force nearest 5 fluid particles (only check within z+-30, x+-30, y+-30)
            dx = fp[:,0] - ibm_pts[si,0]
            dy = fp[:,1] - ibm_pts[si,1]
            dz = fp[:,2] - ibm_pts[si,2]
            r2 = dx*dx + dy*dy + dz*dz
            # Limit search
            near = r2 < 900  # 30m radius
            if near.sum() == 0:
                dists.append(999)
                near_speeds.append(0)
                continue
            r2_near = r2[near]
            fs_near = fs[near]
            top5 = np.argsort(r2_near)[:5]
            dists.append(np.sqrt(r2_near[top5[0]]))
            near_speeds.append(fs_near[top5].mean())

        print(f"  [{z_lo:>4.0f},{z_hi:>4.0f}) | {label:>8} | {type_mask.sum():>6} | "
              f"{np.mean(dists):>16.2f} | {np.mean(near_speeds):>16.2f} | {np.max(near_speeds):>16.2f}")

# ========================================
# Issue 2 구체 분석: 건물 내부 유체 입자 (footprint 안, z < IBM_Z_max)
# 이 입자들의 속도가 0이어야 하는데 왜 non-zero인지
# ========================================
print(f"\n{'='*80}")
print("건물 내부 유체 입자 속도 분포 (z-layer별)")
print(f"{'='*80}")

in_footprint = (fp[:,0] > IBM_X[0]) & (fp[:,0] < IBM_X[1]) & \
               (fp[:,1] > IBM_Y[0]) & (fp[:,1] < IBM_Y[1]) & \
               (fp[:,2] < IBM_Z[1])

fp_in = fp[in_footprint]
fv_in = fv[in_footprint]
fs_in = fs[in_footprint]

print(f"건물 footprint 내 유체: N={in_footprint.sum()}")
z_layers2 = np.arange(0, 35, 2.5)
print(f"\n{'z_center':>8} | {'N':>6} | {'|u|_mean':>10} | {'|u|_max':>10} | {'ux_mean':>10} | {'uz_mean':>10}")
for i in range(len(z_layers2)-1):
    zm = (fp_in[:,2] >= z_layers2[i]) & (fp_in[:,2] < z_layers2[i+1])
    if zm.sum() > 0:
        zc = (z_layers2[i]+z_layers2[i+1])/2
        print(f"  {zc:>6.1f} | {zm.sum():>6} | {fs_in[zm].mean():>10.3f} | {fs_in[zm].max():>10.3f} | "
              f"{fv_in[zm,0].mean():>10.3f} | {fv_in[zm,2].mean():>10.3f}")

# 건물 내부에서 속도가 높은 입자의 위치 (edge 근처인지?)
high_inside = fs_in > 5.0
if high_inside.any():
    hp = fp_in[high_inside]
    print(f"\n  건물 내부 |u|>5 m/s 입자: N={high_inside.sum()}")
    # edge distance
    dx_min = np.minimum(hp[:,0]-IBM_X[0], IBM_X[1]-hp[:,0])
    dy_min = np.minimum(hp[:,1]-IBM_Y[0], IBM_Y[1]-hp[:,1])
    dz_top = IBM_Z[1] - hp[:,2]
    print(f"    x edge 최소거리: mean={dx_min.mean():.1f}, min={dx_min.min():.1f}")
    print(f"    y edge 최소거리: mean={dy_min.mean():.1f}, min={dy_min.min():.1f}")
    print(f"    z top 거리:      mean={dz_top.mean():.1f}, min={dz_top.min():.1f}")
    print(f"    Edge(<5m) 비율: x={((dx_min<5).sum()/len(hp))*100:.0f}%, y={((dy_min<5).sum()/len(hp))*100:.0f}%, z_top={((dz_top<5).sum()/len(hp))*100:.0f}%")

# ========================================
# 시간 발전: 고속도 입자 수 추이
# ========================================
import glob, os
print(f"\n{'='*80}")
print("시간별 비물리적 고속도 유체 입자 수 추이")
print(f"{'='*80}")

files = sorted(glob.glob("Leipzig_small_point/fluid_condition_*stp.vtk"),
               key=lambda f: int(re.search(r"_(\d+)stp", f).group(1)))

print(f"\n{'step':>6} | {'>20m/s':>8} | {'>50m/s':>8} | {'|u|_max':>10} | {'z=1.25층 max':>12} | {'z=1.25층 mean':>12}")
for fpath in files:
    step = int(re.search(r"_(\d+)stp", fpath).group(1))
    if step % 200 != 0 and step != 0:
        continue
    p, f = parse_vtk(fpath)
    v = f["velocity"]
    pt = f["p_type"].astype(int)
    sp = np.sqrt(v[:,0]**2 + v[:,1]**2 + v[:,2]**2)
    fl = pt == 1
    sp_fl = sp[fl]
    z_fl = p[fl, 2]

    z_low_mask = z_fl < 2.0
    z_low_sp = sp_fl[z_low_mask]

    n20 = (sp_fl > 20).sum()
    n50 = (sp_fl > 50).sum()
    umax = sp_fl.max()
    zmax = z_low_sp.max() if z_low_mask.any() else 0
    zmean = z_low_sp.mean() if z_low_mask.any() else 0

    print(f"  {step:>5} | {n20:>8} | {n50:>8} | {umax:>10.1f} | {zmax:>12.1f} | {zmean:>12.2f}")
