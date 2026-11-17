#!/usr/bin/env python3
"""Analyze IBM particle geometry near blowup locations."""
import numpy as np
import re

DTYPE_MAP = {
    "float":  (">f4", 4), "double": (">f8", 8),
    "int":    (">i4", 4), "unsigned_int": (">u4", 4),
    "short":  (">i2", 2), "char":   (">i1", 1),
}

def parse_vtk(filepath):
    with open(filepath, "rb") as f:
        content = f.read()
    m = re.search(rb"POINTS\s+(\d+)\s+(\w+)\n", content)
    if not m: return None, {}
    n_pts = int(m.group(1))
    dt, sz = DTYPE_MAP.get(m.group(2).decode("ascii"), (">f4", 4))
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

# Load building file
pts_b, fields_b = parse_vtk("Leipzig_small_point/building_0stp.vtk")
print(f"Building particles: {len(pts_b)}")
print(f"  z range: [{pts_b[:,2].min():.1f}, {pts_b[:,2].max():.1f}]")
print(f"  z unique (bottom 10): {sorted(np.unique(np.round(pts_b[:,2], 1)))[:10]}")

# Blowup locations (from analysis)
corners = [
    ("Corner A", 125, 65),
    ("Corner B", 125, 135),
    ("Corner C", 292, 142),
    ("Corner D", 137, 32),
]

print(f"\n=== IBM particle coverage near blowup corners ===")
for name, cx, cy in corners:
    r = 15  # search radius
    near = (np.abs(pts_b[:,0] - cx) < r) & (np.abs(pts_b[:,1] - cy) < r)
    if near.any():
        z_vals = sorted(np.unique(np.round(pts_b[near, 2], 1)))
        n = near.sum()
        print(f"\n{name} (~{cx},{cy}):")
        print(f"  N_IBM={n}, x: [{pts_b[near,0].min():.1f}, {pts_b[near,0].max():.1f}]")
        print(f"  y: [{pts_b[near,1].min():.1f}, {pts_b[near,1].max():.1f}]")
        print(f"  z layers: {z_vals[:15]}")
        print(f"  z_min={z_vals[0]:.1f}")
    else:
        print(f"\n{name} (~{cx},{cy}): No IBM particles within r={r}m")

# Load initial fluid file to check boundary particles near corners
pts_f, fields_f = parse_vtk("Leipzig_small_point/fluid_condition_0stp.vtk")
ptype = fields_f.get("p_type", np.zeros(len(pts_f))).astype(int)

print(f"\n=== Particle type distribution near corners (r=10m) ===")
for name, cx, cy in corners:
    r = 10
    near = (np.abs(pts_f[:,0] - cx) < r) & (np.abs(pts_f[:,1] - cy) < r) & (pts_f[:,2] < 10)
    if near.any():
        near_pts = pts_f[near]
        near_pt = ptype[near]
        print(f"\n{name} (~{cx},{cy}), z<10:")
        for pt_val in sorted(np.unique(near_pt)):
            mask = near_pt == pt_val
            z_vals = near_pts[mask, 2]
            print(f"  p_type={pt_val}: N={mask.sum()}, z: [{z_vals.min():.1f}, {z_vals.max():.1f}]")

# Check: is there IBM coverage at z=0 (ground-building junction)?
print(f"\n=== IBM particles at z < 3 (ground-building junction) ===")
low_ibm = pts_b[:,2] < 3
if low_ibm.any():
    z_low = sorted(np.unique(np.round(pts_b[low_ibm, 2], 1)))
    print(f"  N={low_ibm.sum()}, z values: {z_low}")
else:
    print(f"  NONE — IBM starts above z=3!")

# Check the latest step for velocity profile near a corner
pts_last, fields_last = parse_vtk("Leipzig_small_point/fluid_condition_1500stp.vtk")
vel_last = fields_last.get("velocity")
pt_last = fields_last.get("p_type").astype(int)
if vel_last is not None:
    print(f"\n=== Velocity profile near Corner B (125,135) at step 1500 ===")
    cx, cy = 125, 135
    fluid = pt_last == 1
    near = fluid & (np.abs(pts_last[:,0] - cx) < 5) & (np.abs(pts_last[:,1] - cy) < 5)
    if near.any():
        zz = pts_last[near, 2]
        spd = np.sqrt(vel_last[near,0]**2 + vel_last[near,1]**2 + vel_last[near,2]**2)
        ux = vel_last[near, 0]
        uz = vel_last[near, 2]
        order = np.argsort(zz)
        print(f"  {'z':>6} | {'|u|':>8} | {'ux':>8} | {'uz':>8}")
        for idx in order[:20]:
            print(f"  {zz[idx]:>6.1f} | {spd[idx]:>8.2f} | {ux[idx]:>8.2f} | {uz[idx]:>8.2f}")
