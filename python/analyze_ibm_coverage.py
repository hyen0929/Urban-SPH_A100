#!/usr/bin/env python3
"""Check IBM coverage for blowup particles at step 1500."""
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

# Load building and fluid
pts_b, _ = parse_vtk("Leipzig_small_point/building_0stp.vtk")
pts_f, fields_f = parse_vtk("Leipzig_small_point/fluid_condition_1500stp.vtk")
vel = fields_f["velocity"]
ptype = fields_f["p_type"].astype(int)

fluid = ptype == 1
fpts = pts_f[fluid]
fvel = vel[fluid]
speed = np.sqrt(fvel[:,0]**2 + fvel[:,1]**2 + fvel[:,2]**2)

# For each high-velocity fluid particle, find distance to nearest IBM particle
high = speed > 15
hi_pts = fpts[high]
hi_speed = speed[high]

# Parameters
h = 1.5 * 2.5  # h_coeff * dx = 3.75m
IBM_h = 1.0 * h  # IBM_length * h = 3.75m (current setting)
search_range = 2.0 * IBM_h  # = 7.5m

print(f"Parameters: h={h:.2f}m, IBM_h={IBM_h:.2f}m, search_range={search_range:.2f}m")
print(f"High-velocity particles (|u|>15): {high.sum()}")
print(f"\n{'|u|':>8} | {'x':>6} | {'y':>6} | {'z':>5} | {'dist_IBM':>8} | {'N_IBM_near':>10} | {'in_range':>8}")
print("-" * 80)

order = np.argsort(hi_speed)[::-1]
for idx in order[:30]:
    p = hi_pts[idx]
    s = hi_speed[idx]
    # Distance to nearest IBM particle
    dx = pts_b[:,0] - p[0]
    dy = pts_b[:,1] - p[1]
    dz = pts_b[:,2] - p[2]
    dist = np.sqrt(dx*dx + dy*dy + dz*dz)
    min_dist = dist.min()
    n_ibm_near = (dist < search_range).sum()
    in_range = "YES" if min_dist < search_range else "NO"
    print(f"{s:>8.1f} | {p[0]:>6.1f} | {p[1]:>6.1f} | {p[2]:>5.1f} | {min_dist:>8.2f} | {n_ibm_near:>10} | {in_range:>8}")
