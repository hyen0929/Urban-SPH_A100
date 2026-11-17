#!/usr/bin/env python3
"""Check IBM particle velocities in step 200."""
import numpy as np
import re

DTYPE_MAP = {
    "float":  (">f4", 4), "double": (">f8", 8),
    "int":    (">i4", 4), "unsigned_int": (">u4", 4),
}

def parse_vtk(filepath):
    with open(filepath, "rb") as f:
        content = f.read()
    m = re.search(rb"POINTS\s+(\d+)\s+(\w+)\n", content)
    n_pts = int(m.group(1))
    dt, sz = DTYPE_MAP.get(m.group(2).decode("ascii"), (">f4", 4))
    pts = np.frombuffer(content, dtype=dt, count=n_pts*3, offset=m.end()).reshape(n_pts,3).astype(np.float64)

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

pts, fields = parse_vtk("Leipzig_small_point/fluid_200stp.vtk")
vel = fields["velocity"]
ptype = fields["p_type"].astype(int)

# IBM particles
ibm = ptype == 1000
ibm_pts = pts[ibm]
ibm_vel = vel[ibm]
ibm_speed = np.sqrt(ibm_vel[:,0]**2 + ibm_vel[:,1]**2 + ibm_vel[:,2]**2)

print(f"IBM particles: N={ibm.sum()}")
print(f"  |u|: min={ibm_speed.min():.4f}, mean={ibm_speed.mean():.4f}, max={ibm_speed.max():.4f}")
print(f"  ux: min={ibm_vel[:,0].min():.4f}, mean={ibm_vel[:,0].mean():.4f}, max={ibm_vel[:,0].max():.4f}")
print(f"  uy: min={ibm_vel[:,1].min():.4f}, mean={ibm_vel[:,1].mean():.4f}, max={ibm_vel[:,1].max():.4f}")
print(f"  uz: min={ibm_vel[:,2].min():.4f}, mean={ibm_vel[:,2].mean():.4f}, max={ibm_vel[:,2].max():.4f}")
print(f"  z: min={ibm_pts[:,2].min():.1f}, max={ibm_pts[:,2].max():.1f}")

# Show IBM velocity by z-layer
print(f"\n  IBM velocity profile by z:")
z_vals = np.unique(np.round(ibm_pts[:,2], 1))
for z in z_vals[:15]:
    mask = np.abs(ibm_pts[:,2] - z) < 0.3
    if mask.any():
        sp = ibm_speed[mask]
        ux = ibm_vel[mask, 0]
        print(f"    z={z:>6.1f}: N={mask.sum():>5}, |u|={sp.mean():.3f}, ux={ux.mean():.3f}")
