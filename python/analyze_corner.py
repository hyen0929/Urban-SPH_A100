#!/usr/bin/env python3
"""Analyze velocity growth at building corner in fluid_condition VTK files."""
import numpy as np
import os, sys, re, glob

DTYPE_MAP = {
    "float":  (">f4", 4), "double": (">f8", 8),
    "int":    (">i4", 4), "unsigned_int": (">u4", 4),
    "short":  (">i2", 2), "char":   (">i1", 1),
}

def parse_vtk(filepath):
    with open(filepath, "rb") as f:
        content = f.read()
    m = re.search(rb"POINTS\s+(\d+)\s+(\w+)\n", content)
    if not m:
        return None, {}
    n_pts = int(m.group(1))
    dt, sz = DTYPE_MAP.get(m.group(2).decode("ascii"), (">f4", 4))
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

def get_step(fpath):
    m = re.search(r"_(\d+)stp", fpath)
    return int(m.group(1)) if m else -1

vtk_dir = sys.argv[1] if len(sys.argv) > 1 else "Leipzig_small_point"
files = sorted(glob.glob(os.path.join(vtk_dir, "fluid_condition_*stp.vtk")), key=get_step)

print(f"{'step':>6} | {'N':>6} | {'|u|max':>8} | {'|u|mean':>8} | {'#>15':>5} | {'top_ux':>8} | {'top_uz':>8} | {'top_x':>6} | {'top_y':>6} | {'top_z':>6} | {'fb_mag':>10}")
print("-" * 120)

for fpath in files:
    step = get_step(fpath)
    pts, fields = parse_vtk(fpath)
    if pts is None: continue
    vel = fields.get("velocity")
    ptype = fields.get("p_type")
    if vel is None or ptype is None: continue
    ptype = ptype.astype(int)

    # fluid only
    fluid = ptype == 1
    if not fluid.any(): continue
    fvel = vel[fluid]
    fpts = pts[fluid]
    speed = np.sqrt(fvel[:,0]**2 + fvel[:,1]**2 + fvel[:,2]**2)
    above_lim = (speed > 15).sum()

    # IBM fb force on fluid
    fb_info = ""
    if "fbx" in fields or "fb_x" in fields:
        fbx_key = "fbx" if "fbx" in fields else "fb_x"
        fb = fields[fbx_key]
        if fb.ndim == 2:
            fb_mag = np.sqrt(fb[fluid,0]**2 + fb[fluid,1]**2 + fb[fluid,2]**2)
        else:
            fb_mag = np.abs(fb[fluid])
        fb_info = f"{fb_mag.max():.4f}"
    else:
        fb_info = "N/A"

    idx_max = np.argmax(speed)
    p = fpts[idx_max]
    v = fvel[idx_max]
    print(f"{step:>6} | {fluid.sum():>6} | {speed.max():>8.2f} | {speed.mean():>8.3f} | {above_lim:>5} | {v[0]:>8.2f} | {v[2]:>8.2f} | {p[0]:>6.1f} | {p[1]:>6.1f} | {p[2]:>6.1f} | {fb_info:>10}")

# Detailed analysis of the last file
if files:
    step = get_step(files[-1])
    pts, fields = parse_vtk(files[-1])
    vel = fields.get("velocity")
    ptype = fields.get("p_type").astype(int)
    fluid = ptype == 1
    fvel = vel[fluid]
    fpts = pts[fluid]
    speed = np.sqrt(fvel[:,0]**2 + fvel[:,1]**2 + fvel[:,2]**2)

    print(f"\n=== Step {step}: Particles with |u| > 12 ===")
    mask_high = speed > 12
    if mask_high.any():
        hi_speed = speed[mask_high]
        hi_pts = fpts[mask_high]
        hi_vel = fvel[mask_high]
        order = np.argsort(hi_speed)[::-1]
        print(f"  {'|u|':>8} | {'ux':>8} | {'uy':>8} | {'uz':>8} | {'x':>8} | {'y':>8} | {'z':>6}")
        for idx in order[:20]:
            s, v, p = hi_speed[idx], hi_vel[idx], hi_pts[idx]
            print(f"  {s:>8.2f} | {v[0]:>8.2f} | {v[1]:>8.2f} | {v[2]:>8.2f} | {p[0]:>8.1f} | {p[1]:>8.1f} | {p[2]:>6.1f}")

    # check IBM particles
    ibm = ptype == 1000
    if ibm.any():
        ibm_vel = vel[ibm]
        ibm_speed = np.sqrt(ibm_vel[:,0]**2 + ibm_vel[:,1]**2 + ibm_vel[:,2]**2)
        print(f"\n  IBM particles: N={ibm.sum()}, |u|max={ibm_speed.max():.4f}")

    # Show all field names
    print(f"\n  Available fields: {list(fields.keys())}")
