#!/usr/bin/env python3
"""Deeper analysis: which particles blow up, and what's the IBM force magnitude."""
import numpy as np
import os, sys, re, glob

DTYPE_MAP = {
    "float":  (">f4", 4),
    "double": (">f8", 8),
    "int":    (">i4", 4),
    "unsigned_int": (">u4", 4),
    "short":  (">i2", 2),
    "char":   (">i1", 1),
}

def parse_vtk_fluid_full(filepath):
    with open(filepath, "rb") as f:
        content = f.read()
    m = re.search(rb"POINTS\s+(\d+)\s+(\w+)\n", content)
    if not m:
        return None, {}
    n_pts = int(m.group(1))
    dt_str = m.group(2).decode("ascii")
    pts_offset = m.end()
    dt, sz = DTYPE_MAP.get(dt_str, (">f4", 4))
    pts = np.frombuffer(content, dtype=dt, count=n_pts*3, offset=pts_offset).reshape(n_pts, 3).astype(np.float64)

    fields = {}
    m_field = re.search(rb"FIELD FieldData\s+(\d+)\n", content)
    if not m_field:
        return pts, fields
    n_fields = int(m_field.group(1))
    offset = m_field.end()
    for _ in range(n_fields):
        newline = content.index(b"\n", offset)
        header = content[offset:newline].decode("ascii").split()
        fname, ncomp, ntuple, dtype_s = header[0], int(header[1]), int(header[2]), header[3]
        offset = newline + 1
        dt2, itemsize = DTYPE_MAP.get(dtype_s, (">f4", 4))
        nbytes = ncomp * ntuple * itemsize
        arr = np.frombuffer(content, dtype=dt2, count=ncomp*ntuple, offset=offset).astype(np.float64)
        if ncomp > 1:
            arr = arr.reshape(ntuple, ncomp)
        fields[fname] = arr
        offset += nbytes
    return pts, fields

def get_step(fpath):
    m = re.search(r"fluid_(\d+)stp", fpath)
    return int(m.group(1)) if m else -1

def analyze(vtk_dir):
    # Focus on steps around blowup
    target_steps = [200, 500, 800, 1000, 1100, 1200, 1300, 1400, 1500, 1800, 2000, 2100, 2200]
    files = sorted(glob.glob(os.path.join(vtk_dir, "fluid_[0-9]*stp.vtk")), key=get_step)
    files = [f for f in files if get_step(f) in target_steps]

    for fpath in files:
        step = get_step(fpath)
        pts, fields = parse_vtk_fluid_full(fpath)
        if pts is None:
            continue

        print(f"\n{'='*80}")
        print(f"STEP {step}")
        print(f"{'='*80}")

        vel = fields.get("velocity")
        ptype = fields.get("p_type")
        pres = fields.get("pressure")

        if vel is None or ptype is None:
            print("  Missing velocity or p_type field")
            continue

        ptype = ptype.astype(int)
        speed = np.sqrt(vel[:,0]**2 + vel[:,1]**2 + vel[:,2]**2)

        # p_type distribution
        unique_pt, counts_pt = np.unique(ptype, return_counts=True)
        print(f"  p_type distribution: {dict(zip(unique_pt, counts_pt))}")

        # For each p_type, show velocity stats
        for pt_val in unique_pt:
            mask = ptype == pt_val
            sp = speed[mask]
            nan_count = np.isnan(sp).sum()
            if nan_count == len(sp):
                print(f"  p_type={pt_val}: N={mask.sum()}, ALL NaN")
                continue
            valid_sp = sp[~np.isnan(sp)]
            if len(valid_sp) == 0:
                continue
            z_vals = pts[mask, 2]
            print(f"  p_type={pt_val}: N={mask.sum()}, |u|: mean={np.mean(valid_sp):.3f}, "
                  f"max={np.max(valid_sp):.3f}, NaN={nan_count}, "
                  f"z: [{np.min(z_vals):.1f}, {np.max(z_vals):.1f}]")

        # Find top-10 highest velocity FLUID particles (p_type=1)
        fluid_mask = ptype == 1
        if fluid_mask.any():
            fluid_speed = speed[fluid_mask]
            fluid_pts = pts[fluid_mask]
            fluid_vel = vel[fluid_mask]
            valid = ~np.isnan(fluid_speed)
            if valid.any():
                valid_speed = fluid_speed[valid]
                if np.max(valid_speed) > 20:
                    top_idx = np.argsort(valid_speed)[-10:][::-1]
                    print(f"\n  Top-10 fastest FLUID particles:")
                    print(f"  {'|u|':>8} | {'ux':>8} | {'uy':>8} | {'uz':>8} | {'x':>8} | {'y':>8} | {'z':>8}")
                    valid_pts = fluid_pts[valid]
                    valid_vel = fluid_vel[valid]
                    for idx in top_idx:
                        p = valid_pts[idx]
                        v = valid_vel[idx]
                        s = valid_speed[idx]
                        print(f"  {s:>8.2f} | {v[0]:>8.2f} | {v[1]:>8.2f} | {v[2]:>8.2f} | {p[0]:>8.1f} | {p[1]:>8.1f} | {p[2]:>8.1f}")

        # Check IBM particles (p_type >= 1000 in code, but what in VTK?)
        # Show any p_type > 1 to find IBM
        for pt_val in unique_pt:
            if pt_val > 1:
                mask = ptype == pt_val
                sp = speed[mask]
                z_vals = pts[mask, 2]
                print(f"\n  p_type={pt_val} (possible IBM): N={mask.sum()}, "
                      f"|u|: max={np.nanmax(sp):.4f}, mean={np.nanmean(sp):.4f}")

if __name__ == "__main__":
    vtk_dir = sys.argv[1] if len(sys.argv) > 1 else "Leipzig_small_point"
    analyze(vtk_dir)
