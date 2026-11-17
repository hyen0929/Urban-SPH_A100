#!/usr/bin/env python3
"""Analyze velocity fields in fluid VTK files to diagnose blowup."""
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
    """Parse fluid VTK: return pts, and dict of field_name -> array."""
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
    files = sorted(glob.glob(os.path.join(vtk_dir, "fluid_[0-9]*stp.vtk")),
                   key=get_step)
    if not files:
        print(f"No fluid VTK files found in {vtk_dir}")
        return

    print(f"Found {len(files)} fluid VTK files")
    print(f"\n{'step':>6} | {'N_pts':>6} | {'|u|_max':>10} | {'|u|_mean':>10} | {'NaN_u':>6} | {'p_max':>10} | {'p_min':>10} | {'NaN_p':>6} | {'max_loc (x,y,z)':>30}")
    print("-" * 110)

    for fpath in files:
        step = get_step(fpath)
        pts, fields = parse_vtk_fluid_full(fpath)
        if pts is None:
            print(f"{step:>6} | PARSE ERROR")
            continue

        n_pts = len(pts)
        line = f"{step:>6} | {n_pts:>6}"

        if "velocity" in fields:
            vel = fields["velocity"]
            speed = np.sqrt(vel[:,0]**2 + vel[:,1]**2 + vel[:,2]**2)
            nan_u = np.isnan(speed).sum()
            valid = ~np.isnan(speed)
            if valid.any():
                u_max = np.nanmax(speed)
                u_mean = np.nanmean(speed)
                idx_max = np.nanargmax(speed)
                loc = pts[idx_max]
                line += f" | {u_max:>10.3f} | {u_mean:>10.3f} | {nan_u:>6} "
            else:
                line += f" | {'ALL NaN':>10} | {'---':>10} | {nan_u:>6} "
                loc = None
        else:
            line += f" | {'no vel':>10} | {'---':>10} | {'---':>6} "
            loc = None

        if "pressure" in fields:
            pres = fields["pressure"]
            nan_p = np.isnan(pres).sum()
            if (~np.isnan(pres)).any():
                p_max = np.nanmax(pres)
                p_min = np.nanmin(pres)
                line += f"| {p_max:>10.1f} | {p_min:>10.1f} | {nan_p:>6}"
            else:
                line += f"| {'ALL NaN':>10} | {'---':>10} | {nan_p:>6}"
        else:
            line += f"| {'no pres':>10} | {'---':>10} | {'---':>6}"

        if loc is not None:
            line += f" | ({loc[0]:.1f}, {loc[1]:.1f}, {loc[2]:.1f})"

        print(line)

        # If NaN detected or extreme velocity, print p_type distribution
        if "velocity" in fields:
            speed = np.sqrt(fields["velocity"][:,0]**2 + fields["velocity"][:,1]**2 + fields["velocity"][:,2]**2)
            nan_count = np.isnan(speed).sum()
            u_max_val = np.nanmax(speed) if (~np.isnan(speed)).any() else 0
            if nan_count > 0 or u_max_val > 50:
                if "p_type" in fields:
                    pt = fields["p_type"].astype(int)
                    # Show p_type breakdown for high-velocity particles
                    mask_high = speed > 30 if u_max_val > 30 else np.isnan(speed)
                    if mask_high.any():
                        unique, counts = np.unique(pt[mask_high], return_counts=True)
                        ptype_str = ", ".join([f"type{u}:{c}" for u, c in zip(unique, counts)])
                        print(f"       |  -> High-vel p_type: {ptype_str}")
                    # z-range of high velocity particles
                    if (~np.isnan(speed)).any():
                        mask_extreme = speed > 0.8 * u_max_val
                        if mask_extreme.any():
                            z_extreme = pts[mask_extreme, 2]
                            print(f"       |  -> Extreme vel z-range: [{z_extreme.min():.1f}, {z_extreme.max():.1f}]")

if __name__ == "__main__":
    vtk_dir = sys.argv[1] if len(sys.argv) > 1 else "Leipzig_small_point"
    analyze(vtk_dir)
