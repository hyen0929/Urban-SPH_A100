#!/usr/bin/env python3
"""Analyze IBM velocity issues: non-physical high velocities + boundary far-side anomaly."""
import numpy as np
import re, glob, os, sys

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

def analyze_step(fpath):
    step = int(re.search(r"_(\d+)stp", fpath).group(1))
    pts, fields = parse_vtk(fpath)
    vel = fields["velocity"]
    ptype = fields["p_type"].astype(int)
    speed = np.sqrt(vel[:,0]**2 + vel[:,1]**2 + vel[:,2]**2)

    print(f"\n{'='*90}")
    print(f"STEP {step}  (N={len(pts)})")
    print(f"{'='*90}")

    # --- 1. p_type별 속도 통계 ---
    print(f"\n[1] p_type별 속도 통계")
    print(f"  {'p_type':>8} | {'N':>7} | {'|u|_mean':>10} | {'|u|_max':>10} | {'|u|_99%':>10} | {'NaN':>5}")
    for pt in sorted(np.unique(ptype)):
        m = ptype == pt
        sp = speed[m]
        nan_n = np.isnan(sp).sum()
        valid = sp[~np.isnan(sp)]
        if len(valid) == 0:
            print(f"  {pt:>8} | {m.sum():>7} | {'ALL NaN':>10} | {'---':>10} | {'---':>10} | {nan_n:>5}")
            continue
        p99 = np.percentile(valid, 99)
        print(f"  {pt:>8} | {m.sum():>7} | {valid.mean():>10.3f} | {valid.max():>10.3f} | {p99:>10.3f} | {nan_n:>5}")

    # --- 2. 유체 입자 중 비물리적 고속도 ---
    fluid = ptype == 1
    fluid_sp = speed[fluid]
    fluid_pts = pts[fluid]
    fluid_vel = vel[fluid]
    valid_fluid = ~np.isnan(fluid_sp)

    if valid_fluid.any():
        vsp = fluid_sp[valid_fluid]
        p999 = np.percentile(vsp, 99.9)
        p99 = np.percentile(vsp, 99)
        median = np.median(vsp)

        print(f"\n[2] 유체 입자 속도 분포")
        print(f"  median={median:.3f}, 99%={p99:.3f}, 99.9%={p999:.3f}, max={vsp.max():.3f}")

        # 비물리적 고속도 입자 (99.9% 이상)
        extreme_mask = fluid_sp > p999
        extreme_mask_valid = extreme_mask & valid_fluid
        if extreme_mask_valid.any():
            print(f"\n  비물리적 고속도 입자 (>{p999:.1f} m/s): N={extreme_mask_valid.sum()}")
            ep = pts[fluid][extreme_mask[fluid_sp > -1]]  # 전체 fluid mask에서
            # 위치 분석
            extreme_pts = fluid_pts[extreme_mask_valid[valid_fluid]]
            extreme_vel = fluid_vel[extreme_mask_valid[valid_fluid]]
            extreme_sp = fluid_sp[extreme_mask_valid[valid_fluid]]
            print(f"  x: [{extreme_pts[:,0].min():.1f}, {extreme_pts[:,0].max():.1f}]")
            print(f"  y: [{extreme_pts[:,1].min():.1f}, {extreme_pts[:,1].max():.1f}]")
            print(f"  z: [{extreme_pts[:,2].min():.1f}, {extreme_pts[:,2].max():.1f}]")

            # Top 10
            top_idx = np.argsort(extreme_sp)[-min(10, len(extreme_sp)):][::-1]
            print(f"\n  Top-10 fastest fluid particles:")
            print(f"  {'|u|':>8} | {'ux':>8} | {'uy':>8} | {'uz':>8} | {'x':>8} | {'y':>8} | {'z':>8}")
            for idx in top_idx:
                p = extreme_pts[idx]
                v = extreme_vel[idx]
                s = extreme_sp[idx]
                print(f"  {s:>8.2f} | {v[0]:>8.2f} | {v[1]:>8.2f} | {v[2]:>8.2f} | {p[0]:>8.1f} | {p[1]:>8.1f} | {p[2]:>8.1f}")

    # --- 3. IBM 입자 분석: boundary 안쪽(유체와 먼 쪽) 속도 이상 ---
    ibm = ptype >= 1000
    if ibm.any():
        ibm_pts = pts[ibm]
        ibm_vel = vel[ibm]
        ibm_sp = speed[ibm]

        print(f"\n[3] IBM 입자 분석 (N={ibm.sum()})")
        print(f"  |u|: mean={np.nanmean(ibm_sp):.4f}, max={np.nanmax(ibm_sp):.4f}")
        print(f"  위치 범위: x=[{ibm_pts[:,0].min():.1f},{ibm_pts[:,0].max():.1f}], "
              f"y=[{ibm_pts[:,1].min():.1f},{ibm_pts[:,1].max():.1f}], "
              f"z=[{ibm_pts[:,2].min():.1f},{ibm_pts[:,2].max():.1f}]")

        # IBM velocity는 corrector에서 0으로 set되므로 0이어야 함
        if np.nanmax(ibm_sp) > 1e-6:
            print(f"  *** WARNING: IBM 입자 속도가 0이 아님! max={np.nanmax(ibm_sp):.6f}")

    # --- 4. IBM 경계 근처 유체 입자 z-프로필 (boundary far-side 분석) ---
    if ibm.any() and valid_fluid.any():
        # IBM 입자의 z 범위로 boundary 위치 추정
        ibm_z_min = ibm_pts[:, 2].min()
        ibm_z_max = ibm_pts[:, 2].max()
        ibm_x_min = ibm_pts[:, 0].min()
        ibm_x_max = ibm_pts[:, 0].max()

        print(f"\n[4] IBM 경계 근처 유체 속도 z-프로필")
        print(f"  IBM z 범위: [{ibm_z_min:.1f}, {ibm_z_max:.1f}]")

        # IBM 근처 x 범위의 유체만
        fp = fluid_pts[valid_fluid]
        fv = fluid_vel[valid_fluid]
        fs = fluid_sp[valid_fluid]

        near_ibm_x = (fp[:, 0] > ibm_x_min - 50) & (fp[:, 0] < ibm_x_max + 50)

        if near_ibm_x.any():
            fp_near = fp[near_ibm_x]
            fv_near = fv[near_ibm_x]
            fs_near = fs[near_ibm_x]

            # z 방향 프로필 (IBM boundary 주변)
            z_bins = np.arange(max(0, ibm_z_min - 30), ibm_z_max + 60, 5.0)
            print(f"  {'z_center':>8} | {'N':>6} | {'|u|_mean':>10} | {'|u|_max':>10} | {'ux_mean':>10}")
            for i in range(len(z_bins)-1):
                z_mask = (fp_near[:, 2] >= z_bins[i]) & (fp_near[:, 2] < z_bins[i+1])
                if z_mask.sum() > 0:
                    z_center = (z_bins[i] + z_bins[i+1]) / 2
                    sp_z = fs_near[z_mask]
                    ux_z = fv_near[z_mask, 0]
                    print(f"  {z_center:>8.1f} | {z_mask.sum():>6} | {sp_z.mean():>10.3f} | {sp_z.max():>10.3f} | {ux_z.mean():>10.3f}")

    # --- 5. filt(concentration) 분석 - IBM force interpolation 품질 ---
    # concentration이 VTK에 없으므로, 속도 분포로 간접 확인
    # IBM 입자 위치별 주변 유체 속도 gradient 분석
    if ibm.any() and valid_fluid.any():
        print(f"\n[5] IBM 경계 표면별 유체 속도 분석")
        # IBM 입자 x좌표로 건물 표면 추정
        ibm_x_mean = ibm_pts[:, 0].mean()
        ibm_y_mean = ibm_pts[:, 1].mean()

        # 건물 각 면 근처 유체 속도
        dx = 10.0  # 표면에서의 거리
        surfaces = {
            "Left  (x<IBM_xmin)": (fp[:, 0] > ibm_x_min - 2*dx) & (fp[:, 0] < ibm_x_min) & (fp[:, 2] > ibm_z_min) & (fp[:, 2] < ibm_z_max),
            "Right (x>IBM_xmax)": (fp[:, 0] > ibm_x_max) & (fp[:, 0] < ibm_x_max + 2*dx) & (fp[:, 2] > ibm_z_min) & (fp[:, 2] < ibm_z_max),
            "Top   (z>IBM_zmax)": (fp[:, 2] > ibm_z_max) & (fp[:, 2] < ibm_z_max + 2*dx) & (fp[:, 0] > ibm_x_min) & (fp[:, 0] < ibm_x_max),
            "Below (z<IBM_zmin)": (fp[:, 2] > ibm_z_min - 2*dx) & (fp[:, 2] < ibm_z_min) & (fp[:, 0] > ibm_x_min) & (fp[:, 0] < ibm_x_max),
        }
        for name, mask in surfaces.items():
            if mask.sum() > 0:
                sp_s = fs[near_ibm_x][mask[near_ibm_x] if len(mask) == len(fp) else mask]
                # Recalculate with proper indexing
                local_sp = fs_near[mask[:len(fs_near)] if len(mask) >= len(fs_near) else mask]
                print(f"  {name}: N={mask.sum()}")

        # 더 단순한 접근: IBM z 범위 안/밖으로 나눠서 비교
        inside_z = (fp[:, 2] >= ibm_z_min) & (fp[:, 2] <= ibm_z_max)
        above_z = fp[:, 2] > ibm_z_max
        below_z = fp[:, 2] < ibm_z_min

        for label, zmask in [("z < IBM (below)", below_z), ("IBM z range", inside_z), ("z > IBM (above)", above_z)]:
            if zmask.sum() > 0:
                print(f"  {label:>20}: N={zmask.sum():>7}, |u|_mean={fs[zmask].mean():.3f}, |u|_max={fs[zmask].max():.3f}")


if __name__ == "__main__":
    vtk_dir = sys.argv[1] if len(sys.argv) > 1 else "Leipzig_small_point"
    # Analyze specific steps
    target_steps = [500, 1000, 1500, 2000, 2500, 3000]
    files = sorted(glob.glob(os.path.join(vtk_dir, "fluid_condition_*stp.vtk")),
                   key=lambda f: int(re.search(r"_(\d+)stp", f).group(1)))

    if not files:
        print(f"No files found in {vtk_dir}")
        sys.exit(1)

    print(f"Available steps: {[int(re.search(r'_(\d+)stp', f).group(1)) for f in files[:5]]}...{[int(re.search(r'_(\d+)stp', f).group(1)) for f in files[-5:]]}")

    for fpath in files:
        step = int(re.search(r"_(\d+)stp", fpath).group(1))
        if step in target_steps:
            analyze_step(fpath)
