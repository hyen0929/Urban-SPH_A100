#!/usr/bin/env python3
"""
VTK binary FIELD 형식에서 concn peak 위치를 추적하는 스크립트
사용법: python3 track_peak.py <vtk_dir>
"""
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

def parse_vtk_field(filepath):
    """VTK binary FIELD 형식에서 points와 concn 배열을 반환"""
    with open(filepath, "rb") as f:
        content = f.read()

    # --- POINTS 위치 찾기 ---
    m = re.search(rb"POINTS\s+(\d+)\s+(\w+)\n", content)
    if not m:
        return None, None
    n_pts = int(m.group(1))
    pts_offset = m.end()
    pts = np.frombuffer(content, dtype=">f4", count=n_pts*3, offset=pts_offset).reshape(n_pts, 3)

    # --- FIELD FieldData 파싱 ---
    m_field = re.search(rb"FIELD FieldData\s+(\d+)\n", content)
    if not m_field:
        return pts, None
    n_fields = int(m_field.group(1))
    offset = m_field.end()

    for _ in range(n_fields):
        # 필드 헤더: "name ncomp ntuple dtype\n"
        newline = content.index(b"\n", offset)
        header = content[offset:newline].decode("ascii").split()
        fname, ncomp, ntuple, dtype_str = header[0], int(header[1]), int(header[2]), header[3]
        offset = newline + 1

        dtype, itemsize = DTYPE_MAP.get(dtype_str, (">f4", 4))
        nbytes = ncomp * ntuple * itemsize

        if fname == "concn":
            concn = np.frombuffer(content, dtype=dtype, count=ncomp*ntuple, offset=offset)
            return pts, concn.astype(np.float32)

        offset += nbytes

    return pts, None


def track_peaks(vtk_dir, ux=0.3, dt=None, x0=None):
    """
    ux : 유입 속도 [m/s] (default 0.3 for CASE_ID=2002)
    dt : timestep [s] (None이면 데이터에서 자동 추정)
    x0 : 소스 초기 x 위치 (None이면 step=0 peak_x에서 읽음)
    """
    files = sorted(glob.glob(os.path.join(vtk_dir, "fluid_[0-9]*stp.vtk")),
                   key=lambda f: int(re.search(r"fluid_(\d+)stp", f).group(1)))

    results = []
    for fpath in files:
        step = int(re.search(r"fluid_(\d+)stp", fpath).group(1))
        pts, concn = parse_vtk_field(fpath)
        if concn is None:
            continue
        if concn.max() < 1e-6:
            continue
        idx = np.argmax(concn)
        results.append((step, concn[idx], pts[idx, 0], pts[idx, 1], pts[idx, 2]))

    if not results:
        print("No valid data found.")
        return

    # x0 자동 설정
    if x0 is None:
        x0 = results[0][2]  # step=0 의 peak_x

    # dt 자동 추정 (step=0 제외, 2개 이상 데이터 이용)
    if dt is None and len(results) >= 3:
        # 후반부 데이터로 안정적인 dt 추정
        r1, r2 = results[1], results[-1]
        dt = (r2[2] - r1[2]) / (ux * (r2[0] - r1[0]))
        print(f"[Auto dt estimate] dt ≈ {dt:.6f} s/step  (ux={ux} m/s, x0={x0:.4f} m)\n")

    print(f"{'step':>8} | {'max_concn':>12} | {'sim_x':>8} | {'theory_x':>10} | {'error':>8} | {'error%':>7}")
    print("-" * 75)

    for step, cmax, px, py, pz in results:
        theory_x = x0 + ux * dt * step if dt else float('nan')
        err = px - theory_x
        err_pct = (err / (theory_x - x0) * 100) if (dt and (theory_x - x0) != 0) else float('nan')
        print(f"{step:>8} | {cmax:>12.4f} | {px:>8.4f} | {theory_x:>10.4f} | {err:>+8.4f} | {err_pct:>+7.2f}%")


if __name__ == "__main__":
    vtk_dir = sys.argv[1] if len(sys.argv) > 1 else "/home/kdh/Urban-SPH/Transport_advection"
    # 커맨드라인: python3 track_peak.py <dir> [ux] [dt]
    ux = float(sys.argv[2]) if len(sys.argv) > 2 else 0.3
    dt = float(sys.argv[3]) if len(sys.argv) > 3 else None
    track_peaks(vtk_dir, ux=ux, dt=dt)
