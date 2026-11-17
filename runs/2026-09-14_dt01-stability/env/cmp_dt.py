#!/usr/bin/env python3
"""dt=0.01 시험런을 Run A(dt=0.005) 기준곡선과 같은 물리시각에서 대조한다.

  env/cmp_dt.py            표 출력
  env/cmp_dt.py --check    판정만 (통과=0, 의심=1, 실패=2) — 감시 루프용
"""
import csv, os, re, sys, math

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..")
REF  = os.path.join(ROOT, "runs/2026-09-10_leipzig-full-bcfix_2")
TRI  = os.path.join(ROOT, "runs/2026-09-14_dt01-stability")

def ppe(run):
    """time -> (div_avg, div_plain_max, bi_max)"""
    p = os.path.join(run, "out/PPE_diag_explicit.csv")
    out = []
    if not os.path.exists(p): return out
    with open(p, newline="") as f:
        for r in csv.DictReader(f):
            try:
                out.append((float(r["time"]), float(r["div_un1_avg"]),
                            float(r["div_un1_plain_max"]), float(r["bi_max"])))
            except (ValueError, TypeError, KeyError):
                continue          # 정전으로 잘린 마지막 줄 등
    return out

STEP = re.compile(rb"^(\d+)\s+tu_max=([\d.eE+-]+).*?uz_max=([\d.eE+-]+)\s+ftotal_max=([\d.eE+-]+)", re.M)
SUMK = re.compile(rb"^(\d+)\s+tu_max=(?:.*\n)*?.*?\xce\xa3k_turb=([\d.eE+-]+)\s+\xce\xa3\xc2\xbd\|u\|\xc2\xb2=([\d.eE+-]+)", re.M)

def steps(run, dt):
    """time -> (tu_max, uz_max, ftotal_max)"""
    p = os.path.join(run, "run.log")
    if not os.path.exists(p): return []
    raw = open(p, "rb").read()
    return [(int(m.group(1))*dt, float(m.group(2)), float(m.group(3)), float(m.group(4)))
            for m in STEP.finditer(raw)]

def energy(run, dt):
    """time -> (Sk_turb, S0.5|u|^2) — 스텝 줄과 뒤따르는 [domain] 줄을 짝지어 시각을 붙인다"""
    p = os.path.join(run, "run.log")
    if not os.path.exists(p): return []
    raw = open(p, "rb").read()
    return [(int(m.group(1))*dt, float(m.group(2)), float(m.group(3)))
            for m in SUMK.finditer(raw)]

def at(series, t, tol):
    """t 에 가장 가까운 표본 (tol 초 이내)"""
    best, bd = None, tol
    for row in series:
        d = abs(row[0]-t)
        if d <= bd: best, bd = row, d
    return best

def main():
    check = "--check" in sys.argv
    ref_p, tri_p = ppe(REF), ppe(TRI)
    ref_s, tri_s = steps(REF, 0.005), steps(TRI, 0.01)
    ref_e, tri_e = energy(REF, 0.005), energy(TRI, 0.01)

    if not tri_s:
        print("시험런 아직 스텝 없음"); return 1
    tmax = tri_s[-1][0]
    grid = [t for t in (1.0,2.5,5.0,7.5,10.0,12.5,15.0,17.5,20.0,22.5,25.0,27.5,30.0) if t <= tmax + 1e-9]
    if not grid: grid = [tmax]

    worst_avg = worst_max = 0.0
    rows = []
    for t in grid:
        rp, tp = at(ref_p, t, 0.30), at(tri_p, t, 0.30)
        rs, ts = at(ref_s, t, 0.30), at(tri_s, t, 0.30)
        r_avg = rp[1] if rp else float("nan"); t_avg = tp[1] if tp else float("nan")
        r_max = rp[2] if rp else float("nan"); t_max = tp[2] if tp else float("nan")
        ra = t_avg/r_avg if (rp and tp and r_avg) else float("nan")
        rm = t_max/r_max if (rp and tp and r_max) else float("nan")
        if not math.isnan(ra): worst_avg = max(worst_avg, ra)
        if not math.isnan(rm): worst_max = max(worst_max, rm)
        rows.append((t,
                     rs[1] if rs else float("nan"), ts[1] if ts else float("nan"),
                     rs[2] if rs else float("nan"), ts[2] if ts else float("nan"),
                     r_avg, t_avg, ra, r_max, t_max, rm))

    # NaN / 발산 폭주 탐지
    bad = [r for r in tri_p if math.isnan(r[1]) or math.isinf(r[1]) or r[1] > 10.0]

    if not check:
        print(f"{'t[s]':>6} │ {'tu_max A':>8} {'dt01':>7} │ {'uz_max A':>8} {'dt01':>7} │"
              f" {'div_avg A':>10} {'dt01':>10} {'배':>5} │ {'div_pmax A':>10} {'dt01':>10} {'배':>5}")
        print("─"*112)
        for (t,ru,tu,rz,tz,ra_,ta_,ra,rx,tx,rm) in rows:
            f=lambda v,w,p=2: ("{:>%d.%df}"%(w,p)).format(v) if not math.isnan(v) else " "*(w-1)+"-"
            e=lambda v,w: ("{:>%d.2e}"%w).format(v) if not math.isnan(v) else " "*(w-1)+"-"
            print(f"{t:6.1f} │ {f(ru,8)} {f(tu,7)} │ {f(rz,8)} {f(tz,7)} │"
                  f" {e(ra_,10)} {e(ta_,10)} {f(ra,5,2)} │ {e(rx,10)} {e(tx,10)} {f(rm,5,2)}")
        print()
        print(f"시험런 도달: t={tmax:.2f} s / 30.0 s   최악비 div_avg×{worst_avg:.2f}  div_plain_max×{worst_max:.2f}")
        et = tri_e[-1][0] if tri_e else None
        ra_e = at(ref_e, et, 0.30) if et is not None else None
        ta_e = at(tri_e, et, 0.30) if et is not None else None
        if ra_e and ta_e:
            print(f"t={et:.1f} s 기준 — Σk_turb  A={ra_e[1]:.4e}  dt01={ta_e[1]:.4e} (×{ta_e[1]/ra_e[1]:.2f})"
                  f"   Σ½|u|²  A={ra_e[2]:.4e}  dt01={ta_e[2]:.4e} (×{ta_e[2]/ra_e[2]:.3f})")
        if bad: print(f"!! 비정상 div_avg 표본 {len(bad)}개 (NaN/inf/>10)")

    # 판정: div_avg 2배 초과 또는 div_plain_max 3배 초과 = 실패
    if bad: 
        if check: print(f"FAIL NaN/발산 표본 {len(bad)}개")
        return 2
    if worst_avg > 2.0 or worst_max > 3.0:
        if check: print(f"FAIL div_avg×{worst_avg:.2f} (>2) 또는 div_plain_max×{worst_max:.2f} (>3)")
        return 2
    if worst_avg > 1.5 or worst_max > 2.0:
        if check: print(f"WARN div_avg×{worst_avg:.2f}  div_plain_max×{worst_max:.2f}")
        return 1
    if check: print(f"OK div_avg×{worst_avg:.2f}  div_plain_max×{worst_max:.2f}  t={tmax:.1f}s")
    return 0

sys.exit(main())
