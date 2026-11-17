#!/usr/bin/env python3
"""본 런(dt=0.01)을 정전 전 Run A(dt=0.005, t=0~237 s 기록)와 같은 물리시각에서 계속 대조한다.

비교량 (dt 무관하게 맞춘 값):
  - PPE_diag div_un1_avg, div_un1_plain_max : 기록값은 rho*div(u)/dt 이므로 raw*dt 로 비교
  - run.log [domain] Σ½|u|², Σk_turb       : 50 step 마다 찍히는 도메인 합

판정 (Run A 기록이 있는 t<=237 s 에서만):
  FAIL  최근 2 s 중앙값 div_avg 비 > 2.0, 또는 div_plain_max 비 > 3.0, 또는 |Σ½|u|² 비 - 1| > 0.10, 또는 NaN
  WARN  |Σ½|u|² 비 - 1| > 0.05, 또는 Σk_turb 비가 [0.8, 1.25] 밖
  시험런(2026-09-14_dt01-stability, t=0~30 s)에서 div_avg 비는 1.3~1.9, div_plain_max 비는 단발 2.98 까지
  나왔으므로 단발값이 아니라 2 s 이동 중앙값으로 본다.
재개(restart) 뒤 겹쳐 기록된 행은 step 키로 덮어써 마지막 값만 쓴다. 시각은 CSV time 열이 아니라 step*dt.

사용: track_vs_runA.py [--every T] [--check]
  --check : 마지막 줄에 "OK|WARN|FAIL ..." 요약, exit 0/1/2
"""
import csv, math, os, re, sys, argparse, statistics

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.dirname(HERE)
REF = os.path.join(RUN, '..', '2026-09-10_leipzig-full-bcfix_2')
DT_NEW, DT_REF = 0.01, 0.005
WIN = 2.0

SUMK = re.compile(rb"^(\d+)\s+tu_max=(?:.*\n)*?.*?\xce\xa3k_turb=([\d.eE+-]+|nan|inf)\s+"
                  rb"\xce\xa3\xc2\xbd\|u\|\xc2\xb2=([\d.eE+-]+|nan|inf)", re.M)


def load_diag(run, dt):
    # 시각 키는 CSV 의 time 열이 아니라 step*dt 로 만든다. time 열은 float 누적오차가 쌓이고
    # (정전 전 Run A: t=220 s 에서 +0.059 s) %.6g 로 찍혀 t>=100 s 에선 소수 3자리뿐이라,
    # 0.01 반올림으로 두 런을 맞추면 겹치는 시각이 사라져 비교 없이 OK 가 나온다.
    rows = {}
    p = os.path.join(run, 'out', 'PPE_diag_explicit.csv')
    if not os.path.exists(p):
        return rows
    with open(p) as f:
        for r in csv.DictReader(f):
            try:
                rows[round(int(r['step']) * dt, 2)] = (float(r['div_un1_avg']), float(r['div_un1_plain_max']))
            except (ValueError, TypeError, KeyError):
                continue
    return rows


def load_energy(run, dt):
    out = {}
    p = os.path.join(run, 'run.log')
    if not os.path.exists(p):
        return out
    for m in SUMK.finditer(open(p, 'rb').read()):
        out[round(int(m.group(1)) * dt, 2)] = (float(m.group(2)), float(m.group(3)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--every', type=float, default=10.0)
    ap.add_argument('--check', action='store_true')
    a = ap.parse_args()

    dn, dr = load_diag(RUN, DT_NEW), load_diag(REF, DT_REF)
    en, er = load_energy(RUN, DT_NEW), load_energy(REF, DT_REF)
    if not dn:
        print('WAIT 본 런 PPE_diag 없음')
        sys.exit(0)

    t_last = max(dn)
    nan = [t for t, v in dn.items() if not all(map(math.isfinite, v))]
    ts = sorted(t for t in dn if t in dr)
    ra = {t: dn[t][0] * DT_NEW / (dr[t][0] * DT_REF) for t in ts}
    rp = {t: dn[t][1] * DT_NEW / (dr[t][1] * DT_REF) for t in ts}
    te = sorted(t for t in en if t in er)
    rke = {t: en[t][1] / er[t][1] for t in te}
    rk = {t: en[t][0] / er[t][0] for t in te}

    def med(d, t):
        v = [x for s, x in d.items() if t - WIN < s <= t]
        return statistics.median(v) if v else float('nan')

    print(f"{'t[s]':>7} {'avg비(2s중앙)':>13} {'pmax비(2s중앙)':>14} {'KE비':>7} {'k비':>7}")
    last = -1e9
    for t in ts:
        if t - last >= a.every - 1e-9 or t == ts[-1]:
            tk = min(te, key=lambda s: abs(s - t)) if te else None
            ke = rke[tk] if tk is not None and abs(tk - t) < 0.6 else float('nan')
            kk = rk[tk] if tk is not None and abs(tk - t) < 0.6 else float('nan')
            print(f"{t:7.2f} {med(ra, t):13.2f} {med(rp, t):14.2f} {ke:7.3f} {kk:7.3f}")
            last = t

    level, why = 0, []
    if nan:
        level, why = 2, [f'NaN/Inf t={nan[0]}']
    if ts:
        t = ts[-1]
        ma, mp = med(ra, t), med(rp, t)
        if ma > 2.0: level = 2; why.append(f'div_avg비 {ma:.2f}>2')
        if mp > 3.0: level = 2; why.append(f'div_pmax비 {mp:.2f}>3')
    if te:
        t = te[-1]
        dke = abs(rke[t] - 1)
        if dke > 0.10: level = 2; why.append(f'KE비 {rke[t]:.3f}')
        elif dke > 0.05: level = max(level, 1); why.append(f'KE비 {rke[t]:.3f}')
        if not 0.8 <= rk[t] <= 1.25: level = max(level, 1); why.append(f'k비 {rk[t]:.3f}')
    tag = ['OK', 'WARN', 'FAIL'][level]
    ref_note = ' (Run A 대조 구간 밖)' if not ts else (f' (Run A 기록 끝 t={max(dr):.1f} 이후 — 마지막 대조값)' if t_last > max(dr) + 0.5 else '')
    summ = (f"{tag} t={t_last:.2f}"
            + (f" avg비={med(ra, ts[-1]):.2f} pmax비={med(rp, ts[-1]):.2f}" if ts else '')
            + (f" KE비={rke[te[-1]]:.3f} k비={rk[te[-1]]:.3f}" if te else '')
            + ref_note + (' — ' + ', '.join(why) if why else ''))
    print(summ)
    if a.check:
        sys.exit(level)


if __name__ == '__main__':
    main()
