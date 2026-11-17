#!/usr/bin/env python3
"""dt=0.01 시험런 vs Run A(dt=0.005) PPE_diag 를 같은 물리시각에서 대조.

PPE_diag 의 div_un1_* 와 bi_max 는 rho*div(u)/dt 로 기록된다(ISPH_Calc.cuh:499).
dt 가 2배면 같은 속도발산이 raw 값으로는 절반이 되므로, 비교는 값*dt
(= rho*div u, dt 무관) 로 한다.

사용: compare.py [--every T] [--check]
  --check : 판정기준 위반(div_avg 비 > 2, plain_max 비 > 3, NaN)이 있으면 exit 1
"""
import csv, math, sys, os, argparse

HERE = os.path.dirname(os.path.abspath(__file__))
TRIAL = os.path.join(HERE, '..', 'out', 'PPE_diag_explicit.csv')
REF = os.path.join(HERE, '..', '..', '2026-09-10_leipzig-full-bcfix_2', 'out', 'PPE_diag_explicit.csv')
DT_TRIAL, DT_REF = 0.01, 0.005
LIM_AVG, LIM_PMAX = 2.0, 3.0


def load(path):
    rows = {}
    with open(path) as f:
        for r in csv.DictReader(f):
            try:
                t = round(float(r['time']), 2)
                rows[t] = {k: float(v) for k, v in r.items()}
            except (ValueError, TypeError):
                continue
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--every', type=float, default=1.0)
    ap.add_argument('--check', action='store_true')
    a = ap.parse_args()

    tr, rf = load(TRIAL), load(REF)
    ts = sorted(t for t in tr if t in rf)
    bad = []
    print(f"{'t[s]':>6} {'avg_A':>9} {'avg_B':>9} {'ratio':>6} | {'pmax_A':>9} {'pmax_B':>9} {'ratio':>6} | {'bimax_r':>7}")
    last_print = -1e9
    for t in ts:
        A, B = rf[t], tr[t]
        avgA, avgB = A['div_un1_avg'] * DT_REF, B['div_un1_avg'] * DT_TRIAL
        pmA, pmB = A['div_un1_plain_max'] * DT_REF, B['div_un1_plain_max'] * DT_TRIAL
        biA, biB = A['bi_max'] * DT_REF, B['bi_max'] * DT_TRIAL
        ra, rp, rb = avgB / avgA, pmB / pmA, biB / biA
        viol = (not all(map(math.isfinite, (avgB, pmB, biB)))) or ra > LIM_AVG or rp > LIM_PMAX
        if viol:
            bad.append((t, ra, rp))
        if t - last_print >= a.every - 1e-9 or viol or t == ts[-1]:
            print(f"{t:6.2f} {avgA:9.3e} {avgB:9.3e} {ra:6.2f} | {pmA:9.3e} {pmB:9.3e} {rp:6.2f} | {rb:7.2f}{'  <-- 위반' if viol else ''}")
            last_print = t
    if not ts:
        print('겹치는 시각 없음')
    else:
        print(f"\n대조 {len(ts)}점, t=0~{ts[-1]:.2f} s (값은 rho*div u = raw*dt, A=Run A dt0.005, B=시험 dt0.01)")
        print(f"위반 {len(bad)}점" + (f": 첫 위반 t={bad[0][0]:.2f} (avg비 {bad[0][1]:.2f}, pmax비 {bad[0][2]:.2f})" if bad else ''))
    if a.check and bad:
        sys.exit(1)


if __name__ == '__main__':
    main()
