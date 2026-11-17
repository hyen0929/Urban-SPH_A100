#!/usr/bin/env python3
"""fluid_condition 슬래브(y∈[0,2h))의 |u| 분포를 두 런 같은 물리시각에서 대조.
슬래브는 전역 대표성이 없으니(메모: 클립입자 대부분 슬래브 밖) 벌크 유동 A/B 비교용."""
import sys, os, re, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'python'))
from vtk_fields import scan
H = os.path.dirname(os.path.abspath(__file__))
A = os.path.join(H, '..', '..', '2026-09-10_leipzig-full-bcfix_2', 'out')
B = os.path.join(H, '..', 'out')
def stats(p):
    d = dict(scan(p)); u = d['velocity'].reshape(-1, 3); pt = d['p_type']
    m = np.linalg.norm(u[pt == 1], axis=1)
    return np.median(m), np.percentile(m, 99), np.percentile(m, 99.99), m.max(), (m > 10).sum()
print(f"{'t':>4} | {'med_A':>5} {'med_B':>5} | {'p99_A':>5} {'p99_B':>5} | {'p9999A':>6} {'p9999B':>6} | {'max_A':>5} {'max_B':>5} | {'n>10A':>5} {'n>10B':>5}")
for k in range(1, 31):
    fa, fb = f"{A}/fluid_condition_{200*k}stp.vtk", f"{B}/fluid_condition_{100*k}stp.vtk"
    if not (os.path.exists(fa) and os.path.exists(fb)): break
    a, b = stats(fa), stats(fb)
    print(f"{k:4d} | {a[0]:5.2f} {b[0]:5.2f} | {a[1]:5.2f} {b[1]:5.2f} | {a[2]:6.2f} {b[2]:6.2f} | {a[3]:5.2f} {b[3]:5.2f} | {a[4]:5d} {b[4]:5d}")
