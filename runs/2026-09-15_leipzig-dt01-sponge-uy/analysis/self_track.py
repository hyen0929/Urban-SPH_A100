#!/usr/bin/env python3
"""Run A(정전 전) 기록이 끝난 t>237 s 구간 감시 — 이 런 자신의 t=220~237 s 를 기준으로 삼는다.
  FAIL  div_avg(raw) 2 s 중앙값 > 기준 중앙값 x2, 또는 Σ½|u|² 가 기준 평균 대비 ±10 % 밖, 또는 NaN
  WARN  Σ½|u|² ±5 % 밖, 또는 Σk_turb 가 기준 평균 대비 ±25 % 밖
사용: self_track.py  (마지막 줄 "OK|WARN|FAIL ...", exit 0/1/2)
"""
import csv, math, os, re, statistics, sys
RUN = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); DT = 0.01
T0, T1 = 220.0, 237.0
SUMK = re.compile(rb"^(\d+)\s+tu_max=(?:.*\n)*?.*?\xce\xa3k_turb=([\d.eE+-]+|nan|inf)\s+\xce\xa3\xc2\xbd\|u\|\xc2\xb2=([\d.eE+-]+|nan|inf)", re.M)
div = {}
for r in csv.DictReader(open(os.path.join(RUN, 'out', 'PPE_diag_explicit.csv'))):
    try: div[round(int(r['step'])*DT, 2)] = float(r['div_un1_avg'])
    except (ValueError, KeyError): pass
en = {round(int(m.group(1))*DT, 2): (float(m.group(2)), float(m.group(3))) for m in SUMK.finditer(open(os.path.join(RUN, 'run.log'), 'rb').read())}
t_last = max(div)
base_div = statistics.median([v for t, v in div.items() if T0 <= t <= T1])
base_k = statistics.mean([v[0] for t, v in en.items() if T0 <= t <= T1])
base_ke = statistics.mean([v[1] for t, v in en.items() if T0 <= t <= T1])
cur_div = statistics.median([v for t, v in div.items() if t_last - 2 < t <= t_last])
te = max(en); k, ke = en[te]
level, why = 0, []
if any(not math.isfinite(v) for v in div.values()) or not (math.isfinite(k) and math.isfinite(ke)):
    level = 2; why.append('NaN/Inf')
rd, rk, rke = cur_div/base_div, k/base_k, ke/base_ke
if rd > 2: level = 2; why.append(f'div_avg x{rd:.2f}')
if abs(rke-1) > 0.10: level = 2; why.append(f'KE x{rke:.3f}')
elif abs(rke-1) > 0.05: level = max(level, 1); why.append(f'KE x{rke:.3f}')
if abs(rk-1) > 0.25: level = max(level, 1); why.append(f'k x{rk:.3f}')
print(f"{['OK','WARN','FAIL'][level]} t={t_last:.2f} (기준 t={T0:.0f}~{T1:.0f} 자기값) div_avg x{rd:.2f} KE x{rke:.3f} k x{rk:.3f}" + (' — '+', '.join(why) if why else ''))
sys.exit(level)
