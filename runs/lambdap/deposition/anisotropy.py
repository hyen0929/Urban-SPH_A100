#!/usr/bin/env python3
"""λp 가 담지 못하는 축 — 바람 정렬(λf/λp)이 플룸 이방성을 설명하는가.

심사에서 반드시 나올 지적: "λp 는 등방 스칼라인데 플룸은 방향성을 갖는다."
λf(정면 면적 밀도)와 λf/λp 는 이미 산출되어 있고, 침적 플룸의 신장도(elong)와
주축 방향(axis_deg)도 있다. 셋을 맞대어 (1) 이방성 지표가 신장도를 설명하는지,
(2) λp 회귀의 잔차를 설명하는지, (3) 플룸 축이 바람과 정렬되는지를 본다.

바람은 +x 방향이므로 axis_deg 의 바람 편차는 min(a, 180-a) 이다.

출력: anisotropy.csv, 화면 요약
"""
import csv, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MORPH = os.path.normpath(os.path.join(HERE, '..'))
COURT = {3, 6, 10}
R_MAIN = 150
rng = np.random.default_rng(20260921)

met = {}
for r in csv.DictReader(open(os.path.join(HERE, 'dep_metrics.csv'))):
    if r['run'] == 'A':
        met[int(r['src'])] = {k: float(v) for k, v in r.items() if k != 'run'}
mo = {}
for r in csv.DictReader(open(os.path.join(MORPH, 'morph_sources.csv'))):
    if int(float(r['R'])) == R_MAIN:
        mo[int(r['id'])] = {k: (float(v) if k not in ('label',) else v) for k, v in r.items()}

IDS = sorted(met)
KEEP = [i for i in IDS if i not in COURT]


def rank(a): return np.argsort(np.argsort(np.asarray(a, float))).astype(float)
def sp(a, b): return float(np.corrcoef(rank(a), rank(b))[0, 1])


def pcorr(a, b, c):
    """c 를 통제한 a–b 편상관 (순위 기준)."""
    ra, rb, rc = rank(a), rank(b), rank(c)
    ea = ra - np.polyval(np.polyfit(rc, ra, 1), rc)
    eb = rb - np.polyval(np.polyfit(rc, rb, 1), rc)
    return float(np.corrcoef(ea, eb)[0, 1])


def perm_p(a, b, n=20000):
    r0 = abs(sp(a, b)); b = np.asarray(b, float); c = 0
    for _ in range(n):
        if abs(sp(a, rng.permutation(b))) >= r0 - 1e-12: c += 1
    return (c+1)/(n+1)


lam_p = [mo[i]['lambda_p'] for i in KEEP]
lam_f = [mo[i]['lambda_f'] for i in KEEP]
ratio = [mo[i]['lambda_f']/mo[i]['lambda_p'] for i in KEEP]
elong = [met[i]['elong'] for i in KEEP]
a90 = [np.log10(met[i]['area90_m2']) for i in KEEP]
trap = [met[i]['f_r150'] for i in KEEP]
shift = [met[i]['shift_m'] for i in KEEP]
# 바람(+x)으로부터의 주축 편차 [deg], 0~90
dev = [min(met[i]['axis_deg'], 180.0-met[i]['axis_deg']) for i in KEEP]

print(f'가로변 방출원 {len(KEEP)}개, R = {R_MAIN} m\n')
print('1) 형태 지표끼리')
print(f'   λf  vs λp        r = {np.corrcoef(lam_f, lam_p)[0,1]:+.2f}')
print(f'   λf/λp vs λp      r = {np.corrcoef(ratio, lam_p)[0,1]:+.2f}   (중복이면 ±1 에 가까움)')

print('\n2) 플룸 이방성을 설명하는가')
for nm, y in [('신장도 elong', elong), ('주축 바람편차', dev)]:
    print(f'   {nm:<14} vs λp    ρs = {sp(lam_p, y):+.2f}  (p = {perm_p(lam_p, y):.3f})')
    print(f'   {"":<14} vs λf/λp ρs = {sp(ratio, y):+.2f}  (p = {perm_p(ratio, y):.3f})')

print('\n3) λp 회귀의 잔차를 λf/λp 가 설명하는가 (편상관)')
for nm, y in [('log A90', a90), ('150 m 갇힘', trap), ('중심이동', shift)]:
    print(f'   {nm:<12} λf/λp | λp  = {pcorr(ratio, y, lam_p):+.2f}      '
          f'(통제 전 {sp(ratio, y):+.2f})')

print('\n4) 소스별 값')
print(f'{"src":>4}{"λp":>7}{"λf":>7}{"λf/λp":>8}{"신장도":>8}{"주축°":>7}{"편차°":>7}{"A90":>9}')
for i in IDS:
    d = min(met[i]['axis_deg'], 180.0-met[i]['axis_deg'])
    mk = '*' if i in COURT else ' '
    print(f'{mk}S{i:<3}{mo[i]["lambda_p"]:7.3f}{mo[i]["lambda_f"]:7.3f}'
          f'{mo[i]["lambda_f"]/mo[i]["lambda_p"]:8.2f}{met[i]["elong"]:8.2f}'
          f'{met[i]["axis_deg"]:7.0f}{d:7.0f}{met[i]["area90_m2"]:9.0f}')

rows = [dict(src=i, lambda_p=mo[i]['lambda_p'], lambda_f=mo[i]['lambda_f'],
             ratio=mo[i]['lambda_f']/mo[i]['lambda_p'], elong=met[i]['elong'],
             axis_deg=met[i]['axis_deg'],
             wind_dev_deg=min(met[i]['axis_deg'], 180.0-met[i]['axis_deg']),
             area90_m2=met[i]['area90_m2'], f_r150=met[i]['f_r150'],
             courtyard=int(i in COURT)) for i in IDS]
with open(os.path.join(HERE, 'anisotropy.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print('\nanisotropy.csv 저장')
