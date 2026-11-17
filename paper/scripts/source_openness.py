#!/usr/bin/env python3
"""소스별 방향 의존 트임 — 각도마다 λp 상관이 달라지는 이유를 형상만으로 검증.

도시를 +θ 회전시키고 바람을 +x 로 두는 것 = 도시 고정, 바람을 −θ 로 돌리는 것.
  rot 0°→바람 +x(0°)   rot 90°→−y(270°)   rot 180°→−x(180°)   rot 270°→+y(90°)
소스는 도시와 함께 돌므로 도시 좌표계에서 위치가 고정이다. 바람만 바뀐다.

소스마다 ±30° 부채꼴로 광선을 쏘아 건물에 닿기까지의 거리를 평균한다.
  down : 풍하 트임 (플룸이 뻗을 수 있는 거리)
  up   : 풍상 페치 (유입류가 얼마나 방해 없이 오는가)
'건물'은 방출고 5 m 보다 높은 것만 센다(낮은 건물은 5 m 방출을 막지 못한다).
"""
import csv, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MORPH = os.path.normpath(os.path.join(HERE, '..'))
hf = np.load(os.path.join(MORPH, 'heightfield.npz'))
H, cell, x0, y0 = hf['H'], float(hf['cell']), float(hf['x0']), float(hf['y0'])
B = H > 5.0                                  # 방출고보다 높은 건물만 장애물로 본다
ny, nx = H.shape
S = {int(r['id']): (float(r['x']), float(r['y'])) for r in
     csv.DictReader(open(os.path.join(MORPH, 'sources13.csv')))}
lam = {}
for r in csv.DictReader(open(os.path.join(MORPH, 'lambda_p_by_R.csv'))):
    lam[int(float(r['R_m']))] = {i: float(r[f'S{i}']) for i in range(1, 14)}
COURT = {3, 6, 10}
ROT2WIND = {0: 0.0, 90: 270.0, 180: 180.0, 270: 90.0}   # 회전각 → 도시좌표계 바람 방위
MAXR, NS, SECT = 400.0, 320, 30.0

def ray(sx, sy, deg):
    """(sx,sy) 에서 deg 방향으로 건물에 닿기까지 [m]. 도메인을 벗어나면 MAXR."""
    t = np.radians(deg); dx, dy = np.cos(t), np.sin(t)
    px, py = (sx-x0)/cell, (sy-y0)/cell
    for k in range(1, NS+1):
        s = k*(MAXR/NS)/cell
        qx, qy = int(round(px+dx*s)), int(round(py+dy*s))
        if qx < 0 or qx >= nx or qy < 0 or qy >= ny: return MAXR
        if B[qy, qx]: return k*(MAXR/NS)
    return MAXR

def sector(sx, sy, deg):
    ds = np.arange(-SECT, SECT+1e-9, 5.0)
    return float(np.mean([ray(sx, sy, deg+d) for d in ds]))

ids = sorted(S); keep = np.array([i not in COURT for i in ids])
def sp(a, b):
    return float(np.corrcoef(np.argsort(np.argsort(a)), np.argsort(np.argsort(b)))[0, 1])
x150 = np.array([lam[150][i] for i in ids])

rows = []
print('소스별 풍하 트임 down [m]  (±30° 부채꼴, 장애물 = 높이 5 m 초과 건물)')
print(' src  λp150 | rot0(+x)  rot90(-y)  rot180(-x)  rot270(+y)')
tab = {}
for i in ids:
    sx, sy = S[i]
    d = {rot: sector(sx, sy, w) for rot, w in ROT2WIND.items()}
    u = {rot: sector(sx, sy, (w+180) % 360) for rot, w in ROT2WIND.items()}
    tab[i] = (d, u)
    print(f'  {i:2d} {lam[150][i]:6.3f} | ' + '  '.join(f'{d[r]:8.1f}' for r in (0,90,180,270)))
    for r in (0,90,180,270):
        rows.append(dict(src=i, rot=r, wind_deg=ROT2WIND[r], lambda_p150=lam[150][i],
                         down_m=d[r], up_m=u[r]))

print('\nλp(150) 와 방향별 트임의 순위상관 (중정 제외 n=10)')
print(' 회전각   바람      ρs(λp, 풍하트임)   ρs(λp, 풍상페치)')
for r in (0,90,180,270):
    dv = np.array([tab[i][0][r] for i in ids]); uv = np.array([tab[i][1][r] for i in ids])
    print(f'  {r:4d}°  {ROT2WIND[r]:5.0f}°       {sp(x150[keep], dv[keep]):+.2f}              {sp(x150[keep], uv[keep]):+.2f}')

# 실측 지표와 맞물리는지
mp = {}
for r in csv.DictReader(open(os.path.join(HERE, 'rot_metrics.csv'))):
    mp.setdefault(int(r['deg']), {})[int(r['src'])] = r
print('\n실측 지표와 풍하 트임의 순위상관 (중정 제외 n=10)')
print(' 회전각   ρs(트임, A90)   ρs(트임, 갇힘)   [참고] ρs(λp, A90)')
for r in sorted(mp):
    dv = np.array([tab[i][0][r] for i in ids])
    a90 = np.array([float(mp[r][i]['area90_m2']) for i in ids])
    f150 = np.array([float(mp[r][i]['f_r150']) for i in ids])
    print(f'  {r:4d}°      {sp(dv[keep], a90[keep]):+.2f}            {sp(dv[keep], f150[keep]):+.2f}            {sp(x150[keep], a90[keep]):+.2f}')

with open(os.path.join(HERE, 'source_openness.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print('\n저장: source_openness.csv')
