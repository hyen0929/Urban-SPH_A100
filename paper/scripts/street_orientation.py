#!/usr/bin/env python3
"""가로망 방향성 — 회전 실험 결과(0°·180° 는 닮고 90° 만 다름)를 설명할 수 있는가.

도시를 +θ 회전시키고 바람을 +x 로 두는 것은, 도시를 고정하고 바람을 −θ 로 돌리는 것과 같다.
  rot 0°   → 도시 좌표계에서 바람 +x
  rot 90°  → 바람 −y
  rot 180° → 바람 −x
  rot 270° → 바람 +y
따라서 0°/180° 는 x 축 바람, 90°/270° 는 y 축 바람이다. 두 지표가 x 와 y 에서 다르면
"0/180 이 닮고 90/270 이 다르다"가 예측된다.

지표 1  λf(θ) : 정면적 밀도. 벽면을 바람 직각면에 투영한 면적 / 평면적.
지표 2  L(θ)  : 열린 공간의 평균 자유행로. 임의의 열린 점에서 θ 방향으로 건물에 닿기까지 거리.
"""
import numpy as np, csv, os
HERE = os.path.dirname(os.path.abspath(__file__))
hf = np.load('/home/kdh/Urban-SPH/runs/lambdap/heightfield.npz')
H, cell = hf['H'], float(hf['cell'])          # H[y, x]
B = H > 1.0
ny, nx = H.shape
A_plan = ny*nx*cell*cell

# --- 지표 1: 정면적 ---------------------------------------------------------
Ax = np.abs(np.diff(H, axis=1)).sum()*cell    # x 에 수직인 벽 (동/서향)
Ay = np.abs(np.diff(H, axis=0)).sum()*cell    # y 에 수직인 벽 (남/북향)
print(f'벽면적   x 수직 {Ax:12,.0f} m²    y 수직 {Ay:12,.0f} m²    비 Ax/Ay = {Ax/Ay:.3f}')
print(f'λf       바람 x {Ax/A_plan:.4f}        바람 y {Ay/A_plan:.4f}')

# --- 지표 2: 평균 자유행로 (광선 추적) ---------------------------------------
rng = np.random.default_rng(0)
oy, ox = np.nonzero(~B)
pick = rng.choice(len(oy), size=min(20000, len(oy)), replace=False)
py, px = oy[pick].astype(float), ox[pick].astype(float)
MAXR, NS = 400.0, 160                          # 최대 400 m, 2.5 m 간격
steps = (np.arange(1, NS+1)*(MAXR/NS))/cell

rows = []
for deg in range(0, 180, 5):
    t = np.radians(deg); dx, dy = np.cos(t), np.sin(t)
    hit = np.full(len(px), MAXR)
    done = np.zeros(len(px), bool)
    for s in steps:
        qx = np.rint(px + dx*s).astype(int); qy = np.rint(py + dy*s).astype(int)
        out = (qx < 0)|(qx >= nx)|(qy < 0)|(qy >= ny)
        qx = np.clip(qx, 0, nx-1); qy = np.clip(qy, 0, ny-1)
        newhit = (~done) & (B[qy, qx] | out)
        hit[newhit] = s*cell; done |= newhit
        if done.all(): break
    lf = (Ax*abs(np.cos(t)) + Ay*abs(np.sin(t)))/A_plan
    rows.append(dict(deg=deg, L_mean_m=float(hit.mean()), L_med_m=float(np.median(hit)), lambda_f=float(lf)))

L = np.array([r['L_mean_m'] for r in rows]); D = np.array([r['deg'] for r in rows])
print(f'\n자유행로  최대 {L.max():.1f} m @ {D[L.argmax()]}°    최소 {L.min():.1f} m @ {D[L.argmin()]}°    비 {L.max()/L.min():.2f}')
print('\n  바람(도시좌표)   회전실험      λf       평균자유행로')
for deg, rotname in [(0,'rot 0° / 180°'), (90,'rot 90° / 270°')]:
    r = rows[deg//5]
    print(f'   {deg:3d}°          {rotname:14s} {r["lambda_f"]:.4f}   {r["L_mean_m"]:7.1f} m')
print('\n전체 스윕 (5° 간격)')
for r in rows: print(f'  {r["deg"]:3d}°  λf {r["lambda_f"]:.4f}   L {r["L_mean_m"]:6.1f} m')
with open(os.path.join(HERE,'street_orientation.csv'),'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
