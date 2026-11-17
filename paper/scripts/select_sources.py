#!/usr/bin/env python3
"""소스 재선정 — n=13 의 검정력 부족을 해소하되, 무분별하게 늘리지 않는다.

선정 원칙
  1. 기존 13개는 그대로 둔다 (이전 결과와 연결이 끊기지 않게).
  2. λp 를 고르게 덮는다. 상관 검정력은 예측변수의 분산에 비례하므로,
     λp 값이 한쪽에 몰리면 소스를 늘려도 검정력이 안 는다.
  3. 공간적으로 떨어뜨린다. λp(R) 창이 겹치면 두 소스는 같은 형태를 두 번 세는
     셈(pseudo-replication)이라 유효 표본이 늘지 않는다.
  4. 방출점이 열린 곳이어야 하고, 풍하로 뻗을 자리가 있어야 한다.

  python3 select_sources.py [목표개수] [최소간격m]
"""
import csv, os, sys
import numpy as np

HERE  = os.path.dirname(os.path.abspath(__file__))
MORPH = os.path.normpath(os.path.join(HERE, '..'))
N_TARGET = int(sys.argv[1]) if len(sys.argv) > 1 else 40
D_MIN    = float(sys.argv[2]) if len(sys.argv) > 2 else 150.0
R_MAX    = float(sys.argv[3]) if len(sys.argv) > 3 else 1e9   # 회전 중심 (657.5,0) 으로부터 최대 거리
OPEN_MIN = float(sys.argv[4]) if len(sys.argv) > 4 else 0.90  # 방출운(반경 15 m) 안 열린 면적 최소 비율
R_PRIMARY = 100                      # λp 층화에 쓰는 반경
X_INFLOW_MARGIN = 250.0              # 유입면에서 최소 거리
X_DOWNWIND_ROOM = 380.0              # 풍하로 남겨둘 거리 (도메인 x_max=1315)
DOM_X, DOM_Y = 1315.0, 657.5

hf = np.load(os.path.join(MORPH, 'heightfield.npz'))
H, cell, x0, y0 = hf['H'], float(hf['cell']), float(hf['x0']), float(hf['y0'])
ny, nx = H.shape
B = (H > 1.0).astype(np.float64)
gx = x0 + (np.arange(nx)+0.5)*cell
gy = y0 + (np.arange(ny)+0.5)*cell

def lam_field(R):
    """반경 R 원형창 평면적 밀도 (FFT 합성곱)."""
    r = int(round(R/cell)); yy, xx = np.mgrid[-r:r+1, -r:r+1]
    K = ((xx*xx + yy*yy)*cell*cell <= R*R).astype(np.float64); K /= K.sum()
    Kp = np.zeros_like(B); Kp[:K.shape[0], :K.shape[1]] = K
    Kp = np.roll(np.roll(Kp, -(K.shape[0]//2), 0), -(K.shape[1]//2), 1)
    return np.fft.irfft2(np.fft.rfft2(B)*np.fft.rfft2(Kp), B.shape)

LAM = {R: lam_field(R) for R in (75, 100, 150)}

def enclosure(ix, iy, rad=50.0, nray=36):
    """360° 광선 중 rad 안에서 5 m 초과 건물에 막히는 비율. 1 에 가까우면 중정."""
    tall = H > 5.0; hit = 0
    for k in range(nray):
        t = 2*np.pi*k/nray; dx, dy = np.cos(t), np.sin(t)
        for s in np.arange(1, rad/cell):
            qx, qy = int(round(ix+dx*s)), int(round(iy+dy*s))
            if qx < 0 or qx >= nx or qy < 0 or qy >= ny: break
            if tall[qy, qx]: hit += 1; break
    return hit/nray

# ---- 후보 격자 -------------------------------------------------------------
GX, GY = np.meshgrid(gx, gy)
ok  = (H <= 1.0)                                    # 열린 지면
ok &= (GX >= X_INFLOW_MARGIN) & (GX <= DOM_X - X_DOWNWIND_ROOM)
ok &= (np.abs(GY) <= DOM_Y - 150.0)
# 회전 실험에서 네 방향 모두 유효하려면 회전 중심으로부터의 거리를 묶어야 한다.
# 반지름 r 인 점은 회전해도 x 가 [657.5-r, 657.5+r] 안에 머문다.
ok &= (np.hypot(GX - 657.5, GY) <= R_MAX)
# 방출운이 건물 안으로 새지 않아야 한다. 방출은 커널 기각표집으로 반경 kappa*h_src = 15 m
# 까지 퍼지므로, 그 원 안이 충분히 열려 있어야 한다. 건물 안에 떨어진 입자는 IBM 이
# 유동을 붙들어 그 자리에 멈추고 해당 소스의 통계를 망친다. 기존 13 소스는 전부 0.91 이상.
_R=15.0; _r=int(round(_R/cell))
_yy,_xx=np.mgrid[-_r:_r+1,-_r:_r+1]
_K=((_xx*_xx+_yy*_yy)*cell*cell<=_R*_R).astype(np.float64); _K/=_K.sum()
_Kp=np.zeros_like(B); _Kp[:_K.shape[0],:_K.shape[1]]=_K
_Kp=np.roll(np.roll(_Kp,-(_K.shape[0]//2),0),-(_K.shape[1]//2),1)
OPENF=1.0-np.fft.irfft2(np.fft.rfft2(B)*np.fft.rfft2(_Kp), B.shape)
ok &= (OPENF >= OPEN_MIN)
print(f'방출운 열린비율 >= {OPEN_MIN:.2f} 조건 적용')
print(f'후보 셀 {int(ok.sum()):,} 개 (전체 {nx*ny:,})')

# ---- 기존 13개를 먼저 심는다 ----------------------------------------------
S0 = [(int(r['id']), float(r['x']), float(r['y'])) for r in
      csv.DictReader(open(os.path.join(MORPH, 'sources13.csv')))]
sel = [dict(id=i, x=x, y=y, origin='v1') for i, x, y in S0]

def ij(x, y): return int(round((x-x0)/cell)), int(round((y-y0)/cell))
def lam_at(x, y, R):
    i, j = ij(x, y); return float(LAM[R][min(max(j,0),ny-1), min(max(i,0),nx-1)])

for s in sel:
    s['lam100'] = lam_at(s['x'], s['y'], 100)

# ---- λp 층화 + 최대최소거리 탐욕 선택 --------------------------------------
cy, cx = np.nonzero(ok)
cand_x, cand_y = gx[cx], gy[cy]
cand_l = LAM[R_PRIMARY][cy, cx]
lo, hi = np.percentile(LAM[R_PRIMARY][H <= 1.0], [2, 98])
NBIN = 8
edges = np.linspace(lo, hi, NBIN+1)
have = np.histogram([s['lam100'] for s in sel], bins=edges)[0]
print(f'λp({R_PRIMARY}) 후보 범위 {lo:.3f}~{hi:.3f},  기존 13개 구간별 분포 {have.tolist()}')

sx = np.array([s['x'] for s in sel]); sy = np.array([s['y'] for s in sel])
dist = np.hypot(cand_x[:, None]-sx[None, :], cand_y[:, None]-sy[None, :]).min(1)
nid = 14
while len(sel) < N_TARGET:
    need = np.argsort(have)                       # 가장 비어 있는 구간부터
    placed = False
    for b in need:
        m = (cand_l >= edges[b]) & (cand_l < edges[b+1]) & (dist >= D_MIN)
        if not m.any(): continue
        k = np.nonzero(m)[0][np.argmax(dist[m])]  # 그 구간에서 가장 멀리 떨어진 점
        sel.append(dict(id=nid, x=float(cand_x[k]), y=float(cand_y[k]),
                        origin='v2', lam100=float(cand_l[k])))
        have[b] += 1; nid += 1; placed = True
        dist = np.minimum(dist, np.hypot(cand_x-cand_x[k], cand_y-cand_y[k]))
        break
    if not placed:
        print(f'!! 최소간격 {D_MIN:.0f} m 로는 {len(sel)}개가 한계'); break

# ---- 정리 ------------------------------------------------------------------
for s in sel:
    i, j = ij(s['x'], s['y'])
    for R in (75, 100, 150): s[f'lam{R}'] = lam_at(s['x'], s['y'], R)
    s['enclosure'] = round(enclosure(i, j), 3)
P = np.array([[s['x'], s['y']] for s in sel])
DD = np.hypot(P[:, 0][:, None]-P[:, 0][None, :], P[:, 1][:, None]-P[:, 1][None, :])
np.fill_diagonal(DD, 1e9)
for k, s in enumerate(sel): s['nn_m'] = round(float(DD[k].min()), 1)

L = np.array([s['lam150'] for s in sel])
print(f'\n선정 {len(sel)}개  (기존 13 + 신규 {len(sel)-13})')
print(f'  λp(150) {L.min():.3f}~{L.max():.3f}  표준편차 {L.std():.3f}  (기존 13개만: '
      f'{np.array([s["lam150"] for s in sel[:13]]).std():.3f})')
print(f'  최근접 간격 최소 {min(s["nn_m"] for s in sel):.0f} m  중앙값 {np.median([s["nn_m"] for s in sel]):.0f} m')
print(f'  중정형(enclosure>0.8) {sum(1 for s in sel if s["enclosure"]>0.8)}개')
cols = ['id','x','y','lam75','lam100','lam150','enclosure','nn_m','origin']
with open(os.path.join(HERE,'sources_v2.csv'),'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=cols); w.writeheader()
    for s in sel: w.writerow({c: s[c] for c in cols})
print('저장: sources_v2.csv')
