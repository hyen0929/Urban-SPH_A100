#!/usr/bin/env python3
"""벽면 침적이 있었다면 λp–침적률 관계가 달라졌을까 — 노출 적분으로 추정.

가설(사용자): 모델이 평지 z=0 에만 침적시키므로 λp 가 큰 소스는 가용 지면이 좁아 덜 침적한다.
벽면 침적을 넣으면 λp 가 큰 곳에서 벽 면적이 크므로 관계가 바뀔 것이다.

방법 — 새 계산 없이 노출 적분으로 추정한다.
  1. 높이장으로 각 고도층의 건물 마스크 B_k 를 만들고, B_k 에 4-이웃 팽창을 걸어 **벽에 접한 개방 셀**을
     찾는다(격자 2.5 m 정방이라 셀 중심에서 인접 셀의 벽면까지 1.25 m < 2 m). scipy 거리변환 불필요.
  2. LDM 덤프(1 s 간격)마다 활성 입자를 두 무리로 센다.
       지면 후보 : z < DD_h_dep(2 m)                  → 실제 모델이 침적시키는 무리
       벽면 후보 : z >= 2 m 이고 벽까지 거리 < 2 m     → 모델이 무시하는 무리
  3. DD 커널과 같은 식으로 덤프마다 기대 침적수를 더한다: Σ [1-exp(-Vd(z)·Δt/dz)].
     지면은 dz = 입자의 지면 위 높이(커널과 동일), 벽면은 dz = 셀 중심에서 벽까지 1.25 m.
     지면 Vd = 1/(Ra+Rb) + Vs, 벽면 Vd = 1/(Ra+Rb) (중력침강은 수직면에 작용하지 않음).
  4. 지면 추정치를 실제 침적수와 비교해 추정기 자체를 검증한다.
  5. (지면 실측 + 벽면 추정) 으로 침적률을 다시 만들고 λp 상관을 다시 계산한다.

한계 — 벽면 추정은 **상한**이다.
  (a) 지면 침적은 입자를 제거하므로 남은 활성 입자의 지면 노출이 고갈되지만, 벽면은 모델이 흡수하지
      않으므로 같은 입자가 계속 벽 옆에 남아 반복해서 카운트된다.
  (b) 벽면 Ra 를 지면과 같게 두었다(벽면 경계층은 다르다).
  (c) 1 s 간격 표본이라 그 사이의 짧은 접근을 놓친다.
크기 규모(order of magnitude)로만 읽을 것.

출력: wall_exposure.csv, 화면 요약
"""
import csv, os, re
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.normpath(os.path.join(HERE, '..', '..'))
A = os.path.join(RUNS, '2026-10-20_leipzig-s13-most-wall')
MORPH = os.path.normpath(os.path.join(HERE, '..'))
EVERY, DT_DUMP = 100, 1.0      # step 간격 100 = 1 s
H_DEP, WALL_D = 2.0, 1.25      # 지면 트리거 높이 [m], 셀 중심에서 인접 벽면까지 거리 [m]
COURT = {3, 6, 10}
# --- 모델 상수 (functions_DD.cuh / Parameters.cuh) ---
USTAR, Z0, KAPPA = 0.5, 0.1, 0.4
NU, DB, TAU_P, EPS0, VS = 1.5e-5, 2.8e-11, 9.0e-6, 3.0, 8.9e-5

ns = {'__file__': os.path.join(HERE, 'dep_metrics.py')}
exec(open(os.path.join(HERE, 'dep_metrics.py')).read().split('def load(')[0], ns)
read_vtk = ns['read_vtk']


def Vd(z, settling=True):
    """중립(OB=1e6)에서 DD_calc_Vd 재현. z 는 배열도 가능."""
    z_eff = np.maximum(z, 2*Z0)
    Ra = np.maximum(np.log(z_eff/Z0)/(KAPPA*USTAR), 1.0)
    Sc = NU/DB; E_B = Sc**(-2/3)
    St = TAU_P*USTAR**2/NU; E_IM = (St/(St+0.5))**2
    Rb = 1.0/(USTAR*EPS0*(E_B+E_IM))
    return 1.0/(Ra+Rb) + (VS if settling else 0.0)


def p_dep(z, dz, settling=True, dt=1.0):
    """DD_calc_P_dep 그대로: P = 1 - exp(-Vd*dt/dz). 입자 높이 z 를 쓴다."""
    return 1.0 - np.exp(-Vd(z, settling)*dt/np.maximum(dz, 1e-4))


def main():
    hf = np.load(os.path.join(MORPH, 'heightfield.npz'))
    H, x0, y0, cell = hf['H'], float(hf['x0']), float(hf['y0']), float(hf['cell'])
    zlev = np.arange(1.25, 80.0, 2.5)                       # 입자 고도 bin 중심
    wallnb = np.zeros((len(zlev),) + H.shape, dtype=bool)    # 벽에 접한 개방 셀
    for k, z in enumerate(zlev):
        B = H > z                                            # 이 고도에서 건물인 셀
        d = np.zeros_like(B)
        d[1:, :] |= B[:-1, :]; d[:-1, :] |= B[1:, :]
        d[:, 1:] |= B[:, :-1]; d[:, :-1] |= B[:, 1:]
        wallnb[k] = d & ~B
    print(f'[준비] 고도층 {len(zlev)}개 벽 인접셀 마스크 완료 (격자 {H.shape}, cell {cell} m, '
          f'z=1.25 m 층 벽인접 {int(wallnb[0].sum()):,} 셀)')

    S = {int(r['id']): (float(r['x']), float(r['y'])) for r in
         csv.DictReader(open(os.path.join(MORPH, 'sources13.csv')))}
    steps = sorted(int(re.search(r'_(\d+)stp', f).group(1))
                   for f in os.listdir(os.path.join(A, 'out')) if f.startswith('LDM_'))
    steps = [s for s in steps if s >= 15000 and s % EVERY == 0]

    Eg = {i: 0.0 for i in S}; Ew = {i: 0.0 for i in S}; Ndep = {i: 0 for i in S}; Nrel = {i: 0 for i in S}
    for s in steps:
        p, F = read_vtk(os.path.join(A, 'out', f'LDM_{s}stp.vtk'))
        src = F['p_type'].astype(int) - 100
        dep = np.linalg.norm(F['velocity'], axis=1) == 0.0
        air = ~dep
        ix = np.clip(np.floor((p[:, 0]-x0)/cell).astype(int), 0, H.shape[1]-1)
        iy = np.clip(np.floor((p[:, 1]-y0)/cell).astype(int), 0, H.shape[0]-1)
        kz = np.clip(np.searchsorted(zlev, p[:, 2]) - 0, 0, len(zlev)-1)
        inside = p[:, 2] <= H[iy, ix]                       # 건물 내부(있을 수 없지만 방어)
        atwall = wallnb[kz, iy, ix]
        near_ground = air & (p[:, 2] < H_DEP)
        near_wall = air & ~inside & (p[:, 2] >= H_DEP) & (p[:, 2] < zlev[-1]) & atwall
        # 노출이 아니라 **침적 확률 합**을 누적한다 — DD 커널과 같은 식, 같은 dz 규약.
        pg = np.zeros(len(p)); pw = np.zeros(len(p))
        if near_ground.any():
            zz = p[near_ground, 2]
            pg[near_ground] = p_dep(zz, zz, True, DT_DUMP)      # dz = 지면 위 높이
        if near_wall.any():
            pw[near_wall] = p_dep(p[near_wall, 2], WALL_D, False, DT_DUMP)  # dz = 벽까지 거리
        for i in S:
            m = src == i
            Eg[i] += float(pg[m].sum())
            Ew[i] += float(pw[m].sum())
        if s == steps[-1]:
            for i in S:
                m = src == i; Ndep[i] = int((m & dep).sum()); Nrel[i] = int(m.sum())
    vg, vw = Vd(H_DEP, True), Vd(H_DEP, False)
    print(f'[Vd] 지면 {vg:.3e} m/s, 벽면 {vw:.3e} m/s (Vs={VS:.1e}, 비 {vw/vg:.3f})')

    lam = {}
    for r in csv.DictReader(open(os.path.join(MORPH, 'lambda_p_by_R.csv'))):
        lam[int(float(r['R_m']))] = {i: float(r[f'S{i}']) for i in range(1, 14)}
    rows = []
    for i in sorted(S):
        ng_est = Eg[i]   # 이미 기대 침적수 [개]
        nw_est = Ew[i]
        rows.append(dict(src=i, lambda_p150=lam[150][i], n_released=Nrel[i], n_dep_obs=Ndep[i],
                         E_ground=Eg[i], E_wall=Ew[i], n_ground_est=ng_est, n_wall_est=nw_est,
                         est_over_obs=ng_est/max(Ndep[i], 1), wall_over_ground=nw_est/max(ng_est, 1e-9),
                         f_dep_obs=Ndep[i]/max(Nrel[i], 1),
                         f_dep_with_wall=(Ndep[i] + nw_est*Ndep[i]/max(ng_est, 1e-9))/max(Nrel[i], 1)))
    with open(os.path.join(HERE, 'wall_exposure.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

    def sp(a, b):
        return float(np.corrcoef(np.argsort(np.argsort(a)), np.argsort(np.argsort(b)))[0, 1])
    keep = [r for r in rows if r['src'] not in COURT]
    x = np.array([r['lambda_p150'] for r in keep])
    print(f"\n{'src':>4} {'λp150':>7} {'노출 지면':>10} {'노출 벽면':>10} {'벽/지면':>8} "
          f"{'추정/실측(지면)':>14} {'f_dep':>7} {'벽 포함':>8}")
    for r in rows:
        print(f"{r['src']:>4} {r['lambda_p150']:7.3f} {r['E_ground']:10.3e} {r['E_wall']:10.3e} "
              f"{r['wall_over_ground']:8.2f} {r['est_over_obs']:14.2f} {r['f_dep_obs']:7.3f} {r['f_dep_with_wall']:8.3f}")
    print(f"\nλp 상관 (중정 제외 10점): 지면만 {sp(x, np.array([r['f_dep_obs'] for r in keep])):+.2f} "
          f"→ 벽면 포함 {sp(x, np.array([r['f_dep_with_wall'] for r in keep])):+.2f}")
    print(f"벽/지면 비 대 λp: {sp(x, np.array([r['wall_over_ground'] for r in keep])):+.2f}")
    print(f"지면 추정기 검증: 추정/실측 중앙 {np.median([r['est_over_obs'] for r in rows]):.2f} "
          f"(1.0 이면 노출적분+Vd 만으로 실제 침적을 재현)")


if __name__ == '__main__':
    main()
