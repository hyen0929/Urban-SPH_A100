#!/usr/bin/env python3
"""LDM 출력 → CEDVAL 측정점별·퍼프별 농도 시계열.

측정점 주변 상자(수평 BOX×BOX, 연직 [z-DZ, z+DZ], 지면 측정점은 [0, 2·DZ])에 든
퍼프 입자 수를 센다. 퍼프 하나의 전체 질량을 1 로 두면 C = n / (N_puff · V)  [1/m³].
dosage = Σ C Δt  [s/m³].  무차원화는 비교 스크립트에서 한다.

출력: <run>/analysis/receptor_ts.npz
  t[nt], rec[nr,3] (CEDVAL 좌표), C[nr, npuff, nt], t_release[npuff]
사용: python3 validation/cedval/receptor_ts.py runs/<run> [BOX=5] [DZ=2]
"""
import glob, os, re, sys
import numpy as np

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
RUN = os.path.join(ROOT, sys.argv[1]) if len(sys.argv) > 1 else None
BOX = float(sys.argv[2]) if len(sys.argv) > 2 else 5.0
DZ = float(sys.argv[3]) if len(sys.argv) > 3 else 2.0
X_SHIFT = 400.0
N_PUFF, NP_PUFF, T_FIRST, PUFF_EVERY = 60, 20000, 200.0, 5.0

ns = {'__file__': f'{ROOT}/runs/lambdap/deposition/dep_metrics.py'}
exec(open(ns['__file__']).read().split('def load(')[0], ns)
read_vtk = ns['read_vtk']

REC = np.array([
    (198, -135, 2.02), (198, -90, 2.02), (198, -36, 2.02), (198, 0, 2.02), (198, 36, 2.02), (198, 90, 2.02), (198, 135, 2.02),
    (-54, -126, 2.02), (-54, -90, 2.02), (-54, -54, 2.02), (-54, -18, 2.02), (-54, 0, 2.02), (-54, 18, 2.02), (-54, 36, 2.02),
    (-54, 54, 2.02), (-54, 90, 2.02), (-54, 126, 2.02),
    (-162, -90, 2.02), (-162, -36, 2.02), (-162, 36, 2.02), (-162, 90, 2.02),
    (-234, -54, 2.02), (-234, -36, 2.02), (-234, 36, 2.02), (-234, 54, 2.02),
    (-54, 0, 4.5), (-54, 0, 9), (-54, 0, 13.5), (-54, 0, 18), (-54, 0, 27), (-54, 0, 36), (-54, 0, 54), (-54, 0, 72),
    (198, 0, 9), (198, 0, 18), (198, 0, 27), (198, 0, 36), (198, 0, 54), (198, 0, 90),
    (-162, 0, 18), (-162, 0, 27), (-162, 0, 54), (-234, 0, 18), (-234, 0, 27),
], float)


def main():
    files = sorted(glob.glob(f'{RUN}/out/LDM_*stp.vtk'), key=lambda f: int(re.search(r'_(\d+)stp', f).group(1)))
    steps = np.array([int(re.search(r'_(\d+)stp', f).group(1)) for f in files])
    t = steps * 0.01
    lo = REC.copy(); hi = REC.copy()
    lo[:, :2] -= BOX / 2; hi[:, :2] += BOX / 2
    g = REC[:, 2] < 3
    lo[:, 2] = np.where(g, 0.0, REC[:, 2] - DZ); hi[:, 2] = np.where(g, 2 * DZ, REC[:, 2] + DZ)
    lo[:, 0] += X_SHIFT; hi[:, 0] += X_SHIFT
    V = BOX * BOX * (hi[:, 2] - lo[:, 2])
    C = np.zeros((len(REC), N_PUFF, len(files)), np.float32)
    for k, f in enumerate(files):
        p, F = read_vtk(f)
        pid = F['p_type'].astype(int) - 101
        for r in range(len(REC)):
            m = np.all((p >= lo[r]) & (p < hi[r]), axis=1)
            if m.any():
                C[r, :, k] = np.bincount(pid[m], minlength=N_PUFF)[:N_PUFF] / (NP_PUFF * V[r])
        if k % 50 == 0:
            print(f'{k}/{len(files)} t={t[k]:.0f}', flush=True)
    os.makedirs(f'{RUN}/analysis', exist_ok=True)
    np.savez_compressed(f'{RUN}/analysis/receptor_ts.npz', t=t, rec=REC, C=C,
                        t_release=T_FIRST + PUFF_EVERY * np.arange(N_PUFF), box=BOX, dz=DZ)
    print('저장', f'{RUN}/analysis/receptor_ts.npz')


if __name__ == '__main__':
    main()
