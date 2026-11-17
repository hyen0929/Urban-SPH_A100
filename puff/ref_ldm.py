#!/usr/bin/env python3
"""50소스 평균장 LDM 런을 퍼프 모델 비교 기준으로 정리한다.

각 런 좌표계 그대로(바람 +x). 침적 입자(속도 0)·유출 입자(x>1315)는 공기 중에서 뺀다.
농도는 c/Q [s/m^3] (Q = 소스당 방출률).

  python3 puff/ref_ldm.py 0      → puff/data/ref_rot0.npz
"""
import os, sys, glob
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ldm_io import read_ldm

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUN = {0: 'runs/2026-10-06_leipzig-src50-rot0', 90: 'runs/2026-10-06_leipzig-src50-rot90',
       180: 'runs/2026-10-06_leipzig-src50-rot180', 270: 'runs/2026-10-06_leipzig-src50-rot270'}
INP = 'input/Leipzig_src50_rot{0}/Leipzig_src50_rot{0}_LDM.txt'
DT = 0.01
BIN = 5.0
X0, Y0, NXB, NYB = 0.0, -657.5, 263, 263
ZBANDS = [(2.5, 7.5), (17.5, 22.5)]               # 소스 높이(5 m)와 20 m
SNAP_T = [200.0, 250.0, 300.0]
AVG_T = (250.0, 300.0)


def sources(rot):
    a = np.loadtxt(os.path.join(ROOT, INP.format(rot)), skiprows=1)
    ids = np.unique(a[:, 5]).astype(int)
    S = []
    for s in ids:
        b = a[a[:, 5] == s]
        S.append((s, b[0, 0], b[0, 1], b[0, 2], len(b), b[0, 4], b[:, 10].min(), b[:, 10].max()))
    return np.array(S)    # id, x, y, z, n, mass, t_start, t_end


def one(args):
    path, sid = args
    d = read_ldm(path)
    xyz, pt = d['xyz'], d['p_type']
    dep = np.abs(d['velocity']).sum(1) == 0
    out = xyz[:, 0] > 1315.0
    ndep = np.bincount(np.searchsorted(sid, pt[dep]), minlength=len(sid))
    nrel = np.bincount(np.searchsorted(sid, pt), minlength=len(sid))
    air = ~dep & ~out
    H = np.zeros((len(ZBANDS), len(sid), NXB, NYB), np.float32)
    ix = ((xyz[:, 0] - X0) / BIN).astype(int); iy = ((xyz[:, 1] - Y0) / BIN).astype(int)
    ok = air & (ix >= 0) & (ix < NXB) & (iy >= 0) & (iy < NYB)
    si = np.searchsorted(sid, pt)
    for b, (z0, z1) in enumerate(ZBANDS):
        m = ok & (xyz[:, 2] >= z0) & (xyz[:, 2] < z1)
        np.add.at(H[b], (si[m], ix[m], iy[m]), 1.0)
    return ndep, nrel, H


def main(rot):
    S = sources(rot); sid = S[:, 0].astype(int)
    n_src, mass, t0, t1 = S[0, 4], S[0, 5], S[0, 6], S[0, 7] + 3e-3
    Q = n_src * mass / (t1 - t0)                             # 소스당 방출률 [mass/s]
    files = {int(os.path.basename(f)[4:-7]): f for f in glob.glob(os.path.join(ROOT, RUN[rot], 'out', 'LDM_*stp.vtk'))}
    steps = sorted(s for s in files if s * DT >= t0 - 1e-6)
    with Pool(24) as p:
        R = p.map(one, [(files[s], sid) for s in steps])
    T = np.array(steps) * DT
    dep_frac = np.array([r[0] / n_src for r in R])          # (nt, nsrc) 방출 총량 대비 침적 비율
    rel_frac = np.array([r[1] / n_src for r in R])
    vol = BIN * BIN * np.array([z1 - z0 for z0, z1 in ZBANDS])
    conv = mass / Q / vol[:, None, None, None]               # count → c/Q [s/m^3]
    snaps = np.stack([R[np.argmin(abs(T - t))][2] * conv for t in SNAP_T])
    am = (T >= AVG_T[0] - 1e-6) & (T <= AVG_T[1] + 1e-6)
    avg = np.mean([R[i][2] for i in np.where(am)[0]], axis=0) * conv
    os.makedirs(os.path.join(ROOT, 'puff', 'data'), exist_ok=True)
    np.savez_compressed(os.path.join(ROOT, 'puff', 'data', f'ref_rot{rot}.npz'), S=S, T=T, dep_frac=dep_frac,
                        rel_frac=rel_frac, snaps=snaps, snap_t=SNAP_T, avg=avg, avg_t=AVG_T,
                        zbands=ZBANDS, bin=BIN, x0=X0, y0=Y0, Q=Q)
    print(f'rot{rot}: {len(steps)} snaps, Q={Q:.1f}, dep frac @t_end mean={dep_frac[-1].mean():.3f} '
          f'[{dep_frac[-1].min():.3f},{dep_frac[-1].max():.3f}], released={rel_frac[-1].mean():.3f}')


if __name__ == '__main__':
    main(int(sys.argv[1]))
