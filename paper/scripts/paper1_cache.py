#!/usr/bin/env python3
"""1편 4장 그림용 중간 자료 캐시.

무거운 VTK(전체장 2.1 GB)를 한 번만 읽어 그림에 필요한 것만 npz 로 남긴다.

  Run A (전이)  : runs/2026-10-20_leipzig-s13-most-wall
  Run B (평균장): runs/2026-10-23_leipzig-s13-most-meanflow   (fluid_0stp.vtk = 적용된 시간평균장)

출력: paper/data/10_paper1/cache_*.npz
"""
import os, sys
import numpy as np

ROOT = '/home/kdh/Urban-SPH'
A = f'{ROOT}/runs/2026-10-20_leipzig-s13-most-wall/out'
B = f'{ROOT}/runs/2026-10-23_leipzig-s13-most-meanflow/out'
OUT = f'{ROOT}/paper/data/10_paper1'
os.makedirs(OUT, exist_ok=True)

ns = {'__file__': f'{ROOT}/runs/lambdap/deposition/dep_metrics.py'}
exec(open(ns['__file__']).read().split('def load(')[0], ns)
read_vtk = ns['read_vtk']

hf = np.load(f'{ROOT}/paper/data/00_setup/heightfield.npz')
H, hx0, hy0, hc = hf['H'], float(hf['x0']), float(hf['y0']), float(hf['cell'])


def interior(x, y, z):
    """건물 내부(z <= 그 기둥의 건물고) 여부. 높이장 범위 밖은 바깥."""
    i = np.floor((x - hx0) / hc).astype(int)
    j = np.floor((y - hy0) / hc).astype(int)
    ok = (i >= 0) & (i < H.shape[1]) & (j >= 0) & (j < H.shape[0])
    h = np.zeros_like(z)
    h[ok] = H[j[ok], i[ok]]
    return z <= h


def plane(p, F, zz, field):
    m = np.abs(p[:, 2] - zz) < 1e-3
    return p[m, 0], p[m, 1], field[m]


def what(step):
    which = sys.argv[1:] or ['B', 'A', 'ldm']
    return step in which


# ── 평균장: 프로파일 + 보행자층·캐노피 중간 평면 ────────────────────────
if what('B'):
    p, F = read_vtk(f'{B}/fluid_0stp.vtk')
    fl = F['p_type'].astype(int) == 1
    p, v, k = p[fl], F['velocity'][fl], F['k_turb'][fl]
    uh = np.hypot(v[:, 0], v[:, 1])
    inn = interior(p[:, 0], p[:, 1], p[:, 2])
    zs = np.unique(np.round(p[:, 2], 3)); zs = zs[(zs > 0) & (zs < 372)]
    prof = []
    for zz in zs:
        m = (np.abs(p[:, 2] - zz) < 1e-3) & ~inn
        if m.sum() < 100: continue
        prof.append((zz, uh[m].mean(), v[m, 0].mean(), k[m].mean(), m.sum()))
    prof = np.array(prof)
    out = {'prof': prof}
    for zz in (1.25, 11.25):
        x, y, f = plane(p, F, zz, np.c_[v, k])
        out[f'z{zz:g}_x'], out[f'z{zz:g}_y'], out[f'z{zz:g}_f'] = x, y, f
    np.savez_compressed(f'{OUT}/cache_meanfield.npz', **out)
    print('평균장 캐시:', len(prof), '면')

# ── 전이 순시장: 캐노피 중간 와도 ─────────────────────────────────────
if what('A'):
    p, F = read_vtk(f'{A}/fluid_30000stp.vtk')
    fl = F['p_type'].astype(int) == 1
    p, w, v = p[fl], F['vorticity'][fl], F['velocity'][fl]
    out = {}
    for zz in (1.25, 11.25):
        m = np.abs(p[:, 2] - zz) < 1e-3
        out[f'z{zz:g}_x'], out[f'z{zz:g}_y'] = p[m, 0], p[m, 1]
        out[f'z{zz:g}_w'] = np.linalg.norm(w[m], axis=1)
        out[f'z{zz:g}_u'] = np.hypot(v[m, 0], v[m, 1])
    np.savez_compressed(f'{OUT}/cache_inst.npz', **out)
    print('순시장 캐시 저장')

# ── LDM 최종 상태 (A, B) ─────────────────────────────────────────────
if what('ldm'):
    out = {}
    for tag, d in (('A', A), ('B', B)):
        p, F = read_vtk(f'{d}/LDM_30000stp.vtk')
        out[f'{tag}_xyz'] = p.astype(np.float32)
        out[f'{tag}_src'] = (F['p_type'].astype(int) - 100).astype(np.int16)
        out[f'{tag}_bt'] = F['buffer_type'].astype(np.int8)
    np.savez_compressed(f'{OUT}/cache_ldm.npz', **out)
    print('LDM 캐시 저장')
