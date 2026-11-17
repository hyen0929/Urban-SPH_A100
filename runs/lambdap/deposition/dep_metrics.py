#!/usr/bin/env python3
"""소스별 침적 지표 산출 + λp 상관 + Run A(transient) vs Run B(mean-flow) 침적 지도 비교.

입력
  Run A : ../../2026-09-15_leipzig-dt01-sponge-uy/out/LDM_<step>stp.vtk
  Run B : ../../2026-09-15_leipzig-dt01-sponge-uy-meanflow/out/LDM_<step>stp.vtk
  형태  : ../lambda_p_by_R.csv, ../morph_sources.csv, ../sources13.csv

침적 판정: LDM 출력에는 buffer_type 이 없다. KERNEL_dry_deposition(functions_DD.cuh:173-178)이
침적 시 z=L1+1e-4(=1e-4), ux=uy=uz=0 으로 고정하므로 **속도가 정확히 0** 인 입자가 침적 입자다
(활성 입자의 속도 중앙값은 0.4 m/s 이상이라 혼동되지 않는다). 지면은 평지 z=0 이고
지붕·벽면 침적은 모델에 없다.

출력: dep_metrics.csv, dep_lambda_corr.csv, ab_compare.csv, grid_*.npz, 그림 여러 장
"""
import csv, json, os, re, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.normpath(os.path.join(HERE, '..', '..'))
A_DIR = os.path.join(RUNS, '2026-09-15_leipzig-dt01-sponge-uy')
B_DIR = os.path.join(RUNS, '2026-09-15_leipzig-dt01-sponge-uy-meanflow')
MORPH = os.path.normpath(os.path.join(HERE, '..'))
STEP = 30000
CELL = 10.0            # 침적 지도 격자 [m]
COURT = {3, 6, 10}     # 폐쇄 중정 소스 (handoff §5.3)


def read_vtk(path):
    buf = open(path, 'rb').read()
    m = re.search(rb"POINTS\s+(\d+)\s+(\w+)\n", buf); n = int(m.group(1))
    dt = {"float": ">f4", "double": ">f8"}[m.group(2).decode()]
    pts = np.frombuffer(buf, dtype=dt, count=3*n, offset=m.end()).reshape(-1, 3).astype(np.float64)
    f = re.search(rb"FIELD\s+FieldData\s+(\d+)\n", buf); pos = f.end(); F = {}
    for _ in range(int(f.group(1))):
        hm = re.compile(rb"([^\s]+)\s+(\d+)\s+(\d+)\s+(\w+)\n").match(buf, pos)
        if not hm: break
        name, nc, nt, ty = hm.group(1).decode(), int(hm.group(2)), int(hm.group(3)), hm.group(4).decode()
        pos = hm.end(); npd = {"float": ">f4", "double": ">f8", "int": ">i4"}[ty]
        a = np.frombuffer(buf, dtype=npd, count=nc*nt, offset=pos).astype(np.float64)
        pos += nc*nt*np.dtype(npd).itemsize
        if buf[pos:pos+1] == b"\n": pos += 1
        F[name] = a.reshape(nt, nc) if nc > 1 else a
    return pts, F


def load(run_dir, step=STEP):
    p, F = read_vtk(os.path.join(run_dir, 'out', f'LDM_{step}stp.vtk'))
    src = F['p_type'].astype(int) - 100          # 101~113 → 1~13
    dep = np.linalg.norm(F['velocity'], axis=1) == 0.0
    return p, src, dep


def footprint(xy, src_xy, cell=CELL):
    """격자 면적·50/90% 면적·PCA 신장도·질량중심·중심이동."""
    if len(xy) < 10:
        return {}
    ix = np.floor(xy[:, 0]/cell).astype(int); iy = np.floor(xy[:, 1]/cell).astype(int)
    key = ix.astype(np.int64) * 1_000_000 + iy
    uniq, cnt = np.unique(key, return_counts=True)
    order = np.argsort(cnt)[::-1]; c = cnt[order]; cum = np.cumsum(c)/cnt.sum()
    a_cell = cell*cell
    c50 = int(np.searchsorted(cum, 0.50) + 1); c90 = int(np.searchsorted(cum, 0.90) + 1)
    cen = xy.mean(0)
    d = xy - cen
    cov = d.T @ d / len(d)
    ev, evec = np.linalg.eigh(cov)                # ev 오름차순
    major, minor = np.sqrt(ev[1]), np.sqrt(max(ev[0], 1e-12))
    ang = np.degrees(np.arctan2(evec[1, 1], evec[0, 1])) % 180.0
    shift = cen - src_xy
    r = np.linalg.norm(xy - src_xy, axis=1)
    return dict(area_all_m2=len(uniq)*a_cell, area50_m2=c50*a_cell, area90_m2=c90*a_cell,
                elong=major/minor, major_m=major, minor_m=minor, axis_deg=ang,
                cx=cen[0], cy=cen[1], shift_m=float(np.hypot(*shift)),
                shift_deg=float(np.degrees(np.arctan2(shift[1], shift[0])) % 360.0),
                r_med_m=float(np.median(r)), f_r50=float((r <= 50).mean()),
                f_r150=float((r <= 150).mean()))


def spearman(a, b):
    ra = np.argsort(np.argsort(a)).astype(float); rb = np.argsort(np.argsort(b)).astype(float)
    return float(np.corrcoef(ra, rb)[0, 1])


def main():
    S = {int(r['id']): (float(r['x']), float(r['y'])) for r in csv.DictReader(open(os.path.join(MORPH, 'sources13.csv')))}
    lam = {}
    for r in csv.DictReader(open(os.path.join(MORPH, 'lambda_p_by_R.csv'))):
        R = int(float(r['R_m'])); lam[R] = {i: float(r[f'S{i}']) for i in range(1, 14)}
    morph = {}
    for r in csv.DictReader(open(os.path.join(MORPH, 'morph_sources.csv'))):
        morph.setdefault(int(float(r['R'])), {})[int(r['id'])] = r

    rows, grids = [], {}
    for tag, d in (('A', A_DIR), ('B', B_DIR)):
        p, src, dep = load(d)
        for i in sorted(S):
            m = src == i
            md = m & dep
            n_rel, n_dep = int(m.sum()), int(md.sum())
            rec = dict(run=tag, src=i, n_released=n_rel, n_dep=n_dep, f_dep=n_dep/max(n_rel, 1))
            rec.update(footprint(p[md][:, :2], np.array(S[i])))
            # 공중 잔류 중 도메인 밖으로 나간 것 (유출 x>=1307.5)
            air = m & ~dep
            rec['f_air_out'] = float((p[air][:, 0] >= 1307.5).mean()) if air.sum() else 0.0
            rows.append(rec)
        grids[tag] = (p[dep][:, :2], src[dep])
    with open(os.path.join(HERE, 'dep_metrics.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

    # ---- λp 상관 (Run A, Run B, 중정 제외/포함) ----
    byrun = {t: {r['src']: r for r in rows if r['run'] == t} for t in ('A', 'B')}
    metrics = ['f_dep', 'area90_m2', 'elong', 'shift_m', 'f_r150']
    corr_rows = []
    for R in sorted(lam):
        for t in ('A', 'B'):
            for mt in metrics:
                ids = sorted(S)
                x = np.array([lam[R][i] for i in ids]); y = np.array([byrun[t][i][mt] for i in ids])
                keep = np.array([i not in COURT for i in ids])
                corr_rows.append(dict(R=R, run=t, metric=mt,
                                      pearson_all=float(np.corrcoef(x, y)[0, 1]), spearman_all=spearman(x, y),
                                      pearson_nocourt=float(np.corrcoef(x[keep], y[keep])[0, 1]),
                                      spearman_nocourt=spearman(x[keep], y[keep])))
    with open(os.path.join(HERE, 'dep_lambda_corr.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(corr_rows[0].keys())); w.writeheader(); w.writerows(corr_rows)

    # ---- A vs B 침적 지도 비교 ----
    ab = []
    (xyA, sA), (xyB, sB) = grids['A'], grids['B']
    lim = [min(xyA[:, 0].min(), xyB[:, 0].min()), max(xyA[:, 0].max(), xyB[:, 0].max()),
           min(xyA[:, 1].min(), xyB[:, 1].min()), max(xyA[:, 1].max(), xyB[:, 1].max())]
    nx = int((lim[1]-lim[0])/CELL) + 1; ny = int((lim[3]-lim[2])/CELL) + 1

    def grid(xy):
        ix = np.clip(((xy[:, 0]-lim[0])/CELL).astype(int), 0, nx-1)
        iy = np.clip(((xy[:, 1]-lim[2])/CELL).astype(int), 0, ny-1)
        g = np.zeros((nx, ny)); np.add.at(g, (ix, iy), 1.0); return g
    GA, GB = grid(xyA), grid(xyB)
    np.savez_compressed(os.path.join(HERE, 'grid_dep.npz'), GA=GA, GB=GB, lim=np.array(lim), cell=CELL)
    both = (GA > 0) | (GB > 0)
    ab.append(dict(scope='all', cells=int(both.sum()),
                   pearson=float(np.corrcoef(GA[both], GB[both])[0, 1]),
                   spearman=spearman(GA[both], GB[both]),
                   nA=int(GA.sum()), nB=int(GB.sum())))
    for i in sorted(S):
        ga, gb = grid(xyA[sA == i]), grid(xyB[sB == i])
        b = (ga > 0) | (gb > 0)
        ab.append(dict(scope=f'S{i}', cells=int(b.sum()),
                       pearson=float(np.corrcoef(ga[b], gb[b])[0, 1]) if b.sum() > 2 else float('nan'),
                       spearman=spearman(ga[b], gb[b]) if b.sum() > 2 else float('nan'),
                       nA=int(ga.sum()), nB=int(gb.sum())))
    with open(os.path.join(HERE, 'ab_compare.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(ab[0].keys())); w.writeheader(); w.writerows(ab)
    print('dep_metrics.csv / dep_lambda_corr.csv / ab_compare.csv / grid_dep.npz 생성')


if __name__ == '__main__':
    main()
