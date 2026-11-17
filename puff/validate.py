#!/usr/bin/env python3
"""퍼프 모델 vs 50소스 평균장 LDM 비교 (각 런 좌표계, 바람 +x).

모델 3개:
  plume   : 단순 가우시안 플룸 (균일풍 +x, Briggs 도시 D 등급, 지면 반사, 전선 x<U(t-t0) 까지)
  upuff   : 건물 없는 균일장 퍼프 (같은 퍼프 코드, 바람·난류는 소스 높이 영역 평균)
  spuff   : SPH 평균장 퍼프

  python3 puff/validate.py 0 [beta] [vd]
"""
import os, sys, time, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from puffmodel import Field, simulate

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, 'puff', 'data')


def metrics(ref, mod, floor):
    """ref, mod: (Ns, nx, ny). 소스별 지표. floor 보다 큰 칸(둘 중 하나라도)만."""
    out = []
    for r, m in zip(ref, mod):
        sel = (r > floor) | (m > floor)
        rr, mm = np.maximum(r[sel], floor), np.maximum(m[sel], floor)
        ratio = mm / rr
        fac2 = np.mean((ratio >= 0.5) & (ratio <= 2))
        lr = np.corrcoef(np.log10(rr), np.log10(mm))[0, 1]
        # 짙은 영역(각자 최대값의 10% 이상)의 겹침
        a, b = r > 0.1 * r.max(), m > 0.1 * m.max()
        fms = (a & b).sum() / max((a | b).sum(), 1)
        out.append((fac2, lr, fms))
    return np.array(out)


def plume_field(src, xc, yc, z, U, times, t0):
    X, Y = np.meshgrid(xc, yc, indexing='ij')
    out = np.zeros((len(src), len(xc), len(yc)))
    for k, (sx, sy, sz) in enumerate(src):
        dxs = X - sx
        acc = 0
        for t in times:
            x = np.where((dxs > 1.0) & (dxs < U * (t - t0)), dxs, np.nan)
            sy_ = 0.16 * x / np.sqrt(1 + 0.0004 * x); sz_ = 0.14 * x / np.sqrt(1 + 0.0003 * x)
            c = 1 / (2 * np.pi * U * sy_ * sz_) * np.exp(-(Y - sy) ** 2 / (2 * sy_ ** 2)) * \
                (np.exp(-(z - sz) ** 2 / (2 * sz_ ** 2)) + np.exp(-(z + sz) ** 2 / (2 * sz_ ** 2)))
            acc = acc + np.nan_to_num(c)
        out[k] = acc / len(times)
    return out


def main(rot, beta=0.5, vd=0.0, models=('plume', 'upuff', 'spuff'), tag=''):
    ref = np.load(os.path.join(D, f'ref_rot{rot}.npz'))
    S = ref['S']; src = S[:, 1:4]; t0, t1 = S[0, 6], S[0, 7] + 3e-3
    b = float(ref['bin'])
    xc = ref['x0'] + b * (np.arange(ref['avg'].shape[2]) + 0.5)
    yc = ref['y0'] + b * (np.arange(ref['avg'].shape[3]) + 0.5)
    zlev = [0.5 * (z0 + z1) for z0, z1 in ref['zbands']]
    a0, a1 = ref['avg_t']
    avg_times = np.arange(a0, a1 + 1e-6, 1.0)
    out_times = sorted(set(avg_times.tolist()) | set(ref['snap_t'].tolist()))

    fld = Field(os.path.join(D, f'field_rot{rot}.npz'))
    # 균일장 퍼프용: 소스 높이 영역(z 2.5~10 m, 건물 밖) 평균
    kz = slice(1, 4)
    m = ~fld.solid[:, :, kz]
    U_ref = float(fld.F[:, :, kz, 0][m].mean()); TLu = float(np.median(fld.F[:, :, kz, 3][m]))
    sgu = float(np.median(fld.F[:, :, kz, 4][m]))
    print(f'rot{rot}: uniform ref U={U_ref:.2f} m/s TL={TLu:.1f} s sig={sgu:.2f} m/s; beta={beta} vd={vd}')

    res = {}
    for name in models:
        tic = time.time()
        if name == 'plume':
            Cavg = np.stack([plume_field(src, xc, yc, z, U_ref, avg_times, t0) for z in zlev])
            dep = None
        else:
            f = fld if name == 'spuff' else Field(uniform=dict(u=U_ref, TL=TLu, sig=sgu))
            C, Dd = simulate(f, src, t0, t1, a1, out_times, (xc, yc), zlev, beta=beta, vd=vd, seed=1)
            sel = [out_times.index(t) for t in avg_times]
            Cavg = C[sel].mean(0)
            dep = Dd[-1]
            np.save(os.path.join(D, f'val_{name}_rot{rot}{tag}.npy'), Cavg)
        res[name] = Cavg
        line = f'  {name:6s} ({time.time() - tic:5.1f}s)'
        if dep is not None:
            line += f' dep mean={dep.mean():.3f} (LDM {ref["dep_frac"][-1].mean():.3f}), corr per-src={np.corrcoef(dep, ref["dep_frac"][-1])[0, 1]:.2f}'
        print(line)
        for iz, z in enumerate(zlev):
            r = ref['avg'][iz]
            floor = 1e-3 * np.median(r.reshape(len(r), -1).max(1))
            mt = metrics(r, Cavg[iz], floor)
            print(f'     z={z:4.0f} m  FAC2={mt[:, 0].mean():.3f}  logR={np.nanmean(mt[:, 1]):.3f}  '
                  f'overlap(>10%max)={mt[:, 2].mean():.3f}   [median over src: FAC2 {np.median(mt[:, 0]):.3f}]')
    return res


if __name__ == '__main__':
    rot = int(sys.argv[1]); beta = float(sys.argv[2]) if len(sys.argv) > 2 else 0.5
    vd = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0
    models = sys.argv[4].split(',') if len(sys.argv) > 4 else ('plume', 'upuff', 'spuff')
    main(rot, beta, vd, models, tag=f'_b{beta}_vd{vd}')
