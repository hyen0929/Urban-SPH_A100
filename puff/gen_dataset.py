#!/usr/bin/env python3
"""드론 소스추적 학습용 데이터 생성 (SPH 평균장 퍼프, Leipzig, 중립, 풍향 4개).

케이스 = (풍향, 소스 위치). 모든 출력은 **실제 좌표계**(rot0 런 좌표, 건물 고정)로 되돌려 저장한다.
  런 좌표 → 실제 좌표: p_true = C + R(-θ)(p_run - C),  C = (657.5, 0)
  격자 배열은 np.rot90(k=-θ/90, axes=(x,y)) 로 정확히 사상된다(5 m 칸 중심 131 = C).

출력 (out_dir/):
  static.npz            건물 높이장 H(실제 좌표), 각 풍향의 평균 바람·난류 단면 (z 레벨별)
  rot{θ}_b{n}.npz       케이스 묶음: src(실제 좌표), wind_from(°), c[case, t, z, x, y] (float16, c/Q [s/m^3])

  python3 puff/gen_dataset.py --n 50 --rot 0 90 180 270 --out puff/data/dataset_v0
"""
import os, sys, argparse, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from puffmodel import Field, simulate

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, 'puff', 'data')
C = np.array([657.5, 0.0])
BIN, NB = 5.0, 263
XC = 2.5 + BIN * np.arange(NB)
YC = -655.0 + BIN * np.arange(NB)
BETA, VD = 0.8, 0.16              # validate.py 로 50소스 LDM 에 맞춘 값 (2026-11-17)


def to_true(p, th):
    a = np.radians(-th); c, s = np.cos(a), np.sin(a)
    q = p[..., :2] - C
    out = p.copy()
    out[..., 0] = C[0] + c * q[..., 0] - s * q[..., 1]
    out[..., 1] = C[1] + s * q[..., 0] + c * q[..., 1]
    return out


def to_run(p, th):
    return to_true(p, -th)


def grid_to_true(a, th):
    """(..., x, y) 배열을 실제 좌표로."""
    return np.rot90(a, k=-(th // 90), axes=(-2, -1))


def wind_from_deg(th):
    """런 좌표의 +x 바람이 실제 좌표에서 어느 방위에서 불어오는가 (기상 관례: 0=북, 90=동)."""
    a = np.radians(-th)
    u, v = np.cos(a), np.sin(a)                 # 실제 좌표 바람 벡터 (x=동, y=북 가정)
    return float(np.degrees(np.arctan2(-u, -v)) % 360)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=50, help='풍향당 케이스 수')
    ap.add_argument('--batch', type=int, default=50)
    ap.add_argument('--rot', type=int, nargs='+', default=[0, 90, 180, 270])
    ap.add_argument('--z', type=float, nargs='+', default=[5.0, 20.0], help='저장 높이 [m]')
    ap.add_argument('--src_z', type=float, nargs=2, default=[5.0, 5.0])
    ap.add_argument('--t_end', type=float, default=300.0, help='방출 시작 후 기록 끝 [s]')
    ap.add_argument('--every', type=float, default=5.0, help='저장 간격 [s]')
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--out', default=os.path.join(D, 'dataset_v0'))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    static = {}
    times = np.arange(a.every, a.t_end + 1e-6, a.every)
    for th in a.rot:
        fld = Field(os.path.join(D, f'field_rot{th}.npz'))
        # 평균 바람·난류 단면 (실제 좌표, 5 m 칸)
        ii = np.rint((XC - fld.o[0]) / fld.dx).astype(int); jj = np.rint((YC - fld.o[1]) / fld.dx).astype(int)
        for z in a.z:
            k = int(np.rint((z - fld.o[2]) / fld.dx))
            F = fld.F[np.ix_(ii, jj, [k])][:, :, 0]                  # x, y, 5
            ang = np.radians(-th); c, s = np.cos(ang), np.sin(ang)
            u = c * F[..., 0] - s * F[..., 1]; v = s * F[..., 0] + c * F[..., 1]
            static[f'rot{th}_z{z:g}'] = grid_to_true(np.stack([u, v, F[..., 2], F[..., 4]]), th).astype(np.float32)
        if th == 0:
            static['H'] = fld.H
        static[f'rot{th}_wind_from'] = wind_from_deg(th)
        rng = np.random.default_rng(a.seed * 1000 + th)
        src_true = []
        while len(src_true) < a.n:     # 실제 좌표에서 균일하게 뽑고, 런 좌표에서 주변 10 m 가 비었는지 검사
            r = 410.0 * np.sqrt(rng.random()); ph = 2 * np.pi * rng.random()
            p = np.array([C[0] + r * np.cos(ph), C[1] + r * np.sin(ph), rng.uniform(*a.src_z)])
            ring = to_run(p + np.array([[dx, dy, 0] for dx in (-10, 0, 10) for dy in (-10, 0, 10)]), th)
            if not fld.is_solid(ring).any():
                src_true.append(p)
        src_true = np.array(src_true)
        for bi in range(0, a.n, a.batch):
            tic = time.time()
            st = src_true[bi:bi + a.batch]
            Cc, Dd = simulate(fld, to_run(st, th), 0.0, a.t_end, a.t_end, times, (XC, YC), a.z,
                              beta=BETA, vd=VD, seed=a.seed * 7919 + th * 13 + bi)
            # (t, z, src, x, y) → (src, t, z, x, y), 실제 좌표로
            cc = grid_to_true(np.transpose(Cc, (2, 0, 1, 3, 4)), th).astype(np.float16)
            fn = os.path.join(a.out, f'rot{th}_b{bi // a.batch:03d}.npz')
            np.savez_compressed(fn, src=st, wind_from=wind_from_deg(th), rot=th, t=times, z=np.array(a.z),
                                xc=XC, yc=YC, c=cc, dep=Dd[-1])
            print(f'rot{th} batch {bi // a.batch}: {len(st)} cases {time.time() - tic:.0f}s → {fn} '
                  f'({os.path.getsize(fn) / 1e6:.1f} MB)', flush=True)
    np.savez_compressed(os.path.join(a.out, 'static.npz'), xc=XC, yc=YC, z=np.array(a.z), **static)


if __name__ == '__main__':
    main()
