#!/usr/bin/env python3
"""Leipzig 평균장(mean_field.bin)을 퍼프 모델용 정규 격자로 바꾼다.

mean_field.bin 에는 좌표가 없다. 레코드 i = 같은 런의 전체장 VTK(fluid_<step>stp.vtk)
i 번째 점이다(Eulerian SPH 라 위치 고정, 정렬 순서 동일). 풍향 런마다 IBM 이 회전돼
정렬 순서가 다르므로 반드시 **같은 런의** VTK 를 쓴다.

좌표계는 각 런 그대로(바람 +x, 도시가 회전됨). 실제 풍향으로 되돌리는 건 생성 단계에서 한다.

  python3 puff/prep_fields.py            # 4 방향 전부 → puff/data/field_rot{0,90,180,270}.npz
"""
import os, re, sys, struct
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUNS = {0:   'runs/2026-09-15_leipzig-dt01-sponge-uy',
        90:  'runs/2026-09-23_leipzig-rot90',
        180: 'runs/2026-09-23_leipzig-rot180',
        270: 'runs/2026-09-26_leipzig-rot270'}
VTK_STEP = 30000
DX = 2.5
ZMAX = 200.0          # 이 높이까지만 보관 (건물 최고 115 m)
OUT = os.path.join(ROOT, 'puff', 'data')


def read_vtk_points_ptype(path):
    with open(path, 'rb') as f:
        h = f.read(300)
        m = re.search(rb'POINTS\s+(\d+)\s+float\n', h)
        n = int(m.group(1))
        f.seek(m.end())
        xyz = np.fromfile(f, dtype='>f4', count=3 * n).reshape(n, 3).astype(np.float32)
        t = f.read(200)
        m2 = re.search(rb'p_type\s+1\s+(\d+)\s+int\n', t)
        assert m2 and int(m2.group(1)) == n, 'p_type 필드가 첫 필드가 아님'
        f.seek(f.tell() - len(t) + m2.end())
        pt = np.fromfile(f, dtype='>i4', count=n)
    return xyz, pt


def read_mean(path):
    with open(path, 'rb') as f:
        hdr = f.read(48)
        assert hdr[:14] == b'SOPHIA_MEAN_v1'
        n, cnt = struct.unpack_from('<iI', hdr, 16)
        us, L = struct.unpack_from('<ff', hdr, 24)
        rec = np.fromfile(f, dtype=np.dtype([('ux', '<f4'), ('uy', '<f4'), ('uz', '<f4'),
                                             ('k', '<f4'), ('e', '<f4')]), count=n)
    return rec, dict(n=n, avg_count=cnt, us=us, L=L)


def prep(rot):
    run = os.path.join(ROOT, RUNS[rot])
    xyz, pt = read_vtk_points_ptype(os.path.join(run, 'out', f'fluid_{VTK_STEP}stp.vtk'))
    rec, info = read_mean(os.path.join(run, 'out', 'mean_field.bin'))
    assert len(rec) == len(xyz), (len(rec), len(xyz))
    print(f'rot{rot}: n={len(xyz)} avg_count={info["avg_count"]} us={info["us"]:.3f}')

    fl = (pt == 1) & (rec['ux'] != 0) & (xyz[:, 2] <= ZMAX)
    p = xyz[fl]
    x0, y0, z0 = p[:, 0].min(), p[:, 1].min(), p[:, 2].min()
    ijk = np.rint((p - [x0, y0, z0]) / DX).astype(int)
    err = np.abs(ijk * DX + [x0, y0, z0] - p).max()
    assert err < 0.05, f'격자에서 벗어난 점 {err}'
    nx, ny, nz = ijk.max(0) + 1
    print(f'  grid {nx}x{ny}x{nz}, origin ({x0:.2f},{y0:.2f},{z0:.2f}), fluid pts {fl.sum()}')

    # LDM(KERNEL_time_update_LDM) 과 같은 정의: sigma = sqrt(2k/3), T_L = 2 sigma^2 / (C0 eps), C0 = 1.76
    U = np.zeros((4, nx, ny, nz), np.float32)          # ux, uy, uz, T_L
    sig = np.zeros((nx, ny, nz), np.float32)           # sigma_u
    have = np.zeros((nx, ny, nz), bool)
    r = rec[fl]
    i, j, k = ijk.T
    U[0, i, j, k], U[1, i, j, k], U[2, i, j, k] = r['ux'], r['uy'], r['uz']
    s2 = 2.0 / 3.0 * np.maximum(r['k'], 0)
    U[3, i, j, k] = np.minimum(2.0 * s2 / 1.76 / np.maximum(r['e'], 1e-10), 300.0)
    sig[i, j, k] = np.sqrt(s2)
    have[i, j, k] = True

    # 건물: IBM 마커의 기둥 최대 높이 (내부 유체는 살아 있으므로 높이장으로 막는다)
    b = xyz[pt >= 1000]
    bi = np.floor((b[:, 0] - x0) / DX + 0.5).astype(int); bj = np.floor((b[:, 1] - y0) / DX + 0.5).astype(int)
    ok = (bi >= 0) & (bi < nx) & (bj >= 0) & (bj < ny)
    H = np.zeros((nx, ny), np.float32)
    np.maximum.at(H, (bi[ok], bj[ok]), b[ok, 2])
    zc = z0 + DX * np.arange(nz)
    solid = (zc[None, None, :] <= H[:, :, None]) | ~have
    print(f'  building cols {np.count_nonzero(H > 0)}, max H {H.max():.1f}, solid frac {solid.mean():.3f}')

    os.makedirs(OUT, exist_ok=True)
    np.savez(os.path.join(OUT, f'field_rot{rot}.npz'), U=U, sig=sig, solid=solid, H=H,
             origin=np.array([x0, y0, z0]), dx=DX, us=info['us'], rot=rot)


if __name__ == '__main__':
    for rot in ([int(a) for a in sys.argv[1:]] or [0, 90, 180, 270]):
        prep(rot)
