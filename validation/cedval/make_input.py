#!/usr/bin/env python3
"""CEDVAL Complexity 2 (BL3-2-1, 엇갈림 큐브 배열) 입력 생성.

Leipzig_surface 와 같은 규칙을 따른다.
  격자 2.5 m, 지면 p_type 0 네 층(z -8.75 … -1.25), 유체 p_type 1, 천장 p_type -1 세 층,
  y 주기(버퍼 4 면), 유입·유출 버퍼 4 면. 건물은 표면 마커(p_type 1000, 간격 ~1 m,
  질량은 격자 입자와 같게 — Leipzig 와 동일). 건물 내부 격자점은 유체로 둔다(Leipzig 와 동일).

좌표: SPH x = CEDVAL x + X_SHIFT. 바람은 +x. CEDVAL 원점은 배열 중심.
방출점은 자료 머리말에 x=+278.55 로 적혀 있으나 도착 시간 순서로 보아 x=-278.55 (배열 바로 앞).

출력: input/<CASE>/<CASE>.txt, <CASE>_LDM.txt, <CASE>_solv.txt, 그리고 receptors.csv
사용: python3 validation/cedval/make_input.py [CASE]
  CEDVAL_C2       : dx 2.5, 배열 전체 폭 (CASE_ID 1020)
  CEDVAL_C2_dx15  : dx 1.5, y 주기 432 m(건물 12 열 = 지붕 패턴 2 주기)로 배열을 무한 반복,
                    상·하류와 높이를 줄인 도메인 (CASE_ID 1021). 벽·지붕이 격자점 사이 가운데 오게 맞춤
"""
import os, sys
import numpy as np

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
CASE = sys.argv[1] if len(sys.argv) > 1 else 'CEDVAL_C2'
OUT = os.path.join(ROOT, 'input', CASE)
HERE = os.path.dirname(os.path.abspath(__file__))

PRESET = {
    # DX, X_SHIFT(SPH x = CEDVAL x + X_SHIFT), NX, y 격자, 주기 반폭(None=주기 반복 안 함), Z_TOP, 마커 간격, 퍼프
    'CEDVAL_C2': dict(DX=2.5, X_SHIFT=400.0, NX=327, Y=np.arange(-367.5, 367.5 + 1e-9, 2.5), YP=None,
                      Z_TOP=216.25, DS=1.0, T_FIRST=200.0, N_PUFF=60, T_END=1000.0),
    'CEDVAL_C2_dx15': dict(DX=1.5, X_SHIFT=341.25, NX=447, Y=0.75 + 1.5 * np.arange(-148, 148), YP=216.0,
                           Z_TOP=195.75, DS=0.6, T_FIRST=150.0, N_PUFF=60, T_END=1000.0),
}[CASE]
DX = PRESET['DX']
X_SHIFT = PRESET['X_SHIFT']
NX = PRESET['NX']
Z_GROUND = [-DX * (k + 0.5) for k in range(4)][::-1]
Z_TOP = PRESET['Z_TOP']                # 천장 최상면 (= L_ceiling)
NCEIL = 3
RHO, T0, H = 1.225, 300.0, 1.5 * DX
M = RHO * DX**3

# 건물
B, HB, HR = 18.0, 18.0, 27.0         # 한 변, 벽 높이, 27 m 지붕/용마루
PITCH = 36.0
DS = PRESET['DS']                    # 표면 마커 간격 [m] (격자 간격의 ~0.4)

# LDM 퍼프
SRC = (-278.55, 0.0, 1.0)            # CEDVAL 좌표, z 는 첫 유체층 아래
N_PUFF, PUFF_EVERY, T_FIRST = PRESET['N_PUFF'], 5.0, PRESET['T_FIRST']
NP_PUFF = 20000
PUFF_DUR = 0.5                       # 방출 지속 [s] (LDM 갱신 10 step = 0.1 s)
T_END = PRESET['T_END']


def roof_pattern():
    rows = [l.strip() for l in open(os.path.join(HERE, 'setup', 'roof_pattern.txt'))
            if l.strip() and not l.startswith('#')]
    assert len(rows) == 15
    bl = []
    yp = PRESET['YP']
    for r, row in enumerate(rows):
        xc = 252.0 - PITCH * r                       # 위쪽 행 = 하류(+x)
        y0 = 252.0 if len(row) == 15 else 234.0
        if yp is None:
            for k, t in enumerate(row):
                bl.append((xc, y0 - PITCH * k, t))
            continue
        # 주기 반복: 지붕 패턴은 열 6 개 주기 → 2·yp = 12 열이면 이음매 없음
        assert (2 * yp) % (6 * PITCH) == 0
        k0 = 0                                       # 이 행의 가운데 근처 열을 기준으로 패턴을 이어 붙인다
        for kk in range(-20, 21):
            yc = y0 - PITCH * kk
            if abs(yc) > yp + B:                     # 주기 경계를 걸친 건물(|yc|=yp)은 양쪽 모두 둔다
                continue
            bl.append((xc, yc, row[(kk - k0) % 6]))
    return bl


def grid(a, b, n):
    """구간 [a,b] 를 간격 ~DS 로 나눈 셀 중심 n 개."""
    return a + (np.arange(n) + 0.5) * (b - a) / n


def rect(p0, u, v, lu, lv):
    """p0 에서 u,v 방향으로 lu×lv 사각면 위 마커."""
    nu, nv = max(1, round(lu / DS)), max(1, round(lv / DS))
    su, sv = np.meshgrid(grid(0, lu, nu), grid(0, lv, nv), indexing='ij')
    return p0 + su.reshape(-1, 1) * u + sv.reshape(-1, 1) * v


def building(xc, yc, t):
    ex, ey, ez = np.eye(3)
    x0, y0 = xc - B / 2, yc - B / 2
    hw = HR if t == 'X' else HB
    P = [rect(np.array([x0, y0, 0]), ey, ez, B, hw),           # 상류면 (x0)
         rect(np.array([x0 + B, y0, 0]), ey, ez, B, hw),       # 하류면
         rect(np.array([x0, y0, 0]), ex, ez, B, hw),           # y0 면
         rect(np.array([x0, y0 + B, 0]), ex, ez, B, hw)]       # y0+B 면
    if t in 'FX':
        P.append(rect(np.array([x0, y0, hw]), ex, ey, B, B))
        return np.vstack(P)
    # 경사지붕: 'H' 용마루가 y 방향(도면 가로선), 'V' 용마루가 x 방향
    along, across = (ey, ex) if t == 'H' else (ex, ey)
    rise = HR - HB
    slope = np.hypot(B / 2, rise)
    for sgn in (+1, -1):
        # 처마(z=HB)에서 용마루(z=HR)로 올라가는 면
        eave = np.array([xc, yc, HB]) - sgn * across * (B / 2) - along * (B / 2)
        dirup = sgn * across * (B / 2) / slope + ez * rise / slope
        P.append(rect(eave, along, dirup, B, slope))
    # 박공(삼각형) 두 개: 용마루 양 끝 면
    for s in (0.0, B):
        base = np.array([xc, yc, HB]) - along * (B / 2) + along * s - across * (B / 2)
        tri = rect(base, across, ez, B, rise)
        a = tri - (np.array([xc, yc, HB]) - along * (B / 2) + along * s)
        dc = np.abs(a @ across)                      # 용마루 축에서의 거리
        keep = (tri[:, 2] - HB) <= rise * (1 - dc / (B / 2)) + 1e-9
        P.append(tri[keep])
    return np.vstack(P)


def main():
    os.makedirs(OUT, exist_ok=True)
    bl = roof_pattern()
    ibm = np.vstack([building(x, y, t) for x, y, t in bl])
    ibm[:, 0] += X_SHIFT

    if PRESET['YP'] is not None:                     # 버퍼 4 면 바깥 마커는 버린다
        ibm = ibm[np.abs(ibm[:, 1]) < PRESET['YP'] + 4 * DX]
    xs = np.arange(NX) * DX
    ys = PRESET['Y']
    zf = np.arange(DX / 2, Z_TOP - NCEIL * DX + 1e-9, DX)
    zc = Z_TOP - DX * np.arange(NCEIL)[::-1]
    layers = [(z, 0) for z in Z_GROUND] + [(z, 1) for z in zf] + [(z, -1) for z in zc]
    X, Y = np.meshgrid(xs, ys, indexing='ij'); X, Y = X.ravel(), Y.ravel()

    hdr = '1\t2\t3\t4\t7\t8\t9\t10\t12\t26\t27\n'
    fmt = '%.6e\t%.6e\t%.6e\t0.000000e+00\t%.6e\t%d\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n'
    n = {0: 0, 1: 0, -1: 0, 1000: len(ibm)}
    with open(os.path.join(OUT, f'{CASE}.txt'), 'w') as f:
        f.write(hdr)
        for z, pt in layers:
            for x, y in zip(X, Y):
                f.write(fmt % (x, y, z, M, pt, H, T0, RHO, H, M))
            n[pt] += len(X)
        for x, y, z in ibm:
            f.write(fmt % (x, y, z, M, 1000, H, T0, RHO, H, M))
    print('입자:', n, '합', sum(n.values()))

    # LDM: 퍼프마다 p_type 101, 102, …
    mL = DX**3
    with open(os.path.join(OUT, f'{CASE}_LDM.txt'), 'w') as f:
        f.write('1\t2\t3\t4\t7\t8\t9\t12\t26\t27\t28\n')
        sx, sy, sz = SRC[0] + X_SHIFT, SRC[1], SRC[2]
        for k in range(N_PUFF):
            t0 = T_FIRST + k * PUFF_EVERY
            for j in range(NP_PUFF):
                f.write('%.6e\t%.6e\t%.6e\t0.000000e+00\t%.6e\t%d\t%.6e\t1.000000e+00\t%.6e\t%.6e\t%.6e\n'
                        % (sx, sy, sz, mL, 101 + k, H, H, mL, t0 + PUFF_DUR * j / NP_PUFF))
    print('LDM:', N_PUFF * NP_PUFF, f'(퍼프 {N_PUFF} × {NP_PUFF})')

    # solv: Leipzig 4장 런과 같고 끝 시간만 다름
    src = open(os.path.join(ROOT, 'input', 'Leipzig_s10_wall', 'Leipzig_s10_wall_solv.txt')).read()
    src = src.replace('end-time (sec):\t\t\t\t\t\t\t\t\t\t  300.0', f'end-time (sec):\t\t\t\t\t\t\t\t\t\t  {T_END:.1f}')
    assert f'{T_END:.1f}' in src
    src = src.replace('7   concentration p_type velocity e_turb k_turb concn buffer_type', '1   p_type')
    assert '1   p_type' in src
    open(os.path.join(OUT, f'{CASE}_solv.txt'), 'w').write(src)

    # 건물 표 (그림·후처리용)
    with open(os.path.join(OUT, 'buildings.csv'), 'w') as f:
        f.write('x_cedval,y,type\n')
        for x, y, t in bl:
            f.write(f'{x},{y},{t}\n')
    print('범위 x', xs[0], xs[-1], 'y', ys[0], ys[-1], 'z', layers[0][0], layers[-1][0], 'h', H,
          '| 건물', len(bl), '| 마커', len(ibm))


if __name__ == '__main__':
    main()
