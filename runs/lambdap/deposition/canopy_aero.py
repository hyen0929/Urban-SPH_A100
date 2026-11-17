#!/usr/bin/env python3
"""도시 캐노피 공기역학 — 150 s 평균장에서 z0(거칠기 길이)·d(변위 높이)를 구하고
형태학적 추정식(Macdonald et al. 1998)과 비교한다.

입력
  mean_field.bin  : Run A 의 t=150~299 s 평균 (ux,uy,uz,k,e), 좌표 없음.
                    레코드 i = 덤프 시점 dev_P1 의 i 번째 입자 = 같은 런 전체장 VTK 의 i 번째 점.
                    (Eulerian SPH 라 위치 고정, NNPS 정렬은 count==0 한 번뿐)
  fluid_30000stp.vtk : 같은 런의 좌표·p_type (정렬 순서 동일)
  ../heightfield.npz : 2.5 m 격자 건물 높이장 H(x,y) — 건물 내부 입자 제외에 사용

출력: canopy_profile.csv, canopy_aero.json, 화면 요약
"""
import json, os, struct
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.normpath(os.path.join(HERE, '..', '..'))
A = os.path.join(RUNS, '2026-09-15_leipzig-dt01-sponge-uy')
MORPH = os.path.normpath(os.path.join(HERE, '..'))
KAPPA = 0.4

ns = {'__file__': os.path.join(HERE, 'dep_metrics.py')}
exec(open(os.path.join(HERE, 'dep_metrics.py')).read().split('def load(')[0], ns)
read_vtk = ns['read_vtk']


def read_mean(path):
    with open(path, 'rb') as f:
        hdr = f.read(48)
        n, ac = struct.unpack_from('<Ii', hdr, 16)
        us, ob = struct.unpack_from('<ff', hdr, 24)
        rec = np.fromfile(f, dtype=np.dtype([('ux', '<f4'), ('uy', '<f4'), ('uz', '<f4'),
                                             ('k', '<f4'), ('e', '<f4')]), count=n)
    return rec, n, ac, us


def main():
    rec, n, ac, us = read_mean(os.path.join(A, 'out', 'mean_field.bin'))
    p, F = read_vtk(os.path.join(A, 'out', 'fluid_30000stp.vtk'))
    assert len(p) == n, f'입자 수 불일치 {len(p)} vs {n}'
    fluid = F['p_type'] == 1
    nz = rec['ux'] != 0
    inst = F['velocity'][:, 0]
    print(f"[정렬 확인] 평균장 비-0 {nz.sum():,} / VTK 유체 {int(fluid.sum()):,} "
          f"| 교집합 {int((nz & fluid).sum()):,} ({100*(nz & fluid).sum()/max(nz.sum(),1):.2f}% of 비-0)")
    r = np.corrcoef(rec['ux'][nz & fluid], inst[nz & fluid])[0, 1]
    print(f"[정렬 확인] 평균 ux 대 순간 ux 상관 {r:.3f} (순서가 어긋나면 ~0 이 나온다)")

    # --- 건물 내부 입자 제외 (높이장) ---
    hf = np.load(os.path.join(MORPH, 'heightfield.npz'))
    H, x0, y0, cell = hf['H'], float(hf['x0']), float(hf['y0']), float(hf['cell'])
    box = hf['builtup_box']
    # H 는 [y, x] 순서다 (shape 363x368, y 363행). morph_metrics.py 의 frontal/wall 도 같은 축.
    ix = np.floor((p[:, 0]-x0)/cell).astype(int); iy = np.floor((p[:, 1]-y0)/cell).astype(int)
    inside_grid = (iy >= 0) & (iy < H.shape[0]) & (ix >= 0) & (ix < H.shape[1])
    hz = np.zeros(len(p)); hz[inside_grid] = H[iy[inside_grid], ix[inside_grid]]
    in_building = p[:, 2] <= hz + 0.1
    inbox = (p[:, 0] >= box[0]) & (p[:, 0] <= box[1]) & (p[:, 1] >= box[2]) & (p[:, 1] <= box[3])
    use = fluid & nz & inbox & ~in_building
    print(f"[표본] 시가지 박스 내 유체 {int((fluid&inbox).sum()):,} → 건물 내부 제외 후 {int(use.sum()):,}")

    # --- 수평 평균 프로파일 ---
    zb = np.arange(0, 200.1, 2.5)
    zc, U, V, W, TKE, N = [], [], [], [], [], []
    for a, b in zip(zb[:-1], zb[1:]):
        m = use & (p[:, 2] >= a) & (p[:, 2] < b)
        if m.sum() < 500: continue
        zc.append((a+b)/2); N.append(int(m.sum()))
        U.append(float(rec['ux'][m].mean())); V.append(float(rec['uy'][m].mean()))
        W.append(float(rec['uz'][m].mean())); TKE.append(float(rec['k'][m].mean()))
    zc, U, V, W, TKE, N = map(np.array, (zc, U, V, W, TKE, N))
    with open(os.path.join(HERE, 'canopy_profile.csv'), 'w') as f:
        f.write('z_m,U_mean,V_mean,W_mean,k_turb,n_particles\n')
        for row in zip(zc, U, V, W, TKE, N):
            f.write('%.2f,%.5f,%.5f,%.5f,%.5f,%d\n' % row)

    # --- 대수법칙 적합: U(z) = (u*/κ) ln((z-d)/z0), 캐노피 위 구간 ---
    Hmean = float(H[H > 1].mean()); Hmax = float(H.max())
    lo, hi = 2.0*Hmean, 8.0*Hmean
    best = None
    for d in np.arange(0.0, 0.95*Hmean, 0.1):
        m = (zc > max(lo, d + 2.0)) & (zc < hi)
        if m.sum() < 5: continue
        X = np.log(zc[m] - d)
        b, a = np.polyfit(X, U[m], 1)            # U = b*ln(z-d) + a
        res = U[m] - (b*X + a)
        rmse = float(np.sqrt((res**2).mean()))
        if b <= 0: continue
        ustar = b*KAPPA; z0 = float(np.exp(-a/b))
        if best is None or rmse < best['rmse']:
            best = dict(d=float(d), z0=z0, ustar=float(ustar), rmse=rmse, n=int(m.sum()),
                        fit_lo=float(max(lo, d+2.0)), fit_hi=float(hi))

    # --- 형태학적 추정 (Macdonald 1998), 시가지 박스 전체 ---
    inb = np.zeros_like(H, dtype=bool)
    gy = y0 + (np.arange(H.shape[0])+.5)*cell; gx = x0 + (np.arange(H.shape[1])+.5)*cell
    GY, GX = np.meshgrid(gy, gx, indexing='ij')
    inb = (GX >= box[0]) & (GX <= box[1]) & (GY >= box[2]) & (GY <= box[3])
    built = inb & (H > 1)
    lam_p = float(built.sum()/inb.sum())
    # frontal_step_map 은 셀당 면적 [m²] 를 돌려준다(계단높이 x cell) — 추가로 cell 을 곱하지 말 것
    lam_f = float(hf['frontal'][inb].sum()/(inb.sum()*cell*cell))
    Hb = float(H[built].mean())
    Amac, beta, Cd = 4.43, 1.0, 1.2
    d_mac = Hb*(1 + Amac**(-lam_p)*(lam_p - 1))
    z0_mac = Hb*(1 - d_mac/Hb)*np.exp(-(0.5*beta*Cd/KAPPA**2*(1 - d_mac/Hb)*lam_f)**-0.5)

    out = dict(mean_field=dict(records=int(n), avg_count=int(ac), us_ini=float(us)),
               sample=dict(particles=int(use.sum()), box=[float(v) for v in box]),
               morphology=dict(lambda_p=lam_p, lambda_f=lam_f, H_built=Hb, H_max=Hmax),
               fit=best, macdonald=dict(d=float(d_mac), z0=float(z0_mac)))
    json.dump(out, open(os.path.join(HERE, 'canopy_aero.json'), 'w'), indent=2, ensure_ascii=False)

    print(f"\n[형태] 시가지 박스 λp={lam_p:.3f}  λf={lam_f:.3f}  H_built={Hb:.1f} m (최대 {Hmax:.0f} m)")
    print(f"[프로파일] z={zc[0]:.1f}~{zc[-1]:.1f} m, {len(zc)}구간")
    print(f"[대수법칙 적합] 구간 {best['fit_lo']:.0f}~{best['fit_hi']:.0f} m ({best['n']}점): "
          f"d={best['d']:.1f} m, z0={best['z0']:.3f} m, u*={best['ustar']:.3f} m/s, RMSE={best['rmse']:.4f} m/s")
    print(f"[Macdonald 1998] d={d_mac:.1f} m, z0={z0_mac:.3f} m  "
          f"→ 비 d {best['d']/d_mac:.2f}, z0 {best['z0']/z0_mac:.2f}")
    print(f"[참고] 모델 입력 us_ini={us:.2f} m/s, MOST z0 설정값과 비교할 것")


if __name__ == '__main__':
    main()
