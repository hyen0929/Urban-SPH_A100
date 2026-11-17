import sys, numpy as np
sys.path.insert(0,'/tmp/claude-1014/-home-kdh-Urban-SPH/90fe7215-f210-4eca-bb99-e4c2226c7eef/scratchpad')
from vtkread import read
STEPS=list(range(14000,15001,100))
R={'neutral':'runs/2026-09-15_leipzig-dt01-sponge-uy/out/fluid_condition_%dstp.vtk',
   'stable' :'runs/2026-09-21_leipzig-most-stable/out/fluid_condition_%dstp.vtk'}

# 건물 위치 (y=0 슬랩의 IBM 마커 x 범위)
n,pts,fl=read('runs/2026-09-21_leipzig-most-stable/out/fluid_condition_15000stp.vtk')
xb=np.array(pts[0::3]); yb=np.array(pts[1::3]); zb=np.array(pts[2::3]); ptb=np.array(fl['p_type'])
ib=(ptb>=1000)&(np.abs(yb)<1.25)
bx=np.sort(np.unique(np.round(xb[ib]/2.5)*2.5))
print(f"y~0 근방 IBM 마커 {ib.sum()}개, x {bx.min():.0f}~{bx.max():.0f} m")

out={}
for lab,tpl in R.items():
    U=[]
    for s in STEPS:
        n,pts,fl=read(tpl%s)
        y=np.array(pts[1::3]); z=np.array(pts[2::3]); pt=np.array(fl['p_type'])
        m=(pt==1)&(np.abs(y)<1e-3)&(np.abs(z-1.25)<1e-3)
        v=np.array(fl['velocity'])
        if lab not in out: X=np.array(pts[0::3])[m]
        U.append(v[2::3][m])
    out[lab]=(X,np.array(U))
    print(f"{lab}: z=1.25 평면 {m.sum():,}개")

X=out['neutral'][0]
un=out['neutral'][1]; us=out['stable'][1]

print("\n=== z=1.25 평면, 스냅샷별 평균 uz (부호 일관성) ===")
for lab,U in (('neutral',un),('stable',us)):
    per=U.mean(axis=1)
    print(f"  {lab:8s}: "+" ".join(f"{p:+.3f}" for p in per))
    print(f"            음 {sum(p<0 for p in per)}/{len(per)},  평균 {per.mean():+.4f} m/s")

print("\n=== 공간 구조: 건물 유무별 (x 구간) ===")
in_city=(X>=bx.min())&(X<=bx.max())
for nm,sel in (('건물역 (x %.0f~%.0f)'%(bx.min(),bx.max()),in_city),
               ('개방역 (상·하류)',~in_city)):
    a=un[:,sel].mean(); b=us[:,sel].mean()
    fa=(un[:,sel]<0).mean()*100; fb=(us[:,sel]<0).mean()*100
    print(f"  {nm:26s} n={sel.sum():5d} | 중립 {a:+.4f} ({fa:.1f}%)  안정 {b:+.4f} ({fb:.1f}%)")

print("\n=== 균일한가, 조직됐는가 ===")
ms=us.mean(axis=0); mn=un.mean(axis=0)
print(f"  안정 z=1.25 열별 평균 uz: 중앙 {np.median(ms):+.4f}  p10 {np.percentile(ms,10):+.4f}  p90 {np.percentile(ms,90):+.4f}  표준편차 {ms.std():.4f}")
print(f"  중립                   : 중앙 {np.median(mn):+.4f}  p10 {np.percentile(mn,10):+.4f}  p90 {np.percentile(mn,90):+.4f}  표준편차 {mn.std():.4f}")
print(f"  안정에서 x열의 {100*(ms<0).mean():.1f}% 가 하강 (중립 {100*(mn<0).mean():.1f}%)")
np.savez('runs/2026-09-21_leipzig-most-stable/analysis/uz_z125.npz', x=X, uz_n=mn, uz_s=ms, bx_min=bx.min(), bx_max=bx.max())
