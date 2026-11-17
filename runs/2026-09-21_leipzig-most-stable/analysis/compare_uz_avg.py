import sys, numpy as np
sys.path.insert(0,'/tmp/claude-1014/-home-kdh-Urban-SPH/90fe7215-f210-4eca-bb99-e4c2226c7eef/scratchpad')
from vtkread import read

STEPS = list(range(14000, 15001, 100))      # t=140~150 s, 11 스냅샷
RUNS = {
 'neutral': 'runs/2026-09-15_leipzig-dt01-sponge-uy/out/fluid_condition_%dstp.vtk',
 'stable' : 'runs/2026-09-21_leipzig-most-stable/out/fluid_condition_%dstp.vtk',
}
BANDS=[(0,2),(2,5),(5,10),(10,20),(20,30),(30,50),(50,100)]

acc={}
for lab,tpl in RUNS.items():
    UZ=[]; UH=[]; Z=None; X=None
    for s in STEPS:
        n,pts,fl=read(tpl%s)
        y=np.array(pts[1::3]); pt=np.array(fl['p_type'])
        m=(pt==1)&(np.abs(y)<1e-3)
        v=np.array(fl['velocity'])
        if Z is None:
            Z=np.array(pts[2::3])[m]; X=np.array(pts[0::3])[m]
        UZ.append(v[2::3][m]); UH.append(np.hypot(v[0::3][m],v[1::3][m]))
    acc[lab]=dict(z=Z,x=X,uz=np.array(UZ),uh=np.array(UH))
    print(f"{lab}: {len(STEPS)} 스냅샷 x {m.sum():,} 입자")

print("\n" + "="*84)
print(f"{'높이 [m]':>10} | {'시간평균 uz [m/s]':^26} | {'uz 표준편차':^21} | {'하강비율 [%]':^14}")
print(f"{'':>10} | {'중립':>12}{'안정':>13} | {'중립':>10}{'안정':>11} | {'중립':>6}{'안정':>8}")
print("="*84)
for lo,hi in BANDS:
    o={}
    for lab in RUNS:
        d=acc[lab]; s=(d['z']>=lo)&(d['z']<hi)
        u=d['uz'][:,s]
        o[lab]=(u.mean(), u.std(), (u<0).mean()*100)
    a,b=o['neutral'],o['stable']
    print(f"{lo:4d}~{hi:<5d}| {a[0]:+12.4f}{b[0]:+13.4f} | {a[1]:10.4f}{b[1]:11.4f} | {a[2]:6.1f}{b[2]:8.1f}")
print("="*84)

# 시각별 z<5 m 평균 uz — 부호가 일관적인지
print("\n=== z<5 m 시간평균 uz, 스냅샷별 (부호 일관성) ===")
for lab in RUNS:
    d=acc[lab]; s=d['z']<5
    per=[d['uz'][i,s].mean() for i in range(len(STEPS))]
    neg=sum(1 for p in per if p<0)
    print(f"  {lab:8s}: " + " ".join(f"{p:+.3f}" for p in per))
    print(f"            음(하강) {neg}/{len(per)} 스냅샷,  전체평균 {np.mean(per):+.4f} m/s")

np.savez('runs/2026-09-21_leipzig-most-stable/analysis/uz_slab.npz',
         z=acc['neutral']['z'], x=acc['neutral']['x'],
         uz_n=acc['neutral']['uz'].mean(axis=0), uz_s=acc['stable']['uz'].mean(axis=0),
         uh_n=acc['neutral']['uh'].mean(axis=0), uh_s=acc['stable']['uh'].mean(axis=0))
print("\n저장: analysis/uz_slab.npz")
