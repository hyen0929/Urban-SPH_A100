import sys, numpy as np
sys.path.insert(0,'/tmp/claude-1014/-home-kdh-Urban-SPH/90fe7215-f210-4eca-bb99-e4c2226c7eef/scratchpad')
from vtkread import read
STEPS=list(range(14000,15001,100))
R={'neutral':'runs/2026-09-15_leipzig-dt01-sponge-uy/out/fluid_condition_%dstp.vtk',
   'stable' :'runs/2026-09-21_leipzig-most-stable/out/fluid_condition_%dstp.vtk'}
prof={}
for lab,tpl in R.items():
    S={}
    for s in STEPS:
        n,pts,fl=read(tpl%s)
        y=np.array(pts[1::3]); z=np.array(pts[2::3]); pt=np.array(fl['p_type'])
        m=(pt==1)&(np.abs(y)<1e-3)
        zz=z[m]; uz=np.array(fl['velocity'])[2::3][m]
        key=np.round(zz*4)/4
        for k in np.unique(key):
            S.setdefault(k,[]).append(uz[key==k].mean())
    prof[lab]={k:np.mean(v) for k,v in S.items()}
    print(f"{lab}: z 평면 {len(S)}개")

zs=sorted(k for k in prof['neutral'] if 0<k<40)
print("\n=== 평면별 시간평균 uz [m/s] (y=0, t=140~150 s) ===")
print(f"{'z [m]':>8} | {'중립':>10} | {'안정':>10} | {'차이':>10}")
print("-"*46)
prev=None
flips=0
for k in zs:
    a=prof['neutral'][k]; b=prof['stable'][k]
    mark=""
    if prev is not None and np.sign(b)!=np.sign(prev) and abs(b)>0.005 and abs(prev)>0.005:
        mark=" <- 부호전환"; flips+=1
    prev=b
    print(f"{k:8.2f} | {a:+10.4f} | {b:+10.4f} | {b-a:+10.4f}{mark}")
print("-"*46)
print(f"안정 런 z<40 m 구간 부호전환 {flips}회 / 평면 {len(zs)}개")
print(f"  (매 평면 교대 = 체커보드, 1~2회 = 물리적 층 구조)")
