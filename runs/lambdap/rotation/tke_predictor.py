#!/usr/bin/env python3
"""국소 난류에너지가 침적량을 풍향에 무관하게 예측하는가.

λp 는 회전 불변이라 풍향 정보를 담지 못한다. 반면 국소 k 는 그 방향의 유동에서
계산되므로 풍향을 반영한다. 네 방향 모두에서 k 가 f_dep 을 예측하면,
"형태가 아니라 유동량을 쓰라"는 처방이 성립한다.

k 는 각 방향 Run A 의 t=300 s 전체장에서 소스 반경 50 m, z<10 m 유체 입자 평균.
"""
import csv, os
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
DEP=os.path.normpath(os.path.join(HERE,'..','deposition'))
MORPH=os.path.normpath(os.path.join(HERE,'..'))
RUNS=os.path.normpath(os.path.join(HERE,'..','..'))
ns={'__file__':os.path.join(DEP,'dep_metrics.py')}
exec(open(ns['__file__']).read().split('def main(')[0], ns)
read_vtk,load,sp=ns['read_vtk'],ns['load'],ns['spearman']
CX=657.5; COURT={3,6,10}
RUN={0:'2026-09-15_leipzig-dt01-sponge-uy',90:'2026-09-23_leipzig-rot90',
     180:'2026-09-23_leipzig-rot180',270:'2026-09-26_leipzig-rot270'}
S={int(r['id']):(float(r['x']),float(r['y'])) for r in csv.DictReader(open(f'{MORPH}/sources13.csv'))}
lam={int(float(r['R_m'])):{i:float(r[f'S{i}']) for i in range(1,14)}
     for r in csv.DictReader(open(f'{MORPH}/lambda_p_by_R.csv'))}
ids=sorted(S); keep=np.array([i not in COURT for i in ids])
x75=np.array([lam[75][i] for i in ids])
def rot(x,y,d):
    t=np.radians(d);c,s=np.cos(t),np.sin(t);dx,dy=x-CX,y
    return CX+dx*c-dy*s, dx*s+dy*c

rows=[]; K={}; FD={}
for deg,rd in RUN.items():
    p,src,dep=load(f'{RUNS}/{rd}')
    FD[deg]=np.array([float(((src==i)&dep).sum()/max((src==i).sum(),1)) for i in ids])
    print(f'{deg:4d}° 전체장 읽는 중...', flush=True)
    q,F=read_vtk(f'{RUNS}/{rd}/out/fluid_30000stp.vtk')
    m=(F['p_type'].astype(int)==1)&(q[:,2]>0)&(q[:,2]<10.0)
    qq,kk=q[m][:,:2],F['k_turb'][m]
    K[deg]=[]
    for i in ids:
        sx,sy=rot(*S[i],deg)
        sel=(np.abs(qq[:,0]-sx)<50)&(np.abs(qq[:,1]-sy)<50)
        sel&=((qq[:,0]-sx)**2+(qq[:,1]-sy)**2)<2500
        K[deg].append(float(kk[sel].mean()) if sel.sum() else np.nan)
    K[deg]=np.array(K[deg])
    for j,i in enumerate(ids):
        rows.append(dict(deg=deg,src=i,k50=K[deg][j],f_dep=FD[deg][j],lam75=lam[75][i]))

print(f'\nf_dep 예측력 비교 (중정 제외 n={keep.sum()})')
print(f"{'각도':>5} {'ρs(λp75, f_dep)':>18} {'ρs(k50, f_dep)':>16}   k50 범위")
a,b=[],[]
for deg in RUN:
    ra=sp(x75[keep],FD[deg][keep]); rb=sp(K[deg][keep],FD[deg][keep])
    a.append(ra); b.append(rb)
    print(f'{deg:4d}° {ra:18.2f} {rb:16.2f}   {K[deg][keep].min():.3f}~{K[deg][keep].max():.3f}')
print(f"{'폭':>5} {max(a)-min(a):18.2f} {max(b)-min(b):16.2f}")
print(f"{'평균':>5} {np.mean(a):18.2f} {np.mean(b):16.2f}")

# k 는 풍향에 따라 실제로 변하는가 (λp 는 불변)
print('\nk50 의 풍향 의존성 — 소스별 (최대/최소)')
r=np.array([K[d] for d in RUN])
print(f'  소스별 k50 변동비 중앙값 {np.median(r.max(0)/r.min(0)):.2f}배   최대 {np.nanmax(r.max(0)/r.min(0)):.2f}배')
with open(os.path.join(HERE,'tke_predictor.csv'),'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print('저장: tke_predictor.csv')
