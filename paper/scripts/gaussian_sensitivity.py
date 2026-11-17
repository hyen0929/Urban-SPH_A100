#!/usr/bin/env python3
"""운영 가우시안 모델의 설정 민감도 — 어떤 설정에서도 소스 간 구분이 되는가.

축 셋: 풍속(기상탑 외삽 4.41 / 캐노피 실측 2.65 m/s)
       침적속도(LES 식 0.031 / 운영 통상 0.001 / 중간 0.005 m/s)
       안정도(Briggs 도시·전원 × C, D, E)
요점: 평지 가우시안은 소스 위치를 구분할 수단이 없다(도메인 가장자리 효과만 남는다).
      LES 의 소스 간 분산과 비교한다.
"""
import csv, os, itertools
import numpy as np
D=os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..'))
REPO=os.path.normpath(os.path.join(D,'..'))
ns={'__file__':f'{REPO}/runs/lambdap/deposition/dep_metrics.py'}
exec(open(ns['__file__']).read().split('def main(')[0], ns)
read,sp=ns['read_vtk'],ns['spearman']
H_SRC=5.0; T0,T1,TEND=150.0,300.0,300.0; CELL=10.0; LIM=[0.0,1320.0,-660.0,660.0]
CX=657.5; ANG=[0,90,180,270]
SIG={  # Briggs (x in m)
 ('urban','C'):(lambda x:0.22*x*(1+0.0004*x)**-0.5, lambda x:0.20*x),
 ('urban','D'):(lambda x:0.16*x*(1+0.0004*x)**-0.5, lambda x:0.14*x*(1+0.0003*x)**-0.5),
 ('urban','E'):(lambda x:0.11*x*(1+0.0004*x)**-0.5, lambda x:0.08*x*(1+0.0015*x)**-0.5),
 ('rural','C'):(lambda x:0.11*x*(1+0.0001*x)**-0.5, lambda x:0.08*x*(1+0.0002*x)**-0.5),
 ('rural','D'):(lambda x:0.08*x*(1+0.0001*x)**-0.5, lambda x:0.06*x*(1+0.0015*x)**-0.5),
 ('rural','E'):(lambda x:0.06*x*(1+0.0001*x)**-0.5, lambda x:0.03*x*(1+0.0003*x)**-1.0)}
WINDS=[('기상탑 4.41',4.41),('캐노피 2.65',2.65)]
VDS=[('LES식 0.031',0.0313),('중간 0.005',0.005),('운영통상 0.001',0.001)]

def gdep(sx,sy,key,U,VD):
    fy,fz=SIG[key]
    gx=np.arange(LIM[0],LIM[1],CELL)+CELL/2; gy=np.arange(LIM[2],LIM[3],CELL)+CELL/2
    X,Y=np.meshgrid(gx-sx,gy-sy,indexing='ij'); out=np.zeros_like(X); ok=X>1.0
    xs=np.maximum(gx-sx,1.0)
    integ=np.cumsum(VD/(np.sqrt(2*np.pi)*U*fz(xs)))*np.gradient(xs)[0]
    FD=np.broadcast_to(np.exp(-integ)[:,None],X.shape)[ok]
    sy_,sz_=fy(X[ok]),fz(X[ok])
    texp=np.clip(TEND-T0-X[ok]/U,0.0,T1-T0)
    Cg=FD/((T1-T0)*np.pi*U*sy_*sz_)*np.exp(-Y[ok]**2/(2*sy_**2))*np.exp(-H_SRC**2/(2*sz_**2))
    out[ok]=VD*Cg*texp*CELL*CELL
    return out.sum()

S=[(int(r['id']),float(r['x']),float(r['y'])) for r in csv.DictReader(open(f'{D}/data/00_setup/sources_n50.csv'))]
def rot(x,y,d):
    t=np.radians(d);c,s=np.cos(t),np.sin(t);dx,dy=x-CX,y
    return CX+dx*c-dy*s, dx*s+dy*c
LES={}
for a in ANG:
    p,F=read(f'{REPO}/runs/2026-10-06_leipzig-src50-rot{a}/out/LDM_30000stp.vtk')
    pt=F['p_type'].astype(int)-100; dep=np.linalg.norm(F['velocity'],axis=1)==0
    LES[a]=np.array([float(((pt==i)&dep).sum()/max((pt==i).sum(),1)) for i,_,_ in S])
lesall=np.concatenate([LES[a] for a in ANG])
print(f'LES  소스당 침적률 평균 {100*lesall.mean():.2f}%  소스 간 변동계수 CV = {lesall.std()/lesall.mean():.3f}')
print(f'     소스 간 최대/최소 = {lesall.max()/lesall.min():.1f}배\n')
rows=[]
print(f"{'분산계수':>10} {'풍속':>12} {'Vd':>14} | {'침적률':>8} {'G/LES':>7} {'소스CV':>8} {'ρs(LES,G)':>10}")
for (kind,cls),(wl,U),(vl,VD) in itertools.product(SIG,WINDS,VDS):
    g=[]; 
    for a in ANG:
        g.append(np.array([gdep(*rot(x,y,a),(kind,cls),U,VD) for _,x,y in S]))
    gall=np.concatenate(g)
    cv=gall.std()/gall.mean()
    r=sp(lesall,gall)
    rows.append(dict(sigma=f'{kind}-{cls}',wind=wl,vd=vl,dep_pct=round(100*gall.mean(),2),
                     ratio=round(gall.mean()/lesall.mean(),3),cv=round(cv,4),rho=round(r,3)))
    print(f'{kind+"-"+cls:>10} {wl:>12} {vl:>14} | {100*gall.mean():7.2f}% {gall.mean()/lesall.mean():7.2f} {cv:8.4f} {r:+10.2f}')
with open(f'{D}/data/08_gaussian/gaussian_sensitivity.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
cvs=[r['cv'] for r in rows]; rr=[abs(r['rho']) for r in rows]; ra=[r['ratio'] for r in rows]
print(f'\n전체 {len(rows)} 설정')
print(f'  G/LES 배율      {min(ra):.2f} ~ {max(ra):.2f}   (설정에 따라 {max(ra)/min(ra):.0f}배 변한다)')
print(f'  가우시안 소스 CV {min(cvs):.4f} ~ {max(cvs):.4f}   vs LES CV {lesall.std()/lesall.mean():.3f}')
print(f'  |ρs(LES, G)|     {min(rr):.2f} ~ {max(rr):.2f}')
