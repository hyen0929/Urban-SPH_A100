#!/usr/bin/env python3
"""기하 예측의 잔차를 유동량으로 설명할 수 있는가.

분배비(벽/지면)는 기하(벽면적/지면면적)로 ρs=+0.76 까지 설명된다. 남는 것은
소스별 3~5배 편차이고 풍향쌍 간 상관이 +0.10 이라 '풍향이 만드는 것'이다.
그 잔차 = log(sim) - log(geo) 를 방향 의존 유동량으로 회귀한다.
"""
import csv, os, itertools
import numpy as np
D=os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..'))
REPO=os.path.normpath(os.path.join(D,'..'))
ns={'__file__':f'{REPO}/runs/lambdap/deposition/dep_metrics.py'}
exec(open(ns['__file__']).read().split('def main(')[0], ns)
read,sp=ns['read_vtk'],ns['spearman']
ANG=[0,90,180,270]; CX=657.5; R=75.0; RK=50.0
hf=np.load(f'{D}/data/00_setup/heightfield.npz')
H,cell,x0,y0=hf['H'],float(hf['cell']),float(hf['x0']),float(hf['y0'])
ny,nx=H.shape; B=H>1.0; TALL=H>5.0
gx=x0+(np.arange(nx)+0.5)*cell; gy=y0+(np.arange(ny)+0.5)*cell
GX,GY=np.meshgrid(gx,gy)
WA=np.zeros_like(H)
dxx=np.abs(np.diff(H,axis=1))*cell; WA[:,:-1]+=dxx/2; WA[:,1:]+=dxx/2
dyy=np.abs(np.diff(H,axis=0))*cell; WA[:-1,:]+=dyy/2; WA[1:,:]+=dyy/2
S=[(int(r['id']),float(r['x']),float(r['y'])) for r in csv.DictReader(open(f'{D}/data/00_setup/sources_n50.csv'))]
ids=[s[0] for s in S]
def rot(x,y,d):
    t=np.radians(d);c,s=np.cos(t),np.sin(t);dx,dy=x-CX,y
    return CX+dx*c-dy*s, dx*s+dy*c
# --- 기하 예측 (방향 무관) ---
geo=[]
for _,x,y in S:
    m=((GX-x)**2+(GY-y)**2)<=R*R
    geo.append(WA[m].sum()/max((~B[m]).sum()*cell*cell,1e-9))
geo=np.array(geo)
# --- 방향 의존 기하: 풍상 페치 / 풍하 트임 (±30° 부채꼴) ---
MAXR,NS,SECT=400.0,320,30.0
def ray(sx,sy,deg):
    t=np.radians(deg);dx,dy=np.cos(t),np.sin(t)
    px,py=(sx-x0)/cell,(sy-y0)/cell
    for k in range(1,NS+1):
        s=k*(MAXR/NS)/cell
        qx,qy=int(round(px+dx*s)),int(round(py+dy*s))
        if qx<0 or qx>=nx or qy<0 or qy>=ny: return MAXR
        if TALL[qy,qx]: return k*(MAXR/NS)
    return MAXR
def sector(sx,sy,deg):
    return float(np.mean([ray(sx,sy,deg+d) for d in np.arange(-SECT,SECT+1e-9,5.0)]))
ROT2WIND={0:0.0,90:270.0,180:180.0,270:90.0}
# --- 시뮬 분배비 + 유동량 ---
sim={}; K50={}; U50={}; KW={}; UP={}; DN={}
for a in ANG:
    _,F=read(f'{REPO}/runs/2026-10-07_leipzig-src50w-rot{a}/out/LDM_30000stp.vtk')
    pt=F['p_type'].astype(int)-100; b=F['buffer_type'].astype(int)
    g=np.array([float(((pt==i)&(b==1)).sum()/max((pt==i).sum(),1)) for i in ids])
    w=np.array([float(((pt==i)&(b==2)).sum()/max((pt==i).sum(),1)) for i in ids])
    sim[a]=w/np.maximum(g,1e-9)
    print(f'{a:4d}° 전체장 읽는 중...',flush=True)
    q,FF=read(f'{REPO}/runs/'+{0:'2026-09-15_leipzig-dt01-sponge-uy',90:'2026-09-23_leipzig-rot90',
              180:'2026-09-23_leipzig-rot180',270:'2026-09-26_leipzig-rot270'}[a]+'/out/fluid_30000stp.vtk')
    m=(FF['p_type'].astype(int)==1)&(q[:,2]>0)&(q[:,2]<10.0)
    P=q[m][:,:2]; KK=FF['k_turb'][m]; UU=np.hypot(FF['velocity'][m][:,0],FF['velocity'][m][:,1])
    # 건물 인접 여부 (해당 방향 좌표계에서)
    jx=np.clip(np.floor((P[:,0]-x0)/cell).astype(int),0,nx-1)
    jy=np.clip(np.floor((P[:,1]-y0)/cell).astype(int),0,ny-1)
    near=np.zeros(len(P),bool)
    for ddx in (-2,-1,0,1,2):
        for ddy in (-2,-1,0,1,2):
            near |= TALL[np.clip(jy+ddy,0,ny-1), np.clip(jx+ddx,0,nx-1)]
    k5=[];u5=[];kw=[];up=[];dn=[]
    for sid,x,y in S:
        rx,ry=rot(x,y,a)
        sel=((P[:,0]-rx)**2+(P[:,1]-ry)**2)<RK*RK
        k5.append(float(KK[sel].mean()) if sel.sum() else np.nan)
        u5.append(float(UU[sel].mean()) if sel.sum() else np.nan)
        sw=sel&near
        kw.append(float(KK[sw].mean()) if sw.sum() else np.nan)
        wdeg=ROT2WIND[a]
        up.append(sector(x,y,(wdeg+180)%360)); dn.append(sector(x,y,wdeg))
    K50[a]=np.array(k5); U50[a]=np.array(u5); KW[a]=np.array(kw)
    UP[a]=np.array(up); DN[a]=np.array(dn)

res={a: np.log(np.maximum(sim[a],1e-6))-np.log(np.maximum(geo,1e-6)) for a in ANG}
print('\n잔차 = log(시뮬 분배비) − log(기하 분배비)')
print(f"{'각도':>5} {'중앙값':>8} {'표준편차':>9}")
for a in ANG: print(f'{a:4d}° {np.median(res[a]):+8.3f} {res[a].std():9.3f}')
cc=[sp(res[p],res[q]) for p,q in itertools.combinations(ANG,2)]
print(f'풍향쌍 간 잔차 상관 평균 {np.mean(cc):+.3f}  → 잔차는 형상 고정이 아니다')

print('\n잔차를 설명하는가 — ρs(예측자, 잔차)')
PRED=[('국소 k (r<50 m, z<10 m)',K50),('국소 |u_h|',U50),
      ('건물인접 k',KW),('풍상 페치',UP),('풍하 트임',DN)]
print(f"{'예측자':>24} | " + '  '.join(f'{a:>6}°' for a in ANG) + f"  {'평균':>7}")
for lab,P_ in PRED:
    v=[sp(P_[a][~np.isnan(P_[a])], res[a][~np.isnan(P_[a])]) for a in ANG]
    print(f'{lab:>24} | ' + '  '.join(f'{x:+7.2f}' for x in v) + f'  {np.mean(v):+7.2f}')
