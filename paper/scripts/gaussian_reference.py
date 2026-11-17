#!/usr/bin/env python3
"""운영 가우시안 모델 대 건물해상 LES — 평지 가정이 도시에서 얼마나 틀리는가.

운영모델 사용자가 실제로 하는 대로 모사한다:
  · 풍속은 기상탑 값(미교란 멱법칙 프로파일 @ 방출고) — 캐노피 감속을 모른다
  · 분산계수는 Briggs (전원/도시, 중립 D)
  · 지형 평평, 건물 없음, 풍향 하나
유한 방출(150~300 s)·유한 추적(300 s)을 반영한 정상플룸 근사:
  노출시간 = clip(T - t0 - x/u, 0, t1-t0)
침적 = Vd * C_ground * 노출시간.  Vd 는 LES 와 같은 식(Wesely/Zhang, u*=0.5, z0=0.1, 중립).
"""
import csv, os
import numpy as np
D=os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..'))
REPO=os.path.normpath(os.path.join(D,'..'))
ns={'__file__':f'{REPO}/runs/lambdap/deposition/dep_metrics.py'}
exec(open(ns['__file__']).read().split('def main(')[0], ns)
read,sp=ns['read_vtk'],ns['spearman']

# ---- 공통 설정 (LES 와 동일) ----
H_SRC=5.0; T0,T1,TEND=150.0,300.0,300.0; CELL=10.0
U_REF, Z_REF, ALPHA = 10.0, 371.25, 0.19
U = U_REF*(H_SRC/Z_REF)**ALPHA                 # 기상탑 외삽 풍속 @ 5 m
KAP,Z0,USTAR = 0.4,0.1,0.5
# Vd (LES 의 DD 식과 동일, 중립)
nu, DB = 1.5e-5, 2.77e-11
Sc=nu/DB; E_B=Sc**(-2/3); tau_p=9e-6
St=tau_p*USTAR**2/nu; E_IM=(St/(St+0.5))**2
Ra=np.log(max(H_SRC,2*Z0)/Z0)/(KAP*USTAR); Ra=max(Ra,1.0)
Rb=1.0/(USTAR*3.0*(E_B+E_IM)+1e-20)
VD=1.0/(Ra+Rb)+8.9e-5
print(f'기상탑 풍속 U(5 m) = {U:.2f} m/s,  Vd = {VD:.3e} m/s  (Ra={Ra:.1f}, Rb={Rb:.1f})')

def sigmas(x, kind):
    x=np.maximum(x,1.0)
    if kind=='urban':   return 0.16*x*(1+0.0004*x)**-0.5, 0.14*x*(1+0.0003*x)**-0.5
    else:               return 0.08*x*(1+0.0001*x)**-0.5, 0.06*x*(1+0.0015*x)**-0.5

def depletion(xg, kind):
    """선원 고갈 Q(x)/Q0 = exp(-int Vd/(sqrt(2pi) u sz) dx).  표준 source-depletion."""
    sy_,sz_=sigmas(xg,kind)
    integ=np.cumsum(VD/(np.sqrt(2*np.pi)*U*sz_))*np.gradient(xg)[0]
    return np.exp(-integ)

def gaussian_dep(sx, sy, kind, lim, cell=CELL):
    """소스 (sx,sy) 에서 바람 +x. 격자별 침적 비율(총 방출량 대비, 합이 <=1)."""
    gx=np.arange(lim[0],lim[1],cell)+cell/2; gy=np.arange(lim[2],lim[3],cell)+cell/2
    X,Y=np.meshgrid(gx-sx, gy-sy, indexing='ij')
    out=np.zeros_like(X); ok=X>1.0
    xs=np.maximum(gx-sx,1.0); fdep=depletion(xs,kind)          # x 에만 의존
    FD=np.broadcast_to(fdep[:,None],X.shape)[ok]
    sy_,sz_=sigmas(X[ok],kind)
    texp=np.clip(TEND-T0-X[ok]/U, 0.0, T1-T0)                  # 이 지점의 노출 시간 [s]
    # Q = 1/(T1-T0) [방출률], 총 방출량 1 로 정규화
    Cg=FD/( (T1-T0)*np.pi*U*sy_*sz_ )*np.exp(-Y[ok]**2/(2*sy_**2))*np.exp(-H_SRC**2/(2*sz_**2))
    out[ok]=VD*Cg*texp*cell*cell
    return out, gx, gy

# ---- LES 침적 지도 (50 소스 × 4 풍향, 벽면 off) ----
CX=657.5; ANG=[0,90,180,270]
S=[(int(r['id']),float(r['x']),float(r['y'])) for r in csv.DictReader(open(f'{D}/data/00_setup/sources_n50.csv'))]
ids=[s[0] for s in S]
def rot(x,y,d):
    t=np.radians(d);c,s=np.cos(t),np.sin(t);dx,dy=x-CX,y
    return CX+dx*c-dy*s, dx*s+dy*c
LIM=[0.0,1320.0,-660.0,660.0]
rows=[]
print(f"\n{'각도':>5} {'소스당 침적률':>24}   {'공간 상관(log)':>22}")
print(f"{'':>5} {'LES':>8}{'도시G':>8}{'전원G':>8}   {'도시G':>10}{'전원G':>10}")
agg={k:[] for k in ('les','urb','rur','cu','cr')}
for a in ANG:
    p,F=read(f'{REPO}/runs/2026-10-06_leipzig-src50-rot{a}/out/LDM_30000stp.vtk')
    pt=F['p_type'].astype(int)-100; v=np.linalg.norm(F['velocity'],axis=1); dep=v==0
    les_f=[];urb_f=[];rur_f=[];cu=[];cr=[]
    for sid,x,y in S:
        m=pt==sid; n=m.sum()
        rx,ry=rot(x,y,a)
        d=p[m&dep][:,:2]
        les_f.append(float((m&dep).sum()/max(n,1)))
        Gu,gx,gy=gaussian_dep(rx,ry,'urban',LIM); Gr,_,_=gaussian_dep(rx,ry,'rural',LIM)
        urb_f.append(float(Gu.sum())); rur_f.append(float(Gr.sum()))
        # LES 를 같은 격자에 올려 공간 패턴 비교
        ix=np.clip(((d[:,0]-LIM[0])/CELL).astype(int),0,Gu.shape[0]-1)
        iy=np.clip(((d[:,1]-LIM[2])/CELL).astype(int),0,Gu.shape[1]-1)
        L=np.zeros_like(Gu); np.add.at(L,(ix,iy),1.0/max(n,1))
        sel=(L>0)|(Gu>1e-12)
        cu.append(sp(np.log10(L[sel]+1e-9),np.log10(Gu[sel]+1e-9)))
        cr.append(sp(np.log10(L[sel]+1e-9),np.log10(Gr[sel]+1e-9)))
        rows.append(dict(deg=a,src=sid,les=les_f[-1],gauss_urban=urb_f[-1],gauss_rural=rur_f[-1],
                         ratio_urban=urb_f[-1]/max(les_f[-1],1e-9),ratio_rural=rur_f[-1]/max(les_f[-1],1e-9),
                         sp_urban=cu[-1],sp_rural=cr[-1]))
    for k,vv in (('les',les_f),('urb',urb_f),('rur',rur_f),('cu',cu),('cr',cr)): agg[k].append(np.mean(vv))
    print(f'{a:4d}° {np.mean(les_f)*100:7.2f}%{np.mean(urb_f)*100:7.2f}%{np.mean(rur_f)*100:7.2f}%   '
          f'{np.mean(cu):10.2f}{np.mean(cr):10.2f}')
with open(f'{D}/data/08_gaussian/gaussian_vs_les.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
L=np.array([r['les'] for r in rows]); Uu=np.array([r['gauss_urban'] for r in rows]); Rr=np.array([r['gauss_rural'] for r in rows])
print(f'\n전체 (50소스 × 4풍향 = {len(rows)} 사례)')
print(f'  침적률  LES {100*L.mean():.2f}%   도시가우시안 {100*Uu.mean():.2f}%   전원가우시안 {100*Rr.mean():.2f}%')
print(f'  배율    도시G/LES 중앙값 {np.median(Uu/L):.2f}  (10~90퍼센타일 {np.percentile(Uu/L,10):.2f}~{np.percentile(Uu/L,90):.2f})')
print(f'          전원G/LES 중앙값 {np.median(Rr/L):.2f}  (10~90퍼센타일 {np.percentile(Rr/L,10):.2f}~{np.percentile(Rr/L,90):.2f})')
print(f'  소스 순위 재현  ρs(LES, 도시G) {sp(L,Uu):+.2f}   ρs(LES, 전원G) {sp(L,Rr):+.2f}')
