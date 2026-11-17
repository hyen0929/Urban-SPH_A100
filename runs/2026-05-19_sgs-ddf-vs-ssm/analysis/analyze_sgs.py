#!/usr/bin/env python3
# DDF vs SSM SGS-model comparison — CaseAS building-free convective ABL
import struct, re, os, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE='/home/kdh/Urban-SPH'
OUT=os.path.join(BASE,'CaseAS_compare'); os.makedirs(OUT,exist_ok=True)
RUNS={'DDF':'CaseAS_DDF','SSM':'CaseAS_SSM'}
DT=1e-3

def read_vtk(path):
    d=open(path,'rb').read()
    mp=re.search(rb'POINTS[ \t]+(\d+)[ \t]+float\n',d)
    n=int(mp.group(1)); o=mp.end()
    pos=np.frombuffer(d[o:o+12*n],dtype='>f4').reshape(-1,3).astype(np.float64)
    out={'pos':pos,'n':n}
    for name,nc,typ in [('temp',1,'float'),('velocity',3,'float'),
                        ('k_turb',1,'float'),('e_turb',1,'float'),('p_type',1,'int')]:
        pat=name.encode()+rb'[ \t]+'+str(nc).encode()+rb'[ \t]+(\d+)[ \t]+'+typ.encode()+rb'\n'
        m=re.search(pat,d)
        if not m: out[name]=None; continue
        cnt=int(m.group(1)); j=m.end()
        dty='>f4' if typ=='float' else '>i4'
        a=np.frombuffer(d[j:j+cnt*nc*4],dtype=dty).astype(np.float64)
        out[name]=a.reshape(-1,nc) if nc>1 else a
    return out

# z bins
ZED=np.arange(-10,131,5.0); ZC=0.5*(ZED[:-1]+ZED[1:])

def profiles(v):
    """horizontally-averaged profiles over fluid particles."""
    z=v['pos'][:,2]; pt=v['p_type']
    fl=(pt==1)|(pt==-1)
    T=v['temp']; vel=v['velocity']; k=v['k_turb']; e=v['e_turb']
    w=vel[:,2]; spd=np.sqrt((vel**2).sum(1))
    idx=np.digitize(z,ZED)-1
    nb=len(ZC)
    Tp=np.full(nb,np.nan); kp=np.full(nb,np.nan); ep=np.full(nb,np.nan)
    wflux=np.full(nb,np.nan); wvar=np.full(nb,np.nan); spdp=np.full(nb,np.nan)
    for b in range(nb):
        m=fl&(idx==b)
        if m.sum()<20: continue
        Tb=T[m]; wb=w[m]
        Tp[b]=Tb.mean(); kp[b]=k[m].mean(); ep[b]=e[m].mean()
        spdp[b]=spd[m].mean()
        wp=wb-wb.mean(); tp=Tb-Tb.mean()
        wflux[b]=(wp*tp).mean()      # resolved kinematic heat flux <w'T'>
        wvar[b]=(wp*wp).mean()        # vertical velocity variance
    return dict(T=Tp,k=kp,e=ep,wflux=wflux,wvar=wvar,spd=spdp)

# ---- load profile snapshots ----
PSTEPS=[0,2500,5000,10000,15000,20000,25000,30000]
data={r:{} for r in RUNS}
for r,dirn in RUNS.items():
    for s in PSTEPS:
        f=os.path.join(BASE,dirn,f'fluid_condition_{s}stp.vtk')
        if os.path.exists(f): data[r][s]=profiles(read_vtk(f))
print('profile snapshots loaded')

# ---- time series (every 2000 steps) ----
TSTEPS=list(range(0,30001,2000))
ts={r:{'k':[],'wflux':[],'wvar':[],'t':[]} for r in RUNS}
for r,dirn in RUNS.items():
    for s in TSTEPS:
        f=os.path.join(BASE,dirn,f'fluid_condition_{s}stp.vtk')
        if not os.path.exists(f): continue
        p=profiles(read_vtk(f))
        ts[r]['t'].append(s*DT)
        # average over fluid column 0-100 m
        sel=(ZC>=0)&(ZC<=100)
        ts[r]['k'].append(np.nanmean(p['k'][sel]))
        ts[r]['wflux'].append(np.nanmean(p['wflux'][sel]))
        ts[r]['wvar'].append(np.nanmean(p['wvar'][sel]))
print('time series loaded')

C={'DDF':'tab:red','SSM':'tab:blue'}

# === Fig 1: theta(z) profiles ===
fig,ax=plt.subplots(1,2,figsize=(11,6),sharey=True)
for j,s in enumerate([0,5000,10000,20000,30000]):
    for r in RUNS:
        if s not in data[r]: continue
        ls='-' if r=='DDF' else '--'
        ax[0 if r=='DDF' else 1].plot(data[r][s]['T'],ZC,ls,label=f't={s*DT:.0f}s')
for k,r in enumerate(['DDF','SSM']):
    ax[k].set_title(f'{r}: potential temp profile  θ(z)')
    ax[k].set_xlabel('θ [K]'); ax[k].grid(alpha=.3); ax[k].legend(fontsize=8)
ax[0].set_ylabel('z [m]')
plt.tight_layout(); plt.savefig(f'{OUT}/fig1_theta_profiles.png',dpi=110); plt.close()

# === Fig 2: theta(z) overlaid + difference ===
fig,ax=plt.subplots(1,2,figsize=(11,6),sharey=True)
for s in [5000,10000,20000,30000]:
    if s not in data['DDF'] or s not in data['SSM']: continue
    ax[0].plot(data['DDF'][s]['T'],ZC,'-',color=plt.cm.viridis(s/30000),label=f'DDF t={s*DT:.0f}s')
    ax[0].plot(data['SSM'][s]['T'],ZC,'--',color=plt.cm.viridis(s/30000))
    ax[1].plot(data['DDF'][s]['T']-data['SSM'][s]['T'],ZC,'-',
               color=plt.cm.viridis(s/30000),label=f't={s*DT:.0f}s')
ax[0].set_title('θ(z): DDF (solid) vs SSM (dashed)')
ax[0].set_xlabel('θ [K]'); ax[0].set_ylabel('z [m]'); ax[0].grid(alpha=.3); ax[0].legend(fontsize=7)
ax[1].set_title('Δθ(z) = θ(DDF) − θ(SSM)')
ax[1].set_xlabel('Δθ [K]'); ax[1].grid(alpha=.3); ax[1].legend(fontsize=8)
ax[1].axvline(0,color='k',lw=.8)
plt.tight_layout(); plt.savefig(f'{OUT}/fig2_theta_diff.png',dpi=110); plt.close()

# === Fig 3: SGS TKE & dissipation ===
fig,ax=plt.subplots(1,2,figsize=(11,6),sharey=True)
for s in [10000,20000,30000]:
    for r in RUNS:
        if s not in data[r]: continue
        ls='-' if r=='DDF' else '--'
        ax[0].plot(data[r][s]['k'],ZC,ls,color=plt.cm.plasma(s/30000),
                   label=f'{r} t={s*DT:.0f}s')
        ax[1].plot(data[r][s]['e'],ZC,ls,color=plt.cm.plasma(s/30000))
ax[0].set_title('SGS TKE  k_turb(z)'); ax[0].set_xlabel('k [m²/s²]')
ax[0].set_ylabel('z [m]'); ax[0].grid(alpha=.3); ax[0].legend(fontsize=7)
ax[1].set_title('SGS dissipation  e_turb(z)'); ax[1].set_xlabel('ε [m²/s³]')
ax[1].grid(alpha=.3)
plt.tight_layout(); plt.savefig(f'{OUT}/fig3_sgs_tke.png',dpi=110); plt.close()

# === Fig 4: resolved vertical heat flux & w-variance ===
fig,ax=plt.subplots(1,2,figsize=(11,6),sharey=True)
for s in [10000,20000,30000]:
    for r in RUNS:
        if s not in data[r]: continue
        ls='-' if r=='DDF' else '--'
        ax[0].plot(data[r][s]['wflux'],ZC,ls,color=plt.cm.plasma(s/30000),
                   label=f'{r} t={s*DT:.0f}s')
        ax[1].plot(data[r][s]['wvar'],ZC,ls,color=plt.cm.plasma(s/30000))
ax[0].set_title("resolved heat flux  <w'θ'>(z)"); ax[0].set_xlabel("<w'θ'> [K·m/s]")
ax[0].set_ylabel('z [m]'); ax[0].grid(alpha=.3); ax[0].legend(fontsize=7); ax[0].axvline(0,color='k',lw=.8)
ax[1].set_title("vertical velocity variance  <w'²>(z)"); ax[1].set_xlabel("<w'²> [m²/s²]")
ax[1].grid(alpha=.3)
plt.tight_layout(); plt.savefig(f'{OUT}/fig4_heatflux.png',dpi=110); plt.close()

# === Fig 5: time series ===
fig,ax=plt.subplots(1,3,figsize=(15,4.5))
for r in RUNS:
    ax[0].plot(ts[r]['t'],ts[r]['k'],'o-',color=C[r],label=r)
    ax[1].plot(ts[r]['t'],ts[r]['wflux'],'o-',color=C[r],label=r)
    ax[2].plot(ts[r]['t'],ts[r]['wvar'],'o-',color=C[r],label=r)
ax[0].set_title('column-mean SGS TKE (0-100m)'); ax[0].set_ylabel('k [m²/s²]')
ax[1].set_title("column-mean resolved <w'θ'>"); ax[1].set_ylabel("<w'θ'> [K·m/s]")
ax[2].set_title("column-mean <w'²>"); ax[2].set_ylabel("<w'²> [m²/s²]")
for a in ax: a.set_xlabel('time [s]'); a.grid(alpha=.3); a.legend()
plt.tight_layout(); plt.savefig(f'{OUT}/fig5_timeseries.png',dpi=110); plt.close()

# === numeric summary ===
print('\n==== SUMMARY: DDF vs SSM ====')
for s in [10000,20000,30000]:
    if s not in data['DDF'] or s not in data['SSM']: continue
    dT=data['DDF'][s]['T']-data['SSM'][s]['T']
    sel=(ZC>=0)&(ZC<=110)
    print(f' t={s*DT:.0f}s: max|Δθ|={np.nanmax(np.abs(dT[sel])):.3f}K  '
          f'mean k DDF={np.nanmean(data["DDF"][s]["k"][sel]):.4f} '
          f'SSM={np.nanmean(data["SSM"][s]["k"][sel]):.4f}  '
          f'mean<wθ> DDF={np.nanmean(data["DDF"][s]["wflux"][sel]):.4f} '
          f'SSM={np.nanmean(data["SSM"][s]["wflux"][sel]):.4f}')
print('\nplots written to',OUT)
