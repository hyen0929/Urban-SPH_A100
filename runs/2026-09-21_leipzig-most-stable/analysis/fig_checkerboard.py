"""y=0 대조가 아니라 도메인 전체 수평평균 연직속도 프로파일.
중립(부력 off) 대 안정(부력 on) — 2dz 격자 진동의 유무를 보인다."""
import numpy as np, struct, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

RUNS = [('중립  (MOST off, 부력 off)', 'runs/2026-09-15_leipzig-dt01-sponge-uy/out/restart.bin', '#2a78d6', 'o'),
        ('안정  (MOST on,  부력 on)',  'runs/2026-09-21_leipzig-most-stable/out/restart.bin',     '#eb6834', 's')]
ZMAX = 40.0

def profile(path):
    raw = open(path,'rb').read(136)
    sz,_,_,_,n2,_,_,_,step = struct.unpack('<9i', raw[16:52]); t,=struct.unpack('<f', raw[52:56])
    mm = np.memmap(path, dtype=np.uint8, mode='r', offset=136, shape=(n2,sz))
    pt = mm[:,:16].copy().view(np.int32)[:,3]
    fl = mm[:,16:48].copy().view(np.float32)
    f  = pt==1
    z, uz = fl[:,2][f], fl[:,5][f]
    key = np.round(z*4).astype(np.int64)
    zs, us = [], []
    for k in np.unique(key):
        zz = k/4.0
        if not (0 < zz < ZMAX): continue
        zs.append(zz); us.append(uz[key==k].mean())
    return np.array(zs), np.array(us), step, t

def pair_avg(z,u):
    n=(len(u)//2)*2
    return z[:n].reshape(-1,2).mean(1), u[:n].reshape(-1,2).mean(1)

plt.rcParams.update({'font.family':'Noto Sans CJK JP','font.size':9,
    'axes.labelsize':9.5,'axes.titlesize':10,'xtick.labelsize':8.5,'ytick.labelsize':8.5,
    'legend.fontsize':8.5,'axes.edgecolor':'#9a9a97','axes.linewidth':0.8,
    'xtick.color':'#52514e','ytick.color':'#52514e','figure.dpi':150,
    'savefig.facecolor':'#fcfcfb','figure.facecolor':'#fcfcfb','axes.facecolor':'#fcfcfb'})
INK, MUTED = '#0b0b0b', '#8a8a86'

fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.6), sharey=True)
data=[]
for lab, path, col, mk in RUNS:
    z,u,step,t = profile(path); data.append((lab,z,u,col,mk,step,t))

for ax, mode in zip(axes, ('raw','filt')):
    ax.axvline(0, color=MUTED, lw=0.8, zorder=1)
    ax.grid(True, color='#e3e3df', lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    for lab,z,u,col,mk,step,t in data:
        zz,uu = (z,u) if mode=='raw' else pair_avg(z,u)
        ax.plot(uu, zz, color=col, lw=2.0, marker=mk, ms=5.5, mfc=col, mec='#fcfcfb',
                mew=0.9, zorder=3, label=lab, clip_on=False)
    ax.set_xlabel('수평평균 $u_z$  [m/s]', color=INK)
    ax.set_ylim(0, ZMAX); ax.set_xlim(-0.09, 0.09)
    for s in ('top','right'): ax.spines[s].set_visible(False)

axes[0].set_ylabel('높이 $z$  [m]', color=INK)
axes[0].set_title('(a) 격자 평면별 — 안정 런은 평면마다 부호가 교대', color=INK, pad=9, loc='left')
axes[1].set_title('(b) 인접 평면 쌍평균 — 2Δz 모드 제거 후', color=INK, pad=9, loc='left')
axes[0].legend(loc='upper left', frameon=False, labelcolor=INK, handlelength=2.4,
                bbox_to_anchor=(0.02, 0.99))

axes[0].annotate('2Δz 진동\n진폭 0.033 m/s', xy=(-0.072,1.25), xytext=(-0.086,6.8),
    color=INK, fontsize=8.2, ha='left',
    arrowprops=dict(arrowstyle='-|>', color=MUTED, lw=0.9, shrinkA=2, shrinkB=3))
axes[1].annotate('걸러내면 양쪽 모두\n|$u_z$| ≲ 0.03 m/s', xy=(0.005,20), xytext=(-0.082,30),
    color=INK, fontsize=8.2, ha='left')

st = data[1]
fig.suptitle('부력항이 여기시킨 2Δz 연직 격자 진동  (도메인 전체, 평면당 277,729 입자)',
             color=INK, fontsize=10.5, x=0.012, ha='left', y=0.985)
fig.text(0.012, 0.015,
    f"중립 step {data[0][5]} (t={data[0][6]:.0f}s) · 안정 step {st[5]} (t={st[6]:.0f}s)   |   "
    f"부호전환 1/15 대 13/15,  진동/실제신호 0.1배 대 5.3배",
    color='#52514e', fontsize=7.8, ha='left')
fig.tight_layout(rect=[0,0.035,1,0.945])
out='runs/2026-09-21_leipzig-most-stable/analysis/f_checkerboard.png'
fig.savefig(out, bbox_inches='tight')
print("저장:", out)
np.savez('runs/2026-09-21_leipzig-most-stable/analysis/uz_profile_full.npz',
         z_n=data[0][1],uz_n=data[0][2],z_s=data[1][1],uz_s=data[1][2])
