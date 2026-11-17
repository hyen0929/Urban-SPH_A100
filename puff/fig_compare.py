#!/usr/bin/env python3
"""LDM vs 가우시안 플룸 vs SPH 퍼프 시간평균 농도(z=5 m) 그림. validate.py 를 먼저 돌릴 것."""
import os, sys
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from validate import plume_field

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, 'puff', 'data')
rot, tag = int(sys.argv[1]), sys.argv[2] if len(sys.argv) > 2 else '_b0.8_vd0.16'
pick = [int(s) for s in sys.argv[3].split(',')] if len(sys.argv) > 3 else [0, 12, 25, 40]

ref = np.load(os.path.join(D, f'ref_rot{rot}.npz'))
fld = np.load(os.path.join(D, f'field_rot{rot}.npz'))
S = ref['S']; b = float(ref['bin'])
nx, ny = ref['avg'].shape[2:]
xc = ref['x0'] + b * (np.arange(nx) + 0.5); yc = ref['y0'] + b * (np.arange(ny) + 0.5)
sp = np.load(os.path.join(D, f'val_spuff_rot{rot}{tag}.npy'))[0]
avg_t = np.arange(ref['avg_t'][0], ref['avg_t'][1] + 1e-6, 1.0)
H = fld['H']; o = fld['origin']; dx = float(fld['dx'])
ext_b = [o[0] - dx / 2, o[0] + dx * (H.shape[0] - .5), o[1] - dx / 2, o[1] + dx * (H.shape[1] - .5)]
bm = np.ma.masked_where(H.T <= 0, np.ones_like(H.T))

fig, ax = plt.subplots(len(pick), 3, figsize=(15, 4.6 * len(pick)))
norm = LogNorm(1e-6, 1e-2)
for r, k in enumerate(pick):
    sx, sy, sz = S[k, 1:4]
    pl = plume_field(S[k:k + 1, 1:4], xc, yc, 5.0, 2.77, avg_t, S[0, 6])[0]
    for c, (name, F) in enumerate([('LDM (reference)', ref['avg'][0, k]), ('Gaussian plume', pl), ('SPH-field puff', sp[k])]):
        a = ax[r, c]
        a.imshow(bm, origin='lower', extent=ext_b, cmap='Greys', vmin=0, vmax=1.6)
        im = a.imshow(np.ma.masked_less(F.T, 1e-6), origin='lower', norm=norm, cmap='viridis',
                      extent=[xc[0] - b / 2, xc[-1] + b / 2, yc[0] - b / 2, yc[-1] + b / 2])
        a.plot(sx, sy, 'r*', ms=12)
        a.set_xlim(sx - 60, min(sx + 640, 1315)); a.set_ylim(sy - 250, sy + 250); a.set_aspect('equal')
        a.set_title(f'src {int(S[k, 0])}: {name}', fontsize=10)
fig.colorbar(im, ax=ax, shrink=0.4, label='c/Q [s/m$^3$], z=5 m, mean t=250–300 s')
os.makedirs(os.path.join(ROOT, 'puff', 'fig'), exist_ok=True)
out = os.path.join(ROOT, 'puff', 'fig', f'compare_rot{rot}.png')
plt.savefig(out, dpi=80, bbox_inches='tight'); print(out)
