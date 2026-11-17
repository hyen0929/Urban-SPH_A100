#!/usr/bin/env python3
"""
Flag fluid particles inside buildings/terrain as p_type=99 using a column-max
heightfield built from the surface (IBM, p_type>=1000) particles.

    interior(x,y,z)  <=>  z <= z_top(x,y),   z_top = max z of surface particles in the (x,y) bin

Needs only the surface particle positions — no STL, no normals, no watertight
mesh, no flood fill. Overlapping meshes and open bottoms are harmless (max is
idempotent; only the top sheet matters). Surface particles need not lie on a
lattice; bins are centred on the fluid lattice.

The only geometry this gets wrong is a true overhang (air between two surface
sheets in the same column) — it is treated as solid. Such columns are counted
and reported ("overhang candidates"); floor faces (z < ground+floor_z) are
excluded from that test, otherwise every closed-bottom building trips it.

Do NOT enable hole filling on well-sampled surfaces: on Leipzig_surface a
5-of-8 neighbour fill cut concave corners/courtyards (10 real holes out of 623).

Usage:
    python flag_interior_heightfield.py input/<CASE>/<CASE>.txt [options]
      --out PATH        output file (default: <input>_int.txt; never overwrites the input)
      --bin B           (x,y) bin width [m] (default: fluid spacing h/h_coeff)
      --floor-z DZ      exclude surface points below ground+DZ from the overhang test (default 0.5)
      --gap G           vertical hole size that flags an overhang candidate (default 2*bin)
      --hole-fill N     iterations of 5-of-8 neighbour fill of empty bins (default 0; see above)
      --vtk             also write <out>.interior.vtk with the flagged particles
      --reset-99        re-evaluate rows already p_type=99 (default: keep them)
      --dry-run         diagnostics only, write nothing

The CUDA side already maps p_type=99 -> buffer_type=99 at init (functions_NNPS.cuh).
Requires: numpy, awk (the rewrite pass is done by awk for speed on multi-GB inputs).
"""
import argparse, os, subprocess, sys, time
from collections import defaultdict
import numpy as np

H_COEFF = 1.5   # Parameters.cuh: dx = h / h_coeff
SAMPLE_EVERY = 1000   # fluid rows used for the awk-vs-Python consistency check


def read_header(path):
    with open(path) as f:
        cols = f.readline().split()
    pos = {c: i for i, c in enumerate(cols)}          # column index (as in file) -> field position
    for c in ('1', '2', '3', '8', '9'):
        if c not in pos:
            sys.exit(f"header lacks column {c} (need 1=x 2=y 3=z 8=p_type 9=h)")
    return pos


def candidate_expr(fp, reset):
    """awk test for rows that get (re)classified: fluid, plus old p_type=99 when --reset-99."""
    return f'(${fp}==1 || ${fp}==99)' if reset else f'${fp}==1'


def pass1_extract(path, pos, tmp_surface, reset):
    """awk: surface rows -> tmp_surface (x y z), plus fluid anchor/spacing and counts."""
    fx, fy, fz, fp, fh = (pos[c] + 1 for c in ('1', '2', '3', '8', '9'))
    cand = candidate_expr(fp, reset)
    # every SAMPLE_EVERY-th candidate row is also returned: pass 2 must flag exactly the same
    # subset that the Python binning flags (guards the awk/Python hand-off).
    prog = (
        f'NR>1{{c[${fp}]++}} '
        f'NR>1 && ${fp}>=1000{{print ${fx}"\\t"${fy}"\\t"${fz} > "{tmp_surface}"}} '
        f'NR>1 && {cand}{{nf++; if(nf%{SAMPLE_EVERY}==1) print "SAMPLE",${fx},${fy},${fz}; '
        f'if(!a){{print "ANCHOR",${fx},${fy},${fz},${fh}; a=1}}}} '
        f'END{{for(k in c) print "COUNT",k,c[k]}}'
    )
    out = subprocess.run(['awk', '-F', '\t', prog, path], capture_output=True, text=True, check=True,
                         env=dict(os.environ, LC_ALL='C')).stdout
    anchor, counts, sample = None, {}, []
    for line in out.splitlines():
        t = line.split()
        if t[0] == 'ANCHOR':
            anchor = tuple(map(float, t[1:]))
        elif t[0] == 'COUNT':
            counts[int(t[1])] = int(t[2])
        elif t[0] == 'SAMPLE':
            sample.append(tuple(map(float, t[1:])))
    return anchor, counts, np.array(sample).reshape(-1, 3)


def nn_spacing(P, n_sample=5000, cell=2.0):
    key = np.floor(P / cell).astype(int)
    H = defaultdict(list)
    for i, k in enumerate(map(tuple, key)):
        H[k].append(i)
    rng = np.random.default_rng(0)
    idx = rng.choice(len(P), min(n_sample, len(P)), replace=False)
    d = []
    for i in idx:
        kx, ky, kz = key[i]
        cand = [j for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1)
                for j in H.get((kx + a, ky + b, kz + c), []) if j != i]
        d.append(np.sqrt(((P[cand] - P[i]) ** 2).sum(1)).min() if cand else np.inf)
    d = np.array(d)
    return d[np.isfinite(d)]


def build_ztop(S, ox, oy, b, hole_fill):
    si = np.floor((S[:, 0] - ox) / b).astype(int)
    sj = np.floor((S[:, 1] - oy) / b).astype(int)
    nx, ny = si.max() + 2, sj.max() + 2
    ztop = np.full((nx, ny), -np.inf)
    np.maximum.at(ztop, (si, sj), S[:, 2])
    filled = 0
    for _ in range(hole_fill):
        P = np.pad(ztop, 1, constant_values=-np.inf)
        nb = np.stack([P[1 + a:nx + 1 + a, 1 + c:ny + 1 + c]
                       for a in (-1, 0, 1) for c in (-1, 0, 1) if (a, c) != (0, 0)])
        fill = ~np.isfinite(ztop) & (np.isfinite(nb).sum(0) >= 5)
        if not fill.any():
            break
        ztop[fill] = nb.max(0)[fill]
        filled += int(fill.sum())
    return ztop, si, sj, filled


def overhang_candidates(S, si, sj, gap, z_floor):
    k = S[:, 2] >= z_floor
    o = np.lexsort((S[k, 2], sj[k], si[k]))
    col = (si[k].astype(np.int64) * 10 ** 7 + sj[k])[o]
    z = S[k, 2][o]
    same = col[1:] == col[:-1]
    big = (np.diff(z) > gap) & same
    return np.unique(col[1:][big]), z[:-1][big], z[1:][big]


def pass2_rewrite(path, out, pos, table_path, ox, oy, b, reset):
    fx, fy, fz, fp = (pos[c] + 1 for c in ('1', '2', '3', '8'))
    cand = candidate_expr(fp, reset)
    # awk int() truncates toward zero -> explicit floor for particles outside the surface bbox
    prog = (
        f'function fl(v){{return (v<0 && v!=int(v)) ? int(v)-1 : int(v)}} '
        f'BEGIN{{while((getline l < "{table_path}")>0){{split(l,a,"\\t"); zt[a[1]","a[2]]=a[3]+0}}}} '
        f'NR==1{{print; next}} '
        f'{cand}{{nf++; i=fl((${fx}-ox)/b); j=fl((${fy}-oy)/b); k=i","j; '
        f'if((k in zt) && (${fz}+0 <= zt[k]+1e-6)){{${fp}=99; n++; if(nf%{SAMPLE_EVERY}==1) ns++}} '
        f'else {{${fp}=1}}}} '
        f'{{print}} END{{print n+0, ns+0 > "/dev/stderr"}}'
    )
    # plain floats only: repr(np.float64) is "np.float64(...)" on numpy>=2, which awk reads as 0
    ox, oy, b = float(ox), float(oy), float(b)
    with open(out, 'w') as fo:
        r = subprocess.run(['awk', '-F', '\t', '-v', 'OFS=\t', '-v', f'ox={ox:.10g}', '-v', f'oy={oy:.10g}',
                            '-v', f'b={b:.10g}', prog, path], stdout=fo, stderr=subprocess.PIPE, text=True,
                           check=True, env=dict(os.environ, LC_ALL='C'))
    n, ns = map(int, r.stderr.strip().splitlines()[-1].split())
    return n, ns


def write_vtk(out_txt, vtk_path, pos):
    fx, fy, fz, fp = (pos[c] + 1 for c in ('1', '2', '3', '8'))
    prog = f'NR>1 && ${fp}==99{{print ${fx}"\\t"${fy}"\\t"${fz}}}'
    r = subprocess.run(['awk', '-F', '\t', prog, out_txt], capture_output=True, text=True, check=True,
                       env=dict(os.environ, LC_ALL='C'))
    P = np.array(r.stdout.split(), dtype=np.float32).reshape(-1, 3)
    n = len(P)
    with open(vtk_path, 'wb') as f:
        f.write(b"# vtk DataFile Version 3.0\ninterior fluid (p_type=99)\nBINARY\nDATASET POLYDATA\n")
        f.write(f"POINTS {n} float\n".encode()); f.write(P.astype('>f4').tobytes())
        f.write(f"\nVERTICES {n} {2 * n}\n".encode())
        v = np.empty((n, 2), dtype='>i4'); v[:, 0] = 1; v[:, 1] = np.arange(n); f.write(v.tobytes())
        f.write(f"\nPOINT_DATA {n}\nSCALARS p_type int 1\nLOOKUP_TABLE default\n".encode())
        f.write(np.full(n, 99, dtype='>i4').tobytes()); f.write(b"\n")
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('input')
    ap.add_argument('--out')
    ap.add_argument('--bin', type=float)
    ap.add_argument('--floor-z', type=float, default=0.5)
    ap.add_argument('--gap', type=float)
    ap.add_argument('--hole-fill', type=int, default=0)
    ap.add_argument('--vtk', action='store_true')
    ap.add_argument('--reset-99', action='store_true',
                    help='re-evaluate rows already marked p_type=99 (old flags not interior by this rule go back to 1)')
    ap.add_argument('--dry-run', action='store_true')
    a = ap.parse_args()

    src = a.input
    out = a.out or os.path.splitext(src)[0] + '_int.txt'
    if os.path.abspath(out) == os.path.abspath(src):
        sys.exit("refusing to overwrite the input; give --out")
    tmp = out + '.surface.tmp'
    pos = read_header(src)
    t0 = time.time()

    print(f"[1/3] scanning {src} ...", flush=True)
    anchor, counts, sample = pass1_extract(src, pos, tmp, a.reset_99)
    if counts.get(99, 0) and not a.reset_99:
        print(f"      note: {counts[99]} rows already p_type=99 are kept as-is (use --reset-99 to re-evaluate them)")
    print("      p_type counts:", dict(sorted(counts.items())))
    if anchor is None or not any(k >= 1000 for k in counts):
        sys.exit("no fluid (p_type=1) or no surface (p_type>=1000) particles — nothing to do")
    x0, y0, z0, h0 = anchor
    dx = h0 / H_COEFF
    b = a.bin or dx
    S = np.loadtxt(tmp)
    if S.ndim == 1:
        S = S[None, :]
    print(f"      surface particles: {len(S)}  x[{S[:, 0].min():.2f},{S[:, 0].max():.2f}] "
          f"y[{S[:, 1].min():.2f},{S[:, 1].max():.2f}] z[{S[:, 2].min():.2f},{S[:, 2].max():.2f}]")
    nn = nn_spacing(S[:, :3])
    s_med = float(np.median(nn)) if len(nn) else float('nan')
    print(f"      fluid spacing dx={dx:.3f} (h={h0}), bin b={b:.3f}, surface NN spacing median {s_med:.3f} "
          f"(90% {np.percentile(nn, 90):.3f}, max {nn.max():.3f})")
    if s_med > b:
        print(f"      WARNING: surface spacing ({s_med:.2f}) > bin ({b:.2f}): roofs may have empty bins. "
              f"Use --bin >= {1.5 * s_med:.2f} or --hole-fill 1 (check corners afterwards).")

    # fluid-centred bins: bin edges at fluid_coord ± b/2, origin below all surface points
    ox = x0 - b / 2 - b * np.ceil((x0 - S[:, 0].min()) / b + 1)
    oy = y0 - b / 2 - b * np.ceil((y0 - S[:, 1].min()) / b + 1)

    print("[2/3] column-max heightfield ...", flush=True)
    ztop, si, sj, filled = build_ztop(S, ox, oy, b, a.hole_fill)
    nbin = int(np.isfinite(ztop).sum())
    print(f"      bins with surface: {nbin}   hole-filled bins: {filled}")
    gap = a.gap or 2 * b
    z_floor = S[:, 2].min() + a.floor_z
    cols, lo, hi = overhang_candidates(S, si, sj, gap, z_floor)
    print(f"      overhang candidates (vertical hole > {gap:.1f} m, points below z={z_floor:.2f} ignored): "
          f"{len(cols)} columns ({100 * len(cols) / max(nbin, 1):.1f}% of bins)")
    if len(cols):
        print(f"        -> treated as SOLID (conservative). hole lower-z median {np.median(lo):.1f} m, "
              f"hole size median {np.median(hi - lo):.1f} m. Inspect if this is large.")
    if a.dry_run:
        os.remove(tmp); print(f"dry run, nothing written ({time.time() - t0:.0f}s)"); return

    table = out + '.ztop.tmp'
    ii, jj = np.nonzero(np.isfinite(ztop))
    with open(table, 'w') as f:
        for i, j in zip(ii, jj):
            f.write(f"{i}\t{j}\t{ztop[i, j]:.6f}\n")

    # Python-side flags on the sampled fluid rows (same bins) — must equal awk's count on those rows
    fi = np.floor((sample[:, 0] - ox) / b).astype(int)
    fj = np.floor((sample[:, 1] - oy) / b).astype(int)
    inb = (fi >= 0) & (fi < ztop.shape[0]) & (fj >= 0) & (fj < ztop.shape[1])
    zs = np.full(len(sample), -np.inf); zs[inb] = ztop[fi[inb], fj[inb]]
    n_py = int((sample[:, 2] <= zs + 1e-6).sum())

    print(f"[3/3] writing {out} ...", flush=True)
    n99, n_awk = pass2_rewrite(src, out, pos, table, ox, oy, b, a.reset_99)
    nfl = counts.get(1, 0)
    print(f"      flagged p_type=99: {n99}  ({100 * n99 / max(nfl, 1):.2f}% of {nfl} fluid particles)")
    print(f"      self-check on {len(sample)} sampled fluid rows: awk {n_awk} vs python {n_py}", flush=True)
    if n_awk != n_py:
        os.remove(out)
        sys.exit("FAILED: awk and Python binning disagree — output removed, do not use")
    if a.vtk:
        vtk = os.path.splitext(out)[0] + '.interior.vtk'
        print(f"      vtk: {vtk} ({write_vtk(out, vtk, pos)} points)")
    for t in (tmp, table):
        os.remove(t)
    print(f"done ({time.time() - t0:.0f}s)")


if __name__ == '__main__':
    main()
