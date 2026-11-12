"""
Forces while stacking the slices (one printed plate per z) of a design: the axial force on the plate being added,
and the axial force on every individual block (can a block be pulled out of its pocket?), for several orders.

Blocks within one plate share a z centre, so by symmetry they put no axial force on each other: filling a plate
only loads the pocket walls in-plane. Axial (pull-out) loads come only from other plates. So every state is a sum
of slice-pair terms: f[i, j, s] = axial force on each block of slice i from all of slice j when j is moved s mm
along z from its final place relative to i. Those are computed once (pole-face surface-charge integral, as in
forces_12cm.py) and summed for each assembly state.

Signs: "push" > 0 = the incoming plate is pushed away from the stack (you must press it in); < 0 = pulled in.
A block "pulled out" means the force on it points away from its own plate's mid-plane toward a neighbouring
plate, i.e. it wants to leave its pocket through that face; that is what an open-faced pocket must resist.

    python utilities/assembly_12cm.py designs/12cm_bore200/placements.csv
"""
import argparse
import os
import pickle
from concurrent.futures import ProcessPoolExecutor

import magpylib as magpy
import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation as R

MU0 = 4e-7 * np.pi
DIM = np.array([12.7, 12.7, 25.4])
GAPS = np.array([0, 1, 2, 3, 5, 10, 20, 50])          # mm beyond final position (stepwise orders)
SEP = np.array([0, 1, 2, 3, 5, 7, 10, 15, 20, 30, 40, 50, 60, 80, 100, 130, 160, 200, 250, 300, 400, 600, 1000.])
ROD_EXTRA = 2 * (10 + 20)   # mm beyond the plate stack: 10 mm end plate + 20 mm washer/nut/spare thread, each end
BR = 1300.0
NQ = 5
_G = {}


def face_points():
    gy, wy = np.polynomial.legendre.leggauss(NQ)
    gz, wz = np.polynomial.legendre.leggauss(2 * NQ)
    Y, Z = np.meshgrid(gy * DIM[1] / 2, gz * DIM[2] / 2, indexing='ij')
    W = np.outer(wy * DIM[1] / 2, wz * DIM[2] / 2).ravel()
    loc = np.vstack([np.c_[np.full(Y.size, s * (DIM[0] / 2 + 1e-3)), Y.ravel(), Z.ravel()] for s in (1, -1)])
    return loc, np.r_[W, -W]                                    # signed weights, mm^2


def _init(sl):
    _G['sl'] = sl
    _G['loc'], _G['sw'] = face_points()


def pair(args):
    """Axial force (N) on each block of slice i from slice j shifted by s mm in z."""
    i, j, s = args
    (pm, qm), (ps, qs) = _G['sl'][i], _G['sl'][j]
    loc, sw = _G['loc'], _G['sw']
    pts = R.from_quat(np.repeat(qm, len(loc), 0)).apply(np.tile(loc, (len(pm), 1))) + np.repeat(pm, len(loc), 0)
    m = len(pts)
    ps = ps + [0, 0, s]
    Bz = np.zeros(m)
    step = max(1, 700_000 // m)
    for a in range(0, len(ps), step):
        k = min(step, len(ps) - a)
        B = magpy.getB('Cuboid', observers=np.tile(pts, (k, 1)), position=np.repeat(ps[a:a + k], m, 0),
                       orientation=R.from_quat(np.repeat(qs[a:a + k], m, 0)), magnetization=(BR, 0, 0),
                       dimension=DIM)
        Bz += B[:, 2].reshape(k, m).sum(0)
    f = BR * 1e-3 / MU0 * (np.tile(sw, len(pm)) * 1e-6) * (Bz * 1e-3)
    return (i, j, s), f.reshape(len(pm), -1).sum(1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('csv')
    ap.add_argument('--jobs', type=int, default=5)
    a = ap.parse_args()
    df = pd.read_csv(a.csv)
    df['zk'] = df.z.round(1)
    zs = np.sort(df.zk.unique())
    nz = len(zs)
    sl = []
    for z in zs:
        g = df[df.zk == z]
        sl.append((g[['x', 'y', 'z']].to_numpy(),
                   R.from_euler('z', g.rot_z_deg.to_numpy()[:, None], degrees=True).as_quat()))

    # separation grid: plate j moved s mm further from plate i than in the finished magnet (never closer)
    tasks = [(i, j, float(np.sign(j - i) * g)) for i in range(nz) for j in range(nz) if i != j for g in SEP]
    cache = a.csv.rsplit('/', 1)[0] + '/assembly_pairs.pkl'     # slice-pair forces, reused on reruns
    if os.path.exists(cache) and os.path.getmtime(cache) > os.path.getmtime(a.csv):
        F = pickle.load(open(cache, 'rb'))
    else:
        with ProcessPoolExecutor(a.jobs, initializer=_init, initargs=(sl,)) as ex:
            F = dict(ex.map(pair, tasks, chunksize=8))
        pickle.dump(F, open(cache, 'wb'))
    T = {(i, j): np.array([F[(i, j, float(np.sign(j - i) * g))] for g in SEP])
         for i in range(nz) for j in range(nz) if i != j}           # (len(SEP), blocks of i)

    def fpair(i, j, sep):
        """Per-block Fz on plate i from plate j, j separated `sep` >= 0 mm beyond final (linear interp)."""
        sep = min(max(sep, 0.0), SEP[-1])
        k = min(np.searchsorted(SEP, sep, side='right') - 1, len(SEP) - 2)
        t = (sep - SEP[k]) / (SEP[k + 1] - SEP[k])
        return (1 - t) * T[(i, j)][k] + t * T[(i, j)][k + 1]

    def blocks_at(x):
        """Per-block Fz for every plate in dict x (plate -> extra z displacement from final, mm)."""
        return {i: sum(fpair(i, j, (x[j] - x[i]) * np.sign(j - i)) for j in x if j != i) for i in x}

    def state(placed, moving, d, g):
        """Per-block Fz for every slice in placed+moving; moving slices sit d*g mm from final."""
        return blocks_at({k: 0.0 for k in placed} | {k: float(d * g) for k in moving})

    # assembly orders: list of steps (placed slices, moving slices, direction the moving set comes from)
    mid = nz // 2
    orders = {
        'A  one end to the other': [(list(range(k)), [k], +1) for k in range(1, nz)],
        'B  each half from the midplane outward, then join': (
            [(list(range(mid, k)), [k], +1) for k in range(mid + 1, nz)]
            + [(list(range(k + 1, mid)), [k], -1) for k in range(mid - 2, -1, -1)]
            + [(list(range(mid)), list(range(mid, nz)), +1)]),
        'C  each half from its end inward, then join': (
            [(list(range(k + 1, nz)), [k], -1) for k in range(nz - 2, mid - 1, -1)]
            + [(list(range(k)), [k], +1) for k in range(1, mid)]
            + [(list(range(mid)), list(range(mid, nz)), +1)]),
    }

    kg = 1 / 9.81
    print(f'{len(df)} blocks, {nz} plates at z = {", ".join(f"{z:g}" for z in zs)} mm; Br {BR:g} mT')
    print(f'extra gap beyond final position (mm): {", ".join(f"{g:g}" for g in GAPS)}\n')
    worst_block = (0, None)
    for name, steps in orders.items():
        print(f'== order {name}')
        print(f'   {"step":42s} {"plate push N vs gap":>58s}   peak push    worst block pull-out N (gap, plate z)')
        for placed, moving, d in steps:
            push, blk = [], []
            for g in GAPS:
                st = state(placed, moving, d, g)
                push.append(d * sum(st[k].sum() for k in moving))
                # worst block: force away from its own plate toward a neighbour = |Fz| (pockets open on both faces)
                b = max(((np.abs(st[k]).max(), k) for k in st))
                blk.append(b)
            push = np.array(push)
            if len(moving) == 1:
                label = f'add z={zs[moving[0]]:+.1f} onto {len(placed)} plate(s)'
            else:
                label = f'join halves ({len(placed)} + {len(moving)} plates)'
            ip = int(np.argmax(np.abs(push)))
            ib = int(np.argmax([b[0] for b in blk]))
            print(f'   {label:42s} ' + ''.join(f'{v:7.0f}' for v in push)
                  + f'   {push[ip]:+6.0f} ({push[ip] * kg:+5.1f} kgf)'
                  f'   {blk[ib][0]:5.1f} (gap {GAPS[ib]:g}, z={zs[blk[ib][1]]:+.1f})')
            if blk[ib][0] > worst_block[0]:
                worst_block = (blk[ib][0], f'{name.split()[0]}: {label}, gap {GAPS[ib]:g} mm, block in plate '
                                           f'z={zs[blk[ib][1]]:+.1f}')
        print()

    final = state(list(range(nz)), [], 1, 0)
    fin = max(np.abs(final[k]).max() for k in final)
    print(f'worst single-block axial force anywhere in assembly: {worst_block[0]:.1f} N '
          f'({worst_block[0] * kg:.2f} kgf) -- {worst_block[1]}')
    print(f'worst single-block axial force in the finished magnet: {fin:.1f} N  (block weight 0.30 N)')

    # signed per-block axial force through order A (start at the -z end, each new plate comes in from +z):
    # "back" = toward -z = toward the plates already built; "fwd" = toward +z = toward the open end / next plate
    for title, steps in [('A, built from the -z end (new plates arrive from +z)', orders['A  one end to the other'])]:
        lo = {i: np.zeros(len(sl[i][0])) for i in range(nz)}
        hi = {i: np.zeros(len(sl[i][0])) for i in range(nz)}
        for placed, moving, d in [([0], [], 1)] + steps:
            for g in GAPS:
                for i, f in state(placed, moving, d, g).items():
                    lo[i], hi[i] = np.minimum(lo[i], f), np.maximum(hi[i], f)
        print(f'\nsigned per-block axial force, order {title}, over every step and gap (N):')
        print('  plate z  blocks | max back (-z)  blocks ever back | max fwd (+z)  blocks ever fwd | '
              'finished magnet: back / fwd, blocks back')
        tot_b = tot_f = 0
        for i in range(nz):
            fb, ff = final[i].min(), final[i].max()
            nb, nf = int((lo[i] < -0.5).sum()), int((hi[i] > 0.5).sum())
            tot_b, tot_f = tot_b + nb, tot_f + nf
            print(f'  {zs[i]:+7.1f} {len(lo[i]):6d} | {-lo[i].min():12.1f}  {nb:16d} | {hi[i].max():11.1f}  {nf:15d} | '
                  f'{-min(fb, 0):5.1f} / {max(ff, 0):5.1f}, {int((final[i] < 0).sum()):3d} of {len(final[i])}')
        print(f'  blocks that at some point are pushed back (>0.5 N toward -z): {tot_b} of {len(df)}; '
              f'forward (>0.5 N toward +z): {tot_f} of {len(df)}')
        print('  (building from the +z end instead is the mirror image: swap back and fwd)')

    # ---- all plates threaded onto long rods, spread out, then the end nuts are run in together ----
    # Plates slide freely (no friction). Unknowns: the 11 gaps g_i >= 0 beyond final (sum = E) and the end-nut
    # force N. Plate k's contact force from the plate below is C_k = N + sum_{j<k} Fz_j (C_0 = N, the -z nut);
    # each interface is either open (g > 0, C = 0) or touching (g = 0, C >= 0): Fischer-Burmeister residuals.
    # The +z nut carries the same N because the magnetic forces on all plates sum to zero.
    from scipy.optimize import least_squares
    span = df.z.max() - df.z.min() + DIM[2]
    SC = 10.0                                                      # N per mm, residual scaling

    def contacts(g, N):
        x = dict(zip(range(nz), np.r_[0, np.cumsum(g)] - g.sum() / 2))
        st = blocks_at(x)
        Fz = np.array([st[k].sum() for k in range(nz)])
        return st, N + np.cumsum(Fz)[:-1]                          # C_1..C_{nz-1}

    def fb(a_, b_):
        return a_ + b_ - np.sqrt(a_ ** 2 + b_ ** 2)

    print(f'\n== all plates on long rods, closed together by the end nuts (frictionless sliding plates)')
    print(f'closed magnet span {span:.1f} mm; rod length = E + span + {ROD_EXTRA} mm '
          f'(10 mm end plates, 20 mm washer/nut/spare thread per end)')
    print('   E mm  rod mm | end-nut force N (kgf)  per rod (6) | gaps beyond final min/max, # touching | '
          'worst block N: toward mid / toward end')
    v = None
    rows = []
    for E in [1000, 800, 600, 500, 400, 300, 250, 200, 150, 100, 70, 50, 30, 20, 10, 5, 2, 0]:
        if E > 0:
            if v is None:
                v = np.r_[np.full(nz - 1, E / (nz - 1)), 0.0]
            else:
                v = np.r_[v[:-1] * E / max(v[:-1].sum(), 1e-9), v[-1]]

            def resid(v_):
                g, N = v_[:-1], v_[-1]
                _, C = contacts(g, N)
                return np.r_[fb(g, C / SC), g.sum() - E]

            v = least_squares(resid, v, bounds=(np.r_[np.zeros(nz - 1), -np.inf], np.inf),
                              xtol=1e-9, ftol=1e-10, max_nfev=400).x
            g, N = v[:-1], v[-1]
            st, C = contacts(g, N)
            err = np.abs(resid(v)).max() * SC
        else:
            g = np.zeros(nz - 1)
            st, C0 = contacts(g, 0.0)
            N = max(0.0, -C0.min())                                # least clamp that keeps every interface closed
            C = C0 + N
            err = 0.0
        tm = max(max(-st[k].min() if zs[k] > 0 else st[k].max(), 0) for k in range(nz))
        te = max(max(st[k].max() if zs[k] > 0 else -st[k].min(), 0) for k in range(nz))
        rows.append((E, N))
        print(f'  {E:5.0f} {E + span + ROD_EXTRA:6.0f} | {N:8.0f} ({N / 9.81:6.1f})  {N / 6:8.0f} N   | '
              f'{g.min():8.1f} / {g.max():6.1f}, {int((g < 0.05).sum()):2d}          | {tm:6.1f} / {te:6.1f}'
              + (f'   (residual {err:.1f} N)' if err > 2 else ''))
    for lim in (20, 50, 100):
        ok = [E for E, f in rows if f <= lim]
        if ok:
            E = min(ok)
            print(f'  end force <= {lim:3d} N ({lim / 9.81:4.1f} kgf) from E = {E:.0f} mm: rod >= {E + span + ROD_EXTRA:.0f} mm')


if __name__ == '__main__':
    main()
