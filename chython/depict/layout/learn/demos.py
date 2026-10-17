# -*- coding: utf-8 -*-
#
#  Copyright 2026 Ramil Nugmanov <nougmanoff@protonmail.com>
#  This file is part of chython.
#
#  chython is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
"""N1 demonstrations: a beam search from the tree embedding toward a liked drawing, proving it reachable.

| part | rule |
| --- | --- |
| distance | the RMSD, in `L`, of the post-stepped plane to the liked one after the best turn, mirror and shift, \
the least over the first `AUTOMORPHISMS` automorphisms |
| search | a beam of `BEAM` over `Env(everywhere=True)` steps, up to `STEPS` and `BUDGET` seconds, nearest first; each \
component alone |
| best | the nearest state with every ring system drawn |
| output | per component `(distance, keys, component ranks)`; `python -m chython.depict.layout.learn.demos` writes \
`Trajectory`s, skipping sources the run already holds |
| primitive | `primitive_demonstration`: the same beam over `PrimitiveEnv` steps from the cached tree embedding, \
whole-drawing turns left out, one beam entry per drawing: first `_lattice` moves alone, which reach any honeycomb \
drawing, then every move from the best of those; then turns while they improve the score, since the distance ignores \
orientation. \
Its records are `mcts` self-play records with a one-hot visit share, so they train the policy as search targets |

    python -m chython.depict.layout.learn.demos --out RUN [--poses P] [--corpus C] [--workers N] [--limit N]
"""
from argparse import ArgumentParser
from math import atan2, cos, pi, sin, sqrt
from multiprocessing import Pool
from time import perf_counter

from ....core import smiles as read_smiles
from ..env import Env
from ..lattice import L
from ..primitive import SCENE
from ..quality import _automorphisms
from ..relax import fair, regularize, seat, wheels
from ..rings import ring_blocks
from .replay import Trajectory, key, read, write


__all__ = ['aligned', 'post', 'images', 'pending', 'demonstrate', 'primitive_demonstration', 'main']

BEAM = 8
STEPS = 40
BUDGET = 120.           # seconds of search per component; the best state so far stands
AUTOMORPHISMS = 8
LATTICE = (pi / 3, 2 * pi / 3, pi)    # the honeycomb's rotation targets


def aligned(a, b):
    """The RMSD, in `L`, of the point lists `a` and `b` after the best turn, mirror and shift of `b`."""
    n = len(a)
    ca = (sum(p[0] for p in a) / n, sum(p[1] for p in a) / n)
    cb = (sum(p[0] for p in b) / n, sum(p[1] for p in b) / n)
    pa = [(x - ca[0], y - ca[1]) for x, y in a]
    best = None
    for m in (1, -1):
        pb = [(x - cb[0], m * (y - cb[1])) for x, y in b]
        dot = sum(x * u + y * v for (x, y), (u, v) in zip(pa, pb))
        crs = sum(y * u - x * v for (x, y), (u, v) in zip(pa, pb))
        t = atan2(crs, dot)
        c, s = cos(t), sin(t)
        r = sum((x - c * u + s * v) ** 2 + (y - s * u - c * v) ** 2 for (x, y), (u, v) in zip(pa, pb))
        best = r if best is None else min(best, r)
    return sqrt(max(best, 0.) / n) / L


def post(env):
    """`env`'s plane after the ring post-steps, a copy."""
    plane = dict(env.plane)
    rings = [r for _, rs in env.blocks for r in rs]
    rank = env.scope.rank
    regularize(env.mol, plane, rings, rank)
    wheels(env.mol, plane, rings, rank)
    fair(env.mol, plane, rings, rank)
    seat(env.mol, plane, rings, rank)
    return plane


def images(mol):
    """The first `AUTOMORPHISMS` automorphisms of `mol`, or the identity alone where the query refuses: an isotope
    atom is one such."""
    try:
        return list(_automorphisms(mol, AUTOMORPHISMS)) or [{x: x for x in mol}]
    except ValueError:
        return [{x: x for x in mol}]


def pending(run_dir, jobs):
    """The `jobs`, `(source, smiles, plane)`, whose source `run_dir` holds no trajectory of yet."""
    done = {t.source for t in read(run_dir)}
    return [j for j in jobs if j[0] not in done]


def demonstrate(mol, target, beam=BEAM, steps=STEPS, budget=BUDGET):
    """`[(distance, keys, component ranks)]`, per component of `mol` in rank order, of the nearest drawing the search
    reaches to `target`, `{atom: (x, y)}`."""
    rank = mol.canonical_order()
    maps = images(mol)
    out = []
    for c in sorted((set(c) for c in mol.connected_components), key=lambda c: min(rank[x] for x in c)):
        atoms = sorted(c, key=rank.__getitem__)
        goals = [[target[m[x]] for x in atoms] for m in maps]

        def dist(env):
            p = post(env)
            mine = [p[x] for x in atoms]
            return min(aligned(mine, g) for g in goals)
        env = Env(mol, c, rank, [b for b in ring_blocks(mol) if b[0] <= c], everywhere=True)
        states = [(dist(env), [], env)]
        best = None if -1 in env.drawn else states[0]
        deadline = perf_counter() + budget
        for _ in range(steps):
            if perf_counter() > deadline:
                break
            succ = []
            for d, keys, e in states:
                for s in e.legal():
                    if s.move.kind == 'stop':
                        continue
                    f = e.fork()
                    f.step(s)
                    succ.append((dist(f), keys + [key(e, s)], f))
            if not succ:
                break
            succ.sort(key=lambda t: (t[0], len(t[1]), t[1]))
            states = succ[:beam]
            drawn = [t for t in states if -1 not in t[2].drawn]
            if drawn and (best is None or drawn[0][0] < best[0] - 1e-9):
                best = drawn[0]
            elif best is not None and not any(-1 in t[2].drawn for t in states):
                break
        if best is None:
            best = states[0]
        out.append((best[0], best[1], tuple(sorted(rank[x] for x in c))))
    return out


def _lattice(move):
    """A move that keeps a drawing on the honeycomb: a flip, or a rotation to 60, 120 or 180 degrees."""
    return move.kind == 'flip' or move.kind == 'rotate' and any(abs(abs(move.arg) - a) < 1e-9 for a in LATTICE)


def _signature(env, atoms):
    return tuple(round(c, 3) for x in atoms for c in env.plane[x])


def primitive_demonstration(initial, goals, source=None, beam=BEAM, steps=60, budget=BUDGET):
    """`(distance, records)` of the nearest drawing a beam over `initial`'s primitive steps reaches to any of `goals`,
    point lists in `initial.scope.atoms` order; one record per state of the route, its `stop` included."""
    from .mcts import step_key
    atoms = initial.scope.atoms

    def dist(env):
        mine = [env.plane[x] for x in atoms]
        return min(aligned(mine, g) for g in goals)
    deadline = perf_counter() + budget

    def search(first, allowed, steps):
        states = [first]
        seen = {_signature(first[2], atoms)}
        best = first
        for _ in range(steps):
            if perf_counter() > deadline:
                break
            succ = []
            for d, route, e in states:
                for s in e.legal():
                    if s.move.kind in ('stop', 'turn') or not allowed(s.move):
                        continue
                    f = e.fork()
                    f.step(s)
                    succ.append((dist(f), route + ((e, s),), f))
            succ.sort(key=lambda t: (t[0], len(t[1]), [step_key(e, s) for e, s in t[1]]))
            states = []
            for t in succ:                                  # one entry per drawing: routes to it are equal
                sign = _signature(t[2], atoms)
                if sign not in seen:
                    seen.add(sign)
                    states.append(t)
                    if len(states) == beam:
                        break
            if not states:
                break
            if states[0][0] < best[0] - 1e-9:
                best = states[0]
        return best
    start = initial.env()
    lattice = search((dist(start), (), start), _lattice, steps)
    best = search(lattice, lambda m: True, steps - len(lattice[1]))
    distance, route, env = best
    route = list(route)
    for _ in range(len(SCENE)):
        turns = [s for s in env.legal() if s.move.kind == 'turn' and s.rank < env.rank]
        if not turns:
            break
        s = min(turns, key=lambda s: (s.rank, env.tie(s.move)))
        f = env.fork()
        f.step(s)
        route.append((env, s))
        env = f
    if env.legal() and env.legal()[-1].move.kind == 'stop':     # a drawing with a defect is not stopped
        route.append((env, env.legal()[-1]))
    records = [{'source': source, 'plane': [list(e.plane[x]) for x in atoms], 'steps': e.steps,
                'visits': {repr(step_key(e, s)): 1.}} for e, s in route]
    return distance, records


def _work(job):
    source, text, plane = job
    t = perf_counter()
    try:
        mol = read_smiles(text)
        found = demonstrate(mol, dict(zip(mol, map(tuple, plane))))
    except Exception as e:                                  # one record's failure is a result, not a crash
        return source, text, f'{type(e).__name__}: {e}', perf_counter() - t
    return source, text, found, perf_counter() - t


def main(argv=None):
    from ..scorecard import CORPUS, POSES, _read_poses, read_corpus
    ap = ArgumentParser(prog='python -m chython.depict.layout.learn.demos')
    ap.add_argument('--out', required=True)
    ap.add_argument('--poses', default=str(POSES))
    ap.add_argument('--corpus', default=str(CORPUS))
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--limit', type=int)
    a = ap.parse_args(argv)
    rows = {r['id']: r for r in read_corpus(a.corpus)}
    poses = _read_poses(a.poses)
    jobs = [(i, rows[i]['smiles'], p) for i, p in sorted(poses.items())
            if i in rows and '>' not in rows[i]['smiles']]
    if a.limit:
        jobs = jobs[:a.limit]
    jobs = pending(a.out, jobs)
    reached = total = 0
    with Pool(a.workers, maxtasksperchild=10) as pool:
        for source, text, found, took in pool.imap_unordered(_work, jobs):
            if isinstance(found, str):
                print(f'{source}\terror\t{found}', flush=True)
                continue
            write(a.out, [Trajectory(source, text, keys, {'everywhere': True, 'component': list(comp), 'distance': d,
                                                          'seconds': round(took, 2)}) for d, keys, comp in found])
            total += len(found)
            reached += sum(d < .15 for d, _, _ in found)
            print(f'{source}\t{len(found)}\t' + ','.join(f'{d:.3f}' for d, _, _ in found) + f'\t{took:.1f}', flush=True)
    print(f'reached within 0.15 L: {reached} of {total} components')


if __name__ == '__main__':
    main()
