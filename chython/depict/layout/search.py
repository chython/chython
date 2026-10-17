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
"""The controller's search: a beam over edits of the tree embedding, ranked by the score alone.

| step | rule |
| --- | --- |
| embed | each ring system, in rank order, takes one of its `variants`; a masked variant is dropped, so is one leaving \
more cis/trans broken than another, and a system with none left keeps its tree drawing |
| repair | the frontier's flips, rotations and stretches, and each other variant of a drawn system on the frontier, \
once per system; bends only for a state no other move betters |
| beam | the `BEAM` best distinct planes, by broken cis/trans units, then score, ties by `move_key` |
| stop | no successor betters the best state, or `2 × systems + STEPS` moves |
| post | `regularize`, `wheels`, `fair`, `seat`, then `pose`: `level`, `horizontal`, `straight` and `upright` |
| log | `depict:budget` when the move budget ends the search, `depict:open-ring` for a closure left off `L` |

    python -m chython.depict.layout.search stats [--corpus P] [--out F] [--limit N] [--greedy]

writes per molecule its atoms, ring systems, moves, candidates per step, budget, seconds and defects, and prints the
medians per stratum.
"""
from argparse import ArgumentParser
from statistics import median
from time import perf_counter
from typing import NamedTuple

from ...core import AutomorphismBudgetExceeded, LogRecord, ReactionContainer, smiles
from .edits import Move, apply, broken, candidates, embed, masked, move_key, variants
from .env import Env
from .greedy import best_step, drive
from .pose import GAP, horizontal, level, straight, upright
from .quality import defects, rescaled
from .relax import _ct_broken, fair, regularize, seat, wheels
from .rings import ring_blocks
from .score import Scope, add, delta, frontier, full, key, open_closures
from .tree import tree_embedding


__all__ = ['State', 'search', 'controller_layout', 'peptide_search']

BEAM = 4
STEPS = 16


class State(NamedTuple):
    plane: dict
    score: tuple
    ct: int                 # stated cis/trans units drawn against their parity
    drawn: tuple            # per ring system: its variant index, -1 while undrawn
    redrawn: frozenset      # ring systems drawn a second time
    moves: tuple            # the moves taken, in order


def _hash(scope, plane):
    return tuple(round(c, 4) for x in scope.atoms for c in plane[x])


def _successors(scope, state, moves, blocks, drawings, root):
    """The unmasked successors of `state`; of one ring system's redraws only those leaving fewest cis/trans broken."""
    out = []
    least = {}                                              # ring system -> fewest broken units among its redraws
    for m in moves:
        if m.kind == 'embed':
            moved = embed(scope, state.plane, blocks[m.pivot][0], drawings[m.pivot][m.other], root)
        else:
            moved = apply(scope, state.plane, m)
        if masked(scope, state.plane, moved):
            continue
        ct = state.ct - _ct_broken(scope.mol, state.plane, scope.frames, set(moved)) + broken(scope, state.plane, moved)
        if m.kind == 'embed':
            least[m.pivot] = min(least.get(m.pivot, ct), ct)
        out.append((add(state.score, delta(scope, state.plane, moved)), m, state, moved, ct))
    return [t for t in out if t[1].kind != 'embed' or t[4] == least[t[1].pivot]]


def _rank(ct, score):
    """States compare by broken cis/trans units first, then by score."""
    return ct, key(score)


def _keep(scope, succ, beam):
    """The `BEAM` best distinct successors as states."""
    order = {id(s): k for k, s in enumerate(beam)}
    succ.sort(key=lambda t: (_rank(t[4], t[0]), order[id(t[2])], move_key(scope, t[1])))
    out, seen = [], set()
    for score, m, s, moved, ct in succ:
        plane = {**s.plane, **moved}
        h = _hash(scope, plane)
        if h in seen:
            continue
        seen.add(h)
        drawn, redrawn = s.drawn, s.redrawn
        if m.kind == 'embed':
            if drawn[m.pivot] >= 0:
                redrawn = redrawn | {m.pivot}
            drawn = drawn[:m.pivot] + (m.other,) + drawn[m.pivot + 1:]
        out.append(State(plane, score, ct, drawn, redrawn, s.moves + (m,)))
        if len(out) == BEAM:
            break
    return out


def search(scope, blocks, root, plane, stats=None):
    """The best `State` reached from `plane`; `stats` collects `moves`, `candidates` per step and `budget`."""
    drawings = [variants(scope, b) for b in blocks]
    ct = _ct_broken(scope.mol, plane, scope.frames, set(plane))
    beam = [State(dict(plane), full(scope, plane), ct, (-1,) * len(blocks), frozenset(), ())]
    cap = 2 * len(blocks) + STEPS
    counts = []
    for i, block in enumerate(blocks):
        moves = [Move('embed', i, j, 0., frozenset()) for j in range(len(drawings[i]))]
        succ = [t for s in beam for t in _successors(scope, s, moves, blocks, drawings, root)]
        if succ:
            counts.append(len(succ))
            beam = _keep(scope, succ, beam)
    best = min(beam, key=lambda s: _rank(s.ct, s.score))
    budget = False
    while True:
        if len(best.moves) >= cap:
            budget = True
            break
        succ = []
        for s in beam:
            front = frontier(scope, s.plane)
            if not front:
                continue
            moves = candidates(scope, s.plane, front)
            for i, block in enumerate(blocks):
                if s.drawn[i] >= 0 and i not in s.redrawn and block[0] & front:
                    moves += [Move('embed', i, j, 0., frozenset()) for j in range(len(drawings[i])) if j != s.drawn[i]]
            gen = _successors(scope, s, moves, blocks, drawings, root)
            if not any(_rank(t[4], t[0]) < _rank(s.ct, s.score) for t in gen):
                gen += _successors(scope, s, candidates(scope, s.plane, front, violate=True), blocks, drawings, root)
            succ += gen
        if not succ:
            break
        counts.append(len(succ))
        top = _keep(scope, succ, beam)
        if _rank(top[0].ct, top[0].score) >= _rank(best.ct, best.score):
            break
        beam, best = top, top[0]
    if stats is not None:
        stats.update(atoms=len(scope.atoms), systems=len(blocks), moves=len(best.moves), candidates=counts,
                     budget=budget, undrawn=sum(not d for d in drawings))
    return best


def controller_layout(mol, log=None, stats=None, choose=None):
    """`{atom: (x, y)}` for every atom of `mol`, mean bond near `L`; `stats`, a list, gets one dict per component.
    With `choose`, each component is driven greedily through an `Env` by `choose(env, steps)` instead of the beam; a
    chooser whose `everywhere` attribute is true is offered edits on every atom, and one whose `primitive` attribute is
    true drives a `primitive.PrimitiveEnv` from the tree embedding, with no post-step."""
    if not len(mol):
        return {}
    rank = mol.canonical_order()
    blocks = ring_blocks(mol)
    comps = sorted((set(c) for c in mol.connected_components), key=lambda c: min(rank[x] for x in c))
    try:
        orbit = mol.automorphism_orbits()
    except AutomorphismBudgetExceeded:
        orbit = None
    plane = {}
    x0 = 0.
    for c in comps:
        own = [b for b in blocks if b[0] <= c]
        st = {}
        if getattr(choose, 'primitive', False):         # the tree start and primitive steps, no post-step
            from .primitive import Initial
            env = drive(Initial(mol, c, rank).env(), choose)
            st.update(atoms=len(c), systems=len(own), moves=env.steps, budget=env.steps >= env.cap, candidates=[],
                      undrawn=0)
            pts = dict(env.plane)
            xs = [p[0] for p in pts.values()]
            ys = [p[1] for p in pts.values()]
            dx = x0 - min(xs)
            dy = -(min(ys) + max(ys)) / 2
            for a, (x, y) in pts.items():
                plane[a] = (x + dx, y + dy)
            x0 = max(xs) + dx + GAP
            if stats is not None:
                stats.append(st)
            continue
        if choose is None:
            tree = tree_embedding(Scope(mol, c, rank, []))
            scope = Scope(mol, c, rank, tree.closures)
            pts = dict(search(scope, own, tree.root, tree.plane, st).plane)
        else:
            env = drive(Env(mol, c, rank, own, everywhere=getattr(choose, 'everywhere', False)), choose)
            scope, pts = env.scope, dict(env.plane)
            st.update(atoms=len(c), systems=len(own), moves=env.steps, budget=env.steps >= env.cap,
                      candidates=[], undrawn=sum(d == -2 for d in env.drawn))
        rings = [r for _, rs in own for r in rs]
        regularize(mol, pts, rings, rank)
        wheels(mol, pts, rings, rank)
        fair(mol, pts, rings, rank)
        seat(mol, pts, rings, rank)
        level(pts)
        horizontal(pts, rings)
        if not own:
            straight(mol, pts, rank)
        upright(mol, c, pts, orbit)
        order = tuple(sorted(c, key=rank.__getitem__))
        if log is not None:
            if st['budget']:
                log.append(LogRecord('depict:budget', order, 'the controller stopped at its move budget'))
            if gone := open_closures(scope, pts):
                atoms = tuple(sorted({x for b in gone for x in b}, key=rank.__getitem__))
                log.append(LogRecord('depict:open-ring', atoms, 'a ring closure bond is drawn off L'))
        if stats is not None:
            stats.append(st)
        xs = [p[0] for p in pts.values()]
        ys = [p[1] for p in pts.values()]
        dx = x0 - min(xs)
        dy = -(min(ys) + max(ys)) / 2
        for a, (x, y) in pts.items():
            plane[a] = (x + dx, y + dy)
        x0 = max(xs) + dx + GAP
    return plane


def _stats(corpus, out, limit=None, greedy=False):
    from .scorecard import CORPUS, read_corpus
    rows = read_corpus(corpus or CORPUS)
    if limit:
        rows = rows[:limit]
    per = {}
    with open(out, 'w') as f:
        f.write('id\tstratum\tatoms\tsystems\tmoves\tcand_med\tcand_max\tbudget\tseconds\tdefects\topen\n')
        for row in rows:
            mol = smiles(row['smiles'])
            if isinstance(mol, ReactionContainer):
                continue
            st, log = [], []
            t = perf_counter()
            try:
                plane = controller_layout(mol, log, st, choose=best_step if greedy else None)
            except Exception as e:                          # a failure is a result, not a crash
                f.write(f"{row['id']}\t{row['stratum']}\terror\t{type(e).__name__}: {e}\n")
                continue
            sec = perf_counter() - t
            cands = [n for s in st for n in s['candidates']] or [0]
            bad = sum(v for k, v in defects(mol, rescaled(mol, plane)).items() if k != 'ovl61')
            moves = sum(s['moves'] for s in st)
            rec = (len(mol), sum(s['systems'] for s in st), moves, median(cands), max(cands),
                   any(s['budget'] for s in st), round(sec, 3), bad,
                   sum(r.rule == 'depict:open-ring' for r in log))
            f.write('\t'.join(map(str, (row['id'], row['stratum'], *rec))) + '\n')
            f.flush()
            per.setdefault(row['stratum'], []).append(rec)
    print('| stratum | rows | moves med | cand med | cand max med | budget % | s med | defect rows % | open % |')
    print('| --- | --- | --- | --- | --- | --- | --- | --- | --- |')
    for k, rs in sorted(per.items()):
        n = len(rs)
        print(f'| {k} | {n} | {median(r[2] for r in rs)} | {median(r[3] for r in rs)} | {median(r[4] for r in rs)} | '
              f'{100 * sum(r[5] for r in rs) / n:.1f} | {median(r[6] for r in rs)} | '
              f'{100 * sum(r[7] > 0 for r in rs) / n:.1f} | {100 * sum(r[8] > 0 for r in rs) / n:.1f} |')


def main(argv=None):
    ap = ArgumentParser(prog='python -m chython.depict.layout.search')
    sub = ap.add_subparsers(dest='cmd', required=True)
    s = sub.add_parser('stats')
    s.add_argument('--corpus')
    s.add_argument('--out', default='controller_stats.tsv')
    s.add_argument('--limit', type=int)
    s.add_argument('--greedy', action='store_true')
    a = ap.parse_args(argv)
    _stats(a.corpus, a.out, a.limit, a.greedy)


if __name__ == '__main__':
    main()


def _peptide_score(ctx, inv):
    """`(plane, rank key)` of the backbone with the inverted vertices `inv`, or `(None, None)` when it fits no slot:
    collisions and broken backbone rings, then lost links (teleports), then stretch, then height in half-bonds, then
    inverted vertices and vertices off their ring's side."""
    from .peptide import L as PL, draw, frame
    plane, (bad, lost, stretch) = draw(ctx, inv)
    if plane is None:
        return None, None
    differ, same = frame(ctx)[:2]
    ys = [plane[x][1] for x in ctx.peptide.backbone]
    height = round((max(ys) - min(ys)) / (PL / 2))
    broken = sum((a in inv) == (b in inv) for a, b in differ)
    off = sum((a in inv) != (b in inv) for a, b in same)
    return plane, (bad + broken, lost, round(stretch, 9), height, len(inv) + off)


def _pin_moves(ctx, pins, relax, inv, near):
    """The search's moves from `(pins, relax)`: a step or a kept vertex near each stretched anchor in `near`, each
    cross-link's two hosts pinned to one side, then the other, each `same` pair the solve holds let go, and all of
    them at once."""
    from .peptide import frame, hosts, side
    n = len(ctx.peptide.backbone)
    out = []
    for v in sorted(near):
        for j in sorted(range(max(1, v - 3), min(n - 2, v + 3) + 1), key=lambda j: (abs(j - v), j)):
            if j + 1 <= n - 2:
                out.append(({**pins, j: True, j + 1: True}, relax))
            out.append(({**pins, j: False}, relax))
            if j == n - 2:
                out.append(({**pins, j: True}, relax))
    for _, _, ka, kb in hosts(ctx):
        for s in (1, -1):
            new = dict(pins)
            for k in (ka, kb):
                if 1 <= k <= n - 2:
                    new[k] = side(k, False) != s
            out.append((new, relax))
    held = frame(ctx)[1] - set(relax)
    for p in sorted(held):
        out.append((pins, relax | {p}))
    if len(held) > 1:
        out.append((pins, relax | held))
    return out


def peptide_search(ctx):
    """The inverted backbone vertices of a peptide: a beam of `BEAM` over pins and let-go `same` pairs, each state
    completed by `peptide.solve`, from no pin; ranked by `_peptide_score`, ties by the sorted vertices; None when no
    drawing fits."""
    from .peptide import draw, solve
    def state(pins, relax):
        inv = solve(ctx, pins, relax)
        plane, score = _peptide_score(ctx, inv)
        return None if plane is None else (score, sorted(inv), inv, pins, relax)
    best = state({}, frozenset())
    if best is None:
        return None
    beam, seen = [best], {best[2]}
    for _ in range(len(ctx.peptide.backbone) // 2 + STEPS):
        succ = []
        for _, _, inv, pins, relax in beam:
            near = ctx.done[inv][2]
            for p, r in _pin_moves(ctx, pins, relax, inv, near):
                t = state(p, frozenset(r))
                if t is not None and t[2] not in seen:
                    seen.add(t[2])
                    succ.append(t)
        if not succ:
            break
        succ.sort(key=lambda t: (t[0], t[1]))
        if (succ[0][0], succ[0][1]) >= (best[0], best[1]):
            break
        beam, best = succ[:BEAM], succ[0]
    draw(ctx, best[2])
    return best[2]
