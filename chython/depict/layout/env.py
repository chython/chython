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
"""The controller's environment: one component's state, its legal steps, and `step`.

| part | rule |
| --- | --- |
| state | the plane, its score and `TERMS`, broken cis/trans units, the variant drawn per ring system, steps taken |
| start | the tree embedding |
| order | while a ring system is undrawn, only its variants; then frontier edits, arc moves, redraws and `stop` |
| everywhere | edits on every atom, not only the frontier: a policy's taste acts where nothing is wrong |
| legal | an edit `masked` refuses is never offered; of one system's variants only those leaving fewest units broken |
| bend | offered only when no other legal step betters the state |
| defects | an edit adding a defect is never offered, a ring system's variant excepted; nor is `stop` while an edit \
can remove one |
| budget | past `2 × systems + STEPS` steps only `stop` |
| stop | legal once every system is drawn or has no legal variant, and no edit removes a defect |

A `Step` carries its move, the atoms it puts where, the rank the state would take, and its exact `delta_terms`: the
features a policy reads. Ranks compare by broken units first, then by `key` of the score.
"""
from typing import NamedTuple

from .edits import KINDS, Move, apply, arc_moves, broken, candidates, embed, masked, move_key, variants
from .relax import _ct_broken, _frames
from .score import TERMS, Scope, add, combine, delta_terms, frontier, key, terms
from .tree import tree_embedding


__all__ = ['STOP', 'Step', 'Env', 'PeptideEnv', 'STEPS']

STEPS = 16
STOP = Move('stop', -1, -1, 0., frozenset())


class Step(NamedTuple):
    move: Move
    moved: dict             # atom -> point; empty for `stop`
    rank: tuple             # (broken units, key of score) of the state it leads to
    terms: tuple            # its `delta_terms`
    ct: int                 # broken units after it
    state: dict = None      # what else an environment needs to take it: a peptide's pins


class Env:
    """One connected component of `mol` from its tree embedding to `stop`."""
    __slots__ = ('mol', 'scope', 'root', 'blocks', 'drawings', 'plane', 'terms', 'ct', 'drawn', 'redrawn', 'steps',
                 'cap', 'done', 'moves', 'everywhere', '_legal')

    def __init__(self, mol, atoms, rank, blocks, everywhere=False):
        self.mol = mol
        self.everywhere = everywhere
        tree = tree_embedding(Scope(mol, atoms, rank, []))
        self.scope = Scope(mol, atoms, rank, tree.closures)
        self.root = tree.root
        self.blocks = blocks
        self.drawings = [variants(self.scope, b) for b in blocks]
        self.plane = dict(tree.plane)
        self.terms = terms(self.scope, self.plane)
        self.ct = _ct_broken(mol, self.plane, self.scope.frames, set(self.plane))
        self.drawn = [-1 if d else -2 for d in self.drawings]     # -2: a system with no variant keeps its tree drawing
        self.redrawn = set()
        self.steps = 0
        self.cap = 2 * len(blocks) + STEPS
        self.done = False
        self.moves = []
        self._legal = None

    @property
    def score(self):
        return combine(self.terms)

    @property
    def rank(self):
        return self.ct, key(self.score)

    def tie(self, move):
        """The order of two equally ranked moves, `stop` last."""
        return (len(KINDS),) if move.kind == 'stop' else move_key(self.scope, move)

    def _steps(self, moves):
        out, least = [], {}
        for m in moves:
            if m.kind == 'embed':
                moved = embed(self.scope, self.plane, self.blocks[m.pivot][0], self.drawings[m.pivot][m.other],
                              self.root)
            else:
                moved = apply(self.scope, self.plane, m)
            if masked(self.scope, self.plane, moved):
                continue
            ct = self.ct - _ct_broken(self.mol, self.plane, self.scope.frames, set(moved)) + \
                broken(self.scope, self.plane, moved)
            d = delta_terms(self.scope, self.plane, moved)
            if m.kind == 'embed':
                least[m.pivot] = min(least.get(m.pivot, ct), ct)
            out.append(Step(m, moved, (ct, key(add(self.score, combine(d)))), d, ct))
        return [s for s in out if s.move.kind != 'embed' or s.ct == least[s.move.pivot]]

    def legal(self):
        """The legal steps of the state, `stop` last once it is legal; computed once per state."""
        if self._legal is not None:
            return self._legal
        if self.done:
            self._legal = []
            return self._legal
        while -1 in self.drawn:
            i = self.drawn.index(-1)
            out = self._steps([Move('embed', i, j, 0., frozenset()) for j in range(len(self.drawings[i]))])
            if out:
                self._legal = out
                return out
            self.drawn[i] = -2                                  # every variant masked: the tree drawing stays
        stop = Step(STOP, {}, self.rank, (0,) * len(TERMS), self.ct)
        if self.steps >= self.cap:
            self._legal = [stop]
            return self._legal
        front = set(self.scope.atoms) if self.everywhere else frontier(self.scope, self.plane)
        moves = candidates(self.scope, self.plane, front) + arc_moves(self.scope, self.plane, front) if front else []
        for i, (atoms, _) in enumerate(self.blocks):
            if self.drawn[i] >= 0 and i not in self.redrawn and atoms & front:
                moves += [Move('embed', i, j, 0., frozenset()) for j in range(len(self.drawings[i]))
                          if j != self.drawn[i]]
        out = [s for s in self._steps(moves) if s.move.kind == 'embed' or s.terms[0] <= 0]
        if front and not any(s.rank < self.rank for s in out):
            out += [s for s in self._steps(candidates(self.scope, self.plane, front, violate=True)) if s.terms[0] <= 0]
        if any(s.terms[0] < 0 for s in out):                # a defect that can be repaired must be
            self._legal = out
        else:
            self._legal = out + [stop]
        return self._legal

    def fork(self):
        """A copy whose steps leave this one as it is: the plane, terms, counts and moves copied, the molecule, scope
        and variants shared, no legal steps cached."""
        f = Env.__new__(Env)
        f.mol, f.scope, f.root, f.blocks, f.drawings = self.mol, self.scope, self.root, self.blocks, self.drawings
        f.plane, f.terms, f.ct = dict(self.plane), self.terms, self.ct
        f.drawn, f.redrawn, f.steps, f.cap = list(self.drawn), set(self.redrawn), self.steps, self.cap
        f.done, f.moves, f.everywhere, f._legal = self.done, list(self.moves), self.everywhere, None
        return f

    def step(self, s):
        """Take the legal step `s`."""
        m = s.move
        if m.kind == 'stop':
            self.done = True
        else:
            self.plane.update(s.moved)
            self.terms = tuple(a + b for a, b in zip(self.terms, s.terms))
            self.ct = s.ct
            if m.kind == 'embed':
                if self.drawn[m.pivot] >= 0:
                    self.redrawn.add(m.pivot)
                self.drawn[m.pivot] = m.other
            self.steps += 1
        self.moves.append(m)
        self._legal = None


class PeptideEnv:
    """A peptide backbone from no pin to `stop`: each step pins vertices, as `search._pin_moves` offers them, and
    `peptide.solve` completes the flattest backbone around the pins. A pin drawing more stated cis/trans units against
    their parity than the state does is never offered. A step's `terms` is the change of the five places of
    `search._peptide_score`, its `moved` the atoms the new drawing puts elsewhere, its `state` the pins and the
    drawing. `scope`, `blocks`, `drawn`, `terms` and `ct` let `features` read it as it reads an `Env`."""
    __slots__ = ('ctx', 'mol', 'scope', 'blocks', 'drawn', 'frames', 'pins', 'relax', 'inv', 'plane', 'score', 'ct',
                 'steps', 'cap', 'done', 'moves', 'seen', '_legal')

    def __init__(self, ctx):
        from .peptide import solve
        from .search import STEPS as BUDGET, _peptide_score
        self.ctx = ctx
        self.mol = ctx.mol
        self.scope = Scope(ctx.mol, set(ctx.mol), ctx.rank, [])
        self.blocks, self.drawn = [], []
        self.frames = _frames(ctx.mol, {x: None for x in ctx.mol})
        self.pins, self.relax = {}, frozenset()
        self.inv = solve(ctx, self.pins, self.relax)
        self.plane, self.score = _peptide_score(ctx, self.inv)
        self.ct = 0 if self.plane is None else _ct_broken(self.mol, self.plane, self.frames, set(self.mol))
        self.steps = 0
        self.cap = len(ctx.peptide.backbone) // 2 + BUDGET
        self.done = self.plane is None
        self.moves = []
        self.seen = {self.inv}
        self._legal = None

    @property
    def terms(self):
        return self.score + (0,) * (len(TERMS) - len(self.score))     # read as an `Env`'s by `features`

    @property
    def rank(self):
        return self.ct, self.score

    def tie(self, move):
        return move.arg

    def legal(self):
        """The pins one step away whose drawing fits, breaks no further cis/trans and differs from every one seen,
        then `stop`."""
        from .peptide import solve
        from .search import _peptide_score, _pin_moves
        if self._legal is not None:
            return self._legal
        if self.done:
            self._legal = []
            return self._legal
        out = []
        if self.steps < self.cap:
            near = self.ctx.done[self.inv][2]
            for q, (pins, relax) in enumerate(_pin_moves(self.ctx, self.pins, self.relax, self.inv, near)):
                relax = frozenset(relax)
                inv = solve(self.ctx, pins, relax)
                if inv in self.seen:
                    continue
                plane, score = _peptide_score(self.ctx, inv)
                if plane is None:
                    continue
                ct = _ct_broken(self.mol, plane, self.frames, set(self.mol))
                if ct > self.ct:
                    continue
                moved = {x: p for x, p in plane.items() if self.plane.get(x) != p}
                m = Move('pin', -1, -1, float(q), frozenset(moved))
                d = tuple(b - a for a, b in zip(self.score, score)) + (0,) * (len(TERMS) - len(score))
                out.append(Step(m, moved, (ct, score), d, ct,
                                {'pins': pins, 'relax': relax, 'inv': inv, 'plane': plane, 'score': score}))
        self._legal = out + [Step(STOP, {}, self.rank, (0,) * len(TERMS), self.ct)]
        return self._legal

    def step(self, s):
        """Take the legal step `s`."""
        if s.move.kind == 'stop':
            self.done = True
        else:
            self.pins, self.relax, self.inv, self.plane, self.score = \
                (s.state[k] for k in ('pins', 'relax', 'inv', 'plane', 'score'))
            self.ct = s.ct
            self.seen.add(self.inv)
            self.steps += 1
        self.moves.append(s.move)
        self._legal = None
