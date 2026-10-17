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
"""Targets that reward outcomes, not routes: at a state, every legal step whose drawing comes as near the goal as the
best one does is correct.

| part | rule |
| --- | --- |
| distance | `demos.aligned` of a drawing, post-stepped, to the goal, the least over the goal's automorphism images |
| target set | of the legal steps within `TIE` of the nearest step's distance, those removing the most defects; `stop` \
when no step comes nearer than the state is |
| same | two drawings within `SAME` of each other are one outcome, whatever their routes |
| cost | the final distance, `LAMBDA` per step and `DEFECT` per defect left: the shorter of two routes to one outcome \
costs less, and an overlap costs more than any distance |
"""
from .demos import aligned, post


__all__ = ['TIE', 'SAME', 'LAMBDA', 'DEFECT', 'Goal', 'target_set', 'cost']

TIE = .02               # in L: steps this near the best one's distance are equally good
SAME = .15              # in L: two drawings this near are one outcome
LAMBDA = .01            # in L per step
DEFECT = 1.             # in L per defect left in the final drawing: more than any distance a drawing gains


class _View:
    """An environment's frame with another plane: what `post` reads."""
    __slots__ = ('mol', 'plane', 'blocks', 'scope')

    def __init__(self, env, plane):
        self.mol, self.blocks, self.scope, self.plane = env.mol, env.blocks, env.scope, plane


class Goal:
    """A liked drawing of one component, as its automorphism images in the component's rank order."""
    __slots__ = ('atoms', 'images')

    def __init__(self, atoms, images):
        self.atoms = atoms
        self.images = images

    def distance(self, env, plane=None):
        p = post(env if plane is None else _View(env, {**env.plane, **plane}))
        mine = [p[x] for x in self.atoms]
        return min(aligned(mine, g) for g in self.images)


def target_set(env, legal, goal):
    """`(mask, distances)`: per step of `legal`, whether it is a target, and its drawing's distance to `goal`."""
    here = goal.distance(env)
    dist = [here if s.move.kind == 'stop' else goal.distance(env, s.moved) for s in legal]
    moves = [d for d, s in zip(dist, legal) if s.move.kind != 'stop']
    best = min(moves) if moves else here
    if best >= here - 1e-9 and any(s.move.kind == 'stop' for s in legal) and -1 not in getattr(env, 'drawn', ()):
        return [s.move.kind == 'stop' for s in legal], dist
    near = [s.move.kind != 'stop' and d <= best + TIE for d, s in zip(dist, legal)]
    least = min(s.terms[0] for s, n in zip(legal, near) if n)      # of the equally near, those repairing most
    return [n and s.terms[0] == least for s, n in zip(legal, near)], dist


def cost(distance, steps, defects=0):
    """The cost of a route ending `distance` from the goal after `steps` steps with `defects` left."""
    return distance + LAMBDA * steps + DEFECT * defects
