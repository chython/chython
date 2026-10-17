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
"""Tree-search self-play over the primitive environment, AlphaZero's way.

| part | rule |
| --- | --- |
| node | a forked `PrimitiveEnv`, its legal steps with the policy's priors, and per step a visit count and a value sum |
| select | PUCT: the step of highest Q + `C` × P × √N / (1 + n), ties to the first in `legal` order |
| leaf | the policy's value; a `stop` leaf is the true reward of its drawing |
| reward | `reward(env)` of the finished drawing, as `outcome` composes it |
| outcome | `TASTE` × tanh(judge / 2) − `penalty` − `STEP` per step − `LIKED` × min(distance, 2) / 2 where a liked \
drawing is the goal; within [−7, 0.5], so a value head can regress it |
| penalty | `HARD` + `DEFECT` per defect and broken cis/trans unit, up to 10, when there is any: below every clean \
drawing whatever the judge, the distance and the steps say; then `SOFT` × tanh(closure in `L` + `CONVENTION` \
per point) |
| move | sampled from the visit counts for the first `explore` moves, then the most visited; the subtree is kept |
| record | per move the source, the plane in rank order, the steps taken, the visit share per `step_key`, and the \
episode's final reward as `outcome` |
"""
import random
from math import sqrt, tanh

import torch

from ..score import combine, terms
from ..relax import _ct_broken
from .tensors import state, steps


__all__ = ['C', 'DEFECT', 'CONVENTION', 'penalty', 'outcome', 'step_key', 'Search', 'episode', 'rebuild']

C = 1.5
HARD = 3.5
DEFECT = .1
CONVENTION = .2
SOFT = .5
TASTE = .5
STEP = .002
LIKED = 1.
EXPLORE = 8


def penalty(env):
    """What a finished drawing's metrics cost: defects and broken cis/trans first, closure, then conventions."""
    defects, closure, conventions = combine(env.terms)
    hard = getattr(env, 'defects', defects) + env.ct
    return (HARD + DEFECT * min(hard, 10) if hard else 0.) + SOFT * tanh(closure + CONVENTION * conventions)


def outcome(env, judge=0., distance=None):
    """The reward of a finished drawing: the judge's `judge` score, `penalty`, its steps, and its `distance` in `L` to
    a liked drawing when there is one."""
    out = TASTE * tanh(max(-50., min(50., judge)) / 2) - penalty(env) - STEP * env.steps
    if distance is not None:
        out -= LIKED * min(distance, 2.) / 2
    return out


def step_key(env, step):
    """A legal step's key, the same in any atom order."""
    m = step.move
    if m.kind == 'stop':
        return ('stop',)
    rank = env.scope.rank
    return m.kind, rank.get(m.pivot, -1), rank.get(m.other, -1), round(m.arg, 9)


class _Node:
    __slots__ = ('env', 'legal', 'prior', 'visits', 'value', 'children', 'final')

    def __init__(self, env):
        self.env = env
        self.legal = env.legal()
        self.prior = self.visits = self.value = None
        self.children = {}
        self.final = None                               # the reward of a finished drawing


class Search:
    """PUCT search with `model`'s priors and value, `sims` simulations per move."""
    def __init__(self, model, reward, sims=32, c=C, device='cpu'):
        self.model, self.reward, self.sims, self.c, self.device = model.to(device).eval(), reward, sims, c, device

    def _evaluate(self, node):
        if not node.legal:
            node.final = self.reward(node.env)
            return node.final
        with torch.no_grad():
            logits, value = self.model({k: v.to(self.device) for k, v in state(node.env).items()},
                                       {k: v.to(self.device) for k, v in steps(node.env, node.legal).items()})
        node.prior = torch.softmax(logits, 0).tolist()
        node.visits = [0] * len(node.legal)
        node.value = [0.] * len(node.legal)
        return float(value)

    def _select(self, node):
        total = sqrt(sum(node.visits) + 1)
        best = None
        for i, (p, n, w) in enumerate(zip(node.prior, node.visits, node.value)):
            u = (w / n if n else 0.) + self.c * p * total / (1 + n)
            if best is None or u > best[0] + 1e-12:
                best = (u, i)
        return best[1]

    def _child(self, node, i):
        if i not in node.children:
            env = node.env.fork()
            env.step(node.legal[i])
            child = _Node(env)
            node.children[i] = child
            return child, True
        return node.children[i], False

    def run(self, root):
        if root.prior is None and root.final is None:
            self._evaluate(root)
        for _ in range(self.sims):
            node, path = root, []
            while True:
                if node.final is not None:
                    v = node.final
                    break
                i = self._select(node)
                path.append((node, i))
                node, new = self._child(node, i)
                if new:
                    v = self._evaluate(node)
                    break
            for n, i in path:
                n.visits[i] += 1
                n.value[i] += v
        return root


def _plane(env):
    return [list(env.plane[x]) for x in env.scope.atoms]


def episode(source, initial, search, seed=0, explore=EXPLORE):
    """`(records, final env)` of one self-played drawing of `initial`."""
    rng = random.Random(seed)
    root = _Node(initial.env())
    records = []
    while root.final is None and root.legal:
        search.run(root)
        total = sum(root.visits)
        share = [n / total for n in root.visits] if total else [1. / len(root.legal)] * len(root.legal)
        records.append({'source': source, 'plane': _plane(root.env), 'steps': root.env.steps,
                        'visits': {step_key(root.env, s): p for s, p in zip(root.legal, share) if p}})
        if len(records) <= explore:
            i = rng.choices(range(len(share)), weights=share)[0]
        else:
            top = max(share)
            i = min((i for i, p in enumerate(share) if p >= top - 1e-12),
                    key=lambda i: root.env.tie(root.legal[i].move))
        root, _ = search._child(root, i)
        if root.final is None and not root.legal:
            root.final = search.reward(root.env)
    final = root.env
    outcome = search.reward(final) if root.final is None else root.final
    for r in records:
        r['outcome'] = outcome
        r['visits'] = {repr(k): p for k, p in r['visits'].items()}
    return records, final


def rebuild(initial, record):
    """The environment of a `record`'s state."""
    env = initial.env()
    env.plane = {x: tuple(p) for x, p in zip(env.scope.atoms, record['plane'])}
    env.steps = record['steps']
    env.terms = terms(env.scope, env.plane)
    env.ct = _ct_broken(env.mol, env.plane, env.scope.frames, set(env.plane))
    return env
