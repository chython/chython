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
"""Targets reward outcomes, not routes."""
import pytest

from chython.core import smiles
from chython.depict.layout.env import Env
from chython.depict.layout.greedy import best_step, drive
from chython.depict.layout.learn.demos import post
from chython.depict.layout.learn.targets import LAMBDA, TIE, cost, target_set
from chython.depict.layout.learn.train import goal_of
from chython.depict.layout.rings import ring_blocks


def _goal(mol):
    env = drive(Env(mol, set(mol), mol.canonical_order(), ring_blocks(mol), everywhere=True), best_step)
    return goal_of(mol, post(env), set(mol)), env


def test_a_route_one_step_longer_to_one_drawing_costs_exactly_one_step_more():
    assert cost(.05, 4) - cost(.05, 3) == pytest.approx(LAMBDA)


def test_two_steps_to_mirror_drawings_of_a_symmetric_molecule_are_both_targets():
    mol = smiles('CC(C)(C)c1ccc(C(C)(C)C)cc1')           # two equivalent tert-butyls
    goal, _ = _goal(mol)
    env = Env(mol, set(mol), mol.canonical_order(), ring_blocks(mol), everywhere=True)
    legal = env.legal()                                 # the ring system's variants
    mask, dist = target_set(env, legal, goal)
    best = min(dist)
    assert sum(mask) >= 2                               # a variant and its mirror reach one outcome
    assert all(m == (d <= best + TIE) for m, d in zip(mask, dist))


def test_stop_is_the_target_once_no_step_comes_nearer():
    mol = smiles('CCCCCC')
    goal, done = _goal(mol)
    env = Env(mol, set(mol), mol.canonical_order(), [], everywhere=True)
    for m in done.moves[:-1]:
        env.step(next(s for s in env.legal() if s.move == m))
    mask, _ = target_set(env, env.legal(), goal)
    assert [s.move.kind for s, t in zip(env.legal(), mask) if t] == ['stop']
