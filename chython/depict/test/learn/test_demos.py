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
"""N1 demonstrations: reaching a liked drawing from the tree embedding."""
from math import cos, sin

import pytest

from chython.core import smiles
from chython.depict.layout.env import Env
from chython.depict.layout.greedy import best_step, drive
from chython.depict.layout.learn.demos import aligned, demonstrate, post
from chython.depict.layout.learn.replay import replay
from chython.depict.layout.rings import ring_blocks


def test_aligned_distance_ignores_a_turn_a_mirror_and_a_shift():
    a = [(0., 0.), (1., 0.), (1.5, .8), (3., 1.)]
    for t, m in ((.7, 1), (2.1, -1)):
        b = [(cos(t) * x - sin(t) * m * y + 5, sin(t) * x + cos(t) * m * y - 2) for x, y in a]
        assert aligned(a, b) == pytest.approx(0, abs=1e-9)
    assert aligned(a, [(0., 0.), (1., 0.), (1.5, -.8), (2., 2.)]) > .1


@pytest.mark.parametrize('text', ['CC(=O)Oc1ccccc1C(=O)O', 'CCCCC1CCC(CC1)c1ccccc1'])
def test_a_reachable_drawing_is_found_and_replays_to_itself(text):
    mol = smiles(text)
    rank = mol.canonical_order()
    goal = drive(Env(mol, set(mol), rank, ring_blocks(mol), everywhere=True), best_step)
    target = post(goal)
    found = demonstrate(mol, target)
    assert len(found) == 1
    d, keys, comp = found[0]
    assert d < .15
    again = post(replay(mol, keys, everywhere=True, component=comp))
    atoms = sorted(mol, key=rank.__getitem__)
    from chython.depict.layout.learn.demos import aligned as al
    assert al([again[x] for x in atoms], [target[x] for x in atoms]) == pytest.approx(d, abs=1e-6)


def test_each_component_is_its_own_demonstration():
    mol = smiles('[Na+].[O-]C(=O)c1ccccc1')
    rank = mol.canonical_order()
    target = {}
    for c in mol.connected_components:
        env = drive(Env(mol, set(c), rank, [b for b in ring_blocks(mol) if b[0] <= set(c)], everywhere=True),
                    best_step)
        target.update(post(env))
    found = demonstrate(mol, target)
    assert len(found) == 2 and all(d < .15 for d, _, _ in found)


def test_an_isotope_does_not_stop_a_demonstration():
    mol = smiles('[2H]C([2H])([2H])c1ccccc1')
    goal = post(drive(Env(mol, set(mol), mol.canonical_order(), ring_blocks(mol), everywhere=True), best_step))
    (d, _, _), = demonstrate(mol, goal)
    assert d < .15


def test_a_rerun_skips_what_is_done(tmp_path):
    from chython.depict.layout.learn.demos import pending
    from chython.depict.layout.learn.replay import Trajectory, write
    write(tmp_path, [Trajectory('a', 'C', [], {'component': [0], 'distance': 0.})])
    assert pending(tmp_path, [('a', 'C', []), ('b', 'CC', [])]) == [('b', 'CC', [])]


def test_a_huge_record_keeps_within_its_budget():
    from time import perf_counter
    mol = smiles('C' * 150 + 'c1ccccc1' + 'C' * 150)
    goal = {x: (float(i), 0.) for i, x in enumerate(mol)}
    t = perf_counter()
    found = demonstrate(mol, goal, budget=2.)
    assert perf_counter() - t < 20 and len(found) == 1


def test_a_primitive_demonstration_reaches_a_reachable_drawing_and_replays():
    from chython.depict.layout.learn.demos import primitive_demonstration
    from chython.depict.layout.learn.mcts import rebuild, step_key
    from chython.depict.layout.primitive import Initial
    mol = smiles('CCCCC1CCC(CC1)c1ccccc1')
    init = Initial(mol, set(mol))
    goal = drive(init.env(), best_step).plane
    distance, records = primitive_demonstration(init, [[goal[x] for x in init.scope.atoms]])
    assert distance < .15 and records and records[-1]['visits'] == {repr(('stop',)): 1.}
    for r in records:
        (k, p), = r['visits'].items()
        assert p == 1. and k in {repr(step_key(rebuild(init, r), s)) for s in rebuild(init, r).legal()}
