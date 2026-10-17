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
"""Tree-search self-play over the primitive environment. Skipped without torch."""
import pytest

torch = pytest.importorskip('torch')

from chython.core import smiles  # noqa: E402
from chython.depict.layout.learn.mcts import Search, episode, outcome, penalty, rebuild, step_key  # noqa: E402
from chython.depict.layout.learn.policy import Policy  # noqa: E402
from chython.depict.layout.primitive import Initial  # noqa: E402
from chython.depict.layout.score import TERMS  # noqa: E402


def _closure(env):
    return env.terms[TERMS.index('closure')]


def test_an_episode_records_visit_distributions_and_its_outcome():
    torch.manual_seed(0)
    init = Initial(smiles('C1CCCCC1C'), set(smiles('C1CCCCC1C')))
    search = Search(Policy(width=16, layers=2, heads=2), lambda env: -penalty(env), sims=8)
    records, final = episode('test:a', init, search, seed=0)
    assert final.done and records
    for r in records:
        assert abs(sum(r['visits'].values()) - 1) < 1e-9 and r['outcome'] == records[-1]['outcome']


def test_an_open_ring_costs_a_quarter_of_its_closure_error_or_more():
    mol = smiles('C1CCCCCC1')                          # cycloheptane: the honeycomb tree leaves it open
    env = Initial(mol, set(mol)).env()
    assert env.terms[0] == 0 and penalty(env) >= _closure(env) / 4


def test_a_recorded_state_rebuilds_with_the_same_legal_steps():
    torch.manual_seed(0)
    mol = smiles('CC(=O)Oc1ccccc1C(=O)O')
    init = Initial(mol, set(mol))
    search = Search(Policy(width=16, layers=2, heads=2), lambda env: -penalty(env), sims=4)
    records, _ = episode('test:c', init, search, seed=0)
    r = records[len(records) // 2]
    env = rebuild(init, r)
    keys = {repr(step_key(env, s)) for s in env.legal()}
    assert set(r['visits']) <= keys


def test_a_drawing_with_a_defect_scores_below_every_clean_one_whatever_the_judge():
    mol = smiles('C1CCCCCCCCCCC1')                     # a macrocycle: no template, the tree start overlaps
    init = Initial(mol, set(mol))
    bad = init.env()
    assert bad.terms[0] > 0
    from chython.depict.layout.greedy import best_step, drive
    good = drive(init.env(), best_step)
    assert good.terms[0] == 0
    assert outcome(bad, judge=1e6) < outcome(good, judge=-1e6, distance=1e6)


def test_the_outcome_is_bounded():
    mol = smiles('C1CCCCCCCCCCCCCC1')
    init = Initial(mol, set(mol))
    env = init.env()
    for j, d in ((1e9, 0.), (-1e9, 1e9), (0., None)):
        assert -7 <= outcome(env, judge=j, distance=d) <= .5
