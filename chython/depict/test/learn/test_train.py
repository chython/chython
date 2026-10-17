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
"""N2 training end to end on a tiny set. Skipped without torch."""
import pytest

torch = pytest.importorskip('torch')

from chython.core import smiles  # noqa: E402
from chython.depict.layout.env import Env  # noqa: E402
from chython.depict.layout.greedy import best_step, drive  # noqa: E402
from chython.depict.layout.learn.demos import demonstrate, post  # noqa: E402
from chython.depict.layout.learn.policy import Policy  # noqa: E402
from chython.depict.layout.learn.replay import Trajectory  # noqa: E402
from chython.depict.layout.learn.targets import cost  # noqa: E402
from chython.depict.layout.learn.train import chooser, examples, goal_of, imitate, refine  # noqa: E402
from chython.depict.layout.rings import ring_blocks  # noqa: E402


TEXTS = ('CC(=O)Oc1ccccc1C(=O)O', 'CCCCC1CCC(CC1)c1ccccc1')


def _set():
    out = []
    for text in TEXTS:
        mol = smiles(text)
        goal = post(drive(Env(mol, set(mol), mol.canonical_order(), ring_blocks(mol), everywhere=True), best_step))
        (d, keys, comp), = demonstrate(mol, goal)
        out.append((mol, goal, Trajectory(text, text, keys, {'everywhere': True, 'component': list(comp),
                                                              'distance': d})))
    return out


def test_imitation_learns_its_targets_and_drives_to_a_stop():
    torch.manual_seed(0)
    data = [e for mol, goal, t in _set() for e in examples(t, goal)]
    assert data and all(bool(m.any()) for _, _, m, _ in data)
    model = Policy()
    losses = imitate(model, data, epochs=30, patience=30)['train']
    assert losses[-1] < losses[0]
    mol = smiles(TEXTS[0])
    env = drive(Env(mol, set(mol), mol.canonical_order(), ring_blocks(mol), everywhere=True), chooser(model))
    assert env.done and env.steps <= env.cap


def test_refinement_runs_and_rewards_the_end_not_the_route():
    torch.manual_seed(0)
    starts = []
    for mol, goal, t in _set()[:1]:
        atoms = set(mol)
        g = goal_of(mol, goal, atoms)
        starts.append((lambda m=mol: Env(m, set(m), m.canonical_order(), ring_blocks(m), everywhere=True), g))
    rewards = refine(Policy(), starts, rollouts=1)
    assert len(rewards) == 1 and rewards[0] <= 0


def test_the_evaluation_measures_both_choosers_on_one_goal():
    from chython.depict.layout.learn.evaluate import compare
    mol, goal, _ = _set()[0]
    mine = compare(mol, goal, chooser(Policy()))
    base = compare(mol, goal, None)
    assert len(mine) == len(base) == 1 and len(base[0]) == 4
    assert base[0][0] < .15 and base[0][2] == 0                 # the goal is the model-free greedy drawing itself
    assert base[0][3] == 0                                      # and it has no overlap, crossing or atom on a bond


def test_a_policy_chooser_lays_out_with_edits_everywhere():
    from chython.depict.layout.search import controller_layout
    seen = []
    choose = chooser(Policy())

    def spy(env, legal):
        seen.append(env.everywhere)
        return choose(env, legal)
    spy.everywhere = choose.everywhere
    assert choose.everywhere is True
    controller_layout(smiles(TEXTS[0]), choose=spy)
    assert seen and all(seen)


def test_a_defect_left_at_the_end_costs_more_than_any_distance_gain():
    from chython.depict.layout.learn.targets import DEFECT
    assert cost(.0, 3, defects=1) > cost(.9, 3) and DEFECT >= 1.


def test_of_equally_near_steps_the_one_repairing_more_is_the_target():
    from chython.depict.layout.learn.targets import target_set
    mol = smiles('CCCCCCCCCC')
    env = Env(mol, set(mol), mol.canonical_order(), [], everywhere=True)
    env.plane[10] = env.plane[8]
    from chython.depict.layout.score import terms
    env.terms = terms(env.scope, env.plane)
    env._legal = None

    class Flat:                                         # every drawing equally near: only repair tells them apart
        def distance(self, env, plane=None):
            return 0.
    legal = env.legal()
    mask, _ = target_set(env, legal, Flat())
    least = min(s.terms[0] for s in legal)
    assert any(mask) and all(s.terms[0] == least for s, m in zip(legal, mask) if m)


def test_training_stops_at_a_plateau_and_keeps_its_best_weights():
    from chython.depict.layout.learn.train import split, imitate
    torch.manual_seed(0)
    data = [e for mol, goal, t in _set() for e in examples(t, goal)]
    model = Policy(width=16, layers=2, heads=2)
    history = imitate(model, data, epochs=400, valid=data[:2], patience=5)
    assert len(history['train']) < 400 and history['best'] == min(history['valid'])
    train, valid = split(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j'] * 10)
    assert set(train).isdisjoint(valid) and valid


def test_a_primitive_policy_directory_loads_as_a_primitive_chooser(tmp_path):
    import json
    from chython.depict.layout.learn.train import load
    from chython.depict.layout.search import controller_layout
    model = Policy(width=16, layers=2, heads=2)
    torch.save(model.state_dict(), tmp_path / 'policy.pt')
    (tmp_path / 'manifest.json').write_text(json.dumps({'size': model.size, 'kind': 'primitive'}))
    choose = chooser(load(tmp_path))
    seen = []

    def spy(env, legal):
        seen.append(type(env).__name__)
        return choose(env, legal)
    spy.primitive = getattr(choose, 'primitive', False)
    mol = smiles(TEXTS[0])
    plane = controller_layout(mol, choose=spy)
    assert spy.primitive is True and set(plane) == set(mol) and set(seen) == {'PrimitiveEnv'}
