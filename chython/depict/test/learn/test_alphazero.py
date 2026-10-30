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
"""The AlphaZero rounds: tester, learner, judge round, orchestration. Skipped without torch."""
import json

import pytest

torch = pytest.importorskip('torch')

from chython.core import smiles  # noqa: E402
from chython.depict.layout.learn import alphazero as az  # noqa: E402


def test_descriptors_are_fixed_width_and_order_blind():
    mol = smiles('CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C')
    a = az.descriptors(mol)
    assert len(a) == az.DESCRIPTORS and a == az.descriptors(smiles(format(mol, 'r')))


def test_the_tester_samples_hard_molecules_more_and_never_only_them():
    t = az.Tester()
    pool = [('easy', smiles('CCO')), ('hard', smiles('CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C'))]
    t.fit([(az.descriptors(m), 0. if s == 'easy' else -10.) for s, m in pool] * 20, epochs=200)
    picks = t.pick(pool, 400, seed=0)
    assert picks.count('hard') > picks.count('easy') > 0


def test_a_tiny_run_plays_learns_judges_and_resumes(tmp_path):
    corpus = tmp_path / 'corpus.tsv'
    corpus.write_text('id\tstratum\tsmiles\tcore\npub:a\tdrug\tCC(=O)Oc1ccccc1C(=O)O\t\npub:b\tdrug\tC1CCCCC1CCO\t\n')
    poses = tmp_path / 'poses.tsv'
    from chython.depict.layout.greedy import best_step, drive
    from chython.depict.layout.primitive import Initial
    lines = []
    for i, text in (('pub:a', 'CC(=O)Oc1ccccc1C(=O)O'), ('pub:b', 'C1CCCCC1CCO')):
        mol = smiles(text)
        env = drive(Initial(mol, set(mol)).env(), best_step)
        lines.append(f'{i}\t{json.dumps([list(env.plane[x]) for x in mol])}')
    poses.write_text('\n'.join(lines) + '\n')
    cfg = az.Config(corpus=str(corpus), poses=str(poses), dislikes=str(tmp_path / 'none.tsv'),
                    molecules=2, sims=4, workers=0, width=16, layers=2, heads=2, epochs=3, patience=2, window=2)
    run = tmp_path / 'run'
    az.rounds(run, cfg, count=2)
    state = json.loads((run / 'state.json').read_text())
    assert state['round'] == 2
    assert (run / 'round-0000' / 'policy.pt').exists() and (run / 'demos.jsonl').read_text()
    metrics = json.loads((run / 'round-0002' / 'metrics.json').read_text())
    assert {'selfplay', 'policy', 'judge', 'tester', 'eval'} <= set(metrics)
    az.rounds(run, cfg, count=1)                       # resumes after round 2
    assert json.loads((run / 'state.json').read_text())['round'] == 3
