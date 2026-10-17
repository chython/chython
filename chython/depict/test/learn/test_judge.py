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
"""The judge: a preference network over finished drawings. Skipped without torch."""
import pytest

torch = pytest.importorskip('torch')

from chython.core import smiles  # noqa: E402
from chython.depict.layout.greedy import best_step, drive  # noqa: E402
from chython.depict.layout.learn.judge import Judge, judge_score, pairs_of, perturbed, train_judge  # noqa: E402
from chython.depict.layout.primitive import Initial  # noqa: E402


TEXTS = ('CC(=O)Oc1ccccc1C(=O)O', 'CCCCC1CCC(CC1)c1ccccc1', 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O', 'OC1CCCCC1N')


def _liked(text):
    mol = smiles(text)
    init = Initial(mol, set(mol))
    return init, drive(init.env(), best_step).plane


def test_a_perturbed_drawing_differs_and_keeps_every_atom():
    init, liked = _liked(TEXTS[0])
    bad = perturbed(init, liked, steps=3, seed=1)
    assert set(bad) == set(liked) and bad != liked


def test_the_judge_learns_to_prefer_liked_drawings():
    torch.manual_seed(0)
    groups = []
    for text in TEXTS:
        init, liked = _liked(text)
        groups.append((init, liked, [perturbed(init, liked, steps=4, seed=s) for s in range(4)]))
    data = [p for g in groups for p in pairs_of(*g)]
    model = Judge(width=16, layers=2, heads=2)
    history = train_judge(model, data, epochs=200, valid=data, patience=20)
    assert history['accuracy'][-1] >= .9
    init, liked = _liked(TEXTS[0])
    assert judge_score(model, init, liked) > judge_score(model, init, perturbed(init, liked, steps=4, seed=99))
