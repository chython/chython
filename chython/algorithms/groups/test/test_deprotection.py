# -*- coding: utf-8 -*-
#
#  Copyright 2022-2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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
import pytest
from chython import smiles
from chython.algorithms.groups._protective import rules


_test_cases = [(name, *rule) for name, rule in rules.items()]


@pytest.mark.parametrize('name,q,keep,add,test_smi,cleaved_smi,decoys', _test_cases, ids=[x[0] for x in _test_cases])
def test_deprotection(name, q, keep, add, test_smi, cleaved_smi, decoys):
    t = smiles(test_smi)
    t.canonicalize()
    # query matches protected form
    assert q < t, f'query does not match protected form'
    # deprotection produces expected product
    a = smiles(cleaved_smi)
    a.canonicalize()
    t.remove_protection(name)
    assert t == a, f'deprotection gave {t}, expected {a}'
    # decoys are not matched
    for d in decoys:
        d = smiles(d)
        d.canonicalize()
        assert not q < d, f'query matches decoy {d}'
