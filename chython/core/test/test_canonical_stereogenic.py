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
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
"""A parity on a unit that is not stereogenic is not identity: `structure_stereogenic_view`.

| parity on                                   | `==`, `hash`, `canonical_bytes`, `str()` |
| ------------------------------------------- | ---------------------------------------- |
| a stereogenic unit                          | read                                     |
| a perceived unit `mark_stereogenic` refused | ignored; the stored parity is untouched  |
"""
from pytest import mark

from chython.core._core import read_smiles, write_smiles


def _orders(m, count=10):
    return {str(read_smiles(format(m, 'r'))) for _ in range(count)}


@mark.parametrize('stated, plain', [
    ('C[C@H](C)C', 'CC(C)C'),                               # isobutane
    ('F[C@H](Cl)Cl', 'FC(Cl)Cl'),                           # dichlorofluoromethane
    ('C[C@@H](C)[C@H](O)C', 'CC(C)[C@H](O)C'),              # 3-methylbutan-2-ol, one real centre
])
def test_a_refused_parity_is_not_identity(stated, plain):
    a, b = read_smiles(stated), read_smiles(plain)
    assert a == b and hash(a) == hash(b) and a.canonical_bytes == b.canonical_bytes
    assert str(a) == str(b) and _orders(a) == {str(b)}


def test_the_stored_parity_is_untouched():
    m = read_smiles('C[C@H](C)C')
    assert str(m) == 'C(C)(C)C'
    assert [u['parity'] for u in m.stereo_units() if u['parity']] == [2]


def test_a_canonical_string_reads_back_equal():
    """(2E,5E)-hepta-2,5-dien-4-one (1-chloroethylidene)hydrazone: the C=N on the ketone carbon has
    two identical propenyl arms, so its stated parity is not stereogenic."""
    m = read_smiles(r'C/C=C/C(=N/N=C(\Cl)/C)/C=C/C')
    m.canonicalize()
    text = str(m)
    back = read_smiles(text)
    assert back == m and str(back) == text and hash(back) == hash(m)
    assert _orders(m) == {text}


def test_a_real_centre_still_separates_enantiomers():
    a, b = read_smiles('C[C@H](O)CC'), read_smiles('C[C@@H](O)CC')
    assert a != b and write_smiles(a) != write_smiles(b)
