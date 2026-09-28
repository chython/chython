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
"""`thiele()` candidates and `kekule(canonical=True)` do not follow the atom order.

| what                         | fixed by                                                        |
|------------------------------|-----------------------------------------------------------------|
| `thiele()`'s candidate rings | every relevant cycle of at most seven atoms, not the ring basis |
| `kekule(canonical=True)`     | the search walks atoms and neighbours by canonical rank         |

`thiele()` still reads the Kekule form it is handed, so a porphyrinoid's aromatic form follows the form
plain `kekule()` wrote; C60 and PCBM come out in one aromatic form from any form.
"""
import random

from pytest import mark, raises

from chython.core import read_smiles as smiles


C60 = ('c12c3c4c5c1c1c6c7c2c2c8c3c3c9c4c4c%10c5c5c1c1c6c6c%11c7c2c2c7c8c3c3c8c9c4c4c9c%10c5c5c1c1c6c6c%11'
       'c2c2c7c3c3c8c4c4c9c5c1c1c6c2c3c41')
#: PCBM, [6,6]-phenyl-C61-butyric acid methyl ester, in one Kekule form.
PCBM = ('C1=2C3=C4C=5C=6C7=C8C=9C=5C5=C3C=3C=%10C%11%12C(C%13=CC=CC=C%13)(CCCC(OC)=O)C%115C=9C=5C=9C=%11'
        'C%13=C%14C=%15C=%16C%17=C%18C=%15C=%11C=5C%12=C%18C=%10C=5C(C1=3)=C1C(C%17=5)=C3C5=C%10C=%11C=%12'
        'C(C=6C4=C4C=2C1=C5C4=%12)=C1C=%11C(=C%13C(=C71)C=98)C%14=C%10C=%163')
#: phthalocyanine, CAS 574-93-6.
PHTHALOCYANINE = 'c1ccc2c(c1)c1nc2nc2[nH]c(nc3nc(nc4[nH]c(n1)c1ccccc41)c1ccccc31)c1ccccc21'


def orders(text, count):
    base = smiles(text)
    rng = random.getstate()
    random.seed(20260928)
    try:
        for _ in range(count):
            yield smiles(format(base, 'r'))
    finally:
        random.setstate(rng)


@mark.parametrize('text', [C60, PCBM], ids=['C60', 'PCBM'])
def test_a_fullerene_aromatises_to_one_form_from_any_order_and_any_kekule_form(text):
    forms = set()
    for molecule in orders(text, 20):
        molecule.kekule()
        molecule.thiele()
        forms.add(molecule.canonical_bytes)
    assert len(forms) == 1


def test_c60_has_all_thirty_two_faces_as_candidates():
    for molecule in orders(C60, 20):
        molecule.kekule()
        molecule.thiele()
        assert sum(bond.order == 4 for bond in molecule.bonds()) == 90


@mark.parametrize('text', [C60, PCBM, PHTHALOCYANINE, 'c1ccc2ccccc2c1', 'c1cc2ccc3cccc4ccc(c1)c2c34'],
                  ids=['C60', 'PCBM', 'phthalocyanine', 'naphthalene', 'pyrene'])
def test_the_canonical_kekule_form_is_one_form(text):
    forms = set()
    for molecule in orders(text, 20):
        if not any(bond.order == 4 for bond in molecule.bonds()):
            molecule.thiele()
        result = molecule.kekule(canonical=True)
        assert not result.unresolved
        forms.add(molecule.canonical_bytes)
    assert len(forms) == 1


def test_canonical_is_keyword_only():
    with raises(TypeError):
        smiles('c1ccccc1').kekule(None, None, True)
