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
"""`canonicalize()` returns one representation per compound, whatever the atom order and Kekule form.

Each drawing is a random atom order of the compound, written aromatic or re-kekulised in that order --
the second gives a different Kekule form per order wherever the compound has several.  The fullerenes
and porphyrinoids are where the stored ring basis, the kekuliser's search order and a tie between
mobile-hydrogen placements each followed the atom order; the benzenoids are controls.
"""
import random

from pytest import mark

from .. import canonicalize   # the import injects the method
from ...core import read_smiles as smiles


C60 = ('c12c3c4c5c1c1c6c7c2c2c8c3c3c9c4c4c%10c5c5c1c1c6c6c%11c7c2c2c7c8c3c3c8c9c4c4c9c%10c5c5c1c1c6c6c%11'
       'c2c2c7c3c3c8c4c4c9c5c1c1c6c2c3c41')
#: PCBM, [6,6]-phenyl-C61-butyric acid methyl ester, in one Kekule form.
PCBM = ('C1=2C3=C4C=5C=6C7=C8C=9C=5C5=C3C=3C=%10C%11%12C(C%13=CC=CC=C%13)(CCCC(OC)=O)C%115C=9C=5C=9C=%11'
        'C%13=C%14C=%15C=%16C%17=C%18C=%15C=%11C=5C%12=C%18C=%10C=5C(C1=3)=C1C(C%17=5)=C3C5=C%10C=%11C=%12'
        'C(C=6C4=C4C=2C1=C5C4=%12)=C1C=%11C(=C%13C(=C71)C=98)C%14=C%10C=%163')
#: phthalocyanine, CAS 574-93-6, and 5,10,15,20-tetraphenylporphyrin, CAS 917-23-7.
PHTHALOCYANINE = 'c1ccc2c(c1)c1nc2nc2[nH]c(nc3nc(nc4[nH]c(n1)c1ccccc41)c1ccccc31)c1ccccc21'
TPP = 'c1ccc(cc1)-c1c2ccc(n2)c(-c2ccccc2)c2ccc([nH]2)c(-c2ccccc2)c2ccc(n2)c(-c2ccccc2)c2ccc1[nH]2'
PORPHINE = 'c1cc2cc3ccc(cc4ccc(cc5ccc(cc1n2)[nH]5)n4)[nH]3'

CASES = {
    'C60': C60,
    'PCBM': PCBM,
    'phthalocyanine': PHTHALOCYANINE,
    'tetraphenylporphyrin': TPP,
    'porphine': PORPHINE,
    'naphthalene': 'c1ccc2ccccc2c1',
    'pyrene': 'c1cc2ccc3cccc4ccc(c1)c2c34',
    'triphenylene': 'c1ccc2c(c1)c1ccccc1c1ccccc21',
    'biphenylene': 'c1ccc2c(c1)-c1ccccc1-2',
}


def drawings(text, count):
    base = smiles(text)
    rng = random.getstate()
    random.seed(20260928)
    try:
        for i in range(count):
            order = smiles(format(base, 'r'))
            if i % 2:
                order.kekule()
                order = smiles(format(order, 'r'))
            yield order
    finally:
        random.setstate(rng)


@mark.parametrize('name', CASES)
def test_twenty_atom_orders_give_one_representation(name):
    forms = set()
    for molecule in drawings(CASES[name], 20):
        molecule.canonicalize()
        forms.add((str(molecule), molecule.canonical_bytes))
    assert len(forms) == 1, sorted(forms)


@mark.parametrize('name', ['C60', 'PCBM', 'phthalocyanine', 'tetraphenylporphyrin', 'porphine'])
def test_the_canonical_string_reads_back_to_itself(name):
    for molecule in drawings(CASES[name], 4):
        molecule.canonicalize()
        text = str(molecule)
        again = smiles(text)
        again.canonicalize()
        assert again.canonical_bytes == molecule.canonical_bytes
        assert str(again) == text


def test_every_face_of_c60_is_aromatic_from_every_order():
    for molecule in drawings(C60, 20):
        molecule.canonicalize()
        assert all(bond.order == 4 for bond in molecule.bonds())


@mark.parametrize('text', ['CC1=CC2=CC=C3C=CC=C4C=CC(=C1)[C@]2(C)[C@@]34C',
                           'CC=1C=C2C=CC3=CC=CC4=CC=C(C=1)[C@]2(C)[C@@]34C'])
def test_a_perimeter_thiele_does_not_take_keeps_its_bonds_and_its_configuration(text):
    """2,15,16-trimethyl-15,16-dihydropyrene: the 14-atom perimeter has two matchings and no aromatic form.

    Rewriting one matching onto the other on the same atoms keeps the stated parities while moving the
    bonds they are read against, which gives another stereoisomer.  So nothing is rewritten.
    """
    molecule = smiles(text)
    before = molecule.copy()
    molecule.canonicalize()
    assert molecule == before, 'canonicalize() changed the compound'
    assert not any(r.rule == 'kekule-form:canonical' for r in molecule.log)
