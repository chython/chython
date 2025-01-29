# -*- coding: utf-8 -*-
#
#  Copyright 2025 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2025 Tagir Akhmetshin <tagirshin@gmail.com>
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
from chython import smiles
from chython.exceptions import InvalidAromaticRing
from pytest import raises


def test_kekule_basic():
    # Test basic aromatic ring conversion
    mol = smiles('c1ccccc1')  # benzene
    assert mol.kekule()  # should return True for aromatic rings
    assert mol == smiles('C1=CC=CC=C1')


def test_kekule_pyridine():
    # Test pyridine and its derivatives
    mol = smiles('n1ccccc1')  # pyridine
    assert mol.kekule()
    assert mol == smiles('N1=CC=CC=C1')
    assert mol.atom(1).implicit_hydrogens == 0

    # Test protonated pyridine
    mol = smiles('[nH+]1ccccc1')
    assert mol.kekule()
    assert mol == smiles('[NH+]1=CC=CC=C1')


def test_kekule_pyrrole():
    # Test pyrrole and its derivatives
    mol = smiles('[nH]1cccc1')  # pyrrole
    assert mol.kekule()
    assert mol == smiles('N1C=CC=C1')
    assert mol.atom(1).implicit_hydrogens == 1

    mol = smiles('n1cccc1')
    assert mol.kekule()
    assert mol == smiles('N1C=CC=C1')
    assert mol.atom(1).implicit_hydrogens == 1

    # Test N-methylpyrrole
    mol = smiles('Cn1cccc1')
    assert mol.kekule()
    assert mol == smiles('CN1C=CC=C1')
    assert mol.atom(2).implicit_hydrogens == 0


def test_kekule_furan_thiophene():
    # Test oxygen and sulfur containing aromatics
    mol = smiles('o1cccc1')
    assert mol.kekule()
    assert mol == smiles('O1C=CC=C1')
    assert mol.atom(1).implicit_hydrogens == 0

    mol = smiles('s1cccc1')
    assert mol.kekule()
    assert mol == smiles('S1C=CC=C1')
    assert mol.atom(1).implicit_hydrogens == 0


def test_kekule_complex_systems():
    # Test fused ring systems
    mol = smiles('c1ccc2ccccc2c1')
    assert mol.kekule()
    assert mol == smiles('C1=CC2=CC=CC=C2C=C1') or mol == smiles('C1=CC2=C(C=C1)C=CC=C2')

    # Test indole
    mol = smiles('c1ccc2[nH]ccc2c1')
    assert mol.kekule()
    assert mol == smiles('N1C=CC2=C1C=CC=C2') or mol == smiles('N1C=CC2=CC=CC=C12')


def test_kekule_enumeration():
    mol = smiles('Cc1ccccc1C')
    forms = list(mol.enumerate_kekule())
    assert len(forms) == 2  # benzene has 2 Kekule forms
    assert smiles('CC1=C(C)C=CC=C1') in forms
    assert smiles('CC1=CC=CC=C1C') in forms


def test_kekule_invalid_structures():
    # Test invalid aromatic structures
    with raises(InvalidAromaticRing):
        mol = smiles('c1cccc1')  # 5-membered carbon ring (invalid aromatic)
        mol.kekule()

    with raises(InvalidAromaticRing):
        mol = smiles('c1ccc2c1c3ccccc3cc2')  # acenaphthalene (invalid aromatic form)
        mol.kekule()

    with raises(InvalidAromaticRing):
        mol = smiles('c1cccc1C(=O)c1cccc1')  # cyclopentadiene with carbonyl (invalid aromatic)
        mol.kekule()


def test_kekule_charged_species():
    # Test charged aromatic species
    mol = smiles('[n+]1ccccc1')
    assert mol.kekule()
    assert mol == smiles('C=1[NH+]=CC=CC=1')

    mol = smiles('[cH-]1cccc1')
    assert mol.kekule()
    assert mol == smiles('C=1C=C[CH-]C=1')


def test_kekule_multiple_rings():
    # Test molecules with multiple aromatic rings
    mol = smiles('c1ccccc1c2ccccc2')
    assert mol.kekule()
    assert mol == smiles('C1=CC=C(C=C1)C1=CC=CC=C1')


def test_kekule_heteroatoms():
    # Test various heteroatoms in aromatic rings
    mol = smiles('c1cncn1')  # two nitrogens
    assert mol.kekule()
    assert mol == smiles('N1C=CN=C1')

    mol = smiles('o1cncc1')  # oxygen and nitrogen
    assert mol.kekule()
    assert mol == smiles('C1=COC=N1')


def test_kekule_buffer_size():
    # Test buffer size parameter for complex heterocycles
    mol = smiles('c1ccc2[nH]ccc2c1')  # indole
    assert mol.kekule(buffer_size=1)  # small buffer

    mol = smiles('c1ccc2[nH]ccc2c1')  # fresh indole instance
    assert mol.kekule(buffer_size=10)  # large buffer


def test_kekule_radical_species():
    mol = smiles('[c]1ccccc1')
    assert mol.kekule()
    assert mol == smiles('C=1C=CC=[C]C=1 |^1:4|')


def test_kekule_quinones():
    # Test quinone-like structures
    mol = smiles('O=c1ccc(=O)cc1')
    assert mol.kekule()
    assert mol == smiles('C1=CC(C=CC1=O)=O')
