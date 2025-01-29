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


def test_basic_thiele():
    # Test basic aromatic systems
    mol = smiles('C1=CC=CC=C1')  # benzene in Kekule form
    assert mol.thiele()

    # Check that all bonds in the ring are aromatic (order 4)
    for n, m, bond in mol.bonds():
        assert bond.order == 4


def test_pyrrole_thiele():
    # Test pyrrole-like systems
    mol = smiles('N1C=CC=C1')  # pyrrole in Kekule form
    assert mol.thiele()

    # Check that all bonds in the ring are aromatic
    for n, m, bond in mol.bonds():
        assert bond.order == 4


def test_furan_thiele():
    # Test furan-like systems
    mol = smiles('O1C=CC=C1')  # furan in Kekule form
    assert mol.thiele()

    # Check that all bonds in the ring are aromatic
    for n, m, bond in mol.bonds():
        assert bond.order == 4


def test_thiophene_thiele():
    # Test thiophene-like systems
    mol = smiles('S1C=CC=C1')  # thiophene in Kekule form
    assert mol.thiele()

    # Check that all bonds in the ring are aromatic
    for n, m, bond in mol.bonds():
        assert bond.order == 4


def test_condensed_rings():
    # Test condensed ring systems
    mol = smiles('C1=CC=C2C=CC=CC2=C1')  # naphthalene in Kekule form
    assert mol.thiele()

    # Check that all bonds in both rings are aromatic
    for n, m, bond in mol.bonds():
        assert bond.order == 4


def test_tautomer_fix():
    # Test tautomer fixing in condensed rings
    mol = smiles('N1C=CC2=NC=CC2=C1')  # before fix
    assert mol.thiele(fix_tautomers=True)

    for n, m, bond in mol.bonds():
        assert bond.order == 4

    assert mol.atom(1).implicit_hydrogens == 0
    assert mol.atom(5).implicit_hydrogens == 1

    mol = smiles('N1C=CC2=NC=CC2=C1')  # before fix
    assert mol.thiele(fix_tautomers=False)

    for n, m, bond in mol.bonds():
        assert bond.order == 4

    assert mol.atom(1).implicit_hydrogens == 1
    assert mol.atom(5).implicit_hydrogens == 0


def test_quinone_exclusion():
    # Test that quinone-like structures are not aromatized
    mol = smiles('O=C1C=CC(=O)C=C1')  # para-benzoquinone
    assert not mol.thiele()  # should return False
    assert mol == smiles('O=C1C=CC(=O)C=C1')


def test_invalid_rings():
    # Test that invalid rings are not aromatized
    mol = smiles('C1=CC=C1')  # cyclobutadiene - too small
    assert not mol.thiele()  # should return False

    mol = smiles('C1=CC=CC=CC=CC=C1')  # 10-membered ring - too large
    assert not mol.thiele()  # should return False


def test_charged_systems():
    # Test charged aromatic systems
    mol = smiles('[NH+]1C=CC=CC=1')  # pyridinium in Kekule form
    assert mol.thiele()

    # Check that charge is preserved and ring is aromatic
    assert mol.atom(1).charge == 1
    assert mol.atom(1).implicit_hydrogens == 1

    for n, m, bond in mol.bonds():
        assert bond.order == 4  # all bonds should be aromatic


def test_freak_rules():
    # Test special cases handled by freak rules
    mol = smiles('N1C=CN2C=CC=C12')  # special N-fused system in Kekule form
    assert mol.thiele()

    # Check that all bonds are aromatic
    for n, m, bond in mol.bonds():
        assert bond.order == 4


def test_tetracyclic_systems():
    # Test 4-membered ring in condensed systems
    mol = smiles('C1=CC2=C(C=CC=C2)C1')  # benzocyclobutene
    assert mol.thiele()
    assert mol == smiles('C1C=Cc2ccccc12')


def test_phosphorus_rings():
    # Test phosphorus-containing aromatic rings
    mol = smiles('P1C=CC=CC=1')  # phosphabenzene in Kekule form
    assert mol.thiele()

    # Check that all bonds in the ring are aromatic
    for n, m, bond in mol.bonds():
        assert bond.order == 4


def test_boron_rings():
    # Test boron-containing aromatic rings
    mol = smiles('B1C=CC=C1')  # borole in Kekule form
    assert mol.thiele()

    # Check that all bonds in the ring are aromatic
    for n, m, bond in mol.bonds():
        assert bond.order == 4


def test_seven_membered_rings():
    # Test 7-membered rings with different heteroatoms
    # Only boron-containing 7-membered rings should be aromatic
    mol = smiles('B1C=CC=CC=C1')  # 7-membered ring with boron
    assert mol.thiele()

    mol = smiles('N1C=CC=CC=C1')  # 7-membered ring with nitrogen
    assert not mol.thiele()


def test_ferrocene_like():
    # Test negatively charged carbon systems (ferrocene-like)
    mol = smiles('[CH-]1C=CC=C1')  # cyclopentadienyl anion
    assert mol.thiele()

    for n, m, bond in mol.bonds():
        assert bond.order == 4

    assert int(mol) == -1
    assert mol.atom(1).charge == -1


def test_multiple_components():
    # Test systems with multiple aromatic components
    mol = smiles('C1=CC=CC=C1.C1=CC=CC=C1')  # two benzene molecules in Kekule form
    assert mol.thiele()

    # Check that all bonds in both components are aromatic
    for n, m, bond in mol.bonds():
        assert bond.order == 4
    assert mol.connected_components_count == 2


def test_complex_fused_systems():
    # Test complex fused ring systems with multiple heteroatoms

    # Benzothiazole (simplified)
    mol = smiles('C1=CC=C2SC=NC2=C1')  # benzothiazole in Kekule form
    assert mol.thiele()
    for n, m, bond in mol.bonds():
        assert bond.order == 4

    # Thienopyridine (simplified)
    mol = smiles('C1=CC=NC2=CSC=C12')  # thienopyridine in Kekule form
    assert mol.thiele()
    for n, m, bond in mol.bonds():
        assert bond.order == 4


def test_complex_charged_systems():
    # Test charged aromatic systems

    # Basic charged systems that should work
    # Pyridinium
    mol = smiles('[NH+]1C=CC=CC=1')  # pyridinium in Kekule form
    assert mol.thiele()
    for n, m, bond in mol.bonds():
        assert bond.order == 4
    assert mol.atom(1).charge == 1
    assert mol.atom(1).implicit_hydrogens == 1

    # Complex charged systems that are not yet supported
    # N-methylpyridinium (currently not aromatized properly)
    mol = smiles('C[N+]1=CC=CC=C1')  # N-methylpyridinium in Kekule form
    assert mol.thiele()
    assert mol.atom(2).charge == 1


def test_complex_heterocycles():
    # Test heterocyclic systems

    # Basic heterocycles that should work
    # Benzimidazole
    mol = smiles('C1=CC=C2NC=NC2=C1')  # in Kekule form
    assert mol.thiele()
    for n, m, bond in mol.bonds():
        assert bond.order == 4

    # Quinoxaline
    mol = smiles('C1=CC=C2N=CC=NC2=C1')  # in Kekule form
    assert mol.thiele()
    for n, m, bond in mol.bonds():
        assert bond.order == 4

    # Benzimidazole fused to thiophene (works)
    mol = smiles('C1=CC2=C(C=C1)N=CN2C3=CC=CS3')  # in Kekule form
    assert mol.thiele()

    for n, m, bond in mol.bonds():
        # Check if both atoms of the bond are in the same ring
        assert bond.order == 4 or (bond.order == 1 and n in (9, 10) and m in (9, 10))


def test_complex_bridged_systems():
    # Complex bridged system with multiple heteroatoms (works)
    mol = smiles('C1=CC2=C(C=C1)N=C3C(=C2)C=CC4=C3N=CS4')  # in Kekule form
    assert mol.thiele()
    for n, m, bond in mol.bonds():
        assert bond.order == 4

    # Bridged system with N and S (works)
    mol = smiles('C1=CC2=C(C=C1)SC3=C(N=CC=C3)C=C2')  # in Kekule form
    assert mol.thiele()
    assert mol == smiles('S1c2ccccc2C=Cc2[n]cccc12')
