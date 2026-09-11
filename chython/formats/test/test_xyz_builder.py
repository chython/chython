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
"""`chython.formats.xyz.build_molecule`: frames to a container of atoms, coordinates and no bond.

The builder states what the format states.  What it does not do is perceive: bonds come from
`perceive_bonds()` and orders from `saturate()`, both called by the caller, and the tests here assert
the molecule arrives without them.
"""
from ...core import read_smiles
from ..xyz import build_molecule, xyz


#: Water, as an XYZ record: O-H 0.958 A, H-O-H 104.5 deg.
WATER = """3
water
O    0.000000    0.000000    0.000000
H    0.757000    0.586000    0.000000
H   -0.757000    0.586000    0.000000
"""

#: Two frames of one molecule -- the shape a trajectory or a relaxation writes.
TRAJECTORY = WATER + """3
water, stretched
O    0.000000    0.000000    0.000000
H    0.857000    0.586000    0.000000
H   -0.857000    0.586000    0.000000
"""


def test_the_atoms_are_the_frames_atoms_in_file_order():
    mol = build_molecule(xyz(WATER)[0])
    assert [a.atomic_symbol for a in mol.atoms()] == ['O', 'H', 'H']


def test_the_coordinates_land_as_the_first_model():
    mol = build_molecule(xyz(WATER)[0])
    assert mol.has_3d
    assert len(mol.conformers) == 1
    assert mol.conformer(0).xyz_of(2) == (.757, .586, .0)


def test_nothing_is_bonded():
    """The format states no bond, so neither does the molecule; `perceive_bonds()` is the next call."""
    mol = build_molecule(xyz(WATER)[0])
    assert not list(mol.bonds())


def test_every_atom_states_zero_implicit_hydrogens():
    """An XYZ record states EVERY atom, hydrogens included, so an absent hydrogen is absent."""
    mol = build_molecule(xyz(WATER)[0])
    assert [mol.implicit_h_of(n) for n in mol.atom_numbers] == [0, 0, 0]


def test_a_deuterium_keeps_the_isotope_the_reader_read():
    mol = build_molecule(xyz('2\nHD\nH 0. 0. 0.\nD 0. 0. 0.74\n')[0])
    assert [(a.atomic_symbol, a.isotope) for a in mol.atoms()] == [('H', 0), ('H', 2)]


def test_every_frame_becomes_a_model():
    """A list of frames is one molecule with one model per frame, in the order they were read."""
    mol = build_molecule(xyz(TRAJECTORY))
    assert len(mol.conformers) == 2
    assert [c.ext_index for c in mol.conformers] == [0, 1]
    assert mol.conformer(1).xyz_of(2) == (.857, .586, .0)


def test_a_single_frame_and_a_list_of_one_agree():
    assert bytes(build_molecule(xyz(WATER)[0])) == bytes(build_molecule(xyz(WATER)))


def test_a_frame_that_disagrees_with_the_first_is_logged_and_skipped():
    """The frames of one file may hold different molecules; the first is the one that is built."""
    log = []
    mol = build_molecule(xyz(WATER + '1\nlone argon\nAr 0. 0. 0.\n'), log=log)
    assert len(mol.conformers) == 1
    assert any(r.rule == 'xyz:frame-atom-count' for r in log)


def test_an_unreadable_symbol_becomes_a_marker_and_says_so():
    """Element 0 keeps the atom and its coordinate; dropping the row would lose the geometry."""
    log = []
    mol = build_molecule(xyz('2\nunknown\nQq 0. 0. 0.\nH 0. 0. 1.\n')[0], log=log)
    assert [a.atomic_symbol for a in mol.atoms()] == ['R', 'H']
    assert any(r.rule == 'xyz:symbol-not-an-element' for r in log)


def test_the_readers_own_findings_reach_the_molecule():
    """A frame's log is about the record this molecule now IS, so `mol.log` is where it ends up."""
    mol = build_molecule(xyz('1\ncase\nCL 0. 0. 0.\n')[0])
    assert [a.atomic_symbol for a in mol.atoms()] == ['Cl']
    assert any(r.rule == 'xyz:symbol-case-corrected' for r in mol.log)


def test_the_title_is_not_read_as_a_charge():
    """`charge=-1` in a comment names no atom, and nothing here places a charge on a guess."""
    log = []
    mol = build_molecule(xyz('1\ncharge=-1\nCl 0. 0. 0.\n')[0], log=log)
    assert mol.charge_of(1) == 0
    assert any(r.rule == 'xyz:title-charge-not-applied' for r in log)


def test_no_frames_is_an_empty_molecule():
    mol = build_molecule([])
    assert not len(mol)
    assert not mol.has_3d


def test_the_whole_pipeline_gives_the_molecule_the_file_meant():
    """Read, perceive, saturate, implicify: the four explicit calls, on an acetonitrile geometry.

    C-C 1.458 A, C#N 1.157 A, C-H 1.087 A -- so the triple bond is forced by the hydrogen counts and
    not by the distance, which perception never reads as an order.  The hydrogens are folded in last
    because an XYZ record states them as atoms and a chemist's acetonitrile does not.
    """
    from ...chemistry import implicify_hydrogens, perceive_bonds, saturate

    text = """6
acetonitrile
C    0.000000    0.000000    0.000000
C    0.000000    0.000000    1.458000
N    0.000000    0.000000    2.615000
H    1.024000    0.000000   -0.362000
H   -0.512000   -0.887000   -0.362000
H   -0.512000    0.887000   -0.362000
"""
    mol = build_molecule(xyz(text)[0])
    assert perceive_bonds(mol)
    assert saturate(mol)
    assert {int(b) for b in mol.bonds()} == {1, 3}
    assert implicify_hydrogens(mol)
    # As a STRUCTURE and never as a SMILES string: the writer emits arena order, which is file order.
    assert mol == read_smiles('CC#N')
    assert mol.has_3d, 'folding the hydrogens in dropped the model the file stated'
