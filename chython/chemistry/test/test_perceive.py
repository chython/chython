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
"""`perceive_bonds`: connectivity from a stored model, on public experimental geometries.

Every geometry below is built from published bond lengths and angles for the named compound, so a
threshold that bonds a hydrogen-bonded pair or misses a peroxide fails here rather than on a file.
"""
import pytest

from .._perceive import perceive_bonds
from ...core import INFO, LOST, MoleculeContainer


#: Water: O-H 0.958 A, H-O-H 104.5 deg -- so H...H is 1.514 A, the shortest nonbonded pair a small
#: molecule offers and the one a loose threshold turns into a bond.
WATER = (('O', .0, .0, .0), ('H', .757, .586, .0), ('H', -.757, .586, .0))

#: Hydrogen peroxide: O-O 1.475 A, O-H 0.950 A, O-O-H 94.8 deg, dihedral 111.5 deg.  The O-O bond is
#: the longest first-row single bond a tight threshold misses.
PEROXIDE = (('O', .0, .0, .0), ('O', 1.475, .0, .0),
            ('H', -.0794, .9467, .0), ('H', 1.5545, -.3469, .8808))

#: Methane: C-H 1.087 A, tetrahedral -- H...H 1.775 A.
METHANE = (('C', .0, .0, .0), ('H', .6276, .6276, .6276), ('H', -.6276, -.6276, .6276),
           ('H', -.6276, .6276, -.6276), ('H', .6276, -.6276, -.6276))

#: The water dimer: O...O 2.95 A with a nearly linear O-H...O, so H...O is 1.99 A.  TWO molecules,
#: and a threshold that reads a hydrogen bond as a bond answers one.
WATER_DIMER = (('O', .0, .0, .0), ('H', .958, .0, .0), ('H', -.2397, .9276, .0),
               ('O', 2.95, .0, .0), ('H', 3.5033, .7833, .0), ('H', 3.5033, -.7833, .0))

#: Cyclobutadiene: a D2h rectangle, C=C 1.344 A and C-C 1.441 A, C-H 1.083 A along each diagonal.
#: Its transannular C...C is 1.970 A -- the TIGHTEST nonbonded pair a neutral organic molecule offers,
#: and the case a threshold loose enough to reach F2's long F-F bond turns into a bicyclobutane.
CYCLOBUTADIENE = (('C', .672, .7205, .0), ('C', -.672, .7205, .0),
                  ('C', -.672, -.7205, .0), ('C', .672, -.7205, .0),
                  ('H', 1.4107, 1.5124, .0), ('H', -1.4107, 1.5124, .0),
                  ('H', -1.4107, -1.5124, .0), ('H', 1.4107, -1.5124, .0))

#: Benzene: C-C 1.397 A, C-H 1.084 A, so the ring radius is 1.397 and the hydrogens sit at 2.481.
#: The meta C...C pair is 2.420 A, the second-shortest nonbonded pair in the file.
BENZENE = (('C', 1.397, .0, .0), ('C', .6985, 1.2098, .0), ('C', -.6985, 1.2098, .0),
           ('C', -1.397, .0, .0), ('C', -.6985, -1.2098, .0), ('C', .6985, -1.2098, .0),
           ('H', 2.481, .0, .0), ('H', 1.2405, 2.1486, .0), ('H', -1.2405, 2.1486, .0),
           ('H', -2.481, .0, .0), ('H', -1.2405, -2.1486, .0), ('H', 1.2405, -2.1486, .0))


def _placed(geometry, *more):
    """A molecule holding `geometry`'s atoms and coordinates and NO bond, plus a model per extra."""
    mol = MoleculeContainer()
    for element, x, y, z in geometry:
        n = mol.add_atom(element, implicit_h=0)
        mol.set_xyz(n, x, y, z)
    numbers = mol.atom_numbers            # read before the scope: an open edit answers no query
    for extra in more:
        with mol.edit():
            model = mol.add_conformer()
            for n, (_, x, y, z) in zip(numbers, extra):
                mol.set_xyz(n, x, y, z, model=model)
    return mol


def _pairs(mol):
    return {frozenset((b.n, b.m)) for b in mol.bonds()}


def test_water_gets_its_two_bonds():
    mol = _placed(WATER)
    assert perceive_bonds(mol)
    assert _pairs(mol) == {frozenset((1, 2)), frozenset((1, 3))}


def test_the_shortest_nonbonded_pair_in_water_is_not_a_bond():
    """H...H at 1.514 A.  Two hydrogens on one oxygen are 1,3 and no threshold may join them."""
    mol = _placed(WATER)
    perceive_bonds(mol)
    assert frozenset((2, 3)) not in _pairs(mol)


def test_the_peroxide_oxygens_are_bonded():
    """O-O at 1.475 A, the case that rules out the calculated radii: their sum is 0.96 A."""
    mol = _placed(PEROXIDE)
    perceive_bonds(mol)
    assert frozenset((1, 2)) in _pairs(mol)
    assert len(_pairs(mol)) == 3


def test_methane_is_one_carbon_and_four_bonds():
    mol = _placed(METHANE)
    perceive_bonds(mol)
    assert len(_pairs(mol)) == 4


def test_a_hydrogen_bond_is_not_a_bond():
    """The water dimer stays TWO molecules: H...O at 1.99 A is a contact and not a bond."""
    mol = _placed(WATER_DIMER)
    perceive_bonds(mol)
    assert len(_pairs(mol)) == 4
    assert len(mol.connected_components) == 2


def test_benzene_perceives_twelve_bonds_and_no_cross_ring_pair():
    mol = _placed(BENZENE)
    perceive_bonds(mol)
    assert len(_pairs(mol)) == 12
    assert frozenset((1, 3)) not in _pairs(mol)     # meta C...C, 2.420 A


def test_a_four_membered_ring_keeps_its_four_bonds():
    """Cyclobutadiene is a ring of four and not a bicyclobutane: the 1.970 A diagonals are contacts."""
    mol = _placed(CYCLOBUTADIENE)
    perceive_bonds(mol)
    assert len(_pairs(mol)) == 8                    # four ring bonds and four C-H
    assert frozenset((1, 3)) not in _pairs(mol)
    assert frozenset((2, 4)) not in _pairs(mol)


def test_every_perceived_bond_is_single():
    """Order is not this pass's question: `saturate()` raises what the valence rules force."""
    mol = _placed(PEROXIDE)
    perceive_bonds(mol)
    assert {int(b) for b in mol.bonds()} == {1}


def test_a_stated_bond_is_left_alone():
    """A double bond already in the molecule keeps its order, and is not added twice."""
    mol = _placed(WATER)
    with mol.edit():
        mol.add_bond(1, 2, 2)
    assert perceive_bonds(mol)
    orders = {frozenset((b.n, b.m)): int(b) for b in mol.bonds()}
    assert orders == {frozenset((1, 2)): 2, frozenset((1, 3)): 1}


def test_a_molecule_whose_bonds_are_all_there_is_unchanged():
    mol = _placed(WATER)
    perceive_bonds(mol)
    before = bytes(mol)
    assert not perceive_bonds(mol)
    assert bytes(mol) == before


def test_two_atoms_too_far_apart_get_nothing():
    mol = _placed((('C', .0, .0, .0), ('C', 4.0, .0, .0)))
    assert not perceive_bonds(mol)
    assert not _pairs(mol)


def test_the_multiplier_is_the_knob():
    """One C...C pair at 2.2 A: nonbonded at the default, bonded once the multiplier reaches it."""
    mol = _placed((('C', .0, .0, .0), ('C', 2.2, .0, .0)))
    assert not perceive_bonds(mol)
    assert perceive_bonds(mol, radius_multiplier=1.5)
    assert _pairs(mol) == {frozenset((1, 2))}


def test_the_model_is_chosen_and_defaults_to_the_first():
    """Model 0 is water and model 1 pulls one hydrogen 3 A away, so the two answers differ."""
    stretched = (('O', .0, .0, .0), ('H', .757, .586, .0), ('H', -3.0, .586, .0))
    mol = _placed(WATER, stretched)
    assert perceive_bonds(mol, model=1)
    assert _pairs(mol) == {frozenset((1, 2))}

    other = _placed(WATER, stretched)
    perceive_bonds(other)
    assert len(_pairs(other)) == 2


def test_a_model_that_is_not_there_raises():
    """The container's own answer, naming the model; nothing is perceived against an invented one."""
    mol = _placed(WATER)
    with pytest.raises(IndexError):
        perceive_bonds(mol, model=3)


def test_a_molecule_with_no_geometry_raises():
    from ...core import read_smiles

    mol = read_smiles('CCO')
    assert not mol.has_3d
    with pytest.raises(IndexError):
        perceive_bonds(mol)


def test_the_r_marker_gets_no_bond_and_the_log_says_so():
    """Element 0 has no covalent radius, so nothing is bonded to it and the shortfall is recorded."""
    mol = _placed((('C', .0, .0, .0), ('O', 1.43, .0, .0), ('R', -1.5, .0, .0)))
    perceive_bonds(mol)
    assert _pairs(mol) == {frozenset((1, 2))}
    records = [r for r in mol.log if r.rule == 'perceive:no-radius']
    assert len(records) == 1
    assert records[0].severity == LOST
    assert '1' in records[0].message


def test_what_it_did_is_recorded_under_its_own_stage():
    mol = _placed(WATER)
    perceive_bonds(mol)
    records = [r for r in mol.log if r.rule == 'perceive:bonds']
    assert len(records) == 1
    assert records[0].severity == INFO
    assert records[0].stage == 'perceive_bonds'
    assert '2' in records[0].message


def test_it_takes_no_log_argument():
    """A pass writes to `molecule.log` and takes no destination; see `chython.core.recording`."""
    from inspect import signature

    assert 'log' not in signature(perceive_bonds).parameters
