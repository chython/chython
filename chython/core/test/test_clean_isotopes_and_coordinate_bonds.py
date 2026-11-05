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
"""`clean_isotopes()` and `remove_coordinate_bonds()` -- the two standardization entry
points that need no rule table, and so are container methods rather than passes in `chemistry`.

THE DIVIDING LINE IS KNOWLEDGE, NOT MUTABILITY.  `standardize()` asks what a drawing meant and
answers out of 101 rows of chemical knowledge, so its body lives in `chython.chemistry` and reaches
the container through the registration hook.  These two ask nothing: one drops a field, the other
deletes every bond of one order.  `clean_stereo` is here for the same reason.

WHY THEY ARE TESTED TOGETHER.  They are the two halves of one question -- what does a mutation owe
the stereo state it leaves behind -- and the answers differ, which is the whole content of both
docstrings.  A bond is part of its anchor's FRAME, so deleting one is re-based by the journal's apply
and `remove_coordinate_bonds` owes nothing.  An isotope is in no frame, so a parity whose only
justification was the label survives the drop and starts lying; `clean_isotopes` therefore calls
`validate_stereo()`, the one implementation of that question this tree has.
"""
from chython.core import read_smiles


# ----------------------------------------------------------------------------------------------
# clean_isotopes
# ----------------------------------------------------------------------------------------------
def test_every_label_goes_and_the_answer_is_whether_one_was_there():
    mol = read_smiles('[13CH4]')
    assert mol.clean_isotopes() is True
    assert mol.smiles == 'C'
    assert mol.clean_isotopes() is False                 # idempotent, and says so


def test_deuterium_and_tritium_go_too():
    """A caller who wanted the heavy hydrogens kept has not asked for this.  Keeping an isotopic
    hydrogen as an explicit atom is `implicify_hydrogens`' business, where dropping the label would
    lose a fact nobody asked to lose; here losing it is the request."""
    mol = read_smiles('[2H]C([3H])O')
    assert mol.clean_isotopes()
    assert read_smiles(mol.smiles) == read_smiles('[H]C([H])O')
    assert not [a for a in mol.atoms() if a.isotope]


def test_a_molecule_with_no_isotope_is_a_pure_read():
    """No clone, no journal, no generation bump -- the arena a `copy()` shares is still shared."""
    mol = read_smiles('CCO')
    twin = mol.copy()
    assert mol.clean_isotopes() is False
    assert mol.shares_arena_with(twin)


def test_a_parity_the_label_justified_goes_with_the_label():
    """The one place a V3 mutation cannot leave stereo to `_apply`.

    The two methyls of `C[C@H](F)[13CH3]` differ only by the label, so the centre exists only while
    the label does -- and an isotope is in no FRAME, so the apply's re-basing does not see it.  Left
    alone, the molecule would keep writing `[C@H]` on an atom that has no configuration.
    """
    mol = read_smiles('C[C@H](F)[13CH3]')
    assert '@' in mol.smiles
    assert mol.clean_isotopes()
    assert '@' not in mol.smiles, mol.smiles
    assert mol.validate_stereo() == []                   # nothing left for a caller to clear


def test_a_centre_that_stands_without_the_label_keeps_its_parity():
    """The other side of the same test, and why the clear cannot simply be "drop every parity".  The
    label sits ON the centre here, and the four substituents are still four different things."""
    mol = read_smiles('[13C@@H](C)(F)Cl')
    assert mol.clean_isotopes()
    assert '@' in mol.smiles, mol.smiles
    assert read_smiles(mol.smiles) == read_smiles('[C@@H](C)(F)Cl')


# ----------------------------------------------------------------------------------------------
# remove_coordinate_bonds
# ----------------------------------------------------------------------------------------------
def test_the_coordination_sphere_comes_apart_and_the_metal_is_left_stranded():
    """The complement of what `standardize()` does: the metal-organic rules CREATE these bonds, and
    this is the call for a caller who wants them gone.  Stranding the metal is the point."""
    mol = read_smiles('[O+]#[C-]~[Fe](~[C-]#[O+])~[C-]#[O+]')
    assert mol.remove_coordinate_bonds() == 3
    assert read_smiles(mol.smiles) == read_smiles('[Fe].[C-]#[O+].[C-]#[O+].[C-]#[O+]')
    assert mol.remove_coordinate_bonds() == 0            # idempotent


def test_a_molecule_with_no_dative_bond_is_a_pure_read():
    mol = read_smiles('CCO')
    twin = mol.copy()
    assert mol.remove_coordinate_bonds() == 0
    assert mol.shares_arena_with(twin)


def test_a_hydrogen_held_only_by_coordination_keeps_its_bonds_by_default():
    """`keep_stranded_hydrogens=True` protects a hydrogen with NO covalent bond at all, because
    deleting its contacts leaves a disconnected `[H]` that names nothing."""
    kept = read_smiles('[H]~[Fe]~[H]')
    assert kept.remove_coordinate_bonds() == 0
    assert read_smiles(kept.smiles) == read_smiles('[H]~[Fe]~[H]')

    freed = read_smiles('[H]~[Fe]~[H]')
    assert freed.remove_coordinate_bonds(keep_stranded_hydrogens=False) == 2
    assert len(freed.split()) == 3                       # an iron and two loose hydrogens, as asked


def test_an_ordinary_hydrogen_bond_is_not_protected():
    """The guard is about being stranded, not about being a hydrogen: this donor keeps the covalent
    bond it came with, so its contact is deleted like any other."""
    mol = read_smiles('O([H])~O')
    assert mol.remove_coordinate_bonds() == 1
    assert read_smiles(mol.smiles) == read_smiles('[H]O.O')


def test_a_borane_keeps_its_bridging_hydrogen():
    mol = read_smiles('C[B]([H])([H])~[H]')
    assert mol.remove_coordinate_bonds() == 0
    assert len(mol.split()) == 1                         # nothing came off
