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
"""OWNERS AGAINST ANCHOR: the two keys a stereo group is spelled with.

Memory keys a group at the unit's ANCHOR slot, the way `SEG_PARITY` keys a parity.  A wire keys it on
the unit's OWNERS -- the atoms the configuration is named on -- because a pach entry, a V3000
collection and a CXSMILES member list are all written before the reader has a unit table to consult.
`stereo_unit_owners` is the one derivation from the first key to the second, and
`structure_owner_pair_anchor` is its inverse; these tests pin both directions and the one place they
are not inverse, which is a bare atom number naming two different units.
"""
import pytest
from chython.core._core import read_smiles as smiles


def test_a_tetrahedral_centre_owns_itself():
    mol = smiles('C[C@H](N)C(=O)O')                     # alanine
    unit = mol.chiral_atoms()[2]
    assert unit['owners'] == 2
    assert unit['anchor'] == 2
    assert mol.stereo_group_anchor_of(2) == 2


def test_a_cis_trans_axis_owns_its_two_terminals():
    mol = smiles('F/C=C/F')                             # 1,2-difluoroethene
    (pair, unit), = mol.chiral_bonds().items()
    assert pair == (2, 3)
    assert unit['owners'] == (2, 3)
    assert unit['anchor'] == 2
    # either spelling of the pair resolves to the one anchor
    assert mol.stereo_group_anchor_of((2, 3)) == 2
    assert mol.stereo_group_anchor_of((3, 2)) == 2


def test_a_cumulene_owns_two_atoms_that_are_not_bonded():
    mol = smiles('CC=C=C=CC')                           # hexa-2,3,4-triene
    (pair, unit), = mol.chiral_bonds().items()
    assert pair == (2, 5)
    assert unit['owners'] == (2, 5)
    assert 5 not in mol.neighbors_of(2)
    assert mol.stereo_group_anchor_of((2, 5)) == 2


def test_an_allene_owns_its_terminals_and_anchors_between_them():
    mol = smiles('CC=C=CC')                             # penta-2,3-diene
    unit = mol.chiral_atoms()[3]
    assert unit['owners'] == (2, 4)
    assert unit['anchor'] == 3
    # THE ANCHOR IS NEITHER OWNER: the pair spelling resolves through the table, not by arithmetic
    assert mol.stereo_group_anchor_of((2, 4)) == 3
    # and the midpoint spelling resolves to itself
    assert mol.stereo_group_anchor_of(3) == 3


def test_an_atropisomer_owns_both_pivots_however_it_is_anchored():
    mol = smiles('C12=CC=CC=CC=C1C.CC1=CC=CC=CC=C12')
    axes = [(p, u) for p, u in mol.chiral_bonds().items() if u['kind'] == 3]
    assert len(axes) == 1
    pair, unit = axes[0]
    assert pair == (1, 18) and unit['owners'] == (1, 18)
    assert unit['anchor'] == 18
    assert mol.stereo_group_anchor_of((1, 18)) == 18


def test_a_bare_atom_number_names_the_unit_anchored_there():
    """AN ANCHOR OUTRANKS AN OWNERSHIP, which is why an axis has to be spelled as a pair.

    Atom 1 anchors a cis/trans unit AND is an owner of the atropisomer axis anchored at 18.  A bare
    `1` therefore names the cis/trans unit -- the one whose byte lives at slot 1 -- and the axis is
    reachable only as `(1, 18)`.
    """
    mol = smiles('C12=CC=CC=CC=C1C.CC1=CC=CC=CC=C12')
    assert mol.stereo_group_anchor_of(1) == 1
    assert mol.stereo_group_anchor_of((1, 18)) == 18


def test_an_owner_that_anchors_nothing_resolves_to_its_unit():
    mol = smiles('F/C=C/F')
    assert mol.stereo_group_anchor_of(3) == 2           # the higher terminal owns, anchors nothing


def test_an_element_that_owns_no_unit_answers_none():
    mol = smiles('C(=O)(N)N')                           # urea: no candidate anywhere
    assert mol.stereo_units() == []
    assert mol.stereo_group_anchor_of(1) is None
    assert mol.stereo_group_anchor_of((1, 2)) is None


def test_a_candidate_that_is_not_stereogenic_still_answers_its_own_anchor():
    """`set_stereo_group` stores a byte at any atom's slot without consulting stereogenicity, so this
    read may not be narrower than that write: ethanol's two carbons are kind-0 candidates, and a byte
    written on either would live at that atom's own slot.
    """
    mol = smiles('CCO')
    assert [(u['kind'], u['anchor'], u['stereogenic']) for u in mol.stereo_units()] \
        == [(0, 1, False), (0, 2, False)]
    assert mol.stereo_group_anchor_of(1) == 1
    assert mol.stereo_group_anchor_of(2) == 2
    assert mol.stereo_group_anchor_of(3) is None        # the hydroxyl O is no kind's candidate


def test_stereo_group_anchor_of_refuses_a_spelling_it_cannot_read():
    mol = smiles('CCO')
    for bad in ((1,), (1, 2, 3), 'C', 1.0):
        try:
            mol.stereo_group_anchor_of(bad)
        except (TypeError, ValueError):
            continue
        raise AssertionError('accepted %r' % (bad,))


def test_an_unknown_atom_number_raises():
    mol = smiles('CCO')
    for bad in (99, (1, 99)):
        try:
            mol.stereo_group_anchor_of(bad)
        except KeyError:
            continue
        raise AssertionError('accepted %r' % (bad,))


def test_a_pair_naming_one_atom_twice_is_not_a_spelling():
    """`set_bond_stereo_group` raises `KeyError` for `n == m`, so a read that quietly answered the
    bare-int question instead would resolve a pair the write cannot store.
    """
    mol = smiles('F/C=C/F')
    assert mol.stereo_group_anchor_of(3) == 2
    with pytest.raises(ValueError):
        mol.stereo_group_anchor_of((3, 3))
