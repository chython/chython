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
"""ONE SLOT, ONE NAMESPACE: a stereo group lives at its unit's anchor, as a parity does.

`SEG_STEREO_GROUPS` is a byte per atom slot and a unit's byte sits at `u.anchor` -- so an axis's group
is stored on one atom, not on a bond, and an AND 1 stated on a centre and an AND 1 stated on an axis
are ONE collection.  A wire still spells an axis by its owners, which is what `set_stereo_group`
accepts and `stereo_groups()` answers with.
"""
import struct

from chython.core._core import read_smiles as smiles


STEREO_UNSPECIFIED, STEREO_ABS, STEREO_OR, STEREO_AND = 0, 1, 2, 3

# The segment table walk `test_cip_storage.py` uses, and the same three record offsets.
_SEG_CSR_EDGE = 2
_HALFEDGE_RECORD = 8
_HALFEDGE_FLAGS = 6      # uint16_t, and bits 5-12 are HE_FLAGS_LEGACY_GROUP
# `sg_pack(3, 1) << 5`: AND 1 is 0xc1, so the legacy field states a real collection.
_LEGACY_GROUP_BITS = 0x1820


def _forge_half_edge_group_bits(raw):
    """OR `_LEGACY_GROUP_BITS` into CSR half-edge records 0 and 1, in place, and return the buffer.

    In the `F/C=C/F` layout these CSR records are the two halves of the first bond.
    """
    edges_at = struct.unpack_from('<I', raw, 24 + 8 * _SEG_CSR_EDGE)[0]
    for k in (0, 1):
        at = edges_at + _HALFEDGE_RECORD * k + _HALFEDGE_FLAGS
        flags = struct.unpack_from('<H', raw, at)[0]
        struct.pack_into('<H', raw, at, flags | _LEGACY_GROUP_BITS)
    return raw


def test_one_id_space_holds_a_centre_and_an_axis_together():
    mol = smiles('C[C@H](O)/C=C/F')
    mol.set_stereo_group(2, STEREO_AND, 1)              # the centre
    mol.set_stereo_group((4, 5), STEREO_AND, 1)         # the axis, by its owners
    assert mol.stereo_groups() == {(STEREO_AND, 1): [2, (4, 5)]}
    assert mol.bond_stereo_groups() == {(STEREO_AND, 1): [(4, 5)]}


def test_an_axis_group_is_stored_at_the_anchor_atom():
    mol = smiles('F/C=C/F')
    mol.set_stereo_group((2, 3), STEREO_OR, 2)
    # the byte is on atom 2, the anchor -- the same slot a parity would use
    assert mol.stereo_group_of(2) == (STEREO_OR, 2)
    assert mol.stereo_group_anchor_of((2, 3)) == 2
    # and the pair spelling reads it back
    assert mol.stereo_group_of((2, 3)) == (STEREO_OR, 2)
    assert mol.stereo_group_of((3, 2)) == (STEREO_OR, 2)


def test_an_allene_group_is_stored_at_the_midpoint_and_reads_back_as_a_pair():
    mol = smiles('CC=C=CC')                             # penta-2,3-diene
    mol.set_stereo_group((2, 4), STEREO_AND, 3)
    assert mol.stereo_group_of(3) == (STEREO_AND, 3)    # the midpoint carries the byte
    assert mol.stereo_groups() == {(STEREO_AND, 3): [(2, 4)]}
    # the midpoint spelling is accepted too, and lands on the same slot
    mol.set_stereo_group(3, STEREO_OR, 1)
    assert mol.stereo_groups() == {(STEREO_OR, 1): [(2, 4)]}


def test_a_cumulene_group_takes_the_terminal_pair_that_is_not_a_bond():
    mol = smiles('CC=C=C=CC')                           # hexa-2,3,4-triene
    mol.set_stereo_group((2, 5), STEREO_AND, 1)
    assert 5 not in mol.neighbors_of(2)
    assert mol.stereo_groups() == {(STEREO_AND, 1): [(2, 5)]}


def test_a_chain_bond_spelling_reaches_the_axis_it_belongs_to():
    """FORGIVING IN: any honest spelling of the axis, including an interior chain bond."""
    mol = smiles('CC=C=C=CC')
    mol.set_bond_stereo_group(3, 4, STEREO_AND, 1)      # the middle bond of the chain
    assert mol.stereo_groups() == {(STEREO_AND, 1): [(2, 5)]}


def test_set_bond_stereo_group_is_the_pair_spelling_of_set_stereo_group():
    a = smiles('F/C=C/F')
    b = smiles('F/C=C/F')
    a.set_bond_stereo_group(2, 3, STEREO_OR, 1)
    b.set_stereo_group((2, 3), STEREO_OR, 1)
    assert a.stereo_groups() == b.stereo_groups() == {(STEREO_OR, 1): [(2, 3)]}


def test_a_pair_that_owns_no_unit_is_kept_as_a_label_and_logged():
    """DEGRADED, NOT DROPPED.  `(1, 2)` IS a bond in ethanol, so the setter accepts it; no unit owns
    the pair, so the seal keeps the collection on the lower atom and loses only the axis spelling.
    """
    mol = smiles('CCO')
    mol.set_bond_stereo_group(1, 2, STEREO_AND, 1)
    assert mol.stereo_groups() == {(STEREO_AND, 1): [1]}
    assert mol.bond_stereo_groups() == {}
    assert [r.rule for r in mol.log.by_stage('edit')] == ['edit:stereo-group-not-an-axis']


def test_an_atom_that_anchors_no_unit_still_carries_its_byte():
    """A BARE ATOM IS NEVER RESOLVED AWAY.  `set_stereo_group(n, ...)` writes the slot `n` resolves to
    whatever is there, so a file that states a collection on an atom the perception finds no unit at
    keeps it -- input is stored, not judged.

    Ethanol's oxygen is that atom: its two carbons are both kind-0 CANDIDATES, so only atom 3 anchors
    and owns nothing (`stereo_group_anchor_of(3) is None`, measured).
    """
    mol = smiles('CCO')
    mol.set_stereo_group(3, STEREO_AND, 1)
    assert mol.stereo_groups() == {(STEREO_AND, 1): [3]}
    assert mol.bond_stereo_groups() == {}


def test_a_group_survives_a_round_trip_through_bytes():
    mol = smiles('CC=C=CC')
    mol.set_stereo_group((2, 4), STEREO_AND, 3)
    from chython.core._core import MoleculeContainer
    back = MoleculeContainer.from_bytes(mol.to_bytes())
    assert back.stereo_groups() == {(STEREO_AND, 3): [(2, 4)]}


def test_deleting_the_anchor_takes_the_group_with_it():
    """A group dies with the SLOT it is stored at, exactly as a parity does -- not with the axis's
    stereogenicity.  Deleting owner atom 2 is deleting the anchor.
    """
    mol = smiles('F/C=C/F')
    mol.set_stereo_group((2, 3), STEREO_AND, 1)
    with mol.edit() as e:
        e.delete_atom(2)
    assert mol.stereo_groups() == {}


def test_deleting_a_substituent_leaves_the_label_on_the_anchor():
    """THE COMPLEMENT, and the one that proves where the byte lives.  Deleting an F destroys the unit
    -- `stereo_units()` is empty afterwards, measured -- while atom 2's slot survives, so the
    collection survives as a bare label on it.  A group is stated input; losing a unit does not
    withdraw a statement.
    """
    mol = smiles('F/C=C/F')
    mol.set_stereo_group((2, 3), STEREO_AND, 1)
    with mol.edit() as e:
        e.delete_atom(1)
    assert mol.stereo_groups() == {(STEREO_AND, 1): [2]}
    assert mol.bond_stereo_groups() == {}


def test_a_legacy_buffer_drops_its_half_edge_group_bytes_with_one_record():
    """Bits 5-12 of `halfedge_t.flags` are reserved and ignored.  A buffer that sets them loads,
    and the record says how many bonds carried them -- here one, since CSR records 0 and 1 are
    the two halves of one bond."""
    from chython.core._core import MoleculeContainer
    mol = smiles('F/C=C/F')
    raw = bytearray(mol.to_bytes())
    forged = _forge_half_edge_group_bits(raw)
    back = MoleculeContainer.from_bytes(bytes(forged))
    assert back.stereo_groups() == {}
    record, = back.log.by_stage('read')
    assert record.rule == 'container:bond-group-bits-dropped'
    assert '1 bond' in record


def test_an_edit_that_re_anchors_an_axis_moves_the_group_byte_with_it():
    """A CARRIED BYTE FOLLOWS ITS UNIT.  The apply copies group bytes by SLOT and then derives the
    fresh unit table, which re-anchors an atropisomer axis to whichever pivot is free (ruling F45).
    Left where it was, the byte reads back as a bare label on an atom and `bond_stereo_groups()` no
    longer holds the axis -- silently, since the byte is still stored.

    Breaking the ring double bond at pivot 1 frees that pivot, so the axis moves from atom 18 to atom 1
    and the collection has to move with it.  Nothing is lost, so nothing is logged.
    """
    mol = smiles('C12=CC=CC=CC=C1C.CC1=CC=CC=CC=C12')
    mol.set_stereo_group((1, 18), STEREO_AND, 1)
    assert mol.stereo_group_anchor_of((1, 18)) == 18
    with mol.edit() as e:
        e.set_order(1, 2, 1)
    assert mol.stereo_group_anchor_of((1, 18)) == 1, 'the axis re-anchored, which is the premise'
    assert mol.stereo_groups() == {(STEREO_AND, 1): [(1, 18)]}
    assert mol.bond_stereo_groups() == {(STEREO_AND, 1): [(1, 18)]}
    assert mol.stereo_group_of((1, 18)) == (STEREO_AND, 1)
    assert list(mol.log) == []


def test_a_re_anchored_axis_whose_new_anchor_is_taken_keeps_its_label_and_logs():
    """THE COLLISION ARM.  One byte per slot, so a unit that re-anchors onto an atom already stating a
    collection has nowhere to go: the byte stays where it is, as a label on that atom, and the record
    says the axis spelling is what the edit lost.

    Both collections here are stated before the edit -- AND 1 on the axis, anchored at pivot 18, and
    OR 1 on the ring double bond anchored at pivot 1 -- and breaking that double bond re-anchors the
    axis onto the slot the OR 1 byte occupies.
    """
    mol = smiles('C12=CC=CC=CC=C1C.CC1=CC=CC=CC=C12')
    mol.set_stereo_group((1, 18), STEREO_AND, 1)
    mol.set_stereo_group(1, STEREO_OR, 1)
    assert mol.stereo_groups() == {(STEREO_OR, 1): [(1, 2)], (STEREO_AND, 1): [(1, 18)]}
    with mol.edit() as e:
        e.set_order(1, 2, 1)
    assert mol.stereo_groups() == {(STEREO_OR, 1): [(1, 18)], (STEREO_AND, 1): [18]}
    assert mol.bond_stereo_groups() == {(STEREO_OR, 1): [(1, 18)]}
    record, = mol.log.by_stage('edit')
    assert record.rule == 'edit:stereo-group-anchor-taken'
    assert 'atoms 1 and 18' in record and 'atom 18' in record


def test_an_edit_that_touches_no_stereo_unit_moves_no_group_byte():
    """THE MIGRATION IS NOT A SWEEP.  A byte at a unit's anchor stays when the unit's anchor does
    not change: only a byte whose unit re-anchored is migrated.

    `F/C=C/F` with AND 1 on the cis/trans axis (2, 3) is the minimal case: the axis's owners (2, 3)
    survive the element edit and the anchor stays at atom 2, so the migration pass runs, resolves
    (2, 3) back to the same anchor, and declines to move the byte.
    """
    mol = smiles('F/C=C/F')
    mol.set_stereo_group((2, 3), STEREO_AND, 1)
    with mol.edit() as e:
        e.set_element(1, 'Cl')
    assert mol.stereo_groups() == {(STEREO_AND, 1): [(2, 3)]}
    assert list(mol.log.by_stage('edit')) == []
