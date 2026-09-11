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
"""CIP descriptor STORAGE.  Nothing here assigns one -- the arena records what an input stated.

Three properties are worth more than the rest of this file put together, and each has its own section:

* a descriptor SURVIVES remap, copy, a pack round trip, a property-only edit, and `kekule`/`thiele`;
* a descriptor is DROPPED, and logged, by an edit that changes which atoms exist, which are bonded, or
  what a bond's order is -- because a descriptor is a ranking and a ranking reads the whole molecule;
* a descriptor is part of the BYTES and not part of the CANONICAL FORM, so it round-trips through
  `to_bytes` while two molecules that differ only in a descriptor stay equal and stay one hash.

There is no assignment algorithm, so no test here asserts that a descriptor is CORRECT.  A wrong 'R'
goes in and comes back out; that is the arena's posture everywhere and it is deliberate.
"""

import struct

import pytest

from chython.core import MoleculeContainer


# every value the two domains accept, which is also the whole encodable range
ATOM_DESCRIPTORS = ('R', 'S', 'r', 's', 'M', 'P', 'm', 'p')
BOND_DESCRIPTORS = ('E', 'Z', 'M', 'P')

# byte layout, the same constants `test_pack.py` reads the buffer with
_SEG_ATOMS = 0
_SEG_CSR_EDGE = 2
_ATOM_RECORD = 24
_ATOM_RESERVED = 20      # uint32_t, and the low nibble is the atom's CIP code
_HALFEDGE_RECORD = 8
_HALFEDGE_FLAGS = 6      # uint16_t, and bits 2-4 are the bond's CIP code


def _propene():
    """CH3-CH=CH-F: one double bond to carry a bond descriptor, one carbon to carry an atom one.

    Hydrogen counts are STATED rather than omitted, because `thiele` needs them and an atom whose count
    nobody gave cannot be classified as aromatic.  That is not incidental to this file: the
    `kekule`/`thiele` survival tests are the reason the fixture is built this way.
    """
    mol = MoleculeContainer()
    with mol.edit() as e:
        c1 = e.add_atom('C', implicit_h=3)
        c2 = e.add_atom('C', implicit_h=1)
        c3 = e.add_atom('C', implicit_h=1)
        f = e.add_atom('F', implicit_h=0)
        e.add_bond(c1, c2)
        e.add_bond(c2, c3, 2)
        e.add_bond(c3, f)
    return mol, c1, c2, c3, f


def _benzene_with_a_methyl():
    """Toluene, spelled Kekule, so `thiele` has something to change and `kekule` has it back."""
    mol = MoleculeContainer()
    with mol.edit() as e:
        ring = [e.add_atom('C', implicit_h=1) for _ in range(6)]
        for (a, b), order in zip(zip(ring, ring[1:] + ring[:1]), (2, 1, 2, 1, 2, 1)):
            e.add_bond(a, b, order)
        methyl = e.add_atom('C', implicit_h=3)
        e.add_bond(ring[0], methyl)
    return mol, ring, methyl


def _round_trip(mol):
    return MoleculeContainer.from_bytes(mol.to_bytes())


# ------------------------------------------------------------------------------------------------
# the domains: every value, and every neighbouring value that is NOT one
# ------------------------------------------------------------------------------------------------

@pytest.mark.parametrize('descriptor', ATOM_DESCRIPTORS)
def test_every_atom_descriptor_goes_in_and_comes_back_as_itself(descriptor):
    mol, c1, c2, c3, f = _propene()
    mol.set_atom_cip(c2, descriptor)
    assert mol.atom_cip_of(c2) == descriptor
    assert mol.atom_cips() == {c2: descriptor}


@pytest.mark.parametrize('descriptor', BOND_DESCRIPTORS)
def test_every_bond_descriptor_goes_in_and_comes_back_as_itself(descriptor):
    mol, c1, c2, c3, f = _propene()
    mol.set_bond_cip(c2, c3, descriptor)
    assert mol.bond_cip_of(c2, c3) == descriptor
    assert mol.bond_cips() == {(c2, c3): descriptor}


def test_an_atom_with_no_descriptor_answers_none():
    """`None` is the answer and not a KeyError: "no descriptor" is a value this storage holds."""
    mol, c1, c2, c3, f = _propene()
    assert mol.atom_cip_of(c2) is None
    assert mol.bond_cip_of(c2, c3) is None
    assert mol.atom_cips() == {}
    assert mol.bond_cips() == {}


def test_none_clears_a_descriptor_and_leaves_no_trace():
    mol, c1, c2, c3, f = _propene()
    before = mol.to_bytes()
    mol.set_atom_cip(c2, 'R')
    mol.set_bond_cip(c2, c3, 'E')
    mol.set_atom_cip(c2, None)
    mol.set_bond_cip(c2, c3, None)
    assert mol.atom_cip_of(c2) is None
    assert mol.bond_cip_of(c2, c3) is None
    # BYTE-FOR-BYTE back to where it started.  Clearing that left a bit set somewhere would pass the
    # two reads above and still make a cleared molecule a different key from one never labelled.
    assert mol.to_bytes() == before


def test_lowercase_is_a_different_descriptor_and_not_a_spelling_of_the_uppercase():
    """r/s are the pseudo-asymmetric descriptors: a different determination about a different centre.

    Anything that upper-cases on the way in has lost information, so the two must store as two codes
    and 'e' must be a refusal rather than a courtesy read of 'E'.
    """
    mol, c1, c2, c3, f = _propene()
    mol.set_atom_cip(c2, 'R')
    upper = mol.to_bytes()
    mol.set_atom_cip(c2, 'r')
    assert mol.atom_cip_of(c2) == 'r'
    assert mol.to_bytes() != upper, "'r' and 'R' stored as the same code"

    with pytest.raises(ValueError, match='case is significant'):
        mol.set_atom_cip(c2, 'e')
    with pytest.raises(ValueError, match='case is significant'):
        mol.set_bond_cip(c2, c3, 'z')
    # and the refusal did not half-apply: the molecule still holds what it held
    assert mol.atom_cip_of(c2) == 'r'


def test_a_descriptor_from_the_other_domain_is_refused_in_both_directions():
    """'E' on an atom and 'R' on a bond are mistakes, and M/P being valid on both is why.

    One merged table would encode an atom's 'M' and a bond's 'M' as the same code, and then a caller
    who sent a bond descriptor to an atom would get silence instead of an error.  Two tables make it a
    ValueError, and the message says which domain was checked.
    """
    mol, c1, c2, c3, f = _propene()
    with pytest.raises(ValueError, match='for an atom'):
        mol.set_atom_cip(c2, 'E')
    with pytest.raises(ValueError, match='for a bond'):
        mol.set_bond_cip(c2, c3, 'R')
    # M and P are the two letters both domains accept, and they are accepted on both
    mol.set_atom_cip(c2, 'M')
    mol.set_bond_cip(c2, c3, 'M')
    assert mol.atom_cip_of(c2) == 'M'
    assert mol.bond_cip_of(c2, c3) == 'M'


def test_a_descriptor_that_is_not_a_str_is_refused():
    mol, c1, c2, c3, f = _propene()
    with pytest.raises(TypeError, match='must be a str or None'):
        mol.set_atom_cip(c2, b'R')
    with pytest.raises(TypeError, match='must be a str or None'):
        mol.set_atom_cip(c2, 1)
    with pytest.raises(TypeError, match='must be a str or None'):
        mol.set_bond_cip(c2, c3, 2)


def test_a_bond_descriptor_needs_two_distinct_atoms():
    mol, c1, c2, c3, f = _propene()
    with pytest.raises(ValueError, match='two distinct atoms'):
        mol.set_bond_cip(c2, c2, 'E')


def test_a_descriptor_on_an_atom_that_does_not_exist_is_refused():
    mol, c1, c2, c3, f = _propene()
    with pytest.raises(KeyError):
        mol.set_atom_cip(9999, 'R')
    with pytest.raises(KeyError):
        mol.set_bond_cip(c2, 9999, 'E')


def test_a_descriptor_on_a_pair_that_is_not_bonded_is_refused_at_read():
    """The setter takes two live atoms; it is the READ that knows whether they are bonded.

    Worth pinning because the two are asymmetric on purpose: the setter runs while the journal is open
    and the bond may not exist yet, so refusing there would refuse the reader's own ordering.
    """
    mol, c1, c2, c3, f = _propene()
    with pytest.raises(KeyError):
        mol.bond_cip_of(c1, f)


def test_a_descriptor_named_before_its_atom_is_added_fails_loudly():
    """A CIP op naming an atom the journal has not added yet must raise, not land on a stray slot.

    `set_atom_cip` resolves through `_work_index` at replay, and the id it is given here is the one
    `add_atom` is ABOUT to hand out -- the single value most likely to be guessed by a caller building
    a molecule and its descriptors in one pass, which is exactly what a parser does.
    """
    mol = MoleculeContainer()
    with pytest.raises(KeyError):
        with mol.edit() as e:
            first = e.add_atom('C')
            e.set_atom_cip(first + 1, 'R')     # the id the NEXT add_atom would return
            e.add_atom('C')
    # and the refused scope left nothing behind
    assert mol.atom_count == 0


# ------------------------------------------------------------------------------------------------
# a bond descriptor is one statement about one bond, from either end
# ------------------------------------------------------------------------------------------------

def test_a_bond_descriptor_reads_the_same_from_either_end():
    """Not directional, unlike a wedge -- which is why both half-edges are written from one call site.

    The half-edges are stored twice and a one-sided write would pass a read from the low-numbered end
    and answer None from the other, which is the shape of bug that survives a whole test suite.
    """
    mol, c1, c2, c3, f = _propene()
    mol.set_bond_cip(c3, c2, 'Z')            # stated from the HIGH end
    assert mol.bond_cip_of(c2, c3) == 'Z'
    assert mol.bond_cip_of(c3, c2) == 'Z'
    assert mol.bond_cips() == {(c2, c3): 'Z'}, 'bond_cips must list each bond once, canonically'
    # and stating it from the other end is the same statement, not a second one
    mol.set_bond_cip(c2, c3, 'E')
    assert mol.bond_cip_of(c3, c2) == 'E'
    assert mol.bond_cips() == {(c2, c3): 'E'}


def test_a_bond_set_twice_in_one_scope_keeps_the_last_word():
    """And spends one slot doing it: the replay overwrites in place, and the region is sized by bonds.

    Without the in-place overwrite a scope that set every bond twice would run past the end of a region
    sized for one entry per bond.
    """
    mol, c1, c2, c3, f = _propene()
    with mol.edit() as e:
        e.set_bond_cip(c2, c3, 'E')
        e.set_bond_cip(c3, c2, 'Z')          # same bond, other end, later record
        e.set_bond_cip(c2, c3, 'M')
    assert mol.bond_cips() == {(c2, c3): 'M'}


def test_every_bond_labelled_twice_in_one_scope_stays_inside_its_region():
    """The bound the in-place overwrite protects, exercised at the width where it would be exceeded."""
    mol, ring, methyl = _benzene_with_a_methyl()
    pairs = list(zip(ring, ring[1:] + ring[:1])) + [(ring[0], methyl)]
    with mol.edit() as e:
        for a, b in pairs:
            e.set_bond_cip(a, b, 'E')
        for a, b in pairs:
            e.set_bond_cip(b, a, 'Z')        # the reversed pair must find the existing entry
    assert len(mol.bond_cips()) == len(pairs)
    assert set(mol.bond_cips().values()) == {'Z'}


# ------------------------------------------------------------------------------------------------
# SURVIVAL: remap, copy, pack, property-only edits, kekule/thiele
# ------------------------------------------------------------------------------------------------

def test_remap_moves_the_descriptors_with_the_ids():
    """Only labels move, so every descriptor must answer under the new label and none under the old."""
    mol, c1, c2, c3, f = _propene()
    mol.set_atom_cip(c2, 'r')
    mol.set_atom_cip(c3, 'S')
    mol.set_bond_cip(c2, c3, 'M')
    mol.remap({c2: 200, c3: 300})

    assert mol.atom_cip_of(200) == 'r'
    assert mol.atom_cip_of(300) == 'S'
    assert mol.bond_cip_of(200, 300) == 'M'
    assert mol.bond_cip_of(300, 200) == 'M'
    assert mol.atom_cips() == {200: 'r', 300: 'S'}
    assert mol.bond_cips() == {(200, 300): 'M'}
    assert mol.cip_log == (), 'a relabelling is not a change to the molecule'


def test_copy_carries_the_descriptors():
    mol, c1, c2, c3, f = _propene()
    mol.set_atom_cip(c2, 'p')
    mol.set_bond_cip(c2, c3, 'P')
    clone = mol.copy()
    assert clone.atom_cips() == {c2: 'p'}
    assert clone.bond_cips() == {(c2, c3): 'P'}
    # free, and the reason is worth stating: a copy shares the arena, and the descriptors live in it
    assert clone.shares_arena_with(mol)


def test_a_copy_starts_with_an_empty_drop_log():
    """The descriptors travel with the arena; the drop history does not, and that is deliberate.

    `cip_log` records what THIS handle's edits lost, the same standing `sgroup_log` has.  A copy has
    made no edits, so it has lost nothing, and a copied log would attribute one handle's losses to
    another.  Read the log from the container you edited.
    """
    mol, c1, c2, c3, f = _propene()
    mol.set_atom_cip(c2, 'R')
    with mol.edit() as e:
        e.delete_atom(f)
    assert mol.cip_log, 'the delete should have logged'
    assert mol.copy().cip_log == ()


def test_a_pack_round_trip_keeps_every_descriptor():
    mol, ring, methyl = _benzene_with_a_methyl()
    for atom, descriptor in zip(ring, ATOM_DESCRIPTORS):
        mol.set_atom_cip(atom, descriptor)
    for (a, b), descriptor in zip(zip(ring, ring[1:]), BOND_DESCRIPTORS):
        mol.set_bond_cip(a, b, descriptor)

    back = _round_trip(mol)
    assert back.atom_cips() == mol.atom_cips()
    assert back.bond_cips() == mol.bond_cips()
    # both half-edges came back, not just the canonical one the packer wrote
    for a, b in zip(ring, ring[1:]):
        assert back.bond_cip_of(b, a) == mol.bond_cip_of(a, b)
    assert back.to_bytes() == mol.to_bytes()


@pytest.mark.parametrize('descriptor', ATOM_DESCRIPTORS)
def test_each_atom_descriptor_survives_the_bytes_individually(descriptor):
    """Parametrised rather than one molecule carrying all eight, so a code that packs into the wrong
    bit is attributed to the value that did it instead of failing one assertion in a heap."""
    mol, c1, c2, c3, f = _propene()
    mol.set_atom_cip(c2, descriptor)
    assert _round_trip(mol).atom_cip_of(c2) == descriptor


@pytest.mark.parametrize('descriptor', BOND_DESCRIPTORS)
def test_each_bond_descriptor_survives_the_bytes_individually(descriptor):
    mol, c1, c2, c3, f = _propene()
    mol.set_bond_cip(c2, c3, descriptor)
    back = _round_trip(mol)
    assert back.bond_cip_of(c2, c3) == descriptor
    assert back.bond_cip_of(c3, c2) == descriptor


@pytest.mark.parametrize('edit', ['charge', 'isotope', 'radical', 'map_number', 'hydrogens',
                                  'stereo', 'xy', 'wedge', 'stereo_group'])
def test_a_property_only_edit_keeps_the_descriptors(edit):
    """None of these changes which atoms exist, which are bonded, or what a bond's order is.

    ISOTOPE IS THE INTERESTING ONE and it is here on purpose.  CIP Rule 2 ranks by mass, so an isotope
    edit CAN change a computed descriptor -- and it still does not drop, because this layer did not
    compute the stored one.  An assignment algorithm that read a stored descriptor as an input rather
    than recomputing would be wrong for a reason no drop rule here could repair; see RULES.
    """
    mol, c1, c2, c3, f = _propene()
    mol.set_atom_cip(c2, 'S')
    mol.set_bond_cip(c2, c3, 'Z')
    with mol.edit() as e:
        if edit == 'charge':
            e.set_charge(c1, 1)
        elif edit == 'isotope':
            e.set_isotope(c2, 13)
        elif edit == 'radical':
            e.set_radical(c1, True)
        elif edit == 'map_number':
            e.set_map_number(c2, 7)
        elif edit == 'hydrogens':
            e.set_hydrogens(c1, 2)
        elif edit == 'stereo':
            e.set_stereo(c2, True)
        elif edit == 'xy':
            e.set_xy(c2, 1.5, -2.25)
        elif edit == 'wedge':
            e.set_wedge(c2, c1, 1)
        else:
            e.set_stereo_group(c2, 1, 1)

    assert mol.atom_cips() == {c2: 'S'}, f'{edit} dropped the atom descriptor'
    assert mol.bond_cips() == {(c2, c3): 'Z'}, f'{edit} dropped the bond descriptor'
    assert mol.cip_log == ()


def test_kekule_keeps_the_descriptors():
    """One of the two operations in the library allowed to change a representation, and so exempt.

    A descriptor is an assertion the INPUT made about the MOLECULE, and a molecule spelled aromatic and
    the same molecule spelled Kekule are one molecule.  The rule keys on the OPERATION rather than on
    the field, because these two reach the journal as ordinary order changes and cannot be told from a
    caller's own `set_order` any other way.
    """
    mol, ring, methyl = _benzene_with_a_methyl()
    assert mol.thiele().changed, 'fixture must actually aromatise, or this test proves nothing'
    mol.set_atom_cip(ring[0], 'R')
    mol.set_bond_cip(ring[0], ring[1], 'M')
    assert mol.order_of(ring[0], ring[1]) == 4

    result = mol.kekule()
    assert result.changed, 'kekule must actually change the orders here'
    assert mol.order_of(ring[0], ring[1]) != 4
    assert mol.atom_cips() == {ring[0]: 'R'}
    assert mol.bond_cips() == {(ring[0], ring[1]): 'M'}
    assert mol.cip_log == ()


def test_thiele_keeps_the_descriptors():
    mol, ring, methyl = _benzene_with_a_methyl()
    mol.set_atom_cip(ring[0], 'S')
    mol.set_bond_cip(ring[0], methyl, 'P')
    assert mol.order_of(ring[0], ring[1]) == 2

    result = mol.thiele()
    assert result.changed, 'thiele must actually change the orders here'
    assert mol.order_of(ring[0], ring[1]) == 4
    assert mol.atom_cips() == {ring[0]: 'S'}
    assert mol.bond_cips() == {(ring[0], methyl): 'P'}
    assert mol.cip_log == ()


def test_a_hand_written_order_change_is_not_exempt_even_next_to_a_kekule():
    """The exemption belongs to those two functions and not to the op they emit.

    If the flag ever leaked -- set once and not cleared, or set around a scope a caller can reach --
    every `set_order` in the library would silently start preserving descriptors.  This is the test
    that notices, and it runs the caller's edit right after a real `kekule` so a flag left standing
    would be caught rather than merely absent.
    """
    mol, ring, methyl = _benzene_with_a_methyl()
    mol.thiele()
    mol.set_atom_cip(ring[0], 'R')
    mol.kekule()
    assert mol.atom_cips() == {ring[0]: 'R'}

    with mol.edit() as e:
        e.set_order(ring[0], methyl, 2)
    assert mol.atom_cips() == {}
    assert any('dropped' in line for line in mol.cip_log)


# ------------------------------------------------------------------------------------------------
# THE DROP: one test per op, because a dead arm in an if/elif chain is not a warning anywhere
# ------------------------------------------------------------------------------------------------

def _labelled_pair():
    mol, c1, c2, c3, f = _propene()
    mol.set_atom_cip(c2, 'R')
    mol.set_bond_cip(c2, c3, 'E')
    return mol, c1, c2, c3, f


def test_adding_an_atom_drops_the_descriptors():
    """Even an atom bonded to nothing: a ranking reads the whole molecule, and this arm is easy to
    leave dead.

    A staleness flag set by one `elif` naming four ops is where that happens -- two of the four are
    matched by earlier arms of the same chain that count them, so adds invalidate nothing and nothing
    complains.  One test per op is the only thing that catches that class of defect.
    """
    mol, c1, c2, c3, f = _labelled_pair()
    with mol.edit() as e:
        e.add_atom('N')
    assert mol.atom_cips() == {}
    assert mol.bond_cips() == {}
    assert mol.cip_log == ('1 bond CIP descriptor(s) dropped: the molecule changed',
                           '1 atom CIP descriptor(s) dropped: the molecule changed')


def test_adding_a_bond_drops_the_descriptors():
    mol, c1, c2, c3, f = _labelled_pair()
    with mol.edit() as e:
        n = e.add_atom('N')
    mol.set_atom_cip(c2, 'R')
    mol.set_bond_cip(c2, c3, 'E')
    with mol.edit() as e:
        e.add_bond(c1, n)
    assert mol.atom_cips() == {}
    assert mol.bond_cips() == {}
    assert mol.cip_log[-2:] == ('1 bond CIP descriptor(s) dropped: the molecule changed',
                                '1 atom CIP descriptor(s) dropped: the molecule changed')


def test_deleting_an_atom_drops_the_descriptors_even_far_from_the_centre():
    """The methyl carbon is two bonds from the labelled centre and its removal still invalidates.

    Being loudly conservative is the right side to err on: a dropped descriptor is logged and can be
    recomputed, while a carried wrong 'R' is a different molecule to a chemist and nothing downstream
    can tell.
    """
    mol, c1, c2, c3, f = _labelled_pair()
    with mol.edit() as e:
        e.delete_atom(c1)
    assert mol.atom_cips() == {}
    assert mol.bond_cips() == {}
    assert len(mol.cip_log) == 2


def test_deleting_a_bond_drops_the_descriptors():
    mol, c1, c2, c3, f = _labelled_pair()
    with mol.edit() as e:
        e.delete_bond(c3, f)
    assert mol.atom_cips() == {}
    assert mol.bond_cips() == {}
    assert len(mol.cip_log) == 2


def test_changing_a_bond_order_by_hand_drops_the_descriptors():
    mol, c1, c2, c3, f = _labelled_pair()
    with mol.edit() as e:
        e.set_order(c1, c2, 2)
    assert mol.atom_cips() == {}
    assert mol.bond_cips() == {}
    assert len(mol.cip_log) == 2


def test_the_log_names_how_many_of_each_kind_went():
    """Counts, not just a flag: a log line saying "something was dropped" cannot be acted on."""
    mol, ring, methyl = _benzene_with_a_methyl()
    for atom in ring[:3]:
        mol.set_atom_cip(atom, 'R')
    mol.set_bond_cip(ring[0], ring[1], 'E')
    mol.set_bond_cip(ring[2], ring[3], 'Z')
    with mol.edit() as e:
        e.delete_atom(methyl)
    assert mol.cip_log == ('2 bond CIP descriptor(s) dropped: the molecule changed',
                           '3 atom CIP descriptor(s) dropped: the molecule changed')


def test_an_edit_that_drops_nothing_logs_nothing():
    """An unlabelled molecule must not accumulate log noise on every structural edit."""
    mol, c1, c2, c3, f = _propene()
    with mol.edit() as e:
        e.delete_atom(f)
        e.add_atom('Cl')
    assert mol.cip_log == ()


def test_a_drop_is_only_recoverable_from_the_log():
    """Storage cannot tell a never-labelled atom from a dropped one -- both hold code 0.

    So the pair of molecules below are byte-identical, and the ONLY surviving difference is the log.
    That is why the log exists, and why this test asserts on the bytes rather than on a reader.
    """
    labelled, c1, c2, c3, f = _labelled_pair()
    with labelled.edit() as e:
        e.delete_atom(f)

    plain, p1, p2, p3, pf = _propene()
    with plain.edit() as e:
        e.delete_atom(pf)

    assert labelled.to_bytes() == plain.to_bytes()
    assert labelled.cip_log and plain.cip_log == ()


def test_a_descriptor_stated_in_the_same_scope_as_the_edit_still_wins():
    """The order a parser needs: build the molecule and state its descriptors in ONE scope.

    The drop runs after the seed and before the replay, so a descriptor journalled in the same scope as
    the structural edit is applied on top of the cleared state instead of being wiped by it.  Reversing
    those two would make every descriptor a reader states unreachable, which is the whole path.
    """
    mol = MoleculeContainer()
    with mol.edit() as e:
        c1 = e.add_atom('C', implicit_h=3)
        c2 = e.add_atom('C', implicit_h=1)
        c3 = e.add_atom('C', implicit_h=1)
        e.add_bond(c1, c2)
        e.add_bond(c2, c3, 2)
        e.set_atom_cip(c2, 'S')
        e.set_bond_cip(c2, c3, 'Z')
    assert mol.atom_cips() == {c2: 'S'}
    assert mol.bond_cips() == {(c2, c3): 'Z'}
    assert mol.cip_log == (), 'nothing was lost: the descriptors are about the molecule just built'


def test_the_scope_that_states_a_descriptor_wins_wherever_in_the_scope_it_states_it():
    """The drop clears the PRE-SCOPE descriptors, and every descriptor the scope states survives.

    So the position of the statement inside the scope does not matter -- stated before the delete or
    after it, it is applied on top of the cleared state either way.  That is worth pinning rather than
    leaving to the implementation: a rule that read the journal record by record would make a parser's
    output depend on where in its own scope a descriptor happened to land, and a parser that emits
    atoms, bonds and descriptors in file order has no control over that.

    The two descriptors here differ only in when they were stated, and only the older one goes.
    """
    mol, c1, c2, c3, f = _propene()
    mol.set_atom_cip(c2, 'R')                # stated in an earlier scope
    with mol.edit() as e:
        e.set_atom_cip(c3, 'S')              # stated BEFORE the change
        e.delete_atom(f)
    assert mol.atom_cips() == {c3: 'S'}
    assert mol.cip_log == ('1 atom CIP descriptor(s) dropped: the molecule changed',)


# ------------------------------------------------------------------------------------------------
# TWO IDENTITIES: in the bytes, out of the canonical form
# ------------------------------------------------------------------------------------------------

def test_an_unlabelled_molecule_holds_the_bytes_it_held_before_cip_existed():
    """The regression that catches a widened record or a stolen default.

    Stated as a property of the bytes rather than as a stored blob, because a blob would also fail for
    every unrelated format change and could not say which happened.  The property is exact: CIP took
    the low nibble of a `reserved` word that was already serialised and three spare bits of a flags
    word that was already serialised, so an unlabelled molecule's buffer is byte-for-byte the one the
    previous build wrote -- every atom's reserved word is zero and every half-edge's CIP field is zero.
    """
    mol, ring, methyl = _benzene_with_a_methyl()
    data = mol.to_bytes()
    atoms_at = struct.unpack_from('<I', data, 24 + 8 * _SEG_ATOMS)[0]
    edges_at = struct.unpack_from('<I', data, 24 + 8 * _SEG_CSR_EDGE)[0]

    for i in range(mol.atom_count):
        reserved = struct.unpack_from('<I', data, atoms_at + _ATOM_RECORD * i + _ATOM_RESERVED)[0]
        assert reserved == 0, f'atom {i} carries reserved bits an older build would not have written'
    for k in range(2 * mol.bond_count):
        flags = struct.unpack_from('<H', data, edges_at + _HALFEDGE_RECORD * k + _HALFEDGE_FLAGS)[0]
        assert flags & 0x1c == 0, f'half-edge {k} carries CIP bits'


def test_a_labelled_molecule_and_its_unlabelled_twin_have_different_bytes():
    """`to_bytes` IS an identity in this version, so a descriptor must be inside it.

    An atom descriptor and a bond descriptor are asserted separately: they live in two different
    records, and a change that serialised one of them and not the other would pass a single assertion.
    """
    labelled, c1, c2, c3, f = _propene()
    plain, p1, p2, p3, pf = _propene()
    assert labelled.to_bytes() == plain.to_bytes(), 'the twins must start identical'

    labelled.set_atom_cip(c2, 'R')
    atom_only = labelled.to_bytes()
    assert atom_only != plain.to_bytes()

    labelled.set_atom_cip(c2, None)
    labelled.set_bond_cip(c2, c3, 'E')
    bond_only = labelled.to_bytes()
    assert bond_only != plain.to_bytes()
    assert bond_only != atom_only


def test_two_molecules_differing_only_in_a_descriptor_stay_equal_and_stay_one_hash():
    """The other identity, and CIP must stay OUT of it.

    A descriptor is a statement ABOUT a molecule and not a part of what makes it that molecule, so the
    canonical form must not read it: `==` and `hash` answer the same for the labelled and unlabelled
    twins even though their bytes differ.  A dict keyed on molecules must not grow a second entry
    because someone annotated one of them.
    """
    labelled, c1, c2, c3, f = _propene()
    plain, p1, p2, p3, pf = _propene()
    labelled.set_atom_cip(c2, 'R')
    labelled.set_bond_cip(c2, c3, 'E')

    assert labelled.to_bytes() != plain.to_bytes()
    assert labelled == plain
    assert hash(labelled) == hash(plain)
    assert len({labelled, plain}) == 1
    assert labelled.smiles == plain.smiles
    assert labelled.atoms_order == plain.atoms_order


# ------------------------------------------------------------------------------------------------
# the reads that must refuse
# ------------------------------------------------------------------------------------------------

@pytest.mark.parametrize('read', ['atom_cip_of', 'bond_cip_of', 'atom_cips', 'bond_cips'])
def test_reading_a_descriptor_with_pending_edits_is_refused(read):
    """The arena still holds the pre-scope state, so these would answer from stale data."""
    mol, c1, c2, c3, f = _propene()
    mol.set_atom_cip(c2, 'R')
    with pytest.raises(RuntimeError, match='pending edits'):
        with mol.edit() as e:
            e.set_charge(c1, 1)
            if read == 'atom_cip_of':
                mol.atom_cip_of(c2)
            elif read == 'bond_cip_of':
                mol.bond_cip_of(c2, c3)
            elif read == 'atom_cips':
                mol.atom_cips()
            else:
                mol.bond_cips()


def test_a_substructure_is_descriptor_free_and_says_so_in_its_docstring():
    """A cut destroys the frame a descriptor was stated in, the same way it destroys a parity.

    THE LOSS IS NOT LOGGED HERE and that is the documented exception: an edit does not ask to lose a
    descriptor, so an edit reports it, while a cut is a caller asking for a smaller molecule and
    `substructure`'s docstring is where the answer lives.
    """
    mol, c1, c2, c3, f = _propene()
    mol.set_atom_cip(c2, 'R')
    mol.set_bond_cip(c2, c3, 'E')
    cut = mol.substructure([c1, c2, c3])
    assert cut.atom_cips() == {}
    assert cut.bond_cips() == {}
    assert 'CIP' in MoleculeContainer.substructure.__doc__


def test_a_union_drops_the_descriptors_and_logs_it():
    """Conservative, and logged, which is the pairing the drop rule promises everywhere else.

    A union does not touch either component's connectivity, so a future ranking might well survive it.
    Storage does not get to make that call: the rule is stated on the OPERATIONS in the journal, a
    union adds atoms and bonds, and the exemption list is closed by RULES rather than extended here.
    """
    mol, c1, c2, c3, f = _propene()
    mol.set_atom_cip(c2, 'R')
    other, o1, o2, o3, of = _propene()
    joined = mol.union(other)
    assert joined.atom_cips() == {}
    assert any('dropped' in line for line in joined.cip_log)
