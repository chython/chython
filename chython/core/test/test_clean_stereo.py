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
"""`clean_stereo()` -- every kind of stereo state gone, unconditionally.

THE DIFFERENCE FROM `validate_stereo`, which is the reason this method exists at all.
`validate_stereo` asks a question -- "which stated parities can this constitution justify" -- and
clears only the answers it does not like.  This asks nothing.  A caller reaching for it has decided
that whatever the molecule says about configuration is not to be trusted or not to be kept, and the
only correct outcome is a flat molecule.  `clean_stereo` is chython 2's name for it, kept so that no
consumer has to write a partial version of its own.

FOUR KINDS OF STATE, AND WIPING ONLY THE FIRST IS THE TRAP.  Configuration lives in four places:
parities on the atom, wedges on the edge, ABS/AND/OR group membership per atom, and CIP descriptors on
atoms and bonds.  Leave the wedges and the CTfile writer draws them again and the CTfile reader
derives the parities back from them, so the wipe does not survive one round trip -- pinned in
`chython/formats/ctfile/test/test_fidelity.py`, which is where the reader and the writer both are.
Leave a group and an AND membership names a configuration that no longer exists.  Leave a stored
`(R)` and an external consumer -- which is who stored CIP is FOR -- reads a descriptor off an atom
with no parity.
"""
from chython.core import MoleculeContainer


def _flat_butane():
    """CC(Cl)C(Cl)C with coordinates: two tetrahedral centres and nothing stated about them."""
    m = MoleculeContainer()
    xy = [(0.0, 0.0), (0.9, 0.5), (0.9, 1.5), (1.8, 0.0), (1.8, -1.0), (2.7, 0.5)]
    with m.edit():
        sids = [m.add_atom(e, implicit_h=h) for e, h in
                (('C', 3), ('C', 1), ('Cl', 0), ('C', 1), ('Cl', 0), ('C', 3))]
        for i, j in ((0, 1), (1, 2), (1, 3), (3, 4), (3, 5)):
            m.add_bond(sids[i], sids[j], 1)
        for s, (x, y) in zip(sids, xy):
            m.set_xy(s, x, y)
    return m, sids


def _all_four_kinds():
    """`_flat_butane` carrying all four kinds of stereo state at once.

    Both centres get an ODD parity (2), and that is not decoration: word IV bit 6 is set for parity 2
    and only for parity 2, so the round-trip test can detect a clear that missed the byte.
    """
    m, sids = _flat_butane()
    with m.edit():
        m.set_parity(sids[1], 2)
        m.set_parity(sids[3], 2)
        m.set_wedge(sids[1], sids[2], 1)              # narrow at the centre, up
        m.set_stereo_group(sids[1], 3, 1)             # AND1 -- the racemate a CTfile spells
        m.set_stereo_group(sids[3], 3, 1)
        m.set_atom_cip(sids[1], 'R')
        m.set_bond_cip(sids[1], sids[3], 'E')
    return m, sids


# ----------------------------------------------------------------------------------------------
# the wipe
# ----------------------------------------------------------------------------------------------

def test_every_kind_of_stereo_state_is_cleared():
    m, sids = _all_four_kinds()
    # the fixture really does carry all four before the call; without this the test below would pass
    # on a molecule that never had a wedge or a group to lose
    assert [m.parity_of(s) for s in (sids[1], sids[3])] == [2, 2]
    assert m.wedges() and m.stereo_groups() and m.atom_cips() and m.bond_cips()

    report = m.clean_stereo()

    assert [m.parity_of(s) for s in m.atom_numbers] == [0] * 6
    assert [m.stereo_of(s) for s in m.atom_numbers] == [False] * 6
    assert m.wedges() == []
    assert m.stereo_groups() == {}
    assert m.atom_cips() == {}
    assert m.bond_cips() == {}
    assert report == {'parities': [sids[1], sids[3]],
                      'wedges': [(sids[1], sids[2], 1)],
                      'stereo_groups': {(3, 1): [sids[1], sids[3]]},
                      'atom_cips': {sids[1]: 'R'},
                      'bond_cips': {(sids[1], sids[3]): 'E'}}


def test_the_report_is_what_the_five_readers_said_before_the_wipe():
    """The report's shape is not invented here: each value is one reader's own answer, verbatim.

    That is the whole design of the return value -- `validate_stereo` reports one kind of state and a
    flat list of stable ids says all there is to say about it, while this touches five readers, and a
    union list would tell a caller that atom 2 "had something" without saying what.  So the report is
    keyed by reader, and a key is absent when its reader was empty.
    """
    m, sids = _all_four_kinds()
    before = {'parities': sorted(s for s in m.atom_numbers if m.parity_of(s)),
              'wedges': m.wedges(),
              'stereo_groups': m.stereo_groups(),
              'atom_cips': m.atom_cips(),
              'bond_cips': m.bond_cips()}
    assert m.clean_stereo() == before


def test_a_second_call_reports_nothing_and_is_a_pure_read():
    m, sids = _all_four_kinds()
    assert m.clean_stereo()
    other = m.copy()
    assert m.clean_stereo() == {}, 'idempotent: there is nothing left to clear'
    assert m.shares_arena_with(other), 'an empty report is a pure read: no clone, no _gen bump'
    assert m.generation == other.generation


def test_a_molecule_with_no_stereo_reports_nothing():
    m, sids = _flat_butane()
    other = m.copy()
    assert m.clean_stereo() == {}
    assert m.shares_arena_with(other)
    assert m.generation == other.generation


def test_an_empty_molecule_survives_it():
    m = MoleculeContainer()
    assert m.clean_stereo() == {}


# ----------------------------------------------------------------------------------------------
# clone on write: `copy()` shares the arena outright, so the wipe must not reach into it
# ----------------------------------------------------------------------------------------------

def test_a_shared_arena_keeps_its_stereo_when_the_copy_is_cleaned():
    m, sids = _all_four_kinds()
    other = m.copy()
    assert m.shares_arena_with(other)

    other.clean_stereo()

    assert not m.shares_arena_with(other)
    assert [m.parity_of(s) for s in (sids[1], sids[3])] == [2, 2], \
        'the wipe must go into a clone; writing SEG_PARITY in place strips every sharer'
    assert m.wedges() == [(sids[1], sids[2], 1)]
    assert m.stereo_groups() == {(3, 1): [sids[1], sids[3]]}
    assert m.atom_cips() == {sids[1]: 'R'}
    assert m.bond_cips() == {(sids[1], sids[3]): 'E'}


def test_a_shared_arena_keeps_its_stereo_when_the_original_is_cleaned():
    m, sids = _all_four_kinds()
    other = m.copy()

    m.clean_stereo()

    assert [other.parity_of(s) for s in (sids[1], sids[3])] == [2, 2]
    assert other.wedges() == [(sids[1], sids[2], 1)]
    assert other.stereo_groups() == {(3, 1): [sids[1], sids[3]]}
    assert other.atom_cips() == {sids[1]: 'R'}
    assert other.bond_cips() == {(sids[1], sids[3]): 'E'}


def test_the_stale_unit_table_does_not_travel_in_the_clone():
    """`structure_clone` copies derived segments verbatim, marks included, so the table is retired.

    Same guarantee, same measurement and the same sign as `validate_stereo`'s
    `test_a_reported_parity_is_cleared`: the marks in the copied table were computed against the
    parities this call removes, so the clone must SHRINK by exactly one table and the next reader
    must put it back to the byte.  Delete the invalidate and the shrink assertion fails with the two
    numbers equal.
    """
    m, sids = _all_four_kinds()
    m.stereo_units()                       # the table exists, marked, before the wipe
    grown = m.total_len

    assert m.clean_stereo()
    assert m.total_len < grown, 'the stale table was retired in the clone'
    m.stereo_units()
    assert m.total_len == grown, 'and the next reader derived a fresh one of the same size'


# ----------------------------------------------------------------------------------------------
# what must NOT be touched
# ----------------------------------------------------------------------------------------------

def test_coordinates_are_not_dropped():
    """A layout is not a configuration.  Ruled: `clean_stereo` wipes stereo and leaves the drawing.

    The temptation runs the other way -- the parity a CTfile states IS derived from the coordinates,
    so dropping them would make the wipe unrecoverable-by-construction.  It would also destroy the
    only thing a depiction has to work with, on a molecule the caller asked to flatten and not to
    forget, and `clean2d()` is the call that replaces a layout.
    """
    m, sids = _all_four_kinds()
    before = [m.xy_of(s) for s in m.atom_numbers]
    m.clean_stereo()
    assert m.has_coordinates
    assert [m.xy_of(s) for s in m.atom_numbers] == before


def test_the_constitution_is_untouched():
    m, sids = _all_four_kinds()

    def snapshot():
        return (m.atom_count, m.bond_count, sorted(m.atom_numbers),
                [m.element_of(s) for s in m.atom_numbers],
                [m.implicit_h_of(s) for s in m.atom_numbers],
                sorted((b.n, b.m, b.order) for b in m.bonds()))

    before = snapshot()
    m.clean_stereo()
    assert snapshot() == before


# ----------------------------------------------------------------------------------------------
# the parity byte
# ----------------------------------------------------------------------------------------------

def test_a_cleared_parity_does_not_come_back_through_a_round_trip():
    """`clean_stereo` zeroes the byte through `structure_clear_parities`, and `to_bytes` carries the
    segment, so the wipe persists across a round trip.
    """
    m, sids = _all_four_kinds()
    assert m.clean_stereo()
    again = MoleculeContainer.from_bytes(m.to_bytes())
    assert [again.parity_of(s) for s in (sids[1], sids[3])] == [0, 0], \
        'the segment carries the cleared byte, so a zero comes back a zero'
    assert [again.stereo_of(s) for s in (sids[1], sids[3])] == [False, False]
    assert again.wedges() == [] and again.stereo_groups() == {}
    assert again.atom_cips() == {} and again.bond_cips() == {}


def test_the_feature_words_match_a_round_trip_after_the_wipe():
    """Ruling F78: a parity writer outside `rebuild_derived` maintains feature word IV itself.

    Word IV screens the parity VALUE bit at its bit 6 and `features_of()` hands the words to Python
    verbatim, so a wipe that forgot `refresh_parity_features` would leave the words stating the signs
    the arena no longer holds.  A `from_bytes` rebuild derives them from scratch, so it is the oracle.
    The fixture's parities are odd for the reason stated on it -- an even one never sets bit 6 and
    would pass with the maintenance removed.
    """
    m, sids = _all_four_kinds()
    assert m.clean_stereo()
    fresh = MoleculeContainer.from_bytes(m.to_bytes())
    assert [m.features_of(s) for s in m.atom_numbers] == [fresh.features_of(s) for s in m.atom_numbers]
    assert m._union_feature_words == fresh._union_feature_words
