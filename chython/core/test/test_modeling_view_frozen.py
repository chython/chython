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
"""`modeling_view()` against the answer recorded before the array kernel replaced its body.

A RECORDED FIXTURE AND NOT A SECOND IMPLEMENTATION.  Keeping the dict-of-dicts union alive as the
kernel's oracle would be a second derivation of one thing, and its drift would be visible only to the
test comparing the two to each other.
"""
from pytest import mark, skip

from .modeling_view_corpus import RECORDS, load


FROZEN = load()


# DECLARED CONVENTION CHANGES.  A record here is pinned by a named test below instead, because the
# recorded answer states the older convention.
#
# `colliding_map_numbers` has two atoms claiming map number 1: the dict union let the second overwrite
# the first, and the bond between them then closed on itself, so the recorded state carries a degree
# from a self-loop.  The union holds one row per map number and no self-loop, so the second atom
# contributes only its entry in `collisions`.
#
# The other three carry an atom with no map number, which the recorded answer counts and leaves out.
# Leaving it out takes the degree of every neighbour that stayed down with it -- the recorded
# tetrafluoroammonium nitrogen has four bonds and a degree of 0 -- so an unmapped atom now holds a union
# row on the side it came from, keyed negatively.  `unmapped` is unchanged and still recorded, which is
# what `test_what_could_not_be_placed_...` asserts for all four.
DEVIATIONS = {'colliding_map_numbers', 'esterification_partly_mapped', 'unmapped_hydrogenation',
              'tetrafluoroammonium_unknown_h'}


@mark.parametrize('name', sorted(FROZEN))
def test_the_state_quintuples_and_their_order_are_what_was_recorded(name):
    """Order too: a tokenizer takes its atom order from `states`, so a dict compare is too weak."""
    if name in DEVIATIONS:
        skip('declared convention change; see DEVIATIONS')
    recorded = [tuple(row) for row in FROZEN[name]['states']]
    got = [(n, *state) for n, state in RECORDS[name].modeling_view().states.items()]
    assert got == recorded


@mark.parametrize('name', sorted(FROZEN))
def test_every_union_bond_and_its_two_orders_are_what_was_recorded(name):
    if name in DEVIATIONS:
        skip('declared convention change; see DEVIATIONS')
    recorded = {(n, m): (before, after) for n, m, before, after in FROZEN[name]['union_bonds']}
    assert RECORDS[name].modeling_view().union_bonds == recorded


@mark.parametrize('name', sorted(FROZEN))
def test_what_could_not_be_placed_is_reported_as_it_was_recorded(name):
    view = RECORDS[name].modeling_view()
    assert view.unmapped == FROZEN[name]['unmapped']
    assert view.collisions == {side: tuple(numbers)
                               for side, numbers in FROZEN[name]['collisions'].items()}


def test_the_corpus_reaches_the_three_branches_it_was_written_for():
    """A fixture whose sparse records came out empty for the wrong reason pins nothing.

    Read of the RECORDING and not of the current answer: `unmapped_hydrogenation` recorded no atom
    because the convention of the day left an unmapped one out, which is the branch it was written to
    reach and is now `DEVIATIONS`.
    """
    assert FROZEN['unmapped_hydrogenation']['states'] == [], \
        'the unmapped record recorded a union atom; it was not written against that convention'
    assert FROZEN['colliding_map_numbers']['collisions']['reactants'], \
        'the colliding record recorded no collision; it no longer reaches that branch'
    assert any(row[2] == 15 for row in FROZEN['tetrafluoroammonium_unknown_h']['states']), \
        'the H_UNKNOWN record recorded no sentinel; its nitrogen now has a derivable count'


def test_a_colliding_map_number_keeps_its_first_claim_and_invents_no_self_loop():
    """The convention DEVIATIONS names, stated as an assertion rather than left as an absence.

    Map number 1 is claimed twice; its only bond joins the two claimants, so the union cannot place it.
    The product side is untouched by the collision and keeps its bond.
    """
    view = RECORDS['colliding_map_numbers'].modeling_view()
    assert sorted(view.states) == [1, 2, 3]
    assert view.states[1] == (6, 3, 0, 3, 0), 'no bond survives, so no degree on either side'
    assert view.states[2] == (6, 3, 1, 3, 1), 'a product-only pair, unaffected by the collision'
    assert view.union_bonds == {(2, 3): (0, 1)}
    assert view.collisions['reactants'] == (1,)


def test_an_unmapped_atom_holds_a_row_and_leaves_its_neighbour_a_degree():
    """The other convention DEVIATIONS names, on the record whose recording shows the cost.

    The nitrogen of tetrafluoroammonium has four bonds and the record numbers none of the fluorines, so
    the recorded answer gives it a degree of 0 on both sides while the row says it carries four bonds
    worth of hydrogens.  Each fluorine now holds a row, so the degree is 4.  The price is on the same
    record: the two sides are one structure copied, and with nothing pairing the fluorines the union
    reads four bonds breaking and four forming.  `unmapped` counts all eight.
    """
    view = RECORDS['tetrafluoroammonium_unknown_h'].modeling_view()
    assert view.states[1] == (7, 15, 4, 15, 4), 'four bonds, so four heavy neighbours per side'
    assert [n for n in view.states if n < 0] == [-1, -2, -3, -4, -5, -6, -7, -8]
    assert sorted(view.union_bonds.values()) == [(0, 1)] * 4 + [(1, 0)] * 4
    assert view.unmapped == {'reactants': 4, 'products': 4}


def test_a_partly_mapped_record_places_what_the_numbered_part_gains_and_loses():
    """The esterification whose methanol carries no number: the mapped part still reads correctly.

    Map 4, the acid's hydroxyl oxygen, is on both sides and gains the arriving methyl -- one heavy
    neighbour before, two after.  Recorded, it had one after, the bond to an unmapped atom being no
    bond at all.
    """
    view = RECORDS['esterification_partly_mapped'].modeling_view()
    assert view.states[4] == (8, 1, 1, 0, 2)
    assert view.union_bonds[(-3, 4)] == (0, 1), 'the arriving methyl bonds the ester oxygen'
    assert view.union_bonds[(-2, -1)] == (1, 0), 'and the methanol the record does not follow leaves'


def test_a_record_with_no_mapping_at_all_is_its_two_sides_side_by_side():
    """Nothing is conserved because nothing in the record says anything is, and `unmapped` says so."""
    view = RECORDS['unmapped_hydrogenation'].modeling_view()
    assert list(view.states) == [-1, -2, -3, -4, -5, -6]
    assert view.union_bonds == {(-2, -1): (2, 0), (-4, -3): (1, 0), (-6, -5): (0, 1)}
    assert view.unmapped == {'reactants': 4, 'products': 2}
