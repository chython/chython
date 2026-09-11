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


# A DECLARED CONVENTION CHANGE, one record wide.  `colliding_map_numbers` has two atoms claiming map
# number 1: the dict union let the second overwrite the first, and the bond between them then closed on
# itself, so the recorded state carries a degree from a self-loop.  The union holds one row per map
# number and no self-loop, so the second atom contributes only its entry in `collisions` -- which the
# record is in the corpus to pin, and which `test_what_could_not_be_placed_...` still asserts.
DEVIATIONS = {'colliding_map_numbers'}


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
    recorded = {(n, m): (before, after) for n, m, before, after in FROZEN[name]['union_bonds']}
    assert RECORDS[name].modeling_view().union_bonds == recorded


@mark.parametrize('name', sorted(FROZEN))
def test_what_could_not_be_placed_is_reported_as_it_was_recorded(name):
    view = RECORDS[name].modeling_view()
    assert view.unmapped == FROZEN[name]['unmapped']
    assert view.collisions == {side: tuple(numbers)
                               for side, numbers in FROZEN[name]['collisions'].items()}


def test_the_corpus_reaches_the_three_branches_it_was_written_for():
    """A fixture whose sparse records came out empty for the wrong reason pins nothing."""
    assert FROZEN['unmapped_hydrogenation']['states'] == [], \
        'an unmapped reaction must contribute no union atom'
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
