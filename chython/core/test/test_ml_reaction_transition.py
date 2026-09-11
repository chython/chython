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
"""`ReactionContainer.transition_view`: the union of two sides, and what each side reports.

The state rules are `ReactionModelingView`'s, because task 9 rebuilds that view on this kernel.
"""
from pytest import raises

from chython.core import MoleculeContainer, TensorEncoding, read_reaction_smiles

from .modeling_view_corpus import RECORDS


def _by_map(view):
    """Union columns keyed by map number, for a record whose map numbers are unique."""
    return {int(mn): (int(view.elements[i]), int(view.h_before[i]), int(view.n_before[i]),
                      int(view.h_after[i]), int(view.n_after[i]))
            for i, mn in enumerate(view.map_numbers.tolist())}


def _bonds_by_map(view):
    mapping = view.map_numbers.tolist()
    return {tuple(sorted((int(mapping[i]), int(mapping[j])))): (int(b), int(a))
            for (i, j), b, a in zip(view.bonds.tolist(), view.bond_before, view.bond_after)}


def test_a_fully_mapped_substitution_reports_both_sides_of_every_atom():
    view = read_reaction_smiles('[CH3:1][Br:2].[OH-:3]>>[CH3:1][OH:3].[Br-:2]').transition_view()
    assert _by_map(view) == {
        1: (6, 3, 1, 3, 1),      # carbon: three hydrogens and one heavy neighbour on each side
        2: (35, 0, 1, 0, 0),     # bromine: bonded before, a free ion after
        3: (8, 1, 0, 1, 1),      # oxygen: a free hydroxide before, bonded after
    }


def test_a_broken_bond_and_a_formed_one_are_the_reaction_centre():
    view = read_reaction_smiles('[CH3:1][Br:2].[OH-:3]>>[CH3:1][OH:3].[Br-:2]').transition_view()
    assert _bonds_by_map(view) == {(1, 2): (1, 0), (1, 3): (0, 1)}


def test_an_order_change_keeps_one_bond_with_two_orders():
    view = read_reaction_smiles('[CH2:1]=[CH2:2]>>[CH3:1][CH3:2]').transition_view()
    assert _bonds_by_map(view) == {(1, 2): (2, 1)}


def test_a_broken_bond_changes_the_degree_on_one_side_only():
    """The mutation task 7 could not catch: n_before and n_after must come from their own sides."""
    view = read_reaction_smiles('[CH3:1][CH3:2]>>[CH3:1].[CH3:2]').transition_view()
    states = _by_map(view)
    assert states[1][2] == 1, 'one heavy neighbour before'
    assert states[1][4] == 0, 'none after'


def test_the_union_order_is_reactant_atoms_then_product_only_atoms():
    """Components stay contiguous, which is what a distance block over the union needs."""
    view = read_reaction_smiles('[CH3:1][C:2](=[O:3])[OH:9]'
                                '>>[CH3:1][C:2](=[O:3])[O:4][CH3:5].[OH2:9]').transition_view()
    assert view.map_numbers.tolist() == [1, 2, 3, 9, 4, 5]


def test_a_reactant_only_fragment_counts_only_its_own_neighbours_after():
    """A leaving group's after-degree is its degree within the leaving fragment.

    Atom 4 (ester oxygen) is reactant-only: its bonds to atom 2 (retained carbonyl C) break, so
    n_after = 1 (only atom 5, the tert-butyl C, is beside it in the leaving fragment).  Atom 5
    (tert-butyl C) has 4 reactant-only neighbours -- 4, 6, 7, 8 -- ALL of which left with it,
    so n_after = 4.  The ester oxygen (4) left WITH atom 5, not behind it.
    """
    view = read_reaction_smiles(
        '[CH3:1][C:2](=[O:3])[O:4][C:5]([CH3:6])([CH3:7])[CH3:8].[OH2:9]'
        '>>[CH3:1][C:2](=[O:3])[OH:9]').transition_view()
    states = _by_map(view)
    assert states[4][2] == 2, 'the ester oxygen has two heavy neighbours before'
    assert states[4][4] == 1, 'after, only the tert-butyl carbon is still beside it'
    assert states[5][2] == 4
    assert states[5][4] == 4, 'all four reactant-only neighbours (including the ester oxygen) left with it'


def test_a_product_only_fragment_is_the_mirror():
    """Atom 5 (tert-butyl C) has 4 product-only neighbours, so n_before = 4."""
    view = read_reaction_smiles(
        '[CH3:1][C:2](=[O:3])[OH:9]'
        '>>[CH3:1][C:2](=[O:3])[O:4][C:5]([CH3:6])([CH3:7])[CH3:8].[OH2:9]').transition_view()
    states = _by_map(view)
    assert states[4][4] == 2
    assert states[4][2] == 1
    assert states[5][4] == 4
    assert states[5][2] == 4


def test_an_agent_contributes_no_atom_and_no_bond():
    """Agents are excluded by not being handed to the kernel, not by a filter inside it."""
    view = read_reaction_smiles('[CH3:1][C:2](=[O:3])[OH:4].[NH2:5][CH3:6]'
                                '>[CH3:10][CH2:11][OH:12]'
                                '>[CH3:1][C:2](=[O:3])[NH:5][CH3:6].[OH2:4]').transition_view()
    assert sorted(view.map_numbers.tolist()) == [1, 2, 3, 4, 5, 6]


def test_an_unmapped_atom_is_counted_and_left_out_of_the_union():
    view = read_reaction_smiles('[CH3:1][C:2](=[O:3])[OH:4].CO'
                                '>>[CH3:1][C:2](=[O:3])[O:4]C.O').transition_view()
    assert sorted(view.map_numbers.tolist()) == [1, 2, 3, 4]
    assert view.unmapped == {'reactants': 2, 'products': 2}


def test_a_colliding_map_number_is_reported_and_not_refused():
    view = read_reaction_smiles('[CH3:1][CH3:1]>>[CH3:2][CH3:3]').transition_view()
    assert view.collisions == {'reactants': (1,), 'products': ()}
    assert view.map_numbers.tolist() == [1, 2, 3]


def test_a_bond_whose_endpoint_is_unmapped_is_not_a_union_bond():
    """Half a bond cannot be placed in the union, and inventing an endpoint would be worse."""
    view = read_reaction_smiles('[CH3:1]CO>>[CH3:1]C=O').transition_view()
    assert view.bonds.shape == (0, 2)
    assert view.map_numbers.tolist() == [1]


def test_a_reaction_with_no_mapping_at_all_gives_an_empty_union():
    view = read_reaction_smiles('C=C.[H][H]>>CC').transition_view()
    assert view.elements.shape == (0,)
    assert view.bonds.shape == (0, 2)
    assert view.distances.shape == (0, 0)
    assert view.unmapped == {'reactants': 4, 'products': 2}


def test_an_empty_product_side_is_a_record():
    """Both atoms are reactant-only; they leave together, so each counts the other after."""
    view = read_reaction_smiles('[CH3:1][CH3:2]>>').transition_view()
    assert view.map_numbers.tolist() == [1, 2]
    assert _by_map(view)[1] == (6, 3, 1, 3, 1), \
        'both atoms leave together; each counts the other as its reactant-only neighbour'


def test_the_distance_block_spans_the_union_and_not_one_side():
    """A bond formed on the product side shortens a path; that is what the union graph is for."""
    view = read_reaction_smiles('[CH3:1][Br:2].[OH-:3]>>[CH3:1][OH:3].[Br-:2]').transition_view()
    mapping = view.map_numbers.tolist()
    i, j = mapping.index(2), mapping.index(3)
    assert view.distances[i][j] == 2, 'bromine to oxygen, through the carbon, in the union'


def test_a_map_number_at_the_top_of_its_range_costs_nothing_extra():
    """The table is sized by the largest map number, and only touched slots are initialized."""
    view = read_reaction_smiles('[CH3:9999][CH3:2]>>[CH3:9999][CH3:2]').transition_view()
    assert sorted(view.map_numbers.tolist()) == [2, 9999]


def test_the_encoding_reaches_the_reaction_path_too():
    enc = TensorEncoding(element_shift=2, hydrogen_shift=1, neighbor_shift=2, distance_shift=2,
                         disconnected=1, width=8, pad=0, pad_diagonal=1)
    view = read_reaction_smiles('[CH3:1][Br:2].[OH-:3]>>[CH3:1][OH:3].[Br-:2]').transition_view(enc)
    assert view.elements.tolist() == [8, 37, 10, 0, 0, 0, 0, 0]
    assert view.distances.shape == (8, 8)
    assert view.distances[7].tolist() == [0, 0, 0, 0, 0, 0, 0, 1]


def test_a_union_larger_than_width_raises_and_names_its_atom_count():
    rxn = read_reaction_smiles('[CH3:1][C:2](=[O:3])[OH:4].[CH3:5][OH:6]'
                               '>>[CH3:1][C:2](=[O:3])[O:4][CH3:5].[OH2:6]')
    with raises(ValueError, match='6 atoms'):
        rxn.transition_view(TensorEncoding(width=4))


def test_a_collided_atom_contributes_no_bond_and_no_degree():
    """The second claim on a map number contributes its entry in `collisions` and nothing else.

    The middle atom repeats map number 1.  Both of its bonds have it as an endpoint, so neither can be
    placed in the union -- inventing a bond between map numbers 1 and 3, which no side draws, would be
    worse than leaving the record's two atoms unconnected.
    """
    view = read_reaction_smiles('[CH3:1][CH2:1][CH3:3]>>[CH4:1].[CH4:3]').transition_view()
    assert view.collisions == {'reactants': (1,), 'products': ()}
    assert view.map_numbers.tolist() == [1, 3]
    assert view.bonds.shape == (0, 2)
    assert view.n_before.tolist() == [0, 0]


def test_the_bond_arrays_are_sized_by_the_merged_count_not_the_capacity():
    """A merged bond consumes one union slot, and the arrays must be exactly that long."""
    view = read_reaction_smiles('[CH2:1]=[CH2:2]>>[CH3:1][CH3:2]').transition_view()
    assert view.bonds.shape == (1, 2)
    assert view.bond_before.shape == (1,) and view.bond_after.shape == (1,)


def test_every_record_in_the_corpus_produces_a_view():
    """The corpus is the frozen-fixture corpus, so task 9 has a working kernel under every record."""
    for name, rxn in RECORDS.items():
        view = rxn.transition_view()
        n = view.elements.shape[0]
        assert view.distances.shape == (n, n), name
        assert view.bonds.shape[0] == view.bond_before.shape[0], name
        assert view.map_numbers.shape[0] == n, name


def test_a_leaving_atom_below_its_staying_neighbour_still_counts_only_its_own_side():
    """Map 4 (ester oxygen, union index 3, side=1) bonds map 2 (carbonyl C, index 4, side=3).

    Bond n=3 < m=4 has a reactant-only atom at the lower index: side[n]==1 is True, but
    side[n]==1 and side[m]==1 is False -- so n_after[3] must remain 1 (only map-5 neighbour).
    A mutation that weakens the guard to `side[n]==1` alone inflates n_after[3] to 2.
    """
    view = read_reaction_smiles(
        '[CH3:5]([CH3:6])([CH3:7])[O:4][C:2](=[O:3])[CH3:1].[OH2:9]'
        '>>[CH3:1][C:2](=[O:3])[OH:9]').transition_view()
    assert view.map_numbers.tolist() == [5, 6, 7, 4, 2, 3, 1, 9]
    assert view.n_before.tolist() == [3, 1, 1, 2, 3, 1, 1, 0]
    assert view.n_after.tolist() == [3, 1, 1, 1, 3, 1, 1, 1]
    states = _by_map(view)
    assert states[4][4] == 1, 'ester oxygen leaves without carbonyl C; n_after is 1, not 2'


def test_a_collision_on_the_product_side_is_reported_as_its_own_side():
    """A map number claimed twice in the product is listed in collisions[products], not reactants."""
    view = read_reaction_smiles('[CH3:1][CH3:2]>>[CH3:3][CH3:3]').transition_view()
    assert view.collisions == {'reactants': (), 'products': (3,)}
    assert view.map_numbers.tolist() == [1, 2, 3]
    assert _bonds_by_map(view) == {(1, 2): (1, 0)}, 'the product bond has a rejected endpoint'


def test_a_rejected_product_atom_is_not_appended_as_a_second_union_row():
    """A map number claimed twice in the product contributes one row, not two.

    [CH3:2] takes the union slot; the colliding [CH2:2] is counted in `collisions` and nothing else.
    """
    view = read_reaction_smiles('[CH3:1]>>[CH3:2][CH2:2]').transition_view()
    assert view.collisions == {'reactants': (), 'products': (2,)}
    assert view.map_numbers.tolist() == [1, 2]
    assert _by_map(view)[2][3] == 3, 'h_after comes from the first accepted product atom'


def test_a_rejected_product_atom_does_not_overwrite_the_row_that_kept_the_claim():
    """Map 1 is on both sides, so its union row already exists when the duplicate claim arrives.

    `[CH4:1]` sets h_after[0] = 4; the colliding `[CH3:1]` must not pull it down to 3.  This is the
    branch `pseen` protects that a product-only collision cannot reach -- there, the rejected atom
    would be appended rather than overwrite an accepted row.
    """
    view = read_reaction_smiles('[CH4:1]>>[CH4:1].[CH3:1][CH3:2]').transition_view()
    assert view.collisions == {'reactants': (), 'products': (1,)}
    assert view.map_numbers.tolist() == [1, 2]
    assert view.h_after.tolist() == [4, 3], 'the rejected CH3 must not pull row 0 down to 3'


def test_transition_view_refuses_a_reactant_inside_an_open_edit_scope():
    """The arena still holds the pre-scope state; a view from it is from a molecule that may no longer exist."""
    rxn = read_reaction_smiles('[CH3:1][Br:2].[OH-:3]>>[CH3:1][OH:3].[Br-:2]')
    reactant = list(rxn.reactants)[0]
    with raises(RuntimeError, match='pending edits'):
        with reactant.edit() as e:
            e.delete_atom(next(iter(reactant.atoms())).n)
            rxn.transition_view()


def test_modeling_view_refuses_a_reactant_inside_an_open_edit_scope():
    """modeling_view is a dict assembly over the same kernel; the guard must cover it too."""
    rxn = read_reaction_smiles('[CH3:1][Br:2].[OH-:3]>>[CH3:1][OH:3].[Br-:2]')
    reactant = list(rxn.reactants)[0]
    with raises(RuntimeError, match='pending edits'):
        with reactant.edit() as e:
            e.delete_atom(next(iter(reactant.atoms())).n)
            rxn.modeling_view()


def test_transition_view_refuses_a_product_inside_an_open_edit_scope():
    """The guard walks both sides; one that walked only the reactants would pass the two tests above."""
    rxn = read_reaction_smiles('[CH3:1][Br:2].[OH-:3]>>[CH3:1][OH:3].[Br-:2]')
    product = list(rxn.products)[0]
    with raises(RuntimeError, match='pending edits'):
        with product.edit() as e:
            e.delete_atom(next(iter(product.atoms())).n)
            rxn.transition_view()
