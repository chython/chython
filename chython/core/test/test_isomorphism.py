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
import pytest
from chython.core import MoleculeContainer, QueryContainer


def mol(elements, bonds):
    """A molecule from a list of atomic numbers and (i, j, order) triples, 0-based."""
    m = MoleculeContainer()
    ids = [m.add_atom(e) for e in elements]
    for i, j, order in bonds:
        m.add_bond(ids[i], ids[j], order)
    return m, ids


def query(elements, bonds):
    """A query from the same shape: exact elements, exact bond orders, nothing else."""
    q = QueryContainer()
    ids = []
    for e in elements:
        sid = q.add_atom()
        q.atom_primitive(sid, 'element', e)
        ids.append(sid)
    for i, j, order in bonds:
        q.add_bond(ids[i], ids[j])
        q.bond_primitive(ids[i], ids[j], 'bond_order', order)
    return q, ids


def ring_query(orders):
    """A carbon ring query of len(orders) atoms whose i-th bond carries orders[i]."""
    q = QueryContainer()
    ids = [q.add_atom() for _ in orders]
    for sid in ids:
        q.atom_primitive(sid, 'element', 6)
    for i, order in enumerate(orders):
        q.add_bond(ids[i], ids[(i + 1) % len(orders)])
        q.bond_primitive(ids[i], ids[(i + 1) % len(orders)], 'bond_order', order)
    return q


def one_atom(*terms):
    """A one-atom query whose primitives are ANDed: one_atom(('element', 6), ('ring_size', 5))."""
    q = QueryContainer()
    sid = q.add_atom()
    for n, (name, value) in enumerate(terms):
        if n:
            q.atom_operator(sid, 'and_low')
        q.atom_primitive(sid, name, value)
    return q, sid


def spiro45decane():
    """Spiro[4.5]decane: a 5-ring and a 6-ring sharing exactly atom 0.

    Its relevant cycles are the two rings themselves -- their sum is a figure-eight, not a
    cycle -- so atom 0 is the only atom in both a 5- and a 6-membered ring.
    """
    return mol([6] * 10,
               [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1), (5, 0, 1),
                (0, 6, 1), (6, 7, 1), (7, 8, 1), (8, 9, 1), (9, 0, 1)])


def test_a_single_atom_query_finds_every_matching_atom():
    m, _ = mol([6, 6, 8], [(0, 1, 1), (1, 2, 1)])
    q, _ = query([6], [])
    assert q.count(m) == 2


def test_a_single_atom_query_finds_nothing_when_the_element_is_absent():
    m, _ = mol([6, 6], [(0, 1, 1)])
    q, _ = query([7], [])
    assert q.count(m) == 0
    assert not q.is_substructure(m)


def test_a_two_atom_query_matches_both_directions_of_a_symmetric_bond():
    m, _ = mol([6, 6], [(0, 1, 1)])
    q, _ = query([6, 6], [(0, 1, 1)])
    assert q.count(m) == 2, 'C-C matches ethane twice, once per orientation'


def test_a_directed_query_matches_once():
    m, _ = mol([6, 8], [(0, 1, 1)])
    q, _ = query([6, 8], [(0, 1, 1)])
    assert q.count(m) == 1


def test_bond_order_is_enforced():
    m, _ = mol([6, 6], [(0, 1, 1)])
    q, _ = query([6, 6], [(0, 1, 2)])
    assert q.count(m) == 0


def test_a_chain_query_matches_inside_a_longer_chain():
    m, _ = mol([6, 6, 6, 6], [(0, 1, 1), (1, 2, 1), (2, 3, 1)])
    q, _ = query([6, 6, 6], [(0, 1, 1), (1, 2, 1)])
    assert q.count(m) == 4, 'two positions, two orientations'


def test_injectivity_forbids_reusing_an_atom():
    # a 3-atom path query must not match a 2-atom molecule by walking back
    m, _ = mol([6, 6], [(0, 1, 1)])
    q, _ = query([6, 6, 6], [(0, 1, 1), (1, 2, 1)])
    assert q.count(m) == 0


def test_a_branch_query_needs_a_branching_atom():
    m, _ = mol([6, 6, 6, 6], [(0, 1, 1), (1, 2, 1), (1, 3, 1)])
    q, _ = query([6, 6, 6], [(0, 1, 1), (1, 2, 1)])
    assert q.count(m) == 6, 'three ways to pick the centre pair, doubled by orientation'
    q2, ids = query([6, 6, 6, 6], [(0, 1, 1), (1, 2, 1), (1, 3, 1)])
    assert q2.count(m) == 6, 'the three leaves permute'


def test_get_mapping_yields_stable_id_dicts():
    m, mids = mol([6, 8], [(0, 1, 1)])
    q, qids = query([6, 8], [(0, 1, 1)])
    mappings = list(q.get_mapping(m))
    assert mappings == [{qids[0]: mids[0], qids[1]: mids[1]}]


def test_get_raw_mapping_yields_one_tuple_per_embedding():
    # Deliberately NOT the ordering test: for C-O position 0 is the carbon either way, so this
    # assertion cannot tell DFS order from declaration order.  That is the test below.
    m, mids = mol([6, 8], [(0, 1, 1)])
    q, _ = query([6, 8], [(0, 1, 1)])
    assert list(q.get_raw_mapping(m)) == [(mids[0], mids[1])]


def test_raw_mapping_slots_are_in_dfs_order_not_declaration_order():
    """The one case where the two orders differ, plus the accessor that makes slots readable.

    Fails against a kernel that walks positions in declaration order, and against a
    query_numbers that returns declaration order (which would make every get_raw_mapping tuple
    silently mis-keyed for exactly the queries where it matters).
    """
    m, mids = mol([6, 6, 8, 6], [(0, 1, 1), (1, 2, 1), (2, 3, 1)])
    q, qids = query([6, 8, 6], [(0, 1, 1), (1, 2, 1)])
    # declared carbon-first, but the rare oxygen roots the DFS, so slot 0 is the OXYGEN
    assert q.query_numbers() == (qids[1], qids[0], qids[2])
    assert sorted(q.get_raw_mapping(m)) == [(mids[2], mids[1], mids[3]),
                                            (mids[2], mids[3], mids[1])]
    # and the accessor is the key that interprets those tuples: slot 0 really is the query oxygen,
    # which really did match the molecule's only oxygen
    slots = q.query_numbers()
    for row in q.get_raw_mapping(m):
        assert dict(zip(slots, row))[qids[1]] == mids[2]


def test_le_is_is_substructure():
    m, _ = mol([6, 8], [(0, 1, 1)])
    q, _ = query([6, 8], [(0, 1, 1)])
    other, _ = mol([6, 6], [(0, 1, 1)])
    assert q <= m
    assert not q <= other
    assert q.__le__('not a molecule') is NotImplemented


def test_get_mapping_over_a_disconnected_query():
    # two independent carbons, unconstrained relative to each other
    m, _ = mol([6, 6], [(0, 1, 1)])
    q = QueryContainer()
    for _ in range(2):
        sid = q.add_atom()
        q.atom_primitive(sid, 'element', 6)
    assert q.count(m) == 2, 'the two query atoms map to the two molecule atoms, both ways'


def test_a_disconnected_query_still_forbids_sharing_an_atom():
    m, _ = mol([6], [])
    q = QueryContainer()
    for _ in range(2):
        sid = q.add_atom()
        q.atom_primitive(sid, 'element', 6)
    assert q.count(m) == 0


def test_the_rare_element_root_does_not_change_the_answer():
    """The same graph, two declaration orders, one answer.

    The oxygen roots the DFS whichever way the query is spelled, so both spellings must agree.
    Fails against any kernel whose seed depends on declaration order rather than on the arena's
    chosen root -- a single spelling cannot see that, which is why both are built here.
    """
    m, _ = mol([6, 6, 8, 6], [(0, 1, 1), (1, 2, 1), (2, 3, 1)])
    carbon_first, _ = query([6, 8, 6], [(0, 1, 1), (1, 2, 1)])
    oxygen_first, _ = query([8, 6, 6], [(0, 1, 1), (0, 2, 1)])
    assert carbon_first.count(m) == 2
    assert oxygen_first.count(m) == carbon_first.count(m)


def test_a_ring_closure_query_now_matches():
    q, ids = query([6, 6, 6], [(0, 1, 1), (1, 2, 1)])
    q.add_bond(ids[2], ids[0])
    q.bond_primitive(ids[2], ids[0], 'bond_order', 1)
    m, _ = mol([6, 6, 6], [(0, 1, 1), (1, 2, 1), (2, 0, 1)])
    assert q.count(m) == 6


def test_matching_an_empty_molecule_finds_nothing():
    m = MoleculeContainer()
    q, _ = query([6], [])
    assert q.count(m) == 0
    assert not q.is_substructure(m)


def test_is_substructure_agrees_with_count_on_a_large_molecule():
    """'Stops at the first hit' is true (_query_container.pxi calls matcher_next once) but
    unassertable from Python without a step counter, and an assertion of it holds either way.

    What is checkable is that the one-shot path and the exhaustive path never disagree, in both
    directions -- an early-exit search that returned True on an empty candidate scan would pass a
    True-only test.
    """
    m, _ = mol([6] * 20, [(i, i + 1, 1) for i in range(19)])
    present, _ = query([6, 6], [(0, 1, 1)])
    absent, _ = query([7, 6], [(0, 1, 1)])
    assert present.is_substructure(m) is True and present.count(m) == 38
    assert absent.is_substructure(m) is False and absent.count(m) == 0


# ---------------------------------------------------------------------------------------------
# The half-edge word.  A non-root position must be tested against edge_words[k], never against
# the candidate's aggregate feature word 0 -- the aggregate ORs every incident bond's order bit
# together, so a folded box that forbids single rejects a carbonyl carbon that carries a methyl.
# ---------------------------------------------------------------------------------------------


def test_a_non_root_position_tests_the_half_edge_and_not_the_aggregate_word():
    # acetone; O=C-C must find both methyls.  Testing the carbonyl carbon's aggregate word 0
    # against the folded box (which forbids the single-bond bit) yields zero embeddings.
    m, _ = mol([6, 6, 8, 6], [(0, 1, 1), (1, 2, 2), (1, 3, 1)])
    q, _ = query([8, 6, 6], [(0, 1, 2), (1, 2, 1)])
    assert q.count(m) == 2


def test_the_half_edge_word_walks_an_asymmetric_ketone():
    # butan-2-one: O=C-C-C fits the ethyl side only, so the count is odd and orientation-free
    m, _ = mol([6, 6, 8, 6, 6], [(0, 1, 1), (1, 2, 2), (1, 3, 1), (3, 4, 1)])
    q, _ = query([8, 6, 6, 6], [(0, 1, 2), (1, 2, 1), (2, 3, 1)])
    assert q.count(m) == 1


def test_a_double_bond_query_does_not_match_the_single_bonded_neighbour():
    # the same acetone: O=C=C is impossible, and C-C=O picks the carbonyl bond only
    m, _ = mol([6, 6, 8, 6], [(0, 1, 1), (1, 2, 2), (1, 3, 1)])
    assert query([6, 6, 8], [(0, 1, 1), (1, 2, 2)])[0].count(m) == 2
    assert query([6, 6, 8], [(0, 1, 2), (1, 2, 2)])[0].count(m) == 0


def test_a_second_component_root_gets_the_aggregate_word_not_a_half_edge():
    """The aggregate/half-edge choice is root-ness, NOT `depth == 0` -- and only a second
    component can tell those apart.

    Position 1 here is a root (it must get the aggregate feature word) whose depth is not 0.  Key
    that choice on `position == 0` or `depth == 0` and every single-component test in this file
    still passes, while this count silently drops from 2 to 1: position 1's carbon-bucket cursor
    gets read out of edge_words instead, and one of those half-edges points at the oxygen, whose
    element bits the [C] box forbids.  The molecule needs a heteroatom for exactly that reason --
    over C-C both spurious half-edge words carry carbon bits and the bug hides.
    """
    m, mids = mol([8, 6, 6], [(0, 1, 2), (1, 2, 1)])          # O=C-C
    q = QueryContainer()
    o = q.add_atom()
    q.atom_primitive(o, 'element', 8)
    c = q.add_atom()
    q.atom_primitive(c, 'element', 6)                          # a second component: [O].[C]
    assert q.count(m) == 2
    # assert WHICH atoms, so a kernel that seeds the oxygen root from the wrong bucket also fails
    assert {frozenset(d.items()) for d in q.get_mapping(m)} == {
        frozenset({(o, mids[0]), (c, mids[1])}),
        frozenset({(o, mids[0]), (c, mids[2])}),
    }


# ---------------------------------------------------------------------------------------------
# The box disjunction and the any-list quantifier
# ---------------------------------------------------------------------------------------------


def test_a_root_disjunction_admits_every_box_not_just_the_first():
    m, _ = mol([6, 7, 8], [(0, 1, 1), (1, 2, 1)])
    q = QueryContainer()
    sid = q.add_atom()
    q.atom_primitive(sid, 'element', 6)
    q.atom_operator(sid, 'or')
    q.atom_primitive(sid, 'element', 7)
    assert q.count(m) == 2, 'carbon from the first box, nitrogen from the second, never oxygen'


def test_a_folded_bond_keeps_every_box_of_a_disjunction():
    # C-[C,N] over C-C-N: the terminal carbon has one carbon neighbour, the middle carbon has
    # a carbon and a nitrogen.  Testing only the first box gives 2, only the second gives 1.
    m, _ = mol([6, 6, 7], [(0, 1, 1), (1, 2, 1)])
    q = QueryContainer()
    a = q.add_atom()
    q.atom_primitive(a, 'element', 6)
    b = q.add_atom()
    q.atom_primitive(b, 'element', 6)
    q.atom_operator(b, 'or')
    q.atom_primitive(b, 'element', 7)
    q.add_bond(a, b)
    q.bond_primitive(a, b, 'bond_order', 1)
    assert q.count(m) == 3


def test_every_any_entry_must_be_satisfied():
    m, _ = spiro45decane()
    assert one_atom(('element', 6), ('ring_size', 5))[0].count(m) == 5
    assert one_atom(('element', 6), ('ring_size', 6))[0].count(m) == 6
    assert one_atom(('element', 6), ('ring_size', 5), ('ring_size', 6))[0].count(m) == 1, \
        'only the spiro atom lies on both rings'


def test_an_any_entry_no_atom_satisfies_rejects_everything():
    m, _ = mol([6] * 6, [(i, i + 1, 1) for i in range(5)])   # hexane: acyclic
    assert one_atom(('element', 6), ('ring_size', 6))[0].count(m) == 0


# ---------------------------------------------------------------------------------------------
# Injectivity where the atom-count screen cannot short-circuit it
# ---------------------------------------------------------------------------------------------


def test_two_children_of_one_parent_may_not_share_an_atom():
    # C-C plus an isolated carbon: three atoms, so the size screen lets the search run.  The
    # 3-atom path roots at its middle atom, whose two children would both take atom 1 without
    # the injectivity check.
    m, _ = mol([6, 6, 6], [(0, 1, 1)])
    q, _ = query([6, 6, 6], [(0, 1, 1), (1, 2, 1)])
    assert q.count(m) == 0


def test_two_component_roots_may_not_share_the_only_matching_atom():
    m, _ = mol([6, 8], [(0, 1, 1)])
    q = QueryContainer()
    for _ in range(2):
        sid = q.add_atom()
        q.atom_primitive(sid, 'element', 6)
    assert q.count(m) == 0, 'two roots, one carbon'


# ---------------------------------------------------------------------------------------------
# Heavy elements: the root bucket and the word-1 identity test
# ---------------------------------------------------------------------------------------------


def test_a_heavy_element_root_seeds_from_its_own_bucket():
    m, _ = mol([79, 6, 78], [(0, 1, 1), (1, 2, 1)])      # Au-C-Pt
    assert query([79], [])[0].count(m) == 1
    assert query([78], [])[0].count(m) == 1
    assert query([80], [])[0].count(m) == 0
    assert query([79, 6, 78], [(0, 1, 1), (1, 2, 1)])[0].count(m) == 1


# ---------------------------------------------------------------------------------------------
# Guards and lifetime
# ---------------------------------------------------------------------------------------------


def test_a_query_with_a_component_group_is_matched():
    # Task 13: grouped queries are now matched.  A connected query with one group just matches.
    q, ids = query([6, 6], [(0, 1, 1)])
    q.set_group(ids[0], 0)
    m, _ = mol([6, 6], [(0, 1, 1)])
    assert q.count(m) == 2


def test_mutating_the_query_reseals_before_the_next_match():
    m, _ = mol([6, 6, 8], [(0, 1, 1), (1, 2, 1)])
    q, ids = query([6], [])
    assert q.count(m) == 2
    extra = q.add_atom()
    q.atom_primitive(extra, 'element', 8)
    q.add_bond(ids[0], extra)
    q.bond_primitive(ids[0], extra, 'bond_order', 1)
    assert q.count(m) == 1, 'the second seal sees C-O, not C'


def test_get_mapping_survives_a_molecule_edit_mid_iteration():
    m, mids = mol([6, 6], [(0, 1, 1)])
    q, _ = query([6, 6], [(0, 1, 1)])
    it = q.get_mapping(m)
    first = next(it)
    m.delete_atom(mids[1])
    rest = list(it)
    assert len(rest) == 1, 'the generator finishes against the arena it started on'
    assert first != rest[0]


def test_abandoning_a_generator_mid_iteration_leaves_later_searches_intact():
    """Renamed from a 'frees its matcher' claim this body cannot check.

    The leak has no observable consequence in Python: every generator owns its own matcher_t, so
    dropping `finally: matcher_free` costs memory and changes no result.  Asserting on RSS would be
    a flaky test, so the leak was checked out of band instead (200_000 abandoned generators moved
    RSS by 0 bytes) and this pins what IS observable -- that abandoning a partly-consumed search
    corrupts neither the query nor any later search.  Fails against any implementation that cached
    matcher state on the container rather than per generator.
    """
    m, _ = mol([6] * 8, [(i, i + 1, 1) for i in range(7)])
    q, _ = query([6, 6], [(0, 1, 1)])
    for _ in range(200):
        it = q.get_mapping(m)
        next(it)
        del it
    assert q.count(m) == 14, '7 bonds, each matched in both directions'
    assert len(list(q.get_mapping(m))) == 14, 'a full generator still runs to exhaustion'


# ---------------------------------------------------------------------------------------------
# Ring closures (Task 11).  A closure is a query bond both of whose endpoints are already
# mapped when the DFS reaches it; the kernel must find the molecule half-edge and test it.
# ---------------------------------------------------------------------------------------------


def test_a_ring_query_matches_a_ring():
    m, _ = mol([6, 6, 6], [(0, 1, 1), (1, 2, 1), (2, 0, 1)])
    q, ids = query([6, 6, 6], [(0, 1, 1), (1, 2, 1)])
    q.add_bond(ids[2], ids[0])
    q.bond_primitive(ids[2], ids[0], 'bond_order', 1)
    assert q.count(m) == 6, 'cyclopropane has six automorphisms'


def test_a_ring_query_does_not_match_a_chain():
    m, _ = mol([6, 6, 6], [(0, 1, 1), (1, 2, 1)])
    q, ids = query([6, 6, 6], [(0, 1, 1), (1, 2, 1)])
    q.add_bond(ids[2], ids[0])
    q.bond_primitive(ids[2], ids[0], 'bond_order', 1)
    assert q.count(m) == 0


def test_a_chain_query_does_match_a_ring():
    m, _ = mol([6, 6, 6], [(0, 1, 1), (1, 2, 1), (2, 0, 1)])
    q, _ = query([6, 6, 6], [(0, 1, 1), (1, 2, 1)])
    assert q.count(m) == 6, 'a path query does not forbid the extra bond'


def test_a_closure_enforces_its_bond_order():
    m, _ = mol([6, 6, 6], [(0, 1, 1), (1, 2, 1), (2, 0, 1)])
    q, ids = query([6, 6, 6], [(0, 1, 1), (1, 2, 1)])
    q.add_bond(ids[2], ids[0])
    q.bond_primitive(ids[2], ids[0], 'bond_order', 2)
    assert q.count(m) == 0


def test_a_six_ring_query_matches_benzene_kekule():
    orders = [2, 1, 2, 1, 2, 1]
    m = MoleculeContainer()
    ids = [m.add_atom(6) for _ in range(6)]
    for i, o in enumerate(orders):
        m.add_bond(ids[i], ids[(i + 1) % 6], o)
    q = QueryContainer()
    qids = [q.add_atom() for _ in range(6)]
    for i in range(6):
        q.atom_primitive(qids[i], 'element', 6)
    for i, o in enumerate(orders):
        q.add_bond(qids[i], qids[(i + 1) % 6])
        q.bond_primitive(qids[i], qids[(i + 1) % 6], 'bond_order', o)
    # the alternating labels kill the odd rotations and the vertex-centred reflections, leaving
    # three rotations and three reflections through bond midpoints
    assert q.count(m) == 6


def test_two_closures_on_one_atom():
    # bicyclo[1.1.0]: 4 atoms, 5 bonds -- three tree bonds, two closures
    m, _ = mol([6, 6, 6, 6], [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 0, 1), (0, 2, 1)])
    q, ids = query([6, 6, 6, 6], [(0, 1, 1), (1, 2, 1), (2, 3, 1)])
    q.add_bond(ids[3], ids[0])
    q.bond_primitive(ids[3], ids[0], 'bond_order', 1)
    q.add_bond(ids[0], ids[2])
    q.bond_primitive(ids[0], ids[2], 'bond_order', 1)
    assert q.count(m) == 4, 'the two bridge atoms swap, as do the two apex atoms'


def test_a_tree_query_has_no_closures():
    # a check on seal, not the kernel: every bond of an acyclic query is a tree bond
    q, _ = query([6, 6, 6, 6], [(0, 1, 1), (1, 2, 1), (1, 3, 1)])
    assert q.closure_count() == 0


def test_a_closure_bond_disjunction_matches_via_second_box():
    # The closure bond has TWO boxes that survive boxes_merge because they differ in two
    # independent spans:
    #   box 0: acyclic AND single  (not-ring; neg forbids ring_plain, ring_arom, double, ...)
    #   box 1: in-ring AND double  (ring;    neg forbids not_ring, single, ...)
    # The molecule ring bond is in-ring and double, so box 0 REJECTS and box 1 ADMITS.
    # An implementation that reads only boxes[qb.box_begin] (box 0) returns count=0.
    # Pins the `for b in range(qb.box_count)` disjunction inside closures_admit.
    m, _ = mol([6, 6, 6, 6], [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 0, 2)])
    q, ids = query([6, 6, 6, 6], [(0, 1, 1), (1, 2, 1), (2, 3, 1)])
    q.add_bond(ids[3], ids[0])
    # box 0: NOT in ring AND single
    q.bond_primitive(ids[3], ids[0], 'bond_ring', 1, True)  # negated=True → not in ring
    q.bond_operator(ids[3], ids[0], 'and_high')
    q.bond_primitive(ids[3], ids[0], 'bond_order', 1)
    # OR
    q.bond_operator(ids[3], ids[0], 'or')
    # box 1: in ring AND double
    q.bond_primitive(ids[3], ids[0], 'bond_ring', 1)        # in ring
    q.bond_operator(ids[3], ids[0], 'and_high')
    q.bond_primitive(ids[3], ids[0], 'bond_order', 2)
    assert q.count(m) > 0


def test_two_closures_with_different_bond_orders_pin_bond_index():
    # Query: N-A-B-C chain (bslots 0-2, all single, tree bonds, box_count=0) plus three
    # back edges: C->N (bslot 3, single), C->A (bslot 4, double), B->N (bslot 5, single).
    # query_seal confirms: bslots 0,1,2 are tree bonds (box_count=0); bslots 3,4,5 are
    # closure bonds (box_count=1 each). The DFS order is N->A->B->C, so:
    #   pos2 (B) owns closure {to_index=0, bond_index=5}  (B->N, single)
    #   pos3 (C) owns closure {to_index=0, bond_index=3}  (C->N, single)
    #   pos3 (C) owns closure {to_index=1, bond_index=4}  (C->A, double)
    # Wrong implementation `qb = bonds + c` reads bonds[0] for every c=0 closure and
    # bonds[1] for c=1.  Both bonds[0] and bonds[1] are tree bonds (box_count=0), so the
    # disjunction loop never runs, ok stays False, and count=0.  Correct code reads
    # bonds[5], bonds[3], bonds[4] via closures[...].bond_index.  The C->A double and
    # C->N single constraints are different: a molecule with C->N(double)/C->A(single)
    # is rejected even by correct code, confirming bond_index is load-bearing.
    # N (element 7) is the DFS root, ensuring bslot 0 is always a tree bond.
    q, ids = query([7, 6, 6, 6], [(0, 1, 1), (1, 2, 1), (2, 3, 1)])
    q.add_bond(ids[3], ids[0])
    q.bond_primitive(ids[3], ids[0], 'bond_order', 1)
    q.add_bond(ids[3], ids[1])
    q.bond_primitive(ids[3], ids[1], 'bond_order', 2)
    q.add_bond(ids[2], ids[0])
    q.bond_primitive(ids[2], ids[0], 'bond_order', 1)
    m, _ = mol([7, 6, 6, 6], [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 0, 1), (3, 1, 2), (2, 0, 1)])
    assert q.count(m) > 0


# ---------------------------------------------------------------------------------------------
# The automorphism filter (Task 12).  A symmetric query reports one embedding per site rather
# than one per symmetry: the group is computed at seal, from the query's UNFOLDED atom terms and
# its per-edge terms, and the kernel keeps only the lexicographically smallest member of each
# orbit.
# ---------------------------------------------------------------------------------------------


def test_the_filter_collapses_a_symmetric_pattern():
    m, _ = mol([6, 6, 6], [(0, 1, 1), (1, 2, 1), (2, 0, 1)])
    q, ids = query([6, 6, 6], [(0, 1, 1), (1, 2, 1)])
    q.add_bond(ids[2], ids[0])
    q.bond_primitive(ids[2], ids[0], 'bond_order', 1)
    assert q.count(m) == 6
    assert q.count(m, automorphism_filter=True) == 1


def test_the_filter_leaves_an_asymmetric_pattern_alone():
    m, _ = mol([6, 8], [(0, 1, 1)])
    q, _ = query([6, 8], [(0, 1, 1)])
    assert q.count(m, automorphism_filter=True) == 1


def test_the_filter_keeps_distinct_sites():
    # C-C in butane: three distinct bonds, each counted once
    m, _ = mol([6, 6, 6, 6], [(0, 1, 1), (1, 2, 1), (2, 3, 1)])
    q, _ = query([6, 6], [(0, 1, 1)])
    assert q.count(m) == 6
    assert q.count(m, automorphism_filter=True) == 3


def test_the_filter_on_a_star_pattern():
    m, _ = mol([6, 6, 6, 6], [(0, 1, 1), (1, 2, 1), (1, 3, 1)])
    q, ids = query([6, 6, 6, 6], [(0, 1, 1), (1, 2, 1), (1, 3, 1)])
    assert q.count(m) == 6
    assert q.count(m, automorphism_filter=True) == 1, 'the three leaves are interchangeable'


def test_different_elements_break_the_symmetry():
    m, _ = mol([6, 6, 8], [(0, 1, 1), (1, 2, 1), (2, 0, 1)])
    q, ids = query([6, 6, 8], [(0, 1, 1), (1, 2, 1)])
    q.add_bond(ids[2], ids[0])
    q.bond_primitive(ids[2], ids[0], 'bond_order', 1)
    assert q.count(m, automorphism_filter=True) == 1
    assert q.count(m) == 2, 'only the two carbons swap'


def test_the_group_is_computed_once():
    q, ids = query([6, 6, 6], [(0, 1, 1), (1, 2, 1)])
    q.add_bond(ids[2], ids[0])
    q.bond_primitive(ids[2], ids[0], 'bond_order', 1)
    m, _ = mol([6, 6, 6], [(0, 1, 1), (1, 2, 1), (2, 0, 1)])
    q.count(m, automorphism_filter=True)
    assert q.automorphism_count() == 5, 'six automorphisms, identity excluded'
    q.count(m, automorphism_filter=True)
    assert q.automorphism_generation() == 1


def test_an_asymmetric_query_stores_no_permutations():
    q, _ = query([6, 8], [(0, 1, 1)])
    m, _ = mol([6, 8], [(0, 1, 1)])
    q.count(m, automorphism_filter=True)
    assert q.automorphism_count() == 0


def test_get_mapping_honours_the_filter():
    m, _ = mol([6, 6], [(0, 1, 1)])
    q, _ = query([6, 6], [(0, 1, 1)])
    assert len(list(q.get_mapping(m))) == 2
    assert len(list(q.get_mapping(m, automorphism_filter=True))) == 1


def test_get_raw_mapping_and_is_substructure_honour_the_filter():
    m, _ = mol([6, 6], [(0, 1, 1)])
    q, _ = query([6, 6], [(0, 1, 1)])
    assert len(list(q.get_raw_mapping(m))) == 2
    assert len(list(q.get_raw_mapping(m, automorphism_filter=True))) == 1
    # the filter never turns a hit into a miss: the orbit representative always survives
    assert q.is_substructure(m, automorphism_filter=True) is True
    absent, _ = query([7, 7], [(0, 1, 1)])
    assert absent.is_substructure(m, automorphism_filter=True) is False


def test_a_mutation_recomputes_the_group():
    """automorphism_generation resets with the seal it belongs to.

    Two cases on purpose: the first seal's group must be gone after the edit, and the second
    seal's must be the group of the EDITED query -- an implementation that kept the old rows
    would still report a generation of 1 while filtering by a stale permutation.
    """
    ethane, _ = mol([6, 6], [(0, 1, 1)])
    q, ids = query([6, 6], [(0, 1, 1)])
    assert q.count(ethane, automorphism_filter=True) == 1
    assert q.automorphism_generation() == 1
    assert q.automorphism_count() == 1
    extra = q.add_atom()
    q.atom_primitive(extra, 'element', 8)
    q.add_bond(ids[1], extra)
    q.bond_primitive(ids[1], extra, 'bond_order', 1)
    assert q.automorphism_generation() == 0, 'the edit dropped the group with the arena'
    methanol, _ = mol([6, 6, 8], [(0, 1, 1), (1, 2, 1)])
    assert q.count(methanol, automorphism_filter=True) == 1
    assert q.automorphism_count() == 0, 'C-C-O has no symmetry'
    assert q.automorphism_generation() == 1


def test_the_group_comes_from_the_atom_terms_not_the_folded_boxes():
    """query_seal folds each non-root position's tree bond into the atom it leads to, so a root's
    box set and its child's are never equal even when the two query atoms are identical.  A
    partition built from the SEALED boxes therefore puts every root in a class of its own and
    finds no symmetry at all -- C-C would report zero automorphisms and the filter would be a
    no-op.  Two cases: the symmetric query must collapse, the asymmetric one must not.
    """
    ethane, _ = mol([6, 6], [(0, 1, 1)])
    symmetric, _ = query([6, 6], [(0, 1, 1)])
    assert symmetric.automorphism_count() == 1, 'the two carbons swap'
    assert symmetric.count(ethane) == 2
    assert symmetric.count(ethane, automorphism_filter=True) == 1
    methanol, _ = mol([6, 8], [(0, 1, 1)])
    asymmetric, _ = query([6, 8], [(0, 1, 1)])
    assert asymmetric.automorphism_count() == 0
    assert asymmetric.count(methanol, automorphism_filter=True) == 1


def test_bond_orders_are_part_of_the_symmetry():
    """A path of four carbons is symmetric under reversal by adjacency alone; whether the
    reversal is really an automorphism is decided by the bond terms.  Two cases, because a group
    computed from adjacency alone reports the same row for both.
    """
    symmetric, _ = query([6, 6, 6, 6], [(0, 1, 1), (1, 2, 2), (2, 3, 1)])    # C-C=C-C
    broken, _ = query([6, 6, 6, 6], [(0, 1, 1), (1, 2, 2), (2, 3, 2)])       # C-C=C=C
    assert symmetric.automorphism_count() == 1, 'the reversal'
    assert broken.automorphism_count() == 0, 'reversal would swap a single bond for a double one'
    # hexa-2,4-diene: two overlapping C-C=C-C sites, each found in both orientations
    m, _ = mol([6] * 6, [(0, 1, 1), (1, 2, 2), (2, 3, 1), (3, 4, 2), (4, 5, 1)])
    assert symmetric.count(m) == 4
    assert symmetric.count(m, automorphism_filter=True) == 2


def test_ring_bond_orders_shrink_the_automorphism_group():
    """Cyclohexane's query carries all twelve dihedral symmetries; the Kekule query's alternating
    orders kill the odd rotations and the vertex-centred reflections, leaving six.  A group built
    from adjacency alone reports eleven rows for both.
    """
    plain = ring_query([1] * 6)
    kekule = ring_query([2, 1, 2, 1, 2, 1])
    assert plain.automorphism_count() == 11, 'six rotations and six reflections, less identity'
    assert kekule.automorphism_count() == 5, 'three rotations and three reflections, less identity'
    m = MoleculeContainer()
    ids = [m.add_atom(6) for _ in range(6)]
    for i, order in enumerate([2, 1, 2, 1, 2, 1]):
        m.add_bond(ids[i], ids[(i + 1) % 6], order)
    assert kekule.count(m) == 6
    assert kekule.count(m, automorphism_filter=True) == 1


def alternating_ring_query(decorate):
    """A six-carbon ring, every bond a plain single; `decorate(q, a, b)` adds primitives to the
    three alternating bonds.  1-WL cannot see the alternation -- every slot has two neighbours and
    sees one decorated and one plain bond -- so the initial partition and every refinement round
    leave all six slots in one class.  Whatever `decorate` does is therefore visible to the group
    only through _wterm_equal's byte compare of the bond terms.
    """
    q = QueryContainer()
    ids = [q.add_atom() for _ in range(6)]
    for sid in ids:
        q.atom_primitive(sid, 'element', 6)
    for i in range(6):
        a, b = ids[i], ids[(i + 1) % 6]
        q.add_bond(a, b)
        q.bond_primitive(a, b, 'bond_order', 1)
        if i % 2:
            decorate(q, a, b)
    return q


def cyclohexane():
    m = MoleculeContainer()
    ids = [m.add_atom(6) for _ in range(6)]
    for i in range(6):
        m.add_bond(ids[i], ids[(i + 1) % 6], 1)
    return m


def test_a_second_box_on_a_bond_term_is_part_of_the_symmetry():
    """Alternate a plain single bond with single-or-aromatic.  The two terms agree in their first
    box down to the byte -- the disjunction's `-` alternative compiles to exactly the plain box --
    and differ only in that one of them has a second box.  A verification that stops after box
    zero calls them equal, restores the odd rotations and reflections, and over-collapses the
    match count.
    """
    def single_or_aromatic(q, a, b):
        q.bond_operator(a, b, 'or')
        q.bond_primitive(a, b, 'bond_aromatic')

    q = alternating_ring_query(single_or_aromatic)
    assert q.automorphism_count() == 5, 'three rotations and three reflections, less identity'
    m = cyclohexane()
    assert q.count(m) == 12
    assert q.count(m, automorphism_filter=True) == 2


def test_an_any_list_on_a_bond_term_is_part_of_the_symmetry():
    """Alternate a plain single bond with single-and-in-a-six-ring.  A ring size compiles to an
    any list entry and nothing else -- box_fill_defaults writes no default over the ring size span
    -- so the two terms are byte-identical in neg[0..3] and differ only in the any list.  A
    verification that compares neg alone calls them equal and over-collapses the match count.
    """
    def in_a_six_ring(q, a, b):
        q.bond_operator(a, b, 'and_high')
        q.bond_primitive(a, b, 'ring_size', 6)

    q = alternating_ring_query(in_a_six_ring)
    assert q.automorphism_count() == 5, 'three rotations and three reflections, less identity'
    m = cyclohexane()
    assert q.count(m) == 12
    assert q.count(m, automorphism_filter=True) == 2


def test_a_component_group_is_part_of_the_canonical_form():
    """Two identical carbons in two components swap freely; put them in different reaction
    groups and they do not, because the component-group constraint is not symmetric in them.
    """
    ungrouped = QueryContainer()
    for _ in range(2):
        ungrouped.atom_primitive(ungrouped.add_atom(), 'element', 6)
    assert ungrouped.automorphism_count() == 1
    grouped = QueryContainer()
    for group in (0, 1):
        sid = grouped.add_atom()
        grouped.atom_primitive(sid, 'element', 6)
        grouped.set_group(sid, group)
    assert grouped.automorphism_count() == 0, 'group 0 and group 1 are not interchangeable'


def test_map_numbers_and_the_masked_flag_stay_out_of_the_canonical_form():
    """A map number names an atom for the reactor and the masked flag protects it from deletion.
    Neither constrains what the atom matches, so neither may break a symmetry: if they did, a
    template author would silently lose duplicate suppression by numbering their atoms.  Propane
    holds two C-C sites, each found in both directions.
    """
    m, _ = mol([6, 6, 6], [(0, 1, 1), (1, 2, 1)])
    plain, _ = query([6, 6], [(0, 1, 1)])
    assert plain.automorphism_count() == 1
    assert plain.count(m) == 4
    assert plain.count(m, automorphism_filter=True) == 2

    mapped, mapped_ids = query([6, 6], [(0, 1, 1)])
    mapped.set_map_number(mapped_ids[0], 1)
    mapped.set_map_number(mapped_ids[1], 2)
    assert mapped.automorphism_count() == 1, 'the map numbers differ, the constraints do not'
    assert mapped.count(m, automorphism_filter=True) == 2

    masked, masked_ids = query([6, 6], [(0, 1, 1)])
    masked.set_masked(masked_ids[0])
    assert masked.automorphism_count() == 1, 'one atom is masked, the other is not'
    assert masked.count(m, automorphism_filter=True) == 2


def test_the_filter_and_a_disconnected_query():
    # [C].[C] against propane: three unordered pairs of carbons, six ordered ones
    m, _ = mol([6, 6, 6], [(0, 1, 1), (1, 2, 1)])
    q = QueryContainer()
    for _ in range(2):
        q.atom_primitive(q.add_atom(), 'element', 6)
    assert q.count(m) == 6
    assert q.count(m, automorphism_filter=True) == 3


def test_cumulated_and_unmatched_hybridization_are_matchable():
    # derive_scalars gives allene's central carbon z=5 (two doubles, no triple).  The demand has to
    # be expressible, which takes all six bits of the span: a narrower hybridization mask leaves z5
    # and z6 setting no bit in it, so they slip past *every* z demand instead of failing all but their
    # own.
    m, ids = mol([6, 6, 6], [(0, 1, 2), (1, 2, 2)])                 # C=C=C
    assert m.hybridization_of(ids[1]) == 5

    # the two terminal carbons carry one double each, so they are z2 -- every z demand must select
    # exactly its own atoms and no others
    assert [one_atom(('element', 6), ('hybridization', z))[0].count(m) for z in range(1, 7)] == \
        [0, 2, 0, 0, 1, 0]


def test_a_sulfone_sulfur_is_not_sp2():
    # S(=O)(=O) is z5 by the same rule, which is why a sulfone/sulfonamide template written as
    # [S;z2] silently matches nothing.
    m, ids = mol([16, 8, 8, 6, 6], [(0, 1, 2), (0, 2, 2), (0, 3, 1), (0, 4, 1)])
    assert m.hybridization_of(ids[0]) == 5
    assert one_atom(('element', 16), ('hybridization', 2))[0].count(m) == 0
    assert one_atom(('element', 16), ('hybridization', 5))[0].count(m) == 1


# ---------------------------------------------------------------------------
# Task 13: component grouping
# ---------------------------------------------------------------------------

def two_carbons_grouped(group_a, group_b):
    """Two single-atom query components, each carbon, with the given group numbers."""
    q = QueryContainer()
    a = q.add_atom()
    q.atom_primitive(a, 'element', 6)
    b = q.add_atom()
    q.atom_primitive(b, 'element', 6)
    if group_a is not None:
        q.set_group(a, group_a)
    if group_b is not None:
        q.set_group(b, group_b)
    return q


def two_fragments():
    """C-C and O-O in one molecule: two components, four atoms."""
    return mol([6, 6, 8, 8], [(0, 1, 1), (2, 3, 1)])


def test_ungrouped_components_are_unconstrained():
    m, _ = mol([6, 6, 6], [(0, 1, 1)])          # C-C and a lone C: two components
    q = two_carbons_grouped(None, None)
    assert q.count(m) == 6, 'any ordered pair of distinct carbons'


def test_the_same_group_forces_the_same_component():
    m, _ = mol([6, 6, 6], [(0, 1, 1)])
    q = two_carbons_grouped(0, 0)
    assert q.count(m) == 2, 'only the two carbons of the C-C fragment'


def test_different_groups_force_different_components():
    m, _ = mol([6, 6, 6], [(0, 1, 1)])
    q = two_carbons_grouped(0, 1)
    assert q.count(m) == 4, 'one from the pair, one from the lone atom, either order'


def test_the_same_group_across_three_components():
    m, _ = mol([6, 6, 6, 6], [(0, 1, 1), (1, 2, 1)])   # a C3 chain and a lone C
    q = QueryContainer()
    for _ in range(3):
        sid = q.add_atom()
        q.atom_primitive(sid, 'element', 6)
        q.set_group(sid, 0)
    assert q.count(m) == 6, 'all three must sit in the C3 chain'


def test_a_group_of_one_is_still_a_constraint_against_the_others():
    m, _ = two_fragments()
    q = QueryContainer()
    a = q.add_atom()
    q.atom_primitive(a, 'element', 6)
    q.set_group(a, 0)
    b = q.add_atom()
    q.atom_primitive(b, 'element', 6)
    q.set_group(b, 1)
    assert q.count(m) == 0, 'both carbons live in the same molecule component'


def test_a_grouped_and_an_ungrouped_component_coexist():
    m, _ = mol([6, 6, 6], [(0, 1, 1)])
    q = two_carbons_grouped(0, None)
    assert q.count(m) == 6, 'the ungrouped atom is free, including inside group 0'


def test_grouping_a_connected_query_changes_nothing():
    m, _ = mol([6, 6], [(0, 1, 1)])
    q, ids = query([6, 6], [(0, 1, 1)])
    q.set_group(ids[0], 0)
    q.set_group(ids[1], 0)
    assert q.count(m) == 2


def test_groups_compose_with_the_automorphism_filter():
    m, _ = mol([6, 6, 6], [(0, 1, 1)])
    q = two_carbons_grouped(0, 0)
    assert q.count(m, automorphism_filter=True) == 1


def test_a_component_spanning_two_groups_is_an_error():
    q, ids = query([6, 6], [(0, 1, 1)])
    q.set_group(ids[0], 0)
    q.set_group(ids[1], 1)
    m, _ = mol([6, 6], [(0, 1, 1)])
    with pytest.raises(ValueError, match='two groups'):
        q.count(m)


def test_grouped_search_does_not_mutate_molecule_serialisation():
    """A grouped-query search must not modify the molecule's arena (F1 / ruling 1).

    ensure_component_labels appends SEG_COMPONENT_LABEL to the structure arena, which updates the
    segment table entries inside the header.  The header is part of the persistent prefix, so
    to_bytes() returns different bytes before vs. after the match.  With matcher-owned labels the
    arena is never touched, and the two serialisations must be identical.

    This test fails if labels are borrowed from the arena (ensure_component_labels path) and passes
    with the matcher-owned-memory fix.
    """
    m, _ = mol([6, 6, 8, 8], [(0, 1, 1), (2, 3, 1)])  # two components
    q = two_carbons_grouped(0, 0)
    before = m.to_bytes()
    q.count(m)
    after = m.to_bytes()
    assert before == after, 'grouped search must not alter the molecule serialisation'


# ---------------------------------------------------------------------------
# Task 14: the element-demand screen
# ---------------------------------------------------------------------------
# The screen is a sound lower bound: False means no embedding is possible,
# True means "maybe".  The count() and may_match() assertions together pin
# both directions: False is never returned when a match exists (soundness),
# and True is never returned when we know no match can exist.
#
# Ruling 1 note: may_match() calls query_may_match() directly, so the
# `is False` / `is True` assertions cover the screen function itself.  The
# call site in matcher_init is a performance optimisation only (ruling 5):
# removing it leaves every count correct and is therefore a known blind spot.
# Mutation testing confirms this; see the task-14 report.


def test_the_screen_rejects_a_smaller_molecule():
    m, _ = mol([6, 6], [(0, 1, 1)])
    q, _ = query([6, 6, 6], [(0, 1, 1), (1, 2, 1)])
    assert q.may_match(m) is False
    assert q.count(m) == 0


def test_the_screen_rejects_a_missing_element():
    m, _ = mol([6, 6, 6], [(0, 1, 1), (1, 2, 1)])
    q, _ = query([6, 7], [(0, 1, 1)])
    assert q.may_match(m) is False


def test_the_screen_counts_multiplicity():
    m, _ = mol([7, 6, 7], [(0, 1, 1), (1, 2, 1)])
    q = QueryContainer()
    ids = [q.add_atom() for _ in range(3)]
    for sid in ids:
        q.atom_primitive(sid, 'element', 7)
    assert q.may_match(m) is False, 'three nitrogens demanded, two available'


def test_the_screen_admits_an_exact_count():
    m, _ = mol([7, 6, 7], [(0, 1, 1), (1, 2, 1)])
    q = QueryContainer()
    for _ in range(2):
        sid = q.add_atom()
        q.atom_primitive(sid, 'element', 7)
    assert q.may_match(m) is True
    assert q.count(m) == 2


def test_an_element_list_atom_is_not_counted_by_the_screen():
    m, _ = mol([6, 6], [(0, 1, 1)])
    q = QueryContainer()
    a = q.add_atom()
    q.atom_primitive(a, 'element', 6)
    q.atom_operator(a, 'or')
    q.atom_primitive(a, 'element', 7)
    assert q.may_match(m) is True
    assert q.count(m) == 2


def test_the_screen_admits_when_no_descriptor_is_demanded():
    # a query that says nothing about degree must not be screened out by degree: this is the
    # shape Task 15 emits, and it is why propane stays a substructure of isobutane
    m, _ = mol([6, 6, 6, 6], [(0, 1, 1), (1, 2, 1), (1, 3, 1)])
    q, _ = query([6, 6, 6], [(0, 1, 1), (1, 2, 1)])
    assert q.may_match(m) is True
    assert q.count(m) == 6


def test_the_screen_uses_a_descriptor_the_query_does_demand():
    # [C;D2] genuinely does not match isobutane -- no carbon has exactly two heavy neighbours,
    # so rejecting it in the screen is sharpening, not a lost match
    m, _ = mol([6, 6, 6, 6], [(0, 1, 1), (1, 2, 1), (1, 3, 1)])
    q = QueryContainer()
    a = q.add_atom()
    q.atom_primitive(a, 'element', 6)
    q.atom_operator(a, 'and_low')
    q.atom_primitive(a, 'degree', 2)
    assert q.may_match(m) is False
    assert q.count(m) == 0


def test_an_or_over_a_span_demands_nothing_from_the_screen():
    m, _ = mol([6, 6, 6, 6], [(0, 1, 1), (1, 2, 1), (1, 3, 1)])
    q = QueryContainer()
    a = q.add_atom()
    q.atom_primitive(a, 'element', 6)
    q.atom_operator(a, 'and_low')
    q.atom_primitive(a, 'degree', 2)
    q.atom_operator(a, 'or')
    q.atom_primitive(a, 'degree', 3)
    assert q.may_match(m) is True, 'D2,D3 pins no single bit, so it cannot screen'
    assert q.count(m) == 1, 'only the central carbon has degree 3'


def test_the_screen_rejects_a_charge_the_molecule_lacks():
    m, _ = mol([6, 6], [(0, 1, 1)])
    q = QueryContainer()
    a = q.add_atom()
    q.atom_primitive(a, 'element', 6)
    q.atom_operator(a, 'and_low')
    q.atom_primitive(a, 'charge', 1)
    assert q.may_match(m) is False
    assert q.count(m) == 0


def test_the_screen_rejects_an_atom_free_molecule():
    m = MoleculeContainer()
    q, _ = query([6], [])
    assert q.may_match(m) is False


def test_query_search_methods_reject_none_rather_than_crashing():
    q, _ = query([6], [])
    with pytest.raises(TypeError):
        q.may_match(None)
    with pytest.raises(TypeError):
        q.count(None)
    with pytest.raises(TypeError):
        q.is_substructure(None)
    with pytest.raises(TypeError):
        list(q.get_mapping(None))
    with pytest.raises(TypeError):
        list(q.get_raw_mapping(None))
