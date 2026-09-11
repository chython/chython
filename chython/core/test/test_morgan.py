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
from chython.core import MoleculeContainer


def build(bonds, n=None, elements=None, charges=None):
    """bonds: list of (i, j, order) over 0-based positions. Returns (molecule, ids)."""
    if n is None:
        n = max(max(i, j) for i, j, _ in bonds) + 1 if bonds else 0
    m = MoleculeContainer()
    with m.edit():
        ids = [m.add_atom(6 if elements is None else elements[k],
                          charge=0 if charges is None else charges[k]) for k in range(n)]
        for i, j, o in bonds:
            m.add_bond(ids[i], ids[j], o)
    return m, ids


def cycle(orders, elements=None):
    n = len(orders)
    return build([(i, (i + 1) % n, orders[i]) for i in range(n)], elements=elements)


def classes(m, ids):
    """The symmetry partition as a set of frozensets of 0-based positions."""
    order = m.atoms_order
    groups = {}
    for pos, sid in enumerate(ids):
        groups.setdefault(order[sid], set()).add(pos)
    return {frozenset(g) for g in groups.values()}


# --- degenerate inputs ---

def test_empty_molecule_has_empty_order():
    m = MoleculeContainer()
    assert m.atoms_order == {}
    assert m.atoms_order_classes == 0


def test_single_atom_is_one_class():
    m, ids = build([], n=1)
    assert m.atoms_order == {ids[0]: 1}
    assert m.atoms_order_classes == 1


def test_two_isolated_identical_atoms_share_a_class():
    m, ids = build([], n=2)
    assert m.atoms_order[ids[0]] == m.atoms_order[ids[1]]


# --- ranks are 1-based and dense ---

def test_ranks_are_one_based_and_contiguous():
    m, ids = build([(0, 1, 1), (1, 2, 1), (2, 3, 2)], elements=[6, 6, 7, 8])
    ranks = sorted(m.atoms_order.values())
    assert ranks == list(range(1, len(ranks) + 1))
    assert m.atoms_order_classes == len(ranks)


def test_class_count_matches_distinct_ranks():
    m, ids = cycle([2, 1, 2, 1, 2, 1])
    assert m.atoms_order_classes == len({*m.atoms_order.values()})


def test_order_is_keyed_by_stable_id_not_index():
    m, ids = build([(0, 1, 1), (1, 2, 1)], elements=[8, 6, 8])
    assert set(m.atoms_order) == set(ids)
    with m.edit():
        m.delete_atom(ids[0])
    assert set(m.atoms_order) == {ids[1], ids[2]}


# --- symmetry the refinement must find ---

def test_cyclohexane_is_a_single_class():
    m, ids = cycle([1] * 6)
    assert classes(m, ids) == {frozenset(range(6))}


def test_kekule_benzene_is_a_single_class():
    # alternating orders, but every atom still sees one single and one double bond
    m, ids = cycle([2, 1, 2, 1, 2, 1])
    assert classes(m, ids) == {frozenset(range(6))}


def test_cyclopropane_is_a_single_class():
    m, ids = cycle([1, 1, 1])
    assert m.atoms_order_classes == 1


def test_pentane_is_symmetric_about_its_middle():
    m, ids = build([(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1)])
    assert classes(m, ids) == {frozenset({0, 4}), frozenset({1, 3}), frozenset({2})}


def test_isobutane_methyls_are_equivalent():
    m, ids = build([(0, 1, 1), (0, 2, 1), (0, 3, 1)])
    assert classes(m, ids) == {frozenset({0}), frozenset({1, 2, 3})}


def test_neopentane_methyls_are_equivalent():
    m, ids = build([(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1)])
    assert classes(m, ids) == {frozenset({0}), frozenset({1, 2, 3, 4})}


def test_two_ethane_fragments_are_one_class():
    m, ids = build([(0, 1, 1), (2, 3, 1)])
    assert m.atoms_order_classes == 1


def test_disconnected_components_of_different_size_split():
    # two methanes and one ethane: the lone atoms cannot match the bonded ones
    m, ids = build([(2, 3, 1)], n=4)
    assert classes(m, ids) == {frozenset({0, 1}), frozenset({2, 3})}


# --- asymmetry the refinement must not collapse ---

def test_ethanol_has_no_symmetry():
    m, ids = build([(0, 1, 1), (1, 2, 1)], elements=[6, 6, 8])
    assert m.atoms_order_classes == 3


def test_toluene_kekule_ring_is_desymmetrised_by_bond_orders():
    # A fixed Kekule form is genuinely less symmetric than the aromatic ring it stands for:
    # the ring carbon double-bonded to the ipso carbon is not equivalent to the one single-
    # bonded to it. Seven classes is correct here; the aromatic form would give five.
    m, ids = cycle([2, 1, 2, 1, 2, 1])
    with m.edit():
        methyl = m.add_atom(6)
        m.add_bond(ids[0], methyl, 1)
    assert m.atoms_order_classes == 7


def test_propene_terminal_carbons_differ():
    m, ids = build([(0, 1, 2), (1, 2, 1)])
    assert m.atoms_order_classes == 3


def test_bond_order_alone_separates_atoms():
    # butane vs 2-butene skeletons: same graph, different orders, different partitions
    single, ids_s = build([(0, 1, 1), (1, 2, 1), (2, 3, 1)])
    inner, ids_i = build([(0, 1, 1), (1, 2, 2), (2, 3, 1)])
    assert classes(single, ids_s) == classes(inner, ids_i) == {frozenset({0, 3}), frozenset({1, 2})}
    outer, ids_o = build([(0, 1, 2), (1, 2, 1), (2, 3, 1)])
    assert outer.atoms_order_classes == 4


# --- each invariant field is actually read ---
#
# The "ranks before" tests here hold only because every atom in those molecules is already in its
# own class at round 0, so refinement never runs and the round-0 ordering survives. Once a
# hashing round happens the class numbering is arbitrary; assert partitions, not order.

def test_element_separates_otherwise_identical_atoms():
    m, ids = build([(0, 1, 1), (0, 2, 1)], elements=[6, 7, 8])
    assert m.atoms_order_classes == 3


def test_carbon_ranks_before_nitrogen():
    # ordering by the invariant, not by its hash: element is the most significant field
    m, ids = build([(0, 1, 1)], elements=[7, 6])
    assert m.atoms_order[ids[1]] < m.atoms_order[ids[0]]


def test_charge_separates_otherwise_identical_atoms():
    m, ids = build([], n=2, charges=[0, 1])
    assert m.atoms_order_classes == 2


def test_negative_charge_ranks_before_neutral():
    # the charge field is biased so -4 .. +8 stays unsigned and ordered
    m, ids = build([], n=2, charges=[-1, 0])
    assert m.atoms_order[ids[0]] < m.atoms_order[ids[1]]


def test_isotope_separates_otherwise_identical_atoms():
    m = MoleculeContainer()
    with m.edit():
        a = m.add_atom(6)
        b = m.add_atom(6, isotope=13)
        m.add_bond(a, b, 1)
    assert m.atoms_order[a] != m.atoms_order[b]


def test_radical_separates_otherwise_identical_atoms():
    m = MoleculeContainer()
    with m.edit():
        a = m.add_atom(6)
        b = m.add_atom(6, radical=True)
        m.add_bond(a, b, 1)
    assert m.atoms_order[a] != m.atoms_order[b]


def test_ring_and_chain_atoms_never_share_a_class():
    # Cyclopropane and propane in one molecule. The in_ring bit separates them from round 0, and
    # classes only split, so they stay separate. Refinement alone would find this too -- ring
    # membership is largely re-derivable from topology, which is why this asserts the partition
    # rather than claiming the bit is the sole cause.
    m, ids = build([(0, 1, 1), (1, 2, 1), (2, 0, 1), (3, 4, 1), (4, 5, 1)])
    order = m.atoms_order
    assert [m.in_ring_of(s) for s in ids] == [True, True, True, False, False, False]
    assert not {order[s] for s in ids[:3]} & {order[s] for s in ids[3:]}


def test_implicit_hydrogens_separate_otherwise_identical_atoms():
    # core states hydrogen counts rather than deriving them, so set them explicitly
    m = MoleculeContainer()
    with m.edit():
        a = m.add_atom(6, implicit_h=3)
        b = m.add_atom(6, implicit_h=2)
        m.add_bond(a, b, 1)
    assert m.atoms_order[a] != m.atoms_order[b]


def test_fewer_hydrogens_rank_first():
    m = MoleculeContainer()
    with m.edit():
        a = m.add_atom(6, implicit_h=1)
        b = m.add_atom(6, implicit_h=3)
        m.add_bond(a, b, 1)
    assert m.atoms_order[a] < m.atoms_order[b]


# --- refinement must propagate beyond immediate neighbours ---

def test_refinement_propagates_along_a_chain():
    # heptane: symmetry only resolves after information has walked three bonds inward
    m, ids = build([(i, i + 1, 1) for i in range(6)])
    assert classes(m, ids) == {frozenset({0, 6}), frozenset({1, 5}), frozenset({2, 4}),
                               frozenset({3})}


def test_distant_substituent_splits_a_symmetric_looking_pair():
    # two branches identical for two bonds, then diverging: neighbour-only keying merges
    # positions 1 and 4, refinement must not
    m, ids = build([(0, 1, 1), (1, 2, 1), (2, 3, 1),
                    (0, 4, 1), (4, 5, 1), (5, 6, 2)])
    order = m.atoms_order
    assert order[ids[1]] != order[ids[4]]


def test_naphthalene_kekule_partition():
    # C2h symmetry once the double bonds are fixed: five orbits of two
    m, ids = build([(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 0, 1),
                    (4, 6, 1), (6, 7, 2), (7, 8, 1), (8, 9, 2), (9, 3, 1)])
    part = classes(m, ids)
    assert len(part) == 5
    assert all(len(g) == 2 for g in part)


def test_c60_is_vertex_transitive():
    # the truncated icosahedron has one vertex orbit, so a correct refinement must collapse all
    # 60 carbons -- and must not be fooled into splitting them by the pentagon/hexagon faces
    from .test_rings_c60 import C60_EDGES
    m = MoleculeContainer()
    with m.edit():
        ids = [m.add_atom(6) for _ in range(60)]
        for i, j in C60_EDGES:
            m.add_bond(ids[i], ids[j], 1)
    assert m.atoms_order_classes == 1


# --- refined_order: the seeded kernel stereo will use ---

def test_uniform_seed_ignores_atom_records():
    # every atom starting in one class, so only topology and bond orders can split them: the
    # heteroatom in the middle of a chain becomes indistinguishable from a carbon there
    m, ids = build([(0, 1, 1), (1, 2, 1)], elements=[6, 8, 6])
    seed = dict.fromkeys(ids, 0)
    order = m.refined_order(seed)
    assert order[ids[0]] == order[ids[2]]
    assert len({*order.values()}) == 2


def test_seed_split_survives_refinement():
    # cyclohexane is one class, but seeding one atom apart must keep it apart -- and must
    # propagate, since its neighbours now see a distinct rank
    m, ids = cycle([1] * 6)
    assert m.atoms_order_classes == 1
    seed = dict.fromkeys(ids, 0)
    seed[ids[0]] = 1
    order = m.refined_order(seed)
    assert order[ids[0]] != order[ids[1]]
    assert order[ids[1]] == order[ids[5]]     # the two neighbours stay equivalent
    assert order[ids[2]] == order[ids[4]]
    assert len({*order.values()}) == 4        # seeded, ortho, meta, para


def test_seed_labels_need_not_be_dense_or_one_based():
    m, ids = build([(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1)])
    sparse = m.refined_order({s: v for s, v in zip(ids, (700, 40, 40, 40, 700))})
    dense = m.refined_order({s: v for s, v in zip(ids, (2, 1, 1, 1, 2))})
    assert sparse == dense
    assert sorted({*sparse.values()}) == list(range(1, len({*sparse.values()}) + 1))


def test_seed_reproducing_the_default_reproduces_the_default():
    m, ids = build([(0, 1, 1), (1, 2, 1), (2, 3, 2), (0, 4, 1)], elements=[6, 6, 7, 8, 6])
    assert m.refined_order(m.atoms_order) == m.atoms_order


def test_refined_order_rejects_a_missing_atom():
    m, ids = build([(0, 1, 1), (1, 2, 1)])
    with pytest.raises(KeyError):
        m.refined_order({ids[0]: 1, ids[1]: 1})


def test_refined_order_rejects_a_negative_label():
    m, ids = build([(0, 1, 1)])
    with pytest.raises(OverflowError):
        m.refined_order({ids[0]: -1, ids[1]: 0})


def test_refined_order_on_empty_molecule():
    assert MoleculeContainer().refined_order({}) == {}


def test_refined_order_does_not_touch_the_cache():
    m, ids = cycle([1] * 6)
    cached = m.atoms_order
    seed = dict.fromkeys(ids, 0)
    seed[ids[0]] = 1
    assert m.refined_order(seed) != cached
    assert m.atoms_order is cached


# --- caching ---

def test_order_is_cached_between_reads():
    m, ids = cycle([1] * 6)
    assert m.atoms_order is m.atoms_order


def test_edit_invalidates_the_cache():
    m, ids = build([(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1)])
    before = m.atoms_order
    assert before[ids[0]] == before[ids[4]]
    with m.edit():
        m.add_bond(ids[4], m.add_atom(8), 1)
    after = m.atoms_order
    assert after is not before
    assert after[ids[0]] != after[ids[4]]


def test_cache_survives_an_edit_scope_that_changes_nothing():
    m, ids = cycle([1] * 6)
    before = m.atoms_order
    with m.edit():
        pass
    assert m.atoms_order == before


# --- ferrocene: order-8 bonds participate in the invariant but not in ring perception ---

def test_ferrocene_rings_are_equivalent():
    # two cyclopentadienyl rings dative-bonded to one iron. The refinement reads the order-8
    # bonds like any other, so both rings collapse to one class and the iron stands alone.
    bonds = []
    for base in (0, 5):
        orders = (2, 1, 2, 1, 1)
        for k in range(5):
            bonds.append((base + k, base + (k + 1) % 5, orders[k]))
    for k in range(10):
        bonds.append((k, 10, 8))
    m, ids = build(bonds, elements=[6] * 10 + [26])
    order = m.atoms_order
    assert order[ids[10]] not in {order[s] for s in ids[:10]}
    assert {order[s] for s in ids[:5]} == {order[s] for s in ids[5:10]}
