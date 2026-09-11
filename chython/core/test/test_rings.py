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


def build(bonds, n=None, elements=None):
    """bonds: list of (i, j, order) over 0-based positions. Returns (molecule, ids)."""
    if n is None:
        n = max(max(i, j) for i, j, _ in bonds) + 1
    m = MoleculeContainer()
    with m.edit():
        ids = [m.add_atom(6 if elements is None else elements[k]) for k in range(n)]
        for i, j, o in bonds:
            m.add_bond(ids[i], ids[j], o)
    return m, ids


def ring_flags(m, ids):
    return [m.in_ring_of(s) for s in ids]


def test_chain_has_no_ring_atoms_or_bonds():
    m, ids = build([(0, 1, 1), (1, 2, 1)])
    assert ring_flags(m, ids) == [False, False, False]
    assert m.bond_in_ring(ids[0], ids[1]) is False


def test_cyclohexane_is_all_ring():
    m, ids = build([(i, (i + 1) % 6, 1) for i in range(6)])
    assert ring_flags(m, ids) == [True] * 6
    assert all(m.bond_in_ring(ids[i], ids[(i + 1) % 6]) for i in range(6))


def test_toluene_methyl_bond_is_a_bridge():
    bonds = [(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 0, 1),
             (0, 6, 1)]
    m, ids = build(bonds)
    assert ring_flags(m, ids) == [True] * 6 + [False]
    assert m.bond_in_ring(ids[0], ids[6]) is False
    assert m.bond_in_ring(ids[0], ids[1]) is True


def test_biphenyl_linking_bond_is_a_bridge():
    ring_a = [(i, (i + 1) % 6, 1) for i in range(6)]
    ring_b = [(6 + i, 6 + (i + 1) % 6, 1) for i in range(6)]
    m, ids = build(ring_a + ring_b + [(0, 6, 1)])
    assert ring_flags(m, ids) == [True] * 12
    assert m.bond_in_ring(ids[0], ids[6]) is False


def test_naphthalene_fusion_bond_is_not_a_bridge():
    # 0-1-2-3-4-5-0 fused to 4-5 via 5-6-7-8-9-4
    bonds = [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1), (5, 0, 1),
             (4, 9, 1), (9, 8, 1), (8, 7, 1), (7, 6, 1), (6, 5, 1)]
    m, ids = build(bonds)
    assert ring_flags(m, ids) == [True] * 10
    assert m.bond_in_ring(ids[4], ids[5]) is True


def test_spiro_atom_is_in_ring_and_both_rings_are_bridgeless():
    # cyclobutane 0-1-2-3-0 spiro-fused at 0 to cyclopentane 0-4-5-6-7-0
    bonds = [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 0, 1),
             (0, 4, 1), (4, 5, 1), (5, 6, 1), (6, 7, 1), (7, 0, 1)]
    m, ids = build(bonds)
    assert ring_flags(m, ids) == [True] * 8
    assert all(m.bond_in_ring(ids[i], ids[j]) for i, j, _ in bonds)


def test_disconnected_components_are_handled_independently():
    m, ids = build([(0, 1, 1), (1, 2, 1), (2, 0, 1), (3, 4, 1)])
    assert ring_flags(m, ids) == [True, True, True, False, False]


def test_isolated_atom_is_not_in_ring():
    m = MoleculeContainer()
    a = m.add_atom(6)
    assert m.in_ring_of(a) is False


def test_long_chain_does_not_overflow_the_stack():
    n = 20000
    m, ids = build([(i, i + 1, 1) for i in range(n - 1)], n=n)
    assert m.in_ring_of(ids[0]) is False
    assert m.in_ring_of(ids[n - 1]) is False
    assert not any(m.in_ring_of(s) for s in ids)


def test_bond_and_atom_views_expose_in_ring():
    m, ids = build([(i, (i + 1) % 6, 1) for i in range(6)] + [(0, 6, 1)])
    assert m.atom(ids[6]).in_ring is False
    assert m.atom(ids[0]).in_ring is True
    assert m.bond(ids[0], ids[6]).in_ring is False
    assert m.bond(ids[0], ids[1]).in_ring is True


def test_opening_a_ring_clears_in_ring_on_every_atom():
    m, ids = build([(i, (i + 1) % 6, 1) for i in range(6)])
    assert ring_flags(m, ids) == [True] * 6
    m.delete_bond(ids[0], ids[1])
    # atom_t.flags is carried forward across the fold, so the pass must write an
    # authoritative value rather than only setting the bit when it finds a ring bond
    assert ring_flags(m, ids) == [False] * 6
    assert m.atom(ids[0]).in_ring is False


def test_deleting_a_ring_atom_clears_in_ring_on_the_survivors():
    m, ids = build([(i, (i + 1) % 6, 1) for i in range(6)])
    m.delete_atom(ids[0])
    assert ring_flags(m, ids[1:]) == [False] * 5


def test_bond_in_ring_is_symmetric():
    m, ids = build([(i, (i + 1) % 6, 1) for i in range(6)] + [(0, 6, 1)])
    for a, b in [(0, 1), (0, 6)]:
        assert m.bond_in_ring(ids[a], ids[b]) is m.bond_in_ring(ids[b], ids[a])


def test_half_edge_flag_values_are_stable():
    from chython.core import _core as _structure
    assert _structure.HE_IN_RING == 1
    assert _structure.HE_AROMATIC == 2


def ring_sizes(m):
    return sorted(len(r) for r in m.rings)


def test_no_rings_in_a_chain():
    m, ids = build([(0, 1, 1), (1, 2, 1)])
    assert m.rings == []
    assert m.rings_count == 0


def test_benzene_has_one_ring():
    m, ids = build([(i, (i + 1) % 6, 1) for i in range(6)])
    assert ring_sizes(m) == [6]
    assert set(m.rings[0]) == set(ids)


def test_naphthalene_two_six_rings_not_a_ten_ring():
    bonds = [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1), (5, 0, 1),
             (4, 9, 1), (9, 8, 1), (8, 7, 1), (7, 6, 1), (6, 5, 1)]
    m, ids = build(bonds)
    assert ring_sizes(m) == [6, 6]


def test_anthracene_reports_only_six_rings():
    # three linearly fused six-rings, 14 carbons
    bonds = [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1), (5, 0, 1),
             (2, 6, 1), (6, 7, 1), (7, 8, 1), (8, 9, 1), (9, 3, 1),
             (7, 10, 1), (10, 11, 1), (11, 12, 1), (12, 13, 1), (13, 8, 1)]
    m, ids = build(bonds)
    assert ring_sizes(m) == [6, 6, 6]


def test_azulene_five_and_seven():
    # bicyclo[5.3.0], 10 carbons: 5-ring 0-1-2-3-4, 7-ring 0-4-5-6-7-8-9
    bonds = [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 0, 1),
             (4, 5, 1), (5, 6, 1), (6, 7, 1), (7, 8, 1), (8, 9, 1), (9, 0, 1)]
    m, ids = build(bonds)
    assert ring_sizes(m) == [5, 7]


def test_indane_five_and_six():
    bonds = [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1), (5, 0, 1),
             (0, 6, 1), (6, 7, 1), (7, 8, 1), (8, 1, 1)]
    m, ids = build(bonds)
    assert ring_sizes(m) == [5, 6]


def test_norbornane_two_five_rings_and_no_six_ring():
    # bicyclo[2.2.1]heptane: bridgeheads 0 and 3
    bonds = [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1), (5, 0, 1),
             (0, 6, 1), (6, 3, 1)]
    m, ids = build(bonds)
    # the six-ring is the GF(2) sum of the two five-rings, both strictly
    # shorter than it, so it is not relevant
    assert ring_sizes(m) == [5, 5]
    assert m.rings_count == 2


def test_bicyclo222octane_reports_two_six_rings_at_rank_two():
    # bridgeheads 0 and 4, three two-carbon bridges. All three 6-rings are relevant, but
    # `rings` is a minimum cycle basis, so it holds two of them and their GF(2) sum is the
    # third. The per-atom descriptors still see every relevant cycle.
    bonds = [(0, 1, 1), (1, 2, 1), (2, 4, 1),
             (0, 3, 1), (3, 5, 1), (5, 4, 1),
             (0, 6, 1), (6, 7, 1), (7, 4, 1)]
    m, ids = build(bonds)
    assert ring_sizes(m) == [6, 6]
    assert m.rings_count == 2          # circuit rank is 9 - 8 + 1 = 2
    assert m.ring_count_of(ids[0]) == 3          # bridgehead: all three relevant 6-rings
    assert m.ring_sizes_of(ids[1]) == frozenset({6})


def test_adamantane_basis_is_three_of_its_four_relevant_six_rings():
    # bridgeheads 0-3, one CH2 (4-9) bridging each of the six bridgehead pairs
    bonds = [(0, 4, 1), (4, 1, 1), (0, 5, 1), (5, 2, 1), (0, 6, 1), (6, 3, 1),
             (1, 7, 1), (7, 2, 1), (1, 8, 1), (8, 3, 1), (2, 9, 1), (9, 3, 1)]
    m, ids = build(bonds)
    assert m.rings_count == 3          # circuit rank is 12 - 10 + 1 = 3
    assert ring_sizes(m) == [6, 6, 6]


def test_cubane_basis_is_five_of_its_six_faces():
    # the sixth face is the GF(2) sum of the other five, so a basis cannot hold it
    bonds = [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 0, 1),
             (4, 5, 1), (5, 6, 1), (6, 7, 1), (7, 4, 1),
             (0, 4, 1), (1, 5, 1), (2, 6, 1), (3, 7, 1)]
    m, ids = build(bonds)
    assert m.rings_count == 5          # circuit rank is 12 - 8 + 1 = 5
    assert ring_sizes(m) == [4] * 5


def test_spiro34octane_two_rings():
    bonds = [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 0, 1),
             (0, 4, 1), (4, 5, 1), (5, 6, 1), (6, 7, 1), (7, 0, 1)]
    m, ids = build(bonds)
    assert ring_sizes(m) == [4, 5]


def test_biphenyl_two_rings_not_a_macrocycle():
    ring_a = [(i, (i + 1) % 6, 1) for i in range(6)]
    ring_b = [(6 + i, 6 + (i + 1) % 6, 1) for i in range(6)]
    m, ids = build(ring_a + ring_b + [(0, 6, 1)])
    assert ring_sizes(m) == [6, 6]


def test_macrocycle_is_reported_at_full_size():
    n = 30
    m, ids = build([(i, (i + 1) % n, 1) for i in range(n)])
    assert ring_sizes(m) == [30]
    assert all(m.in_ring_of(s) for s in ids)


def test_rings_are_reported_in_stable_id_terms():
    m, ids = build([(i, (i + 1) % 6, 1) for i in range(6)])
    assert set(m.rings[0]) == set(ids)


def test_k33_basis_is_four_of_its_nine_relevant_four_rings():
    # girth 4, and C(3,2) x C(3,2) = 9 four-cycles.  No shorter cycle exists, so none of the
    # nine can be a GF(2) sum of strictly shorter cycles: all nine are relevant by definition.
    # Circuit rank is only 9 - 6 + 1 = 4, so the basis holds four -- but every vertex still
    # counts all six relevant 4-rings it lies on, which is what the descriptors are for.
    bonds = [(a, b, 1) for a in (0, 1, 2) for b in (3, 4, 5)]
    m, ids = build(bonds)
    assert m.rings_count == 4
    assert ring_sizes(m) == [4] * 4
    for s in ids:
        assert m.ring_count_of(s) == 6      # C(2,1) x C(2,1) x ... : 6 of the 9 touch each vertex
        assert m.ring_sizes_of(s) == frozenset({4})


def test_tricyclic_cage_basis_is_three_of_its_five_relevant_four_rings():
    # cyclobutane C1-C6-C2-C5 with one CH2 bridging C5..C6 and another bridging C1..C2:
    # five 4-rings, all of girth size, all relevant.  Circuit rank is 8 - 6 + 1 = 3.
    bonds = [(1, 6, 1), (0, 5, 1), (2, 3, 1), (0, 6, 1),
             (1, 3, 1), (2, 6, 1), (2, 5, 1), (1, 5, 1)]
    m, ids = build(bonds)
    assert m.rings_count == 3
    assert ring_sizes(m) == [4] * 3


def test_empty_molecule_has_no_rings():
    m = MoleculeContainer()
    assert m.rings == []
    assert m.rings_count == 0


def test_single_atom_has_no_rings():
    m = MoleculeContainer()
    a = m.add_atom(6)
    assert m.rings == []
    assert m.rings_count == 0
    assert m.in_ring_of(a) is False


def test_rings_survive_a_second_edit():
    m, ids = build([(i, (i + 1) % 6, 1) for i in range(6)])
    assert m.rings_count == 1
    with m.edit():
        m.add_atom(6)
    assert m.rings_count == 1
    assert ring_sizes(m) == [6]


def test_large_macrocycle_does_not_overflow_the_path_stack():
    # the shortest-path DAG depth is ~n/2 here; a recursive enumerator raises RecursionError
    n = 2000
    m, ids = build([(i, (i + 1) % n, 1) for i in range(n)])
    assert m.rings_count == 1
    assert ring_sizes(m) == [n]


def test_chain_atoms_have_no_ring_descriptors():
    m, ids = build([(0, 1, 1), (1, 2, 1)])
    for s in ids:
        assert m.ring_count_of(s) == 0
        assert m.ring_sizes_of(s) == frozenset()
        assert m.macrocycle_of(s) is False
        assert m.in_ring_of(s) == (m.ring_count_of(s) > 0)
    assert not m.shares_ring(ids[0], ids[1])
    # a chain atom shares no ring with itself
    assert not m.shares_ring(ids[0], ids[0])


def test_benzene_every_atom_in_one_six_ring():
    m, ids = build([(i, (i + 1) % 6, 1) for i in range(6)])
    for s in ids:
        assert m.ring_count_of(s) == 1
        assert m.ring_sizes_of(s) == frozenset({6})
        assert m.ring_sizes_word_of(s) == 1 << 6
        assert m.in_ring_of(s) == (m.ring_count_of(s) > 0)
    assert m.shares_ring(ids[0], ids[3])
    # a ring atom shares its ring with itself
    assert m.shares_ring(ids[0], ids[0])


def test_naphthalene_fusion_atoms_carry_count_two():
    bonds = [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1), (5, 0, 1),
             (4, 9, 1), (9, 8, 1), (8, 7, 1), (7, 6, 1), (6, 5, 1)]
    m, ids = build(bonds)
    counts = [m.ring_count_of(s) for s in ids]
    assert counts == [1, 1, 1, 1, 2, 2, 1, 1, 1, 1]
    for s in ids:
        assert m.ring_sizes_of(s) == frozenset({6})
        assert m.in_ring_of(s) == (m.ring_count_of(s) > 0)


def test_anthracene_atoms_see_only_six_rings():
    bonds = [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1), (5, 0, 1),
             (2, 6, 1), (6, 7, 1), (7, 8, 1), (8, 9, 1), (9, 3, 1),
             (7, 10, 1), (10, 11, 1), (11, 12, 1), (12, 13, 1), (13, 8, 1)]
    m, ids = build(bonds)
    # no atom of anthracene belongs to a ring of any size other than six --
    # not to the 10-ring, not to the 14-ring
    for s in ids:
        assert m.ring_sizes_of(s) == frozenset({6})
        assert m.in_ring_of(s) == (m.ring_count_of(s) > 0)
    counts = [m.ring_count_of(s) for s in ids]
    assert counts == [1, 1, 2, 2, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1]


def test_azulene_fusion_atoms_carry_both_sizes():
    bonds = [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 0, 1),
             (4, 5, 1), (5, 6, 1), (6, 7, 1), (7, 8, 1), (8, 9, 1), (9, 0, 1)]
    m, ids = build(bonds)
    assert m.ring_sizes_of(ids[0]) == frozenset({5, 7})
    assert m.ring_sizes_of(ids[4]) == frozenset({5, 7})
    assert m.ring_sizes_of(ids[2]) == frozenset({5})
    assert m.ring_sizes_of(ids[6]) == frozenset({7})
    assert m.ring_count_of(ids[0]) == 2
    assert m.ring_count_of(ids[6]) == 1


def test_spiro_atom_carries_both_sizes_but_arms_do_not_share_a_ring():
    bonds = [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 0, 1),
             (0, 4, 1), (4, 5, 1), (5, 6, 1), (6, 7, 1), (7, 0, 1)]
    m, ids = build(bonds)
    assert m.ring_sizes_of(ids[0]) == frozenset({4, 5})
    assert m.ring_count_of(ids[0]) == 2
    assert m.shares_ring(ids[0], ids[2])
    assert m.shares_ring(ids[0], ids[5])
    assert not m.shares_ring(ids[2], ids[5])


def test_norbornane_bridge_atoms_are_in_both_five_rings():
    bonds = [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1), (5, 0, 1),
             (0, 6, 1), (6, 3, 1)]
    m, ids = build(bonds)
    counts = [m.ring_count_of(s) for s in ids]
    assert counts == [2, 1, 1, 2, 1, 1, 2]
    for s in ids:
        assert m.ring_sizes_of(s) == frozenset({5})
        assert m.in_ring_of(s) == (m.ring_count_of(s) > 0)


def test_adamantane_bridgeheads_are_in_three_rings():
    bonds = [(0, 4, 1), (4, 1, 1), (0, 5, 1), (5, 2, 1), (0, 6, 1), (6, 3, 1),
             (1, 7, 1), (7, 2, 1), (1, 8, 1), (8, 3, 1), (2, 9, 1), (9, 3, 1)]
    m, ids = build(bonds)
    counts = [m.ring_count_of(s) for s in ids]
    assert counts == [3, 3, 3, 3, 2, 2, 2, 2, 2, 2]
    for s in ids:
        assert m.ring_sizes_of(s) == frozenset({6})
        assert m.macrocycle_of(s) is False


def test_cubane_every_vertex_in_three_faces():
    bonds = [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 0, 1),
             (4, 5, 1), (5, 6, 1), (6, 7, 1), (7, 4, 1),
             (0, 4, 1), (1, 5, 1), (2, 6, 1), (3, 7, 1)]
    m, ids = build(bonds)
    for s in ids:
        assert m.ring_count_of(s) == 3
        assert m.ring_sizes_of(s) == frozenset({4})
        assert m.in_ring_of(s) == (m.ring_count_of(s) > 0)
    # every pair of cubane vertices lies on a common face except the four
    # body diagonals
    assert m.shares_ring(ids[0], ids[6]) is False
    assert m.shares_ring(ids[0], ids[5]) is True


def test_macrocycle_is_reported_without_an_exact_size():
    n = 30
    m, ids = build([(i, (i + 1) % n, 1) for i in range(n)])
    for s in ids:
        assert m.ring_count_of(s) == 1
        assert m.ring_sizes_of(s) == frozenset()   # 30 is out of the 3-24 range
        assert m.macrocycle_of(s) is True
        assert m.ring_sizes_word_of(s) == 1


def test_macrocycle_flag_starts_exactly_where_exact_sizes_stop():
    # 24 is the largest size ring_sizes can name; everything above it is reported as a
    # macrocycle instead, and `rings` is where the exact size comes from
    for n in (24, 25, 32, 33, 48, 49, 60):
        m, ids = build([(i, (i + 1) % n, 1) for i in range(n)])
        if n == 24:
            assert m.ring_sizes_of(ids[0]) == frozenset({24})
            assert m.macrocycle_of(ids[0]) is False
        else:
            assert m.ring_sizes_of(ids[0]) == frozenset()
            assert m.macrocycle_of(ids[0]) is True
        assert m.atom(ids[0]).macrocycle == m.macrocycle_of(ids[0])
        assert len(m.rings[0]) == n


def test_ring_bitmap_survives_a_property_edit():
    bonds = [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1), (5, 0, 1)]
    m, ids = build(bonds)
    m.set_charge(ids[0], 1)          # an atomic mutation folds a fresh arena
    assert m.ring_count_of(ids[0]) == 1
    assert m.shares_ring(ids[0], ids[3])


def test_shares_ring_rejects_unknown_ids():
    m, ids = build([(i, (i + 1) % 6, 1) for i in range(6)])
    with pytest.raises(KeyError):
        m.shares_ring(ids[0], 99999)


def test_opening_a_ring_clears_the_ring_descriptors():
    m, ids = build([(i, (i + 1) % 6, 1) for i in range(6)])
    assert m.ring_sizes_of(ids[0]) == frozenset({6})
    m.delete_bond(ids[0], ids[1])
    # atom_t is memcpy'd forward across the fold, so ring_sizes and ring_counts
    # must be written authoritatively rather than accumulated
    for s in ids:
        assert m.ring_sizes_of(s) == frozenset()
        assert m.ring_count_of(s) == 0
    assert not m.shares_ring(ids[0], ids[3])


def test_shrinking_a_ring_replaces_the_size_rather_than_adding_one():
    # cyclohexane -> cyclopentane: drop one ring CH2 and close the gap
    m, ids = build([(i, (i + 1) % 6, 1) for i in range(6)])
    with m.edit():
        m.delete_atom(ids[5])
        m.add_bond(ids[4], ids[0], 1)
    assert m.ring_sizes_of(ids[0]) == frozenset({5})
    assert m.ring_count_of(ids[0]) == 1


def test_ring_bitmap_spans_two_words():
    # 10x10 grid: 81 unit squares, so words = ceil(81 / 64) = 2. Every molecule
    # elsewhere in this suite has <= 32 rings and therefore exercises words == 1
    # only, which cannot distinguish the word index from the bit index.
    coord = {}
    bonds = []
    idx = 0
    for r in range(10):
        for c in range(10):
            coord[(r, c)] = idx
            idx += 1
    for r in range(10):
        for c in range(10):
            if c + 1 < 10:
                bonds.append((coord[(r, c)], coord[(r, c + 1)], 1))
            if r + 1 < 10:
                bonds.append((coord[(r, c)], coord[(r + 1, c)], 1))
    m, ids = build(bonds)
    assert m.atom_count == 100
    assert m.bond_count == 180
    assert m.rings_count == 81

    for s in ids:
        assert m.ring_sizes_of(s) == frozenset({4})
        assert m.macrocycle_of(s) is False
        assert m.in_ring_of(s) == (m.ring_count_of(s) > 0)

    # a grid vertex belongs to one square per quadrant it has both neighbours for:
    # corners 1, other edge vertices 2, interior 4
    for r in range(10):
        for c in range(10):
            on_r_edge = r == 0 or r == 9
            on_c_edge = c == 0 or c == 9
            if on_r_edge and on_c_edge:
                expected = 1
            elif on_r_edge or on_c_edge:
                expected = 2
            else:
                expected = 4
            assert m.ring_count_of(ids[coord[(r, c)]]) == expected

    # two atoms diagonally opposite in one square share exactly that square, and no
    # two atoms in different squares share anything. Both must hold for squares whose
    # bit index lands in the SECOND word -- that is the assertion the transposed
    # index dies on.
    assert m.shares_ring(ids[coord[(9, 9)]], ids[coord[(8, 8)]])
    assert not m.shares_ring(ids[coord[(0, 0)]], ids[coord[(9, 9)]])
    assert not m.shares_ring(ids[coord[(0, 0)]], ids[coord[(0, 3)]])


def test_ring_count_saturates_at_255():
    # K25: 2300 relevant triangles, each atom in C(24, 2) = 276 of them, so the low byte of
    # ring_counts must clamp rather than wrap (276 & 0xff would be 20). `rings` is a basis, so
    # it holds only the circuit rank 300 - 25 + 1 = 276 of them; the per-atom count does not
    # come from the basis.
    bonds = []
    for i in range(25):
        for j in range(i + 1, 25):
            bonds.append((i, j, 1))
    m, ids = build(bonds)
    assert m.bond_count == 300
    assert m.rings_count == 276
    for s in ids:
        assert m.ring_count_of(s) == 255
        assert m.ring_sizes_of(s) == frozenset({3})


def test_a_cyclophane_does_not_materialise_its_exponential_relevant_set():
    # Twelve para-disubstituted benzenes closed into a macrocycle. Each benzene offers two
    # equal-length arms, so there are 2**12 = 4096 distinct 48-membered relevant cycles, all
    # generated by ONE Vismara prototype. Enumerating them is what an earlier revision did:
    # on a 20-benzene version of this motif -- an ordinary macrocyclic aryl sulfone -- it spent
    # 36 seconds and 423 MB. Prototypes are polynomial; the cycles they stand for are not.
    #
    # rings_count is the guard: the basis has 84 - 72 + 1 = 13 rings, and any revision that
    # goes back to storing the relevant set reports 4096 + 12 here instead.
    k = 12
    bonds = []
    for r in range(k):
        base = 6 * r
        for i in range(6):
            bonds.append((base + i, base + (i + 1) % 6, 1))
        bonds.append((base + 3, 6 * ((r + 1) % k), 1))   # para link to the next ring
    m, ids = build(bonds)
    assert m.atom_count == 72
    assert m.bond_count == 84
    assert m.rings_count == 13
    assert sorted(len(x) for x in m.rings) == [6] * 12 + [48]

    for s in ids:
        assert m.ring_sizes_of(s) == frozenset({6})
        assert m.macrocycle_of(s) is True        # the 48-ring is past the exact-size range
        assert m.ring_count_of(s) >= 1
        assert m.in_ring_of(s)
    for i, j, _ in bonds:
        assert m.bond_in_ring(ids[i], ids[j])

    # the macrocycle prototype spans the whole ring system, so atoms six benzenes apart do
    # share a relevant cycle
    assert m.shares_ring(ids[0], ids[36])


def test_a_cage_whose_smallest_ring_set_is_not_a_cycle_basis():
    # A pentacyclic cage with relevant sizes [3, 3, 4, 4, 4, 5, 5, 5, 5] at circuit rank 5.
    # The three smallest rings after the two triangles are the three 4-rings, and taking all
    # five of [3, 3, 4, 4, 4] gives GF(2) rank 4 -- not a basis at all -- while leaving bond
    # (2, 3) covered by no ring. A minimum cycle basis has to reach for one of the 5-rings, which is
    # the whole reason this graph is here: taking the smallest rings by size is not enough.
    bonds = [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1), (5, 6, 1),
             (6, 7, 1), (7, 0, 1), (0, 2, 1), (1, 7, 1), (3, 6, 1), (4, 7, 1)]
    m, ids = build(bonds)
    assert m.rings_count == 5                 # circuit rank is 12 - 8 + 1 = 5
    assert ring_sizes(m) == [3, 3, 4, 4, 5]   # weight 19: the minimum, not [3, 3, 4, 4, 4]
    for i, j, _ in bonds:
        assert m.bond_in_ring(ids[i], ids[j])

    # per-atom sizes come from the relevant prototypes, so they are unaffected by which
    # 5-ring the basis happened to pick
    assert [sorted(m.ring_sizes_of(s)) for s in ids] == [
        [3, 5], [3, 5], [3, 5], [4, 5], [4, 5], [4], [4, 5], [3, 4, 5]]


# The basis property itself, checked against the definition rather than against expected sizes.
# Every helper below reimplements the arithmetic independently of core.

def _mu(adj):
    seen, comps = set(), 0
    for start in adj:
        if start in seen:
            continue
        comps += 1
        stack = [start]
        seen.add(start)
        while stack:
            v = stack.pop()
            for w in adj[v]:
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
    return sum(len(v) for v in adj.values()) // 2 - len(adj) + comps


def _gf2_insert(basis, vector):
    """Reduce `vector` against `basis`; append and return True if it raised the rank."""
    for b in basis:
        vector = min(vector, vector ^ b)
    if vector:
        basis.append(vector)
        basis.sort(reverse=True)
        return True
    return False


def _all_cycles(adj):
    """Every simple cycle of `adj` exactly once, as a node list. Exponential; small graphs only."""
    order = {v: i for i, v in enumerate(sorted(adj))}
    out = []
    for start in adj:
        stack = [(start, [start], {start})]
        while stack:
            v, path, seen = stack.pop()
            for w in adj[v]:
                if w == start:
                    # one canonical rotation and direction per cycle
                    if len(path) > 2 and order[path[1]] < order[path[-1]] \
                            and order[start] == min(order[x] for x in path):
                        out.append(list(path))
                elif w not in seen and order[w] > order[start]:
                    stack.append((w, path + [w], seen | {w}))
    return out


def assert_minimum_cycle_basis(bonds, brute_force=True, elements=None):
    m, ids = build(bonds, elements=elements)
    position = {s: i for i, s in enumerate(ids)}
    adj = {}
    for i, j, _ in bonds:
        adj.setdefault(i, set()).add(j)
        adj.setdefault(j, set()).add(i)
    edge_id = {}
    for i, j, _ in bonds:
        edge_id[frozenset((i, j))] = len(edge_id)

    def vector(ring):
        out = 0
        for k in range(len(ring)):
            out |= 1 << edge_id[frozenset((ring[k], ring[k - 1]))]
        return out

    rings = [[position[s] for s in r] for r in m.rings]
    mu = _mu(adj)

    assert m.rings_count == mu
    assert len(rings) == mu, 'basis has the wrong number of members'

    for ring in rings:                     # each one is a real simple cycle of the graph
        assert len(set(ring)) == len(ring) >= 3
        for k in range(len(ring)):
            assert ring[k - 1] in adj[ring[k]]

    basis = []                             # GF(2) independent over the bond set, hence a basis
    for ring in rings:
        assert _gf2_insert(basis, vector(ring)), 'basis members are linearly dependent'

    covered = set()                        # and spanning: no cycle bond is left uncovered
    for ring in rings:
        for k in range(len(ring)):
            covered.add(frozenset((ring[k], ring[k - 1])))
    for i, j, _ in bonds:
        if m.bond_in_ring(ids[i], ids[j]):
            assert frozenset((i, j)) in covered, f'ring bond {i}-{j} is in no basis member'

    if brute_force:
        # greedy over every cycle shortest-first is a minimum cycle basis, so its weight is
        # the bound to beat
        reference, weight = [], 0
        for cycle in sorted(_all_cycles(adj), key=len):
            if _gf2_insert(reference, vector(cycle)):
                weight += len(cycle)
        assert len(reference) == mu
        assert sum(len(r) for r in rings) == weight, 'basis is independent but not minimum'
    return m, ids


def _cycle(n, offset=0):
    return [(offset + i, offset + (i + 1) % n, 1) for i in range(n)]


def _complete(n):
    return [(i, j, 1) for i in range(n) for j in range(i + 1, n)]


def _bipartite(a, b):
    return [(i, a + j, 1) for i in range(a) for j in range(b)]


BASIS_CASES = {
    'benzene': _cycle(6),
    'cyclopropane': _cycle(3),
    'macrocycle-30': _cycle(30),
    'spiro[4.4]': _cycle(5) + [(0, 5, 1), (5, 6, 1), (6, 7, 1), (7, 8, 1), (8, 0, 1)],
    'two-rings-one-bridge': _cycle(5) + _cycle(6, 10) + [(0, 10, 1)],
    'disjoint-benzenes': _cycle(6) + _cycle(6, 6),
    'benzene-with-a-tail': _cycle(6) + [(0, 100, 1), (100, 101, 1), (101, 102, 1)],
    'norbornane': [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1), (5, 0, 1),
                   (2, 6, 1), (6, 5, 1)],
    'bicyclo[2.2.2]octane': [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1), (5, 0, 1),
                             (0, 6, 1), (6, 7, 1), (7, 3, 1)],
    'adamantane': [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1), (5, 0, 1), (1, 6, 1),
                   (6, 7, 1), (7, 8, 1), (8, 3, 1), (5, 9, 1), (9, 7, 1)],
    'tetrahedrane': _complete(4),
    'prismane': [(0, 1, 1), (1, 2, 1), (2, 0, 1), (3, 4, 1), (4, 5, 1), (5, 3, 1),
                 (0, 3, 1), (1, 4, 1), (2, 5, 1)],
    'cubane': [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 0, 1), (4, 5, 1), (5, 6, 1), (6, 7, 1),
               (7, 4, 1), (0, 4, 1), (1, 5, 1), (2, 6, 1), (3, 7, 1)],
    'pentacyclic-cage': [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1), (5, 6, 1),
                         (6, 7, 1), (7, 0, 1), (0, 2, 1), (1, 7, 1), (3, 6, 1), (4, 7, 1)],
    'K5': _complete(5),
    'K6': _complete(6),
    'K3,3': _bipartite(3, 3),
    'K4,4': _bipartite(4, 4),
    'petersen': [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 0, 1),
                 (5, 7, 1), (7, 9, 1), (9, 6, 1), (6, 8, 1), (8, 5, 1),
                 (0, 5, 1), (1, 6, 1), (2, 7, 1), (3, 8, 1), (4, 9, 1)],
}


@pytest.mark.parametrize('name', sorted(BASIS_CASES))
def test_rings_are_a_minimum_cycle_basis(name):
    assert_minimum_cycle_basis(BASIS_CASES[name])


def test_dodecahedron_is_a_basis():
    # 20 vertices, rank 11: brute-force cycle enumeration is too slow, the basis checks are not
    assert_minimum_cycle_basis(
        [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 0, 1),
         (0, 5, 1), (1, 6, 1), (2, 7, 1), (3, 8, 1), (4, 9, 1),
         (5, 10, 1), (5, 14, 1), (6, 10, 1), (6, 11, 1), (7, 11, 1), (7, 12, 1), (8, 12, 1),
         (8, 13, 1), (9, 13, 1), (9, 14, 1),
         (10, 15, 1), (11, 16, 1), (12, 17, 1), (13, 18, 1), (14, 19, 1),
         (15, 16, 1), (16, 17, 1), (17, 18, 1), (18, 19, 1), (19, 15, 1)],
        brute_force=False)


def test_a_wheel_loses_its_rim_to_the_minimum_basis():
    # Ferrocene's shape with the Fe-C bonds spelled as ordinary single bonds. Each Fe-Cp
    # fragment is a wheel: five triangles weigh 15, four triangles plus the Cp five-ring weigh
    # 17. So the rim is in no minimum cycle basis -- and being in none, it is not a relevant
    # cycle either, so it is absent from the per-atom descriptors too. No ring algorithm can
    # recover it. This is the whole reason dative bonds have to be excluded before perception.
    bonds = _cycle(5) + _cycle(5, 5) + [(10, i, 1) for i in range(10)]
    m, ids = assert_minimum_cycle_basis(bonds, brute_force=False, elements=[6] * 10 + [26])
    assert sorted(len(r) for r in m.rings) == [3] * 10
    assert m.ring_sizes_of(ids[0]) == frozenset({3})


def test_ferrocene_is_two_five_rings_with_the_iron_outside_them():
    # The same skeleton with the Fe-C bonds as order 8, which is what they are. mark_bridges
    # does not admit them, so the rim survives as the basis and the iron is in no ring at all.
    m, ids = build(_cycle(5) + _cycle(5, 5) + [(10, i, 8) for i in range(10)],
                   elements=[6] * 10 + [26])
    assert m.rings_count == 2
    assert sorted(len(r) for r in m.rings) == [5, 5]
    assert m.sssr == m.rings

    fe = ids[10]
    assert m.in_ring_of(fe) is False
    assert m.ring_count_of(fe) == 0
    assert m.ring_sizes_of(fe) == frozenset()
    for i in range(10):
        assert m.ring_sizes_of(ids[i]) == frozenset({5})
        assert m.in_ring_of(ids[i]) is True
        assert m.bond_in_ring(fe, ids[i]) is False
    for i in range(5):
        assert m.bond_in_ring(ids[i], ids[(i + 1) % 5]) is True
    assert m.shares_ring(ids[0], ids[5]) is False   # the two Cp rings share nothing


def test_a_dative_bond_cannot_close_a_ring():
    # cyclohexane where one bond is dative: no ring, and the ring flags say so everywhere
    bonds = [(i, (i + 1) % 6, 1) for i in range(6)]
    bonds[5] = (5, 0, 8)
    m, ids = build(bonds)
    assert m.rings_count == 0
    assert m.rings == []
    for i, j, _ in bonds:
        assert m.bond_in_ring(ids[i], ids[j]) is False
    for s in ids:
        assert m.in_ring_of(s) is False
        assert m.ring_sizes_of(s) == frozenset()
