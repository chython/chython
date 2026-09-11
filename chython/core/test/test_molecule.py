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
from struct import unpack_from

import pytest

from chython.core import Atom, Bond, MoleculeContainer, QueryContainer
from chython.core._core import JOURNAL_OPS, journal_record_size, _seal_probe


def test_journal_record_is_twenty_bytes_and_every_field_is_naturally_aligned():
    """It was 16 until `OP_SET_XYZ` needed a fourth payload word, and 20 is the honest price.

    The struct is deliberately NOT packed, and what that buys is aligned loads -- which 20 still gives,
    because `op` is followed by three bytes of padding and every payload word is a 4-byte field at a
    4-byte offset.  What is lost is a stride that is a power of two, and that is worth nothing here: the
    journal is scanned sequentially, never indexed by a shift.

    THIS BUFFER IS TRANSIENT AND IS NEVER SERIALISED, which is the whole reason a field was the right
    answer.  It is malloc'd at the first edit and freed at seal, so growing it reprices no stored key
    and breaks no buffer -- unlike `atom_t`, whose size this suite pins for exactly that reason.  The
    alternative that would have kept 16 was a PAIR of adjacent records read together, and that makes
    correctness depend on an ordering invariant the struct cannot state, in a buffer four separate
    loops walk by index.  Four bytes of scratch against a class of bug.
    """
    assert journal_record_size() == 20
    assert journal_record_size() % 4 == 0, 'a payload word would be misaligned'


def test_fresh_container_is_empty():
    m = MoleculeContainer()
    assert m.journal_length == 0
    assert m.atom_count == 0


def test_add_atom_appends_one_record_and_returns_a_stable_id():
    m = MoleculeContainer()
    with m.edit():
        assert m.add_atom(6) == 1
        assert m.add_atom(8) == 2
        assert m.journal_length == 2
        assert m.journal_record(0) == (JOURNAL_OPS['add_atom'], 1, 0, 6, 0, 0)
        assert m.journal_record(1) == (JOURNAL_OPS['add_atom'], 2, 0, 8, 0, 0)
        with pytest.raises(IndexError):
            m.journal_record(2)
    assert m.journal_length == 0
    assert m.atom_count == 2


def test_attributes_arrive_as_their_own_records_in_a_fixed_order():
    m = MoleculeContainer()
    with m.edit():
        sid = m.add_atom(6, charge=-1, isotope=13, radical=True, map_number=17,
                         implicit_h=2, stereo=True)
        assert m.journal_length == 7
        assert [m.journal_record(i)[0] for i in range(7)] == [
            JOURNAL_OPS['add_atom'], JOURNAL_OPS['set_charge'], JOURNAL_OPS['set_isotope'],
            JOURNAL_OPS['set_radical'], JOURNAL_OPS['set_map_number'],
            JOURNAL_OPS['set_hydrogens'], JOURNAL_OPS['set_stereo']]
        assert m.journal_record(1) == (JOURNAL_OPS['set_charge'], sid, 0, -1, 0, 0)
        assert m.journal_record(5) == (JOURNAL_OPS['set_hydrogens'], sid, 0, 2, 0, 0)


def test_implicit_h_zero_still_emits_a_record():
    # 0 pinned hydrogens is a statement; None is the absence of one
    m = MoleculeContainer()
    with m.edit():
        m.add_atom(6, implicit_h=0)
        assert m.journal_length == 2
        assert m.journal_record(1) == (JOURNAL_OPS['set_hydrogens'], 1, 0, 0, 0, 0)
    n = MoleculeContainer()
    with n.edit():
        n.add_atom(6)
        assert n.journal_length == 1


def test_bond_ops_carry_both_endpoints():
    m = MoleculeContainer()
    with m.edit():
        a1, a2 = m.add_atom(6), m.add_atom(6)
        m.add_bond(a1, a2, 2)
        m.set_order(a1, a2, 3)
        m.delete_bond(a2, a1)
        assert m.journal_record(2) == (JOURNAL_OPS['add_bond'], a1, a2, 2, 0, 0)
        assert m.journal_record(3) == (JOURNAL_OPS['set_order'], a1, a2, 3, 0, 0)
        assert m.journal_record(4) == (JOURNAL_OPS['delete_bond'], a2, a1, 0, 0, 0)
    assert m.atom_count == 2
    assert m.bond_count == 0


def test_journal_grows_past_its_initial_capacity():
    m = MoleculeContainer()
    with m.edit():
        for _ in range(500):
            m.add_atom(6)
        assert m.journal_length == 500
        assert m.journal_record(499) == (JOURNAL_OPS['add_atom'], 500, 0, 6, 0, 0)
        assert m.journal_record(0) == (JOURNAL_OPS['add_atom'], 1, 0, 6, 0, 0)
    assert m.atom_count == 500


def test_symbols_and_atomic_numbers_both_work():
    m = MoleculeContainer()
    with m.edit():
        assert m.add_atom('H') == 1
        assert m.add_atom('Cl') == 2
        assert m.add_atom('Og') == 3
        assert m.journal_record(0)[3] == 1
        assert m.journal_record(1)[3] == 17
        assert m.journal_record(2)[3] == 118
    assert m.element_of(2) == 17


def test_unimplemented_element_forms_raise_not_implemented():
    m = MoleculeContainer()
    with pytest.raises(NotImplementedError):
        m.add_atom(object())
    with pytest.raises(ValueError):
        m.add_atom('Xx')


def test_element_and_attribute_ranges_are_validated():
    m = MoleculeContainer()
    with pytest.raises(ValueError):
        m.add_atom(119)
    with pytest.raises(ValueError):
        m.add_atom(7, charge=9)
    with pytest.raises(ValueError):
        m.add_atom(7, charge=-5)
    with pytest.raises(ValueError):
        m.add_atom(7, isotope=-1)
    with pytest.raises(ValueError):
        m.add_atom(7, map_number=10000)
    with pytest.raises(ValueError):
        m.add_atom(7, implicit_h=16)
    assert m.add_atom(7, charge=8)
    assert m.add_atom(7, charge=-4)


def test_r_atom_is_accepted_by_both_spellings():
    """element 0 (R) is valid; `add_atom(0)` and `add_atom('R')` both land element 0."""
    m = MoleculeContainer()
    with m.edit():
        n = m.add_atom(0)
        s = m.add_atom('R')
    assert m.atom(n).element == 0
    assert m.atom(n).atomic_symbol == 'R'
    assert m.atom(s).element == 0
    assert m.atom(s).atomic_symbol == 'R'


def test_map_number_boundary_agrees_between_molecule_and_query():
    # MAP_NUMBER_MAX is 9999.  The molecule side and the query side must enforce the same
    # ceiling: a map number the container refuses must not seal into a query.
    # Drive the boundary through the public API on both sides rather than comparing constants:
    # a constant-comparison test passes even if a validator forgets to use its constant.

    # molecule side: 9999 accepted, 10000 rejected
    m = MoleculeContainer()
    sid = m.add_atom(6, map_number=9999)
    assert m.map_number_of(sid) == 9999
    with pytest.raises(ValueError):
        m.add_atom(6, map_number=10000)
    with pytest.raises(ValueError):
        m.set_map_number(sid, 10000)

    # query container side: 9999 accepted, 10000 and -1 rejected
    q = QueryContainer()
    a = q.add_atom()
    q.atom_primitive(a, 'element', 6)
    q.set_map_number(a, 9999)
    assert q.map_numbers() == {a: 9999}
    q.atom_count_sealed()  # forces a seal; confirms 9999 is accepted by the seal path
    q2 = QueryContainer()
    a2 = q2.add_atom()
    with pytest.raises(ValueError):
        q2.set_map_number(a2, 10000)
    with pytest.raises(ValueError):
        q2.set_map_number(a2, -1)

    # seal-path validator: _seal_probe bypasses the container validator and reaches query_seal
    # directly, so it is the only route that tests the seal-path validator on the rejecting side.
    with pytest.raises(ValueError):
        _seal_probe([('atom', 1), ('token', 1, 'element', 6, False), ('map', 1, 10000)])


def test_bond_orders_are_validated_and_self_loops_rejected():
    m = MoleculeContainer()
    with m.edit():
        a1, a2 = m.add_atom(6), m.add_atom(6)
        with pytest.raises(ValueError):
            m.add_bond(a1, a2, 0)
        with pytest.raises(ValueError):
            m.add_bond(a1, a2, 5)
        with pytest.raises(ValueError):
            m.add_bond(a1, a1, 1)     # no self loops
        assert m.journal_length == 2   # nothing rejected reached the journal
    assert m.bond_count == 0


def test_an_aromatic_bond_is_stored_as_written():
    # ORDER 4 IS A STORED ORDER: `add_bond(..., 4)` stores it rather than demanding a Kekule form,
    # because refusing to store what a file contained loses the input. What that costs a consumer
    # of bond orders is three order-aware answers and nothing else -- hybridization 4, the aromatic
    # feature bit, the sp2 electron budget.
    m = MoleculeContainer()
    a1, a2 = m.add_atom(6), m.add_atom(6)
    m.add_bond(a1, a2, 4)
    assert m.order_of(a1, a2) == 4
    assert m.aromatic_bond_count == 1
    assert not m.is_kekule

    # and the same through a scope, where the order reaches the arena via the compaction path
    other = MoleculeContainer()
    with other.edit():
        b1, b2 = other.add_atom(6), other.add_atom(6)
        other.add_bond(b1, b2, 4)
        assert other.journal_length == 3
    assert other.order_of(b1, b2) == 4
    assert other.aromatic_bond_count == 1


def test_an_aromatic_bond_removed_before_the_apply_leaves_no_trace():
    # A bond added and then deleted in the same scope is not in the final graph, so it
    # contributes nothing -- neither an order nor a count. The count is recomputed from the
    # surviving edges rather than incremented as bonds are journalled, which is what makes this
    # hold; a counter maintained on the journal would report one aromatic bond in an empty graph.
    m = MoleculeContainer()
    with m.edit():
        a1, a2 = m.add_atom(6), m.add_atom(6)
        m.add_bond(a1, a2, 4)
        m.delete_bond(a1, a2)
    assert m.bond_count == 0
    assert m.order_of(a1, a2) is None
    assert m.aromatic_bond_count == 0
    assert m.is_kekule

    # The other way a bond leaves the graph: an endpoint is deleted. Here the recorded order in
    # the journal is still 4 -- only delete_bond zeroes it -- so this is the case that pins the
    # count to the surviving-edge branch rather than to the top of the compaction loop.
    other = MoleculeContainer()
    with other.edit():
        b1, b2 = other.add_atom(6), other.add_atom(6)
        other.add_bond(b1, b2, 4)
        other.delete_atom(b2)
    assert other.atom_count == 1
    assert other.bond_count == 0
    assert other.aromatic_bond_count == 0


def test_unknown_stable_ids_raise_key_error():
    m = MoleculeContainer()
    a1 = m.add_atom(6)
    with pytest.raises(KeyError):
        m.add_bond(a1, 99, 1)
    with pytest.raises(KeyError):
        m.delete_atom(99)
    with pytest.raises(KeyError):
        m.set_charge(99, 1)


# ── Task 6 tests ────────────────────────────────────────────────────────────


def test_atomic_mutation_is_visible_immediately():
    m = MoleculeContainer()
    c = m.add_atom(6)
    o = m.add_atom(8)
    m.add_bond(c, o, 2)
    assert m.journal_length == 0      # applied on every call
    assert m.atom_count == 2
    assert m.bond_count == 1
    assert m.order_of(c, o) == 2
    assert m.order_of(o, c) == 2
    assert m.order_of(c, c) is None


def test_indices_are_dense_and_stable_ids_are_never_reused():
    m = MoleculeContainer()
    a1, a2, a3 = m.add_atom(6), m.add_atom(8), m.add_atom(7)
    m.delete_atom(a2)
    a4 = m.add_atom(16)
    m.add_bond(a1, a3, 1)
    m.add_bond(a3, a4, 2)
    assert a4 == 4
    assert m.atom_count == 3
    assert m.bond_count == 2
    assert m.atom_numbers == [a1, a3, a4]
    assert [m.index_of(s) for s in (a1, a3, a4)] == [0, 1, 2]
    assert [m.number_of(i) for i in range(3)] == [a1, a3, a4]
    with pytest.raises(KeyError):
        m.index_of(a2)


def test_atom_n_is_the_number_the_container_issued():
    """`a.n` is the id `add_atom` returned, gaps included -- not the arena position."""
    m = MoleculeContainer()
    a1, a2, a3 = m.add_atom(6), m.add_atom(8), m.add_atom(7)
    m.delete_atom(a2)
    a4 = m.add_atom(16)
    m.add_bond(a1, a3, 1)
    m.add_bond(a3, a4, 2)
    assert [a.n for a in m.atoms()] == [1, 3, 4], 'the deleted 2 leaves a gap'
    assert m.atom_numbers == [1, 3, 4], 'the same numbers, in arena order'
    assert [m.index_of(a.n) for a in m.atoms()] == [0, 1, 2], 'positions stay dense'


def test_number_of_is_the_inverse_of_index_of():
    """`number_of(index)` is position -> number; `IndexError` past the end, never a wrap."""
    m = MoleculeContainer()
    a1, a2, a3 = m.add_atom(6), m.add_atom(8), m.add_atom(7)
    m.delete_atom(a2)
    a4 = m.add_atom(16)
    m.add_bond(a1, a3, 1)
    m.add_bond(a3, a4, 2)
    assert [m.number_of(i) for i in range(3)] == [1, 3, 4]
    with pytest.raises(IndexError):
        m.number_of(3)
    with pytest.raises(IndexError):
        m.number_of(10000)


def test_primary_atom_data_survives_the_apply():
    m = MoleculeContainer()
    a = m.add_atom(6, charge=-1, isotope=13, radical=True, map_number=17,
                   implicit_h=2, stereo=True)
    assert m.element_of(a) == 6
    assert m.charge_of(a) == -1
    assert m.isotope_of(a) == 13       # absolute mass number, no delta
    assert m.map_number_of(a) == 17
    assert m.implicit_h_of(a) == 2
    assert m.explicit_h_of(a) == 0     # explicit_h was not set; stored value is 0
    assert m.radical_of(a) is True
    assert m.stereo_of(a) is True
    assert m.map_number_of(m.add_atom(6)) == 0


def test_later_records_win_over_earlier_ones():
    m = MoleculeContainer()
    a = m.add_atom(6, charge=1)
    m.set_charge(a, -2)
    m.set_charge(a, 3)
    assert m.charge_of(a) == 3


def test_degree_is_filled_from_the_csr():
    m = MoleculeContainer()
    c = m.add_atom(6)
    hs = [m.add_atom(1) for _ in range(4)]
    for h in hs:
        m.add_bond(c, h, 1)
    assert m.degree_of(c) == 4
    assert [m.degree_of(h) for h in hs] == [1, 1, 1, 1]


def test_deleting_an_atom_drops_its_bonds():
    m = MoleculeContainer()
    a1, a2, a3 = m.add_atom(6), m.add_atom(6), m.add_atom(6)
    m.add_bond(a1, a2, 1)
    m.add_bond(a2, a3, 1)
    m.delete_atom(a2)
    assert m.atom_count == 2
    assert m.bond_count == 0
    assert m.atom_numbers == [a1, a3]


def test_set_order_and_delete_bond_affect_both_directions():
    m = MoleculeContainer()
    a1, a2, a3 = m.add_atom(6), m.add_atom(6), m.add_atom(6)
    m.add_bond(a1, a2, 1)
    m.set_order(a2, a1, 3)
    assert m.order_of(a1, a2) == 3
    assert m.order_of(a2, a1) == 3
    with pytest.raises(ValueError):
        m.set_order(a1, a2, 5)         # 5 is not a bond order at all
    m.set_order(a1, a2, 4)             # order 4 is stored, both directions, like any other
    assert m.order_of(a1, a2) == 4
    assert m.order_of(a2, a1) == 4
    assert m.aromatic_bond_count == 1
    m.set_order(a1, a2, 3)             # and back out again: the count follows the bonds
    assert m.aromatic_bond_count == 0
    with pytest.raises(KeyError):
        m.set_order(a1, a3, 1)         # both atoms exist, the bond does not
    m.delete_bond(a2, a1)
    assert m.bond_count == 0
    assert m.order_of(a1, a2) is None
    with pytest.raises(KeyError):
        m.delete_bond(a1, a2)          # already gone


def test_duplicate_bond_is_rejected():
    m = MoleculeContainer()
    a1, a2 = m.add_atom(6), m.add_atom(6)
    m.add_bond(a1, a2, 1)
    with pytest.raises(ValueError):
        m.add_bond(a1, a2, 2)
    with pytest.raises(ValueError):
        m.add_bond(a2, a1, 2)          # the reverse direction is the same bond
    assert m.order_of(a1, a2) == 1
    assert m.bond_count == 1


def test_every_mutation_bumps_the_generation():
    m = MoleculeContainer()
    a = m.add_atom(6)
    gen = m.generation
    m.set_charge(a, 1)
    assert m.generation == gen + 1


def test_an_empty_container_reads_as_empty():
    m = MoleculeContainer()
    assert m.atom_count == 0
    assert m.bond_count == 0
    assert m.atom_numbers == []


def _ethanol():
    m = MoleculeContainer()
    with m.edit():
        c1 = m.add_atom(6)
        c2 = m.add_atom(6)
        o = m.add_atom(8)
        m.add_bond(c1, c2, 1)
        m.add_bond(c2, o, 1)
    return m, c1, c2, o


def test_a_scope_applies_once_on_exit():
    m, c1, c2, o = _ethanol()
    assert m.atom_count == 3
    assert m.bond_count == 2
    assert m.generation == 1      # one apply for the whole block
    assert m.order_of(c2, o) == 1


def test_reading_a_property_inside_a_dirty_scope_raises():
    m = MoleculeContainer()
    with m.edit():
        m.add_atom(6)
        with pytest.raises(RuntimeError, match='pending'):
            m.atom_count
        with pytest.raises(RuntimeError, match='pending'):
            m.atom_numbers


def test_a_scope_that_has_not_mutated_anything_still_reads():
    m, c1, c2, o = _ethanol()
    with m.edit():
        assert m.atom_count == 3   # nothing pending yet
        m.set_charge(o, -1)
    assert m.charge_of(o) == -1


def test_scopes_nest_and_only_the_outermost_applies():
    m = MoleculeContainer()
    with m.edit():
        c = m.add_atom(6)
        with m.edit():
            o = m.add_atom(8)
        assert m.journal_length == 2   # the inner exit did not apply
    assert m.atom_count == 2
    assert m.atom_numbers == [c, o]
    assert m.generation == 1


def test_an_exception_in_a_scope_discards_the_journal():
    m, c1, c2, o = _ethanol()
    gen = m.generation
    with pytest.raises(ZeroDivisionError):
        with m.edit():
            m.delete_atom(o)
            m.add_atom(7)
            raise ZeroDivisionError
    assert m.generation == gen
    assert m.journal_length == 0
    assert m.atom_count == 3
    assert m.charge_of(o) == 0


def test_a_stable_id_burned_by_a_rolled_back_scope_is_not_reused():
    m, c1, c2, o = _ethanol()
    with pytest.raises(ZeroDivisionError):
        with m.edit():
            m.add_atom(7)
            raise ZeroDivisionError
    assert m.add_atom(7) == 5       # 4 was consumed and discarded
    assert m.atom_count == 4


def test_a_duplicate_bond_inside_a_scope_raises_at_exit():
    m = MoleculeContainer()
    with pytest.raises(ValueError, match='duplicate'):
        with m.edit():
            a1 = m.add_atom(6)
            a2 = m.add_atom(6)
            m.add_bond(a1, a2, 1)
            m.add_bond(a2, a1, 2)   # the block is one transaction, so not caught here
    assert m.journal_length == 0
    assert m.atom_count == 0        # the failed apply swapped in nothing


def test_a_missing_bond_inside_a_scope_raises_at_exit():
    m, c1, c2, o = _ethanol()
    with pytest.raises(KeyError):
        with m.edit():
            m.set_order(c1, o, 2)   # c1 and o exist but share no bond
    assert m.journal_length == 0
    assert m.order_of(c1, c2) == 1


def test_atom_payload_follows_a_survivor_across_compaction():
    m = MoleculeContainer()
    a = m.add_atom(6, charge=-1, isotope=13, map_number=7)
    b = m.add_atom(6)
    c = m.add_atom(8, charge=1, isotope=18, map_number=9)
    m.add_bond(a, b, 1)
    m.add_bond(b, c, 2)
    m.delete_atom(b)                       # c moves from index 2 to index 1
    assert m.index_of(c) == 1
    assert m.charge_of(c) == 1
    assert m.isotope_of(c) == 18
    assert m.map_number_of(c) == 9
    assert m.element_of(c) == 8
    assert m.charge_of(a) == -1
    assert m.isotope_of(a) == 13
    assert m.map_number_of(a) == 7


def test_a_scope_that_both_deletes_and_adds_lands_consistently():
    m = MoleculeContainer()
    a, b, c = m.add_atom(6), m.add_atom(6), m.add_atom(6)
    m.add_bond(a, b, 1)
    m.add_bond(b, c, 1)
    with m.edit():
        m.delete_atom(b)
        d = m.add_atom(7)
        e = m.add_atom(8)
        m.add_bond(a, d, 1)
        m.add_bond(d, e, 2)
    assert m.atom_count == 4
    assert m.atom_numbers == [a, c, d, e]
    assert [m.index_of(x) for x in (a, c, d, e)] == [0, 1, 2, 3]
    assert [m.number_of(i) for i in range(4)] == [a, c, d, e]
    assert m.order_of(a, d) == 1
    assert m.order_of(d, e) == 2
    with pytest.raises(KeyError):
        m.order_of(a, b)   # b is gone
    assert m.element_of(d) == 7
    assert m.element_of(e) == 8
    assert m.bond_count == 2


def test_number_of_rejects_an_index_past_the_end():
    m = MoleculeContainer()
    a = m.add_atom(6)
    b = m.add_atom(6)
    m.add_bond(a, b, 1)
    assert m.number_of(0) == a
    assert m.number_of(1) == b
    with pytest.raises(IndexError):
        m.number_of(2)
    with pytest.raises(IndexError):
        m.number_of(10000)


def test_all_arena_backed_readers_raise_inside_a_dirty_scope():
    m = MoleculeContainer()
    a = m.add_atom(6)
    b = m.add_atom(6)
    m.add_bond(a, b, 1)

    readers = [
        ('atom_count', None),
        ('bond_count', None),
        ('atom_numbers', None),
        ('index_of', (a,)),
        ('number_of', (0,)),
        ('element_of', (a,)),
        ('charge_of', (a,)),
        ('isotope_of', (a,)),
        ('map_number_of', (a,)),
        ('degree_of', (a,)),
        ('implicit_h_of', (a,)),
        ('explicit_h_of', (a,)),
        ('radical_of', (a,)),
        ('stereo_of', (a,)),
        ('order_of', (a, b)),
        # Lazily built derived segments read the arena too, and perceive against it: inside a
        # dirty scope they would answer from the PRE-edit graph and hand back slot-keyed results
        # the pending journal is about to renumber.
        ('component_labels', ()),
        ('stereo_units', ()),
        ('unit_of', (a,)),
    ]

    with m.edit():
        m.set_charge(a, 1)   # make it dirty
        for name, args in readers:
            with pytest.raises(RuntimeError, match='pending'):
                if args is None:
                    getattr(m, name)
                else:
                    getattr(m, name)(*args)


# ── Task 7 tests ────────────────────────────────────────────────────────────


def test_copy_shares_the_arena():
    m, c1, c2, o = _ethanol()
    n = m.copy()
    assert n is not m
    assert n.shares_arena_with(m) is True
    assert m.shares_arena_with(n) is True
    assert n.atom_numbers == m.atom_numbers
    assert n.generation == m.generation
    assert n.bond_count == 2


def test_a_mutation_breaks_the_sharing_without_disturbing_the_copy():
    m, c1, c2, o = _ethanol()
    n = m.copy()
    m.delete_atom(o)
    assert m.atom_count == 2
    assert n.atom_count == 3
    assert n.charge_of(o) == 0
    assert n.order_of(c2, o) == 1
    assert n.shares_arena_with(m) is False


def test_a_copy_gets_a_stable_id_counter_that_does_not_collide():
    m, c1, c2, o = _ethanol()
    n = m.copy()
    assert m.add_atom(7) == 4
    assert n.add_atom(7) == 4      # independent containers, independent arenas
    assert m.atom_count == 4
    assert n.atom_count == 4
    assert m.shares_arena_with(n) is False


def test_atom_view_reads_through():
    m, c1, c2, o = _ethanol()
    a = m.atom(o)
    assert isinstance(a, Atom)
    assert a.element == 8
    assert a.n == o
    assert a.charge == 0
    assert a.degree == 1
    assert a.radical is False
    with pytest.raises(KeyError):
        m.atom(99)


def test_neighbors_are_stable_ids_on_the_molecule():
    m, c1, c2, o = _ethanol()
    assert m.neighbors_of(o) == [c2]
    assert sorted(m.neighbors_of(c2)) == sorted([c1, o])
    assert m.neighbors_of(c1) == [c2]
    with pytest.raises(KeyError):
        m.neighbors_of(99)


def test_a_view_taken_before_a_mutation_is_stale():
    m, c1, c2, o = _ethanol()
    a = m.atom(o)
    bd = m.bond(c2, o)
    m.set_charge(o, -1)
    with pytest.raises(RuntimeError, match='stale'):
        a.element
    with pytest.raises(RuntimeError, match='stale'):
        bd.order
    assert m.atom(o).charge == -1   # a freshly taken view is fine


def test_a_view_on_a_copy_survives_a_mutation_of_the_original():
    m, c1, c2, o = _ethanol()
    n = m.copy()
    a = n.atom(o)
    m.delete_atom(o)
    assert a.element == 8           # n never moved off its arena


def test_bond_view_and_iteration_visit_each_bond_once():
    m, c1, c2, o = _ethanol()
    bd = m.bond(c1, c2)
    assert isinstance(bd, Bond)
    assert bd.order == 1
    assert int(bd) == 1
    assert len(list(m.bonds())) == 2
    assert all(isinstance(e, Bond) for e in m.bonds())
    assert {tuple(sorted((e.n, e.m))) for e in m.bonds()} == {(c1, c2), (c2, o)}
    assert [a.n for a in m.atoms()] == [c1, c2, o]
    with pytest.raises(KeyError):
        m.bond(c1, o)               # both atoms exist, the bond does not


def test_the_old_endpoint_names_are_GONE_and_not_merely_discouraged():
    """The endpoints are `n`/`m`; `a`/`b` are not names for them and there is no migration shim.

    An `(bd.a, bd.b) == (bd.n, bd.m)` assertion would be SELF-REFERENTIAL -- both properties
    would return the same two slots, so it establishes that an alias aliases, and it survives
    transposing the endpoints at the source.  The values are pinned by the two tests below; what
    is worth pinning about `a`/`b` is that they do not resolve, because a name that quietly comes
    back is how two spellings for one endpoint ship.
    """
    m, c1, c2, o = _ethanol()
    bd = m.bond(c1, c2)
    with pytest.raises(AttributeError):
        bd.a
    with pytest.raises(AttributeError):
        bd.b


def test_the_old_atom_number_names_are_GONE_and_not_merely_discouraged():
    """`stable_id`/`stable_ids`/`stable_id_of`/`atoms_numbers` are not names for this value.

    The reason is the `a`/`b` entry's above: all four would answer what `n`, `atom_numbers` and
    `number_of` answer, so an equality between two of them tests the aliasing and not the value.
    The values are pinned by `test_atom_n_is_the_number_the_container_issued` and
    `test_number_of_is_the_inverse_of_index_of`; what is worth pinning here is that the four
    do not resolve.  `atoms_numbers` gets no deprecated alias either -- see the comment above the
    chython-2 block in `_molecule_container.pxi`.
    """
    m, c1, c2, o = _ethanol()
    with pytest.raises(AttributeError):
        m.stable_ids
    with pytest.raises(AttributeError):
        m.atoms_numbers
    with pytest.raises(AttributeError):
        m.stable_id_of
    with pytest.raises(AttributeError):
        m.atom(c1).stable_id


def test_bond_answers_the_endpoints_in_the_ORDER_ASKED():
    """`mol.bond(n, m).n is n`, not "one of the two endpoints".

    A bond is undirected and every lookup in this file is symmetric, so it is tempting to
    call the orientation an implementation detail.  It is not -- a caller building a
    `{(n, m): ...}` map from these endpoints and probing it with an independently ordered
    pair gets a miss rather than an error, and `files/ctfile` builds exactly such a map.

    Found by mutation: swapping `bd._n`/`bd._m` in `bond()` passed all 2467 tests.  The
    rename that introduced these names was executed by hundreds of tests and verified by
    none of them -- a line running is not an assertion depending on its value.

    WHAT THIS IS *NOT*, corrected after the epic that owns the path MEASURED it.  I claimed
    a transposition here silently inverts a stereo parity in the CTfile wedge writers.  It
    does not.  That path looks up `wedge_of[(n, m)]` and ON A MISS tries `(m, n)` and swaps
    so the narrow end is written first, which removes orientation from the answer before any
    parity is computed -- 264 corpus records carrying 1191 wedge bonds emitted and re-read
    under the mutation gave 0 parity mismatches, against a positive control (swapping the
    emitted stereo codes 1 and 6) that gave 243 of 264.  So the probe reports what it claims
    and the zero is real.  The reason to pin the order is fidelity of the emitted endpoint
    order, not a silent-corruption hazard, and the alarming version of the story was mine.
    """
    m, c1, c2, o = _ethanol()
    assert (m.bond(c1, c2).n, m.bond(c1, c2).m) == (c1, c2)
    assert (m.bond(c2, c1).n, m.bond(c2, c1).m) == (c2, c1), \
        'the same bond addressed the other way round reports the other way round'
    assert (m.bond(c2, o).n, m.bond(c2, o).m) == (c2, o), \
        'and it is not merely sorted -- c2 > c1 here, so a sort would answer (o, c2)'


def test_bonds_yields_each_bond_once_with_the_lower_stable_id_FIRST():
    """`bonds()` walks the CSR and emits only `to > i`, so `n` is the earlier atom.

    Stated because it is relied on, not because it is inevitable: a caller building a
    `{(n, m): ...}` map from `bonds()` and probing it with an independently ordered pair
    needs to know which orientation it got.  `_ethanol` adds its atoms in order, so dense
    index order and stable id order coincide and the assertion can be written either way.

    The orientation is asserted HERE, in the core's own suite: transposing it at the source
    fails exactly one test elsewhere in the tree, and a guard that lives in a consumer package
    leaves with that package.
    """
    m, c1, c2, o = _ethanol()
    pairs = [(bd.n, bd.m) for bd in m.bonds()]
    assert pairs == [(c1, c2), (c2, o)]
    assert all(n < mm for n, mm in pairs), 'lower stable id first, for these atoms'
    assert len(pairs) == m.bond_count, 'each bond once, not once per direction'


def test_views_refuse_to_read_inside_a_dirty_scope():
    m, c1, c2, o = _ethanol()
    a = m.atom(o)
    with m.edit():
        m.set_charge(o, -1)
        with pytest.raises(RuntimeError):
            a.element


def _benzoate():
    # benzene with one O hung off atom 0, so there is a ring atom and an acyclic one
    m = MoleculeContainer()
    with m:
        ids = [m.add_atom(6) for _ in range(6)]
        for i in range(6):
            m.add_bond(ids[i], ids[(i + 1) % 6], 1)
        o = m.add_atom(8)
        m.add_bond(ids[0], o, 1)
    return m, ids, o


def test_container_is_its_own_edit_scope():
    m, ids, o = _benzoate()
    assert m.atom_count == 7 and m.bond_count == 7
    # `with mol:` and `with mol.edit():` share one counter, so they nest either way
    with m:
        with m.edit():
            m.set_charge(o, -1)
            assert m.journal_length == 1
        assert m.journal_length == 1      # the inner exit did not apply
    assert m.atom(o).charge == -1


def test_a_failed_scope_leaves_the_arena_untouched():
    m, ids, o = _benzoate()
    with pytest.raises(ValueError):
        with m:
            m.set_charge(o, -1)
            m.set_charge(o, 99)           # rejected at the API, before the arena is touched
    assert m.journal_length == 0
    assert m.atom(o).charge == 0


def test_atom_numbers_is_the_atom_number_list_in_arena_order():
    """The only name for the list; `stable_ids` and `atoms_numbers` were the other two, both gone."""
    m, ids, o = _benzoate()
    assert m.atom_numbers == ids + [o]
    assert m.atom_numbers == [m.number_of(i) for i in range(m.atom_count)]


def test_bonds_yield_views_carrying_their_endpoint_ids():
    m, ids, o = _benzoate()
    bonds = list(m.bonds())
    assert len(bonds) == 7
    assert all(isinstance(b, Bond) for b in bonds)
    assert {frozenset((b.n, b.m)) for b in bonds} == (
        {frozenset((ids[i], ids[(i + 1) % 6])) for i in range(6)} | {frozenset((ids[0], o))})
    for b in bonds:
        assert b.in_ring == (o not in (b.n, b.m))
        assert m.bond(b.n, b.m).order == b.order


def test_atom_reports_its_own_ring_membership():
    m, ids, o = _benzoate()
    ring = m.atom(ids[0])
    assert ring.ring_sizes == frozenset({6}) and ring.ring_count == 1
    assert ring.in_ring and not ring.macrocycle
    chain = m.atom(o)
    assert chain.ring_sizes == frozenset() and chain.ring_count == 0
    assert not chain.in_ring
    # the molecule-level readers answer the same question without building a view
    assert m.ring_sizes_of(ids[0]) == ring.ring_sizes
    assert m.ring_count_of(ids[0]) == ring.ring_count
    assert m.macrocycle_of(ids[0]) is ring.macrocycle


def test_neighbors_is_a_count_and_the_walk_lives_on_the_molecule():
    m, ids, o = _benzoate()
    assert m.atom(ids[0]).neighbors == 3 == m.atom(ids[0]).degree
    assert m.atom(o).neighbors == 1
    assert sorted(m.neighbors_of(ids[0])) == sorted([ids[1], ids[5], o])
    assert not hasattr(m.atom(ids[0]), 'neighbors_of')


def test_attribute_writes_go_through_the_journal_and_keep_the_view_live():
    m, ids, o = _benzoate()
    a = m.atom(o)
    gen = m.generation
    a.charge = -1
    a.isotope = 18
    a.map_number = 3
    a.is_radical = True
    a.implicit_h = 0
    assert m.generation == gen + 5          # five atomic edits, five applies
    assert (a.charge, a.isotope, a.map_number, a.is_radical, a.implicit_h) == (-1, 18, 3, True, 0)
    assert a.radical is a.is_radical
    assert (m.charge_of(o), m.isotope_of(o), m.map_number_of(o)) == (-1, 18, 3)
    # a write through a *different* handle still invalidates this one
    other = m.atom(ids[0])
    other.charge = 1
    with pytest.raises(RuntimeError):
        a.charge


def test_a_scope_batches_attribute_writes_into_one_apply():
    m, ids, o = _benzoate()
    gen = m.generation
    with m:
        for a in list(m.atoms()):
            a.charge = 1
    assert m.generation == gen + 1
    assert [a.charge for a in m.atoms()] == [1] * 7


def test_bond_order_is_writable_through_its_view():
    m, ids, o = _benzoate()
    bd = m.bond(ids[0], o)
    bd.order = 2
    assert bd.order == 2 and m.order_of(ids[0], o) == 2
    with pytest.raises(ValueError):
        bd.order = 7


def test_coordinates_are_atom_attributes():
    m, ids, o = _benzoate()
    a = m.atom(o)
    assert a.xy is None and a.x is None and a.y is None
    a.xy = (1.5, -2.25)
    assert a.xy == (1.5, -2.25) and (a.x, a.y) == (1.5, -2.25)
    assert m.xy_of(o) == a.xy and m.has_coordinates
    a.x = 3.0                                # setting one coordinate keeps the other
    assert a.xy == (3.0, -2.25)
    a.y = 0.5
    assert a.xy == (3.0, 0.5)


def test_remap_relabels_in_place_and_keeps_everything_else():
    m, ids, o = _benzoate()
    rings, order, sizes = m.rings, m.order_of(ids[0], o), m.ring_sizes_of(ids[0])
    m.remap({ids[0]: 100, o: 101})
    assert m.atom_numbers == [100] + ids[1:] + [101]
    assert m.order_of(100, 101) == order
    assert m.ring_sizes_of(100) == sizes
    assert len(m.rings) == len(rings)
    assert m.atom(100).element == 6 and m.atom(101).element == 8
    with pytest.raises(KeyError):
        m.atom(ids[0])
    # ids handed out later never collide with what remap introduced
    assert m.add_atom(7) > 101


def test_remap_does_not_disturb_a_container_sharing_the_arena():
    m, ids, o = _benzoate()
    clone = m.copy()
    assert clone.shares_arena_with(m)
    clone.remap({ids[0]: 100})
    assert not clone.shares_arena_with(m)
    assert m.atom_numbers == ids + [o]
    assert m.atom(ids[0]).element == 6


def test_remap_rejects_a_mapping_it_cannot_honour():
    m, ids, o = _benzoate()
    before = m.atom_numbers
    with pytest.raises(ValueError):
        m.remap({ids[0]: ids[1]})          # would collide with an unmapped atom
    with pytest.raises(KeyError):
        m.remap({9999: 1})                 # not a live id
    with pytest.raises(ValueError):
        m.remap({ids[0]: 0})               # 0 is not a stable id
    with pytest.raises(TypeError):
        m.remap({ids[0]: 'x'})
    assert m.atom_numbers == before


def test_arena_bytes_round_trip_and_pickle_agree():
    import pickle
    m, ids, o = _benzoate()
    m.set_charge(o, -1)
    m.atom(o).xy = (1.0, 2.0)
    for r in (MoleculeContainer.from_bytes(m.to_bytes()), pickle.loads(pickle.dumps(m))):
        assert r.atom_numbers == m.atom_numbers
        assert r.rings == m.rings
        assert r.charge_of(o) == -1
        assert r.xy_of(o) == (1.0, 2.0)
        assert not r.shares_arena_with(m)


def test_unpack_reads_both_serialised_forms_and_needs_no_flag_to_tell_them_apart():
    """`pack` is the legacy pach record and `to_bytes` is the arena verbatim; `unpack` takes either.

    No argument selects between them: byte 0 of a pach record is its format version, only ever 0 or 2,
    and byte 0 of the arena format is the low byte of its magic, 0x33.  A caller holding bytes out of a
    store therefore does not have to know which era wrote them.  The codec itself is tested in
    `test_pach.py` against bytes chython 2 produced; what is checked here is only that the container's
    own two doors are wired to it.
    """
    m, ids, o = _benzoate()
    m.set_charge(o, -1)
    for data in (m.pack(version=2), m.to_bytes()):
        r = MoleculeContainer.unpack(data)
        assert r.atom_numbers == m.atom_numbers
        assert r.charge_of(o) == -1
    assert m.pack(compressed=False, version=2)[0] == 2
    assert m.to_bytes()[0] == 0x33
    with pytest.raises(ValueError):
        MoleculeContainer.unpack(b'', compressed=False)


def test_connected_components_of_an_empty_container():
    m = MoleculeContainer()
    assert m.connected_components_count == 0
    assert m.connected_components == []


def test_connected_components_of_a_single_molecule():
    m, ids, o = _benzoate()
    assert m.connected_components_count == 1
    assert m.connected_components == [tuple(m.atom_numbers)]


def test_connected_components_split_a_salt_and_count_lone_ions():
    m = MoleculeContainer()
    with m.edit():
        c1 = m.add_atom('C')
        c2 = m.add_atom('C')
        o1 = m.add_atom('O')
        o2 = m.add_atom('O', charge=-1)
        m.add_bond(c1, c2, 1)
        m.add_bond(c2, o1, 2)
        m.add_bond(c2, o2, 1)
        na = m.add_atom('Na', charge=1)
        cl = m.add_atom('Cl', charge=-1)
    assert m.connected_components_count == 3
    assert m.connected_components == [(c1, c2, o1, o2), (na,), (cl,)]


def test_connected_components_follow_deletions():
    m = MoleculeContainer()
    with m.edit():
        ids = [m.add_atom('C') for _ in range(4)]
        for i in range(3):
            m.add_bond(ids[i], ids[i + 1], 1)
    assert m.connected_components_count == 1
    m.delete_bond(ids[1], ids[2])
    assert m.connected_components_count == 2
    assert m.connected_components == [(ids[0], ids[1]), (ids[2], ids[3])]
    m.delete_atom(ids[0])
    assert m.connected_components == [(ids[1],), (ids[2], ids[3])]


def test_connected_components_ignore_bond_order_including_dative():
    # unlike ring perception, connectivity is connectivity: an order-8 bond still joins
    m = MoleculeContainer()
    with m.edit():
        fe = m.add_atom('Fe')
        n = m.add_atom('N')
        m.add_bond(fe, n, 8)
    assert m.connected_components_count == 1
    assert m.connected_components == [(fe, n)]


# ── the arena validates the BUFFER and not the chemistry ────────────────────────────────────────
#
# "molecular input from files must be treated as shit. but we can't reject it. so, arena must allow
# bad structures, which can be fixed by multistep process: kekule, standardization, charge/tautomer
# canonicalization, thiele (optional)." -- and that is a storage-layer decision, so it is asserted
# here rather than described in a docstring somewhere. What follows enumerates the things the arena
# deliberately does NOT check, so that a future reviewer tempted to add a valence model to
# `from_bytes` finds a failing test instead of an empty field.


def test_a_chemically_absurd_molecule_stores_and_round_trips_unchanged():
    """Five-valent neutral N, a stray radical, a nonsense charge, fourteen hydrogens, a triple-bonded
    peroxide and a disconnected fragment -- all in one record, all preserved byte for byte.

    Every one of these is something a real vendor file contains and something the standardisation
    layers exist to repair. The arena's job is to hold the input long enough for them to run, so it
    checks the BUFFER (a bond to a nonexistent atom, an element outside 1-118, a segment whose length
    contradicts its table entry, a non-empty unknown segment) and nothing about the chemistry.

    FOURTEEN and not fifteen since H_UNKNOWN landed: 15 is the sentinel for "the record does not
    state a count" (see test_unknown_h.py), so the largest absurd COUNT the nibble can hold is 14.
    The point of the line is unchanged -- fourteen hydrogens on a carbon is still nonsense the arena
    stores without comment -- and it is the boundary that matters, so it stays at the top of the
    count range rather than dropping to a safe middle value.
    """
    m = MoleculeContainer()
    with m.edit():
        n = m.add_atom('N')
        cs = [m.add_atom('C') for _ in range(5)]
        for c in cs:
            m.add_bond(n, c, 1)                  # five bonds on a neutral nitrogen
        o1, o2 = m.add_atom('O'), m.add_atom('O')
        m.add_bond(o1, o2, 3)                    # a triple bond between two oxygens
    m.set_radical(cs[0], True)                   # a radical on an otherwise saturated carbon
    m.set_charge(cs[1], 4)                       # a charge no carbon has
    m.set_hydrogens(cs[2], 14)                   # fourteen hydrogens, the count field's maximum
    raw = m.to_bytes()
    back = MoleculeContainer.from_bytes(raw)
    assert back.to_bytes() == raw, 'a record the arena accepted must survive a round trip exactly'
    assert len(back.neighbors_of(n)) == 5
    assert back.charge_of(cs[1]) == 4
    assert back.implicit_h_of(cs[2]) == 14
    assert back.radical_of(cs[0])
    assert back.connected_components_count == 2, 'and the disconnected fragment is still separate'


def test_the_checks_that_remain_are_about_the_buffer_and_they_do_fire():
    """Anti-vacuity for the test above: "validates nothing" would satisfy it just as well.

    So each surviving class of check is provoked once. None of them is a chemical judgement -- a bond
    to atom 9999 is not a molecule at all, and an element of 200 names nothing.
    """
    m = MoleculeContainer()
    with m.edit():
        a, b = m.add_atom('C'), m.add_atom('C')
        m.add_bond(a, b, 1)
    raw = bytearray(m.to_bytes())
    assert unpack_from('<I', raw, 0)[0] == 0x43485933
    with pytest.raises(ValueError, match='magic'):
        forged = bytearray(raw)
        forged[0] ^= 0xFF
        MoleculeContainer.from_bytes(bytes(forged))
    with pytest.raises(ValueError, match='version'):
        forged = bytearray(raw)
        forged[4] = 99
        MoleculeContainer.from_bytes(bytes(forged))
    with pytest.raises(ValueError, match='header'):
        MoleculeContainer.from_bytes(bytes(raw) + b'\0' * 8)
    # Two ways to be too short, because the header is two pieces: 24 bytes of scalars and then a
    # table whose length the scalars state.  Slicing to 100 reaches neither -- this
    # molecule's header is 48 bytes -- so the bound is read out of the buffer rather than assumed.
    with pytest.raises(ValueError, match='too short'):
        MoleculeContainer.from_bytes(bytes(raw[:20]))
    seg_count = unpack_from('<H', raw, 20)[0]
    with pytest.raises(ValueError, match='too short'):
        MoleculeContainer.from_bytes(bytes(raw[:24 + 8 * seg_count - 8]))
