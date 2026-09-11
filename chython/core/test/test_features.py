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
from chython.core._core import sig_mask


def chain(elements, orders=None):
    """Linear chain of the given elements; single bonds unless orders given."""
    m = MoleculeContainer()
    ids = [m.add_atom(e) for e in elements]
    if orders is None:
        orders = [1] * (len(ids) - 1)
    for i, o in enumerate(orders):
        m.add_bond(ids[i], ids[i + 1], o)
    return m, ids


def ring(elements, orders=None):
    m = MoleculeContainer()
    ids = [m.add_atom(e) for e in elements]
    n = len(ids)
    if orders is None:
        orders = [1] * n
    for i, o in enumerate(orders):
        m.add_bond(ids[i], ids[(i + 1) % n], o)
    return m, ids


def test_union_feature_words_is_four_words():
    m, ids = chain([6, 6, 6])
    sig = m._union_feature_words
    assert len(sig) == 4
    assert all(isinstance(w, int) and 0 <= w < 1 << 64 for w in sig)


def test_union_feature_words_is_the_or_of_atom_words():
    m, ids = chain([6, 6, 8])
    expected = [0, 0, 0, 0]
    for s in ids:
        for k, w in enumerate(m.features_of(s)):
            expected[k] |= w
    assert m._union_feature_words == tuple(expected)


def test_light_element_bit_is_fiftyseven_minus_atomic_number():
    m, ids = chain([6, 7, 8])
    assert m.features_of(ids[0])[0] >> (57 - 6) & 1
    assert m.features_of(ids[1])[0] >> (57 - 7) & 1
    assert m.features_of(ids[2])[0] >> (57 - 8) & 1
    # and the transfer bit stays clear for light elements
    assert not m.features_of(ids[0])[0] & 1


def test_heavy_element_sets_the_transfer_bit_and_word_two():
    m, ids = chain([92, 8])                      # uranyl-ish U-O fragment
    w1, w2, _, _ = m.features_of(ids[0])
    assert w1 & 1                                # transfer
    assert w1 >> 1 & 0xffffffffffffff == 0       # no light element bit
    assert w2 >> (92 - 57) & 1


def test_barium_is_the_last_light_element():
    m, ids = chain([56, 17])
    w1, w2, _, _ = m.features_of(ids[0])
    assert not w1 & 1
    assert w1 >> 1 & 1                           # 57 - 56 == 1
    m2, ids2 = chain([57, 17])                   # lanthanum crosses over
    assert m2.features_of(ids2[0])[0] & 1
    assert m2.features_of(ids2[0])[1] & 1        # 57 - 57 == 0


def test_radical_occupies_word_two_top_bits():
    m = MoleculeContainer()
    a = m.add_atom(6, radical=True)
    b = m.add_atom(6)
    m.add_bond(a, b, 1)
    assert m.features_of(a)[1] >> 63 & 1
    assert not m.features_of(a)[1] >> 62 & 1
    assert m.features_of(b)[1] >> 62 & 1
    assert not m.features_of(b)[1] >> 63 & 1


def test_bond_bits_reflect_orders_and_ring_membership():
    m, ids = ring([6] * 6, [1, 2, 1, 2, 1, 2])
    w1 = m.features_of(ids[0])[0]
    # this ring is written Kekule, so every ring bond is a plain ring bond (bit 62); bit 57 belongs
    # to a stored order 4 and nothing derives it from a bond pattern
    assert not w1 >> 57 & 1                      # no aromatic bond: none was written
    assert not w1 >> 58 & 1                      # no acyclic bond
    assert w1 >> 59 & 1                          # single present
    assert w1 >> 60 & 1                          # double present
    assert not w1 >> 61 & 1                      # no triple

    assert w1 >> 62 & 1                          # plain ring bond
    assert not w1 >> 63 & 1                      # and no dative

    m2, ids2 = chain([6, 6, 6], [3, 1])
    w = m2.features_of(ids2[1])[0]
    assert w >> 58 & 1                           # acyclic
    assert not w >> 57 & 1                       # no ring+aromatic bond
    assert w >> 61 & 1                           # triple
    assert w >> 59 & 1                           # single


def test_dative_bonds_occupy_the_top_order_bit():
    m, ids = chain([7, 5], [8])                  # amine-borane dative bond
    assert m.features_of(ids[0])[0] >> 63 & 1
    assert not m.features_of(ids[0])[0] >> 59 & 1


def test_counts_land_in_the_declared_word_three_fields():
    # CH3-O-CH3 with the hydrogen counts stated, since nothing derives them
    m = MoleculeContainer()
    c1 = m.add_atom(6, implicit_h=3)
    o = m.add_atom(8, implicit_h=0)
    c2 = m.add_atom(6, implicit_h=3)
    m.add_bond(c1, o, 1)
    m.add_bond(o, c2, 1)
    w = m.features_of(o)[2]
    assert w >> (0 + 0) & 1                      # ether O: both neighbours are C
    assert w >> (9 + 2) & 1                      # degree 2
    assert w >> (17 + 0) & 1                     # 0 implicit H
    assert w >> (22 + 0) & 1                     # 0 explicit H
    assert w >> (27 + 0) & 1                     # 0 total H
    c = m.features_of(c1)[2]
    assert c >> (0 + 1) & 1                      # methyl carbon: x == 1 (the O)
    assert c >> (9 + 1) & 1                      # degree 1
    assert c >> (17 + 3) & 1                     # 3 implicit H
    assert c >> (27 + 3) & 1                     # 3 total H


def test_word_three_counts_skip_dative_bonds_and_the_stored_scalars_do_not():
    """RULING, and the one place in the tree where two counts of the same thing are both correct.

    Word 3 feeds the `D` and `x` query primitives, and a coordination contact is NOT a substituent:
    trimethylamine donating its lone pair to an iron has three substituents in every sense a rule
    cares about, which is why `derive_scalars` already calls it sp3 and counts no dative bond
    towards hybridization.  So `D` and `x` agree with `z` -- `D3`, not `D4`, and `x0`, not `x1`.

    `degree_of` and `heteroatoms_of` KEEP counting it, because they are structural: degree is the
    CSR row length, `_pach.pxi` writes and reads it as that row length, and `_stereo.pxi` reads it
    as connectivity.  Both numbers are asserted here so that "these disagree" is a documented
    property rather than something a later reader takes for a bug and unifies.
    """
    m = MoleculeContainer()
    n = m.add_atom(7, implicit_h=0)
    iron = m.add_atom(26)
    for _ in range(3):
        m.add_bond(n, m.add_atom(6, implicit_h=3), 1)
    m.add_bond(n, iron, 8)

    w = m.features_of(n)[2]
    assert w >> (9 + 3) & 1 and not w >> (9 + 4) & 1     # the query sees D3
    assert w >> (0 + 0) & 1 and not w >> (0 + 1) & 1     # and x0
    assert m.degree_of(n) == 4                           # the structure still has four bonds
    assert m.heteroatoms_of(n) == 1                      # and the iron is still a heteroatom
    assert m.hybridization_of(n) == 1                    # `z` was already excluding it

    # and from the acceptor's side, where the dative bond is the atom's ONLY bond: D0, not D1
    w = m.features_of(iron)[2]
    assert w >> (9 + 0) & 1 and not w >> (9 + 1) & 1
    assert w >> (0 + 0) & 1 and not w >> (0 + 1) & 1
    assert m.degree_of(iron) == 1


def test_count_fields_saturate_rather_than_overflow():
    m = MoleculeContainer()
    centre = m.add_atom(16)                     # sulfur hub, degree 8
    for _ in range(8):
        leaf = m.add_atom(9)
        m.add_bond(centre, leaf, 1)
    w = m.features_of(centre)[2]
    assert w >> (0 + 8) & 1                      # heteroatoms cap at 8+
    assert w >> (9 + 7) & 1                      # degree cap at 7+

    # saturation: 9 heteroatom neighbours land on the same cap bit as 8
    m2 = MoleculeContainer()
    c2 = m2.add_atom(16)
    for _ in range(9):
        m2.add_bond(c2, m2.add_atom(9), 1)
    w2 = m2.features_of(c2)[2]
    assert w2 >> (0 + 8) & 1                     # still at the 8+ cap


def test_charge_field_spans_minus_four_to_plus_eight():
    for charge, offset in ((-4, 0), (-1, 3), (0, 4), (3, 7), (8, 12)):
        m = MoleculeContainer()
        a = m.add_atom(7)
        b = m.add_atom(6)
        m.add_bond(a, b, 1)
        m.set_charge(a, charge)
        assert m.features_of(a)[2] >> (33 + offset) & 1


def test_charge_outside_the_field_cannot_be_built():
    m = MoleculeContainer()
    a = m.add_atom(6)
    with pytest.raises(ValueError):
        m.set_charge(a, -8)


def test_isotope_is_a_delta_from_the_mdl_isotope():
    m = MoleculeContainer()
    a = m.add_atom(6, isotope=13)               # MDL isotope of C is 12
    b = m.add_atom(6)
    m.add_bond(a, b, 1)
    assert m.features_of(a)[2] >> (46 + 8 + 1) & 1
    assert not m.features_of(a)[2] >> 63 & 1     # isotope is set
    assert m.features_of(b)[2] >> 63 & 1         # isotope unset
    assert not m.features_of(b)[2] >> (46 + 8) & 1


def test_isotope_delta_saturates_at_both_ends():
    m = MoleculeContainer()
    a = m.add_atom(6, isotope=24)               # delta +12, saturates to +8
    b = m.add_atom(6, isotope=1)                # delta -11, saturates to -8
    m.add_bond(a, b, 1)
    assert m.features_of(a)[2] >> (46 + 8 + 8) & 1
    assert m.features_of(b)[2] >> (46 + 8 - 8) & 1
    # a real, unsaturated delta for contrast
    m2 = MoleculeContainer()
    t = m2.add_atom(1, isotope=3)               # tritium: delta +2
    c = m2.add_atom(6)
    m2.add_bond(t, c, 1)
    assert m2.features_of(t)[2] >> (46 + 8 + 2) & 1


def test_hybridization_and_stereo_bits():
    m, ids = chain([6, 6, 6], [2, 1])
    assert m.features_of(ids[0])[3] >> (2 - 1) & 1     # sp2
    assert m.features_of(ids[2])[3] >> (1 - 1) & 1     # sp3
    assert not m.features_of(ids[0])[3] >> 6 & 1       # no stereo label


def test_ring_size_field_is_the_atom_ring_word_shifted():
    m, ids = ring([6] * 7)
    for s in ids:
        assert m.features_of(s)[3] >> 22 & 0x1ffffff == m.ring_sizes_word_of(s)
        assert m.features_of(s)[3] >> (22 + 7) & 1      # exact size seven
        assert m.features_of(s)[3] >> (47 + 1) & 1      # one ring


def test_acyclic_atoms_report_zero_rings():
    m, ids = chain([6, 6, 6])
    for s in ids:
        assert m.features_of(s)[3] >> 22 & 0x1ffffff == 0
        assert m.features_of(s)[3] >> 47 & 1            # ring count zero
        assert m.features_of(s)[3] >> 56 & 1            # aromatic count zero


def test_reserved_bits_are_zero():
    # The free span moved when Task 11 took word IV bits 7-8 for the parity-configured span: bits
    # 9-21 are what is left between it and ring_sizes at bit 22.  Bit 8 IS SET on every atom here --
    # none of them carries a configured parity -- which is the span's "not configured" half.
    m, ids = ring([6] * 6, [1, 2, 1, 2, 1, 2])
    for s in ids:
        assert m.features_of(s)[3] >> 9 & 0x1fff == 0
        assert m.features_of(s)[3] >> 7 & 3 == 2, 'no parity configured, so the span says so'


def test_features_are_recomputed_after_an_edit():
    m, ids = ring([6] * 6)
    assert m.features_of(ids[0])[3] >> 22 & 0x1ffffff == 1 << 6
    before = m._union_feature_words
    m.delete_bond(ids[0], ids[1])
    # SEG_FEATURES is rebuilt on every fold, so a chain must not still look cyclic
    assert m.features_of(ids[0])[3] >> 22 & 0x1ffffff == 0
    assert m.features_of(ids[0])[3] >> 47 & 1
    assert m._union_feature_words != before


def test_screen_admits_a_true_subgraph():
    # a cyclohexane ring is a subgraph of decalin; the screen must not reject it
    sub, _ = ring([6] * 6)
    sup = MoleculeContainer()
    ids = [sup.add_atom(6) for _ in range(10)]
    for i in range(5):
        sup.add_bond(ids[i], ids[(i + 1) % 6], 1)
    sup.add_bond(ids[5], ids[0], 1)
    for a, b in ((4, 6), (6, 7), (7, 8), (8, 9), (9, 5)):
        sup.add_bond(ids[a], ids[b], 1)
    assert sup.may_contain(sub)


def test_screen_rejects_a_missing_element():
    sub, _ = chain([6, 17])
    sup, _ = chain([6, 6, 6])
    assert not sup.may_contain(sub)


def test_screen_is_reflexive():
    m, _ = ring([6, 6, 6, 7, 6, 6])
    assert m.may_contain(m)


def test_screen_admits_every_relation_chython_two_calls_a_substructure():
    # Each pair is a substructure relation under chython 2's `<`. A screen that rejects
    # any of them makes its callers skip a search that would have matched, which is a
    # silently wrong answer rather than a slow one. These four pairs are exactly the
    # fields that do not survive embedding: degree, ring membership of a bond,
    # hybridization, and heteroatom count.
    propane, _ = chain([6, 6, 6])
    isobutane = MoleculeContainer()
    with isobutane.edit():
        centre = isobutane.add_atom(6)
        for _ in range(3):
            isobutane.add_bond(centre, isobutane.add_atom(6), 1)
    assert isobutane.may_contain(propane)

    ethane, _ = chain([6, 6])
    cyclohexane, _ = ring([6] * 6)
    assert cyclohexane.may_contain(ethane)

    butadiene = MoleculeContainer()
    with butadiene.edit():
        bd = [butadiene.add_atom(6) for _ in range(4)]
        butadiene.add_bond(bd[0], bd[1], 2)
        butadiene.add_bond(bd[1], bd[2], 1)
        butadiene.add_bond(bd[2], bd[3], 2)
    assert butadiene.may_contain(ethane)

    methanol, _ = chain([6, 8])
    methanediol, _ = chain([8, 6, 8])
    assert methanediol.may_contain(methanol)


def test_screen_admits_a_ring_the_relevant_cycles_do_not_report():
    # Cyclohexane IS a subgraph of norbornane, but norbornane's relevant cycles are its
    # two 5-rings, so no norbornane atom ever reports ring size 6. Ring descriptors are
    # therefore unsound to screen on even as pure graph theory, not merely because the
    # local environment grows.
    cyclohexane, _ = ring([6] * 6)
    norbornane = MoleculeContainer()
    with norbornane.edit():
        a = [norbornane.add_atom(6) for _ in range(7)]
        for i in range(5):
            norbornane.add_bond(a[i], a[i + 1], 1)
        norbornane.add_bond(a[5], a[0], 1)
        norbornane.add_bond(a[0], a[6], 1)      # the one-carbon bridge
        norbornane.add_bond(a[6], a[3], 1)
    assert 6 not in {s for i in norbornane.atom_numbers for s in norbornane.ring_sizes_of(i)}
    assert norbornane.may_contain(cyclohexane)


def test_screen_still_rejects_on_the_fields_it_keeps():
    # The mask must not collapse the screen to "always True". Every field left in it is
    # one that substructure matching compares exactly, so each of these is a real reject.
    ethane, _ = chain([6, 6])
    cyclohexane, _ = ring([6] * 6)

    chloromethane, _ = chain([6, 17])
    assert not ethane.may_contain(chloromethane)                      # element

    ethene = MoleculeContainer()
    with ethene.edit():
        e = [ethene.add_atom(6) for _ in range(2)]
        ethene.add_bond(e[0], e[1], 2)
    assert not ethane.may_contain(ethene)                             # bond order

    cation = MoleculeContainer()
    with cation.edit():
        c = [cation.add_atom(6, charge=1), cation.add_atom(6)]
        cation.add_bond(c[0], c[1], 1)
    assert not ethane.may_contain(cation)                             # charge

    heavy = MoleculeContainer()
    with heavy.edit():
        h = [heavy.add_atom(6, isotope=13), heavy.add_atom(6)]
        heavy.add_bond(h[0], h[1], 1)
    assert not ethane.may_contain(heavy)                              # isotope

    radical = MoleculeContainer()
    with radical.edit():
        r = [radical.add_atom(6, radical=True), radical.add_atom(6)]
        radical.add_bond(r[0], r[1], 1)
    assert not ethane.may_contain(radical)                            # radical

    # the topology triple is outside the mask: the screen is permissive about ring bonds
    assert ethane.may_contain(cyclohexane)                            # screen can't reject on ring bits
    assert cyclohexane.may_contain(ethane)


def test_screen_rejects_none_rather_than_crashing():
    m, _ = chain([6, 6])
    with pytest.raises(TypeError):
        m.may_contain(None)


def test_signature_mask_drops_exactly_the_environment_fields():
    w1, w2, w3, w4 = sig_mask()
    assert w1 == 0xB9FFFFFFFFFFFFFF           # bond topology triple (bits 57, 58, 62) all dropped
    assert not w1 >> 57 & 1                   # ring-arom bit dropped: aromaticity changes under embedding
    assert not w1 >> 58 & 1                   # not-in-ring bit dropped: ring membership changes
    assert not w1 >> 62 & 1                   # ring-plain bit dropped: ring membership changes
    assert w1 >> 63 & 1                       # bond order 8 (dative) stays in the mask
    assert w2 == (1 << 64) - 1                # elements and radical are exact
    assert w3 == ((1 << 64) - 1) ^ ((1 << 33) - 1)
    assert w3 >> 33 & 1 and not w3 & 1        # charge kept, heteroatom count dropped
    assert w4 == 0                            # nothing in word IV survives embedding


# The isotope and MDL reference tables are generated C arrays, so their tests live in
# test_element_tables.py, next to the generator that produces them -- including the two that use
# chython 2 as the oracle.


def test_element_buckets_group_atoms_by_element():
    m, ids = chain([6, 8, 6, 7, 6])
    assert m.atoms_of_element(6) == (ids[0], ids[2], ids[4])
    assert m.atoms_of_element(8) == (ids[1],)
    assert m.atoms_of_element(7) == (ids[3],)


def test_absent_elements_have_empty_buckets():
    m, ids = chain([6, 6])
    assert m.atoms_of_element(8) == ()
    assert m.atoms_of_element(1) == ()
    assert m.atoms_of_element(118) == ()


def test_element_counts_sum_to_the_atom_count():
    m, ids = chain([6, 8, 6, 7, 6])
    assert m.element_counts == {6: 3, 7: 1, 8: 1}
    assert sum(m.element_counts.values()) == m.atom_count


def test_every_atom_appears_in_exactly_one_bucket():
    m, ids = chain([6, 7, 8, 9, 16, 17, 35, 53, 5, 15])
    seen = []
    for number in range(1, 119):
        seen.extend(m.atoms_of_element(number))
    assert sorted(seen) == sorted(m.atom_numbers)


def test_bucket_for_the_last_element_needs_no_special_case():
    m, ids = chain([118, 6])
    assert m.atoms_of_element(118) == (ids[0],)
    assert m.atoms_of_element(6) == (ids[1],)


def test_buckets_are_in_ascending_index_order():
    m, ids = chain([6] * 12)
    assert m.atoms_of_element(6) == tuple(ids)


def test_element_buckets_reject_a_number_out_of_range():
    m, ids = chain([6, 6])
    with pytest.raises(ValueError):
        m.atoms_of_element(0)
    with pytest.raises(ValueError):
        m.atoms_of_element(119)


def test_element_buckets_are_rebuilt_after_an_edit():
    m, ids = chain([6, 6, 6])
    with m.edit():
        m.add_bond(ids[0], m.add_atom(8), 1)
    assert len(m.atoms_of_element(8)) == 1
    assert m.atoms_of_element(6) == (ids[0], ids[1], ids[2])


def test_element_buckets_survive_a_pack_round_trip():
    m, ids = chain([6, 8, 6, 7, 6])
    back = MoleculeContainer.from_bytes(m.to_bytes())
    assert back.element_counts == m.element_counts
    assert back.atoms_of_element(6) == m.atoms_of_element(6)


def test_empty_molecule_has_no_buckets():
    m = MoleculeContainer()
    assert m.element_counts == {}
    assert m.atoms_of_element(6) == ()


def test_bond_topology_bits_are_mutually_exclusive_per_bond():
    # propane: two acyclic single bonds. Every atom sees only the not-in-ring bit.
    m, ids = chain([6, 6, 6])
    for s in ids:
        w0 = m.features_of(s)[0]
        assert w0 & 1 << 58, 'acyclic bond must set bit 58'
        assert not w0 & 1 << 57, 'acyclic bond must not set the aromatic-ring bit'
        assert not w0 & 1 << 62, 'acyclic bond must not set the plain-ring bit'


def test_plain_ring_bond_sets_bit_62_only():
    # cyclohexane: six single ring bonds, no aromaticity
    m, ids = ring([6] * 6)
    for s in ids:
        w0 = m.features_of(s)[0]
        assert w0 & 1 << 62
        assert not w0 & 1 << 57
        assert not w0 & 1 << 58


def test_aromatic_ring_bond_sets_bit_57_only():
    # benzene as the file drew it: six order-4 bonds. This test was an xfail until order 4 became
    # a stored order -- it asked for bit 57 from a KEKULE ring, on the theory that a perception
    # pass would mark the bonds aromatic. It gets the bit from the ORDER instead, which is the
    # difference between deriving a chemical judgement and recording what the input said.
    m, ids = ring([6] * 6, [4] * 6)
    assert m.aromatic_bond_count == 6
    for s in ids:
        w0 = m.features_of(s)[0]
        assert w0 & 1 << 57, 'aromatic ring bond must set bit 57'
        assert not w0 & 1 << 62, 'an aromatic bond is not a plain ring bond'
        assert not w0 & 1 << 58


def test_the_same_ring_written_kekule_answers_differently_and_that_is_the_point():
    # The could-have-failed half of the test above: bit 57 is not something every benzene has.
    # A Kekule benzene is a DIFFERENT molecule to the feature words -- plain ring bonds, orders 1
    # and 2 -- and a query for an aromatic bond does not match it. Nothing in the core silently
    # bridges the two spellings; `kekule()` and `thiele()` are the only crossing.
    arom, _ = ring([6] * 6, [4] * 6)
    kek, ids = ring([6] * 6, [2, 1, 2, 1, 2, 1])
    assert kek.aromatic_bond_count == 0
    for s in ids:
        w0 = kek.features_of(s)[0]
        assert not w0 & 1 << 57
        assert w0 & 1 << 62, 'a Kekule ring bond is a plain ring bond'
    assert arom._union_feature_words != kek._union_feature_words


def test_an_atom_may_carry_two_topology_bits_from_two_bonds():
    # toluene, ring written aromatic: the methyl carbon is acyclic, the ring carbon it hangs off
    # sees both. The union word ORs every incident bond, so one atom carries two of a span that is
    # one-hot PER BOND -- which is why a folded query box must be given a half-edge word and not
    # this one (see `atom_admits`).
    m = MoleculeContainer()
    ids = [m.add_atom(6) for _ in range(7)]
    for i in range(6):
        m.add_bond(ids[i], ids[(i + 1) % 6], 4)
    m.add_bond(ids[0], ids[6], 1)
    w0 = m.features_of(ids[0])[0]
    assert w0 & 1 << 57, 'two aromatic ring bonds'
    assert w0 & 1 << 58, 'one acyclic bond'
    assert not w0 & 1 << 62
    assert m.features_of(ids[6])[0] & 1 << 58


def test_sig_mask_excludes_every_topology_bit():
    # ring membership and aromaticity do not survive substructure embedding, so none of the
    # three topology bits may take part in the screen.
    #
    # Bit 57 stays out even though it now fires exactly on order 4, and the reason is the one this
    # test was written for rather than an oversight: it is the ONLY separator between order 4 and
    # order 8 in this word (the two share bit 63 -- word 0 is 64/64 spent), so admitting it to the
    # screen would let a query demanding a coordination bond screen out a molecule whose only
    # order-8-bit-carrying bonds are aromatic, and vice versa. Screening is a cheap prefilter that
    # must never reject a true match; the exact separation happens in the box, where both bits are
    # available together.
    assert sig_mask()[0] == 0xB9FFFFFFFFFFFFFF
    for bit in (57, 58, 62):
        assert not sig_mask()[0] & 1 << bit
    assert sig_mask()[0] & 1 << 63, 'bond order 8 is compared, so it stays in the mask'


def test_a_nonaromatic_ring_double_bond_is_not_aromatic():
    # cyclohexene: the ring C=C is sp2 at both ends but not aromatic. Deriving aromaticity
    # from hybridization would light bit 57 here and make [C;a] match cyclohexene.
    m, ids = ring([6] * 6, orders=[2, 1, 1, 1, 1, 1])
    for s in ids:
        assert not m.features_of(s)[0] & 1 << 57


# ---------------------------------------------------------------------------
# Task 3: SEG_EDGE_WORD — one u64 per half-edge, target element + bond bits
# ---------------------------------------------------------------------------


def test_edge_word_carries_the_target_element():
    # C-O: the halfedge out of the carbon must describe oxygen (bit 57 - 8 == 49)
    m, ids = chain([6, 8])
    (word,) = m.edge_words_of(ids[0])
    assert word & 1 << 49, 'the halfedge out of C must carry O in the element span'
    assert not word & 1 << 51, 'it must not carry the source element'
    (word,) = m.edge_words_of(ids[1])
    assert word & 1 << 51, 'the halfedge out of O must carry C'


def test_edge_word_carries_this_bonds_order_only():
    m, ids = chain([6, 6, 6], [1, 2])
    orders = sorted(w & 0xB800000000000000 for w in m.edge_words_of(ids[1]))
    assert orders == sorted([1 << 59, 1 << 60]), 'one order bit per halfedge, not the union'


def test_edge_word_topology_matches_the_bond():
    m, ids = ring([6] * 6)
    for word in m.edge_words_of(ids[0]):
        assert word & 1 << 62, 'plain ring bond'
        assert not word & 1 << 58
    m, ids = chain([6, 6])
    for word in m.edge_words_of(ids[0]):
        assert word & 1 << 58, 'acyclic bond'


def test_edge_word_marks_heavy_targets_with_bit_zero():
    # uranium is element 92: word 0 carries only the "heavy" marker, bit 0
    m, ids = chain([6, 92])
    (word,) = m.edge_words_of(ids[0])
    assert word & 1, 'a target above element 56 sets bit 0'
    assert not word & 0x01FFFFFFFFFFFFFE, 'and no other element bit'


def test_single_atom_has_no_edge_words():
    # a single atom has no bonds, so edge_words_of must return an empty list
    m = MoleculeContainer()
    a = m.add_atom(6)
    assert m.edge_words_of(a) == []


def test_edge_words_survive_a_pack_round_trip():
    m, ids = chain([6, 8, 7])
    back = MoleculeContainer.from_bytes(m.to_bytes())
    assert back.edge_words_of(back.atom_numbers[0]) == m.edge_words_of(ids[0])


def test_edge_word_barium_is_at_the_light_heavy_boundary():
    # Barium (56) is the last light element: encoded at bit 57-56=1.  Lanthanum (57) is
    # the first heavy element: it sets the transfer bit (bit 0) not any element span bit.
    # This is the narrowest boundary in the element encoding -- one bit apart.
    m_ba, ids_ba = chain([6, 56])
    (word,) = m_ba.edge_words_of(ids_ba[0])
    assert word >> 1 & 1, 'halfedge to Ba (element 56) must set bit 1 (57 - 56)'
    assert not word & 1, 'halfedge to Ba must not set the heavy-element transfer bit'
    m_la, ids_la = chain([6, 57])
    (word,) = m_la.edge_words_of(ids_la[0])
    assert word & 1, 'halfedge to La (element 57) must set the heavy-element transfer bit'
    assert not word >> 1 & 1, 'halfedge to La must not set bit 1 (the Ba bit)'


# ---------------------------------------------------------------------------
# Task 4: SEG_COMPONENT_LABEL — lazy connected-component labels
# ---------------------------------------------------------------------------


def test_component_labels_are_dense_and_zero_based():
    m = MoleculeContainer()
    a, b = m.add_atom(6), m.add_atom(6)
    m.add_bond(a, b, 1)
    c = m.add_atom(11)          # a lone sodium: its own component
    labels = m.component_labels()
    assert labels[a] == labels[b]
    assert labels[c] != labels[a]
    assert sorted(set(labels.values())) == [0, 1]


def test_component_labels_are_idempotent():
    m, ids = chain([6, 6, 6])
    assert m.component_labels() == m.component_labels()
    assert set(m.component_labels().values()) == {0}


def test_component_labels_follow_an_edit():
    m, ids = chain([6, 6, 6])
    assert set(m.component_labels().values()) == {0}
    m.delete_bond(ids[0], ids[1])
    assert len(set(m.component_labels().values())) == 2
