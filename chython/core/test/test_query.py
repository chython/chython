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
from chython.core import QueryContainer
from chython.core._core import (_compile_probe, _prim_probe, _query_alloc_probe,
                                _query_header_size, _query_record_sizes, _query_segment_ids,
                                _seal_probe)


def _prim_probe_with_element(name, value, negated, element):
    """A primitive applied to a box that already settled on `element`."""
    from chython.core._core import _prim_probe_seeded
    return _prim_probe_seeded(name, value, negated, element)


def test_record_sizes_are_the_documented_ones():
    sizes = _query_record_sizes()
    assert sizes['qatom_t'] == 24, 'a query atom is the same 24 bytes as an arena atom'
    assert sizes['qbox_t'] == 40
    assert sizes['qany_t'] == 16
    assert sizes['qbond_t'] == 8
    assert sizes['qclosure_t'] == 8
    assert sizes['qcomp_t'] == 12
    # position + readiness + four refs + sign + n_refs + one spare byte, packed
    assert sizes['qstereo_t'] == 28


def test_header_has_one_entry_per_segment():
    # QSEG_STEREO makes the count 11, so the header is one segment_t longer than the ten segments.
    # stereo_count is carved out of `reserved`, which is why the scalar half is still 12 uint32.
    assert _query_segment_ids()['QSEG_COUNT'] == 11
    # 12 uint32 (10 scalars + 2 reserved) + 4 uint64 signature + 11 segments * 8 bytes
    assert _query_header_size() == 12 * 4 + 4 * 8 + 11 * 8 == 168


def test_alloc_lays_segments_out_in_order_without_overlap():
    info = _query_alloc_probe(3, 2, 5, 2, 1, 1, 1)
    assert info['magic'] == 0x43485951
    assert info['version'] == 1
    assert info['atom_count'] == 3
    offsets = info['segments']
    assert offsets[0][0] == _query_header_size(), 'atoms start right after the header'
    prev_end = _query_header_size()
    for offset, length in offsets:
        assert offset >= prev_end, 'segments must not overlap'
        assert offset % 8 == 0, 'every segment is 8-byte aligned'
        prev_end = offset + length
    assert prev_end <= info['total_len']


def test_element_demand_is_always_present_and_zeroed():
    info = _query_alloc_probe(1, 0, 1, 0, 0, 0, 1)
    assert info['element_demand'] == [0] * 120


def test_empty_segments_read_from_the_zero_page():
    # no bonds, no closures: reading those segments must be safe and yield zeros
    info = _query_alloc_probe(1, 0, 1, 0, 0, 0, 1)
    assert info['closure_reads'] == [0, 0]


def test_element_forbids_every_other_element_bit():
    from chython.core._core import _prim_probe
    r = _prim_probe('element', 6, False)
    assert not r['neg'][0] & 1 << 51, 'carbon is bit 57 - 6'
    assert r['neg'][0] & 0x01FFFFFFFFFFFFFF == 0x01FFFFFFFFFFFFFF ^ (1 << 51)
    assert r['neg'][1] & 0x3FFFFFFFFFFFFFFF == 0x3FFFFFFFFFFFFFFF, 'and every heavy element'
    assert r['touched'][0] & 0x01FFFFFFFFFFFFFF


def test_heavy_element_keeps_bit_zero_and_one_word_one_bit():
    from chython.core._core import _prim_probe
    r = _prim_probe('element', 92, False)          # uranium
    assert not r['neg'][0] & 1, 'the heavy marker must stay allowed'
    assert r['neg'][0] & 0x01FFFFFFFFFFFFFE == 0x01FFFFFFFFFFFFFE
    assert not r['neg'][1] & 1 << 35, '92 - 57 == 35'


def test_negated_element_forbids_only_that_element():
    from chython.core._core import _prim_probe
    r = _prim_probe('element', 6, True)
    assert r['neg'][0] & 0x01FFFFFFFFFFFFFF == 1 << 51


def test_any_touches_nothing():
    from chython.core._core import _prim_probe
    r = _prim_probe('any', 0, False)
    assert r['neg'] == (0, 0, 0, 0)
    assert r['touched'] == (0, 0, 0, 0)


def test_metal_element_set_has_expected_popcount():
    from chython.core._core import _prim_probe
    # 93 metals total: 34 light (e<=56) + 1 heavy-marker bit in word 0 = 35 allowed in w0
    # 59 heavy metals (e>56, excluding At/Rn/Og) in word 1
    # Verify by checking neg mask popcount (forbidden = non-metal bits)
    r = _prim_probe('metal', 0, False)
    W0_ELEMENT_SPAN = 0x01FFFFFFFFFFFFFF
    W1_ELEMENT_SPAN = 0x3FFFFFFFFFFFFFFF
    allowed_w0 = W0_ELEMENT_SPAN & ~r['neg'][0]
    allowed_w1 = W1_ELEMENT_SPAN & ~r['neg'][1]
    assert bin(allowed_w0).count('1') == 35, '34 light metals + 1 heavy-marker bit'
    assert bin(allowed_w1).count('1') == 59, '59 heavy metals (excl. At, Rn, Og)'
    # spot checks: Na (11) is allowed in word 0 (bit 57-11=46)
    assert not r['neg'][0] & 1 << 46, 'sodium must be allowed'
    # carbon (6) is forbidden in word 0 (bit 57-6=51)
    assert r['neg'][0] & 1 << 51, 'carbon must be forbidden'
    # lead (82) is allowed in word 1 (bit 82-57=25)
    assert not r['neg'][1] & 1 << 25, 'lead must be allowed in word 1'


def test_charge_is_biased_by_four_from_bit_33():
    from chython.core._core import _prim_probe
    assert not _prim_probe('charge', 0, False)['neg'][2] & 1 << 37
    assert not _prim_probe('charge', 1, False)['neg'][2] & 1 << 38
    assert not _prim_probe('charge', -1, False)['neg'][2] & 1 << 36
    r = _prim_probe('charge', 0, False)
    assert r['neg'][2] & 0x3FFE00000000 == 0x3FFE00000000 ^ (1 << 37)


def test_charge_saturates_at_the_span_edges():
    from chython.core._core import _prim_probe
    assert _prim_probe('charge', 8, False)['neg'][2] == _prim_probe('charge', 12, False)['neg'][2]
    assert _prim_probe('charge', -4, False)['neg'][2] == _prim_probe('charge', -9, False)['neg'][2]


def test_radical_is_a_one_hot_pair():
    from chython.core._core import _prim_probe
    assert _prim_probe('radical', 0, False)['neg'][1] & 0xC000000000000000 == 1 << 62
    assert _prim_probe('radical', 0, True)['neg'][1] & 0xC000000000000000 == 1 << 63


def test_degree_saturates_at_seven():
    from chython.core._core import _prim_probe
    assert _prim_probe('degree', 7, False)['neg'][2] == _prim_probe('degree', 9, False)['neg'][2]
    assert not _prim_probe('degree', 2, False)['neg'][2] & 1 << 11
    # Saturation puts the demanded bit at position 16; all other span bits must be forbidden
    assert _prim_probe('degree', 7, False)['neg'][2] & 0x1FE00 == 0x1FE00 ^ (1 << 16)
    assert _prim_probe('degree', 7, False)['touched'][2] == 0x1FE00


def test_degree_span_covers_bit_sixteen_not_bit_eight():
    from chython.core._core import _prim_probe
    r = _prim_probe('degree', 6, False)
    assert r['neg'][2] & (1 << 16), 'degree 7 must be excluded by D6'
    assert not r['neg'][2] & (1 << 8), 'bit 8 belongs to the heteroatom span and must be untouched'
    assert r['touched'][2] == 0x1FE00


def test_heteroatoms_span_excludes_only_other_bits():
    from chython.core._core import _prim_probe
    r = _prim_probe('heteroatoms', 2, False)
    assert not r['neg'][2] & (1 << 2), 'x2 demanded bit must be allowed'
    assert r['neg'][2] & 0x1FF == 0x1FF ^ (1 << 2)


@pytest.mark.parametrize('name,value,word,span', (
    ('heteroatoms', 2, 2, 0x1FF),
    ('degree', 2, 2, 0x1FE00),
    ('implicit_h', 2, 2, 0x3E0000),
    ('total_h', 2, 2, 0x1F8000000),
    ('charge', 0, 2, 0x3FFE00000000),
    ('hybridization', 1, 3, 0x3F),
    ('ring_count', 1, 3, 0xFF800000000000),
))
def test_one_hot_primitives_touch_exactly_their_span(name, value, word, span):
    from chython.core._core import _prim_probe
    r = _prim_probe(name, value, False)
    assert r['touched'][word] == span, f'{name}: touched {r["touched"][word]:#x} != span {span:#x}'
    assert r['neg'][word] & ~span == 0, f'{name}: neg {r["neg"][word]:#x} has bits outside span {span:#x}'


def test_hydrogen_counts_use_separate_spans():
    from chython.core._core import _prim_probe
    implicit = _prim_probe('implicit_h', 1, False)
    total = _prim_probe('total_h', 1, False)
    assert implicit['touched'][2] == 0x3E0000
    assert total['touched'][2] == 0x1F8000000
    assert not implicit['neg'][2] & 1 << 18
    assert not total['neg'][2] & 1 << 28
    # implicit_h saturates at 4: bit 21 is the highest in the span
    assert _prim_probe('implicit_h', 4, False)['neg'][2] == _prim_probe('implicit_h', 7, False)['neg'][2]
    assert not _prim_probe('implicit_h', 4, False)['neg'][2] & (1 << 21)
    assert _prim_probe('implicit_h', 4, False)['neg'][2] & 0x3E0000 == 0x3E0000 ^ (1 << 21)
    # total_h saturates at 5: bit 32 is the highest in the span
    assert _prim_probe('total_h', 5, False)['neg'][2] == _prim_probe('total_h', 9, False)['neg'][2]
    assert not _prim_probe('total_h', 5, False)['neg'][2] & (1 << 32)
    assert _prim_probe('total_h', 5, False)['neg'][2] & 0x1F8000000 == 0x1F8000000 ^ (1 << 32)


def test_isotope_is_a_delta_from_the_mdl_common_isotope():
    # 13C: MDL common carbon is 12, delta 1, and _bit_of(1, -8, 8) is 9 -> bit 46 + 9
    r = _prim_probe_with_element('isotope', 13, False, 6)
    assert not r['neg'][2] & 1 << 55
    assert r['neg'][2] & 0xFFFFC00000000000 == 0xFFFFC00000000000 ^ (1 << 55)


def test_the_common_isotope_is_the_middle_of_the_span():
    r = _prim_probe_with_element('isotope', 12, False, 6)
    assert not r['neg'][2] & 1 << 54, 'delta 0 -> bit 46 + 8'


def test_no_isotope_is_its_own_primitive_at_the_top_of_the_span():
    from chython.core._core import _prim_probe
    r = _prim_probe('no_isotope', 0, False)
    assert r['neg'][2] & 0xFFFFC00000000000 == 0xFFFFC00000000000 ^ (1 << 63)


def test_isotope_without_a_settled_element_is_rejected():
    from chython.core._core import _prim_probe
    with pytest.raises(ValueError, match='isotope'):
        _prim_probe('isotope', 13, False)


def test_hybridization_and_aromatic_share_a_span():
    from chython.core._core import _prim_probe
    assert _prim_probe('hybridization', 4, False)['neg'][3] & 0x3F == 0x3F ^ (1 << 3)
    assert _prim_probe('hybridization', 1, False)['neg'][3] & 0x3F == 0x3F ^ 1


def test_the_hybridization_span_covers_all_six_states():
    # derive_scalars emits 5 (cumulated doubles) and 6 (unmatched combination) alongside 1-4, so the
    # span is six bits wide, not four.  A four-bit span would let every z5/z6 atom slip past every
    # z demand, since no bit in the mask would be set in its feature word.
    from chython.core._core import _prim_probe
    for z in range(1, 7):
        assert _prim_probe('hybridization', z, False)['neg'][3] & 0x3F == 0x3F ^ (1 << (z - 1))
    with pytest.raises(ValueError, match='hybridization value 7 is out of range'):
        _prim_probe('hybridization', 7, False)
    with pytest.raises(ValueError, match='hybridization value 0 is out of range'):
        _prim_probe('hybridization', 0, False)


def test_positive_ring_size_becomes_an_any_entry():
    from chython.core._core import _prim_probe
    r = _prim_probe('ring_size', 5, False)
    assert r['neg'][3] & 0x7FFFFFC00000 == 0, 'a multi-hot span cannot be forbidden into'
    assert r['any'] == ((3, 1 << 27),), 'bit 22 + 5'
    assert r['touched'][3] & 0x7FFFFFC00000


def test_negated_ring_size_is_a_plain_forbidden_bit():
    from chython.core._core import _prim_probe
    r = _prim_probe('ring_size', 5, True)
    assert r['any'] == ()
    assert r['neg'][3] & 0x7FFFFFC00000 == 1 << 27


def test_large_ring_sizes_fall_into_buckets():
    from chython.core._core import _prim_probe
    assert _prim_probe('ring_size', 30, False)['any'] == ((3, 1 << 22),)
    assert _prim_probe('ring_size', 40, False)['any'] == ((3, 1 << 23),)
    assert _prim_probe('ring_size', 60, False)['any'] == ((3, 1 << 24),)


def test_bare_ring_count_means_at_least_one_ring():
    from chython.core._core import _prim_probe
    r = _prim_probe('ring_count', 0, True)          # '!R0' is 'R': not acyclic
    assert r['neg'][3] & 0xFF800000000000 == 1 << 47
    acyclic = _prim_probe('ring_count', 0, False)   # '!R'
    assert acyclic['neg'][3] & 0xFF800000000000 == 0xFF800000000000 ^ (1 << 47)


def test_bond_order_span():
    from chython.core._core import _prim_probe
    order_span = 0xB800000000000000
    single = _prim_probe('bond_order', 1, False)['neg'][0] & order_span
    assert single == order_span ^ (1 << 59), 'a single bond forbids the rest of the span'
    dative = _prim_probe('bond_order', 8, False)['neg'][0] & order_span
    assert dative == order_span ^ (1 << 63)


def test_aromatic_bond_forbids_both_non_aromatic_topologies():
    from chython.core._core import _prim_probe
    r = _prim_probe('bond_aromatic', 0, False)
    assert r['neg'][0] & 0x4600000000000000 == (1 << 58) | (1 << 62)


def test_ring_bond_forbids_only_the_acyclic_bit():
    from chython.core._core import _prim_probe
    r = _prim_probe('bond_ring', 0, False)
    assert r['neg'][0] & 0x4600000000000000 == 1 << 58
    assert _prim_probe('bond_ring', 0, True)['neg'][0] & 0x4600000000000000 == \
        (1 << 57) | (1 << 62)


def test_a_stereo_primitive_puts_no_SIGN_bit_in_a_box():
    # A stored parity is a statement in the molecule's ruling-F26 frame and a query's is in its own, so
    # feature word IV bit 6 cannot screen the sign -- it would reject true matches whose two frames
    # differ by an odd permutation.  The sign travels through QSEG_STEREO instead, and bit 6 has to
    # come back untouched, or boxes_merge and box_unsatisfiable would both be reading a bit that
    # means nothing to them.
    #
    # The one demand that IS frame-free is bits 7-8, "a parity is configured" against
    # "none is".  Stated as the forbidding of bit 8 over that span, which is how every one-hot demand
    # is stated, and it is what puts bit 7 into the query signature.
    from chython.core._core import _prim_probe
    for value, sign in ((1, 1), (2, 2)):
        r = _prim_probe('stereo', value, False)
        assert r['neg'] == (0, 0, 0, 1 << 8), 'not configured is forbidden; the sign is not a bit'
        assert r['touched'] == (0, 0, 0, 3 << 7)
        assert not r['neg'][3] & 1 << 6 and not r['touched'][3] & 1 << 6, 'bit 6 stays out'
        assert r['any'] == ()
        # ruling F87: what it does write for the SIGN is the box's own sign field, which is not a
        # feature bit and is read by the kernel rather than by any mask AND
        assert r['sign'] == sign


def test_a_stereo_primitive_validates_its_value_and_refuses_negation():
    from chython.core._core import _prim_probe
    with pytest.raises(ValueError, match='stereo value'):
        _prim_probe('stereo', 3, False)
    with pytest.raises(ValueError, match='stereo value'):
        _prim_probe('stereo', 0, False)
    with pytest.raises(ValueError, match='cannot be negated'):
        _prim_probe('stereo', 1, True)


# ---------------------------------------------------------------------------
# Task 7: compile_term tests
# ---------------------------------------------------------------------------

def P(name, value=0, negated=False):
    return ('prim', name, value, negated)


OR = ('or',)
AND_LOW = ('and_low',)
AND_HIGH = ('and_high',)


def test_a_single_primitive_is_one_box():
    boxes = _compile_probe([P('element', 6)])
    assert len(boxes) == 1


def test_comma_makes_two_boxes():
    boxes = _compile_probe([P('degree', 2), OR, P('degree', 3)], merge=False)
    assert len(boxes) == 2
    # each box allows exactly one degree
    span = 0x1FE00
    assert boxes[0]['neg'][2] & span == span ^ (1 << 11)
    assert boxes[1]['neg'][2] & span == span ^ (1 << 12)


def test_semicolon_ands_into_one_box():
    boxes = _compile_probe([P('element', 6), AND_LOW, P('degree', 2)])
    assert len(boxes) == 1
    assert boxes[0]['neg'][0] & 0x01FFFFFFFFFFFFFF == 0x01FFFFFFFFFFFFFF ^ (1 << 51)
    assert boxes[0]['neg'][2] & 0x1FE00 == 0x1FE00 ^ (1 << 11)


def test_semicolon_distributes_over_comma():
    # [C;D2,D3] -> two boxes, both carbon
    boxes = _compile_probe([P('element', 6), AND_LOW, P('degree', 2), OR, P('degree', 3)],
                           merge=False)
    assert len(boxes) == 2
    carbon = 0x01FFFFFFFFFFFFFF ^ (1 << 51)
    assert all(b['neg'][0] & 0x01FFFFFFFFFFFFFF == carbon for b in boxes)
    assert {b['neg'][2] & 0x1FE00 for b in boxes} == {0x1FE00 ^ (1 << 11), 0x1FE00 ^ (1 << 12)}


def test_two_semicolon_groups_cross_multiply():
    # [C,N;D2,D3] -> 4 boxes
    boxes = _compile_probe([P('element', 6), OR, P('element', 7), AND_LOW,
                            P('degree', 2), OR, P('degree', 3)], merge=False)
    assert len(boxes) == 4
    pairs = {(b['neg'][0] & 0x01FFFFFFFFFFFFFF, b['neg'][2] & 0x1FE00) for b in boxes}
    assert len(pairs) == 4, 'every combination appears exactly once'


def test_and_high_binds_tighter_than_comma():
    # C&D2,N -> (carbon AND degree 2) OR nitrogen: two boxes, only the first constrains degree
    boxes = _compile_probe([P('element', 6), AND_HIGH, P('degree', 2), OR, P('element', 7)])
    assert len(boxes) == 2
    assert boxes[0]['neg'][2] & 0x1FE00 == 0x1FE00 ^ (1 << 11)
    assert boxes[1]['neg'][2] & 0x1FE00 == 0, 'the nitrogen branch says nothing about degree'


def test_charge_defaults_to_neutral_when_untouched():
    boxes = _compile_probe([P('element', 6)])
    charge_span = 0x3FFE00000000
    assert boxes[0]['neg'][2] & charge_span == charge_span ^ (1 << 37), 'default charge is 0'


def test_an_explicit_charge_suppresses_the_default():
    boxes = _compile_probe([P('element', 6), AND_LOW, P('charge', 1)])
    charge_span = 0x3FFE00000000
    assert boxes[0]['neg'][2] & charge_span == charge_span ^ (1 << 38)


def test_the_default_is_per_box_not_per_atom():
    # [C;+,D2] -> the '+' box keeps charge +1, the 'D2' box gets the neutral default
    boxes = _compile_probe([P('element', 6), AND_LOW, P('charge', 1), OR, P('degree', 2)])
    charge_span = 0x3FFE00000000
    charged = [b for b in boxes if b['neg'][2] & charge_span == charge_span ^ (1 << 38)]
    neutral = [b for b in boxes if b['neg'][2] & charge_span == charge_span ^ (1 << 37)]
    assert len(charged) == 1 and len(neutral) == 1


def test_a_negated_charge_also_suppresses_the_default():
    # [C;!+1] must allow 0 and -1 and +2, not just 0
    boxes = _compile_probe([P('element', 6), AND_LOW, P('charge', 1, True)])
    charge_span = 0x3FFE00000000
    assert boxes[0]['neg'][2] & charge_span == 1 << 38, 'only +1 is forbidden'


def test_radical_defaults_to_not_a_radical():
    boxes = _compile_probe([P('element', 6)])
    default_neg1 = boxes[0]['neg'][1]
    assert default_neg1 & 0xC000000000000000 == 1 << 63
    # the default demands exactly what an explicitly negated radical primitive demands
    assert (default_neg1 & 0xC000000000000000 ==
            _prim_probe('radical', 0, True)['neg'][1] & 0xC000000000000000)


def test_counts_and_topology_default_to_unconstrained():
    boxes = _compile_probe([P('element', 6)])
    assert boxes[0]['neg'][2] & 0x1FE00 == 0, 'degree is free'
    assert boxes[0]['neg'][2] & 0x3E0000 == 0, 'implicit H is free'
    assert boxes[0]['neg'][3] & 0xFF800000000000 == 0, 'ring count is free'
    assert boxes[0]['neg'][3] & 0x3F == 0, 'hybridization is free'


def test_boxes_differing_in_one_span_merge():
    # [C;D2,D3] can collapse to one box: same everything, union of two degree bits
    boxes = _compile_probe([P('element', 6), AND_LOW, P('degree', 2), OR, P('degree', 3)])
    assert len(boxes) == 1, 'the two degree boxes merge into one'
    span = 0x1FE00
    assert boxes[0]['neg'][2] & span == span ^ ((1 << 11) | (1 << 12))


def test_a_sign_difference_blocks_a_merge_that_would_otherwise_happen():
    """Ruling F87, at the merge: the sign is not a feature bit, so no span loop would notice it.

    `[C;D2@,D3]` differs from the merging `[C;D2,D3]` above in nothing a SPAN_MASK covers -- only in
    which of the two disjuncts demands a configuration.  Merging them would put the `@` demand on
    the D3 arm, which is exactly the '[C;@,D3] refuses a parity-less carbon' defect.
    """
    boxes = _compile_probe([P('element', 6), AND_LOW, P('degree', 2), AND_HIGH, P('stereo', 1),
                            OR, P('degree', 3)])
    assert len(boxes) == 2, 'the merge the sign-free version performs must not happen here'
    assert sorted(b['sign'] for b in boxes) == [0, 1], 'and each box keeps its own demand'


def test_two_signs_anded_into_one_box_compile_to_a_contradiction_and_are_kept():
    """`[C;@;@@]` is unsatisfiable and still constructible.

    Pruning the box would leave the term empty, and an empty term is a ValueError -- which ruling F87
    forbids here for the same reason F77 case 2 gives: a query that cannot be satisfied is not a
    construction error.  The kernel refuses it at match time instead
    (test_a_centre_cannot_satisfy_both_signs_at_once in test_stereo_query.py).
    """
    boxes = _compile_probe([P('element', 6), AND_LOW, P('stereo', 1), AND_LOW, P('stereo', 2)])
    assert len(boxes) == 1
    assert boxes[0]['sign'] == 3, 'QSIGN_CW | QSIGN_CCW on one box: no configuration satisfies it'


def test_boxes_differing_in_two_spans_do_not_merge():
    boxes = _compile_probe([P('element', 6), AND_HIGH, P('degree', 2), OR,
                            P('element', 7), AND_HIGH, P('degree', 3)])
    assert len(boxes) == 2


def test_the_spec_example_compiles_to_exactly_two_boxes():
    # [C,N&+;D4]: (C) OR (N AND +1), then AND D4 -- two boxes, both degree 4
    boxes = _compile_probe([P('element', 6), OR, P('element', 7), AND_HIGH, P('charge', 1),
                            AND_LOW, P('degree', 4)])
    assert len(boxes) == 2
    assert all(b['neg'][2] & 0x1FE00 == 0x1FE00 ^ (1 << 13) for b in boxes), 'D4 -> bit 9 + 4'
    charge_span = 0x3FFE00000000
    charges = {b['neg'][2] & charge_span for b in boxes}
    assert charges == {charge_span ^ (1 << 37), charge_span ^ (1 << 38)}, 'C neutral, N cationic'


def test_saturation_survives_the_pipeline():
    # D7 is the top of the degree span, so it also admits degree 9
    boxes = _compile_probe([P('element', 6), AND_LOW, P('degree', 7)])
    assert boxes[0]['neg'][2] & 0x1FE00 == 0x1FE00 ^ (1 << 16)


def test_an_unsatisfiable_box_is_dropped():
    # [C;N,D2] -> the C-and-N box is impossible; the C-and-D2 box survives
    boxes = _compile_probe([P('element', 6), AND_LOW, P('element', 7), OR, P('degree', 2)])
    assert len(boxes) == 1
    assert boxes[0]['neg'][2] & 0x1FE00 == 0x1FE00 ^ (1 << 11)


def test_a_wholly_unsatisfiable_term_is_an_error():
    with pytest.raises(ValueError, match='never match'):
        _compile_probe([P('element', 6), AND_LOW, P('element', 7)])


def test_ring_fusion_is_expressible():
    # [c;r5;r6] -- two positive multi-hot demands in one box
    boxes = _compile_probe([P('element', 6), AND_LOW, P('hybridization', 4), AND_LOW,
                            P('ring_size', 5), AND_LOW, P('ring_size', 6)])
    assert len(boxes) == 1
    assert boxes[0]['any'] == ((3, 1 << 27), (3, 1 << 28))


def test_any_lists_survive_the_cross_product():
    boxes = _compile_probe([P('ring_size', 5), OR, P('ring_size', 6), AND_LOW,
                            P('hybridization', 4)])
    assert len(boxes) == 2
    assert {b['any'] for b in boxes} == {((3, 1 << 27),), ((3, 1 << 28),)}


def test_negated_element_compiles():
    # [!C] — only carbon's bit is forbidden; all other elements are allowed
    boxes = _compile_probe([P('element', 6, True)])
    assert len(boxes) == 1


def test_metal_query_compiles():
    # [M] — non-metals forbidden, metals allowed
    boxes = _compile_probe([P('metal', 0, False)])
    assert len(boxes) == 1


def test_negated_metal_compiles():
    # [!M] — metals forbidden, non-metals allowed
    boxes = _compile_probe([P('metal', 0, True)])
    assert len(boxes) == 1


def test_element_compatible_with_metal_compiles():
    # [Fe;M] — Fe is a metal, so its bit survives both constraints
    boxes = _compile_probe([P('element', 26), AND_LOW, P('metal', 0, False)])
    assert len(boxes) == 1


def test_element_and_negated_different_element_compiles():
    # [C;!N] — carbon's bit is allowed; forbidding nitrogen changes nothing for carbon
    boxes = _compile_probe([P('element', 6), AND_LOW, P('element', 7, True)])
    assert len(boxes) == 1


def test_element_self_negated_is_an_error():
    # [C;!C] — every element bit ends up forbidden
    with pytest.raises(ValueError, match='never match'):
        _compile_probe([P('element', 6), AND_LOW, P('element', 6, True)])


def test_heavy_element_compiles():
    # [U] — uranium is Z 92, so its identity lives in word 1 bit 35 and word 0 keeps only the
    # heavy-marker bit 0 allowed.  Nothing is contradictory here.
    boxes = _compile_probe([P('element', 92)])
    assert len(boxes) == 1


def test_a_negated_heavy_element_spares_the_other_heavy_elements():
    # [!U] must forbid uranium's word-1 identity bit and nothing else.  Word 0 bit 0 is the
    # heavy marker that every element above 56 sets, so forbidding it there would silently
    # reject thorium, lanthanum and the rest of the heavy block.
    box = _prim_probe('element', 92, True)
    assert box['neg'][0] & 1 == 0, 'the shared heavy marker must stay allowed'
    assert box['neg'][1] & (1 << 35) != 0, "uranium's own identity bit is forbidden"
    assert box['neg'][1] & (1 << 33) == 0, "thorium's identity bit must stay allowed"
    assert box['neg'][1] == 1 << 35, 'no other element bit is touched'


def test_heavy_element_self_negated_is_an_error():
    # [U;!U] — this is the only case that distinguishes the two halves of the satisfiability
    # test.  Word 0's heavy-marker bit stays *allowed* (a heavy element could still match word 0),
    # and it is word 1 losing its last heavy identity that empties the set.  A predicate joining
    # the two halves with AND instead of OR compiles this and matches nothing at runtime;
    # [C;!C] forbids both, so it cannot tell the two spellings apart.
    with pytest.raises(ValueError, match='never match'):
        _compile_probe([P('element', 92), AND_LOW, P('element', 92, True)])


def test_heavy_and_light_element_is_an_error():
    # [C;U] — carbon forbids the heavy marker, uranium forbids every light bit
    with pytest.raises(ValueError, match='never match'):
        _compile_probe([P('element', 6), AND_LOW, P('element', 92)])


def test_negated_ring_sizes_are_not_merged():
    # [C;D2,D3;!r5,!r6] — the cross product produces 4 boxes pre-merge; verify with merge=False.
    # With merge=True the degree dimension collapses to 2 boxes (D23×!r5, D23×!r6), but the
    # uncovered-bits guard must keep those two separate — their ring-size bits differ in the
    # region of word 3 that is outside SPAN_COVERED.
    boxes = _compile_probe([P('element', 6), AND_LOW,
                            P('degree', 2), OR, P('degree', 3), AND_LOW,
                            P('ring_size', 5, True), OR, P('ring_size', 6, True)],
                           merge=False)
    assert len(boxes) == 4, 'four degree×negated-ring-size combinations pre-merge'
    # Now verify the merged form: 2 boxes, each carrying exactly one ring-size bit
    boxes_merged = _compile_probe([P('element', 6), AND_LOW,
                                   P('degree', 2), OR, P('degree', 3), AND_LOW,
                                   P('ring_size', 5, True), OR, P('ring_size', 6, True)])
    assert len(boxes_merged) == 2, 'degree merges within each ring group; ring-size distinction survives'
    ring_bits = {b['neg'][3] & 0x7FFFFFC00000 for b in boxes_merged}
    assert ring_bits == {1 << 27, 1 << 28}, \
        'both !r5 and !r6 bit patterns must appear as separate boxes'


def test_positive_ring_sizes_do_not_merge():
    # [C;D2,D3;r5,r6] — positive ring sizes live in the any lists (not neg), so the
    # uncovered-bits guard cannot help here; the any-list guard holds the line.  The two
    # surviving boxes have identical neg in all four words, so deleting that guard makes
    # diff_cnt == 0 dedup them into one box and turns an AND of demands into an OR.
    # Pre-merge: 4 boxes. Post-merge: 2 (degree merges within each ring group;
    # any-list guard prevents those 2 from collapsing to 1).
    boxes_pre = _compile_probe([P('element', 6), AND_LOW,
                                P('degree', 2), OR, P('degree', 3), AND_LOW,
                                P('ring_size', 5), OR, P('ring_size', 6)],
                               merge=False)
    assert len(boxes_pre) == 4, 'four boxes pre-merge'
    boxes = _compile_probe([P('element', 6), AND_LOW,
                            P('degree', 2), OR, P('degree', 3), AND_LOW,
                            P('ring_size', 5), OR, P('ring_size', 6)])
    assert len(boxes) == 2, 'degree merges; any-list guard blocks the final 2→1'


def test_a_trailing_operator_is_an_error():
    with pytest.raises(ValueError, match='malformed'):
        _compile_probe([P('element', 6), AND_LOW])


def test_two_operators_in_a_row_are_an_error():
    with pytest.raises(ValueError, match='malformed'):
        _compile_probe([P('element', 6), OR, AND_LOW, P('element', 7)])


def test_an_empty_term_is_an_error():
    with pytest.raises(ValueError, match='malformed'):
        _compile_probe([])


# ---------------------------------------------------------------------------
# Seal — Task 8
# ---------------------------------------------------------------------------

def atoms_chain(n, element=6):
    """A linear chain of n atoms, stable ids 1..n, all carbon, single bonds."""
    ops = []
    for i in range(1, n + 1):
        ops.append(('atom', i))
        ops.append(('token', i, 'element', element, False))
    for i in range(1, n):
        ops.append(('bond', i, i + 1))
    return ops


def test_a_single_atom_seals_to_one_component_one_root():
    info = _seal_probe([('atom', 1), ('token', 1, 'element', 6, False)])
    assert info['atom_count'] == 1
    assert info['bond_count'] == 0
    assert info['components'] == [(0, 1, -1)]
    assert info['roots'] == [0]
    assert info['back'] == [0]


def test_a_chain_orders_atoms_so_every_atom_touches_an_earlier_one():
    info = _seal_probe(atoms_chain(4))
    assert info['atom_count'] == 4
    for position in range(1, 4):
        assert info['back'][position] < position
    assert info['closure_count'] == 0, 'a tree has no closures'


def test_two_fragments_become_two_components():
    ops = atoms_chain(2)
    ops += [('atom', 3), ('token', 3, 'element', 8, False)]
    info = _seal_probe(ops)
    assert info['component_count'] == 2
    assert [c[:2] for c in info['components']] == [(0, 2), (2, 3)]
    assert len(info['roots']) == 2


def test_the_rarest_element_wins_the_root():
    # C-C-O: oxygen is the most constrained atom, so the DFS seeds from oxygen
    ops = [('atom', 1), ('token', 1, 'element', 6, False),
           ('atom', 2), ('token', 2, 'element', 6, False),
           ('atom', 3), ('token', 3, 'element', 8, False),
           ('bond', 1, 2), ('bond', 2, 3)]
    info = _seal_probe(ops)
    assert info['order'][0] == 2, 'slot 2 is the oxygen'


def test_an_unconstrained_element_never_wins_the_root():
    # [A]-C: the any-atom has infinite element demand, so carbon roots
    ops = [('atom', 1), ('token', 1, 'any', 0, False),
           ('atom', 2), ('token', 2, 'element', 6, False),
           ('bond', 1, 2)]
    info = _seal_probe(ops)
    assert info['order'][0] == 1


def test_a_ring_produces_exactly_one_closure():
    ops = atoms_chain(6)
    ops.append(('bond', 6, 1))
    info = _seal_probe(ops)
    assert info['closure_count'] == 1
    to_position, _ = info['closures'][0]
    assert to_position == 0, 'the closure points back at the root'


def test_two_fused_rings_produce_two_closures():
    # naphthalene skeleton: 10 atoms, 11 bonds, 10 tree bonds
    ops = atoms_chain(10)
    ops.append(('bond', 10, 1))
    ops.append(('bond', 5, 10))
    info = _seal_probe(ops)
    assert info['bond_count'] == 11
    assert info['closure_count'] == 2


def test_a_tree_bond_folds_its_boxes_into_the_atom():
    ops = [('atom', 1), ('token', 1, 'element', 6, False),
           ('atom', 2), ('token', 2, 'element', 8, False),
           ('bond', 1, 2), ('btoken', 1, 2, 'bond_order', 2, False)]
    info = _seal_probe(ops)
    assert info['bond_box_counts'] == [0], 'the tree bond has no boxes left of its own'
    # the second position's box now forbids everything but a double bond
    order_span = 0xB800000000000000
    position = 1
    box = info['boxes'][position][0]
    assert box['neg'][0] & order_span == order_span ^ (1 << 60)


def test_a_closure_bond_keeps_its_boxes():
    ops = atoms_chain(6)
    ops.append(('bond', 6, 1))
    ops.append(('btoken', 6, 1, 'bond_order', 2, False))
    info = _seal_probe(ops)
    assert sum(info['bond_box_counts']) > 0, 'the closure bond keeps its own boxes'


def test_folding_multiplies_box_counts():
    # [C,N] with a -,= bond: 2 atom boxes x 2 bond boxes = 4 pre-merge, and (C, single) /
    # (C, double) differ in exactly the order span, so they merge -- as do the two N boxes.
    ops = [('atom', 1), ('token', 1, 'element', 6, False),
           ('atom', 2), ('token', 2, 'element', 6, False), ('op', 2, 'or'),
           ('token', 2, 'element', 7, False),
           ('bond', 1, 2), ('btoken', 1, 2, 'bond_order', 1, False), ('bop', 1, 2, 'or'),
           ('btoken', 1, 2, 'bond_order', 2, False)]
    info = _seal_probe(ops)
    assert len(info['boxes'][1]) == 2, 'the two bond orders merge, the two elements do not'


def test_an_implicit_bond_means_single():
    ops = [('atom', 1), ('token', 1, 'element', 6, False),
           ('atom', 2), ('token', 2, 'element', 6, False), ('bond', 1, 2)]
    info = _seal_probe(ops)
    order_span = 0xB800000000000000
    assert info['boxes'][1][0]['neg'][0] & order_span == order_span ^ (1 << 59)


def test_element_demand_counts_atoms_per_element():
    ops = [('atom', 1), ('token', 1, 'element', 6, False),
           ('atom', 2), ('token', 2, 'element', 6, False),
           ('atom', 3), ('token', 3, 'element', 8, False),
           ('bond', 1, 2), ('bond', 2, 3)]
    demand = _seal_probe(ops)['element_demand']
    assert demand[6] == 2
    assert demand[8] == 1
    assert demand[7] == 0


def test_an_unconstrained_atom_demands_no_element():
    ops = [('atom', 1), ('token', 1, 'any', 0, False)]
    assert _seal_probe(ops)['element_demand'] == [0] * 120


def test_groups_land_on_components():
    ops = atoms_chain(2)
    ops += [('atom', 3), ('token', 3, 'element', 8, False)]
    ops += [('group', 1, 0), ('group', 2, 0), ('group', 3, 0)]
    info = _seal_probe(ops)
    assert [c[2] for c in info['components']] == [0, 0]
    assert info['flags'] & 1, 'QFLAG_HAS_GROUP'


def test_masked_and_map_survive_the_reorder():
    ops = atoms_chain(3)
    ops += [('masked', 3), ('map', 3, 42)]
    info = _seal_probe(ops)
    position = info['order'].index(2)          # slot 2 is stable id 3
    assert info['masked'][position] is True
    assert info['map_numbers'][position] == 42
    assert info['flags'] & 4, 'QFLAG_HAS_MASKED'


def test_a_bond_to_an_unknown_atom_is_an_error():
    with pytest.raises(ValueError, match='unknown atom'):
        _seal_probe([('atom', 1), ('token', 1, 'element', 6, False), ('bond', 1, 9)])


def test_a_self_loop_is_an_error():
    with pytest.raises(ValueError, match='self'):
        _seal_probe([('atom', 1), ('token', 1, 'element', 6, False), ('bond', 1, 1)])


def test_a_duplicate_bond_is_an_error():
    ops = atoms_chain(2)
    ops.append(('bond', 2, 1))
    with pytest.raises(ValueError, match='duplicate'):
        _seal_probe(ops)


def test_an_atom_without_primitives_is_an_error():
    with pytest.raises(ValueError, match='no primitives'):
        _seal_probe([('atom', 1)])


def test_the_offending_atom_is_named_in_a_compile_error():
    with pytest.raises(ValueError, match='atom 1'):
        _seal_probe([('atom', 1), ('token', 1, 'element', 6, False), ('op', 1, 'and_low'),
                     ('token', 1, 'element', 7, False)])


def test_a_bond_token_for_a_pair_that_is_not_a_bond_is_an_error():
    # Without this check the constraint the caller wrote would vanish with no error at all.
    ops = atoms_chain(3)
    ops.append(('btoken', 1, 3, 'bond_order', 2, False))
    with pytest.raises(ValueError, match='unknown bond'):
        _seal_probe(ops)


def test_a_bond_operator_for_a_pair_that_is_not_a_bond_is_an_error():
    ops = atoms_chain(3)
    ops.append(('bop', 1, 3, 'or'))
    with pytest.raises(ValueError, match='unknown bond'):
        _seal_probe(ops)


def test_a_component_spanning_two_groups_is_an_error():
    # A group is per-atom in the journal but per-component in the arena, so two atoms of one
    # component carrying different groups has no representation.
    ops = atoms_chain(2)
    ops += [('group', 1, 0), ('group', 2, 1)]
    with pytest.raises(ValueError, match='two groups'):
        _seal_probe(ops)


def test_only_the_closure_bond_carries_boxes():
    ops = atoms_chain(6)
    ops.append(('bond', 6, 1))
    ops.append(('btoken', 6, 1, 'bond_order', 2, False))
    info = _seal_probe(ops)
    # bond slot 5 is the (6, 1) bond: the only one the DFS did not use as a tree edge
    assert info['bond_box_counts'] == [0, 0, 0, 0, 0, 1]
    _, bond_slot = info['closures'][0]
    assert bond_slot == 5, 'the closure names the bond whose boxes survived'


def test_an_unknown_operator_name_raises_key_error():
    with pytest.raises(KeyError):
        _seal_probe([('atom', 1), ('token', 1, 'element', 6, False), ('op', 1, 'nand')])


def test_a_closure_bond_box_is_the_unfolded_bond_term():
    # The closure's boxes must be the bond's own -- an edge-word test, not the atom fold.  If
    # the emit pass wrote atom boxes into the bond-box segment, the element span would be set
    # here and the order span would carry the fold's three forbidden bits instead of one.
    ops = atoms_chain(6)
    ops.append(('bond', 6, 1))
    ops.append(('btoken', 6, 1, 'bond_order', 2, False))
    info = _seal_probe(ops)
    order_span = 0xB800000000000000
    element_span = 0x01FFFFFFFFFFFFFF
    box = info['bond_boxes'][5][0]
    assert box['neg'][0] & order_span == order_span ^ (1 << 60)
    assert box['neg'][0] & element_span == 0, 'a bond box constrains no element'


def test_a_lone_stereo_primitive_seals_with_a_frame_no_target_can_satisfy():
    # A stereo atom with no neighbours seals -- an unsatisfiable query is not a construction error --
    # and records n_refs = 0, which the kernel refuses at match time (ruling F77 case 2).  Readiness
    # is the atom's own position, because a frame with no directions is complete as soon as the
    # anchor is bound.  'kind' is SU_TETRA and 'demand' 0: a sign's own state is in the box it came
    # from, and only a geometry's is in the record.
    info = _seal_probe([('atom', 1), ('token', 1, 'stereo', 1, False)])
    assert info['flags'] & 2, 'QFLAG_HAS_STEREO'
    assert info['stereo'] == [{'position': 0, 'readiness': 0, 'refs': (0xFFFFFFFF,) * 4, 'sign': 1,
                               'n_refs': 0, 'kind': 0, 'demand': 0}]


def test_the_stereo_record_holds_the_querys_own_f26_order():
    # C(F)(Cl)(Br) with '@@' on the carbon, the fluorine added LAST so that ascending query slot
    # and the plan's own order cannot coincide.  refs must be the neighbours by ascending slot --
    # creation order -- mapped through pos_of, and readiness the greatest of those positions and
    # the anchor's.
    ops = [('atom', 1), ('token', 1, 'element', 6, False), ('op', 1, 'and_low'),
           ('token', 1, 'stereo', 2, False),
           ('atom', 2), ('token', 2, 'element', 17, False),
           ('atom', 3), ('token', 3, 'element', 35, False),
           ('atom', 4), ('token', 4, 'element', 9, False),
           ('bond', 1, 2), ('btoken', 1, 2, 'bond_order', 1, False),
           ('bond', 1, 3), ('btoken', 1, 3, 'bond_order', 1, False),
           ('bond', 1, 4), ('btoken', 1, 4, 'bond_order', 1, False)]
    info = _seal_probe(ops)
    assert len(info['stereo']) == 1
    rec = info['stereo'][0]
    assert rec['sign'] == 2 and rec['n_refs'] == 3
    # info['order'][p] is the query slot at plan position p, so this inverts it.
    pos_of_slot = {slot: p for p, slot in enumerate(info['order'])}
    assert rec['position'] == pos_of_slot[0]
    assert rec['refs'][:3] == (pos_of_slot[1], pos_of_slot[2], pos_of_slot[3]), \
        'slots 1, 2, 3 are Cl, Br, F in creation order'
    assert rec['refs'][3] == 0xFFFFFFFF, 'the unnamed direction needs no entry'
    assert rec['readiness'] == max(rec['position'], *rec['refs'][:3])
    assert rec['readiness'] > rec['position'], 'the frame completes after the anchor is bound'


def test_a_two_element_atom_demands_neither_element():
    # [C,N]-C: only the plain carbon reaches the histogram.  QSEG_ELEMENT_DEMAND has to be a
    # sound lower bound, so an atom counts toward slot e only when EVERY box of its disjunction
    # allows exactly e.  An implementation that took the first box, or read wbox_t.element_single,
    # would score demand[6] == 2 -- and the screen would then reject a molecule holding one
    # carbon and one nitrogen, a false negative with no error anywhere.
    ops = [('atom', 1), ('token', 1, 'element', 6, False), ('op', 1, 'or'),
           ('token', 1, 'element', 7, False),
           ('atom', 2), ('token', 2, 'element', 6, False),
           ('bond', 1, 2)]
    demand = _seal_probe(ops)['element_demand']
    assert demand[6] == 1, '[C,N] guarantees no carbon, so only the plain carbon counts'
    assert demand[7] == 0


def test_a_heavy_element_lands_in_its_own_demand_slot():
    # Word 0 bit 56 is the light/heavy boundary and it has already produced one Critical in this
    # epic.  Uranium's identity lives in word 1 bit 35, so the histogram has to read it from
    # there and add 57 back; getting that wrong lands the atom in some other slot, or none.
    demand = _seal_probe([('atom', 1), ('token', 1, 'element', 92, False)])['element_demand']
    assert demand[92] == 1
    assert sum(demand) == 1, 'uranium demands uranium and nothing else'


def test_a_heavy_element_is_no_more_constrained_than_a_light_one():
    # C-[U]: both atoms allow exactly one element and both have degree 1, so the whole tie-break
    # chain falls through to the lower slot.  An element-cardinality function blind to the heavy
    # span would score [U] as allowing zero elements and hand it the root.
    ops = [('atom', 1), ('token', 1, 'element', 6, False),
           ('atom', 2), ('token', 2, 'element', 92, False),
           ('bond', 1, 2)]
    assert _seal_probe(ops)['order'][0] == 0


# ---------------------------------------------------------------------------
# The automorphism group, computed at seal from the pre-fold scratch
# ---------------------------------------------------------------------------

QFLAG_ASYMMETRIC = 8
QFLAG_PARTIAL_AUTOMORPHISM = 16


def test_the_automorphism_segment_is_sized_by_alloc():
    # Sized like every other segment rather than appended to a sealed arena: three rows of four
    # positions is 48 bytes, and nothing after it may overlap.
    info = _query_alloc_probe(4, 3, 4, 0, 0, 0, 1, 3)
    offset, length = info['segments'][_query_segment_ids()['QSEG_AUTOMORPHISM']]
    assert length == 3 * 4 * 4
    assert offset + length <= info['total_len']
    assert info['automorphism_count'] == 3


def test_a_symmetric_query_stores_its_rows_as_position_permutations():
    # C-C: one row, the swap.  Rows index DFS positions, so the row must be a permutation of
    # range(atom_count) -- storing slots instead would still read as a permutation here, which is
    # why the ring test below uses a query whose slot order and position order differ.
    ops = atoms_chain(2)
    info = _seal_probe(ops)
    assert info['automorphism_count'] == 1
    assert info['automorphisms'] == [(1, 0)]
    assert not info['flags'] & QFLAG_ASYMMETRIC


def test_an_asymmetric_query_is_flagged_and_stores_nothing():
    ops = [('atom', 1), ('token', 1, 'element', 6, False),
           ('atom', 2), ('token', 2, 'element', 8, False),
           ('bond', 1, 2)]
    info = _seal_probe(ops)
    assert info['automorphism_count'] == 0
    assert info['automorphisms'] == []
    assert info['flags'] & QFLAG_ASYMMETRIC, 'the refinement proved every class a singleton'


def test_every_stored_row_is_a_permutation_of_positions_and_never_the_identity():
    # cyclopropane: the full S3, five rows.  The oxygen-free ring makes the DFS pick slot 0 as
    # root, so slot order and position order agree; what this pins is that the identity is
    # excluded and that no row repeats or leaves a position out.
    ops = atoms_chain(3)
    ops.append(('bond', 3, 1))
    info = _seal_probe(ops)
    assert info['automorphism_count'] == 5
    identity = tuple(range(3))
    assert len(set(info['automorphisms'])) == 5
    for row in info['automorphisms']:
        assert sorted(row) == list(identity)
        assert row != identity


def test_rows_are_position_permutations_not_slot_permutations():
    """C-C-O with the closure C-O: slot 2 is the oxygen and it roots the DFS, so position order
    is (2, 0, 1) or (2, 1, 0) -- never the identity on slots.  The group is the swap of the two
    carbons, which is slots {0, 1} and therefore positions {1, 2}: a row of (0, 2, 1).  Storing
    the slot permutation would write (1, 0, 2) instead, and mapping_is_canonical -- which indexes
    m.mapping by POSITION -- would then compare the oxygen against a carbon.
    """
    ops = [('atom', 1), ('token', 1, 'element', 6, False),
           ('atom', 2), ('token', 2, 'element', 6, False),
           ('atom', 3), ('token', 3, 'element', 8, False),
           ('bond', 1, 2), ('bond', 2, 3), ('bond', 3, 1)]
    info = _seal_probe(ops)
    assert info['order'][0] == 2, 'the oxygen roots the DFS'
    assert info['automorphisms'] == [(0, 2, 1)]


def test_the_row_cap_marks_the_group_partial():
    # Seven interchangeable disconnected carbons have 5039 non-identity automorphisms; the arena
    # stores Q_AUTOMORPHISM_MAX_ROWS of them and says so.  A partial group filters less, so this
    # is the safe failure mode.
    ops = []
    for i in range(1, 8):
        ops.append(('atom', i))
        ops.append(('token', i, 'element', 6, False))
    info = _seal_probe(ops)
    assert info['automorphism_count'] == 1024
    assert info['flags'] & QFLAG_PARTIAL_AUTOMORPHISM
    assert not info['flags'] & QFLAG_ASYMMETRIC
    for row in info['automorphisms']:
        assert sorted(row) == list(range(7))


def test_a_one_atom_query_has_a_trivial_group():
    info = _seal_probe([('atom', 1), ('token', 1, 'element', 6, False)])
    assert info['automorphism_count'] == 0
    assert info['flags'] & QFLAG_ASYMMETRIC


# ---------------------------------------------------------------------------
# Task 9: QueryContainer tests
# ---------------------------------------------------------------------------

def test_building_a_two_atom_query():
    q = QueryContainer()
    a = q.add_atom()
    b = q.add_atom()
    q.atom_primitive(a, 'element', 6)
    q.atom_primitive(b, 'element', 8)
    q.add_bond(a, b)
    assert len(q) == 2
    assert q.atom_count == 2
    assert q.bond_count == 1


def test_stable_ids_start_at_one_and_increment():
    q = QueryContainer()
    assert q.add_atom() == 1
    assert q.add_atom() == 2


def test_an_empty_query_refuses_to_seal():
    q = QueryContainer()
    with pytest.raises(ValueError, match='empty query'):
        q.atom_count_sealed()


def test_operators_interleave_with_primitives():
    # Two atoms connected by a bond: [C;D2,D3]-[O]. Both should compile to 1 box each.
    # D2 and D3 merge into one box (the main invariant), O is a single-box atom.
    # Using two atoms exercises box_counts() over all positions (not just the first).
    q = QueryContainer()
    a = q.add_atom()
    q.atom_primitive(a, 'element', 6)
    q.atom_operator(a, 'and_low')
    q.atom_primitive(a, 'degree', 2)
    q.atom_operator(a, 'or')
    q.atom_primitive(a, 'degree', 3)
    b = q.add_atom()
    q.atom_primitive(b, 'element', 8)
    q.add_bond(a, b)
    counts = q.box_counts()
    assert len(counts) == 2, 'two atoms, two positions'
    assert all(c == 1 for c in counts), 'D2 and D3 merge into one box; [O] is also one box'


def test_a_bad_primitive_name_is_rejected_at_append_time():
    q = QueryContainer()
    a = q.add_atom()
    with pytest.raises(KeyError):
        q.atom_primitive(a, 'nonsense', 1)


def test_a_bond_between_unknown_atoms_is_rejected_at_append_time():
    q = QueryContainer()
    q.add_atom()
    with pytest.raises(ValueError, match='unknown atom'):
        q.add_bond(1, 7)


def test_edit_discards_on_exception():
    q = QueryContainer()
    a = q.add_atom()
    q.atom_primitive(a, 'element', 6)
    with pytest.raises(RuntimeError):
        with q.edit():
            b = q.add_atom()
            q.atom_primitive(b, 'element', 8)
            q.add_bond(a, b)
            raise RuntimeError('nope')
    assert q.atom_count == 1
    assert q.bond_count == 0


def test_edit_commits_on_success():
    q = QueryContainer()
    a = q.add_atom()
    q.atom_primitive(a, 'element', 6)
    with q.edit():
        b = q.add_atom()
        q.atom_primitive(b, 'element', 8)
        q.add_bond(a, b)
    assert q.atom_count == 2


def test_sealing_twice_reuses_the_arena():
    q = QueryContainer()
    a = q.add_atom()
    q.atom_primitive(a, 'element', 6)
    assert q.seal_generation() == 0
    q.box_counts()
    assert q.seal_generation() == 1
    q.box_counts()
    assert q.seal_generation() == 1, 'the second call reuses the sealed arena'


def test_a_mutation_invalidates_the_seal():
    q = QueryContainer()
    a = q.add_atom()
    q.atom_primitive(a, 'element', 6)
    q.box_counts()
    b = q.add_atom()
    q.atom_primitive(b, 'element', 8)
    q.box_counts()
    assert q.seal_generation() == 2


def test_masked_and_map_number_round_trip():
    # Two atoms: only atom `a` is masked and has a map number.  Both map_numbers() and
    # masked_atoms() must iterate ALL atoms, so this two-atom case distinguishes "every"
    # from "any" -- a wrong implementation that only checks the first or last atom would
    # either include `b` or miss `a`.
    q = QueryContainer()
    a = q.add_atom()
    q.atom_primitive(a, 'element', 6)
    q.set_masked(a)
    q.set_map_number(a, 7)
    b = q.add_atom()
    q.atom_primitive(b, 'element', 8)
    # b is not masked and has no map number; it must not appear in either result
    assert q.map_numbers() == {a: 7}
    assert q.masked_atoms() == frozenset({a})


def test_the_element_wildcard_flags_are_derived_from_the_boxes_not_declared():
    # There is no `set_wildcard` and there must not be one: QATOM_ANY_ELEMENT and
    # QATOM_METAL_ELEMENT are computed at seal from the compiled term, so a query built through
    # the journal -- no reader anywhere near it -- carries them exactly as `[A]`/`[M]` would.
    # Fails against: a provenance flag journalled by the SMARTS reader, which would leave every
    # programmatically built wildcard invisible to `wildcard_atoms()`.
    q = QueryContainer()
    any_atom = q.add_atom()
    # PRIM_ANY, the ELEMENT wildcard the reader emits for `[A]` and for `*`.  It touches no span at
    # all, which is why its value is ignored and passed as 0 -- and why it is not `'any_charge'`,
    # PRIM_ANY_CHARGE, which touches the charge span to withdraw its default.
    q.atom_primitive(any_atom, 'any', 0)
    metal = q.add_atom()
    q.atom_primitive(metal, 'metal', 0)         # PRIM_METAL, `[M]`: value ignored likewise
    named = q.add_atom()
    q.atom_primitive(named, 'element', 6)
    counted = q.add_atom()
    q.atom_primitive(counted, 'degree', 2)      # names no element either, and says so
    assert q.wildcard_atoms() == {any_atom: 'any', metal: 'metal', counted: 'any'}


def test_a_masked_wildcard_reports_as_both():
    # The two live in one uint16 flags field, so a mask must not clear the wildcard bit or the
    # other way about.  `[M;M]` -- a metal, masked -- is the shape the reactor's context atoms take.
    q = QueryContainer()
    a = q.add_atom()
    q.atom_primitive(a, 'metal', 0)
    q.set_masked(a)
    assert q.masked_atoms() == frozenset({a})
    assert q.wildcard_atoms() == {a: 'metal'}


def test_a_stereo_primitive_seals_through_the_container():
    # The stereo primitive is in PRIM_NAMES so the append succeeds (no KeyError).
    # To actually reach prim_apply(PRIM_STEREO), the token stream must be valid:
    # two consecutive OPC_PRIM without an operator is malformed and fires a ValueError
    # before prim_apply runs.  The and_high operator is the implicit juxtaposition that
    # SMARTS uses for [C@@], so [C;and_high;stereo] is the correct minimal test.
    q = QueryContainer()
    a = q.add_atom()
    q.atom_primitive(a, 'element', 6)
    q.atom_operator(a, 'and_high')
    q.atom_primitive(a, 'stereo', 1)
    assert q.box_counts() == [1], 'the sign is not a box: [C@] compiles to the boxes of [C]'
    assert q.seal_generation() == 1


def test_seal_generation_does_not_increment_on_a_failed_seal():
    # A failing seal is needed, and a stereo primitive is not one -- it seals happily.  '[C;!C]' --
    # carbon and not carbon -- is: compile_term prunes its only box and raises.
    q = QueryContainer()
    a = q.add_atom()
    q.atom_primitive(a, 'element', 6)
    q.atom_operator(a, 'and_low')
    q.atom_primitive(a, 'element', 6, negated=True)
    with pytest.raises(ValueError, match='can never match'):
        q.box_counts()
    assert q.seal_generation() == 0
    with pytest.raises(ValueError, match='can never match'):
        q.box_counts()
    assert q.seal_generation() == 0, 'a failed seal leaves nothing behind to reuse'


def test_a_stereo_query_reports_its_automorphism_group_as_partial():
    """A stereo primitive is invisible to the automorphism search, so query_seal drops the rows.

    `[C@](F)(F)Cl`: exchanging the two fluorines is a graph automorphism and an ODD permutation of
    the centre's directions, so it inverts the very thing the primitive states.  mapping_is_canonical
    only ever REJECTS embeddings, so keeping such a row would drop a real match; discarding every row
    only over-reports duplicates.  Without this assertion the rows could be re-admitted with the
    suite green, which is why it is pinned on the flag and the count rather than on a match result.
    """
    ops = [('atom', 1), ('token', 1, 'element', 6, False),
           ('atom', 2), ('token', 2, 'element', 9, False),
           ('atom', 3), ('token', 3, 'element', 9, False),
           ('atom', 4), ('token', 4, 'element', 17, False),
           ('bond', 1, 2), ('btoken', 1, 2, 'bond_order', 1, False),
           ('bond', 1, 3), ('btoken', 1, 3, 'bond_order', 1, False),
           ('bond', 1, 4), ('btoken', 1, 4, 'bond_order', 1, False)]
    plain = _seal_probe(ops)
    assert plain['automorphism_count'] == 1, 'the two fluorines really do exchange'
    assert not plain['flags'] & 16, 'and the group is reported in full'

    stereo = _seal_probe(ops + [('op', 1, 'and_low'), ('token', 1, 'stereo', 1, False)])
    assert stereo['automorphism_count'] == 0
    assert stereo['automorphisms'] == []
    assert stereo['flags'] & 16, 'QFLAG_PARTIAL_AUTOMORPHISM: under-reported, not proven trivial'
    assert not stereo['flags'] & 8, 'and not QFLAG_ASYMMETRIC'


# ---------------------------------------------------------------------------
# Task 9 fix-round 1 tests
# ---------------------------------------------------------------------------

def test_box_counts_returns_per_position_values():
    # [C,N]-[O]: atom a has 2 boxes (C or N), atom b has 1 box (O).
    # Oxygen (1 element allowed) is rarer than C,N (2 elements allowed) so the DFS roots
    # at oxygen: position 0 = O (1 box), position 1 = C,N (2 boxes).
    # Fails against: atoms[0].box_count reused for every position (would give [1, 1]).
    # Also fails against reversed-order iteration (would give [2, 1]).
    q = QueryContainer()
    a = q.add_atom()
    q.atom_primitive(a, 'element', 6)
    q.atom_operator(a, 'or')
    q.atom_primitive(a, 'element', 7)
    b = q.add_atom()
    q.atom_primitive(b, 'element', 8)
    q.add_bond(a, b)
    assert q.box_counts() == [1, 2]


def test_atom_count_sealed_returns_correct_count():
    # return 0 passes all prior tests; this pins a non-zero return value.
    q = QueryContainer()
    a = q.add_atom()
    b = q.add_atom()
    q.atom_primitive(a, 'element', 6)
    q.atom_primitive(b, 'element', 8)
    q.add_bond(a, b)
    assert q.atom_count_sealed() == 2


def test_bond_primitive_lands_on_the_bond_not_the_atom():
    # Wrong impl: op.op = QOP_ATOM_TOKEN at bond_primitive line 277.
    # With QOP_ATOM_TOKEN, the bond_order primitive becomes an extra OPC_PRIM token for atom a,
    # creating a malformed two-OPC_PRIM-in-a-row sequence that raises ValueError at seal.
    # With correct QOP_BOND_TOKEN, the bond order folds into atom b's box: both atoms have
    # exactly 1 box.  The seal succeeding and returning [1, 1] is the observable.
    q = QueryContainer()
    a = q.add_atom()
    q.atom_primitive(a, 'element', 6)
    b = q.add_atom()
    q.atom_primitive(b, 'element', 8)
    q.add_bond(a, b)
    q.bond_primitive(a, b, 'bond_order', 2)
    assert q.box_counts() == [1, 1]


def test_bond_primitive_negated():
    # Negated bond primitive: forbids double bonds.  Seals successfully.
    q = QueryContainer()
    a = q.add_atom()
    q.atom_primitive(a, 'element', 6)
    b = q.add_atom()
    q.atom_primitive(b, 'element', 8)
    q.add_bond(a, b)
    q.bond_primitive(a, b, 'bond_order', 2, negated=True)
    assert q.box_counts() == [1, 1]


def test_bond_operator_lands_on_the_bond_not_the_atom():
    # Wrong impl: QOP_ATOM_TOKEN for bond_operator.  With that bug the OPC_OR token is
    # appended to atom a's token stream: atom a then has the sequence
    # [OPC_PRIM(element), OPC_PRIM(bond_order1), OPC_OR, OPC_PRIM(bond_order2)].
    # The two consecutive OPC_PRIM tokens without an operator between them are malformed
    # and query_seal raises ValueError('malformed query term').
    # With correct QOP_BOND_TOKEN the bond tokens are routed to the bond segment and the
    # atom token streams are both clean singletons → box_counts() returns [1, 1].
    q = QueryContainer()
    a = q.add_atom()
    q.atom_primitive(a, 'element', 6)
    b = q.add_atom()
    q.atom_primitive(b, 'element', 8)
    q.add_bond(a, b)
    q.bond_primitive(a, b, 'bond_order', 1)
    q.bond_operator(a, b, 'or')
    q.bond_primitive(a, b, 'bond_order', 2)
    assert q.box_counts() == [1, 1]


def test_set_group_writes_group_opcode_not_masked_or_map():
    # set_group is unobservable via query inspection (group is per-component in the arena
    # and no getter exists by design).  Pin it negatively: the op must not write QOP_SET_MASKED
    # or QOP_SET_MAP, which would corrupt masked_atoms() / map_numbers().
    # Fails against: set_group accidentally writing QOP_SET_MASKED (masked_atoms would be
    # non-empty) or QOP_SET_MAP (map_numbers would be non-empty).
    q = QueryContainer()
    a = q.add_atom()
    q.atom_primitive(a, 'element', 6)
    q.set_group(a, 5)
    assert q.masked_atoms() == frozenset()
    assert q.map_numbers() == {}


def test_atom_primitive_negated():
    # negated=True was not exercised through the container.  [!C] forbids carbon.
    # Seals to 1 box.  Fails against: negated flag silently dropped (would still produce
    # 1 box, but the box would allow carbon instead of forbidding it).
    # We verify via box_counts that the seal succeeds and produces exactly 1 box.
    q = QueryContainer()
    a = q.add_atom()
    q.atom_primitive(a, 'element', 6, negated=True)
    assert q.box_counts() == [1]


def test_edit_discards_restores_next_id():
    # Pins self._q._next_id = self._saved_next_id in _QueryEditScope.__exit__.
    # Without it: after discarding an edit that called add_atom, _next_id stays advanced.
    # Pre-edit, a = add_atom() consumed id 1 (_next_id becomes 2).  The edit scope captures
    # _saved_next_id = 2.  Inside, add_atom() consumes id 2 (_next_id becomes 3).  After
    # rollback, _next_id must be restored to 2 so the next add_atom() returns 2, not 3.
    q = QueryContainer()
    a = q.add_atom()
    q.atom_primitive(a, 'element', 6)
    with pytest.raises(RuntimeError):
        with q.edit():
            q.add_atom()
            raise RuntimeError()
    assert q.add_atom() == 2   # wrong impl (no restore): would return 3


def test_edit_discards_invalidates_stale_sealed_arena():
    # Pins self._q._invalidate() in _QueryEditScope.__exit__.
    # Without it: _query still holds the two-atom arena after rollback.  Any code that
    # calls box_counts() post-rollback would silently match against the discarded atom.
    q = QueryContainer()
    a = q.add_atom()
    q.atom_primitive(a, 'element', 6)
    q.box_counts()   # seals one-atom query (_seal_generation becomes 1)
    with pytest.raises(RuntimeError):
        with q.edit():
            b = q.add_atom()
            q.atom_primitive(b, 'element', 8)
            q.box_counts()   # seals two-atom arena, _query now holds it
            raise RuntimeError()
    # wrong impl (no _invalidate): _query still holds two-atom arena → [1, 1]
    assert q.box_counts() == [1]


def test_edit_discard_makes_empty_query_raise_on_next_seal():
    # Companion to test_edit_discards_restores_next_id for the fully-empty case.
    # A fresh container discards its only add_atom: _next_id must be restored to 1,
    # so sealed() raises ValueError (the _next_id == 1 guard fires).
    # Without _next_id restore: _next_id stays 2, query_seal gets a zero-op journal,
    # and atom_count_sealed() returns 0 through an except 0 signature with no exception.
    q = QueryContainer()
    with pytest.raises(RuntimeError):
        with q.edit():
            q.add_atom()
            raise RuntimeError()
    with pytest.raises(ValueError, match='empty query'):
        q.atom_count_sealed()


def test_nested_edit_outer_rollback_discards_both_levels():
    # _scope_depth must ensure only the outermost scope rolls back.
    # Without _scope_depth the inner scope would roll back on exception and
    # the outer scope would see scope_depth=0 with no exception, committing nothing.
    q = QueryContainer()
    with pytest.raises(RuntimeError):
        with q.edit():            # outer: saves len=0, next_id=1, depth→1
            a = q.add_atom()
            q.atom_primitive(a, 'element', 6)
            with q.edit():        # inner: saves len=2, next_id=2, depth→2
                b = q.add_atom()
                q.atom_primitive(b, 'element', 8)
            # inner exits cleanly: depth→1, no action
            raise RuntimeError()
        # outer exits with exception at depth=0: rollback to len=0, next_id=1
    assert q.atom_count == 0
    assert q.add_atom() == 1    # _next_id restored to 1, not left at 3


def test_nested_edit_inner_exception_caught_between_scopes_keeps_ops():
    # An exception caught *between* the inner and outer scope leaves the inner ops in
    # the journal.  The outer scope exits cleanly and commits.
    q = QueryContainer()
    with q.edit():                # outer: depth→1
        a = q.add_atom()
        q.atom_primitive(a, 'element', 6)
        try:
            with q.edit():        # inner: depth→2
                b = q.add_atom()
                q.atom_primitive(b, 'element', 8)
                raise RuntimeError()
            # inner exits with exception: depth→1, scope_depth != 0 → no rollback
        except RuntimeError:
            pass
        # b's ops are still in the journal; outer scope exits cleanly
    # outer exits cleanly: depth→0, no rollback → both atoms committed
    assert q.atom_count == 2
