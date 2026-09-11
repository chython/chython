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
"""H_UNKNOWN: the implicit-hydrogen count that a record does not state.

The state exists because real files omit the number and no rule recovers it: an MDL record may carry
an aromatic atom in a ring that will not kekulise, and a SMILES atom may be an element and charge no
valence rule covers.  The alternative is to turn both into a zero, which is a different molecule -- an
atom with no hydrogens is methane's neighbour, an atom with no RECORDED hydrogens is an unanswered
question.

NO POPULATION COUNT IS QUOTED HERE, deliberately.  "How many atoms hit a table miss" and "how many
atoms end up with no answer" measure different events and yield different numbers on one corpus, so a
frequency argues for nothing: a representation must hold "no answer" whether or not today's corpus
reaches it, and the tests below assert behaviour and never a frequency.

What these tests pin, in the order the value travels:

  storage      the sentinel lives in the implicit nibble and survives to_bytes/from_bytes, copy,
               substructure and union without becoming a count
  surface      writers say H_UNKNOWN, readers answer None, and None on the writer side still means
               "leave it at zero" -- the three-way distinction in add_atom's docstring
  arithmetic   float() leaves the mass light rather than adding fifteen hydrogens, and
               unknown_h_count is how a caller learns the mass is a lower bound
  screening    both hydrogen spans of feature word III go FULL, which in this encoding refuses
               every h demand rather than satisfying all of them -- positive and negated alike
  perception   no stereo unit is perceived on such an atom WHERE THE MISSING NUMBER COULD HAVE
               CHANGED THE FRAME -- and one IS perceived where the named directions already fill it,
               because there no value the count could take would be admitted anyway
"""
import pytest
from chython.core import MoleculeContainer, QueryContainer, H_IMPLICIT_MAX, H_UNKNOWN
from chython.core._core import kekule


# The two feature-word III spans this file asserts on, from SPAN_MASK in _query_boxes.pxi.
IMPLICIT_H_SPAN = 0x00000000003E0000   # word 2 bits 17-21
EXPLICIT_H_SPAN = 0x0000000007C00000   # word 2 bits 22-26
TOTAL_H_SPAN = 0x00000001F8000000      # word 2 bits 27-32


def build(atoms, bonds=()):
    """atoms is a list of (element, implicit_h); implicit_h is passed through verbatim, so None
    means "do not write a count" exactly as it does on `add_atom`."""
    m = MoleculeContainer()
    ids = []
    with m.edit():
        for element, h in atoms:
            ids.append(m.add_atom(element, implicit_h=h))
        for a, b, *rest in bonds:
            m.add_bond(ids[a], ids[b], rest[0] if rest else 1)
    return m, ids


def h_query(value, negated=False, element=6):
    q = QueryContainer()
    with q.edit():
        s = q.add_atom()
        q.atom_primitive(s, 'element', element, False)
        q.atom_operator(s, 'and_low')
        q.atom_primitive(s, 'implicit_h', value, negated)
    return q


def matched_atoms(q, m):
    return sorted(next(iter(mapping.values())) for mapping in q.get_mapping(m))


# ---------------------------------------------------------------------------------------------
# the constant


def test_unknown_h_is_the_nibble_value_and_not_a_flag():
    """15, because no atom carries fifteen implicit hydrogens -- the sentinel costs no storage.

    Asserted rather than assumed because the number is written into files by the MDL and SMILES
    readers and read back by this package; a silent change to 14 would turn a real count into a
    sentinel in every record already on disk.
    """
    assert H_UNKNOWN == 15


def test_the_bound_is_exported_beside_the_sentinel_and_is_one_below_it():
    """A reader validating a count it parsed needs the number 14, so it is on the surface too.

    Together, because they are one fact: the count stops where the sentinel starts.  Exported
    because the alternative is what actually happened -- the CTfile reader wrote its own literal
    `H_MAX = 15` from the nibble's width and admitted the sentinel as a count on three write paths.
    A bound restated as a literal is a bound that drifts, and this one drifting turns a stated
    count into "nobody knows".
    """
    assert H_IMPLICIT_MAX == 14
    assert H_UNKNOWN == H_IMPLICIT_MAX + 1


def test_fourteen_is_still_a_count():
    """The boundary.  14 is the largest count the nibble can hold beside the sentinel, and it is
    a count -- absurd chemistry, but the arena stores records, not judgements."""
    m, (a,) = build([('C', 14)])
    assert m.implicit_h_of(a) == 14
    assert m.unknown_h_count == 0


def test_no_bound_on_an_implicit_count_admits_the_sentinel():
    """A validator that stops at the nibble's WIDTH accepts 15 as a count, and that is the one
    thing the third state cannot survive.

    The nibble is four bits wide and the explicit half really does use all sixteen values, so
    "0..15" is a true statement about the layout and a false one about an implicit count.  Every
    place that bounds a caller-stated implicit count therefore stops at 14, and the two entry
    points a caller can reach are checked here together because they were written apart:
    `set_hydrogens`, where 15 IS the sentinel and 16 is out of range, and `kekule`'s `stated_h`,
    where 15 is out of range outright -- a caller who means "the input said nothing" has
    AROM_H_UNSTATED (-1) there and does not need the arena's spelling.
    """
    m, (a,) = build([('C', None)])
    with pytest.raises(ValueError, match='0..14'):
        m.set_hydrogens(a, 16)
    m.set_hydrogens(a, H_UNKNOWN)                  # 15 is accepted, AS THE SENTINEL
    assert m.implicit_h_of(a) is None

    ring = MoleculeContainer()
    ids = [ring.add_atom('C') for _ in range(6)]
    aromatic = [(ids[i], ids[(i + 1) % 6]) for i in range(6)]
    for x, y in aromatic:
        ring.add_bond(x, y, 1)
    with pytest.raises(ValueError, match=r'stated_h\[%d\] = 15 is outside 0\.\.14' % ids[0]):
        kekule(ring, aromatic, {ids[0]: 15})
    assert kekule(ring, aromatic, {ids[0]: 1}).changed, 'a real count still kekulises'


# ---------------------------------------------------------------------------------------------
# the surface: three kinds of value in, TWO kinds out, and `None` is the same statement as the
# sentinel rather than a third one


def test_the_three_writer_values_are_only_TWO_different_statements():
    """`3` states a count; `H_UNKNOWN` and `None` both state that nobody counted.

    THERE IS NO THIRD STATEMENT, AND ITS ABSENCE IS THE POINT.  `implicit_h=None` does not store a
    zero: a caller who omits an argument has not made a statement, and charging them with "zero
    hydrogens" is the same "unset is indistinguishable from stated zero" trap as `H_IMPLICIT_MAX` and
    `SGROUP_NO_INDEX`.  A caller who means zero has a way to say it, and it is `implicit_h=0` -- which
    is why 0 is asserted here alongside the other two.
    """
    m, (stated, zero, unknown, omitted) = build([('C', 3), ('C', 0), ('C', H_UNKNOWN), ('O', None)])
    assert m.implicit_h_of(stated) == 3
    assert m.implicit_h_of(zero) == 0, 'a stated zero is still a count and still reads back as one'
    assert m.implicit_h_of(unknown) is None
    assert m.implicit_h_of(omitted) is None, 'omission is not a count'
    # And the two spellings of "nobody counted" are indistinguishable ON PURPOSE -- one stored value,
    # so no consumer has to handle two flavours of unknown.
    assert m.implicit_h_of(omitted) == m.implicit_h_of(unknown)


def test_readers_answer_none_and_never_the_number():
    """None on every read path, so an untaught caller gets a TypeError on the first arithmetic
    rather than a plausible fifteen that propagates into a formula."""
    m, (a,) = build([('C', H_UNKNOWN)])
    assert m.implicit_h_of(a) is None
    assert m.total_h_of(a) is None
    assert m.atom(a).implicit_h is None
    assert m.atom(a).total_h is None


def test_explicit_hydrogens_stay_a_number():
    """The asymmetry the design rests on: an explicit hydrogen is an atom someone drew, so its
    count is known even when the implicit one is not.  `explicit_h_of` is the way to ask for the
    half of the total that IS known."""
    m, (c, h) = build([('C', H_UNKNOWN), ('H', 0)], [(0, 1)])
    assert m.explicit_h_of(c) == 1
    assert m.implicit_h_of(c) is None
    assert m.total_h_of(c) is None       # a sum with an unknown term


def test_set_hydrogens_writes_the_sentinel_and_can_write_over_it():
    """Both directions.  Recording "unknown" is a write, and so is later learning the answer --
    a reader that finds the count in a later record block must be able to replace the sentinel."""
    m, (a,) = build([('C', 3)])
    m.set_hydrogens(a, H_UNKNOWN)
    assert m.implicit_h_of(a) is None
    m.set_hydrogens(a, 2)
    assert m.implicit_h_of(a) == 2


def test_atom_view_setter_writes_the_sentinel():
    m, (a,) = build([('C', 3)])
    m.atom(a).implicit_h = H_UNKNOWN
    assert m.atom(a).implicit_h is None


@pytest.mark.parametrize('bad', [-1, 16, 20, 255])
def test_counts_stop_at_fourteen(bad):
    """A count above 14 is refused, so the sentinel cannot be reached by an off-by-one in a
    parser's arithmetic -- only by naming it."""
    m = MoleculeContainer()
    with pytest.raises(ValueError, match='implicit_h must be in 0..14'):
        m.add_atom('C', implicit_h=bad)


@pytest.mark.parametrize('bad', [-1, 16, 20])
def test_set_hydrogens_has_the_same_domain(bad):
    m, (a,) = build([('C', 0)])
    with pytest.raises(ValueError, match='implicit_h must be in 0..14'):
        m.set_hydrogens(a, bad)


def test_unknown_h_count_is_zero_on_a_complete_record():
    """The property is the one test a caller needs before trusting a derived number, so it has to
    be quiet on the normal case."""
    m, _ = build([('C', 3), ('C', 2), ('O', 1)], [(0, 1), (1, 2)])
    assert m.unknown_h_count == 0


def test_unknown_h_count_counts_atoms_and_not_hydrogens():
    m, _ = build([('C', H_UNKNOWN), ('C', H_UNKNOWN), ('O', 1)], [(0, 1), (1, 2)])
    assert m.unknown_h_count == 2


# ---------------------------------------------------------------------------------------------
# storage: the sentinel survives every path that copies an atom record


def test_bytes_round_trip_preserves_the_sentinel():
    """The requirement the MDL reader asked for in exactly these words: a value `implicit_h_of`
    can return that survives pack/unpack.  It travels as the nibble it is -- no load-time
    validation narrows it, and `rebuild_derived` re-derives the feature words from it."""
    m, ids = build([('C', 3), ('C', H_UNKNOWN), ('O', 0)], [(0, 1), (1, 2)])
    r = MoleculeContainer.from_bytes(m.to_bytes())
    assert [r.implicit_h_of(x) for x in r.atom_numbers] == [3, None, 0]
    assert r.unknown_h_count == 1


def test_copy_preserves_the_sentinel():
    m, (a, b) = build([('C', 3), ('C', H_UNKNOWN)], [(0, 1)])
    c = m.copy()
    assert c.implicit_h_of(a) == 3
    assert c.implicit_h_of(b) is None


def test_substructure_does_not_invent_a_count():
    """An atom nobody had a count for must not acquire one by being cut out of a bigger molecule.
    The copy path passes the RAW nibble for this reason."""
    m, (a, b, c) = build([('C', 3), ('C', H_UNKNOWN), ('O', 0)], [(0, 1), (1, 2)])
    sub = m.substructure([a, b])
    assert sorted(x for x in (sub.implicit_h_of(i) for i in sub.atom_numbers) if x is not None) == [3]
    assert sub.unknown_h_count == 1


def test_union_preserves_both_sides():
    m, _ = build([('C', H_UNKNOWN)])
    u = m.union(m.copy())
    assert u.unknown_h_count == 2


# ---------------------------------------------------------------------------------------------
# arithmetic


def test_mass_is_light_and_says_nothing():
    """chython 2's behaviour, kept deliberately: `float()` is not the place to raise, and there is
    no better number than "no hydrogens here".  The mass of the unknown-H molecule is therefore a
    LOWER BOUND, and `unknown_h_count` is how a caller finds that out."""
    stated, _ = build([('C', 4)])
    unknown, _ = build([('C', H_UNKNOWN)])
    assert float(stated) == pytest.approx(16.043, abs=1e-3)
    assert float(unknown) == pytest.approx(12.011, abs=1e-3)
    assert float(unknown) < float(stated)
    assert unknown.unknown_h_count == 1


def test_mass_does_not_add_fifteen_hydrogens():
    """The failure this whole commit exists to prevent, stated as a number: reading the sentinel as
    a count puts 15 protons on the atom and reports methane as 27 daltons."""
    unknown, _ = build([('C', H_UNKNOWN)])
    assert float(unknown) < 13.0


# ---------------------------------------------------------------------------------------------
# screening


def test_both_hydrogen_spans_go_full():
    """FULL SPANS MEAN "NO h DEMAND MATCHES", not "all of them do".  The kernel's atom test is
    `word & box.neg` and a box states `h2` by forbidding the rest of its span, so an atom carrying
    the whole span is refused by any box that touched it."""
    m, (stated, unknown) = build([('C', 3), ('C', H_UNKNOWN)], [(0, 1)])
    fs = m.features_of(stated)
    fu = m.features_of(unknown)
    assert fs[2] & IMPLICIT_H_SPAN != IMPLICIT_H_SPAN     # one-hot: exactly one bit
    assert fu[2] & IMPLICIT_H_SPAN == IMPLICIT_H_SPAN
    assert fu[2] & TOTAL_H_SPAN == TOTAL_H_SPAN


def test_the_explicit_span_stays_one_hot():
    """The known half of the record keeps screening exactly.  An unknown implicit count does not
    make the drawn hydrogens unknown."""
    m, (c, _) = build([('C', H_UNKNOWN), ('H', 0)], [(0, 1)])
    span = m.features_of(c)[2] & EXPLICIT_H_SPAN
    assert span and span & (span - 1) == 0               # exactly one bit set


def test_no_h_demand_matches_in_either_direction():
    """DIVERGES FROM chython 2 ON THE NEGATED FORM, on purpose.  There `h` compared against
    `Element.implicit_hydrogens`, so `None != 2` was True and `[C;!h2]` matched an atom with no
    count at all -- a match granted by the absence of data.  Here it does not match.
    """
    m, (stated, unknown) = build([('C', 3), ('C', H_UNKNOWN)], [(0, 1)])
    assert matched_atoms(h_query(3), m) == [stated]
    assert matched_atoms(h_query(3, negated=True), m) == []
    assert matched_atoms(h_query(0), m) == []
    assert matched_atoms(h_query(0, negated=True), m) == [stated]


def test_round_trip_re_derives_the_same_words():
    """The feature words are derived, so from_bytes recomputes them -- and must reach the same
    full spans, or a query would answer differently before and after a save."""
    m, _ = build([('C', 3), ('C', H_UNKNOWN), ('O', 0)], [(0, 1), (1, 2)])
    r = MoleculeContainer.from_bytes(m.to_bytes())
    assert ([r.features_of(x) for x in r.atom_numbers]
            == [m.features_of(x) for x in m.atom_numbers])


# ---------------------------------------------------------------------------------------------
# perception


def test_no_tetrahedral_unit_where_the_missing_count_could_have_mattered():
    """CHFClBr is a stereocentre; the same skeleton with no recorded hydrogen count is not one that
    can be configured.  Three heavy neighbours plus a hydrogen is a centre, three heavy neighbours
    and nothing is not, and the missing number is precisely which."""
    stated, _ = build([('C', 1), ('F', 0), ('Cl', 0), ('Br', 0)], [(0, 1), (0, 2), (0, 3)])
    unknown, _ = build([('C', H_UNKNOWN), ('F', 0), ('Cl', 0), ('Br', 0)],
                       [(0, 1), (0, 2), (0, 3)])
    assert len(stated.stereo_units()) == 1
    assert unknown.stereo_units() == []


def test_a_fully_substituted_centre_is_perceived_despite_the_sentinel():
    """CFClBrI: four heavy neighbours leave no room for a hydrogen, so the missing count cannot
    change the frame and the unit must be built.

    A refusal here does not surface as a missing unit -- it surfaces as a STATED PARITY THAT CANNOT BE
    WRITTEN.  A query format hands over exactly this shape: a fully substituted centre with a
    configuration and no hydrogen count anywhere.  Refuse perception and the writer holds a parity with
    no unit under it, so the caller gets a SMILES with no `@` in it and no explanation.  The sentinel
    reads as zero HERE ONLY, and only because every other value it could have taken is refused by the
    four-direction test regardless.

    `implicit_h_of` is asserted alongside on purpose: the unit exists and the count is still not
    known.  Perceiving the frame may not be allowed to turn into a derivation of the number.
    """
    m, ids = build([('C', H_UNKNOWN), ('F', 0), ('Cl', 0), ('Br', 0), ('I', 0)],
                   [(0, 1), (0, 2), (0, 3), (0, 4)])
    units = m.stereo_units()
    assert len(units) == 1
    assert units[0]['kind'] == 0                        # SU_TETRA
    assert units[0]['n_refs'] == 4
    assert units[0]['unnamed_mask'] == 0, 'no unnamed direction may be invented for the sentinel'
    assert m.implicit_h_of(ids[0]) is None, 'the frame is known; the count is still not'

    m.set_parity(ids[0], 1)
    assert m.parity_of(ids[0]) == 1
    assert MoleculeContainer.from_bytes(m.to_bytes()).parity_of(ids[0]) == 1


def test_a_fully_substituted_cumulene_terminal_is_perceived_despite_the_sentinel():
    """The same argument two positions over: a terminal has two in-plane directions, and when both
    are named the missing count cannot decide anything.  A tetrasubstituted allene with the sentinel
    on one terminal keeps its axis."""
    m, _ = build([('F', 0), ('Cl', 0), ('C', H_UNKNOWN), ('C', 0), ('C', 0), ('Br', 0), ('I', 0)],
                 [(0, 2), (1, 2), (2, 3, 2), (3, 4, 2), (4, 5), (4, 6)])
    units = m.stereo_units()
    assert [u['kind'] for u in units] == [2]            # SU_ALLENE
    assert units[0]['unnamed_mask'] == 0


def test_a_fully_substituted_nitrogen_is_not_called_protic_by_the_sentinel():
    """The second half of the same defect, one stage down: `_anchor_is_protic` asks whether a group
    15/16 anchor has a hydrogen to invert through, and the sentinel is not a count there either.

    Reading it raw makes `implicit_h != 0` true and the centre is denied as protic -- so narrowing
    perception alone would have produced a UNIT THAT IS NEVER STEREOGENIC, which is a subtler wrong
    answer than no unit at all.  Four named directions on the nitrogen mean no unnamed slot, and an
    implicit hydrogen is by definition unnamed, so the frame itself settles the question without
    reading the nibble."""
    m, ids = build([('N', H_UNKNOWN), ('C', 3), ('F', 0), ('Cl', 0), ('Br', 0)],
                   [(0, 1), (0, 2), (0, 3), (0, 4)])
    m.set_charge(ids[0], 1)
    anchor = [u for u in m.stereo_units() if u['anchor'] == 1]
    assert len(anchor) == 1
    assert anchor[0]['stereogenic'] is True, 'the sentinel was read as a hydrogen count'

    protic, _ = build([('N', 1), ('C', 3), ('F', 0), ('Cl', 0)], [(0, 1), (0, 2), (0, 3)])
    assert [u['stereogenic'] for u in protic.stereo_units() if u['anchor'] == 1] == [False]


def test_a_sulfur_lone_pair_still_leaves_room_and_is_still_refused():
    """The narrowing is by DIRECTION COUNT and not by heavy-atom degree, which matters exactly once:
    a sulfur's lone pair is a direction nothing names.  Four heavy neighbours on sulfur plus the pair
    is five, so this record is refused -- and it is refused by the four-direction test rather than by
    the sentinel, which is why treating the sentinel as zero here is safe."""
    m, _ = build([('S', H_UNKNOWN), ('F', 0), ('Cl', 0), ('Br', 0), ('I', 0)],
                 [(0, 1), (0, 2), (0, 3), (0, 4)])
    assert m.stereo_units() == []


def test_no_cis_trans_unit_when_a_terminal_h_is_unknown():
    """But-2-ene has a cis/trans unit; the same skeleton with an unrecorded count on one
    double-bond carbon does not.  One hydrogen and one methyl is a genuine E/Z pair, two hydrogens
    is `=CH2` and has no isomer, and the count is what would say which."""
    stated, _ = build([('C', 3), ('C', 1), ('C', 1), ('C', 3)], [(0, 1), (1, 2, 2), (2, 3)])
    unknown, _ = build([('C', 3), ('C', H_UNKNOWN), ('C', 1), ('C', 3)],
                       [(0, 1), (1, 2, 2), (2, 3)])
    kinds = [u['kind'] for u in stated.stereo_units()]
    assert 1 in kinds                                   # SU_CIS_TRANS present
    assert 1 not in [u['kind'] for u in unknown.stereo_units()]


def test_the_methyl_units_are_untouched_by_the_cut():
    """The refusal is per atom and not per molecule: but-2-ene's two methyl carbons keep their
    (non-stereogenic) tetrahedral units when the middle carbon's count goes missing."""
    unknown, _ = build([('C', 3), ('C', H_UNKNOWN), ('C', 1), ('C', 3)],
                       [(0, 1), (1, 2, 2), (2, 3)])
    assert [u['kind'] for u in unknown.stereo_units()] == [0, 0]


def test_symmetry_keeps_unknown_apart_from_a_count():
    """An automorphism may not map an atom whose hydrogens were recorded onto one whose were not.
    Propane's two methyls are interchangeable; make one of them unrecorded and they are not, so the
    two molecules cannot have the same canonical form."""
    both, _ = build([('C', 3), ('C', 2), ('C', 3)], [(0, 1), (1, 2)])
    one, _ = build([('C', 3), ('C', 2), ('C', H_UNKNOWN)], [(0, 1), (1, 2)])
    assert both != one
    assert both.canonical_bytes != one.canonical_bytes


def test_two_unknowns_are_still_symmetric():
    """The other direction: unknown maps onto unknown, so a molecule with the count missing from
    both ends keeps its symmetry and its canonical form does not depend on which end is which."""
    a, _ = build([('C', H_UNKNOWN), ('C', 2), ('C', H_UNKNOWN)], [(0, 1), (1, 2)])
    b, _ = build([('C', H_UNKNOWN), ('C', 2), ('C', H_UNKNOWN)], [(2, 1), (1, 0)])
    assert a == b
