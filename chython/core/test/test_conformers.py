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
"""`SEG_CONFORMERS` -- the third coordinate, at one model.

Reads the core directly and never `chython`, per the isolation ratchet.
"""
from struct import pack, pack_into, unpack_from

import pytest

from chython.core import MoleculeContainer
from chython.core import _core


# One water, placed. Public, tiny, and its three atoms are enough to tell a dropped column from a
# shifted one -- which two atoms would not be.
WATER = ((8, 0.0, 0.0, 0.0), (1, 0.9572, 0.0, 0.0), (1, -0.2400, 0.9266, 0.0))


def _water(with_xyz=True, with_xy=False):
    m = MoleculeContainer()
    with m.edit():
        o = m.add_atom(8)
        h1 = m.add_atom(1)
        h2 = m.add_atom(1)
        m.add_bond(o, h1, 1)
        m.add_bond(o, h2, 1)
    with m.edit():
        for n, (_, x, y, z) in zip((1, 2, 3), WATER):
            if with_xyz:
                m.set_xyz(n, x, y, z)
            if with_xy:
                m.set_xy(n, x, y)
    return m


def _segment(data, seg):
    """`(offset, length)` of table entry `seg`, or None when the buffer's table stops short."""
    seg_count = unpack_from('<H', data, 20)[0]
    if seg >= seg_count:
        return None
    return unpack_from('<II', data, 24 + 8 * seg)


# ---------------------------------------------------------------------------------------------
# the round trip, both ways round
# ---------------------------------------------------------------------------------------------

def test_a_molecule_with_no_conformer_says_so_and_carries_no_segment():
    """M == 0 is the common case and must cost nothing -- not a zeroed segment, not a header entry.

    The point of a separate segment rather than a wider `xy_t` is exactly this: a molecule with no
    third coordinate has the bytes it always had. `_segment` returning None here is that claim in the
    header itself, and not merely a length of zero.
    """
    m = _water(with_xyz=False)
    assert m.has_3d is False
    assert m.xyz_of(1) is None
    assert _segment(m.to_bytes(), _core.SEG_CONFORMERS) is None


def test_coordinates_survive_to_bytes_and_back_exactly():
    """Exactly, and not approximately: the wire format is fixed point on the same grid as the input.

    `XY_SCALE` is 10000, so a coordinate written with four decimals round-trips bit for bit. The
    equality here is therefore a real assertion and not a tolerance in disguise -- a float comparison
    that needed a tolerance would mean the scale had changed under it.
    """
    m = _water()
    assert m.has_3d is True
    back = MoleculeContainer.from_bytes(m.to_bytes())
    assert back.has_3d is True
    for n, (_, x, y, z) in zip((1, 2, 3), WATER):
        assert m.xyz_of(n) == (x, y, z)
        assert back.xyz_of(n) == (x, y, z)
    assert back.to_bytes() == m.to_bytes()


def test_the_conformer_segment_is_the_whole_difference_in_the_bytes():
    """A 3D water's buffer differs from a flat one's by one segment and its table entry, and no more.

    This is byte-neutrality stated from the consumer's side. Every other segment's LENGTH is unchanged
    and the atom records are identical, which is what makes 3.1's multi-model growth purely additive:
    nothing the buffer already held has moved.
    """
    flat = _water(with_xyz=False).to_bytes()
    solid = _water().to_bytes()
    conf = _segment(solid, _core.SEG_CONFORMERS)
    assert conf is not None
    # 3 atoms * 12 bytes + an 8-byte header + one 4-byte record = 48, already 8-aligned.
    assert conf[1] == 48
    assert len(solid) - len(flat) == 48 + 8 * (_core.SEG_CONFORMERS + 1
                                               - unpack_from('<H', flat, 20)[0])
    # LENGTHS, and an absent entry counts as zero. Spending id 8 pushes `seg_count` to 9, so entries
    # this molecule does not use EXIST as `(0, 0)` empties rather than the table stopping short of
    # them -- documented behaviour ("interior empties still occupy their slot"), not a segment
    # appearing. What must not change is a byte of PAYLOAD, so lengths are the comparison, and the
    # table's own growth is accounted for in the size assertion above.
    for seg in range(_core.SEG_CONFORMERS):
        a, b = _segment(flat, seg), _segment(solid, seg)
        assert (0 if a is None else a[1]) == (0 if b is None else b[1]), \
            f'segment {seg} changed length'


def test_the_record_starts_life_saying_the_file_named_no_model():
    """`ext_index` is stamped with `CONF_NO_INDEX` and not left as the memset's zero.

    Zero is a REAL model number -- a PDB `MODEL 0` is representable and nothing forbids it -- so a
    fresh conformer that reported 0 would be claiming a file said something it never said. The
    sentinel is the top u32 value for that reason, and this test reads it out of the serialised bytes
    rather than through an accessor, because there is no accessor for it in the slice.
    """
    data = _water().to_bytes()
    off, _ = _segment(data, _core.SEG_CONFORMERS)
    count, reserved = unpack_from('<II', data, off)
    assert count == 1
    assert reserved == 0
    ext_index, = unpack_from('<I', data, off + 8)
    assert ext_index == _core.CONF_NO_INDEX
    assert _core.CONFORMER_RECORD_SIZE == 4, 'the record is that one field and nothing after it'


# ---------------------------------------------------------------------------------------------
# in the bytes, out of the identity
# ---------------------------------------------------------------------------------------------

def test_geometry_is_in_the_bytes_and_out_of_the_identity():
    """`to_bytes()` sees the conformer; `==`, `hash` and the canonical form do not.

    The rule is RULES.md 1.5's for a CIP descriptor, applied to geometry: the canonical form answers
    "which molecule is this", and two conformers of one molecule are one molecule. `to_bytes()`
    answers "which stored record is this", and they are two records. Getting this backwards either way
    is expensive -- an identity that saw geometry would make every conformer a new compound, and a
    `to_bytes()` that did not would silently discard one of the two on a round trip through a store.
    """
    solid, flat = _water(), _water(with_xyz=False)
    assert solid == flat
    assert hash(solid) == hash(flat)
    assert solid.to_bytes() != flat.to_bytes()
    # and the canonical numbering is untouched, which is the thing `==` is actually built on
    assert solid.atoms_order == flat.atoms_order


def test_two_conformers_of_one_molecule_are_equal_and_differently_stored():
    m1, m2 = _water(), _water()
    with m2.edit():
        m2.set_xyz(1, 5.0, 5.0, 5.0)
    assert m1 == m2 and hash(m1) == hash(m2)
    assert m1.to_bytes() != m2.to_bytes()


# ---------------------------------------------------------------------------------------------
# 2D and 3D are independent
# ---------------------------------------------------------------------------------------------

def test_the_two_coordinate_segments_do_not_touch_each_other():
    """`set_xyz` leaves `SEG_XY` byte-identical and `set_xy` leaves the conformer byte-identical.

    Design D2, and the reason it matters is `clean2d()`: a 3D file fills both segments, and a layout
    engine that rewrote the depiction would destroy the geometry if the two shared storage. Asserted
    as a byte comparison of the segments themselves, because a comparison through the accessors would
    pass even if one segment were a view of the other.
    """
    both = _water(with_xyz=True, with_xy=True)
    data = both.to_bytes()
    xy = _segment(data, _core.SEG_XY)
    conf = _segment(data, _core.SEG_CONFORMERS)
    assert xy is not None and conf is not None
    xy_bytes = data[xy[0]:xy[0] + xy[1]]
    conf_bytes = data[conf[0]:conf[0] + conf[1]]

    with both.edit():
        both.set_xyz(1, -7.5, -7.5, -7.5)
    after = both.to_bytes()
    assert after[xy[0]:xy[0] + xy[1]] == xy_bytes, 'set_xyz moved the depiction'
    assert after[conf[0]:conf[0] + conf[1]] != conf_bytes

    with both.edit():
        both.set_xy(1, 3.0, 4.0)
    after2 = both.to_bytes()
    assert after2[conf[0]:conf[0] + conf[1]] == after[conf[0]:conf[0] + conf[1]], \
        'set_xy moved the geometry'
    assert both.xy_of(1) == (3.0, 4.0)
    assert both.xyz_of(1) == (-7.5, -7.5, -7.5)


def test_a_three_dimensional_molecule_need_not_have_a_depiction():
    m = _water()
    assert m.has_3d is True
    assert m.has_coordinates is False
    assert m.xy_of(1) is None


# ---------------------------------------------------------------------------------------------
# edits
# ---------------------------------------------------------------------------------------------

def test_deleting_an_atom_drops_that_atom_s_column_and_no_other():
    """The one property a per-atom segment can get wrong, and the way it gets it wrong is a shift.

    Geometry follows the atom's slot through the same `newidx` that moves the atom record, so a
    surviving atom keeps ITS coordinate rather than its late neighbour's. Deleting the FIRST atom is
    the case that catches an off-by-one; deleting the last would not.
    """
    m = _water()
    with m.edit():
        m.delete_atom(1)
    assert m.has_3d is True
    assert list(m) == [2, 3]
    assert m.xyz_of(2) == WATER[1][1:]
    assert m.xyz_of(3) == WATER[2][1:]


def test_an_edit_that_touches_no_atom_carries_the_geometry_through_unchanged():
    """A bond edit rewrites no coordinate -- design D7's second half.

    Every atom stays live and `newidx[i] == i`, so the carry is a copy. Worth a test because the
    conformer travels through the same fresh-allocation path as everything else: an `_apply` that
    forgot to request the segment would silently return a flat molecule from a bond order change.
    """
    m = _water()
    before = [m.xyz_of(n) for n in (1, 2, 3)]
    with m.edit():
        m.delete_bond(1, 3)
    assert m.has_3d is True
    assert [m.xyz_of(n) for n in (1, 2, 3)] == before


def test_an_added_atom_sits_at_the_origin_and_this_is_the_slice_s_known_defect():
    """(0, 0, 0) and not None, because the segment is one dense column per model.

    There is no honest coordinate for an atom the caller placed by connectivity alone, and the per-atom
    validity bitmap that would let `xyz_of` answer None for it was declined for the slice (design D7).
    Pinned so the behaviour is a documented answer rather than a discovery: a caller building a
    molecule atom by atom on a 3D template gets an atom at the origin, and must place it.
    """
    m = _water()
    with m.edit():
        n = m.add_atom(1)
        m.add_bond(2, n, 1)
    assert m.xyz_of(n) == (0.0, 0.0, 0.0)
    assert m.xyz_of(1) == WATER[0][1:], 'the existing atoms kept their places'


def test_copy_and_substructure_and_union_all_carry_the_geometry():
    m = _water()
    assert m.copy().xyz_of(2) == WATER[1][1:]
    sub = m.substructure([1, 2])
    assert sub.has_3d is True and sub.xyz_of(2) == WATER[1][1:]
    # a union with a flat partner: the geometry it had survives, the partner has none to contribute
    other = MoleculeContainer()
    with other.edit():
        other.add_atom(7)
    joined = m.union(other)
    assert joined.xyz_of(1) == WATER[0][1:]
    assert joined.xyz_of(4) == (0.0, 0.0, 0.0)


def test_setting_an_sgroup_on_a_three_dimensional_molecule_does_not_flatten_it():
    """`structure_respan` has to carry the conformer, and nothing else in the suite would notice.

    An S-group write reallocates the buffer through a path of its own rather than through `_apply`'s
    fresh allocation, so a segment it does not know about is dropped in silence. Found by reading
    rather than by a failure, which is why the test exists at all.
    """
    m = _water()
    m.set_aliases({1: 'OH2'})
    assert m.aliases == {1: b'OH2'}
    assert m.has_3d is True
    assert m.xyz_of(1) == WATER[0][1:]
    assert m.xyz_of(3) == WATER[2][1:]


# ---------------------------------------------------------------------------------------------
# the domain, and the refusals
# ---------------------------------------------------------------------------------------------

def test_a_coordinate_outside_the_slot_s_range_is_refused():
    """The bound is the SLOT's and not any file's field -- RULES.md 1.4.

    `int32_t` at `XY_SCALE` reaches +/-214748.0, which is far wider than a V3000 line's ten-column
    `F10.4`. A format that wants the narrower one enforces it at its own boundary; the container
    refuses only what it cannot store.
    """
    m = _water(with_xyz=False)
    with m.edit():
        m.set_xyz(1, 214748.0, -214748.0, 0.0)      # the edges are IN range
        for bad in ((214749.0, 0.0, 0.0), (0.0, -214749.0, 0.0), (0.0, 0.0, 1e9)):
            with pytest.raises(ValueError, match='fixed point range'):
                m.set_xyz(1, *bad)
    assert m.xyz_of(1) == (214748.0, -214748.0, 0.0)


@pytest.mark.parametrize('field,value,message', [
    ('count', 0, 'conformer'),
    ('count', 0x10000, 'conformer'),
    ('reserved', 1, 'conformer'),
])
def test_a_malformed_conformer_header_is_refused(field, value, message):
    """The header is checked before a byte of the payload is trusted.

    `count == 0` is refused rather than read as "no models", because a segment that is PRESENT and
    empty is a writer that got confused -- an absent segment is how "no models" is spelled.
    `count > CONF_MAX_MODELS` cites the exported bound rather than the field's width: the field is
    `uint32_t` and would admit four billion models, whose payload overruns the buffer limit long
    before that.
    """
    data = bytearray(_water().to_bytes())
    off, _ = _segment(data, _core.SEG_CONFORMERS)
    pack_into('<I', data, off + (0 if field == 'count' else 4), value)
    with pytest.raises(ValueError, match=message):
        MoleculeContainer.from_bytes(bytes(data))


def test_a_declared_length_that_disagrees_with_the_model_count_is_refused():
    """An EQUALITY check, and that is the payoff for storing `count` in the payload.

    A segment whose length is merely a floor leaves a blind spot: a buffer declaring more bytes than
    it needs passes, and the surplus is unreachable but trusted. Here the length is a function of
    `count` and `atom_count`, both of which the buffer states, so the reader can demand the exact
    number and reject anything else. Shortening the DECLARED length is the test -- emptying the
    segment entirely would only exercise the `count == 0` arm above.
    """
    data = bytearray(_water().to_bytes())
    off, length = _segment(data, _core.SEG_CONFORMERS)
    pack_into('<I', data, 24 + 8 * _core.SEG_CONFORMERS + 4, length - 8)
    with pytest.raises(ValueError, match='conformer'):
        MoleculeContainer.from_bytes(bytes(data))


def test_the_journal_op_is_on_the_public_op_map():
    """`set_xyz` is in `JOURNAL_OPS`, which is how a test names an op without restating its number.

    `drop_conformer` is the highest op: `_apply` range-checks against `OP_HIGHEST`, and adding an op
    without raising that bound left the arm that would handle it raising `NotImplementedError`.
    A number pinned in two places is pinned in neither.
    """
    from chython.core._core import JOURNAL_OPS
    assert 'set_xyz' in JOURNAL_OPS
    assert 'add_conformer' in JOURNAL_OPS and 'drop_conformer' in JOURNAL_OPS
    assert JOURNAL_OPS['drop_conformer'] == max(JOURNAL_OPS.values())


def test_the_journal_record_reports_the_fourth_payload_word():
    m = _water(with_xyz=False)
    with m.edit():
        m.set_xyz(1, 1.0, -2.0, 3.5)
        op, a, b, v, w, model = m.journal_record(m.journal_length - 1)
        from chython.core._core import JOURNAL_OPS
        assert op == JOURNAL_OPS['set_xyz']
        assert a == 1
        # `b` is unsigned in the record and the coordinate is signed, so -2.0 arrives as its
        # two's complement -- the same thing `set_xy` already does with x.
        assert b == 10000
        assert v == -20000
        assert w == 35000
        assert model == 0


def test_the_journal_carries_the_model_in_its_padding():
    """The model index is 16 bits after `op`, so the record stays twenty bytes and no fifth payload
    word appears."""
    from chython.core._core import JOURNAL_OPS, journal_record_size
    m = _water()
    with m.edit():
        second = m.add_conformer(ext_index=9)
        m.set_xyz(1, 1., 2., 3., model=second)
        assert m.journal_record(0) == (JOURNAL_OPS['add_conformer'], 9, 0, 0, 0, 1)
        assert m.journal_record(1)[5] == 1
    assert journal_record_size() == 20


# ---------------------------------------------------------------------------------------------
# N models
# ---------------------------------------------------------------------------------------------

def test_a_second_model_is_added_and_read_back():
    m = _water()
    numbers = m.atom_numbers
    with m.edit():
        second = m.add_conformer(ext_index=7)
        for n, xyz in zip(numbers, [(1., 1., 1.), (2., 2., 2.), (3., 3., 3.)]):
            m.set_xyz(n, *xyz, model=second)
    assert second == 1
    assert len(m.conformers) == 2
    assert m.conformer(1).xyz_of(numbers[1]) == (2., 2., 2.)
    assert m.conformer(1).ext_index == 7
    assert m.conformer(0).xyz_of(numbers[1]) == m.xyz_of(numbers[1])
    assert m.xyz_of(numbers[1]) != (2., 2., 2.), 'xyz_of stays model 0'


def test_an_added_model_starts_at_the_origin():
    """One dense column per model with no per-atom validity bitmap, so a model nothing wrote reads
    the origin -- design D7 at N, and the same answer an atom added to a 3D molecule gets."""
    m = _water()
    with m.edit():
        m.add_conformer()
    assert m.conformer(1).coordinates == [(0., 0., 0.)] * 3
    assert m.conformer(1).ext_index is None


def test_a_dropped_model_compacts_and_its_neighbour_keeps_its_number():
    m = _water()
    numbers = m.atom_numbers
    with m.edit():
        for ext in (11, 22):
            model = m.add_conformer(ext_index=ext)
            for n in numbers:
                m.set_xyz(n, float(ext), 0., 0., model=model)
    with m.edit():
        m.drop_conformer(1)
    assert [c.ext_index for c in m.conformers] == [None, 22]
    assert m.conformer(1).xyz_of(numbers[0]) == (22., 0., 0.)


def test_dropping_every_model_leaves_no_segment():
    """The molecule reaches the state a never-placed one is in, down to the header.

    Not a present-and-empty entry: the seal spends segment id 8 only when there is a model, so the
    table stops short of it and the bytes are a flat molecule's -- which is what
    `test_a_molecule_with_no_conformer_says_so_and_carries_no_segment` asserts from the other side.
    """
    m = _water()
    with m.edit():
        m.drop_conformer(0)
    assert not m.has_3d
    assert m.conformers == ()
    assert m.xyz_of(m.atom_numbers[0]) is None
    assert _segment(m.to_bytes(), _core.SEG_CONFORMERS) is None


def test_a_session_either_adds_or_drops():
    """A drop shifts the models above it while an add counts from the unshifted source count, so a
    returned index would mean one thing before the seal and another after."""
    m = _water()
    with m.edit():
        m.add_conformer()
        with pytest.raises(ValueError, match='already added'):
            m.drop_conformer(0)
    m2 = _water()
    first = m2.atom_numbers[0]
    with m2.edit():
        m2.drop_conformer(0)
        with pytest.raises(ValueError, match='already dropped'):
            m2.add_conformer()
        with pytest.raises(ValueError, match='dropped a conformer'):
            m2.set_xyz(first, 1., 1., 1.)


def test_setting_a_model_that_does_not_exist_is_refused():
    m = _water()
    with pytest.raises(ValueError, match='does not exist'):
        m.set_xyz(m.atom_numbers[0], 1., 1., 1., model=3)


def test_the_first_model_of_a_flat_molecule_is_created_by_setting_it():
    m = _water(with_xyz=False, with_xy=True)
    assert not m.has_3d
    numbers = m.atom_numbers
    with m.edit():
        for n in numbers:
            m.set_xyz(n, 1., 2., 3.)
    assert len(m.conformers) == 1
    assert m.conformer(0).ext_index is None


def test_the_implicitly_created_first_model_counts_as_an_add():
    """A session that fills model 0 by setting it and then adds gets index 1, not 0 again.

    The implicit first model is journalled as an add rather than derived at the seal: `add_conformer`
    counts from this session's adds, so a derived one would hand back 0 and overwrite what the `set_xyz`
    calls just placed.  [mutant: derive the first model from a bare `set_xyz` at the seal]
    """
    m = _water(with_xyz=False, with_xy=True)
    numbers = m.atom_numbers
    with m.edit():
        for n in numbers:
            m.set_xyz(n, 1., 2., 3.)
        assert m.add_conformer(ext_index=7) == 1
        for n in numbers:
            m.set_xyz(n, 4., 5., 6., model=1)
    assert len(m.conformers) == 2
    assert [c.ext_index for c in m.conformers] == [None, 7]
    assert m.conformer(0).xyz_of(numbers[0]) == (1., 2., 3.)
    assert m.conformer(1).xyz_of(numbers[0]) == (4., 5., 6.)


def test_an_external_index_is_stored_verbatim_and_the_top_value_is_refused():
    m = _water()
    with m.edit():
        m.add_conformer(ext_index=0)
    assert m.conformer(1).ext_index == 0, 'zero is a MODEL number a file can state'
    with pytest.raises(ValueError, match='outside'):
        with m.edit():
            m.add_conformer(ext_index=_core.CONF_NO_INDEX)
    assert _core.CONF_EXT_INDEX_MAX == _core.CONF_NO_INDEX - 1


def test_n_models_survive_a_byte_round_trip():
    m = _water()
    numbers = m.atom_numbers
    with m.edit():
        for ext in (3, 4, 5):
            model = m.add_conformer(ext_index=ext)
            for i, n in enumerate(numbers):
                m.set_xyz(n, float(ext), float(i), -1., model=model)
    raw = m.to_bytes()
    back = MoleculeContainer.from_bytes(raw)
    assert back.to_bytes() == raw
    assert [c.ext_index for c in back.conformers] == [None, 3, 4, 5]
    assert back.conformer(2).coordinates == m.conformer(2).coordinates


# ---------------------------------------------------------------------------------------------
# the `Conformer` view
# ---------------------------------------------------------------------------------------------

def test_a_conformer_view_answers_for_its_own_model():
    m = _water()
    numbers = m.atom_numbers
    with m.edit():
        model = m.add_conformer(ext_index=4)
        for n in numbers:
            m.set_xyz(n, 5., 6., 7., model=model)
    first, second = m.conformers
    assert (first.index, second.index) == (0, 1)
    assert first.ext_index is None and second.ext_index == 4
    assert second.coordinates == [(5., 6., 7.)] * 3
    assert second.xyz_of(numbers[0]) == (5., 6., 7.)
    assert second == m.conformer(1) and second != first
    assert hash(second) == hash(m.conformer(1))
    with pytest.raises(IndexError):
        m.conformer(2)
    with pytest.raises(KeyError):
        first.xyz_of(999)


def test_a_conformer_view_goes_stale_with_the_molecule():
    """The same rule `Atom` and `Bond` hold: a borrowed handle raises rather than read an arena an
    edit has moved out from under it."""
    m = _water()
    view = m.conformer(0)
    with m.edit():
        m.add_conformer()
    with pytest.raises(RuntimeError, match='stale Conformer view'):
        view.coordinates
    m2 = _water()
    gone = m2.conformer(0)
    with m2.edit():
        m2.drop_conformer(0)
    with pytest.raises(RuntimeError, match='stale Conformer view'):
        gone.ext_index


def test_conformers_stay_out_of_the_molecule_identity():
    """Extends `test_geometry_is_in_the_bytes_and_out_of_the_identity` to N: the canonical form is
    about the chemistry, and a model count is not part of it."""
    flat = _water(with_xyz=False, with_xy=True)
    many = _water(with_xyz=False, with_xy=True)
    numbers = many.atom_numbers
    with many.edit():
        for _ in range(3):
            model = many.add_conformer()
            for n in numbers:
                many.set_xyz(n, 1., 1., 1., model=model)
    assert len(many.conformers) == 3
    assert flat == many
    assert hash(flat) == hash(many)


# ---------------------------------------------------------------------------------------------
# the carry-through paths
# ---------------------------------------------------------------------------------------------

def _two_model_water():
    m = _water()
    numbers = m.atom_numbers
    with m.edit():
        model = m.add_conformer(ext_index=2)
        for i, n in enumerate(numbers):
            m.set_xyz(n, float(i), 9., 9., model=model)
    return m


def test_copy_and_substructure_carry_every_model():
    m = _two_model_water()
    same = m.copy()
    assert [c.ext_index for c in same.conformers] == [None, 2]
    assert same.conformer(1).coordinates == m.conformer(1).coordinates
    keep = m.atom_numbers[:2]
    sub = m.substructure(keep)
    assert len(sub.conformers) == 2
    assert sub.conformer(1).ext_index == 2
    assert [sub.conformer(1).xyz_of(n) for n in sub.atom_numbers] == \
        [m.conformer(1).xyz_of(n) for n in keep]


def test_union_carries_the_wider_side_and_the_narrow_one_sits_at_the_origin():
    """A model index one side does not have leaves that side's atoms at the origin, which is the same
    answer a union with a molecule that has no geometry already gives the other way round."""
    many = _two_model_water()
    one = _water()
    joined = many | one
    assert len(joined.conformers) == 2
    assert [c.ext_index for c in joined.conformers] == [None, 2]
    other_atoms = list(joined.atom_numbers)[len(many.atom_numbers):]
    assert [joined.conformer(1).xyz_of(n) for n in other_atoms] == [(0., 0., 0.)] * 3
    flipped = one | many
    assert len(flipped.conformers) == 2
    assert flipped.conformer(1).ext_index == 2
    assert [flipped.conformer(1).xyz_of(n) for n in list(flipped.atom_numbers)[:3]] == \
        [(0., 0., 0.)] * 3


# ---------------------------------------------------------------------------------------------
# a version-5 buffer, whose conformer record is four words
# ---------------------------------------------------------------------------------------------

def _forge_v5(data, ext_index=_core.CONF_NO_INDEX):
    """This build's buffer rewritten as a version-5 one.

    The version byte, and the conformer table widened back to four words per model: `ext_index`,
    then the three words version 5 carried and this build does not model, given nonzero values so a
    migration that read them as coordinates would be caught rather than pass by luck.
    """
    off, length = _segment(data, _core.SEG_CONFORMERS)
    models, _ = unpack_from('<II', data, off)
    atoms = unpack_from('<I', data, 8)[0]
    xyz_at = off + 8 + 4 * models

    seg = bytearray(data[off:off + 8])
    for model in range(models):
        seg += pack('<I', ext_index) + pack('<Iq', 1, -1)
    seg += data[xyz_at:xyz_at + models * atoms * 12]
    while len(seg) % 8:
        seg += b'\0'

    grown = len(seg) - length
    out = bytearray(data[:off]) + seg + bytearray(data[off + length:])
    out[4:6] = pack('<H', _core.STRUCT_VERSION_V5)
    out[16:20] = pack('<I', unpack_from('<I', data, 16)[0] + grown)
    # Every segment starting after this one moves by the widening, and this one's length grows.
    for seg_id in range(unpack_from('<H', data, 20)[0]):
        entry = 24 + 8 * seg_id
        s_off, s_len = unpack_from('<II', data, entry)
        if not s_len:
            continue
        elif s_off == off:
            out[entry + 4:entry + 8] = pack('<I', s_len + grown)
        elif s_off > off:
            out[entry:entry + 4] = pack('<I', s_off + grown)
    return bytes(out)


def test_a_version_five_buffer_with_conformers_migrates():
    """The coordinates and `ext_index` survive; the three discarded words do not become geometry."""
    m = _water()
    back = MoleculeContainer.from_bytes(_forge_v5(m.to_bytes(), ext_index=7))

    assert back.has_3d
    for n in (1, 2, 3):
        assert back.xyz_of(n) == m.xyz_of(n)
    # THE RECORD READ OUT OF THE MIGRATED BYTES, there being no accessor for `ext_index` in the slice.
    data = back.to_bytes()
    off, length = _segment(data, _core.SEG_CONFORMERS)
    assert length == 48, "the migrated segment is laid out at this build's record width"
    assert unpack_from('<I', data, off + 8)[0] == 7, "the file's own model number travels"


def test_a_version_five_buffer_with_conformers_reports_the_narrowing():
    """`structure_from_bytes` has no log destination, so `from_bytes` is what says the words went."""
    back = MoleculeContainer.from_bytes(_forge_v5(_water().to_bytes()))

    assert [r for r in back.log if r.rule == 'container:conformer-narrowed'], \
        f'the narrowing must be on the molecule; log: {back.log}'


def test_a_version_five_buffer_with_no_conformers_migrates_by_the_version_byte():
    """The common case: no conformer segment, so nothing is re-laid and the layout is untouched."""
    m = _water(with_xyz=False, with_xy=True)
    raw = bytearray(m.to_bytes())
    raw[4:6] = pack('<H', _core.STRUCT_VERSION_V5)
    back = MoleculeContainer.from_bytes(bytes(raw))

    assert not back.has_3d
    for n in (1, 2, 3):
        assert back.xy_of(n) == m.xy_of(n)
    assert not [r for r in back.log if r.rule == 'container:conformer-narrowed']
