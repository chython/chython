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
import struct

import pytest

from chython.core import STEREO_AND, WEDGE_UP, MoleculeContainer
from chython.core import _core

# Segment table starts at byte 24; each entry is 8 bytes (offset uint32 + length uint32).
# The segment indices below match the cdef enum in _structure.pxd.
_SEG_ATOMS = 0
_SEG_CSR_PTR = 1
_SEG_CSR_EDGE = 2
_SEG_XY = 3
_SEG_STEREO_GROUPS = 4
# Each atom_t record is 24 bytes (pinned by test_struct_sizes_are_locked).
_ATOM_RECORD_SIZE = 24


def _chain(elements):
    """Linear chain of atoms with bond order 1 between consecutive pairs."""
    m = MoleculeContainer()
    ids = [m.add_atom(e) for e in elements]
    for i in range(len(ids) - 1):
        m.add_bond(ids[i], ids[i + 1], 1)
    return m, ids


def loaded():
    """A molecule exercising every persistent field at once."""
    m = MoleculeContainer()
    ring = [m.add_atom(6) for _ in range(6)]
    for i in range(6):
        m.add_bond(ring[i], ring[(i + 1) % 6], 1 if i % 2 else 2)
    n = m.add_atom(7, charge=1, map_number=17)
    o = m.add_atom(8, charge=-1, isotope=18)
    c = m.add_atom(6, isotope=13, radical=True, implicit_h=1)
    m.add_bond(ring[0], n, 1)
    m.add_bond(n, o, 1)
    m.add_bond(ring[3], c, 1)
    for k, s in enumerate(ring):
        m.set_xy(s, 1.5 + k, -2.25 - k)
    m.set_xy(n, 0.0001, -0.0001)
    m.set_wedge(ring[0], n, 1)
    m.set_stereo_group(c, STEREO_AND, 3)
    return m, ring + [n, o, c]


def test_pack_length_is_the_persistent_prefix():
    m, ids = loaded()
    data = m.to_bytes()
    assert isinstance(data, bytes)
    assert len(data) == m.persistent_len
    assert m.persistent_len < m.total_len      # derived segments are excluded


def test_pack_is_deterministic():
    m, ids = loaded()
    assert m.to_bytes() == m.to_bytes()


def test_pack_survives_a_copy():
    m, ids = loaded()
    assert m.copy().to_bytes() == m.to_bytes()


def test_round_trip_preserves_every_persistent_field():
    m, ids = loaded()
    back = MoleculeContainer.from_bytes(m.to_bytes())

    assert back.atom_count == m.atom_count
    assert back.bond_count == m.bond_count
    assert back.atom_numbers == m.atom_numbers
    for s in m.atom_numbers:
        assert back.element_of(s) == m.element_of(s)
        assert back.charge_of(s) == m.charge_of(s)
        assert back.isotope_of(s) == m.isotope_of(s)
        assert back.map_number_of(s) == m.map_number_of(s)
        assert back.radical_of(s) == m.radical_of(s)
        assert back.implicit_h_of(s) == m.implicit_h_of(s)
        assert back.stereo_group_of(s) == m.stereo_group_of(s)
        assert back.xy_of(s) == m.xy_of(s)
    assert sorted(back.wedges()) == sorted(m.wedges())
    for a in m.atom_numbers:
        for b in m.atom_numbers:
            assert back.order_of(a, b) == m.order_of(a, b)


def test_round_trip_rebuilds_derived_layers():
    m, ids = loaded()
    back = MoleculeContainer.from_bytes(m.to_bytes())

    assert back._union_feature_words == m._union_feature_words
    assert back.rings_count == m.rings_count
    assert sorted(sorted(r) for r in back.rings) == sorted(sorted(r) for r in m.rings)
    for s in m.atom_numbers:
        assert back.in_ring_of(s) == m.in_ring_of(s)
        assert back.ring_count_of(s) == m.ring_count_of(s)
        assert back.ring_sizes_of(s) == m.ring_sizes_of(s)
        assert back.degree_of(s) == m.degree_of(s)
        assert back.heteroatoms_of(s) == m.heteroatoms_of(s)
        assert back.hybridization_of(s) == m.hybridization_of(s)
        assert back.total_h_of(s) == m.total_h_of(s)
        assert back.features_of(s) == m.features_of(s)


def test_round_trip_rebuilds_a_multi_word_ring_bitmap():
    # loaded() has 2 rings, so its ring bitmap is one uint64_t per atom and cannot
    # tell the word index from the bit index. Rebuilding after unpack is where a
    # words mismatch is genuinely plausible -- the ring count is recomputed from a
    # buffer rather than carried in the same pass that wrote it -- so the round trip
    # needs at least one molecule with more than 64 relevant cycles.
    #
    # 10x10 grid: 100 atoms, 180 bonds, 81 unit squares, words = ceil(81 / 64) = 2.
    coord = {}
    m = MoleculeContainer()
    with m.edit():
        idx = 0
        for r in range(10):
            for c in range(10):
                coord[(r, c)] = m.add_atom(6)
                idx += 1
        for r in range(10):
            for c in range(10):
                if c + 1 < 10:
                    m.add_bond(coord[(r, c)], coord[(r, c + 1)], 1)
                if r + 1 < 10:
                    m.add_bond(coord[(r, c)], coord[(r + 1, c)], 1)
    assert m.rings_count == 81

    back = MoleculeContainer.from_bytes(m.to_bytes())
    assert back.rings_count == 81
    assert sorted(sorted(r) for r in back.rings) == sorted(sorted(r) for r in m.rings)
    for s in m.atom_numbers:
        assert back.ring_count_of(s) == m.ring_count_of(s)
        assert back.ring_sizes_of(s) == m.ring_sizes_of(s)

    # Rings are ordered by their lowest atom index, so the square whose top-left
    # corner is (tr, tc) is ring tr * 9 + tc. Atom (9, 9) is a grid corner and so
    # belongs to exactly one square -- ring 80, which lives in the SECOND word.
    # Atom (1, 7) belongs to rings 6, 7, 15 and 16 and to nothing else.
    #
    # This pair is the discriminator: it must NOT share a ring. Under an in-bounds
    # word-index bug -- `r >> 7` instead of `r >> 6`, which stays inside the
    # allocation and therefore passes every count and size assertion above --
    # ring 80 folds onto word 0 bit 16 (80 & 63 == 16) and collides with ring 16,
    # so shares_ring would answer True.
    assert not back.shares_ring(coord[(9, 9)], coord[(1, 7)])
    # the true case, also resolving through word 1: the diagonal of ring 80
    assert back.shares_ring(coord[(9, 9)], coord[(8, 8)])


def _wedged_grid():
    """A 10x10 carbon grid (100 atoms, 180 bonds, 81 rings) with one wedged bond.

    Big enough that `rebuild_derived`'s derived segments push the arena past its current
    allocation, so PyMem_Realloc has a real chance of MOVING the buffer -- which is the
    precondition for the defect the test below guards.  The wedge is what makes it the
    right shape: `from_bytes` reads the wedge segment before the rebuild (Ruling F54's
    legacy discriminator) and the stable ids after it.
    """
    coord = {}
    m = MoleculeContainer()
    with m.edit():
        for r in range(10):
            for c in range(10):
                coord[(r, c)] = m.add_atom(6)
        for r in range(10):
            for c in range(10):
                if c + 1 < 10:
                    m.add_bond(coord[(r, c)], coord[(r, c + 1)], 1)
                if r + 1 < 10:
                    m.add_bond(coord[(r, c)], coord[(r + 1, c)], 1)
        m.set_wedge(coord[(0, 0)], coord[(0, 1)], WEDGE_UP)
    return m, coord


# How many times the round trip below is repeated inside the one test.  The assertions are
# deterministic; their DETECTION is not -- see the docstring.  This count is LOAD-BEARING and
# must not be trimmed: detection is a threshold effect on how far the arena has grown, not a
# series of independent trials, so the rate does not decay gracefully as the count drops.
# Measured against a build with the fix reverted: 1 round trip 0/1040, 2 round trips 0%,
# 4 -> 1%, 8 -> 40%, 16 -> 82-100%, 32 -> 99%.  Sixteen is the first count that detects
# reliably, and the whole test costs under 5 ms.
_ARENA_ROUND_TRIPS = 16


def test_wedge_round_trip_keeps_stable_ids():
    """Regression test for a stale arena pointer in `from_bytes` (Ruling F60).

    `from_bytes` fetches the atom array before the legacy wedge normalisation, then calls
    `rebuild_derived`, which appends derived segments and REALLOCATES (and may move) the
    arena.  Between round 2 of task 5 and its fix, the `_numbers` list was then built by
    reading `n` through the pre-rebuild pointer, i.e. out of freed memory -- silently
    wrong stable ids, and hence a wrong `rings`, a wrong `_index_of` and everything keyed by
    them, on `from_bytes` calls for a molecule this size.

    Each assertion here is deterministic, but whether the defect is VISIBLE on any one
    round trip is not: it needs the realloc to actually move the buffer AND the freed block
    to have been reused before it is read.  Both depend on allocator state, so the rate varies
    widely between processes and is not a fixed per-call probability -- measured with the
    re-fetch reverted, over builds made from clean archives: 1177 / 2000 in-process (59%
    overall, but 14% to 99% depending on the process), 39 / 60 across fresh pytest processes
    (65%), and 307 / 400 on a fresh-molecule harness (77%).  With the re-fetch in place,
    **0 / 2200**.  So a green run of this test is evidence, not proof: do not read one as
    showing that no arena pointer is stale -- read the comment at `rebuild_derived` instead.
    """
    m, coord = _wedged_grid()
    expected = m.atom_numbers
    assert len(expected) == 100 and len(set(expected)) == 100
    expected_rings = sorted(sorted(r) for r in m.rings)
    data = m.to_bytes()

    for attempt in range(_ARENA_ROUND_TRIPS):
        back = MoleculeContainer.from_bytes(data)
        assert back.atom_numbers == expected, \
            'stable ids read through a stale arena pointer (attempt %d)' % attempt
        # the ring reader maps atom indices through _numbers, so a corrupt id lands here too
        assert sorted(sorted(r) for r in back.rings) == expected_rings, attempt
        # and the wedge -- read before the rebuild -- must survive it
        assert back.wedge_of(coord[(0, 0)], coord[(0, 1)]) == WEDGE_UP


def test_unpack_recomputes_forged_atom_descriptors():
    m, ids = loaded()
    buf = bytearray(m.to_bytes())
    # the header is self-describing: the segment table starts at byte 24 and
    # SEG_ATOMS is entry 0, so its offset is the first uint32 there
    atoms_at = struct.unpack_from('<I', buf, 24)[0]
    # atom_t is a *packed* struct, so the layout is the cumulative field sum:
    # element 0, charge 1, hydrogens 2, flags 3, isotope 4, map_number 6,
    # n 8, degree 12, heteroatoms 13, ring_sizes 14, ring_counts 18,
    # reserved 20 -- 24 bytes total, which test_struct_sizes_are_locked pins
    for i in range(len(ids)):
        rec = atoms_at + 24 * i
        buf[rec + 12] = 99                                   # degree
        buf[rec + 13] = 99                                    # heteroatoms
        struct.pack_into('<I', buf, rec + 14, 0xffffffff)      # every ring-size bit
        struct.pack_into('<H', buf, rec + 18, 200)            # 200 rings
    back = MoleculeContainer.from_bytes(bytes(buf))
    for s in m.atom_numbers:
        assert back.degree_of(s) == m.degree_of(s)
        assert back.heteroatoms_of(s) == m.heteroatoms_of(s)
        assert back.ring_sizes_of(s) == m.ring_sizes_of(s)
        assert back.ring_count_of(s) == m.ring_count_of(s)


def test_unpack_clears_forged_ring_flags():
    m, ids = loaded()
    ring, n, o, c = ids[:6], ids[6], ids[7], ids[8]
    buf = bytearray(m.to_bytes())
    atoms_at = struct.unpack_from('<I', buf, 24)[0]
    edges_at = struct.unpack_from('<I', buf, 24 + 8 * 2)[0]   # SEG_CSR_EDGE is entry 2
    for i in range(len(ids)):
        buf[atoms_at + 24 * i + 3] |= 0x04                   # atom_t.flags bit 2, in_ring
    # halfedge_t is packed too: to 0, order 4, wedge 5, flags 6 -- 8 bytes
    for k in range(2 * m.bond_count):
        at = edges_at + 8 * k + 6
        struct.pack_into('<H', buf, at, struct.unpack_from('<H', buf, at)[0] | 1)
    back = MoleculeContainer.from_bytes(bytes(buf))
    # every bridge comes back acyclic even though the buffer claimed otherwise
    assert back.in_ring_of(n) is False
    assert back.in_ring_of(o) is False
    assert back.in_ring_of(c) is False
    assert back.bond_in_ring(ring[0], n) is False
    assert back.bond_in_ring(n, o) is False
    assert back.bond_in_ring(ring[3], c) is False
    # and the real ring survives -- the pass distinguishes, it does not blanket-clear
    assert all(back.in_ring_of(s) for s in ring)
    assert back.bond_in_ring(ring[0], ring[1]) is True


def test_round_trip_is_idempotent():
    m, ids = loaded()
    once = MoleculeContainer.from_bytes(m.to_bytes())
    assert once.to_bytes() == m.to_bytes()
    assert MoleculeContainer.from_bytes(once.to_bytes())._union_feature_words == m._union_feature_words


def test_stated_hydrogens_and_the_UNSTATED_SENTINEL_survive_the_round_trip():
    """A stated count and an unstated one are two different stored values, and both cross the bytes.

    The unstated atom asserted `== 0` here until the `add_atom` default moved to `H_UNKNOWN`, which
    is why this test is worth having at all: the thing being round-tripped is a DISTINCTION, and a
    zero it shared with a stated zero was a distinction the bytes were not really carrying.  Now the
    nibble holds 15 for the one and 1 for the other, and a pack that dropped the sentinel would show
    up as an unstated atom coming back with a count.
    """
    m = MoleculeContainer()
    stated = m.add_atom(7, implicit_h=1)
    unstated = m.add_atom(7)                 # nobody said anything about this one
    h = m.add_atom(1)
    m.add_bond(stated, unstated, 1)
    m.add_bond(unstated, h, 1)
    back = MoleculeContainer.from_bytes(m.to_bytes())
    assert back.implicit_h_of(stated) == 1
    assert back.implicit_h_of(unstated) is None
    # explicit H is derived — unpack recomputed it from the restored bonds
    assert back.explicit_h_of(unstated) == 1
    assert back.explicit_h_of(stated) == 0
    # the nibble rode along in the atom record, so folding a fresh arena still tells the stated
    # count apart from the unstated one
    back.set_map_number(stated, 0)   # any mutation re-derives from the primary record
    assert back.implicit_h_of(stated) == 1
    assert back.implicit_h_of(unstated) is None


def test_stable_ids_survive_a_deletion_gap():
    m = MoleculeContainer()
    a, b, c = m.add_atom(6), m.add_atom(6), m.add_atom(6)
    m.add_bond(a, b, 1)
    m.add_bond(b, c, 1)
    m.delete_atom(b)
    assert m.atom_numbers == [a, c]
    assert MoleculeContainer.from_bytes(m.to_bytes()).atom_numbers == [a, c]


def test_molecule_without_coordinates_round_trips():
    m = MoleculeContainer()
    a, b = m.add_atom(6), m.add_atom(8)
    m.add_bond(a, b, 2)
    assert not m.has_coordinates
    back = MoleculeContainer.from_bytes(m.to_bytes())
    assert not back.has_coordinates
    assert back.xy_of(a) is None
    assert back.order_of(a, b) == 2


def test_empty_molecule_round_trips():
    m = MoleculeContainer()
    back = MoleculeContainer.from_bytes(m.to_bytes())
    assert back.atom_count == 0
    assert back.atom_numbers == []
    assert back.rings_count == 0


def test_persistent_view_matches_pack():
    m, ids = loaded()
    view = m.persistent_view
    assert view.readonly
    assert bytes(view) == m.to_bytes()
    assert len(view) == m.persistent_len


def test_unpack_rejects_empty_input():
    with pytest.raises(ValueError, match='too short'):
        MoleculeContainer.from_bytes(b'')


def test_unpack_rejects_a_truncated_header():
    """Truncated below the 24 fixed bytes, where the segment table's own length is stated.

    Not 64 bytes: `loaded()` has coordinates but no stereo groups and no S-groups, so it spends
    four table entries for a 56-byte header, and a 64-byte slice is a truncated PAYLOAD caught by
    the length-versus-`persistent_len` check instead.  Truncation inside the table is covered
    separately below.
    """
    m, ids = loaded()
    with pytest.raises(ValueError, match='too short'):
        MoleculeContainer.from_bytes(m.to_bytes()[:20])


def test_unpack_rejects_a_bad_magic():
    m, ids = loaded()
    data = bytearray(m.to_bytes())
    data[0] ^= 0xff
    with pytest.raises(ValueError, match='magic'):
        MoleculeContainer.from_bytes(bytes(data))


def test_unpack_rejects_a_future_version():
    m, ids = loaded()
    data = bytearray(m.to_bytes())
    data[4] = 99
    with pytest.raises(ValueError, match='version'):
        MoleculeContainer.from_bytes(bytes(data))


def _header_len(data):
    """24 fixed bytes plus the table entries this buffer actually declares."""
    return 24 + 8 * struct.unpack_from('<H', data, 20)[0]


def test_the_segment_table_stops_after_the_last_segment_the_molecule_uses():
    """A molecule pays for the entries it uses and not for the thirteen the struct reserves.

    `seg_count` is written as one past the highest non-empty persistent id, so the header is 48
    bytes for a molecule with no coordinates, 56 with them.  A fixed 128-byte header costs a 13-atom
    aspirin 80 bytes of pure padding on a 704-byte record, which is what the sizing buys.  The floor
    is three entries because SEG_ATOMS, SEG_CSR_PTR and SEG_CSR_EDGE are dereferenced
    unconditionally by every reader.
    """
    m, ids = _chain('CCO')
    assert not m.has_coordinates
    data = m.to_bytes()
    assert struct.unpack_from('<H', data, 20)[0] == 3
    assert _header_len(data) == 48
    # the atom payload begins where the table ends, with no gap
    assert struct.unpack_from('<I', data, 24 + 8 * _SEG_ATOMS)[0] == 48

    for sid in ids:
        m.set_xy(sid, 1.0, 2.0)
    data = m.to_bytes()
    assert struct.unpack_from('<H', data, 20)[0] == 4
    assert _header_len(data) == 56
    assert MoleculeContainer.from_bytes(data).xy_of(ids[0]) == (1.0, 2.0)


def test_an_interior_empty_entry_keeps_its_slot():
    """Trailing empties are dropped; an interior one is not, because ids are positional.

    A molecule with stereo groups (id 4) and no coordinates (id 3) must still spend five entries
    with the fourth empty -- renumbering to close the hole would make id 4 mean id 3 to every
    reader, including one reading a buffer written by an older build.  This is the case that
    distinguishes "truncate the tail" from "compact the table", and only the first is legal.
    """
    m, ids = _chain('CCO')
    m.set_stereo_group(ids[0], STEREO_AND, 1)
    assert not m.has_coordinates
    data = m.to_bytes()
    assert struct.unpack_from('<H', data, 20)[0] == 5
    assert struct.unpack_from('<II', data, 24 + 8 * _SEG_XY)[1] == 0, 'the hole must stay a hole'
    assert struct.unpack_from('<II', data, 24 + 8 * _SEG_STEREO_GROUPS)[1] != 0

    back = MoleculeContainer.from_bytes(data)
    assert not back.has_coordinates
    assert back.stereo_group_of(ids[0]) == (STEREO_AND, 1)


def test_unpack_rejects_a_table_too_short_for_the_csr_segments():
    """A hand-made buffer may not claim fewer than three entries.

    Without this the entries for SEG_CSR_PTR and SEG_CSR_EDGE would be read out of the atom
    payload and used as unvalidated offsets -- the validation loop skips entries past `seg_count`,
    so nothing else would look at them.  The writer never emits fewer than three; this is what
    stops a forged buffer from saying it did.
    """
    m, ids = loaded()
    for claimed in (1, 2):
        data = bytearray(m.to_bytes())
        struct.pack_into('<H', data, 20, claimed)
        with pytest.raises(ValueError, match='segment table'):
            MoleculeContainer.from_bytes(bytes(data))
    data = bytearray(m.to_bytes())
    struct.pack_into('<H', data, 20, 0)
    with pytest.raises(ValueError, match='empty segment table'):
        MoleculeContainer.from_bytes(bytes(data))


def test_unpack_rejects_a_buffer_shorter_than_the_table_it_declares():
    """`seg_count` is data, so it can name a table that does not fit -- checked before it is read."""
    m, ids = loaded()
    data = m.to_bytes()
    with pytest.raises(ValueError, match='too short'):
        MoleculeContainer.from_bytes(data[:_header_len(data) - 8])


def test_unpack_rejects_a_table_claiming_a_segment_this_build_cannot_model():
    """Forward compatibility as a refusal: an entry above the persistent block must be EMPTY.

    A buffer from a later release that puts real data in the first entry above the persistent block
    carries information this build cannot model, and writing the molecule back would drop it silently.
    The table is grown here by raising `seg_count` past SEG_PERSISTENT_COUNT and giving the new entry a
    length, which is what such a buffer would look like from here.

    EVERY NUMBER HERE IS DERIVED FROM `_persistent_segment_count()`, and the two that were literals --
    a `grown = 9` and an `entry 8` -- came due the moment `SEG_CONFORMERS` was appended: entry 8 became
    a segment this build models perfectly well, so the test was asserting a refusal of something no
    longer refusable and failed for a change that broke nothing.  The fact being pinned is "the FIRST
    unmodellable entry", which is a relation to the size of the persistent block and never a number.
    """
    m, ids = loaded()
    data = bytearray(m.to_bytes())
    unknown = _core._persistent_segment_count()     # the lowest id this build cannot model
    grown = unknown + 1                             # a table one entry longer than that
    head = _header_len(data)
    struct.pack_into('<H', data, 20, grown)
    # the table now overlaps the payload, so rebuild the buffer with the header widened
    data = bytearray(data[:head]) + bytearray(8 * (grown - (head - 24) // 8)) + bytearray(data[head:])
    struct.pack_into('<I', data, 16, len(data))                     # persistent_len
    for seg in range(grown):
        offset = struct.unpack_from('<I', data, 24 + 8 * seg)[0]
        if offset:
            struct.pack_into('<I', data, 24 + 8 * seg, offset + 8 * (grown - (head - 24) // 8))
    struct.pack_into('<II', data, 24 + 8 * unknown, 0, 8)           # that entry: non-empty
    with pytest.raises(ValueError, match=f'unknown segment {unknown}'):
        MoleculeContainer.from_bytes(bytes(data))


def test_unpack_rejects_a_length_mismatch():
    m, ids = loaded()
    with pytest.raises(ValueError, match='length'):
        MoleculeContainer.from_bytes(m.to_bytes() + b'\x00' * 8)
    with pytest.raises(ValueError, match='length'):
        MoleculeContainer.from_bytes(m.to_bytes()[:-8])


def test_unpack_rejects_a_segment_out_of_bounds():
    m, ids = loaded()
    data = bytearray(m.to_bytes())
    data[24] = 0xff                     # first segment offset, low byte
    data[25] = 0xff
    data[26] = 0xff
    with pytest.raises(ValueError, match='segment'):
        MoleculeContainer.from_bytes(bytes(data))


def test_unpack_rejects_a_dirty_topology_flag():
    m, ids = loaded()
    data = bytearray(m.to_bytes())
    data[6] |= 2                        # FLAG_TOPOLOGY_DIRTY in header.flags
    with pytest.raises(ValueError, match='dirty'):
        MoleculeContainer.from_bytes(bytes(data))


def test_unpack_rejects_a_wrong_atom_count():
    m, ids = loaded()
    data = bytearray(m.to_bytes())
    data[8] = (data[8] + 1) & 0xff      # atom_count, low byte
    with pytest.raises(ValueError, match='segment'):
        MoleculeContainer.from_bytes(bytes(data))


def test_unpack_rejects_broken_stable_ids():
    # atom_t is 24 bytes with `n` at offset 8. All three branches are here because
    # each one is undetectable later: a zero or duplicate id leaves the index dict short,
    # and 0xFFFFFFFF wraps unpack's `_next_id = high + 1` to 0 so the next add_atom starts
    # reissuing ids that the packed buffer already used.
    m, ids = loaded()
    atoms_at = struct.unpack_from('<I', m.to_bytes(), 24)[0]        # SEG_ATOMS offset

    data = bytearray(m.to_bytes())
    struct.pack_into('<I', data, atoms_at + 8, 0)               # atom 0
    with pytest.raises(ValueError, match='never issued'):
        MoleculeContainer.from_bytes(bytes(data))

    data = bytearray(m.to_bytes())
    struct.pack_into('<I', data, atoms_at + 8, 0xFFFFFFFF)
    with pytest.raises(ValueError, match='reserved stable id'):
        MoleculeContainer.from_bytes(bytes(data))

    data = bytearray(m.to_bytes())
    first = struct.unpack_from('<I', data, atoms_at + 8)[0]
    struct.pack_into('<I', data, atoms_at + 24 + 8, first)      # atom 1 aliases atom 0
    with pytest.raises(ValueError, match='appears twice'):
        MoleculeContainer.from_bytes(bytes(data))


def test_unpack_rejects_a_bond_count_that_would_wrap():
    # bond_count is a uint32 at header offset 12. 0x80000000 is the value that makes
    # `2 * bond_count` wrap to 0 in 32-bit arithmetic -- which would zero out both the
    # edge-segment size check and the half-edge walk, adopting the buffer while its
    # header claimed 2.1 billion bonds. The size_t casts in structure_from_bytes are
    # what turn this into a rejection.
    m, ids = loaded()
    data = bytearray(m.to_bytes())
    struct.pack_into('<I', data, 12, 0x80000000)
    with pytest.raises(ValueError, match='csr edge segment is too small'):
        MoleculeContainer.from_bytes(bytes(data))


def test_unpack_rejects_an_out_of_range_element():
    m, ids = loaded()
    atoms_at = struct.unpack_from('<I', m.to_bytes(), 24)[0]   # SEG_ATOMS offset
    # element 0 is now R (valid); 119 and 255 are truly out of range
    for forged in (119, 255):
        data = bytearray(m.to_bytes())
        data[atoms_at] = forged                            # atom 0, element byte
        with pytest.raises(ValueError, match='element'):
            MoleculeContainer.from_bytes(bytes(data))


def test_unpacked_molecule_is_mutable_again():
    m, ids = loaded()
    back = MoleculeContainer.from_bytes(m.to_bytes())
    back.set_charge(ids[0], 1)
    assert back.charge_of(ids[0]) == 1
    assert back.rings_count == 1
    # unpack restored the stable id counter past every id in the packed arena
    assert back.add_atom(6) not in ids


# ---------------------------------------------------------------------------
# Fix round 1: Item 1 — CSR symmetry (heap out-of-bounds write in mark_bridges)
# ---------------------------------------------------------------------------

def test_unpack_rejects_a_half_edge_with_no_twin():
    # One edited `to` byte in a legitimate pack() reaches an out-of-bounds write in
    # mark_bridges, which stores through the failed-twin-search sentinel ptr[child + 1].
    m = MoleculeContainer()
    with m.edit():
        ids = [m.add_atom(6) for _ in range(6)]
        m.add_bond(ids[0], ids[1], 1)
        m.add_bond(ids[0], ids[5], 1)
        m.add_bond(ids[1], ids[5], 1)
    data = bytearray(m.to_bytes())
    edge_offset = struct.unpack_from('<I', data, 24 + 8 * _SEG_CSR_EDGE)[0]
    struct.pack_into('<I', data, edge_offset + 8 * 5, 0)
    with pytest.raises(ValueError):
        MoleculeContainer.from_bytes(bytes(data))


def test_unpack_rejects_a_self_loop_half_edge():
    m, ids = _chain([6, 6, 6])
    data = bytearray(m.to_bytes())
    edge_offset = struct.unpack_from('<I', data, 24 + 8 * _SEG_CSR_EDGE)[0]
    struct.pack_into('<I', data, edge_offset, 0)          # atom 0's half-edge points at itself
    with pytest.raises(ValueError):
        MoleculeContainer.from_bytes(bytes(data))


def test_unpack_rejects_half_edges_out_of_order():
    # csr_build always leaves each atom's range sorted by `to`; an unsorted range means the
    # buffer did not come from csr_build, and the cursor twin search assumes the order.
    m = MoleculeContainer()
    with m.edit():
        ids = [m.add_atom(6) for _ in range(4)]
        m.add_bond(ids[0], ids[1], 1)
        m.add_bond(ids[0], ids[2], 1)
        m.add_bond(ids[0], ids[3], 1)
    data = bytearray(m.to_bytes())
    edge_offset = struct.unpack_from('<I', data, 24 + 8 * _SEG_CSR_EDGE)[0]
    first = bytes(data[edge_offset:edge_offset + 8])
    second = bytes(data[edge_offset + 8:edge_offset + 16])
    data[edge_offset:edge_offset + 8] = second
    data[edge_offset + 8:edge_offset + 16] = first
    with pytest.raises(ValueError):
        MoleculeContainer.from_bytes(bytes(data))


def test_unpack_rejects_half_edges_that_disagree_on_order():
    # A bond has one order, so both halves must carry it. csr_build writes the same value
    # into both and no mutator can separate them, but an edited buffer can -- and then
    # order_of(a, b) and order_of(b, a) answer differently, with hybridization, the feature
    # words and the union row all derived from a graph that is not a graph.
    m, ids = _chain([6, 6, 6])
    data = bytearray(m.to_bytes())
    edge_offset = struct.unpack_from('<I', data, 24 + 8 * _SEG_CSR_EDGE)[0]
    struct.pack_into('<B', data, edge_offset + 4, 2)      # atom 0's half only
    with pytest.raises(ValueError, match='disagree on order'):
        MoleculeContainer.from_bytes(bytes(data))


def test_unpack_keeps_a_molecule_whose_bonds_have_genuinely_mixed_orders():
    # The order-symmetry check compares a half-edge against its twin, not against its
    # neighbours, so a molecule with several different bond orders must still round trip.
    m = MoleculeContainer()
    with m.edit():
        ids = [m.add_atom(6) for _ in range(4)]
        m.add_bond(ids[0], ids[1], 1)
        m.add_bond(ids[1], ids[2], 2)
        m.add_bond(ids[2], ids[3], 3)
    back = MoleculeContainer.from_bytes(m.to_bytes())
    assert [back.order_of(a, b) for a, b in ((ids[0], ids[1]), (ids[1], ids[2]),
                                             (ids[2], ids[3]))] == [1, 2, 3]
    assert [back.order_of(b, a) for a, b in ((ids[0], ids[1]), (ids[1], ids[2]),
                                             (ids[2], ids[3]))] == [1, 2, 3]


# ---------------------------------------------------------------------------
# Fix round 1: Item 2 — optional segment size checks
# ---------------------------------------------------------------------------

def test_unpack_rejects_a_coordinate_segment_too_small_for_its_atoms():
    m, ids = _chain([6, 6, 6])
    for i in ids:
        m.set_xy(i, 1.0, 2.0)
    data = bytearray(m.to_bytes())
    struct.pack_into('<I', data, 24 + 8 * _SEG_XY + 4, 8)   # claim room for one atom
    with pytest.raises(ValueError):
        MoleculeContainer.from_bytes(bytes(data))


# ---------------------------------------------------------------------------
# Fix round 1: Item 3 — stable-id space exhaustion
# ---------------------------------------------------------------------------

def test_stable_id_space_exhaustion_raises_instead_of_wrapping():
    # Forging the largest issuable id leaves _next_id one short of the wrap. The next
    # add_atom must refuse rather than reissue 0 and alias an atom already in the buffer.
    m, ids = _chain([6, 6])
    data = bytearray(m.to_bytes())
    atom_offset = struct.unpack_from('<I', data, 24 + 8 * _SEG_ATOMS)[0]
    struct.pack_into('<I', data, atom_offset + 8, 0xFFFFFFFE)
    back = MoleculeContainer.from_bytes(bytes(data))
    # edit() returns False from __exit__ on exception, so it propagates unchanged
    with pytest.raises(OverflowError):
        with back.edit():
            back.add_atom(6)


# ---------------------------------------------------------------------------
# Fix round 1: Item 4 — fill_ring_descriptors must be unconditional
# ---------------------------------------------------------------------------

def test_unpack_clears_forged_ring_data_on_an_acyclic_molecule():
    # The cyclic forgery test cannot pin this: with cycles present, fill_ring_descriptors
    # runs either way. Only an acyclic molecule distinguishes "cleared unconditionally"
    # from "cleared when there was something to fill", and the acyclic case is the one a
    # foreign buffer exploits.
    m, ids = _chain([6, 6, 6])
    clean = m._union_feature_words
    data = bytearray(m.to_bytes())
    atom_offset = struct.unpack_from('<I', data, 24 + 8 * _SEG_ATOMS)[0]
    for i in range(3):
        base = atom_offset + i * _ATOM_RECORD_SIZE
        data[base + 3] |= 0x04                              # in_ring flag
        struct.pack_into('<I', data, base + 14, 0xFFFFFFFF)  # ring_sizes
        struct.pack_into('<H', data, base + 18, 200)         # ring_counts
    back = MoleculeContainer.from_bytes(bytes(data))
    for i in back.atom_numbers:
        assert back.ring_sizes_of(i) == frozenset()
    assert back._union_feature_words == clean


# ---------------------------------------------------------------------------
# Fix round 1: Item 5 — write-path field ranges and HE_AROMATIC
# ---------------------------------------------------------------------------

def test_unpack_rejects_fields_the_write_path_cannot_produce():
    # Each forgery is a value add_atom/add_bond/set_wedge already refuse. A validator that
    # admits more than the write path produces admits states no other code is written for.
    # HE_AROMATIC (flags bit 1 = 2): mark_bridges only clears HE_IN_RING (bit 0), so a
    # forged aromatic flag would survive into feature word I bit 62 and break invariants.
    m, ids = _chain([6, 6, 6])
    base_packed = m.to_bytes()
    atom_offset = struct.unpack_from('<I', bytearray(base_packed), 24 + 8 * _SEG_ATOMS)[0]
    edge_offset = struct.unpack_from('<I', bytearray(base_packed), 24 + 8 * _SEG_CSR_EDGE)[0]

    def forged(offset, fmt, value):
        data = bytearray(base_packed)
        struct.pack_into(fmt, data, offset, value)
        return bytes(data)

    cases = [(edge_offset + 4, '<B', 0),                       # bond order 0
             (edge_offset + 4, '<B', 200),                     # bond order 200
             (edge_offset + 5, '<B', 200),                     # wedge 200
             (edge_offset + 6, '<H', 2),                       # HE_AROMATIC
             (atom_offset + 1, '<b', -128),                    # charge
             (atom_offset + 6, '<H', 65535),                   # map number
             # `reserved`'s low nibble is the atom's CIP code, so 1 is a LEGAL value there ('R') and
             # is absent from this list.  What is refused is a code with no descriptor (9 is one past
             # the last), and any bit above the defined nibble.
             (atom_offset + 20, '<I', 9),                       # CIP code past the domain
             (atom_offset + 20, '<I', 0x10),                    # reserved bit above the CIP nibble
             (atom_offset + 20, '<I', 0x80000000)]              # reserved top bit
    for offset, fmt, value in cases:
        with pytest.raises(ValueError):
            MoleculeContainer.from_bytes(forged(offset, fmt, value))


# ---------------------------------------------------------------------------
# segment alignment
# ---------------------------------------------------------------------------

def test_unpack_rejects_a_misaligned_segment_offset():
    m, ids = _chain([6, 6, 6])
    data = bytearray(m.to_bytes())
    offset = struct.unpack_from('<I', data, 24 + 8 * _SEG_ATOMS)[0]
    struct.pack_into('<I', data, 24 + 8 * _SEG_ATOMS, offset + 1)
    with pytest.raises(ValueError):
        MoleculeContainer.from_bytes(bytes(data))


def test_unpack_rejects_a_segment_length_that_is_not_8_aligned():
    """A length is checked on the READER, not argued from the writer's align8.

    `structure_respan` memcpy's the DESTINATION's align8 length out of the source segment, so a
    source whose stated length is short of a multiple of 8 is read past its end.  A forged length is
    a FORMAT violation and the reader refuses it, which is the one class of refusal a reader makes.

    Three segments, because the check is in the loop over the whole persistent table rather than an
    arm for one id: it must fire wherever the forged entry is.
    """
    m, ids = _chain([6, 6, 6])
    data = bytearray(m.to_bytes())
    for seg in (_SEG_ATOMS, _SEG_CSR_PTR, _SEG_CSR_EDGE):
        d = bytearray(data)
        length = struct.unpack_from('<I', d, 24 + 8 * seg + 4)[0]
        assert length and not length & 7, 'segment %d must be present and aligned to begin with' % seg
        struct.pack_into('<I', d, 24 + 8 * seg + 4, length - 1)
        with pytest.raises(ValueError, match='not 8-aligned'):
            MoleculeContainer.from_bytes(bytes(d))


# ---------------------------------------------------------------------------
# reserved flag bits
# ---------------------------------------------------------------------------

def test_a_version_5_buffer_may_not_set_the_reserved_flag_bits():
    """Bits 1 and 7 of `atom_t.flags` are reserved from version 5 on, and a reader says so.

    The same rule the segment table has: a field this build does not model must be empty, because
    writing the molecule back would silently drop whatever it meant.  Version 3 and version 4 buffers
    are exempt -- those bits ARE their parity, and `from_bytes` adopts them.
    """
    m, ids = _chain([6, 6, 6])
    data = bytearray(m.to_bytes())
    atoms_at = struct.unpack_from('<I', data, 24 + 8 * _SEG_ATOMS)[0]
    for bits in (0x02, 0x80, 0x82):
        d = bytearray(data)
        d[atoms_at + 24 * 0 + 3] |= bits
        with pytest.raises(ValueError, match='reserved flag bits'):
            MoleculeContainer.from_bytes(bytes(d))
