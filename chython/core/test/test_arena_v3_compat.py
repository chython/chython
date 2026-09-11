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
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
"""A current build must read the frozen version-3 buffers, and the evidence is real bytes rather than an argument.

`v3_fixtures.py` holds seven records serialised by the last v3 build (75f6c60), each captured twice
-- `cold`, straight out of the builder, and `warm`, after three reads had built derived segments and
written their table entries into the persistent prefix. Both forms exist in the wild, so both must
read, and the warm form is the interesting one: v3's table entries 5..11 hold DERIVED offsets that
point past `persistent_len`, whereas v4's ids 5..8 are persistent segments. A v4 reader that
validated a warm v3 buffer against its own persistent id range would reject it for carrying an
S-group segment it does not have.

Alongside the bytes each fixture carries `answers` -- thirty-four values the v3 build computed from
that record, from element counts through canonical order to `validate_stereo`. The replay below uses
the generator's own `snapshot()` function, so the comparison is against the same code that produced
the frozen values and a drift in either direction shows up as a key mismatch rather than as silence.

Exactly ONE of those values is exempt, named as a `(fixture, key)` pair in `MOVED` and justified by a
test of its own: closing the mirror-automorphism defect in the canonical search moved
`a_bigger_one`'s canonical order to the other of two labellings related by a reflection of its
macrocycle. The bytes are untouched and so is every other answer, `atoms_order` included -- what
moved is which member of a tie an extremal search takes, which is what the fix was for.

WHAT MAKES THIS EVIDENCE RATHER THAN DECORATION (ruling F102). Two of the tests here exist only to
show the suite could have failed: `test_a_single_flipped_payload_byte_does_not_go_unnoticed` sweeps
every byte of every persistent payload and requires each flip to be REJECTED or to change an answer,
with an explicit, justified list of the byte positions that are blind -- and
`test_the_replay_can_detect_a_wrong_answer` corrupts an expected value and watches the comparison
catch it.
"""
from struct import unpack_from

from pytest import mark, raises

from chython.core import MoleculeContainer
from chython.core._core import read_smiles, write_smiles, STRUCT_VERSION, _parity_bytes

from .gen_v3_fixtures import snapshot
from .v3_fixtures import V3_FIXTURES


ALL = sorted(V3_FIXTURES)
BOTH = [(name, warmth) for name in ALL for warmth in ('cold', 'warm')]

# THE ONE FROZEN ANSWER A LATER BUILD IS ALLOWED TO MOVE, LISTED RATHER THAN TOLERATED.
# `a_bigger_one` is a twenty-membered ring carrying twenty methyls and twenty configured centres, and
# the v3 build's canonical order for it was one of TWO extremal labellings related by a reflection of
# the macrocycle -- the mirror automorphism defect. Closing that defect (a parity tail on the leaf
# certificate, an orbit prune refined by parity) made the search choose the other one. The bytes did
# not move and neither did any other answer, including `atoms_order`: the CONSTITUTIONAL colouring is
# identical, so what changed is which member of a tie the search takes. It is named here as a pair
# rather than excluded as a key, so a drift in `canonical_order` on any of the other six fixtures --
# none of which has a symmetry that could excuse one -- is still a failure, and
# `test_the_one_moved_answer_is_a_REFLECTION_and_nothing_else_moved` both justifies this entry and
# fails if the answer ever comes back, so the entry cannot outlive its reason.
MOVED = {('a_bigger_one', 'canonical_order')}

V3_HEADER_LEN = 128
V3_SEG_COUNT = 13
ATOM_RECORD = 24
HALFEDGE = 8
SEG_ATOMS = 0        # persistent id 0, the same in both versions
ATOM_FLAGS = 3       # flags at byte 3 of a packed record
RESERVED_FLAGS = 0x82   # bits 1 and 7: reserved now, a parity in v3 and v4


def _reserved_bits_masked(records):
    """The SOURCE atoms segment with every record's reserved flag bits cleared.

    Used to compare the output against the source: the output's reserved bits are asserted zero
    directly (ingest must clear them), and this mask makes the source comparable by zeroing the bits
    the source may have set as a v3/v4 parity.  Masked rather than skipped so the other 23 bytes of
    every record still have to survive exactly.
    """
    out = bytearray(records)
    for i in range(ATOM_FLAGS, len(out), ATOM_RECORD):
        out[i] &= ~RESERVED_FLAGS & 0xFF
    return bytes(out)


def _table(data, count=V3_SEG_COUNT):
    """The segment table as [(offset, length)], parsed by hand from the bytes.

    Deliberately not asked of the module under test: the point of these tests is to check the reader
    against an independent account of the layout, and a helper that called into the reader would
    agree with it by construction.  `count` defaults to v3's fixed thirteen entries; a current buffer
    states its own in bytes 20-21.
    """
    return [unpack_from('<II', data, 24 + 8 * i) for i in range(count)]


@mark.parametrize('name,warmth', BOTH)
def test_the_fixtures_really_are_v3_bytes(name, warmth):
    """Anti-vacuity for the whole module: what follows proves nothing unless these are v3 buffers.

    If someone reran the generator on a v4 build -- which its docstring forbids and this catches --
    every other test here would still pass while measuring a v4 reader against v4 bytes.
    """
    data = V3_FIXTURES[name][warmth]
    assert unpack_from('<I', data, 0)[0] == 0x43485933, 'not an arena buffer at all'
    assert unpack_from('<H', data, 4)[0] == 3, 'the fixtures are no longer v3 bytes'


@mark.parametrize('name', ALL)
def test_the_warm_form_really_differs_from_the_cold_one(name):
    """The other half of the anti-vacuity: `warm` must be a DIFFERENT byte string from `cold`.

    Were they equal, the fourteen cases below would be seven cases run twice and the warm-buffer
    hazard -- derived table entries naming offsets a v4 reader would misread as S-group segments --
    would never be exercised. They differ because of the v3 defect this branch fixed: a read wrote
    into the persistent prefix. So this assertion is a record of that defect, and it is frozen; it
    can never be satisfied by a post-v3 build, only by the bytes in the fixture file.
    """
    cold, warm = V3_FIXTURES[name]['cold'], V3_FIXTURES[name]['warm']
    assert cold != warm, 'the fixture was captured without the reads that make a warm buffer'
    assert len(cold) == len(warm), 'v3 wrote the derived segments past persistent_len, so the ' \
                                   'serialised prefix was the same length either way'
    assert cold[V3_HEADER_LEN:] == warm[V3_HEADER_LEN:], \
        'the v3 defect was confined to the header; a payload difference means something else'


@mark.parametrize('name,warmth', BOTH)
def test_a_v3_record_reproduces_every_answer_the_v3_build_gave(name, warmth):
    """The headline: thirty-four answers per record, recomputed on this build, byte-for-byte equal.

    This is the demonstration the brief asked for in place of an argument. It covers constitution
    (elements, charges, bonds, hydrogens), the derived scalars (hybridization, heteroatom and degree
    counts), ring perception, components, canonical order and the whole stereo surface including
    `validate_stereo`, which MUTATES -- so the replay exercises the clone-and-clear path on a v3
    buffer too, not only the read paths.

    One answer is exempt and it is exempt BY NAME: see `MOVED` above. Everything else, on every
    fixture, is compared.
    """
    expected = V3_FIXTURES[name]['answers']
    got = snapshot(MoleculeContainer.from_bytes(V3_FIXTURES[name][warmth]))
    wrong = sorted(k for k in expected if got[k] != expected[k] and (name, k) not in MOVED)
    assert not wrong, 'this build disagrees with the v3 build on %r' % (wrong,)
    assert len(expected) >= 30, 'the frozen answer set has shrunk; this compares almost nothing'
    assert set(got) == set(expected), 'snapshot() changed shape, so the fixtures are stale'


@mark.parametrize('warmth', ('cold', 'warm'))
def test_the_one_moved_answer_is_a_REFLECTION_and_nothing_else_moved(warmth):
    """The justification for `MOVED`, and the thing that deletes it if the answer ever comes back.

    An exemption list is a liability unless something proves the exemption is earned, so this test
    makes the case in four measurements rather than in prose.

    IT REALLY MOVED. Asserted first, because if the new order agreed with the v3 one again the entry in
    `MOVED` would be dead weight silently hiding a future drift, and this is the assertion that says
    so out loud.

    NOTHING ELSE MOVED. `canonical_order` is the only differing key out of the thirty-four, and
    `atoms_order` -- the constitutional refinement the canonical search starts from -- is IDENTICAL.
    That pair is what separates a tie-break moving from a perception bug: had ring perception, valence
    or the parity store drifted, the refinement would have drifted with them.

    THE TWO ORDERS ARE RELATED BY A GRAPH AUTOMORPHISM OF ORDER TWO. Compose the v3 labelling with the
    inverse of the new one and the result is a bijection that preserves element, charge, hydrogen count
    and the whole adjacency relation, and squares to the identity: a REFLECTION of the macrocycle. So
    both labellings are extremal labellings of one graph -- the v3 build was not wrong about the graph,
    it just had no way to choose between two mirror candidates, which is the defect by name. An
    arbitrary renumbering would fail this assertion, so it is not a formality.

    AND THE NEW ONE IS THE STABLE ONE. The molecule's canonical SMILES is a fixed point of write-read-
    write, and the molecule read back out of it is `==` to the original with the same
    `canonical_bytes`. Reading a string is a genuinely different presentation -- slots come from string
    position -- so this is the invariant the fix was specified against, checked on the one record in
    this file big enough to have exercised the defect.
    """
    name = 'a_bigger_one'
    expected = V3_FIXTURES[name]['answers']
    m = MoleculeContainer.from_bytes(V3_FIXTURES[name][warmth])
    got = snapshot(m)

    assert got['canonical_order'] != expected['canonical_order'], \
        'the new order agrees with the v3 one again -- delete this test and the MOVED entry'
    assert sorted(k for k in expected if got[k] != expected[k]) == ['canonical_order']
    assert got['atoms_order'] == expected['atoms_order'], \
        'the constitutional refinement moved too, so this is not a tie-break'

    old, new = expected['canonical_order'], got['canonical_order']
    assert sorted(new.values()) == list(range(m.atom_count)), 'not a labelling at all'
    old_at = {position: atom for atom, position in old.items()}
    new_at = {position: atom for atom, position in new.items()}
    sigma = {old_at[i]: new_at[i] for i in range(m.atom_count)}
    assert sorted(sigma) == sorted(sigma.values()), 'not a bijection of the atoms'
    for a, b in sigma.items():
        assert (m.element_of(a), m.charge_of(a), m.implicit_h_of(a)) == \
               (m.element_of(b), m.charge_of(b), m.implicit_h_of(b)), a
        assert {sigma[c] for c in m.neighbors_of(a)} == set(m.neighbors_of(b)), a
    assert any(a != b for a, b in sigma.items()), 'the identity, so the orders were equal'
    assert all(sigma[sigma[a]] == a for a in sigma), 'not an involution, so not a reflection'

    first = write_smiles(m)
    back = read_smiles(first)
    assert write_smiles(back) == first
    assert back.canonical_bytes == m.canonical_bytes
    assert back == m


def test_the_replay_can_detect_a_wrong_answer():
    """Show the comparison failing (ruling F102).

    The test above is a loop that asserts an empty list. If `snapshot()` returned a constant, or the
    keys silently stopped matching, it would pass on every record. So corrupt one expected value and
    require the same comparison to name that key and no other.
    """
    expected = dict(V3_FIXTURES['with_everything']['answers'])
    expected['canonical_order'] = tuple(reversed(expected['canonical_order']))
    got = snapshot(MoleculeContainer.from_bytes(V3_FIXTURES['with_everything']['cold']))
    wrong = sorted(k for k in expected if got[k] != expected[k])
    assert wrong == ['canonical_order']


@mark.parametrize('name,warmth', BOTH)
def test_reserialising_a_v3_buffer_carries_every_payload_byte(name, warmth):
    """A v3 buffer read and written again keeps every v3 payload byte, at its own segment's length.

    This is the concrete pay-off of a design decision that could have gone the other way. Removing
    the derived segments freed v3's `total_len` field, and the four bytes it occupied are exactly what
    `seg_count` needed -- so the table still begins at offset 24 and persistent ids 0-4 mean the same
    thing in both versions. Nothing in a v3 payload has to be repacked to be read now, which is why
    the conversion cannot corrupt a record it misunderstands: there is nothing to misunderstand about
    a payload copied byte for byte.

    The record can nonetheless GROW by one segment, and one class of record does: a v3 buffer states
    its parities in `atom_t.flags`, ingest moves them into SEG_PARITY, and a segment the incoming
    buffer has no room for means a fresh block. So the five v3 payloads are compared at the offsets
    each header names -- a claim that holds whether or not the record grew -- and the growth itself is
    asserted to be exactly the records that state a parity.

    The move is a move and not a copy: ingest clears flags bits 1 and 7, which from version 5 are
    reserved and refused on ingest, so leaving them set would make this molecule's own `to_bytes`
    unreadable by its own `from_bytes`. Those two bits per atom record are therefore the one part of a
    v3 payload allowed to differ, and the atoms segment is compared with them masked while every other
    byte of every record must survive exactly.
    """
    data = V3_FIXTURES[name][warmth]
    mol = MoleculeContainer.from_bytes(data)
    out = mol.to_bytes()
    assert unpack_from('<H', out, 4)[0] == STRUCT_VERSION, 'the record was not upgraded'
    stated = any(mol.parity_of(n) for n in mol.atom_numbers)
    src = _table(data)
    dst = _table(out, unpack_from('<H', out, 20)[0])
    for seg in range(5):
        off, length = src[seg]
        new_off, new_length = dst[seg]
        assert new_length == length, 'persistent segment %d changed length' % seg
        if seg == SEG_ATOMS:
            got = out[new_off:new_off + length]
            assert all(got[i] & RESERVED_FLAGS == 0 for i in range(ATOM_FLAGS, length, ATOM_RECORD)), \
                'ingest left a reserved flag bit set, so this buffer is unreadable by from_bytes'
            assert got == _reserved_bits_masked(data[off:off + length]), \
                'the atoms segment changed outside flags bits 1 and 7'
        else:
            assert out[new_off:new_off + length] == data[off:off + length], \
                'persistent segment %d changed content' % seg
        if not stated:
            assert src[seg] == dst[seg], 'segment %d table entry moved without a parity segment' % seg
    assert bool(_parity_bytes(mol)) == stated, 'the flag parities did not reach SEG_PARITY'


@mark.parametrize('name', ALL)
def test_a_warm_v3_buffer_relabelled_v4_is_refused(name):
    """The version dispatch must really branch, not merely record which number it saw.

    A warm v3 buffer's bytes 20-21 held v3's `total_len`; on the v4 arm they become `seg_count`.
    That value is large (typically hundreds) so the hoisted header-length check fires:
    `'packed molecule is too short for its N-entry segment table'`.  If the dispatch accepted,
    a v3 buffer would be read under v4's rules and the tests above would measure one code path
    while claiming to measure two.
    """
    data = bytearray(V3_FIXTURES[name]['warm'])
    data[4] = 4
    with raises(ValueError, match=r'too short for its'):
        MoleculeContainer.from_bytes(bytes(data))


@mark.parametrize('version', (0, 1, 2, 7, 255))
def test_an_unrecognised_version_is_refused(version):
    """Neither older-than-v3 nor newer-than-v6 may be guessed at.

    Forward rejection matters as much as backward acceptance: a v7 buffer handed to this build must
    fail rather than be read as v6 and quietly misinterpreted, because the whole reason a version
    field exists is that the reader cannot know what it does not know.
    """
    data = bytearray(V3_FIXTURES['with_everything']['cold'])
    data[4:6] = version.to_bytes(2, 'little')
    with raises(ValueError, match='version'):
        MoleculeContainer.from_bytes(bytes(data))


def _blind_offsets(data):
    """Payload byte positions a flip cannot change an answer at, and why each one is there.

    Three kinds, all of them consequences of the reader being STRICTER than the format:

      * DERIVED FIELDS INSIDE A PERSISTENT SEGMENT. `atom_t` carries degree, heteroatom count, ring
        sizes and ring counts, and `halfedge_t` carries the in-ring flag; all of them are recomputed
        by `rebuild_derived` on every read, so a forged value is overwritten rather than trusted.
        That is the desired behaviour -- a packed buffer's derived fields are not evidence -- and
        finding these positions blind is a check on it, not a gap in it.
      * `map_number`'s low byte. In range, and no answer in the frozen set reads a map number,
        because none of these records carries one. A record that did would detect it.
      * ALIGNMENT PADDING between segments and after the last one. Not data at anybody's request.

    Everything else -- every element, charge, hydrogen count, isotope, stable id, atom flag byte,
    CSR index, bond order, wedge and stereo group byte -- must be load-bearing.

    THE HALFEDGE FLAG BYTE IS NOT BLIND, and its absence from this list is the sharpest single piece
    of evidence that storing order 4 tightened validation rather than loosening it.
    `structure_from_bytes` requires `HE_AROMATIC` to agree with `order == 4` and rejects every other
    flag bit, so flipping byte 6 of a half-edge either sets HE_AROMATIC on an order-1 bond (refused)
    or sets a reserved bit (refused). Byte 7, the high half of the flags field, is refused for the
    second reason alone.
    """
    table = _table(data)
    atoms_off, atoms_len = table[0]
    ptr_off, ptr_len = table[1]
    sg_off, sg_len = table[4]

    blind = set()
    for i in range(atoms_len // ATOM_RECORD):
        base = atoms_off + i * ATOM_RECORD
        blind.add(base + 6)                                   # map_number, low byte
        blind.update(range(base + 12, base + 20))             # degree..ring_counts, recomputed
    # csr_ptr holds atom_count + 1 entries of 4 bytes and is then padded up to 8
    used = 4 * (unpack_from('<I', data, 8)[0] + 1)
    blind.update(range(ptr_off + used, ptr_off + ptr_len))
    # SEG_STEREO_GROUPS is one byte per atom, padded up to 8
    if sg_len:
        blind.update(range(sg_off + unpack_from('<I', data, 8)[0], sg_off + sg_len))
    # the trailing alignment of whichever persistent segment is laid out last
    last = max((off + ln for off, ln in table[:5] if ln), default=0)
    blind.update(range(last, unpack_from('<I', data, 16)[0]))
    return blind


@mark.parametrize('name', ALL)
def test_a_single_flipped_payload_byte_does_not_go_unnoticed(name):
    """Sweep every persistent payload byte. Each flip must be rejected or must change an answer.

    This is what turns the replay above from a corpus into evidence (ruling F102): it establishes
    that the frozen answers actually DEPEND on the bytes, so a reader that ignored the buffer and
    answered from thin air could not pass. Anything blind is enumerated and justified by
    `_blind_offsets`, and the count of blind positions is required to stay a minority of the payload
    -- if a change made half the record inert, that assertion is where it surfaces.
    """
    data = V3_FIXTURES[name]['cold']
    expected = V3_FIXTURES[name]['answers']
    # persistent_len is v3 header bytes 16..19; the derived tail past it is not serialised state
    payload_end = unpack_from('<I', data, 16)[0]
    blind = _blind_offsets(data)
    rejected = changed = 0
    unnoticed = []
    for i in range(V3_HEADER_LEN, payload_end):
        forged = bytearray(data)
        forged[i] ^= 0xFF
        try:
            got = snapshot(MoleculeContainer.from_bytes(bytes(forged)))
        except Exception:
            rejected += 1
            continue
        if any(got[k] != expected[k] for k in expected):
            changed += 1
        elif i not in blind:
            unnoticed.append(i)
    assert not unnoticed, 'flipping payload byte(s) %r changed nothing and is not a known-blind ' \
                          'position' % (unnoticed,)
    assert rejected, 'not one forged byte was rejected, so validation is doing nothing'
    assert changed, 'not one forged byte changed an answer, so the answers ignore the buffer'
    assert len(blind) * 2 < payload_end - V3_HEADER_LEN, \
        'most of the record is inert: %d of %d payload bytes' % (len(blind),
                                                                 payload_end - V3_HEADER_LEN)
