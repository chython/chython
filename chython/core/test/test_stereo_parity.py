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
from itertools import permutations

import pytest

from chython.core import MoleculeContainer, _core, WEDGE_DOWN, read_smiles


def _atom_records(data):
    """Where the atom records begin and end in a packed buffer, read out of its segment table.

    The forging tests below hunt for a flags byte by scanning bytes, so they must know where the
    atom payload starts.  That is NOT a constant: `seg_count` is data, the header is
    `24 + 8 * seg_count` bytes, and a molecule with no coordinates and no S-groups stops its table
    at three entries -- a 48-byte header, not 128.  Scanning from a hard-coded 128 would start
    somewhere inside the third atom record and forge the wrong atom, which is a test that passes
    for the wrong reason rather than one that fails.
    """
    seg_count = struct.unpack_from('<H', data, 20)[0]
    assert 24 + 8 * seg_count <= len(data), 'segment table does not fit the buffer'
    offset, length = struct.unpack_from('<II', data, 24)   # SEG_ATOMS is id 0, the first entry
    return offset, offset + length


ATOM_RECORD = 24     # sizeof(atom_t)
ATOM_FLAGS = 3       # its `flags` byte, at byte 3 of a packed record


def _forged_v4_legacy(m, index):
    """`m`'s buffer relabelled version 4, with atom `index`'s flags forged to bit 1 alone.

    A version-4 buffer states a parity in `atom_t.flags`, bit 1 the value and bit 7 the configured
    bit.  Bit 1 alone is a writer that spent bit 1 on a drawing flag instead, which no writer in this
    build produces -- a wedge writes neither parity bit -- so the pattern has to be forged.  The
    version byte is forged with it because ingest reads the flags only for a buffer older than
    STRUCT_VERSION, and `m` must therefore carry no parity segment of its own: a table of three
    entries is a layout version 4 shares, a tenth entry is not.
    """
    data = bytearray(m.to_bytes())
    assert struct.unpack_from('<H', data, 20)[0] == 3, 'this buffer has a segment version 4 lacks'
    offset, end = _atom_records(data)
    at = offset + ATOM_RECORD * index + ATOM_FLAGS
    assert at < end, 'atom %d is past the atom records' % index
    data[at] = (data[at] | 0x02) & ~0x80
    struct.pack_into('<H', data, 4, 4)
    return bytes(data)


# AN ANCHOR'S HYDROGEN COUNT IS STATED, A SUBSTITUENT'S IS USUALLY NOT, and the difference in the
# fixtures below is deliberate rather than sloppy.  `_stereo.pxi` refuses a unit whose ANCHOR has an
# unknown count -- three heavy neighbours plus an unknown hydrogen is a stereocentre or is not, and
# the missing number is precisely which -- and asks nothing about the substituents.  So the oxime
# nitrogens and the sulfoxide sulfur below say `0`, which is what `=N-` and `>S=O` actually carry;
# spare methyls and the OH oxygen keep `None`, because the fixture genuinely does not care and
# `None` is now how that is spelled.  These `0`s USED TO BE `None` and passed on the arena's stored
# zero; when the default moved to `H_UNKNOWN` the anchors lost their direction count and thirteen
# tests here went from asserting on a parity to asserting on an empty unit list.


def _chiral_methane():
    """CFClBr with parity=1: refs are (F, Cl, Br, None) in CSR ascending order.

    The core never derives an implicit hydrogen count (test_stereo_units.py:46-49), so the
    carbon's one implicit H is stated explicitly.  Without implicit_h=1 the C would have only
    three directions and no stereo unit would be perceived.  F, Cl and Br get implicit_h=None
    (the default) so h_pinned is not set on them for no reason.
    """
    m = MoleculeContainer()
    with m.edit():
        sids = [m.add_atom(e, implicit_h=(1 if e == 'C' else None))
                for e in ('C', 'F', 'Cl', 'Br')]
        for s in sids[1:]:
            m.add_bond(sids[0], s, 1)
        m.set_parity(sids[0], 1)
    return m, sids


def _but_2_ene():
    """But-2-ene with parity=1: refs are (sids[0], None, sids[3], None), unnamed_mask=0b1010.

    The double-bond carbons each carry one implicit H; their methyl neighbours need no explicit H
    (they have only one heavy neighbour, the double-bond carbon, so would get three implicit H
    automatically -- but here we leave implicit_h=None to avoid adding tetrahedral units for them).
    Atom order: C0(methyl), C1(=CH-, anchor), C2(=CH-), C3(methyl).
    Bonds: C0-C1 single, C1=C2 double, C2-C3 single.
    Anchor is sids[1] (lower double-bond terminal).
    """
    m = MoleculeContainer()
    with m.edit():
        sids = [m.add_atom('C', implicit_h=h)
                for h in (None, 1, 1, None)]
        m.add_bond(sids[0], sids[1], 1)
        m.add_bond(sids[1], sids[2], 2)
        m.add_bond(sids[2], sids[3], 1)
        m.set_parity(sids[1], 1)
    return m, sids


def _dichlorobut_2_ene():
    """2,3-Dichlorobut-2-ene with parity=1: all four refs are named (unnamed_mask=0).

    Atom order: C0(CCl=), Cl1, C2(=CCl), Cl3, C4(methyl on C0 side), C5(methyl on C2 side).
    Bonds: C0-Cl1, C0=C2, C0-C4, C2-Cl3, C2-C5.
    Anchor is sids[0] (lower double-bond terminal).
    refs = (sids[1], sids[4], sids[3], sids[5]) — CSR ascending order within each pair.
    """
    m = MoleculeContainer()
    with m.edit():
        sids = [m.add_atom(e) for e in ('C', 'Cl', 'C', 'Cl', 'C', 'C')]
        m.add_bond(sids[0], sids[1], 1)   # C0-Cl1
        m.add_bond(sids[0], sids[2], 2)   # C0=C2
        m.add_bond(sids[0], sids[4], 1)   # C0-C4(methyl)
        m.add_bond(sids[2], sids[3], 1)   # C2-Cl3
        m.add_bond(sids[2], sids[5], 1)   # C2-C5(methyl)
        m.set_parity(sids[0], 1)
    return m, sids


def _acetaldoxime():
    """Acetaldoxime CH3-CH=N-OH with parity=1: unnamed_mask=0b0010.

    Slot 1 is the carbon's implicit H (real unnamed direction, mask bit 1 set).
    Slot 3 is the nitrogen lone pair (NOT a direction, mask bit 3 clear = PINNED).
    Atom order: C0(methyl), C1(=CH-), N2(=N-), O3(OH).
    Bonds: C0-C1 single, C1=N2 double, N2-O3 single.
    Anchor is sids[1].  refs = (sids[0], None, sids[3], None), unnamed_mask=0b0010.
    """
    m = MoleculeContainer()
    with m.edit():
        sids = [m.add_atom(e, implicit_h=h)
                for e, h in (('C', None), ('C', 1), ('N', 0), ('O', None))]
        m.add_bond(sids[0], sids[1], 1)
        m.add_bond(sids[1], sids[2], 2)
        m.add_bond(sids[2], sids[3], 1)
        m.set_parity(sids[1], 1)
    return m, sids


# ---------------------------------------------------------------------------
# STANDING REQUIREMENT for the bond-kind fixtures above and below.
#
# A bond-kind fixture set must VARY THE NUMBER OF UNNAMED SLOTS PER PAIR, because that is the
# dimension `translate_stereo` branches on.  All three cases must be present:
#
#   0 unnamed slots            -- _dichlorobut_2_ene (both terminals disubstituted)
#   1 per pair (2 in total)    -- _but_2_ene, _acetaldoxime, _acetaldoxime_on
#   exactly 1 in total         -- _chlorobut_2_ene / _chlorobut_2_ene_reversed
#
# The third is the commonest real E/Z shape there is -- one terminal disubstituted, the other
# bearing a hydrogen -- and it was absent from this file until fix round 4.  Two of this epic's
# defects have now hidden behind one-sided fixtures: first behind which pair held the wildcard,
# then behind how many unnamed slots each pair held.  A fixture set that is symmetric in the
# dimension the code branches on cannot see a bug in that branch.
# ---------------------------------------------------------------------------

def _chlorobut_2_ene():
    """(Z)-2-Chlorobut-2-ene CH3-C(Cl)=CH-CH3, anchored at the CHLORINATED terminal.

    The one bond-kind shape the suite lacked until fix round 4: exactly ONE SU_NO_REF slot in
    the whole unit, so the two pairs hold DIFFERENT numbers of unnamed slots.  Pair 0 (the
    anchor's) is fully named; pair 1 holds the other terminal's methyl and its implicit H.

    Atom order: C0(=C(Cl)-, anchor), Cl1, C2(methyl on C0), C3(=CH-), C4(methyl on C3).
    Bonds: C0-Cl1, C0-C2, C0=C3, C3-C4.  C3 carries the one implicit H.
    refs = (Cl1, C2, C4, None), unnamed_mask = 0b1000 (slot 3 = C3's implicit H, a real
    unnamed direction; there is no pinned slot in this unit).
    """
    m = MoleculeContainer()
    with m.edit():
        sids = [m.add_atom(e, implicit_h=h)
                for e, h in (('C', None), ('Cl', None), ('C', None), ('C', 1), ('C', None))]
        m.add_bond(sids[0], sids[1], 1)   # C0-Cl1
        m.add_bond(sids[0], sids[2], 1)   # C0-C2(methyl)
        m.add_bond(sids[0], sids[3], 2)   # C0=C3
        m.add_bond(sids[3], sids[4], 1)   # C3-C4(methyl)
        m.set_parity(sids[0], 1)
    return m, sids


def _chlorobut_2_ene_reversed():
    """The same (Z)-2-chlorobut-2-ene spelled from the CH terminal, which is then the anchor.

    Same one-unnamed-slot shape as `_chlorobut_2_ene`, with the unnamed slot in the OTHER pair:
    pair 0 (the anchor's) holds the methyl and the implicit H, pair 1 is fully named.  Both
    spellings are needed for the same reason round 2 needed both oxime spellings -- the pair
    that holds the wildcard must not be a constant of the fixture set.

    Atom order: C0(=CH-, anchor), C1(methyl on C0), C2(=C(Cl)-), Cl3, C4(methyl on C2).
    Bonds: C0-C1, C0=C2, C2-Cl3, C2-C4.  C0 carries the one implicit H.
    refs = (C1, None, Cl3, C4), unnamed_mask = 0b0010 (slot 1 = C0's implicit H).
    """
    m = MoleculeContainer()
    with m.edit():
        sids = [m.add_atom(e, implicit_h=h)
                for e, h in (('C', 1), ('C', None), ('C', None), ('Cl', None), ('C', None))]
        m.add_bond(sids[0], sids[1], 1)   # C0-C1(methyl)
        m.add_bond(sids[0], sids[2], 2)   # C0=C2
        m.add_bond(sids[2], sids[3], 1)   # C2-Cl3
        m.add_bond(sids[2], sids[4], 1)   # C2-C4(methyl)
        m.set_parity(sids[0], 1)
    return m, sids


# ---------------------------------------------------------------------------
# Atom-kind (SU_TETRA) tests -- these are the original 13 tests, trimmed and fixed.
# ---------------------------------------------------------------------------

def test_identity_order_returns_the_stored_parity():
    m, sids = _chiral_methane()
    assert m.translate_stereo(sids[0], (sids[1], sids[2], sids[3], None)) == 1


def test_one_swap_flips_the_parity():
    m, sids = _chiral_methane()
    assert m.translate_stereo(sids[0], (sids[2], sids[1], sids[3], None)) == 2


def test_two_swaps_restore_the_parity():
    m, sids = _chiral_methane()
    assert m.translate_stereo(sids[0], (sids[2], sids[3], sids[1], None)) == 1


def test_every_permutation_matches_its_inversion_count():
    m, sids = _chiral_methane()
    refs = (sids[1], sids[2], sids[3], None)
    for perm in permutations(range(4)):
        inversions = sum(1 for i in range(4) for j in range(i + 1, 4) if perm[i] > perm[j])
        want = tuple(refs[k] for k in perm)
        expected = 1 if inversions % 2 == 0 else 2
        assert m.translate_stereo(sids[0], want) == expected, perm


def test_unset_parity_translates_to_unset():
    m = MoleculeContainer()
    with m.edit():
        sids = [m.add_atom(e, implicit_h=(1 if e == 'C' else None))
                for e in ('C', 'F', 'Cl', 'Br')]
        for s in sids[1:]:
            m.add_bond(sids[0], s, 1)
    # Check identity (pp=0) and a non-identity order (pp=1): unset parity must stay 0.
    # The early-return guard in translate_parity is the contract; the pp!=0 case is the seatbelt.
    for want in [
        (sids[1], sids[2], sids[3], None),   # identity: pp = 0
        (sids[2], sids[1], sids[3], None),   # one swap:  pp = 1
        (sids[3], sids[1], sids[2], None),   # odd perm
    ]:
        assert m.translate_stereo(sids[0], want) == 0, want


def test_implicit_hydrogen_sorts_last_in_the_stored_refs():
    m, sids = _chiral_methane()
    assert m.unit_of(sids[0])['refs'][3] is None


def test_set_parity_round_trips_through_parity_of():
    m, sids = _chiral_methane()
    assert m.unit_of(sids[0])['parity'] == 1
    assert m.parity_of(sids[0]) == 1
    with m.edit():
        m.set_parity(sids[0], 2)
    assert m.unit_of(sids[0])['parity'] == 2
    assert m.parity_of(sids[0]) == 2


def test_parity_is_unset_until_stated():
    m = MoleculeContainer()
    with m.edit():
        sids = [m.add_atom(e, implicit_h=(1 if e == 'C' else None))
                for e in ('C', 'F', 'Cl', 'Br')]
        for s in sids[1:]:
            m.add_bond(sids[0], s, 1)
    assert m.parity_of(sids[0]) == 0, 'a drawn centre with no wedge is unset, not even'


def test_clearing_a_parity_returns_it_to_unset():
    m, sids = _chiral_methane()
    with m.edit():
        m.set_parity(sids[0], 0)
    assert m.parity_of(sids[0]) == 0


def test_legacy_set_stereo_still_reads_back_as_a_bool():
    # set_stereo / stereo_of is exercised on its own terms; no external dependency.
    m, sids = _chiral_methane()
    with m.edit():
        m.set_stereo(sids[0], True)
    assert m.stereo_of(sids[0]) is True
    assert m.parity_of(sids[0]) == 2, 'a legacy write must land configured, not unset'
    with m.edit():
        m.set_stereo(sids[0], False)
    assert m.stereo_of(sids[0]) is False
    assert m.parity_of(sids[0]) == 1


def test_set_parity_rejects_out_of_range():
    m, sids = _chiral_methane()
    with pytest.raises(ValueError):
        with m.edit():
            m.set_parity(sids[0], 3)


def test_translate_rejects_a_non_permutation():
    m, sids = _chiral_methane()
    with pytest.raises(ValueError):
        m.translate_stereo(sids[0], (sids[1], sids[1], sids[3], None))


def test_translate_rejects_an_unanchored_atom():
    m, sids = _chiral_methane()
    with pytest.raises(KeyError):
        m.translate_stereo(sids[1], (sids[0], None, None, None))


# ---------------------------------------------------------------------------
# Guard tests (M7, M8, M9) -- unit record write test (Item 3).
# ---------------------------------------------------------------------------

def test_translate_rejects_wrong_none_count():
    """M7 guard: a caller who passes two Nones for a unit with one unnamed direction raises."""
    m, sids = _chiral_methane()
    # refs = (F, Cl, Br, None): only one real unnamed direction
    with pytest.raises(ValueError):
        m.translate_stereo(sids[0], (sids[1], sids[2], None, None))


def test_translate_rejects_wrong_order_length():
    """M8 guard: order must have exactly n_refs elements."""
    m, sids = _chiral_methane()
    with pytest.raises(ValueError):
        m.translate_stereo(sids[0], (sids[1], sids[2], sids[3]))


def test_translate_rejects_foreign_stable_id():
    """M9 guard: a stable id not in the molecule raises ValueError (not KeyError)."""
    m, sids = _chiral_methane()
    with pytest.raises(ValueError):
        m.translate_stereo(sids[0], (9999, sids[2], sids[3], None))


def test_translate_stereo_does_not_write_unit_record():
    """Item 3: translate_stereo must not mutate u.parity; _stereo_emit is the only writer.

    _unit_parity_raw_probe reads u.parity directly from the derived segment.  It must be 0
    before and after any translate_stereo call.  Re-introducing the u.parity write in
    translate_stereo makes this assertion fail.
    """
    m, sids = _chiral_methane()
    assert _core._unit_parity_raw_probe(m, sids[0]) == 0, 'perception must not write u.parity'
    result = m.translate_stereo(sids[0], (sids[1], sids[2], sids[3], None))
    assert result == 1
    assert _core._unit_parity_raw_probe(m, sids[0]) == 0, 'translate_stereo must not write u.parity'


# ---------------------------------------------------------------------------
# Bit-1-alone in an older buffer, and the wedge that discriminates the two readings.
# ---------------------------------------------------------------------------

def _forgeable_chiral_methane(wedge):
    """CHFClBr with NO parity, optionally wedged from C to F, and C's stable id.

    No parity, so the buffer stops its table at three entries and `_forged_v4_legacy` can relabel it.
    The wedge is what the discriminator reads, and the atom it names is C: the narrow end.
    """
    m = MoleculeContainer()
    with m.edit():
        c = m.add_atom('C', implicit_h=1)
        f = m.add_atom('F')
        cl = m.add_atom('Cl')
        br = m.add_atom('Br')
        m.add_bond(c, f, 1)
        m.add_bond(c, cl, 1)
        m.add_bond(c, br, 1)
        if wedge:
            m.set_wedge(c, f, WEDGE_DOWN)
    return m, c


def test_legacy_buffer_no_wedge_normalises_to_parity_2():
    """Bit 1 alone on an atom with no wedge is a stored parity direction: promoted to configured-odd.

    Bit 7 is set, so `parity_of` is 2 (odd) -- an even original is not recoverable from that encoding
    -- and the promoted parity is then adopted into SEG_PARITY like any other one the buffer states.
    """
    m, c = _forgeable_chiral_methane(False)
    m2 = MoleculeContainer.from_bytes(_forged_v4_legacy(m, m.index_of(c)))
    assert m2.parity_of(c) == 2, 'bit-1-alone without wedge must normalise to parity 2'
    assert _core._parity_bytes(m2)[m2.index_of(c)] == 2, 'and the promoted parity must reach the segment'


def test_legacy_buffer_wedge_narrow_end_normalises_to_no_parity():
    """Bit 1 alone on a wedge's narrow end was set by the wedge: cleared rather than promoted.

    The atom is found through `halfedge_t.wedge` in SEG_CSR_EDGE, part of the persistent prefix, so
    the discriminator is readable before `rebuild_derived`.  Nothing is configured, so the record
    states no parity and adoption lays out no segment for it.
    """
    m, c = _forgeable_chiral_methane(True)
    m2 = MoleculeContainer.from_bytes(_forged_v4_legacy(m, m.index_of(c)))
    assert m2.parity_of(c) == 0, 'narrow-end wedge atom must have parity_of == 0 after normalisation'
    assert m2.stereo_of(c) is False, 'narrow-end wedge atom must have stereo_of False after normalisation'
    assert _core._parity_bytes(m2) == b'', 'nothing is configured, so nothing is adopted'


def test_a_legacy_buffer_s_conformer_survives_the_adoption():
    """The adoption reallocates, so every persistent payload is copied by hand -- including the 3D one.

    No frozen fixture carries both a conformer and a stated parity, so the buffer is forged: a
    current-version buffer's tenth table entry is dropped, which shortens the header by 8 bytes and
    moves every payload offset with it.  A dropped `SEG_CONFORMERS` copy would leave the segment
    allocated and zeroed, so the assertion is on the coordinates and not on `has_3d`.
    """
    m = read_smiles('C[C@H](N)C(=O)O')
    with m.edit():
        for i, n in enumerate(m.atom_numbers, 1):
            m.set_xyz(n, 0.1 * i, 0.2 * i, 0.3 * i)
    assert m.has_3d is True

    data = bytearray(m.to_bytes())
    seg_count = struct.unpack_from('<H', data, 20)[0]
    assert seg_count == 10, 'this buffer does not end its table at the parity segment'
    header_len = 24 + 8 * seg_count
    par_off, par_len = struct.unpack_from('<II', data, 24 + 8 * 9)
    assert par_len, 'the fixture states no parity'

    out = bytearray(data[:header_len - 8] + data[header_len:par_off])
    struct.pack_into('<H', out, 4, 4)                 # version
    struct.pack_into('<H', out, 20, 9)                # seg_count
    struct.pack_into('<I', out, 16, len(out))         # persistent_len -- from_bytes checks it
    for i in range(9):
        offset, length = struct.unpack_from('<II', data, 24 + 8 * i)
        struct.pack_into('<II', out, 24 + 8 * i, offset - 8 if length else 0, length)
    at = struct.unpack_from('<II', out, 24)[0] + ATOM_RECORD * 1 + ATOM_FLAGS
    out[at] |= 0x82                                   # configured, odd

    # AND THE OLD CONFORMER RECORD IS FOUR WORDS WIDE.  A version-4 buffer's record is
    # `CONFORMER_RECORD_V5` bytes where this build's is one word, so the table is widened back with
    # the three dropped words left zero -- otherwise the ingest length check, which is an equality,
    # refuses the forgery before the adoption under test can run.  The segment is last in the
    # truncated buffer, so widening it only grows the tail.
    c_off, c_len = struct.unpack_from('<II', out, 24 + 8 * 8)
    assert c_len, 'the fixture carries no conformer'
    assert struct.unpack_from('<I', out, c_off)[0] == 1, 'one model, so one record to widen'
    # The xyz block is taken at its OWN length rather than to the end of the buffer: the segment's
    # length is padded to an 8-boundary, and carrying that padding across the widening would produce
    # a length the old stride's arithmetic does not come to.
    xyz_len = struct.unpack_from('<I', out, 8)[0] * 12
    widened = (out[c_off:c_off + 8] + out[c_off + 8:c_off + 12] + bytes(12)
               + out[c_off + 12:c_off + 12 + xyz_len])
    while len(widened) % 8:
        widened += b'\0'
    out = out[:c_off] + widened
    struct.pack_into('<I', out, 24 + 8 * 8 + 4, len(widened))
    struct.pack_into('<I', out, 16, len(out))

    back = MoleculeContainer.from_bytes(bytes(out))
    assert _core._parity_bytes(back) == _core._parity_bytes(m)
    assert back.xyz_of(1) == m.xyz_of(1)
    assert back.xyz_of(6) == m.xyz_of(6) != (0.0, 0.0, 0.0)
    assert str(back) == str(m)


# ---------------------------------------------------------------------------
# Bond-kind (SU_BOND) tests -- Items 1 and 2.
# ---------------------------------------------------------------------------

def test_but_2_ene_unit_structure():
    """Verify the but-2-ene unit structure before testing translation."""
    m, sids = _but_2_ene()
    units = m.stereo_units()
    assert len(units) == 1
    u = units[0]
    assert u['kind'] == 1 and u['anchor'] == sids[1]
    assert u['refs'] == (sids[0], None, sids[3], None)
    assert u['unnamed_mask'] == 0b1010


def test_but_2_ene_identity_returns_stored_parity():
    m, sids = _but_2_ene()
    assert m.translate_stereo(sids[1], (sids[0], None, sids[3], None)) == 1


def test_but_2_ene_pair_exchange_does_not_flip():
    """The wholesale pair exchange (C1,None,C0,None) is even -- parity must be unflipped."""
    m, sids = _but_2_ene()
    assert m.translate_stereo(sids[1], (sids[3], None, sids[0], None)) == 1


def test_but_2_ene_within_pair_0_swap_flips():
    m, sids = _but_2_ene()
    assert m.translate_stereo(sids[1], (None, sids[0], sids[3], None)) == 2


def test_but_2_ene_within_pair_1_swap_flips():
    m, sids = _but_2_ene()
    assert m.translate_stereo(sids[1], (sids[0], None, None, sids[3])) == 2


def test_but_2_ene_cross_pair_raises():
    """An order that mixes pair 0 and pair 1 directions must raise ValueError."""
    m, sids = _but_2_ene()
    with pytest.raises(ValueError):
        m.translate_stereo(sids[1], (sids[0], sids[3], None, None))


def test_dichlorobut_2_ene_unit_structure():
    """Verify the 2,3-dichlorobut-2-ene unit structure."""
    m, sids = _dichlorobut_2_ene()
    units = m.stereo_units()
    assert len(units) == 1
    u = units[0]
    assert u['kind'] == 1 and u['anchor'] == sids[0]
    assert u['unnamed_mask'] == 0


def test_dichlorobut_2_ene_identity_returns_stored_parity():
    m, sids = _dichlorobut_2_ene()
    u = m.unit_of(sids[0])
    refs = u['refs']
    assert m.translate_stereo(sids[0], refs) == 1


def test_dichlorobut_2_ene_pair_exchange_does_not_flip():
    m, sids = _dichlorobut_2_ene()
    u = m.unit_of(sids[0])
    refs = u['refs']  # (sids[1], sids[4], sids[3], sids[5]) or similar
    # Pair exchange: swap pair 0 and pair 1 wholesale
    exchanged = refs[2:4] + refs[0:2]
    assert m.translate_stereo(sids[0], exchanged) == 1


def test_dichlorobut_2_ene_within_pair_swap_flips():
    m, sids = _dichlorobut_2_ene()
    u = m.unit_of(sids[0])
    refs = u['refs']
    # Swap within pair 0
    swapped = (refs[1], refs[0], refs[2], refs[3])
    assert m.translate_stereo(sids[0], swapped) == 2


def test_dichlorobut_2_ene_cross_pair_raises():
    m, sids = _dichlorobut_2_ene()
    u = m.unit_of(sids[0])
    refs = u['refs']
    # Cross-pair: refs[0] from pair 0 and refs[2] from pair 1 mixed into first two positions
    mixed = (refs[0], refs[2], refs[1], refs[3])
    with pytest.raises(ValueError):
        m.translate_stereo(sids[0], mixed)


def test_acetaldoxime_unit_structure():
    """Verify acetaldoxime's mask: bit 1 set (real direction), bit 3 clear (pinned)."""
    m, sids = _acetaldoxime()
    units = m.stereo_units()
    assert len(units) == 1
    u = units[0]
    assert u['kind'] == 1 and u['anchor'] == sids[1]
    assert u['refs'] == (sids[0], None, sids[3], None)
    assert u['unnamed_mask'] == 0b0010, 'slot 1 is a real direction; slot 3 is pinned'
    assert u['unnamed_mask'] & 0b1000 == 0, 'bit 3 must be clear: nitrogen lone pair is pinned'


def test_acetaldoxime_identity_returns_stored_parity():
    m, sids = _acetaldoxime()
    assert m.translate_stereo(sids[1], (sids[0], None, sids[3], None)) == 1


def test_acetaldoxime_within_pair_swap_flips():
    """Slot 1 is a real unnamed direction and can be swapped within pair 0."""
    m, sids = _acetaldoxime()
    assert m.translate_stereo(sids[1], (None, sids[0], sids[3], None)) == 2


def test_acetaldoxime_pair_exchange_does_not_flip():
    """Wholesale pair exchange is even; parity must be unflipped.

    CH3-C=N-O: pair exchange gives (O, None, CH3, None).  Pair exchange is (0 2)(1 3) which
    is two transpositions = even, so parity stays 1.
    """
    m, sids = _acetaldoxime()
    assert m.translate_stereo(sids[1], (sids[3], None, sids[0], None)) == 1


def test_acetaldoxime_cross_pair_raises():
    """An order mixing named atoms from different pairs must raise ValueError."""
    m, sids = _acetaldoxime()
    with pytest.raises(ValueError):
        m.translate_stereo(sids[1], (sids[0], sids[3], None, None))


def test_acetaldoxime_moving_pinned_slot_raises():
    """Any order that attempts to swap the pinned other-slot (position 3) must raise ValueError.

    Slot 3 is the N lone pair (pinned: unnamed_mask bit 3 clear).  Moving it to any non-None
    position, or swapping within its pair, is forbidden by Ruling F55.
    """
    m, sids = _acetaldoxime()
    with pytest.raises(ValueError):
        # want[2]=None, want[3]=O: O goes to pair-1 in swapped order, but pair-1's other-slot
        # is the N lone pair (pinned) -- frozen within-pair order raises.
        m.translate_stereo(sids[1], (sids[0], None, None, sids[3]))


# ---------------------------------------------------------------------------
# Acetaldoxime O-N=C-CH3 (reversed spelling) -- Fix Round 2, Item 1.
# The pinned slot (lone pair) is at refs index 1 (within P0), not index 3.
# ---------------------------------------------------------------------------

def _acetaldoxime_on():
    """O-N=C-CH3 acetaldoxime with parity=1: unnamed_mask=0b1000.

    Atom order: O(0), N(1), C(2), CH3(3).  Bonds: O-N single, N=C double, C-CH3 single.
    N is the anchor (lower double-bond terminal by atom number).
    P0 = N's substituents: (O=sids[0], None=lone_pair_PINNED).
    P1 = C's substituents: (CH3=sids[3], None=H_implicit_real_direction).
    refs = (sids[0], None, sids[3], None), unnamed_mask=0b1000.
    Bit 3 SET (slot 3 = H_implicit, real unnamed direction).
    Bit 1 CLEAR (slot 1 = N lone pair, PINNED -- not a direction, Ruling F55).

    This is the same compound as _acetaldoxime() spelled from the O side.  In round 1 the
    algorithm returned a silently wrong parity for this spelling's pair exchange because the
    absolute-index pin (`perm[j]=j`) collided with per-pair placement.  Fix Round 2 replaces
    the entire bond-kind path with a direct pair decomposition (Ruling F55 / F56) that handles
    both spellings correctly.
    """
    m = MoleculeContainer()
    with m.edit():
        sids = [m.add_atom(e, implicit_h=h)
                for e, h in (('O', None), ('N', 0), ('C', 1), ('C', None))]
        m.add_bond(sids[0], sids[1], 1)
        m.add_bond(sids[1], sids[2], 2)
        m.add_bond(sids[2], sids[3], 1)
        m.set_parity(sids[1], 1)
    return m, sids


def test_acetaldoxime_on_unit_structure():
    """Verify the O-N=C-CH3 unit: unnamed_mask=0b1000 (bit 1 clear = N lone pair is PINNED)."""
    m, sids = _acetaldoxime_on()
    units = m.stereo_units()
    assert len(units) == 1
    u = units[0]
    assert u['kind'] == 1 and u['anchor'] == sids[1]
    assert u['refs'] == (sids[0], None, sids[3], None)
    assert u['unnamed_mask'] == 0b1000, 'bit 3 set = H implicit is a real direction; bit 1 clear = lone pair is pinned'
    assert u['unnamed_mask'] & 0b0010 == 0, 'bit 1 must be clear: N lone pair is pinned, not a direction'


def test_acetaldoxime_on_identity_returns_stored_parity():
    m, sids = _acetaldoxime_on()
    assert m.translate_stereo(sids[1], (sids[0], None, sids[3], None)) == 1


def test_acetaldoxime_on_pair_exchange_does_not_flip():
    """Wholesale pair exchange (CH3,None,O,None) is even; parity must stay 1.

    An algorithm counting the exchange as one swap answers 2 here.
    """
    m, sids = _acetaldoxime_on()
    assert m.translate_stereo(sids[1], (sids[3], None, sids[0], None)) == 1


def test_acetaldoxime_on_within_pair_1_swap_flips():
    """Swap within pair 1 (CH3 / H_implicit): one transposition flips parity."""
    m, sids = _acetaldoxime_on()
    assert m.translate_stereo(sids[1], (sids[0], None, None, sids[3])) == 2


def test_acetaldoxime_on_within_pair_0_pinned_swap_raises():
    """Pair 0's other-slot is a pinned N lone pair; swapping within pair 0 must raise ValueError."""
    m, sids = _acetaldoxime_on()
    with pytest.raises(ValueError):
        m.translate_stereo(sids[1], (None, sids[0], sids[3], None))


def test_acetaldoxime_on_cross_pair_raises():
    """An order mixing named atoms from different pairs must raise ValueError."""
    m, sids = _acetaldoxime_on()
    with pytest.raises(ValueError):
        m.translate_stereo(sids[1], (sids[0], sids[3], None, None))


# ---------------------------------------------------------------------------
# Fix Round 3 — phase-1 validation on both kinds (Item 1, N7).
#
# Before Fix Round 3 the uniqueness and None-count checks sat inside the
# atom-kind arm only; bond-kind orders with duplicated atoms or wrong None
# counts were accepted and returned a parity.  The tests below verify that
# phase 1 fires for bond kinds in the same way it does for atom kinds, and
# thus that both paths agree on what a valid permutation is.
#
# WHICH BRANCH EACH TEST BELOW ACTUALLY PINS (corrected in fix round 4).
# Three branches are the SOLE reason their test's order is rejected, so a
# bare pytest.raises pins them -- disable the branch and only that test fails:
#
#   uniqueness  -- test_dichlorobut2ene_duplicate_pair0_named_atom_raises,
#                  test_dichlorobut2ene_duplicate_other_slot_atom_raises,
#                  test_dichlorobut2ene_duplicate_cross_pair_named_atom_raises,
#                  test_chiral_methane_duplicate_named_atom_raises
#   None count   -- test_chiral_methane_too_many_nones_raises (atom kind only)
#   cross-pair   -- test_dichlorobut_2_ene_cross_pair_raises
#
# Two branches are REACHABLE BUT NEVER THE SOLE VIOLATION, measured over 2,229
# orders on nine bond fixtures: None correspondence, and the None-count check
# on BOND kinds.  Every order that reaches them is also rejected by a later
# check (the other pair's cross-pair test, or None correspondence itself), so
# disabling one changes the message and not the outcome.  Those two tests --
# test_chlorobut_2_ene_none_does_not_correspond_raises and
# test_but_2_ene_too_many_nones_raises -- therefore use match= on purpose.
# That is not a design smell about the branches; it is a fact about check
# ORDERING, and a match= is the only thing that can tell the two apart.
#
# test_but_2_ene_all_none_in_first_half_raises pins nothing on its own: its
# order violates src_pair, cross-pair and None correspondence at once, and it
# is kept for the src_pair message it produces today.  It is documented as
# such rather than advertised as disjoint.
# ---------------------------------------------------------------------------

# -- Bond kind: duplicate named atom (phase 1 uniqueness) --

def test_dichlorobut2ene_duplicate_pair0_named_atom_raises():
    """Both chlorines named twice drops both methyls: phase 1 uniqueness raises.

    refs = (Cl1, C4, Cl3, C5).  Order (Cl1, Cl1, Cl3, Cl3) names Cl1 twice --
    caught by phase 1's uniqueness check on the second Cl1 occurrence.
    Corresponds to the reviewer's (2,2,4,4) example.
    """
    m, sids = _dichlorobut_2_ene()
    with pytest.raises(ValueError):
        # sids[1]=Cl1, sids[3]=Cl3; sids[4]=C4 and sids[5]=C5 are dropped
        m.translate_stereo(sids[0], (sids[1], sids[1], sids[3], sids[3]))


def test_dichlorobut2ene_duplicate_other_slot_atom_raises():
    """A pair's other-slot atom named twice: phase 1 uniqueness raises, and only it.

    refs = (Cl1, C4, Cl3, C5), so C5 is pair 1's other-slot (offset 1) atom.  The order
    (Cl1, C4, C5, C5) names C5 twice and drops Cl3, and uniqueness is its SOLE violation:
    want[0:2] resolves to pair 0 cleanly, both halves stay inside their pair, there is no
    None anywhere and no pinned slot, so with uniqueness disabled the function returns a
    parity for an order that names three of four directions.  That is what makes this test
    pin the branch it names.

    (The neighbouring order (5,5,4,6) = (C4, C4, Cl3, C5) is asserted below as well, but it does
    NOT pin uniqueness: it names no pair's slot-0 atom in want[0:2], so the src_pair raise catches
    it instead.)
    """
    m, sids = _dichlorobut_2_ene()
    with pytest.raises(ValueError):
        # sids[1]=Cl1, sids[4]=C4, sids[5]=C5 twice; sids[3]=Cl3 dropped
        m.translate_stereo(sids[0], (sids[1], sids[4], sids[5], sids[5]))
    with pytest.raises(ValueError):
        # the reviewer's (5,5,4,6): rejected, but by uniqueness OR src_pair
        m.translate_stereo(sids[0], (sids[4], sids[4], sids[3], sids[5]))


def test_dichlorobut2ene_duplicate_cross_pair_named_atom_raises():
    """Cl3 named twice (once per pair position): phase 1 uniqueness raises.

    refs = (Cl1, C4, Cl3, C5).  Order (Cl1, C4, Cl3, Cl3) names Cl3 twice.
    Corresponds to the reviewer's (2,5,4,4) example.
    """
    m, sids = _dichlorobut_2_ene()
    with pytest.raises(ValueError):
        # sids[1]=Cl1, sids[4]=C4, sids[3]=Cl3 twice; sids[5]=C5 dropped
        m.translate_stereo(sids[0], (sids[1], sids[4], sids[3], sids[3]))


def test_but_2_ene_duplicate_named_atom_raises():
    """C0 named twice on but-2-ene: phase 1 uniqueness raises.

    refs = (C0, None, C3, None).  Order (C0, C0, C3, None) names C0 twice.
    Corresponds to the reviewer's (C0,C0,C3,None) example.
    """
    m, sids = _but_2_ene()
    with pytest.raises(ValueError):
        m.translate_stereo(sids[1], (sids[0], sids[0], sids[3], None))


# -- Bond kind: wrong None count (phase 1 None-count) --

def test_but_2_ene_too_many_nones_raises():
    """Three Nones where refs has two: phase 1's None count raises, and says so.

    refs = (C0, None, C3, None) has 2 SU_NO_REF.  Order (C0, None, None, None) has 3 Nones.

    The `match=` is load-bearing, not decoration.  On BOND kinds the None-count check is
    reachable-first but never the SOLE violation -- every bond-kind order with a wrong None
    count also misplaces a None within some pair, so None correspondence would reject it one
    step later.  Measured over 2,229 orders on nine bond fixtures, no bond-kind order exists
    whose only violation is the None count, so no bare pytest.raises can pin this branch on a
    bond kind and matching the message is the only way.  The atom kind is where the check IS
    sole (49 such orders on CFClBr): test_chiral_methane_too_many_nones_raises pins it there.
    """
    m, sids = _but_2_ene()
    with pytest.raises(ValueError, match='None count does not match'):
        m.translate_stereo(sids[1], (sids[0], None, None, None))


# -- Bond kind: src_pair raise (Item 2, N11) --

def test_but_2_ene_all_none_in_first_half_raises():
    """want[0:2] all None -- no named atom -- raises at src_pair assignment.

    refs = (C0, None, C3, None).  Order (None, None, C0, C3) passes phase 1
    (uniqueness: C0 and C3 each appear once; None count: 2 in want, 2 in refs).
    It then reaches phase 2 where the src_pair loop finds no named atom in
    want[0:2] and raises ValueError.

    This test pins NO branch on its own: the same order violates src_pair, the k=0 cross-pair
    check and None correspondence, so disabling any one of the three leaves it passing.  It is
    kept for what it asserts today -- that the order is rejected rather than answered -- and
    not as coverage of the src_pair raise.  No order rejected only by src_pair exists in any
    of this file's fixtures (measured over 2,229 orders on nine of them).
    """
    m, sids = _but_2_ene()
    with pytest.raises(ValueError):
        m.translate_stereo(sids[1], (None, None, sids[0], sids[3]))


# -- Atom kind: same shapes as bond-kind tests above, showing both paths agree --

def test_chiral_methane_duplicate_named_atom_raises():
    """F named twice, Cl dropped: phase 1 uniqueness raises for atom kind.

    refs = (F, Cl, Br, None).  Order (F, F, Br, None) names F twice.
    Shows the atom-kind path is identical to the bond-kind path for phase 1.
    """
    m, sids = _chiral_methane()
    with pytest.raises(ValueError):
        m.translate_stereo(sids[0], (sids[1], sids[1], sids[3], None))


def test_chiral_methane_too_many_nones_raises():
    """Two Nones where refs has one: phase 1 None-count raises for atom kind.

    refs = (F, Cl, Br, None) has 1 SU_NO_REF.  Order (F, Cl, None, None) has 2 Nones.
    """
    m, sids = _chiral_methane()
    with pytest.raises(ValueError):
        m.translate_stereo(sids[0], (sids[1], sids[2], None, None))


# ---------------------------------------------------------------------------
# Fix Round 4 — (Z)-2-chlorobut-2-ene: the one-unnamed-slot bond shape (N13, Ruling F59).
#
# Until this round every bond fixture was either fully named (0 unnamed slots) or symmetric
# (1 per pair).  Nothing covered the commonest E/Z shape of all -- one terminal disubstituted,
# the other bearing a hydrogen -- and that gap hid a live rejection branch behind a comment
# claiming it was unreachable.  See the STANDING REQUIREMENT beside the fixtures.
# ---------------------------------------------------------------------------

def test_chlorobut_2_ene_unit_structure():
    """The Cl-first spelling: one SU_NO_REF in the unit, in pair 1, and it is a real direction."""
    m, sids = _chlorobut_2_ene()
    units = m.stereo_units()
    assert len(units) == 1
    u = units[0]
    assert u['kind'] == 1 and u['anchor'] == sids[0]
    assert u['refs'] == (sids[1], sids[2], sids[4], None)
    assert u['unnamed_mask'] == 0b1000, 'slot 3 is the other terminal implicit H, a real direction'
    assert sum(1 for x in u['refs'] if x is None) == 1, \
        'exactly one unnamed slot: the shape no other bond fixture in this file has'


def test_chlorobut_2_ene_reversed_unit_structure():
    """The CH-first spelling: the same one unnamed slot, now in pair 0."""
    m, sids = _chlorobut_2_ene_reversed()
    units = m.stereo_units()
    assert len(units) == 1
    u = units[0]
    assert u['kind'] == 1 and u['anchor'] == sids[0]
    assert u['refs'] == (sids[1], None, sids[3], sids[4])
    assert u['unnamed_mask'] == 0b0010, 'slot 1 is the anchor terminal implicit H'
    assert sum(1 for x in u['refs'] if x is None) == 1


def test_chlorobut_2_ene_identity_returns_stored_parity():
    m, sids = _chlorobut_2_ene()
    assert m.translate_stereo(sids[0], (sids[1], sids[2], sids[4], None)) == 1
    m2, sids2 = _chlorobut_2_ene_reversed()
    assert m2.translate_stereo(sids2[0], (sids2[1], None, sids2[3], sids2[4])) == 1


def test_chlorobut_2_ene_within_pair_swap_flips():
    """One within-pair transposition is odd on either spelling; both slots are free to move."""
    m, sids = _chlorobut_2_ene()
    # swap pair 0 (Cl1 / C2), both named
    assert m.translate_stereo(sids[0], (sids[2], sids[1], sids[4], None)) == 2
    # swap pair 1 (C4 / the implicit H) -- legal, the slot is a real direction, not pinned
    assert m.translate_stereo(sids[0], (sids[1], sids[2], None, sids[4])) == 2
    m2, sids2 = _chlorobut_2_ene_reversed()
    # swap pair 0 (C1 / the implicit H)
    assert m2.translate_stereo(sids2[0], (None, sids2[1], sids2[3], sids2[4])) == 2


def test_chlorobut_2_ene_pair_exchange_does_not_flip():
    """The wholesale pair exchange is (0 2)(1 3), even, so it contributes nothing (F56)."""
    m, sids = _chlorobut_2_ene()
    assert m.translate_stereo(sids[0], (sids[4], None, sids[1], sids[2])) == 1
    m2, sids2 = _chlorobut_2_ene_reversed()
    assert m2.translate_stereo(sids2[0], (sids2[3], sids2[4], sids2[1], None)) == 1


def test_chlorobut_2_ene_none_does_not_correspond_raises():
    """A None placed where its pair holds a NAMED slot: the None-correspondence branch (F59).

    This is the test the None-correspondence branch never had, and this fixture is why: the
    branch is reachable exactly when the two pairs hold DIFFERENT numbers of SU_NO_REF slots.

    Cl-first spelling, refs = (Cl1, C2, C4, None), order (Cl1, None, C4, C2):
      phase 1 passes -- Cl1, C4 and C2 each appear once, and the None count is 1 = 1;
      src_pair resolves to pair 0 on want[0] = Cl1;
      the k=0 cross-pair check passes -- Cl1 is in pair 0 and want[1] is None, so it is skipped;
      then swap[0] = 0 (want[0] IS the named atom), so want[1]'s None must land on refs[1] --
      which is C2, a named atom.  That is the violation, and it fires at base 0.
    The reversed spelling reaches the same branch at base 2.

    `match=` is load-bearing: the branch is reachable but never the SOLE violation.  Here the
    misplaced None also strands C2 outside its pair, which the k=1 cross-pair check rejects one
    iteration later -- so with None correspondence disabled the call still raises, just with a
    different message.  Matching the message is what makes disabling the branch fail this test.
    """
    m, sids = _chlorobut_2_ene()
    with pytest.raises(ValueError, match='does not correspond to a'):
        m.translate_stereo(sids[0], (sids[1], None, sids[4], sids[2]))

    m2, sids2 = _chlorobut_2_ene_reversed()
    with pytest.raises(ValueError, match='does not correspond to a'):
        m2.translate_stereo(sids2[0], (sids2[3], None, sids2[1], sids2[4]))


# ---------------------------------------------------------------------------
# Task 8: parity re-basing in the journal apply.
#
# The stored bit means "the arrangement of THESE directions, in THIS order, has this
# handedness", and both halves of that sentence are properties of the arena the bit was
# written against.  The apply is the only writer, so the apply is where the bit is either
# re-expressed in the new order or dropped -- see the fragment comment on `rebase_parity`
# in `_stereo.pxi` for which of the two, and why.
#
# THE TWO SIDES OF THE DROP LINE, which these tests are split along:
#   * a frame that EXISTED AND WAS DESTROYED -- the atom anchored a unit, the edit left it
#     anchoring none or anchoring one whose directions are not the old ones.  The bit dies:
#     keeping it would silently re-read a sign against a frame nobody wrote it in.
#   * a frame that HAS NOT YET EXISTED -- the atom never anchored a unit at all.  The bit
#     lives: it was never interpreted against anything, and the molecule may still be
#     completed into something it means.  `test_a_parity_whose_frame_never_existed_survives`
#     is that side, and it is why the harvest snapshots configured UNITS and not configured
#     ATOMS.
# ---------------------------------------------------------------------------

def _named_refs(m, sid):
    """The unit's refs as element numbers, so they are comparable across relabelings."""
    unit = m.unit_of(sid)
    return tuple(None if r is None else m.element_of(r) for r in unit['refs'])


def _chiral_methane_with_a_spare_iodine():
    """`_chiral_methane` preceded by an unbonded iodine, so the anchor is NOT slot 0.

    The reachability fixture for tetrahedral re-basing.  `add_atom` appends, so a freshly
    created atom can never sort below an existing neighbour and can never permute the
    anchor's row -- the only way in is `add_bond` to an atom ALREADY in the molecule at a
    lower slot, which is what the spare iodine is here to be.
    """
    m = MoleculeContainer()
    with m.edit():
        sids = [m.add_atom(e, implicit_h=(1 if e == 'C' else None))
                for e in ('I', 'C', 'F', 'Cl', 'Br')]
        for s in sids[2:]:
            m.add_bond(sids[1], s, 1)
        m.set_parity(sids[1], 1)
    return m, sids


def _but_2_ene_with_a_spare_chlorine():
    """`_but_2_ene` preceded by an unbonded chlorine, for the same reason as the iodine above.

    Atom order: Cl0 (unbonded), C1(=CH-, anchor), C2(=CH-), C3(methyl on C1), C4(methyl on C2).
    The anchor terminal is slot 1, so Cl0 is the one atom that can enter its row in FRONT of
    the methyl.
    """
    m = MoleculeContainer()
    with m.edit():
        sids = [m.add_atom(e, implicit_h=h)
                for e, h in (('Cl', None), ('C', 1), ('C', 1), ('C', None), ('C', None))]
        m.add_bond(sids[1], sids[2], 2)
        m.add_bond(sids[1], sids[3], 1)
        m.add_bond(sids[2], sids[4], 1)
        m.set_parity(sids[1], 1)
    return m, sids


def _penta_2_3_diene_with_a_spare_chlorine():
    """CH3-CH=C=CH-CH3 preceded by an unbonded chlorine; the allene reachability fixture.

    Atom order: Cl0 (unbonded), C1(methyl), C2(=CH-, terminal), C3(centre, anchor),
    C4(=CH-, terminal), C5(methyl).  Chain length 3 is odd, so the unit is axial and anchored
    on C3; refs are the two terminals' pairs, C2's leading because C2 is the lower terminal.
    """
    m = MoleculeContainer()
    with m.edit():
        sids = [m.add_atom(e, implicit_h=h)
                for e, h in (('Cl', None), ('C', None), ('C', 1),
                             ('C', None), ('C', 1), ('C', None))]
        m.add_bond(sids[1], sids[2], 1)
        m.add_bond(sids[2], sids[3], 2)
        m.add_bond(sids[3], sids[4], 2)
        m.add_bond(sids[4], sids[5], 1)
        m.set_parity(sids[3], 1)
    return m, sids


# --- remap: ruling F64, the invariant that replaces the brief's `rebase_parities` ----------

def test_remap_leaves_every_stored_parity_bit_untouched():
    """Ruling F64: `remap` moves LABELS, not slots, so there is nothing to re-base.

    Parity is stored against the anchor's CSR row and a CSR row is keyed by SLOT.  `remap`
    clones the arena and writes new stable ids into the same slots, so every row keeps its
    order however the labels are permuted -- here they are reversed, F/Cl/Br taking
    descending ids, which is the mapping that WOULD flip the bit if the bit were stored
    against label order.  The brief asked for a `rebase_parities` for this path; it has no
    caller because this test's answer is "unchanged".
    """
    m, sids = _chiral_methane()
    stored_before = m.unit_of(sids[0])['parity']
    assert stored_before == 1, 'the fixture must arrive configured'
    m2 = m.copy()
    m2.remap({sids[0]: 100, sids[1]: 103, sids[2]: 102, sids[3]: 101})
    assert m2.unit_of(100)['parity'] == stored_before
    assert m2.parity_of(100) == stored_before


def test_remap_preserves_the_configuration():
    m, sids = _chiral_methane()
    before = m.translate_stereo(sids[0], (sids[1], sids[2], sids[3], None))
    m2 = m.copy()
    m2.remap({sids[0]: 100, sids[1]: 103, sids[2]: 102, sids[3]: 101})
    # the SAME geometric arrangement: F, Cl, Br under their new labels, in the same order
    after = m2.translate_stereo(100, (103, 102, 101, None))
    assert after == before, 'the same geometric arrangement must keep the same parity'
    # and the row itself did not move -- the refs are the same elements in the same places
    assert _named_refs(m2, 100) == _named_refs(m, sids[0])


# --- the drop rule: a frame that existed and was destroyed ---------------------------------

def test_a_fourth_heavy_neighbour_destroys_the_unit_and_clears_the_bit():
    """Five directions is not a unit, so the sign has nothing left to be a sign of.

    The core never derives implicit hydrogens, so a new bond does not consume the stated one:
    three heavy plus one implicit H is four directions, and the fourth heavy makes five.  The
    brief called this a re-basing case; it is a drop case, and `add_atom` appending is why it
    could not have been the other thing -- a fresh iodine lands last and permutes nothing.
    """
    m, sids = _chiral_methane()
    assert m.unit_of(sids[0]) is not None, 'the fixture must arrive as a unit'
    assert m.parity_of(sids[0]) == 1
    with m.edit():
        i = m.add_atom('I')
        m.add_bond(sids[0], i, 1)
    assert m.unit_of(sids[0]) is None
    assert m.parity_of(sids[0]) == 0, 'the bit must not survive its frame'
    assert m.stereo_of(sids[0]) == 0
    # and it is gone from the ATOM RECORD, not merely from this view: the flags are what
    # `to_bytes` serialises, so a bit still set here would travel into every buffer written from
    # this molecule and come back as a configuration nothing in the graph supports
    assert MoleculeContainer.from_bytes(m.to_bytes()).parity_of(sids[0]) == 0


def test_deleting_a_neighbour_drops_the_unit_and_clears_the_bit():
    m, sids = _chiral_methane()
    assert m.parity_of(sids[0]) == 1
    with m.edit():
        m.delete_atom(sids[3])
    assert m.unit_of(sids[0]) is None
    assert m.parity_of(sids[0]) == 0
    assert m.stereo_of(sids[0]) == 0


def test_deleting_a_bond_drops_the_unit_and_clears_the_bit():
    """The same frame destruction with the atom left in place: three directions, no unit."""
    m, sids = _chiral_methane()
    with m.edit():
        m.delete_bond(sids[0], sids[3])
    assert m.unit_of(sids[0]) is None
    assert m.parity_of(sids[0]) == 0


def test_growing_the_cumulene_chain_drops_the_centre_parity():
    """An allene's centre stops anchoring anything when the chain turns even.

    Chain length 3 is axial and anchored on the centre; extending it to 4 makes the unit
    cis/trans-like and anchors it on the lower TERMINAL, so the centre anchors nothing and
    its bit dies with the frame.
    """
    m, sids = _penta_2_3_diene_with_a_spare_chlorine()
    assert m.unit_of(sids[3])['kind'] == 2, 'the fixture must arrive axial'
    assert m.parity_of(sids[3]) == 1
    with m.edit():
        m.set_order(sids[4], sids[5], 2)      # C4=C5: the chain is now four atoms
        m.set_hydrogens(sids[4], 0)           # ...and C4 is an interior atom, not a terminal
    assert m.unit_of(sids[3]) is None
    assert m.parity_of(sids[3]) == 0


# --- the re-basing rule: the frame survived with its directions permuted -------------------

def test_substituting_the_implicit_hydrogen_rebases_a_tetrahedral_sign():
    """The one tetrahedral shape that survives an edit with its row permuted.

    Four directions before (F, Cl, Br, implicit H) and four after (I, F, Cl, Br), because the
    same edit that adds the bond states the hydrogen away.  The iodine was already in the
    molecule at slot 0, so it enters the anchor's row in FRONT of the fluorine and every named
    direction shifts one place: the cycle (3 0 1 2), three transpositions, odd, so the stored
    bit must flip.

    The second assertion is the one that says WHY: read in the OLD positional order -- the
    three heavy directions where they were, and the iodine standing in the place the implicit
    hydrogen occupied -- the answer is the parity that was stored before the edit.  The
    arrangement in space did not move; only its spelling did.
    """
    m, sids = _chiral_methane_with_a_spare_iodine()
    i, c, f, cl, br = sids
    assert m.unit_of(c)['refs'] == (f, cl, br, None)
    assert m.parity_of(c) == 1
    with m.edit():
        m.add_bond(c, i, 1)
        m.set_hydrogens(c, 0)
    assert m.unit_of(c)['refs'] == (i, f, cl, br), 'the row gained a front element'
    assert m.parity_of(c) == 2, 'an odd permutation of the directions flips the stored bit'
    assert m.translate_stereo(c, (f, cl, br, i)) == 1, 'the old frame still reads the old sign'


def test_substituting_the_implicit_hydrogen_rebases_a_cis_trans_sign():
    """The bond-kind re-basing case: one within-pair transposition, so the bit flips.

    The anchor terminal's pair is (methyl, implicit H); bonding the already-present chlorine
    -- slot 0, below the methyl -- and stating the hydrogen away makes it (Cl, methyl).  The
    methyl moves from offset 0 to offset 1 of its own pair, which is one transposition and odd.
    """
    m, sids = _but_2_ene_with_a_spare_chlorine()
    chlorine, anchor, far, near_me, far_me = sids
    assert m.unit_of(anchor)['refs'] == (near_me, None, far_me, None)
    assert m.unit_of(anchor)['unnamed_mask'] == 0b1010
    assert m.parity_of(anchor) == 1
    with m.edit():
        m.add_bond(anchor, chlorine, 1)
        m.set_hydrogens(anchor, 0)
    unit = m.unit_of(anchor)
    assert unit is not None, 'a disubstituted terminal is still a cis/trans terminal'
    assert unit['refs'] == (chlorine, near_me, far_me, None)
    assert m.parity_of(anchor) == 2, 'the within-pair transposition flips the stored bit'
    assert m.translate_stereo(anchor, (near_me, chlorine, far_me, None)) == 1


def test_substituting_the_implicit_hydrogen_rebases_an_allene_sign():
    """The same within-pair transposition on an axial unit, whose anchor is the chain centre."""
    m, sids = _penta_2_3_diene_with_a_spare_chlorine()
    chlorine, near_me, near, centre, far, far_me = sids
    assert m.unit_of(centre)['refs'] == (near_me, None, far_me, None)
    assert m.parity_of(centre) == 1
    with m.edit():
        m.add_bond(near, chlorine, 1)
        m.set_hydrogens(near, 0)
    unit = m.unit_of(centre)
    assert unit is not None, 'the axis survives a disubstituted terminal'
    assert unit['kind'] == 2 and unit['refs'] == (chlorine, near_me, far_me, None)
    assert m.parity_of(centre) == 2
    assert m.translate_stereo(centre, (near_me, chlorine, far_me, None)) == 1


def test_explicitating_the_hydrogen_leaves_the_sign_alone():
    """Ruling F26's promise, measured: drawing the hydrogen re-bases nothing.

    A hydrogen direction sorts after every heavy one whether it is named or not, so the
    implicit H at slot 3 becomes the drawn H at slot 3 -- the identity permutation -- and the
    stored bit is neither flipped nor dropped.  This is the case a strict "the direction set
    must be identical" rule would silently lose, and losing it would cost the parity of every
    molecule that gets its hydrogens drawn.
    """
    m, sids = _chiral_methane()
    with m.edit():
        h = m.add_atom('H')
        m.add_bond(sids[0], h, 1)
        m.set_hydrogens(sids[0], 0)
    unit = m.unit_of(sids[0])
    assert unit is not None
    assert unit['refs'] == (sids[1], sids[2], sids[3], h)
    assert unit['unnamed_mask'] == 0, 'the direction is named now'
    assert m.parity_of(sids[0]) == 1, 'the same order, so the same bit'
    assert m.translate_stereo(sids[0], (sids[1], sids[2], sids[3], h)) == 1


def test_explicitating_both_cis_trans_hydrogens_leaves_the_sign_alone():
    """Ruling F69's other half: ONE substitution IN EACH PAIR is two in the unit and must be kept.

    Both terminals of but-2-ene carry an implicit hydrogen, and drawing both of them in one edit
    leaves two old positions with nothing of their own to correspond to -- one in each pair.  A
    leftover budget counted over the WHOLE UNIT would see two and drop the sign, which ruling F26
    forbids; counted PER DIRECTION LIST it is one and one, each the identity within its own pair, so
    the bit stands.  This is the case that makes the per-list granularity load-bearing, and it is the
    exact converse of `test_rebase_refuses_two_arrivals_into_one_cis_trans_terminal`, where the same
    two leftovers fall in ONE pair and the sign dies.
    """
    m, sids = _but_2_ene()
    anchor = sids[1]
    assert m.unit_of(anchor)['refs'] == (sids[0], None, sids[3], None)
    assert m.unit_of(anchor)['unnamed_mask'] == 0b1010, 'one implicit H on each terminal'
    assert m.parity_of(anchor) == 1
    with m.edit():
        near_h = m.add_atom('H')
        far_h = m.add_atom('H')
        m.add_bond(anchor, near_h, 1)
        m.add_bond(sids[2], far_h, 1)
        m.set_hydrogens(anchor, 0)
        m.set_hydrogens(sids[2], 0)
    unit = m.unit_of(anchor)
    assert unit is not None, 'both terminals still carry two directions'
    assert unit['refs'] == (sids[0], near_h, sids[3], far_h)
    assert unit['unnamed_mask'] == 0, 'every direction is named now'
    assert m.parity_of(anchor) == 1, 'the identity within each pair, so the same bit'
    assert m.translate_stereo(anchor, (sids[0], near_h, sids[3], far_h)) == 1


# --- re-basing across the PINNED / UNNAMED flavours ----------------------------------------
#
# This file's STANDING REQUIREMENT above applies to the re-basing fixtures too, and it was breached
# by the first round of them: every one carried its unnamed directions in the same flavour, so the
# EMPTY-versus-UNNAMED law `rebase_parity` now enforces (ruling F69) was untested.  The three tests
# below cross the flavours on real molecules -- an oxime nitrogen's lone pair, which is a slot with NO
# direction in it, an oxime carbon's implicit hydrogen, which is a real unnamed direction in the very
# same record, and a sulfoxide sulfur's lone pair, which is an unnamed direction on an ATOM kind.
#
# AND IT BOUNDS WHAT THESE CAN REACH, structurally: NO PERCEIVED DIRECTION LIST HOLDS BOTH FLAVOURS.  A
# SU_TETRA record's four slots are all real directions -- a neighbour, an implicit hydrogen, or the one
# sulfur lone pair -- so an EMPTY slot cannot occur there at all; EMPTY occurs only as the second slot of
# a bond-kind pair whose terminal has one direction, and `_terminal_pair` refuses a pair whose BOTH slots
# are nameless.  At most one nameless slot per list means the correspondence is forced whichever way the
# mask reads, so on reachable frames the mask confirms rather than decides -- which is exactly why the law
# needs a FORGED frame to be tested at all (`test_rebase_obeys_the_empty_versus_unnamed_law`), and why it
# is kept: it becomes decisive the day a kind with two nameless slots in one list is perceived.


def _acetaldoxime_with_a_spare_methyl():
    """`_acetaldoxime` preceded by an unbonded methyl, so a new N substituent sorts BEFORE the O.

    Atom order: C0 (spare methyl, unbonded), C1 (methyl), C2 (=CH-, anchor), N3, O4 (OH).
    refs = (C1, None, O4, None), unnamed_mask = 0b0010: slot 1 is C2's implicit hydrogen, a real
    unnamed direction, and slot 3 is the nitrogen's lone pair, which `_terminal_pair` does not count
    as a direction at all.  The spare methyl is at slot 0 because CSR order is what decides where a
    new N substituent lands in the nitrogen's own pair.
    """
    m = MoleculeContainer()
    with m.edit():
        sids = [m.add_atom(e, implicit_h=h)
                for e, h in (('C', 3), ('C', 3), ('C', 1), ('N', 0), ('O', 1))]
        m.add_bond(sids[1], sids[2], 1)
        m.add_bond(sids[2], sids[3], 2)
        m.add_bond(sids[3], sids[4], 1)
        m.set_parity(sids[2], 1)
    return m, sids


def _methanesulfinyl_fixture():
    """Dimethyl sulfoxide preceded by an unbonded chlorine; the sulfur-lone-pair fixture.

    Atom order: Cl0 (unbonded), S1, O2, C3 (methyl), C4 (methyl).  The sulfur reaches four directions
    as two sigma, one pi and ONE LONE PAIR, so refs = (O2, C3, C4, None) with unnamed_mask = 0b1000 --
    an SU_TETRA record whose unnamed direction is a lone pair rather than an implicit hydrogen, which
    is the flavour no re-basing fixture carried before.
    """
    m = MoleculeContainer()
    with m.edit():
        sids = [m.add_atom(e, implicit_h=h)
                for e, h in (('Cl', None), ('S', 0), ('O', None), ('C', 3), ('C', 3))]
        m.add_bond(sids[1], sids[2], 2)
        m.add_bond(sids[1], sids[3], 1)
        m.add_bond(sids[1], sids[4], 1)
        m.set_parity(sids[1], 1)
    return m, sids


def test_displacing_an_sp2_nitrogen_lone_pair_rebases_the_sign():
    """An empty slot filled by a real substituent: the pair's own transposition, so the bit flips.

    Methylating acetaldoxime's nitrogen gives the nitrone CH3-CH=N(+)(CH3)-O(-).  The nitrogen's pair
    was `(O, <no direction>)`; it becomes `(CH3, O)`, because the spare methyl sits at a lower slot
    than the oxygen, so the oxygen moves from offset 0 to offset 1 of its own pair.  One within-pair
    transposition, odd, so the stored bit flips.

    The anchor's OWN pair is untouched and its implicit hydrogen has to stay where it is: it is an
    UNNAMED direction and the only unnamed slot of pair 0, and a correspondence that let it drift into
    the nitrogen's freed position would compute a different answer.
    """
    m, sids = _acetaldoxime_with_a_spare_methyl()
    spare, methyl, anchor, n, o = sids
    unit = m.unit_of(anchor)
    assert unit is not None and unit['kind'] == 1
    assert unit['refs'] == (methyl, None, o, None)
    assert unit['unnamed_mask'] == 0b0010, 'slot 1 is a direction, slot 3 is not'
    assert m.parity_of(anchor) == 1
    with m.edit():
        m.add_bond(n, spare, 1)
        m.set_charge(n, 1)
        m.set_charge(o, -1)
        m.set_hydrogens(o, 0)
    after = m.unit_of(anchor)
    assert after is not None, 'a disubstituted nitrogen is still a cis/trans terminal'
    assert after['refs'] == (methyl, None, spare, o), 'the lone pair slot is a substituent now'
    assert after['unnamed_mask'] == 0b0010, "and the carbon's hydrogen is still unnamed"
    assert m.parity_of(anchor) == 2, 'one within-pair transposition flips the stored bit'
    assert m.translate_stereo(anchor, (methyl, None, o, spare)) == 1, \
        'read in the old positional order the arrangement is the sign that was stored'


def test_an_oxime_empty_slot_corresponds_to_the_empty_slot_not_to_the_drawn_hydrogen():
    """Both flavours in one record, and each keeps to its own kind: the sign is untouched.

    Drawing acetaldoxime's CH hydrogen names slot 1 and leaves slot 3 -- the nitrogen's lone pair --
    still no direction at all.  The empty slot must correspond to the empty slot; a correspondence
    that handed it the newly drawn hydrogen instead would produce a permutation across the pair
    boundary and drop a sign ruling F26 promises to keep.
    """
    m, sids = _acetaldoxime()
    methyl, anchor, n, o = sids
    assert m.unit_of(anchor)['refs'] == (methyl, None, o, None)
    assert m.unit_of(anchor)['unnamed_mask'] == 0b0010
    assert m.parity_of(anchor) == 1
    with m.edit():
        h = m.add_atom('H')
        m.add_bond(anchor, h, 1)
        m.set_hydrogens(anchor, 0)
    after = m.unit_of(anchor)
    assert after is not None
    assert after['refs'] == (methyl, h, o, None), 'slot 1 named, slot 3 still not a direction'
    assert after['unnamed_mask'] == 0, 'nothing unnamed is left; slot 3 was never unnamed'
    assert m.parity_of(anchor) == 1, 'the identity permutation, so the same bit'
    assert m.translate_stereo(anchor, (methyl, h, o, None)) == 1


def test_substituting_a_neighbour_across_a_sulfur_lone_pair_rebases_the_sign():
    """An SU_TETRA re-base whose fourth direction is a LONE PAIR, not an implicit hydrogen.

    Dimethyl sulfoxide becomes methanesulfinyl chloride: one methyl leaves and the already-present
    chlorine arrives at slot 0, in FRONT of the oxygen, so every named direction shifts.  The lone
    pair is an unnamed direction and stays in slot 3 -- it corresponds to itself, by flavour, and the
    substitution is confined to the position the methyl vacated.  The row `(O, CH3, CH3, pair)`
    becomes `(Cl, O, CH3, pair)`, one transposition, odd, so the bit flips.
    """
    m, sids = _methanesulfinyl_fixture()
    chlorine, s, o, leaving, staying = sids
    unit = m.unit_of(s)
    assert unit is not None and unit['kind'] == 0
    assert unit['refs'] == (o, leaving, staying, None)
    assert unit['unnamed_mask'] == 0b1000, 'the fourth direction is the sulfur lone pair'
    assert m.parity_of(s) == 1
    with m.edit():
        m.delete_bond(s, leaving)
        m.add_bond(s, chlorine, 1)
    after = m.unit_of(s)
    assert after is not None, 'two sigma, one pi and the pair is still four directions'
    assert after['refs'] == (chlorine, o, staying, None)
    assert after['unnamed_mask'] == 0b1000, 'and the pair is still the unnamed one'
    assert m.parity_of(s) == 2, 'one transposition of the row flips the stored bit'
    assert m.translate_stereo(s, (o, chlorine, staying, None)) == 1


def test_rebase_obeys_the_empty_versus_unnamed_law():
    """The mask is CONSULTED: the same refs with two different flavour maps give two answers.

    Toluene's methyl carbon is a real record with three unnamed directions, `(ring, -, -, -)` and mask
    `0b1110`.  Read the same refs back as a frame whose slots 1 and 3 were EMPTY and only slot 2 a
    real direction, and the two empty slots have nothing of their own to correspond to -- an empty
    slot may only answer to an empty slot, and this record has none.  That is two leftovers in one
    direction list and the sign dies, where the honest mask re-bases it through the identity.

    A FORGED frame, and it has to be: a perceived SU_TETRA record has no empty slot at all and
    `_terminal_pair` refuses a bond-kind pair whose both slots are nameless, so no perceived direction
    list holds both flavours (the note above this block has the argument).  It is here because the law is
    what the code implements, and a law that is only enforced where no caller can see it is a law that
    will be deleted by the next reader -- and because the day a kind with two nameless slots in one list
    is perceived, this stops being defence in depth and starts deciding signs.
    """
    m, sids = _build(['C'] * 7, [3, 0, 1, 1, 1, 1, 1],
                     [(0, 1, 1), (1, 2, 2), (2, 3, 1), (3, 4, 2), (4, 5, 1), (5, 6, 2), (6, 1, 1)])
    unit = m.unit_of(sids[0])
    assert unit is not None, 'the fixture must arrive as a unit'
    assert unit['refs'] == (sids[1], None, None, None)
    assert unit['unnamed_mask'] == 0b1110, 'three implicit hydrogens, all real directions'
    frame = (sids[1], None, None, None)
    assert _core._rebase_parity_probe(m, sids[0], 0, 1, frame, 0b1110) == 1
    assert _core._rebase_parity_probe(m, sids[0], 0, 1, frame, 0b0100) == -1


# --- what the apply must NOT touch ---------------------------------------------------------

def test_an_edit_elsewhere_leaves_the_configuration_alone():
    """An edit that does not reach the anchor's frame must change neither bit nor reading."""
    m, sids = _chiral_methane()
    before = m.translate_stereo(sids[0], (sids[1], sids[2], sids[3], None))
    with m.edit():
        a = m.add_atom('C', implicit_h=3)
        b = m.add_atom('O', implicit_h=1)
        m.add_bond(a, b, 1)
    assert m.parity_of(sids[0]) == before
    assert m.translate_stereo(sids[0], (sids[1], sids[2], sids[3], None)) == before


def test_a_parity_whose_frame_never_existed_survives_the_apply():
    """A bit on an atom that never anchored a unit is data, not a stale sign.

    The drop rule is scoped to the units the apply HARVESTED, deliberately: a global "clear
    every parity whose atom anchors no unit" sweep would delete this datum at apply time, and
    a container mid-edit legitimately carries a parity nothing justifies yet.  Reporting and
    clearing it is `validate_stereo`'s job, on the consumer's demand -- not the apply's.
    """
    m = MoleculeContainer()
    with m.edit():
        c = m.add_atom('C')
        m.set_parity(c, 1)
    assert m.unit_of(c) is None, 'one direction is not a frame'
    assert m.parity_of(c) == 1, 'the bit was never interpreted, so nothing invalidated it'
    # and a second, unrelated edit does not invent a reason to clear it either
    with m.edit():
        m.add_atom('O')
    assert m.parity_of(c) == 1


def test_a_parity_stated_in_the_edit_that_destroys_the_frame_is_the_callers_own():
    """Ruling F73: the journal's own write outranks the DROP as well as any re-basing.

    The harvest skips an anchor the journal states a parity for, and that skip is not narrower than it
    looks -- it removes the unit from the snapshot entirely, so the apply neither re-bases the old
    value nor clears it.  THAT IS THE DECISION, not a leak.  A caller who states a parity inside the
    very edit that destroys the frame is stating a sign against the molecule the edit produces, not
    asking for the old one to be carried over; the apply has no standing to overrule it and no way to
    tell "meant it" from "forgot", so it leaves the datum alone.

    What is left behind is ruling F66's second case -- a configured bit whose frame has not yet
    existed -- and it is READABLE and it SURVIVES SERIALISATION, asserted below so that the next
    reader meets it as a decision rather than discovering it as a surprise.  Deciding whether such a
    bit is justified, reporting it, and clearing it if the consumer asks, is `validate_stereo`'s remit.
    The apply must never clear it eagerly: the same sweep that would tidy this case away
    also deletes the legitimate mid-edit parity in
    `test_a_parity_whose_frame_never_existed_survives_the_apply`.
    """
    m, sids = _chiral_methane()
    c = sids[0]
    assert m.unit_of(c) is not None and m.parity_of(c) == 1
    with m.edit():
        i = m.add_atom('I')
        m.add_bond(c, i, 1)     # five directions: the frame is destroyed
        m.set_parity(c, 2)      # ...and the caller states the OTHER value anyway
    assert m.unit_of(c) is None, 'five directions is not a unit'
    # 2 is neither the harvested 1 nor the dropped 0, so this distinguishes all three outcomes
    assert m.parity_of(c) == 2, "the caller's own write is not the apply's to drop"
    assert m.stereo_of(c) == 1, 'and it reached the atom flags, not just this view'
    # ...so it travels: this is exactly the bit `validate_stereo` is for
    assert MoleculeContainer.from_bytes(m.to_bytes()).parity_of(c) == 2


def test_an_explicit_parity_in_the_same_edit_wins_over_the_rebase():
    """The journal's own write is the caller speaking; a re-base must not overwrite it."""
    m, sids = _chiral_methane_with_a_spare_iodine()
    i, c = sids[0], sids[1]
    with m.edit():
        m.add_bond(c, i, 1)
        m.set_hydrogens(c, 0)
        m.set_parity(c, 1)          # the re-base alone would have made this 2
    assert m.parity_of(c) == 1


def test_the_rebase_needs_one_batched_edit():
    """Split the same two operations across two applies and the sign is gone, correctly.

    After the first apply the centre has five directions and anchors nothing, so its frame is
    destroyed and the bit dies there -- the second apply has nothing left to re-base.  This is
    a property of the edits, not a defect: an edit that means to preserve a configuration has
    to leave the molecule a frame to preserve it against at every apply boundary.
    """
    m, sids = _chiral_methane_with_a_spare_iodine()
    i, c = sids[0], sids[1]
    m.add_bond(c, i, 1)
    assert m.parity_of(c) == 0
    m.set_hydrogens(c, 0)
    assert m.unit_of(c) is not None, 'the frame is back'
    assert m.parity_of(c) == 0, 'but the sign it held is not, and is not invented'


# --- the table is derived, so a round trip must rebuild it from the CSR and SEG_PARITY -----

def test_configuration_survives_a_serialisation_round_trip():
    """to_bytes/from_bytes, not pack/unpack -- `MoleculeContainer.pack` and
    `MoleculeContainer.unpack` raise NotImplementedError on this branch (they are reserved for
    the chython 2 pach format).

    SEG_STEREO_UNIT is a DERIVED segment and is not serialised, so `from_bytes` has to rebuild
    the whole unit table from the persistent CSR plus the anchor's parity byte in SEG_PARITY --
    which is the same derivation the apply's replay depends on.
    """
    m, sids = _chiral_methane()
    before = m.translate_stereo(sids[0], (sids[1], sids[2], sids[3], None))
    m2 = MoleculeContainer.from_bytes(m.to_bytes())
    assert m2.unit_of(sids[0]) is not None, 'the table must be rebuilt, not carried'
    assert m2.translate_stereo(sids[0], (sids[1], sids[2], sids[3], None)) == before


# --- the fourth kind: an atropisomer's directions cannot be permuted -----------------------
# An atropisomer pivot's two directions are both NAMED ring atoms, so there is no unnamed slot for
# a new neighbour to displace -- and the pivot's direction set IS its two non-axis bonds, so
# changing that set changes its degree away from 3 or takes a bond out of a ring, either of which
# ends the axis.  What is left is the PAIR EXCHANGE that ruling F45's anchor relocation performs,
# and the two tests below measure both halves: no permutation across an edit that keeps the anchor,
# and a drop when the anchor itself moves.

_BIPHENYL_BONDS = [(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 0, 1),
                   (6, 7, 2), (7, 8, 1), (8, 9, 2), (9, 10, 1), (10, 11, 2), (11, 6, 1),
                   (0, 6, 1)]
_BIPHENYL_H = [0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1]
# Kekule cyclooctatetraene twice: ring 1 is atoms 0-7 and ring 2 atoms 8-15, each alternating from
# its own lowest atom.  Eight is the smallest ring the small-ring cut does not reach, so each pivot
# carries a cis/trans unit of its own -- which is what makes the axis relocate.
_COT_ONE = [(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 6, 1), (6, 7, 2), (7, 0, 1)]
_COT_TWO = [(8, 9, 2), (9, 10, 1), (10, 11, 2), (11, 12, 1), (12, 13, 2), (13, 14, 1),
            (14, 15, 2), (15, 8, 1)]


def _build(elements, hydrogens, bonds):
    m = MoleculeContainer()
    with m.edit():
        sids = [m.add_atom(e, implicit_h=h) for e, h in zip(elements, hydrogens)]
        for i, j, o in bonds:
            m.add_bond(sids[i], sids[j], o)
    return m, sids


def _chlorofluorobiphenyl():
    """2-chloro-2'-fluorobiphenyl: an atropisomer axis anchored on its lower pivot, atom 0."""
    m, sids = _build(['C'] * 12 + ['Cl', 'F'], _BIPHENYL_H + [0, 0],
                     _BIPHENYL_BONDS + [(1, 12, 1), (7, 13, 1)])
    m.set_parity(sids[0], 1)
    return m, sids


def test_an_atropisomer_keeps_its_ref_order_and_sign_across_an_edit():
    """Measured: nothing an edit can do permutes an atropisomer's directions.

    The axis' four directions are the two pivots' ring neighbours, all four named.  An edit that
    changes any of them changes a pivot's degree away from 3 or unrings a bond, and then there is no
    axis to re-base; an edit that does not touch them leaves CSR order alone, because slot
    compaction is monotone and no operation in the core permutes slots.  So the sign is carried
    through unchanged -- neither re-based nor dropped.
    """
    m, sids = _chlorofluorobiphenyl()
    unit = m.unit_of(sids[0])
    assert unit is not None and unit['kind'] == 3
    assert unit['refs'] == (sids[1], sids[5], sids[7], sids[11])
    assert unit['unnamed_mask'] == 0, 'all four directions are named ring atoms'
    with m.edit():                       # an unrelated fragment: methanol, touching nothing
        a = m.add_atom('C', implicit_h=3)
        b = m.add_atom('O', implicit_h=1)
        m.add_bond(a, b, 1)
    after = m.unit_of(sids[0])
    assert after is not None and after['refs'] == (sids[1], sids[5], sids[7], sids[11])
    assert m.parity_of(sids[0]) == 1


def test_relocating_an_atropisomer_anchor_drops_the_sign():
    """The one permutation an axis has -- the pair exchange -- comes with the anchor moving, and a
    harvest keyed on the ANCHOR drops it rather than following it.

    An ortho-dichloro bi(cyclooctatetraenyl) anchors its axis on the UPPER pivot, atom 15, because
    the lower pivot's own ring double bond claimed it first (ruling F45).  Moving that double bond
    one bond round ring 1 frees atom 0, the axis returns to the lower pivot, and the two ref pairs
    trade places -- the same molecule, the same axis, described against a different atom.

    THE SIGN IS DROPPED, and that is a known cost rather than a correctness bug: the snapshot names
    its anchor, the far pivot is not among the refs, and finding the relocated unit would mean
    re-walking the axis.  Dropping is the safe direction -- the alternative is a sign read against
    pairs that traded places.  A caller who relocates an anchor must re-state the configuration.
    """
    m, sids = _build(['C'] * 16 + ['Cl', 'Cl'],
                     [0, 0] + [1] * 6 + [1] * 6 + [0, 0] + [0, 0],
                     _COT_ONE + _COT_TWO + [(0, 15, 1), (1, 16, 1), (14, 17, 1)])
    axis = m.unit_of(sids[15])
    assert axis is not None and axis['kind'] == 3
    assert axis['refs'] == (sids[8], sids[14], sids[1], sids[7]), 'the anchor pair leads'
    m.set_parity(sids[15], 1)
    assert m.parity_of(sids[15]) == 1
    with m.edit():                        # 0=1 becomes 1=2: atom 0 is no longer a cis/trans terminal
        m.set_order(sids[0], sids[1], 1)
        m.set_order(sids[1], sids[2], 2)
    moved = m.unit_of(sids[0])
    assert moved is not None and moved['kind'] == 3, 'the axis is still there, on the lower pivot'
    assert moved['refs'] == (sids[1], sids[7], sids[8], sids[14]), 'the pairs traded places'
    assert m.unit_of(sids[15]) is None, 'and atom 15 anchors nothing now'
    assert m.parity_of(sids[15]) == 0, 'so its sign is dropped, not moved'
    assert m.parity_of(sids[0]) == 0, 'and never invented on the new anchor'


# --- the two refusals a molecule cannot stage, measured through the probe ------------------

def test_rebase_refuses_a_kind_change_under_an_unmoved_anchor():
    """A sign about an axis is not a sign about a centre, even on the same atom.

    No edit reaches this: a kind change also changes the anchor's directions, so the leftover
    budget refuses one step earlier.  The guard is still live -- the pair arithmetic
    below it is only meaningful for the kind that was measured -- so it is measured here with the
    frame supplied by hand.
    """
    m, sids = _chiral_methane()
    frame = (sids[1], sids[2], sids[3], None)
    assert m.unit_of(sids[0])['unnamed_mask'] == 0b1000, 'slot 3 is the implicit hydrogen'
    # the same frame read as the kind it actually is: carried through unchanged
    assert _core._rebase_parity_probe(m, sids[0], 0, 1, frame, 0b1000) == 1
    # ...and read as a cis/trans frame: dropped
    assert _core._rebase_parity_probe(m, sids[0], 1, 1, frame, 0b1000) == -1


def test_rebase_refuses_a_bond_frame_whose_atoms_changed_ends():
    """Ruling F75: an old direction found in the OTHER pair MIGRATED, and no arithmetic saves that.

    2,3-dichlorobut-2-ene has all four directions named, which is the only shape in which a frame can
    cross the pair boundary without also losing a direction.  Both crossing frames below are forged, for
    the same reason: `rebase_parity` looks the new unit up by the anchor slot the old frame came from,
    and by ruling F26 and `_stereo_emit` pair 0 is that anchor's own end in both records, so a named
    atom that changed pairs cannot have been merely re-listed -- it is bonded to the other end of the
    double bond now, which is a different molecule and not a permutation of the old frame.

    That covers the WHOLESALE EXCHANGE as well as the partial mix, and the exchange is the one the
    reachable test `test_a_substituent_that_changes_ends_drops_the_sign` measures through real edits: an
    earlier ruling accepted it as even (`(0 2)(1 3)` is two transpositions, and the arithmetic is not
    what was wrong) and that KEPT a sign across four broken bonds.

    `translate_stereo` accepts the very same exchange as even
    (`test_dichlorobut_2_ene_pair_exchange_does_not_flip`) and is right to, and the divergence is
    deliberate: its frame is a CALLER'S ordering, where listing the two pairs the other way round is a
    re-ordering of one geometry, while both frames here are perceived, canonical and keyed to one anchor,
    where the pair order is nobody's to choose.  Different questions, so different answers.
    """
    m, sids = _dichlorobut_2_ene()
    assert m.unit_of(sids[0])['refs'] == (sids[1], sids[4], sids[3], sids[5])
    assert m.unit_of(sids[0])['unnamed_mask'] == 0, 'all four directions are named'
    # the frame as it stands: no permutation, no flip
    assert _core._rebase_parity_probe(m, sids[0], 1, 1,
                                      (sids[1], sids[4], sids[3], sids[5]), 0) == 1
    # one within-pair transposition: the sign flips, and this is the arithmetic the reachable
    # re-basing tests exercise through a real edit
    assert _core._rebase_parity_probe(m, sids[0], 1, 1,
                                      (sids[4], sids[1], sids[3], sids[5]), 0) == 2
    # old pair 0 is (Cl1, Cl3), whose atoms are in DIFFERENT new pairs: dropped
    assert _core._rebase_parity_probe(m, sids[0], 1, 1,
                                      (sids[1], sids[3], sids[4], sids[5]), 0) == -1
    # ...and the wholesale exchange, where every atom of each pair is in the other one: also dropped,
    # for both stored values, because there is no answer rather than an answer that happens to be even
    assert _core._rebase_parity_probe(m, sids[0], 1, 1,
                                      (sids[3], sids[5], sids[1], sids[4]), 0) == -1
    assert _core._rebase_parity_probe(m, sids[0], 1, 2,
                                      (sids[3], sids[5], sids[1], sids[4]), 0) == -1
    # exchanged with a within-pair transposition on top of it: still a migration, still dropped
    assert _core._rebase_parity_probe(m, sids[0], 1, 1,
                                      (sids[5], sids[3], sids[1], sids[4]), 0) == -1
    assert _core._rebase_parity_probe(m, sids[0], 1, 1,
                                      (sids[5], sids[3], sids[4], sids[1]), 0) == -1


def _bromochlorofluoroiodoethene(spare=False):
    """1-bromo-1-chloro-2-fluoro-2-iodoethene, parity 1: the frame with FOUR named directions.

    Atom order: C0 (the Cl/Br end and the anchor, lower slot), C1, Cl2, Br3, F4, I5, and with
    `spare=True` an unbonded At6 to substitute with.  refs = (Cl2, Br3, F4, I5), mask 0 -- every slot
    named, which is what makes a migration across the double bond visible as a permutation instead of
    dying on the leftover budget first.  No `implicit_h` anywhere: both carbons are already at three
    sigma neighbours plus the pi bond.
    """
    m = MoleculeContainer()
    with m.edit():
        els = ['C', 'C', 'Cl', 'Br', 'F', 'I'] + (['At'] if spare else [])
        sids = [m.add_atom(e) for e in els]
        m.add_bond(sids[0], sids[1], 2)
        m.add_bond(sids[0], sids[2], 1)
        m.add_bond(sids[0], sids[3], 1)
        m.add_bond(sids[1], sids[4], 1)
        m.add_bond(sids[1], sids[5], 1)
        m.set_parity(sids[0], 1)
    return m, sids


def test_a_substituent_that_changes_ends_drops_the_sign():
    """Ruling F75, and it is REACHABLE -- four real batched edits, no probe.

    Pair 0 of the record is the anchor's own end and the anchor does not move here, so a named direction
    that turns up in the other pair is bonded to the other carbon now: the four bonds were broken and
    remade, and the old sign says nothing about the result.  All four spellings of that are dropped, and
    the first two are the ones an earlier ruling KEPT -- accepting the wholesale exchange as an even
    permutation, which the arithmetic supports and the chemistry does not.
    """
    for case, edit, want_refs in (
            # every substituent moves to the opposite end: reads as the wholesale pair exchange
            ('the wholesale exchange',
             lambda m, s: (m.delete_bond(s[0], s[2]), m.delete_bond(s[0], s[3]),
                           m.delete_bond(s[1], s[4]), m.delete_bond(s[1], s[5]),
                           m.add_bond(s[0], s[4], 1), m.add_bond(s[0], s[5], 1),
                           m.add_bond(s[1], s[2], 1), m.add_bond(s[1], s[3], 1)),
             lambda s: (s[4], s[5], s[2], s[3])),
            # ...and two atoms of one pair going to different new pairs: the partial mix
            ('one atom from each pair crossing',
             lambda m, s: (m.delete_bond(s[0], s[3]), m.delete_bond(s[1], s[5]),
                           m.add_bond(s[0], s[5], 1), m.add_bond(s[1], s[3], 1)),
             lambda s: (s[2], s[5], s[3], s[4])),
            # one atom crosses and an implicit hydrogen takes its place: one leftover per list, so the
            # budget would allow it and only the crossing refusal does not
            ('one crossing with a hydrogen behind it',
             lambda m, s: (m.delete_bond(s[0], s[3]), m.delete_bond(s[1], s[5]),
                           m.add_bond(s[1], s[3], 1), m.set_hydrogens(s[0], 1)),
             lambda s: (s[2], None, s[3], s[4])),
    ):
        m, sids = _bromochlorofluoroiodoethene()
        assert m.unit_of(sids[0])['refs'] == (sids[2], sids[3], sids[4], sids[5]), case
        assert m.parity_of(sids[0]) == 1, case
        with m.edit():
            edit(m, sids)
        unit = m.unit_of(sids[0])
        assert unit is not None, f'{case}: the double bond and both terminals survive'
        assert unit['refs'] == want_refs(sids), case
        assert m.parity_of(sids[0]) == 0, f'{case}: nothing left to re-base along'

    # the exchange with a substitution on top of it -- one leftover in each list, so this one passes
    # the per-list budget and is refused purely for having changed ends
    m, sids = _bromochlorofluoroiodoethene(spare=True)
    with m.edit():
        m.delete_bond(sids[0], sids[2])
        m.delete_bond(sids[0], sids[3])
        m.delete_bond(sids[1], sids[4])
        m.delete_bond(sids[1], sids[5])
        m.add_bond(sids[0], sids[4], 1)
        m.add_bond(sids[0], sids[5], 1)
        m.add_bond(sids[1], sids[2], 1)
        m.add_bond(sids[1], sids[6], 1)
        m.delete_atom(sids[3])
    unit = m.unit_of(sids[0])
    assert unit is not None, 'both terminals still carry two directions'
    assert unit['refs'] == (sids[4], sids[5], sids[2], sids[6])
    assert m.parity_of(sids[0]) == 0, 'exchanged AND substituted is still exchanged'


def test_rebase_probe_refuses_an_out_of_range_kind_or_parity():
    """The Python-visible door validates rather than passing rubbish into the arithmetic.

    Unguarded, `parity=7` comes back out as 7, which is not a parity at all.  The probe is test-only,
    but the convention in this file is to guard, and a probe that launders bad input is a probe that
    can certify a defect.
    """
    m, sids = _chiral_methane()
    frame = (sids[1], sids[2], sids[3], None)
    with pytest.raises(ValueError, match='parity must be'):
        _core._rebase_parity_probe(m, sids[0], 0, 7, frame, 0b1000)
    with pytest.raises(ValueError, match='kind must be'):
        _core._rebase_parity_probe(m, sids[0], 9, 1, frame, 0b1000)
    with pytest.raises(ValueError, match='old_unnamed_mask must be'):
        _core._rebase_parity_probe(m, sids[0], 0, 1, frame, 0b10000)
    # ...and 0 is a legal parity, meaning "nothing configured", which re-bases to itself
    assert _core._rebase_parity_probe(m, sids[0], 0, 0, frame, 0b1000) == 0


def test_rebase_takes_one_substituted_direction_and_refuses_two():
    """One direction replaced is the hydrogen case with a heavier atom in it; two is a guess.

    The vanished direction leaves its POSITION to whatever the new refs hold there, which for the
    first frame below is the bromine that used to sit behind it -- an identity permutation, so the
    sign stands.  The second frame has two positions with nothing to correspond to, and no
    arithmetic decides which of the two survivors took which; the sign dies instead.
    """
    m, sids = _chiral_methane()
    with m.edit():
        spare = m.add_atom('I')          # in the molecule, but not a neighbour of the centre
        other = m.add_atom('At')
    assert _core._rebase_parity_probe(m, sids[0], 0, 1,
                                      (sids[1], sids[2], spare, None), 0b1000) == 1
    assert _core._rebase_parity_probe(m, sids[0], 0, 1,
                                      (sids[1], spare, other, None), 0b1000) == -1


def test_rebase_refuses_two_arrivals_into_one_tetrahedral_row():
    """Ruling F69's Major: one vanished NAMED direction plus the unnamed one consumed is TWO.

    `(F, Cl, Br, implicitH)` becoming `(Cl, Br, I, At)` has two old positions with nothing to
    correspond to -- the departed fluorine and the stated-away hydrogen -- and two new positions the
    old frame never mentioned.  Assigning the two arrivals to the two free positions one way gives a
    parity of 1 and the other way 2, so there is no answer to compute and the sign dies.  A rule that
    counted only vanished NAMED directions saw one and kept the guess.

    Reachable in one batched edit, so it is measured through the edit rather than the probe.
    """
    m = MoleculeContainer()
    with m.edit():
        # two spare atoms in FRONT of the centre, so the arrivals can also permute the row
        sids = [m.add_atom(e, implicit_h=(1 if e == 'C' else None))
                for e in ('I', 'At', 'C', 'F', 'Cl', 'Br')]
        for s in sids[3:]:
            m.add_bond(sids[2], s, 1)
        m.set_parity(sids[2], 1)
    c = sids[2]
    assert m.unit_of(c)['refs'] == (sids[3], sids[4], sids[5], None)
    assert m.unit_of(c)['unnamed_mask'] == 0b1000, 'slot 3 is the centre\'s implicit hydrogen'
    with m.edit():
        m.delete_bond(c, sids[3])        # the fluorine leaves
        m.add_bond(c, sids[0], 1)        # the iodine arrives
        m.add_bond(c, sids[1], 1)        # ...and so does the astatine
        m.set_hydrogens(c, 0)            # ...into the hydrogen's place
    unit = m.unit_of(c)
    assert unit is not None, 'four heavy directions is still a unit'
    assert unit['refs'] == (sids[0], sids[1], sids[4], sids[5])
    assert unit['unnamed_mask'] == 0, 'and nothing unnamed is left to correspond with'
    assert m.parity_of(c) == 0, 'two arrivals into one row is a guess, so the sign dies'


def test_rebase_refuses_two_arrivals_into_one_cis_trans_terminal():
    """The same Major on a bond kind: `(methyl, implicitH)` becoming `(Cl, Br)` is two at once.

    The anchor terminal's whole pair is replaced in one edit -- the methyl leaves, a chlorine and a
    bromine arrive, one of them into the implicit hydrogen's position -- so the terminal's two
    directions have no correspondence to the old two.  A unit-wide count of vanished named directions
    reports only 1 here, which is not enough to drop the bit.
    """
    m = MoleculeContainer()
    with m.edit():
        # Cl0 and Br1 are spare; C2 is the anchor terminal, C3 the far one, C4 its methyl, C5 the
        # anchor's own methyl.
        sids = [m.add_atom(e, implicit_h=h)
                for e, h in (('Cl', None), ('Br', None), ('C', 1), ('C', 1),
                             ('C', None), ('C', None))]
        m.add_bond(sids[2], sids[3], 2)
        m.add_bond(sids[3], sids[4], 1)
        m.add_bond(sids[2], sids[5], 1)
        m.set_parity(sids[2], 1)
    anchor = sids[2]
    assert m.unit_of(anchor)['refs'] == (sids[5], None, sids[4], None)
    assert m.unit_of(anchor)['unnamed_mask'] == 0b1010, 'one implicit H on each terminal'
    with m.edit():
        m.delete_bond(anchor, sids[5])   # the methyl leaves
        m.add_bond(anchor, sids[0], 1)   # chlorine arrives
        m.add_bond(anchor, sids[1], 1)   # ...and bromine, into the hydrogen's position
        m.set_hydrogens(anchor, 0)
    unit = m.unit_of(anchor)
    assert unit is not None, 'a disubstituted terminal is still a terminal'
    assert unit['refs'] == (sids[0], sids[1], sids[4], None)
    assert m.parity_of(anchor) == 0, 'two arrivals into one terminal is a guess, so the sign dies'


def test_two_deleted_directions_drop_the_sign():
    """Two `RB_GONE` positions: the reachable spelling of the same refusal, via `delete_atom`.

    `delete_bond` leaves the departed neighbour in the molecule and the re-basing sees it as an atom
    that stopped being a direction; `delete_atom` takes the atom away entirely and the replay marks
    the position `RB_GONE` instead.  Both are leftovers and both count against the budget, but only
    this spelling exercises the `RB_GONE` arm -- and it was untested.
    """
    m = MoleculeContainer()
    with m.edit():
        sids = [m.add_atom(e, implicit_h=(1 if e == 'C' else None))
                for e in ('I', 'At', 'C', 'F', 'Cl', 'Br')]
        for s in sids[3:]:
            m.add_bond(sids[2], s, 1)
        m.set_parity(sids[2], 1)
    c = sids[2]
    assert m.unit_of(c)['refs'] == (sids[3], sids[4], sids[5], None)
    with m.edit():
        m.delete_atom(sids[3])           # the fluorine's ATOM is gone: RB_GONE
        m.delete_atom(sids[4])           # ...and the chlorine's: a second RB_GONE
        m.add_bond(c, sids[0], 1)
        m.add_bond(c, sids[1], 1)
    unit = m.unit_of(c)
    assert unit is not None, 'four directions again, so still a unit'
    assert unit['refs'] == (sids[0], sids[1], sids[5], None)
    assert m.parity_of(c) == 0, 'two deleted directions leave nothing to re-base along'


# --- one direction substituted: reachable, and the mirror of the hydrogen cases ------------

def test_implicitating_the_hydrogen_leaves_the_sign_alone():
    """Ruling F26's promise in the other direction: undrawing a hydrogen re-bases nothing.

    The drawn H is deleted and the count stated back in the same edit, so the unit keeps four
    directions and the H's position is now unnamed -- the identity permutation.  A rule that
    dropped a sign whenever a named direction left would lose the configuration of every molecule
    that gets its hydrogens undrawn, which is the same molecules `to_bytes` is asked to shrink.
    """
    m = MoleculeContainer()
    with m.edit():
        sids = [m.add_atom(e, implicit_h=(0 if e == 'C' else None))
                for e in ('C', 'F', 'Cl', 'Br')]
        h = m.add_atom('H')
        for s in sids[1:] + [h]:
            m.add_bond(sids[0], s, 1)
        m.set_parity(sids[0], 1)
    assert m.unit_of(sids[0])['refs'] == (sids[1], sids[2], sids[3], h)
    with m.edit():
        m.delete_atom(h)
        m.set_hydrogens(sids[0], 1)
    assert m.unit_of(sids[0])['refs'] == (sids[1], sids[2], sids[3], None)
    assert m.parity_of(sids[0]) == 1


def test_substituting_one_neighbour_rebases_the_sign():
    """A chlorine replaced by an iodine that sorts to the front of the row: one transposition.

    The iodine inherits the chlorine's POSITION in the frame -- that is the substitution rule --
    and then the row it lands in is `(I, F, Br, H)`, where the frame's own order reads
    `(F, I, Br, H)`.  One transposition, odd, so the stored bit flips; and the second assertion
    reads the new unit back in that frame order and gets the parity that was stored before.
    """
    m, sids = _chiral_methane_with_a_spare_iodine()
    i, c, f, cl, br = sids
    with m.edit():
        m.delete_bond(c, cl)
        m.add_bond(c, i, 1)
    assert m.unit_of(c)['refs'] == (i, f, br, None)
    assert m.parity_of(c) == 2
    assert m.translate_stereo(c, (f, i, br, None)) == 1


def test_undoing_a_substitution_restores_the_sign():
    """Two edits that cancel as graph operations must cancel as parity operations.

    An independent check on the arithmetic: the forward edit is measured to flip the bit and the
    reverse edit is a different permutation (a 4-cycle rather than the 3-cycle's inverse spelling),
    so getting back to 1 is the two agreeing rather than one of them being applied twice.
    """
    m, sids = _chiral_methane_with_a_spare_iodine()
    i, c = sids[0], sids[1]
    with m.edit():
        m.add_bond(c, i, 1)
        m.set_hydrogens(c, 0)
    assert m.parity_of(c) == 2
    with m.edit():
        m.delete_bond(c, i)
        m.set_hydrogens(c, 1)
    assert m.unit_of(c)['refs'] == (sids[2], sids[3], sids[4], None)
    assert m.parity_of(c) == 1, 'the sign came back, so the two re-basings are inverses'


def test_substituting_two_neighbours_at_once_drops_the_sign():
    """Reachable, and the drop side of the substitution rule."""
    m = MoleculeContainer()
    with m.edit():
        # two spare atoms in FRONT, so a substitution can also permute the row
        sids = [m.add_atom(e, implicit_h=(1 if e == 'C' else None))
                for e in ('I', 'Br', 'C', 'F', 'Cl', 'Br')]
        for s in sids[3:]:
            m.add_bond(sids[2], s, 1)
        m.set_parity(sids[2], 1)
    c = sids[2]
    assert m.unit_of(c)['refs'] == (sids[3], sids[4], sids[5], None)
    with m.edit():
        m.delete_bond(c, sids[3])
        m.delete_bond(c, sids[4])
        m.add_bond(c, sids[0], 1)
        m.add_bond(c, sids[1], 1)
    assert m.unit_of(c) is not None, 'still four directions, so still a unit'
    assert m.unit_of(c)['refs'] == (sids[0], sids[1], sids[5], None)
    assert m.parity_of(c) == 0, 'but no correspondence to re-base along'


# ---------------------------------------------------------------------------
# the parity segment
# ---------------------------------------------------------------------------

def test_an_edit_that_rebases_a_parity_writes_the_segment():
    """The re-base runs after the seal and against the sealed arena, so it is a second writer."""
    from chython.core._core import _parity_bytes

    m, sids = _chiral_methane_with_a_spare_iodine()
    i, c = sids[0], sids[1]
    with m.edit():
        m.add_bond(c, i, 1)
        m.set_hydrogens(c, 0)
    assert m.parity_of(c) == 2, 'the re-base flipped the sign'
    par = _parity_bytes(m)
    for n in m.atom_numbers:
        assert par[m.index_of(n)] == m.parity_of(n)


def test_a_parity_the_smiles_reader_stated_lands_in_the_segment():
    """The reader is a writer like any other, and its parities are storage too."""
    from chython.core._core import _parity_bytes, read_smiles

    m = read_smiles('C[C@H](N)C(=O)O')
    par = _parity_bytes(m)
    assert par, 'a molecule with a stated parity carries the segment'
    for n in m.atom_numbers:
        assert par[m.index_of(n)] == m.parity_of(n)


def test_a_molecule_with_no_stated_parity_carries_no_segment():
    """Absent is unset: the byte is not paid by the molecules that have nothing to say."""
    from chython.core._core import _parity_bytes, read_smiles

    assert _parity_bytes(read_smiles('CC(=O)Oc1ccccc1C(=O)O')) == b''


def test_request_parity_lays_out_the_segment_with_nothing_stated():
    """The door a post-seal writer needs: the segment exists and every byte is 0."""
    from chython.core._core import _parity_bytes, read_smiles

    m = read_smiles('CCO')
    assert _parity_bytes(m) == b''
    with m.edit() as e:
        e.request_parity()
    assert _parity_bytes(m) == bytes(len(m))


def test_set_parity_asks_for_the_segment_itself():
    """`request_parity` is only ever needed by a writer that states its parity after the seal."""
    from chython.core._core import _parity_bytes, read_smiles

    m = read_smiles('CCO')
    with m.edit():
        m.set_parity(m.atom_numbers[1], 2)
    par = _parity_bytes(m)
    assert par == b'\x00\x02\x00'
    assert par[1] == m.parity_of(m.atom_numbers[1])


def test_an_adopted_parity_is_dropped_from_the_segment_by_an_edit_that_drops_its_frame():
    """A version-4 buffer's parity is adopted at ingest, so the edit path finds it in the segment.

    The re-base is a writer like the replay is, and a dropped frame has to clear the byte -- a
    nonzero byte without a frame is the parity the molecule no longer holds.
    """
    from chython.core._core import _parity_bytes

    from .v4_fixtures import V4_SGROUP_STEREO_BYTES

    m = MoleculeContainer.from_bytes(V4_SGROUP_STEREO_BYTES)
    assert m.parity_of(2) == 2
    assert _parity_bytes(m)[m.index_of(2)] == 2, 'the flag parity was not adopted'
    with m.edit():
        m.delete_bond(2, 3)
    assert m.parity_of(2) == 0, 'the frame is gone, so the sign is'
    assert _parity_bytes(m)[m.index_of(2)] == 0, 'and so is the byte'


def test_an_adopted_segment_agrees_with_the_v4_records_atom_by_atom():
    """Every SEG_PARITY byte matches the parity the corresponding version-4 record stated.

    Adoption reads bit 7 and bit 1 from each v4 atom record and writes the three-state byte: 0 when
    bit 7 is clear, 2 if bit 1 is also set, 1 otherwise.  `request_parity` inside an edit must not
    overwrite the adopted values.
    """
    from chython.core._core import _parity_bytes

    from .v4_fixtures import V4_SGROUP_STEREO_BYTES

    start, _end = _atom_records(V4_SGROUP_STEREO_BYTES)
    m = MoleculeContainer.from_bytes(V4_SGROUP_STEREO_BYTES)
    with m.edit() as e:
        e.request_parity()
    par = _parity_bytes(m)
    assert par, 'the segment is gone'
    for n in m.atom_numbers:
        i = m.index_of(n)
        flags = V4_SGROUP_STEREO_BYTES[start + i * ATOM_RECORD + ATOM_FLAGS]
        v4_par = (2 if flags & 0x02 else 1) if flags & 0x80 else 0
        assert par[i] == v4_par, \
            'atom %d: segment byte is %d, v4 record stated %d' % (n, par[i], v4_par)
    assert par[m.index_of(2)] == 2, \
        'atom 2 stated parity 2 in the v4 record; it must reach the segment'


def test_every_writer_fills_the_parity_segment():
    """Three SMILES doors a parity comes in by: an `@` centre, a `/` double bond, and both at once.

    The pach v2 decoder is the fourth, in `test_pach.py` beside the corpus it needs.  All four write
    into a SEALED arena, so each has to name the segment before the seal, and each is a separate place
    to forget.
    """
    from chython.core._core import _parity_bytes

    for smi in ('C[C@H](N)C(=O)O', 'F/C=C/F', 'C[C@@H](O)/C=C\\C'):
        mol = read_smiles(smi)
        par = _parity_bytes(mol)
        assert par, '%s states a configuration and carries no parity segment' % smi
        for n in mol.atom_numbers:
            assert par[list(mol.atom_numbers).index(n)] == mol.parity_of(n), smi


def test_clearing_stereo_clears_the_parity_segment():
    """Both clear paths, and neither may leave a byte behind.

    The byte IS the molecule's stereo -- a molecule whose stereo was explicitly dropped would answer
    with the configuration it dropped.
    """
    from chython.core._core import _parity_bytes

    mol = read_smiles('C[C@H](N)C(=O)O')
    mol.clean_stereo()
    assert _parity_bytes(mol) in (b'', bytes(len(mol))), 'clean_stereo left a parity byte'

    # `validate_stereo` drops only the configurations it refuses, so the input has to state one it
    # can realize and one it cannot: atom 4 of 3-amino-2-methylbutan-2-ol carries two methyls, and
    # no geometry answers for it.  Both directions are asserted -- the refused byte goes to 0, the
    # realizable byte stays -- because a clear that took the whole segment would satisfy either one
    # alone, which is the shape a molecule with a single centre cannot tell apart.
    mol = read_smiles('C[C@H](N)[C@](C)(C)O')
    numbers = list(mol.atom_numbers)
    before = _parity_bytes(mol)
    assert before[numbers.index(2)] == 2 and before[numbers.index(4)] == 2, \
        'both centres state a parity going in'
    assert mol.validate_stereo() == [4], 'atom 4 is the one configuration nothing can realize'
    after = _parity_bytes(mol)
    assert len(after) == len(before), 'the segment must stay, holding the parity that survived'
    assert after[numbers.index(4)] == 0, 'validate_stereo left the refused parity byte behind'
    assert after[numbers.index(2)] == 2, 'validate_stereo took a byte it had no business taking'
    assert (mol.parity_of(4), mol.parity_of(2)) == (0, 2), 'and the accessor reads those same bytes'
