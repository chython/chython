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
"""The legacy pach codec, held against bytes chython 2 wrote.

THIS FILE IS AN ACCEPTANCE TEST AGAINST A FROZEN ORACLE, not a unit test of the decoder's internals.
Every expected value comes from `pach_corpus.py`'s three corpora, whose answers were produced by
chython 2's own unpacker before chython 2 left this tree.  Where chython 2's writer lost something,
the answers lost it too, so what is asserted here is that V3 recovers what the FORMAT carried -- not
a fidelity the format never had.  See `pach_corpus.py` for provenance and for why the corpora may
not be regenerated from this tree.

WHAT IS ASSERTED EXACTLY, and what is only measured:

  ASSERTED   Every record in all three corpora decodes, and its atoms (number, element, isotope,
             charge, radical, implicit hydrogens, degree), its bond set with orders, and its record
             length agree with the frozen answers.  Aromatic input stays aromatic.  Nothing in the
             codec calls kekule(), thiele() or standardize().  Garbage never crashes the decoder.
             The codec is idempotent: decode(encode(decode(x))) encodes to the same bytes.
  MEASURED   Byte-identity of a re-encoded record against the original.  It is NOT 100% and cannot
             be, for two enumerated reasons the module docstring of `_pach.pxi` states: chython 2
             wrote a neighbour list in insertion order where the arena's is index-ascending, and
             chython 2's float16 conversion truncates the mantissa where a rounding writer does not.
             The test pins the rate so a regression is visible, and prints it.

COORDINATES ARE COMPARED WITH A TOLERANCE, and the tolerance is a finding rather than a convenience.
The arena stores a display coordinate as a x10000 fixed-point integer; pach stores a float16.  Half
of 1e-4 is the largest error the arena can introduce, so `1.2001953125` (a float16) comes back as
`1.2002`.  That is the arena being unable to hold what the old format carried, and it is reported.
"""
import zlib

import pytest

from chython.core import MoleculeContainer, pach_dump, pach_load, pach_record_length, read_smiles
from .pach_corpus import V0_NATIVE_PATH, V0_PATH, V2_PATH, load_corpus


# The corpora are read once per session: five thousand records through gzip and json is a second of
# work and every test below wants all of them.
_CACHE = {}


def corpus(path):
    if path not in _CACHE:
        _CACHE[path] = load_corpus(path)
    return _CACHE[path]


def v2():
    return corpus(V2_PATH)


def v0():
    return corpus(V0_PATH)


def v0_native():
    return corpus(V0_NATIVE_PATH)


def every_record():
    """(corpus name, record name, bytes, answers) over all three corpora."""
    for tag, records in (('v2', v2()), ('v0', v0()), ('v0-native', v0_native())):
        for name, data, answers in records:
            yield tag, name, data, answers


# The fixed-point quantum of the arena's coordinate store, plus a hair for the division.
XY_TOLERANCE = 5.01e-5


def decoded_atoms(mol):
    """The molecule's atoms in the answers' own shape, minus the coordinates and the parity.

    Coordinates are compared separately because they need a tolerance; everything else here is exact.

    THE PARITY IS NOT IN THIS COMPARISON AT ALL, and the reason is a real difference between the two
    codebases rather than a tolerance.  chython 2 kept a cis/trans sign on the central BOND and left
    both terminal atoms' own stereo field None, so the answers' atom rows say "no configuration" for
    every atom of every double bond that has one.  The arena keeps that sign on the terminal atom it
    anchored the unit at, so `parity_of` says the opposite -- correctly.  Even the reduced question
    "is there a configuration here at all" therefore has two different right answers, and a column
    that must disagree on every cis/trans record is not a check.  The stereo tests below carry the
    parity instead, calibrated atom by atom against the V3 SMILES reader, which is a comparison in one
    frame between two independent readers of the same string.
    """
    out = []
    for n in mol.atom_numbers:
        isotope = mol.isotope_of(n)
        out.append([n, mol.element_of(n), isotope if isotope else None, mol.charge_of(n),
                    1 if mol.radical_of(n) else 0, mol.implicit_h_of(n), mol.degree_of(n)])
    return out


def answer_atoms(answers):
    """The frozen answers with the coordinates and the stereo field removed, as above."""
    out = []
    for a in answers['atoms']:
        row = list(a)
        del row[8]
        del row[7]
        del row[6]
        out.append(row)
    return out


def decoded_bonds(mol):
    out = set()
    for n in mol.atom_numbers:
        for m in mol.neighbors_of(n):
            out.add((min(n, m), max(n, m), mol.order_of(n, m)))
    return sorted(out)


def answer_bonds(answers):
    return sorted(tuple(b) for b in answers['bonds'])


# ------------------------------------------------------------------------------------------------
# The corpora themselves.  A silently truncated fixture would make every test below vacuous, so the
# counts are pinned here and nowhere else.
# ------------------------------------------------------------------------------------------------

def test_the_three_corpora_are_present_and_the_expected_size():
    assert len(v2()) == 2492
    assert len(v0()) == 2471
    assert len(v0_native()) == 236


def test_every_record_declares_a_version_this_codec_claims_to_read():
    for tag, name, data, _ in every_record():
        assert data[0] in (0, 2), '%s:%s declares version %d' % (tag, name, data[0])


# ------------------------------------------------------------------------------------------------
# Decoding.  The whole point of the branch: read what is stored.
# ------------------------------------------------------------------------------------------------

def test_every_record_of_every_corpus_decodes_without_a_problem():
    failed = []
    for tag, name, data, _ in every_record():
        mol, problems = pach_load(data, compressed=False)
        if mol is None or problems:
            failed.append((tag, name, problems))
    assert not failed, failed[:20]


def test_every_atom_field_agrees_with_chython_twos_own_answers():
    bad = []
    for tag, name, data, answers in every_record():
        mol, _ = pach_load(data, compressed=False)
        if mol is None:
            bad.append((tag, name, 'undecodable'))
            continue
        got, want = decoded_atoms(mol), answer_atoms(answers)
        if got != want:
            for g, w in zip(got, want):
                if g != w:
                    bad.append((tag, name, g, w))
                    break
    assert not bad, bad[:20]


def test_every_bond_agrees_with_chython_twos_own_answers():
    bad = []
    for tag, name, data, answers in every_record():
        mol, _ = pach_load(data, compressed=False)
        if mol is None:
            bad.append((tag, name, 'undecodable'))
            continue
        if decoded_bonds(mol) != answer_bonds(answers):
            bad.append((tag, name, decoded_bonds(mol)[:6], answer_bonds(answers)[:6]))
    assert not bad, bad[:10]


def test_coordinates_agree_within_the_arenas_fixed_point_quantum():
    bad = []
    for tag, name, data, answers in every_record():
        mol, _ = pach_load(data, compressed=False)
        if mol is None:
            continue
        for row in answers['atoms']:
            n, x, y = row[0], float(row[7]), float(row[8])
            if not mol.has_coordinates:
                if x or y:
                    bad.append((tag, name, n, 'coordinates dropped', x, y))
                continue
            gx, gy = mol.xy_of(n)
            if abs(gx - x) > XY_TOLERANCE or abs(gy - y) > XY_TOLERANCE:
                bad.append((tag, name, n, (gx, gy), (x, y)))
    assert not bad, bad[:20]


def test_the_record_length_a_reader_of_a_stream_needs_agrees_with_chython_two():
    bad = []
    for tag, name, data, answers in every_record():
        got = pach_record_length(data, compressed=False)
        if got != answers['size']:
            bad.append((tag, name, got, answers['size']))
    assert not bad, bad[:10]


def test_a_record_followed_by_trailing_bytes_still_decodes():
    """`pach_record_length` exists so a caller can walk a concatenated stream; the decoder must
    therefore ignore whatever follows the record it was given rather than reject the buffer."""
    name, data, answers = v2()[0]
    mol, problems = pach_load(data + b'garbage after the record', compressed=False)
    assert mol is not None and not problems, (name, problems)
    assert decoded_bonds(mol) == answer_bonds(answers)


def test_the_compressed_door_and_the_raw_door_agree():
    for name, data, _ in v2()[:200]:
        raw, _ = pach_load(data, compressed=False)
        packed, _ = pach_load(zlib.compress(data, 9), compressed=True)
        assert raw is not None and packed is not None, name
        assert decoded_atoms(raw) == decoded_atoms(packed), name
        assert decoded_bonds(raw) == decoded_bonds(packed), name


def test_unpack_is_the_versioned_front_door_and_reads_both_formats():
    """One method, two formats, dispatched on the first byte.

    A pach record's first byte is its format version, 0 or 2.  The arena's first byte is part of its
    magic and is neither.  So `unpack` can answer for both without a flag, and a caller holding a
    stored key does not have to know which decade it came from.
    """
    _, data, answers = v2()[0]
    from_pach = MoleculeContainer.unpack(data, compressed=False)
    assert decoded_bonds(from_pach) == answer_bonds(answers)
    from_arena = MoleculeContainer.unpack(from_pach.to_bytes(), compressed=False)
    assert decoded_bonds(from_arena) == answer_bonds(answers)
    assert from_arena.to_bytes() == from_pach.to_bytes()


# ------------------------------------------------------------------------------------------------
# IO IS NOT A MUTATOR OF REPRESENTATION.  Both directions.
# ------------------------------------------------------------------------------------------------

def test_an_aromatic_record_stays_aromatic_and_a_kekule_one_stays_kekule():
    """The corpus carries the same ring written both ways, deliberately.  A decoder that repaired
    either spelling into the other would pass every other test in this file.
    """
    seen = 0
    for name, data, answers in v2():
        if not name.startswith('ring_spelling:'):
            continue
        seen += 1
        mol, _ = pach_load(data, compressed=False)
        assert mol is not None, name
        want = {(b[0], b[1]): b[2] for b in answers['bonds']}
        aromatic_in = sum(1 for v in want.values() if v == 4)
        aromatic_out = sum(1 for n, m, o in decoded_bonds(mol) if o == 4)
        assert aromatic_in == aromatic_out, name
        if name.endswith(':aromatic'):
            assert aromatic_out, '%s lost its aromatic bonds' % name
        elif name.endswith(':kekule'):
            assert not aromatic_out, '%s gained aromatic bonds' % name
    assert seen >= 30, 'the ring-spelling half of the corpus is missing'


def test_the_codec_calls_no_repair_anywhere():
    """A structural check, because a repair call is invisible in a round trip that happens to be a
    fixed point of the repair.  `_pach.pxi` may not name kekule, thiele or standardize at all.
    """
    from pathlib import Path
    source = (Path(__file__).parent.parent / '_pach.pxi').read_text(encoding='utf-8')
    code = '\n'.join(line for line in source.split('\n')
                     if not line.lstrip().startswith('#'))
    for forbidden in ('kekule(', 'thiele(', 'standardize(', 'canonicalize('):
        assert forbidden not in code, 'the codec calls %s' % forbidden


def test_an_unknown_implicit_hydrogen_count_survives_as_unknown():
    """chython 1.42 left aromatic ring atoms' implicit hydrogen count unstated and the sentinel 7 is
    in the native v0 records.  None is not zero, in either direction.
    """
    seen = 0
    for name, data, answers in v0_native():
        mol, _ = pach_load(data, compressed=False)
        for row in answers['atoms']:
            if row[5] is None:
                seen += 1
                assert mol.implicit_h_of(row[0]) is None, (name, row[0])
    assert seen, 'the native v0 corpus no longer carries the unknown-hydrogen sentinel'


# ------------------------------------------------------------------------------------------------
# Encoding, and the round trip.
# ------------------------------------------------------------------------------------------------

def test_the_codec_is_idempotent_on_every_record():
    """decode -> encode -> decode -> encode is a fixed point.

    This is the round-trip property that CAN be 100%, and it is the one that matters for stored data
    written from V3 onwards: the second encoding has nothing left to normalise.
    """
    bad = []
    for tag, name, data, _ in every_record():
        mol, _ = pach_load(data, compressed=False)
        if mol is None:
            continue
        once = pach_dump(mol, compressed=False, version=2)
        again, problems = pach_load(once, compressed=False)
        if again is None:
            bad.append((tag, name, 'second decode failed', problems))
            continue
        twice = pach_dump(again, compressed=False, version=2)
        if once != twice:
            bad.append((tag, name, 'not a fixed point'))
    assert not bad, bad[:20]


def test_a_re_encoded_record_carries_the_same_answers():
    """Byte-identity is not achievable in general; MEANING-identity is, and this is the assertion
    that says the writer and the reader agree about the same molecule.
    """
    bad = []
    for tag, name, data, answers in every_record():
        mol, _ = pach_load(data, compressed=False)
        if mol is None:
            continue
        again, _ = pach_load(pach_dump(mol, compressed=False, version=2), compressed=False)
        if again is None:
            bad.append((tag, name, 'undecodable after encode'))
            continue
        if decoded_atoms(again) != answer_atoms(answers):
            bad.append((tag, name, 'atoms'))
        elif decoded_bonds(again) != answer_bonds(answers):
            bad.append((tag, name, 'bonds'))
    assert not bad, bad[:20]


def test_byte_identity_of_a_re_encoded_record_is_measured_and_explained(capsys):
    """The number, and the two reasons it is not the whole corpus.

    (1) chython 2 wrote each atom's neighbour list in its `_bonds` dict insertion order.  The arena's
        CSR is index-ascending.  For an acyclic molecule read from SMILES the two coincide; for a
        ring-closure neighbour they do not, and the connection table -- and with it the order block,
        which is consumed in connection-table order -- differs.
    (2) chython 2's `double_to_float16` TRUNCATES the mantissa (`<unsigned short> f`).  A V3 writer
        that truncated too would be enshrining a V2 defect in the V3 writer, which is forbidden, so
        it rounds; the two differ on the last mantissa bit of about half of all coordinates.

    A regression that broke the writer wholesale would drive this to zero, so the floor is asserted.
    """
    identical = total = 0
    for _, _, data, _ in every_record():
        mol, _ = pach_load(data, compressed=False)
        if mol is None:
            continue
        total += 1
        if pach_dump(mol, compressed=False, version=2) == data:
            identical += 1
    with capsys.disabled():
        print('\n  pach re-encode byte-identical: %d / %d (%.1f%%)'
              % (identical, total, 100. * identical / total))
    assert identical > total // 4, 'the writer stopped reproducing chython 2 bytes at all'


def test_packed_size_against_chython_two_is_measured(capsys):
    """pach against the arena, compressed and not, on the same molecules."""
    pach_raw = pach_zlib = arena = 0
    for _, data, _ in v2():
        mol, _ = pach_load(data, compressed=False)
        if mol is None:
            continue
        pach_raw += len(data)
        pach_zlib += len(zlib.compress(data, 9))
        arena += len(mol.to_bytes())
    with capsys.disabled():
        print('\n  pach raw %d, pach zlib %d, arena %d bytes over the v2 corpus'
              % (pach_raw, pach_zlib, arena))
    assert pach_raw and arena


def test_pack_and_unpack_are_the_containers_own_doors():
    _, data, answers = v2()[0]
    mol = MoleculeContainer.unpack(data, compressed=False)
    assert mol.pack(compressed=False, version=2) == pach_dump(mol, compressed=False, version=2)
    assert zlib.decompress(mol.pack(version=2)) == mol.pack(compressed=False, version=2)
    assert answer_bonds(answers) == decoded_bonds(MoleculeContainer.unpack(mol.pack(version=2)))


# ------------------------------------------------------------------------------------------------
# STEREO.  The one field whose meaning is not a bit but a bit RELATIVE TO AN ORDER, and the two
# orders are different, so the corpus carries the source string as the calibration channel.
# ------------------------------------------------------------------------------------------------

def stereo_records():
    """Records that came from a SMILES string and carry at least one configured parity."""
    for name, data, answers in v2():
        if 'smiles' not in answers:
            continue
        if any(row[6] != -1 for row in answers['atoms']) or answers['ct']:
            yield name, data, answers


def test_the_corpus_has_a_stereo_population_worth_calibrating_against():
    assert sum(1 for _ in stereo_records()) >= 40


def test_every_decoded_parity_says_what_the_smiles_reader_says(capsys):
    """The acceptance test for the stereo half of the codec.

    A parity is a statement about ONE ordering of an atom's directions.  chython 2's ordering was its
    `_bonds` insertion order with hydrogens excluded; the arena's is ruling F26's.  A bit alone
    therefore cannot say which arena parity is right, and the corpus answers cannot say either --
    they are chython 2's bit in chython 2's frame.  What CAN say is the source string: both readers
    number atoms 1..n in token order, so re-reading it with the V3 SMILES reader gives an independent
    answer for the same atom in the SAME frame the decoder must produce.

    Disagreements are reported per atom, not counted, because a single sign flip on a single centre
    is the whole failure mode this test exists to catch.

    ONLY ATOMS THE RECORD ACTUALLY CARRIES A CONFIGURATION FOR ARE COMPARED.  Where the record says
    nothing and the V3 reader says something, the two readers disagree about whether the atom is a
    stereocentre at all, and a decoder cannot invent what was never written: `C/C=C/[C@H](O)/C=C/C` is
    the case in the corpus -- its central carbon's two arms are constitutionally identical and differ
    only in the configuration of their double bonds, chython 2 does not perceive it and stored -1, and
    the V3 reader does.  Those atoms are counted and printed as a finding, not asserted on.
    """
    bad = []
    compared = 0
    unperceived = []
    for name, data, answers in stereo_records():
        mol, problems = pach_load(data, compressed=False)
        assert mol is not None, (name, problems)
        reference = read_smiles(answers['smiles'])
        if reference.atom_numbers != mol.atom_numbers:
            # Not a failure of the codec: the two readers disagree about the graph, so there is no
            # atom-by-atom comparison to make.  Reported by the count, which the assertion below
            # keeps honest.
            continue
        # What the RECORD carries for each atom: its own stereo column, plus membership of a cis/trans
        # entry, whose sign the record keeps on a bond and the arena on a terminal.
        carried = {row[0] for row in answers['atoms'] if row[6] != -1}
        for entry in answers['ct']:
            carried.add(entry[0])
            carried.add(entry[1])
        for n in mol.atom_numbers:
            want, got = reference.parity_of(n), mol.parity_of(n)
            if want == 0 and got == 0:
                continue
            if n not in carried:
                unperceived.append((name, answers['smiles'], n, want))
                continue
            compared += 1
            if want != got:
                bad.append((name, answers['smiles'], n, got, want))
    with capsys.disabled():
        print('\n  parities calibrated against the SMILES reader: %d, and %d centre(s) the V3 reader '
              'perceives and no record carries' % (compared, len(unperceived)))
    assert compared >= 40, 'the calibration channel is empty; the corpus lost its source strings'
    assert not bad, bad[:20]


def test_a_configured_parity_survives_a_pach_round_trip():
    """Encode then decode: the parity must come back, in the arena's frame, unchanged.  This holds
    even where the frames differ, because the writer states the sign in the frame it writes.
    """
    bad = []
    for name, data, _ in stereo_records():
        mol, _ = pach_load(data, compressed=False)
        again, problems = pach_load(pach_dump(mol, compressed=False, version=2), compressed=False)
        assert again is not None, (name, problems)
        for n in mol.atom_numbers:
            if mol.parity_of(n) != again.parity_of(n):
                bad.append((name, n, mol.parity_of(n), again.parity_of(n)))
    assert not bad, bad[:20]


def test_a_cis_trans_sign_survives_a_pach_round_trip():
    bad = []
    seen = 0
    for name, data, answers in v2():
        if not answers['ct']:
            continue
        seen += 1
        mol, _ = pach_load(data, compressed=False)
        again, _ = pach_load(pach_dump(mol, compressed=False, version=2), compressed=False)
        for n in mol.atom_numbers:
            unit = mol.unit_of(n)
            if unit is None or unit['parity'] == 0:
                continue
            other = again.unit_of(n)
            if other is None or other['parity'] != unit['parity']:
                bad.append((name, n, unit['kind'], unit['parity'],
                            None if other is None else other['parity']))
    assert seen, 'the corpus lost its cis/trans block'
    assert not bad, bad[:20]


# ------------------------------------------------------------------------------------------------
# INPUT BY DEFAULT IS GARBAGE.  A corrupt record may not take down a loop over forty thousand of
# them, so the decoder never raises; it returns what it could build and what was wrong with it.
# ------------------------------------------------------------------------------------------------

def test_truncation_at_every_length_of_a_real_record_never_crashes():
    _, data, _ = v2()[40]
    assert len(data) > 40
    for cut in range(len(data)):
        mol, problems = pach_load(data[:cut], compressed=False)
        assert mol is not None or problems, 'silence at cut %d' % cut


def test_truncation_of_every_record_in_the_corpus_never_crashes():
    """Not a duplicate of the above: one record has one atom block layout, and the failure modes
    live in the interaction of the counts with the buffer's length.
    """
    for _, _, data, _ in every_record():
        for cut in (1, 2, 3, 4, 5, 9, 13, len(data) // 2, len(data) - 1):
            if cut < 0 or cut >= len(data):
                continue
            mol, problems = pach_load(data[:cut], compressed=False)
            assert mol is not None or problems


def test_a_single_byte_flip_anywhere_never_crashes():
    _, data, _ = v2()[7]
    for i in range(len(data)):
        for bit in (0x01, 0x10, 0x80):
            broken = bytearray(data)
            broken[i] ^= bit
            mol, problems = pach_load(bytes(broken), compressed=False)
            assert mol is not None or problems


def test_the_empty_buffer_and_the_unknown_version_are_reported_not_raised():
    for data in (b'', b'\x00', b'\x02', b'\x01\x00\x00\x00', b'\x7f' * 64, bytes(64)):
        mol, problems = pach_load(data, compressed=False)
        assert mol is None or not problems or isinstance(problems, tuple)


def test_data_that_is_not_zlib_is_reported_not_raised():
    mol, problems = pach_load(b'not compressed at all', compressed=True)
    assert mol is None and problems


def test_the_answer_boundary_refuses_where_the_decoder_reports():
    """`pach_load` is the loop-safe door and never raises.  `unpack` is the answer boundary: a caller
    who asked for a molecule and cannot have one is told so.  Both doors, one decoder.
    """
    with pytest.raises(ValueError):
        MoleculeContainer.unpack(b'\x01\x00\x00\x00', compressed=False)
    with pytest.raises(ValueError):
        MoleculeContainer.unpack(b'', compressed=False)


#: `_dense_record`'s answer, built once: the two tests below decode the same 140 KB buffer.
_DENSE = []


def _dense_record():
    """A well-formed v2 record for a near-15-regular graph on the format's 4095 atoms.

    Nothing about the bytes is damaged; the graph is what `perceive_rings` refuses.  Degree is a
    nibble and an atom number is 12 bits, so v2 cannot state the complete graph `test_pach3.py` uses,
    and a random regular graph is the densest cycle space the format's own ceilings allow: 30672
    bonds, 140377 bytes, refused in about half a second.  Built by the configuration model -- fifteen
    slots per atom, shuffled once, paired -- so degree is 15 by construction rather than by rejection,
    and seeded so the record is the same on every run.
    """
    from random import Random

    if not _DENSE:
        n = 4095
        slots = [i for i in range(1, n + 1) for _ in range(15)]
        Random(20260907).shuffle(slots)
        edges = {(min(a, b), max(a, b)) for a, b in zip(slots[::2], slots[1::2]) if a != b}
        atoms = [(i, 0, 0, 6, b'\x00\x00', b'\x00\x00', 0xe0) for i in range(1, n + 1)]
        _DENSE.append(synthetic(atoms, [(a, b, 1) for a, b in sorted(edges)]))
    return _DENSE[0]


def test_a_graph_the_ring_perception_refuses_is_no_molecule_and_a_sentence():
    """The legacy decoder answers a derivation refusal the same way the v3/v4 one does.

    `rebuild_derived` ends every builder and `perceive_rings` is an answer boundary with a
    relevant-cycle prototype limit; the decode path is not one, so it reports.
    """
    mol, problems = pach_load(_dense_record(), compressed=False)
    assert mol is None
    assert any('could not be derived' in p and 'prototype limit' in p for p in problems)


def test_unpack_raises_where_that_record_reports():
    with pytest.raises(ValueError, match='prototype limit'):
        MoleculeContainer.unpack(_dense_record(), compressed=False)


# ------------------------------------------------------------------------------------------------
# chython 2 defects, converted into V3 acceptance tests.  None of them is fixed in chython 2 and
# none is reproduced by the V3 writer; each is a statement about what the DECODER does with bytes
# that already exist in stored data.
# ------------------------------------------------------------------------------------------------

def synthetic(atoms, bonds, cis_trans=(), version=2):
    """A pach record built by hand, so a field can hold a value chython 2's writer never wrote.

    `atoms` is [(number, stereo_nibble, isotope_shift, atomic_number, x_bytes, y_bytes, hcr), ...]
    with the neighbour count filled in from `bonds`; `bonds` is [(number, number, order), ...].
    """
    degree = {a[0]: 0 for a in atoms}
    for n, m, _ in bonds:
        degree[n] += 1
        degree[m] += 1
    out = bytearray()
    out.append(version)
    out.append(len(atoms) >> 4)
    out.append(((len(atoms) & 0x0f) << 4) | (len(cis_trans) >> 8))
    out.append(len(cis_trans) & 0xff)
    for n, stereo, isotope, element, xb, yb, hcr in atoms:
        out.append(n >> 4)
        out.append(((n & 0x0f) << 4) | degree[n])
        out.append(stereo | (isotope >> 1))
        out.append(((isotope & 1) << 7) | element)
        out += xb
        out += yb
        out.append(hcr)
    # Connection table in atom order, and the order block in the order the table is consumed.
    table, orders, seen = [], [], set()
    adjacency = {a[0]: [] for a in atoms}
    order_of = {}
    for n, m, o in bonds:
        adjacency[n].append(m)
        adjacency[m].append(n)
        order_of[(min(n, m), max(n, m))] = o
    for a in atoms:
        seen.add(a[0])
        for m in adjacency[a[0]]:
            table.append(m)
            if m not in seen:
                orders.append(order_of[(min(a[0], m), max(a[0], m))] - 1)
    for i in range(0, len(table), 2):
        n = table[i]
        m = table[i + 1] if i + 1 < len(table) else 0
        out.append(n >> 4)
        out.append(((n & 0x0f) << 4) | (m >> 8))
        out.append(m & 0xff)
    bits = ''.join('{:03b}'.format(o) for o in orders)
    if version == 0:
        while len(bits) % 15:
            bits += '0'
        for i in range(0, len(bits), 15):
            # THE PAD BIT IS THE TOP BIT of the 16-bit group, not the bottom one.  `_unpack_v0v2.pyx`
            # reads the group's first order out of the high nibble of the first byte (`a >> 4`, four
            # bits and one of them the pad) and its last out of the low three bits of the second, so
            # the free bit is bit 15 and the five values are right-aligned.  Writing it at the bottom
            # instead shifted every order of every synthetic v0 record by one position, which made the
            # unmasked-shift test below assert against a record no writer could produce.
            chunk = '0' + bits[i:i + 15]
            out.append(int(chunk[:8], 2))
            out.append(int(chunk[8:], 2))
    else:
        while len(bits) % 8:
            bits += '0'
        for i in range(0, len(bits), 8):
            out.append(int(bits[i:i + 8], 2))
    for n, m, s in cis_trans:
        out.append(n >> 4)
        out.append(((n & 0x0f) << 4) | (m >> 8))
        out.append(m & 0xff)
        out.append(1 if s else 0)
    return bytes(out)


def ethane(hcr_a=0xe0, hcr_b=0xe0, element=6, isotope=0, stereo=0, order=1, version=2):
    return synthetic([(1, stereo, isotope, element, b'\x00\x00', b'\x00\x00', hcr_a),
                      (2, 0, 0, 6, b'\x00\x00', b'\x00\x00', hcr_b)],
                     [(1, 2, order)], version=version)


def test_the_synthetic_builder_agrees_with_chython_twos_writer():
    """The hand builder is only evidence if it produces what the real writer produces.  Ethane with
    no coordinates and no hydrogen counts is in the corpus under a name, so compare against it.
    """
    mol, problems = pach_load(ethane(), compressed=False)
    assert mol is not None and not problems, problems
    assert mol.atom_numbers == [1, 2]
    assert decoded_bonds(mol) == [(1, 2, 1)]
    assert mol.implicit_h_of(1) is None and mol.implicit_h_of(2) is None


def test_an_implicit_hydrogen_count_above_six_cannot_be_written_and_is_not_pretended_otherwise():
    """chython 2's field is 3 bits with 7 reserved for "unknown", so counts 7..14 have no spelling.
    Its writer `<unsigned char> py_nan_int << 5` silently truncates 8 to 0.  The V3 WRITER refuses
    instead -- it does not enshrine the truncation -- and names the atom.
    """
    mol = MoleculeContainer()
    n = mol.add_atom(5, implicit_h=8)
    m = mol.add_atom(5)
    mol.add_bond(n, m, 1)
    with pytest.raises(ValueError, match='implicit'):
        pach_dump(mol, compressed=False, version=2)


def test_the_hydrogen_sentinel_is_a_sentinel_and_seven_is_not_a_count():
    mol, _ = pach_load(ethane(hcr_a=(7 << 5) | (4 << 1)), compressed=False)
    assert mol.implicit_h_of(1) is None
    mol, _ = pach_load(ethane(hcr_a=(6 << 5) | (4 << 1)), compressed=False)
    assert mol.implicit_h_of(1) == 6


def test_a_charge_the_field_can_hold_but_chython_two_could_not_is_decoded_and_reported():
    """The charge field is `charge + 4` in four bits, so it admits +5..+11 -- values chython 2's own
    `Element` refuses to hold, and which its writer therefore never wrote.  A record carrying one is
    corrupt, and the decoder's job is to say so rather than to raise inside somebody's loop.
    """
    mol, problems = pach_load(ethane(hcr_a=0xe0 | (13 << 1)), compressed=False)
    assert mol is not None
    assert any('charge' in p for p in problems), problems
    assert mol.charge_of(1) <= 8


def test_an_isotope_shift_out_of_the_five_bit_window_is_reported():
    """Shift 1..31 spells common-15..common+15 and 0 spells "unset".  chython 2's writer computes
    `isotope - common_isotope` with no range check, so an exotic isotope wraps into the wrong element
    mass or into 0.  The decoder reads what is there and reports an impossible mass number.
    """
    mol, problems = pach_load(ethane(element=1, isotope=1), compressed=False)
    assert mol is not None
    # Hydrogen's MDL reference mass is 1, so shift 1 asks for mass number 1 - 15 = -14.
    assert any('isotope' in p for p in problems), problems
    assert mol.isotope_of(1) == 0


def test_an_isotope_inside_the_window_round_trips():
    mol, _ = pach_load(ethane(element=6, isotope=17), compressed=False)
    assert mol.isotope_of(1) == 13
    again, _ = pach_load(pach_dump(mol, compressed=False, version=2), compressed=False)
    assert again.isotope_of(1) == 13


def test_a_bond_order_the_arena_does_not_have_is_reported_not_invented():
    """Three bits admit 0..7, so orders 1..8; the arena holds 1, 2, 3, 4 and 8.  Orders 5, 6 and 7
    are unrepresentable and chython 2's writer never wrote one -- but the v0 decoder can MANUFACTURE
    one, see below, so the case is reachable from stored bytes.
    """
    mol, problems = pach_load(ethane(order=6), compressed=False)
    assert mol is not None
    assert any('order' in p for p in problems), problems
    assert mol.order_of(1, 2) == 8


def test_the_v0_decoders_unmasked_shift_is_not_reproduced():
    """`_unpack_v0v2.pyx` reads the first order of a v0 pair as `a >> 4` with no mask, so a set pad
    bit -- bit 7 of the first byte, which nothing in the format defines -- yields an order of 9..16.
    chython 2 then built a Bond of order 9, which no chython 2 rule admits.  The V3 decoder masks to
    three bits, and this test is the statement that it does.
    """
    data = bytearray(ethane(order=2, version=0))
    # The order block is the last two bytes: five 3-bit orders and one pad bit.
    data[-2] |= 0x80
    mol, problems = pach_load(bytes(data), compressed=False)
    assert mol is not None, problems
    assert mol.order_of(1, 2) == 2, 'the pad bit leaked into the order'


def test_an_atomic_number_outside_the_periodic_table_is_refused_for_that_record_only():
    mol, problems = pach_load(ethane(element=0), compressed=False)
    assert mol is None and problems
    mol, problems = pach_load(ethane(element=119), compressed=False)
    assert mol is None and problems
    # ...and the next record in the loop is unaffected, which is the whole point.
    mol, problems = pach_load(ethane(), compressed=False)
    assert mol is not None and not problems


def test_an_asymmetric_connection_table_is_reported_where_chython_two_raised_keyerror():
    """chython 2's decoder does `py_bonds[py_m][py_n]` for a neighbour it has already seen, so a
    table where A lists B but B does not list A raises KeyError out of the C extension and takes the
    caller's loop with it.  The V3 decoder drops the half-bond and says so.
    """
    data = bytearray(ethane())
    # Atom 2's neighbour entry is the second 12-bit number of the connection table; point it at
    # atom 2 itself so that atom 1's claim is unreciprocated.
    table = 4 + 9 * 2
    data[table + 1] = (data[table + 1] & 0xf0) | 0
    data[table + 2] = 2
    mol, problems = pach_load(bytes(data), compressed=False)
    assert mol is not None, problems
    assert problems


def test_a_table_that_declares_more_pairs_than_its_own_bond_count_does_not_overrun():
    """A record whose degree sum implies four bonds and whose table names nine pairs.

    The bond count is a FUNCTION OF THE HEADER -- half the sum of the neighbour counts -- and a damaged
    header can make it disagree with the number of pairs the table goes on to name.  Sizing the decoder's
    edit buffer by that count wrote past its end for a record like this one and corrupted the heap; the
    crash surfaced thousands of records later, in the loop the decoder exists to keep alive.  So the
    buffer is sized by the number of entries that can be examined instead, and this is the record that
    settles it: a nine-armed star whose arms have had their neighbour counts patched to zero, which
    leaves every one of the nine pairs unreciprocated as well.
    """
    star = [(1, 0, 0, 6, b'\x00\x00', b'\x00\x00', 0xe0)]
    bonds = []
    for n in range(2, 11):
        star.append((n, 0, 0, 6, b'\x00\x00', b'\x00\x00', 0xe0))
        bonds.append((1, n, 1))
    data = bytearray(synthetic(star, bonds))
    for i in range(1, 10):                        # the arms, whose declared degree becomes 0
        at = 4 + 9 * i + 1
        data[at] &= 0xf0
    mol, problems = pach_load(bytes(data), compressed=False)
    assert mol is not None, problems
    assert len(mol) == 10
    assert sum(mol.degree_of(n) for n in mol.atom_numbers) == 0
    assert sum(1 for p in problems if 'not named back' in p) == 9, problems
    # and the recovered molecule is a usable one, not a half-built arena
    assert mol.rings_count == 0 and str(mol)


def test_a_record_whose_graph_has_no_usable_stereo_unit_table_reports_the_loss_and_stays_readable():
    """A configuration cannot be placed when the GRAPH the record describes has no anchor to put it on.

    Six carbons in a ring, two of the bonds cumulated, and the middle sp carbon declaring two implicit
    hydrogens -- a valence no writer would produce and exactly what a damaged record states.  The core's
    stereo derivation refuses such a graph outright ("two stereo units claim one anchor atom"), and it
    refuses it identically when the same graph is built with `add_atom`/`add_bond`, so this is a property
    of the graph and not of the decoder that happened to produce it.  The decoder therefore asks for the
    unit table inside a guard: the record still yields a usable molecule, the configuration it states is
    named as dropped, and nothing propagates into the caller's loop.
    """
    hcr = {0: 0x08, 1: 0x28, 2: 0x48, None: 0xe8}
    atoms, bonds = [], []
    for i, h in enumerate((0, 2, 1, None, None, None)):
        atoms.append((i + 1, 0xc0 if i == 0 else 0, 0, 6, b'\x00\x00', b'\x00\x00', hcr[h]))
        bonds.append((i + 1, (i + 1) % 6 + 1, 2 if i < 2 else 1))
    data = synthetic(atoms, bonds)

    mol, problems = pach_load(data, compressed=False)
    assert mol is not None, problems
    assert len(mol) == 6 and len(decoded_bonds(mol)) == 6
    assert any('stereo unit table' in p for p in problems), problems
    assert all(mol.parity_of(n) == 0 for n in mol.atom_numbers)

    # and the same graph built through the container API refuses in the same place, which is why the
    # decoder reports rather than repairs.  The refusal is a finding about the graph, not about pach.
    rebuilt = MoleculeContainer()
    ids = [rebuilt.add_atom(6, implicit_h=h) for h in (0, 2, 1, None, None, None)]
    for i in range(6):
        rebuilt.add_bond(ids[i], ids[(i + 1) % 6], 2 if i < 2 else 1)
    with pytest.raises(RuntimeError, match='anchor'):
        str(rebuilt)


def test_a_duplicate_atom_number_is_refused_for_that_record_only():
    data = synthetic([(1, 0, 0, 6, b'\x00\x00', b'\x00\x00', 0xe0),
                      (1, 0, 0, 6, b'\x00\x00', b'\x00\x00', 0xe0)], [])
    mol, problems = pach_load(data, compressed=False)
    assert mol is None and any('duplicate' in p for p in problems), problems


def test_the_tetrahedron_and_allene_nibbles_are_both_read_as_chython_two_read_them():
    """chython 2's decoder collapses the two 2-bit fields with a catch-all `else: True`, so the
    tetrahedron/allene distinction the WRITER made is not in the answers.  The V3 decoder cannot
    recover a distinction the reference never had, so it does the same thing: the sign is read, and
    which unit it belongs to is decided by the graph.
    """
    for nibble, expect in ((0x80, 1), (0xc0, 2), (0x20, 1), (0x30, 2)):
        mol, problems = pach_load(ethane(stereo=nibble), compressed=False)
        assert mol is not None, problems
        # Ethane has no stereo unit, so the sign has nowhere to go and that is reported rather than
        # written to an atom that cannot hold it.
        assert mol.parity_of(1) == 0
        assert problems, (nibble, expect)


def test_a_zero_length_pack_of_a_molecule_with_no_bonds_round_trips():
    """chython 2 refused to pack a molecule with no bonds (`check=True` raises on empty `_bonds`),
    which made a lone sodium cation unstorable.  The format itself has no such restriction: the
    order block is simply zero bytes long.  V3 writes it.
    """
    mol = MoleculeContainer()
    mol.add_atom(11, charge=1)
    data = pach_dump(mol, compressed=False, version=2)
    again, problems = pach_load(data, compressed=False)
    assert again is not None and not problems, problems
    assert again.atom_numbers == mol.atom_numbers
    assert again.charge_of(mol.atom_numbers[0]) == 1


# ------------------------------------------------------------------------------------------------
# What the format cannot carry.  Every one of these is data the arena holds and pach has no slot
# for, and the writer's contract is that it says so instead of dropping it.
# ------------------------------------------------------------------------------------------------

def small():
    mol = MoleculeContainer()
    a = mol.add_atom(6)
    b = mol.add_atom(8)
    mol.add_bond(a, b, 1)
    return mol, a, b


def test_a_map_number_has_no_slot_and_is_not_dropped_silently():
    mol, a, _ = small()
    mol.set_map_number(a, 7)
    with pytest.raises(ValueError, match='map_number'):
        pach_dump(mol, compressed=False, version=2)
    assert pach_dump(mol, compressed=False, drop=['map_number'], version=2)
    assert pach_dump(mol, compressed=False, drop='*', version=2)


def test_a_title_has_no_slot_and_is_not_dropped_silently():
    mol, _, _ = small()
    mol.set_title('a public compound')
    with pytest.raises(ValueError, match='title'):
        pach_dump(mol, compressed=False, version=2)
    assert pach_dump(mol, compressed=False, drop=['title'], version=2)


def test_an_unknown_drop_name_is_refused_rather_than_ignored():
    mol, a, _ = small()
    mol.set_map_number(a, 7)
    with pytest.raises(ValueError, match='drop'):
        pach_dump(mol, compressed=False, drop=['mapnumbers'])


def test_an_atom_number_the_twelve_bit_field_cannot_hold_is_refused():
    mol = MoleculeContainer()
    ids = [mol.add_atom(6) for _ in range(3)]
    mol.add_bond(ids[0], ids[1], 1)
    mol.add_bond(ids[1], ids[2], 1)
    mol.remap({ids[2]: 4096})
    with pytest.raises(ValueError, match='4095'):
        pach_dump(mol, compressed=False, version=2)


def test_a_degree_above_fifteen_is_refused():
    mol = MoleculeContainer()
    centre = mol.add_atom(1)
    for _ in range(16):
        mol.add_bond(centre, mol.add_atom(6), 1)
    with pytest.raises(ValueError, match='neighbo'):
        pach_dump(mol, compressed=False, version=2)


def test_pack_refuses_a_non_empty_meta():
    """Record metadata is the eighth thing pach has no field for, so it is refused by name."""
    mol = read_smiles('CCO')
    mol.meta['boiling_point'] = '78.37'
    with pytest.raises(ValueError, match='metadata key'):
        mol.pack(version=2)
    assert MoleculeContainer.unpack(mol.pack(drop=['meta'], version=2)) == mol


def test_pack_is_silent_about_an_untouched_meta():
    """Reading `mol.meta` creates the dict; an empty one is not a loss."""
    mol = read_smiles('CCO')
    assert mol.meta == {}
    mol.pack(version=2)


def test_a_reaction_waiver_reaches_its_components():
    from chython.core import ReactionContainer

    mol = read_smiles('CCO')
    mol.meta['k'] = 'v'
    rxn = ReactionContainer([mol], [read_smiles('CC=O')])
    with pytest.raises(ValueError, match='metadata key'):
        rxn.pack()
    rxn.pack(drop=['meta'])


def test_a_v2_record_lands_its_parities_in_the_segment():
    """The legacy pach decoder writes after the seal too, so it names the segment at build time."""
    from chython.core._core import _parity_bytes

    checked = 0
    for name, data, answers in v2():
        mol, _ = pach_load(data, compressed=False)
        if any(mol.parity_of(n) for n in mol.atom_numbers):
            par = _parity_bytes(mol)
            assert par, 'a decoded v2 record states a parity and carries no segment: %s' % name
            for i, n in enumerate(mol.atom_numbers):
                assert par[i] == mol.parity_of(n), name
            checked += 1
    assert checked, 'the v2 corpus states no parity; the test is not exercising the path'
