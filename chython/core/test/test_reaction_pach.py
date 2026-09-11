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
"""The reaction-level pach codec: four header bytes over the molecule codec.

WHAT THE FORMAT IS, and it is four bytes and no more:

    byte 0        0x05 for the current format (version 5), or 0x01 for the legacy one (version 1)
    byte 1        reactant count,  uint8
    byte 2        AGENT count,     uint8    <- the middle field, not the last one
    byte 3        product count,   uint8
    the rest      the molecules' pach records, uncompressed, concatenated in the order
                  reactants -> agents -> products, i.e. `ReactionContainer.molecules()`

At version 5 the body records are version 3 or version 4 molecule records (third-generation, with a
map-number block); at version 1 they are version 0 or version 2 (one number field per atom, no map
block).  Bytes 1-3 and the concatenated body structure are identical between the two.

There is no length field anywhere: a molecule record's length is a function of its own header, so a
reader walks the stream with `pach_record_length`.  That is why a wrong count and a truncated buffer
are the same class of failure and both are tested below.

WHY THE MAP NUMBER TESTS ARE THE POINT OF THIS FILE.  chython 2 had no separate map-number field on
an atom: the atom's NUMBER was its mapping.  Version 1's reaction writer moved map numbers into pach's
12-bit atom-number field, and version 1's reader takes them out of it again.  Version 5 has a map
block in each molecule record, so neither step happens.  A codec that lost the mapping would pass
every test about atoms and bonds and still be useless, because a mapping is the only reason a reaction
gets packed rather than written as a string.

THE FIXTURE IS THE SPECIFICATION, not the round trip.  `reaction_pach_v2_corpus.bin.gz` was written by
an installed chython 2.24 -- see `reaction_pach_corpus.py` for provenance.  A writer and a reader that
agree with each other and not with chython 2 pass every round-trip assertion in this file and still
fail the only requirement that matters, so the fixture tests are the ones to read first.

PUBLIC COMPOUNDS ONLY, here and in the fixture: textbook amidation, Suzuki, nitration, esterification.
"""
import gzip
import random
import zlib
from struct import pack as struct_pack

import pytest

from chython.core import (MoleculeContainer, ReactionContainer, read_smiles, reaction_pach_dump,
                          reaction_pach_load, reaction)
from .reaction_pach_corpus import V2_PATH, load_corpus


# --------------------------------------------------------------------------------------------------
# the reactions the round-trip tests use.  All public: an amide coupling, a Suzuki, a nitration.
# --------------------------------------------------------------------------------------------------

AMIDATION = ('[CH3:1][C:2](=[O:3])[OH:4].[NH2:5][CH3:6]'
             '>[CH3:10][CH2:11][OH:12]'
             '>[CH3:1][C:2](=[O:3])[NH:5][CH3:6].[OH2:4]')


def _reaction(spec):
    """Build a reaction from a `reactants>agents>products` string with the V3 SMILES reader.

    Written out rather than reached for through `chython.formats`, because `chython/core/` may not
    import anything above itself -- see `test_no_chython_two_imports.py`.
    """
    sides = spec.split('>')
    assert len(sides) == 3
    out = []
    for side in sides:
        out.append([read_smiles(s) for s in side.split('.') if s])
    return ReactionContainer(out[0], out[2], out[1])


def _snapshot(rxn):
    """Everything a round trip has to preserve, as plain data.

    Keyed by MAP NUMBER where there is one, because that is the label the format carries and the one
    that has to survive; by stable id where there is not, since an unmapped molecule's atoms still
    have to come back with the same elements on the same graph.  A VERSION 1 record's number field
    holds the stable id and its reader reports that as the mapping, so the two keys coincide there;
    a version 5 record carries the field and an unmapped molecule comes back with 0, so the stable id
    is the key on both sides.  Either way the snapshot compares.
    """
    out = []
    for side in (rxn.reactants, rxn.agents, rxn.products):
        molecules = []
        for mol in side:
            key = {a.n: (a.map_number or a.n) for a in mol.atoms()}
            atoms = {key[a.n]: (a.element, a.charge, a.isotope, a.is_radical, a.implicit_h,
                                a.stereo) for a in mol.atoms()}
            assert len(atoms) == len(mol), 'the snapshot key collided; it cannot compare anything'
            bonds = sorted(tuple(sorted((key[b.n], key[b.m]))) + (b.order,) for b in mol.bonds())
            molecules.append((atoms, bonds))
        out.append(molecules)
    return out


# --------------------------------------------------------------------------------------------------
# round trips
# --------------------------------------------------------------------------------------------------

def test_a_mapped_three_sided_reaction_round_trips_completely():
    rxn = _reaction(AMIDATION)
    back = ReactionContainer.unpack(rxn.pack())
    assert (len(back.reactants), len(back.agents), len(back.products)) == (2, 1, 2)
    assert _snapshot(back) == _snapshot(rxn)


def test_map_numbers_survive_and_are_the_ones_that_were_written():
    """The assertion the whole codec exists for, spelled out rather than folded into a snapshot."""
    rxn = _reaction(AMIDATION)
    back = ReactionContainer.unpack(rxn.pack())
    assert [sorted(a.map_number for a in m.atoms()) for m in back.reactants] == [[1, 2, 3, 4], [5, 6]]
    assert [sorted(a.map_number for a in m.atoms()) for m in back.agents] == [[10, 11, 12]]
    assert [sorted(a.map_number for a in m.atoms()) for m in back.products] == [[1, 2, 3, 5, 6], [4]]
    # and the mapping still relates the two sides: the acid's carbonyl carbon is the amide's
    acid = {a.map_number: a for m in back.reactants for a in m.atoms()}
    amide = {a.map_number: a for m in back.products for a in m.atoms()}
    assert acid[2].element == amide[2].element == 6
    assert acid[5].element == amide[5].element == 7


def test_charges_isotopes_and_radicals_survive():
    rxn = _reaction('[13CH3:1][C:2](=[O:3])[O-:4].[Na+:5]>>[13CH3:1][C:2](=[O:3])[O-:4].[Na+:5]')
    back = ReactionContainer.unpack(rxn.pack())
    atoms = {a.map_number: a for m in back.reactants for a in m.atoms()}
    assert atoms[1].isotope == 13
    assert atoms[4].charge == -1
    assert atoms[5].charge == 1
    radical = _reaction('[CH3:1][CH2:2][O:3][O:4]>>[CH3:1][CH2:2][O:3].[O:4]')
    with radical.products[1].edit():
        radical.products[1].set_radical(next(iter(radical.products[1].atom_numbers)), True)
    back = ReactionContainer.unpack(radical.pack())
    assert next(iter(back.products[1].atoms())).is_radical


def test_a_reaction_carrying_stereo_round_trips_through_the_wrapper():
    """(S)-lactic acid esterified: a tetrahedral parity has to come back on the right atom."""
    rxn = _reaction('[CH3:1][C@H:2]([OH:3])[C:4](=[O:5])[OH:6].[CH3:7][OH:8]'
                    '>>[CH3:1][C@H:2]([OH:3])[C:4](=[O:5])[O:6][CH3:7].[OH2:8]')
    before = _snapshot(rxn)
    back = ReactionContainer.unpack(rxn.pack())
    assert _snapshot(back) == before
    stereo = {a.map_number: a.stereo for m in back.reactants for a in m.atoms()}
    assert stereo[2] is not None
    # and the parity is the one the source string states, not merely "some parity"
    assert stereo[2] == {a.map_number: a.stereo for m in rxn.reactants for a in m.atoms()}[2]


def test_an_aromatic_ring_stays_aromatic_across_the_wrapper():
    """Nitration of toluene.  The codec repairs nothing in either direction, so aromatic bonds go in
    aromatic and come back aromatic -- no kekulisation on the way through."""
    rxn = _reaction('[cH:1]1[cH:2][cH:3][c:4]([CH3:5])[cH:6][cH:7]1'
                    '>>[c:1]1([N+:8](=[O:9])[O-:10])[cH:2][cH:3][c:4]([CH3:5])[cH:6][cH:7]1')
    back = ReactionContainer.unpack(rxn.pack())
    carbons = {1, 2, 3, 4, 6, 7}          # the ring; 5 is the methyl and 8..10 the nitro group
    for side in (back.reactants, back.products):
        maps = {a.n: a.map_number for a in side[0].atoms()}
        ring = [b.order for b in side[0].bonds() if maps[b.n] in carbons and maps[b.m] in carbons]
        assert len(ring) == 6 and set(ring) == {4}, ring
    assert _snapshot(back) == _snapshot(rxn)


@pytest.mark.parametrize('spec, counts', [
    ('CCO.CC(=O)O>>CC(=O)OCC.O', (2, 0, 2)),          # no agents
    ('CCO>CC(=O)O>', (1, 1, 0)),                      # no products
    ('>CC(=O)O>CCO', (0, 1, 1)),                      # no reactants
    ('>>', (0, 0, 0)),                                # nothing at all
])
def test_empty_sides_round_trip(spec, counts):
    """AND ONE OF THESE PINS A chython 2 DIVERGENCE RATHER THAN REPRODUCING IT.  chython 2's writer
    wrote `(1, 1, 0, 0)` for a record with no products and its own reader then read that buffer as
    `CCO>>CCO`: `molecules[-products:]` with `products == 0` is `molecules[0:]`, the whole list, so
    the reactant reappeared as a product.  Not patched there -- V3 reads the header it was given."""
    rxn = _reaction(spec)
    packed = rxn.pack()
    # the header spells the sides in the SAME order as `counts`: reactants, agents, products
    assert zlib.decompress(packed)[1:4] == bytes(counts)
    back = ReactionContainer.unpack(packed)
    assert (len(back.reactants), len(back.agents), len(back.products)) == counts
    assert _snapshot(back) == _snapshot(rxn)


def test_packing_does_not_mutate_the_reaction_it_was_handed():
    """Version 1's writer relabels a molecule onto its map numbers to spell the record, and it does
    that on a COPY.  A serialiser that renumbered its caller's molecules would be a mutating getter;
    version 5 has nothing to relabel and is packed here so neither door can grow one."""
    rxn = _reaction(AMIDATION)
    before = [(m.atom_numbers, {a.n: a.map_number for a in m.atoms()}) for m in rxn.molecules()]
    rxn.pack()
    rxn.pack(version=1)
    after = [(m.atom_numbers, {a.n: a.map_number for a in m.atoms()}) for m in rxn.molecules()]
    assert after == before


def test_the_codec_is_idempotent():
    """A second round trip changes nothing, which is what makes a stored buffer a stable key."""
    rxn = _reaction(AMIDATION)
    once = rxn.pack(compressed=False)
    twice = ReactionContainer.unpack(once, compressed=False).pack(compressed=False)
    assert once == twice


# --------------------------------------------------------------------------------------------------
# the layout, asserted directly
# --------------------------------------------------------------------------------------------------

def test_the_header_is_version_reactants_agents_products():
    rxn = _reaction(AMIDATION)
    raw = rxn.pack(compressed=False, version=1)
    assert raw[0] == 1
    assert (raw[1], raw[2], raw[3]) == (2, 1, 2)
    # and the body is the molecule records, in `molecules()` order, byte for byte
    shift = 4
    for mol in rxn.molecules():
        n = len(mol)
        assert raw[shift] == 2, 'each body record is a pach version 2 molecule'
        assert ((raw[shift + 1] << 4) | (raw[shift + 2] >> 4)) == n
        shift += _record_length(raw, shift)
    assert shift == len(raw), 'the records exactly fill the buffer'


def _record_length(raw, shift):
    from chython.core import pach_record_length
    return pach_record_length(bytes(raw[shift:]), compressed=False)


def test_pack_is_zlib_compressed_by_default_and_raw_on_request():
    rxn = _reaction(AMIDATION)
    raw = rxn.pack(compressed=False)
    assert zlib.decompress(rxn.pack()) == raw
    assert raw[0] == 5


def test_unpack_sniffs_compression_and_can_be_told():
    rxn = _reaction(AMIDATION)
    raw = rxn.pack(compressed=False)
    packed = rxn.pack()
    assert len(ReactionContainer.unpack(raw)) == 5
    assert len(ReactionContainer.unpack(packed)) == 5
    assert len(ReactionContainer.unpack(raw, compressed=False)) == 5
    assert len(ReactionContainer.unpack(packed, compressed=True)) == 5
    with pytest.raises(ValueError):
        ReactionContainer.unpack(raw, compressed=True)
    with pytest.raises(ValueError):
        ReactionContainer.unpack(packed, compressed=False)


def test_pack_len_reports_atom_counts_side_by_side():
    rxn = _reaction(AMIDATION)
    assert ReactionContainer.pack_len(rxn.pack()) == ((4, 2), (3,), (5, 1))
    # and it does not confuse an empty side for a full one, which chython 2's version did
    empty = _reaction('CCO>CC(=O)O>')
    assert ReactionContainer.pack_len(empty.pack()) == ((3,), (4,), ())


# --------------------------------------------------------------------------------------------------
# garbage.  Input is garbage by default: every one of these raises a clear Python exception, and the
# process survives all of them.
# --------------------------------------------------------------------------------------------------

def test_a_truncated_buffer_raises_and_does_not_return_half_a_reaction():
    raw = _reaction(AMIDATION).pack(compressed=False)
    for cut in range(1, len(raw)):
        with pytest.raises(ValueError):
            ReactionContainer.unpack(raw[:cut], compressed=False)


def test_a_wrong_header_byte_raises():
    raw = bytearray(_reaction(AMIDATION).pack(compressed=False))
    for wrong in (0, 2, 3, 0x33, 0xff):
        raw[0] = wrong
        with pytest.raises(ValueError) as err:
            ReactionContainer.unpack(bytes(raw), compressed=False)
        assert 'version' in str(err.value) or 'not a reaction' in str(err.value)


def test_a_count_that_overruns_the_data_raises():
    raw = bytearray(_reaction(AMIDATION).pack(compressed=False))
    raw[1] = 9                     # nine reactants, five molecules present
    with pytest.raises(ValueError):
        ReactionContainer.unpack(bytes(raw), compressed=False)
    raw[1], raw[3] = 2, 200        # and on the far side of the stream
    with pytest.raises(ValueError):
        ReactionContainer.unpack(bytes(raw), compressed=False)


def test_the_empty_buffer_and_random_noise_raise():
    with pytest.raises(ValueError):
        ReactionContainer.unpack(b'')
    with pytest.raises(ValueError):
        ReactionContainer.unpack(b'\x01')
    with pytest.raises(ValueError):
        ReactionContainer.unpack(bytes(range(64)))


def test_reaction_pach_load_never_raises_on_any_of_them():
    """The loop-safe door, mirroring `pach_load`: a store of forty thousand records must not be
    stopped by one of them, so this reports instead of raising."""
    raw = _reaction(AMIDATION).pack(compressed=False)
    cases = [b'', b'\x01', b'\x01\x02\x00\x02', bytes(range(64)), raw[:20],
             b'\x07' + raw[1:], zlib.compress(b'nonsense'), raw]
    for case in cases:
        rxn, problems = reaction_pach_load(case)
        assert isinstance(problems, list)
        assert all(isinstance(p, str) for p in problems)
        if rxn is None:
            assert problems, case
        else:
            assert isinstance(rxn, ReactionContainer)
    # the last case is the good one and it comes back clean
    rxn, problems = reaction_pach_load(raw)
    assert problems == [] and len(rxn) == 5


def test_a_component_whose_graph_the_ring_perception_refuses_is_no_reaction_and_a_sentence():
    """The molecule-level report reaches the reaction door unchanged.

    `reaction_pach_load` inherits `pach_load`'s answer for each component, so a well-formed record
    whose graph trips `perceive_rings`' resource limits costs the reaction and not the process.  The
    forger lives beside the molecule test because the record is a molecule record.
    """
    from .test_pach3 import complete_graph_record

    raw = bytes([5, 1, 0, 0]) + complete_graph_record()
    rxn, problems = reaction_pach_load(raw, compressed=False)
    assert rxn is None
    assert any('could not be derived' in p and 'prototype limit' in p for p in problems)


def test_a_seeded_fuzz_of_the_whole_buffer_never_crashes_the_interpreter():
    """4000 mutants of a real record: bytes replaced at random and the tail cut off at random.

    THE POINT IS THE ABSENCE OF A SEGFAULT, not any particular answer.  Everything under here is
    reading a C extension's memory off lengths the buffer itself declares, so a wrong length is a
    wrong pointer; the only proof that it cannot be is to try.  Seeded, so a failure is reproducible.
    A larger run -- 60000 mutants -- was done by hand and is not committed, because a test that takes
    a minute to say nothing new is a test people learn to skip.
    """
    raw = _reaction(AMIDATION).pack(compressed=False)
    rand = random.Random(20260903)
    read = refused = 0
    for _ in range(4000):
        mutant = bytearray(raw)
        for _ in range(rand.randint(1, 6)):
            mutant[rand.randrange(len(mutant))] = rand.randrange(256)
        if rand.random() < 0.3:
            del mutant[rand.randrange(len(mutant)):]
        rxn, problems = reaction_pach_load(bytes(mutant), compressed=False)
        if rxn is None:
            assert problems
            refused += 1
        else:
            # whatever it read, it read the WHOLE declared reaction and not part of one
            assert len(rxn) == mutant[1] + mutant[2] + mutant[3]
            read += 1
    assert read and refused, 'the fuzz degenerated: %d read, %d refused' % (read, refused)


def test_every_single_byte_flip_in_a_header_either_reads_or_raises():
    """A property test over the four header bytes: 1024 buffers, no crash, no partial reaction."""
    raw = _reaction(AMIDATION).pack(compressed=False)
    for i in range(4):
        for value in range(256):
            mutant = bytearray(raw)
            mutant[i] = value
            rxn, problems = reaction_pach_load(bytes(mutant), compressed=False)
            if rxn is not None:
                assert len(rxn) == mutant[1] + mutant[2] + mutant[3]


# --------------------------------------------------------------------------------------------------
# more than 255 molecules on a side
# --------------------------------------------------------------------------------------------------

def test_more_than_255_molecules_on_a_side_is_refused_by_name():
    """chython 2 did `bytearray((1, len(reactants), ...))`, which raises `ValueError: byte must be in
    range(0, 256)` -- so it never silently truncated and no stored buffer can hold such a record.
    V3 therefore refuses too, and says which side and how many rather than talking about bytes."""
    water = read_smiles('O')
    for side, name in ((0, 'reactants'), (1, 'products'), (2, 'agents')):
        sides = [[], [], []]
        sides[side] = [water.copy() for _ in range(256)]
        sides[(side + 1) % 3] = [water.copy()]
        rxn = ReactionContainer(sides[0], sides[1], sides[2])
        with pytest.raises(ValueError) as err:
            rxn.pack()
        assert name in str(err.value) and '256' in str(err.value)
    # 255 is fine, and reads back
    rxn = ReactionContainer([water.copy() for _ in range(255)], [water.copy()])
    back = ReactionContainer.unpack(rxn.pack())
    assert (len(back.reactants), len(back.products)) == (255, 1)


# --------------------------------------------------------------------------------------------------
# the writer refuses to lose things quietly, exactly as the molecule writer does
# --------------------------------------------------------------------------------------------------

def test_the_writer_refuses_reaction_meta_and_title_and_takes_a_waiver():
    rxn = _reaction(AMIDATION)
    rxn.meta['SOURCE'] = 'a textbook'
    with pytest.raises(ValueError) as err:
        rxn.pack()
    assert 'meta' in str(err.value)
    assert ReactionContainer.unpack(rxn.pack(drop=['meta'])).meta == {}

    rxn = _reaction(AMIDATION)
    rxn.set_title(b'acetamide from acetic acid')
    with pytest.raises(ValueError) as err:
        rxn.pack()
    assert 'title' in str(err.value)
    assert ReactionContainer.unpack(rxn.pack(drop=['title'])).title == ''


def test_a_molecule_level_refusal_reaches_the_caller_and_names_the_field():
    rxn = _reaction(AMIDATION)
    rxn.reactants[0].set_title(b'acetic acid')
    with pytest.raises(ValueError) as err:
        rxn.pack()
    assert 'title' in str(err.value)
    assert len(ReactionContainer.unpack(rxn.pack(drop=['title']))) == 5


def test_an_unrecognised_drop_name_is_refused_rather_than_ignored():
    rxn = _reaction(AMIDATION)
    with pytest.raises(ValueError) as err:
        rxn.pack(drop=['mapping'])
    assert 'mapping' in str(err.value)
    assert len(rxn.pack(drop='*')) > 4


def test_a_partially_mapped_molecule_is_refused_rather_than_half_written():
    """A molecule with a mapping on some atoms and not others has no honest spelling in a format with
    one number field per atom, so it is refused by name.  `drop=['map_number']` writes it without the
    mapping, which is a loss the caller asked for."""
    rxn = _reaction('[CH3:1][C:2](=[O:3])O.CN>>CC')
    with pytest.raises(ValueError) as err:
        rxn.pack(version=1)
    assert 'map_number' in str(err.value)
    back = ReactionContainer.unpack(rxn.pack(drop=['map_number'], version=1))
    assert len(back) == 3


def test_a_map_number_above_the_formats_12_bit_field_is_refused():
    rxn = _reaction('[CH3:1][OH:4096]>>[CH3:1][OH:4096]')
    with pytest.raises(ValueError) as err:
        rxn.pack(version=1)
    assert '4096' in str(err.value) or '4095' in str(err.value)


def test_two_atoms_sharing_a_map_number_inside_one_molecule_are_refused():
    mol = read_smiles('[CH3:1][OH:1]')
    rxn = ReactionContainer([mol], [mol.copy()])
    with pytest.raises(ValueError) as err:
        rxn.pack(version=1)
    assert 'map_number' in str(err.value)


# --------------------------------------------------------------------------------------------------
# the molecule door says what a reaction buffer is instead of guessing
# --------------------------------------------------------------------------------------------------

@pytest.mark.parametrize('compressed', [True, False])
def test_molecule_unpack_names_a_reaction_record_instead_of_failing_obscurely(compressed):
    raw = _reaction(AMIDATION).pack(compressed=compressed)
    with pytest.raises(ValueError) as err:
        MoleculeContainer.unpack(raw)
    assert 'reaction' in str(err.value).lower()


# --------------------------------------------------------------------------------------------------
# THE FIXTURE.  Bytes an installed chython 2.24 wrote, and the answers chython 2.24's own unpacker
# gave for them.  This is the only test in the file that can fail when the writer and the reader
# agree with each other and with nothing else.
# --------------------------------------------------------------------------------------------------

def _corpus():
    return load_corpus(V2_PATH)


def test_the_fixture_is_present_and_was_written_by_chython_two():
    records = _corpus()
    assert len(records) >= 8
    for record in records:
        assert record['data'][0] == 1, record['name']


def test_every_chython_two_reaction_record_decodes_to_the_answers_chython_two_gave():
    records = _corpus()
    for record in records:
        rxn, problems = reaction_pach_load(record['data'], compressed=False)
        assert rxn is not None, (record['name'], problems)
        assert problems == [], (record['name'], problems)
        answers = record['answers']
        assert [len(rxn.reactants), len(rxn.agents), len(rxn.products)] == answers['counts'], \
            record['name']
        got = []
        for mol in rxn.molecules():
            maps = {a.n: a.map_number for a in mol.atoms()}
            atoms = sorted([a.map_number, a.element, a.isotope or None, a.charge,
                            int(a.is_radical), a.implicit_h, a.degree] for a in mol.atoms())
            bonds = sorted([min(maps[b.n], maps[b.m]), max(maps[b.n], maps[b.m]), b.order]
                           for b in mol.bonds())
            got.append({'atoms': atoms, 'bonds': bonds})
        assert got == answers['molecules'], record['name']


def test_the_v3_writer_reproduces_chython_twos_reaction_layout():
    """The reaction LAYER byte for byte: the four header bytes and the record boundaries.

    Not the whole buffer, and the reason is stated in `test_pach.py`: chython 2 wrote a neighbour
    list in insertion order where the arena's is index-ascending, so the molecule records are not
    byte-identical for every molecule and cannot be made so without reproducing that. What the
    reaction layer owns -- the header and where each record starts -- is identical for every record,
    and the molecule-level identity rate is measured and pinned below.
    """
    records = _corpus()
    for record in records:
        rxn, problems = reaction_pach_load(record['data'], compressed=False)
        assert not problems, record['name']
        mine = rxn.pack(compressed=False, version=1)
        assert mine[:4] == record['data'][:4], record['name']
        assert _boundaries(mine) == _boundaries(record['data']), record['name']


def _boundaries(raw):
    """The offsets each molecule record starts at, and the atom count each one declares."""
    from chython.core import pach_record_length
    out = []
    shift = 4
    for _ in range(raw[1] + raw[2] + raw[3]):
        out.append((shift, (raw[shift + 1] << 4) | (raw[shift + 2] >> 4)))
        shift += pach_record_length(raw[shift:], compressed=False)
    return out


def test_the_measured_byte_identity_rate_against_chython_two(capsys):
    """MEASURED, not asserted at 100%: see the docstring above.  Pinned so a regression is visible.

    8 of the 10 records re-encode byte for byte.  The two that do not are `suzuki_mapped` and
    `toluene_nitration`, and both are the MOLECULE-layer divergence `test_pach.py` already names:
    chython 2 wrote an atom's neighbour list in `_bonds` insertion order, the arena's is
    index-ascending, and a substituted ring is where the two first disagree.  The record LENGTHS and
    the record boundaries are identical -- only the order of entries inside the connection table
    differs, which is why the layout test above passes for all ten.
    """
    records = _corpus()
    identical = []
    for record in records:
        rxn, _ = reaction_pach_load(record['data'], compressed=False)
        if rxn.pack(compressed=False, version=1) == record['data']:
            identical.append(record['name'])
        else:
            assert len(rxn.pack(compressed=False, version=1)) == len(record['data']), record['name']
    with capsys.disabled():
        print('\n  reaction pach re-encode byte-identical: %d / %d' % (len(identical), len(records)))
    assert len(identical) >= 8, 'the writer stopped reproducing chython 2 bytes: %s' % identical


def test_pack_len_agrees_with_chython_twos_own_pack_len():
    """The atom counts, from the same bytes, without decoding -- against chython 2's answer for them.

    chython 2's `pack_len` walked the stream with arithmetic of its own rather than a length function,
    so this is a second, independent statement that the record boundaries are where V3 puts them.
    """
    for record in _corpus():
        expected = tuple(tuple(side) for side in record['answers']['atom_counts'])
        assert ReactionContainer.pack_len(record['data'], compressed=False) == expected, record['name']


def test_the_unmapped_fixture_record_comes_back_with_the_numbers_the_record_held():
    """THE AMBIGUITY, PINNED.  `esterification_unmapped` was parsed by chython 2 from a string with no
    mapping in it; chython 2 numbered its atoms 1..10 across the whole record and wrote those numbers
    into the format's one number field.  V3 reports them as map numbers, because the format cannot say
    which of the two they were and losing a real mapping is the worse error of the two.  If this
    assertion is ever deliberately changed, the module docstring of `reaction.py` changes with it."""
    record = next(r for r in _corpus() if r['name'] == 'esterification_unmapped')
    rxn, problems = reaction_pach_load(record['data'], compressed=False)
    assert not problems
    numbers = [sorted(a.map_number for a in m.atoms()) for m in rxn.molecules()]
    assert numbers == [[1, 2, 3, 4], [5, 6], [7, 8, 9, 10, 11], [12]], numbers
    assert all(a.map_number == a.n for m in rxn.molecules() for a in m.atoms())


def test_the_fixture_carries_a_mapped_reaction_and_the_mapping_comes_back():
    """The fixture is not merely decodable: at least one record is a mapped reaction whose map
    numbers relate the two sides, and those numbers are what the reader reports."""
    found = 0
    for record in _corpus():
        rxn, _ = reaction_pach_load(record['data'], compressed=False)
        left = {a.map_number for m in rxn.reactants for a in m.atoms()}
        right = {a.map_number for m in rxn.products for a in m.atoms()}
        if left and right and len(left & right) > 1:
            found += 1
    assert found >= 4, 'the fixture has no mapped reaction in it and cannot pin the mapping'


def test_the_gzipped_container_is_readable_without_the_helper():
    """A sanity check on the artefact itself, so a corrupted commit fails here and not in ten
    assertions that blame the codec."""
    with gzip.open(V2_PATH, 'rb') as f:
        blob = f.read()
    count = int.from_bytes(blob[:4], 'little')
    assert count == len(_corpus())
    assert blob[:4] == struct_pack('<I', count)


# --------------------------------------------------------------------------------------------------
# version 5: molecules carry their own map numbers
# --------------------------------------------------------------------------------------------------

def test_a_reaction_packs_as_version_5_by_default():
    rxn = _reaction(AMIDATION)
    raw = rxn.pack(compressed=False)
    assert raw[0] == 5
    assert (raw[1], raw[2], raw[3]) == (2, 1, 2)
    shift = 4
    for mol in rxn.molecules():
        assert raw[shift] in (3, 4), 'each body record is a third-generation molecule'
        assert (raw[shift + 2] | (raw[shift + 3] << 8)) == len(mol)
        shift += _record_length(raw, shift)
    assert shift == len(raw), 'the records exactly fill the buffer'


def test_version_5_carries_the_mapping_on_the_molecules_themselves():
    rxn = _reaction(AMIDATION)
    back = ReactionContainer.unpack(rxn.pack())
    assert _snapshot(back) == _snapshot(rxn)


def test_version_5_writes_a_partially_mapped_molecule_that_version_1_refuses():
    rxn = ReactionContainer(reactants=(read_smiles('[CH3:1][OH:2]'), read_smiles('[CH3:3]O')),
                            products=(read_smiles('[CH3:1][NH2:4]'),))
    with pytest.raises(ValueError, match='PARTIALLY mapped'):
        rxn.pack(version=1)
    back = ReactionContainer.unpack(rxn.pack())
    assert [[back_mol.map_number_of(k) for k in back_mol.atom_numbers]
            for back_mol in back.molecules()] == [[1, 2], [3, 0], [1, 4]]


def test_an_unmapped_version_5_reaction_comes_back_unmapped():
    """The one ambiguity version 1 has: a record it wrote cannot tell "mapped 1..N" from "unmapped"."""
    rxn = ReactionContainer(reactants=(read_smiles('CCO'),), products=(read_smiles('CC=O'),))
    back = ReactionContainer.unpack(rxn.pack())
    for mol in back.molecules():
        assert all(mol.map_number_of(k) == 0 for k in mol.atom_numbers)


def test_version_5_carries_a_stereo_group_that_version_1_cannot():
    from chython.core import STEREO_AND
    alanine = read_smiles('C[C@H](N)C(=O)O')
    alanine.set_stereo_group(alanine.atom_numbers[1], STEREO_AND, 1)
    rxn = ReactionContainer(reactants=(read_smiles('CCO'),), products=(alanine,))
    with pytest.raises(ValueError, match='stereo_groups'):
        rxn.pack(version=1)                            # the version 2 record has no field for one
    back = ReactionContainer.unpack(rxn.pack())
    product = back.products[0]
    assert product.stereo_group_of(product.atom_numbers[1]) == (STEREO_AND, 1)


def test_version_1_is_still_written_on_request():
    rxn = _reaction(AMIDATION)
    raw = rxn.pack(compressed=False, version=1)
    assert raw[0] == 1 and raw[4] == 2
    assert _snapshot(ReactionContainer.unpack(raw, compressed=False)) == _snapshot(rxn)


def test_an_unwritable_reaction_version_is_refused_by_name():
    with pytest.raises(ValueError, match='not a writable reaction pach version'):
        _reaction(AMIDATION).pack(version=2)


def test_pack_len_reads_a_version_5_record():
    rxn = _reaction(AMIDATION)
    assert ReactionContainer.pack_len(rxn.pack()) == ((4, 2), (3,), (5, 1))


def test_an_unknown_reaction_version_byte_is_reported():
    raw = bytearray(_reaction(AMIDATION).pack(compressed=False))
    raw[0] = 7
    rxn, problems = reaction_pach_load(bytes(raw), compressed=False)
    assert rxn is None
    assert any('reaction pach version' in p for p in problems)


def test_a_version_5_record_with_a_missing_molecule_returns_no_reaction():
    raw = _reaction(AMIDATION).pack(compressed=False)
    rxn, problems = reaction_pach_load(raw[:-4], compressed=False)
    assert rxn is None and problems


def test_version_5_carries_coordinates_and_wedges_through_a_reaction():
    """Version 5's body records are the molecule codec's own, so a molecule that has a depiction keeps it.
    `_snapshot` compares neither field, and version 1 stored no coordinates at all, so this is the only
    place either is asserted at the reaction layer -- and the only test reaching `_pach_molecule_atom_count`'s
    version 3 branch, since a record built from SMILES is coordinate-free version 4."""
    mol1 = read_smiles('CCO')
    n1 = mol1.atom_numbers
    with mol1.edit() as e:
        e.set_xy(n1[0], 0.0, 0.0)
        e.set_xy(n1[1], 1.5, 0.0)
        e.set_xy(n1[2], 0.0, -2.25)
        e.set_wedge(n1[0], n1[1], 1)
    mol2 = read_smiles('CCO')                               # no coordinates -> version 4 body record
    mol3 = read_smiles('CCN')
    n3 = mol3.atom_numbers
    with mol3.edit() as e:
        e.set_xy(n3[0], 0.0, 0.0)
        e.set_xy(n3[1], 1.5, 0.0)
        e.set_xy(n3[2], 0.0, -2.25)
    rxn = ReactionContainer(reactants=(mol1, mol2), products=(mol3,))
    assert ReactionContainer.pack_len(rxn.pack()) == ((3, 3), (), (3,))
    raw = rxn.pack(compressed=False)
    shift = 4
    body_versions = []
    for _ in rxn.molecules():
        body_versions.append(raw[shift])
        shift += _record_length(raw, shift)
    assert body_versions == [3, 4, 3], 'mol1 and mol3 have coordinates (v3); mol2 does not (v4)'
    back = ReactionContainer.unpack(rxn.pack())
    r1, r2, p = back.reactants[0], back.reactants[1], back.products[0]
    assert r1.xy_of(r1.atom_numbers[0]) == (0.0, 0.0)
    assert r1.xy_of(r1.atom_numbers[1]) == (1.5, 0.0)
    assert r1.xy_of(r1.atom_numbers[2]) == (0.0, -2.25)
    assert not r2.has_coordinates
    assert p.xy_of(p.atom_numbers[1]) == (1.5, 0.0)
    assert r1.wedges() == mol1.wedges()


def test_version_5_stores_a_mapping_version_1_has_no_field_for():
    """One number field per atom cannot hold a number above 4095 or the same number twice, and version 1
    refuses both by name.  Version 5's map block holds a uint16 per atom with no uniqueness rule."""
    rxn_over = ReactionContainer(reactants=(read_smiles('[CH3:1][OH:4096]'),),
                                 products=(read_smiles('[CH3:1][OH:4096]'),))
    with pytest.raises(ValueError, match='4096'):
        rxn_over.pack(version=1)
    back = ReactionContainer.unpack(rxn_over.pack())
    assert [[back_mol.map_number_of(k) for k in back_mol.atom_numbers]
            for back_mol in back.molecules()] == [[1, 4096], [1, 4096]]

    rxn_dup = ReactionContainer(reactants=(read_smiles('[CH3:1][OH:1]'),),
                                products=(read_smiles('[CH3:1][OH:1]'),))
    with pytest.raises(ValueError, match='map_number'):
        rxn_dup.pack(version=1)
    back2 = ReactionContainer.unpack(rxn_dup.pack())
    assert [[back_mol.map_number_of(k) for k in back_mol.atom_numbers]
            for back_mol in back2.molecules()] == [[1, 1], [1, 1]]

    rxn = _reaction(AMIDATION)
    back3 = ReactionContainer.unpack(rxn.pack(drop=['map_number']))
    assert [[back_mol.map_number_of(k) for k in back_mol.atom_numbers]
            for back_mol in back3.molecules()] == [[0, 0, 0, 0], [0, 0], [0, 0, 0],
                                                   [0, 0, 0, 0, 0], [0]]


def test_a_version_in_the_writable_set_with_no_writer_is_refused(monkeypatch):
    """The frozenset and the writer's dispatch move together.  A version added to the set alone would
    otherwise be written as version 5's bytes under its own header byte -- a record with no reader."""
    monkeypatch.setattr(reaction, '_PACH_REACTION_VERSIONS', frozenset({1, 5, 6}))
    with pytest.raises(ValueError, match='in the writable version set but has no writer'):
        _reaction(AMIDATION).pack(version=6)


def test_a_boolean_version_is_refused_by_name():
    """True.__class__ is bool, not int, so it must not pass as version 1."""
    with pytest.raises(ValueError, match='not a writable reaction pach version'):
        _reaction(AMIDATION).pack(version=True)
