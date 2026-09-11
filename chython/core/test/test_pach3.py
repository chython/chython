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
"""Versions 3 and 4 of the pach record.

ASSERTED, not measured: every byte layout here is the one `docs/superpowers/specs/
2026-09-07-pach-v3-design.md` states, so a test that disagrees with the spec is a defect in the test.
The v0/v2 corpora in `test_pach.py` are the other kind -- there the format is the subject and chython
2's own writer produced the bytes.
"""
import tracemalloc

from pytest import mark, raises

from chython.core import H_UNKNOWN, MoleculeContainer, pach_load, pach_record_length, read_smiles
from .pach3_corpus import BUILDERS, V3_PATH, V4_PATH, answers, drawn, load_corpus
from .pach3_corpus import _BIPHENYL_BONDS, _BIPHENYL_H, _BUTANE_XY, _chlorofluorobiphenyl, _drawn_butane
from .pach_corpus import V0_NATIVE_PATH, V0_PATH, V2_PATH


def _v3_header(version=4, flags=0, atoms=0, bonds=0, stereo=0, sgroups=0):
    """A 12 byte v3/v4 header as bytes, every count stated."""
    out = bytearray(12)
    out[0] = version
    out[1] = flags
    out[2:4] = atoms.to_bytes(2, 'little')
    out[4:6] = bonds.to_bytes(2, 'little')
    out[6:8] = stereo.to_bytes(2, 'little')
    out[8:10] = sgroups.to_bytes(2, 'little')
    return bytes(out)


def test_length_of_a_version_4_header_is_arithmetic():
    # three atoms at 3 bytes, two bonds at 5, one stereo record at 9, one sgroup entry at 3
    raw = _v3_header(4, atoms=3, bonds=2, stereo=1, sgroups=1) + bytes(9 + 10 + 9 + 3)
    assert pach_record_length(raw, compressed=False) == 12 + 9 + 10 + 9 + 3


def test_length_of_a_version_3_header_counts_nine_bytes_per_atom():
    raw = _v3_header(3, atoms=3, bonds=2) + bytes(27 + 10)
    assert pach_record_length(raw, compressed=False) == 12 + 27 + 10


def test_length_counts_the_map_block_when_the_flag_is_set():
    raw = _v3_header(4, flags=1, atoms=3, bonds=2) + bytes(9 + 10 + 6)
    assert pach_record_length(raw, compressed=False) == 12 + 9 + 10 + 6


def test_a_truncated_version_4_header_is_refused_by_name():
    with raises(ValueError, match='12 byte header'):
        pach_record_length(b'\x04\x00\x03\x00', compressed=False)


def test_an_unknown_version_byte_is_refused_by_name():
    with raises(ValueError, match='not a pach version'):
        pach_record_length(b'\x07' + bytes(11), compressed=False)


def test_the_legacy_doors_still_answer_version_2():
    raw = read_smiles('CCO').pack(compressed=False, version=2)
    assert raw[0] == 2
    assert pach_record_length(raw, compressed=False) == len(raw)


def test_a_lone_sodium_cation_is_a_header_and_three_bytes():
    raw = read_smiles('[Na+]').pack(compressed=False, version=4)
    assert raw == _v3_header(4, atoms=1) + bytes([11, 0x80, 0x50])
    # element 11, no isotope and no radical, h_pinned set (bracket atom), 0 implicit H,
    # charge +1 as 5 in the high nibble


def test_the_isotope_field_is_a_shift_of_thirty_two():
    raw = read_smiles('[13CH4]').pack(compressed=False, version=4)
    # carbon's MDL reference is 12, so 13 spells 13 - 12 + 32 = 33; bracket atom so h_pinned is set:
    # byte 1 = 0x80 | 33 = 0xa1
    assert raw[12:15] == bytes([6, 0xa1, 0x44])


def test_h_unknown_is_the_arena_nibble_verbatim():
    mol = read_smiles('[SeH4]')
    mol.set_hydrogens(mol.atom_numbers[0], H_UNKNOWN)
    raw = mol.pack(compressed=False, version=4)
    assert raw[14] & 0x0f == 15


def test_a_pinned_hydrogen_count_sets_bit_seven():
    raw = read_smiles('[CH3-]').pack(compressed=False, version=4)
    assert raw[13] & 0x80


def test_a_radical_sets_bit_six():
    raw = read_smiles('[CH3] |^1:0|').pack(compressed=False, version=4)
    assert raw[13] & 0x40


def test_an_r_marker_sets_bit_seven_of_the_element_byte():
    raw = read_smiles('[R7]').pack(compressed=False, version=4)
    assert raw[12] == 0x87
    raw = read_smiles('[R]').pack(compressed=False, version=4)
    assert raw[12] == 0x80


def test_a_title_is_refused_by_name():
    mol = read_smiles('[Na+]')
    mol.set_title('sodium')
    with raises(ValueError, match="drop=\\['title'\\]"):
        mol.pack(version=4)
    assert mol.pack(compressed=False, version=4, drop=['title'])[12] == 11


def test_version_3_on_a_coordinate_free_molecule_writes_version_4():
    raw = read_smiles('[Na+]').pack(compressed=False, version=3)
    assert raw == _v3_header(4, atoms=1) + bytes([11, 0x80, 0x50])


def test_the_declared_length_is_the_buffer_length():
    """`pach_record_length` reads the header and nothing else, so a header that under-declares its own
    record walks a caller stepping through a store into the middle of the next one."""
    for text in ('[Na+]', 'CCO', 'N[C@@H](C)C(=O)O'):
        raw = read_smiles(text).pack(compressed=False)
        assert pach_record_length(raw, compressed=False) == len(raw), text


def test_the_bond_block_is_five_bytes_per_bond():
    raw = read_smiles('CCO').pack(compressed=False, version=4)
    assert raw[:12] == _v3_header(4, atoms=3, bonds=2)
    bonds = raw[12 + 9:]
    assert bonds == bytes([0, 0, 1, 0, 0x01,        # atom 0 -- atom 1, single, no wedge
                           1, 0, 2, 0, 0x01])


def test_an_aromatic_order_is_the_whole_of_aromaticity():
    raw = read_smiles('c1ccccc1').pack(compressed=False, version=4)
    orders = {raw[12 + 6 * 3 + 5 * k + 4] for k in range(6)}
    assert orders == {0x04}


def test_a_dative_bond_is_order_eight():
    raw = read_smiles('[NH3]->[BH3]').pack(compressed=False, version=4)
    assert raw[12 + 2 * 3 + 4] == 0x08


def _drawn_amino_propanol():
    """CC(N)O with a drawing, so a wedge has coordinates to mean something against."""
    mol = read_smiles('CC(N)O')
    n = mol.atom_numbers
    with mol.edit() as e:
        e.set_xy(n[0], 0.0, 0.0)
        e.set_xy(n[1], 1.0, 0.0)
        e.set_xy(n[2], 1.5, 1.0)
        e.set_xy(n[3], 2.0, 0.0)
    return mol, n


def test_a_wedge_is_written_at_its_own_end():
    mol, n = _drawn_amino_propanol()
    with mol.edit() as e:
        e.set_wedge(n[1], n[2], 1)
    raw = mol.pack(compressed=False, version=3)
    assert raw[0] == 3
    block = raw[12 + 4 * 9:]
    found = [block[5 * k:5 * k + 5] for k in range(3)]
    assert bytes([1, 0, 2, 0, 0x11]) in found


def test_a_wedge_whose_narrow_end_is_the_higher_slot_is_written_first():
    mol, n = _drawn_amino_propanol()
    with mol.edit() as e:
        e.set_wedge(n[2], n[1], 1)
    raw = mol.pack(compressed=False, version=3)
    block = raw[12 + 4 * 9:]
    found = [block[5 * k:5 * k + 5] for k in range(3)]
    assert bytes([2, 0, 1, 0, 0x11]) in found
    assert bytes([1, 0, 2, 0, 0x11]) not in found


def test_a_version_4_record_writes_no_wedge_and_says_so():
    mol = read_smiles('CC(N)O')
    n = mol.atom_numbers
    with mol.edit() as e:
        e.set_wedge(n[1], n[2], 1)
    raw = mol.pack(compressed=False, version=4)
    block = raw[12 + 4 * 3:]
    assert all(block[5 * k + 4] >> 4 == 0 for k in range(3))
    assert any(r.rule == 'pach:wedge-lost' for r in mol.log)


def test_a_drawn_molecule_packs_as_version_3():
    raw = _drawn_butane().pack(compressed=False)
    assert raw[0] == 3
    assert pach_record_length(raw, compressed=False) == len(raw) == 12 + 4 * 9 + 3 * 5


def test_a_coordinate_is_an_exact_int24_at_ten_thousand():
    raw = _drawn_butane().pack(compressed=False)
    assert raw[12 + 9 + 3:12 + 9 + 6] == (15000).to_bytes(3, 'little')       # atom 1 x = 1.5
    assert raw[12 + 9 + 6:12 + 9 + 9] == (0).to_bytes(3, 'little')           # atom 1 y = 0.0, written
    assert raw[12 + 2 * 9 + 6:12 + 2 * 9 + 9] == (-12500 & 0xffffff).to_bytes(3, 'little')


def test_a_string_read_molecule_packs_as_version_4():
    assert read_smiles('CCO').pack(compressed=False)[0] == 4


def test_dropping_coordinates_selects_version_4():
    raw = _drawn_butane().pack(compressed=False, drop=['coordinates'])
    assert raw[0] == 4
    assert len(raw) == 12 + 4 * 3 + 3 * 5


def test_a_coordinate_beyond_the_int24_range_is_refused_by_name():
    mol = read_smiles('CCO')
    n = mol.atom_numbers
    with mol.edit() as e:
        e.set_xy(n[0], 0.0, 0.0)
        e.set_xy(n[1], 900.0, 0.0)
        e.set_xy(n[2], 0.0, 0.0)
    with raises(ValueError, match='838.8607'):
        mol.pack(version=3)


# ----- decoder tests: the atom block read back -----

def test_a_one_atom_record_round_trips():
    mol = read_smiles('[13CH3-]')
    back, problems = pach_load(mol.pack(compressed=False), compressed=False)
    assert problems == []
    assert back == mol


def test_an_r_marker_round_trips_with_its_index():
    mol = read_smiles('[R7]')
    back, problems = pach_load(mol.pack(compressed=False), compressed=False)
    assert problems == []
    assert back.atom(back.atom_numbers[0]).r_index == 7


def test_h_unknown_round_trips_as_h_unknown():
    mol = read_smiles('[SeH4]')
    mol.set_hydrogens(mol.atom_numbers[0], H_UNKNOWN)
    back, problems = pach_load(mol.pack(compressed=False), compressed=False)
    assert problems == []
    assert back.implicit_h_of(back.atom_numbers[0]) is None


def test_the_atom_block_being_short_returns_no_molecule_and_says_so():
    raw = _v3_header(4, atoms=3) + bytes(6)
    mol, problems = pach_load(raw, compressed=False)
    assert mol is None
    assert any('atom block' in p for p in problems)


def test_a_reserved_flag_bit_is_reported_and_ignored():
    raw = _v3_header(4, flags=0x02, atoms=1) + bytes([6, 0, 0x44])
    mol, problems = pach_load(raw, compressed=False)
    assert mol is not None
    assert any('flags' in p for p in problems)


def test_an_element_byte_outside_the_domain_reads_as_a_bare_r():
    raw = _v3_header(4, atoms=1) + bytes([0, 0, 0x40])
    mol, problems = pach_load(raw, compressed=False)
    assert mol.element_of(mol.atom_numbers[0]) == 0    # element 0 is the R marker
    assert any('element byte' in p for p in problems)


def test_a_charge_beyond_the_arena_is_clamped_and_reported():
    raw = _v3_header(4, atoms=1) + bytes([6, 0, 0xf0])       # charge nibble 15 == +11
    mol, problems = pach_load(raw, compressed=False)
    assert mol.charge_of(mol.atom_numbers[0]) == 8
    assert any('charge' in p for p in problems)


def test_unpack_answers_a_version_4_record():
    mol = read_smiles('[13CH3-]')
    assert MoleculeContainer.unpack(mol.pack()) == mol


def test_a_buffer_shorter_than_the_header_returns_no_molecule_and_says_so():
    raw = _v3_header(4, atoms=1)[:8]        # 8 bytes, not the required 12
    mol, problems = pach_load(raw, compressed=False)
    assert mol is None
    assert any('at least a 12 byte header' in p for p in problems)


def test_an_r_index_above_the_maximum_is_reported_and_decoded_as_a_bare_r():
    # R_INDEX_MAX is 99; bit 7 set with index 100 = byte 0x80 | 100 = 0xe4
    raw = _v3_header(4, atoms=1) + bytes([0xe4, 0, 0x40])
    mol, problems = pach_load(raw, compressed=False)
    assert mol.element_of(mol.atom_numbers[0]) == 0      # bare R, element 0
    assert mol.atom(mol.atom_numbers[0]).r_index == 0    # index cleared
    assert any('states R index' in p for p in problems)


# ----- decoder tests: the bond block and the wedges -----

def test_ethanol_round_trips_through_version_4():
    mol = read_smiles('CCO')
    back, problems = pach_load(mol.pack(compressed=False), compressed=False)
    assert problems == []
    assert back == mol


def test_a_drawn_molecule_round_trips_byte_for_byte():
    first = _drawn_butane().pack(compressed=False)
    back, problems = pach_load(first, compressed=False)
    assert problems == []
    assert back.pack(compressed=False) == first


def test_an_aromatic_ring_round_trips_unkekulised():
    mol = read_smiles('c1ccccc1')
    back, problems = pach_load(mol.pack(compressed=False), compressed=False)
    assert problems == []
    assert back.aromatic_bond_count == 6
    assert back == mol


def test_a_wedge_round_trips_on_the_narrow_end():
    mol = _drawn_butane()
    n = mol.atom_numbers
    with mol.edit() as e:
        e.set_wedge(n[1], n[2], 2)
    back, problems = pach_load(mol.pack(compressed=False), compressed=False)
    assert problems == []
    assert back.wedge_of(2, 3) == 2
    assert back.wedge_of(3, 2) == 0


def test_a_bond_naming_an_atom_that_is_not_there_is_dropped_and_reported():
    raw = _v3_header(4, atoms=2, bonds=2) + bytes([6, 0, 0x43, 6, 0, 0x43]) \
        + bytes([0, 0, 1, 0, 0x01]) + bytes([0, 0, 9, 0, 0x01])
    mol, problems = pach_load(raw, compressed=False)
    assert mol.bond_count == 1
    assert any('atom index' in p for p in problems)


def test_a_bond_order_of_zero_reads_as_single_and_is_reported():
    raw = _v3_header(4, atoms=2, bonds=1) + bytes([6, 0, 0x43, 6, 0, 0x43]) \
        + bytes([0, 0, 1, 0, 0x00])
    mol, problems = pach_load(raw, compressed=False)
    assert mol.order_of(1, 2) == 1
    assert any('order' in p for p in problems)


def test_a_wedge_in_a_version_4_record_is_read_as_none_and_reported():
    raw = _v3_header(4, atoms=2, bonds=1) + bytes([6, 0, 0x43, 6, 0, 0x43]) \
        + bytes([0, 0, 1, 0, 0x11])
    mol, problems = pach_load(raw, compressed=False)
    assert mol.wedge_of(1, 2) == 0
    assert any('wedge' in p for p in problems)


def test_a_truncated_bond_block_keeps_the_bonds_it_has():
    raw = read_smiles('CCO').pack(compressed=False)
    mol, problems = pach_load(raw[:-5], compressed=False)
    assert mol.bond_count == 1
    assert problems == ['the header declares 2 bond(s) and the buffer holds 1; the rest of the block '
                        'was not read']


def test_a_self_loop_bond_is_dropped_and_reported():
    raw = _v3_header(4, atoms=2, bonds=1) + bytes([6, 0, 0x43, 6, 0, 0x43]) \
        + bytes([0, 0, 0, 0, 0x01])
    mol, problems = pach_load(raw, compressed=False)
    assert mol.bond_count == 0
    assert any('itself' in p for p in problems)


def test_a_duplicate_edge_is_dropped_and_reported():
    raw = _v3_header(4, atoms=2, bonds=2) + bytes([6, 0, 0x43, 6, 0, 0x43]) \
        + bytes([0, 0, 1, 0, 0x01]) + bytes([0, 0, 1, 0, 0x01])
    mol, problems = pach_load(raw, compressed=False)
    assert mol.bond_count == 1
    assert any('repeat' in p for p in problems)


# ----- encoder tests: the stereo block -----

def _record(raw, atoms, bonds, index=0):
    """The `index`-th stereo record of `raw`, nine bytes."""
    stride = 9 if raw[0] == 3 else 3
    at = 12 + atoms * stride + bonds * 5 + index * 9
    return raw[at:at + 9]


def _slots(rec):
    """A record's four slots and its kind, as the tuple a test states literally.  Slots are atom
    indices, so they are the positions in `mol.atom_numbers` rather than the stable ids `refs` uses."""
    return (int.from_bytes(rec[0:2], 'little'), int.from_bytes(rec[2:4], 'little'),
            int.from_bytes(rec[4:6], 'little'), int.from_bytes(rec[6:8], 'little'), rec[8] & 0x07)


def _sign(rec):
    """A record's parity bit: 0 even, 1 odd.  Byte 8's high nibble is reserved and checked here, so
    every test that reads a sign pins it as zero."""
    assert rec[8] & 0xf0 == 0, 'byte 8 bits 4-7 are reserved and must be written zero'
    return (rec[8] >> 3) & 1


def test_a_tetrahedral_centre_is_one_record_of_centre_and_three_directions():
    mol = read_smiles('N[C@@H](C)C(=O)O')
    raw = mol.pack(compressed=False)
    assert raw[6:8] == (1).to_bytes(2, 'little')
    rec = _record(raw, 6, 5)
    # the alpha carbon, then N, C and C in refs order; the implicit H is the implied fourth direction
    assert _slots(rec) == (1, 0, 2, 3, 0)
    # the same three directions as stable ids -- the frame the sign is the parity in
    assert _sign(rec) == mol.translate_stereo(2, (1, 3, 4, None)) - 1


def test_a_double_bond_is_one_record_of_two_ends_and_one_direction_each():
    mol = read_smiles('C/C=C/C')
    raw = mol.pack(compressed=False)
    assert raw[6:8] == (1).to_bytes(2, 'little')
    rec = _record(raw, 4, 3)
    # end 1 with its own direction 0, then end 2 with its own direction 3
    assert _slots(rec) == (1, 0, 2, 3, 1)
    assert _sign(rec) == mol.translate_stereo(2, (1, None, 4, None)) - 1


def test_an_allene_stores_its_ends_and_not_its_centre():
    # 1,3-dibromo-1,3-difluoroallene; `[C@]` on the centre is OpenSMILES extended tetrahedral, which
    # the arena calls SU_ALLENE.  Atom order: F0 C1 Br2 C3(centre) C4 F5 Br6.
    mol = read_smiles('FC(Br)=[C@]=C(F)Br')
    rec = _record(mol.pack(compressed=False), 7, 6)
    # the chain ends 1 and 4 and never the centre 3, each followed by a direction OF ITS OWN: F0 is
    # bonded to C1 and F5 to C4, which is the pairing the encoder picks by looking for the bond.
    assert _slots(rec) == (1, 0, 4, 5, 2)
    assert _sign(rec) == mol.translate_stereo(4, (1, 3, 6, 7)) - 1


def test_an_atropisomer_stores_both_pivots():
    mol = _chlorofluorobiphenyl()
    rec = _record(mol.pack(compressed=False), 14, 15)
    # the two pivots of the biaryl axis, each followed by one of its own ring neighbours
    assert _slots(rec) == (0, 1, 6, 7, 3)
    assert _sign(rec) == mol.translate_stereo(1, (2, 6, 8, 12)) - 1


def test_the_two_parities_of_one_centre_differ_only_in_the_sign_bit():
    left = read_smiles('N[C@@H](C)C(=O)O').pack(compressed=False)
    right = read_smiles('N[C@H](C)C(=O)O').pack(compressed=False)
    assert left[:-1] == right[:-1]
    assert left[-1] ^ right[-1] == 0x08
    assert (_sign(_record(left, 6, 5)), _sign(_record(right, 6, 5))) == (0, 1)


def test_a_version_3_record_puts_the_stereo_block_after_nine_bytes_an_atom():
    mol = read_smiles('N[C@@H](C)C(=O)O')
    with mol.edit() as e:
        for i, n in enumerate(mol.atom_numbers):
            e.set_xy(n, float(i), 0.0)
    raw = mol.pack(compressed=False, version=3)
    assert raw[0] == 3
    # `_record` finds the block at 12 + 6 * 9 + 5 * 5 = 91, so the slots it reads prove the offset
    assert _slots(_record(raw, 6, 5)) == (1, 0, 2, 3, 0)
    assert pach_record_length(raw, compressed=False) == len(raw) == 100


def test_a_tetrahedral_centre_with_two_named_directions_is_refused_by_name():
    # A sulfonium: its lone pair is one direction and its hydrogen another, so only two directions have
    # an atom of their own and a record with three slots for them has nothing to put in the third.
    mol = read_smiles('C[S@@H+]CC')
    with raises(ValueError, match='names 2 direction'):
        mol.pack()
    raw = mol.pack(compressed=False, drop=['stereo'])
    assert raw[6:8] == b'\x00\x00'
    assert len(raw) == 12 + 4 * 3 + 3 * 5


def test_dropping_stereo_writes_no_stereo_block():
    raw = read_smiles('N[C@@H](C)C(=O)O').pack(compressed=False, drop=['stereo'])
    assert raw[6:8] == b'\x00\x00'
    assert len(raw) == 12 + 6 * 3 + 5 * 5


_STEREO_CORPUS = ['N[C@@H](C)C(=O)O',                            # alanine
                  '[C@@H](N)(C)C(=O)O',                          # the same centre written anchor first
                  'N[C@@]([H])(C)C(=O)O',                        # ... and with its hydrogen named
                  'C/C=C/C', 'C/C=C\\C',                         # both parities of one double bond
                  'F/C=C/C=C/F',                                 # two records in one molecule
                  'C(/F)=C(\\F)Cl',                              # one end named twice, one once
                  'FC(Br)=[C@]=C(F)Br',                          # extended tetrahedral
                  'O[C@H]1CC[C@@H](O)CC1',                       # cyclohexane-1,4-diol
                  'C[C@H](O)[C@@H](N)C(=O)O',                    # threonine
                  'OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O']  # glucose, five centres


def test_the_frame_is_refs_order_on_every_configured_unit():
    """The harness for `_pach3.pxi`'s claim that the record's direction frame is `refs` order.

    Per configured unit: no direction list puts its unnamed direction before its named one, so the frame
    is `refs` as it stands; slot1 is a direction of slot0 and slot3 of slot2; and the sign is the parity
    read in that frame, which `translate_stereo` answers by its own path.
    """
    kinds = set()
    for mol in [read_smiles(s) for s in _STEREO_CORPUS] + [_chlorofluorobiphenyl()]:
        raw = mol.pack(compressed=False, version=4)
        numbers = mol.atom_numbers
        slot_of = {n: i for i, n in enumerate(numbers)}
        units = [u for u in mol.stereo_units() if u['parity']]
        assert raw[6:8] == len(units).to_bytes(2, 'little')
        for index, u in enumerate(units):
            refs, kind, anchor = u['refs'], u['kind'], u['anchor']
            rec = _record(raw, mol.atom_count, mol.bond_count, index)
            slot0, slot1, slot2, slot3, written = _slots(rec)
            assert written == kind
            if kind == 0:
                named = tuple(r for r in refs if r is not None)
                assert refs[:len(named)] == named, 'an unnamed direction ahead of a named one'
                frame = named + (None,) * (4 - len(named))
                assert slot0 == slot_of[anchor]
                assert (slot1, slot2, slot3) == tuple(slot_of[r] for r in frame[:3])
                assert all(mol.order_of(anchor, r) is not None for r in frame[:3])
            else:
                assert refs[0] is not None and refs[2] is not None, 'a list whose unnamed direction leads'
                frame = refs
                assert (slot1, slot3) == (slot_of[frame[0]], slot_of[frame[2]])
                assert mol.order_of(numbers[slot0], frame[0]) is not None
                assert mol.order_of(numbers[slot2], frame[2]) is not None
            assert _sign(rec) == mol.translate_stereo(anchor, frame) - 1
            kinds.add(kind)
    assert kinds == {0, 1, 2, 3}, 'the corpus stopped covering a stereo kind'


# ----- decoder tests: the stereo block read back -----

def test_a_tetrahedral_parity_survives_the_round_trip():
    mol = read_smiles('N[C@@H](C)C(=O)O')
    back, problems = pach_load(mol.pack(compressed=False), compressed=False)
    assert problems == []
    assert back == mol


def test_the_two_enantiomers_stay_two_molecules():
    left = read_smiles('N[C@@H](C)C(=O)O')
    right = read_smiles('N[C@H](C)C(=O)O')
    back_left, _ = pach_load(left.pack(compressed=False), compressed=False)
    back_right, _ = pach_load(right.pack(compressed=False), compressed=False)
    assert back_left == left and back_right == right
    assert back_left != back_right


def test_flipping_the_sign_bit_decodes_the_mirror_image():
    """The frame algebra's sharpest statement: the record's sign bit is the parity IN THE RECORD'S
    frame, so flipping it and nothing else is exactly the enantiomer."""
    raw = bytearray(read_smiles('N[C@@H](C)C(=O)O').pack(compressed=False))
    raw[-1] ^= 0x08
    back, problems = pach_load(bytes(raw), compressed=False)
    assert problems == []
    assert back == read_smiles('N[C@H](C)C(=O)O')


def test_both_cis_trans_configurations_survive():
    for text in ('C/C=C/C', 'C/C=C\\C'):
        mol = read_smiles(text)
        back, problems = pach_load(mol.pack(compressed=False), compressed=False)
        assert problems == []
        assert back == mol, text


def test_both_allene_configurations_survive():
    for text in ('FC(Br)=[C@]=C(F)Br', 'FC(Br)=[C@@]=C(F)Br'):
        mol = read_smiles(text)
        back, problems = pach_load(mol.pack(compressed=False), compressed=False)
        assert problems == []
        assert back == mol, text


def test_an_atropisomer_survives_both_ways():
    for parity in (1, 2):
        mol = _chlorofluorobiphenyl()
        mol.set_parity(mol.atom_numbers[0], parity)
        back, problems = pach_load(mol.pack(compressed=False), compressed=False)
        assert problems == []
        assert back.parity_of(back.atom_numbers[0]) == parity


def test_an_implicit_hydrogen_double_bond_keeps_its_configuration():
    """(E)- and (Z)-2-butenoic acid: each terminal orders one named direction and one implicit
    hydrogen, so two of the four frame slots are `SU_NO_REF` and carry no atom to be matched on."""
    for text in ('C/C=C/C(=O)O', 'C/C=C\\C(=O)O'):
        mol = read_smiles(text)
        back, problems = pach_load(mol.pack(compressed=False), compressed=False)
        assert problems == []
        assert back == mol, text


def test_every_configured_unit_in_the_corpus_survives_the_round_trip():
    """The encoder's own corpus, read back.  All four kinds, and every one of them by equality rather
    than by parity byte, so a configuration landing on the wrong anchor fails here too."""
    for mol in [read_smiles(s) for s in _STEREO_CORPUS] + [_chlorofluorobiphenyl()]:
        back, problems = pach_load(mol.pack(compressed=False, version=4), compressed=False)
        assert problems == []
        assert back == mol


def test_a_permuted_frame_is_read_in_the_frame_the_record_states():
    """The sign is the parity in the record's OWN direction order, so the same three directions named
    in a different order state the other sign for one molecule.

    Deliberate, because the encoder cannot produce such a record: ruling F26 makes the writer's frame
    `refs` order on every configured unit perception emits, so the permutation is the identity there and
    a decoder that computed it and then ignored it would pass every round trip above.  Transposing
    alanine's slot1 and slot2 is one exchange, which `translate_stereo` prices at the other parity, so
    the swapped record reads as the enantiomer with its sign unchanged and as the original with its sign
    flipped.
    """
    left = read_smiles('N[C@@H](C)C(=O)O')
    right = read_smiles('N[C@H](C)C(=O)O')
    assert left.translate_stereo(2, (1, 3, 4, None)) != left.translate_stereo(2, (3, 1, 4, None))
    raw = read_smiles('N[C@@H](C)C(=O)O').pack(compressed=False)
    swapped = bytearray(raw)
    swapped[-7:-5], swapped[-5:-3] = raw[-5:-3], raw[-7:-5]
    back, problems = pach_load(bytes(swapped), compressed=False)
    assert problems == []
    assert back == right
    swapped[-1] ^= 0x08
    back, problems = pach_load(bytes(swapped), compressed=False)
    assert problems == []
    assert back == left


def _far_owner_first(raw):
    """The last stereo record with its two (owner, direction) pairs exchanged: the same configuration
    stated from the other owner."""
    out = bytearray(raw)
    out[-9:-5], out[-5:-1] = raw[-5:-1], raw[-9:-5]
    return bytes(out)


def test_a_record_stating_the_far_owner_first_reads_the_same_configuration():
    """Exchanging the two direction LISTS is an even permutation, so the sign does not move.

    Unreachable from the encoder, which writes the unit's anchor at slot0; it is the branch ruling F45's
    anchor relocation needs, and the one `SU_NO_REF`'s anonymity makes easy to get wrong.  Both of
    `C/C=C/C`'s lists carry an implicit hydrogen, so matching the two unnamed frame slots against the
    unit's by POSITION rather than by list adds a transposition and inverts every such record.
    """
    for text in ('C/C=C/C', 'C/C=C\\C', 'C(/F)=C(\\F)Cl', 'FC(Br)=[C@]=C(F)Br'):
        mol = read_smiles(text)
        back, problems = pach_load(_far_owner_first(mol.pack(compressed=False)), compressed=False)
        assert problems == []
        assert back == mol, text


def test_a_stereo_record_naming_an_absent_atom_is_dropped_and_reported():
    raw = bytearray(read_smiles('N[C@@H](C)C(=O)O').pack(compressed=False))
    raw[-9:-7] = (900).to_bytes(2, 'little')            # slot0, past the atom block
    mol, problems = pach_load(bytes(raw), compressed=False)
    assert mol.parity_of(mol.atom_numbers[1]) == 0
    assert any('atom index' in p for p in problems)


def test_a_stereo_record_resolving_to_no_unit_is_dropped_and_reported():
    """A hand-built tetrahedral record on ethanol's hydroxyl oxygen, which anchors no unit.

    Not on a carbon: the lookup is against the UNMARKED unit table (ruling F70), which is perception
    without the stereogenicity filter, so ethanol's two carbons do each anchor a tetrahedral record and
    only the oxygen anchors none.
    """
    raw = _v3_header(4, atoms=3, bonds=2, stereo=1) \
        + bytes([6, 0, 0x43, 6, 0, 0x42, 8, 0, 0x41]) \
        + bytes([0, 0, 1, 0, 0x01]) + bytes([1, 0, 2, 0, 0x01]) \
        + (2).to_bytes(2, 'little') + (1).to_bytes(2, 'little') \
        + (0).to_bytes(2, 'little') + (1).to_bytes(2, 'little') + bytes([0x00])
    mol, problems = pach_load(raw, compressed=False)
    assert mol is not None and mol.bond_count == 2
    assert any('no stereo unit' in p for p in problems)


def test_a_stereo_record_naming_a_direction_the_centre_lacks_is_dropped_and_reported():
    """The same hand-built record on ethanol's first carbon, which anchors a unit ordering one named
    direction, so two of the three the record states are not the centre's at all."""
    raw = _v3_header(4, atoms=3, bonds=2, stereo=1) \
        + bytes([6, 0, 0x43, 6, 0, 0x42, 8, 0, 0x41]) \
        + bytes([0, 0, 1, 0, 0x01]) + bytes([1, 0, 2, 0, 0x01]) \
        + (0).to_bytes(2, 'little') + (1).to_bytes(2, 'little') \
        + (2).to_bytes(2, 'little') + (1).to_bytes(2, 'little') + bytes([0x00])
    mol, problems = pach_load(raw, compressed=False)
    assert mol is not None and mol.parity_of(mol.atom_numbers[0]) == 0
    assert any('not one of its four' in p for p in problems)


def test_a_stereo_record_naming_a_direction_the_axis_lacks_is_dropped_and_reported():
    raw = bytearray(read_smiles('C/C=C/C').pack(compressed=False))
    raw[-7:-5] = (3).to_bytes(2, 'little')             # slot1, the FAR terminal's own substituent
    mol, problems = pach_load(bytes(raw), compressed=False)
    assert mol.parity_of(mol.atom_numbers[1]) == 0
    assert any('are not the ones' in p for p in problems)


def test_an_unknown_stereo_kind_is_dropped_and_reported():
    raw = bytearray(read_smiles('N[C@@H](C)C(=O)O').pack(compressed=False))
    raw[-1] = 0x06                                      # kind 6, and no kind 6 exists
    mol, problems = pach_load(bytes(raw), compressed=False)
    assert mol.parity_of(mol.atom_numbers[1]) == 0
    assert any('kind 6' in p for p in problems)


def test_a_second_record_for_one_unit_is_dropped_and_reported():
    raw = read_smiles('N[C@@H](C)C(=O)O').pack(compressed=False)
    doubled = bytearray(raw)
    doubled[6:8] = (2).to_bytes(2, 'little')
    doubled += raw[-9:]
    mol, problems = pach_load(bytes(doubled), compressed=False)
    assert mol == read_smiles('N[C@@H](C)C(=O)O')
    assert any('already configured' in p for p in problems)


def test_a_truncated_stereo_block_drops_the_one_record_it_cannot_read():
    raw = read_smiles('N[C@@H](C)C(=O)O').pack(compressed=False)
    mol, problems = pach_load(raw[:-4], compressed=False)
    assert mol is not None
    assert mol.parity_of(mol.atom_numbers[1]) == 0
    assert problems == ['the header declares 1 stereo configuration(s) and the buffer holds 0; the rest '
                        'of the block was not read']


def test_a_truncated_stereo_block_keeps_the_records_it_has():
    """Two records, four bytes short: the first is whole and is applied, the second is not read.

    The clip is `ns_declared = ns_have`, and only a block holding more than one record tells that apart
    from dropping the block -- so the assertion is the SURVIVING configuration, not the report.
    """
    raw = read_smiles('F/C=C/C=C/F').pack(compressed=False)
    mol, problems = pach_load(raw[:-4], compressed=False)
    assert [mol.parity_of(n) for n in mol.atom_numbers] == [0, 1, 0, 0, 0, 0]
    assert mol == read_smiles('F/C=C/C=CF')
    assert problems == ['the header declares 2 stereo configuration(s) and the buffer holds 1; the rest '
                        'of the block was not read']


def test_a_reserved_bit_of_the_stereo_flag_byte_is_reported():
    mol = read_smiles('N[C@@H](C)C(=O)O')
    raw = bytearray(mol.pack(compressed=False))
    assert raw[6:8] == bytearray(b'\x01\x00')            # one stereo record, or the offset is wrong
    raw[12 + len(mol) * 3 + mol.bond_count * 5 + 8] |= 0xf0
    back, problems = pach_load(bytes(raw), compressed=False)
    assert back == mol                                  # ignored, not a reason to drop
    assert any('reserved and must be 0' in p for p in problems)


def test_a_stereo_count_a_truncated_bond_block_never_reached_is_reported():
    """The stereo block starts where the header's bond count puts it, which here is past the buffer's
    end, so the count is reported and NOTHING IS READ AT IT.

    The full list is the assertion: a reader taking its offsets from what it managed to read instead
    fabricates records off the end of the buffer, and the list runs to hundreds of entries.
    """
    raw = read_smiles('N[C@@H](C)C(=O)O').pack(compressed=False)
    mol, problems = pach_load(raw[:-11], compressed=False)     # two bytes into the last bond record
    assert mol is not None and mol.bond_count == 4
    assert problems == [
        'the header declares 5 bond(s) and the buffer holds 4; the rest of the block was not read',
        'the header declares 1 stereo configuration(s) and the buffer holds 0; the rest of the block '
        'was not read']


def test_a_stereo_count_with_no_atoms_to_name_is_reported():
    mol, problems = pach_load(_v3_header(4, atoms=0, bonds=0, stereo=1) + bytes(9), compressed=False)
    assert mol is not None and mol.atom_count == 0
    assert any('no atoms for them to name' in p for p in problems)


def test_a_bond_count_with_no_atoms_to_name_is_reported():
    """Four counted blocks and one rule, so the bond count is named there too, and first."""
    mol, problems = pach_load(_v3_header(4, atoms=0, bonds=5, stereo=1) + bytes(34), compressed=False)
    assert mol is not None and mol.atom_count == 0
    assert problems == ['the header declares 5 bond(s) and no atoms for them to name; none were read',
                        'the header declares 1 stereo configuration(s) and no atoms for them to name; '
                        'none were read']


# ----- encoder and decoder tests: the enhanced-stereo block -----

def test_an_and_group_round_trips_with_its_index():
    from chython.core import STEREO_AND
    mol = read_smiles('N[C@@H](C)C(=O)O')
    mol.set_stereo_group(mol.atom_numbers[1], STEREO_AND, 3)
    raw = mol.pack(compressed=False)
    assert raw[8:10] == (1).to_bytes(2, 'little')
    assert raw[-3:] == (1).to_bytes(2, 'little') + bytes([0xc3])     # kind 3 << 6 | group 3
    back, problems = pach_load(raw, compressed=False)
    assert problems == []
    assert back.stereo_group_of(back.atom_numbers[1]) == (STEREO_AND, 3)


def test_an_abs_group_needs_no_index():
    from chython.core import STEREO_ABS
    mol = read_smiles('N[C@@H](C)C(=O)O')
    mol.set_stereo_group(mol.atom_numbers[1], STEREO_ABS)
    raw = mol.pack(compressed=False)
    assert raw[-1] == 0x40
    back, problems = pach_load(raw, compressed=False)
    assert problems == []
    assert back.stereo_group_of(back.atom_numbers[1]) == (STEREO_ABS, 0)


def test_a_group_on_an_atom_owning_no_unit_survives():
    """`set_stereo_group` accepts any atom, so the block must carry an entry a unit record could not."""
    from chython.core import STEREO_OR
    mol = read_smiles('CCO')
    mol.set_stereo_group(mol.atom_numbers[2], STEREO_OR, 7)
    back, problems = pach_load(mol.pack(compressed=False), compressed=False)
    assert problems == []
    assert back.stereo_group_of(back.atom_numbers[2]) == (STEREO_OR, 7)


def test_every_grouped_atom_gets_its_own_entry_in_atom_order():
    from chython.core import STEREO_ABS, STEREO_AND, STEREO_OR
    mol = read_smiles('CC(N)C(O)C')
    mol.set_stereo_group(mol.atom_numbers[4], STEREO_OR, 2)
    mol.set_stereo_group(mol.atom_numbers[1], STEREO_AND, 9)
    mol.set_stereo_group(mol.atom_numbers[3], STEREO_ABS)
    raw = mol.pack(compressed=False)
    assert raw[8:10] == (3).to_bytes(2, 'little')
    assert raw[-9:] == bytes([1, 0, 0xc9, 3, 0, 0x40, 4, 0, 0x82])
    back, problems = pach_load(raw, compressed=False)
    assert problems == []
    assert back.stereo_groups() == mol.stereo_groups()


def test_the_declared_length_covers_the_group_block():
    """The allocation's upper bound and the `buf[:at]` prefix both grew by the block, so the header's
    arithmetic and the buffer still agree."""
    from chython.core import STEREO_AND
    mol = read_smiles('N[C@@H](C)C(=O)O')
    mol.set_stereo_group(mol.atom_numbers[1], STEREO_AND, 3)
    raw = mol.pack(compressed=False)
    assert pach_record_length(raw, compressed=False) == len(raw) == 12 + 6 * 3 + 5 * 5 + 9 + 3


def test_a_molecule_with_no_groups_writes_no_block():
    raw = read_smiles('CCO').pack(compressed=False)
    assert raw[8:10] == b'\x00\x00'
    assert len(raw) == 12 + 3 * 3 + 2 * 5


def test_dropping_stereo_drops_the_groups_too():
    from chython.core import STEREO_AND
    mol = read_smiles('N[C@@H](C)C(=O)O')
    mol.set_stereo_group(mol.atom_numbers[1], STEREO_AND, 3)
    raw = mol.pack(compressed=False, drop=['stereo'])
    assert raw[6:8] == b'\x00\x00' and raw[8:10] == b'\x00\x00'
    assert len(raw) == 12 + 6 * 3 + 5 * 5


def test_a_group_entry_naming_an_absent_atom_is_dropped_and_reported():
    raw = _v3_header(4, atoms=3, bonds=2, sgroups=1) \
        + bytes([6, 0, 0x43, 6, 0, 0x42, 8, 0, 0x41]) \
        + bytes([0, 0, 1, 0, 0x01]) + bytes([1, 0, 2, 0, 0x01]) \
        + (9).to_bytes(2, 'little') + bytes([0xc3])
    mol, problems = pach_load(raw, compressed=False)
    assert mol is not None and not mol.has_stereo_groups
    assert problems == ['enhanced-stereo entry 0 names atom index 9 and this record has 3 atom(s); the '
                        'entry was dropped']


def test_an_or_entry_with_no_group_index_is_dropped_and_reported():
    raw = _v3_header(4, atoms=1, sgroups=1) + bytes([6, 0, 0x44]) \
        + (0).to_bytes(2, 'little') + bytes([0x80])          # kind 2, group 0
    mol, problems = pach_load(raw, compressed=False)
    assert mol.stereo_group_of(mol.atom_numbers[0]) == (0, 0)
    assert problems == ['enhanced-stereo entry 0 is or and states no group index, and 1 to 63 is what '
                        'one takes; the entry was dropped']


def test_an_abs_entry_carrying_a_group_index_keeps_the_kind_and_reports():
    raw = _v3_header(4, atoms=1, sgroups=1) + bytes([6, 0, 0x44]) \
        + (0).to_bytes(2, 'little') + bytes([0x45])          # kind 1, group 5
    from chython.core import STEREO_ABS
    mol, problems = pach_load(raw, compressed=False)
    assert mol.stereo_group_of(mol.atom_numbers[0]) == (STEREO_ABS, 0)
    assert problems == ['enhanced-stereo entry 0 is abs and carries group index 5, which only or and '
                        'and take; read as abs with none']


def test_an_entry_stating_no_kind_and_no_group_is_the_default_and_not_damage():
    """A writer that emits a zero byte has stated the default, which is what an absent entry states."""
    raw = _v3_header(4, atoms=1, sgroups=1) + bytes([6, 0, 0x44]) \
        + (0).to_bytes(2, 'little') + bytes([0x00])
    mol, problems = pach_load(raw, compressed=False)
    assert problems == []
    assert mol.stereo_group_of(mol.atom_numbers[0]) == (0, 0) and not mol.has_stereo_groups


def test_an_entry_stating_a_group_index_with_no_kind_is_dropped_and_reported():
    raw = _v3_header(4, atoms=1, sgroups=1) + bytes([6, 0, 0x44]) \
        + (0).to_bytes(2, 'little') + bytes([0x05])          # kind 0, group 5
    mol, problems = pach_load(raw, compressed=False)
    assert not mol.has_stereo_groups
    assert any('states no kind and group index 5' in p for p in problems)


def test_a_second_entry_for_one_atom_is_dropped_and_reported():
    from chython.core import STEREO_AND
    raw = _v3_header(4, atoms=1, sgroups=2) + bytes([6, 0, 0x44]) \
        + (0).to_bytes(2, 'little') + bytes([0xc3]) \
        + (0).to_bytes(2, 'little') + bytes([0x81])
    mol, problems = pach_load(raw, compressed=False)
    assert mol.stereo_group_of(mol.atom_numbers[0]) == (STEREO_AND, 3)
    assert any('repeats atom index 0' in p for p in problems)


def test_a_truncated_group_block_keeps_the_entries_it_has():
    from chython.core import STEREO_AND
    raw = _v3_header(4, atoms=2, sgroups=2) + bytes([6, 0, 0x44, 6, 0, 0x44]) \
        + (0).to_bytes(2, 'little') + bytes([0xc3]) + (1).to_bytes(2, 'little')
    mol, problems = pach_load(raw, compressed=False)
    assert mol.stereo_group_of(mol.atom_numbers[0]) == (STEREO_AND, 3)
    assert mol.stereo_group_of(mol.atom_numbers[1]) == (0, 0)
    assert problems == ['the header declares 2 enhanced-stereo entr(ies) and the buffer holds 1; the '
                        'rest of the block was not read']


def test_a_group_count_with_no_atoms_to_name_is_reported():
    mol, problems = pach_load(_v3_header(4, atoms=0, sgroups=1) + bytes(3), compressed=False)
    assert mol is not None and mol.atom_count == 0
    assert any('enhanced-stereo entr(ies) and no atoms' in p for p in problems)


def test_a_group_block_shortfall_costs_its_own_entry_and_not_the_stereo_record_ahead_of_it():
    """The stereo record's nine declared bytes are all in the buffer, so the configuration is read; the
    entry behind it is not, and that is the whole report.

    Both cuts answer the same, which is the rule stated: a block starts where the header's counts put
    it, so how the tail was lost does not move a block that arrived whole.  What the parity assertion
    gates is the stereo block reserving the group block's declared bytes -- a reserve reads 0 of its 1
    declared configuration here and discards a record that is in the buffer whole.
    """
    from chython.core import STEREO_AND
    mol = read_smiles('N[C@@H](C)C(=O)O')
    mol.set_stereo_group(mol.atom_numbers[0], STEREO_AND, 3)
    raw = mol.pack(compressed=False)
    assert raw[-3:] == bytes([0, 0, 0xc3])
    for cut in (raw[:-1], raw[:-3]):
        back, problems = pach_load(cut, compressed=False)
        assert back.parity_of(back.atom_numbers[1]) == mol.parity_of(mol.atom_numbers[1]) != 0
        assert not back.has_stereo_groups
        assert problems == ['the header declares 1 enhanced-stereo entr(ies) and the buffer holds 0; '
                            'the rest of the block was not read']


def test_a_truncated_bond_block_leaves_one_shortfall_sentence_per_block_behind_it():
    """One clip rule for four blocks, so a bond block that ended the buffer states the same shortfall
    three times over.

    `not mol.has_stereo_groups` is the load-bearing line and this is the tree's only gate on the
    declared offsets: reading the group block from where the bond loop stopped instead of from
    `groups_at` fabricates an entry out of the truncated bond record's own bytes.
    """
    from chython.core import STEREO_AND
    mol = read_smiles('N[C@@H](C)C(=O)O')
    mol.set_stereo_group(mol.atom_numbers[1], STEREO_AND, 3)
    mol, problems = pach_load(mol.pack(compressed=False)[:53], compressed=False)
    assert mol.bond_count == 4 and not mol.has_stereo_groups
    assert problems == [
        'the header declares 5 bond(s) and the buffer holds 4; the rest of the block was not read',
        'the header declares 1 stereo configuration(s) and the buffer holds 0; the rest of the block '
        'was not read',
        'the header declares 1 enhanced-stereo entr(ies) and the buffer holds 0; the rest of the block '
        'was not read']


def test_a_group_count_that_agrees_with_the_record_stops_before_the_map_block():
    """A group block stops at its DECLARED count, so a header that agrees with what was written ends
    the block before the map bytes and a byte missing from the map block costs a map number and not the
    entry ahead of it.

    The other side of that count: the header's counts define the record's layout, so the same two atoms
    with `sgroups=2` read the map bytes as a second entry, put both atoms in AND group 3, take their map
    numbers out of whatever follows and report nothing.  A concatenated store is a record followed by
    more bytes, and telling one from the other is what the count is for.
    """
    from chython.core import STEREO_AND
    raw = _v3_header(4, flags=1, atoms=2, sgroups=1) + bytes([6, 0, 0x44, 6, 0, 0x44]) \
        + bytes([0, 0, 0xc3]) + bytes([1, 0, 0xc3, 0])
    assert len(raw) == 25 == pach_record_length(raw, compressed=False)
    over = bytearray(raw + bytes([7, 0, 9, 0]))
    over[8:10] = (2).to_bytes(2, 'little')
    mol, problems = pach_load(bytes(over), compressed=False)
    assert problems == []      # the count defines the layout, so this is the format's answer
    assert [mol.stereo_group_of(k) for k in mol.atom_numbers] == [(STEREO_AND, 3), (STEREO_AND, 3)]
    assert [mol.map_number_of(k) for k in mol.atom_numbers] == [1792, 2304]
    mol, problems = pach_load(raw, compressed=False)
    assert problems == []
    assert mol.stereo_group_of(mol.atom_numbers[0]) == (STEREO_AND, 3)
    assert mol.stereo_group_of(mol.atom_numbers[1]) == (0, 0)
    assert [mol.map_number_of(k) for k in mol.atom_numbers] == [1, 195]
    mol, problems = pach_load(raw[:-1], compressed=False)
    assert mol.stereo_group_of(mol.atom_numbers[0]) == (STEREO_AND, 3)
    assert [mol.map_number_of(k) for k in mol.atom_numbers] == [1, 0]
    assert problems == ['the header declares 2 map number(s) and the buffer holds 1; the rest of the '
                        'block was not read']


def test_dropping_stereo_groups_alone_keeps_the_configurations():
    """`drop=['stereo_groups']` is the name the version 0 and 2 writers' refusal tells a caller to
    pass, so the version 3 and 4 writer honours it as well as the blanket `drop=['stereo']`."""
    from chython.core import STEREO_AND
    mol = read_smiles('N[C@@H](C)C(=O)O')
    mol.set_stereo_group(mol.atom_numbers[1], STEREO_AND, 3)
    raw = mol.pack(compressed=False, drop=['stereo_groups'])
    assert raw[6:8] == (1).to_bytes(2, 'little') and raw[8:10] == b'\x00\x00'
    assert len(raw) == 12 + 6 * 3 + 5 * 5 + 9
    back, problems = pach_load(raw, compressed=False)
    assert problems == [] and not back.has_stereo_groups
    assert back.parity_of(back.atom_numbers[1]) == mol.parity_of(mol.atom_numbers[1]) != 0


# ----- encoder and decoder tests: the map block -----

def test_a_mapped_molecule_round_trips():
    mol = read_smiles('[CH3:1][OH:2]')
    raw = mol.pack(compressed=False)
    assert raw[1] & 0x01
    assert len(raw) == 12 + 2 * 3 + 1 * 5 + 2 * 2
    assert raw[-4:] == (1).to_bytes(2, 'little') + (2).to_bytes(2, 'little')
    back, problems = pach_load(raw, compressed=False)
    assert problems == []
    assert [back.map_number_of(k) for k in back.atom_numbers] == [1, 2]


def test_a_partially_mapped_molecule_round_trips():
    mol = read_smiles('[CH3:1]O')
    back, problems = pach_load(mol.pack(compressed=False), compressed=False)
    assert problems == []
    assert [back.map_number_of(k) for k in back.atom_numbers] == [1, 0]


def test_an_unmapped_molecule_writes_no_map_block():
    raw = read_smiles('CCO').pack(compressed=False)
    assert raw[1] & 0x01 == 0
    assert len(raw) == 12 + 3 * 3 + 2 * 5


def test_dropping_map_numbers_clears_the_flag():
    raw = read_smiles('[CH3:1][OH:2]').pack(compressed=False, drop=['map_number'])
    assert raw[1] & 0x01 == 0
    assert len(raw) == 12 + 2 * 3 + 1 * 5


def test_a_map_number_beyond_the_arena_is_read_as_none_and_reported():
    raw = _v3_header(4, flags=1, atoms=1) + bytes([6, 0, 0x44]) + (60000).to_bytes(2, 'little')
    mol, problems = pach_load(raw, compressed=False)
    assert mol.map_number_of(mol.atom_numbers[0]) == 0
    assert any('map number' in p for p in problems)


def test_a_truncated_map_block_keeps_the_numbers_it_has():
    raw = read_smiles('[CH3:1][OH:2]').pack(compressed=False)
    mol, problems = pach_load(raw[:-2], compressed=False)
    assert [mol.map_number_of(k) for k in mol.atom_numbers] == [1, 0]
    assert problems == ['the header declares 2 map number(s) and the buffer holds 1; the rest of the '
                        'block was not read']


def test_a_clipped_stereo_count_still_lays_out_the_parity_segment():
    """`want_parity` is the header's DECLARED count and not what resolves, because the arena's
    persistent block is laid out before `_pach3_apply_stereo` has a graph to resolve against.

    Stated as a size because that is the only place it shows: the alanine's clipped record carries the
    segment its header asked for and reads the same length as the whole one, where the same molecule
    drawn without stereo carries no segment and reads shorter by it.
    """
    raw = read_smiles('N[C@@H](C)C(=O)O').pack(compressed=False)
    clipped, problems = pach_load(raw[:-4], compressed=False)
    whole, _ = pach_load(raw, compressed=False)
    flat, _ = pach_load(read_smiles('NC(C)C(=O)O').pack(compressed=False), compressed=False)
    assert any('stereo configuration(s) and the buffer holds 0' in p for p in problems)
    assert clipped.parity_of(clipped.atom_numbers[1]) == 0
    assert len(clipped.to_bytes()) == len(whole.to_bytes()) > len(flat.to_bytes())


def test_an_over_declared_bond_count_shifts_the_blocks_behind_it():
    """The header's counts define the layout, so a bond count one too high moves the stereo block by a
    record and each block states its own shortfall.

    The fabricated fifth bond is read out of the stereo record's bytes -- it is the cost of trusting the
    header, and the alternative is a block that reserves what follows it and so discards readable records
    behind a single missing byte.
    """
    mol = read_smiles('N[C@@H](C)C(=O)O')
    raw = bytearray(mol.pack(compressed=False))
    raw[4:6] = (6).to_bytes(2, 'little')                  # declare one bond more than is written
    del raw[12 + 6 * 3:12 + 6 * 3 + 5]                    # and drop bond #0, so the fabrication is no repeat
    back, problems = pach_load(bytes(raw), compressed=False)
    assert back.bond_count == 5
    assert problems == ['the header declares 6 bond(s) and the buffer holds 5; the rest of the block was '
                        'not read',
                        'the header declares 1 stereo configuration(s) and the buffer holds 0; the rest '
                        'of the block was not read']


def test_one_missing_byte_costs_one_record_and_not_the_blocks_behind_it():
    """The reason a block does not reserve the bytes the blocks behind it declare: one byte off the end of
    an 82 byte record costs the one group entry it falls in and nothing else."""
    from chython.core import STEREO_OR
    mol = read_smiles('N[C@@H](C)C(=O)O')
    for k in mol.atom_numbers:
        mol.set_stereo_group(k, STEREO_OR, 1)
    raw = mol.pack(compressed=False)
    back, problems = pach_load(raw[:-1], compressed=False)
    assert back.parity_of(back.atom_numbers[1]) == 1       # the stereo block is whole and was read
    assert back.has_stereo_groups
    assert len([k for k in back.atom_numbers if back.stereo_group_of(k) != (0, 0)]) == 5
    assert problems == ['the header declares 6 enhanced-stereo entr(ies) and the buffer holds 5; the '
                        'rest of the block was not read']


def _rich_record():
    """One record exercising every block: coordinates, bonds, a centre, a group, map numbers."""
    from chython.core import STEREO_AND
    mol = read_smiles('[NH2:1][C@@H:2]([CH3:3])[C:4](=[O:5])[OH:6]')
    for i, k in enumerate(mol.atom_numbers):
        mol.set_xy(k, i * 1.5, 0.0)
    mol.set_stereo_group(mol.atom_numbers[1], STEREO_AND, 1)
    return mol.pack(compressed=False)


def test_a_zero_atom_record_is_an_empty_molecule():
    mol, problems = pach_load(_v3_header(4), compressed=False)
    assert problems == []
    assert len(mol) == 0


def test_the_reserved_header_bytes_are_reported_and_ignored():
    raw = bytearray(_v3_header(4, atoms=1) + bytes([6, 0, 0x44]))
    raw[10] = 0x01
    mol, problems = pach_load(bytes(raw), compressed=False)
    assert len(mol) == 1
    assert any('reserved' in p for p in problems)


def test_trailing_bytes_belong_to_the_next_record():
    raw = _rich_record()
    first, problems = pach_load(raw + raw, compressed=False)
    assert problems == []
    assert pach_record_length(raw + raw, compressed=False) == len(raw)
    assert first == pach_load(raw, compressed=False)[0]


def test_no_truncation_of_a_rich_record_raises_and_a_prefix_answers_iff_the_atoms_fit():
    """The atom block is the one all-or-nothing part, so a prefix answers with a molecule exactly when
    the atom block fits in it, and never silently."""
    raw = _rich_record()
    head = 12 + int.from_bytes(raw[2:4], 'little') * (9 if raw[0] == 3 else 3)
    for cut in range(len(raw) + 1):
        mol, problems = pach_load(raw[:cut], compressed=False)
        assert mol is not None or problems, cut
        assert (mol is not None) == (cut >= head), cut
        assert problems or cut == len(raw), cut


def test_no_single_byte_corruption_raises_and_no_payload_byte_costs_the_molecule():
    """A payload byte cannot make the record unreadable: only the header's own fields decide whether
    there is a molecule at all, so every mutation from byte 12 on answers with one."""
    raw = _rich_record()
    for i in range(len(raw)):
        for value in (0x00, 0x7f, 0xff):
            mutant = bytearray(raw)
            mutant[i] = value
            mol, problems = pach_load(bytes(mutant), compressed=False)
            assert mol is not None or problems, (i, value)
            assert mol is not None or i < 12, (i, value)


def test_a_declared_count_does_not_size_the_scratch_block():
    """The bond regions are sized from what the buffer holds, so a 15 byte record cannot ask for three
    quarters of a megabyte of scratch on the read path."""
    raw = _v3_header(4, atoms=1, bonds=65535, stereo=65535, sgroups=65535) + bytes([6, 0, 0x44])
    pach_load(raw, compressed=False)                       # warm the import-time allocations
    tracemalloc.start()
    pach_load(raw, compressed=False)
    peak = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()
    assert peak < 64 * 1024, peak                          # 2,433 measured; 854,393 from the declared count


def test_unpack_raises_where_pach_load_reports():
    raw = _v3_header(4, atoms=3)                              # the atom block is not there at all
    mol, problems = pach_load(raw, compressed=False)
    assert mol is None and problems
    with raises(ValueError, match='atom block'):
        MoleculeContainer.unpack(raw, compressed=False)


def test_the_atom_count_ceiling_is_refused_by_name():
    mol = MoleculeContainer()
    with mol.edit() as e:
        for _ in range(65536):
            e.add_atom('C')
    with raises(ValueError, match='atom count is a 16 bit field'):
        mol.pack(compressed=False)


def test_the_bond_count_ceiling_is_refused_by_name():
    mol = MoleculeContainer()
    with mol.edit() as e:
        ids = [e.add_atom('C') for _ in range(65535)]
        for a, b in zip(ids, ids[1:]):
            e.add_bond(a, b, 1)
        e.add_bond(ids[0], ids[1000], 1)                      # two closures past the chain's 65534
        e.add_bond(ids[5], ids[2000], 1)
    with raises(ValueError, match='bond count is a 16 bit field'):
        mol.pack(compressed=False)


def complete_graph_record(n=150):
    """A well-formed version 4 record for the complete graph on `n` carbons.

    No container can produce it -- `edit()`'s seal derives the same rings the decode does, so the
    molecule cannot be built to be packed -- and there is nothing damaged about the bytes: every count
    agrees with every block.  `n=150` is 11175 bonds and 56337 bytes, and reaches `perceive_rings`'
    relevant-cycle prototype limit in a fraction of a second.  The wall-clock deadline is the other
    resource refusal and is not usable in a test.
    """
    bonds = [(a, b) for a in range(n) for b in range(a + 1, n)]
    out = bytearray(_v3_header(4, atoms=n, bonds=len(bonds)))
    for _ in range(n):
        out += bytes([6, 0x80, 0x40])                          # carbon, count pinned at 0, neutral
    for a, b in bonds:
        out += a.to_bytes(2, 'little') + b.to_bytes(2, 'little') + b'\x01'
    return bytes(out)


def test_a_graph_the_ring_perception_refuses_is_no_molecule_and_a_sentence():
    """The decode path is not an answer boundary: `perceive_rings` raises where `pach_load` reports.

    `rebuild_derived` ends every builder, and a record can state a graph its resource limits refuse
    without stating one damaged byte.  Answering `None` with the reason is what the door's contract
    already promises; raising is what the input-posture rule forbids.
    """
    mol, problems = pach_load(complete_graph_record(), compressed=False)
    assert mol is None
    assert any('could not be derived' in p and 'prototype limit' in p for p in problems)


def test_unpack_raises_where_that_record_reports():
    """The same bytes at the answer boundary, which has no way to say "unknown"."""
    with raises(ValueError, match='prototype limit'):
        MoleculeContainer.unpack(complete_graph_record(), compressed=False)


def test_dropping_coordinates_beats_an_explicit_version_three():
    """`drop=` is the waiver that wins: it selects version 4 whatever `version=` asked for.

    Otherwise the loop does not close -- a coordinate outside the int24 range is refused with advice to
    pass `drop=['coordinates']`, and that call would come back version 3 with the drawing intact.
    """
    mol, _ = _drawn_amino_propanol()
    raw = mol.pack(compressed=False, version=3, drop=['coordinates'])
    assert raw[0] == 4
    back, problems = pach_load(raw, compressed=False)
    assert problems == []
    assert not back.has_coordinates
    assert back == mol


def test_the_refusals_advice_is_a_call_that_works():
    mol, n = _drawn_amino_propanol()
    with mol.edit() as e:
        e.set_xy(n[0], 1000.0, 0.0)
    with raises(ValueError, match="drop=\\['coordinates'\\]"):
        mol.pack(compressed=False, version=3)
    assert mol.pack(compressed=False, version=3, drop=['coordinates'])[0] == 4


# ----- frozen corpus tests: the pinned version 3 and version 4 records -----

def _fields(mol, with_xy):
    """The molecule as plain data, for comparing one version's decode against another's."""
    index = {k: i for i, k in enumerate(mol.atom_numbers)}
    out = [[a.element, a.r_index, a.isotope, a.charge, int(a.is_radical), a.implicit_h]
           for a in mol.atoms()]
    bonds = sorted((index[b.n], index[b.m], b.order) for b in mol.bonds())
    stereo = sorted((u['kind'], index[u['anchor']],
                     tuple(sorted(index[r] for r in u['refs'] if r is not None)), u['parity'])
                    for u in mol.stereo_units() if u['parity'])
    xy = [list(a.xy) if a.xy is not None else None for a in mol.atoms()] if with_xy else None
    return (out, bonds, stereo, xy)


@mark.parametrize('path,version', [(V3_PATH, 3), (V4_PATH, 4)])
def test_the_fixture_is_present_and_holds_every_builder(path, version):
    records = load_corpus(path)
    assert [name for name, _, _ in records] == [name for name, _ in BUILDERS]
    assert {data[0] for _, data, _ in records} == {version}


@mark.parametrize('path,version', [(V3_PATH, 3), (V4_PATH, 4)])
def test_the_writer_reproduces_every_pinned_record(path, version):
    built = dict(BUILDERS)
    for name, data, _ in load_corpus(path):
        mol = built[name]() if version == 4 else drawn(built[name]())
        assert mol.pack(compressed=False, version=version) == data, name


@mark.parametrize('path', [V3_PATH, V4_PATH])
def test_every_record_decodes_to_the_answers_it_was_written_with(path):
    for name, data, expected in load_corpus(path):
        mol, problems = pach_load(data, compressed=False)
        assert problems == [], name
        assert answers(mol, path is V3_PATH) == expected, name


@mark.parametrize('path', [V3_PATH, V4_PATH])
def test_decoding_and_re_encoding_is_byte_identical(path):
    for name, data, _ in load_corpus(path):
        mol, _ = pach_load(data, compressed=False)
        assert mol.pack(compressed=False, version=data[0]) == data, name


@mark.parametrize('path', [V3_PATH, V4_PATH, V0_PATH, V0_NATIVE_PATH, V2_PATH])
def test_pach_record_length_agrees_with_the_true_length(path):
    for name, data, _ in load_corpus(path):
        assert pach_record_length(data, compressed=False) == len(data), name


@mark.parametrize('path', [V0_PATH, V0_NATIVE_PATH, V2_PATH])
def test_every_legacy_record_re_encodes_as_version_3_and_4_and_agrees(path):
    """Each legacy record decoded, written as both new versions, decoded again, compared.

    Not `problems == []`: a legacy corpus holds records whose own writer left a field underivable, and
    a report about one of those is the reader working.  `mol is not None` is the bar.
    """
    _expected_refused = {V0_PATH: ['xy:14', 'xy:15', 'xy:16', 'xy:17', 'xy:18'],
                         V2_PATH: ['xy:14', 'xy:15', 'xy:16', 'xy:17', 'xy:18'],
                         V0_NATIVE_PATH: []}
    refused = []
    for name, data, _ in load_corpus(path):
        mol, _ = pach_load(data, compressed=False)
        assert mol is not None, name
        try:
            three, problems = pach_load(mol.pack(compressed=False, version=3), compressed=False)
        except ValueError as err:
            # a version 0/2 coordinate can sit outside version 3's field; the chemistry still travels
            assert 'pach coordinate field' in str(err), name
            refused.append(name)
        else:
            assert problems == [], name
            assert _fields(three, True) == _fields(mol, True), name
        four, problems = pach_load(mol.pack(compressed=False, version=4), compressed=False)
        assert problems == [], name
        assert _fields(four, False) == _fields(mol, False), name
    assert refused == _expected_refused[path]


def test_no_configured_direction_list_has_two_unnamed_directions():
    """Why four slots are complete: a stored configuration always has one implied direction at most.

    `_list_has_two_unnamed` refuses to call a list with two unnamed directions stereogenic, so every
    configuration the arena holds names all but one direction per list and every slot holds a real atom
    index -- which is what makes the record's frame total and sentinels unnecessary.
    """
    built = dict(BUILDERS)
    checked = 0
    for name, _ in BUILDERS:
        for unit in built[name]().stereo_units():
            if not unit['parity']:
                continue
            refs = unit['refs']
            lists = [refs] if unit['kind'] == 0 else [refs[:2], refs[2:]]
            for one in lists:
                assert sum(1 for r in one if r is None) <= 1, (name, unit)
            checked += 1
    assert checked, 'no builder states a configuration; the sweep is not exercising the path'


def test_a_two_unnamed_direction_list_is_not_stereogenic():
    """Ethanol's CH2: two heavy neighbours and two hydrogens, so two of its four directions are unnamed
    -- and the arena does not call it stereogenic, so there is no configuration for the record to fail
    to express.  This is the other half of the sweep above: four slots are complete BECAUSE a list that
    needs two implied directions is never configured."""
    mol = read_smiles('CCO')
    unit, = [u for u in mol.stereo_units() if u['anchor'] == mol.atom_numbers[1]]
    assert sum(1 for r in unit['refs'] if r is None) == 2
    assert not unit['stereogenic'] and not unit['parity']
