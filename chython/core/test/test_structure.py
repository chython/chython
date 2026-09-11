# -*- coding: utf-8 -*-
import pytest
from chython.core import _core as _structure


def test_struct_sizes_are_locked():
    # sizeof(StructureHeader), which is the header's MAXIMUM and not its usual size: a buffer's
    # header is 24 + 8 * seg_count and most molecules stop the table at three entries.
    assert _structure._header_size() == 128
    assert _structure._atom_record_size() == 24
    assert _structure._halfedge_size() == 8


def test_the_sgroup_index_bound_is_one_below_its_sentinel():
    """`SGROUP_NO_INDEX` spends the top u16 value, so the greatest REAL index is 0xFFFE.

    Written before any S-group validator exists, which is the whole point.  The identical
    shape -- `H_UNKNOWN` against `H_IMPLICIT_MAX` -- was diagnosed only AFTER a reader had
    derived its bound from the nibble's WIDTH and so admitted the sentinel as a count, on
    three write paths, the worst of them a valence clamp that put a computed number onto the
    value meaning "nobody could compute it".  That family has now cost three commits across
    two epics.

    Asserted as the RELATION and not as two numbers: two literals are two things that can
    drift, and the fact being pinned is that the bound sits one below the sentinel, whatever
    either happens to be.  Any validator for `index`, `ext_index` or `parent` cites
    SGROUP_INDEX_MAX; none of them may cite the field's width.
    """
    assert _structure.SGROUP_INDEX_MAX == _structure.SGROUP_NO_INDEX - 1
    assert _structure.SGROUP_NO_INDEX == 0xFFFF, 'the sentinel is the top u16 value'
    assert _structure.SGROUP_INDEX_MAX == 0xFFFE, \
        'and the bound is NOT 0xFFFF -- a width is not a bound'
    # Read as module attributes ON PURPOSE: `globals().update(_segment_ids())` is what publishes
    # them, so this is the surface a caller outside the core actually sees.  §6.3 -- an unexported
    # domain does not get asked about, it gets re-derived wrongly, and this bound has a history.


def test_segment_ids_are_dense_and_unique():
    """Dense, and split into a persistent block and a derived block IN THAT ORDER.

    The order is load-bearing rather than tidy.  `structure_append` refuses any id below
    SEG_PERSISTENT_COUNT and `structure_alloc_full` lays out only the ids below it, so "persistent"
    and "derived" are decided by a single comparison against one number.  Were the two interleaved,
    every one of those tests would have to become a switch.  The persistent PREFIX is also frozen
    across releases -- a v3 buffer's ids 0..4 mean in v4 what they meant in v3 -- while the derived
    ids may be renumbered freely, because nothing serialises them.
    """
    persistent = [_structure.SEG_ATOMS, _structure.SEG_CSR_PTR, _structure.SEG_CSR_EDGE,
                  _structure.SEG_XY, _structure.SEG_STEREO_GROUPS, _structure.SEG_SGROUP_RECORD,
                  _structure.SEG_SGROUP_INDEX, _structure.SEG_OPAQUE_BLOB,
                  _structure.SEG_CONFORMERS, _structure.SEG_PARITY]
    derived = [_structure.SEG_RING_BITS, _structure.SEG_RELEVANT_RINGS, _structure.SEG_FEATURES,
               _structure.SEG_ELEMENT_INDEX, _structure.SEG_EDGE_WORD,
               _structure.SEG_COMPONENT_LABEL, _structure.SEG_STEREO_UNIT]
    assert sorted(persistent + derived) == list(range(len(persistent) + len(derived)))
    assert _structure._segment_count() == len(persistent) + len(derived)
    assert _structure._persistent_segment_count() == len(persistent)
    assert max(persistent) < min(derived), 'the persistent ids must be a prefix, not interleaved'


def test_parity_is_a_persistent_segment():
    """The stated parity is storage, not a cache: id 9, inside the persistent block, and in the table.

    Ten persistent ids of the thirteen SEG_TABLE_MAX admits, so three remain for a later release.
    """
    from chython.core._core import SEG_PARITY, SEG_PERSISTENT_COUNT, SEG_TABLE_MAX, SEG_MASK_PARITY

    assert SEG_PARITY == 9
    assert SEG_PARITY < SEG_PERSISTENT_COUNT == 10
    assert SEG_PERSISTENT_COUNT <= SEG_TABLE_MAX == 13
    assert SEG_MASK_PARITY == 4                      # SEG_MASK_XY 1 and SEG_MASK_STEREO 2 are taken


def test_the_arena_states_version_six():
    """The version a fresh buffer stamps, and the three older ones a reader accepts."""
    from chython.core._core import (STRUCT_VERSION, STRUCT_VERSION_V5, STRUCT_VERSION_V4,
                                    STRUCT_VERSION_V3)

    assert (STRUCT_VERSION, STRUCT_VERSION_V5, STRUCT_VERSION_V4, STRUCT_VERSION_V3) == (6, 5, 4, 3)


def test_the_segment_table_has_room_left_and_the_header_did_not_move():
    """13 table entries in v4 and in v5, and 13 of them come to a header the size of v3's.

    Deleting v3's `total_len` field freed exactly the four bytes `seg_count` needed, so the table
    still begins at offset 24 and 13 entries still come to 128 bytes.  That is why a v3 buffer needs
    no payload relocation to be read here -- every persistent offset it names is still correct.  Five
    spare entries was the headroom the format change bought; `SEG_CONFORMERS` and `SEG_PARITY` have
    spent two of them and three remain, and the point of `seg_count` is that running out costs a
    larger header rather than another format version.

    A buffer this build writes does not usually spend 13.  `seg_count` is one past the highest entry
    the molecule uses, so 128 is the header's ceiling -- what this test locks is the ceiling and the
    offset the table starts at, both of which the v3 read path depends on.
    """
    assert _structure._header_size() == 128
    assert _structure._segment_table_max() == 13
    assert _structure._persistent_segment_count() < _structure._segment_table_max(), \
        'the table is full again and a new persistent segment cannot be added'


def test_the_conformer_struct_sizes_are_locked():
    """The two structs `SEG_CONFORMERS` is made of, and the third coordinate is NOT in `xy_t`.

    `xy_t` staying 8 bytes is the load-bearing half.  Widening it to hold z would have been the
    smaller change to write and would have grown every purely 2D molecule's buffer by four bytes per
    atom -- and since `to_bytes()` is a molecule identity, that reprices every stored key for a
    molecule that has no third coordinate at all.  A separate segment costs nothing to a molecule
    that does not use it, which is the reason the crystals design recommended one.
    """
    assert _structure._xy_size() == 8, 'xy_t must not have grown a z'
    assert _structure._xyz_size() == 12
    assert _structure._conformer_record_size() == 4


def test_the_conformer_model_bound_is_one_below_nothing_and_the_no_index_sentinel_is_the_top():
    """Two domains, declared once each, and they are NOT the same shape -- which is the point.

    `CONF_NO_INDEX` spends the top u32 value because `ext_index` must round-trip a file's own model
    number VERBATIM including zero: a PDB `MODEL 0` is representable and nothing forbids it, so zero
    cannot double as "no file number" without making the two indistinguishable.  That is
    `SGROUP_NO_INDEX`'s argument applied to a wider field.

    `CONF_MAX_MODELS` is a DIFFERENT KIND OF NUMBER and deliberately not `CONF_NO_INDEX - 1`.  It
    bounds how many models a molecule may hold, and a `uint32_t` count would admit four billion of
    them -- `M * atom_count * 12` bytes overruns the 4 GiB buffer limit long before that, so the
    width is not the bound.  §6.3: the bound is real, is exported, and a format module refusing a
    trajectory cites it rather than re-deriving one.
    """
    assert _structure.CONF_NO_INDEX == 0xFFFFFFFF, 'the sentinel is the top u32 value'
    assert _structure.CONF_MAX_MODELS == 0xFFFF
    assert _structure.CONF_MAX_MODELS < _structure.CONF_NO_INDEX, \
        'a model count must not be able to reach the value meaning "no file number"'


def test_allocation_reports_magic_and_counts():
    info = _structure._alloc_probe(5, 4, False)
    assert info['magic'] == 0x43485933  # 'CHY3' -- the magic names the project, not the version
    assert info['version'] == _structure.STRUCT_VERSION
    assert info['atom_count'] == 5
    assert info['bond_count'] == 4
    # NOT `>= 128`: the header is 24 + 8 * seg_count and this molecule spends three entries, so it
    # is 48 bytes.  128 is the size of the C struct and the size of the header only for a molecule
    # that uses the last persistent segment.
    assert info['total_len'] >= 24 + 8 * info['seg_count']
    assert info['total_len'] % 8 == 0


def test_wide_index_flag_is_recorded():
    assert _structure._alloc_probe(1, 0, True)['flags'] == 1
    assert _structure._alloc_probe(1, 0, False)['flags'] == 0


def test_empty_molecule_allocates_a_valid_csr_pointer_slot():
    info = _structure._alloc_probe(0, 0, False)
    assert info['atom_count'] == 0
    assert info['bond_count'] == 0
    # A header plus the one csr_ptr entry an atom-less molecule still owns -- the point of the test.
    assert info['total_len'] > 24 + 8 * info['seg_count']
    assert info['total_len'] % 8 == 0
    assert info['persistent_len'] <= info['total_len']


def test_absent_segment_reads_as_zero():
    probe = _structure._zero_page_probe(3, 2)
    assert probe['has_xy'] is False
    assert probe['xy_reads'] == [0, 0, 0, 0, 0, 0]
    assert probe['has_atoms'] is True


def test_atom_flag_and_nibble_packing():
    probe = _structure._atom_field_probe()
    assert probe['implicit_h'] == 3
    assert probe['explicit_h'] == 2
    assert probe['radical'] is True
    assert probe['h_pinned'] is True
    assert probe['in_ring'] is True
    assert probe['hybridization'] == 5
    assert probe['rings_count'] == 2
    assert probe['aromatic_ring_count'] == 1
    # nothing bled into a neighbouring field
    assert probe['element'] == 6
    assert probe['charge'] == -1
    assert probe['isotope'] == 13
    assert probe['map_number'] == 4095
    assert probe['n'] == 7


def test_every_bit_of_the_atom_flags_byte_is_allocated_or_reserved():
    """The flags byte is either allocated or reserved, and this is the test that says so.

    Three single-bit flags hold 0x01, 0x04, 0x40 and hybridization holds bits 3-5 (0x38).
    Bits 1 and 7 are reserved -- every writer leaves them 0, and a version-5 buffer that sets
    one is refused.  A new flag may not take a reserved bit; doing so is a format change that
    reprices every key ever produced by `to_bytes()`.

    What this catches is a new flag quietly taking a bit that is already spoken for: the probe sets
    all active fields at once, so an overlap makes one of them read back wrong.  The last assertion
    is the one worth having: it is the only place a test says a writer leaves the reserved bits alone.
    """
    probe = _structure._atom_field_probe()
    single_bit_flags = 0x45               # radical, in_ring, h_pinned
    hybridization_field = 0x38            # bits 3-5
    reserved_flags = 0x82                 # bits 1 and 7
    assert single_bit_flags & hybridization_field == 0, 'a flag overlaps the hybridization field'
    assert single_bit_flags & reserved_flags == 0, 'a flag took a reserved bit'
    assert single_bit_flags | hybridization_field | reserved_flags == 0xFF, \
        'a bit of the flags byte is neither allocated nor reserved -- say what it is for'
    assert probe['flags'] == single_bit_flags | (5 << 3)
    assert probe['flags'] & reserved_flags == 0, 'a writer set a reserved bit'


def test_h_nibbles_are_isolated_at_the_boundary():
    probe = _structure._atom_h_nibble_probe(15, 15)
    assert probe['implicit_h'] == 15
    assert probe['explicit_h'] == 15
    assert probe['raw'] == 0xff

    probe = _structure._atom_h_nibble_probe(0, 15)
    assert probe['implicit_h'] == 0
    assert probe['explicit_h'] == 15

    probe = _structure._atom_h_nibble_probe(15, 0)
    assert probe['implicit_h'] == 15
    assert probe['explicit_h'] == 0


def test_hybridization_is_range_checked():
    with pytest.raises(ValueError):
        _structure._atom_set_hybridization_probe(7)
    with pytest.raises(ValueError):
        _structure._atom_set_hybridization_probe(0)


def test_csr_build_mirrors_every_bond():
    # propane skeleton: 0-1, 1-2 single; plus 0-2 double to force a ring
    probe = _structure._csr_probe(3, [(0, 1, 1), (1, 2, 1), (0, 2, 2)])
    assert probe['ptr'] == [0, 2, 4, 6]
    assert probe['neighbors'] == {0: [(1, 1), (2, 2)], 1: [(0, 1), (2, 1)],
                                  2: [(0, 2), (1, 1)]}
    assert probe['canonical'] == [(0, 1, 1), (0, 2, 2), (1, 2, 1)]


def test_csr_find_returns_none_for_unbonded():
    probe = _structure._csr_find_probe(4, [(0, 1, 1), (2, 3, 1)])
    assert probe[(0, 1)] == 1
    assert probe[(1, 0)] == 1
    assert probe[(0, 2)] is None
    assert probe[(0, 3)] is None


def test_csr_handles_isolated_atoms():
    probe = _structure._csr_probe(3, [(0, 2, 1)])
    assert probe['ptr'] == [0, 1, 1, 2]
    assert probe['neighbors'] == {0: [(2, 1)], 1: [], 2: [(0, 1)]}


def test_csr_with_no_bonds_never_writes_the_zero_page():
    probe = _structure._csr_probe(3, [])
    assert probe['ptr'] == [0, 0, 0, 0]
    assert probe['neighbors'] == {0: [], 1: [], 2: []}
    assert probe['canonical'] == []


def test_the_parity_segment_is_laid_out_only_when_asked():
    """One byte per atom, 8-aligned, and absent unless the mask asks -- SEG_STEREO_GROUPS' own shape."""
    from chython.core._core import SEG_MASK_PARITY, SEG_PARITY, _alloc_probe

    plain = _alloc_probe(13, 13, False, 0)
    assert plain['lengths'][SEG_PARITY] == 0
    assert plain['seg_count'] <= SEG_PARITY

    asked = _alloc_probe(13, 13, False, SEG_MASK_PARITY)
    assert asked['lengths'][SEG_PARITY] == 16            # align8(13)
    assert asked['seg_count'] == SEG_PARITY + 1
    assert asked['offsets'][SEG_PARITY] % 8 == 0


def test_a_parity_byte_outside_the_domain_is_refused():
    """The byte is one field with three values; a fourth is a buffer this build cannot model."""
    from chython.core._core import MoleculeContainer, SEG_PARITY, _segment_span, read_smiles

    mol = read_smiles('C[C@H](N)C(=O)O')
    raw = bytearray(mol.to_bytes())
    off, ln = _segment_span(mol, SEG_PARITY)
    assert ln, 'a molecule with a stated parity must carry the segment'
    raw[off] = 3
    with pytest.raises(ValueError, match='parity byte 0 states 3'):
        MoleculeContainer.from_bytes(bytes(raw))


def test_the_arena_names_the_version_whose_conformer_record_was_wider():
    """A version-5 buffer's conformer record is 16 bytes wide, so the old stride is a named constant.

    The one site that reads it -- `structure_from_bytes`, checking the segment's declared length --
    checks an EQUALITY, so it must compute the length at the stride the buffer itself states.  A
    literal 16 there would be a second copy of a layout the struct no longer states.
    """
    assert _structure.STRUCT_VERSION_V5 == 5
    assert _structure.CONFORMER_RECORD_V5 == 16
    assert _structure.STRUCT_VERSION_V4 < _structure.STRUCT_VERSION_V5


def test_the_segment_length_takes_the_record_stride():
    """One copy of the layout arithmetic, parameterised by the record width rather than per version.

    2 models of 5 atoms is chosen so both strides land on an 8-boundary already: the function pads to
    align8, so a difference asserted at arbitrary counts would measure the padding rather than the
    stride.  Zero models is zero bytes at any stride -- an absent segment and no models are one state.
    """
    header = _structure._conformer_header_size()
    xyz = 2 * 5 * _structure._xyz_size()
    assert _structure._conformer_seg_len_probe(2, 5, 4) == header + 2 * 4 + xyz
    assert _structure._conformer_seg_len_probe(2, 5, _structure.CONFORMER_RECORD_V5) == \
        header + 2 * _structure.CONFORMER_RECORD_V5 + xyz
    assert _structure._conformer_seg_len_probe(0, 5, _structure.CONFORMER_RECORD_SIZE) == 0
    assert _structure.CONFORMER_RECORD_SIZE == _structure._conformer_record_size(), \
        'the exported width and the struct must be the same number'
