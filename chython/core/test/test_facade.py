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
"""`smiles()`, `pach()`, `unpach()` and `unpack`: the core's bidirectional short doors.

The direction is the argument's type, and each keyword serves one direction and is ignored by the
other -- the contract `ctfile`'s `mol()` already holds.  What the codecs themselves guarantee is
`test_pach.py`, `test_pach3.py` and `test_reaction_pach.py`; what is asserted here is the dispatch,
the two error policies, and that chython 2's three dropped keywords are refused rather than ignored.
"""
import pytest

from .. import (MoleculeContainer, QueryContainer, ReactionContainer, pach, pach_dump, pach_load,
                read_smarts, read_smiles, reaction_pach_dump, smiles, unpach, unpack,
                write_reaction_smiles, write_smiles)


#: Header byte 1 bit 4 is undefined, so setting it is damage the decoder reports and reads past --
#: the one shape that exercises "a structure AND a complaint", which the two error policies differ on.
_UNDEFINED_FLAG = 0x10


def _molecule():
    return read_smiles('C/C=C/C')


def _reaction():
    return read_smiles('[CH3:1][OH:2]>>[CH3:1][NH2:3]')


def _damaged_molecule_record():
    """A version 3 record that decodes to the whole molecule and one complaint."""
    raw = bytearray(_molecule().pack(compressed=False))
    raw[1] ^= _UNDEFINED_FLAG
    mol, problems = pach_load(bytes(raw), compressed=False)
    assert mol is not None and problems, 'the fixture no longer states recoverable damage'
    return bytes(raw)


def _damaged_reaction_record():
    """The same damage in the first molecule record of a reaction record."""
    raw = bytearray(_reaction().pack(compressed=False))
    raw[4 + 1] ^= _UNDEFINED_FLAG
    return bytes(raw)


# ------------------------------------------------------------------------------------- smiles()

def test_a_string_reads_and_a_molecule_writes():
    assert smiles(smiles('CCO')) == write_smiles(read_smiles('CCO'))
    assert isinstance(smiles('CCO'), MoleculeContainer)
    assert isinstance(smiles(smiles('CCO')), str)


def test_an_arrow_makes_it_a_reaction_in_both_directions():
    rxn = smiles('[CH3:1][OH:2]>>[CH3:1][NH2:3]')
    assert isinstance(rxn, ReactionContainer)
    assert smiles(rxn) == write_reaction_smiles(rxn) == rxn.smiles


def test_bytes_read_too_because_read_smiles_takes_them():
    assert smiles(b'CCO').atom_count == 3


def test_a_query_has_no_write_direction_and_says_so():
    query = read_smarts('[C;a]')
    assert isinstance(query, QueryContainer)
    with pytest.raises(TypeError, match='no SMILES form'):
        smiles(query)


def test_the_log_is_positional_because_a_file_loop_spells_it_that_way():
    log = []
    # a ring closure whose two labels state different orders: one line, on both destinations
    molecule = smiles('C-1CCCCC=1', log)
    assert [x.rule for x in log] == ['smiles:ring-bond-order-conflict']
    assert [x.rule for x in molecule.log] == [x.rule for x in log]


def test_the_spec_selects_the_same_string_format_does():
    molecule = _molecule()
    for spec in ('', 'a', 'h', '!s'):
        assert smiles(molecule, spec=spec) == format(molecule, spec)


def test_each_keyword_serves_one_direction_and_the_other_ignores_it():
    """`spec=` on import and `log=` on export are accepted and ignored, as `mol()`'s `version=` is on
    import: a caller reading and writing in one loop passes one keyword set to both calls."""
    assert smiles('CCO', spec='a').atom_count == 3
    log = []
    assert smiles(_molecule(), log) == format(_molecule(), '')
    assert log == [], 'the writer has nothing to report and must not invent a line'


# ---------------------------------------------------------------------------- pach(), both ways

def test_a_container_writes_and_the_bytes_read_back():
    for structure in (_molecule(), _reaction()):
        record = pach(structure)
        assert isinstance(record, bytes)
        assert str(unpach(record)) == str(structure)


def test_the_short_door_writes_exactly_what_pack_writes():
    molecule, rxn = _molecule(), _reaction()
    assert pach(molecule) == molecule.pack() == molecule.pach() == pach_dump(molecule)
    assert pach(rxn) == rxn.pack() == rxn.pach() == reaction_pach_dump(rxn)


def test_the_current_layout_is_what_no_version_writes():
    """Version 3 with coordinates and 4 without, which is `pack()`'s rule and not a second one."""
    flat = _molecule()
    assert not flat.has_coordinates
    assert pach(flat, compressed=False)[0] == 4
    drawn = _molecule()
    drawn.clean2d()
    assert pach(drawn, compressed=False)[0] == 3
    assert pach(drawn, compressed=False, version=4)[0] == 4


def test_a_stated_version_is_obeyed_in_both_eras():
    molecule, rxn = _molecule(), _reaction()
    assert pach(molecule, compressed=False, version=2, drop='*')[0] == 2
    assert pach(rxn, compressed=False, version=1, drop='*')[0] == 1
    assert unpach(pach(molecule, version=2, drop='*')).atom_count == 4
    assert isinstance(unpach(pach(rxn, version=1, drop='*')), ReactionContainer)


def test_compressed_is_the_export_default_and_false_writes_the_raw_record():
    molecule = _molecule()
    assert pach(molecule, compressed=False)[0] in (3, 4)
    assert pach(molecule)[0] not in (0, 1, 2, 3, 4, 5, 0x33), 'a zlib header, not a version byte'
    assert pach(molecule) == pach(molecule, compressed=None)


def test_drop_reaches_the_encoder_rather_than_being_swallowed():
    molecule = _molecule()
    molecule.set_title(b'x')
    with pytest.raises(ValueError, match='title'):
        pach(molecule)
    assert unpach(pach(molecule, drop=['title'])).title == ''


def test_the_import_direction_ignores_the_export_keywords():
    record = pach(_molecule())
    assert pach(record, version=2, drop='*').atom_count == 4


def test_a_query_has_no_record():
    with pytest.raises(TypeError, match='no record for a query'):
        pach(read_smarts('[C;a]'))


# ------------------------------------------------------------------- what byte 0 dispatches on

def test_all_three_eras_and_both_shapes_come_back_through_one_door():
    molecule, rxn = _molecule(), _reaction()
    cases = {3: molecule.pack(), 4: molecule.pack(version=4), 2: molecule.pack(version=2, drop='*'),
             5: rxn.pack(), 1: rxn.pack(version=1, drop='*')}
    for version, record in cases.items():
        obj = unpach(record)
        assert isinstance(obj, ReactionContainer if version in (1, 5) else MoleculeContainer), version
    # and the arena, whose first byte is the magic's low byte rather than a version
    arena = molecule.to_bytes()
    assert arena[0] == 0x33
    assert unpach(arena).canonical_bytes == molecule.canonical_bytes


def test_raw_and_compressed_are_both_read_with_nothing_declared():
    molecule = _molecule()
    for record in (molecule.pack(), molecule.pack(compressed=False)):
        assert unpach(record).atom_count == 4


def test_compressed_states_it_and_a_buffer_that_disagrees_is_refused():
    molecule = _molecule()
    with pytest.raises(ValueError, match='compressed=True'):
        unpach(molecule.pack(compressed=False), compressed=True)
    with pytest.raises(ValueError, match='compressed=False'):
        unpach(molecule.pack(), compressed=False)


def test_a_byte_zero_no_era_claims_is_named_rather_than_guessed_at():
    log = []
    assert unpach(b'\x07\x00\x00\x00', compressed=False, log=log) is None
    assert 'neither a pach version' in log[0]


def test_a_damaged_molecule_record_is_not_retried_as_a_reaction():
    """chython 2 tried the molecule door and fell through to the reaction one on `ValueError`, so a
    molecule record it could not read was reported as an unreadable reaction.  The version byte says
    which door was meant, so the complaint is about the record the caller actually stored."""
    truncated = _molecule().pack(compressed=False)[:6]
    log = []
    assert unpach(truncated, compressed=False, log=log) is None
    assert log and not any('reaction' in x for x in log), log


# -------------------------------------------------------------------- log-or-raise, both shapes

def test_without_a_log_a_recoverable_complaint_is_still_an_error():
    """The answer boundary: a caller who asked for a structure and can only be given a damaged one is
    told, and told everything -- including what the decoder recovered from."""
    with pytest.raises(ValueError, match='damaged'):
        unpach(_damaged_molecule_record(), compressed=False)
    with pytest.raises(ValueError, match='damaged'):
        unpach(_damaged_reaction_record(), compressed=False)


def test_a_log_takes_the_complaints_and_the_structure_comes_back():
    log = []
    molecule = unpach(_damaged_molecule_record(), compressed=False, log=log)
    assert molecule is not None and molecule.atom_count == 4
    assert log and 'flags' in log[0]

    log = []
    rxn = unpach(_damaged_reaction_record(), compressed=False, log=log)
    assert rxn is not None and str(rxn) == str(_reaction())
    assert log and 'flags' in log[0]


def test_an_unreadable_buffer_is_none_with_a_log_and_a_raise_without():
    log = []
    assert unpach(b'', log=log) is None
    assert log
    with pytest.raises(ValueError, match='not a readable pach record'):
        unpach(b'')


def test_the_loop_safe_door_does_not_end_the_loop():
    records = [_molecule().pack(), b'not a record at all', _reaction().pack()]
    log, out = [], []
    for record in records:
        out.append(unpach(record, log=log))
    assert [x is None for x in out] == [False, True, False]
    assert log


# ------------------------------------------------------------- the chython 2 names and signatures

def test_unpack_is_unpach_and_not_a_second_wrapper():
    assert unpack is unpach


def test_data_is_positional_only_as_chython_2_had_it():
    with pytest.raises(TypeError):
        unpach(data=_molecule().pack())


def test_the_method_aliases_are_the_methods():
    molecule, rxn = _molecule(), _reaction()
    assert molecule.pach() == molecule.pack()
    assert rxn.pach() == rxn.pack()
    assert MoleculeContainer.unpach(molecule.pack()).canonical_bytes == molecule.canonical_bytes
    assert str(ReactionContainer.unpach(rxn.pack())) == str(rxn)


@pytest.mark.parametrize('dead', ['check', 'order', 'skip_labels_calculation'])
def test_chython_2s_dropped_keywords_are_refused_rather_than_ignored(dead):
    """A silently accepted `check=False` would promise a refusal was waived and let the encoder raise
    anyway; a silently accepted `order=` would promise an atom order pach has never carried."""
    molecule, rxn = _molecule(), _reaction()
    for call in (lambda **kw: pach(molecule, **kw), lambda **kw: molecule.pach(**kw),
                 lambda **kw: molecule.pack(**kw), lambda **kw: rxn.pach(**kw),
                 lambda **kw: rxn.pack(**kw), lambda **kw: unpach(molecule.pack(), **kw)):
        with pytest.raises(TypeError):
            call(**{dead: True})


def test_a_container_handed_to_the_import_half_names_the_other_door():
    for structure in (_molecule(), _reaction(), read_smarts('[C;a]')):
        with pytest.raises(TypeError, match='pach\\(\\) writes one'):
            unpach(structure)
