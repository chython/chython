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
"""Legacy PDB reader tests.

``[mutant: ...]`` names the implementation line whose alteration makes the assertion fail.  Fixtures
under ``/test/`` are hand-written and column-exact; inline decks come from :func:`_atom`, the format's
column layout written once, and the fixtures are the independent check on it."""
from pathlib import Path

import pytest

from ..pdb import PDBRecord, mmcif, pdb, read_pdb
from ..pdb._records import normalize_element


_DATA = Path(__file__).resolve().parent.parent.parent.parent / 'test'


def _read(name):
    log = []
    return pdb((_DATA / name).read_bytes().decode(), log=log), log


def _read_text(text):
    """The same pair for a deck built inline, so a test asserting on the log needs no local list."""
    log = []
    return pdb(text, log=log), log


def _atom(serial, name, res, chain, seq, x=0., y=0., z=0., element='', charge='', alt=' ', ins=' ',
          occupancy=1., b=10., tag='ATOM'):
    """One ``ATOM``/``HETATM`` line, by column.

    1-6 record, 7-11 serial, 13-16 name, 17 altLoc, 18-20 resName, 22 chainID, 23-26 resSeq,
    27 iCode, 31-38/39-46/47-54 x/y/z, 55-60 occupancy, 61-66 tempFactor, 77-78 element,
    79-80 charge.
    """
    head = f'{tag:<6}{serial:>5} {name:<4}{alt:1}{res:>3} {chain:1}{seq:>4}{ins:1}   '
    assert len(head) == 30
    return f'{head}{x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{b:6.2f}          {element:>2}{charge:>2}'


def _conect(serial, *partners):
    """``CONECT``: serial at 7-11, partners at 12-16, 17-21, 22-26, 27-31."""
    return f'CONECT{serial:>5}' + ''.join(f'{p:>5}' for p in partners)


def _link(name1, res1, chain1, seq1, name2, res2, chain2, seq2):
    """``LINK``: 13-16/43-46 atom names, 18-20/48-50 residues, 22/52 chains, 23-26/53-56 sequences."""
    first = f'{"LINK":<6}      {name1:<4} {res1:>3} {chain1:1}{seq1:>4} '
    assert len(first) == 27
    second = f'{" " * 15}{name2:<4} {res2:>3} {chain2:1}{seq2:>4} '
    assert len(second) == 30
    return first + second


# the element, and only from 77-78


def test_the_element_comes_from_columns_77_78_and_never_from_the_atom_name():
    """A haem pyrrole nitrogen is named ``NA`` and is nitrogen, not sodium.

    Deriving the element from the atom name reads the four haem nitrogens ``NA``, ``NB``, ``NC``,
    ``ND`` as sodium.  [mutant: reading the element from `_field(line, 13, 16)`]
    """
    records = pdb('\n'.join((
        _atom(1, 'NA', 'HEM', 'A', 1, element=' N', tag='HETATM'),
        _atom(2, 'NB', 'HEM', 'A', 1, element=' N', tag='HETATM'),
        _atom(3, 'CA', 'HEM', 'A', 1, element=' C', tag='HETATM'),
        _atom(4, 'FE', 'HEM', 'A', 1, element='FE', tag='HETATM'))))
    assert [a.element for a in records[0].atoms] == ['N', 'N', 'C', 'Fe']
    assert [a.atom_name for a in records[0].atoms] == ['NA', 'NB', 'CA', 'FE']


@pytest.mark.parametrize('column,element', [('FE', 'Fe'), ('ZN', 'Zn'), ('CL', 'Cl'), (' C', 'C'),
                                            ('Fe', 'Fe'), ('SE', 'Se'), (' D', 'H'), (' T', 'H'),
                                            ('fe', 'Fe'), ('zn', 'Zn'), ('d', 'H'),
                                            ('Fe2+', 'Fe'), ('XX', None), ('', None)])
def test_the_element_column_is_case_folded_in_both_directions(column, element):
    """The field is right-justified and upper case, so ``FE`` is iron and not fluorine-einsteinium.

    The fold works both ways: a writer emitting ``fe`` still means iron.  ``D`` and ``T`` are hydrogen
    isotopes, which the format writes as elements.  [mutant: `token.upper()` in `normalize_element`]
    """
    found, isotope, message = normalize_element(column)
    assert found == element
    assert (isotope != 0) is (column.strip().upper() in ('D', 'T'))
    # a field that resolves exactly is not news; one that needed repair or failed outright is
    assert (message is None) is (column.strip().upper() in ('FE', 'ZN', 'CL', 'C', 'SE', ''))


def test_an_atom_with_no_element_column_keeps_everything_else_and_is_counted():
    """Pre-1990s files have no element column, and the atoms are still atoms.

    They keep coordinates, residue and serial, and lose only the fact the file does not state.  The
    count is one log line, not one per atom.  [mutant: the `nothing in the element columns` branch]
    """
    records, log = _read('pdb_damaged.pdb')
    nameless = [a for a in records[0].atoms if a.element is None]
    assert len(nameless) == 3
    assert nameless[0].atom_name == 'C2' and nameless[0].x == 1.5
    assert any(str(m).startswith('atom: 2 line(s) hold an atom with nothing in the element columns 77-78')
               for m in log), log


def test_an_element_column_holding_something_else_names_the_token():
    """``XX`` is not an element; the atom keeps its coordinates and the token is quoted in the log.

    The quoted token is what tells a caller which writer produced the file.  [mutant:
    `_Damage.element`]
    """
    _, log = _read('pdb_damaged.pdb')
    assert any(str(m).startswith('atom: element columns 77-78 hold something that is not an element '
                                 "symbol: 'XX' (1)") for m in log), log


def test_the_charge_column_reads_both_spellings():
    """The format spells it ``2+`` in columns 79-80; writers also spell it ``+2``.

    The zinc in the fixture carries ``2+``.  [mutant: the `_charge` swap]
    """
    records, _ = _read('pdb_ssbond_link.pdb')
    zinc = records[0].atoms[4]
    assert (zinc.element, zinc.charge) == ('Zn', 2)


# bonds the file states


def test_conect_records_are_the_bonds_and_a_reciprocal_pair_is_not_news():
    """``CONECT`` is read, and the reciprocal statement the format asks for is not damage.

    Every well-formed file states each bond twice, once from each atom.  [mutant: the `allowance` of 2
    for conect]
    """
    records, log = _read('pdb_conect.pdb')
    record = records[0]
    assert len(record.bonds) == 2
    assert {(record.atoms[b.a].atom_name, record.atoms[b.b].atom_name) for b in record.bonds} == \
        {('C1', 'C2'), ('C2', 'O')}
    assert not any('restate a pair' in m for m in log), log


def test_a_pair_stated_a_third_time_is_named():
    """Past reciprocity, a repeated pair means a double bond to some writers and a duplicate to others.

    Neither reading is the format's, so the pair is one single bond and the repetition is news.
    [mutant: the `repeated` sum]
    """
    _, log = _read('pdb_damaged.pdb')
    assert any(str(m).startswith('bond: 1 stated bond(s) restate a pair beyond the reciprocal CONECT')
               for m in log), log


def test_no_bond_in_the_file_carries_an_order_and_the_log_says_so():
    """``CONECT``, ``SSBOND`` and ``LINK`` state no bond order, in any writer.

    So every bond this reader builds is a single bond that may not be one, told once per record.
    ``stated_order`` carries the same fact per bond, which is what kekulisation needs.  [mutant: the
    `none of which states an order` line]
    """
    records, log = _read('pdb_conect.pdb')
    assert all(b.order == 1 and not b.stated_order for b in records[0].bonds)
    assert any(str(m).startswith('bond: 2 bond(s) come from CONECT, SSBOND or LINK records, none of which '
                                 'states an order') for m in log), log


def test_a_file_with_no_connectivity_yields_its_atoms_and_names_the_unbonded_count():
    """Zero bonds is an answer; silence about zero bonds is the defect.

    Most legacy entries state ``CONECT`` only for their ligands, so a protein read from one has no
    backbone connectivity at all.  [mutant: the `states no connectivity` branch]
    """
    records, log = _read('pdb_no_conect.pdb')
    assert len(records) == 1 and len(records[0].atoms) == 3 and not records[0].bonds
    assert any('bond: the file states no connectivity for this record; 3 atom(s) are unbonded' in x for x in log)


def test_a_partly_bonded_record_names_the_atoms_no_bond_touches():
    """The water in this fixture is bonded to nothing, and one line says exactly that.

    [mutant: the `touched by no stated bond` branch]
    """
    records, log = _read('pdb_conect.pdb')
    assert records[0].unbonded_count() == 1
    assert any('bond: 1 atom(s) are touched by no stated bond' in x for x in log)


def test_ssbond_and_link_are_bonds():
    """A disulfide and a metal contact, both stated as records and both built.

    ``LINK`` names both atoms, so it is the one legacy record that can state a metal-ligand contact.
    The fixture geometry is deliberately unreal.  [mutant: `_ssbond`, `_link`]
    """
    records, _ = _read('pdb_ssbond_link.pdb')
    record = records[0]
    assert len(record.bonds) == 2
    by_source = {b.source: (record.atoms[b.a].atom_name, record.atoms[b.b].atom_name)
                 for b in record.bonds}
    assert by_source['ssbond'] == ('SG', 'SG')
    assert by_source['link'] == ('SG', 'ZN')
    disulfide = next(b for b in record.bonds if b.source == 'ssbond')
    assert record.atoms[disulfide.a].residue_seq != record.atoms[disulfide.b].residue_seq


def test_a_symmetry_image_partner_gets_no_bond_and_the_operator_is_named():
    """``SSBOND`` under operator 2555 joins an image that is not among the coordinates in the file.

    The message names the offending operator and not the identity one beside it.  [mutant: `_symmetry`
    collecting the non-identity operators]
    """
    records, log = _read('pdb_damaged.pdb')
    assert not any(b.source == 'ssbond' for b in records[0].bonds)
    assert any(str(m).startswith('unsupported: SSBOND on line 12 joins a symmetry image under operator '
                                 '2555') for m in log), log


def test_a_conect_naming_an_absent_serial_is_named():
    """A ``CONECT`` to serial 99 in a file with five atoms loses its bond, and says which serial.

    One line for the whole record, naming the count and the first serial -- which is what a reader of
    the file greps for.  [mutant: the `names atom serial` branch]
    """
    _, log = _read('pdb_damaged.pdb')
    assert any(str(m) == 'bond: 1 CONECT record(s) name an atom serial that is not in this record (99 '
               'first); no bond built for those'
               for m in log), log


def test_alternate_conformers_are_bonded_by_serial_and_not_by_name():
    """``CONECT`` names serials, which identify a conformer exactly, so there is nothing to infer.

    Unlike mmCIF's ``_chem_comp_bond``, which names atoms by name and needs the cross-conformer pairs
    excluded, a ``CONECT`` pair is already unambiguous -- so the reader builds exactly what the file
    states, including a bond from the shared ``CA`` to both conformers of ``CB``.  [mutant: `by_serial`
    keyed on the atom name]
    """
    records, _ = _read('pdb_altloc.pdb')
    record = records[0]
    pairs = {(record.atoms[b.a].atom_name, record.atoms[b.a].alt_loc,
              record.atoms[b.b].atom_name, record.atoms[b.b].alt_loc) for b in record.bonds}
    assert ('CB', 'A', 'OG', 'A') in pairs
    assert ('CB', 'B', 'OG', 'B') in pairs
    assert ('CB', 'A', 'OG', 'B') not in pairs
    assert ('CB', 'B', 'OG', 'A') not in pairs
    assert ('CA', None, 'CB', 'A') in pairs and ('CA', None, 'CB', 'B') in pairs


def test_an_ssbond_naming_an_atom_in_two_conformers_says_which_it_took():
    """``SSBOND`` and ``LINK`` have no altLoc field, so a residue with two conformers is ambiguous.

    The reader takes the first conformer and says so.  [mutant: the `ambiguous` counter]
    """
    log = []
    records = pdb('\n'.join((
        _atom(1, 'CB', 'CYS', 'A', 10, x=1.9, alt='A', element=' C'),
        _atom(2, 'SG', 'CYS', 'A', 10, x=3.7, alt='A', element=' S'),
        _atom(3, 'CB', 'CYS', 'A', 10, x=1.8, alt='B', element=' C'),
        _atom(4, 'SG', 'CYS', 'A', 10, x=3.6, alt='B', element=' S'),
        _atom(5, 'SG', 'CYS', 'A', 20, x=6.0, element=' S'),
        'SSBOND   1 CYS A   10    CYS A   20                          1555   1555  2.03')),
        log=log)
    bond = next(b for b in records[0].bonds if b.source == 'ssbond')
    assert bond.a == 1 and bond.b == 4               # the first conformer, and the unambiguous atom
    assert any(str(m).startswith('bond: 1 SSBOND record(s) name an atom that this record holds in more '
                                 'than one alternate conformer') for m in log), log


# models, records, whole files


def test_each_model_becomes_its_own_record_and_the_bonds_apply_to_all_of_them():
    """``CONECT`` sits after the last ``ENDMDL`` and applies to every model.

    So this reader materialises its records before yielding them, where the mmCIF reader streams: the
    connectivity is not known until the file ends.  [mutant: `_resolve` called from the `MODEL` branch]
    """
    records, log = _read('pdb_two_models.pdb')
    assert [r.model for r in records] == [1, 2]
    assert [len(r.atoms) for r in records] == [2, 2]
    assert [len(r.bonds) for r in records] == [1, 1]
    assert records[0].atoms[0].x != records[1].atoms[0].x
    assert sum(1 for m in log if str(m).startswith('record: the file states 2 models')) == 2


def test_a_file_level_message_reaches_the_caller_once_and_every_record_always():
    """A record has to be readable on its own; the caller's flat log must not repeat a file-level line.

    Both halves pull opposite ways, so both are asserted.  [mutant: `log.extend(file_log)` moved inside
    the record loop]
    """
    records, log = _read('pdb_two_models.pdb')
    assert sum(1 for m in log if str(m).startswith('unsupported: 2 record(s)')) == 1
    for record in records:
        assert any(str(m).startswith('unsupported: 2 record(s)') for m in record.log)
        assert any(str(m).startswith('record: the file states 2 models') for m in record.log)


def test_a_record_type_is_either_read_or_named():
    """``CRYST1``, ``TER`` and ``REMARK`` are in the file and not in the reader, so they are named.

    One aggregated line with a count per type.  ``END`` and ``MASTER`` are the deck's own bookkeeping
    and carry nothing to lose, so they are neither read nor named.  [mutant: the `unread` dict, or
    `_STRUCTURAL`]
    """
    _, log = _read('pdb_conect.pdb')
    assert any(str(m) == 'unsupported: 7 record(s) of 3 type(s) are not modelled: CRYST1 (1), REMARK (5), '
               'TER (1)' for m in log), log
    assert not any('END' in m and 'not modelled' in m for m in log), log


def test_the_entry_id_and_title_are_read():
    """``HEADER`` columns 63-66 and every ``TITLE`` continuation line, joined.

    The entry id is four characters at a fixed offset inside a mostly free-text record.  [mutant:
    `_field(line, 63, 66)`]
    """
    records, _ = _read('pdb_conect.pdb')
    assert records[0].entry_id == 'TST1'
    assert records[0].title == 'ETHANOL AND ONE WATER'


def test_water_is_identifiable_without_an_entity_table():
    """Legacy PDB has no ``_entity``, so water is the two component ids wwPDB issues for it.

    Reading a residue name against a registry of two is reading the file, not perceiving anything.
    [mutant: `_WATER`]
    """
    records, _ = _read('pdb_conect.pdb')
    assert [a.is_water for a in records[0].atoms] == [False, False, False, True]


def test_reading_from_a_path_matches_reading_from_text():
    """`read_pdb` takes a path and yields lazily; `pdb` takes the text.  A `str` is never a path.

    [mutant: `_iter_lines`]
    """
    log = []
    records = list(read_pdb(_DATA / 'pdb_conect.pdb', log=log))
    text_records, text_log = _read('pdb_conect.pdb')
    assert [a.atom_name for a in records[0].atoms] == [a.atom_name for a in text_records[0].atoms]
    assert log == text_log


def test_a_file_with_no_atoms_at_all_still_reports():
    """A deck with a header and nothing else is not an exception; it is a log line.

    [mutant: the `states no ATOM or HETATM records` branch]
    """
    log = []
    records = pdb('HEADER    NOTHING AT ALL\nEND\n', log=log)
    assert records == []
    assert any('record: the file states no ATOM or HETATM records; no atoms read' in x for x in log)


# damage, none of it fatal


def test_a_record_is_not_a_container():
    """The reader stops at a neutral record: three coordinates and every annotation, no chemistry.

    [mutant: any container construction in `read_pdb`]
    """
    records, _ = _read('pdb_conect.pdb')
    assert isinstance(records[0], PDBRecord)
    assert (records[0].atoms[3].x, records[0].atoms[3].y, records[0].atoms[3].z) == (8., 8., 8.)
    assert records[0].atoms[3].residue_name == 'HOH'
    assert records[0].atoms[3].chain == 'W'


def test_the_line_endings_do_not_matter():
    """A deck punched on one platform and read on another arrives CRLF.

    A lone ``\\r`` makes columns 79-80 unreadable and every stripped field end in it.  All three ways in
    are checked: ``str.splitlines()`` and text-mode files drop it, while an iterable of lines that kept
    their endings hands it to the field reader, where the per-field ``.strip()`` protects it.  Applied
    here rather than committed as a fixture because ``core.autocrlf`` is ``input``.  [mutant: the
    `.strip()` in `_field`]
    """
    raw = (_DATA / 'pdb_damaged.pdb').read_text(encoding='utf-8').replace('\r\n', '\n').replace('\n', '\r\n')
    assert raw.count('\r\n') > 10
    for source in (raw, list(raw.split('\n')), raw.splitlines()):
        records = pdb(source) if isinstance(source, str) else list(read_pdb(source))
        assert records[0].atoms[0].element == 'C'
        assert records[0].atoms[0].residue_name == 'LIGX'
        assert records[0].atoms[0].chain == 'A'
        assert records[0].title == 'EVERY MALFORMATION BELOW IS DELIBERATE'


def test_a_four_character_residue_name_is_read_and_the_chain_stays_in_column_22():
    """A writer with a four-character component id uses column 21, which the format leaves blank.

    Reading only columns 18-20 renames the residue silently; chasing the non-blank run past column 21
    eats the chain id.  So the field widens by exactly the one blank column and stops.  [mutant:
    `_residue_name`]
    """
    records, log = _read('pdb_damaged.pdb')
    wide = records[0].atoms[0]
    assert (wide.residue_name, wide.chain, wide.residue_seq) == ('LIGX', 'A', 501)
    assert any(str(m).startswith('atom: 1 line(s) hold a four-character residue name in columns 18-21')
               for m in log), log
    # the three-character names in the same file are unaffected, and so are their chains
    assert [(a.residue_name, a.chain) for a in records[0].atoms[1:]] == [('LIG', 'A')] * 4


def test_an_occupancy_outside_its_range_and_a_b_factor_too_wide_are_stored_as_stated():
    """Both are wrong and both are the file's; the reader stores them and says so.

    A B factor of 1234.5 does not fit six columns at the two decimals the format asks for, so the
    writer wrote something the format cannot express.  [mutant: `range_messages`, `_B_FACTOR_LIMIT`]
    """
    records, log = _read('pdb_damaged.pdb')
    atom = records[0].atoms[3]
    assert atom.occupancy == -0.5 and atom.b_factor == 1234.5
    assert any(str(m).startswith('atom: 1 line(s) hold an occupancy outside 0.0-1.0') for m in log), log
    assert any(str(m).startswith('atom: 1 line(s) hold a B factor over 999.99') for m in log), log


def test_a_truncated_line_keeps_the_fields_it_has():
    """A line ending inside its coordinates keeps x; a line ending before its serial is no atom.

    [mutant: the two `len(line.rstrip())` tests]
    """
    records, log = _read('pdb_damaged.pdb')
    truncated = records[0].atoms[4]
    assert (truncated.x, truncated.y, truncated.z) == (3., None, None)
    assert truncated.atom_name == 'N1' and truncated.residue_seq == 501
    assert len(records[0].atoms) == 5           # the line ending before its serial is not one of them
    assert any(str(m).startswith('atom: 1 line(s) hold an ATOM/HETATM line that ends inside or before its '
                                 'coordinate fields') for m in log), log
    assert any(str(m).startswith('atom: 1 line(s) hold an ATOM/HETATM line that ends before its serial '
                                 'field') for m in log), log


def test_a_repeated_serial_is_named_and_the_first_atom_wins():
    """Two atoms with serial 3 make every ``CONECT`` naming 3 ambiguous.

    The reader resolves to the first and says so; there is nothing in the file to prefer either.
    [mutant: the `duplicated` line]
    """
    records, log = _read('pdb_damaged.pdb')
    assert [a.serial for a in records[0].atoms] == [1, 2, 3, 3, 5]
    assert any(str(m).startswith('atom: 1 atom serial(s) occur more than once (3 first)') for m in log), log


def test_every_malformation_in_the_damaged_fixture_is_read_and_logged():
    """The input posture in one file: nothing is refused, everything is named.

    Asserted as prefixes rather than a count, so a new line does not break it.  [mutant: any single
    damage branch]
    """
    records, log = _read('pdb_damaged.pdb')
    assert len(records) == 1 and len(records[0].atoms) == 5 and len(records[0].bonds) == 1
    for expected in ('unsupported: SSBOND on line 12 joins a symmetry image',
                     'atom: 1 line(s) hold an occupancy outside 0.0-1.0',
                     'atom: 1 line(s) hold a four-character residue name',
                     'atom: 2 line(s) hold an atom with nothing in the element columns 77-78',
                     'atom: 1 line(s) hold a B factor over 999.99',
                     'atom: 1 line(s) hold an ATOM/HETATM line that ends inside or before',
                     'atom: 1 line(s) hold an ATOM/HETATM line that ends before its serial',
                     'atom: element columns 77-78 hold something that is not an element symbol',
                     'unsupported: 9 record(s) of 1 type(s) are not modelled: REMARK (9)',
                     'atom: 1 atom serial(s) occur more than once',
                     'bond: 1 CONECT record(s) name an atom serial that is not in this record',
                     'bond: 1 stated bond(s) restate a pair',
                     'bond: 1 bond(s) come from CONECT, SSBOND or LINK records',
                     'bond: 3 atom(s) are touched by no stated bond'):
        assert any(str(m).startswith(expected) for m in log), (expected, log)


def test_a_damage_counter_names_the_line_it_first_saw():
    """Aggregated damage still has to point somewhere, or a caller cannot find it.

    The counters aggregate, so each carries the line it first saw.  [mutant: `_Damage.first`]
    """
    _, log = _read('pdb_damaged.pdb')
    assert any('(first on line 13)' in m for m in log), log
    assert any('(first on line 18)' in m for m in log), log


# the columns, filled to their edges


def test_every_field_written_to_its_last_column_is_read_from_its_own_columns():
    """A deck whose numeric fields leave no blank column, which is what an off-by-one hides behind.

    Every other fixture writes ``   1.234`` and ``  1.00``, so a field read one column off still strips
    to the same text.  This one leaves no padding: eight-column coordinates with a negative and a
    four-digit integer part, an occupancy filling all six columns, a B factor at the format's limit, a
    four-digit sequence number, a non-blank insertion code, a four-character atom name, a two-character
    *negative* charge (a positive one reads the same with its sign column dropped), a four-digit model
    serial and a title reaching column 80.  [mutant: any column boundary in `_atom`]
    """
    records, log = _read('pdb_wide_columns.pdb')
    assert len(records) == 1
    record = records[0]
    assert record.entry_id == '9XYZ' and record.model == 1234
    assert record.title == ('LEUCINE AND ASPARTATE, EVERY NUMERIC FIELD WRITTEN TO ITS LAST '
                            'COLUMN.')

    first, second = record.atoms
    assert (first.serial, first.atom_name, first.alt_loc, first.residue_name) == (99999, 'HG12', 'A',
                                                                                  'LEU')
    assert (first.chain, first.residue_seq, first.ins_code) == ('B', 9999, 'A')
    assert (first.x, first.y, first.z) == (-999.999, 1234.567, -123.456)
    assert (first.occupancy, first.b_factor) == (1.0, 123.45)
    assert (first.element, first.charge) == ('H', 0)

    assert (second.serial, second.atom_name, second.residue_name) == (99998, 'OD1', 'ASP')
    assert (second.x, second.y, second.z) == (-888.888, -777.777, -666.666)
    assert (second.occupancy, second.b_factor) == (0.1235, 999.99)
    assert (second.element, second.charge) == ('O', -1)
    # nothing above is damage: the B factor is *at* the six-column limit and the occupancy a fraction
    assert [str(x) for x in log] == ['bond: the file states no connectivity for this record; 2 atom(s) are unbonded']


def test_a_link_reads_its_second_atom_name_from_all_four_of_its_columns():
    """``LINK`` names two atoms, and the second field is the one with no other test over it.

    A four-character name is what makes the field's last column significant: with a three-character
    name the same bond is built whether the reader stops at column 45 or 46.  [mutant: the LINK second
    atom-name columns]
    """
    deck = [_atom(1, 'HG12', 'LEU', 'A', 1, element='H'),
            _atom(2, 'HD21', 'ASN', 'A', 2, element='H'),
            _link('HG12', 'LEU', 'A', 1, 'HD21', 'ASN', 'A', 2)]
    records, log = _read_text('\n'.join(deck))
    assert len(records[0].bonds) == 1
    bond = records[0].bonds[0]
    assert (bond.a, bond.b) == (0, 1) and bond.source == 'link'
    assert not any(str(m).startswith('bond: 1 LINK record(s) name a residue or atom') for m in log), log


@pytest.mark.parametrize('column', [32, 80])
def test_the_obsolete_conect_fields_are_named_from_their_first_column_and_their_last(column):
    """``CONECT`` past column 31 held hydrogen bonds and salt bridges until 2011.

    The field spans columns 32-80, so a reader looking one column in from either end misses a writer
    that used only the first or only the last.  [mutant: the obsolete-field columns]
    """
    line = _conect(1, 2).ljust(column - 1) + '3'
    records, log = _read_text('\n'.join([_atom(1, 'C1', 'LIG', 'A', 1, element='C'),
                                         _atom(2, 'C2', 'LIG', 'A', 1, element='C'),
                                         line]))
    assert len(records[0].bonds) == 1
    assert any(str(m).startswith('unsupported: CONECT on line 3 carries the obsolete hydrogen-bond and '
                                 'salt-bridge fields past column 31') for m in log), log


# aggregated damage, and the serial


def test_absent_conect_serials_are_counted_into_one_line_naming_the_first():
    """A file that lost a chain names its every serial from the CONECT records that survived it.

    The count is the report and the first serial is what points into the file.  [mutant: the
    `missing['conect']` counter, the report branch that names it]
    """
    deck = [_atom(1, 'C1', 'LIG', 'A', 1, element='C'),
            _conect(1, 91), _conect(1, 92), _conect(1, 93)]
    records, log = _read_text('\n'.join(deck))
    assert not records[0].bonds
    conect = [m for m in log if 'CONECT record(s) name an atom serial' in m]
    assert [str(x) for x in conect] == ['bond: 3 CONECT record(s) name an atom serial that is not in this record (91 '
                                        'first); no bond built for those'], log


def test_a_duplicate_serial_is_reported_by_a_file_that_states_no_bonds_at_all():
    """The serial feeds bond resolution only, so a file with no ``CONECT`` never used it -- and is
    still a broken file.

    Reporting it only when something depended on it makes the line a property of the reader's work
    rather than of the file.  [mutant: the early return before the duplicate-serial report]
    """
    records, log = _read_text('\n'.join([_atom(1, 'C1', 'LIG', 'A', 1, element='C'),
                                         _atom(1, 'C2', 'LIG', 'A', 1, element='C')]))
    assert [a.serial for a in records[0].atoms] == [1, 1]
    assert any(str(m) == 'atom: 1 atom serial(s) occur more than once (1 first); a CONECT naming one is '
               'resolved to the atom that came first' for m in log), log


# a bond is what the file said


def test_a_distance_decides_nothing_in_either_direction():
    """Two atoms half an Angstrom apart with no ``CONECT``, and two 500 Angstroms apart with one.

    Every other fixture states the bonds its geometry implies, so only this deck fails the moment a
    coordinate reaches a bond decision -- as a cutoff that adds a bond or a check that drops one.
    There is no covalent radius and no distance function in this package.  [mutant: any distance test
    inside the bond builders]
    """
    records, log = _read_text('\n'.join([
        _atom(1, 'C1', 'LIG', 'A', 1, x=0., y=0., z=0., element='C'),
        _atom(2, 'C2', 'LIG', 'A', 1, x=0.5, y=0., z=0., element='C'),
        _atom(3, 'C3', 'LIG', 'A', 1, x=500., y=0., z=0., element='C'),
        _conect(1, 3), _conect(3, 1)]))
    bonds = {bond.key for bond in records[0].bonds}
    assert bonds == {(0, 2)}                      # the stated one, at 500 A; not the touching pair
    assert any(str(m).startswith('bond: 1 atom(s) are touched by no stated bond') for m in log), log


# cross-reader occupancy invariant


def test_five_atoms_with_bad_occupancy_produce_one_aggregate_line():
    """Five atoms carrying out-of-range occupancy count as one writer defect, not five findings.

    ``range_messages`` keys each pair on a kind string, and ``_atom`` must feed that kind to
    ``damage.hit`` so the counter reaches five and ``_report`` emits one line.  [mutant: the
    `damage.hit` call in the occupancy loop inside `_atom`]
    """
    deck = '\n'.join(_atom(i, 'C1', 'LIG', 'A', 1, element='C', occupancy=2.0) for i in range(1, 6))
    _, log = _read_text(deck)
    occ_lines = [m for m in log if 'occupancy' in m and 'outside 0.0-1.0' in m]
    assert len(occ_lines) == 1, log
    assert str(occ_lines[0]).startswith('atom: 5 line(s) hold an occupancy outside 0.0-1.0'), log


def test_atoms_with_valid_occupancy_log_nothing_about_occupancy():
    """In-range occupancy on every atom must not produce any occupancy log entry.

    The armed negative: without it, nothing stops ``damage.hit`` becoming unconditional.  [mutant:
    removing `not 0.0 <= occupancy <= 1.0` in `range_messages`]
    """
    deck = '\n'.join(_atom(i, 'C1', 'LIG', 'A', 1, element='C', occupancy=0.8) for i in range(1, 3))
    _, log = _read_text(deck)
    assert not any('occupancy' in m and 'outside' in m for m in log), log


def test_legacy_and_mmcif_each_produce_one_occupancy_line_for_the_same_bad_value():
    """Both readers share ``range_messages`` and must both aggregate by its kind string.

    A caller counting findings otherwise gets a different answer per format: five legacy atoms at
    occupancy 2.0 giving five lines against mmCIF's one.  [mutant: removing `damage.hit` in the
    occupancy loop inside `_legacy._atom`]
    """
    # legacy: five atoms, each with occupancy 2.0
    legacy_deck = '\n'.join(_atom(i, 'C1', 'LIG', 'A', 1, element='C', occupancy=2.0)
                            for i in range(1, 6))
    legacy_log = []
    pdb(legacy_deck, log=legacy_log)

    # mmCIF: five rows, same occupancy, minimal loop
    mmcif_block = (
        'data_PAIR\n'
        'loop_\n'
        '_atom_site.id\n'
        '_atom_site.type_symbol\n'
        '_atom_site.occupancy\n'
        + '\n'.join(f'{i} C 2.0' for i in range(1, 6))
        + '\n'
    )
    mmcif_log = []
    mmcif(mmcif_block, log=mmcif_log)

    legacy_occ = [m for m in legacy_log if 'outside 0.0-1.0' in m]
    mmcif_occ = [m for m in mmcif_log if 'outside 0.0-1.0' in m]
    assert len(legacy_occ) == 1, legacy_log
    assert len(mmcif_occ) == 1, mmcif_log


# parse-failure aggregation


def test_five_atoms_with_bad_x_coordinate_produce_one_aggregate_line():
    """Five atoms with an unparseable x field count as one writer defect, not five findings.

    ``_number`` routes through ``damage.hit`` when called from ``_atom``, so the counter reaches five
    and ``_report`` emits one line.  The bad x field is spliced into columns 31-38 (0-indexed 30-37),
    where ``_atom`` writes ``x:8.3f`` and ``_legacy._number`` reads it.  [mutant: reverting `_number`
    to `sink.append` unconditionally in `_legacy.py`]
    """
    good = _atom(1, 'C', 'LIG', 'A', 1, element='C')
    bad_x = good[:30] + '   BADC ' + good[38:]
    deck = '\n'.join(bad_x for _ in range(5))
    _, log = _read_text(deck)
    x_lines = [m for m in log if 'x coordinate' in m and 'not a number' in m]
    assert len(x_lines) == 1, log
    assert str(x_lines[0]).startswith('atom: 5 line(s) hold an x coordinate that is not a number'), log


def test_atoms_with_valid_numeric_fields_log_nothing_about_parse_failures():
    """Well-formed coordinates, sequence numbers and charges produce no parse-failure log entry.

    The armed negative for the three ``damage.hit`` calls in ``_number``, ``_integer`` and ``_charge``.
    [mutant: removing the ``try/except`` guard in ``_number`` in `_legacy.py`]
    """
    deck = '\n'.join(_atom(i, 'C', 'LIG', 'A', 1, element='C') for i in range(1, 4))
    _, log = _read_text(deck)
    assert not any('not a number' in m or 'not an integer' in m or 'not a charge' in m
                   for m in log), log


def test_two_different_broken_fields_produce_two_distinct_log_lines():
    """A bad x coordinate and a bad sequence number are different defects and get different counters.

    Distinct ``what`` strings in ``_number`` and ``_integer`` give distinct kind keys in ``_Damage``,
    so each gets its own aggregate line.  The bad x field is in columns 31-38 (0-indexed 30-37); the
    bad sequence number in columns 23-26 (0-indexed 22-25).  [mutant: making ``_number`` and
    ``_integer`` share the same ``what`` key in ``damage.hit``]
    """
    good = _atom(1, 'C', 'LIG', 'A', 1, element='C')
    bad_x = good[:30] + '   BADC ' + good[38:]
    bad_seq = good[:22] + 'ABCD' + good[26:]
    _, log = _read_text(bad_x + '\n' + bad_seq)
    parse_lines = [m for m in log if 'not a number' in m or 'not an integer' in m]
    assert len(parse_lines) == 2, log
    assert any('x coordinate' in m for m in parse_lines), log
    assert any('sequence number' in m for m in parse_lines), log
