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
"""MOL2 reader tests.

Fixtures under /test/: mol2_simple.mol2 (benzene, C.ar types), mol2_two_records.mol2,
mol2_no_charge.mol2 (no charge column), mol2_bad_type.mol2 (unknown SYBYL types).
"""

from pathlib import Path

import pytest
from pytest import raises

from ..mol2 import Mol2ParseError, mol2_mol, read_mol2


_TEST_DIR = Path(__file__).resolve().parent.parent.parent.parent / 'test'


# mol2_mol helpers

def _mol(text, log=None):
    """Parse text as one MOL2 record, returning the molecule."""
    if log is None:
        log = []
    return mol2_mol(text, log=log)


def _record(*lines):
    """Assemble a minimal MOL2 record string from content lines."""
    return '\n'.join(['@<TRIPOS>MOLECULE'] + list(lines))


# basic reads

def test_benzene_aromatic_bonds_stored_as_order_4():
    """Aromatic bonds (MOL2 type 'ar') are stored as bond order 4, not kekulized.

    [mutant: change `'ar': 4` in _BOND_ORDER to `'ar': None`]
    """
    mol = _mol((_TEST_DIR / 'mol2_simple.mol2').read_text(encoding='utf-8'))
    assert mol.atom_count == 6
    assert mol.bond_count == 6
    assert mol.aromatic_bond_count == 6


def test_benzene_title_preserved():
    """The MOLECULE name line becomes the molecule's title."""
    mol = _mol((_TEST_DIR / 'mol2_simple.mol2').read_text(encoding='utf-8'))
    assert mol.title == 'benzene'


def test_benzene_coordinates_set():
    """3D coordinates are stored; xy_of returns non-zero values."""
    mol = _mol((_TEST_DIR / 'mol2_simple.mol2').read_text(encoding='utf-8'))
    xs = [mol.xy_of(sid)[0] for sid in mol.atom_numbers]
    assert any(abs(x) > 0.01 for x in xs)


def test_benzene_no_log():
    """A well-formed record with recognised types produces an empty log."""
    log = []
    _mol((_TEST_DIR / 'mol2_simple.mol2').read_text(encoding='utf-8'), log=log)
    assert not log, f'unexpected log entries: {log}'


def test_sp3_carbon_element():
    """A C.3 atom becomes a carbon atom."""
    text = _record(
        'methane', ' 1 0 0 0 0', 'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        '      1 C1  0.0 0.0 0.0 C.3  1 LIG 0.000',
        '@<TRIPOS>BOND',
    )
    mol = _mol(text)
    assert mol.atom_count == 1
    assert mol.element_of(list(mol.atom_numbers)[0]) == 6   # atomic number for C


def test_sp2_nitrogen_element():
    """A N.2 atom becomes a nitrogen atom."""
    text = _record(
        'test', ' 1 0 0 0 0', 'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        '      1 N1  0.0 0.0 0.0 N.2  1 LIG 0.000',
        '@<TRIPOS>BOND',
    )
    mol = _mol(text)
    assert mol.element_of(list(mol.atom_numbers)[0]) == 7   # N


def test_aromatic_nitrogen():
    """N.ar becomes a nitrogen; the aromatic flag is on the bond, not the atom type."""
    text = _record(
        'test', ' 1 0 0 0 0', 'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        '      1 N1  0.0 0.0 0.0 N.ar  1 LIG 0.000',
        '@<TRIPOS>BOND',
    )
    mol = _mol(text)
    assert mol.element_of(list(mol.atom_numbers)[0]) == 7


# multi-record

def test_two_record_file_yields_two_molecules():
    """Multi-record files yield one molecule per @<TRIPOS>MOLECULE block.

    [mutant: stop after the first record]
    """
    path = _TEST_DIR / 'mol2_two_records.mol2'
    results = list(read_mol2(path))
    assert len(results) == 2


def test_two_record_file_titles():
    """Each record in a multi-record file gets its own title."""
    path = _TEST_DIR / 'mol2_two_records.mol2'
    results = list(read_mol2(path))
    titles = [mol.title for mol, _ in results]
    assert titles == ['methane', 'formaldehyde']


def test_two_record_file_atom_counts():
    """Multi-record file: each record's atom count is independent."""
    path = _TEST_DIR / 'mol2_two_records.mol2'
    results = list(read_mol2(path))
    counts = [mol.atom_count for mol, _ in results]
    assert counts == [1, 3]


def test_read_mol2_yields_log_per_record():
    """read_mol2 yields (molecule, log) pairs; each log is a fresh list."""
    path = _TEST_DIR / 'mol2_two_records.mol2'
    results = list(read_mol2(path))
    assert len(results) == 2
    mols = [m for m, _ in results]
    logs = [rec for _, rec in results]
    assert all(isinstance(rec, list) for rec in logs)
    assert logs[0] is not logs[1]


def test_read_mol2_from_string():
    """read_mol2 accepts a string as well as a path."""
    text = (_TEST_DIR / 'mol2_two_records.mol2').read_text(encoding='utf-8')
    results = list(read_mol2(text))
    assert len(results) == 2


def test_read_mol2_from_stream():
    """read_mol2 accepts an open text file handle."""
    path = _TEST_DIR / 'mol2_two_records.mol2'
    with open(path, encoding='utf-8') as fh:
        results = list(read_mol2(fh))
    assert len(results) == 2


# charge column absent

def test_no_charge_column_atoms_stored():
    """A record without the charge column stores atoms with charge 0.

    The charge column is optional in MOL2; its absence is not an error.

    [mutant: raise on a missing 9th field]
    """
    path = _TEST_DIR / 'mol2_no_charge.mol2'
    log = []
    results = list(read_mol2(path, log_factory=list))
    mol, log = results[0]
    assert mol.atom_count == 9
    charge_lines = [rec for rec in log if 'charge' in str(rec).lower()]
    assert not charge_lines, f'unexpected charge log: {charge_lines}'


def test_no_charge_column_produces_no_charge_log():
    """An absent charge column is normal and produces no charge-related log entries."""
    path = _TEST_DIR / 'mol2_no_charge.mol2'
    results = list(read_mol2(path, log_factory=list))
    _, log = results[0]
    charge_lines = [rec for rec in log if 'charge' in str(rec).lower()]
    assert not charge_lines, f'unexpected charge log: {charge_lines}'


# malformed atom types

def test_bare_element_accepted():
    """A SYBYL type that is just an element symbol (no dot) is read correctly.

    Bare types like 'Br', 'Cl', 'Na' appear when the writer does not know the hybridization.
    """
    path = _TEST_DIR / 'mol2_bad_type.mol2'
    results = list(read_mol2(path, log_factory=list))
    mol, _ = results[0]
    # 'C' and 'Br' should both be stored; 'N.bogus' logs but N is stored too; 'O.3' is clean
    assert mol.atom_count == 4


def test_bare_element_no_error_log():
    """A bare element type is recognised; no error log line is emitted for it."""
    path = _TEST_DIR / 'mol2_bad_type.mol2'
    results = list(read_mol2(path, log_factory=list))
    _, log = results[0]
    # 'C' and 'Br' should not produce log lines; only 'N.bogus' should
    c_lines = [rec for rec in log if "'C'" in rec and 'not recognised' in rec]
    br_lines = [rec for rec in log if "'Br'" in rec and 'not recognised' in rec]
    assert not c_lines, c_lines
    assert not br_lines, br_lines


def test_unknown_sybyl_tag_is_unsupported_not_a_broken_file():
    """An unrecognised type tag (N.bogus) is our table's limit, so the line carries the
    `unsupported:` prefix rather than blaming the file for a type Tripos may well document.

    [mutant: `_resolve_type` returns the prefix element with no note]
    """
    path = _TEST_DIR / 'mol2_bad_type.mol2'
    results = list(read_mol2(path, log_factory=list))
    _, log = results[0]
    tag_lines = [rec for rec in log if 'N.bogus' in rec]
    assert tag_lines, f'expected a log line naming the unknown tag, got: {log}'
    for line in tag_lines:
        assert str(line).startswith('unsupported'), line
    assert any('bogus' in rec and 'not interpreted' in rec for rec in tag_lines), tag_lines


# partial charges

def test_partial_charges_not_stored():
    """When charge_type is GASTEIGER the float charge column is not used as formal charge, and
    an 'unsupported:' line says why.

    [mutant: drop GASTEIGER from `_PARTIAL_CHARGE_TYPES`]
    """
    text = _record(
        'test', ' 1 0 0 0 0', 'SMALL', 'GASTEIGER', '',
        '@<TRIPOS>ATOM',
        '      1 N1  0.0 0.0 0.0 N.3  1 LIG  0.345',
        '@<TRIPOS>BOND',
    )
    log = []
    mol = _mol(text, log=log)
    sids = list(mol.atom_numbers)
    assert mol.charge_of(sids[0]) == 0, 'partial charge must not be used as formal charge'
    assert any(str(rec).startswith('unsupported') for rec in log), \
        f'expected unsupported log for GASTEIGER, got: {log}'


def test_no_charges_formal_charge_rounded():
    """Integer formal charges (stored as floats) are rounded to the nearest integer.

    A nitrogen with charge=-1.0 in a NO_CHARGES file is a formal negative charge.
    """
    text = _record(
        'test', ' 1 0 0 0 0', 'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        '      1 N1  0.0 0.0 0.0 N.4  1 LIG  -1.0',
        '@<TRIPOS>BOND',
    )
    log = []
    mol = _mol(text, log=log)
    sid = list(mol.atom_numbers)[0]
    assert mol.charge_of(sid) == -1, f'charge {mol.charge_of(sid)} != -1'


# bond types

def test_amide_bond_stored_as_single_with_unsupported_log():
    """Bond type 'am' (amide) is not a distinct order in this library; stored as single with an
    'unsupported:' log line.

    [mutant: drop the 'unsupported' prefix from the amide log line]
    """
    text = _record(
        'test', ' 2 1 0 0 0', 'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        '      1 C1  0.0 0.0 0.0 C.2  1 LIG 0.0',
        '      2 N1  1.3 0.0 0.0 N.am 1 LIG 0.0',
        '@<TRIPOS>BOND',
        '      1 1 2 am',
    )
    log = []
    mol = _mol(text, log=log)
    assert mol.bond_count == 1
    sids = list(mol.atom_numbers)
    assert mol.order_of(sids[0], sids[1]) == 1
    assert any(str(rec).startswith('unsupported') for rec in log), f'expected unsupported log, got: {log}'


def test_dummy_bond_skipped_with_unsupported_log():
    """Bond type 'du' (dummy) is not modelled; the bond is skipped.

    [mutant: give 'du' order 1 in _BOND_ORDER]
    """
    text = _record(
        'test', ' 2 1 0 0 0', 'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        '      1 C1  0.0 0.0 0.0 C.3  1 LIG 0.0',
        '      2 C2  1.5 0.0 0.0 C.3  1 LIG 0.0',
        '@<TRIPOS>BOND',
        '      1 1 2 du',
    )
    log = []
    mol = _mol(text, log=log)
    assert mol.bond_count == 0
    assert any(str(rec).startswith('unsupported') for rec in log), \
        f'expected unsupported log for du, got: {log}'


def test_nc_bond_skipped_with_unsupported_log():
    """Bond type 'nc' (not connected) is not modelled; the bond is skipped."""
    text = _record(
        'test', ' 2 1 0 0 0', 'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        '      1 C1  0.0 0.0 0.0 C.3  1 LIG 0.0',
        '      2 C2  1.5 0.0 0.0 C.3  1 LIG 0.0',
        '@<TRIPOS>BOND',
        '      1 1 2 nc',
    )
    log = []
    mol = _mol(text, log=log)
    assert mol.bond_count == 0
    assert any(str(rec).startswith('unsupported') for rec in log), \
        f'expected unsupported log for nc, got: {log}'


def test_unknown_bond_type_logged_and_stored_as_single():
    """A completely unknown bond type string is logged and the bond is stored as single."""
    text = _record(
        'test', ' 2 1 0 0 0', 'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        '      1 C1  0.0 0.0 0.0 C.3  1 LIG 0.0',
        '      2 C2  1.5 0.0 0.0 C.3  1 LIG 0.0',
        '@<TRIPOS>BOND',
        '      1 1 2 xz',
    )
    log = []
    mol = _mol(text, log=log)
    assert mol.bond_count == 1
    sids = list(mol.atom_numbers)
    assert mol.order_of(sids[0], sids[1]) == 1
    assert log, 'expected at least one log entry for unknown bond type'


# malformed ATOM block

def test_atom_block_too_few_fields_raises():
    """An ATOM line with fewer than 6 fields raises Mol2ParseError: the atom's identity
    (element, coordinates) cannot be read at all, unlike a chemically-broken atom which is
    stored-and-logged.
    """
    text = _record(
        'test', ' 1 0 0 0 0', 'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        '      1 C1  0.0 0.0',     # only 4 fields
        '@<TRIPOS>BOND',
    )
    with raises(Mol2ParseError):
        _mol(text)


def test_malformed_coordinates_stored_as_zero():
    """Non-numeric coordinate fields are replaced by 0.0 and logged.

    [mutant: raise instead of logging on the coordinate error path]
    """
    text = _record(
        'test', ' 1 0 0 0 0', 'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        '      1 C1  0.0 NOTANUMBER 0.0 C.3  1 LIG 0.0',
        '@<TRIPOS>BOND',
    )
    log = []
    mol = _mol(text, log=log)
    assert mol.atom_count == 1
    assert any('coordinate' in rec for rec in log), f'expected coordinate log, got: {log}'


def test_atom_count_mismatch_logged():
    """A discrepancy between the MOLECULE count and the actual ATOM block is logged; the reader
    uses the block count, not the header count.

    [mutant: raise instead of logging]
    """
    text = _record(
        'test', ' 5 0 0 0 0', 'SMALL', 'NO_CHARGES', '',  # claims 5 atoms
        '@<TRIPOS>ATOM',
        '      1 C1  0.0 0.0 0.0 C.3  1 LIG 0.0',          # only 1 present
        '@<TRIPOS>BOND',
    )
    log = []
    mol = _mol(text, log=log)
    assert mol.atom_count == 1
    assert any('claims' in rec or 'ATOM block' in rec for rec in log), \
        f'expected count mismatch log, got: {log}'


def test_bond_referencing_missing_atom_is_logged():
    """A bond referencing an atom ID not present in the ATOM block is logged and skipped.

    [mutant: skip silently in the id_to_index check, without logging]
    """
    text = _record(
        'test', ' 1 1 0 0 0', 'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        '      1 C1  0.0 0.0 0.0 C.3  1 LIG 0.0',
        '@<TRIPOS>BOND',
        '      1 1 99 1',   # atom 99 does not exist
    )
    log = []
    mol = _mol(text, log=log)
    assert mol.bond_count == 0
    assert any('not in the ATOM block' in rec or 'atom id 99' in rec for rec in log), \
        f'expected missing-atom log, got: {log}'


def test_self_loop_bond_logged():
    """A self-loop bond (same atom both ends) is logged and dropped."""
    text = _record(
        'test', ' 1 1 0 0 0', 'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        '      1 C1  0.0 0.0 0.0 C.3  1 LIG 0.0',
        '@<TRIPOS>BOND',
        '      1 1 1 1',    # self-loop
    )
    log = []
    mol = _mol(text, log=log)
    assert mol.bond_count == 0
    assert any('self-loop' in rec for rec in log), f'expected self-loop log, got: {log}'


def test_duplicate_bond_logged():
    """A duplicated bond between the same atoms is logged and the second occurrence dropped."""
    text = _record(
        'test', ' 2 2 0 0 0', 'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        '      1 C1  0.0 0.0 0.0 C.3  1 LIG 0.0',
        '      2 C2  1.5 0.0 0.0 C.3  1 LIG 0.0',
        '@<TRIPOS>BOND',
        '      1 1 2 1',
        '      2 1 2 1',   # duplicate
    )
    log = []
    mol = _mol(text, log=log)
    assert mol.bond_count == 1
    assert any('duplicate' in rec for rec in log), f'expected duplicate log, got: {log}'


# pseudo-atoms

def test_lone_pair_skipped_with_unsupported_log():
    """An LP (lone pair) pseudo-atom has no nucleus and must not become a graph atom.

    LP is a legitimate MOL2 construct this library does not model, so the drop carries the
    `unsupported:` prefix rather than accusing the file.

    [mutant: let LP through, or drop the 'unsupported' prefix]
    """
    text = _record(
        'test', ' 2 0 0 0 0', 'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        '      1 C1  0.0 0.0 0.0 C.3  1 LIG 0.0',
        '      2 LP1 0.5 0.5 0.0 LP   1 LIG 0.0',
        '@<TRIPOS>BOND',
    )
    log = []
    mol = _mol(text, log=log)
    assert mol.atom_count == 1   # LP was not stored
    assert any(str(rec).startswith('unsupported') and 'LP' in rec for rec in log), \
        f'expected unsupported log for LP, got: {log}'


def test_dummy_atom_skipped_with_unsupported_log():
    """A Du (dummy atom) is a legitimate MOL2 pseudo-atom; its drop is `unsupported:` prefixed."""
    text = _record(
        'test', ' 2 0 0 0 0', 'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        '      1 C1  0.0 0.0 0.0 C.3  1 LIG 0.0',
        '      2 DU1 1.0 0.0 0.0 Du   1 LIG 0.0',
        '@<TRIPOS>BOND',
    )
    log = []
    mol = _mol(text, log=log)
    assert mol.atom_count == 1
    assert any(str(rec).startswith('unsupported') and 'Du' in rec for rec in log), \
        f'expected unsupported log for Du, got: {log}'


def test_pseudo_atom_does_not_cause_false_count_mismatch():
    """A record with pseudo-atoms in the ATOM block must not trigger the count-mismatch log: the
    header claims 2 atoms and the block has 2 lines, and not storing the LP is our limitation.

    [mutant: compare the header against `len(atoms)` instead of `total_atom_lines`]
    """
    text = _record(
        'test', ' 2 0 0 0 0', 'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        '      1 C1  0.0 0.0 0.0 C.3  1 LIG 0.0',
        '      2 LP1 0.5 0.5 0.0 LP   1 LIG 0.0',
        '@<TRIPOS>BOND',
    )
    log = []
    _mol(text, log=log)
    count_lines = [rec for rec in log if 'header claims' in rec and 'ATOM' in rec]
    assert not count_lines, \
        f'false count-mismatch log when header was correct: {count_lines}'


def test_genuine_count_mismatch_still_logged():
    """A header claiming more atoms than the block contains -- a truncated record -- is reported.

    [mutant: make `total_atom_lines` always equal `num_atoms`]
    """
    text = _record(
        'test', ' 5 0 0 0 0', 'SMALL', 'NO_CHARGES', '',   # claims 5
        '@<TRIPOS>ATOM',
        '      1 C1  0.0 0.0 0.0 C.3  1 LIG 0.0',           # only 1 line
        '@<TRIPOS>BOND',
    )
    log = []
    mol = _mol(text, log=log)
    assert mol.atom_count == 1
    assert any('claims' in rec and 'ATOM' in rec for rec in log), \
        f'expected count mismatch log for genuinely short block, got: {log}'


def test_bond_to_pseudo_atom_is_unsupported_not_dangling():
    """A bond from a real atom to an LP/Du endpoint is `unsupported:`, not a dangling reference:
    the atom id WAS in the ATOM block and we chose not to store it.

    [mutant: use the plain `bond` prefix for pseudo-atom endpoints in _parse_bonds]
    """
    text = _record(
        'test', ' 2 1 0 0 0', 'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        '      1 O1  0.0 0.0 0.0 O.3  1 LIG 0.0',
        '      2 LP1 0.5 0.0 0.0 LP   1 LIG 0.0',
        '@<TRIPOS>BOND',
        '      1 1 2 1',
    )
    log = []
    mol = _mol(text, log=log)
    assert mol.atom_count == 1
    assert mol.bond_count == 0
    pseudo_bond_lines = [rec for rec in log if 'bond' in rec and '2' in rec]
    assert pseudo_bond_lines, f'no bond log at all, got: {log}'
    for line in pseudo_bond_lines:
        assert str(line).startswith('unsupported'), \
            f'bond to pseudo-atom must be unsupported, not: {line!r}'


# SUBSTRUCTURE block

def test_substructure_block_present_record_reads_cleanly():
    """A SUBSTRUCTURE block is ignored gracefully; the record is not rejected.

    [mutant: raise on SUBSTRUCTURE in the section splitter]
    """
    text = _record(
        'test', ' 1 0 1 0 0', 'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        '      1 C1  0.0 0.0 0.0 C.3  1 LIG 0.0',
        '@<TRIPOS>BOND',
        '@<TRIPOS>SUBSTRUCTURE',
        '     1 LIG       1 TEMP              0 ****  ****    0 ROOT',
    )
    log = []
    mol = _mol(text, log=log)
    assert mol.atom_count == 1


def test_atoms_spanning_two_substructures_fires_one_log_line():
    """When atoms belong to more than one (subst_id, subst_name) pair, exactly one 'unsupported:'
    line fires, naming the number of substructures.  The two columns appear in practically every
    record real tools write, so only a record that genuinely spans them is worth a line; the exact
    text is asserted because the count and the column names are what make it actionable.

    [mutant: remove the `if len(substructures) > 1:` guard]
    """
    text = _record(
        'test', ' 2 0 0 0 0', 'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        '      1 C1  0.0 0.0 0.0 C.3  1 LIG  0.0',
        '      2 N1  1.5 0.0 0.0 N.3  2 RES  0.0',
        '@<TRIPOS>BOND',
    )
    log = []
    mol = _mol(text, log=log)
    assert mol.atom_count == 2
    subst_lines = [rec for rec in log if 'substructures' in rec]
    assert len(subst_lines) == 1
    assert str(subst_lines[0]) == (
        'unsupported: the ATOM block assigns its atoms to 2 substructures '
        '(the subst_id and subst_name columns); chython stores no residue annotation'
    )


def test_status_bit_column_fires_one_log_line():
    """A tenth field on any ATOM line is the optional status_bit column, whose values chython
    does not model; exactly one 'unsupported:' line fires.  The exact text is asserted because
    'status_bit' is the term a caller would search for.

    [mutant: remove the `if status_bits:` guard]
    """
    text = _record(
        'test', ' 1 0 0 0 0', 'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        '      1 C1  0.0 0.0 0.0 C.3  1 LIG  0.0  ****',
        '@<TRIPOS>BOND',
    )
    log = []
    mol = _mol(text, log=log)
    assert mol.atom_count == 1
    status_lines = [rec for rec in log if 'status_bit' in rec]
    assert len(status_lines) == 1
    assert str(status_lines[0]) == 'unsupported: the ATOM block states status_bit values; chython stores none'


def test_single_substructure_no_status_bits_logs_nothing_about_them():
    """A record whose atoms all share one (subst_id, subst_name) and carry no status_bit column
    logs nothing about either: an unconditional line carries no information and destroys
    ``any(rec.startswith('unsupported'))`` as a screen for callers.

    [mutant: change `> 1` to `>= 1`]
    """
    text = _record(
        'test', ' 2 1 0 0 0', 'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        '      1 C1  0.0 0.0 0.0 C.3  1 LIG  0.0',
        '      2 N1  1.5 0.0 0.0 N.3  1 LIG  0.0',
        '@<TRIPOS>BOND',
        '      1 1 2 1',
    )
    log = []
    _mol(text, log=log)
    noisy = [rec for rec in log if 'substructure' in rec or 'status_bit' in rec]
    assert not noisy, f'ordinary record produced substructure/status_bit log: {noisy}'


# Windows line endings

def test_windows_crlf_accepted():
    """Records with \\r\\n line endings are parsed identically to Unix endings.

    [mutant: `rstrip('\\n')` instead of `rstrip('\\r\\n')`, leaving \\r in field values]
    """
    lines = [
        '@<TRIPOS>MOLECULE\r\n',
        'test\r\n',
        ' 1 0 0 0 0\r\n',
        'SMALL\r\n',
        'NO_CHARGES\r\n',
        '\r\n',
        '@<TRIPOS>ATOM\r\n',
        '      1 C1  0.0 0.0 0.0 C.3  1 LIG 0.0\r\n',
        '@<TRIPOS>BOND\r\n',
    ]
    from io import StringIO
    results = list(read_mol2(StringIO(''.join(lines))))
    assert len(results) == 1
    mol, _ = results[0]
    assert mol.atom_count == 1


# the third coordinate

def test_z_coordinates_are_stored_and_no_longer_logged_as_a_loss():
    """MOL2 is a 3D format: every record states a z, and the stated z is stored rather than read
    and dropped with an `unsupported:` line.  Both segments are filled -- `has_coordinates` for
    the depiction, `has_3d` for the geometry.
    """
    text = _record(
        'test', ' 1 0 0 0 0', 'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        '      1 C1  1.0 2.0 3.5 C.3  1 LIG 0.0',
        '@<TRIPOS>BOND',
    )
    log = []
    m = _mol(text, log=log)
    assert m.has_3d is True
    assert m.xyz_of(1) == (1.0, 2.0, 3.5)
    assert m.has_coordinates is True
    assert m.xy_of(1) == (1.0, 2.0)
    assert not any('z' in rec and rec.startswith('unsupported') for rec in log), \
        f'the z loss is no longer a loss, so it must not be logged as one: {log}'


def test_a_flat_mol2_record_stores_no_geometry():
    """z == 0 for every atom is a record with no geometry, and a conformer would claim one.

    `solid` is a whole-record test rather than a per-atom one: `has_3d` must answer False here, or
    a caller cannot tell a placed structure from a flat drawing that arrived in a 3D format.
    """
    text = _record(
        'test', ' 1 0 0 0 0', 'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        '      1 C1  1.0 2.0 0.0 C.3  1 LIG 0.0',
        '@<TRIPOS>BOND',
    )
    m = _mol(text)
    assert m.has_coordinates is True
    assert m.has_3d is False
    assert m.xyz_of(1) is None


# mol2_mol convenience

def test_mol2_mol_accepts_string_with_molecule_tag():
    """mol2_mol accepts a string that includes the @<TRIPOS>MOLECULE tag."""
    text = (
        '@<TRIPOS>MOLECULE\n'
        'water\n'
        ' 1 0 0 0 0\n'
        'SMALL\nNO_CHARGES\n\n'
        '@<TRIPOS>ATOM\n'
        '      1 O1  0.0 0.0 0.0 O.3  1 LIG 0.0\n'
        '@<TRIPOS>BOND\n'
    )
    mol = mol2_mol(text)
    assert mol.atom_count == 1
    assert mol.element_of(list(mol.atom_numbers)[0]) == 8   # O


def test_mol2_mol_accepts_list_of_lines():
    """mol2_mol accepts a list of line strings as well as a plain string."""
    lines = [
        '@<TRIPOS>MOLECULE',
        'water',
        ' 1 0 0 0 0',
        'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        '      1 O1  0.0 0.0 0.0 O.3  1 LIG 0.0',
        '@<TRIPOS>BOND',
    ]
    mol = mol2_mol(lines)
    assert mol.atom_count == 1


# sybyl_types table

def test_sybyl_types_table_loads():
    """The table loads without error and contains the expected carbon types."""
    from chython.chemistry._tables import sybyl_types
    table = sybyl_types()
    assert 'C.3' in table
    assert 'C.ar' in table
    assert 'N.4' in table
    assert 'O.co2' in table
    assert 'S.o2' in table


def test_sybyl_types_pseudo_atoms_have_empty_element():
    """Pseudo-atom types (LP, Du, etc.) have an empty element string in the table."""
    from chython.chemistry._tables import sybyl_types
    table = sybyl_types()
    for t in ('LP', 'Du', 'Du.C', 'Any', 'Hal', 'Het', 'Hev'):
        assert table[t].element == '', f'{t} should have empty element, got {table[t].element!r}'


def test_sybyl_types_carbon_sp3_hybridization():
    """C.3 maps to hybridization 1 (sp3)."""
    from chython.chemistry._tables import sybyl_types
    assert sybyl_types()['C.3'].hybridization == 1


def test_sybyl_types_aromatic_hybridization():
    """C.ar and N.ar map to hybridization 4 (aromatic)."""
    from chython.chemistry._tables import sybyl_types
    table = sybyl_types()
    assert table['C.ar'].hybridization == 4
    assert table['N.ar'].hybridization == 4


def test_sybyl_types_sulfone_cumulated():
    """S.o2 and S.O2 map to hybridization 5 (cumulated, chython's z5)."""
    from chython.chemistry._tables import sybyl_types
    table = sybyl_types()
    assert table['S.o2'].hybridization == 5
    assert table['S.O2'].hybridization == 5


# packaging gate

def test_sybyl_types_tsv_ships_in_package():
    """The TSV is reachable through importlib.resources, which is the failure mode
    test_packaging.py cannot catch in a source checkout.
    """
    from importlib.resources import files
    text = files('chython.chemistry').joinpath('tables/sybyl_types.tsv').read_text(encoding='utf-8')
    assert 'C.3' in text
    assert 'C.ar' in text


# atom identity

def test_a_duplicate_atom_id_is_reported_and_binds_no_bonds():
    """Two ATOM lines stating the same id: the first claim wins and the second is reachable by
    nobody.  Rebinding the id to the later atom moves every bond that names it; which atom the
    bond meant is not ours to decide, so the second is reported as unnameable instead.
    """
    text = _record(
        'test', ' 3 1 0 0 0', 'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        '      1 C1  0.0 0.0 0.0 C.3  1 LIG 0.0',
        '      1 O1  1.5 0.0 0.0 O.3  1 LIG 0.0',   # the same id again
        '      3 N1  3.0 0.0 0.0 N.3  1 LIG 0.0',
        '@<TRIPOS>BOND',
        '      1 1 3 1',
    )
    log = []
    mol = _mol(text, log=log)
    assert mol.atom_count == 3
    sids = list(mol.atom_numbers)
    # The bond named id 1, which is the carbon: the first atom, not the oxygen.
    assert mol.element_of(sids[0]) == 6
    assert list(mol.neighbors_of(sids[0])) == [sids[2]]
    assert list(mol.neighbors_of(sids[1])) == []
    assert any('already stated' in rec for rec in log), f'expected a duplicate-id log, got: {log}'


def test_an_unreadable_atom_id_is_not_replaced_by_an_invented_one():
    """An id nobody can read becomes no id, not a line number.

    A line number lives in the same space as a stated id, so an invented value can equal a real
    one further down the block and quietly take its bonds.
    """
    text = _record(
        'test', ' 3 1 0 0 0', 'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        '      X C1  0.0 0.0 0.0 C.3  1 LIG 0.0',   # unreadable id, on line 1
        '      1 O1  1.5 0.0 0.0 O.3  1 LIG 0.0',   # a real id 1 -- what an invented 1 would shadow
        '      3 N1  3.0 0.0 0.0 N.3  1 LIG 0.0',
        '@<TRIPOS>BOND',
        '      1 1 3 1',
    )
    log = []
    mol = _mol(text, log=log)
    assert mol.atom_count == 3
    sids = list(mol.atom_numbers)
    # The bond named 1, which is the oxygen; the carbon whose id was unreadable holds no bond.
    assert list(mol.neighbors_of(sids[0])) == []
    assert sorted(mol.neighbors_of(sids[1])) == [sids[2]]
    assert any('not an integer' in rec for rec in log), f'expected an unreadable-id log, got: {log}'


def test_every_atom_line_is_named_the_same_way():
    """One identity scheme in the log: the file line and the id the file stated on it.  The four
    possible spellings (stated id, line index within the section, position among the atoms kept,
    stable id) disagree for a block holding a pseudo-atom, and a user matching a log line to a
    text editor needs one convention.
    """
    text = _record(
        'test', ' 4 0 0 0 0', 'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        '     11 C1  0.0 0.0 0.0 C.3   1 LIG 0.0',
        '     12 LP1 0.5 0.5 0.5 LP    1 LIG 0.0',
        '     13 C2  x   0.0 0.0 C.3   1 LIG 0.0',
        '     14 Q1  2.0 0.0 0.0 Zz.9  1 LIG 0.0',
        '@<TRIPOS>BOND',
    )
    log = []
    _mol(text, log=log)
    assert len(log) == 3, log
    for lineno, file_id in ((2, 12), (3, 13), (4, 14)):
        assert any(f'atom line {lineno} (id {file_id})' in rec for rec in log), (lineno, file_id, log)


# counts

def test_a_header_of_zero_against_a_populated_block_is_reported():
    """A stated 0 is a claim, not the absence of one, so both directions of the mismatch are
    reported: treating 0 as 'the writer said nothing' silences a writer who said something false.
    """
    text = _record(
        'test', ' 0 0 0 0 0', 'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        '      1 C1  0.0 0.0 0.0 C.3  1 LIG 0.0',
        '      2 C2  1.5 0.0 0.0 C.3  1 LIG 0.0',
        '@<TRIPOS>BOND',
        '      1 1 2 1',
    )
    log = []
    mol = _mol(text, log=log)
    assert mol.atom_count == 2 and mol.bond_count == 1
    assert any('claims 0 atoms' in rec for rec in log), log
    assert any('claims 0 bonds' in rec for rec in log), log


def test_a_self_consistent_record_holding_an_unstorable_bond_is_not_called_broken():
    """Both count checks compare the header to the number of block LINES: a bond we chose not to
    store is our limitation, and comparing against what was stored accuses a file that counted its
    own lines correctly.  Any legal `du` or `nc` bond and any pseudo-atom endpoint hits this.
    """
    text = _record(
        'test', ' 3 2 0 0 0', 'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        '      1 C1  0.0 0.0 0.0 C.3  1 LIG 0.0',
        '      2 C2  1.5 0.0 0.0 C.3  1 LIG 0.0',
        '      3 LP1 0.5 0.9 0.0 LP   1 LIG 0.0',
        '@<TRIPOS>BOND',
        '      1 1 2 1',
        '      2 1 3 1',
    )
    log = []
    mol = _mol(text, log=log)
    assert mol.bond_count == 1
    assert not [rec for rec in log if 'header claims' in rec], f'a self-consistent file was accused: {log}'
    assert all(str(rec).startswith('unsupported') for rec in log), log


# sections

def test_every_section_the_reader_does_not_claim_is_named():
    """A discarded section is named, because 'a section was ignored' cannot be looked up in the
    file.
    """
    text = _record(
        'test', ' 1 0 0 0 1', 'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        '      1 C1  0.0 0.0 0.0 C.3  1 LIG 0.0',
        '@<TRIPOS>BOND',
        '@<TRIPOS>CRYSIN',
        '  10.0 10.0 10.0 90.0 90.0 90.0 1 1',
        '@<TRIPOS>SET',
        'SET1 STATIC ATOMS <user> **** Comment',
        '@<TRIPOS>COMMENT',
        'anything at all',
    )
    log = []
    _mol(text, log=log)
    for name in ('CRYSIN', 'SET', 'COMMENT'):
        assert any(str(rec).startswith('unsupported') and name in rec for rec in log), (name, log)


def test_the_section_that_holds_formal_charges_names_that_consequence():
    """UNITY writers put FORMAL CHARGES in UNITY_ATOM_ATTR, so the line names that section and
    says the atoms kept the ATOM block's column instead.
    """
    text = _record(
        'test', ' 1 0 0 0 0', 'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        '      1 N1  0.0 0.0 0.0 N.4  1 LIG 0.0',
        '@<TRIPOS>BOND',
        '@<TRIPOS>UNITY_ATOM_ATTR',
        '1 1',
        'charge 1',
    )
    log = []
    mol = _mol(text, log=log)
    assert mol.charge_of(list(mol.atom_numbers)[0]) == 0
    named = [rec for rec in log if 'UNITY_ATOM_ATTR' in rec]
    assert named and all(str(rec).startswith('unsupported') for rec in named), log
    assert any('formal charge' in rec for rec in named), named


def test_a_blank_line_in_the_molecule_section_does_not_shift_its_fields():
    """The MOLECULE section is positional: name, counts, mol_type, charge_type.  Filtering blank
    lines out re-indexes it, so a record whose name line is empty -- legal, and common from
    converters -- loses its charge type and reads the counts line as the name.
    """
    text = _record(
        '', ' 1 0 0 0 0', 'SMALL', 'GASTEIGER', '',
        '@<TRIPOS>ATOM',
        '      1 N1  0.0 0.0 0.0 N.3  1 LIG -0.51',
        '@<TRIPOS>BOND',
    )
    log = []
    mol = _mol(text, log=log)
    assert mol.title == ''
    assert mol.charge_of(list(mol.atom_numbers)[0]) == 0, 'a partial charge became a formal one'
    assert any(str(rec).startswith('unsupported') and 'GASTEIGER' in rec for rec in log), log


# charges

def test_the_partial_charge_types_are_all_eleven():
    """Every Tripos charge type except NO_CHARGES carries continuous charges.  USER_CHARGES is
    what most real writers emit for partial charges and DICT_CHARGES came from a dictionary
    lookup; rounding either into a formal charge invents chemistry on every atom.
    """
    from ..mol2 import _PARTIAL_CHARGE_TYPES
    spec = {'DEL_RE', 'GASTEIGER', 'GAST_HUCK', 'HUCKEL', 'PULLMAN', 'GAUSS80_CHARGES',
            'AMPAC_CHARGES', 'MULLIKEN_CHARGES', 'DICT_CHARGES', 'MMFF94_CHARGES', 'USER_CHARGES'}
    assert _PARTIAL_CHARGE_TYPES == spec


@pytest.mark.parametrize('charge_type', ['USER_CHARGES', 'DICT_CHARGES', 'Gasteiger'])
def test_a_declared_partial_charge_type_is_honoured_however_it_is_spelled(charge_type):
    """A lower-cased spelling is a writer's habit, not a statement that the charges are formal."""
    text = _record(
        'test', ' 1 0 0 0 0', 'SMALL', charge_type, '',
        '@<TRIPOS>ATOM',
        '      1 N1  0.0 0.0 0.0 N.3  1 LIG -0.51',
        '@<TRIPOS>BOND',
    )
    log = []
    mol = _mol(text, log=log)
    assert mol.charge_of(list(mol.atom_numbers)[0]) == 0
    assert any(str(rec).startswith('unsupported') and charge_type in rec for rec in log), log


def test_reading_the_charge_column_with_no_declared_type_is_reported():
    """The charge type line is optional and the column is read as formal charges when it is missing.
    That is a decision about the file's meaning, so it goes in the log rather than being assumed.
    """
    text = '\n'.join(['@<TRIPOS>MOLECULE', 'test', ' 1 0 0 0 0', '',
                      '@<TRIPOS>ATOM',
                      '      1 N1  0.0 0.0 0.0 N.4  1 LIG 1.0',
                      '@<TRIPOS>BOND'])
    log = []
    mol = _mol(text, log=log)
    assert mol.charge_of(list(mol.atom_numbers)[0]) == 1
    assert any('no charge type' in rec for rec in log), log


# bond types

def test_an_unknown_order_bond_is_reported_like_its_three_siblings():
    """`un` means "order unknown" and there is no such order here, so single is a guess and the
    guess is reported, like the other three un-representable Tripos bond types.
    """
    text = _record(
        'test', ' 2 1 0 0 0', 'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        '      1 C1  0.0 0.0 0.0 C.3  1 LIG 0.0',
        '      2 C2  1.5 0.0 0.0 C.3  1 LIG 0.0',
        '@<TRIPOS>BOND',
        '      1 1 2 un',
    )
    log = []
    mol = _mol(text, log=log)
    sids = list(mol.atom_numbers)
    assert mol.order_of(sids[0], sids[1]) == 1
    assert any(str(rec).startswith('unsupported') and '"un"' in rec for rec in log), log


def test_a_bond_type_the_format_does_not_define_is_the_files_fault_not_ours():
    """Writers emit '4' for aromatic although Tripos spells it 'ar'.  Reading it as aromatic is the
    right recovery and it is not an `unsupported:` case: the file is wrong and we coped.
    """
    text = _record(
        'test', ' 2 1 0 0 0', 'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        '      1 C1  0.0 0.0 0.0 C.ar 1 LIG 0.0',
        '      2 C2  1.5 0.0 0.0 C.ar 1 LIG 0.0',
        '@<TRIPOS>BOND',
        '      1 1 2 4',
    )
    log = []
    mol = _mol(text, log=log)
    sids = list(mol.atom_numbers)
    assert mol.order_of(sids[0], sids[1]) == 4
    type_lines = [rec for rec in log if '"4"' in rec]
    assert type_lines and not any(str(rec).startswith('unsupported') for rec in type_lines), log


def test_an_unstorable_bond_between_an_already_bonded_pair_keeps_its_own_reason():
    """The bond type is resolved before the duplicate test, because a `du` bond is not a duplicate of
    anything -- it is a construct we do not model, and reporting it as a duplicate loses the only line
    that explains why the file's bond count and ours differ.
    """
    text = _record(
        'test', ' 2 2 0 0 0', 'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        '      1 C1  0.0 0.0 0.0 C.3  1 LIG 0.0',
        '      2 C2  1.5 0.0 0.0 C.3  1 LIG 0.0',
        '@<TRIPOS>BOND',
        '      1 1 2 1',
        '      2 1 2 du',
    )
    log = []
    mol = _mol(text, log=log)
    assert mol.bond_count == 1
    assert any(str(rec).startswith('unsupported') and '"du"' in rec for rec in log), log
    assert not any('duplicate' in rec for rec in log), log


def test_a_bond_to_an_unrecognised_type_atom_is_not_called_a_dangling_reference():
    """The atom id WAS in the ATOM block; we could not store the atom.  Saying the id is absent
    accuses the file of an error it did not make, and both lines are our limitation, not its.
    """
    text = _record(
        'test', ' 2 1 0 0 0', 'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        '      1 C1  0.0 0.0 0.0 C.3   1 LIG 0.0',
        '      2 Q1  1.5 0.0 0.0 Zz.9  1 LIG 0.0',
        '@<TRIPOS>BOND',
        '      1 1 2 1',
    )
    log = []
    mol = _mol(text, log=log)
    assert mol.atom_count == 1 and mol.bond_count == 0
    assert all(str(rec).startswith('unsupported') for rec in log), log
    assert not any('not in the ATOM block' in rec for rec in log), log


# SYBYL types

@pytest.mark.parametrize('sybyl_type,element', [('Cr.th', 24), ('Cr.oh', 24), ('Co.oh', 27),
                                                ('Ru.oh', 44)])
def test_a_metal_type_states_a_geometry_we_have_no_code_for(sybyl_type, element):
    """The four types Tripos documents whose suffix is a coordination geometry.  The file is fine
    and the hybridization code set is the limitation, so the line is `unsupported:`.
    """
    text = _record(
        'test', ' 1 0 0 0 0', 'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        f'      1 M1  0.0 0.0 0.0 {sybyl_type} 1 LIG 0.0',
        '@<TRIPOS>BOND',
    )
    log = []
    mol = _mol(text, log=log)
    assert mol.element_of(list(mol.atom_numbers)[0]) == element
    assert any(str(rec).startswith('unsupported') and sybyl_type in rec for rec in log), log


def test_geometry_only_types_matches_tsv_qualifying_rows():
    """_GEOMETRY_ONLY_TYPES must hold exactly the TSV rows carrying a real non-hydrogen element
    and hybridization 0.

    That is how the TSV encodes a coordination-geometry type; the dict maps those same types to
    their geometry name for the log message.  Divergence either under-reports a type or names one
    that no longer exists.
    """
    from ..mol2 import _GEOMETRY_ONLY_TYPES
    from chython.chemistry._tables import sybyl_types
    table = sybyl_types()
    expected = {
        t for t, row in table.items()
        if row.element and row.element != 'H' and row.hybridization == 0
    }
    assert set(_GEOMETRY_ONLY_TYPES) == expected


def test_the_type_resolver_carries_its_own_reason():
    """Three outcomes and three reasons, returned rather than reconstructed by the caller: a
    `hybridization == 0 and element != 'H'` guess is wrong for the metal geometry rows above.
    """
    from ..mol2 import _resolve_type
    assert _resolve_type('C.3') == ('C', 1, None)
    assert _resolve_type('H') == ('H', 0, None)
    element, hybridization, note = _resolve_type('Cr.oh')
    assert (element, hybridization) == ('Cr', 0) and 'octahedral' in note
    element, hybridization, note = _resolve_type('N.bogus')
    assert (element, hybridization) == ('N', 0) and 'not interpreted' in note
    element, hybridization, note = _resolve_type('Zz.9')
    assert element == '' and 'not an element symbol' in note   # '' is unstorable, None is a pseudo-atom
    assert _resolve_type('LP') == (None, 0, None)


def test_a_stated_hybridization_that_the_bonds_contradict_is_reported():
    """Hybridization is derived from bonds and no file's claim is stored, so the claim can only be
    a check -- the one thing a SYBYL type says that the bonds do not.
    """
    text = _record(
        'test', ' 2 1 0 0 0', 'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        '      1 C1  0.0 0.0 0.0 C.ar 1 LIG 0.0',
        '      2 C2  1.5 0.0 0.0 C.ar 1 LIG 0.0',
        '@<TRIPOS>BOND',
        '      1 1 2 1',      # a single bond between two atoms typed aromatic
    )
    log = []
    mol = _mol(text, log=log)
    assert all(mol.hybridization_of(sid) == 1 for sid in mol.atom_numbers)
    assert len([rec for rec in log if 'the bonds are used' in rec]) == 2, log


def test_a_stated_hybridization_the_bonds_agree_with_says_nothing():
    """Benzene's C.ar atoms hold aromatic bonds, so there is nothing to report -- a line per atom
    of every well-formed record would destroy the log's value as a screen.
    """
    log = []
    _mol((_TEST_DIR / 'mol2_simple.mol2').read_text(encoding='utf-8'), log=log)
    assert not log, log


# indented tags

def test_an_indented_tag_is_a_tag_everywhere():
    """One lstrip-aware tag test, shared by the record splitter, the section splitter and the
    string entry point: an indented `@<TRIPOS>MOLECULE` must not become the record's title, nor an
    indented `@<TRIPOS>ATOM` an empty block.
    """
    text = '\n'.join(['   @<TRIPOS>MOLECULE', 'test', ' 1 0 0 0 0', 'SMALL', 'GASTEIGER', '',
                      '  @<TRIPOS>ATOM',
                      '      1 N1  0.0 0.0 0.0 N.3  1 LIG -0.51',
                      '\t@<TRIPOS>BOND'])
    log = []
    mol = _mol(text, log=log)
    assert mol.title == 'test'
    assert mol.atom_count == 1
    assert mol.charge_of(list(mol.atom_numbers)[0]) == 0
    assert not [rec for rec in log if 'header claims' in rec], log


# entry points

def test_a_multi_record_string_returns_the_first_record_and_says_so():
    """A multi-record string yields its first record, matching `mol()`, and the log line names
    `read_mol2()` as the call that yields them all.
    """
    text = '\n'.join([_record('first', ' 1 0 0 0 0', 'SMALL', 'NO_CHARGES', '',
                              '@<TRIPOS>ATOM',
                              '      1 C1  0.0 0.0 0.0 C.3  1 LIG 0.0',
                              '@<TRIPOS>BOND'),
                      _record('second', ' 1 0 0 0 0', 'SMALL', 'NO_CHARGES', '',
                              '@<TRIPOS>ATOM',
                              '      1 O1  0.0 0.0 0.0 O.3  1 LIG 0.0',
                              '@<TRIPOS>BOND')])
    log = []
    mol = _mol(text, log=log)
    assert mol.title == 'first'
    assert mol.element_of(list(mol.atom_numbers)[0]) == 6
    assert any('read_mol2' in rec for rec in log), f'the discarded record was not named: {log}'


def test_a_repeated_section_is_reported_rather_than_overwritten():
    """The section split is a list, not a dict: a dict cannot hold two sections of one name, so a
    second ATOM block would silently replace the first.
    """
    text = _record(
        'test', ' 2 0 0 0 0', 'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        '      1 C1  0.0 0.0 0.0 C.3  1 LIG 0.0',
        '@<TRIPOS>ATOM',
        '      2 O1  1.5 0.0 0.0 O.3  1 LIG 0.0',
    )
    log = []
    mol = _mol(text, log=log)
    assert mol.atom_count == 2
    assert sorted(mol.element_of(sid) for sid in mol.atom_numbers) == [6, 8]
    assert any('ATOM' in rec and 'second' in rec for rec in log), log


def test_a_misspelled_path_fails_at_the_call_site_however_it_is_spelled():
    """A `str` with no `@<TRIPOS>` in it is a path, and the open must not be deferred into the
    generator: that turns a typo into 'this file has no molecules'.
    """
    with raises(FileNotFoundError):
        read_mol2('does_not_exist.mol2')
    with raises(FileNotFoundError):
        read_mol2(Path('does_not_exist.mol2'))


def test_a_record_that_cannot_be_parsed_does_not_cost_the_tail_of_the_file():
    """An unparseable record is filed as a `FailedRecord` and reading continues, so one bad line
    does not cost the rest of the file.  Same shape as the CTfile side, so a caller reading both
    learns one vocabulary.
    """
    from ..ctfile import FailedRecord

    def rec(name, atom):
        return _record(name, ' 1 0 0 0 0', 'SMALL', 'NO_CHARGES', '',
                       '@<TRIPOS>ATOM', atom, '@<TRIPOS>BOND')

    text = '\n'.join([rec('good_one', '      1 C1 0.0 0.0 0.0 C.3 1 LIG 0.0'),
                      rec('broken', '      1 C1 0.0 0.0'),
                      rec('good_two', '      1 O1 0.0 0.0 0.0 O.3 1 LIG 0.0')])
    results = list(read_mol2(text))
    assert len(results) == 3
    assert [getattr(m, 'title', None) for m, _ in results] == ['good_one', None, 'good_two']
    failed, failed_log = results[1]
    assert isinstance(failed, FailedRecord)
    assert failed.position == 1
    assert 'broken' in failed.text
    assert any(str(rec).startswith('record:') for rec in failed_log), failed_log


# the molecule's own log


def test_a_records_damage_reaches_the_molecule_with_nothing_passed_in():
    """`mol.log` is the destination, so a caller who passed no `log=` still has the damage report.

    [mutant: drop the `mol.log.absorb('read', own)` in `_parse_record`]
    """
    text = _record(
        'test', ' 2 1 0 0 0', 'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        '      1 C1  0.0 0.0 0.0 C.3  1 LIG 0.0',
        '      2 C2  1.5 0.0 0.0 C.3  1 LIG 0.0',
        '@<TRIPOS>BOND',
        '      1 1 2 xx',
    )
    mol = mol2_mol(text)                              # no log= anywhere
    assert any(rec.rule == 'mol2:unknown-bond-type' for rec in mol.log), mol.log
    assert all(rec.stage == 'read' for rec in mol.log), mol.log
    assert mol.log.repaired(), mol.log                 # the stage and severity fields are readable


def test_the_caller_list_and_the_molecule_log_hold_the_same_lines():
    """The `log=` list is a second copy of the same records, not the storage.

    [mutant: return the record's own list from `_parse_record` without extending the caller's]
    """
    text = _record(
        'test', ' 1 0 0 0 0', 'SMALL', 'NO_CHARGES', '',
        '@<TRIPOS>ATOM',
        '      1 N1  0.0 0.0 0.0 N.9  1 LIG 0.0',
    )
    log: list = []
    mol = mol2_mol(text, log=log)
    assert [str(rec) for rec in log] == [str(rec) for rec in mol.log], (log, mol.log)


def test_each_molecule_of_a_multi_record_file_holds_only_its_own_lines():
    """Record 2's damage is on molecule 2 and on no other molecule: never pooled.

    [mutant: absorb the caller's flat list instead of the record's own]
    """
    from ..mol2 import mol2

    def rec(name, atom, bond_type):
        return _record(name, ' 2 1 0 0 0', 'SMALL', 'NO_CHARGES', '',
                       '@<TRIPOS>ATOM',
                       f'      1 {atom}  0.0 0.0 0.0 C.3  1 LIG 0.0',
                       '      2 C2  1.5 0.0 0.0 C.3  1 LIG 0.0',
                       '@<TRIPOS>BOND',
                       f'      1 1 2 {bond_type}')

    log: list = []
    molecules = mol2('\n'.join([rec('clean', 'C1', '1'),
                                rec('damaged', 'C1', 'xx'),
                                rec('clean_again', 'C1', '1')]), log=log)
    assert [m.title for m in molecules] == ['clean', 'damaged', 'clean_again']
    assert not molecules[0].log and not molecules[2].log, [m.log for m in molecules]
    assert [rec.rule for rec in molecules[1].log] == ['mol2:unknown-bond-type']
    assert len(log) == 1                               # and the caller's flat list holds it once


def test_the_first_of_several_records_is_told_it_was_the_first():
    """`mol2_mol` on a multi-record string reports the choice on the molecule it returned.

    [mutant: log `mol2:multiple-records` to the caller's list only]
    """
    def rec(name):
        return _record(name, ' 1 0 0 0 0', 'SMALL', 'NO_CHARGES', '',
                       '@<TRIPOS>ATOM',
                       '      1 C1  0.0 0.0 0.0 C.3  1 LIG 0.0')

    mol = mol2_mol('\n'.join([rec('first'), rec('second')]))
    assert mol.title == 'first'
    assert [rec.rule for rec in mol.log] == ['mol2:multiple-records']
