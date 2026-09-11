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
"""R-atom read and write: the four V2000 spellings and the V3000 R#/RGROUPS dialect."""
from chython.core import read_smiles
from .._v2000 import emit_v2000, parse_v2000
from .._v3000 import emit_v3000, parse_v3000
from .test_v2000 import _ETHANOL


def _with_symbol(symbol, extra=()):
    """Ethanol with its first atom's symbol column replaced, plus any extra property lines."""
    lines = list(_ETHANOL)
    lines[4] = f'    0.0000    0.0000    0.0000 {symbol:<3s} 0  0  0  0  0  0  0  0  0  0  0  0'
    if extra:
        lines[-1:-1] = list(extra)          # before `M  END`
    mol, store, log = parse_v2000(lines, []).build()
    return mol, store, log


def _the_r(mol):
    return next(a for a in mol.atoms() if a.is_r)


def test_r_hash_with_rgp():
    mol, _, _ = _with_symbol('R#', ('M  RGP  1   1   3',))
    assert _the_r(mol).r_index == 3


def test_r_digits_in_the_symbol_column():
    mol, _, _ = _with_symbol('R7')
    assert _the_r(mol).r_index == 7


def test_bare_r():
    mol, _, _ = _with_symbol('R')
    assert _the_r(mol).r_index == 0


def test_star():
    mol, _, _ = _with_symbol('*')
    assert _the_r(mol).r_index == 0


def test_a_case_folded_marker_reads_the_same_as_an_upper_case_one():
    # This column is fixed-width text and writers fold whole records, so `r7` is `R7`.  All six
    # spellings answer the same marker and the same index.
    for spelling, index in (('r', 0), ('r#', 0), ('r7', 7), ('R', 0), ('R#', 0), ('R7', 7)):
        mol, _, _ = _with_symbol(spelling)
        assert _the_r(mol).r_index == index, spelling


def test_the_elements_that_start_with_r_are_still_elements():
    for spelling in ('Rb', 'rb', 'RU', 'Rn', 'Re', 'Ra', 'Rf', 'Rg'):
        mol, _, _ = _with_symbol(spelling)
        assert not next(iter(mol.atoms())).is_r, spelling


def test_r_hash_without_rgp_reads_as_a_plain_r():
    # A file that names R# and never states the group: an underivable index is 0, not a refusal.
    mol, _, _ = _with_symbol('R#')
    assert _the_r(mol).r_index == 0


def test_the_record_is_whole_and_the_marker_is_not_an_alias():
    # The other two atoms and both bonds survive, and the marker is an element rather than a label:
    # nothing goes in `aliases`, and the atom is not the placeholder carbon a free-text symbol gets.
    mol, store, _ = _with_symbol('R1')
    assert len([*mol.atom_numbers]) == 3
    assert len(list(mol.bonds())) == 2
    assert mol.aliases == {}
    assert store.aliases == {}


def test_the_neighbour_of_an_r_keeps_the_hydrogens_a_carbon_neighbour_would_leave():
    # Ethanol's first carbon has one carbon neighbour and three hydrogens; with the marker in its
    # place the second carbon's count must not move.
    marked, _, _ = _with_symbol('R')
    plain, _, _ = _with_symbol('C')
    second = [*marked.atom_numbers][1]
    assert marked.implicit_h_of(second) == plain.implicit_h_of([*plain.atom_numbers][1])


def test_rgp_naming_an_absent_atom_is_logged_not_fatal():
    mol, _, log = _with_symbol('R#', ('M  RGP  1   9   3',))
    assert len([*mol.atom_numbers]) == 3
    assert _the_r(mol).r_index == 0
    assert any('RGP' in x for x in log), log


def test_rgp_on_an_atom_that_is_not_a_marker_is_logged_not_fatal():
    # `Ctab.build` applies an index only to element 0, so an `M  RGP` aimed elsewhere is dropped
    # there and logged rather than reaching `set_r_index`.
    mol, _, log = _with_symbol('R', ('M  RGP  1   2   3',))
    assert len([*mol.atom_numbers]) == 3
    assert any('RGP' in x for x in log), log


def test_an_index_past_the_cap_is_logged_not_fatal():
    mol, _, log = _with_symbol('R#', ('M  RGP  1   1 100',))
    assert _the_r(mol).r_index == 0
    assert any('RGP' in x for x in log), log


def _write(mol):
    lines, _ = emit_v2000(mol)
    return lines


def _reread(mol):
    """`mol` written as V2000 and read straight back."""
    back, _, _ = parse_v2000(_write(mol), []).build()
    return back


def _atom_line(lines, needle):
    """The atom line holding `needle`, or None -- an absent line then fails an assert, not the test."""
    return next((line for line in lines if needle in line), None)


def test_an_unindexed_r_writes_a_bare_r_and_no_rgp():
    lines = _write(read_smiles('[R]c1ccccc1'))
    assert _atom_line(lines, ' R   0') is not None, lines
    assert not any(line.startswith('M  RGP') for line in lines)


def test_an_indexed_r_writes_r_hash_plus_rgp():
    lines = _write(read_smiles('[R3]c1ccccc1'))
    assert _atom_line(lines, ' R#  0') is not None, lines
    # The marker is atom 1 of the record and its group is 3, in `M  CHG`'s count-then-pairs columns.
    assert 'M  RGP  1   1   3' in lines


def test_the_highest_index_writes_and_reads_back():
    lines = _write(read_smiles('[R99]C'))
    assert _atom_line(lines, ' R#  0') is not None, lines
    assert 'M  RGP  1   1  99' in lines
    assert next(a for a in _reread(read_smiles('[R99]C')).atoms() if a.is_r).r_index == 99


def test_the_r_line_states_no_valence():
    # `vvv` 0 is "no marking".  A marker holds a bond, so its total valence is not zero, and 15 --
    # which states zero out loud -- would contradict the bond block.
    line = _atom_line(_write(read_smiles('[R]C')), ' R   0')
    assert line[48:51] == '  0'


def test_the_atom_line_keeps_its_width():
    # `R#` is the same three characters as any element symbol, so no field right of it shifts.
    lines_r = _write(read_smiles('[R7]C'))
    lines_c = _write(read_smiles('CC'))
    r_line = _atom_line(lines_r, ' R#  0')
    c_line = _atom_line(lines_c, ' C   0')
    assert r_line is not None, lines_r
    assert c_line is not None, lines_c
    assert len(r_line) == len(c_line)


def test_round_trip_keeps_the_index():
    mol = _reread(read_smiles('[R3]c1ccccc1'))
    assert next(a for a in mol.atoms() if a.is_r).r_index == 3
    assert mol.atom_count == 7


def test_round_trip_keeps_two_different_indices_apart():
    mol = _reread(read_smiles('[R1]CCC[R2]'))
    assert sorted(a.r_index for a in mol.atoms() if a.is_r) == [1, 2]


def test_round_trip_keeps_an_unindexed_marker_unindexed():
    mol = _reread(read_smiles('[R]CC[R]'))
    markers = [a for a in mol.atoms() if a.is_r]
    assert len(markers) == 2
    assert [a.r_index for a in markers] == [0, 0]


def test_nine_markers_wrap_the_rgp_block_at_eight():
    # Eight entries per line is the format's limit, not a wrapping preference.
    mol = read_smiles('C(' + ')('.join(f'[R{i}]' for i in range(1, 10)) + ')C')
    rgp = [line for line in _write(mol) if line.startswith('M  RGP')]
    assert len(rgp) == 2
    assert rgp[0].startswith('M  RGP  8')
    assert rgp[1].startswith('M  RGP  1')
    assert sorted(a.r_index for a in _reread(mol).atoms() if a.is_r) == list(range(1, 10))


# V3000

def _v3000_lines(mol):
    lines, _ = emit_v3000(mol)
    return lines


def _v3000_atom_line(lines, position):
    """The `M  V30` atom line for the atom at 1-based `position`, without its tag."""
    body = iter(lines)
    for line in body:
        if line.strip().endswith('BEGIN ATOM'):
            break
    atoms = []
    for line in body:
        if line.strip().endswith('END ATOM'):
            break
        atoms.append(line.split('V30 ', 1)[1])
    return atoms[position - 1]


def _read_v3000(text_lines):
    mol, store, log = parse_v3000(text_lines, []).build()
    return mol, store, log


def _v3000_record(atom_field):
    """A two-atom V3000 record whose first atom line carries `atom_field` after the index."""
    return ['', '', '', '  0  0  0  0  0  0            999 V3000',
            'M  V30 BEGIN CTAB',
            'M  V30 COUNTS 2 1 0 0 0',
            'M  V30 BEGIN ATOM',
            f'M  V30 1 {atom_field}',
            'M  V30 2 C 1.5 0.0 0.0 0',
            'M  V30 END ATOM',
            'M  V30 BEGIN BOND',
            'M  V30 1 1 1 2',
            'M  V30 END BOND',
            'M  V30 END CTAB',
            'M  END']


def test_v3000_reads_r_hash_with_rgroups():
    mol, _, _ = _read_v3000(_v3000_record('R# 0.0 0.0 0.0 0 RGROUPS=(1 5)'))
    assert next(a for a in mol.atoms() if a.is_r).r_index == 5


def test_v3000_reads_r_hash_without_rgroups_as_an_unindexed_marker():
    mol, _, _ = _read_v3000(_v3000_record('R# 0.0 0.0 0.0 0'))
    assert next(a for a in mol.atoms() if a.is_r).r_index == 0


def test_v3000_reads_a_bare_r():
    mol, _, _ = _read_v3000(_v3000_record('R 0.0 0.0 0.0 0'))
    assert next(a for a in mol.atoms() if a.is_r).r_index == 0


def test_v3000_reads_the_index_in_the_type_token():
    mol, _, _ = _read_v3000(_v3000_record('R7 0.0 0.0 0.0 0'))
    assert next(a for a in mol.atoms() if a.is_r).r_index == 7


def test_v3000_reads_star_as_a_marker():
    mol, _, _ = _read_v3000(_v3000_record('* 0.0 0.0 0.0 0'))
    assert next(a for a in mol.atoms() if a.is_r).r_index == 0


def test_v3000_a_case_folded_type_token_reads_the_same():
    # `_SYMBOL_FOLD` exists because real files case-fold whole records; the marker folds with them.
    for token, expected in (('r', 0), ('r#', 0), ('r7', 7), ('R', 0), ('R#', 0), ('R7', 7)):
        mol, _, _ = _read_v3000(_v3000_record(f'{token} 0.0 0.0 0.0 0'))
        marker = next(a for a in mol.atoms() if a.is_r)
        assert marker.r_index == expected, token


def test_v3000_the_elements_that_start_with_r_are_still_elements():
    for token in ('Rb', 'rb', 'RU', 'Rn', 'Re', 'Ra', 'Rf', 'Rg'):
        mol, _, _ = _read_v3000(_v3000_record(f'{token} 0.0 0.0 0.0 0'))
        assert not any(a.is_r for a in mol.atoms()), token


def test_v3000_rgroups_naming_more_than_one_group_keeps_the_first_and_logs():
    # `RGROUPS=(2 4 6)` is an Rgroup *member of two groups*, which one atom cannot represent here.
    mol, _, log = _read_v3000(_v3000_record('R# 0.0 0.0 0.0 0 RGROUPS=(2 4 6)'))
    assert next(a for a in mol.atoms() if a.is_r).r_index == 4
    assert any('RGROUPS' in x for x in log), log


def test_v3000_an_index_past_the_cap_is_logged_not_fatal():
    from chython.core import R_INDEX_MAX

    mol, _, log = _read_v3000(_v3000_record(f'R# 0.0 0.0 0.0 0 RGROUPS=(1 {R_INDEX_MAX + 1})'))
    assert mol.atom_count == 2
    assert next(a for a in mol.atoms() if a.is_r).r_index == 0
    assert any('RGROUPS' in x for x in log), log


def test_v3000_writes_r_hash_with_rgroups():
    line = _v3000_atom_line(_v3000_lines(read_smiles('[R5]C')), 1)
    assert line.split()[1] == 'R#'
    assert 'RGROUPS=(1 5)' in line


def test_v3000_writes_a_bare_r_with_no_rgroups():
    line = _v3000_atom_line(_v3000_lines(read_smiles('[R]C')), 1)
    assert line.split()[1] == 'R'
    assert 'RGROUPS' not in line


def test_v3000_writes_no_valence_for_a_marker():
    # A marker holds a bond, so `VAL=-1` -- V3000's "zero valence" -- would be false.
    assert 'VAL=' not in _v3000_atom_line(_v3000_lines(read_smiles('[R]C')), 1)


def test_v3000_round_trips_the_index():
    mol, _, _ = _read_v3000(_v3000_lines(read_smiles('[R5]c1ccccc1')))
    assert next(a for a in mol.atoms() if a.is_r).r_index == 5
    assert mol.atom_count == 7


def test_v3000_round_trips_two_indices_and_an_unindexed_marker():
    mol, _, _ = _read_v3000(_v3000_lines(read_smiles('[R1]CC([R])C[R2]')))
    assert sorted(a.r_index for a in mol.atoms() if a.is_r) == [0, 1, 2]


def test_v3000_rgroups_overrides_the_type_token_and_says_so():
    mol, _, log = _read_v3000(_v3000_record('R7 0.0 0.0 0.0 0 RGROUPS=(1 5)'))
    assert next(a for a in mol.atoms() if a.is_r).r_index == 5
    assert any('RGROUPS' in x for x in log), log
