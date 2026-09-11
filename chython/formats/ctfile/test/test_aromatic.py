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
"""Bond type 4 -- a delocalised bond -- is read and written by both CTAB versions.

Every assertion is against the literal columns of the emitted line, not against the writer having
returned without raising.  The second half covers the ``vvv`` total-valence column, which the writer
declines for an atom holding an aromatic bond: a face-value drawn sum counts such a bond as 4, giving
a pyrrole nitrogen 9 against a real total valence of 3.  The S-group carries the count instead.
"""

from pytest import mark

from ....core import read_smiles
from .._hydrogens import MRV_IMPLICIT_H, valence_for_write
from .._v2000 import emit_v2000, parse_v2000
from .._v3000 import emit_v3000, parse_v3000


#: The V2000 atom line's `vvv` field -- the stated total valence, 0 meaning "not stated".
_VALENCE = (48, 51)


def _v2000_record(symbols, bonds, title='aromatic'):
    """A V2000 record from ``['C', 'N', ...]`` and ``[(a, b, type), ...]``, 1-based bonds."""
    lines = [title, '  test', '',
             f'{len(symbols):3d}{len(bonds):3d}  0  0  0  0            999 V2000']
    for i, symbol in enumerate(symbols):
        lines.append(f'{float(i):10.4f}    0.0000    0.0000 {symbol:<3} 0  0  0  0  0  0'
                     f'  0  0  0  0  0')
    for a, b, order in bonds:
        lines.append(f'{a:3d}{b:3d}{order:3d}  0  0  0  0')
    lines.append('M  END')
    return lines


def _read_v2000(lines):
    return parse_v2000(lines, []).build()


def _bond_block_v2000(lines):
    """``[(a, b, type), ...]`` read back out of an emitted V2000 bond block, by column.

    Columns and not a regex: the type must land in characters 6 to 9, which is all a fixed-column
    reader looks at.
    """
    count = int(lines[3][0:3]), int(lines[3][3:6])
    atoms, bonds = count
    out = []
    for line in lines[4 + atoms:4 + atoms + bonds]:
        out.append((int(line[0:3]), int(line[3:6]), int(line[6:9])))
    return out


def _bond_block_v3000(lines):
    """``[(a, b, type), ...]`` from an emitted V3000 bond block, in its own free format."""
    out = []
    inside = False
    for line in lines:
        body = line[7:] if line.startswith('M  V30 ') else ''
        if body == 'BEGIN BOND':
            inside = True
        elif body == 'END BOND':
            inside = False
        elif inside and body:
            _, order, a, b = body.split()[:4]
            out.append((int(a), int(b), int(order)))
    return out


def _state(mol):
    """Everything about a molecule this round trip must preserve, in position terms.

    Structural and not a SMILES comparison: agreeing SMILES strings are not evidence of the same
    molecule.
    """
    sids = list(mol.atom_numbers)
    index = {s: i for i, s in enumerate(sids)}
    atoms = [(mol.element_of(s), mol.charge_of(s), mol.isotope_of(s), bool(mol.radical_of(s)),
              mol.implicit_h_of(s)) for s in sids]
    bonds = sorted((min(index[b.n], index[b.m]), max(index[b.n], index[b.m]), b.order)
                   for b in mol.bonds())
    return atoms, bonds


#: Benzene as six type-4 bonds.
_BENZENE = (['C'] * 6, [(1, 2, 4), (2, 3, 4), (3, 4, 4), (4, 5, 4), (5, 6, 4), (6, 1, 4)])

#: Styrene: an aromatic ring and an exocyclic Kekule double bond, in one bond block.
_MIXED = (['C'] * 8,
          [(1, 2, 4), (2, 3, 4), (3, 4, 4), (4, 5, 4), (5, 6, 4), (6, 1, 4), (1, 7, 1), (7, 8, 2)])


# the round trip

def test_a_v2000_record_of_type_four_bonds_writes_type_four_bonds_back():
    """Asserted against the bond type in its own three columns, because "no exception" is satisfied
    by a writer that emits nothing.  The write log must be empty too: benzene has nothing to report,
    every carbon's count being derivable, so neither write channel fires.
    """
    mol, store, log = _read_v2000(_v2000_record(*_BENZENE))
    assert log == [], 'the reader takes an aromatic record without comment; that is the premise'
    assert mol.aromatic_bond_count == 6

    lines, write_log = emit_v2000(mol, store)
    assert write_log == [], write_log
    types = [t for _, _, t in _bond_block_v2000(lines)]
    assert types == [4] * 6, lines
    assert '  1  2  4  0  0  0  0' in lines, lines


def test_a_v3000_record_of_type_four_bonds_writes_type_four_bonds_back():
    """The same for V3000."""
    mol, store, _ = _read_v2000(_v2000_record(*_BENZENE))
    lines, _ = emit_v3000(mol, store)
    types = [t for _, _, t in _bond_block_v3000(lines)]
    assert types == [4] * 6, lines
    assert 'M  V30 1 4 1 2' in lines, lines


def test_read_write_read_keeps_the_six_aromatic_bonds():
    """The molecule that comes back out is the molecule that went in, compared as a structure."""
    mol, store, _ = _read_v2000(_v2000_record(*_BENZENE))
    lines, _ = emit_v2000(mol, store)
    again, _, log = _read_v2000(lines)
    assert _state(again) == _state(mol)
    assert again.aromatic_bond_count == 6
    assert log == [], log


# the caller's explicit route

def test_kekule_first_still_writes_alternating_orders():
    """``kekule()`` is the caller's route to an alternating file, and the writer does not
    second-guess it: it writes the representation the molecule is in, never 4 for every ring bond.
    """
    mol, store, _ = _read_v2000(_v2000_record(*_BENZENE))
    mol.kekule()
    assert mol.aromatic_bond_count == 0

    lines, _ = emit_v2000(mol, store)
    types = sorted(t for _, _, t in _bond_block_v2000(lines))
    assert types == [1, 1, 1, 2, 2, 2], lines
    assert 4 not in types

    lines3, _ = emit_v3000(mol, store)
    assert sorted(t for _, _, t in _bond_block_v3000(lines3)) == [1, 1, 1, 2, 2, 2], lines3


def test_a_mixed_molecule_writes_both_kinds_of_bond_in_one_block():
    """Styrene: an aromatic ring and a Kekule double bond side by side -- the case a naive
    order-to-type table gets wrong by picking one representation per record.
    """
    mol, store, _ = _read_v2000(_v2000_record(*_MIXED))
    lines, _ = emit_v2000(mol, store)
    types = sorted(t for _, _, t in _bond_block_v2000(lines))
    assert types == [1, 2, 4, 4, 4, 4, 4, 4], lines

    lines3, _ = emit_v3000(mol, store)
    assert sorted(t for _, _, t in _bond_block_v3000(lines3)) == [1, 2, 4, 4, 4, 4, 4, 4], lines3


# the version pair

@mark.parametrize('fixture', [_BENZENE, _MIXED])
def test_both_writers_agree_with_each_other_and_with_the_input(fixture):
    """One molecule through both writers: both state type 4, and both read back the same structure.

    A pair-check, because changing one half of the version pair and not the other is the signature
    failure shape here.
    """
    mol, store, _ = _read_v2000(_v2000_record(*fixture))
    two, _ = emit_v2000(mol, store)
    three, _ = emit_v3000(mol, store)

    assert sorted(t for _, _, t in _bond_block_v2000(two)) \
        == sorted(t for _, _, t in _bond_block_v3000(three))
    assert 4 in [t for _, _, t in _bond_block_v2000(two)]

    from_two, _, _ = _read_v2000(two)
    from_three, _, _ = parse_v3000(three, []).build()
    assert _state(from_two) == _state(mol)
    assert _state(from_three) == _state(mol)


# the valence column, declined

def _pyrrole_nitrogen(mol):
    return next(s for s in mol.atom_numbers if mol.element_of(s) == 7)


def test_no_valence_is_stated_for_an_atom_holding_an_aromatic_bond():
    """``vvv`` stays 0 for a pyrrole nitrogen.  Column and function both asserted.

    The face-value drawn sum would give 9 -- two aromatic bonds at 4 apiece plus one hydrogen --
    against a real total valence of 3, beside an ``MRV_IMPLICIT_H`` datum stating one hydrogen.
    """
    mol = read_smiles('c1cc[nH]c1')
    nitrogen = _pyrrole_nitrogen(mol)
    assert mol.implicit_h_of(nitrogen) == 1, 'the fixture has to carry a count worth stating'
    assert valence_for_write(mol, nitrogen) is None

    lines, _ = emit_v2000(mol)
    atom_lines = lines[4:4 + len(list(mol.atom_numbers))]
    line = next(x for x in atom_lines if x[31:34].strip() == 'N')
    assert line[_VALENCE[0]:_VALENCE[1]] == '  0', f'a valence was stated after all: {line!r}'
    assert all(x[_VALENCE[0]:_VALENCE[1]] == '  0' for x in atom_lines), atom_lines

    lines3, _ = emit_v3000(mol)
    assert not any('VAL=' in x for x in lines3), lines3


def test_a_bonded_atom_with_no_aromatic_bond_does_get_its_valence_stated():
    """The positive the test above needs: without it, ``if mol.degree_of(sid): return None`` --
    decline for any bonded atom -- passes the whole package.

    Methanesulfinyl hydride: the sulfur draws 1 + 2 and carries one hydrogen, total 4, which the
    valence rules do not reproduce, so both channels fire and agree.
    """
    mol = read_smiles('C[SH]=O')
    sulfur = next(s for s in mol.atom_numbers if mol.element_of(s) == 16)
    assert mol.degree_of(sulfur) == 2, 'the atom has to be bonded or the gate is not exercised'
    assert valence_for_write(mol, sulfur) == 4

    lines, log = emit_v2000(mol)
    atom_lines = lines[4:4 + len(list(mol.atom_numbers))]
    line = next(x for x in atom_lines if x[31:34].strip() == 'S')
    assert line[_VALENCE[0]:_VALENCE[1]] == '  4', f'no valence was stated: {line!r}'

    lines3, _ = emit_v3000(mol)
    assert any('VAL=4' in x for x in lines3), lines3

    # The aromatic caveat must not appear on a molecule holding no aromatic bond.
    assert len(log) == 1, log
    assert 'aromatic' not in log[0], log
    assert 'and in the valence field' in log[0], log


def test_the_hydrogen_count_still_goes_out_and_still_comes_back():
    """Silence in ``vvv`` costs nothing: the S-group is the channel that survives a read, and
    without this the test above is indistinguishable from dropping the count.
    """
    mol = read_smiles('c1cc[nH]c1')
    # Bound *and* asserted: the checks below match substrings of the file and the log, which a
    # molecule with no nitrogen would also satisfy.
    nitrogen = _pyrrole_nitrogen(mol)
    assert mol.implicit_h_of(nitrogen) == 1, 'the fixture states no count, so there is nothing to carry'

    lines, log = emit_v2000(mol)
    assert any(MRV_IMPLICIT_H in x for x in lines), lines
    assert any('IMPL_H1' in x for x in lines), lines
    assert any('in the valence field for the 0 of them' in x for x in log), log

    again, _, read_log = _read_v2000(lines)
    back = _pyrrole_nitrogen(again)
    assert again.implicit_h_of(back) == 1, read_log
    assert read_log == [], read_log

    lines3, _ = emit_v3000(mol)
    assert any('IMPL_H1' in x for x in lines3), lines3
    from_three, _, log3 = parse_v3000(lines3, []).build()
    assert from_three.implicit_h_of(_pyrrole_nitrogen(from_three)) == 1, log3


# the count that is not known, either version

#: A pyrrole ring drawn as five type-4 bonds and nothing else -- no valence, no data S-group, so no
#: sentinel had to be planted.  The nitrogen's count is genuinely undetermined: it either donates its
#: lone pair and carries a hydrogen or takes a ring double bond and carries none, and only the ring
#: decides.
_PYRROLE_SKELETON = (['C', 'C', 'C', 'N', 'C'],
                     [(1, 2, 4), (2, 3, 4), (3, 4, 4), (4, 5, 4), (5, 1, 4)])


def test_an_unknown_count_on_an_aromatic_atom_is_written_as_silence_and_read_back_as_unknown():
    """Both write channels stay silent, for two different reasons.  ``vvv`` is declined because a
    total valence is not computable here; ``MRV_IMPLICIT_H`` because there is no count to state, and
    ``IMPL_H0`` would turn "nobody said" into "there are none".  Either fallback would round-trip as a
    definite count and look like agreement.  Asserted on the emitted text and on the molecule that
    comes back, since only the second half proves the absence was read as one.
    """
    mol, store, log = _read_v2000(_v2000_record(*_PYRROLE_SKELETON))
    nitrogen = _pyrrole_nitrogen(mol)
    assert mol.implicit_h_of(nitrogen) is None, 'the fixture has to arrive with the count unknown'
    assert mol.unknown_h_count == 1
    assert any('only the ring decides' in x for x in log), log

    lines, write_log = emit_v2000(mol, store)
    atom_lines = lines[4:4 + len(list(mol.atom_numbers))]
    assert all(x[_VALENCE[0]:_VALENCE[1]] == '  0' for x in atom_lines), atom_lines
    assert not any(MRV_IMPLICIT_H in x for x in lines), lines
    assert not any('IMPL_H' in x for x in lines), lines
    assert any('no known implicit hydrogen count' in x for x in write_log), write_log

    again, _, read_log = _read_v2000(lines)
    back = _pyrrole_nitrogen(again)
    assert again.implicit_h_of(back) is None, f'the unknown came back as a number: {read_log}'
    assert again.unknown_h_count == 1
    assert _state(again) == _state(mol)

    lines3, _ = emit_v3000(mol, store)
    assert not any('VAL=' in x for x in lines3), lines3
    assert not any('IMPL_H' in x for x in lines3), lines3
    from_three, _, log3 = parse_v3000(lines3, []).build()
    assert from_three.implicit_h_of(_pyrrole_nitrogen(from_three)) is None, log3
    assert _state(from_three) == _state(mol)


def test_kekulising_that_same_record_first_settles_the_count_without_stating_it():
    """The other end of the same fixture, so the silence above is not a dead end.

    After ``kekule()`` the nitrogen's class is settled: one hydrogen.  The file still states no count,
    but for the opposite reason -- the valence rules reproduce it from the bonds drawn -- so the
    molecule that comes back carries the count rather than the sentinel.
    """
    mol, store, _ = _read_v2000(_v2000_record(*_PYRROLE_SKELETON))
    assert not mol.kekule().unresolved
    nitrogen = _pyrrole_nitrogen(mol)
    assert mol.implicit_h_of(nitrogen) == 1
    assert mol.unknown_h_count == 0

    lines, write_log = emit_v2000(mol, store)
    assert not any('IMPL_H' in x for x in lines), lines
    assert valence_for_write(mol, nitrogen) is None, 'the rules reproduce it, so it is not stated'
    assert not any('no known implicit hydrogen count' in x for x in write_log), write_log

    again, _, read_log = _read_v2000(lines)
    assert again.implicit_h_of(_pyrrole_nitrogen(again)) == 1, read_log
    assert again.unknown_h_count == 0, read_log
    assert sorted(b.order for b in again.bonds()) == [1, 1, 1, 2, 2], 'and it is a Kekule file now'
