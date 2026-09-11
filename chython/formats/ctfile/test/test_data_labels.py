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
"""`add_data_sgroup()`: a CTfile DAT label attached in one call.

The record, its references and its FIELDDISP anchor are already modelled -- these tests are about the
one-call door and what it computes for a caller who states only atoms.
"""

from pytest import raises

from .._facade import mol
from .._sgroup import FIELDDISP_TAIL, add_data_sgroup, data_sgroups


# butane with 2D coordinates, so an anchor can be computed rather than stated.
_BUTANE = '\n'.join(['butane', '', '',
                     '  4  3  0  0  0  0            999 V2000',
                     '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
                     '    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
                     '    2.0000    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
                     '    3.0000    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
                     '  1  2  1  0  0  0  0',
                     '  2  3  1  0  0  0  0',
                     '  3  4  1  0  0  0  0',
                     'M  END'])


def test_an_atom_label_survives_a_v3000_round_trip():
    m = mol(_BUTANE)
    n = m.atom_numbers[1]
    add_data_sgroup(m, 'StereoLabel', '(R)', atoms=[n])

    again = mol(mol(m, version=3000))
    record, = data_sgroups(again, 'StereoLabel')
    assert record.field_data == '(R)'
    assert record.atoms == [again.atom_numbers[1]]
    assert record.disp[:2] == (1.0, 0.0)


def test_an_atom_label_survives_a_v2000_round_trip():
    """The same record through `M  STY`/`M  SAL`/`M  SDT`/`M  SDD`/`M  SED`."""
    m = mol(_BUTANE)
    add_data_sgroup(m, 'StereoLabel', '(R)', atoms=[m.atom_numbers[1]])
    written = mol(m, version=2000)
    assert 'M  SDT   1 StereoLabel' in written
    record, = data_sgroups(mol(written), 'StereoLabel')
    assert record.field_data == '(R)'


def test_a_bond_label_keeps_its_bond_and_anchors_between_its_atoms():
    m = mol(_BUTANE)
    a, b = m.atom_numbers[1], m.atom_numbers[2]
    record = add_data_sgroup(m, 'StereoLabel', '(Z)', atoms=[a, b], bonds=[(a, b)])
    assert record.bonds == [(a, b)]
    assert record.disp[:2] == (1.5, 0.5)

    again = mol(mol(m, version=3000))
    back, = data_sgroups(again, 'StereoLabel')
    assert len(back.bonds) == 1


def test_the_anchor_carries_the_fixed_display_tail():
    m = mol(_BUTANE)
    record = add_data_sgroup(m, 'StereoLabel', '(R)', atoms=[m.atom_numbers[0]])
    assert record.disp == (0.0, 0.0, FIELDDISP_TAIL)


def test_a_stated_anchor_wins_and_a_stated_tail_is_kept():
    m = mol(_BUTANE)
    stated = add_data_sgroup(m, 'X', 'v', atoms=[m.atom_numbers[0]], disp=(9.0, -9.0))
    assert stated.disp == (9.0, -9.0, FIELDDISP_TAIL)
    own = add_data_sgroup(m, 'Y', 'v', atoms=[m.atom_numbers[0]], disp=(1.0, 2.0, '  DAU'))
    assert own.disp == (1.0, 2.0, '  DAU')


def test_disp_false_writes_no_anchor():
    m = mol(_BUTANE)
    record = add_data_sgroup(m, 'X', 'v', atoms=[m.atom_numbers[0]], disp=False)
    assert record.disp is None
    assert 'M  SDD' not in mol(m)


def test_no_coordinates_means_no_anchor_and_a_log_line():
    """The anchor is a display position; a molecule with no drawing has none to give."""
    from chython.core import read_smiles

    m = read_smiles('CCO')
    log = []
    record = add_data_sgroup(m, 'X', 'v', atoms=[m.atom_numbers[0]], log=log)
    assert record.disp is None
    assert any('FIELDDISP not written' in x for x in log), log


def test_a_second_label_does_not_replace_the_first():
    """The core's `set_sgroups` replaces the whole set; this verb appends to it."""
    m = mol(_BUTANE)
    add_data_sgroup(m, 'StereoLabel', '(R)', atoms=[m.atom_numbers[1]])
    add_data_sgroup(m, 'StereoLabel', '(S)', atoms=[m.atom_numbers[2]])
    assert sorted(r.field_data for r in data_sgroups(m, 'StereoLabel')) == ['(R)', '(S)']
    assert len({r.index for r in data_sgroups(m)}) == 2


def test_a_multi_value_datum_keeps_its_order():
    m = mol(_BUTANE)
    record = add_data_sgroup(m, 'NOTE', ['first', 'second'], atoms=[m.atom_numbers[0]])
    assert record.data == [b'first', b'second']
    assert record.field_data == 'first\nsecond'


def test_an_atom_not_in_the_molecule_is_refused_by_number():
    m = mol(_BUTANE)
    with raises(ValueError, match='atom 99'):
        add_data_sgroup(m, 'X', 'v', atoms=[99])


def test_a_bond_that_does_not_exist_is_refused_by_its_endpoints():
    m = mol(_BUTANE)
    a, d = m.atom_numbers[0], m.atom_numbers[3]
    with raises(ValueError, match='bond'):
        add_data_sgroup(m, 'X', 'v', atoms=[a, d], bonds=[(a, d)])


def test_data_sgroups_returns_every_dat_record_when_no_name_is_given():
    m = mol(_BUTANE)
    add_data_sgroup(m, 'A', '1', atoms=[m.atom_numbers[0]])
    add_data_sgroup(m, 'B', '2', atoms=[m.atom_numbers[1]])
    assert sorted(r.name for r in data_sgroups(m)) == ['A', 'B']
    assert [r.name for r in data_sgroups(m, 'A')] == ['A']


def test_the_stereo_label_job_runs_on_chython_alone():
    """The job an RDKit script is usually written for: a descriptor drawn beside each centre as a DAT
    S-group, plus an SD data field, written V3000 and read back.

    The descriptor letters come from the caller -- chython stores CIP and computes none (`set_atom_cip`
    is storage only), so a labelling job supplies them, typically from a column of its input.
    """
    from io import StringIO

    from .._stream import ESDFWrite, SDFRead

    m = mol(_BUTANE)
    a, b = m.atom_numbers[1], m.atom_numbers[2]
    with m.edit() as e:
        e.set_atom_cip(a, 'R')
        e.set_bond_cip(a, b, 'Z')

    # what the RDKit version called CreateMolDataSubstanceGroup + SetAtoms/SetBonds
    add_data_sgroup(m, 'StereoLabel', f'({m.atom(a).cip})', atoms=[a])
    add_data_sgroup(m, 'StereoLabel', '(Z)', atoms=[a, b], bonds=[(a, b)])
    # and the hand-written `>  <KEY>` framing plus `$$$$`
    m.meta['StereoDescriptors'] = f'{a}R'

    buf = StringIO()
    with ESDFWrite(buf) as f:      # V3000; `SDFWrite` takes no version argument
        f.write(m)

    with SDFRead(StringIO(buf.getvalue())) as r:
        back = next(iter(r))

    labels = sorted(x.field_data for x in data_sgroups(back, 'StereoLabel'))
    assert labels == ['(R)', '(Z)']
    assert back.meta['StereoDescriptors'] == f'{a}R'
    assert all(x.disp is not None for x in data_sgroups(back, 'StereoLabel'))
