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
"""`mol()` and `rxn()`: one callable per format, both directions.

The direction follows the argument's type -- a container in means export, a string in means import.
`version=` applies only to export and is accepted and ignored on import.
"""

import pytest
from pytest import raises

from chython.core import MoleculeContainer, read_smiles
from chython.formats.ctfile import MalformedCtfile, mol, needs_v3000


_ETHANOL_V2000 = '\n'.join(['ethanol', '', '',
                            '  3  2  0  0  0  0            999 V2000',
                            '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
                            '    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
                            '    2.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0',
                            '  1  2  1  0  0  0  0',
                            '  2  3  1  0  0  0  0',
                            'M  END'])


def test_mol_imports_a_string():
    m = mol(_ETHANOL_V2000)
    assert isinstance(m, MoleculeContainer)
    assert m.atom_count == 3
    assert m.title == 'ethanol'


def test_mol_tolerates_the_record_separator():
    """A record copied out of an SDF brings its ``$$$$`` along; that is not an error."""
    assert mol(_ETHANOL_V2000 + '\n$$$$\n').atom_count == 3


def test_mol_tolerates_crlf():
    assert mol(_ETHANOL_V2000.replace('\n', '\r\n')).atom_count == 3


def test_mol_exports_a_molecule_as_v2000_by_default():
    out = mol(mol(_ETHANOL_V2000))
    assert isinstance(out, str)
    assert 'V2000' in out.split('\n')[3]
    assert out.rstrip().endswith('M  END')


def test_mol_round_trips_through_both_versions():
    original = mol(_ETHANOL_V2000)
    for version in (2000, 3000):
        again = mol(mol(original, version=version))
        assert again.atom_count == original.atom_count
        assert [b.order for b in again.bonds()] == [b.order for b in original.bonds()]


def test_mol_reads_v3000_with_no_keyword():
    """The version is sniffed, so import takes no version argument at all."""
    v3000 = mol(mol(_ETHANOL_V2000), version=3000)
    assert 'V3000' in v3000.split('\n')[3]
    assert mol(v3000).atom_count == 3


def test_needs_v3000_is_none_for_an_ordinary_molecule():
    assert needs_v3000(mol(_ETHANOL_V2000)) is None


@pytest.mark.parametrize('coord', [100000.0, -10000.0])
def test_needs_v3000_names_an_over_wide_coordinate(coord):
    # 100000.0 overflows on the positive side; -10000.0 overflows on the negative side (sign eats a
    # character). Both are storable in the arena but too wide for V2000's 10-character column.
    m = mol(_ETHANOL_V2000)
    with m.edit() as e:
        e.set_xy(next(iter(m.atom_numbers)), coord, 0.0)
    reason = needs_v3000(m)
    assert reason is not None and 'coordinate' in reason, reason


def test_a_charge_is_not_a_reason_to_escalate():
    """There is no charge check in emit_v2000: the ccc column takes 0 and M  CHG takes the truth."""
    m = mol(_ETHANOL_V2000)
    with m.edit() as e:
        e.set_charge(next(iter(m.atom_numbers)), 4)
    assert needs_v3000(m) is None
    assert 'M  CHG' in mol(m)


def test_explicit_v2000_does_not_escalate():
    """`version=2000` is an instruction.  An over-wide coordinate raises rather than switching."""
    m = mol(_ETHANOL_V2000)
    with m.edit() as e:
        e.set_xy(next(iter(m.atom_numbers)), 100000.0, 0.0)
    with raises(MalformedCtfile, match='V2000 10-character column'):
        mol(m, version=2000)
    assert 'V3000' in mol(m, version=3000).split('\n')[3]


def test_auto_escalates_and_logs_the_reason():
    m = mol(_ETHANOL_V2000)
    with m.edit() as e:
        e.set_xy(next(iter(m.atom_numbers)), 100000.0, 0.0)
    log = []
    out = mol(m, log=log)
    assert 'V3000' in out.split('\n')[3]
    assert any('V3000' in x and 'coordinate' in x for x in log), log


def test_a_bad_version_number_says_what_is_accepted():
    with raises(ValueError, match='2000') as exc_info:
        mol(mol(_ETHANOL_V2000), version=4000)
    assert not isinstance(exc_info.value, MalformedCtfile), \
        'a bad version= argument is a caller error, not a malformed file; must not be catchable as MalformedCtfile'


def test_mol_round_trips_data_fields():
    """INVERTED: this asserted a log line saying the data fields were dropped, which they no longer
    are -- the molecule holds them, so `mol(mol(text))` keeps them and the report would be false."""
    text = _ETHANOL_V2000 + '\n> <ACTIVITY>\n5.0\n\n'
    log = []
    m = mol(text, log=log)
    assert isinstance(m, MoleculeContainer) and m.meta == {'ACTIVITY': '5.0'}
    assert not any('not returned by mol()' in x for x in log), log
    assert mol(mol(m)).meta == {'ACTIVITY': '5.0'}, 'out on the write and back in on the read'


# Each test below asserts the negative case from the same threshold, so one edit cannot move a bound
# without breaking both halves.

@pytest.mark.parametrize('kind', [3, 2])
def test_needs_v3000_names_an_and_or_stereo_group(kind):
    """AND (kind 3) and OR (kind 2) stereo groups have no V2000 spelling.

    Why needs_v3000 is a predicate and not a try/except: emit_v2000 does not raise here, it logs and
    writes the record with the collections absent, so a catch-and-retry dispatcher loses the stereo.
    """
    m = read_smiles('ClC(Br)(F)I')
    c = next(sid for sid in m.atom_numbers if m.atom(sid).element == 6)
    m.set_stereo_group(c, kind, 1)
    reason = needs_v3000(m)
    assert reason is not None and 'stereo' in reason.lower(), reason
    log = []
    out = mol(m, log=log)
    assert 'V3000' in out.split('\n')[3]
    assert any('V3000' in x for x in log), log


def test_abs_stereo_does_not_escalate():
    """ABS (kind 1) is a single known enantiomer; V2000's chiral flag can express it."""
    m = read_smiles('ClC(Br)(F)I')
    c = next(sid for sid in m.atom_numbers if m.atom(sid).element == 6)
    m.set_stereo_group(c, 1)  # ABS, group 0 (the only valid value for kind 1)
    assert needs_v3000(m) is None


def test_needs_v3000_for_too_many_atoms():
    """V2000 can express at most 999 atoms in the 3-character count field."""
    m999 = MoleculeContainer()
    with m999.edit() as e:
        ids = [e.add_atom(6) for _ in range(999)]
        e.add_bond(ids[0], ids[1], 1)
    assert needs_v3000(m999) is None  # 999 atoms fits

    m1000 = MoleculeContainer()
    with m1000.edit() as e:
        ids = [e.add_atom(6) for _ in range(1000)]
        e.add_bond(ids[0], ids[1], 1)
    reason = needs_v3000(m1000)
    assert reason is not None and 'atom' in reason, reason  # 1000 atoms does not


_RXN_TEXT = '\n'.join(['$RXN', 'oxidation', '', '',
                       '  1  1',
                       '$MOL'] + _ETHANOL_V2000.split('\n') + ['$MOL'] + _ETHANOL_V2000.split('\n'))


def test_rxn_imports_a_string():
    from chython.core.reaction import ReactionContainer
    from chython.formats.ctfile import rxn

    r = rxn(_RXN_TEXT)
    assert isinstance(r, ReactionContainer)
    assert len(r.reactants) == 1 and len(r.products) == 1
    assert r.title == 'oxidation'


def test_rxn_exports_v2000_by_default_and_round_trips():
    from chython.formats.ctfile import rxn

    original = rxn(_RXN_TEXT)
    text = rxn(original)
    assert text.split('\n')[0] == '$RXN'
    again = rxn(text)
    assert len(again.reactants) == len(original.reactants)
    assert [m.atom_count for m in again.molecules()] == [m.atom_count for m in original.molecules()]


def test_rxn_round_trips_v3000():
    from chython.formats.ctfile import rxn

    original = rxn(_RXN_TEXT)
    text = rxn(original, version=3000)
    assert text.split('\n')[0] == '$RXN V3000'
    again = rxn(text)
    assert [m.atom_count for m in again.molecules()] == [m.atom_count for m in original.molecules()]


def test_rxn_auto_escalates_when_any_component_needs_it():
    from chython.formats.ctfile import rxn

    r = rxn(_RXN_TEXT)
    target = r.products[0]
    with target.edit() as e:
        # 100000.0 overflows the V2000 10-character coordinate column
        e.set_xy(next(iter(target.atom_numbers)), 100000.0, 0.0)
    log = []
    text = rxn(r, log=log)
    assert text.split('\n')[0] == '$RXN V3000'
    assert any('coordinate' in x for x in log), log


def test_agents_alone_do_not_escalate():
    """An agent is expressible in V2000 by the third-count convention, so it is not a reason."""
    from chython.formats.ctfile import rxn

    # replace the first (and only) '  1  1' counts line with '  1  1  1' to add an agent count
    text = _RXN_TEXT.replace('  1  1', '  1  1  1', 1) + '\n$MOL\n' + _ETHANOL_V2000
    r = rxn(text)
    assert len(r.agents) == 1
    log = []
    out = rxn(r, log=log)
    assert out.split('\n')[0] == '$RXN', out.split('\n')[0]
    assert any('third count' in x for x in log), log


def test_mol_and_rxn_refuse_each_other_by_naming_the_other():
    from pytest import raises

    from chython.formats.ctfile import mol, rxn

    with raises(TypeError, match='rxn'):
        mol(rxn(_RXN_TEXT))
    with raises(TypeError, match='mol'):
        rxn(mol(_ETHANOL_V2000))


# $RFMT leader handling: clean, RFMT-led, and no $RXN at all

def test_rxn_accepts_clean_rxn_input():
    """A string that starts directly with $RXN has no leading lines to skip."""
    from chython.formats.ctfile import rxn

    log = []
    r = rxn(_RXN_TEXT, log=log)
    assert r.title == 'oxidation'
    assert not any('skipped' in x for x in log), log


def test_rxn_strips_rfmt_header_and_logs():
    """Input lifted from an RDfile may start with $RFMT and a datestamp before $RXN; rxn() slices
    from the first $RXN line and logs how many lines it skipped."""
    from chython.formats.ctfile import rxn

    # Real RDfiles use "$RFMT" (the record marker) and "$DATUM" (data line) — neither starts "$RXN"
    rfmt_prefix = '$RFMT\n$DATUM reaction_id 42\n'
    rfmt_led = rfmt_prefix + _RXN_TEXT
    log = []
    r = rxn(rfmt_led, log=log)
    assert r.title == 'oxidation'
    assert any('skipped' in x and '$RXN' in x for x in log), log


def test_rxn_raises_when_no_rxn_line():
    """A string with no $RXN line at all is unparseable and raises MalformedCtfile."""
    from chython.formats.ctfile import MalformedCtfile, rxn

    with raises(MalformedCtfile, match=r'\$RXN'):
        rxn('just some garbage\nno rxn header here\n')


def test_needs_v3000_for_too_many_bonds():
    """V2000 can express at most 999 bonds in the 3-character count field.

    K45 has 45*44/2 = 990 bonds (fits); K46 has 46*45/2 = 1035 bonds (does not).
    """
    m_under = MoleculeContainer()
    with m_under.edit() as e:
        ids = [e.add_atom(6) for _ in range(45)]
        for i in range(45):
            for j in range(i + 1, 45):
                e.add_bond(ids[i], ids[j], 1)
    assert needs_v3000(m_under) is None  # 990 bonds fits

    m_over = MoleculeContainer()
    with m_over.edit() as e:
        ids = [e.add_atom(6) for _ in range(46)]
        for i in range(46):
            for j in range(i + 1, 46):
                e.add_bond(ids[i], ids[j], 1)
    reason = needs_v3000(m_over)
    assert reason is not None and 'bond' in reason, reason  # 1035 bonds does not


def test_mol_writes_no_data_fields_when_told_none():
    """`None` is "I did not say" and `{}` is "I said none", as for title and sgroups."""
    m = mol(_ETHANOL_V2000)
    m.meta['ACTIVITY'] = '5.0'
    assert 'ACTIVITY' not in mol(m, meta={})


def test_mol_writes_the_fields_it_is_given_instead_of_the_molecules_own():
    m = mol(_ETHANOL_V2000)
    m.meta['ACTIVITY'] = '5.0'
    out = mol(m, meta={'SOURCE': 'plan'})
    assert '>  <SOURCE>' in out
    assert 'ACTIVITY' not in out


def test_mol_says_where_the_data_fields_went():
    """A molfile has no data-field section; SD framing with no `$$$$` is worth one line."""
    m = mol(_ETHANOL_V2000)
    m.meta['ACTIVITY'] = '5.0'
    log = []
    out = mol(m, log=log)
    assert '>  <ACTIVITY>' in out and '$$$$' not in out
    assert any('after M  END' in x for x in log), log
