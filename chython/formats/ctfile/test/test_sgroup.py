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
"""S-groups: what a CTfile says about a molecule that is not its constitution.

Fixtures are hand-written because the corpus carries only ``DAT``/``MRV_IMPLICIT_H`` groups.  The
invariant throughout: a bond reference is an endpoint pair, never a bond number.  A CTfile bond
number is a position in the bond block, so it does not survive an edit or a reordered write.
"""

from pytest import raises

from .._errors import UnsupportedCtfile
from .._sgroup import (DISP_MAX, NO_INDEX, SGroup, SGroupStore, format_fielddisp, parse_fielddisp)
from .._v2000 import emit_v2000, parse_v2000


# butane, so there are four atoms and three bonds to refer to.
_BUTANE = ['butane', '  test', '', '  4  3  0  0  0  0            999 V2000',
           '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
           '    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
           '    2.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
           '    3.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
           '  1  2  1  0  0  0  0',
           '  2  3  1  0  0  0  0',
           '  3  4  1  0  0  0  0']


def _record(*properties):
    return _BUTANE + list(properties) + ['M  END']


def _read(*properties):
    ctab = parse_v2000(_record(*properties), [])
    mol, store, log = ctab.build()
    return mol, store, log


def _shape(mol, store):
    """A store in position terms: what a round trip must preserve.

    Positions rather than stable ids or file numbers -- neither the S-group number nor a bond index
    is promised to survive, only which atoms and bonds each record refers to and what it says.
    """
    position = {sid: i for i, sid in enumerate(mol.atom_numbers)}
    out = []
    for r in store.records:
        out.append((
            r.type, r.subtype, r.name,
            tuple(sorted(position[a] for a in r.atoms if a in position)),
            tuple(sorted(position[a] for a in r.patoms if a in position)),
            tuple(sorted(tuple(sorted((position[a], position[b])))
                         for a, b in r.bonds if a in position and b in position)),
            tuple(r.data), r.disp,
            tuple(sorted((k, tuple(v)) for k, v in r.fields.items())),
            tuple(sorted((None if p is None else tuple(sorted((position[p[0]], position[p[1]]))),
                          t) for p, t in r.cstates)),
        ))
    return sorted(out)


def _round_trip(*properties):
    """Read, write, read.  Returns the two shapes and the two stores."""
    mol, store, _ = _read(*properties)
    lines, log = emit_v2000(mol, store)
    ctab2 = parse_v2000(lines, [])
    mol2, store2, _ = ctab2.build()
    return _shape(mol, store), _shape(mol2, store2), store2, lines


# the DAT group

def test_a_data_group_round_trips_with_its_field_name_and_datum():
    # The field name occupies 30 fixed columns and the type the 2 after it.
    sdt = 'M  SDT   1 ' + f'{"BOILING.POINT":<30s}' + 'N'
    before, after, store, _ = _round_trip('M  STY  1   1 DAT', 'M  SAL   1  2   1   2', sdt,
                                          'M  SED   1 -0.5C')
    assert before == after
    assert [r.name for r in store.records] == ['BOILING.POINT']
    assert store.records[0].fields['FIELDTYPE'] == ['N']
    assert store.records[0].field_data == '-0.5C'


def test_a_datum_longer_than_one_line_is_reassembled_and_re_split():
    """``M  SCD`` continues a datum and ``M  SED`` closes it, at 69 characters a line."""
    text = 'A' * 80
    before, after, store, lines = _round_trip('M  STY  1   1 DAT', 'M  SAL   1  1   1',
                                              'M  SDT   1 LONG',
                                              f'M  SCD   1 {text[:69]}',
                                              f'M  SED   1 {text[69:]}')
    assert store.records[0].field_data == text
    assert before == after
    assert sum(1 for x in lines if x.startswith('M  SCD')) == 1, lines


def test_an_unclosed_scd_datum_is_kept_and_reported():
    """The datum is real; only its terminator is missing."""
    mol, store, log = _read('M  STY  1   1 DAT', 'M  SAL   1  1   1', 'M  SDT   1 X',
                            'M  SCD   1 value')
    assert store.records[0].field_data == 'value'
    assert any('not closed' in x for x in log), log


def test_the_display_position_survives_a_version_neutral_round_trip():
    """V2000 ``M  SDD`` and V3000 ``FIELDDISP=`` are byte-for-byte the same layout, so one model
    serves both and a DAT group keeps its anchor across a version change."""
    before, after, store, _ = _round_trip('M  STY  1   1 DAT', 'M  SAL   1  1   1',
                                          'M  SDT   1 X',
                                          'M  SDD   1     1.2300   -4.5600    DAU   ALL  0       0',
                                          'M  SED   1 v')
    assert before == after
    assert store.records[0].disp is not None
    assert store.records[0].disp[:2] == (1.23, -4.56)


def test_an_unparseable_display_line_is_kept_as_text_rather_than_dropped():
    mol, store, log = _read('M  STY  1   1 DAT', 'M  SAL   1  1   1', 'M  SDT   1 X',
                            'M  SDD   1 nonsense', 'M  SED   1 v')
    assert store.records[0].disp is None
    assert store.records[0].fields.get('FIELDDISP') == ['nonsense']


def test_fielddisp_formatting_is_its_own_inverse():
    for text in ('    1.2300   -4.5600    DAU   ALL  0       0',
                 '    0.0000    0.0000    DA    ALL  1       1',
                 # Both ends of the field: the sign costs a column, so ten characters reach
                 # `-9999.9999` at the low end and `99999.9999` at the high one.
                 '99999.9999-9999.9999    DA    ALL  1       1'):
        parsed = parse_fielddisp(text, [])
        assert parsed is not None, text
        assert parse_fielddisp(format_fielddisp(parsed), []) == parsed, text


def test_a_display_anchor_too_wide_for_the_field_is_kept_as_text_and_not_as_a_number():
    """``1e9`` is a number ``float()`` accepts and F10.4 cannot write, so it is refused at the door
    and lands in the same verbatim fallback as any other unparseable anchor."""
    log = []
    assert parse_fielddisp('       1e9       1e9    DA    ALL  1       5', log) is None
    assert any('outside F10.4' in x for x in log), log
    mol, store, _ = _read('M  STY  1   1 DAT', 'M  SAL   1  1   1', 'M  SDT   1 X',
                          'M  SDD   1        1e9       1e9    DA', 'M  SED   1 v')
    assert store.records[0].disp is None
    assert store.records[0].fields.get('FIELDDISP') == ['       1e9       1e9    DA']


def test_a_display_anchor_too_wide_for_the_field_is_clamped_rather_than_shifting_the_columns():
    """The failure this prevents is a wrong *y*, not a wrong x.

    ``f'{1e9:10.4f}'`` is eleven characters, so every column after it starts one early and
    ``1e9, 1e9`` comes back as ``1e9, 1e-05``.  Clamping keeps the layout and says so.
    """
    log = []
    text = format_fielddisp((1e9, 1e9, '    DA    ALL  1       5'), log)
    assert any('clamped' in x for x in log), log
    assert len(text[:20]) == 20 and text[20:] == '    DA    ALL  1       5', repr(text)
    assert parse_fielddisp(text, []) == (DISP_MAX, DISP_MAX, '    DA    ALL  1       5')

    # The unguarded format read back by column: the y field is nine of x's spilled digits.
    unguarded = f'{1e9:10.4f}{1e9:10.4f}'
    assert (float(unguarded[0:10]), float(unguarded[10:20])) == (1e9, 1e-05), \
        'the y coordinate the guard exists to protect'


# the types the corpus lacks

def test_a_superatom_survives_with_its_abbreviation():
    """``SUP`` is how a file says "these four atoms are drawn as Ph"."""
    before, after, store, _ = _round_trip('M  STY  1   1 SUP', 'M  SAL   1  2   1   2',
                                          'M  SMT   1 Et')
    assert before == after
    assert store.records[0].type == 'SUP'
    assert store.records[0].fields['LABEL'] == ['Et']


def test_a_multiple_group_keeps_its_parent_atom_subset():
    """``MUL`` says "this fragment repeats n times"; ``SPA`` says which copy is drawn.  Without
    ``SPA`` every copy becomes structural."""
    before, after, store, _ = _round_trip('M  STY  1   1 MUL', 'M  SAL   1  4   1   2   3   4',
                                          'M  SPA   1  2   1   2', 'M  SMT   1 2')
    assert before == after
    assert store.records[0].type == 'MUL'
    assert len(store.records[0].patoms) == 2
    assert store.records[0].fields['MULT'] == ['2']


def test_a_repeat_unit_keeps_the_bonds_it_is_cut_at():
    """``SRU`` with ``SBL``: the crossing bonds are the content of a polymer repeat unit, and they
    are bond *numbers* in the file."""
    before, after, store, _ = _round_trip('M  STY  1   1 SRU', 'M  SAL   1  2   2   3',
                                          'M  SBL   1  2   1   3', 'M  SST  1   1 HT')
    assert before == after
    assert store.records[0].subtype == 'HT'
    assert len(store.records[0].bonds) == 2


def test_a_bond_reference_is_an_endpoint_pair_and_not_a_bond_number():
    """The bond block is rewritten in the molecule's own order, so the numbers out are not the
    numbers in; the *bonds* referred to must still be the same ones."""
    mol, store, _ = _read('M  STY  1   1 SRU', 'M  SAL   1  2   2   3', 'M  SBL   1  1   1')
    position = {sid: i for i, sid in enumerate(mol.atom_numbers)}
    pairs = {tuple(sorted((position[a], position[b]))) for a, b in store.records[0].bonds}
    assert pairs == {(0, 1)}, 'M  SBL 1 is the first bond line, which joins atoms 1 and 2'
    lines, _ = emit_v2000(mol, store)
    store2 = parse_v2000(lines, []).build()[1]
    mol2 = parse_v2000(lines, []).build()[0]
    position2 = {sid: i for i, sid in enumerate(mol2.atom_numbers)}
    pairs2 = {tuple(sorted((position2[a], position2[b]))) for a, b in store2.records[0].bonds}
    assert pairs2 == pairs


def test_a_cstate_vector_travels_with_its_bond_and_not_with_its_number():
    """``M  SBV``'s first value is a bond number, so keeping the line as text after a renumbering
    points the vector -- a direction in the drawing -- at a different bond."""
    before, after, store, lines = _round_trip('M  STY  1   1 SRU', 'M  SAL   1  2   2   3',
                                              'M  SBV   1   3    1.0000    0.0000')
    assert before == after
    assert len(store.records[0].cstates) == 1
    assert any(x.startswith('M  SBV') for x in lines), lines


def test_a_bond_the_file_declared_backwards_is_still_found():
    """A file may write a bond as ``2  1``, so the stored pair is ``(2, 1)`` while the molecule
    reports ``(1, 2)``.  The file's order is not recoverable, so the writer's bond-number table holds
    both orientations; ``_shape`` normalises with ``sorted`` and is blind to this.  Asserted as "the
    reference survived", since the failure is a dropped reference the writer only reports.
    """
    butane = list(_BUTANE)
    # By content, not offset: atom and bond blocks are both fixed-width, so an off-by-one edits an
    # atom instead.
    butane[butane.index('  1  2  1  0  0  0  0')] = '  2  1  1  0  0  0  0'
    record = butane + ['M  STY  1   1 SRU', 'M  SAL   1  2   1   2', 'M  SBL   1  1   1',
                       'M  SBV   1   1    0.5000    0.5000', 'M  END']
    mol, store, _ = parse_v2000(record, []).build()
    assert store.records[0].bonds == [(2, 1)], 'stored in the order the FILE gave, which is the point'
    assert [(b.n, b.m) for b in mol.bonds()] == [(1, 2), (2, 3), (3, 4)], \
        'while the molecule reports the other order'
    lines, log = emit_v2000(mol, store)
    assert not any('dropped' in x for x in log), log
    assert any(x.startswith('M  SBL') for x in lines), lines
    assert any(x.startswith('M  SBV') for x in lines), lines


def test_a_bond_reference_whose_bond_was_deleted_is_dropped_and_reported():
    """An S-group loses a bond reference only when the bond goes while its atoms stay.  ``M  SBL``
    and ``M  SBV`` resolve on separate lines, so one record naming the same bond twice drives both
    failure branches.  The V3000 twin is in ``test_v3000.py``."""
    mol, store, _ = _read('M  STY  1   1 SRU', 'M  SAL   1  2   1   2', 'M  SBL   1  1   1',
                          'M  SBV   1   1    0.5000    0.5000')
    assert store.records[0].bonds and store.records[0].cstates, 'the fixture lost its references early'
    with mol.edit() as e:
        e.delete_bond(*store.records[0].bonds[0])
    lines, log = emit_v2000(mol, store)
    assert sum('no longer' in x for x in log) == 2, log
    assert not any(x.startswith('M  SBL') or x.startswith('M  SBV') for x in lines), lines
    assert any(x.startswith('M  SAL') for x in lines), 'the atoms are still referenced'


def test_an_atom_alias_is_read_from_the_line_after_its_header():
    """``A  aaa`` then free text.  The text may look like a property line, so the reader takes the
    next line whole rather than pattern-matching it."""
    mol, store, _ = _read('A    2', 'M  CHG  1   1   0')
    assert store.aliases and list(store.aliases.values()) == ['M  CHG  1   1   0']


def test_an_atom_value_line_is_stored_the_same_way_as_an_alias():
    mol, store, _ = _read('V    2 some text')
    assert list(store.aliases.values()) == ['some text']


def test_an_alias_for_an_atom_that_does_not_exist_is_reported():
    mol, store, log = _read('A    9', 'text')
    assert any('out of 1..4' in x for x in log), log


# numbering and recovery

def test_an_sgroup_numbered_zero_is_a_real_sgroup():
    """The "no number" sentinel is 0xFFFF and not 0: ``M  STY``'s 3-character field can hold
    ``  0``, and nothing in CTfile forbids it."""
    mol, store, _ = _read('M  STY  1   0 DAT', 'M  SAL   0  1   1', 'M  SDT   0 X',
                          'M  SED   0 v')
    assert [r.index for r in store.records] == [0]
    assert store.records[0].field_data == 'v'


def test_a_group_used_before_it_is_declared_is_created_and_reported():
    """The lines may arrive in any order, including data ahead of the declaring ``M  STY``."""
    mol, store, log = _read('M  SAL   1  1   1', 'M  SDT   1 X', 'M  SED   1 v',
                            'M  STY  1   1 DAT')
    assert [r.type for r in store.records] == ['DAT']
    assert any('used before' in x for x in log), log


def test_a_group_never_declared_at_all_is_read_as_gen():
    mol, store, log = _read('M  SAL   1  1   1', 'M  SDT   1 X', 'M  SED   1 v')
    assert [r.type for r in store.records] == ['GEN']
    assert any('read as GEN' in x for x in log), log


def test_an_out_of_range_atom_reference_is_dropped_and_the_record_kept():
    """The record still says a field was attached to something."""
    mol, store, log = _read('M  STY  1   1 DAT', 'M  SAL   1  2   1   9', 'M  SDT   1 X',
                            'M  SED   1 v')
    assert len(store.records) == 1
    assert len(store.records[0].atoms) == 1
    assert any('out of 1..4' in x for x in log), log


def test_an_out_of_range_bond_reference_is_dropped_and_the_record_kept():
    mol, store, log = _read('M  STY  1   1 SRU', 'M  SAL   1  1   1', 'M  SBL   1  1   9')
    assert len(store.records) == 1 and not store.records[0].bonds
    assert any('out of 1..3' in x for x in log), log


def test_a_record_whose_atoms_all_vanish_stays_present():
    """``translate`` is where a molecule that lost atoms meets a record referring to them.  An
    empty DAT record still states that a field was attached."""
    sg = SGroup('DAT', index=1)
    sg.atoms.extend([10, 11])
    sg.name = 'X'
    log = []
    out = sg.translate({}, log)
    assert out.name == 'X' and out.atoms == []
    assert any('reference(s) dropped' in x for x in log), log


def test_an_unmodelled_keyword_is_named_in_the_log_and_not_swept_up():
    """A keyword this release does not model is a known gap; an unrecognised one is a surprise, and
    the log must distinguish them."""
    mol, store, log = _read('M  SDS EXP  1   1')
    assert any('M  SDS is not modelled' in x for x in log), log
    mol, store, log = _read('M  ZZZ nonsense')
    assert any('unrecognised property' in x for x in log), log


def test_a_query_atom_list_is_refused_by_name():
    """``M  ALS`` is a query, not a molecule, and the refusal names the reader to reach for."""
    with raises(UnsupportedCtfile, match='query reader'):
        parse_v2000(_record('M  ALS   1  2 F     F   Cl'), [])


def test_rgp_on_a_non_marker_is_logged_not_fatal():
    # `M  RGP` on a carbon has no destination: the assignment is dropped and logged.
    _, _, log = _read('M  RGP  1   1   1')
    assert any('RGP' in x for x in log), log


# the store

def test_the_store_finds_records_by_field_name():
    mol, store, _ = _read('M  STY  2   1 DAT   2 DAT', 'M  SAL   1  1   1', 'M  SDT   1 A',
                          'M  SED   1 1', 'M  SAL   2  1   2', 'M  SDT   2 B', 'M  SED   2 2')
    assert [r.field_data for r in store.by_name('B')] == ['2']
    assert len(store.data_records()) == 2


def test_next_index_does_not_collide_with_a_number_already_in_use():
    store = SGroupStore([SGroup('DAT', index=1), SGroup('DAT', index=3)])
    assert store.next_index() not in (1, 3)


def test_a_record_with_no_index_is_numbered_on_write_and_not_dropped():
    mol, store, _ = _read()
    sg = SGroup('DAT', index=NO_INDEX)
    sg.atoms.append(next(iter(mol.atom_numbers)))
    sg.name = 'X'
    sg.data.append(b'v')
    lines, _ = emit_v2000(mol, SGroupStore([sg]))
    again = parse_v2000(lines, []).build()[1]
    assert [r.name for r in again.records] == ['X']
    assert again.records[0].index != NO_INDEX, 'the writer has to invent a number to refer to it by'


def test_an_invented_number_is_not_one_another_record_already_states():
    """An invented number drawn from the record's position collides with the file's own 1..n
    numbers, and then two groups' ``M  SAL``/``M  SDT`` lines name the same group."""
    mol, store, _ = _read()
    first = next(iter(mol.atom_numbers))
    unnumbered, numbered = SGroup('DAT', index=NO_INDEX), SGroup('DAT', index=1)
    for sg, name in ((unnumbered, 'U'), (numbered, 'N')):
        sg.atoms.append(first)
        sg.name = name
        sg.data.append(b'v')
    lines, _ = emit_v2000(mol, SGroupStore([unnumbered, numbered]))
    again = parse_v2000(lines, []).build()[1]
    assert sorted(r.name for r in again.records) == ['N', 'U'], 'two records, not one shared number'


def test_a_number_too_wide_for_v2000_is_rewritten_and_its_parent_follows_it():
    """The model holds V3000's 65534 and V2000 writes three columns, so reading V3000 and writing
    V2000 reaches the gap.  A wide number would push every column right, so the record is
    renumbered; ``M  SPL`` names a number, so it goes through the same table or the hierarchy
    reattaches elsewhere."""
    mol, store, _ = _read()
    first = next(iter(mol.atom_numbers))
    parent, child = SGroup('SUP', index=1500), SGroup('SUP', index=7)
    child.parent = 1500
    for sg in (parent, child):
        sg.atoms.append(first)
    lines, log = emit_v2000(mol, SGroupStore([parent, child]))
    assert any('does not fit' in x and '1500' in x for x in log), log
    assert not any('1500' in x for x in lines), 'the wide number reached the file'
    again = parse_v2000(lines, []).build()[1]
    written = {r.index for r in again.records}
    assert 1500 not in written and len(written) == 2, written
    child_again = next(r for r in again.records if r.parent != NO_INDEX)
    assert child_again.parent in written - {child_again.index}, \
        'the parent reference points at the group it always pointed at, under its new number'


def test_a_cstate_whose_bond_died_keeps_its_vector_text_instead_of_vanishing():
    """A CSTATE's vector tail lives in a run whose length is derived from the cstate count, so
    compacting a dead entry out re-pairs every surviving vector with the entry before it.  Demoting
    to "unresolved" keeps the count and the file's vector text."""
    sg = SGroup('SRU', index=1)
    sg.atoms = [1, 2]
    sg.cstates = [((1, 2), '1.0 0.0 0.0'), ((2, 3), '0.0 1.0 0.0')]
    log = []
    out = sg.translate({1: 11, 2: 12}, log)
    assert out.cstates == [((11, 12), '1.0 0.0 0.0'), (None, '0.0 1.0 0.0')]
    assert any('dropped' in x for x in log), log
