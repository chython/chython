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
"""V2000 reading and writing, against chython 2 and against the corpus.

The oracle test is the centre of the file: two independent implementations reading the same bytes and
reporting the same molecule.  What the oracle cannot check is a record chython 2 refuses, so the
recoveries have their own tests below with the malformation written out inline.
"""

from pytest import raises

from .._errors import MalformedCtfile, UnsupportedCtfile
from .._hydrogens import implicit_for_atom
from .._sdf import sniff_version
from .._sgroup import SGroupStore
from .._v2000 import V2000_STAMP, emit_v2000, parse_v2000
from .conftest import holds_an_aromatic_bond


def _v2000_only(records):
    """``(index, record)`` for the V2000 records of a file.

    A ``.sdf`` may mix versions -- ``test/implicit.sdf`` opens with a V3000 record.  The index is
    kept so the oracle stays aligned with chython 2's record list.
    """
    return [(i, r) for i, r in enumerate(records) if sniff_version(r, []) == V2000_STAMP]


def _state(mol):
    """Everything about a molecule a round trip must preserve, in position terms."""
    sids = list(mol.atom_numbers)
    index = {s: i for i, s in enumerate(sids)}
    atoms = [(mol.element_of(s), mol.charge_of(s), mol.isotope_of(s), bool(mol.radical_of(s)),
              mol.implicit_h_of(s)) for s in sids]
    bonds = sorted((min(index[b.n], index[b.m]), max(index[b.n], index[b.m]), b.order)
                   for b in mol.bonds())
    return atoms, bonds, [mol.parity_of(s) for s in sids]


def _read(lines):
    ctab = parse_v2000(lines, [])
    return ctab.build()


# the corpus sweep and the oracle

def test_every_record_is_read_or_refused_but_never_crashes(corpus):
    """No record in the corpus raises anything other than a ``CtfileError``.

    A `CtfileError` is this package saying "no, and here is why"; any other exception is a bug -- an
    `IndexError` on a short line, a `KeyError` on a missing atom.
    """
    read = refused = 0
    for name, records in corpus.items():
        for n, record in _v2000_only(records):
            try:
                parse_v2000(record, []).build()
                read += 1
            except UnsupportedCtfile:
                refused += 1
            except MalformedCtfile as e:  # allowed, but only with a message worth reading
                assert str(e), f'{name} record {n}: refusal with no message'
                refused += 1
    assert read > 400, f'only {read} records read; the corpus should yield hundreds'
    assert read + refused > 500


def test_agrees_with_chython2_wherever_chython2_commits(corpus, v2_molecules):
    """Constitution identical to chython 2's, and hydrogen counts identical wherever V2 gives one.

    Honouring V2000 ``vvv`` over the valence rules changes a ferrocene Cp carbon's hydrogen count
    from 1 to 2, and this test fails on it.
    """
    compared = 0
    for name, records in corpus.items():
        v2 = v2_molecules.get(name)
        if v2 is None or len(v2) != len(records):
            continue
        for n, record in _v2000_only(records):
            m2 = v2[n]
            try:
                mol, _, _ = parse_v2000(record, []).build()
            except UnsupportedCtfile:
                continue
            sids = list(mol.atom_numbers)
            nums = list(m2)
            where = f'{name} record {n}'
            assert len(sids) == len(nums), f'{where}: atom count {len(sids)} vs {len(nums)}'
            for sid, num in zip(sids, nums):
                a2, a3 = m2.atom(num), mol.atom(sid)
                assert a3.element == a2.atomic_number, f'{where} atom {num}: element'
                assert a3.charge == a2.charge, f'{where} atom {num}: charge'
                assert bool(a3.is_radical) == bool(a2.is_radical), f'{where} atom {num}: radical'
                assert (a3.isotope or 0) == (a2.isotope or 0), f'{where} atom {num}: isotope'
                if a2.implicit_hydrogens is not None:
                    assert mol.implicit_h_of(sid) == a2.implicit_hydrogens, \
                        f'{where} atom {num}: implicit H'
            ours = {frozenset((sids.index(b.n), sids.index(b.m))): b.order for b in mol.bonds()}
            theirs = {frozenset((nums.index(u), nums.index(v))): b.order for u, v, b in m2.bonds()}
            assert ours == theirs, f'{where}: bonds'
            compared += 1
    assert compared > 400, f'only {compared} records compared against chython 2'


def test_where_chython2_says_unknown_we_say_unknown_or_it_was_derivable(corpus, v2_molecules):
    """An atom V2 reports as unknown must be one of three things here:

    * flagged unknown here too;
    * determined by something the file stated -- a hydrogen count, or a total valence equal to the
      orders drawn, which says nothing is left over for hydrogen;
    * holding an aromatic bond, where the core's classifier settles the count because every Kekule
      form of the ring gives the same one.  V2 committed only for a neutral aromatic carbon.

    A non-aromatic recovery would mean the two valence tables had drifted, and is still caught.
    """
    checked = 0
    for name, records in corpus.items():
        v2 = v2_molecules.get(name)
        if v2 is None or len(v2) != len(records):
            continue
        for n, record in _v2000_only(records):
            m2 = v2[n]
            try:
                ctab = parse_v2000(record, [])
                mol, _, _ = ctab.build()
            except UnsupportedCtfile:
                continue
            unknown = set(ctab.unknown_hydrogens)
            stated = {i for i, a in enumerate(ctab.atoms)
                      if a.stated_h is not None or a.valence is not None}
            for i, (sid, num) in enumerate(zip(mol.atom_numbers, list(m2))):
                if m2.atom(num).implicit_hydrogens is not None:
                    continue
                checked += 1
                aromatic = holds_an_aromatic_bond(mol, sid)
                assert sid in unknown or i in stated or aromatic, (
                    f'{name} record {n} atom {num}: chython 2 says the hydrogen count is unknown, '
                    f'but this reader committed to {mol.implicit_h_of(sid)} with nothing stated in '
                    f'the file and no aromatic bond to settle it')
    assert checked, 'the corpus contains no unknown hydrogen counts; this test proved nothing'


def test_read_write_read_is_identical(corpus):
    """The round trip preserves constitution, hydrogen counts, parities and the unknown set.

    An aromatic record is round-tripped twice, as read and through ``kekule()``, and the unknown set
    is asserted differently for each: written as read it must come back **exactly** the same, since an
    atom holding ``H_UNKNOWN`` gets neither a valence nor an ``MRV_IMPLICIT_H`` group; kekulised first
    it may come back **smaller** but never larger, an unknown appearing meaning a write path invented
    a doubt.
    """
    checked = aromatic = repaired = 0
    for name, records in corpus.items():
        for n, record in _v2000_only(records):
            ctab = parse_v2000(record, [])
            mol, store, _ = ctab.build()
            kekulised = False
            if mol.aromatic_bond_count:
                aromatic += 1
                # Written as read first: the bond block states type 4 back, so the exact invariant
                # applies to an aromatic record as it does to a Kekule one.
                as_read = parse_v2000(emit_v2000(mol, store, title=ctab.title)[0], [])
                same, _, _ = as_read.build()
                assert _state(same) == _state(mol), \
                    f'{name} record {n}: writing the aromatic representation changed the molecule'
                assert set(as_read.unknown_hydrogens) == set(ctab.unknown_hydrogens), \
                    f'{name} record {n}: writing as read changed which counts are unknown'
                result = mol.kekule()
                if result.unresolved:
                    continue  # no Kekule form exists; there is nothing more to test on this record
                kekulised = True
            lines, _ = emit_v2000(mol, store, title=ctab.title)
            ctab2 = parse_v2000(lines, [])
            mol2, _, _ = ctab2.build()
            where = f'{name} record {n}'
            before, after = set(ctab.unknown_hydrogens), set(ctab2.unknown_hydrogens)
            if kekulised:
                assert after <= before, \
                    f'{where}: kekulising and writing invented an unknown count: {after - before}'
                if after < before:
                    repaired += 1
                    continue  # the resolved count is a real number now and need not equal the 0
            else:
                assert before == after, \
                    f'{where}: round trip changed which hydrogen counts are unknown'
            assert _state(mol) == _state(mol2), f'{where}: round trip changed the molecule'
            checked += 1
    assert checked > 400
    assert aromatic, ('no aromatic record reached the writer; either the corpus changed or the '
                      'reader stopped storing bond type 4, and this test stopped covering it')
    assert repaired, ('no record had an unknown count resolved by kekulising, so the branch above is '
                      'untested; the corpus used to contain such records')


# the details

_ETHANOL = """ethanol
  test
comment
  3  2  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END""".split('\n')


def test_the_baseline_record_reads_as_ethanol():
    mol, _, log = _read(_ETHANOL)
    assert [mol.element_of(s) for s in mol.atom_numbers] == [6, 6, 8]
    assert [mol.implicit_h_of(s) for s in mol.atom_numbers] == [3, 2, 1]
    assert not log, f'a clean record should log nothing, got {log}'


def test_counts_line_is_a_hint_not_a_contract():
    """A record whose counts line over-promises still yields the atoms it does contain."""
    lines = list(_ETHANOL)
    lines[3] = '  9  8  0  0  0  0            999 V2000'
    mol, _, log = _read(lines)
    assert len([*mol.atom_numbers]) == 3
    assert any('atom' in x for x in log), f'the shortfall must be logged, got {log}'


def test_a_bond_to_a_nonexistent_atom_is_dropped_with_a_log():
    lines = list(_ETHANOL)
    lines[3] = '  3  3  0  0  0  0            999 V2000'
    lines.insert(9, '  1  9  1  0  0  0  0')
    mol, _, log = _read(lines)
    assert len(list(mol.bonds())) == 2
    assert any('atom number out of' in x and 'dropped' in x for x in log), log


def test_query_bond_types_are_refused_by_name():
    """Bond types 5, 6 and 7 are "single or double" and friends: a constraint, not a bond, which a
    molecule cannot hold.  The refusal names the type so the caller reaches for a query reader.
    """
    for order in (5, 6, 7):
        lines = list(_ETHANOL)
        lines[7] = f'  1  2{order:3d}  0  0  0  0'
        with raises(UnsupportedCtfile, match='query bond'):
            _read(lines)


def test_the_hhh_query_field_is_ignored_with_a_log():
    """``hhh`` states a minimum, not a count."""
    lines = list(_ETHANOL)
    lines[4] = '    0.0000    0.0000    0.0000 C   0  0  0  2  0  0  0  0  0  0  0  0'
    mol, _, log = _read(lines)
    assert mol.implicit_h_of(next(iter(mol.atom_numbers))) == 3
    assert any('minimum' in x for x in log), log


def test_a_charge_beyond_the_ccc_field_comes_from_m_chg():
    """The properties block wins, per the spec, and is the only way a big charge survives."""
    lines = list(_ETHANOL)
    lines.insert(-1, 'M  CHG  1   3  -1')
    mol, _, _ = _read(lines)
    assert [mol.charge_of(s) for s in mol.atom_numbers] == [0, 0, -1]


def test_m_chg_supersedes_the_atom_block_rather_than_adding_to_it():
    lines = list(_ETHANOL)
    lines[6] = '    2.0000    0.0000    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0'
    lines.insert(-1, 'M  CHG  1   3  -1')
    mol, _, _ = _read(lines)
    assert [mol.charge_of(s) for s in mol.atom_numbers] == [0, 0, -1]


def test_m_iso_supersedes_the_mass_difference_column():
    lines = list(_ETHANOL)
    lines.insert(-1, 'M  ISO  1   3  18')
    mol, _, _ = _read(lines)
    assert [mol.isotope_of(s) for s in mol.atom_numbers] == [0, 0, 18]


def test_a_truncated_record_without_m_end_still_reads():
    mol, _, log = _read(_ETHANOL[:-1])
    assert len([*mol.atom_numbers]) == 3
    assert any('M  END' in x for x in log), log


def test_a_free_text_symbol_is_read_as_an_alias_and_not_as_a_lost_record():
    """``Me`` in the symbol column: the reference lists no free text there and readers differ over
    what it means, so the record is kept and the text goes where a display label goes -- the atom's
    alias.  The atom itself is the marker, so nothing reads a borrowed element as chemistry: element
    0 holds no hydrogens, contributes nothing to a formula and has no valence to violate.
    """
    lines = list(_ETHANOL)
    lines[4] = '    0.0000    0.0000    0.0000 Me  0  0  0  0  0  0  0  0  0  0  0  0'
    ctab = parse_v2000(lines, [])
    mol, store, log = ctab.build()
    sid, *_ = mol.atom_numbers
    assert len([*mol.atom_numbers]) == 3, 'the whole record, not two atoms and a renumbering'
    assert len(list(mol.bonds())) == 2
    assert mol.aliases == {sid: b'Me'}
    assert store.aliases == {sid: 'Me'}
    assert mol.element_of(sid) == 0, 'the marker, not a borrowed element'
    assert mol.implicit_h_of(sid) == 0, 'a marker holds no hydrogens'
    assert not ctab.unknown_hydrogens, 'nothing about a marker is undeterminable'
    assert any('names no element' in x for x in log), log


def test_an_explicit_alias_line_outranks_the_symbol_column():
    """``A  1`` states the label outright, and a label read out of the symbol column is a recovery, so
    the statement wins."""
    lines = list(_ETHANOL)
    lines[4] = '    0.0000    0.0000    0.0000 Me  0  0  0  0  0  0  0  0  0  0  0  0'
    lines.insert(-1, 'A  1')
    lines.insert(-1, 'Ethyl')
    mol, _, _ = _read(lines)
    sid, *_ = mol.atom_numbers
    assert mol.aliases == {sid: b'Ethyl'}


def test_a_g_line_is_read_as_the_superatom_label_it_duplicates_and_not_as_an_alias():
    """``G  aaappp`` names both ends of the bond crossing a contracted group's boundary, and both atoms
    are in the atom block.  So it adds nothing to expand: the atoms are drawn and only their DISPLAY is
    contracted.  Reading the text as atom `aaa`'s alias would put a label on an atom whose group is
    already there, and `expand_abbreviations` would then graft a second copy of it.
    """
    lines = list(_ETHANOL)
    lines[-1:-1] = ['G    3  2', 'OMe', 'M  STY  1   1 SUP', 'M  SAL   1  1   3', 'M  SMT   1 OMe']
    mol, _, log = _read(lines)
    assert not mol.aliases, 'the label belongs to the S-group, not to an atom'
    assert [mol.element_of(s) for s in mol.atom_numbers] == [6, 6, 8], 'ethanol, unexpanded'
    stated = [x for x in log if x.rule == 'v2000:group-abbreviation']
    assert len(stated) == 1 and stated[0].severity == 'info', log
    assert 'SUP S-group already states' in stated[0]
    assert not [x for x in log if 'unrecognised-props-line' in x.rule], \
        'the G line and its text line are both consumed, not two lost lines'


def test_a_g_line_with_no_superatom_to_hold_it_is_a_lost_label():
    """The same line without the ``SUP`` record that gives the group's extent: the label names a group
    the file never delimited, so there are no atoms to carry it."""
    lines = list(_ETHANOL)
    lines[-1:-1] = ['G    3  2', 'OMe']
    mol, _, log = _read(lines)
    assert not mol.aliases
    stated = [x for x in log if x.rule == 'v2000:group-abbreviation']
    assert len(stated) == 1 and stated[0].severity == 'lost', log
    assert 'no SUP S-group gives its extent' in stated[0]


def test_a_query_symbol_is_still_refused_by_name():
    """Negative control for the recovery above: ``A`` and ``L`` are listed constraints rather than
    free text, and a constraint has nowhere to go in a molecule.  ``*`` is excluded: it is now the
    marker element 0, not a query symbol."""
    for symbol in ('A  ', 'Q  ', 'L  '):
        lines = list(_ETHANOL)
        lines[4] = f'    0.0000    0.0000    0.0000 {symbol} 0  0  0  0  0  0  0  0  0  0  0  0'
        with raises(UnsupportedCtfile, match='query reader'):
            _read(lines)


def test_a_blank_symbol_is_still_refused():
    """Negative control: free text is a label, and no text at all is not a label."""
    lines = list(_ETHANOL)
    lines[4] = '    0.0000    0.0000    0.0000     0  0  0  0  0  0  0  0  0  0  0  0'
    with raises(MalformedCtfile, match='no element symbol'):
        _read(lines)


def test_a_record_shorter_than_its_header_is_refused_as_malformed():
    """Negative control: some inputs really are not records."""
    with raises(MalformedCtfile, match='header alone needs'):
        parse_v2000(['only', 'two'], [])


def test_a_counts_line_that_is_not_a_counts_line_is_refused():
    lines = list(_ETHANOL)
    lines[3] = 'this is not a counts line at all'
    with raises(MalformedCtfile, match='counts line'):
        parse_v2000(lines, []).build()


def test_more_atoms_than_the_count_field_holds_is_refused_by_the_writer():
    """1000 atoms cannot be written as V2000 -- the count field is 3 characters -- and the message
    names V3000 as the fix.  Truncating the count produces a file that reads back as another molecule.
    """
    from ....core import MoleculeContainer

    mol = MoleculeContainer()
    with mol.edit():
        for _ in range(1000):
            mol.add_atom('C')
    with raises(MalformedCtfile, match='V3000'):
        emit_v2000(mol)


def test_the_writer_states_a_charge_in_both_places():
    """Many readers look at only one of the two, so a charge is written in both."""
    lines = list(_ETHANOL)
    lines.insert(-1, 'M  CHG  1   3  -1')
    mol, store, _ = _read(lines)
    out, _ = emit_v2000(mol, store)
    assert any(x.startswith('M  CHG') for x in out), out
    charge_column = [x[36:39] for x in out[4:7]]
    assert charge_column[2].strip() == '5', f'ccc should carry the -1 code 5, got {charge_column}'


def test_bond_orders_are_written_as_stored_and_not_normalised():
    lines = list(_ETHANOL)
    lines[7] = '  1  2  2  0  0  0  0'
    mol, store, _ = _read(lines)
    out, _ = emit_v2000(mol, store)
    assert out[7].startswith('  1  2  2'), out[7]


def test_an_sgroup_survives_a_round_trip_through_the_writer(corpus):
    """Every S-group in the corpus, re-emitted and re-parsed, comes back the same.

    Position alphabet on both sides: a CTfile S-group number is not promised to survive, only which
    atoms and bonds the record refers to.
    """
    from .._v2000 import _emit_sgroups

    checked = 0
    for name, records in corpus.items():
        for n, record in _v2000_only(records):
            try:
                ctab = parse_v2000(record, [])
            except MalformedCtfile:
                continue
            if not ctab.sgroups:
                continue
            position = {i: i + 1 for i in range(len(ctab.atoms))}
            bond_position = {}
            for i, b in enumerate(ctab.bonds, 1):
                bond_position[(b.a, b.b)] = bond_position[(b.b, b.a)] = i
            props = _emit_sgroups(SGroupStore(ctab.sgroups), position, bond_position, [])
            head = record[:4 + len(ctab.atoms) + len(ctab.bonds)]
            again = parse_v2000(head + props + ['M  END'], [])
            before = [(r.type, r.name, tuple(r.atoms), tuple(r.data), r.disp)
                      for r in ctab.sgroups]
            after = [(r.type, r.name, tuple(r.atoms), tuple(r.data), r.disp)
                     for r in again.sgroups]
            assert sorted(before) == sorted(after), f'{name} record {n}: S-groups changed'
            checked += 1
    assert checked, 'no S-groups in the corpus; this test proved nothing'


_BENZENE = """benzene
  test
comment
  6  6  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.8000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    1.6000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.6000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5000    0.8000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  4  0  0  0  0
  2  3  4  0  0  0  0
  3  4  4  0  0  0  0
  4  5  4  0  0  0  0
  5  6  4  0  0  0  0
  6  1  4  0  0  0  0
M  END""".split('\n')


def test_charge_code_four_means_a_radical_and_not_a_charge():
    """``ccc=4`` is the one member of the charge-code table that is not a charge: neutral doublet.

    The rest reads as arithmetic (1 is +3, 2 is +2, 3 is +1, 5 is -1 ...), so regenerating the table
    from that pattern gets 4 wrong and drops the radical.  Nothing in the corpus carries a 4.

    The write direction is deliberately not symmetric: a radical goes out as ``M  RAD``, never as
    ``ccc=4``, because ``ccc`` holds one fact and a charged radical has two.
    """
    lines = list(_BENZENE[:4]) + [_BENZENE[4][:36] + '  4' + _BENZENE[4][39:], 'M  END']
    lines[3] = '  1  0  0  0  0  0            999 V2000'
    mol, store, _ = _read(lines)
    sid, = mol.atom_numbers
    assert mol.radical_of(sid), 'charge code 4 read as a plain neutral atom; the radical was dropped'
    assert mol.charge_of(sid) == 0, 'code 4 carries no charge'

    out, _ = emit_v2000(mol, store)
    atom_line, = [x for x in out if x.endswith('0  0  0  0  0  0  0  0  0  0')]
    assert atom_line[36:39].strip() != '4', \
        'a radical was written into ccc, which cannot also carry a charge'
    assert any(x.startswith('M  RAD') and x.rstrip().endswith('2') for x in out), out
    mol2, _, _ = _read(out)
    assert _state(mol) == _state(mol2), 'the radical did not survive the round trip'


def test_bond_type_four_is_stored_as_an_aromatic_bond():
    """A file that says aromatic is stored aromatic: the arena holds order 4 as a first-class
    order, and kekulising on the way in would change what the file said."""
    mol, _, log = _read(_BENZENE)
    assert mol.aromatic_bond_count == 6
    assert [mol.order_of(b.n, b.m) for b in mol.bonds()] == [4] * 6
    assert [mol.implicit_h_of(s) for s in mol.atom_numbers] == [1] * 6, \
        'an aromatic CH is derivable without a Kekule form, and chython 2 agrees'
    assert not any('unreliable' in x for x in log), \
        'an aromatic record is now ordinary; a warning per record would be noise'


def test_a_kekule_drawing_stays_kekule():
    """The other direction: nothing here aromatises either."""
    lines = list(_BENZENE)
    for i, order in enumerate((2, 1, 2, 1, 2, 1)):
        lines[10 + i] = lines[10 + i][:6] + f'{order:3d}' + lines[10 + i][9:]
    mol, _, _ = _read(lines)
    assert mol.aromatic_bond_count == 0
    assert sorted(b.order for b in mol.bonds()) == [1, 1, 1, 2, 2, 2]


def test_an_aromatic_molecule_is_written_as_read_and_kekule_first_still_works():
    """An aromatic bond is a normal bond for this format: bond type 4 is what the block has for one.

    The reference lists type 4 among the query bond types, and the convention on top of that is what
    files carry: measured over five toolkits, three write 4 by default for a ring they perceive as
    aromatic and every one of the five reads 4 back as aromatic (see ``CTFILE.md`` section 12).
    ``kekule()`` first still writes 1 and 2, which is the caller's explicit route.
    """
    mol, _, _ = _read(_BENZENE)
    lines, _ = emit_v2000(mol)
    assert sorted(b.order for b in parse_v2000(lines, []).bonds) == [4, 4, 4, 4, 4, 4]

    mol.kekule()
    lines, _ = emit_v2000(mol)
    assert sorted(b.order for b in parse_v2000(lines, []).bonds) == [1, 1, 1, 2, 2, 2]


def test_an_aromatic_heteroatom_whose_count_needs_a_kekule_form_is_marked_not_refused():
    """Pyridine drawn with bond type 4: the nitrogen's count needs a Kekule form nobody supplied,
    so it is genuinely unknown and chython 2 answers ``None`` for it too.

    The record is read anyway: the count is stored as ``H_UNKNOWN``, registered in
    ``unknown_hydrogens`` and logged per atom.  Answering ``None`` rather than 0 is what distinguishes
    "nobody could say" from "this nitrogen has no hydrogen".
    """
    lines = list(_BENZENE)
    lines[4] = '    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0'
    ctab = parse_v2000(lines, [])
    mol, _, log = ctab.build()
    assert mol.aromatic_bond_count == 6, 'stored as drawn; this reader does not kekulise'
    assert len(ctab.unknown_hydrogens) == 1, 'findable, not asserted as zero'
    nitrogen, = ctab.unknown_hydrogens
    assert mol.implicit_h_of(nitrogen) is None, 'the sentinel, not a zero the caller cannot question'
    assert mol.unknown_h_count == 1
    assert any('no Kekule form' in x for x in log), log
    assert any('unknown implicit hydrogen count' in x for x in log), \
        'the summary line too, so one grep answers the question for a whole SDF'


def test_implicit_hydrogens_are_computed_and_not_left_at_zero():
    """Building an atom derives no count and an unset one stores as ``H_UNKNOWN``, so what is
    pinned here is that this reader runs the core's derivation for every atom."""
    mol, _, _ = _read(_ETHANOL)
    for sid in mol.atom_numbers:
        rule, _ = implicit_for_atom(mol, sid)
        assert mol.implicit_h_of(sid) == rule


def test_a_non_utf8_title_is_emitted_as_the_byte_it_was():
    """The CTfile writer takes NO loss on a name line: it emits the surrogate and the stream re-encodes
    it, which is what makes the round trip byte for byte.

    The title byte ``\\xe9`` is the latin-1 e-acute and not valid UTF-8 on its own.  Nothing here
    replaces it, so there is no replacement count to check and no genuine U+FFFD to be mistaken for a
    substituted byte.  The loss is XML's alone; see ``chython/formats/xml/_dialect.py::xml_text``.
    """
    from ....core import read_smiles

    # methane, so every implicit hydrogen count is known and an empty log means the title cost nothing
    mol = read_smiles('C')
    mol.set_title(b'caf\xe9')
    lines, log = emit_v2000(mol)
    assert not log, log
    assert lines[0].encode('utf8', 'surrogateescape') == b'caf\xe9'


def test_build_puts_meta_on_the_molecule():
    ctab = parse_v2000(_ETHANOL, [])
    ctab.meta['k'] = 'v'
    mol, _, _ = ctab.build()
    assert mol.meta == {'k': 'v'}


def test_build_puts_the_log_on_the_molecule_as_well_as_returning_it():
    """A charge code outside 0-7 is read as neutral and reported, so this record has a line to carry.

    The plan reached for an unstamped counts line; that one is `sniff_version`'s, and `_read` calls
    `parse_v2000` directly, so the log would have been empty and the assertion vacuous.
    """
    lines = list(_ETHANOL)
    lines[4] = lines[4].replace('C   0  0', 'C   0  9')
    mol, _, log = _read(lines)
    assert log and [str(x) for x in mol.log] == [str(x) for x in log]
