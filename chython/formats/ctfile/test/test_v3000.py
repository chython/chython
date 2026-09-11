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
"""V3000: keywords instead of columns, and indices instead of positions.

A V3000 atom index is **not** a position -- any positive integer, in any order, with gaps -- so
``INDEX=7`` means the atom whose line began with 7.  The parser keeps a map; the writer regenerates
indices from scratch.  ``test/implicit.sdf`` has exactly one V3000 record, so the corpus contribution
is thin and everything else is written out inline.
"""

from pytest import raises

from .._ctab import STEREO_AND, STEREO_OR
from .._errors import MalformedCtfile, UnsupportedCtfile
from .._v2000 import emit_v2000, parse_v2000
from .._hydrogens import MRV_IMPLICIT_H
from .._sdf import sniff_version, split_records
from .._sgroup import NO_INDEX
from .._v3000 import V3000_STAMP, emit_v3000, parse_v3000
from ....core import read_smiles


def _record(*body, title='t'):
    return [title, '  test', '', '  0  0  0  0  0  0            999 V3000',
            'M  V30 BEGIN CTAB', *[f'M  V30 {x}' for x in body], 'M  V30 END CTAB', 'M  END']


_ETHANOL = _record('COUNTS 3 2 0 0 0',
                   'BEGIN ATOM',
                   '1 C 0 0 0 0',
                   '2 C 1.5 0 0 0',
                   '3 O 3 0 0 0',
                   'END ATOM',
                   'BEGIN BOND',
                   '1 1 1 2',
                   '2 1 2 3',
                   'END BOND')


def _read(lines):
    return parse_v3000(lines, []).build()


def test_the_baseline_record_reads_as_ethanol():
    mol, _, log = _read(_ETHANOL)
    assert [mol.element_of(s) for s in mol.atom_numbers] == [6, 6, 8]
    assert [mol.implicit_h_of(s) for s in mol.atom_numbers] == [3, 2, 1]


def test_an_atom_index_is_not_a_position():
    """Indices 7, 3, 99 in that order: a reader using them as positions gets the bond wrong or
    crashes."""
    mol, _, log = _read(_record('COUNTS 3 2 0 0 0',
                                'BEGIN ATOM', '7 C 0 0 0 0', '3 O 1.5 0 0 0', '99 N 3 0 0 0',
                                'END ATOM',
                                'BEGIN BOND', '1 1 7 3', '2 2 3 99', 'END BOND'))
    elements = [mol.element_of(s) for s in mol.atom_numbers]
    assert elements == [6, 8, 7], 'atoms keep file order, not index order'
    sids = list(mol.atom_numbers)
    assert mol.order_of(sids[0], sids[1]) == 1
    assert mol.order_of(sids[1], sids[2]) == 2


def test_a_repeated_atom_index_is_reported_and_the_second_ignored():
    """Two atoms cannot share an index: every reference to it would be ambiguous."""
    log = []
    ctab = parse_v3000(_record('COUNTS 2 0 0 0 0', 'BEGIN ATOM', '1 C 0 0 0 0', '1 O 1 0 0 0',
                               'END ATOM'), log)
    assert len(ctab.atoms) == 1
    assert any('repeated' in x for x in log), log


def test_a_bond_to_an_unknown_index_is_dropped_with_a_log():
    log = []
    ctab = parse_v3000(_record('COUNTS 2 1 0 0 0', 'BEGIN ATOM', '1 C 0 0 0 0', '2 C 1 0 0 0',
                               'END ATOM', 'BEGIN BOND', '1 1 1 9', 'END BOND'), log)
    assert not ctab.bonds
    assert any('unknown atom index' in x for x in log), log


def _one_bond(type_):
    """An Fe-N record whose single bond carries `type_`.  Iron so that no valence rule interferes."""
    return _record('COUNTS 2 1 0 0 0',
                   'BEGIN ATOM', '1 Fe 0 0 0 0', '2 N 1.5 0 0 0', 'END ATOM',
                   'BEGIN BOND', f'1 {type_} 1 2', 'END BOND')


def _read_with_log(lines, log):
    """:func:`_read` with the parse log and the build log in one list, which is what a caller sees."""
    ctab = parse_v3000(lines, log)
    mol, store, build_log = ctab.build()
    log.extend(build_log)
    return mol, store, log


def test_a_v3000_bond_type_goes_through_the_same_translation_as_a_v2000_one():
    """Both readers call ``order_from_bond_type``, so a file's bond type is translated rather than
    passed into ``CtabBond.order`` raw for the storability guard to mop up."""
    for type_, order in ((1, 1), (2, 2), (3, 3), (4, 4), (8, 8), (9, 8), (10, 8), (12, 1)):
        mol, _, _ = _read_with_log(_one_bond(type_), [])
        a, b = mol.atom_numbers
        assert mol.order_of(a, b) == order, f'type {type_} read as order {mol.order_of(a, b)}'


def test_a_vendors_coordination_bond_and_our_own_type_8_now_reach_the_same_order():
    """Both emitters put ``bond.order`` into the type column, so a chython dative bond goes out as
    type 8 in either version; the spec's own coordination type is 9.  Both must come back as order 8,
    or a conformant vendor file degrades to a single bond while our own output round-trips.
    """
    log8, log9 = [], []
    m8, _, _ = _read_with_log(_one_bond(8), log8)
    m9, _, _ = _read_with_log(_one_bond(9), log9)
    assert m8.order_of(*m8.atom_numbers) == m9.order_of(*m9.atom_numbers) == 8

    # Each says which it was: the two are not the same statement about the file.
    assert any('type 8 read as chython order 8' in x and 'query "any bond"' in x for x in log8), log8
    assert any('coordination bond type 9 read as chython order 8' in x for x in log9), log9
    assert not any('read as single' in x for x in log8 + log9), log8 + log9


def test_writing_a_dative_bond_states_type_8_and_says_nothing_about_it():
    """The asymmetry is deliberate: the reader warns about type 8 because the spec calls it the
    query *any bond* and it cannot know the file's provenance; the writer knows, so it says nothing.

    Asserted as an empty log, and against the emitted column in both versions, since a writer that
    dropped the bond would also have an empty log.  See ``_ctab.order_from_bond_type``.
    """
    mol = read_smiles('[Fe]~N(C)(C)C')
    assert 8 in {b.order for b in mol.bonds()}, 'the fixture carries no dative bond'

    lines, log = emit_v3000(mol)
    assert log == [], log
    assert 'M  V30 1 8 1 2' in lines, lines

    lines2, log2 = emit_v2000(mol)
    assert log2 == [], log2
    assert '  1  2  8  0  0  0  0' in lines2, lines2


def test_a_v3000_query_bond_refuses_the_record_exactly_as_a_v2000_one_does():
    """The file is well formed and says "either of two orders", which a molecule has nowhere to
    put, so this is a refusal at the answer boundary naming the type.  One specification cannot have
    two readings of type 5, so V2000 and V3000 must refuse alike."""
    for type_ in (5, 6, 7):
        with raises(UnsupportedCtfile, match='query bond'):
            parse_v3000(_one_bond(type_), [])


def test_a_hydrogen_bond_is_stored_as_a_contact_and_not_as_a_covalent_bond():
    """V3000's type 10, which V2000 has no spelling for.  chython has no hydrogen-bond order, and
    order 1 would join the valence arithmetic while order 8 is excluded from valence, degree and
    heteroatom counts.  So the contact survives, its kind does not, and ``unsupported: `` says so."""
    log = []
    mol, _, _ = _read_with_log(_one_bond(10), log)
    a, b = mol.atom_numbers
    assert mol.order_of(a, b) == 8
    assert mol.degree_of(a) == 1, 'stored structurally, like every order-8 contact'
    lines = [x for x in log if 'hydrogen bond type 10' in x]
    assert lines and str(lines[0]).startswith('unsupported: '), log
    assert 'not preserved' in lines[0], lines


def test_the_two_ctab_versions_agree_on_every_bond_type_they_share():
    """A V2000 bond line and a V3000 bond entry carrying the same type must produce the same order.
    Both readers call ``order_from_bond_type``, so there is no second spelling to drift.
    """
    for type_ in (1, 2, 3, 4, 8, 9, 12):
        v3, _, _ = _read_with_log(_one_bond(type_), [])
        v2_lines = ['t', '  test', '', '  2  1  0  0  0  0            999 V2000',
                    '    0.0000    0.0000    0.0000 Fe  0  0  0  0  0  0  0  0  0  0  0  0',
                    '    1.5000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0',
                    f'  1  2  {type_}  0  0  0  0', 'M  END']
        v2, _, _ = parse_v2000(v2_lines, []).build()
        assert v3.order_of(*v3.atom_numbers) == v2.order_of(*v2.atom_numbers), f'type {type_}'


def test_the_counts_line_is_compared_but_the_blocks_win():
    """A wrong COUNTS is a writer bug and the blocks are the data: trusting the count truncates a
    file whose atom block is longer than it claims."""
    log = []
    ctab = parse_v3000(_record('COUNTS 9 8 0 0 0', 'BEGIN ATOM', '1 C 0 0 0 0', 'END ATOM'), log)
    assert len(ctab.atoms) == 1
    assert any('COUNTS says 9 atoms' in x for x in log), log


def test_a_ctab_with_no_atom_block_is_refused():
    """Negative control: a record with no atoms is not a molecule with no atoms, and this is the
    shape a V3000 record takes when a V2000 parser reads it."""
    with raises(MalformedCtfile, match='atom block'):
        parse_v3000(_record('COUNTS 0 0 0 0 0'), [])


def test_a_missing_m_end_is_not_fatal():
    log = []
    lines = [x for x in _ETHANOL if not x.startswith('M  END')]
    ctab = parse_v3000(lines, log)
    assert len(ctab.atoms) == 3
    assert any('no M  END' in x for x in log), log


def test_an_unclosed_block_is_read_to_the_end_and_reported():
    log = []
    ctab = parse_v3000(_record('COUNTS 1 0 0 0 0', 'BEGIN ATOM', '1 C 0 0 0 0'), log)
    assert len(ctab.atoms) == 1
    assert any('not closed' in x for x in log), log


def test_blocks_may_arrive_in_any_order():
    """Files with SGROUP before ATOM exist, so index resolution is deferred to the end of the
    CTAB."""
    mol, store, log = _read(_record('COUNTS 2 1 0 0 0',
                                    'BEGIN SGROUP',
                                    '1 DAT 0 ATOMS=(1 2) FIELDNAME=X FIELDDATA=v',
                                    'END SGROUP',
                                    'BEGIN ATOM', '1 C 0 0 0 0', '2 O 1.5 0 0 0', 'END ATOM',
                                    'BEGIN BOND', '1 1 1 2', 'END BOND'))
    assert len(store.records) == 1 and len(store.records[0].atoms) == 1


# atom keywords

def test_charge_mass_and_radical_keywords():
    mol, _, _ = _read(_record('COUNTS 2 0 0 0 0',
                              'BEGIN ATOM', '1 N 0 0 0 0 CHG=1 MASS=15', '2 O 1 0 0 0 RAD=2',
                              'END ATOM'))
    sids = list(mol.atom_numbers)
    assert mol.charge_of(sids[0]) == 1 and mol.isotope_of(sids[0]) == 15
    assert bool(mol.radical_of(sids[1]))


def test_a_diradical_is_stored_as_one_radical_and_reported():
    """The core holds one radical bit, so a singlet or triplet diradical cannot be said in full;
    half the truth is kept, labelled."""
    log = []
    ctab = parse_v3000(_record('COUNTS 1 0 0 0 0', 'BEGIN ATOM', '1 C 0 0 0 0 RAD=3', 'END ATOM'),
                       log)
    assert ctab.atoms[0].radical
    assert any('triplet diradical' in x for x in log), log


def test_val_minus_one_is_the_spec_s_spelling_of_valence_zero():
    """A bare metal ion states VAL=-1, which is a real statement and not the same as VAL absent."""
    ctab = parse_v3000(_record('COUNTS 1 0 0 0 0', 'BEGIN ATOM', '1 Na 0 0 0 0 VAL=-1 CHG=1',
                               'END ATOM'), [])
    assert ctab.atoms[0].valence == 0


def test_a_stated_zero_valence_is_written_back_as_val_minus_one_and_never_as_fifteen():
    """The two versions spell "valence zero" differently, and this is the only test on the seam.

    ``ZERO_VALENCE`` is 15 because that is V2000's ``vvv`` spelling, while 15 in V3000's ``VAL`` is a
    valence of *fifteen*, so the emitter's translation is load-bearing in one direction only.
    """
    from ....core import MoleculeContainer

    # A bare carbon holding no hydrogens: the rules say four, so the writer states the disagreement
    # as a total valence of zero.  Built rather than read, since the rules outrank a stated valence.
    mol = MoleculeContainer()
    with mol.edit():
        sid = mol.add_atom('C')
        mol.set_hydrogens(sid, 0)

    lines, _ = emit_v3000(mol)
    atom_line, = [x for x in lines if ' C ' in x]
    assert 'VAL=-1' in atom_line, atom_line
    assert 'VAL=15' not in atom_line, "V2000's spelling of zero valence leaked into a V3000 file"

    mol2, _, _ = parse_v3000(lines, []).build()
    assert mol2.implicit_h_of(next(iter(mol2.atom_numbers))) == 0, \
        'the zero survives on the MRV_IMPLICIT_H channel, which fires on the same predicate'


def test_val_minus_one_on_a_bond_free_atom_produces_zero_hydrogens_end_to_end():
    """VAL=-1 must carry from the parser to the hydrogen count, not just be stored internally.

    The two adjacent tests each pin half the path and both stay green without the parser's VAL=-1
    handling, because the MRV_IMPLICIT_H channel fires on the same predicate and outranks the
    derivation.  Carbon, because its free-atom default of 4 is far enough from zero that the log's
    "preference to the 4" shows the stated-zero path ran.  The V2000 analog is vvv=15 in
    test_hydrogens.py.
    """
    mol, _, log = _read(_record('COUNTS 1 0 0 0 0',
                                'BEGIN ATOM',
                                '1 C 0 0 0 0 VAL=-1',
                                'END ATOM'))
    assert mol.implicit_h_of(next(iter(mol.atom_numbers))) == 0, \
        'stated zero valence: not methane'
    assert any('no bonds drawn' in x and 'preference to the 4' in x for x in log), log


def test_hcount_is_logged_and_ignored():
    """A query field: "n or more", not a count."""
    log = []
    mol, _, build_log = parse_v3000(_record('COUNTS 1 0 0 0 0',
                                            'BEGIN ATOM', '1 C 0 0 0 0 HCOUNT=2', 'END ATOM'),
                                    log).build()
    assert mol.implicit_h_of(next(iter(mol.atom_numbers))) == 4
    assert any('query field' in x for x in log), log


def test_a_free_text_atom_type_is_read_as_an_alias_and_not_as_a_lost_record():
    """The same recovery as V2000's: a label in the atom-type field keeps its record.  V3000 has no
    alias line of its own, so this is the only channel the text has here.
    """
    ctab = parse_v3000(_record('COUNTS 3 2 0 0 0',
                               'BEGIN ATOM', '1 Me 0 0 0 0', '2 C 1.5 0 0 0', '3 O 3 0 0 0',
                               'END ATOM',
                               'BEGIN BOND', '1 1 1 2', '2 1 2 3', 'END BOND'), [])
    mol, store, log = ctab.build()
    sid, *_ = mol.atom_numbers
    assert len([*mol.atom_numbers]) == 3
    assert mol.aliases == {sid: b'Me'}
    assert mol.element_of(sid) == 0, 'the marker, not a borrowed element'
    assert mol.implicit_h_of(sid) == 0
    assert not ctab.unknown_hydrogens, 'nothing about a marker is undeterminable'
    assert any('names no element' in x for x in log), log


def test_a_query_atom_type_is_still_refused_by_name():
    """Negative control: the listed query types stay refusals, a constraint being a different thing
    from a label."""
    for token in ('A', 'Q', 'M'):
        with raises(UnsupportedCtfile, match='query'):
            _read(_record('COUNTS 1 0 0 0 0', 'BEGIN ATOM', f'1 {token} 0 0 0 0', 'END ATOM'))


def test_an_unknown_atom_keyword_is_named_in_the_log():
    log = []
    parse_v3000(_record('COUNTS 1 0 0 0 0', 'BEGIN ATOM', '1 C 0 0 0 0 ZZZ=1', 'END ATOM'), log)
    assert any('ZZZ' in x for x in log), log


# stereo groups

def test_an_enhanced_stereo_collection_is_read_with_the_right_semantics():
    """STERAC is a racemate -- both enantiomers -- so it is AND; STEREL is one enantiomer of
    unknown configuration, so it is OR.  Inverting them changes what the file says about a sample."""
    ctab = parse_v3000(_record('COUNTS 2 0 0 0 0',
                               'BEGIN ATOM', '1 C 0 0 0 0', '2 C 1 0 0 0', 'END ATOM',
                               'BEGIN COLLECTION',
                               'MDLV30/STERAC1 ATOMS=(1 1)',
                               'MDLV30/STEREL2 ATOMS=(1 2)',
                               'END COLLECTION'), [])
    assert ctab.groups == {0: (STEREO_AND, 1), 1: (STEREO_OR, 2)}


def test_a_collection_referring_to_an_unknown_atom_is_reported():
    log = []
    parse_v3000(_record('COUNTS 1 0 0 0 0', 'BEGIN ATOM', '1 C 0 0 0 0', 'END ATOM',
                        'BEGIN COLLECTION', 'MDLV30/STEABS ATOMS=(1 9)', 'END COLLECTION'), log)
    assert any('unknown atom index' in x for x in log), log


def test_a_group_id_above_the_storable_range_is_renumbered_and_not_dropped():
    """The id is a label -- the partition is the statement -- so an id the arena cannot hold is mapped
    to a free one of its own KIND rather than costing the record its group.  `MDLV30/STERAC1384` occurs
    in the wild.

    Three things the mapping has to get right: every line naming 1384 lands in one group, an id the
    file itself uses is not stolen for the renumbering, and AND and OR number independently -- so
    STEREL64 may become OR 1 while AND 1 is taken.
    """
    log = []
    ctab = parse_v3000(_record('COUNTS 4 0 0 0 0',
                               'BEGIN ATOM', '1 C 0 0 0 0', '2 C 1 0 0 0', '3 C 2 0 0 0', '4 C 3 0 0 0',
                               'END ATOM',
                               'BEGIN COLLECTION',
                               'MDLV30/STERAC1384 ATOMS=(1 1)',
                               'MDLV30/STERAC1 ATOMS=(1 2)',
                               'MDLV30/STERAC1384 ATOMS=(1 3)',
                               'MDLV30/STEREL64 ATOMS=(1 4)',
                               'END COLLECTION'), log)
    assert ctab.groups == {0: (STEREO_AND, 2), 1: (STEREO_AND, 1), 2: (STEREO_AND, 2), 3: (STEREO_OR, 1)}
    assert sum('renumbered' in x for x in log) == 2, log


def test_an_unknown_collection_is_ignored_by_name():
    log = []
    parse_v3000(_record('COUNTS 1 0 0 0 0', 'BEGIN ATOM', '1 C 0 0 0 0', 'END ATOM',
                        'BEGIN COLLECTION', 'MDLV30/HILITE ATOMS=(1 1)', 'END COLLECTION'), log)
    assert any('HILITE' in x for x in log), log


# writing and back

def test_read_write_read_is_identical():
    mol, store, _ = _read(_ETHANOL)
    lines, _ = emit_v3000(mol, store, title='t')
    assert sniff_version(lines, []) == V3000_STAMP
    mol2, _, _ = _read(lines)
    assert [mol2.element_of(s) for s in mol2.atom_numbers] == [6, 6, 8]
    assert [mol2.implicit_h_of(s) for s in mol2.atom_numbers] == [3, 2, 1]


def test_the_writer_regenerates_indices_from_scratch():
    """Indices are positions in the written file and nothing else; the format does not treat the
    original numbers as an identity."""
    mol, store, _ = _read(_record('COUNTS 2 1 0 0 0',
                                  'BEGIN ATOM', '7 C 0 0 0 0', '99 O 1.5 0 0 0', 'END ATOM',
                                  'BEGIN BOND', '4 1 7 99', 'END BOND'))
    lines, _ = emit_v3000(mol, store)
    atoms = [x for x in lines if x.startswith('M  V30 1 C') or x.startswith('M  V30 2 O')]
    assert len(atoms) == 2, lines


def test_a_thousand_atoms_is_not_a_problem_for_v3000():
    """Why the V2000 writer's refusal names V3000: no fixed columns, so it fits here."""
    from ....core import MoleculeContainer

    mol = MoleculeContainer()
    with mol.edit():
        for _ in range(1000):
            mol.add_atom('C')
    lines, _ = emit_v3000(mol)
    assert any('COUNTS 1000' in x for x in lines)


def test_an_unknown_hydrogen_count_suppresses_the_valence_keyword():
    """Same rule as V2000: a count nobody stated is not written as a fact.  The implicit-H data
    S-groups are counted on the COUNTS line, so they are built before it is written."""
    ctab = parse_v3000(_record('COUNTS 6 5 0 0 0',
                               'BEGIN ATOM', '1 C 0 0 0 0', '2 F 1 0 0 0', '3 F 2 0 0 0',
                               '4 F 3 0 0 0', '5 F 4 0 0 0', '6 F 5 0 0 0', 'END ATOM',
                               'BEGIN BOND', '1 1 1 2', '2 1 1 3', '3 1 1 4', '4 1 1 5',
                               '5 1 1 6', 'END BOND'), [])
    mol, store, _ = ctab.build()
    assert ctab.unknown_hydrogens
    lines, _ = emit_v3000(mol, store)
    assert not any('VAL=' in x for x in lines), lines
    assert not any(MRV_IMPLICIT_H in x for x in lines), lines


def test_a_stated_hydrogen_count_is_read_from_the_sgroup_block():
    """``MRV_IMPLICIT_H`` is honoured by both parsers, through the one helper in ``_hydrogens``,
    called after the S-group block is parsed because that is where the statement arrives.

    Pyrrole, drawn aromatic: the nitrogen's count is not derivable without a Kekule form, which is
    the case the group gets written for."""
    lines = _record('COUNTS 5 5 1 0 0',
                    'BEGIN ATOM', '1 N 0 0 0 0', '2 C 1 0 0 0', '3 C 2 0 0 0', '4 C 3 0 0 0',
                    '5 C 4 0 0 0', 'END ATOM',
                    'BEGIN BOND', '1 4 1 2', '2 4 2 3', '3 4 3 4', '4 4 4 5', '5 4 5 1',
                    'END BOND',
                    'BEGIN SGROUP',
                    f'1 DAT 0 ATOMS=(1 1) FIELDNAME={MRV_IMPLICIT_H} FIELDDATA=IMPL_H1',
                    'END SGROUP')
    ctab = parse_v3000(lines, [])
    mol, _, _ = ctab.build()
    assert mol.implicit_h_of(next(iter(mol.atom_numbers))) == 1
    assert not ctab.unknown_hydrogens, 'the file stated the count, so nothing is unknown'


def test_without_that_sgroup_the_same_count_is_unknown():
    """The negative control: the same pyrrole minus the one line stating its nitrogen's count is
    read and the nitrogen marked, which distinguishes *stated 1* from *not known* on one record."""
    lines = _record('COUNTS 5 5 0 0 0',
                    'BEGIN ATOM', '1 N 0 0 0 0', '2 C 1 0 0 0', '3 C 2 0 0 0', '4 C 3 0 0 0',
                    '5 C 4 0 0 0', 'END ATOM',
                    'BEGIN BOND', '1 4 1 2', '2 4 2 3', '3 4 3 4', '4 4 4 5', '5 4 5 1',
                    'END BOND')
    ctab = parse_v3000(lines, [])
    mol, _, log = ctab.build()
    nitrogen = next(iter(mol.atom_numbers))
    assert ctab.unknown_hydrogens == (nitrogen,), ctab.unknown_hydrogens
    assert mol.implicit_h_of(nitrogen) is None, \
        'the sentinel: `None` is not a number that can be mistaken for an answer'
    assert mol.unknown_h_count == 1
    assert any('no Kekule form' in x for x in log), log


def test_a_data_sgroup_round_trips():
    mol, store, _ = _read(_record('COUNTS 2 1 0 0 0',
                                  'BEGIN ATOM', '1 C 0 0 0 0', '2 O 1.5 0 0 0', 'END ATOM',
                                  'BEGIN BOND', '1 1 1 2', 'END BOND',
                                  'BEGIN SGROUP',
                                  '1 DAT 0 ATOMS=(1 1) FIELDNAME=BP FIELDDATA=42',
                                  'END SGROUP'))
    lines, _ = emit_v3000(mol, store)
    _, store2, _ = _read(lines)
    assert [(r.name, r.field_data) for r in store2.records] == [('BP', '42')]


def test_what_the_parser_reported_is_still_in_the_build_log():
    """``Ctab.build`` starts its log from ``Ctab.log``, so a parser diagnostic reaches the caller
    without them having to pass and then read their own list."""
    _, _, log = _read(_record('COUNTS 1 0 0 0 0', 'BEGIN ATOM', '1 C 0 0 0 0 ZZZ=1', 'END ATOM'))
    assert any('ZZZ' in x for x in log), log


def _sgroup_record(*sgroup_lines, bond='1 1 1 2'):
    """`bond` is exposed because a file may declare a bond's endpoints in either order, and which
    order it chose is not recoverable from the molecule -- see the backwards-bond test below."""
    return _record('COUNTS 2 1 0 0 0',
                   'BEGIN ATOM', '1 C 0 0 0 0', '2 O 1.5 0 0 0', 'END ATOM',
                   'BEGIN BOND', bond, 'END BOND',
                   'BEGIN SGROUP', *sgroup_lines, 'END SGROUP')


def test_an_unmodelled_index_valued_keyword_is_dropped_and_not_re_emitted():
    """An index-valued keyword may be translated or dropped, never passed through.

    ``SAP`` names a superatom's attachment atom by the file's numbering and the writer regenerates
    that numbering, so re-emitting the value verbatim moves the attachment point -- silently, and only
    on files whose atoms were not already numbered 1..n in order.
    """
    mol, store, log = _read(_sgroup_record('1 SUP 0 ATOMS=(2 1 2) SAP=(3 1 2 1) LABEL=Ph'))
    assert any('SAP' in x and 'not be re-emitted' in x for x in log), log
    lines, _ = emit_v3000(mol, store)
    assert not any('SAP' in x for x in lines), lines
    # What is modelled still survives, so the drop is targeted.
    assert any('LABEL=Ph' in x for x in lines), lines


def test_the_bond_keyword_of_the_other_sgroup_type_is_dropped_too():
    """``XBONDS`` on a ``DAT`` group is bond indices under a keyword that record type does not
    use.  Merging it into the modelled bond list would re-emit it as ``CBONDS``, so it is dropped."""
    mol, store, log = _read(_sgroup_record('1 DAT 0 ATOMS=(1 1) XBONDS=(1 1) FIELDNAME=BP'))
    assert any('XBONDS' in x and 'not be re-emitted' in x for x in log), log
    lines, _ = emit_v3000(mol, store)
    assert not any('XBONDS' in x for x in lines), lines


def test_a_crossing_bond_is_written_back_as_a_bond_number():
    """``XBONDS`` on an ``SRU`` is the one S-group keyword whose values are bond *positions*, so it
    is the V3000 writer's only consumer of the bond-number table -- every other S-group test here uses
    a record type whose bond list is dropped on the way in.
    """
    mol, store, _ = _read(_sgroup_record('1 SRU 0 ATOMS=(2 1 2) XBONDS=(1 1) CONNECT=HT'))
    assert [tuple(r.bonds) for r in store.records] == [((1, 2),)], 'an endpoint pair, not a number'
    lines, log = emit_v3000(mol, store)
    assert not log, log
    assert any('XBONDS=(1 1)' in x for x in lines), lines


def test_a_cstate_is_written_back_with_its_bond_number():
    """``CSTATE``'s leading value is a bond position too, resolved through the same table by a
    separate line of code."""
    mol, store, _ = _read(_sgroup_record('1 SRU 0 ATOMS=(2 1 2) XBONDS=(1 1) '
                                         'CSTATE=(4 1 0.5 0.5 0)'))
    lines, log = emit_v3000(mol, store)
    assert not log, log
    assert any('CSTATE=(4 1 0.5 0.5 0)' in x for x in lines), lines


def test_a_bond_the_file_declared_backwards_is_still_found():
    """A file may write a bond as ``2 1``, so the stored pair is ``(2, 1)`` while the molecule
    reports ``(1, 2)``.  The file's order is not recoverable, so the writer's bond-number table holds
    both orientations -- which a test normalising the pair before comparing cannot see.  This one
    asserts the reference survived instead.
    """
    for keyword, where in (('XBONDS=(1 1)', 'bonds'), ('CSTATE=(4 1 0.5 0.5 0)', 'cstates')):
        # The two keywords land in different fields and resolve on different lines of the writer,
        # so one record would let either carry the other.
        record = _sgroup_record(f'1 SRU 0 ATOMS=(2 1 2) {keyword}', bond='1 1 2 1')
        mol, store, _ = _read(record)
        stored = store.records[0].bonds if where == 'bonds' \
            else [p for p, _ in store.records[0].cstates]
        assert list(stored) == [(2, 1)], \
            f'{keyword}: the pair is stored in the order the FILE gave, which is the point'
        assert [(b.n, b.m) for b in mol.bonds()] == [(1, 2)], \
            f'{keyword}: while the molecule reports the other order'
        lines, log = emit_v3000(mol, store)
        assert not log, f'{keyword}: {log}'
        assert any(keyword in x for x in lines), f'{keyword}: {lines}'


def test_the_bond_block_names_its_endpoints_in_the_molecules_own_order():
    """The writer promises byte-identical output for two files stating the same thing, and a bond
    line is two of those bytes.  A read-write-read comparison is orientation-blind."""
    mol, _, _ = _read(_ETHANOL)
    written = [x[7:] for x in emit_v3000(mol)[0] if x.startswith('M  V30 ')]
    body = written[written.index('BEGIN BOND') + 1:written.index('END BOND')]
    assert body == ['1 1 1 2', '2 1 2 3'], body


def test_a_parent_reference_follows_the_sgroups_renumbering():
    """``PARENT`` holds an S-group index and the writer renumbers S-groups by position, so it is
    translated like any other index.  The groups are 3 and 7 because with 1 and 2 the bug is
    invisible."""
    mol, store, _ = _read(_sgroup_record('3 SUP 0 ATOMS=(1 1) LABEL=A',
                                         '7 SUP 0 ATOMS=(1 2) LABEL=B PARENT=3'))
    lines, log = emit_v3000(mol, store)
    parents = [x for x in lines if 'PARENT' in x]
    assert len(parents) == 1 and 'PARENT=1' in parents[0], (parents, log)


def test_a_parent_naming_a_group_that_is_gone_is_reported_not_guessed():
    mol, store, _ = _read(_sgroup_record('7 SUP 0 ATOMS=(1 2) LABEL=B PARENT=3'))
    lines, log = emit_v3000(mol, store)
    assert not any('PARENT' in x for x in lines), lines
    assert any('PARENT=3' in x and 'dropped' in x for x in log), log


def test_the_case_of_a_datum_survives_in_both_directions():
    """A ``STEREOLABEL`` datum is a CIP descriptor, and lowercase ``r``/``s`` are the
    pseudo-asymmetric descriptors from the auxiliary rules -- a different kind of centre, not a
    spelling of ``R``/``S``.  So only the V3000 keyword is case-folded, never the value, and not the
    field's own name either.
    """
    for stated in ('R', 'r'):
        mol, store, _ = _read(_sgroup_record(f'1 DAT 1 ATOMS=(1 1) FIELDNAME=STEREOLABEL '
                                             f'FIELDDATA={stated}'))
        assert store.records[0].data == [stated.encode()], stated
        assert store.records[0].name == 'STEREOLABEL'
        lines, _ = emit_v3000(mol, store)
        assert any(f'FIELDDATA={stated}' in x for x in lines), (stated, lines)


def test_an_unmodelled_keyword_keeps_its_values_in_order_but_not_its_place_on_the_line():
    """The declared limit of "rides through verbatim".

    ``fields`` is keyed by keyword, so keywords come back sorted and interleaving across keywords is
    gone, while each keyword's own values keep the file's order -- the half that carries meaning
    (``BRKXYZ`` twice is two brackets, in order).  Sorting keywords is what makes two files stating
    the same S-groups produce byte-identical output.
    """
    mol, store, _ = _read(_sgroup_record('1 SUP 0 ATOMS=(1 1) LABEL=a NATREPLACE=x LABEL=b'))
    assert store.records[0].fields == {'LABEL': ['a', 'b'], 'NATREPLACE': ['x']}
    line = next(x for x in emit_v3000(mol, store)[0] if 'NATREPLACE' in x)
    assert line.index('LABEL=a') < line.index('LABEL=b'), 'a keyword\'s own values keep their order'
    assert line.index('LABEL=b') < line.index('NATREPLACE'), 'and the keywords themselves are sorted'


def test_the_number_the_model_spells_absent_with_is_not_read_as_absent():
    """V3000 states an S-group number as an unbounded integer, so a file may state 65535 -- the
    value this model uses for "unnumbered".  Read as unnumbered, the record loses the number its own
    keywords refer to it by, so the guard is explicit rather than implied by the field's size."""
    mol, store, log = _read(_sgroup_record('65535 DAT 0 ATOMS=(1 1) FIELDNAME=BP'))
    assert len(store.records) == 1
    assert store.records[0].index != NO_INDEX, 'read as "no number"'
    assert any('65535' in x and 'renumbered' in x for x in log), log


def test_a_number_outside_the_domain_does_not_take_one_another_group_states():
    """The replacement is the lowest *free* number, not 1: group 1 is taken here, and handing it out
    twice makes two records answer to the same reference."""
    mol, store, log = _read(_sgroup_record('1 SUP 0 ATOMS=(1 1) LABEL=A',
                                           '70000 SUP 0 ATOMS=(1 2) LABEL=B'))
    indices = [r.index for r in store.records]
    assert len(set(indices)) == 2 and all(0 <= i < NO_INDEX for i in indices), indices
    assert any('70000' in x and 'renumbered' in x for x in log), log


def test_a_parent_the_model_cannot_hold_is_reported_rather_than_read_as_no_parent():
    """``PARENT=65535`` and an omitted ``PARENT`` land on the same stored value, so the check is at
    the parse site rather than in the renumbering pass.  The reference is lost out loud."""
    mol, store, log = _read(_sgroup_record('1 SUP 0 ATOMS=(1 1) LABEL=A PARENT=65535'))
    assert store.records[0].parent == NO_INDEX
    assert any('PARENT=65535' in x and 'dropped' in x for x in log), log


def test_a_bond_reference_whose_bond_was_deleted_is_dropped_and_reported():
    """Live atoms with no bond between them is the only shape that reaches the writer's two
    bond-number lookups, and both keywords name the same bond here so one fixture drives both."""
    mol, store, _ = _read(_sgroup_record('1 SUP 0 ATOMS=(2 1 2) XBONDS=(1 1) '
                                         'CSTATE=(4 1 1.0 0.0 0.0) LABEL=Ph'))
    assert store.records[0].bonds and store.records[0].cstates, 'the fixture lost its references early'
    first, second = list(mol.atom_numbers)[:2]
    with mol.edit() as e:
        e.delete_bond(first, second)
    lines, log = emit_v3000(mol, store)
    assert sum('no longer exists' in x for x in log) == 2, log
    body = [x for x in lines if x.startswith('M  V30 1 SUP')]
    assert len(body) == 1 and 'XBONDS' not in body[0] and 'CSTATE' not in body[0], body
    assert 'ATOMS=(2 1 2)' in body[0], 'the atoms are still there; only the bond references went'


def test_the_one_v3000_record_in_the_corpus_agrees_with_chython2(root, v2_molecules):
    """Thin, and named as thin: ``implicit.sdf`` has exactly one V3000 record, and it is the only
    V3000 in this repository a second implementation has an opinion about."""
    path = root / 'test' / 'implicit.sdf'
    if not path.exists() or 'implicit.sdf' not in v2_molecules:
        return
    with path.open(encoding='utf8', errors='replace') as f:
        records = list(split_records(f))
    v2 = v2_molecules['implicit.sdf']
    assert len(v2) == len(records), 'chython 2 read a different number of records; not an oracle'
    checked = 0
    for n, record in enumerate(records):
        if sniff_version(record, []) != V3000_STAMP:
            continue
        # The one V3000 record is heteroaromatic, its nitrogen count needing a Kekule form; both
        # readers decline it, so the comparison below covers only where V2 commits.
        mol, _, _ = parse_v3000(record, []).build()
        m2 = v2[n]
        nums = list(m2)
        assert len([*mol.atom_numbers]) == len(nums)
        for sid, num in zip(mol.atom_numbers, nums):
            assert mol.element_of(sid) == m2.atom(num).atomic_number
            if m2.atom(num).implicit_hydrogens is not None:
                assert mol.implicit_h_of(sid) == m2.atom(num).implicit_hydrogens
        checked += 1
    assert checked, 'no V3000 record found in implicit.sdf; this test proved nothing'


def test_query_atom_message_does_not_promise_a_query_reader():
    """A message must not name an API that does not exist."""
    from pytest import raises

    from chython.formats.ctfile import UnsupportedCtfile, parse_v3000

    lines = ['query', '', '',
             '  0  0  0     0  0            999 V3000',
             'M  V30 BEGIN CTAB',
             'M  V30 COUNTS 1 0 0 0 0',
             'M  V30 BEGIN ATOM',
             'M  V30 1 A 0 0 0 0',
             'M  V30 END ATOM',
             'M  V30 END CTAB',
             'M  END']
    with raises(UnsupportedCtfile) as info:
        parse_v3000(lines, [])
    message = str(info.value)
    assert 'query reader' not in message, message
    assert 'query atom' in message, message


def test_a_non_utf8_title_is_emitted_as_the_byte_it_was_v3000():
    """Same invariant as V2000, tested separately because the call site is separate: the writer takes no
    loss on a name line, and the genuine-U+FFFD negative is gone with the replacement step that needed
    it.
    """
    from ....core import read_smiles

    # methane, so every implicit hydrogen count is known and an empty log means the title cost nothing
    mol = read_smiles('C')
    mol.set_title(b'caf\xe9')
    lines, log = emit_v3000(mol)
    assert not log, log
    assert lines[0].encode('utf8', 'surrogateescape') == b'caf\xe9'
