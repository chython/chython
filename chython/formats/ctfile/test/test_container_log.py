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
"""Every record a CTfile reader makes about a container it returns is on that container's `log`.

A reader keeps its `log=` list -- a parse has no container to write to until it has produced one, and a
framing decision belongs to the file rather than to any record -- but nothing needs to be passed for the
records to exist, and the stage is `'read'`.

The scoping is what these tests are for: a `log=` list is per record, so record 3's lines are on record
3's container and on nothing else.  A reader that folded one flat list onto every container it produced
would pass a "the log is not empty" assertion and be useless.
"""

from io import StringIO

from chython.formats.ctfile import RDFRead, SDFRead, mol, parse_record, parse_rxn, rxn


#: A record with no version stamp on its counts line: one repair, and the smallest one there is.
_NO_STAMP = ['no stamp', '', '',
             '  1  0  0  0  0  0            999',
             '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
             'M  END']

#: A record whose symbol column holds free text, which is read as a display label with the hydrogen
#: count left unknown -- three lines, none of them from the framing.
_LABELLED = ['labelled', '', '',
             '  2  1  0  0  0  0            999 V2000',
             '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
             '    1.0000    0.0000    0.0000 Xx  0  0  0  0  0  0  0  0  0  0  0  0',
             '  1  2  1  0  0  0  0',
             'M  END']

_CLEAN = ['clean', '', '',
          '  1  0  0  0  0  0            999 V2000',
          '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
          'M  END']


def _sdf(*records):
    return ''.join('\n'.join(x) + '\n$$$$\n' for x in records)


# ------------------------------------------------------------------------------- molfile and SDF

def test_a_repair_is_on_the_molecule_with_nothing_passed_in():
    """`mol()` called with no `log=` leaves the molecule as the only place the record can be."""
    molecule = mol('\n'.join(_NO_STAMP))
    assert [x.rule for x in molecule.log] == ['sdf:no-stamp']
    assert molecule.log[0].stage == 'read'
    assert molecule.log.repaired(), 'reading a stamp that is not there is a repair'


def test_the_facade_writes_to_the_molecule_too():
    assert any(x.rule == 'sdf:no-stamp' for x in mol('\n'.join(_NO_STAMP)).log)


def test_a_parse_records_on_the_molecule_and_on_the_callers_list():
    log = []
    molecule = parse_record(_LABELLED, log)
    assert [str(x) for x in log] == [str(x) for x in molecule.log]
    assert {x.stage for x in molecule.log} == {'read'}


def test_the_read_stage_does_not_leak_into_the_sgroup_view():
    """`sgroup_log` is `log` filtered to `edit:sgroup`; a reader's records must not appear in it."""
    molecule = parse_record(_LABELLED)
    assert molecule.log and molecule.sgroup_log == ()


def test_each_sdf_record_carries_its_own_records_and_no_others():
    with SDFRead(StringIO(_sdf(_CLEAN, _CLEAN, _LABELLED, _CLEAN))) as f:
        molecules = f.read()
    assert [x.title for x in molecules] == ['clean', 'clean', 'labelled', 'clean']
    assert [len(x.log) for x in molecules] == [0, 0, 1, 0], [list(x.log) for x in molecules]
    assert molecules[2].log[0].rule == 'v2000:symbol-as-label'


def test_a_shared_caller_list_is_not_folded_onto_a_later_molecule():
    """The list accumulates for the caller; the molecule gets its own record and no earlier one."""
    log = []
    first = parse_record(_LABELLED, log)
    second = parse_record(_NO_STAMP, log)
    assert len(log) == 2, log
    assert len(first.log) == 1 and [x.rule for x in second.log] == ['sdf:no-stamp']


def test_the_callers_list_is_still_filled_when_one_is_given():
    """The `log=` parameter did not go anywhere: a caller collecting many records still gets them."""
    log = []
    molecule = mol('\n'.join(_LABELLED), log=log)
    assert [str(x) for x in log] == [str(x) for x in molecule.log]


# ------------------------------------------------------------------------------------- reactions

_RXN = ['$RXN', 'labelled reactant', '', '',
        '  1  1',
        '$MOL', *_LABELLED,
        '$MOL', *_NO_STAMP]


def test_a_reaction_records_the_framing_and_names_the_component():
    """The record's own lines land bare; a component's are mirrored under the `subject` naming it.

    Both ends of the pair are asserted: `rxn.log.by_subject('reactants[0]')` and
    `rxn.reactants[0].log` answer the same question, which is the arrangement `ReactionContainer.log`
    documents and the only one in which a stable id in `atoms` names one container.
    """
    reaction = parse_rxn(_RXN)
    assert [x.rule for x in reaction.log.by_subject('reactants[0]')] == \
           [x.rule for x in reaction.reactants[0].log]
    assert [x.rule for x in reaction.log.by_subject('products[0]')] == ['sdf:no-stamp']
    assert {x.stage for x in reaction.log} == {'read'}


def test_no_component_record_reaches_the_reaction_twice():
    """The caller's list carries the role prefix, `reaction.log` carries the subject -- never both."""
    log = []
    reaction = parse_rxn(_RXN, log)
    assert any(str(x).startswith('component 1: ') for x in log), log
    assert not any('component 1: ' in str(x) for x in reaction.log), list(reaction.log)
    assert len(reaction.log) == len(log)


def test_the_component_keeps_its_own_records_unprefixed():
    reaction = parse_rxn(_RXN)
    assert [x.rule for x in reaction.reactants[0].log] == ['v2000:symbol-as-label']
    assert not any('component' in str(x) for x in reaction.reactants[0].log)


def test_the_facade_folds_its_own_repair_onto_the_reaction():
    """`rxn()` slices an RDfile paste back to `$RXN`; the skip count is about the reaction it returns."""
    reaction = rxn('$RFMT\n' + '\n'.join(_RXN))
    assert any(x.rule == 'ctfile:rfmt-skipped' for x in reaction.log), list(reaction.log)


_RXN_V3000 = ['$RXN V3000', 'v3000 with a labelled agent', '', '',
              'M  V30 COUNTS 1 1 1',
              'M  V30 BEGIN REACTANT',
              'M  V30 BEGIN CTAB',
              'M  V30 COUNTS 1 0 0 0 0',
              'M  V30 BEGIN ATOM',
              'M  V30 1 C 0 0 0 0',
              'M  V30 END ATOM',
              'M  V30 END CTAB',
              'M  V30 END REACTANT',
              'M  V30 BEGIN PRODUCT',
              'M  V30 BEGIN CTAB',
              'M  V30 COUNTS 1 0 0 0 0',
              'M  V30 BEGIN ATOM',
              'M  V30 1 O 0 0 0 0',
              'M  V30 END ATOM',
              'M  V30 END CTAB',
              'M  V30 END PRODUCT',
              'M  V30 BEGIN AGENT',
              'M  V30 BEGIN CTAB',
              'M  V30 COUNTS 1 0 0 0 0',
              'M  V30 BEGIN ATOM',
              'M  V30 1 Xx 0 0 0 0',
              'M  V30 END ATOM',
              'M  V30 END CTAB',
              'M  V30 END AGENT',
              'M  END']


def test_a_v3000_reaction_subjects_its_agent():
    """`_located`'s order is reactants, agents, products, so an agent is `agents[0]` and not `[2]`."""
    reaction = parse_rxn(_RXN_V3000)
    assert [x.rule for x in reaction.log.by_subject('agents[0]')][0] == 'v3000:atom-type-as-label'
    assert reaction.log.by_subject('reactants[0]') == []


def test_an_sd_record_holding_a_rxn_reports_on_the_reaction():
    """`SDFRead.log` is the container's own, so a rescue that stopped at the reader would be invisible."""
    with SDFRead(StringIO('\n'.join(_RXN) + '\n$$$$\n')) as f:
        reaction = f.read_record()
    assert any(x.rule == 'sdf:record-holds-rxn' for x in reaction.log), list(reaction.log)
    assert f.log is reaction.log


# ------------------------------------------------------------------------------------------ RDfile

_RDF = ('$RDFILE 1\n'
        '$DATM    09/06/26 12:00\n'
        '$MFMT\n' + '\n'.join(_CLEAN) + '\n'
        '$DTYPE name\n$DATUM first\n'
        '$MFMT\n' + '\n'.join(_LABELLED) + '\n'
        '$MIREG 42\n'
        '$RFMT\n' + '\n'.join(_RXN) + '\n')


def test_each_rdf_record_carries_its_own_records():
    with RDFRead(StringIO(_RDF)) as f:
        records = f.read()
    assert [len(x.log) for x in records] == [0, 2, 2], [list(x.log) for x in records]
    # The registry reference is a metadata line about record 2 and lands with the structure's own.
    assert [x.rule for x in records[1].log][-1] == 'rdf:registry-reference'
    assert records[0].meta == {'name': 'first'} and not records[0].log


def test_the_rdf_reaction_record_is_subjected():
    with RDFRead(StringIO(_RDF)) as f:
        reaction = f.read()[-1]
    assert reaction.log.by_subject('reactants[0]') and reaction.log.by_subject('products[0]')


def test_a_file_level_line_stays_on_the_file_log():
    """A stray line before the first record tag names no record, so no container may claim it."""
    with RDFRead(StringIO('$RDFILE 1\nstray\n$MFMT\n' + '\n'.join(_CLEAN) + '\n')) as f:
        molecule = f.read_record()
    assert [x.rule for x in f.file_log] == ['rdf:pre-record-line']
    assert not molecule.log
