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
"""The bidirectional conformance matrix: one cell per (feature, format, toolkit, direction, seed).

A format is a lingua franca or it is nothing, so both directions are asserted:

* INBOUND -- chython recovers the feature from a file a reference implementation wrote.  Four writers
  produce four dialects of the same format, and a feature chython reads out of Marvin's V2000 but not out
  of CDK's is a finding about a dialect, recorded as one.
* OUTBOUND -- a reference implementation recovers the feature from the file chython wrote.  This is the
  half that decides whether a collaborator can open a chython file at all.

The reference toolkit's own probe is read FIRST off its own output, in both directions.  Without that a
feature the writer never wrote would classify as the reader dropping it, and the failure would land on the
wrong side.

A cell states no expected outcome, so a reader that gets better never fails its own matrix.  Where the
format says nothing and implementations differ, the cell is listed in `DIALECT` with the reason -- a
statement of what each side can represent, never a judgement of an implementation.
"""
from typing import NamedTuple

from pytest import mark, skip

from .oracles import (CDK_JAR, CHYTHON_PROBES, DATA, FEATURE_FORMATS, MOLCONVERT, PROBES, READS,
                      SEEDS, TOOLKITS, Unrepresentable, WRITERS, container_of, is_reaction,
                      molecule_to_inchi, parse_once, read_chython, require, skeleton_of, versions,
                      warm, write_once)


#: A cell's four outcomes, plus the two that are not findings.
RECOVERED = 'recovered'      # the feature survived
ABSENT = 'absent'            # the file or the reader carried it, and the value differs
REFUSED = 'refused'          # the reader would not parse the file at all
NOT_WRITTEN = 'not-written'  # the writer stated no such feature in this file: nothing to recover
DECLARED = 'declared'        # the reader's API states no access to the feature: a declared gap


class Cell(NamedTuple):
    """One measurement.  `direction` is `'IN'` (others write, chython reads) or `'OUT'`."""
    feature: str
    fmt: str
    toolkit: str
    direction: str
    seed: str


def _id(cell):
    return f'{cell.feature}/{cell.fmt}/{cell.toolkit}/{cell.direction}/{cell.seed}'


def _cells(direction):
    """Every cell that exists: the feature can be stated in the format, both sides handle the format,
    and the seed declares the feature."""
    out = []
    for seed in SEEDS:
        for feature in sorted(seed.features):
            for fmt in sorted(FEATURE_FORMATS[feature]):
                if fmt not in READS['chython']:
                    continue
                for toolkit in TOOLKITS:
                    if fmt not in READS[toolkit]:
                        continue
                    writer = 'chython' if direction == 'OUT' else toolkit
                    if (writer, fmt) not in WRITERS:
                        continue
                    out.append(Cell(feature, fmt, toolkit, direction, seed.key))
    return tuple(out)


INBOUND = _cells('IN')
OUTBOUND = _cells('OUT')

#: The formats chython reads and does not write.  They carry atoms and coordinates and no bond block, so
#: the features measured on them are `element_multiset` and `coordinate_count`, and only the inbound half
#: of each cell exists.  `saturate()` is out of scope here: it is a separate explicitly invoked pass in
#: `chython.chemistry`, and importing that from `formats/test/` would reverse the dependency direction.
READ_ONLY = frozenset({'xyz', 'pdb', 'mmcif', 'mol2'})

#: The read-only cells, for the report.  The day a writer for one of these lands, `_cells('OUT')` grows
#: the outbound half by itself and the ratchet below is what says so.
INBOUND_ONLY = tuple(cell for cell in INBOUND if cell.fmt in READ_ONLY)

#: Cells where the format cannot state the feature, or states it in more than one spelling, and the
#: implementations differ.  A claim is a STATEMENT OF CAPABILITY -- what each side can represent and where
#: the format is silent -- never a judgement of an implementation.  A claimed cell that starts
#: round-tripping fails, so no claim outlives the behaviour it describes.
DIALECT = {
    Cell('stereo_groups', 'v2000', 'cdk', 'IN', 'hexanediol'):
        'CTfile V2000 has no enhanced-stereo block: relative stereochemistry is stated by the chiral '
        'flag on the counts line, and CDK 2.12 writes that flag as 0 with the parities in the atom '
        'block.  CDK reads its own file back as one racemic collection over those parities; chython '
        'reports a collection only where a file states one, which in CTfile means a V3000 '
        'MDLV30/STERAC or STEREL block.  The two formats differ in what they can state, and the readers '
        'differ in what they infer from the flag.',
    Cell('stereo_groups', 'sdf', 'cdk', 'IN', 'hexanediol'):
        'The same V2000 CTAB inside an SD frame, so the same chiral-flag reading: see the v2000/cdk '
        'claim.',
    Cell('stereo_groups', 'v2000', 'indigo', 'IN', 'hexanediol'):
        'Indigo 1.45 writes the same V2000 shape -- chiral flag 0, wedges in the bond block -- and reads '
        'it back as `stereocenterType() == 3` (AND) on both centres.  chython reads the parities and '
        'reports no collection, V2000 having no block that states one.',
    Cell('dat_sgroup', 'v3000', 'cdk', 'OUT', 'aspirin'):
        'CDK 2.12\'s MDLV3000Reader logs `Skipping unrecognized SGROUP type: DAT` for the group its own '
        'V2000 reader reads, so the same data S-group survives chython\'s V2000 and not its V3000.',
    Cell('implicit_h', 'v2000', 'cdk', 'OUT', 'pyrrole'):
        'Only the pyrrole nitrogen, and not the bond order: CDK 2.12 derives 1 hydrogen for each of the '
        'four aromatic carbons of chython\'s order-4 ring and 0 for the N.  chython states that N-H in a '
        '`MRV_IMPLICIT_H` data S-group (`M  SDT 1 MRV_IMPLICIT_H` / `M  SED 1 IMPL_H1`), which RDKit '
        '2026.03.4, Indigo 1.45 and `molconvert` all consume.  CDK 2.12\'s MDLV2000Reader reads the group '
        '-- it comes back as `(MRV_IMPLICIT_H, IMPL_H1, 1)` from `getSgroups()` -- and states the count '
        'from its own configuration rather than from the group, that field being a ChemAxon convention '
        'and not a CTfile one.  A Kekule form, written after `kekule()`, carries the count to CDK '
        'without it.',
    Cell('implicit_h', 'v3000', 'cdk', 'OUT', 'pyrrole'):
        'The same nitrogen through V3000, where CDK 2.12 does not read the group at all: '
        '`Skipping unrecognized SGROUP type: DAT`.  See the dat_sgroup/v3000/cdk claim.',
    Cell('implicit_h', 'sdf', 'cdk', 'OUT', 'pyrrole'):
        'The same V2000 CTAB inside an SD frame, same reading as the v2000/cdk claim.',
    Cell('sd_fields', 'cml', 'cdk', 'OUT', 'glycerol'):
        'A molecule property has two spellings in CML and the two writers use one each.  chython writes '
        'the `<propertyList>` of `<property>`/`<scalar>` that `molconvert cml` writes and Marvin reads '
        'back; CDK 2.12\'s CMLCoreModule logs `Ignoring scalar: '
        '/cml/molecule/propertyList/property/scalar/` for it and writes its own as a `<scalar '
        'dictRef="cdk:molecularProperty" title="BATCH_ID">` straight under `<molecule>`.  chython reads '
        'that second spelling -- the inbound cell recovers both fields -- so the difference is in what '
        'the two write, CML stating no one form for a data field.',
    Cell('sd_fields', 'cml', 'indigo', 'OUT', 'glycerol'):
        'Indigo 1.45 states properties on a molecule (`setProperty`/`iterateProperties`) and its CML '
        'reader and writer state none: `loadMolecule` on chython\'s `<propertyList>` yields a molecule '
        'with no property, and `.cml()` on a molecule carrying one writes no property element, which is '
        'why the inbound cell has nothing to recover and skips.',
}

#: Cells that are chython's own reading or writing to close, with the measurement that says what closes
#: them.  Kept as claims rather than as a report file: a claim is rechecked on every run and fails the day
#: it is fixed, which is how the finding reaches whoever fixed it.
HANDOFF = {}

#: Every claim, whichever kind.  A cell may not be in both.
CLAIMS = {**DIALECT, **HANDOFF}


def _stated(cache, toolkit, fmt, feature, text):
    """One side's own probe value for a file, plus the evidence, or a reason it has none."""
    parsed = parse_once(cache, toolkit, fmt, text)
    if parsed.error is not None:
        return None, parsed.evidence, parsed.error
    try:
        return PROBES[toolkit][feature](parsed.obj), parsed.evidence, None
    except Unrepresentable as gap:
        return None, parsed.evidence, f'{DECLARED}: {gap}'


def classify_inbound(cell, written, parses):
    """What chython does with a reference writer's file: `(outcome, evidence)`."""
    try:
        text = write_once(written, cell.toolkit, cell.fmt, _seed(cell))
    except Exception as error:
        return NOT_WRITTEN, f'{cell.toolkit} has no {cell.fmt} spelling for {cell.seed}: {error}'
    stated, call, error = _stated(parses, cell.toolkit, cell.fmt, cell.feature, text)
    if error is not None:
        return NOT_WRITTEN, f'{cell.toolkit} does not state {cell.feature} for its own {cell.fmt}: {error}'
    if not stated:
        return NOT_WRITTEN, f'{cell.toolkit} wrote {cell.fmt} without the feature; nothing to recover'
    log = []
    try:
        got = read_chython(text, cell.fmt, log=log)
    except Exception as error:
        return REFUSED, (f'chython refused: {type(error).__name__}: {error}\n'
                         f'call: {call}\n--- bytes {cell.toolkit} wrote ---\n{text}')
    mine = CHYTHON_PROBES[cell.feature](got)
    if mine == stated:
        return RECOVERED, ''
    return ABSENT, (f'{cell.toolkit} stated {stated!r}, chython read {mine!r}\n'
                    f'chython log: {log}\ncall: {call}\n--- bytes {cell.toolkit} wrote ---\n{text}')


def classify_outbound(cell, written, parses):
    """What a reference reader recovers from chython's file: `(outcome, evidence)`."""
    try:
        text = write_once(written, 'chython', cell.fmt, _seed(cell))
    except Exception as error:
        return NOT_WRITTEN, f'chython has no {cell.fmt} spelling for {cell.seed}: {error}'
    mine = CHYTHON_PROBES[cell.feature](read_chython(text, cell.fmt))
    if not mine:
        return NOT_WRITTEN, f'chython wrote {cell.fmt} without the feature; nothing to recover'
    back, call, error = _stated(parses, cell.toolkit, cell.fmt, cell.feature, text)
    if error is not None:
        if error.startswith(DECLARED):
            return DECLARED, f'{error}\ncall: {call}'
        return REFUSED, (f'{cell.toolkit} refused the file: {error}\n'
                         f'call: {call}\n--- bytes chython wrote ---\n{text}')
    if back == mine:
        return RECOVERED, ''
    return ABSENT, (f'chython wrote {mine!r}, {cell.toolkit} read {back!r}\n'
                    f'call: {call}\n--- bytes chython wrote ---\n{text}')


def _seed(cell):
    for seed in SEEDS:
        if seed.key == cell.seed:
            return seed
    raise KeyError(cell.seed)


def _assert(cell, outcome, evidence, versions):
    if outcome in (NOT_WRITTEN, DECLARED):
        skip(evidence)
    if cell in CLAIMS:
        assert outcome != RECOVERED, (
            f'{_id(cell)} is claimed as {CLAIMS[cell]!r} and now round-trips.  Delete the claim -- a '
            f'stale claim is a cell nobody measures.')
        return
    assert outcome == RECOVERED, (
        f'{_id(cell)} with {cell.toolkit} at {versions[cell.toolkit]}: {outcome}\n{evidence}\n'
        f'Either fix the reader or the writer, or state the cell in DIALECT (the format cannot say it, or '
        f'says it two ways) or in HANDOFF (chython\'s own writing, with the measurement that closes it).')


@mark.parametrize('cell', INBOUND, ids=[_id(c) for c in INBOUND])
def test_chython_recovers_the_feature_a_reference_writer_stated(cell, written, parses,
                                                                oracle_versions):
    require(cell.toolkit)
    _assert(cell, *classify_inbound(cell, written, parses), oracle_versions)


@mark.parametrize('cell', OUTBOUND, ids=[_id(c) for c in OUTBOUND])
def test_a_reference_reader_recovers_the_feature_chython_wrote(cell, written, parses, oracle_versions):
    require(cell.toolkit)
    _assert(cell, *classify_outbound(cell, written, parses), oracle_versions)


@mark.parametrize('cell', OUTBOUND, ids=[_id(c) for c in OUTBOUND])
def test_a_reference_reader_reads_chython_as_the_same_molecule(cell, written, parses, oracle_versions):
    """Constitution before feature: formula, `/c` and `/h` from the reader's own InChI must match
    chython's.  A feature can survive a write the molecule does not, and that is the worse loss."""
    require(cell.toolkit)
    if is_reaction(_seed(cell)):
        skip('a reaction record has no single InChI; the feature assertion carries this cell')
    try:
        text = write_once(written, 'chython', cell.fmt, _seed(cell))
    except Exception as error:
        skip(str(error))
    parsed = parse_once(parses, cell.toolkit, cell.fmt, text)
    if parsed.error is not None:
        skip(f'{cell.toolkit} refused the file; the feature assertion records that: {parsed.error}')
    if parsed.skeleton is None:
        skip(f'{cell.toolkit} states no InChI for {cell.fmt}; the feature assertion still runs')
    mine = skeleton_of(molecule_to_inchi(container_of(read_chython(text, cell.fmt))))
    assert parsed.skeleton == mine, (
        f'{_id(cell)} with {cell.toolkit} at {oracle_versions[cell.toolkit]}: chython wrote {mine}, '
        f'{cell.toolkit} read {parsed.skeleton}\ncall: {parsed.evidence}\n'
        f'--- bytes chython wrote ---\n{text}')


def report_rows():
    """Every cell measured once: `(direction, format, feature, seed, toolkit, outcome, note)`.

    The note is the claim where one is stated and the first line of the evidence otherwise, so a row is
    readable without the bytes behind it.
    """
    written, parses = {}, {}
    warm(written, parses)
    rows = []
    for cells, classify in ((INBOUND, classify_inbound), (OUTBOUND, classify_outbound)):
        for cell in cells:
            if versions()[cell.toolkit] is None:
                outcome, evidence = 'skipped', f'{cell.toolkit} absent'
            else:
                outcome, evidence = classify(cell, written, parses)
            note = CLAIMS.get(cell, evidence.split('\n')[0])
            rows.append((cell.direction, cell.fmt, cell.feature, cell.seed, cell.toolkit, outcome,
                         note))
    return rows


def print_matrix():
    """`python -m chython.formats.test.oracles matrix` -- the measurement, as markdown, on stdout.

    This command IS the harness behind any conformance claim about formats: a coverage claim ships with
    the harness that produced it or it is not made.  Nothing here is committed as a report file, because
    a report file is a claim nobody remeasured.
    """
    from datetime import date

    print(f'# chython format conformance, measured {date.today().isoformat()}\n')
    for name, version in versions().items():
        print(f'{name:8} {version if version is not None else "absent"}')
    if versions()['cdk'] is not None:
        print(f'{"":8} jar {CDK_JAR}')
    if versions()['marvin'] is not None:
        print(f'{"":8} bin {MOLCONVERT}')
    print()

    rows = report_rows()
    outcomes = {}
    for row in rows:
        outcomes[row[5]] = outcomes.get(row[5], 0) + 1
    print('| dir | format | feature | seed | ' + ' | '.join(TOOLKITS) + ' |')
    print('|---' * (4 + len(TOOLKITS)) + '|')
    by_key = {}
    for direction, fmt, feature, seed, toolkit, outcome, _ in rows:
        by_key.setdefault((direction, fmt, feature, seed), {})[toolkit] = outcome
    for (direction, fmt, feature, seed), cells in by_key.items():
        marks = ' | '.join(_MARK.get(cells.get(t), '') for t in TOOLKITS)
        print(f'| {direction} | {fmt} | {feature} | {seed} | {marks} |')
    print()
    for outcome, count in sorted(outcomes.items()):
        print(f'{outcome:12} {count}')
    print(f'\n{len(DIALECT)} dialect claims, {len(HANDOFF)} handed to chython\'s writers:\n')
    for cell, reason in CLAIMS.items():
        kind = 'DIALECT' if cell in DIALECT else 'HANDOFF'
        print(f'* **{_id(cell)}** [{kind}] {reason}')
    return 0


#: The one-character spelling of each outcome in the table.  A blank cell does not exist.
_MARK = {RECOVERED: 'ok', ABSENT: 'differs', REFUSED: 'refused', NOT_WRITTEN: '-', DECLARED: 'no api',
         'skipped': 'skipped'}


def test_the_sd_field_fixture_states_a_value_on_two_lines():
    """The control under the eight `sd_fields` cells.

    A field whose second line is dropped is an `ABSENT` outcome and not a `REFUSED` one -- the frame
    parses either way -- so the only thing that makes those cells measure it is the fixture.  Were it
    quietly reduced to one line, all eight would still pass and none would measure anything.
    """
    from .oracles import SD_FIELDS, seed_of, write_chython

    fields = dict(SD_FIELDS)
    assert '\n' in fields['NOTES'], 'the multi-line SD value is what a lost second line shows up in'
    frame = write_chython(seed_of('glycerol'), 'sdf').splitlines()
    assert ['first line', 'second line'] == frame[frame.index('>  <NOTES>') + 1:][:2]


# ---------------------------------------------------------------------------------- atom lists (query)
# The one axis with no matrix cell.  An atom list is a query construct: CTfile states it as `L` in the
# atom-line symbol field plus `M  ALS`, or as `[C,N]` in a V3000 atom line, and a MoleculeContainer has
# nowhere to put it, in either direction.  So instead of a cell the divergence is recorded here -- the
# four readers hand back four different kinds of object for the same bytes, and none of the four is the
# shape this test asserts.

#: The query the fixture states, as SMARTS.  Written by RDKit rather than by hand: the `M  ALS` tail is
#: fixed-width and inventing its columns manufactures findings about the fixture.
ATOM_LIST_SMARTS = '[#6,#7]CO'

#: What each reader makes of that file, MEASURED at the versions in the report header, `L`/`M  ALS` for
#: v2000 and `[C,N]` for v3000.  A record and not a claim: a reader that changes fails here and is
#: remeasured.  The vocabulary is the KIND of answer -- both list members, an opaque query object, a
#: label, or a refusal -- because that is what differs, and CTfile leaves how to model a query open.
ATOM_LIST_KINDS = {
    ('chython', 'v2000'): 'refused: UnsupportedCtfile, naming a query reader',
    ('chython', 'v3000'): 'refused: UnsupportedCtfile, naming a query reader',
    ('rdkit', 'v2000'): 'MolFromMolBlock: a mol whose atom states both members',
    ('rdkit', 'v3000'): 'MolFromMolBlock: a mol whose atom states both members',
    ('indigo', 'v2000'): 'loadMolecule: refused, atom lists being for queries; '
                         'loadQueryMolecule: both members',
    ('indigo', 'v3000'): 'loadMolecule: refused, atom lists being for queries; '
                         'loadQueryMolecule: both members',
    ('cdk', 'v2000'): 'MDLV2000Reader: a QueryAtom, symbol None',
    ('cdk', 'v3000'): 'MDLV3000Reader: a PseudoAtom, labelled R',
    ('marvin', 'v2000'): 'molconvert -g smarts: both members',
    ('marvin', 'v3000'): 'molconvert -g smarts: both members',
}


def _atom_list_fixture(version):
    """The query written as a CTAB by RDKit, both versions."""
    from rdkit import Chem
    from rdkit.Chem import rdDepictor

    query = Chem.MolFromSmarts(ATOM_LIST_SMARTS)
    rdDepictor.Compute2DCoords(query)
    return Chem.MolToMolBlock(query) if version == 'v2000' else Chem.MolToV3KMolBlock(query)


def _atom_list_kind(toolkit, version, text):
    """The kind of answer *toolkit* gives for the fixture, in the vocabulary of `ATOM_LIST_KINDS`."""
    if toolkit == 'chython':
        from ..ctfile import mol
        from ..ctfile._errors import UnsupportedCtfile
        try:
            mol(text)
        except UnsupportedCtfile as refusal:
            assert 'query' in str(refusal), refusal
            return 'refused: UnsupportedCtfile, naming a query reader'
        return 'read as a structure CTAB'
    if toolkit == 'rdkit':
        from rdkit import Chem
        back = Chem.MolFromMolBlock(text)
        if back is None:
            return 'MolFromMolBlock: refused'
        smarts = Chem.MolToSmarts(back)
        both = '#6' in smarts and '#7' in smarts
        return f'MolFromMolBlock: a mol whose atom states {"both members" if both else "one element"}'
    if toolkit == 'indigo':
        from indigo import Indigo
        session = Indigo()
        try:
            session.loadMolecule(text)
            first = 'loadMolecule: read as a structure'
        except Exception as refusal:
            assert 'quer' in str(refusal), refusal
            first = 'loadMolecule: refused, atom lists being for queries'
        smarts = session.loadQueryMolecule(text).smarts()
        both = '#6' in smarts and '#7' in smarts
        return f'{first}; loadQueryMolecule: {"both members" if both else "one element"}'
    if toolkit == 'cdk':
        from .oracles import cdk
        JClass = cdk()
        name = f'org.openscience.cdk.io.MDL{version.upper()}Reader'
        builder = JClass('org.openscience.cdk.silent.SilentChemObjectBuilder').getInstance()
        reader = JClass(name)(JClass('java.io.StringReader')(text))
        atom = next(iter(reader.read(builder.newAtomContainer()).atoms()))
        kind, symbol = type(atom).__name__.rsplit('.', 1)[-1], str(atom.getSymbol())
        return f'MDL{version.upper()}Reader: a {kind}, ' + \
               (f'labelled {symbol}' if kind == 'PseudoAtom' else f'symbol {symbol}')
    from .oracles import molconvert
    code, out, err, _ = molconvert(text, 'smarts')
    if code or not out.strip():
        return f'molconvert -g smarts: refused, {(err or out).strip()[:80]}'
    both = '#6' in out and '#7' in out
    return f'molconvert -g smarts: {"both members" if both else "one element"}'


@mark.parametrize('version', ('v2000', 'v3000'))
@mark.parametrize('toolkit', ('chython',) + TOOLKITS)
def test_an_atom_list_is_modelled_differently_by_every_reader(toolkit, version):
    """Four object kinds for one file, recorded rather than reconciled.

    CTfile states the atom list; it does not state what a reader builds from one, and these four build a
    query mol, a query-only door, a query atom, a pseudo atom and a refusal.  So nothing here asserts a
    shape -- only that each reader still answers the way the record says, which is what makes a change
    show up as something to remeasure instead of passing unnoticed.
    """
    require('rdkit')  # the fixture's writer, whichever toolkit reads it
    if toolkit != 'chython':
        require(toolkit)
    text = _atom_list_fixture(version)
    assert ATOM_LIST_KINDS[toolkit, version] == _atom_list_kind(toolkit, version, text), \
        f'{toolkit} answers the {version} atom list differently now; remeasure the record.\n{text}'


# ------------------------------------------------------------------------------------ reaction records
# What a reader recovers from a reaction file, not what a mapper produces: Stream 2 changes what the
# reactor writes and nothing here depends on it.

#: The counts line chython's V2000 RXN states for the three-sided seed, MEASURED.  The spec's counts
#: line officially carries two fields; a third with the agent count is what the ecosystem writes, and
#: `emit_rxn` logs the convention.
V2000_AGENT_COUNTS_LINE = '  2  2  1'

#: What each reader reports for that file, MEASURED at the versions in the report header.  All four read
#: the third field as an agent count, so this is a record and not a claim -- a reader that changes fails
#: here and is remeasured rather than being assumed.
V2000_AGENT_SIDES = {'rdkit': (2, 1, 2), 'indigo': (2, 1, 2), 'cdk': (2, 1, 2), 'marvin': (2, 1, 2)}

#: The committed RDF corpus as `(name, records, reactions)`.  Among the four, Marvin writes RDF and none
#: of the other three reads it, so these fixtures are the rest of the RDF inbound half.  `MR.rdf` mixes
#: V2000 and V3000 CTABs in one file, which is why the import side of `mol()` sniffs the version per CTAB.
RDF_FIXTURES = (('MR', 4, 2), ('ions', 1, 1), ('standardize', 6, 6), ('reaction_centerslist', 2, 2))


def _esterification_v2000():
    from ..ctfile import rxn

    from .oracles import build_chython, seed_of

    return rxn(build_chython(seed_of('esterification')), version=2000)


def test_the_reaction_seed_states_two_reactants_one_agent_and_two_products():
    """The baseline the reaction cells rest on, and the order `molecules()` yields, which is
    load-bearing: reactants, then agents, then products."""
    from .oracles import build_chython, seed_of

    reaction = build_chython(seed_of('esterification'))
    assert CHYTHON_PROBES['reaction_sides'](reaction) == (2, 1, 2)
    assert [id(x) for x in reaction.molecules()] == \
           [id(x) for side in (reaction.reactants, reaction.agents, reaction.products) for x in side]


@mark.parametrize('toolkit', TOOLKITS)
def test_the_third_field_of_a_v2000_rxn_counts_line_is_read_as_the_agent_count(toolkit):
    """The agent question, measured instead of assumed.

    V2000's counts line officially carries two fields, so a file with three is a place the spec is
    silent.  Every one of these four reads the third as an agent count, and the record above is what
    fails if that changes.
    """
    require(toolkit)
    text = _esterification_v2000()
    assert V2000_AGENT_COUNTS_LINE in text.splitlines(), 'chython no longer states the third field'
    parsed = parse_once({}, toolkit, 'rxn', text)
    assert parsed.error is None, f'{toolkit} refused the file: {parsed.error}\ncall: {parsed.evidence}'
    assert PROBES[toolkit]['reaction_sides'](parsed.obj) == V2000_AGENT_SIDES[toolkit], \
        f'{toolkit} at {versions()[toolkit]} now reads the counts line differently\n{text}'


@mark.parametrize('name,records,reactions', RDF_FIXTURES)
def test_chython_reads_every_record_of_a_committed_rdf(name, records, reactions):
    """An RDfile record count and how many of them are reactions -- `MR.rdf` carries both kinds."""
    from io import StringIO

    from ..ctfile import RDFRead

    with RDFRead(StringIO((DATA / f'{name}.rdf').read_text(encoding='utf-8'))) as handle:
        read = list(handle)
    assert len(read) == records
    assert sum(1 for x in read if hasattr(x, 'reactants')) == reactions


@mark.parametrize('toolkit', TOOLKITS)
@mark.parametrize('name', [x[0] for x in RDF_FIXTURES])
def test_a_reference_reader_recovers_the_sides_and_mapping_of_a_re_emitted_rdf_record(name, toolkit):
    """The outbound half for the RDF corpus: chython reads a committed record and writes it as RXN, and
    the reference reader must report the same sides and the same map numbers.

    The first reaction of each fixture, one file per fixture: the point is the dialects the corpus
    carries -- `MR.rdf`'s two CTAB versions among them -- not the record count.
    """
    from io import StringIO

    from ..ctfile import RDFRead, rxn

    require(toolkit)
    with RDFRead(StringIO((DATA / f'{name}.rdf').read_text(encoding='utf-8'))) as handle:
        reaction = next(x for x in handle if hasattr(x, 'reactants'))
    text = rxn(reaction, version=3000)
    parsed = parse_once({}, toolkit, 'rxn', text)
    assert parsed.error is None, (f'{toolkit} refused the RXN of {name}.rdf: {parsed.error}\n'
                                  f'call: {parsed.evidence}\n{text}')
    for feature in ('reaction_sides', 'mapping'):
        try:
            back = PROBES[toolkit][feature](parsed.obj)
        except Unrepresentable as gap:
            skip(f'{toolkit} states no access to {feature}: {gap}')
        assert back == CHYTHON_PROBES[feature](reaction), (
            f'{name}.rdf {feature} with {toolkit} at {versions()[toolkit]}\ncall: {parsed.evidence}\n'
            f'--- bytes chython wrote ---\n{text}')


#: PDBx/mmCIF fixtures, one record per `_atom_site.pdbx_PDB_model_num`.  The damaged and unterminated
#: fixtures are not here: they measure refusal and repair, which `test_pdb_builder.py` owns.
MMCIF_FIXTURES = ('mmcif_dipeptide', 'mmcif_disulfide', 'mmcif_ligand_water', 'mmcif_no_bonds',
                  'mmcif_altloc', 'mmcif_two_models')


def atom_site_rows(text):
    """The `_atom_site` loop of a PDBx/mmCIF file, as `{model: [row, ...]}`.

    Ten lines rather than an oracle call, because none of the four reads PDBx/mmCIF: `rdkit.Chem` states
    no mmCIF reader, CDK 2.12's `CIFReader` returns no container for these files (its tags are the
    crystallographic ones), and `molconvert` answers `Cannot recognize format by any of the supported
    molecular file format recognizers`.  So the reference side of an mmCIF cell is the file's own bytes,
    counted here by a different code path than the one under test.
    """
    from shlex import split

    models, tags, in_loop = {}, [], False
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith('#'):
            continue
        if stripped == 'loop_':
            tags, in_loop = [], True
        elif in_loop and stripped.startswith('_'):
            tags.append(stripped.split()[0])
        elif in_loop and tags and tags[0].startswith('_atom_site.'):
            row = dict(zip(tags, split(stripped)))
            models.setdefault(row.get('_atom_site.pdbx_PDB_model_num', '1'), []).append(row)
        elif stripped.startswith('_'):
            in_loop = False
    return models


@mark.parametrize('name', MMCIF_FIXTURES)
def test_chython_reads_every_atom_site_row_of_a_committed_mmcif_fixture(name):
    """The inbound mmCIF cells, fixture-backed rather than oracle-generated.

    A row of the loop is an atom of the record and a `Cartn_x`/`Cartn_y`/`Cartn_z` triple is a
    coordinate, so a dropped row or a dropped model shows up as a count.  Alternate locations are all
    here: the reader stores the rows and `build_molecule` is where one location is chosen.
    """
    from chython.formats.pdb import read_mmcif

    text = (DATA / f'{name}.cif').read_text(encoding='utf-8')
    models = atom_site_rows(text)
    records = list(read_mmcif(text.splitlines()))
    assert len(records) == len(models), (f'{name}.cif states {sorted(models)} models in '
                                         f'`_atom_site.pdbx_PDB_model_num`, chython read {len(records)}')
    for (model, rows), record in zip(sorted(models.items()), records):
        elements = tuple(sorted(row['_atom_site.type_symbol'] for row in rows))
        assert CHYTHON_PROBES['element_multiset'](record) == elements, \
            f'{name}.cif model {model}: the `_atom_site` loop states {elements}'
        assert CHYTHON_PROBES['coordinate_count'](record) == len(rows), \
            f'{name}.cif model {model}: the loop states {len(rows)} rows with Cartn_x/y/z'


def test_the_coordinate_formats_are_inbound_only_until_a_writer_lands():
    """XYZ, PDB, PDBx/mmCIF and MOL2 are read and not written, so half of each cell exists.

    Stated as a ratchet and not as a comment: when a writer for one of them lands, this fails, and the
    outbound cells the matrix grows at the same moment are then measured rather than assumed.
    """
    from .oracles import WRITERS as writers

    have = sorted(fmt for fmt in READ_ONLY if ('chython', fmt) in writers)
    assert not have, f'chython now writes {have}: the outbound half of those cells is live'
    assert INBOUND_ONLY, 'the read-only formats lost their inbound cells'
    assert not [c for c in OUTBOUND if c.fmt in READ_ONLY]


def test_every_dialect_claim_names_a_cell_the_matrix_measures():
    """A claim cannot be added for a cell that does not exist -- that is a claim nobody rechecks."""
    live = set(INBOUND) | set(OUTBOUND)
    assert not set(CLAIMS) - live, sorted(_id(c) for c in set(CLAIMS) - live)
    assert all(CLAIMS.values()), 'every claim needs a written reason'
    both = set(DIALECT) & set(HANDOFF)
    assert not both, f'a cell is one kind of claim or the other: {sorted(_id(c) for c in both)}'
