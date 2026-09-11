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
"""The writer against real molecules, with TWO independent oracles and a corpus-wide invariance sweep.

Every other writer test states a rule and shows one molecule obeying it.  This one states nothing and
asks 5000 public records whether the string that came out means what went in.  The round trip here
goes out through this writer and back through somebody ELSE's parser, which is stronger than a
self-round-trip, not weaker: a writer and a reader that share a misunderstanding agree with each
other and with nobody else.  The closed loop our own reader makes possible lives in
`test_smiles_roundtrip.py` -- a complement to this file and not a replacement, for exactly the reason
in the sentence before this one.

THE ORACLES, in order of how much they are trusted:

1. **RDKit canonical SMILES.**  `MolToSmiles(MolFromSmiles(input)) == MolToSmiles(MolFromSmiles(our
   output))`.  Independent of chython entirely, insensitive to the Kekule-versus-aromatic spelling
   (RDKit re-perceives), and its canonical form IS a fixpoint.  This is the gate.
2. **chython 2's canonical string.**  Secondary, because it is NOT automorphism-invariant: measured on
   this corpus, 11 records where V2 hands back a different string for the same molecule presented in a
   different atom order, and 4 of those are cases where V2 does not even reach a fixpoint.  V2 is an
   oracle, not a contract, and this is the file that measures how good an oracle it is.  It is reached
   as a CO-PROCESS and never imported -- see `oracle.py` and the note above the `oracle` fixture --
   so the core's suite does not require chython 2 to be in this tree, and skips when it is absent.
3. **Ourselves, swept.**  Every record written from six different creation orders; the strings must be
   byte-identical.  This is the epic's central claim and the only one no external tool can check.

CORPORA.  `first_5K.smi` is 4999 public NCI records shipped inside RDKit, so nothing is fetched and
there is no chance of two files with one name; plus the repo's own `test/arenes.sdf` and
`test/heterocycles_charges.smi`, which are dense in the aromatic and charged-heterocycle cases the
fixtures reach for one at a time.  **A corpus is read from `test/` or from a shipped RDKit data file and
from nowhere else** -- a loose file in the repo root is scratch and is never wired in here.  Every
corpus is optional at import time: a missing one skips, it does not fail.

THE NCI CORPUS CARRIES NO STEREO AT ALL -- measured, not assumed: zero `@` and zero `/` in 4999
records.  So everything above it is a CONSTITUTION differential and nothing in it can fail on a
configuration, which is why the second half of this file exists and reaches for two other corpora:
`test/stereo.sdf` (300 records, 964 tetrahedral signs and 11 allene signs) for the tetrahedral and
allene half, and RDKit's own NIBR PubChem example table for the cis/trans half, the repo having no
double-bond geometry anywhere in `test/`.  The stereo half has its own bridge, its own oracle
ordering, and its own measured failure -- see `bridge_stereo` and the mirror-automorphism tests.
"""
from csv import DictReader
from pathlib import Path
from random import Random

from pytest import fixture, mark, skip

from chython.core import H_UNKNOWN, MoleculeContainer
from chython.core._core import smw_traversal, write_smiles
from . import oracle as oracle_module


def _rdkit_root():
    """Where RDKit is installed, or None.  The corpora live INSIDE the package, so this is also the
    check that they exist -- and it is a path nobody has to edit to run the file on another machine."""
    try:
        import rdkit
    except ImportError:
        return None
    return Path(rdkit.__file__).resolve().parent


RDKIT = _rdkit_root()
NCI = None if RDKIT is None else RDKIT / 'Data' / 'NCI' / 'first_5K.smi'
# 444 rows of substructure filters, each with up to five PubChem example structures -- a public table
# whose examples are dense in double-bond geometry, which is the one thing `test/` has none of.
NIBR = None if RDKIT is None else (RDKIT / 'Contrib' / 'NIBRSubstructureFilters' /
                                   'SubstructureFilter_HitTriaging_wPubChemExamples.csv')
REPO = Path(__file__).resolve().parents[3] / 'test'


# THE ORACLE IS A CO-PROCESS AND NOT AN IMPORT.  chython 2 runs in its own interpreter against an
# INSTALLED copy -- see `oracle.py` -- so this file, and therefore the core's whole suite, does not
# have to exist in the same tree as the library it is a differential against.  Nothing here
# reconstructs a V2 molecule: a `Record` is what V2 said about one, plus a handle for asking more.
#
# Both readers stay on V2's side by their UNSHADOWED module paths.  `chython.SDFRead` in this tree is
# the chython 3 CTfile reader (`formats/__init__.py`), so following that name would replace the oracle
# with the subject: the stereo corpus is read so that V2 can be ASKED for a sign the writer is then
# checked against, and a V3 reader answers with the arena's own parity bytes instead, which is the
# thing under test.  The failure that prevents is a green suite comparing the writer to itself.

@fixture(scope='module')
def oracle():
    """One chython 2 interpreter for the whole file. Skips when it is not provisioned."""
    live = oracle_module.session()
    yield live
    live.close()


def rdkit_canonical():
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog('rdApp.*')

    def canonical(text):
        mol = Chem.MolFromSmiles(text, sanitize=True)
        return None if mol is None else Chem.MolToSmiles(mol)
    return canonical


def _bridge_atoms(source, order=None):
    """`(V3 molecule, {V2 atom number: stable id})`, constitution only.

    `implicit_hydrogens is None` in V2 means the valence model could not derive a count -- which is
    `H_UNKNOWN` and NOT zero, and mapping it to zero would hide the only class of difference this file
    finds.  `order` permutes the creation order, which is what the invariance sweep needs.
    """
    mol = MoleculeContainer()
    ids = {}
    atoms = source.atom_rows
    for k in (range(len(atoms)) if order is None else order):
        atom = atoms[k]
        h = atom['h']
        ids[atom['n']] = mol.add_atom(atom['z'], implicit_h=H_UNKNOWN if h is None else h,
                                      charge=atom['charge'], radical=atom['radical'],
                                      isotope=atom['isotope'])
    for bond in source.bond_rows:
        mol.add_bond(ids[bond['n']], ids[bond['m']], bond['order'])
    return mol, ids


def bridge(source, order=None):
    """The constitution bridge the whole first half of this file runs on."""
    return _bridge_atoms(source, order)[0]


@fixture(scope='module')
def nci(oracle):
    """`(text, Record)` per readable NCI record.  Parsed once; V2's parser is the slow part."""
    if NCI is None or not NCI.is_file():
        skip('the RDKit NCI corpus is not installed (looked for %s)' % NCI)
    with NCI.open(encoding='utf-8') as f:
        texts = [line.split()[0] for line in f if line.split()]
    out = oracle.read_smiles(texts)     # V2 cannot read a few, and those have nothing to compare
    assert len(out) > 4900, 'the corpus shrank: %d records' % len(out)
    return out


@fixture(scope='module')
def arenes(oracle):
    """`test/arenes.sdf` plus `test/heterocycles_charges.smi` -- aromatic and charged ring systems."""
    sdf = REPO / 'arenes.sdf'
    smi = REPO / 'heterocycles_charges.smi'
    if not sdf.is_file() or not smi.is_file():
        skip('the repo test data is not present')
    out = [(record.canonical, record)
           for record in oracle.read_sdf({'arenes': sdf})['arenes']]
    with smi.open(encoding='utf-8') as f:
        texts = [line.split()[0] for line in f if line.split()]
    return out + oracle.read_smiles(texts)


# ------------------------------------------------------------------------------------------------
# ORACLE 1.  RDKit, and the five records where it disagrees -- all five for the same reason.
@mark.parametrize('corpus', ['nci', 'arenes'])
def test_rdkits_verdict_on_our_output_is_its_verdict_on_the_input(corpus, request):
    """A PARSEABILITY claim, and the one that catches a token-level defect: a dropped ring closure, an
    unbalanced bracket, a bond token in the wrong place.

    Stated as an EQUALITY of verdicts rather than "RDKit reads everything we write", because RDKit
    rejects 25 of these records on input -- V2 accepts valences RDKit does not, `[N+](=O)(=O)` and an
    aromatic pyrrole-2,5-dione among them -- and demanding a readable output for an unreadable input
    would be demanding that the writer repair its input.  The measurement is that there is not ONE
    record in either direction: none we broke, and none we accidentally fixed.
    """
    canonical = rdkit_canonical()
    broke, fixed, neither = [], [], 0
    for text, source in request.getfixturevalue(corpus):
        out = write_smiles(bridge(source))
        before, after = canonical(text), canonical(out)
        if before is None and after is None:
            neither += 1
        elif after is None:
            broke.append((text, out))
        elif before is None:
            fixed.append((text, out))
    assert broke == []
    assert fixed == []
    assert neither < 30


@mark.parametrize('corpus', ['nci', 'arenes'])
def test_rdkit_agrees_on_the_constitution_except_where_the_hydrogen_count_is_UNSTATED(corpus,
                                                                                      request):
    """The gate: 4986 of 4991 comparable NCI records identical, 78 of 78 arene records, and all five
    differences are one thing.

    THE FIVE ARE NOT WRITER DEFECTS, and the assertion says so structurally rather than by listing
    them: every mismatching record contains an atom whose implicit hydrogen count V2 could not derive,
    and the writer REPORTS every one of those in `unknown_h`.  SMILES cannot spell "unstated" -- inside
    brackets an absent H term means zero, and a bare symbol means "derive it" -- so on such an atom the
    string necessarily states something, and the reader's derivation is its own.  `NC[S](=O)=O` is the
    whole story in one record: the input brackets say zero, we write a bare `S` because nothing else
    forces a bracket, and RDKit derives one hydrogen for it.  chython 2 loses the same hydrogen
    silently; the difference is the report.

    The count is asserted as an upper bound with the structural condition on top, so a NEW mismatch of
    any other kind fails even if the total stays five.
    """
    canonical = rdkit_canonical()
    mismatched = []
    comparable = 0
    for text, source in request.getfixturevalue(corpus):
        before = canonical(text)
        if before is None:
            continue                          # RDKit cannot read the record; not our disagreement
        comparable += 1
        mol = bridge(source)
        out = write_smiles(mol)
        if canonical(out) != before:
            mismatched.append((text, out, smw_traversal(mol)['unknown_h']))
    assert comparable > (4900 if corpus == 'nci' else 70)
    assert len(mismatched) <= (5 if corpus == 'nci' else 0)
    for text, out, unknown in mismatched:
        assert unknown != (), 'a mismatch with every hydrogen count stated: %s -> %s' % (text, out)


# ------------------------------------------------------------------------------------------------
# ORACLE 2.  chython 2, and how good an oracle it turns out to be.
def test_chython_2_reads_every_string_we_produce(nci, oracle):
    """The second parser, and a different set of rules to violate: V2 rejects valences RDKit allows."""
    outs = [write_smiles(bridge(source)) for _, source in nci]
    # one call for 5000 strings: `parse_canonical` answers None where V2 refuses, which is the whole
    # question here, and unlike `parse` it leaves no handle behind for a string we only wanted read
    answers = oracle.call('parse_canonical', texts=outs)
    failed = [(text, out) for (text, _), out, answer in zip(nci, outs, answers) if answer is None]
    assert failed == []


def test_chython_2s_canonical_string_is_not_automorphism_invariant(nci, oracle):
    """MEASURED, and it is why V2 is the SECOND oracle: 11 records, 4 without even a fixpoint.

    Write the record, hand the string back to V2, ask V2 for its canonical form of both molecules --
    and on eleven records the two differ.  Phloroglucinol (`OC1=CC(=CC(=C1)O)O`) is the clearest: the
    two strings are the same Kekule form of the same molecule related by a MIRROR automorphism, both
    are fixpoints of V2's own write-parse-write loop, and V2 has two attractors for one molecule.  Four
    more records do not even reach a fixpoint: `format(v2(format(m)))` differs from `format(m)`.

    So this is not a bug hunt, it is a bound on the oracle -- and the reason the RDKit test above is
    the gate rather than this one.  Every one of the eleven is checked ISOMORPHIC through V2's own
    matcher, which is the part of V2 that does not depend on a string.
    """
    pairs = [(source.handle, write_smiles(bridge(source))) for _, source in nci]
    differing = []
    for (text, _), row in zip(nci, oracle.call('roundtrip', pairs=pairs)):
        assert row['error'] is None, 'V2 could not read our string for %s: %s' % (text, row['error'])
        if row['here'] == row['there']:
            continue
        # "does not reach a fixpoint" is the negation of BOTH being fixpoints, not of either: the
        # allene measurement below needs the stronger statement, so the oracle reports the two flags
        # separately and each caller combines them the way its own claim requires
        oscillates = not (row['here_fixpoint'] and row['there_fixpoint'])
        differing.append((text, row['here'], row['there'], oscillates, row['same']))
    assert len(differing) <= 11
    for text, here, there, oscillates, same in differing:
        assert same, 'V2 says these are different molecules: %s vs %s' % (here, there)
    assert sum(1 for row in differing if row[3]) <= 4


# ------------------------------------------------------------------------------------------------
# ORACLE 3.  Ourselves, swept.  The claim no external tool can check.
@mark.parametrize('corpus', ['nci', 'arenes'])
def test_the_string_does_not_depend_on_the_creation_order(corpus, request):
    """Same molecule, six creation orders, one byte-identical string -- on every record of the corpus.

    THE EPIC'S CENTRAL CLAIM, and the fixtures cannot make it: a hand-built molecule has 5 to 12 atoms
    and a shape somebody chose, while these have up to 100 and shapes nobody chose.  A tie in the
    canonical ranking that the fixtures never reach is exactly what breaks here and only here.

    The seed is fixed so a failure is reproducible, and the orders include the identity so a corpus
    where every molecule is symmetric could not pass by accident.
    """
    rnd = Random(20260903)
    broken = []
    for text, source in request.getfixturevalue(corpus):
        n = len(source.atom_rows)
        strings = set()
        for k in range(6):
            order = list(range(n))
            if k:
                rnd.shuffle(order)
            strings.add(write_smiles(bridge(source, order)))
        if len(strings) != 1:
            broken.append((text, sorted(strings)))
    assert broken == []


def test_the_sweep_can_fail(nci):
    """Ruling F102.  Stored-slot order (`i`) is the writer with canonicalisation switched off.

    Asserted as "most records disagree" rather than all: a molecule with no branches and no symmetry
    can write the same string from two orders by luck, and a corpus this size will contain some.  The
    number is a floor on the sweep's resolution -- if it ever drops, the sweep above has stopped being
    a measurement.
    """
    rnd = Random(20260903)
    varying = 0
    sample = nci[:300]
    for text, source in sample:
        n = len(source.atom_rows)
        strings = set()
        for k in range(4):
            order = list(range(n))
            if k:
                rnd.shuffle(order)
            strings.add(write_smiles(bridge(source, order), 'i'))
        if len(strings) > 1:
            varying += 1
    assert varying > 0.9 * len(sample)


# ------------------------------------------------------------------------------------------------
# THE STEREO HALF.  Everything above this line runs on a corpus with no configuration in it.
#
# THE BRIDGE IS THE HARD PART, and it is worth saying why it is not a copy of the constitution
# bridge with one more field.  A stored parity byte is meaningless on its own: it is the parity of
# the anchor's refs IN THE ORDER THE ARENA HAPPENS TO HOLD THEM, so the creation order decides it
# (ruling F26).  What can be transferred between two libraries is a SIGN IN A STATED FRAME.  So the
# transfer runs in the direction that has an answer in both: V3's unit table names the frame, V2 is
# asked for its sign in THAT frame through `_translate_*_sign`, and the parity byte is then chosen --
# by trying both -- so that `translate_stereo` in that same frame reproduces it.  Nothing reads a
# byte from one library and writes it into the other.
#
# THE POLARITY IS MEASURED, NOT ASSUMED.  Which V2 sign means which V3 parity is a convention, and
# `flip=True` is the opposite convention: `test_the_bridges_polarity_is_MEASURED` runs the whole
# agreement test with it and the agreement collapses (290 -> 159 tetrahedral, 340 -> 1 cis/trans).
# Ruling F102 -- a fixture that cannot fail is not a measurement.
#
# WHAT THE BRIDGE CANNOT CARRY, all reported rather than swallowed:
#  * an atropisomer unit (9 of them, in 6 NIBR records) -- V2 has no model for one;
#  * a V2 sign no V3 unit anchors -- would silently lose a configuration V2 states;
#  * a frame where no parity reproduces the sign.
# A unit V3 anchors and V2 leaves UNSET is not a failure: V2's own string says nothing there either,
# so the comparison is still sound, and on the PubChem corpus this is the common case (3993 unset
# tetrahedral units against 49 configured cis/trans ones).
V2_KINDS = {0: 'tetrahedral', 1: 'cis/trans', 2: 'allene', 3: 'atropisomer'}


def bridge_stereo(source, order=None, flip=False):
    """`(V3 molecule, carried, unset, failures)` -- the constitution plus every sign V2 will state.

    `failures` is a list of strings; a non-empty one means this record cannot be compared, and the
    tests below say so rather than comparing anyway.  `flip` inverts the sign-to-parity convention
    and exists only so that the agreement can be shown to depend on it.
    """
    mol, ids = _bridge_atoms(source, order)
    back = {v: k for k, v in ids.items()}
    carried = unset = 0
    failures = []
    covered = set()

    # PHASE ONE, entirely on the core's side: enumerate the units and name a frame for each.  The
    # frames are what the oracle has to be asked about, and they are all known once the core has
    # enumerated its units -- so the whole record costs ONE round trip rather than one per unit.
    asking = []
    for unit in mol.stereo_units():
        kind = unit['kind']
        anchor = back[unit['anchor']]
        refs = unit['refs']
        covered.add(anchor)
        if kind == 3:
            failures.append('atropisomer at %d: chython 2 has no model for one' % anchor)
            continue
        if kind and (refs[0] is None or refs[2] is None):
            failures.append('%s at %d: a pair with no named slot 0' % (V2_KINDS[kind], anchor))
            continue
        if kind == 0:
            query = [0, anchor, [back[x] for x in refs if x is not None]]
        elif kind == 1:
            # a double bond V2 does not hold a counterpart for is a bond V2's cis/trans model does
            # not know about, so there is no frame to ask about and the unit counts as UNSET.  V2's own
            # lookup raises a KeyError there; it stays an unset here rather than becoming a crash,
            # because a crash would make the corpus unusable instead of reporting one unit V2 is
            # silent on.
            query = ([1, anchor, source.counterpart[anchor], back[refs[0]], back[refs[2]]]
                     if anchor in source.counterpart else None)
        else:
            query = [2, anchor, back[refs[0]], back[refs[2]]]
        asking.append((unit, kind, anchor, refs, query))

    # PHASE TWO: chython 2's sign in each of those frames.  `None` is V2's KeyError -- V2 states
    # nothing there, and neither will we.  Only the frames that exist are sent; the rest are unset
    # without a round trip.
    answers = iter(source.translate([q for *_, q in asking if q is not None]))
    for unit, kind, anchor, refs, query in asking:
        sign = None if query is None else next(answers)
        if sign is None:
            unset += 1
            continue
        want = (1 if sign else 2) if flip else (2 if sign else 1)
        for parity in (1, 2):
            mol.set_parity(unit['anchor'], parity)
            if mol.translate_stereo(unit['anchor'], refs) == want:
                carried += 1
                break
        else:
            mol.set_parity(unit['anchor'], 0)
            failures.append('%s at %d: no parity reproduces sign %r in the named frame'
                            % (V2_KINDS[kind], anchor, sign))
    for atom in source.atom_rows:
        if atom['stereo'] is not None and atom['n'] not in covered:
            failures.append('a chython 2 atom sign at %d that no V3 unit anchors' % atom['n'])
    for bond in source.bond_rows:
        if bond['stereo'] is not None and bond['n'] not in covered and bond['m'] not in covered:
            failures.append('a chython 2 bond sign at %d-%d that no V3 unit anchors'
                            % (bond['n'], bond['m']))
    return mol, carried, unset, failures


@fixture(scope='module')
def stereo_sdf(oracle):
    """`test/stereo.sdf` -- 300 records, 964 tetrahedral signs and 11 allene signs, no cis/trans.

    V2 perceives the configuration from the 2D coordinates and the wedges, so the signs here did not
    come from a SMILES string and cannot have been shaped by anybody's SMILES conventions.
    """
    path = REPO / 'stereo.sdf'
    if not path.is_file():
        skip('the repo test data is not present')
    out = oracle.read_sdf({'stereo': path})['stereo']
    assert len(out) == 300, 'the corpus changed: %d records' % len(out)
    return out


@fixture(scope='module')
def cis_trans(oracle):
    """The NIBR table's PubChem examples that V2 reads AND gives at least one bond sign: 346 records.

    Filtered on `/` or `\\` in the text before parsing, because parsing 2000 records to find 346 is
    the slow way round.
    """
    if NIBR is None or not NIBR.is_file():
        skip('the RDKit NIBR example table is not installed (looked for %s)' % NIBR)
    texts = []
    with NIBR.open(encoding='utf-8') as f:
        for row in DictReader(f):
            for key in ('EX1', 'EX2', 'EX3', 'EX4', 'EX5'):
                text = (row.get(key) or '').strip()
                if text and ('/' in text or '\\' in text):
                    texts.append(text)
    out = [(text, record) for text, record in oracle.read_smiles(texts)
           if any(bond['stereo'] is not None for bond in record.bond_rows)]
    assert len(out) > 300, 'the corpus shrank: %d records' % len(out)
    return out


def _stereo_agreement(records, flip=False):
    """`(agreed, disagreed, skipped)` against RDKit's canonical form of chython 2's OWN string.

    The reference is V2's output and not the input text, which matters here and did not above: the
    bridge can only carry what V2 perceived, so an input stating a configuration V2 drops would make
    the writer answer for V2's perception.  V2's own string is exactly what V2's sign store means.

    A record where the bridge carried NO sign is not counted at all.  It would agree under any
    convention, so counting it inflates both this measurement and the polarity measurement that
    depends on it.
    """
    canonical = rdkit_canonical()
    agreed = []
    disagreed = []
    skipped = 0
    for record in records:
        source = record[1] if isinstance(record, tuple) else record
        reference = canonical(source.canonical)
        if reference is None:
            continue
        mol, carried, _, failures = bridge_stereo(source, flip=flip)
        if failures:
            skipped += 1
            continue
        if not carried:
            continue
        out = write_smiles(mol)
        (agreed if canonical(out) == reference else disagreed).append((source.canonical, out))
    return agreed, disagreed, skipped


# ------------------------------------------------------------------------------------------------
# ORACLE 1 AGAIN, now with the configuration in it.
def test_rdkit_agrees_on_the_CONFIGURATION_and_not_only_the_constitution(stereo_sdf):
    """212 of 212 records that carry a configuration -- 975 signs, 964 tetrahedral and 11 allene.

    This is the fixture the epic's warning is about, and the reason it is a differential rather than
    a handful of hand-built centres: a wrong sign convention, a frame read in the wrong direction, or
    a ring-closure neighbour counted in the wrong place all produce a string that PARSES, and only a
    second implementation of the configuration model notices.  RDKit re-derives the configuration
    from our tokens and its canonical form is a fixpoint, so equality here is agreement about the
    molecule and not about the spelling.

    Six records are skipped and all six for one reason: they are biaryls where V3 finds an
    atropisomer unit and V2 has no atropisomer model, so there is no sign to ask V2 for.  The other
    82 of the 300 carry no configuration V2 will state -- 1858 unset units -- and are not counted,
    because a record with nothing to get wrong is not evidence.
    """
    agreed, disagreed, skipped = _stereo_agreement(stereo_sdf)
    assert disagreed == []
    assert len(agreed) >= 212
    assert skipped <= 6


def test_rdkit_agrees_on_every_cis_trans_double_bond(cis_trans):
    """340 of 346, and the six skipped are the atropisomers.  The `/` and `\\` half of the writer.

    `test/` has no double-bond geometry in any file, so without this corpus the direction tokens
    would be tested only by hand-built fixtures -- and direction tokens are the one part of SMILES
    stereo where the token is not attached to the atom it describes: a `/` belongs to a bond, its
    meaning depends on which end was written first, and a chain of them has to stay consistent
    across ring closures and branches.  That is exactly the kind of thing a corpus finds.
    """
    agreed, disagreed, skipped = _stereo_agreement(cis_trans)
    assert disagreed == []
    assert len(agreed) >= 340
    assert skipped <= 6


@mark.parametrize('corpus,agree_flipped', [('stereo_sdf', 80), ('cis_trans', 1)])
def test_the_bridges_polarity_is_MEASURED(corpus, agree_flipped, request):
    """Ruling F102 for the two tests above: with the opposite convention the agreement collapses.

    212 -> 80 on the SDF and 340 -> 1 on the double bonds.  Neither number is asserted to be zero,
    because inverting every sign of a molecule whose configuration is its own mirror image changes
    nothing, and 80 of the SDF records are like that.  That residue is also why `_stereo_agreement`
    counts no record whose bridge carried nothing: with those in, the flipped run agrees on 162 and the
    honest one on 159, and a measurement that goes UP when the convention is wrong is measuring the
    size of the corpus.
    """
    agreed, disagreed, _ = _stereo_agreement(request.getfixturevalue(corpus), flip=True)
    assert len(disagreed) > 0, 'the polarity does not matter, so this bridge proves nothing'
    assert len(agreed) <= agree_flipped


# ------------------------------------------------------------------------------------------------
# ORACLE 3 AGAIN.  The sweep: the epic's central claim, for double bonds and for tetrahedral centres,
# and no external tool can check either one.
def test_the_cis_trans_string_does_not_depend_on_the_creation_order(cis_trans):
    """340 records, four creation orders each, one string.  Including the design note's own witness.

    `smw_stereo_seed` exists because (2Z,4E)-hexa-2,4-diene wrote two strings over its 720 creation
    orders -- the constitution is symmetric end to end and the configuration is not, so the
    constitutional ranking tied and the slot order broke the tie.  Seeding the refinement with the
    frame-free parity code splits that tie, and this is the corpus-scale evidence that it does.
    """
    rnd = Random(20260903)
    broken = []
    for text, source in cis_trans:
        mol, _, _, failures = bridge_stereo(source)
        if failures:
            continue
        n = len(source.atom_rows)
        strings = set()
        for k in range(4):
            order = list(range(n))
            if k:
                rnd.shuffle(order)
            candidate, _, _, bad = bridge_stereo(source, order)
            strings.add('BRIDGE FAILED' if bad else write_smiles(candidate))
        if len(strings) != 1:
            broken.append((text, sorted(strings)))
    assert broken == []


def test_the_diene_the_seed_was_written_for_is_invariant_over_two_hundred_orders(oracle):
    """The single molecule the seed was written for, swept over 200 creation orders."""
    (_, source), = oracle.read_smiles(['C/C=C/C=C\\C'])
    rnd = Random(20260903)
    strings = set()
    for k in range(200):
        order = list(range(6))
        if k:
            rnd.shuffle(order)
        mol, _, _, failures = bridge_stereo(source, order)
        assert failures == []
        strings.add(write_smiles(mol))
    assert len(strings) == 1


def test_the_tetrahedral_string_does_not_depend_on_the_creation_order(stereo_sdf):
    """THE EPIC'S CENTRAL CLAIM, FOR TETRAHEDRAL CENTRES TOO: 0 of 294 records oscillate.

    The labelling is canonical up to the PARITY-REFINED automorphism group, and it has to be.  An
    element of the CONSTITUTIONAL group which INVERTS parity -- a mirror symmetry of an achiral
    molecule -- carries one extremal labelling to another whose string differs in every stereo token,
    and no colouring can split that, because the colouring is what the automorphism preserves.  Which
    is why the seed (`smw_stereo_seed`, splitting the ties a colouring CAN split) carries the cis/trans
    sweep above and not this one.

    Two mechanisms hold it, both inside the canonical search rather than in the writer.  The leaf
    certificate carries a PARITY TAIL, so two labellings that agree on constitution and disagree on
    configuration are not tied.  And the orbit prune does not take the constitutional orbits on faith:
    it feeds `mol_automorphisms` the parity-refined colouring `cls * 4 + digit`, so an automorphism only
    collapses two candidates when it preserves configuration as well as constitution, and where the
    colouring cannot NAME a configured unit's frame it declines to prune at all.  The digits come from
    `_canon_stereo_digits`, which calls the same
    `_frame_free_parity_seed` that `mol_identity_bytes` uses, so the search's reading of a parity and
    the identity's reading of it cannot drift apart.
    """
    rnd = Random(20260903)
    broken = []
    for source in stereo_sdf:
        mol, _, _, failures = bridge_stereo(source)
        if failures:
            continue
        n = len(source.atom_rows)
        strings = set()
        for k in range(4):
            order = list(range(n))
            if k:
                rnd.shuffle(order)
            candidate, _, _, bad = bridge_stereo(source, order)
            strings.add('BRIDGE FAILED' if bad else write_smiles(candidate))
        if len(strings) != 1:
            broken.append((source.canonical, sorted(strings)))
    assert broken == []


def test_nothing_oscillates_and_a_REGRESSION_would_have_to_be_a_MIRROR_PAIR(stereo_sdf):
    """Zero records oscillate, and the diagnostic that bounds a regression is kept as a tripwire.

    The final assertion is the one with teeth: at a count of zero the loop body does not run.  It stays
    because if the parity tail or the parity-refined orbit prune regresses, this is the test that says
    which kind of regression it is -- the mirror case un-fixed (the variants are still one molecule by
    RDKit's reckoning, so the assertion inside the loop passes while the count assertion fails) or a
    configuration corrupted outright (the inner assertion fails first and its message names the
    strings).  Those are very different bugs and the message should not have to be guessed at.

    chython 2, measured the same way on the same corpus, oscillates on 48 of 300, which is why V2 is
    not the oracle for this behaviour.
    """
    canonical = rdkit_canonical()
    rnd = Random(20260903)
    oscillating = 0
    for source in stereo_sdf:
        mol, _, _, failures = bridge_stereo(source)
        if failures:
            continue
        n = len(source.atom_rows)
        strings = set()
        for k in range(4):
            order = list(range(n))
            if k:
                rnd.shuffle(order)
            candidate, _, _, bad = bridge_stereo(source, order)
            if not bad:
                strings.add(write_smiles(candidate))
        if len(strings) == 1:
            continue
        oscillating += 1
        molecules = set()
        for text in strings:
            molecules.add(canonical(text))
        assert len(molecules) == 1, 'the variants are DIFFERENT MOLECULES: %s' % sorted(strings)
    assert oscillating == 0


def test_the_minimal_witness_is_cis_1_4_dimethylcyclohexane_and_chython_2_still_fails_it(oracle):
    """Eight atoms.  V3 writes ONE string across 200 creation orders; chython 2 writes two.

    cis-1,4-dimethylcyclohexane is the smallest witness in the corpus and it is ACHIRAL -- measured,
    not assumed: RDKit embeds `C[C@H]1CC[C@@H](C)CC1` with both methyls on the same face of the ring,
    and its canonical form is invariant under inverting every sign in the string.  The graph
    automorphism that swaps the two CH2 branches fixes both stereocentres, exchanges two of each one's
    neighbours, and therefore inverts both parities: it maps the molecule to its mirror image, which
    is itself.  So `[C@H]...[C@@H]` and `[C@@H]...[C@H]` are two spellings of one compound and there
    is no CONSTITUTIONAL reason to prefer either.  Which is why a search pruning on constitutional
    orbits alone takes whichever candidate its own slot order reaches first, and why a seed cannot fix
    it: the seed labels are IDENTICAL for the two creation orders that produce the two strings --
    `{methyl: 12, CH: 5, CH2: 8}` in both -- so the tie is not one the colouring failed to split, it is
    one the colouring cannot see.  The parity-refined orbit prune CAN see it: with a digit on each anchor
    the swap does not preserve the colouring, so it does not collapse the two candidates, and the parity
    tail on the certificate decides between them on configuration.

    The second half records what chython 2 does, measured.  2.24, remapped the same 200 ways, writes two
    strings for this molecule and they are each other's mirror spelling: swapping `@` for `@@` and back
    carries one to the other byte for byte, so the skeleton, the branches and the ring closure are
    identical and the whole difference is which enantiomeric labelling its search reaches.  It is
    asserted so that V2's behaviour is a measured fact and not a remembered one; nothing in V3 is shaped
    to match it, and if a future 2.x answers differently this assertion is the thing to delete.
    """
    (_, source), = oracle.read_smiles(['C[C@H]1CC[C@@H](C)CC1'])
    rnd = Random(20260903)
    strings = set()
    for k in range(200):
        order = list(range(8))
        if k:
            rnd.shuffle(order)
        mol, carried, _, failures = bridge_stereo(source, order)
        assert failures == [] and carried == 2
        strings.add(write_smiles(mol))
    assert len(strings) == 1, sorted(strings)
    canonical = rdkit_canonical()
    assert canonical(strings.pop()) == canonical('C[C@H]1CC[C@@H](C)CC1')

    # and chython 2 measured THE SAME WAY we measure ourselves -- renumbered, not re-parsed, so the
    # comparison is between two canonical labellings of one molecule and not between two parsers
    orders = []
    for k in range(200):
        shuffled = list(source.numbers)
        if k:
            rnd.shuffle(shuffled)
        orders.append(shuffled)
    v2_strings = set(oracle.call('remap_canonical', handle=source.handle, orders=orders))
    assert len(v2_strings) == 2, 'chython 2 no longer oscillates here: %s' % sorted(v2_strings)
    first, second = sorted(v2_strings)
    assert first.replace('@@', '\0').replace('@', '@@').replace('\0', '@') == second


# ------------------------------------------------------------------------------------------------
# ALLENES.  Eleven signs, no external oracle, and the one place the writer REFUSES.
def test_the_allene_configuration_is_written_or_REPORTED_and_never_dropped(stereo_sdf):
    """Ten records carry an allene sign: eight get a token, two are refused and named in `lost`.

    The two refusals are cumulated -- `C=C=C=C=C`, an axis of four double bonds -- and SMILES has no
    agreed notation for an axis longer than an allene's, so the writer declines to invent one.  What
    makes that acceptable is the second half of the assertion: the declined unit arrives in
    `smw_traversal(...)['lost']`, so a caller who needs the configuration can see that the string
    does not carry it.  A silent drop here would be indistinguishable from a molecule that never had
    a configuration, and that is the failure mode this test exists to prevent.
    """
    spelled = 0
    refused = 0
    for source in stereo_sdf:
        centres = [atom['n'] for atom in source.atom_rows
                   if atom['stereo'] is not None and atom['n'] in source.allenes]
        if not centres:
            continue
        mol, carried, _, failures = bridge_stereo(source)
        assert failures == [], failures
        out = write_smiles(mol)
        lost = smw_traversal(mol)['lost']
        axis = [unit for unit in mol.stereo_units() if unit['kind'] == 2]
        assert len(axis) == len(centres)
        if lost:
            refused += 1
            assert len(lost) == len(centres)
            assert '@' not in out            # the two refused records have no other centre either
        else:
            spelled += 1
            assert '@' in out
    assert spelled == 8
    assert refused == 2


def test_chython_2_has_TWO_ATTRACTORS_for_an_allene_so_it_cannot_arbitrate(stereo_sdf, oracle):
    """Why the allene sign is not asserted against an oracle: there is no oracle here to assert it
    against, and this test measures the absence rather than leaving it unsaid.

    RDKit's SMILES parser drops extended tetrahedral stereo -- `ClC=C=CCl` is its canonical form of
    both handednesses of 1,3-dichloroallene, and its InChI of either loses the axis too -- so oracle
    1 is blind here.  chython 2 cannot arbitrate either: our string for that record, re-read by V2, is
    a FIXPOINT of V2's own write-parse-write loop, and it is not V2's fixpoint.  Two attractors for one
    molecule, exactly as with phloroglucinol above.

    And the two strings are the same molecule, which can be settled by hand because it is four
    substituents: ours is `ClC([H])=[C@]=C([H])Cl`, substituents in written order (Cl, H | H, Cl);
    V2's is `[H]C(Cl)=[C@@]=C([H])Cl`, (H, Cl | H, Cl).  One transposition apart, one symbol apart,
    same configuration.  So the disagreement is a spelling, and the assertion below is the
    structural one that is actually available: every record we differ on is a record where V2 has
    settled into a second fixpoint of its own.
    """
    pairs = []
    for source in stereo_sdf:
        if not any(atom['stereo'] is not None and atom['n'] in source.allenes
                   for atom in source.atom_rows):
            continue
        mol, _, _, failures = bridge_stereo(source)
        assert failures == []
        pairs.append((source.handle, write_smiles(mol)))

    differing = 0
    for row in oracle.call('roundtrip', pairs=pairs):
        assert row['error'] is None, row['error']
        if row['here'] == row['there']:
            continue
        differing += 1
        # BOTH fixpoints, separately -- the claim is that V2 has settled into a SECOND attractor of
        # its own loop, which is strictly stronger than "at least one of the two strings oscillates"
        assert row['there_fixpoint'], 'not a fixpoint, so V2 could still arbitrate'
        assert row['here_fixpoint']
    assert differing == 4       # two of them the refusals above, two the two-attractor pair


def test_the_bridge_reports_the_atropisomers_it_cannot_carry(cis_trans):
    """Six records, nine units, and a skip that is written down instead of being quietly absent.

    chython 2 has no atropisomer model at all, so there is no sign to ask it for -- which means these
    records cannot be compared, NOT that they are correct.  The test asserts the count so that a
    future bridge that starts carrying them changes a number here rather than silently widening what
    the corpus is taken to prove.
    """
    reported = 0
    records = 0
    for text, source in cis_trans:
        _, _, _, failures = bridge_stereo(source)
        atropisomers = [f for f in failures if f.startswith('atropisomer')]
        if atropisomers:
            records += 1
            reported += len(atropisomers)
    assert (records, reported) == (6, 9)
