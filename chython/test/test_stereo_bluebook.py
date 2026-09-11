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
"""The stereo epic's acceptance gate: the core against an independent authority, on real molecules.

`test/stereo.sdf` holds 300 IUPAC Blue Book stereochemistry examples with the units each one carries
annotated in `STEREOGENIC_UNITS`.  That annotation is the only oracle here: there is no differential
test against another implementation anywhere in this file, and the corpus is read by `chython.formats`.

WHY THIS FILE SITS IN `chython/test/` AND NOT IN `chython/core/test/`, WHICH IS WHERE ITS SUBJECT
LIVES.  It needs two things from two different layers: the core's stereo perception, and a reader for
an MDL corpus.  The dependency direction is `core <- chemistry <- {formats, ...}`, and
`chython/core/test/test_no_chython_two_imports.py` holds the core's own suite to it -- the core must be
shippable and testable alone, so nothing under `chython/core/` may import `chython.formats`.  A gate
that genuinely needs both layers therefore belongs above both, which is this directory.

WHAT CROSSES THE CONSTITUTION BRIDGE, AND WHY THAT LIST IS EXACTLY THIS LONG.  Element,
isotope, charge, radical, implicit hydrogen count and bond order.  **No stereo parity crosses it,
and neither do the 2d coordinates.**  The reader is asked for `ignore_stereo=True` and the rebuild
below copies no parity, so what this gate sees is a constitution the core must find its own units in
-- which is the thing under test.  Reading the file's wedges and then asserting the units they imply
would test the MDL stereo layer instead, and that layer has its own suite.

Per ruling F26 a stored parity byte is FRAME-RELATIVE to atom creation order, so a `_relabel` that
copied `parity_of()` verbatim would produce a *different molecule* and gate 1 would then assert
invariance over subjects that are not the same subject.  Nothing here reads a parity at all --
`test_no_parity_crosses_the_bridge` asserts that positively -- so `translate_stereo` is not needed and
ruling F26 is satisfied by construction.  If this bridge ever grows a stereo write, the parity must move
through `translate_stereo` into the new frame, never by assignment.

TWO DECISIONS TAKEN BEFORE THE CORPUS LOOP WAS WRITTEN, stated here rather than left implicit in a
bare `except`:

1. **Ruling F100's `RuntimeError: two stereo units claim one anchor atom` is EXCLUDED from this
   gate, not xfailed.**  It fires on 10 of 3,000 randomised multi-component records, identically at
   BASE and at HEAD, and is pre-existing (filed against Task 3 / ruling F45).  This gate excludes
   the class by choosing a corpus that does not contain it, and the exclusion is a *measured
   precondition* rather than a hope: `test_no_two_units_contest_one_anchor` asserts over all 300
   records that no two units share an anchor.  The alternative -- a randomised multi-component
   corpus with 10 xfails in it -- would put a known unrelated defect inside the acceptance gate and
   make the gate's own colour depend on it.  A corpus this gate cannot state a verdict on is not
   evidence, so it is not in the corpus.
2. **`stereo_truncated is False` is NEVER asserted corpus-wide without a size bound.**  The flag
   returns at roughly 5,760 atoms on disjoint-copy records through the whole-call node budget, with
   the marks still exactly correct -- a false alarm, not unsoundness.  So every place this file
   asserts the flag asserts an atom-count bound in the same breath
   (`test_no_record_truncates_below_the_node_budget`), and the bound is what makes the assertion
   true rather than luck.

RULING F102 -- A CORPUS IS EVIDENCE ONLY IF IT COULD HAVE FAILED.  300 records of clean output prove
nothing on their own, so this file carries its own positive control:
`test_gate1_reports_a_dirty_result_when_the_relabelling_drops_a_property` breaks `_relabel` on
purpose, one constitutional property at a time, and requires every record that carries that property
to be reported as oscillating.  It is what turns gate 1's zero into a measurement.

That control also fixes a real hole in gate 1 as originally specified.  Gate 1 compares 60
relabelings *to each other*, so a `_relabel` that drops a property UNIFORMLY produces 60 forms that
all agree and the gate passes -- measured: dropping the isotope from `_relabel` alone leaves gate 1
green on all 300 records.  `_forms()` therefore includes the ORIGINAL molecule's canonical string in
the compared set, and with that one line the same sabotage reports 16 oscillating records.

Gates 2 and 5 of the spec's §9 are Task 6's pseudo-asymmetry tests and Task 9's silent-loss tests.
This file covers gates 1, 3 and 4 and asserts nothing twice.
"""
from collections import Counter
from pathlib import Path
from random import Random

import pytest

from chython.core import MoleculeContainer
from chython.core.test.test_stereo_perception import _disjoint


# test/stereo.sdf is tracked in git, so a missing file is a broken checkout and must fail loudly.
# A skipif here would silently retire the epic's entire acceptance gate.
#
# The root is found by looking for `pyproject.toml` rather than by counting `parents[N]`.  A count is
# a second fact about where this file sits, and it is wrong the moment the file moves -- which is
# exactly what happened when this gate came up out of `chython/core/test/`.
def _repo_root():
    for candidate in Path(__file__).resolve().parents:
        if (candidate / 'pyproject.toml').is_file():
            return candidate
    raise RuntimeError('cannot locate the repository root: no pyproject.toml above this file')


SDF = _repo_root() / 'test' / 'stereo.sdf'

# `stereo_unit_t.kind`, as `stereo_units()` reports it.  4 (helical) is RESERVED and never produced.
SU_TETRA, SU_CIS_TRANS, SU_ALLENE, SU_ATROPISOMER, SU_HELICAL = 0, 1, 2, 3, 4

# The largest record in the file, measured.  The node budget that sets `stereo_truncated` is two
# orders of magnitude above this, which is what licenses the flag assertions in this file.
MAX_ATOMS = 72


def _constitution(src):
    """Rebuild a molecule carrying its constitution and nothing else -- no stereo, no coordinates.

    The implicit hydrogen count is part of that constitution: the core never derives one, so
    dropping it would leave every CH stereocentre in the file with three directions and no unit.

    THE TWO "UNSET" NORMALISATIONS, chosen so this function and `_relabel` are the same function of
    the same atom.  `add_atom(isotope=0)` means "unset" and `isotope_of` returns 0 for unset, so
    `a.isotope or 0` here and `isotope_of` there land on the same value.  `add_atom(implicit_h=None)`
    means "unstated", the arena stores `H_UNKNOWN` for it and `implicit_h_of` answers None -- so an
    unstated count (2 atoms in the whole file, both on the S#I bond of VS170, where no valence rule
    applies) survives the round trip as itself.

    THE SECOND ONE IS AN IDENTITY AND NOT A COINCIDENCE: None goes in, None comes out, and VS170's two
    undeterminable counts are still undeterminable at the far end.  Neither normalisation is free: see
    `test_gate1_reports_a_dirty_result_when_the_relabelling_drops_a_property`.
    """
    m = MoleculeContainer()
    with m.edit():
        fresh = {}
        for a in src.atoms():
            fresh[a.n] = m.add_atom(a.element, isotope=a.isotope or 0, charge=a.charge,
                                    radical=a.is_radical, implicit_h=a.implicit_h)
        for b in src.bonds():
            m.add_bond(fresh[b.n], fresh[b.m], b.order)
    return m


@pytest.fixture(scope='module')
def records():
    # `ignore_stereo=True` because this gate is about the core finding its own units in a
    # constitution: reading the file's wedges would hand it the answer.  Read off the reader here,
    # and `mol.meta` holds the same dict -- copied out because the fixture outlives the reader.
    from chython.formats import SDFRead
    out = {}
    with SDFRead(str(SDF), ignore_stereo=True) as f:
        for mol in f:
            out[f.meta['STRUCTURE_ID']] = (_constitution(mol), dict(f.meta))
    assert len(out) == 300, f'expected 300 records, read {len(out)}'
    return out


def _edges(m):
    return {(min(a, b), max(a, b), m.order_of(a, b))
            for a in m.atom_numbers for b in m.neighbors_of(a)}


def _relabel(m, rng, drop=None):
    """Another encoding of `m`: the same constitution built in a shuffled atom order.

    `drop` names a constitutional property to omit, and exists for the ruling F102 positive control
    ONLY.  Every gate calls this with `drop=None`; a gate that passed with a property dropped would
    be a gate whose subject that property does not reach.

    EVERY `drop` WRITES THE PROPERTY'S EMPTY VALUE, and for `implicit_h` that is `0` and not `None`.
    `None` was the empty value while it stored a zero; it now stores `H_UNKNOWN`, which is a DIFFERENT
    value rather than an absent one, and dropping to it moved seven records whose counts were all
    genuinely zero -- they carry nothing, so a sabotage must not be able to disturb them, and the test
    below asserts exactly that by requiring `dirty == carrying` and not `dirty >= carrying`.  An
    erasing sabotage tests sensitivity; a substituting one tests that two values differ, which is a
    weaker claim wearing the same assertion.
    """
    perm = list(m.atom_numbers)
    rng.shuffle(perm)
    out = MoleculeContainer()
    with out.edit():
        fresh = {sid: out.add_atom(m.element_of(sid),
                                   isotope=0 if drop == 'isotope' else m.isotope_of(sid),
                                   charge=0 if drop == 'charge' else m.charge_of(sid),
                                   radical=False if drop == 'radical' else m.radical_of(sid),
                                   implicit_h=0 if drop == 'implicit_h' else m.implicit_h_of(sid))
                 for sid in perm}
        for a, b, order in _edges(m):
            out.add_bond(fresh[a], fresh[b], order)
    return out


def _canonical_string(m):
    """`m`'s constitution written out in `canonical_order()` positions.

    Compared through `canonical_order` and per-atom properties, never through a canonical SMILES
    string: that output oscillates on symmetric stereocentres and is unsound as an identity check.
    The atom tuple carries every property `_relabel` copies, so a relabelling that CHANGED one
    cannot compare equal -- which is the whole point of the tuple and is what the F102 control
    measures.
    """
    pos = m.canonical_order()
    atoms = [None] * len(pos)
    for sid, p in pos.items():
        atoms[p] = (m.element_of(sid), m.isotope_of(sid), m.charge_of(sid),
                    m.radical_of(sid), m.implicit_h_of(sid))
    edges = sorted((min(pos[a], pos[b]), max(pos[a], pos[b]), order)
                   for a, b, order in _edges(m))
    return repr((atoms, edges))


def _forms(m, rng, count=60, drop=None):
    """The distinct canonical forms of `m` over `count` relabelings AND the original.

    THE ORIGINAL IS IN THE SET DELIBERATELY.  Comparing relabelings only to each other cannot see a
    `_relabel` that drops a property uniformly -- measured: with the isotope dropped, all 300 records
    still give exactly one form.  With the original included, the 16 isotope-carrying records give
    two.  One line, and it is the difference between a gate and a tautology.
    """
    return {_canonical_string(m)} | {_canonical_string(_relabel(m, rng, drop)) for _ in range(count)}


def _codes(meta):
    """The annotated unit codes. STEREOGENIC_UNITS holds 'TH', 'CT,TH', 'TH5,CT' and so on;
    BLUE_BOOK_REF holds the P-number and must NOT be used for this.

    Exact token equality after splitting on ',', never substring matching: 'TH' is a substring of
    'TH3' and 'TH5', so a substring test would silently inflate every TH assertion.

    The field is ABSENT on VS002 and VS004, which is what the default handles -- see
    `test_the_two_unannotated_records_have_nothing_to_annotate` for what that absence means.
    """
    return {p.strip() for p in meta.get('STEREOGENIC_UNITS', '').split(',') if p.strip()}


# code -> the unit kind it must produce. Odd-length axes are axial, even ones cis/trans-like.
CODE_KIND = {'TH': SU_TETRA, 'CT': SU_CIS_TRANS, 'CT4': SU_CIS_TRANS,
             'TH3': SU_ALLENE, 'TH5': SU_ALLENE, 'AT': SU_ATROPISOMER}

# the census the suite actually contains, confirmed by reading the file
CODE_COUNT = {'TH': 249, 'CT': 65, 'AT': 7, 'TH3': 8, 'CT4': 5, 'TH5': 2, 'HE': 2}


# --- The corpus is the file we think it is, and the shape both up-front decisions rest on --------

def test_the_corpus_is_the_file_we_think_it_is(records):
    """The gate's own preconditions, measured rather than assumed.

    Three of these are load-bearing for assertions elsewhere in the file: `MAX_ATOMS` licenses every
    `stereo_truncated is False`, the multi-component census names the one record whose units live in
    two components, and the radical count is an admission rather than a claim.

    THE RADICAL FIELD OF `_canonical_string` CARRIES NO COVERAGE HERE.  No record in this file has a
    radical, so `_relabel(drop='radical')` changes nothing and the F102 control below reports zero
    for it -- correctly.  That is a documented gap in what this corpus can exercise, not a property
    the corpus certifies.  A radical is covered by `test_features.py`, not by this file.
    """
    assert len(records) == 300
    assert all(sid.startswith('VS') for sid in records), 'STRUCTURE_ID is the VS-numbered oracle key'

    sizes = {sid: m.atom_count for sid, (m, _) in records.items()}
    assert max(sizes.values()) == MAX_ATOMS, 'the size bound the truncation assertions rest on'

    multi = {sid for sid, (m, _) in records.items() if m.connected_components_count > 1}
    assert multi == {'VS031', 'VS054', 'VS068', 'VS095', 'VS165', 'VS186'}, \
        'six multi-component records; five carry a one-atom counter-ion and only VS186 has units ' \
        'in two components, which is why ruling F103 needs a record this corpus does not contain'

    assert not [sid for sid, (m, _) in records.items()
                if any(m.radical_of(a) for a in m.atom_numbers)], \
        'no radical anywhere: the radical axis is uncovered by this corpus, and said so'
    assert len([sid for sid, (m, _) in records.items()
                if any(m.isotope_of(a) for a in m.atom_numbers)]) == 16
    assert len([sid for sid, (m, _) in records.items()
                if any(m.charge_of(a) for a in m.atom_numbers)]) == 17


def test_no_two_units_contest_one_anchor(records):
    """Ruling F100's exclusion, as a measured precondition instead of a bare `except`.

    Perception's unit SET is order-dependent where two bond-kind units contest one anchor atom, and
    the collision raises `RuntimeError` outright.  It is pre-existing (filed against Task 3 / ruling
    F45), it fires on 10 of 3,000 randomised multi-component records, and this gate excludes the
    class by not containing it.  This test is what makes "does not contain it" a fact: if a future
    change gave any Blue Book record two units on one anchor, this fails by name here rather than
    surfacing as an opaque raise from an unrelated gate.
    """
    contested = {}
    for sid, (m, _) in records.items():
        counts = Counter(u['anchor'] for u in m.stereo_units())
        if any(v > 1 for v in counts.values()):
            contested[sid] = [a for a, v in counts.items() if v > 1]
    assert not contested, f'ruling F100 reaches this corpus after all: {contested}'


def test_no_record_truncates_below_the_node_budget(records):
    """`stereo_truncated is False` -- asserted WITH the size bound that makes it true.

    The flag is a false alarm at scale, not unsoundness: it returns at roughly 5,760 atoms on
    disjoint-copy records through the whole-call node budget while the marks stay exactly correct.
    So the assertion is only meaningful next to the bound, and the bound is asserted here rather
    than trusted.  Never write this one corpus-wide without it.
    """
    assert max(m.atom_count for m, _ in records.values()) == MAX_ATOMS
    truncated = [sid for sid, (m, _) in records.items() if m.stereo_truncated]
    assert not truncated, f'{truncated} truncated at {MAX_ATOMS} atoms or fewer'


def test_no_parity_crosses_the_bridge(records):
    """The bridge's central claim, asserted positively rather than inferred from a clean validator.

    `validate_stereo() == []` in `test_no_record_raises` says nothing was written UNJUSTIFIABLY; this
    says nothing was written at all, which is the stronger statement and the one gate 1's soundness
    depends on.  Per ruling F54 a `parity_of` of 0 means "no parity configured" and never "no wedge
    drawn" -- nothing here reads it as evidence about a drawing.

    The group views come back empty for the same reason, and per Task 11 `{}` from
    `canonical_stereo_groups()` means SEG_STEREO_GROUPS is absent, not "no stereo".  `()` from
    `canonical_stereo_group_ambiguities()` is the strongest of its answers and is what a molecule
    with no stored parity must give.  Neither is joined to `canonical_order()` anywhere.
    """
    for sid, (m, _) in records.items():
        assert all(m.parity_of(a) == 0 for a in m.atom_numbers), sid
        assert all(m.stereo_of(a) is False for a in m.atom_numbers), sid
        assert m.canonical_stereo_groups() == {}, sid
        assert m.canonical_stereo_group_ambiguities() == (), sid


# --- Gate 1: random-relabeling invariance, the spec's headline defect -------------------
#
# EVERY comparison in this gate is chython-against-chython, and that is a hard rule, not a
# convenience. The canonical order takes its extremum over an invariantly-selected subset of the
# target cell rather than the whole cell (spec 6), so it produces a *different* invariant
# representative than a whole-cell search would -- and than nauty or RDKit do. All of those are
# equally canonical; none of them agree numerically. So: never assert against nauty or RDKit
# canonical ranks, never against a hash captured from a build that did a whole-cell search, and
# never store a canonical string as a fixture expecting a later build to reproduce it. Compare
# forms produced by the same build, within one test run. What is guaranteed is that relabelings
# of one molecule agree with each other; what is NOT guaranteed is which representative they
# agree on.


def test_gate1_canonical_order_is_relabeling_invariant(records):
    """60 relabelings of each record, plus the original, must give exactly one canonical form."""
    rng = Random(20260901)
    broken = {}
    for sid, (m, _) in records.items():
        forms = _forms(m, rng)
        if len(forms) != 1:
            broken[sid] = len(forms)
    assert not broken, f'{len(broken)} records still oscillate: {broken}'


@pytest.mark.parametrize('structure_id', ['VS009', 'VS196', 'VS252', 'VS255'])
def test_gate1_the_four_worst_offenders_converge(records, structure_id):
    """The corpus's four worst oscillators, named, so a regression on them cannot hide in an aggregate."""
    rng = Random(1)
    assert len(_forms(records[structure_id][0], rng)) == 1


def _read(m, prop, sid):
    """The one property `_relabel(drop=prop)` omits, read back off the core."""
    if prop == 'isotope':
        return m.isotope_of(sid)
    elif prop == 'charge':
        return m.charge_of(sid)
    elif prop == 'radical':
        return m.radical_of(sid)
    return m.implicit_h_of(sid)


# (property, how many of the 300 records carry it, the named records the sabotage must reach).
# The named sets are the records this file makes OTHER claims about, so a sabotage that missed them
# would be a sabotage that misses the evidence: the three isotope records are gate 4's isotope-decided
# centres, the three charge records are gate 4's onium and N-oxide centres, and VS002/VS004 are the
# unannotated pair whose entire "nothing to annotate" verdict rests on a CH2's hydrogen count.
# VS009 -- gate 1's worst oscillator -- is deliberately NOT among them: it carries every hydrogen
# explicitly and so is not an implicit_h carrier at all.
SABOTAGE = [('isotope', 16, {'VS180', 'VS183', 'VS186'}),
            ('charge', 17, {'VS031', 'VS045', 'VS054'}),
            ('implicit_h', 293, {'VS002', 'VS004', 'VS010'}),
            ('radical', 0, set())]


@pytest.mark.parametrize('prop,carriers,named', SABOTAGE)
def test_gate1_reports_a_dirty_result_when_the_relabelling_drops_a_property(records, prop, carriers,
                                                                            named):
    """RULING F102: gate 1's zero above is evidence only from an instrument that can report dirt.

    `_relabel` is sabotaged one property at a time, and the requirement is EXACT rather than "some
    record fails": the set of records reported as oscillating must be precisely the set of records
    that CARRY the dropped property.  A weaker instrument would report a subset (insensitive on the
    rest) and a broken comparison would report a superset (manufacturing its own difference), so
    equality rules out both directions at once.  Measured: isotope 16, charge 17, implicit_h 293.

    `radical` is parametrised at ZERO on purpose.  No record in this file has one, so the property
    is uncarried, the sabotage is a no-op and the honest answer is an empty set -- an admission that
    `_canonical_string`'s radical field is untested here, kept visible instead of dropped from the
    parametrisation where nobody would notice it was missing.

    Four relabelings, not sixty: a property that reaches `canonical_order` at all separates the
    original from the very first shuffle, and this test is about sensitivity rather than about
    convergence.
    """
    carrying = {sid for sid, (m, _) in records.items()
                if any(_read(m, prop, a) for a in m.atom_numbers)}
    assert len(carrying) == carriers, f'{prop} is carried by {len(carrying)} records, not {carriers}'

    rng = Random(4242)
    dirty = {sid for sid, (m, _) in records.items() if len(_forms(m, rng, 4, drop=prop)) != 1}
    assert dirty == carrying, \
        f'dropping {prop} was reported on {len(dirty)} records but is carried by {len(carrying)}'
    assert named <= dirty, f'{sorted(named - dirty)} carry {prop} and were not reported'


# --- Gate 3: spiro and ring double bonds, the obligatory criterion ----------------------

@pytest.mark.parametrize('structure_id', ['VS063', 'VS078', 'VS120'])
def test_gate3_spiro_and_macrocyclic_cumulenes_are_found(records, structure_id):
    """The obligatory criterion: spiro and ring-double-bond handling is the one behaviour that must not
    regress.

    VS063 is annotated CT4,TH -- an even cumulene in a spiro system, so kind 1, NOT an allene.
    VS078 and VS120 are TH,TH3: macrocyclic allenes carrying a remote tetrahedral centre in
    the same ring, so kind 2.
    """
    m, meta = records[structure_id]
    kinds = {u['kind'] for u in m.stereo_units()}
    want = {CODE_KIND[c] for c in _codes(meta) if c in CODE_KIND}
    assert want, f'{structure_id} carries no code this gate knows: {_codes(meta)}'
    assert want <= kinds, f'expected kinds {sorted(want)}, got {sorted(kinds)}'


@pytest.mark.parametrize('code', [
    'TH3', 'TH5', 'CT4',
    # RULING F85: this fails today and is landed as a strict xfail rather than weakened.  VS023 is
    # annotated AT and the core finds no kind-3 unit on it: the record is a BRIDGED biaryl -- the two
    # aryl rings are joined by the pivot bond AND by a three-atom -CH2-C(=O)-CH2- bridge -- so the
    # axis is itself a ring bond, and `_is_atropisomer_axis` refuses a ring bond outright
    # (`e.flags & HE_IN_RING`).  The other six AT records all pass.  Chemically the annotation is
    # right: a bridge that short locks the axis absolutely, so this is the one class of atropisomer
    # where a ring axis is MORE hindered rather than less.  Owned by Task 3 (`_is_atropisomer_axis`,
    # spec 4.5's heuristic); strict, so the moment it is fixed this flips to a failure and the
    # marker is removed rather than rotting into a permanent exemption.
    pytest.param('AT', marks=pytest.mark.xfail(strict=True, reason=(
        'VS023: a bridged biaryl whose axis is a ring bond, refused by _is_atropisomer_axis\'s '
        'HE_IN_RING clause. Annotated AT by the Blue Book; 6 of 7 AT records pass. '
        'Ruling F85 -- a finding for Task 3, not an assertion to weaken.'))),
])
def test_gate3_annotated_codes_produce_their_unit_kind(records, code):
    """Each axial, even-cumulene and atropisomer annotation must yield the right kind.

    RULING F83: `missing` is empty both when every annotated record passes and when NO record
    carries the code, so the subject is counted here rather than a census test two screens away
    doing it at a distance.  `TH` is deliberately absent from the parametrisation -- VS138 is
    annotated `TH` on a phosphine the core refuses by scope decision, which
    `test_gate4_phosphine_lone_pair_stays_excluded` is about.
    """
    want = CODE_KIND[code]
    carrying = [sid for sid, (_, meta) in records.items() if code in _codes(meta)]
    assert len(carrying) == CODE_COUNT[code] > 0, \
        f'{code} is carried by {len(carrying)} records; the assertion below would be vacuous at 0'

    missing = [sid for sid in carrying
               if not any(u['kind'] == want for u in records[sid][0].stereo_units())]
    assert not missing, f'{code} -> kind {want} missing on {missing}'


@pytest.mark.parametrize('code,expected', sorted(CODE_COUNT.items()))
def test_gate3_the_annotation_census_is_what_we_think(records, code, expected):
    """A guard on the gate itself: if the SDF changes, the tests above silently weaken."""
    got = sum(1 for _, (_, meta) in records.items() if code in _codes(meta))
    assert got == expected


# --- Gate 4: the 26 centres of the spec's 10.2 census -----------------------------------

# (structure id, the anchor's atomic number, a label for the failure message)
# from the spec's 10.2 census: the centre kinds a stereo model can fail to represent at all
UNREPRESENTABLE = [
    ('VS071', 15, 'phosphine oxide'),
    ('VS147', 16, 'sulfoxide'),
    ('VS180', 6, 'isotopic methane'),
    ('VS183', 16, 'isotopic sulfone'),
    ('VS186', 16, 'sulfonimidoyl'),
    ('VS045', 7, 'N-oxide'),
    ('VS031', 7, 'ammonium'),
    ('VS054', 15, 'phosphonium'),
    ('VS104', 14, 'silicon'),
    ('VS130', 16, 'thio-sulfonyl'),
]


@pytest.mark.parametrize('structure_id,element,label', UNREPRESENTABLE)
def test_gate4_previously_unrepresentable_centres_are_candidates(records, structure_id,
                                                                 element, label):
    """`stereo_units()` and not `stereogenic_units()`: gate 4 asks about CANDIDATES.

    Being in that list is being a place that could carry a configuration, which is exactly the
    question "can the model represent the centre at all".  Whether the candidate survives the
    automorphism group is the next question and a different one; the three records where the isotope
    decides it are `test_gate4_the_isotope_decided_centres_are_genuinely_stereogenic`.
    """
    m = records[structure_id][0]
    anchors = [u['anchor'] for u in m.stereo_units() if u['kind'] == SU_TETRA]
    assert any(m.element_of(a) == element for a in anchors), \
        f'{label} centre on {structure_id} is still unrepresentable'


@pytest.mark.parametrize('structure_id,element', [('VS180', 6), ('VS183', 16), ('VS186', 16)])
def test_gate4_the_isotope_decided_centres_are_genuinely_stereogenic(records, structure_id,
                                                                     element):
    """The three records where the ISOTOPE is the only thing that makes the centre real.

    Candidacy alone is insensitive to the isotope here, and that is worth stating rather than
    assuming: bromochloromethane's carbon has four direction SLOTS whether or not one hydrogen is
    deuterium, so VS180 passes gate 4 above even from a bridge that dropped the isotope entirely.
    The verdict does not -- drop the isotope and the two hydrogens become interchangeable, an
    automorphism swaps them oddly, and the centre is refused.  So this test is where the isotope
    becomes load-bearing at the level a chemist cares about.

    The discriminating property is asserted, not narrated: each anchor must have two directions that
    agree on element and DIFFER in isotope, which is the only thing separating these centres from
    their achiral parents.
    """
    m = records[structure_id][0]
    units = [u for u in m.stereogenic_units()
             if u['kind'] == SU_TETRA and m.element_of(u['anchor']) == element]
    assert len(units) == 1, f'{structure_id}: one stereogenic centre on element {element}'

    refs = [r for r in units[0]['refs'] if r is not None]
    by_element = Counter(m.element_of(r) for r in refs)
    twins = [e for e, c in by_element.items() if c > 1]
    assert twins, 'the centre has two like directions, or the isotope decides nothing'
    assert any(len({m.isotope_of(r) for r in refs if m.element_of(r) == e}) > 1 for e in twins), \
        'and they differ in isotope, which is the property under test'


def test_gate4_phosphine_lone_pair_stays_excluded(records):
    """VS138 is a phosphine, and the core deliberately does NOT find the unit the oracle annotates.

    The file annotates VS138 `TH`.  Refusing it is the user's scope decision -- "forget about amines
    and phosphines, it only introduces noise" -- and not a bug: a P(III) lone-pair centre inverts at
    room temperature, so a toolkit that reported it would flood every phosphine ligand with
    stereocentres nobody can isolate.  This gate certifies the decision rather than the annotation.

    RULING F83, twice over.  The subject is checked first, because a negative assertion over units
    that do not exist proves nothing; and the assertion is that NO kind-0 unit is anchored on that
    phosphorus AT ALL, rather than the narrow conjunction `kind == 0 and element == 15 and
    n_refs == 4 and degree == 3` the brief proposed -- that one lets a three-ref tetrahedral unit on
    the same atom through, and tests an implementation detail where the scope decision is the point.
    """
    m = records['VS138'][0]
    phosphorus = [a for a in m.atom_numbers if m.element_of(a) == 15]
    assert len(phosphorus) == 1, 'one phosphorus, which is the subject of the assertion below'
    p = phosphorus[0]
    assert (m.charge_of(p), m.degree_of(p), m.total_h_of(p)) == (0, 3, 0), \
        'neutral, three-coordinate, no hydrogen: a P(III) lone-pair centre and not a phosphonium'
    assert 'TH' in _codes(records['VS138'][1]), 'the oracle does annotate it, which is the point'

    units = m.stereo_units()
    assert units, 'the record does produce candidates, so the negative below is not vacuous'
    assert not [u for u in units if u['anchor'] == p], \
        'no unit of any kind is anchored on the phosphine phosphorus'


def test_gate4_helicity_is_still_not_perceived(records):
    """VS010 and VS011 are [6]helicene of opposite helicity. Kind 4 is RESERVED and never produced.

    RULING F83: `all(u['kind'] != 4 for u in [])` is true of an empty list, and this is exactly that
    case -- both records return `[]`, measured.  A 26-carbon [6]helicene has no tetrahedral centre,
    no isolable cis/trans bond, and no atropisomer axis either, since every one of its pivot bonds
    is a ring-fusion bond and `_is_atropisomer_axis` refuses a ring bond.  Emptiness is the CORRECT
    answer today; helicity is deferred by choice.

    So the measured shape is asserted -- `== []`, not a `!= 4` that empties trivially -- next to a
    skeleton pin, so a bridge that produced an empty MOLECULE could not pass this by accident.  The
    corpus-wide statement that kind 4 never appears is asserted here too, where it is not vacuous.
    """
    for sid in ('VS010', 'VS011'):
        m, meta = records[sid]
        assert _codes(meta) == {'HE'}, f'{sid} is the helicity annotation'
        assert (m.atom_count, m.bond_count, m.rings_count) == (26, 31, 6), \
            f'{sid} is C26 with six fused rings, so the empty answer below is about a real molecule'
        assert {m.element_of(a) for a in m.atom_numbers} == {6}
        assert m.stereo_units() == [], f'{sid} produces no candidate at all, which is correct today'

    assert not [(sid, u) for sid, (m, _) in records.items()
                for u in m.stereo_units() if u['kind'] == SU_HELICAL], \
        'kind 4 is reserved across the whole corpus, where the statement has 300 subjects'


def test_the_two_unannotated_records_have_nothing_to_annotate(records):
    """VS002 and VS004 carry no STEREOGENIC_UNITS field, and here that means "no unit exists".

    The absence is ambiguous in general -- it could equally mean "not annotated" -- so the chemistry
    decides rather than the field.  VS002 is cyclohepta-1,3,5-triene (C7H8, one 7-ring, three ring
    double bonds) and VS004 is 1,2-dihydronaphthalene (C10H10, two rings, four double bonds).  In
    both, every double bond lies in a ring far too small for its two configurations to be separable,
    and every sp3 carbon is a ring CH2 whose two hydrogens are indistinguishable.  There is nothing
    to annotate, so the assertion is on the VERDICT and the skeleton, and it is a real one.

    What is deliberately NOT asserted: that `stereo_units()` is empty.  It is not -- each CH2 is a
    candidate with two unnamed directions in one direction list -- and that is right.  Candidacy is
    constitution; the two unnamed directions are what the verdict then refuses.
    """
    for sid, atoms, bonds, rings, doubles in (('VS002', 7, 7, 1, 3), ('VS004', 10, 11, 2, 4)):
        m, meta = records[sid]
        assert 'STEREOGENIC_UNITS' not in meta, f'{sid} is one of the two unannotated records'
        assert (m.atom_count, m.bond_count, m.rings_count) == (atoms, bonds, rings), sid
        assert {m.element_of(a) for a in m.atom_numbers} == {6}, f'{sid} is a hydrocarbon'
        assert sum(1 for _, _, o in _edges(m) if o == 2) == doubles, sid

        units = m.stereo_units()
        assert units, f'{sid} does produce candidates, so the emptiness below is a verdict'
        for u in units:
            assert u['kind'] == SU_TETRA
            assert m.total_h_of(u['anchor']) == 2, 'a ring CH2: two indistinguishable hydrogens'
            assert bin(u['unnamed_mask']).count('1') == 2, 'two unnamed directions in one list'
        assert m.stereogenic_units() == [], f'{sid} has no stereogenic unit, which is the annotation'


# --- Ruling F103: the invariant this corpus is blind to --------------------------------

def test_the_corpus_multi_component_record_answers_component_by_component(records):
    """VS186, the only Blue Book record whose units live in two components -- certified by name.

    Task 6b restricts each unit's stereogenicity witness search to its anchor's own component by
    pinning the complement, and the pin array is reused across units.  A complement pin left
    standing pins the NEXT unit's own component too, the whole record becomes the identity, the
    identity is even, no witness is found, and the unit is wrongly marked -- silently, with
    `stereo_truncated` still False.

    VS186 is a tetrabutylammonium arylsulfonate: 17 candidates in the ammonium component, none of
    them stereogenic, and one stereogenic sulfur in the sulfonate.  That is the shape the defect
    lives in, and asserting it by IDENTITY rather than by count matters, because the defect invents
    a centre in a DIFFERENT component from the one that has any.  This record does not itself change
    under the leak, which is precisely why ruling F103 also demands the synthetic record below.
    """
    m, _ = records['VS186']
    labels = m.component_labels()
    assert m.connected_components_count == 2
    assert sorted(Counter(labels.values()).values()) == [11, 17]

    per_component = Counter(labels[u['anchor']] for u in m.stereo_units())
    assert sorted(per_component.values()) == [1, 17], 'units on both sides of the salt'

    chiral = m.chiral_atoms()
    assert len(chiral) == 1 and m.chiral_bonds() == {}
    centre = next(iter(chiral))
    assert m.element_of(centre) == 16, 'the sulfonate sulfur, and nothing in the ammonium'
    assert per_component[labels[centre]] == 1, \
        'the marked component holds one unit; the seventeen next door are all refused'
    for u in m.stereo_units():
        assert u['stereogenic'] is (u['anchor'] == centre)


@pytest.mark.parametrize('names,marked_block', [
    (('dimethylcyclohexane', 'methylcyclohexane'), 0),
    (('methylcyclohexane', 'dimethylcyclohexane'), 1),
])
def test_a_refused_unit_in_a_non_first_component_stays_refused(names, marked_block):
    """RULING F103: the acceptance gate's counterpart to Task 6b's mechanism fixture.

    `test/stereo.sdf` cannot express the component-pin leak -- 294 of its 300 records are
    single-component, five of the remaining six carry a one-atom counter-ion, and VS186 does not
    change under the leak -- so 300 clean records certify nothing about the restriction the branch
    now rests on.  Corpus size is not sensitivity.  This is the minimal witness: 1,4-dimethyl-
    cyclohexane plus methylcyclohexane as two disjoint components, which answers TWO chiral atoms
    correctly and THREE under the leak, by inventing one on the methylcyclohexane ring.

    Both orders, because the leak's first victim is the second component processed, so the
    refusal-first order survives it.  The two fragments differ in exactly one thing -- whether the
    ring's arm swap is odd at a second CH and therefore refutes itself -- so one is refused by the
    witness search and the other survives it, and the property under test is the only thing
    separating them.

    Not redundant with `test_stereo_perception.py`'s fixture: that one tests the mechanism, this one
    certifies that the acceptance gate can see it.  The fragments and the builder are imported from
    there rather than rebuilt, so the two cannot drift apart.
    """
    m, blocks = _disjoint(*names)
    where = {sid: (b, i) for b, sids in enumerate(blocks) for i, sid in enumerate(sids)}
    assert m.connected_components_count == 2, 'two components, not one fused record'
    assert m.atom_count == 15, 'far below the node budget, so the flag below is decided, not a guess'
    assert m.stereo_truncated is False, 'these marks are proven, not taken conservatively'

    # BY IDENTITY, and anchor-free per ruling F101: which end of a unit is keyed in SEG_PARITY is a
    # slot-order artifact, so nothing here keys on an anchor or reaches a unit through `unit_of`.
    assert {where[sid] for sid in m.chiral_atoms()} == {(marked_block, 1), (marked_block, 4)}, \
        'both centres on the dimethyl ring, and none on the ring next door'
    assert m.chiral_bonds() == {}, 'neither fragment carries a bond-kind unit'

    other = 1 - marked_block
    assert any(b == other for b, _ in where.values()), 'the refused component is really in there'
    assert not [sid for sid, key in where.items() if key[0] == other and m.is_chiral(sid)], \
        'the methylcyclohexane CH is a candidate and stays refused wherever it sits'


# --- No record may crash or hang -------------------------------------------------------

def test_no_record_raises(records):
    """Every reader runs on all 300 records, and validation is clean by construction: no parity
    crosses the bridge, so nothing can be reported as unjustified. A non-empty report here means
    the bridge grew a stereo write.

    `canonical_stereo_groups()` can raise `AutomorphismBudgetExceeded`, and nothing here catches it:
    swallowing that into an "unknown stereo" fallback would make distinct mixtures compare equal on
    exactly the hard molecules, which is the failure gate 1 exists to catch.  It is called in
    `test_no_parity_crosses_the_bridge`; if a record ever raises there, that is a finding to report
    with the record id, and a tension with ruling F62 for the final review to resolve.
    """
    for sid, (m, _) in records.items():
        assert m.stereo_units() is not None
        assert len(m.canonical_order()) == len(m.atom_numbers)
        assert m.validate_stereo() == [], sid
