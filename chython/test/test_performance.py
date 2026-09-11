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
"""THE RELEASE CLAIM, MECHANIZED.  chython 3 exists for speed, and "faster than RDKit" is a claim about
measurable operations on a stated corpus -- so it lives in the test suite, next to every other claim the
tree makes about itself, rather than in a throwaway script whose numbers nobody can reproduce.

Run it as a benchmark rather than as a gate::

    pytest chython/test/test_performance.py -s

`-s` lets the comparison table through; without it the assertions still run and the table is swallowed.

WHAT IS MEASURED, AND WHY NAIVELY TIMING A LOOP GETS IT WRONG.  Both toolkits memoize derived data on the
molecule object, and they memoize *different* things:

* chython caches the canonical identity (`_identity_cache`, keyed on the edit generation) and the SMILES
  string.  Measured here: a second pass over the same objects costs 0.03 us against 24 us for the first.
* RDKit caches the Crippen atom contributions on the mol (`Crippen.MolLogP`: 0.13 us on a second pass
  against 89 us on the first) and computes ring info lazily on first use.

So there is no per-operation classification a harness can hardcode -- "does this cache?" has a different
answer per toolkit, and a table of answers rots.  Two disciplines are used instead, and neither of them
subtracts one timing from another:

**COLD BATCH** -- the whole corpus is parsed into a list *outside* the timer, then the operation is
applied once to each object *inside* it.  Every object is therefore cold, every call is measured exactly
once, and no call can read a cache that an earlier call wrote.  Repeats rebuild the corpus, so a repeat
is as cold as the first pass.  The predecessor of this file timed `parse N objects + op` and subtracted a
pure-parse baseline; that works, but it spends the whole parse cost as measurement noise, which is fatal
for the operations that cost less than a percent of a parse (mass, TPSA).  Building outside the timer
costs nothing and bounds nothing.

**WHOLE PIPELINE** -- the timer covers a string going in and an answer coming out, with no object handed
in at all.  Used for `SMILES parse` (there is nothing to hand in) and for `dedup key`, where the thing a
caller actually pays for is *string to key*: RDKit's key is a canonical SMILES and chython's is
`canonical_bytes`, and quoting either without its parse would describe an operation nobody performs.

`min` over the repeats, not the mean: interference from the rest of the machine only ever adds time, so
the smallest observation is the closest one to the cost of the code.  **And the repeats of the two
toolkits alternate** -- see `_Bench` for why that is load-bearing rather than tidy.

WHAT IS ASSERTED.  Ratios with a wide margin, never absolute microseconds -- a CI box under load is
several times slower than an idle laptop and a test that fails there is worthless.  Each row's docstring
states the ratio measured on the development machine, and each floor is roughly a third of it, so a row
survives chython losing three quarters of its relative advantage before it fails.  That is loose on
purpose: a benchmark that fails on a loaded machine gets disabled, and a disabled benchmark measures
nothing.

**The two rows chython does not win are asserted as ceilings, not as wins.**  TPSA is genuinely slower
and is known to be (deferred to 3.1); a ceiling there catches a regression while leaving an improvement
-- even one that overtakes RDKit -- free to pass, and the printed table is where the improvement shows.
Pinning a loss as a loss would make fixing it a test failure.

ONE MEASUREMENT NOTE ON MOLECULAR MASS: `float(mol)` is measured rather than `mol.molecular_mass`,
its other spelling, because the property adds a Python-level attribute lookup to 70 ns of arithmetic.
The row is asserted only as parity, since 70 ns is small enough that the comparison is really between
two toolkits' Python call overheads and not between two mass computations.
"""
from importlib.util import find_spec
from pathlib import Path
from time import perf_counter_ns

from pytest import fixture, mark, skip

# `read_smiles` and not `smiles`: the facade's door dispatches on the argument's type, and a
# benchmark of the PARSER must not charge it that isinstance chain per record.
from chython import SDFRead, inchi_library_loaded, molecule_to_inchi, read_smiles as smiles, smarts


# InChI is a separate optional piece of the build -- `build_inchi.py` skips when cmake or the submodule
# are absent -- so its row skips on the same terms as every other `needs_inchi` in the tree.
needs_inchi = mark.skipif(not inchi_library_loaded(), reason='libinchi not loaded')

# numpy is optional too (`chython[ml]`), and the fingerprint row is the one measurement here whose
# CHYTHON side needs it.  In practice the `corpus` fixture skips this file first on such an install --
# RDKit's own wheel requires numpy, so a checkout with one and not the other is not something pip can
# produce -- but the row states its own dependency rather than inheriting a skip from the reference
# toolkit's packaging.  `find_spec` rather than `importorskip`, for the reason `corpus` gives.
needs_numpy = mark.skipif(find_spec('numpy') is None,
                          reason='numpy is not installed; morgan_fingerprint answers an array')


# The root is found by looking for `pyproject.toml`, not by counting `parents[N]`: a count is a second
# fact about where this file sits and it is wrong the moment the file moves.  Same reasoning, same
# spelling, as `test_stereo_bluebook.py`.
def _repo_root():
    for candidate in Path(__file__).resolve().parents:
        if (candidate / 'pyproject.toml').is_file():
            return candidate
    raise RuntimeError('cannot locate the repository root: no pyproject.toml above this file')


#: The corpus, four tracked files chosen so that no single kind of structure dominates it: 300 IUPAC
#: Blue Book stereochemistry examples, 37 polycycles, 73 arenes and 2 peptides.  Aromaticity, ring
#: fusion, stereocentres and long chains are all represented, and the mean is 18 heavy atoms -- which
#: is the size a benchmark of a cheminformatics toolkit should be about, drug-like rather than either
#: a two-atom microbenchmark or a protein.
CORPUS_FILES = ('stereo.sdf', 'cycle.sdf', 'arenes.sdf', 'peptide.sdf')

#: Repeats per measurement.  The expensive rows get fewer: InChI alone costs ~85 us per molecule per
#: toolkit, so five repeats of both sides of that row would be most of the file's runtime for a
#: measurement that is already stable at three.
REPEATS = 5
REPEATS_SLOW = 3

#: Substructure patterns, spelled once per toolkit.  Three of them because a single pattern measures one
#: point in the matcher's behaviour: a rare terminal group that fails fast, a common carbonyl, and an
#: aromatic bond that hits early and often.  The chython spellings carry map numbers because that is how
#: the tree writes a query; RDKit's are the nearest equivalent in Daylight SMARTS.
#:
#: Note the `:` on the aromatic pattern.  A chython SMARTS bond written implicitly matches single bonds
#: ONLY, aromatic bonds not included, so `[C;a][C;a]` would measure a pattern that never matches.
PATTERNS = (('primary amine', '[N;D1;z1;x0:1][C;z1:2]', '[NX3;H2][CX4]'),
            ('carbonyl', '[O;z2;x0:2]=[C;z2:1]', '[OX1]=[CX3]'),
            ('aromatic bond', '[C;a:1]:[C;a:2]', 'c:c'))

#: Rows of the printed table, filled as the tests run: `(operation, rdkit_us, chython_us, floor)`.
_TABLE = []

#: Records in the comparison corpus, so the printed table states what it was measured over.
_MEASURED_OVER = []


@fixture(scope='module')
def chython_corpus():
    """The corpus as canonical SMILES strings, chython's own reading of it and nothing else.

    SMILES rather than the SDF records themselves, so that the two toolkits are handed *the same input*
    for the parse row and equal molecules for every other row.  Canonicalized first because the files
    hold Kekule CTABs: without it chython would be timed on aromatic input and RDKit on Kekule input for
    the aromatize row, which measures the corpus and not the code.

    Separate from `corpus` because the chython-against-chython row -- the `canonical_order()` tail --
    must run on a machine with no RDKit installed, and a fixture that imports RDKit to filter would take
    it down with the comparison rows.
    """
    root = _repo_root() / 'test'
    out = []
    for name in CORPUS_FILES:
        with SDFRead(root / name) as f:
            for m in f:
                m.canonicalize()
                out.append(str(m))
    assert len(out) > 300, f'corpus collapsed to {len(out)} records; a reader is broken, not slow'
    return tuple(out)


@fixture(scope='module')
def corpus(chython_corpus):
    """The subset both toolkits parse -- and the skip guard for every comparison in this file.

    RDKit is an optional extra (`pip install chython[rdkit]`) and must never become anything more: the
    standing ruling is "NO DEPS ON OTHER TOOLKITS", and a benchmark that made the reference toolkit
    mandatory would smuggle one in through the test suite.  The guard is a `skip()` inside the fixture
    rather than a mark on eleven tests, which is the idiom `chython/interop/test/conftest.py` uses for
    the same problem; `find_spec` rather than `importorskip` at module scope, because importing RDKit
    costs a noticeable fraction of a second on every run of the whole suite, including runs that select
    nothing here.

    The three records RDKit declines (a ruthenium cluster, a closo-borane cage, and one record whose
    aromatic ring RDKit will not accept) are dropped rather than tolerated.  A benchmark corpus has to
    be one corpus; a molecule only one side can read cannot appear in a ratio.
    """
    if find_spec('rdkit') is None:
        skip('rdkit is not installed')
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog('rdApp.*')
    out = tuple(s for s in chython_corpus if Chem.MolFromSmiles(s) is not None)
    _MEASURED_OVER.append(out)
    return out


class _Bench:
    """One row of the table: two timing disciplines, INTERLEAVED, in microseconds per molecule.

    THE REPEATS ALTERNATE BETWEEN THE TOOLKITS, and that is not cosmetic.  Timing all of RDKit's
    repeats and then all of chython's takes `min` over two *different* windows of the machine's life,
    so a load spike that lands in the second window is read as a slowdown in whichever toolkit was
    measured there.  Observed on this branch while three other agents were building: chython's parse
    came out at 26.6 us against its true 8.5, RDKit's unaffected at 72, and the row failed for a
    reason that had nothing to do with either toolkit.  Alternating puts both sides in the same
    windows, so `min` picks each one's cleanest pass out of the same weather.
    """
    __slots__ = ('corpus', '_rdkit_parse')

    def __init__(self, corpus):
        from rdkit import Chem

        self.corpus = corpus
        self._rdkit_parse = Chem.MolFromSmiles

    def _cold_pass(self, build, op):
        """One pass: build the whole corpus OUTSIDE the timer, then apply `op` once to each object."""
        objects = [build(x) for x in self.corpus]
        started = perf_counter_ns()
        for o in objects:
            op(o)
        return perf_counter_ns() - started

    def _whole_pass(self, fn):
        """One pass: string in, answer out, everything timed -- there is no object to hand in."""
        started = perf_counter_ns()
        for x in self.corpus:
            fn(x)
        return perf_counter_ns() - started

    def _run(self, name, floor, rdkit_pass, chython_pass, repeats):
        rdkit_ns = chython_ns = None
        for _ in range(repeats):
            spent = rdkit_pass()
            if rdkit_ns is None or spent < rdkit_ns:
                rdkit_ns = spent
            spent = chython_pass()
            if chython_ns is None or spent < chython_ns:
                chython_ns = spent

        scale = len(self.corpus) * 1000
        rdkit_us, chython_us = rdkit_ns / scale, chython_ns / scale
        ratio = rdkit_us / chython_us
        _TABLE.append((name, rdkit_us, chython_us, floor))
        assert ratio >= floor, (f'{name}: chython {chython_us:.2f} us vs RDKit {rdkit_us:.2f} us is '
                                f'{ratio:.2f}x, below the asserted floor of {floor}x')

    def cold(self, name, floor, rdkit_op, chython_op, repeats=REPEATS):
        """Compare `op` on cold objects.  The parsers are baked in: every row uses the same two."""
        self._run(name, floor, lambda: self._cold_pass(self._rdkit_parse, rdkit_op),
                  lambda: self._cold_pass(smiles, chython_op), repeats)

    def whole(self, name, floor, rdkit_fn, chython_fn, repeats=REPEATS):
        """Compare two string-in/answer-out pipelines."""
        self._run(name, floor, lambda: self._whole_pass(rdkit_fn),
                  lambda: self._whole_pass(chython_fn), repeats)


@fixture(scope='module')
def bench(corpus):
    return _Bench(corpus)


@fixture(scope='module', autouse=True)
def _print_table():
    """Print the comparison table once the module's measurements are in.  Visible under `-s`."""
    yield
    if not _TABLE:
        return
    from rdkit import rdBase

    size = len(_MEASURED_OVER[0]) if _MEASURED_OVER else 0
    print(f'\n\nchython 3 against RDKit {rdBase.rdkitVersion}, {size} molecules from '
          f'{", ".join(CORPUS_FILES)}')
    print(f'{"operation":<28} {"RDKit us":>9} {"chython us":>11} {"ratio":>8}   {"asserted":>12}')
    print('-' * 74)
    for name, rdkit_us, chython_us, floor in _TABLE:
        print(f'{name:<28} {rdkit_us:9.2f} {chython_us:11.2f} {rdkit_us / chython_us:7.2f}x   '
              f'{"ratio >= " + format(floor, "g"):>12}')
    print('-' * 74)
    print('us per molecule, min of the repeats.  ratio > 1 means chython is faster.')


# ---------------------------------------------------------------------------------------------------
# The methodology self-check.  It runs first because every number below it is only meaningful if it
# holds.
# ---------------------------------------------------------------------------------------------------

def test_a_second_pass_is_a_cache_read(bench, corpus):
    """The memoization this file is built around is real, and is what forbids reusing objects.

    Asserted as `any`, not as a fact about a named operation: either toolkit is free to drop a cache,
    and this test must then keep measuring the discipline rather than failing for an improvement. What
    it will not survive is somebody "simplifying" `cold()` into building the corpus once -- which is the
    mistake it exists to catch, because that edit would silently turn most rows below into a comparison
    of two dictionary lookups.
    """
    from rdkit import Chem
    from rdkit.Chem import Crippen

    probes = [('chython canonical_bytes', smiles, lambda m: m.canonical_bytes),
              ('chython SMILES write', smiles, str),
              ('RDKit Crippen.MolLogP', Chem.MolFromSmiles, Crippen.MolLogP)]

    speedups = []
    for name, build, op in probes:
        objects = [build(x) for x in corpus]
        started = perf_counter_ns()
        for o in objects:
            op(o)
        first = perf_counter_ns() - started
        started = perf_counter_ns()
        for o in objects:
            op(o)
        second = perf_counter_ns() - started
        speedups.append((name, first / max(second, 1)))

    assert any(s >= 5 for _, s in speedups), \
        f'no probed operation memoizes any more; second-pass speedups were {speedups}'


# ---------------------------------------------------------------------------------------------------
# The rows chython wins.  Every floor is roughly a third of the ratio measured on the development
# machine, so a row survives chython losing three quarters of its relative advantage before it fails --
# which is the headroom a CI box under an unknown load needs.  The floors are deliberately not tight:
# a benchmark that fails on a loaded machine gets disabled, and a disabled benchmark measures nothing.
# ---------------------------------------------------------------------------------------------------

def test_smiles_parse(bench):
    """Measured 7.8x (72.2 us / 9.2 us).  The hottest path in the library and the widest margin."""
    from rdkit import Chem

    bench.whole('SMILES parse', 2.5, Chem.MolFromSmiles, smiles)


@mark.parametrize('label,chython_pattern,rdkit_pattern', PATTERNS, ids=[p[0] for p in PATTERNS])
def test_substructure_match(bench, label, chython_pattern, rdkit_pattern):
    """Measured 8.1x-10.8x (~1.0 us / 0.10-0.13 us), pattern-dependent.

    The floor is the same for all three: a per-pattern floor would be fitting the assertion to the
    measurement, and the point of the row is that the matcher is an order of magnitude ahead wherever
    it is probed.
    """
    from rdkit import Chem

    query = smarts(chython_pattern)
    rdkit_query = Chem.MolFromSmarts(rdkit_pattern)
    assert rdkit_query is not None, f'RDKit rejected the reference pattern {rdkit_pattern!r}'

    bench.cold(f'substructure: {label}', 2.5, lambda m: m.HasSubstructMatch(rdkit_query),
               lambda m: query.is_substructure(m))


def test_aromatize(bench):
    """Measured 6.9x (45.7 us / 6.7 us).  `kekule()` then `thiele()`, the repair pipeline's aromatic half.

    RDKit's equivalent is `Kekulize(clearAromaticFlags=True)` then `SetAromaticity`: the same round trip
    out of and back into the aromatic form, on molecules that arrived aromatic.  Cold batch, because
    both operations mutate -- a second pass would find the work already done.
    """
    from rdkit import Chem

    def chython_op(m):
        m.kekule()
        m.thiele()

    def rdkit_op(m):
        Chem.Kekulize(m, clearAromaticFlags=True)
        Chem.SetAromaticity(m)

    bench.cold('aromatize (kekule+thiele)', 2.0, rdkit_op, chython_op)


def test_crippen_logp(bench):
    """Measured 3.3x (91.8 us / 27.9 us).

    The row that most needs cold batch: RDKit caches its Crippen contributions on the mol, so a reused
    object reports 0.13 us and RDKit appears 200x faster than chython at the same arithmetic.
    """
    from rdkit.Chem import Crippen

    bench.cold('logP (Crippen)', 1.4, Crippen.MolLogP, lambda m: m.crippen_logp)


def test_dedup_key(bench):
    """Measured 2.9x (97.6 us / 34.0 us).  String in, identity key out -- what a deduplicating pass pays.

    Whole pipeline rather than cold batch, and deliberately: both toolkits cache the key on the object,
    but more importantly a caller deduplicating a file has strings and not molecules, so the parse is
    part of the operation.  RDKit's key is its canonical SMILES; chython's is `canonical_bytes`, which
    is bytes and not a string for the reason its docstring gives.
    """
    from rdkit import Chem

    bench.whole('dedup key (string->key)', 1.4, lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s)),
                lambda s: smiles(s).canonical_bytes)


@needs_numpy
def test_morgan_fingerprint(bench):
    """Measured 2.5x (12.5 us / 5.0 us), radius 2 into 1024 bits on both sides.

    This compares cost and says nothing about the bits: RDKit's bits and chython's are not the same set.
    """
    from rdkit.Chem import rdFingerprintGenerator

    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)

    bench.cold('Morgan fingerprint', 1.2, generator.GetFingerprintAsNumPy,
               lambda m: m.morgan_fingerprint())


# ---------------------------------------------------------------------------------------------------
# Parity and the known-slower rows.  Ceilings, not wins: the assertion is that chython has not fallen
# further behind, and an improvement past RDKit passes untouched and shows up in the printed table.
# ---------------------------------------------------------------------------------------------------

@needs_inchi
def test_inchi(bench):
    """Measured 0.96x (93.5 us / 97.7 us): parity, and both sides are the same C library.

    An earlier hand measurement put this at 1.5x in chython's favour.  It is not reproduced here and the
    likely reason is the corpus: that one came from SDF records carrying 2D coordinates, which RDKit
    hands to the InChI library so that it derives double-bond geometry from them, while chython passes
    none.  That is a difference between the two callers rather than between the two toolkits, and it is
    the same asymmetry `interop/test/test_rdkit.py` strips before comparing InChI strings.  A corpus of
    SMILES removes it, and what is left is parity -- which is the honest number, since the work is
    libinchi's either way and only the marshalling is ours.
    """
    from rdkit import Chem

    bench.cold('InChI', 0.6, Chem.MolToInchi, molecule_to_inchi, REPEATS_SLOW)


def test_smiles_write(bench):
    """Measured 1.07x (25.9 us / 24.3 us): parity, and the ceiling says it must stay there.

    Writing is where chython is level with RDKit rather than ahead, which is expected: the canonical
    labelling dominates and both toolkits do the same amount of it.
    """
    from rdkit import Chem

    bench.cold('SMILES write', 0.6, Chem.MolToSmiles, str)


def test_tpsa_is_known_slower(bench):
    """Measured 0.27x -- chython is 3.7x slower here (1.8 us / 6.7 us).  Known, and deferred to 3.1.

    KNOWN-SLOWER AND PINNED AS SUCH.  The ceiling is 10x slower, so the row catches a regression while leaving
    every improvement free to pass, including one that overtakes RDKit; the printed table is where an
    improvement becomes visible.  Asserting the loss itself would make fixing it a test failure, which
    is how a benchmark suite starts defending the thing it measures.

    The cost is the table lookup: chython reads `tables/tpsa.tsv` contributions through the featurizer
    spine, where RDKit has the contributions compiled in.
    """
    from rdkit.Chem import rdMolDescriptors

    bench.cold('TPSA (known slower)', 0.1, rdMolDescriptors.CalcTPSA, lambda m: m.tpsa)


def test_molecular_mass(bench):
    """Measured 2.6x (0.22 us / 0.09 us) for `float(mol)`, asserted only as parity.

    Two reasons for the weak assertion.  Both numbers are well under a microsecond, so most of each is
    the interpreter's cost of one call and not the summation -- a ratio there is not a claim about
    either toolkit's arithmetic.  And the fastest RDKit spelling is chosen (`CalcExactMolWt` is faster
    than `Descriptors.MolWt`, a lambda around it), because a benchmark that quotes the slowest available
    spelling of the reference toolkit's answer is not measuring anything.
    """
    from rdkit.Chem import rdMolDescriptors

    bench.cold('molecular mass', 0.4, rdMolDescriptors.CalcExactMolWt, float)


def test_canonical_order_has_no_pathological_tail(chython_corpus):
    """`canonical_order()` must stay within a small factor of `canonical_bytes` on every molecule.

    A ratio rather than a microsecond figure, so it states something about the algorithm and not about
    the machine: the two share the extremal canonical labelling, so `canonical_order()` can only be a
    relabelled view of work `canonical_bytes` already does cheaply.  A tail of 100x is not the cost of
    canonicalizing a symmetric molecule, it is a search that fails to prune -- which is what it was
    before `_canon_order` began from the parity fold, and the worst ratio is now 1.13x.

    The worst molecule is found, not named: hardcoding today's worst record would let the pathology move
    to another one unnoticed.  Both sides are timed on a freshly parsed molecule because
    `canonical_bytes` memoizes and `canonical_order()` does not, so reusing one object would compare a
    computation against a cache read.  Molecules under 10 us either way are skipped -- their ratio is
    dominated by the timer (the worst is 6.9x, which would fail a 4x bound while saying nothing).
    """
    #: The lower bound on a measurement this test will draw a conclusion from, and the repeat count
    #: that makes each measurement a `min` rather than a single sample.
    floor_us, repeats, bound = 10, 3, 4

    def cost(op, source):
        best = None
        for _ in range(repeats):
            m = smiles(source)
            started = perf_counter_ns()
            op(m)
            spent = perf_counter_ns() - started
            if best is None or spent < best:
                best = spent
        return best / 1000

    worst = None
    measured = 0
    for s in chython_corpus:
        order_us = cost(lambda m: m.canonical_order(), s)
        bytes_us = cost(lambda m: m.canonical_bytes, s)
        if max(order_us, bytes_us) < floor_us:
            continue
        measured += 1
        ratio = order_us / bytes_us
        if worst is None or ratio > worst[0]:
            worst = (ratio, s, order_us, bytes_us)

    assert measured >= 20, (f'only {measured} molecules cost more than {floor_us} us either way; the '
                            f'bound below would be measuring almost nothing')
    ratio, s, order_us, bytes_us = worst
    assert ratio <= bound, (f'canonical_order() costs {order_us:.1f} us against canonical_bytes '
                            f'{bytes_us:.1f} us ({ratio:.0f}x) on {s}')
