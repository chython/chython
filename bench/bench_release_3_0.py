# -*- coding: utf-8 -*-
"""The numbers the 3.0 release notes quote.  Not source; see CLAUDE.md.

    python bench/bench_release_3_0.py

Runs each operation on this tree and on an installed chython 2.24, and prints one row per operation
with both times and their ratio.  chython 2 answers in another interpreter under `-I`, the same
isolation `chython/core/test/oracle.py` uses and for the same reason: the two libraries never share a
`sys.modules`, so neither can be measured against a copy of itself.

WHAT A ROW MEANS.  A pipeline number, not a microbenchmark of one function: `parse + canonicalize`
includes the parse, because that is what a caller pays to get a canonical form out of a string.  Each
operation runs over the whole corpus `REPEATS` times and the FASTEST pass is reported -- a slower pass
measured the machine, not the library.  `standardize` subtracts the parse it needed, so a negative
number there would mean the subtraction is noise and the row is not reportable.

Both interpreters must be the same Python: 2.24 is pure Python where 3.0 is compiled, so a version
difference between the two sides lands entirely on chython 2's side of the ratio.  The header prints
both versions and refuses to print a ratio when they differ.
"""
from json import dumps, loads
from os import environ
from pathlib import Path
from platform import python_version
from re import search
from subprocess import run
from sys import argv, executable, path as _path
from time import perf_counter


ROOT = Path(__file__).resolve().parent.parent

# `sys.path[0]` is `bench/`, so an installed chython answers instead of the tree -- and this machine has
# one, so without this line the 3.0 column would measure chython 2 against chython 2.  The oracle child
# must NOT get this: it runs under `-I` to import the INSTALLED chython 2, which is the whole point.
if '--measure' not in argv[1:] and str(ROOT) not in _path:
    _path.insert(0, str(ROOT))


#: Public compounds, a spread of sizes: nothing internal is measured, and a corpus of ethanol would
#: report the call overhead rather than the algorithms.
COMPOUNDS = [
    'CCO',                                                          # ethanol
    'CC(=O)Nc1ccc(O)cc1',                                           # paracetamol
    'CN1C=NC2=C1C(=O)N(C)C(=O)N2C',                                 # caffeine
    'CC(C)Cc1ccc(cc1)C(C)C(=O)O',                                   # ibuprofen
    'CC(=O)Oc1ccccc1C(=O)O',                                        # aspirin
    'OCC1OC(O)C(O)C(O)C1O',                                         # glucose
    'CN1CCC[C@H]1c1cccnc1',                                         # nicotine
    'CC(C)(C)NC[C@H](O)c1ccc(O)c(CO)c1',                            # salbutamol
    'Clc1ccccc1C1=NCC(=O)Nc2ccc(Cl)cc21',                           # a benzodiazepine
    'CC1=C(C(=O)Nc2ccccc2)S(=O)(=O)c2ccccc21',
    'C[C@]12CC[C@H]3[C@@H](CC[C@@H]4CC(=O)CC[C@]34C)[C@@H]1CCC2=O',  # a steroid skeleton
    'CN(C)CCCN1c2ccccc2CCc2ccccc21',                                # imipramine
    'OC(=O)c1ccccc1Nc1ccccc1',                                      # fenamic acid
    'CC(C)NCC(O)COc1ccccc1OCC=C',                                   # oxprenolol
    'Nc1nc(=O)n([C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)cc1',             # cytidine
    'CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O',      # penicillin G
    'COc1cc2c(cc1OC)C(=O)c1ccccc1C2',
    'c1ccc2c(c1)ccc1c2ccc2c1cccc2',                                 # a fused polyarene
    'O=C(O)[C@@H](N)Cc1c[nH]c2ccccc12',                             # tryptophan
    'COc1ccc2cc(ccc2c1)[C@@H](C)C(=O)O',                            # naproxen
    'CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C',  # cholesterol
]

#: A mixture, because the canonical form is computed per component and a salt or a formulation is
#: where that shows.  One string with this many dots, parsed and canonicalized as one record.
MIXTURE_PARTS = 40

#: How many times each operation walks the corpus.  The fastest pass is the answer, so this buys
#: confidence that some pass ran without the machine interfering rather than an average of interference.
REPEATS = 15

#: The mixture is ONE record, so a pass over it is a single timing and the machine shows through -- at
#: `REPEATS` its best of 15 moved by 2x between runs.  Its own count, high enough that the minimum is
#: reproducible, which is the only form a quotable number comes in.
MIXTURE_REPEATS = 60


def _best(call, corpus, repeats=None) -> float:
    """Microseconds per record for `call` over `corpus`, from the fastest of `repeats` passes."""
    best = None
    for _ in range(repeats or REPEATS):
        started = perf_counter()
        for record in corpus:
            call(record)
        elapsed = perf_counter() - started
        if best is None or elapsed < best:
            best = elapsed
    return best / len(corpus) * 1e6


def _version(module) -> str:
    """`module`'s version, from the tree's `pyproject.toml` when it IS the tree and metadata otherwise.

    Not metadata alone: this machine has chython 2 installed, so in the parent `importlib.metadata`
    answers about the install rather than about the checkout `sys.path` puts first.
    """
    if Path(module.__file__).resolve().is_relative_to(ROOT):
        return search(r"(?m)^version = '([^']+)'", (ROOT / 'pyproject.toml').read_text()).group(1)
    return __import__('importlib.metadata', fromlist=['version']).version('chython')


def measure() -> dict:
    """Every operation, timed against whichever chython this interpreter imports."""
    import chython

    version = _version(chython)
    parse = chython.smiles if version.startswith('2') else chython.read_smiles

    parsed = [parse(s) for s in COMPOUNDS]
    mixture = '.'.join(COMPOUNDS[:4] * (MIXTURE_PARTS // 4))
    # compiled once and asked of every record, which is how a filter is actually run.  Two queries that
    # hit and one that misses, so the row is not a measurement of the early exit alone
    read_query = chython.smarts if version.startswith('2') else chython.read_smarts
    queries = [read_query(q) for q in ('[N;D2]C(=O)', '[C;z1]OC(=O)', '[S;D2][C;a]')]

    def canonical(s):
        molecule = parse(s)
        molecule.canonicalize()
        return molecule.smiles

    def standardized(molecule):
        clone = molecule.copy()
        clone.standardize()
        return clone

    out = {'version': version, 'python': python_version(), 'path': chython.__file__,
           'parse': _best(parse, COMPOUNDS),
           'parse + canonicalize': _best(canonical, COMPOUNDS),
           'linear fingerprint': _best(lambda m: m.linear_fingerprint(), parsed),
           'morgan fingerprint': _best(lambda m: m.morgan_fingerprint(), parsed),
           'three SMARTS asked of a record': _best(lambda m: [q < m for q in queries], parsed)}
    # a COPY is subtracted and not a parse: the pass needs a fresh molecule per call, and subtracting
    # the parse left a residue big enough to move the row by 2x between runs.  It is also the one row
    # whose two sides are not the same work -- 2.24's `standardize` is the whole repair pipeline where
    # 3.0's is one pass of it, which `canonicalize` above runs in full on both sides
    out['standardize (copy subtracted)'] = (_best(standardized, parsed)
                                           - _best(lambda m: m.copy(), parsed))
    out[f'{MIXTURE_PARTS}-component mixture, parse + canonicalize'] = _best(
        lambda s: canonical(s), [mixture], MIXTURE_REPEATS) / 1e3      # milliseconds, one record
    return out


def _oracle() -> str | None:
    """The interpreter with chython 2 installed, or `None` when nothing provisioned one."""
    if (given := environ.get('CHYTHON2_ORACLE')):
        return given
    candidate = Path.home() / '.cache/chython2-oracle/bin/python'
    return str(candidate) if candidate.is_file() else None


def main() -> int:
    if '--measure' in argv[1:]:                      # the child: one JSON object on stdout
        print(dumps(measure()))
        return 0

    three = measure()
    if (oracle := _oracle()) is None:
        print('no chython 2 to compare against; provision one as `oracle.py` documents')
        two = None
    else:
        # `-I` drops cwd, PYTHONPATH and user site together, so the child imports the INSTALLED
        # chython 2 and not this tree -- without it the comparison is this tree against itself
        done = run([oracle, '-I', str(Path(__file__).resolve()), '--measure'],
                   capture_output=True, text=True)
        if done.returncode:
            print(f'the oracle failed:\n{done.stderr.strip()[-2000:]}')
            return 1
        two = loads(done.stdout)

    print(f'chython {three["version"]} on Python {three["python"]}: {three["path"]}')
    if two is not None:
        print(f'chython {two["version"]} on Python {two["python"]}: {two["path"]}')
        comparable = two['python'] == three['python']
        if not comparable:
            print('the two Pythons differ, so no ratio is printed: 2.24 is pure Python where 3.0 is '
                  'compiled, and the version difference would land on 2.24\'s side of it')
    print()
    unit = {f'{MIXTURE_PARTS}-component mixture, parse + canonicalize': 'ms'}
    width = max(len(k) for k in three if k not in ('version', 'python', 'path'))
    print(f'{"":{width}}  {"3.0":>12}  {"2.24":>12}   ratio')
    for key, value in three.items():
        if key in ('version', 'python', 'path'):
            continue
        suffix = unit.get(key, 'us')
        row = f'{key:{width}}  {value:9.2f} {suffix:2}'
        if two is not None:
            other = two[key]
            row += f'  {other:9.2f} {suffix:2}'
            if comparable:
                row += f'  {other / value:6.1f}x' if value > 0 else '        --'
        print(row)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
