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
"""chython 2 as an oracle, reached through a subprocess instead of an import.

WHY A SUBPROCESS AND NOT AN IMPORT.

chython 2 is a good second opinion.  What it must not be is a *dependency of this tree*: a test under
`chython/core/test/` that says `from chython.periodictable import ...` makes the V3 suite fail the
moment chython 2 is not in the working tree, and that coupling appears on no feature list.

`import chython.core` runs `chython/__init__.py`, which pulls in the whole of V2 -- so an in-tree
import of the oracle is not separable from the library under test even in principle.  Running the
oracle in another interpreter, against an installed chython 2, separates them completely: the two
libraries never share a `sys.modules`, and the oracle's version is pinned by whoever provisioned it
rather than by whatever the working tree happens to contain.

`test_no_chython_two_imports.py` enforces the distinction: a subprocess call to another interpreter
is allowed, `import chython.periodictable` is not.  Do not "simplify" this module back into an
import -- that reinstates the blocker it exists to remove.

PROVISIONING.  An isolated virtualenv with chython 2 installed:

    python3.10 -m venv ~/.cache/chython2-oracle
    ~/.cache/chython2-oracle/bin/pip install chython==2.24

Point `CHYTHON2_ORACLE` at another interpreter to override.  When neither exists every test that
needs it skips, because the oracle is an OPTIONAL test dependency: a machine that has not
provisioned it still runs the whole suite and every assertion that states an answer directly.

WHY THIS IS ONE MODULE AND NOT ONE PER CALLER.  A caller that spawns the oracle itself can forget
`-I`, and a missing `-I` does not break a test, it VOIDS one -- `python -c` puts the current directory
on `sys.path`, pytest runs from the repository root, so the child imports `./chython/` and the
differential compares the tree under test to itself and always agrees.
`importlib.metadata.version('chython')` finds `./chython.egg-info` and reports 3.0 for the same
reason, and a worktree with no egg-info passes where the main checkout fails -- so a green worktree
run does not clear it.  A failure mode that silently converts a differential into a tautology earns
exactly one definition, and this is it.  `run` is the only place in the repository that spawns the
oracle, and `ISOLATION` is the flag it always passes.

THE THREE GUARDS, all of them here so that no caller can forget one:

  1. `-I` on every invocation.  Not `-E`, not a scrubbed `env=`, not `cwd=` somewhere else -- `-I`
     drops cwd, PYTHONPATH and the user site directory together, and makes the answer independent of
     the directory pytest was invoked from.
  2. `VERSION` is pinned EXACTLY and asserted.  A differential against a moving target proves
     nothing, and a version bump has to be loud: some of what a newer V2 changed may be the very
     defects a caller has written down as expected divergences.
  3. The child's `chython.__file__` is checked against the repository root.  Guard 2 catches a leak
     only by inference -- a reader seeing `oracle is chython 3.0` has to work out why -- and this one
     says it directly and keeps saying it if the two version numbers ever coincide.

An ABSENT oracle skips; a WRONG one fails loudly.  That asymmetry is deliberate: skipping on a
version mismatch would make a provisioned-but-wrong oracle indistinguishable from an unprovisioned
one, which is how a stale pin survives.  `test_oracle.py` holds the negative control for guard 1 --
it drops `-I` and asserts that guard 3 then fires.

The oracle is a CORRECTNESS ORACLE AND NOT A CONTRACT.  A disagreement is a question, not a
verdict against the core -- where V2 diverges, the divergence is named at the site that works around
it.  Never "fix" V2 to make a comparison pass, and never freeze
a wrong V2 answer as an expected one: a live oracle hands out its wrong answers forever, so an
encoded one turns a temporary defect into a permanent specification.  Mark the case and report it.
"""
from __future__ import annotations

import json
import os
import pathlib
import subprocess
import sys


ENV_VAR = 'CHYTHON2_ORACLE'
DEFAULT_PATH = pathlib.Path.home() / '.cache' / 'chython2-oracle' / 'bin' / 'python'

# PINNED EXACTLY, not a floor: see guard 2 in the module docstring.
VERSION = '2.24'

# The tree under test. `parents[3]` is the repository root from `chython/core/test/oracle.py`.
ROOT = pathlib.Path(__file__).resolve().parents[3]

# Named rather than spelled inline so that `test_oracle.py` can drop it on purpose and demonstrate
# that the guard which depends on it actually fires. A constant nobody can vary is not a guard.
ISOLATION = '-I'

_RESOLVED = ...          # sentinel: not looked up yet
_PROBE = None


PREAMBLE = """\
import json, sys
_payload = json.loads(sys.stdin.read() or 'null')
def _emit(value):
    sys.stdout.write('\\x1e' + json.dumps(value))
"""

# `importlib.metadata` and not `chython.__version__`, which chython 2 does not define; and both
# facts in one run, because they are one question -- "which library did this interpreter import".
_PROBE_SCRIPT = ('import chython, json\n'
                 'from importlib.metadata import version\n'
                 'print(json.dumps([version("chython"), chython.__file__]))\n')


def run(*args, isolation=ISOLATION, **kwargs):
    """`subprocess.run` on the oracle interpreter. THE ONLY PLACE THE ORACLE IS SPAWNED.

    `isolation` exists to be dropped by the negative control in `test_oracle.py` and by nothing
    else.  Every real caller takes the default, which is why the default is the safe one.
    """
    python = require()
    kwargs.setdefault('capture_output', True)
    kwargs.setdefault('text', True)
    kwargs.setdefault('timeout', 600)
    return subprocess.run([str(python), *([isolation] if isolation else []), *args], **kwargs)


def interpreter():
    """The oracle interpreter path, or None -- EXISTENCE ONLY. Resolved once and cached.

    Deliberately does not check the version: that is `verify`, and it FAILS rather than returning
    None, so a wrong oracle cannot masquerade as an absent one.
    """
    global _RESOLVED
    if _RESOLVED is not ...:
        return _RESOLVED
    override = os.environ.get(ENV_VAR)
    candidate = pathlib.Path(override) if override else DEFAULT_PATH
    _RESOLVED = candidate if candidate.is_file() else None
    return _RESOLVED


def probe(isolation=ISOLATION):
    """`(version, chython.__file__)` as the oracle interpreter itself reports them.

    Cached for the default isolation, because it is one fact about the machine; an explicit
    `isolation` bypasses the cache so the negative control measures what it asked for.
    """
    global _PROBE
    if isolation != ISOLATION:
        return _read_probe(run('-c', _PROBE_SCRIPT, isolation=isolation))
    if _PROBE is None:
        _PROBE = _read_probe(run('-c', _PROBE_SCRIPT))
    return _PROBE


def _read_probe(out):
    if out.returncode:
        raise AssertionError('the chython 2 oracle could not import chython '
                             f'(exit {out.returncode}):\n{out.stderr.strip()}')
    version, source = json.loads(out.stdout.strip().splitlines()[-1])
    return version, source


def verify():
    """Guards 3 then 2, once per session. Raises -- an oracle that is present and wrong is not a skip.

    THE ORDER IS LOAD-BEARING and was measured the wrong way round first.  A leak trips BOTH guards:
    the child imports `./chython/` and reads `./chython.egg-info`, so it reports version 3.0.  With
    the version pin checked first, the one failure a reader ever sees for the commonest cause is
    "the oracle is 3.0", which describes a symptom and not the cause -- and the whole reason guard 3
    exists is to say the cause out loud.  So guard 3 goes first and guard 2 is what is left over:
    a genuinely mis-provisioned but genuinely separate chython 2.
    """
    version, source = probe()
    assert not source.startswith(str(ROOT)), (
        f'the oracle imported THE TREE UNDER TEST from {source}: the differential is comparing the '
        f'core to itself and cannot fail. `{ISOLATION}` is missing from an invocation somewhere.')
    assert version == VERSION, (
        f'the chython 2 oracle is {version}, every differential in this tree is pinned to {VERSION}. '
        f'Re-run the comparisons deliberately and move the pin -- a newer V2 is not automatically a '
        f'better oracle, since what it changed may be the divergences the callers expect.')


def require():
    """The interpreter, or skip the test. The oracle is optional; the suite is not."""
    from pytest import skip

    found = interpreter()
    if found is None:
        skip(f'the chython 2 oracle is not provisioned; see {__name__}.__doc__ '
             f'(set {ENV_VAR}, or install chython=={VERSION} into {DEFAULT_PATH.parent.parent})')
    return found


def requires_oracle(function):
    """Decorator form, for callers that want the check without taking a fixture."""
    from functools import wraps

    @wraps(function)
    def wrapper(*args, **kwargs):
        require()
        verify()
        return function(*args, **kwargs)
    return wrapper


def ask(source, payload=None, timeout=600, argv=()):
    """Run `source` in the oracle interpreter and return what it passed to `_emit`.

    `payload` is handed over as `_payload`, already decoded, and `argv` as `sys.argv[1:]`.  The
    answer travels as one JSON document after an ASCII record separator, so anything the oracle
    prints on its own -- and V2 prints deprecation warnings -- cannot be mistaken for the result.
    """
    verify()
    out = run('-c', PREAMBLE + source, *argv, input=json.dumps(payload), timeout=timeout)
    if out.returncode:
        raise AssertionError(f'the chython 2 oracle failed (exit {out.returncode}):\n'
                             f'{out.stderr.strip()}')
    if '\x1e' not in out.stdout:
        raise AssertionError('the oracle produced no answer; it must call `_emit(value)` exactly '
                             f'once. stdout was:\n{out.stdout.strip()}\n{out.stderr.strip()}')
    return json.loads(out.stdout.split('\x1e', 1)[1])


def ask_text(source, timeout=600):
    """`ask` for a caller that wants the oracle's plain stdout rather than a JSON document."""
    verify()
    out = run('-c', source, timeout=timeout)
    assert not out.returncode, (f'the chython 2 oracle failed (exit {out.returncode}):\n'
                                f'{out.stderr.strip()}')
    return out.stdout


# --- the persistent session ----------------------------------------------------------------------
#
# `ask` is enough when the questions are known before the first answer.  The SMILES writer's
# differential is not like that: it asks chython 2 for a stereo sign IN A FRAME THE CORE COMPUTED,
# so the V3 side has to run between one question and the next.  There are about 3,600 such rounds
# and a fresh interpreter costs 220ms, so one-shot calls would turn a 22-second file into a
# fifteen-minute one.
#
# So the oracle can also be a co-process: one chython 2 interpreter, alive for the session, holding
# V2 molecules behind integer handles and answering line-delimited JSON.  A round trip is about
# 100us, which is three orders of magnitude cheaper and makes an interactive differential affordable
# without moving any of chython 2's logic into this tree.
#
# THE MOLECULES STAY ON THEIR SIDE.  Nothing tries to reconstruct a V2 object here -- a handle is an
# integer, and every question about the molecule behind it is a method call over the pipe.  That is
# what keeps this a differential against chython 2 rather than a reimplementation of it.

SERVER = r'''
import json, sys, traceback

HANDLES = {}


def _put(mol):
    HANDLES[len(HANDLES)] = mol
    return len(HANDLES) - 1


def op_version(_):
    from importlib.metadata import version
    return version('chython')


def op_parse(payload):
    """A handle per string, `None` where V2 will not read it -- which is not an error here.

    A record chython 2 cannot parse has nothing to compare, so the caller drops it; raising would
    make one unreadable string in a 5,000-record corpus lose the whole corpus.
    """
    from chython import smiles
    out = []
    for text in payload['texts']:
        try:
            out.append(_put(smiles(text)))
        except Exception:
            out.append(None)
    return out


def _v2_sdf_reader():
    """chython 2's SDF reader BY ITS UNSHADOWED MODULE PATH, and the reason that matters.

    Never `chython.SDFRead`. In chython 3 the formats package registers its own reader over that name
    -- `files/__init__.py` imports it last, deliberately -- so following the package root would hand
    back the reader UNDER TEST. Every comparison in `formats/ctfile/test/` would then pass by identity,
    and the wedge suite's thousand-centre stereo oracle would be reading the arena's own parity bytes
    and calling the agreement evidence.

    That is the `-I` bug by a second route, and here it is impossible rather than merely avoided: this
    code runs inside chython 2's interpreter, where the V3 package does not exist to shadow anything.
    """
    from chython.files.SDFrw import SDFRead
    return SDFRead


def _v2_rdf_reader():
    """chython 2's RDfile reader by its unshadowed module path; see :func:`_v2_sdf_reader` for why."""
    from chython.files.RDFrw import RDFRead
    return RDFRead


def op_read_sdf(payload):
    """Whole SDF by path. `tolerant` turns a file V2 itself cannot read into `None`, not an error --
    such a file is simply not an oracle for itself."""
    SDFRead = _v2_sdf_reader()
    out = {}
    for name, path in payload['paths'].items():
        try:
            with SDFRead(path) as f:
                out[name] = [_put(mol) for mol in f]
        except Exception:
            if not payload.get('tolerant'):
                raise
            out[name] = None
    return out


def op_read_rdf(payload):
    """`{name: [{'reactants': n, 'products': n, 'agents': n, 'atoms': [...]}, ...]}` from chython 2.

    Counts, not handles.  `op_read_sdf` returns handles into this interpreter because the caller
    compares structures; this gate compares side assignment and component order, so numbers cross the
    pipe and nothing else.  A molecule record contributes `None` -- it has no sides to compare.

    `tolerant` maps a file V2 itself cannot read to `None`, exactly as in `op_read_sdf`.  It is NOT
    the default: swallowing every exception would turn a broken venv or a renamed V2 attribute into
    "V2 cannot read this file" on all four fixtures, and the gate would go green having compared
    nothing.
    """
    RDFRead = _v2_rdf_reader()
    out = {}
    for name, path in payload['paths'].items():
        try:
            with RDFRead(path) as f:
                records = []
                for r in f:
                    # `hasattr` and not `isinstance`: it needs no import, and the one import that
                    # would answer it is the package root this file is careful never to follow.
                    if hasattr(r, 'reactants'):
                        records.append({'reactants': len(r.reactants), 'products': len(r.products),
                                        'agents': len(r.reagents),
                                        'atoms': [len(m) for m in (*r.reactants, *r.reagents,
                                                                   *r.products)]})
                    else:
                        records.append(None)
            out[name] = records
        except Exception:
            if not payload.get('tolerant'):
                raise
            out[name] = None
    return out


def op_read_mdl_text(payload):
    """One MDL record given as TEXT rather than as a path, for a hand-written fixture.

    A test that writes down the parity it expects is testing that the implementation has not changed;
    a test that asks chython 2 is testing that the two stacks agree, which is the only statement about
    a sign convention worth making. So the fixtures need a reader that takes a string.
    """
    from io import StringIO
    SDFRead = _v2_sdf_reader()
    out = []
    for text in payload['texts']:
        if not text.endswith('\n'):
            text += '\n'
        try:
            with SDFRead(StringIO(text + '$$$$\n')) as f:
                out.append(_put(next(iter(f))))
        except Exception:
            out.append(None)
    return out


def op_record(payload):
    """Everything about a molecule that does not depend on a question, in one document."""
    out = []
    for handle in payload['handles']:
        mol = HANDLES[handle]
        atoms = []
        for n, atom in mol.atoms():
            atoms.append({'n': n, 'z': atom.atomic_number, 'h': atom.implicit_hydrogens,
                          'charge': atom.charge, 'radical': atom.is_radical,
                          'isotope': atom.isotope or 0, 'stereo': atom.stereo})
        bonds = [{'n': n, 'm': m, 'order': int(bond), 'stereo': bond.stereo}
                 for n, m, bond in mol.bonds()]
        out.append({'atoms': atoms, 'bonds': bonds, 'canonical': format(mol, ''),
                    'numbers': list(mol),
                    'allenes': sorted(mol.stereogenic_allenes),
                    'cis_trans_counterpart': {str(k): v for k, v
                                              in mol._stereo_cis_trans_counterpart.items()}})
    return out


def op_translate(payload):
    """Signs for a batch of frames the CALLER computed. A KeyError is `None`: V2 states nothing.

    One call per bridged record rather than per stereo unit -- the frames for a record are all known
    once the core has enumerated its units, and a record has one to a few units.
    """
    mol = HANDLES[payload['handle']]
    out = []
    for query in payload['queries']:
        kind = query[0]
        try:
            if kind == 0:
                out.append(mol._translate_tetrahedron_sign(query[1], tuple(query[2])))
            elif kind == 1:
                out.append(mol._translate_cis_trans_sign(query[1], query[2], query[3], query[4]))
            else:
                out.append(mol._translate_allene_sign(query[1], query[2], query[3]))
        except KeyError:
            out.append(None)
    return out


def op_parse_canonical(payload):
    """`format(smiles(text))` for each text, or None where V2 will not read it."""
    from chython import smiles
    out = []
    for text in payload['texts']:
        try:
            out.append(format(smiles(text), ''))
        except Exception:
            out.append(None)
    return out


def op_roundtrip(payload):
    """Our string back through V2, against V2's own canonical form of the source molecule.

    Answers the three questions the automorphism-invariance measurement asks together, because they
    are three questions about one parse and splitting them would parse the same string three times:
    what V2 calls our molecule, whether either string is a fixpoint of V2's own loop, and whether
    V2's OWN MATCHER -- the part of V2 that does not depend on a string -- calls them the same graph.
    """
    from chython import smiles
    out = []
    for handle, text in payload['pairs']:
        mol = HANDLES[handle]
        here = format(mol, '')
        try:
            back = smiles(text)
        except Exception as e:
            out.append({'error': repr(e)})
            continue
        there = format(back, '')
        row = {'here': here, 'there': there, 'error': None}
        if here != there:
            # the two fixpoint checks stay separate rather than being OR-ed here: the allene
            # measurement needs to say that OUR string is a fixpoint of V2's loop AND V2's own is,
            # which is a stronger statement than "at least one of them oscillates"
            row['here_fixpoint'] = format(smiles(here), '') == here
            row['there_fixpoint'] = format(smiles(there), '') == there
            row['same'] = len(mol) == len(back) and mol.get_mapping(back) is not None
        out.append(row)
    return out


def op_remap_canonical(payload):
    """V2's own string over a batch of renumberings -- V2 measured the same way we measure ourselves."""
    mol = HANDLES[payload['handle']]
    numbers = list(mol)
    out = []
    for shuffled in payload['orders']:
        other = mol.copy()
        other.remap(dict(zip(numbers, shuffled)))
        out.append(format(other, ''))
    return out


OPS = {name[3:]: value for name, value in list(globals().items()) if name.startswith('op_')}

for line in sys.stdin:
    line = line.strip()
    if not line:
        continue
    request = json.loads(line)
    try:
        answer = {'ok': OPS[request['op']](request)}
    except Exception:
        answer = {'error': traceback.format_exc()}
    sys.stdout.write(json.dumps(answer) + '\n')
    sys.stdout.flush()
'''


class AtomView:
    """One atom as chython 2 described it, under chython 2's own attribute names.

    The names are V2's on purpose: a caller written against a V2 molecule object reads this view
    unchanged, where renaming `implicit_hydrogens` to `h` would mean editing every assertion in it.
    `implicit_hydrogens is None` keeps its V2 meaning: the valence model could not derive a count.
    That is a real answer to be matched, never a gap to be filled with zero.
    """
    __slots__ = ('n', 'atomic_number', 'charge', 'is_radical', 'isotope', 'implicit_hydrogens',
                 'stereo')

    def __init__(self, row):
        self.n = row['n']
        self.atomic_number = row['z']
        self.charge = row['charge']
        self.is_radical = row['radical']
        self.isotope = row['isotope'] or None
        self.implicit_hydrogens = row['h']
        self.stereo = row['stereo']


class BondView:
    """One bond, likewise. `int(bond)` is V2's own spelling for the order and is kept working."""
    __slots__ = ('order', 'stereo')

    def __init__(self, row):
        self.order = row['order']
        self.stereo = row['stereo']

    def __int__(self):
        return self.order


class Record:
    """What chython 2 says about one molecule, plus the handle for asking it something new.

    The constitution, the sign store's own verdict and V2's canonical string are fetched once and
    cached, because they do not depend on any question.  A stereo SIGN does depend on one -- it is a
    sign IN A FRAME THE CALLER COMPUTED -- so those go over the pipe as they are needed, which is the
    whole reason the session exists.

    Two surfaces, deliberately.  `atom_rows`/`bond_rows` are the decoded documents, which is what a
    caller building its own molecule wants.  `atoms()`, `bonds()`, `atom(n)` and iteration mimic a V2
    molecule, which is what the callers that were written against one want.  Neither is a
    reimplementation of chython 2: every one of them is a field V2 filled in.
    """
    __slots__ = ('_session', 'handle', 'atom_rows', 'bond_rows', 'canonical', 'numbers', 'allenes',
                 'counterpart', '_atoms')

    def __init__(self, live, handle, data):
        self._session = live
        self.handle = handle
        self.atom_rows = data['atoms']
        self.bond_rows = data['bonds']
        self.canonical = data['canonical']
        self.numbers = data['numbers']
        self.allenes = set(data['allenes'])
        self.counterpart = {int(k): v for k, v in data['cis_trans_counterpart'].items()}
        self._atoms = {row['n']: AtomView(row) for row in self.atom_rows}

    # --- the V2-molecule surface
    def __iter__(self):
        """V2 iterates a molecule over its atom NUMBERS, in its own order. So does this."""
        return iter(self.numbers)

    def __len__(self):
        return len(self.numbers)

    def atom(self, n):
        return self._atoms[n]

    def atoms(self):
        return ((row['n'], self._atoms[row['n']]) for row in self.atom_rows)

    def bonds(self):
        return ((row['n'], row['m'], BondView(row)) for row in self.bond_rows)

    # --- the questions that need a frame
    def translate(self, queries):
        """chython 2's sign for each frame, `None` where V2's store says nothing about it."""
        if not queries:
            return []
        return self._session.call('translate', handle=self.handle, queries=queries)

    def tetrahedron_sign(self, anchor, env):
        """One tetrahedral frame. `None` when V2 cannot express it -- not an oracle for that centre."""
        sign, = self.translate([[0, anchor, list(env)]])
        return sign

    def cis_trans_sign(self, n, m, a, b):
        sign, = self.translate([[1, n, m, a, b]])
        return sign


class Session:
    """One chython 2 interpreter, alive for as long as the fixture that owns it.

    Use it through the `session` fixture rather than building one per test: the process costs 220ms
    to start and the molecules it holds are worth reusing across the tests that share a corpus.
    """
    def __init__(self, python):
        # `ISOLATION` and not a literal, for the same reason `run` uses it: the co-process is a
        # spawn of the oracle like any other and must not be the one place the guard is spelled by hand
        self._process = subprocess.Popen(
            [str(python), ISOLATION, '-c', SERVER],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
            bufsize=1)

    def call(self, op, **kwargs):
        kwargs['op'] = op
        if self._process.poll() is not None:
            raise AssertionError('the chython 2 oracle process died:\n'
                                 f'{self._process.stderr.read()}')
        self._process.stdin.write(json.dumps(kwargs) + '\n')
        self._process.stdin.flush()
        line = self._process.stdout.readline()
        if not line:
            raise AssertionError('the chython 2 oracle process stopped answering:\n'
                                 f'{self._process.stderr.read()}')
        answer = json.loads(line)
        if 'error' in answer:
            raise AssertionError(f'the chython 2 oracle failed on {op}:\n{answer["error"]}')
        return answer['ok']

    def close(self):
        if self._process.poll() is None:
            self._process.stdin.close()
            try:
                self._process.wait(timeout=30)
            except subprocess.TimeoutExpired:
                self._process.kill()

    # --- the three ways a corpus arrives, each returning `Record`s and each batched
    #
    # BATCHED, not because the pipe is slow, but because V2's PARSER is: reading 5,000 SMILES or 512
    # SDF records dominates everything else in these files, and one call that reads them all keeps
    # that cost where it was instead of adding a round trip per record on top of it.

    def records(self, handles):
        """`Record` per handle, `None` where the handle is `None`. One `record` call for the batch."""
        wanted = [h for h in handles if h is not None]
        got = dict(zip(wanted, self.call('record', handles=wanted))) if wanted else {}
        return [None if h is None else Record(self, h, got[h]) for h in handles]

    def read_smiles(self, texts):
        """`(text, Record)` for each string chython 2 accepts, the rest DROPPED.

        A record V2 cannot parse has nothing to compare, so it is not an error: raising would let one
        unreadable string cost a 5,000-record corpus.
        """
        texts = list(texts)
        handles = self.call('parse', texts=texts)
        return [(text, record) for text, record
                in zip(texts, self.records(handles)) if record is not None]

    def read_sdf(self, paths, tolerant=False):
        """`{name: [Record, ...]}` for `{name: path}`. `tolerant` maps an unreadable file to `None`."""
        answer = self.call('read_sdf', paths={k: str(v) for k, v in paths.items()},
                           tolerant=tolerant)
        flat = [h for handles in answer.values() if handles is not None for h in handles]
        made = iter(self.records(flat))
        return {name: None if handles is None else [next(made) for _ in handles]
                for name, handles in answer.items()}

    def read_rdf(self, paths, tolerant=False):
        """`{name: [{...}, ...]}` for `{name: path}`.  `tolerant` maps an unreadable file to `None`."""
        return self.call('read_rdf', paths={k: str(v) for k, v in paths.items()}, tolerant=tolerant)

    def read_mdl(self, texts):
        """`Record` per MDL record TEXT, `None` where V2 refuses it."""
        return self.records(self.call('read_mdl_text', texts=list(texts)))


def session():
    """A `Session`, or skip. Meant to back a session-scoped pytest fixture.

    `verify` runs first, so a module whose whole point is a differential does not get as far as its
    corpus before finding out that the oracle it is differing against is the tree under test.
    """
    python = require()
    verify()
    return Session(python)


def main():
    """`python -m chython.core.test.oracle` -- report whether the oracle is usable, and say why not."""
    found = interpreter()
    if found is None:
        print(f'no chython 2 oracle: {ENV_VAR} unset and nothing at {DEFAULT_PATH}')
        print(f'provision one:\n    python3.10 -m venv {DEFAULT_PATH.parent.parent}\n'
              f'    {DEFAULT_PATH.parent}/pip install chython=={VERSION}')
        return 1
    print(f'chython 2 oracle: {found}')
    version, source = probe()
    print(f'  version {version} (pinned {VERSION})')
    print(f'  imports {source}')
    try:
        verify()
    except AssertionError as e:
        print(f'UNUSABLE: {e}')
        return 1
    live = Session(found)
    try:
        handle, = live.call('parse', texts=['CCO'])
        print(f'  co-process answers: {live.call("record", handles=[handle])[0]["canonical"]}')
    finally:
        live.close()
    return 0


if __name__ == '__main__':
    sys.exit(main())
