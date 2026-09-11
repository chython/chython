#!/usr/bin/env python3
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
"""The valence rule collection: derive it, compile it, measure it.

Three verbs on one subject, which is why they are one file:

    python chython/core/test/gen_valence_rules.py compile
        chython/core/valence_rules.tsv -> the generated tables in chython/core/_valence.pxi.
        Run this after editing the TSV.  Idempotent; also rewrites the TSV itself if its rows
        are out of scan order, because the file's row order IS part of the semantics (see
        `canonical_order`).

    python chython/core/test/gen_valence_rules.py derive
        chython 2's `periodictable` -> the `common` and `curated` rows of the TSV.  This is how
        those two provenances stay checkable; rows with any other provenance are carried through
        untouched.  Run it when chython 2's tables change.

        chython 2 is reached through `oracle`, which runs an INSTALLED copy in another
        interpreter.  This file imports no chython 2 and neither does its test, which is what lets
        chython 2 be deleted from this repository without the core's suite dying with it -- see
        `test_no_chython_two_imports.py`.  Unprovisioned, `derive` and the one test that calls it
        skip; the TSV is checked in, so nothing else in the collection's coverage is lost.

    python chython/core/test/gen_valence_rules.py coverage [--v2] FILE...
        Parse molecules with the core's own reader and report how many atoms the collection has
        no rule for, by element.  This is the data-driven-development hook: an element at the top
        of that report is one whose rules are missing, and a `mined:<corpus>` row is how it gets
        fixed.  `--v2` reads with chython 2 instead (through the oracle); the two must produce the
        same report and `test_both_readers_measure_the_same_coverage` requires it, which is what
        keeps switching the parser from being a silent corpus swap.

WHY A TSV AND NOT PYTHON.

The rules are data.  Spelled as 118 Python classes with two properties each, a rule change is a code
change, the collection sits behind an import, and no tool that is not Python can read it.  One
tab-separated file can be diffed, sorted, grepped, joined against a mining run's output, and
reviewed by a chemist who does not read Cython.  `test_valence.py` re-derives every `common` and
`curated` row from chython 2 and fails if a single one drifted.

Reading the file needs no chython 2 either: the element names come from `elements.tsv`, so a data
file's parser does not depend on a whole library being importable, and `derive` requires the two
libraries to agree on all 118 before it writes a row keyed by one.

THE FLAT SHAPE, AND WHAT IT COSTS.

chython 2 states a rule as `(charge, radical, implicit, environment)` where `implicit` is a
*maximum*: an alcohol oxygen's `(0, False, 1, ((1, 'C'),))` means "one single-bonded carbon, and
then either one hydrogen or one more bond".  One authored rule therefore describes several
valence states.  The TSV has one row per *state* instead, so that rule is two rows -- bonds 1
with 1 hydrogen, bonds 2 with 0.  The cost is that editing "the alcohol rule" means editing two
adjacent rows.  The gain is that a row is exactly one observation: element, charge, radical,
bond-order sum, hydrogen count, neighbourhood.  That is the shape a corpus miner produces and
the shape the coverage report counts, and "implicit means up to" is a semantic nobody has to
learn twice.

THE COLUMNS.

    element     symbol, as the core spells it
    charge      -4..4
    radical     0 or 1
    bonds       sum of the bond orders to explicit neighbours, hydrogens included
    implicit_h  the hydrogen count this state carries.  Exact, not a maximum
    env         `*` for "any neighbourhood", else space-separated `<glyph><symbol>` tokens with
                `-` `=` `#` for orders 1, 2, 3.  A LOWER BOUND, never an exact match: `-C`
                matches methanol and dimethyl ether alike, and a repeated token means a
                multiplicity, so `=O =O` does not match a sulfoxide
    provenance  `common`   from chython 2's `_common_valences`
                `curated`  from chython 2's `_valences_exceptions`
                `mined:<corpus>` added from data.  Say which corpus; that is the whole point
"""
from __future__ import annotations

import pathlib
import sys
from collections import Counter, defaultdict


# The scan order's field order and the key packing are DEFINED in _valence.pxi (VAL_KEY_*).
# They are repeated here because this script emits the packed keys; nothing checks the two
# copies by inspection, and nothing needs to -- if they disagreed, no lookup in the exhaustive
# sweep in test_valence.py would find its own row.
CHARGE_BIAS = 8
KEY_Z_SHIFT = 16
KEY_CHARGE_SHIFT = 12
KEY_RADICAL_SHIFT = 11

ORDER_GLYPH = {1: '-', 2: '=', 3: '#'}
GLYPH_ORDER = {v: k for k, v in ORDER_GLYPH.items()}

# The dense hydrogen table's extents, mirrored in _valence.pxi as VAL_H_*.  These are the
# collection's measured extents, not a guess, and `assert_hot_invariants` fails the compile if a
# new row leaves them.
H_Z_MAX = 118
H_CHARGE_MIN = -4
H_CHARGE_MAX = 4
H_BONDS_MAX = 8
H_NO_RULE = -1
H_CONSULT = -2

ROOT = pathlib.Path(__file__).resolve().parent.parent

# run as a script, sys.path[0] is this directory and the repo is not on it at all, so `chython`
# would resolve to whatever is installed -- a different checkout's rules, silently.  The
# derivation and the tables it feeds must come from the tree this file is in
if str(ROOT.parent.parent) not in sys.path:
    sys.path.insert(0, str(ROOT.parent.parent))

TSV = ROOT / 'valence_rules.tsv'
PXI = ROOT / '_valence.pxi'
BEGIN = '# --- BEGIN GENERATED TABLES: python chython/core/test/gen_valence_rules.py compile ---'
END = '# --- END GENERATED TABLES ---'

HEADER = ('element', 'charge', 'radical', 'bonds', 'implicit_h', 'env', 'provenance')

PREAMBLE = """\
# The valence rule collection.  THIS FILE IS THE AUTHORITY; the C tables in _valence.pxi are
# compiled from it by `python chython/core/test/gen_valence_rules.py compile`.
#
# It answers a question about a MOLECULE -- is this a valence state chemistry is known to allow,
# and how many hydrogens does it come with.  It is a collection of states that have been
# observed, not a theory: it is neither complete nor ideal, it exists to catch bad input, and it
# improves by adding rows with evidence behind them.  Nothing in chython may REJECT a structure
# because of a verdict from this file.
#
# It is NOT the model that decides how many hydrogens a bracketless SMILES atom implies.  That is
# a question about a notation, its authority is the OpenSMILES specification, and it lives in the
# SMILES layer.  The two disagree on purpose -- a bare `S` with six single-bonded carbons is
# legal SMILES implying zero hydrogens, and has no row here.  Merging them breaks both
# directions; a test in each suite fails if anyone tries.
#
# ROW ORDER IS SEMANTICS, NOT COSMETICS.  Several rows can describe one `(element, charge,
# radical, bonds)` key, and "how many hydrogens" answers with the FIRST of them that matches the
# neighbourhood.  Rows are therefore stored in scan order: by key, and within a key `env=*` rows
# before the rest.  `compile` reorders the file if it is not, and a test fails if a commit lands
# it out of order.  Fifteen keys in the shipped data give different answers under a different
# order, so this is load-bearing -- see test_row_order_within_a_key_is_observable.
#
# WHEN YOU ADD A ROW: an `env=*` row is consulted before every row with an environment at the
# same key, so it will shadow them for the hydrogen count.  If you mean "only when X is present",
# give it an environment.  "Is this legal" is unaffected -- that question accepts any matching
# row, which is why the two are not each other's inverse.
#
# That ordering is also what lets the hydrogen count be precomputed into a dense table instead of
# searched for, so `compile` refuses a file that breaks it, and refuses a row whose bonds, charge
# or element leaves the extents the dense table covers (bonds 0..8, charge -4..4).  Widening those
# is a one-line change in `gen_valence_rules.py`; it is guarded because it is a size change nobody
# would otherwise notice.
#
# ALUMINIUM'S ROWS ARE NOT OUT OF ORDER -- DO NOT "FIX" THEM.  Nine elements carry both a hydride
# ladder and a bare-atom row at one key.  Eight put the hydride first, so `[C]` is methane; Al puts
# the bare atom first, which makes its ladder read as self-contradictory (nothing on it: 0 H, one
# carbon on it: 2 H).  The difference is that Al is a METAL: a lone Al in a connectivity file is the
# metal or its ion, while a substituted Al is an organoaluminium, where the hydride is ordinary and
# load-bearing -- DIBAL-H is a two-coordinate aluminium with one.  A notation that means alane can
# say so and is believed (`[AlH3]`); a file that states nothing is stating a metal.  Both readings
# stay legal either way -- only the derived default differs, which is why this is an ordering
# question at all.  Pinned by test_a_LONE_metal_is_the_metal_and_only_a_SUBSTITUTED_one_takes_hydrides.
#
# The same convention is why As Sb Ga In Tl Sn Pb Bi Po At have NO neutral hydride ladder: their
# bare rows are metals, and a partially substituted one returns *unknown* rather than a guess, which
# is reported and honest.  Filling those in needs a per-element valence default, and for exactly
# these elements there are two -- R2Sn is a stannylene or R2SnH2, R-Tl is Tl(I) or Tl(III) -- so a
# first-matching row there invents a compound.  Evidence first; this is not a gap to close by hand.
#
# Do not hand-edit `common` or `curated` rows: they are derived from chython 2 and a test
# re-derives them.  Re-run `derive`, or add a row with a new provenance.
"""


class Rule:
    """One row: a valence state and the hydrogen count it carries."""
    __slots__ = ('z', 'charge', 'radical', 'bonds', 'h', 'env', 'provenance')

    def __init__(self, z, charge, radical, bonds, h, env, provenance):
        self.z = z
        self.charge = charge
        self.radical = radical
        self.bonds = bonds
        self.h = h
        self.env = env                  # sorted tuple of (order, atomic number), possibly empty
        self.provenance = provenance

    @property
    def key(self):
        return ((self.z << KEY_Z_SHIFT) | ((self.charge + CHARGE_BIAS) << KEY_CHARGE_SHIFT) |
                (self.radical << KEY_RADICAL_SHIFT) | self.bonds)

    def __eq__(self, other):
        return self.row() == other.row()

    def __repr__(self):
        return f'Rule({"|".join(self.row())})'

    def row(self):
        return (symbol(self.z), str(self.charge), '1' if self.radical else '0', str(self.bonds),
                str(self.h), format_env(self.env), self.provenance)


def symbols():
    """Atomic number -> symbol, from the core's own table.

    `elements.tsv` is the authority for the names, so reading the TSV needs no chython 2 import --
    a dependency nobody would look for in a data file's parser.  `test_element_tables.py` checks
    that the two libraries spell all 118 the same way.
    """
    from chython.core._core import element_symbols

    table = element_symbols()
    return {z: table[z] for z in range(1, len(table))}


_SYMBOLS = None
_NUMBERS = None


def symbol(z):
    global _SYMBOLS
    if _SYMBOLS is None:
        _SYMBOLS = symbols()
    return _SYMBOLS[z]


def number(name):
    global _NUMBERS
    if _NUMBERS is None:
        _NUMBERS = {s: z for z, s in symbols().items()}
    if name not in _NUMBERS:
        raise ValueError(f'unknown element symbol {name!r}')
    return _NUMBERS[name]


def format_env(env):
    if not env:
        return '*'
    return ' '.join(f'{ORDER_GLYPH[order]}{symbol(z)}' for order, z in env)


def parse_env(text):
    if text == '*':
        return ()
    out = []
    for token in text.split():
        if token[0] not in GLYPH_ORDER:
            raise ValueError(f'environment token {token!r} must start with one of -=#')
        out.append((GLYPH_ORDER[token[0]], number(token[1:])))
    return tuple(sorted(out))


def canonical_order(rules):
    """Scan order: by key, and within a key `env=*` first, otherwise as authored.

    Both halves are measured rather than assumed.  Sorting by key is free -- lookup is by exact
    key.  The `env=*` partition is a *stable* one, and on the shipped collection it is the
    identity permutation: chython 2 already emits every empty-environment rule at a key before
    every non-empty one, so this reordering changes no answer.  That was checked over all 576
    keys before it was written down, and `test_the_shipped_file_is_already_in_scan_order` keeps
    it true.  A stable partition also matters for its own sake: fourteen of the fifteen
    order-observable keys have `env=*` on both sides of the pair, so an unstable sort within the
    `*` group would silently change the hydrogen count of a bare atom.
    """
    order = {}
    for i, rule in enumerate(rules):
        order[id(rule)] = i
    return sorted(rules, key=lambda r: (r.key, 0 if not r.env else 1, order[id(r)]))


def read_tsv(path=TSV):
    rules = []
    for lineno, line in enumerate(path.read_text(encoding='utf-8').splitlines(), 1):
        line = line.rstrip('\n')
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        fields = line.split('\t')
        if tuple(fields) == HEADER:
            continue
        if len(fields) != len(HEADER):
            raise ValueError(f'{path}:{lineno}: {len(fields)} fields, expected {len(HEADER)}')
        element, charge, radical, bonds, h, env, provenance = fields
        if radical not in ('0', '1'):
            raise ValueError(f'{path}:{lineno}: radical must be 0 or 1, got {radical!r}')
        rules.append(Rule(number(element), int(charge), radical == '1', int(bonds), int(h),
                          parse_env(env), provenance))
    return rules


def write_tsv(rules, path=TSV):
    lines = [PREAMBLE, '\t'.join(HEADER)]
    for rule in rules:
        lines.append('\t'.join(rule.row()))
    path.write_text('\n'.join(lines) + '\n')


# chython 2's two properties per element, expanded into one row per valence state.  This runs IN
# THE ORACLE INTERPRETER, not here -- see `oracle` for why the second opinion is a subprocess.
#
# It reproduces `_compiled_valence_rules` -- including the *order* in which it appends, which is
# observable through the first-match rule -- and then flattens its dict into rows.  The arithmetic
# is chython 2's, quoted rather than reasoned about: only the first common valence grants implicit
# hydrogens, hydrogen itself is excluded so `[H][H]` never grows a third hydrogen, and an element
# whose first common valence is 0 takes the same path as hydrogen.
DERIVE_IN_ORACLE = """
from collections import defaultdict

import chython.periodictable                  # noqa: F401  -- registers the subclasses
from chython.periodictable.base.element import Element

lifted = {}
for cls in Element.__subclasses__():
    lifted[cls.__name__] = cls.atomic_number.fget(None)

rules = defaultdict(list)           # key tuple -> list of (h, env, provenance)
for cls in Element.__subclasses__():
    z = lifted[cls.__name__]
    common = cls._common_valences.fget(None)
    if common and common[0] and z != 1:
        valence = common[0]
        for h in range(valence + 1):
            rules[(z, 0, False, valence - h)].append((h, (), 'common'))
        for valence in common[1:]:
            rules[(z, 0, False, valence)].append((0, (), 'common'))
    else:
        for valence in common:
            rules[(z, 0, False, valence)].append((0, (), 'common'))

    for charge, radical, implicit, environment in cls._valences_exceptions.fget(None):
        env = tuple(sorted((order, lifted[e]) for order, e in environment))
        explicit = sum(order for order, _ in environment)
        if implicit:
            for h in range(implicit + 1):
                rules[(z, charge, radical, explicit + implicit - h)].append((h, env, 'curated'))
        else:
            rules[(z, charge, radical, explicit)].append((0, env, 'curated'))

out = []
for (z, charge, radical, bonds), entries in rules.items():
    for h, env, provenance in entries:
        out.append((z, charge, radical, bonds, h, sorted(env), provenance))
# the symbol map travels with the rows so that the caller can require the two libraries spell all
# 118 elements the same way before it writes any of them into a file keyed by symbol
_emit({'rules': out, 'symbols': {str(v): k for k, v in lifted.items()}})
"""


def derive():
    """chython 2's rules, flattened into one row per valence state. Needs the oracle.

    Skips when the oracle is not provisioned -- the TSV is checked in and every other test in
    `test_valence.py` reads it, so an unprovisioned machine loses this one comparison and keeps
    the rest of the collection's coverage.
    """
    from chython.core.test.oracle import ask

    answer = ask(DERIVE_IN_ORACLE)
    # V2 and the core must agree on all 118 names or the TSV would round-trip through a different
    # element.  Checked here rather than in a test of its own: this is the one place where a
    # disagreement would corrupt data instead of merely failing a comparison.
    theirs = {int(z): name for z, name in answer['symbols'].items()}
    ours = symbols()
    assert theirs == ours, {z: (theirs.get(z), ours.get(z)) for z in set(theirs) | set(ours)
                            if theirs.get(z) != ours.get(z)}
    return [Rule(z, charge, radical, bonds, h, tuple(sorted(tuple(t) for t in env)), provenance)
            for z, charge, radical, bonds, h, env, provenance in answer['rules']]


def emit_array(name, ctype, values, per_line, comment):
    lines = [f'    /* {comment} */', f'    static const {ctype} {name}[{len(values)}] = {{']
    for i in range(0, len(values), per_line):
        chunk = ', '.join(str(v) for v in values[i:i + per_line])
        lines.append(f'    {chunk},' if i + per_line < len(values) else f'    {chunk}')
    lines.append('    };')
    return lines


def hydrogen_tables(rules):
    """The HOT artifact: a dense table that answers "how many hydrogens" without a search.

    Two questions share one collection and only one of them runs on every atom of every parsed
    molecule.  Of the 1036 rows, 1031 answer the hydrogen question with no reference to a
    neighbourhood at all, and the hydrogen question never needs to look at the other five unless the
    atom is a phosphorus or a tin.  Making that path a binary search over 576 keys plus a multiset
    compare is paying the checking question's price on the parsing question's traffic.

    So the hydrogen answer is precomputed for every state the collection covers:

        VAL_H_PAT[(z * 9 + charge + 4) * 2 + radical]  ->  a pattern index
        VAL_H_ROW[pattern * 9 + bonds]                 ->  the hydrogen count, or a sentinel

    One indexed load, one branch.  It is a projection and not a summary: every value in it is the
    answer `val_scan` gives for the same state with an empty environment, and the differential test
    sweeps all of them.  65 distinct patterns behind 2142 slots, so the whole hot artifact is about
    2.7 KB and a molecule touches a handful of cache lines of it.

    THE SENTINELS ARE THREE, NOT TWO.  `-1` is "no row covers this state" -- which is not zero
    hydrogens, and callers depend on telling those apart.  `-2` is "an environment decides here",
    and it is why this table cannot be a summary of the collection: at 99 keys the only rows are
    rows with an environment, and skipping them would answer "no rule" for a sulfone, a nitro
    group and a perchlorate.  On the public NCI 5K, 81,518 of 82,157 non-aromatic atoms are
    answered by the dense load, 484 fall through to the environment scan (466 of them sulfur), and
    155 are "no rule" without touching the table at all.

    WHY THIS IS SOUND, and the one invariant it rests on: the shipped TSV is sorted with `env=*`
    rows before environment rows at the same key, so wherever an `env=*` row exists it is the first
    match and the environment cannot change the answer.  That sort was measured to be the identity
    permutation on this collection -- it reorders nothing -- which is what makes precomputing legal
    rather than a behaviour change.  `assert_hot_invariants` re-checks it on every compile.
    """
    assert_hot_invariants(rules)

    # first env-free row per (key, bonds) -- the same "first match wins" the scanner implements
    star = defaultdict(dict)
    consult = defaultdict(set)
    for rule in rules:
        if rule.env:
            consult[(rule.z, rule.charge, rule.radical)].add(rule.bonds)
        else:
            star[(rule.z, rule.charge, rule.radical)].setdefault(rule.bonds, rule.h)

    patterns = {}
    index = []
    for z in range(H_Z_MAX + 1):
        for charge in range(H_CHARGE_MIN, H_CHARGE_MAX + 1):
            for radical in (False, True):
                key = (z, charge, radical)
                found = star.get(key, {})
                envs = consult.get(key, ())
                pattern = tuple(found.get(bonds, H_CONSULT if bonds in envs else H_NO_RULE)
                                for bonds in range(H_BONDS_MAX + 1))
                if pattern not in patterns:
                    patterns[pattern] = len(patterns)
                index.append(patterns[pattern])
    assert len(patterns) < 256, f'{len(patterns)} patterns will not fit an unsigned char index'

    flat = []
    for pattern in patterns:
        flat.extend(pattern)
    return index, flat, len(patterns)


def assert_hot_invariants(rules):
    """The three facts the dense table is built on.  A failure here means the TSV outgrew it."""
    for rule in rules:
        assert rule.bonds <= H_BONDS_MAX, f'{rule!r} exceeds H_BONDS_MAX = {H_BONDS_MAX}'
        assert H_CHARGE_MIN <= rule.charge <= H_CHARGE_MAX, f'{rule!r} is outside the charge span'
        assert rule.z <= H_Z_MAX, f'{rule!r} is outside the element span'
    seen_env = set()
    for rule in rules:
        if rule.env:
            seen_env.add(rule.key)
        else:
            assert rule.key not in seen_env, (
                f'{rule!r} follows an environment row at the same key.  The dense hydrogen table '
                f'assumes `env=*` rows are scanned first -- run `compile` to restore scan order')


def compile_tables(rules):
    """The generated block: a sorted key index, the rules behind each key, interned environments.

    Environments are interned by content -- 219 distinct neighbourhoods behind 1036 rules -- so
    the flat environment array is 732 entries rather than 1824.  That is not only size: an
    interned block means two rules that demand the same neighbourhood compare their requirement
    against the same bytes, and the table cannot drift between them.

    This is the COLD artifact.  It answers "is this state described, and does anything accept it",
    it carries every row including the 583 that exist only to be checked against, and nothing on
    the parse path reaches it except the 0.6% of atoms whose hydrogen count an environment decides.
    """
    keys = []
    key_off = []
    key_len = []
    rule_h = []
    rule_env_off = []
    rule_env_len = []
    env_flat = []
    interned = {}

    for rule in rules:
        if rule.env not in interned:
            interned[rule.env] = len(env_flat)
            for order, z in rule.env:
                env_flat.append((order << 8) | z)

    for rule in rules:
        if not keys or keys[-1] != rule.key:
            keys.append(rule.key)
            key_off.append(len(rule_h))
            key_len.append(0)
        key_len[-1] += 1
        rule_h.append(rule.h)
        rule_env_off.append(interned[rule.env])
        rule_env_len.append(len(rule.env))

    h_index, h_flat, h_patterns = hydrogen_tables(rules)

    provenance = Counter(r.provenance for r in rules)
    elements = len({r.z for r in rules})
    lines = [
        BEGIN,
        '# Compiled from chython/core/valence_rules.tsv, which is the authority.  Do not edit by',
        '# hand -- run the command in the marker above.  The TSV is in scan order and so is this,',
        '# so the k-th entry here is the k-th row there.',
        '#',
        f'# {len(rules)} rules over {len(keys)} keys on {elements} elements; '
        f'{len(interned)} distinct environments, {len(env_flat)} entries.',
        '# Provenance: ' + ', '.join(f'{n} {p}' for p, n in sorted(provenance.items())) + '.',
        '#',
        f'# The hot artifact is separate and is a projection of the same rows: {h_patterns} '
        f'distinct hydrogen',
        f'# patterns behind {len(h_index)} states.  See `hydrogen_tables`.',
        'cdef extern from *:',
        '    """',
    ]
    lines += emit_array('VAL_KEY', 'unsigned int', keys, 8,
                        'sorted: (z << 16) | ((charge + 8) << 12) | (radical << 11) | bonds')
    lines += emit_array('VAL_KEY_OFF', 'unsigned short', key_off, 12,
                        'first rule of the k-th key')
    lines += emit_array('VAL_KEY_LEN', 'unsigned char', key_len, 16,
                        'how many rules that key has, in scan order')
    lines += emit_array('VAL_H', 'unsigned char', rule_h, 24,
                        'implicit hydrogen count of the k-th rule')
    lines += emit_array('VAL_ENV_OFF', 'unsigned short', rule_env_off, 12,
                        'the k-th rule\'s environment, interned by content')
    lines += emit_array('VAL_ENV_LEN', 'unsigned char', rule_env_len, 24,
                        'how many neighbours it demands; 0 is `env=*`')
    lines += emit_array('VAL_ENV', 'unsigned short', env_flat, 12,
                        '(bond order << 8) | atomic number')
    lines += emit_array('VAL_H_PAT', 'unsigned char', h_index, 16,
                        'hot: pattern of ((z * 9 + charge + 4) * 2 + radical)')
    lines += emit_array('VAL_H_ROW', 'signed char', h_flat, 9,
                        'hot: [pattern][bonds] -> hydrogens, -1 no rule, -2 an environment decides')
    lines += [
        '    """',
        f'    const uint32_t VAL_KEY[{len(keys)}]',
        f'    const uint16_t VAL_KEY_OFF[{len(key_off)}]',
        f'    const uint8_t VAL_KEY_LEN[{len(key_len)}]',
        f'    const uint8_t VAL_H[{len(rule_h)}]',
        f'    const uint16_t VAL_ENV_OFF[{len(rule_env_off)}]',
        f'    const uint8_t VAL_ENV_LEN[{len(rule_env_len)}]',
        f'    const uint16_t VAL_ENV[{len(env_flat)}]',
        f'    const uint8_t VAL_H_PAT[{len(h_index)}]',
        f'    const int8_t VAL_H_ROW[{len(h_flat)}]',
        '',
        f'DEF VAL_KEY_COUNT = {len(keys)}',
        f'DEF VAL_H_PATTERNS = {h_patterns}',
        END,
    ]
    return '\n'.join(lines)


def rewrite_pxi(block, path=PXI):
    text = path.read_text(encoding='utf-8')
    start = text.index(BEGIN)
    stop = text.index(END) + len(END)
    if text[start:stop] == block:
        return False
    path.write_text(text[:start] + block + text[stop:])
    return True


# The two readers, behind one shape: (z, charge, radical, order_sum, implicit_h, env, aromatic) per
# atom, for a whole list of strings at a time.  Two of them and not one because a parser swap has to
# be EVIDENCED rather than announced -- a silent one changes the corpus underneath a coverage number
# that people quote, which is the more insidious of the two failure modes a shadowed name creates.
# `test_both_readers_measure_the_same_coverage` runs the pair over the same strings and requires the
# same counts.
#
# The V2 reader runs in the oracle interpreter, so this file imports no chython 2 -- see
# `oracle`.  A whole corpus per subprocess and not a string per subprocess: the round trip costs
# more than the parse, and at one call per line a 5K corpus would take twenty minutes.
#
# `implicit_h` is the state being CHECKED, so it stays out of the bond-order sum and out of the
# environment.  An explicit H atom is an ordinary neighbour and is in both, which is what chython 2
# does and what the rows written against a hydrogen neighbour need.

def _atom_states_core(lines):
    from chython.core._core import read_smiles

    out = []
    for line in lines:
        try:
            mol = read_smiles(line)
        except Exception:
            out.append(None)
            continue
        states = []
        for atom in mol.atoms():
            env = []
            order_sum = 0
            aromatic = False
            for other in mol.neighbors_of(atom.n):
                order = mol.bond(atom.n, other).order
                if order == 4:
                    aromatic = True
                    break
                if order == 8:
                    continue        # a dative contact carries no electron pair; see _valence.pxi
                order_sum += order
                env.append((order, mol.atom(other).element))
            states.append((atom.element, atom.charge, atom.is_radical, order_sum, atom.implicit_h,
                           env, aromatic))
        out.append(states)
    return out


# by module path, never `from chython import smiles`: the root name became the core's reader, and a
# root import here would quietly make this the other one
STATES_IN_ORACLE = """
from chython.files.daylight.smiles import smiles

out = []
for line in _payload:
    try:
        mol = smiles(line)
    except Exception:
        out.append(None)
        continue
    atoms = mol._atoms
    states = []
    for n, atom in atoms.items():
        env = []
        order_sum = 0
        aromatic = False
        for m, bond in mol._bonds[n].items():
            order = int(bond)
            if order == 4:
                aromatic = True
                break
            if order == 8:
                continue
            order_sum += order
            env.append((order, atoms[m].atomic_number))
        states.append((atom.atomic_number, atom.charge, atom.is_radical, order_sum,
                       atom.implicit_hydrogens or 0, env, aromatic))
    out.append(states)
_emit(out)
"""


def _atom_states_v2(lines):
    from chython.core.test.oracle import ask

    out = []
    for states in ask(STATES_IN_ORACLE, list(lines)):
        if states is None:
            out.append(None)
        else:
            out.append([(z, charge, radical, order_sum, implicit_h,
                         [tuple(t) for t in env], aromatic)
                        for z, charge, radical, order_sum, implicit_h, env, aromatic in states])
    return out


READERS = {'core': _atom_states_core, 'v2': _atom_states_v2}


def coverage(paths, reader='core'):
    """How many atoms the collection has nothing to say about, by element.

    Reads with the core's own reader; the verdicts come from the compiled core, so this measures the
    shipped tables and not the Python they came from.  `reader='v2'` reads with chython 2 instead,
    which is how the switch stays checkable rather than trusted -- see the note above `READERS`.

    Aromatic atoms are counted separately rather than guessed at: there are no aromatic rows, by
    design, so an aromatic atom is not a coverage gap.
    """
    from chython.core._core import valence_check

    states_of = READERS[reader]
    seen = Counter()
    unknown = Counter()
    violation = Counter()
    aromatic = Counter()
    molecules = 0
    lines = []
    for name in paths:
        for line in pathlib.Path(name).read_text(encoding='utf-8').splitlines():
            line = line.split()[0] if line.split() else ''
            if line:
                lines.append(line)
    for states in states_of(lines):
        if states is None:              # a string this reader does not accept is not a coverage gap
            continue
        molecules += 1
        for z, charge, radical, order_sum, implicit_h, env, is_aromatic in states:
            if is_aromatic:
                aromatic[symbol(z)] += 1
                continue
            seen[symbol(z)] += 1
            verdict = valence_check(z, charge, radical, order_sum, implicit_h, env)
            if verdict == 'violation':
                violation[symbol(z)] += 1
            elif verdict == 'unknown':
                unknown[symbol(z)] += 1

    total = sum(seen.values())
    print(f'{molecules} molecules, {total} non-aromatic atoms, '
          f'{sum(aromatic.values())} aromatic atoms skipped')
    print(f'{sum(unknown.values())} no rule, {sum(violation.values())} violate a known rule')
    print(f'{"element":>8} {"atoms":>8} {"no rule":>8} {"violation":>10}')
    for element, n in seen.most_common():
        if unknown[element] or violation[element]:
            print(f'{element:>8} {n:>8} {unknown[element]:>8} {violation[element]:>10}')
    return unknown, violation


def main(argv):
    verb = argv[0] if argv else 'compile'
    if verb == 'derive':
        foreign = []
        if TSV.exists():
            for rule in read_tsv():
                if rule.provenance not in ('common', 'curated'):
                    foreign.append(rule)
        rules = canonical_order(derive() + foreign)
        write_tsv(rules)
        print(f'{TSV}: {len(rules)} rows'
              + (f', {len(foreign)} carried through' if foreign else ''))
    elif verb == 'compile':
        rules = read_tsv()
        ordered = canonical_order(rules)
        if [r.row() for r in ordered] != [r.row() for r in rules]:
            write_tsv(ordered)
            print(f'{TSV}: reordered into scan order')
        if rewrite_pxi(compile_tables(ordered)):
            print(f'{PXI}: {len(ordered)} rules written')
        else:
            print(f'{PXI}: already current')
    elif verb == 'coverage':
        reader = 'core'
        argv = list(argv)
        for name in ('--v2', '--core'):
            if name in argv:
                argv.remove(name)
                reader = name[2:]
        if len(argv) < 2:
            raise SystemExit('coverage needs at least one file of SMILES')
        coverage(argv[1:], reader)
    else:
        raise SystemExit(__doc__)


if __name__ == '__main__':
    # nothing to import first: `symbols` reads the core's table and `derive` registers V2's 118
    # Element subclasses inside the oracle interpreter, where forgetting to would make the
    # derivation silently trivial rather than merely wrong
    main(sys.argv[1:])
