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
"""Derive `tables/standardize_{groups,metals}.tsv` from chython 2, and check they have not drifted.

    gen_standardize_rules.py derive   V2 `_groups.py` / `_metal_organics.py` -> the TSVs, idempotent
    gen_standardize_rules.py check    re-derive into memory and diff; non-zero on disagreement

V2 is deleted, so `derive` cannot run and the TSVs are the authority; the columns are documented in
each TSV's own header comment, which `_preamble` reads back rather than duplicating here.
"""
from __future__ import annotations

import ast
import difflib
import pathlib
import re
import sys


ROOT = pathlib.Path(__file__).resolve().parent.parent          # chython/chemistry
TABLES = ROOT / 'tables'
V2 = ROOT.parent / 'algorithms' / 'standardize'

SOURCES = (
    # tag,       V2 source,                  emitted TSV
    ('groups', V2 / '_groups.py', TABLES / 'standardize_groups.tsv'),
    ('metals', V2 / '_metal_organics.py', TABLES / 'standardize_metals.tsv'),
)

# `after`, `examples` and `why` are hand-measured and NOT derivable: `derive` carries them forward from
# the checked-in TSV by id.  A patch slot is an atom NUMBER, not a position -- V2 numbers a query's
# atoms by taking an explicit `:N` where written and otherwise the next unused integer from 1, so the
# first atom of a pattern that maps three later atoms is 4.  `atom_fix` entry order is significant
# (`metals:01` needs the metal fixed first or the charge arithmetic is wrong).
HEADER = ('id', 'smarts', 'atom_fix', 'bonds_fix', 'tautomer', 'after', 'examples', 'why')
EMPTY = '-'

RADICAL_GLYPH = {None: '-', False: '0', True: '1'}
GLYPH_RADICAL = {'-': None, '0': False, '1': True}

# V2's `z3` and the core's `z3` are different primitives sharing a spelling: V2 saturates and caps at 3,
# the core reports what it found (3 is sp only, 5 two cumulated doubles, 6 anything else).  Copying a V2
# `z3` through narrows the rule to alkynes and nitriles and its repair silently stops happening.  So
# every `z3` needs a ruling: 'A' keeps `z3`, 'C' -> `z5`, 'D' -> `z6`, where arguable the mapping widens.
# Keyed by SMARTS so V2 line drift cannot mis-assign one; the line number only cross-checks
# docs/superpowers/research/2026-09-03-z3-port-mapping.md.
TARGET = {'A': 'z3', 'C': 'z5', 'D': 'z6'}

Z3_MAP = {
    # Group A -- genuine sp, `z3` stays `z3` (14 primitives on 12 lines)
    '[N;D2;z3;+](#[N;D1])[C,N,O;z1;-]': ('A', 290),
    '[N;D2;z3;x2]([N;D2;z1])#[N;D1]': ('A', 314),
    '[N;D2;z3;x1]([N,O,S;D1])#[C;D1,D2]': ('A', 354),
    '[N;D2;z3;x1]([N;D2;z1])#[C;D1,D2]': ('A', 362),
    '[N;D2;z3;x1;+]([N,O,S;D1])#[C;D1;-]': ('A', 370),
    '[N;D2;z3;x1;+]([N;D2;z1])#[C;D1;-]': ('A', 378),
    '[N;D2;z3]([A])#[C;D1]': ('A', 386),
    '[N;D3;z3;x1](#[N;D1])(C)C': ('A', 442),
    '[N;D1;x0;z3]#[C;D2;z3;x2][O;D1]': ('A', 466),            # two primitives
    '[N;D1;x0;z3]#[C;D2;z3;x2][O;D1;-]': ('A', 474),          # two primitives
    '[C;D2;z3;x1]([N,O,S;D1])#[C;D1,D2]': ('A', 864),
    '[C;D2;z3;x1]([N;D2;z1])#[C;D1,D2]': ('A', 872),

    # Group C -- two cumulated doubles and no triple, `z3` -> `z5` (15).  These repair mis-drawn
    # pentavalent nitro / sulfonyl / phosphonate spellings, so they must keep matching them
    '[N;D3;z3;x2](=[O;D1])([O;D1])=C': ('C', 172),
    '[N;D3;z3](=[O;D1])(=[C,N,O])-[A]': ('C', 184),
    '[N;D3;z3](=[N;D3;z2;+])(=[O;D1])[A]': ('C', 207),
    '[N;D3;z3](=[N;D3;z2;+])(=[N;D1,D2;z2])[A]': ('C', 217),
    '[N;D3;z3](=[N;D1,D2;z2])(=[C,N])[A]': ('C', 229),
    '[N;D3;z3](=[O;D1])(=[O;D1])[A;-]': ('C', 251),
    '[N;D2;z3;x2;-](=[O;D1])=[O;D1]': ('C', 274),
    '[N;D2;z3;x2](=[N;D2;z2])=[N;D1]': ('C', 306),
    '[N;D2;z3;x1;+](=[N;D1])=[C;D1,D2;z2;-]': ('C', 395),
    '[P;D4;z3;-](=[O;D1])(=[O;D1])([A])[A]': ('C', 667),
    '[S;D1;-][S;D4;z3](=[O;D1])(=[A])[A]': ('C', 703),
    '[S;D1][S;D4;z3](=[O;D1])(=[A])[A]': ('C', 715),
    '[S;D3;z3;-](=[O;D1])(=[O;D1])[A]': ('C', 785),
    '[S;D3;z3;x3;-]([S;D1;-])(=[O;D1])=[O;D1]': ('C', 809),
    '[S;D4;z3:1]([O;D1:2])(=[N;D1,D2;z2:3])(=[A])[A]': ('C', 855),

    # Group D -- `z3` -> `z6`, the core's catch-all.  Six are double-plus-triple; 797 is three
    # `=O` on a D4 sulfur, which `D4` plus the explicit orders already pins (7)
    '[N;D2;z3](#[N;D1])=[C,N,O]': ('D', 282),
    '[N;D2;z3;x2](#[N;D2;+][A])=[N;D1;-]': ('D', 298),
    '[N;D2;z3;x2](=[N;D2;z2])#[N;D1;-]': ('D', 322),
    '[N;D2;z3;x2](=[N;D1;-])#[N;D1]': ('D', 330),
    '[N;D2;z3;x1](=[N;D1])#[C;D1,D2]': ('D', 338),
    '[N;D2;z3;x1](=[N,O;z2])#[C;D1,D2]': ('D', 346),
    '[S;D4;z3;-](=[O;D1])(=[O;D1])(=[O;D1])[A]': ('D', 797),
}

# Empty on purpose.  It held three "spell the four single bonds out" rewrites working around a core `D`
# that counted dative bonds; the core's `D`, `x` and `z` now all skip order 8, so re-applying them would
# put an adjacency walk back into three of the cheapest rules in the table.
DEGREE_MAP = {}

# 36 primitives on 34 lines in `_groups.py`.  The mapping document's "50 on 47" also counts
# `algorithms/groups/_functional.py`, which is a different port and not this file's input.
Z3_PRIMITIVES = 36
Z3_LINES = 34

# What was merged into what: `{surviving id: ((absorbed id, the pattern it shipped), ...)}`.  The
# patterns are the text the tables actually shipped at 413c735, so the union proof in
# `test_standardize_rules_merges.py` compares against what a caller really had.
#
# A merged row sits at the position of the LAST row it replaced, not the first: a union is at least as
# general as its members and row order is semantics, so placed early it beats a more specific rule that
# was supposed to see the site first (`groups:00`+`groups:11` at position 0 turns `BN(C)(C)C` into
# `B[N+](C)(C)C` instead of the amine-borane adduct).
#
# `metals:11` is a widening and not a union -- see `WIDENED`.
MERGES = {
    'groups:10': (
        ('groups:00', '[P;D4;x0;z1](-[*])(-[*])(-[*])-[*]'),
        ('groups:11', '[N;D4;z1](-[*])(-[*])(-[*])-[*]'),
    ),
    'groups:14': (
        ('groups:15', '[N;D3;z5](=[N;D3;z2;+])(=[O;D1])[A]'),
        ('groups:16', '[N;D3;z5](=[N;D3;z2;+])(=[N;D1,D2;z2])[A]'),
    ),
    'groups:29': (
        ('groups:31', '[N;D2;z3;x1]([N,O,S;D1])#[C;D1,D2]'),
        ('groups:32', '[N;D2;z3;x1]([N;D2;z1])#[C;D1,D2]'),
    ),
    'groups:30': (
        ('groups:33', '[N;D2;z3;x1;+]([N,O,S;D1])#[C;D1;-]'),
        ('groups:34', '[N;D2;z3;x1;+]([N;D2;z1])#[C;D1;-]'),
    ),
    'groups:69': (
        ('groups:78', '[S;D4;z2](=[O;D1])[O;D1]'),
        ('groups:79', '[S;D4;z2](=[O;D1])[N;D1,D2;z1]'),
    ),
    'groups:73': (
        ('groups:83', '[C;D2;z3;x1]([N,O,S;D1])#[C;D1,D2]'),
        ('groups:84', '[C;D2;z3;x1]([N;D2;z1])#[C;D1,D2]'),
    ),
    'groups:78': (
        ('groups:62', '[P;D4;z1;+][O;D1;-]'),
        ('groups:70', '[S,Se,Si;D3;z1;+][O;D1;-]'),
        ('groups:71', '[S;D2,D4;z2;+][O;D1;-]'),
        ('groups:90', '[Cl,Br,I;D2;z1;+][O;D1;-]'),
    ),
    'groups:79': (
        ('groups:72', '[S;D3;z2;+2]([O;D1;-])[O;D1;-]'),
        ('groups:73', '[S;D4;z1;+2]([O;D1;-])[O;D1;-]'),
        ('groups:91', '[Cl,Br,I;D3;z1;+2]([O;D1;-])[O;D1;-]'),
    ),
    'metals:11': (
        ('metals:11', '[M;*;^,!^:1]-1-2-3-4-[C:2]-5-[C:3]-1-[C:4]-2-[C:5]-3-[C:6]-4-5'),
        ('metals:12', '[M;*;^,!^:1]-1-2-3-4-[C:2]-5-[C:3]-1=[C:4]-2-[C:5]-3=[C:6]-4-5'),
    ),
    'metals:16': (
        ('metals:17', '[M;*;^,!^:1]-[P;D4;z1;+:2]-C'),
        ('metals:19', '[M;*;^,!^:1]-[N;z1;+:2]-C'),
    ),
    'metals:17': (
        ('metals:18', '[M;*;^,!^:1]-[P;D4;z1:2]-C'),
        ('metals:20', '[M;*;^,!^:1]-[N;z1:2]-C'),
    ),
}

# The one row deleted outright, and the row that made it redundant.  `[C;D1,D2,D3;z1;+]-[N;D3;z1;x0]`
# is a proper subset of `[C;D1,D2,D3;z1;+]-[N;D3;z1]` standing one row in front of it, with the same
# patch and the same `after`, so no molecule can tell the two tables apart.
DELETED = {'groups:87': 'groups:76'}

# The one row in `MERGES` that is a widening: its two members differed in whether the cyclopentadienyl's
# two ring double bonds were drawn and the merged row says `-,=` on both, so it also matches the
# half-drawn ring -- deliberate, since that is garbage input the row exists to repair.  The union proof
# asserts a superset for this row and equality for the other ten; a second entry needs its own approval.
WIDENED = frozenset(('metals:11',))

# The rows written after the port, which no V2 rule stands behind: all three repair a NEUTRAL
# OVER-VALENT NITROGEN drawn Kekule, the spelling `kekule()` repairs in the reader whenever there is an
# aromatic system to resolve and the ported table covered only for the radical (`groups:34`,
# `groups:35`) and four-coordinate (`groups:10`, `groups:33`) cases.  `derive` does not know them and
# would drop them.  Gated by their own `examples` cells and by `test_standardize_overvalent_nitrogen.py`;
# they never fire on the V2 gate corpus, which is why `test_standardize_groups_port.py` lists them
# UNREACHED.
ADDED = {
    'groups:82': 'an N-oxide drawn as a hydroxyl on a nitrogen that already has a double bond',
    'groups:83': 'the same, half separated already, the oxygen holding the minus',
    'groups:84': 'a three-coordinate sp2 nitrogen with nothing to take a counter-charge: a cation',
}


# The TSV's header comment is the single copy of the column documentation; this reads it back so that
# documenting a column in the file a chemist opens cannot leave the generator describing another table.
def _preamble(path):
    """Every leading `#` line of a TSV, verbatim and newline-terminated."""
    lines = []
    for line in path.read_text(encoding='utf-8').split('\n'):
        if not line.startswith('#'):
            break
        lines.append(line)
    if not lines:
        raise ValueError(f'{path} has no header comment; the columns are documented there')
    return ''.join(f'{line}\n' for line in lines)


PREAMBLE = {tag: _preamble(tsv) for tag, _, tsv in SOURCES}


class Rule:
    """One row: a pattern, the patch it applies, and the claim it makes."""
    __slots__ = ('id', 'smarts', 'atom_fix', 'bonds_fix', 'tautomer', 'after', 'examples', 'why',
                 'lineno')

    def __init__(self, id, smarts, atom_fix, bonds_fix, tautomer, after=(), examples=(), why='',
                 lineno=0):
        self.id = id
        self.smarts = smarts
        self.atom_fix = atom_fix        # list of (slot, charge delta, None | False | True)
        self.bonds_fix = bonds_fix      # list of (a, b, order)
        self.tautomer = tautomer
        self.after = tuple(after)       # ids this rule must follow
        self.examples = tuple(examples)  # one 'IN>>OUT' per alternative of the pattern, both SMILES
        self.why = why
        self.lineno = lineno            # of the `smarts(...)` call in V2; not emitted

    def row(self):
        atom_fix = ';'.join(f'{n}:{d}:{RADICAL_GLYPH[r]}' for n, d, r in self.atom_fix) or EMPTY
        bonds_fix = ';'.join(f'{a}:{b}:{o}' for a, b, o in self.bonds_fix) or EMPTY
        return (self.id, self.smarts, atom_fix, bonds_fix, '1' if self.tautomer else '0',
                ';'.join(self.after) or EMPTY, ';'.join(self.examples) or EMPTY, self.why or EMPTY)


def parse_atom_fix(text):
    if text == EMPTY:
        return []
    out = []
    for entry in text.split(';'):
        slot, delta, radical = entry.split(':')
        if radical not in GLYPH_RADICAL:
            raise ValueError(f'atom_fix radical must be one of -01, got {radical!r}')
        out.append((int(slot), int(delta), GLYPH_RADICAL[radical]))
    return out


def parse_bonds_fix(text):
    if text == EMPTY:
        return []
    out = []
    for entry in text.split(';'):
        a, b, order = entry.split(':')
        out.append((int(a), int(b), int(order)))
    return out


def parse_after(text):
    if text == EMPTY:
        return ()
    return tuple(text.split(';'))


def parse_examples(text):
    """`;`-separated `IN>>OUT`s.  `;` is free as a separator: it is a SMARTS character, and this
    column is SMILES."""
    if text == EMPTY:
        return ()
    return tuple(text.split(';'))


def comment_above(lines, lineno):
    """The contiguous run of `#` lines immediately above `lineno`, as one collapsed string.

    This is why the extraction reads source text rather than compiled rules: V2's ASCII art is where a
    rule says what it is for, and it is thrown away at import.  Blank comment lines pad the art and are
    dropped rather than emitted as empty ` / ` runs.
    """
    block = []
    i = lineno - 2                                  # 0-based index of the line above
    while i >= 0 and lines[i].lstrip().startswith('#'):
        block.append(lines[i].lstrip()[1:])
        i -= 1
    block.reverse()
    parts = [re.sub(r'\s+', ' ', text).strip() for text in block]
    return ' / '.join(p for p in parts if p)


def extract(path):
    """Every `rules.append((...))` in `_rules`, in declaration order, from the source text.

    Names resolve against the assignments preceding the append, which is how V2 spells a rule.  A `q`
    that is not rebound is a rule reusing the previous pattern, and `_groups.py` does that on purpose.
    `_metal_organics.py` appends 3-tuples, so a missing fourth element reads as `is_tautomer` False.
    """
    text = path.read_text(encoding='utf-8')
    lines = text.splitlines()
    tree = ast.parse(text, filename=str(path))

    body = None
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name == '_rules':
            body = node.body
            break
    if body is None:
        raise ValueError(f'{path}: no `_rules` function')

    env = {}                                        # name -> value
    smarts_lineno = {}                              # name -> lineno of its smarts() call
    out = []
    for stmt in body:
        if isinstance(stmt, ast.Assign) and len(stmt.targets) == 1 and \
                isinstance(stmt.targets[0], ast.Name):
            name = stmt.targets[0].id
            value = stmt.value
            if isinstance(value, ast.Call) and isinstance(value.func, ast.Name) and \
                    value.func.id == 'smarts':
                env[name] = ast.literal_eval(value.args[0])
                smarts_lineno[name] = value.lineno
            else:
                try:
                    env[name] = ast.literal_eval(value)
                except (ValueError, TypeError, SyntaxError):
                    env.pop(name, None)             # not a literal, so not rule data
            continue

        if not isinstance(stmt, ast.Expr) or not isinstance(stmt.value, ast.Call):
            continue
        call = stmt.value
        if not (isinstance(call.func, ast.Attribute) and call.func.attr == 'append' and
                isinstance(call.func.value, ast.Name) and call.func.value.id == 'rules'):
            continue
        if len(call.args) != 1 or not isinstance(call.args[0], ast.Tuple):
            raise ValueError(f'{path}:{call.lineno}: rules.append with an unexpected argument')

        fields = []
        for element in call.args[0].elts:
            if isinstance(element, ast.Name):
                fields.append(env[element.id])
            else:
                fields.append(ast.literal_eval(element))
        smarts, atom_fix, bonds_fix = fields[0], fields[1], fields[2]
        tautomer = fields[3] if len(fields) > 3 else False
        lineno = smarts_lineno[call.args[0].elts[0].id]

        out.append(Rule(id=None,
                        smarts=smarts,
                        atom_fix=[(n, d, r) for n, (d, r) in atom_fix.items()],
                        bonds_fix=[tuple(b) for b in bonds_fix],
                        tautomer=bool(tautomer),
                        why=comment_above(lines, lineno),
                        lineno=lineno))
    return out


def translate_z(rules, tag):
    """Apply `Z3_MAP` to every pattern carrying a `z3`, and refuse to guess about one it misses.

    Everything here fails rather than warns except the line numbers, which are only a cross-reference
    to the mapping document and drift harmlessly.  An unmapped `z3` is not harmless -- see `Z3_MAP`.
    """
    if tag == 'metals':
        for rule in rules:
            if 'z3' in rule.smarts:
                raise SystemExit(f'{tag}: unexpected `z3` in {rule.smarts!r}; _metal_organics.py '
                                 f'is documented to have none.  Add it to Z3_MAP with a ruling')
        return {}

    seen = set()
    counts = {'A': 0, 'C': 0, 'D': 0}
    drift = []
    for rule in rules:
        occurrences = rule.smarts.count('z3')
        if not occurrences:
            if rule.smarts in Z3_MAP:
                raise SystemExit(f'{tag}: Z3_MAP has a ruling for {rule.smarts!r}, which has no '
                                 f'`z3` in it')
            continue
        if rule.smarts not in Z3_MAP:
            raise SystemExit(f'{tag}:{rule.lineno}: no `z3` ruling for {rule.smarts!r}.  V2\'s '
                             f'`z3` is not the core\'s -- see the mapping document, then add a row '
                             f'to Z3_MAP.  Refusing to copy it through')
        group, doc_lineno = Z3_MAP[rule.smarts]
        if rule.smarts not in seen:
            # `groups:81`/`groups:82` share one `q`, so a primitive is counted per PATTERN and not per
            # rule -- otherwise the population is 37 and the cross-check against the mapping document,
            # which counts primitives in the source, fails for a spurious reason
            seen.add(rule.smarts)
            counts[group] += occurrences
            if rule.lineno != doc_lineno:
                drift.append((rule.lineno, doc_lineno, rule.smarts))
        rule.smarts = rule.smarts.replace('z3', TARGET[group])

    missing = set(Z3_MAP) - seen
    if missing:
        raise SystemExit(f'{tag}: Z3_MAP rules nothing in the source matched: '
                         + ', '.join(sorted(missing)))

    total = sum(counts.values())
    if total != Z3_PRIMITIVES or len(seen) != Z3_LINES:
        raise SystemExit(f'{tag}: found {total} `z3` primitives on {len(seen)} lines, expected '
                         f'{Z3_PRIMITIVES} on {Z3_LINES}')
    if drift:
        print('!' * 78, file=sys.stderr)
        print(f'!! {len(drift)} `z3` rule(s) have moved since the mapping document was written.',
              file=sys.stderr)
        print('!! The rulings were applied anyway -- Z3_MAP is keyed by SMARTS, not by line -- but',
              file=sys.stderr)
        print('!! docs/superpowers/research/2026-09-03-z3-port-mapping.md now cites wrong lines.',
              file=sys.stderr)
        for actual, expected, smarts in drift:
            print(f'!!   {expected} -> {actual}  {smarts}', file=sys.stderr)
        print('!' * 78, file=sys.stderr)
    return counts


def translate_degree(rules, tag):
    """Apply `DEGREE_MAP`, and refuse to emit a ruling that matched nothing.

    Unlike `z`, a `D` needing a rewrite is indistinguishable from one that does not by looking at the
    string -- the primitive's meaning moved, not its spelling -- so the map is exhaustive by inspection
    and the only mechanical check is that every entry still has a rule to rewrite.
    """
    used = set()
    for rule in rules:
        if rule.smarts in DEGREE_MAP:
            used.add(rule.smarts)
            rule.smarts = DEGREE_MAP[rule.smarts]
    missing = set(DEGREE_MAP) - used
    if missing:
        raise SystemExit(f'{tag}: DEGREE_MAP rewrites nothing in the source matched: '
                         + ', '.join(sorted(missing)) + '.  Either the pattern moved or the ruling '
                         'is stale; re-read it against the source before deleting it')
    return used


def derive(tag, source):
    """The patch columns from V2; `after`, `examples` and `why` read back from the TSV by id.

    Those three are hand-measured, so a re-derivation must not regenerate them.  `why` falls back to
    V2's comment block only for a row the TSV does not have; the other two are left empty for a human,
    because an example nobody executed is worse than no example.
    """
    rules = extract(source)
    counts = translate_z(rules, tag)
    if tag == 'groups':
        translate_degree(rules, tag)
    tsv = next(path for name, _, path in SOURCES if name == tag)
    prose = {r.id: r for r in read_tsv(tsv)} if tsv.exists() else {}
    for i, rule in enumerate(rules):
        rule.id = f'{tag}:{i:02d}'
        kept = prose.get(rule.id)
        if kept is not None:
            rule.after, rule.examples = kept.after, kept.examples
            rule.why = kept.why or rule.why
    return rules, counts


def render(tag, rules):
    lines = [PREAMBLE[tag], '\t'.join(HEADER)]
    for rule in rules:
        row = rule.row()
        for name, cell in zip(HEADER, row):
            if '\t' in cell or '\n' in cell:
                raise SystemExit(f'{rule.id}: column {name} contains a tab or a newline')
            if not cell:
                raise SystemExit(f'{rule.id}: column {name} is empty; write {EMPTY!r} instead')
        lines.append('\t'.join(row))
    return '\n'.join(lines) + '\n'


def read_tsv(path):
    """The TSV back into `Rule`s.  This is what a consumer of the collection uses."""
    rules = []
    for lineno, line in enumerate(path.read_text(encoding='utf-8').splitlines(), 1):
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        fields = line.split('\t')
        if tuple(fields) == HEADER:
            continue
        if len(fields) != len(HEADER):
            raise ValueError(f'{path}:{lineno}: {len(fields)} fields, expected {len(HEADER)}')
        id, smarts, atom_fix, bonds_fix, tautomer, after, examples, why = fields
        if tautomer not in ('0', '1'):
            raise ValueError(f'{path}:{lineno}: tautomer must be 0 or 1, got {tautomer!r}')
        rules.append(Rule(id, smarts, parse_atom_fix(atom_fix), parse_bonds_fix(bonds_fix),
                          tautomer == '1', parse_after(after), parse_examples(examples),
                          '' if why == EMPTY else why))
    return rules


def check():
    """Re-derive into memory and diff.  Returns the number of files that disagree."""
    bad = 0
    for tag, source, tsv in SOURCES:
        rules, _ = derive(tag, source)
        want = render(tag, rules)
        have = tsv.read_text(encoding='utf-8') if tsv.exists() else ''
        if want == have:
            print(f'{tsv.name}: current, {len(rules)} rows')
            continue
        bad += 1
        print(f'{tsv.name}: DIFFERS from what `derive` produces')
        sys.stdout.writelines(difflib.unified_diff(have.splitlines(keepends=True),
                                                   want.splitlines(keepends=True),
                                                   fromfile=f'{tsv.name} (checked in)',
                                                   tofile=f'{tsv.name} (derived)'))
    return bad


def main(argv):
    verb = argv[0] if argv else 'check'
    if verb == 'derive':
        for tag, source, tsv in SOURCES:
            rules, counts = derive(tag, source)
            tsv.write_text(render(tag, rules))
            note = ''
            if counts:
                note = (f'; `z3` -> ' + ', '.join(f'{n} {TARGET[g]}'
                                                  for g, n in sorted(counts.items()) if n))
            print(f'{tsv}: {len(rules)} rows{note}')
    elif verb == 'check':
        raise SystemExit(1 if check() else 0)
    else:
        raise SystemExit(__doc__)


if __name__ == '__main__':
    main(sys.argv[1:])
