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
"""The two standardization rule TSVs, read on their own: row counts, that every cell decodes, that every
patch slot is an atom the pattern declares, that no `z3` survived where the port ruled `z5`/`z6`, and
which rows carry the `tautomer` flag.  The three documentation columns are checked for shape only --
whether an example is true and whether an ordering obligation is real are questions no TSV can answer,
and `test_standardize_rules_examples.py` executes both against the core.
"""
import re

import pytest

from . import gen_standardize_rules as gen


GROUPS = gen.read_tsv(gen.TABLES / 'standardize_groups.tsv')
METALS = gen.read_tsv(gen.TABLES / 'standardize_metals.tsv')
ALL = GROUPS + METALS

V2_PRESENT = (gen.V2 / '_groups.py').exists() and (gen.V2 / '_metal_organics.py').exists()

# atom tokens.  Two-letter symbols first, or `Cl` reads as a carbon followed by a chlorine
CX_SUFFIX = re.compile(r'\s*\|[^|]*\|\s*$')
BARE_ATOM = re.compile(r'Cl|Br|Si|Se|[BCNOPSFIbcnops]')
MAP_NUMBER = re.compile(r':(\d+)$')
# ring closures, bond orders, branches and the `,` of a bond-order alternative.  `^` is here as the
# DATIVE BOND, which the metals table writes between atoms; the `^` of a radical primitive is inside a
# bracket and never reaches this set, because a `[` consumes through its `]`
STRUCTURE = set('-=#:~^/\\,().0123456789')


def atom_numbers(smarts):
    """The atom numbering a query's patch slots refer to, recomputed from the pattern text.

    Walk the atoms in declaration order; an atom's number is its explicit `:N` map where one is
    written, otherwise the lowest integer from 1 that no map claims and no earlier unmapped atom took.
    So the first atom of a pattern that maps three of its later atoms is 4.  Recomputed rather than
    read off a compiled query so that the TSV stays checkable without importing anything.
    """
    text = CX_SUFFIX.sub('', smarts)
    written = []
    i = 0
    while i < len(text):
        if text[i] == '[':
            j = text.index(']', i)
            found = MAP_NUMBER.search(text[i + 1:j])
            written.append(int(found.group(1)) if found else None)
            i = j + 1
            continue
        bare = BARE_ATOM.match(text, i)
        if bare:
            written.append(None)
            i = bare.end()
            continue
        if text[i] in STRUCTURE:
            i += 1
            continue
        raise ValueError(f'{smarts!r}: cannot read {text[i]!r} at {i}')

    taken = {n for n in written if n is not None}
    out = []
    nxt = 1
    for n in written:
        if n is None:
            while nxt in taken:
                nxt += 1
            n = nxt
            taken.add(n)
        out.append(n)
    return out


@pytest.mark.skipif(not V2_PRESENT,
                    reason='chython 2 standardize rules are gone; the TSV is the authority now')
def test_check_mode_passes():
    """The checked-in TSVs are what `derive` produces from the V2 source as it stands."""
    assert gen.check() == 0, 'run `python chython/chemistry/test/gen_standardize_rules.py derive`'


def test_row_counts():
    """85 and 19.  A silent gain or loss of a rule is a behaviour change.  `gen_standardize_rules.MERGES`
    and `DELETED` say which rows collapsed into which, and the three rows past the port -- `groups:82` to
    `groups:84`, the neutral over-valent nitrogen written Kekule -- are gated by
    `test_standardize_overvalent_nitrogen.py`."""
    assert len(GROUPS) == 85
    assert len(METALS) == 19


@pytest.mark.parametrize('rule', ALL, ids=[r.id for r in ALL])
def test_every_row_decodes(rule):
    """Both patch columns round-trip, and the numbers in them are real atoms of the pattern."""
    assert rule.smarts
    row = rule.row()
    assert gen.parse_atom_fix(row[2]) == rule.atom_fix
    assert gen.parse_bonds_fix(row[3]) == rule.bonds_fix

    declared = set(atom_numbers(rule.smarts))
    for slot, delta, radical in rule.atom_fix:
        assert slot in declared, f'{rule.id}: atom_fix slot {slot} is not an atom of {rule.smarts}'
        assert -4 <= delta <= 4, f'{rule.id}: implausible charge delta {delta}'
        assert radical in (None, False, True)
    for a, b, order in rule.bonds_fix:
        assert a in declared, f'{rule.id}: bonds_fix atom {a} is not an atom of {rule.smarts}'
        assert b in declared, f'{rule.id}: bonds_fix atom {b} is not an atom of {rule.smarts}'
        assert a != b, f'{rule.id}: bonds_fix names a self-loop on {a}'
        assert order in (1, 2, 3, 4, 8), f'{rule.id}: unknown bond order {order}'


def test_no_cell_holds_a_tab_or_is_empty():
    for rule in ALL:
        for name, cell in zip(gen.HEADER, rule.row()):
            assert cell, f'{rule.id}: column {name} is empty; `-` is how emptiness is spelled'
            assert '\t' not in cell, f'{rule.id}: column {name} holds a tab'


def test_z3_is_gone_where_the_port_ruled_it_z5_or_z6():
    """A `z3` left on one of these patterns would narrow it to alkynes and nitriles and the rule would
    quietly stop firing.  Checked against the emitted file, not the generator's intent."""
    emitted = {rule.smarts for rule in GROUPS}
    absorbed = {pattern for members in gen.MERGES.values() for _, pattern in members}
    for original, (group, _) in gen.Z3_MAP.items():
        translated = original.replace('z3', gen.TARGET[group])
        if translated in absorbed:
            # a merged row does not carry any member's text verbatim; that the ruling survived the
            # merge is checked against the surviving row in `test_standardize_rules_merges.py`
            continue
        assert translated in emitted, f'no row carries the translation of {original!r}'
        if group != 'A':
            assert 'z3' not in translated, f'{original!r} is group {group} but kept a `z3`'
    assert sum(r.smarts.count('z3') for r in GROUPS) == 11, (
        'group A is 14 primitives, three of them on rows that merged into two')
    assert sum(r.smarts.count('z5') for r in GROUPS) == 15, (
        'group C is 15, one of them duplicated and the duplicate absorbed by a merge')
    assert sum(r.smarts.count('z6') for r in GROUPS) == 7
    for rule in METALS:
        assert 'z3' not in rule.smarts


@pytest.mark.parametrize('tag', [tag for tag, _, _ in gen.SOURCES])
def test_the_files_round_trip_through_the_generator(tag):
    """The checked-in bytes are what `render(read_tsv(...))` produces; V2 is not consulted.

    The file is self-consistent under the only reader and writer this package has: every cell survives
    a decode and re-encode, the header row is `HEADER`, the comment block is intact, nothing is
    separated by anything but a tab.  A hand edit `read_tsv` silently drops fails here.  It also keeps
    the column documentation single-copy, since `gen.PREAMBLE` reads the comment block back out.
    """
    tsv = next(path for name, _, path in gen.SOURCES if name == tag)
    assert gen.render(tag, gen.read_tsv(tsv)) == tsv.read_text(encoding='utf-8'), (
        f'{tsv.name} is not what the generator writes from what it reads.  Either a cell holds '
        'something `read_tsv` cannot represent, or the header row no longer matches `gen.HEADER`.')


def test_no_metal_rule_is_a_tautomer_fix():
    """No metal repair moves a hydrogen, so the column is not a place a metal rule opts in."""
    assert not any(rule.tautomer for rule in METALS)


# The 26 rows whose repair moves a hydrogen from one heavy atom to another, which is what the column
# means -- see the TSV header.  Listed rather than counted because the flag is consumed:
# `standardize(fix_tautomers=False)` withholds exactly these, so a row joining or leaving the set is a
# behaviour change for every caller who passes the flag.
TAUTOMER_ROWS = frozenset((
    'groups:11', 'groups:13', 'groups:27', 'groups:30', 'groups:32', 'groups:39', 'groups:42',
    'groups:45', 'groups:46', 'groups:47', 'groups:48', 'groups:49', 'groups:50', 'groups:51',
    'groups:52', 'groups:53', 'groups:54', 'groups:55', 'groups:56', 'groups:57', 'groups:62',
    'groups:64', 'groups:70', 'groups:71', 'groups:72', 'groups:73',
))

# The two rows the port flags differently from V2, deliberately: both are diazo repairs that move a
# hydrogen from nitrogen to carbon (`A-C#N=NH >> A-[CH]=[N+]=[N-]`), which is what the column means.
DISAGREES_WITH_V2 = frozenset(('groups:27', 'groups:32'))


def test_the_tautomer_column_is_the_documented_set():
    """Which rows carry the flag, named.  A count would not catch a swap."""
    carried = {rule.id for rule in GROUPS if rule.tautomer}
    assert carried == TAUTOMER_ROWS, (
        f'the tautomer column has moved.\n  gained: {sorted(carried - TAUTOMER_ROWS)}\n'
        f'  lost:   {sorted(TAUTOMER_ROWS - carried)}\n\nThe column means "the repair moves a '
        'hydrogen between heavy atoms" and `standardize(fix_tautomers=False)` withholds exactly '
        'these rows, so this is a behaviour change and not a documentation edit.')
    assert DISAGREES_WITH_V2 <= carried
    assert len(carried) == 26, (
        'V2 flags 26 of its 94: the two diazo rows in DISAGREES_WITH_V2 are extra, and two merges '
        'each collapsed a flagged PAIR into one flagged row, so the two counts coincide by accident')


def test_the_documented_duplicate_is_a_duplicate():
    """`groups:71` and `groups:72` are one rule appended twice, on purpose: the pattern matches a
    second, overlapping site only after the first has been patched, so one pass leaves it unrepaired.
    If a future edit makes the two rows differ, that intent has been lost."""
    first = next(r for r in GROUPS if r.id == 'groups:71')
    second = next(r for r in GROUPS if r.id == 'groups:72')
    assert first.smarts == second.smarts
    assert first.atom_fix == second.atom_fix
    assert first.bonds_fix == second.bonds_fix
    assert first.tautomer == second.tautomer


def test_ids_are_declaration_order():
    """Row order is semantics -- several rules are documented as order dependent."""
    assert [r.id for r in GROUPS] == [f'groups:{i:02d}' for i in range(len(GROUPS))]
    assert [r.id for r in METALS] == [f'metals:{i:02d}' for i in range(len(METALS))]


@pytest.mark.parametrize('table', [GROUPS, METALS], ids=['groups', 'metals'])
def test_the_after_column_is_well_formed(table):
    """An obligation names earlier rows of the same table, each once.

    Rules run top to bottom, so "must follow a later row" is unsatisfiable and a cross-table entry is
    meaningless -- the group table runs to completion before the metal table starts.  Whether an
    obligation is real is `test_standardize_rules_examples.py`'s business.
    """
    index = {rule.id: i for i, rule in enumerate(table)}
    for i, rule in enumerate(table):
        assert len(set(rule.after)) == len(rule.after), f'{rule.id}: after repeats an id'
        for earlier in rule.after:
            assert earlier in index, f'{rule.id}: after names {earlier!r}, not a row of this table'
            assert index[earlier] < i, f'{rule.id}: after names {earlier!r}, which is not earlier'


def test_every_row_carries_an_executable_example_per_alternative():
    """`IN>>OUT`, one arrow, both halves present, one per `,` alternative of the pattern.

    Counting them is this file's job because it is a property of the text: a merged row whose second
    branch carries no example has an ungated chemical claim, and no execution can notice a missing case.
    """
    for rule in ALL:
        assert rule.examples, f'{rule.id}: no example.  Every row is gated by its own row'
        for example in rule.examples:
            assert example.count('>>') == 1, f'{rule.id}: {example!r} is not a single `IN>>OUT`'
            source, product = example.split('>>')
            assert source and product, f'{rule.id}: {example!r} has an empty side'
        # a floor and not an equality: not every `,` is a separate claim (`[C;D1,D2]` is one), so a row
        # states at least as many examples as it has merged members
        members = len(gen.MERGES.get(rule.id, (None,)))
        assert len(rule.examples) >= members, (
            f'{rule.id} absorbed {members} rows and carries {len(rule.examples)} examples; each '
            'merged member brought its own and none may be dropped')


def test_every_row_says_why_in_prose():
    """The `why` column, which is also the log message: prose about the chemistry, long enough to say
    what was drawn wrong and never a transcribed comment block."""
    for rule in ALL:
        assert len(rule.why) >= 40, f'{rule.id}: `why` is too short to say what was drawn wrong'
        assert ' / ' not in rule.why, (
            f'{rule.id}: `why` holds a ` / `, the separator V2\'s collapsed comment blocks used.  '
            'This column is prose about the chemistry, not transcribed ASCII art')
