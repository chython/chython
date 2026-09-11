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
"""`tables/maccs.tsv` as a table: the schema, the four kinds and the two dialect traps."""
from re import findall, finditer

from chython.chemistry._tables import (MACCS_KINDS, MACCS_PREDICATES, MACCS_UNSET_KEYS, maccs_rules,
                                       maccs_rules_by_key)


def implicit_junctions(pattern):
    """Every place a bracket atom is followed by another bracket atom with no bond token between them.

    Returns the text found between each such pair -- `''` for `[C;z4][C;z4]`, `'1'` for a ring closure,
    `'('` for a branch -- so a caller can assert the list is empty and see what it found when it is not.
    Ring-closure digits, `(`, `)` and `%` are the only things allowed between two brackets without a
    bond token; anything in that set still leaves the bond implicit, which is the whole point.
    """
    return [m.group('between') for m in finditer(r'\](?P<between>[()%\d]*)\[', pattern)]


def test_the_keys_are_exactly_one_through_one_hundred_and_sixty_six():
    keys = [r.key for r in maccs_rules()]
    assert keys == list(range(1, 167))          # ordered, complete, no duplicate, no gap


def test_every_row_is_well_formed():
    for r in maccs_rules():
        assert r.id == f'maccs:{r.key}'
        assert r.kind in MACCS_KINDS
        assert r.description and r.description != '-'
        if r.kind == 'smarts':
            assert r.query is not None and r.count == 1 and r.predicate == '-'
        elif r.kind == 'count':
            assert r.query is not None and r.count >= 2 and r.predicate == '-'
        elif r.kind == 'predicate':
            assert r.query is None and r.count == 0 and r.predicate in MACCS_PREDICATES
        else:
            assert r.query is None and r.count == 0 and r.predicate == '-'


def test_a_description_saying_more_than_one_is_a_count_row():
    # the published wording "> 1" is a threshold, and a key that reads it as presence is wrong.
    for r in maccs_rules():
        if '> 1' in r.description and r.kind != 'predicate':
            assert r.kind == 'count' and r.count >= 2, r.id


def test_every_bracket_atom_frees_the_charge():
    # an unstated charge means charge ZERO in chython, and a MACCS key never means neutral-only.
    for r in maccs_rules():
        if r.query is None:
            continue
        for bracket in findall(r'\[[^]]*\]', r.pattern):
            assert '*' in bracket or '+' in bracket or '-' in bracket, (r.id, bracket)


def test_no_pattern_leaves_a_bond_implicit():
    """An absent bond matches SINGLE ONLY, so `[C;z4][C;z4]` is a dead pattern and `[C;z4]:[C;z4]` is not.

    Every `%A`, `$A`, `!A` and bare-`A` shorthand in the MACCS legend translates to an EXPLICIT bond
    token, so a pattern with an implicit junction is a transcription slip by construction -- there is no
    MACCS key whose correct reading is "single bond only between two atoms written as `A`".

    `implicit_junctions` looks for the SHAPE and not for a spelling: a gate written as a substring search
    for `']z4]['` cannot fire on any legal chython SMARTS at all, which is worse than no gate.
    """
    for r in maccs_rules():
        if r.query is None:
            continue
        found = implicit_junctions(r.pattern)
        assert not found, (r.id, r.pattern, f'bond left implicit at {found}; an absent bond matches '
                                            f'single only -- write the bond token')


def test_the_implicit_bond_gate_can_fail():
    """Negative control, because a gate green on every input is a gate that has stopped gating.

    Positive cases: two aromatics juxtaposed, the same with a wildcard, a heteroatom pair, and a ring
    written with no bond tokens at all.  Negative cases: the same patterns with the bond stated.
    """
    assert implicit_junctions('[C;z4][C;z4]') == ['']
    assert implicit_junctions('[*;z4][*;z4]') == ['']
    assert implicit_junctions('[O;*][C;*]') == ['']
    assert implicit_junctions('[*]1[*][*][*]1') == ['1', '', '']
    assert implicit_junctions('[C;z4]:[C;z4]') == []
    assert implicit_junctions('[O;*]=[C;*]-[N;*]') == []
    assert implicit_junctions('[*]~;@[*]~;!@[*]') == []


def test_every_predicate_name_is_registered():
    used = {r.predicate for r in maccs_rules() if r.kind == 'predicate'}
    assert used <= set(MACCS_PREDICATES)


def test_the_unset_keys_are_exactly_the_ones_with_no_published_definition():
    """Key 44's published description is the placeholder `OTHER`; there is nothing to transcribe.

    Asserted as an equality in both directions, so neither a new `unset` row nor a quiet promotion of
    key 44 to a guessed pattern can happen without this failing.  A row shipping unset must say why in
    its own description -- a caller reading the table is the person this row exists for.
    """
    rows = maccs_rules_by_key()
    assert tuple(sorted(r.key for r in maccs_rules() if r.kind == 'unset')) == MACCS_UNSET_KEYS
    for key in MACCS_UNSET_KEYS:
        r = rows[key]
        assert r.pattern == '-' and r.predicate == '-' and r.count == 0 and r.query is None
        assert 'unset' in r.description.lower()
