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
"""The row merges, re-proved: a merged row selects exactly the sites its members selected.

`gen_standardize_rules.MERGES`/`DELETED` record which rows collapsed into which; this measures that the
record is true.  Separate from `test_standardize_rules_tsv.py`, as it carries pre-merge patterns as
literals.  Sites compare as PATCHED-ATOM TUPLES: whole embeddings differ in atom count, and bare atom
sets lose the number the patch is keyed by.  Ten merges prove equality; the deliberate widening
`metals:11` proves a superset (see `gen.WIDENED`).
"""
import pytest

from . import gen_standardize_rules as gen
from ._corpus import CORPUS
from .._smarts import compile_smarts
from .._standardize import _pass
from .._tables import groups_rules, metals_rules
from ...core import read_smiles


GROUPS = groups_rules()
METALS = metals_rules()
BY_ID = {rule.id: rule for rule in GROUPS + METALS}

# Every SMILES this package can reach without opening a data file: the shared corpus, plus both halves
# of every example in the two tables -- a merged row's examples are exactly the drawings its members
# were chosen to repair.
_SOURCES = list(CORPUS)
for _rule in GROUPS + METALS:
    for _example in _rule.examples:
        _SOURCES.extend(_example.split('>>'))

MOLECULES = []
for _smiles in _SOURCES:
    try:
        MOLECULES.append(read_smiles(_smiles))
    except Exception:                            # pragma: no cover - a corpus entry, not this test
        pass


def _sites(smarts, patch_numbers, molecule):
    """Every site `smarts` selects in `molecule`, as the tuple of atoms its patch would write.

    `patch_numbers` are taken from the SURVIVING row, which is legitimate only because the merge
    required every member to carry an identical patch.
    """
    query, numbers, _ = compile_smarts(smarts)
    out = set()
    for mapping in query.get_mapping(molecule):
        out.add(tuple(mapping[numbers[n]] for n in patch_numbers))
    return out


def _patch_numbers(rule):
    """The atom numbers the row's patch names, in a fixed order.  Non-empty for every merged row --
    a rule with no patch would be a no-op, and `metals:12` is the only patchless row in either table."""
    numbers = list(dict.fromkeys([n for n, _, _ in rule.atom_fix]
                                 + [n for pair in rule.bonds_fix for n in pair[:2]]))
    assert numbers, f'{rule.id}: no patch, so there is nothing to compare sites by'
    return numbers


@pytest.mark.parametrize('surviving', sorted(gen.MERGES), ids=sorted(gen.MERGES))
def test_a_merged_row_selects_what_its_members_selected(surviving):
    """The union proof, per merged row, over every molecule this package has."""
    rule = BY_ID[surviving]
    numbers = _patch_numbers(rule)
    widened = surviving in gen.WIDENED
    for molecule in MOLECULES:
        merged = _sites(rule.smarts, numbers, molecule)
        union = set()
        for _, pattern in gen.MERGES[surviving]:
            union |= _sites(pattern, numbers, molecule)
        if widened:
            assert union <= merged, (
                f'{surviving} is recorded as a WIDENING of {union - merged} and must still match '
                f'everything its members did on {format(molecule, "")!r}')
        else:
            assert merged == union, (
                f'{surviving} does not select what {[i for i, _ in gen.MERGES[surviving]]} selected '
                f'on {format(molecule, "")!r}: gained {sorted(merged - union)}, '
                f'lost {sorted(union - merged)}')


def test_the_widening_is_the_only_one_and_it_is_the_half_drawn_ring():
    """`metals:11` matches a cyclopentadienyl with ONE of its two ring double bonds drawn.

    The rows it replaced matched the all-single and both-double spellings and nothing between.  A
    half-drawn ring is garbage input of the kind the row exists to repair, which is why the merge was
    allowed to widen here and nowhere else.
    """
    assert gen.WIDENED == frozenset(('metals:11',))
    rule = BY_ID['metals:11']
    numbers = _patch_numbers(rule)
    half = read_smiles('[Fe]1234C5C1=C2C3C45')
    assert _sites(rule.smarts, numbers, half), 'the merged row must match the half-drawn ring'
    for _, pattern in gen.MERGES['metals:11']:
        assert not _sites(pattern, numbers, half), (
            f'{pattern!r} already matched the half-drawn ring, so `metals:11` is not a widening '
            'after all and the record is wrong')


def test_the_deleted_row_was_a_subset_of_the_row_that_replaced_it():
    """`groups:87` is gone because `groups:76` matches everything it did, with the same patch.

    There is no query-into-query substructure test, so the containment is measured over molecules.
    Both halves are asserted: a containment that holds because NEITHER pattern ever matches proves
    nothing.
    """
    assert gen.DELETED == {'groups:87': 'groups:76'}
    surviving = BY_ID['groups:76']
    deleted = '[C;D1,D2,D3;z1;+]-[N;D3;z1;x0]'
    assert surviving.smarts == '[C;D1,D2,D3;z1;+]-[N;D3;z1]', (
        'the surviving row has been edited; the containment below is about the pattern it had when '
        f'{deleted!r} was deleted, not about whatever it says now')
    numbers = _patch_numbers(surviving)
    hits = 0
    for molecule in MOLECULES + [read_smiles(s) for s in ('C[N+](C)=C', 'C[N+](C)(C)C',
                                                          'CN(C)[CH2+]', 'CN(O)[CH2+]')]:
        narrow = _sites(deleted, numbers, molecule)
        wide = _sites(surviving.smarts, numbers, molecule)
        hits += len(narrow)
        assert narrow <= wide, (
            f'{deleted!r} selects {sorted(narrow - wide)} on {format(molecule, "")!r} and '
            f'{surviving.smarts!r} does not, so deleting the narrow row dropped a repair')
    assert hits, f'{deleted!r} matches nothing anywhere, so the containment above is vacuous'
    # the surviving row is STRICTLY wider, which is why the narrow one was the one dropped: an amine
    # with a heteroatom neighbour fails the deleted row's `x0`
    strictly_wider = read_smiles('CN(O)[CH2+]')
    assert not _sites(deleted, numbers, strictly_wider)
    assert _sites(surviving.smarts, numbers, strictly_wider)


def test_no_merged_row_reintroduces_a_z3_the_port_ruled_out():
    """The `z` translation survives the merge, on the alternative that needed it.

    `test_standardize_rules_tsv.py` checks the whole table; here it is checked where a merge could
    have lost it silently -- the eight `Z3_MAP` rulings absorbed into a merged row.
    """
    absorbed = {pattern: surviving for surviving, members in gen.MERGES.items()
                for _, pattern in members}
    checked = 0
    for original, (group, _) in gen.Z3_MAP.items():
        translated = original.replace('z3', gen.TARGET[group])
        if translated not in absorbed:
            continue
        checked += 1
        rule = BY_ID[absorbed[translated]]
        assert gen.TARGET[group] in rule.smarts, (
            f'{absorbed[translated]} absorbed {translated!r} but carries no {gen.TARGET[group]!r}; '
            'the port ruling was lost in the merge')
        if group != 'A':
            assert 'z3' not in rule.smarts, (
                f'{absorbed[translated]} is group {group} and must not carry a `z3`')
    assert checked == 8, f'{checked} absorbed z3 rulings, not the 8 the merge record accounts for'


def test_the_merged_table_is_the_size_the_record_says():
    """104 rows, and the arithmetic from 116 checks out against `MERGES`, `DELETED` and `ADDED`.

    Catches a row added or dropped without a record: the record must explain the whole difference.
    `ADDED` is the only term that is not port arithmetic -- rows written since, each named there and
    each gated by a test of its own.
    """
    absorbed = sum(len(members) - 1 for members in gen.MERGES.values())
    assert len(GROUPS) + len(METALS) == 116 - absorbed - len(gen.DELETED) + len(gen.ADDED) == 104
    assert len(GROUPS) == 85
    assert len(METALS) == 19
    assert set(gen.ADDED) <= {rule.id for rule in GROUPS}, (
        'a row `ADDED` names is not in the table; the record is of rows that were added, so a removal '
        'has to come out of the record too')
