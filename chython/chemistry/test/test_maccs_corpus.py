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
"""The acceptance corpus for `tables/maccs.tsv`: one molecule that must set each key, one that must not.

A mis-transcribed key usually still parses and still matches something, so the must-not-set half is what
catches it.  The corpus rows are the oracle here -- no other implementation is consulted, and a bit that
differs from another toolkit's is a documented difference rather than a defect.
"""
from pytest import mark

from chython.chemistry._maccs import maccs_keys, maccs_match_counts
from chython.chemistry._tables import MACCS_UNSET_KEYS, maccs_corpus, maccs_rules_by_key
from chython.core import read_smiles


CORPUS = [(r.key, r.expectation, r.smiles, r.name) for r in maccs_corpus()]


@mark.parametrize('key,expectation,smi,name', CORPUS)
def test_the_corpus_molecule_sets_or_does_not_set_its_key(key, expectation, smi, name):
    m = read_smiles(smi)
    got = bool(maccs_keys(m)[key])
    want = expectation == 'set'
    assert got == want, (f'key {key} ({maccs_rules_by_key()[key].description!r}) on {name}: '
                         f'expected {"set" if want else "unset"}, matched '
                         f'{maccs_match_counts(m).get(key, 0)} sites')


def test_every_corpus_row_names_a_real_key_and_a_valid_expectation():
    keys = set(maccs_rules_by_key())
    for r in maccs_corpus():
        assert r.key in keys, r.key
        assert r.expectation in ('set', 'unset')
        assert r.name and r.name != '-'
        read_smiles(r.smiles)                # must parse


def test_a_key_with_no_published_definition_has_no_corpus_row():
    """`MACCS_UNSET_KEYS` cannot be exemplified, so a row for one is a mistake, not an omission.

    A `set` row for such a key can never pass, and an `unset` row would pass for the wrong reason --
    every molecule leaves the bit down -- so it would read as evidence for a pattern that does not
    exist.  Refusing both is what keeps the completeness test honest about what it covers.
    """
    offenders = sorted({r.key for r in maccs_corpus()} & set(MACCS_UNSET_KEYS))
    assert not offenders, (f'keys {offenders} have no published definition and ship permanently '
                           f'unset; a corpus row for one asserts nothing')


def test_the_corpus_covers_every_key_with_a_published_definition():
    """All 166 keys minus `MACCS_UNSET_KEYS`, which cannot be exemplified -- see the sibling test."""
    have: dict[int, set[str]] = {}
    for r in maccs_corpus():
        have.setdefault(r.key, set()).add(r.expectation)
    want = [k for k in range(1, 167) if k not in MACCS_UNSET_KEYS]
    missing = [k for k in want if have.get(k) != {'set', 'unset'}]
    assert not missing, f'keys with no set/unset pair: {missing}'
    assert len(want) == 165, 'exactly one key has no published definition'
