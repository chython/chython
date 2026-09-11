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
"""The `examples` and `after` columns of the rule tables, executed -- one gate per `,` alternative.

An example asserts `standardize()` turns IN into OUT, that the row's own rule fired and nothing else.
An `after` obligation is measured by running the pair alone in both orders; a stale one fails here.
"""
import pytest

from .._implicit import calc_implicit
from .._standardize import _pass, standardize
from .._tables import groups_rules, metals_rules
from ...core import read_smiles


GROUPS = groups_rules()
METALS = metals_rules()
ALL = GROUPS + METALS

# `groups:50` is unreachable as shipped and kept rather than deleted -- see its `why`.  Every site it
# matches is also a `groups:47` site and `groups:47` runs first, giving the same product, so the
# shadowing costs nothing.  If an edit makes `groups:50` reachable, this entry stops being true.
SHADOWED = {'groups:50': 'groups:47'}

# `groups:72` is `groups:71` appended twice on purpose, so its example legitimately reports both.
ALSO_FIRES = {'groups:72': frozenset(('groups:71',))}

# Metal rows whose product does not equal a plain re-read of its own SMILES, for a reason unrelated to
# the rule: the repaired metal has no valence row, so `calc_implicit` leaves it `H_UNKNOWN`, while a
# bracket atom in the expectation states 0 hydrogens.  For these the expectation is re-derived the same
# way the product was -- see `_derived`.  Named rather than counted, so a new row taking this path has
# to be looked at first.
NEEDS_DERIVED_HYDROGENS = frozenset((
    'metals:01', 'metals:03', 'metals:04', 'metals:05', 'metals:06', 'metals:11', 'metals:12',
    'metals:14', 'metals:15', 'metals:16', 'metals:18',
))


def _derived(smiles):
    """A molecule read from SMILES with every implicit hydrogen count re-derived.

    Only for `NEEDS_DERIVED_HYDROGENS`.  Legitimate because no rule writes a hydrogen count, so
    deriving both sides compares the repair rather than the expectation's spelling of a consequence.
    """
    molecule = read_smiles(smiles)
    for atom in list(molecule.atoms()):
        calc_implicit(molecule, atom.n)
    return molecule


def _run(smiles, table):
    """One rule table over one molecule, hydrogens recomputed -- `standardize()` with a table of choice.

    `_pass` deliberately does not recompute; `standardize` does, over the atoms it wrote.  Reading a
    product without the recompute gives a stale hydrogen count and a misleading SMILES.
    """
    molecule = read_smiles(smiles)
    for n in sorted(_pass(molecule, table, molecule.log)):
        calc_implicit(molecule, n)
    return molecule


# A row carries one example per `,` alternative of its pattern, so the gate is per example and not per
# row: a merged row whose first alternative works and whose second does not would otherwise pass on the
# strength of the first.  `[P&D4&x0&z1,N&D4&z1]` is a phosphonium claim and an ammonium claim.
EXAMPLES = tuple((rule, example) for rule in ALL for example in rule.examples)
EXAMPLE_IDS = [f'{rule.id}-{rule.examples.index(example)}' for rule, example in EXAMPLES]


@pytest.mark.parametrize('rule,example', EXAMPLES, ids=EXAMPLE_IDS)
def test_the_example_is_what_the_pass_does(rule, example):
    """IN standardizes to OUT.  Every row of both tables carries at least one; this executes them."""
    assert rule.examples, f'{rule.id}: no example.  Every rule is gated by its own row'
    source, expected = example.split('>>')
    molecule = read_smiles(source)
    assert standardize(molecule), f'{rule.id}: nothing changed on {source!r}'
    if rule.id in NEEDS_DERIVED_HYDROGENS:
        assert molecule == _derived(expected)
    else:
        assert molecule == read_smiles(expected), (
            f'{rule.id}: {source!r} standardizes to {format(molecule, "")!r}, not {expected!r}')


@pytest.mark.parametrize('rule,example', EXAMPLES, ids=EXAMPLE_IDS)
def test_the_example_fires_this_rule_and_no_other(rule, example):
    """The row's own id is in the log, and nothing else is.

    An input some earlier rule repairs first would pass the test above while proving nothing about this
    row.
    """
    source = example.split('>>')[0]
    molecule = read_smiles(source)
    log = molecule.log
    standardize(molecule)
    # the reader writes to `mol.log` too, and an input drawn wrong enough to reach a rule is often
    # drawn wrong enough for the reader to have said so first: the stage is what says whose record it is
    fired = {record.rule for record in log if record.stage != 'read'}
    expected = {SHADOWED.get(rule.id, rule.id)} | ALSO_FIRES.get(rule.id, frozenset())
    assert fired == expected, (
        f'{rule.id}: {source!r} fired {sorted(fired)}.  An example has to exercise its own row -- '
        'if an earlier rule gets there first, the example is measuring that rule instead')


def test_the_shadowed_rule_is_shadowed_for_the_documented_reason():
    """`groups:50` is unreachable, and `groups:47` produces the same answer where it would have run --
    which is why leaving the row in the table is harmless rather than a latent behaviour difference.
    """
    for shadowed, shadow in SHADOWED.items():
        rule = next(r for r in GROUPS if r.id == shadowed)
        source, expected = rule.examples[0].split('>>')
        alone = _run(source, (rule,))
        assert alone == read_smiles(expected), (
            f'{shadowed} on its own gives {format(alone, "")!r}, not the {shadow} product '
            f'{expected!r}.  The two rules no longer agree, so the shadowing now changes the answer')
        assert next(r for r in GROUPS if r.id == shadow).query.is_substructure(read_smiles(source))


AFTER = tuple((rule, earlier) for rule in ALL for earlier in rule.after)


@pytest.mark.parametrize('rule,earlier', AFTER, ids=[f'{r.id}-after-{e}' for r, e in AFTER])
def test_every_ordering_obligation_has_a_witness(rule, earlier):
    """Some example in the collection gives a different answer when the pair is swapped.

    The pair is run alone, which isolates the claim: a whole-table swap also moves the later row past
    everything between them.  Measured across every example in both tables, because the witness for "B
    must follow A" is usually a molecule where A wins, which cannot be B's own example.
    """
    before = next(r for r in ALL if r.id == earlier)
    for _, example in EXAMPLES:
        source = example.split('>>')[0]
        molecule = read_smiles(source)
        # The screen is `or` deliberately: the commonest witness is a molecule only one of the pair can
        # match as drawn, the other becoming matchable once the first has patched it.  Requiring both
        # loses two of the ten obligations.  Sound because if neither matches, nothing fires either way.
        if not (rule.query.may_match(molecule) or before.query.may_match(molecule)):
            continue
        if _run(source, (before, rule)) != _run(source, (rule, before)):
            return
    pytest.fail(f'{rule.id} declares it must follow {earlier}, and no example in either table shows '
                f'it: the two rules give the same answer in both orders on all {len(EXAMPLES)} of '
                'them.  '
                'Either the obligation is stale, or the example that witnessed it has been edited')


def test_the_witnesses_v2_annotated_are_still_the_measured_ones():
    """The two orderings inherited as prose annotations, pinned by name in the `after` column.

    `groups:28` must follow `groups:27`, and `groups:76` must follow the `[A-]`-`[C+]` rules
    `groups:58` and `groups:65`.  A check on the measurement rather than on the code.
    """
    after = {rule.id: rule.after for rule in GROUPS}
    assert after['groups:28'] == ('groups:27',)
    assert after['groups:76'] == ('groups:58', 'groups:65')
