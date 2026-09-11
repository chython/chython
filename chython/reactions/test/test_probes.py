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
"""Every row of `reactions.tsv` against its own `probe`.

WHY THE TABLE CARRIES ITS OWN PROBES.  A row that composes cleanly and never fires is invisible: an
absent reaction and an unmatchable one are the same empty result from `react()`.  A `z` or an `x` typed
wrong, or a product side that numbers an atom the group does not have, both land there.
"""
import pytest
from re import search
from .._tables import read_table, reaction_rules
from ...core import read_smiles


RULES = tuple(rule for family in reaction_rules().values() for rule in family)


def _prepared(smiles: str):
    """One record, in the one state a probe is read and compared in: kekulized, then aromatized."""
    molecule = read_smiles(smiles)
    molecule.kekule()
    molecule.thiele()
    return molecule


@pytest.mark.parametrize('rule', RULES, ids=lambda rule: rule.id.replace(':', '_'))
def test_every_row_yields_its_own_probe_product(rule):
    """The one gate a new row cannot pass by accident."""
    reactants, _, expected = rule.probe.partition('>>')
    got = set()
    for template in rule.templates:
        for reaction in template(_prepared(reactants)):
            for product in reaction.products:
                product.kekule()
                product.thiele()
                assert all(atom.implicit_h is not None for atom in product.atoms()), (
                    f'{rule.id} ({rule.name}) leaves an unknown hydrogen count: an all-`:` product side '
                    'states aromaticity the row cannot derive a count from.  Spell the new ring Kekule '
                    'and let `thiele()` aromatize it.')
                got.add(format(product))
    assert format(_prepared(expected)) in got, (
        f'{rule.id} ({rule.name}) does not yield {expected!r} from {reactants!r}; it yielded '
        f'{sorted(got) or "nothing"}.  `compose_smirks(rule.groups, rule.product)` prints the SMIRKS '
        'the row composes to.')


def test_every_row_has_a_probe():
    """The column is not optional, which is what makes the gate above a ratchet rather than a sample."""
    missing = [row['name'] for row in read_table('reactions.tsv') if '>>' not in row['probe']]
    assert not missing, f'reactions.tsv rows whose probe is not `<reactants>>><product>`: {missing}'


def test_a_probe_has_one_arrow_and_one_product_record():
    """Each side of `>>` is ONE record, whose components the reactant side finds wherever they are."""
    for rule in RULES:
        assert rule.probe.count('>>') == 1, f'{rule.id} ({rule.name}) probe is not `<reactants>>><product>`'
        _reactants, _, expected = rule.probe.partition('>>')
        assert '.' not in expected, (
            f'{rule.id} ({rule.name}) names {expected.count(".") + 1} product components; a row states one '
            'product, and a leaving group goes by absence rather than by being spelled')


def test_a_row_that_states_a_configuration_probes_one():
    """A `@~`, `@=` or `&<n>` whose probe cannot show it is a claim the gate above never fires.

    The token acts at the reaction centre, so the probe's own product is the only place it becomes visible:
    a substrate with no centre -- isopropanol for an inverting row, a vinyl Grignard for a retaining one --
    yields the same string with the token and without it.
    """
    unshown = [f'{rule.id} ({rule.name})' for rule in RULES
               if search(r'@~|@=|&\d', rule.product)
               and not search(r'[@/\\]', rule.probe.partition('>>')[2])]
    assert not unshown, ('reactions.tsv rows stating a configuration their probe does not carry: '
                         f'{unshown}. Draw the probe stereodefined, or drop the token.')
