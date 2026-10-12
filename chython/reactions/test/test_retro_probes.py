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
"""Every row of `reactions.tsv` run backwards from its own `probe` product.

The mirror of `test_probes.py`: a row whose reverse query cannot match the product it makes, or whose
spellings rebuild something the forward row rejects, is a row `retro()` can never answer with.
"""
from collections import Counter
import pytest
from .._enumerate import retro
from .._tables import reaction_rules
from ...core import read_smiles


RULES = tuple(rule for family in reaction_rules().values() for rule in family)


def _prepared(smiles):
    molecule = read_smiles(smiles)
    molecule.kekule()
    molecule.thiele()
    return molecule


@pytest.mark.parametrize('rule', RULES, ids=lambda rule: rule.id.replace(':', '_'))
def test_every_row_reverses_its_own_probe(rule):
    _reactants, _, product = rule.probe.partition('>>')
    hits = [out for out in retro(_prepared(product), reaction=rule.name) if out.rule_id == rule.id]
    assert hits, (f'{rule.id} ({rule.name}) gives nothing back from its own probe product {product!r}.  '
                  '`compose_retro(rule.groups, rule.product)` prints what it ran.')


def test_most_rows_rebuild_their_probe_reactants_exactly():
    """A ratchet on how often a spelling IS the probe's reagent.  The rest differ by a counter-ion the
    probe draws and no spelling states, or by a reagent the spelling did not choose (NHS for NHPI)."""
    exact = 0
    for rule in RULES:
        reactants, _, product = rule.probe.partition('>>')
        want = Counter(m.canonical_bytes for m in _prepared(reactants).split())
        for out in retro(_prepared(product), reaction=rule.name):
            if out.rule_id == rule.id and Counter(
                    p.canonical_bytes for m in out.reaction.reactants for p in m.split()) == want:
                exact += 1
                break
    assert exact >= 644, f'{exact} of {len(RULES)} rows rebuild their probe reactants exactly; it was 644'
