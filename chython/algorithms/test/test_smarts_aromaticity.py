# -*- coding: utf-8 -*-
#
#  Copyright 2026 Tagir Akhmetshin <tagirshin@gmail.com>
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
"""Regression tests for case-sensitive aromaticity in SMARTS atom symbols.

Per standard SMARTS semantics:
- ``[C]`` (uppercase) matches **aliphatic** carbon only;
- ``[c]`` (lowercase) matches **aromatic** carbon only.

Historically chython treated the case of the symbol as just a default
hybridization hint, so a rule with ``[C;D3]`` would silently match a
``[c;D3]`` in the target — producing chemically incorrect retros
("methylate an aromatic phenol as if it were an aliphatic alcohol").
"""
from chython import smiles, smarts
from chython.reactor import Reactor


def test_smarts_uppercase_C_does_not_match_aromatic_c():
    # A bare [C] query must not match an aromatic carbon (`c` in
    # phenol).
    q = smarts("[C]")
    target = smiles("c1ccc(O)cc1")
    aromatic_carbons = [n for n, a in target.atoms() if a.atomic_symbol == "C" and a.hybridization == 4]
    assert aromatic_carbons, "test fixture must contain aromatic carbons"
    for m in q.get_mapping(target):
        n = next(iter(m.values()))
        assert target.atom(n).hybridization != 4, (
            f"[C] aliphatic query matched aromatic atom {n}: hybridization={target.atom(n).hybridization}"
        )


def test_smarts_lowercase_c_does_not_match_aliphatic_C():
    # And conversely: lowercase [c] must not match an aliphatic carbon.
    q = smarts("[c]")
    target = smiles("CC")
    mappings = list(q.get_mapping(target))
    assert mappings == [], (
        f"[c] aromatic query must not match aliphatic CC, got {mappings}"
    )


def test_reactor_aliphatic_C_does_not_match_aromatic_phenol():
    # End-to-end: a rule that methylates an aliphatic D3-C-OH must not
    # fire on aromatic phenol-like OH.
    rule_rxn = smarts("[C;D3:1]-[O;D1:2]>>[C;D3:1]-[O;D2:2]-[C;D1:3]")
    reactor = Reactor(
        patterns=tuple(rule_rxn.reactants),
        products=tuple(rule_rxn.products),
        delete_atoms=False,
    )
    # The OH sits on an aromatic carbon. Per SMARTS semantics [C;D3]
    # must not match the aromatic carbon, so the reactor must yield
    # nothing.
    target = smiles("c1ccc(O)cc1")
    assert list(reactor(target)) == []
    # Direct mapping check too — same constraint at the matcher level.
    assert list(rule_rxn.reactants[0].get_mapping(target)) == []


def test_smarts_uppercase_C_still_matches_aliphatic_C():
    # The positive control: [C;D3] in the rule must still match an
    # actual aliphatic D3 carbon. No false negatives from the fix.
    rule_rxn = smarts("[C;D3:1]-[O;D1:2]>>[C;D3:1]-[O;D2:2]-[C;D1:3]")
    target = smiles("CC(C)O")  # isopropanol — aliphatic D3-C with OH
    mappings = list(rule_rxn.reactants[0].get_mapping(target))
    assert len(mappings) >= 1, (
        "aliphatic [C;D3]-[O;D1] must still match isopropanol"
    )
