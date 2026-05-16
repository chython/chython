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
"""Regression tests for aromaticity preservation through ``Reactor``.

Applying a SMARTS rule whose pattern uses aromatic-bond (``:``) tokens
on atoms inside a fused aromatic ring system must not silently
dearomatize neighbouring rings that the rule does not touch.

The pre-fix behaviour was to emit a non-aromatizable bond-order
assignment that even ``thiele()`` could not rescue.
"""
from chython import smiles, smarts
from chython.reactor import Reactor


def _aromatic_atom_count(mol) -> int:
    return sum(1 for _, a in mol.atoms() if a.hybridization == 4)


def _run_first_product(rule_str: str, target_smiles: str):
    rule = smarts(rule_str)
    target = smiles(target_smiles)
    reactor = Reactor(
        patterns=tuple(rule.reactants),
        products=tuple(rule.products),
        delete_atoms=False,
    )
    results = list(reactor(target))
    if not results:
        return target, None
    return target, results[0].products[0]


def test_reactor_preserves_adjacent_aromaticity_quinolinone():
    # Rule oxidizes CH2OH on a quinolinone-type ring to CHO. The fused
    # aromatic system on the quinolinone parent ring must remain
    # aromatic — the touched atoms are aromatic carbons and an
    # aromatic nitrogen.
    rule_str = (
        "[C:1]:[C;D3:2](-[C;D2:3]-[O;D1:4]):[N:5]"
        ">>"
        "[C:1]:[C;D3:2](-[C;D2:3]=[O;D1:4]):[N:5]"
    )
    target_smi = "Cc1noc(C)c1-c2cc3c(c(=O)cc(CO)n3-c4ccc(cc4F)F)cc2"
    target, product = _run_first_product(rule_str, target_smi)
    assert product is not None, "reactor produced no output on quinolinone fixture"
    target_aromatic = _aromatic_atom_count(target)
    product_aromatic = _aromatic_atom_count(product)
    # The rule preserves the count of aromatic atoms (no aromatic atom
    # leaves the aromatic system, the only chemistry is CH2OH -> CHO).
    assert product_aromatic >= target_aromatic, (
        f"aromaticity degraded by reactor: in={target_aromatic}, "
        f"out={product_aromatic} on {product}"
    )


def test_reactor_output_aromaticity_stable_under_thiele():
    # A correct reactor output must already be in a re-aromatized state
    # (or at minimum, a state from which ``thiele()`` can recover the
    # full set of aromatic atoms).
    rule_str = (
        "[C:1]:[C;D3:2](-[C;D2:3]-[O;D1:4]):[N:5]"
        ">>"
        "[C:1]:[C;D3:2](-[C;D2:3]=[O;D1:4]):[N:5]"
    )
    target_smi = "Cc1noc(C)c1-c2cc3c(c(=O)cc(CO)n3-c4ccc(cc4F)F)cc2"
    _, product = _run_first_product(rule_str, target_smi)
    assert product is not None
    before = _aromatic_atom_count(product)
    product.thiele()
    after = _aromatic_atom_count(product)
    assert after >= before, (
        f"thiele lost aromatic atoms: before={before}, after={after}"
    )
    # And the recovered count must reach the input's aromatic count
    # (the rule doesn't touch aromatic atoms).
    target_aromatic = _aromatic_atom_count(smiles(target_smi))
    assert after == target_aromatic, (
        f"reactor output cannot reach full aromaticity even after thiele: "
        f"target_aromatic={target_aromatic}, output_aromatic={after}"
    )


def test_reactor_preserves_aromaticity_purine_like_ring_opening():
    # imidazo[4,5-d]pyrimidine-style bicycle: rule opens the imidazole
    # portion. The pyrimidine ring portion is chemically valid as an
    # aromatic heterocycle in the product; the benzene side chain must
    # also remain aromatic.
    rule_str = (
        "[N;D3:1]:1:[C:2]:[N:3]:[C;D3:4]:[C;D3:5]:2:[C;D3:6]:1"
        ":[N;D2:7]:[C;D2:8]:[N;D2:9]:2"
        ">>"
        "[C;D2:8](=[O;D1:10])-[O;D1:11]."
        "[N;D3:1]:1:[C:2]:[N:3]:[C;D3:4]:[C;D3:5](:[C;D3:6]:1-[N;D1:7])-[N;D1:9]"
    )
    target_smi = "c12[nH]cnc1n(c(=S)[nH]c2=O)-c3cc(CC)ccc3"
    target, product = _run_first_product(rule_str, target_smi)
    assert product is not None, (
        "rule must fire on the purine-like target — otherwise this test "
        "passes vacuously and does not actually check aromaticity"
    )
    # The benzene side chain (6 aromatic carbons) MUST remain aromatic in
    # any chemically valid output even after the imidazole ring is opened.
    # The pyrimidine ring should also remain aromatic. Total: at least 12.
    benzene_atoms_out = _aromatic_atom_count(product)
    assert benzene_atoms_out >= 6, (
        f"benzene side chain dearomatized: aromatic atoms in product = "
        f"{benzene_atoms_out}, expected >= 6. Product: {product}"
    )
