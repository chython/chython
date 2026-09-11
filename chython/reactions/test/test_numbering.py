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
from ...core import ReactionContainer, read_smiles as smiles
from .._numbering import fast_mapping, mapping_agrees


def test_fast_mapping_is_an_isomorphism():
    a = smiles('c1ccc(CO)cc1')
    b = smiles('OCc1ccccc1')
    pairs = fast_mapping(a, b)
    assert pairs is not None
    assert len(pairs) == len(a.atom_numbers)
    assert sorted(pairs.values()) == sorted(b.atom_numbers)
    for n, m in pairs.items():
        assert a.atom(n).element == b.atom(m).element        # `.element` is the atomic NUMBER
        assert a.atom(n).implicit_h == b.atom(m).implicit_h
    for bond in a.bonds():
        assert b.bond(pairs[bond.n], pairs[bond.m]).order == bond.order


def test_fast_mapping_refuses_different_structures():
    assert fast_mapping(smiles('CCO'), smiles('CCC')) is None


def test_fast_mapping_is_stereo_aware():
    # meso versus (R,R) tartaric acid: one graph, two structures.  A graph-only correspondence
    # would answer here; `canonical_order()` must not.
    assert fast_mapping(smiles('O[C@@H](C(=O)O)[C@H](O)C(=O)O'),
                        smiles('O[C@@H](C(=O)O)[C@@H](O)C(=O)O')) is None


def test_mapping_agrees_on_itself():
    rxn = _mapped()
    # _mapped() has two products (brominated toluene + HBr); mapping_agrees scores all products,
    # so the expected count is the total atom count across all products, not just products[0].
    total = sum(len(p.atom_numbers) for p in rxn.products)
    assert mapping_agrees(rxn, rxn) == (total, 0, 0)


def test_mapping_agrees_excuses_an_automorphic_swap():
    # The two ortho carbons of toluene are one orbit, so swapping their numbers is the same answer.
    reference = _mapped()
    produced = reference.copy()
    product = produced.products[0]
    ring = [n for n in product.atom_numbers if product.atom(n).hybridization == 4]
    a, b = ring[1], ring[-1]
    na, nb = product.map_number_of(a), product.map_number_of(b)
    product.set_map_number(a, nb)
    product.set_map_number(b, na)
    agreed, disagreed, missing = mapping_agrees(produced, reference)
    assert disagreed == 0 and missing == 0


def test_mapping_agrees_counts_missing_for_orphaned_map_number():
    # Orphan a product atom by assigning it a map number that no input carries.
    # The atom is absent from `got`, so it must land in `missing`, not `disagreed`.
    reference = _mapped()
    produced = reference.copy()
    product = produced.products[0]
    n = next(iter(product.atom_numbers))
    product.set_map_number(n, 9999)
    agreed, disagreed, missing = mapping_agrees(produced, reference)
    assert missing == 1


def test_mapping_agrees_counts_a_real_disagreement():
    # Swap the map numbers of two reactant atoms that are NOT in the same automorphism orbit
    # (methyl C, sp3, map 1, versus ipso ring C, sp2, map 2 — unambiguously distinct).
    # Every product atom's key survives in both got and want because every map number is still
    # present on some input atom, so nothing is missing; the two swapped atoms disagree.
    reference = _mapped()
    produced = reference.copy()
    reactant = produced.reactants[0]   # toluene [CH3:1][c:2]1[cH:3][cH:4][cH:5][cH:6][cH:7]1
    n1 = next(n for n in reactant.atom_numbers if reactant.map_number_of(n) == 1)
    n2 = next(n for n in reactant.atom_numbers if reactant.map_number_of(n) == 2)
    reactant.set_map_number(n1, 2)
    reactant.set_map_number(n2, 1)
    agreed, disagreed, missing = mapping_agrees(produced, reference)
    total = sum(len(p.atom_numbers) for p in reference.products)  # 9
    assert missing == 0
    assert disagreed == 2
    assert agreed + disagreed == total


def _mapped():
    """A fully mapped two-product reaction: toluene bromination, giving 4-bromotoluene and HBr,
    written with map numbers 1-9."""
    return ReactionContainer([smiles('[CH3:1][c:2]1[cH:3][cH:4][cH:5][cH:6][cH:7]1'),
                              smiles('[Br:8][Br:9]')],
                             [smiles('[CH3:1][c:2]1[cH:3][cH:4][c:5]([Br:8])[cH:6][cH:7]1'),
                              smiles('[BrH:9]')])
