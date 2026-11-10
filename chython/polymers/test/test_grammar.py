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
"""HELM and BILN strings into the monomer graph.  Edges compare without orientation: `a-b` is `b-a`."""

from pytest import mark, raises

from .._biln import parse_biln
from .._errors import BilnError, HelmError
from .._helm import parse_helm


def edges(graph):
    return {frozenset(((e.a, e.ra), (e.b, e.rb))) for e in graph.edges}


def e(a, ra, b, rb):
    return frozenset(((a, ra), (b, rb)))


def symbols(graph):
    return [n.symbol for n in graph.nodes]


def rules(graph):
    return {r.rule for r in graph.log}


def test_a_linear_peptide():
    g = parse_helm('PEPTIDE1{A.G.[Aib]}$$$$')
    assert symbols(g) == ['A', 'G', 'Aib']
    assert [(n.chain, n.position, n.polymer_type) for n in g.nodes] == [
        ('PEPTIDE1', 1, 'PEPTIDE'), ('PEPTIDE1', 2, 'PEPTIDE'), ('PEPTIDE1', 3, 'PEPTIDE')]
    assert edges(g) == {e(0, 2, 1, 1), e(1, 2, 2, 1)}


def test_a_connection_between_chains():
    g = parse_helm('PEPTIDE1{A.C}|PEPTIDE2{G}$PEPTIDE1,PEPTIDE2,2:R3-1:R1$$$')
    assert symbols(g) == ['A', 'C', 'G']
    assert edges(g) == {e(0, 2, 1, 1), e(1, 3, 2, 1)}


def test_a_head_to_tail_ring():
    g = parse_helm('PEPTIDE1{A.G.G.L}$PEPTIDE1,PEPTIDE1,4:R2-1:R1$$$')
    assert edges(g) == {e(0, 2, 1, 1), e(1, 2, 2, 1), e(2, 2, 3, 1), e(3, 2, 0, 1)}


def test_an_inline_smiles_monomer():
    g = parse_helm('PEPTIDE1{A.[C([*:1])C([*:2])=O]}$$$$')
    assert [n.inline for n in g.nodes] == [False, True]
    assert g.nodes[1].symbol == 'C([*:1])C([*:2])=O'


def test_nucleotides_are_sugar_base_phosphate():
    g = parse_helm('RNA1{R(A)P.R(U)}$$$$')
    assert symbols(g) == ['R', 'A', 'P', 'R', 'U']
    assert [n.position for n in g.nodes] == [1, 2, 3, 4, 5]
    assert edges(g) == {e(0, 3, 1, 1), e(0, 2, 2, 1), e(2, 2, 3, 1), e(3, 3, 4, 1)}


def test_an_exact_repeat_expands():
    g = parse_helm("PEPTIDE1{A'3'}$$$$")
    assert symbols(g) == ['A', 'A', 'A'] and edges(g) == {e(0, 2, 1, 1), e(1, 2, 2, 1)}


def test_a_repeat_range_is_opaque_and_logged():
    g = parse_helm("PEPTIDE1{A'2-5'}$$$$")
    assert len(g.opaque) == 1 and g.opaque[0][1] == "A'2-5'"
    assert 'polymers:opaque-unit' in rules(g)


@mark.parametrize('unit', ['(A,G)', '(A+G)', '?', '*'])
def test_a_unit_that_names_no_one_monomer_is_opaque(unit):
    g = parse_helm(f'PEPTIDE1{{A.{unit}.G}}$$$$')
    assert [i for i, _ in g.opaque] == [1] and len(g.nodes) == 3


def test_biln_is_the_same_graph():
    b = parse_biln('A-G-C(1,3).C(1,3)')
    h = parse_helm('PEPTIDE1{A.G.C}|PEPTIDE2{C}$PEPTIDE1,PEPTIDE2,3:R3-1:R3$$$$')
    assert symbols(b) == symbols(h) and edges(b) == edges(h)
    assert [n.chain for n in b.nodes] == [n.chain for n in h.nodes]


def test_biln_takes_bracketed_names_and_several_bonds():
    g = parse_biln('[Aib](1,1)-G-C(1,2)')
    assert symbols(g) == ['Aib', 'G', 'C'] and edges(g) == {e(0, 2, 1, 1), e(1, 2, 2, 1), e(0, 1, 2, 2)}


def test_a_hydrogen_pair_is_not_a_bond():
    g = parse_helm('RNA1{R(A)P}|RNA2{R(U)P}$RNA1,RNA2,2:pair-2:pair$$$')
    assert len(g.edges) == 4 and 'polymers:hydrogen-pair-not-built' in rules(g)


def test_annotations_are_logged():
    g = parse_helm('PEPTIDE1{A.G}"heavy chain"$$${"note":1}$')
    assert symbols(g) == ['A', 'G'] and 'polymers:annotation-not-kept' in rules(g)


@mark.parametrize('text', ['PEPTIDE1{A.G$$$$', 'PEPTIDE1{A}$PEPTIDE1,PEPTIDE9,1:R3-1:R1$$$',
                           'PEPTIDE1{A}$PEPTIDE1,PEPTIDE1,5:R3-1:R1$$$', 'NOPE1{A}$$$$', ''])
def test_a_malformed_helm_is_refused(text):
    with raises(HelmError):
        parse_helm(text)


@mark.parametrize('text', ['A-G-C(1,3)', 'A-(G', ''])
def test_a_malformed_biln_is_refused(text):
    with raises(BilnError):
        parse_biln(text)
