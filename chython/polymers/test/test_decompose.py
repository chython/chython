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
"""A molecule into the monomer graph: cut at library attachment bonds, fewest monomers, stereo exact."""

import json
from time import perf_counter

from pytest import mark, raises

from chython.core import smiles

from ._fixtures import OXYTOCIN_HELM, assemble_text, rules
from .._decompose import decompose
from .._errors import StereoMismatch
from .._library import get_monomer_library, set_monomer_library


def run(mol, log=None):
    return decompose(mol, get_monomer_library(), [] if log is None else log)


def labels(graph):
    """Each edge as ``{(symbol, R), (symbol, R)}``, orientation-free."""
    return sorted(sorted(((graph.nodes[e.a].symbol, e.ra), (graph.nodes[e.b].symbol, e.rb))) for e in graph.edges)


def test_a_peptide_decomposes_to_its_monomers():
    g = run(smiles('C[C@H](N)C(=O)NCC(=O)O'))
    assert sorted(n.symbol for n in g.nodes) == ['A', 'G']
    assert labels(g) == [[('A', 2), ('G', 1)]]


@mark.parametrize('residue', ['N', 'Q'])
def test_a_side_chain_amide_is_not_cut(residue):
    g = run(assemble_text(f'PEPTIDE1{{A.{residue}.G}}$$$$'))
    assert sorted(n.symbol for n in g.nodes) == sorted(['A', residue, 'G'])


def test_a_disulfide_is_a_cross_link():
    g = run(assemble_text(OXYTOCIN_HELM))
    assert sorted(n.symbol for n in g.nodes) == sorted('CYIQNCPLG') + ['am']
    assert [('C', 3), ('C', 3)] in labels(g)


def test_a_head_to_tail_ring():
    g = run(assemble_text('PEPTIDE1{G.A.G}$PEPTIDE1,PEPTIDE1,3:R2-1:R1$$$'))
    assert len(g.nodes) == 3 and len(g.edges) == 3


def test_an_undefined_centre_raises():
    with raises(StereoMismatch, match='A'):
        run(smiles('CC(N)C(=O)NCC(=O)O'))


def test_the_wrong_configuration_with_no_monomer_for_it_raises(tmp_path):
    # the public set holds every Thr, Ile and Hyp diastereomer, so attach one diastereomer of an amino
    # acid it lacks -- a tert-butyl threonine analogue -- and draw another
    p = tmp_path / 'one.json'
    p.write_text(json.dumps([{'symbol': 'tBuThr', 'polymerType': 'PEPTIDE', 'monomerType': 'Backbone',
                              'smiles': 'CC(C)(C)[C@@H](O)[C@H](N[H:1])C([OH:2])=O'}]))
    set_monomer_library(p)
    try:
        assert sorted(n.symbol for n in run(smiles('CC(C)(C)[C@@H](O)[C@H](N)C(=O)NCC(=O)O')).nodes) == \
            ['G', 'tBuThr']
        with raises(StereoMismatch, match='tBuThr'):
            run(smiles('CC(C)(C)[C@H](O)[C@H](N)C(=O)NCC(=O)O'))
    finally:
        set_monomer_library()


def test_a_region_matching_nothing_is_one_inline_monomer():
    # no public monomer attaches at an ester oxygen, so the glycyl ester is one region the library cannot
    # split; the alanine still matches
    log = []
    g = run(smiles('C[C@H](N)C(=O)NCC(=O)OCC1CCCCC1'), log)
    (inline,) = [n for n in g.nodes if n.inline]
    assert '[*:1]' in inline.symbol
    assert [n.symbol for n in g.nodes if not n.inline] == ['A']
    assert 'polymers:inline-monomer' in rules(log)


def test_node_order_does_not_depend_on_atom_order():
    a = run(smiles('C[C@H](N)C(=O)NCC(=O)N[C@@H](CC(C)C)C(=O)O'))
    b = run(smiles('OC(=O)[C@H](CC(C)C)NC(=O)CNC(=O)[C@@H](N)C'))
    assert [n.symbol for n in a.nodes] == [n.symbol for n in b.nodes]
    assert labels(a) == labels(b)


#: A click macrocycle's monomers, public chemistry: ethynylglycine (R3 cap on its stereocentre), 5-azido-
#: norvaline (R3 azide cap) and the 1,4-triazole the click leaves (R1 on carbon, R4 on nitrogen).
CLICK = [{'symbol': 'Egly', 'polymerType': 'PEPTIDE', 'smiles': '[H:1]N[C@@H]([C:3](#C))C([OH:2])=O'},
         {'symbol': 'Anv', 'polymerType': 'PEPTIDE', 'smiles': '[H:1]N[C@@H](CCC[N:3](=[N+]=[N-]))C([OH:2])=O'},
         {'symbol': 'tz', 'polymerType': 'PEPTIDE', 'smiles': '[H:1]c1cn([H:4])nn1'}]


def test_a_click_macrocycle_round_trips(tmp_path):
    # an R3 on an alpha carbon and an alkyl R3 must not license cutting side-chain C-C bonds, nor the
    # triazole's carbon every benzyl bond: F stays F, Anv stays Anv
    p = tmp_path / 'click.json'
    p.write_text(json.dumps(CLICK))
    set_monomer_library(p)
    try:
        text = 'PEPTIDE1{[Egly].F.G.[Anv]}|PEPTIDE2{[tz]}$PEPTIDE1,PEPTIDE2,1:R3-1:R1|PEPTIDE1,PEPTIDE2,4:R3-1:R4$$$'
        g = run(assemble_text(text))
        assert sorted(n.symbol for n in g.nodes) == ['Anv', 'Egly', 'F', 'G', 'tz']
    finally:
        set_monomer_library()


def test_a_small_ring_is_not_cut(tmp_path):
    # with an alkyl R3 in the library, Tic's ring CH2-N bond looks like a junction; a ring of six is not one
    p = tmp_path / 'click.json'
    p.write_text(json.dumps(CLICK))
    set_monomer_library(p)
    try:
        g = run(assemble_text('PEPTIDE1{A.[Tic].G}$$$$'))
        assert sorted(n.symbol for n in g.nodes) == ['A', 'G', 'Tic']
    finally:
        set_monomer_library()


def test_a_long_peptide_with_unknown_regions_decomposes_quickly():
    # unmatched pieces along a long chain made the covering search exponential when it was one problem
    spacer = '[[*:1]C(=O)CCOCCN[*:2]]'
    m = assemble_text('PEPTIDE1{' + '.'.join(['A', 'C', 'D', 'E', 'F', spacer, 'G', 'H', 'K', 'L'] * 4) + '}$$$$')
    start = perf_counter()
    g = run(m)
    assert sum(n.inline for n in g.nodes) == 4 and perf_counter() - start < 2.0
