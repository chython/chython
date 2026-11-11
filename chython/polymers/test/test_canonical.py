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
"""The canonical writers: one molecule, one string, whatever its atom order or HELM spelling."""

from pytest import mark, raises

from chython.core import smiles

from ._fixtures import OXYTOCIN_HELM, assemble_text, same
from .._biln import parse_biln
from .._decompose import decompose
from .._errors import BilnError
from .._library import get_monomer_library
from .._write import write_biln, write_helm


def to_helm(mol):
    return write_helm(decompose(mol, get_monomer_library(), []))


def to_biln(mol):
    return write_biln(decompose(mol, get_monomer_library(), []))


def reordered(mol):
    """The same molecule with its atoms in a random order: a random SMILES read back."""
    return smiles(format(mol, 'r'))


def test_a_linear_peptide_is_written_n_to_c():
    assert to_helm(smiles('C[C@H](N)C(=O)NCC(=O)O')) == 'PEPTIDE1{A.G}$$$$'


@mark.parametrize('text', [OXYTOCIN_HELM, 'PEPTIDE1{A.N.G.[Aib].P}$$$$',
                           'PEPTIDE1{G.A.L.G}$PEPTIDE1,PEPTIDE1,4:R2-1:R1$$$'])
def test_atom_order_does_not_change_the_string(text):
    m = assemble_text(text)
    want = to_helm(m)
    for _ in range(20):
        assert to_helm(reordered(m)) == want


def test_an_inline_monomer_is_spelt_canonically():
    # a piece no monomer matches; whatever order its atoms arrive in, its inline SMILES is one text
    m = smiles('C[C@H](N)C(=O)NCC(=O)OCC1CCC(C(=O)O)CC1')
    want = to_helm(m)
    assert '[*:1]' in want
    for _ in range(20):
        assert to_helm(reordered(m)) == want


def test_a_biaryl_bond_spelt_aromatic_is_a_single_bond():
    # OpenSMILES leaves the biphenyl bond unmarked, which reads aromatic; it is still the Bip side chain
    biaryl = smiles('C[C@H](N)C(=O)N[C@@H](Cc1ccc(cc1)c1ccccc1)C(=O)NCC(=O)O')
    assert any(b.order == 4 and not b.in_ring for b in biaryl.bonds())
    assert to_helm(biaryl) == to_helm(assemble_text('PEPTIDE1{A.[Bip].G}$$$$')) == 'PEPTIDE1{A.[Bip].G}$$$$'


def test_the_string_reads_back_to_the_molecule():
    m = assemble_text(OXYTOCIN_HELM)
    assert same(assemble_text(to_helm(m)), m)


def test_two_spellings_write_one_string():
    a = 'PEPTIDE1{A.G.G.L}$PEPTIDE1,PEPTIDE1,4:R2-1:R1$$$'
    b = 'PEPTIDE1{G.L.A.G}$PEPTIDE1,PEPTIDE1,4:R2-1:R1$$$'      # the same ring, opened elsewhere
    assert to_helm(assemble_text(a)) == to_helm(assemble_text(b))


def test_chain_order_does_not_change_the_string():
    a = 'PEPTIDE1{A.C}|PEPTIDE2{G.C}$PEPTIDE1,PEPTIDE2,2:R3-2:R3$$$'
    b = 'PEPTIDE1{G.C}|PEPTIDE2{A.C}$PEPTIDE2,PEPTIDE1,2:R3-2:R3$$$'
    assert to_helm(assemble_text(a)) == to_helm(assemble_text(b))


def test_a_cross_link_is_a_connection():
    helm = to_helm(assemble_text(OXYTOCIN_HELM))
    chains, connections = helm.split('$')[:2]
    assert chains.startswith('PEPTIDE1{') and connections.count(':R3') == 2


def test_biln_is_written_and_reads_back():
    m = assemble_text(OXYTOCIN_HELM)
    text = to_biln(m)
    assert text.count('(1,3)') == 2
    graph = parse_biln(text)
    assert sorted(n.symbol for n in graph.nodes) == sorted('CYIQNCPLG') + ['am']


def test_biln_refuses_what_it_cannot_spell():
    with raises(BilnError):
        to_biln(smiles('C[C@H](N)C(=O)NCC(=O)OCC1CCCCC1'))
