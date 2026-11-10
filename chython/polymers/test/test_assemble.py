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
"""The monomer graph assembled into a molecule through one SMILES parse."""

from pytest import raises

from chython.core import smiles

from ._fixtures import OXYTOCIN_HELM, OXYTOCIN_SMILES, assemble_text, rules, same
from .._errors import HelmError


def test_a_dipeptide_is_the_molecule():
    assert assemble_text('PEPTIDE1{A.G}$$$$') == smiles('C[C@H](N)C(=O)NCC(=O)O')


def test_stereo_is_kept():
    assert assemble_text('PEPTIDE1{[dA].G}$$$$') == smiles('C[C@@H](N)C(=O)NCC(=O)O')


def test_proline_closes_its_own_ring():
    assert assemble_text('PEPTIDE1{A.P}$$$$') == smiles('C[C@H](N)C(=O)N1CCC[C@H]1C(=O)O')


def test_a_disulfide_ring():
    assert same(assemble_text(OXYTOCIN_HELM), smiles(OXYTOCIN_SMILES))


def test_a_head_to_tail_ring():
    assert assemble_text('PEPTIDE1{G.G.G}$PEPTIDE1,PEPTIDE1,3:R2-1:R1$$$') == smiles('C1C(=O)NCC(=O)NCC(=O)N1')


def test_no_explicit_hydrogen_is_left():
    m = assemble_text('PEPTIDE1{C.A}$$$$')
    assert all(m.atom(a).atomic_symbol != 'H' for a in m.atom_numbers)
    assert not any(m.atom(a).map_number for a in m.atom_numbers)


def test_an_inline_monomer():
    assert assemble_text('PEPTIDE1{A.[[*:1]CC([*:2])=O].G}$$$$') == \
        smiles('C[C@H](N)C(=O)CC(=O)NCC(=O)O')


def test_an_unknown_symbol_is_a_marker_with_its_alias():
    log = []
    m = assemble_text('PEPTIDE1{A.[NoSuch].G}$$$$', log)
    (marker,) = [a for a in m.atom_numbers if m.atom(a).is_r]
    assert m.aliases[marker] in ('NoSuch', b'NoSuch') and len(m.neighbors_of(marker)) == 2
    assert 'polymers:monomer-unknown' in rules(log)


def test_an_attachment_the_monomer_lacks_is_refused():
    with raises(HelmError, match='R3'):
        assemble_text('PEPTIDE1{A.G}|PEPTIDE2{G}$PEPTIDE1,PEPTIDE2,1:R3-1:R1$$$')


def test_more_than_ninety_bonds_assemble():
    m = assemble_text('PEPTIDE1{' + '.'.join('A' * 120) + '}$$$$')
    assert m.atom_count == 120 * 5 + 1          # 120 Ala residues plus the terminal OH oxygen
