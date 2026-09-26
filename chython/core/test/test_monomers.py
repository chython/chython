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
"""`monomers()`: residues, main chain, branches and cross-links of public peptides."""
from dataclasses import FrozenInstanceError

import pytest

from chython.core import smiles
from chython.core.monomers import monomers
from .peptides import CASES, PEPTIDES, explicit, molecule


@pytest.mark.parametrize('name, residues, main, links, h2t', [
    ('leu-enkephalin', 5, 5, 0, False), ('oxytocin', 9, 9, 1, False), ('linaclotide', 14, 14, 3, False),
    ('gramicidin S', 10, 10, 1, True), ('cyclosporin A', 11, 11, 1, True), ('semaglutide', 35, 31, 0, False),
    ('beta/gamma mix', 5, 5, 0, False), ('lactam', 7, 7, 1, False)])
def test_segmentation(name, residues, main, links, h2t):
    p = monomers(molecule(name))
    assert (len(p.residues), len(p.main), len(p.crosslinks), p.head_to_tail) == (residues, main, links, h2t)


def test_kinds_of_the_beta_gamma_mix():
    p = monomers(molecule('beta/gamma mix'))
    assert [p.residues[i].kind for i in p.main] == ['beta', 'beta', 'gamma', 'alpha', 'beta']


@pytest.mark.parametrize('name', PEPTIDES + ('branched', 'lactam', 'Ac3c', 'side-chain branch', 'triazole staple',
                                             'alkene staple'))
def test_backbone_is_a_bonded_path(name):
    mol = molecule(name)
    p = monomers(mol)
    assert len(set(p.backbone)) == len(p.backbone)
    assert all(b in mol.neighbors_of(a) for a, b in zip(p.backbone, p.backbone[1:]))
    for branch in p.branches:
        assert branch.path[0] in mol.neighbors_of(branch.host)
        assert all(b in mol.neighbors_of(a) for a, b in zip(branch.path, branch.path[1:]))


@pytest.mark.parametrize('name', ['triazole staple', 'alkene staple'])
def test_staple_links_on_a_single_bond(name):
    """A small ring or a double bond is claimed whole, so the cross-link falls on a single bond beside it."""
    mol = molecule(name)
    (a, b), = monomers(mol).crosslinks
    assert mol.order_of(a, b) == 1
    assert not any(len(r) <= 7 and a in r and b in r for r in mol.sssr)


def test_side_chain_branch_path_runs_to_the_next_residue():
    mol = molecule('side-chain branch')
    br, = monomers(mol).branches
    assert len(br.chain) == 4 and not any(mol.element_of(x) == 8 for x in br.path)


def test_explicit_hydrogens_leave_the_backbone():
    bb = monomers(molecule('leu-enkephalin')).backbone
    assert monomers(explicit(molecule('leu-enkephalin'))).backbone == bb


def test_head_to_tail_link_is_last():
    mol = molecule('gramicidin S')
    p = monomers(mol)
    c, n = p.crosslinks[-1]
    assert (c, n) == (p.residues[p.main[-1]].exit, p.residues[p.main[0]].entry)
    assert n == p.backbone[0] and c == p.backbone[-1]


def test_branches():
    p = monomers(molecule('semaglutide'))
    assert len(p.branches) == 1 and len(p.branches[0].chain) == 4
    p = monomers(molecule('branched'))
    assert len(p.branches) == 1 and len(p.branches[0].chain) == 3


def test_members_partition_the_molecule():
    mol = molecule('semaglutide')
    p = monomers(mol)
    covered = [x for i in (*p.main, *(j for b in p.branches for j in b.chain)) for x in p.members[i]]
    assert len(covered) == len(set(covered)) == len(mol)


@pytest.mark.parametrize('smi', [CASES['tripeptide'], 'c1ccccc1', 'CC(=O)Oc1ccccc1C(=O)O', 'CC(=O)NC'])
def test_not_a_peptide(smi):
    assert monomers(smiles(smi)) is None


def test_threshold():
    assert monomers(molecule('tetrapeptide')) is not None


def test_library_assigns_symbols_without_moving_boundaries():
    mol = molecule('leu-enkephalin')
    library = {'Phe': smiles('NC(Cc1ccccc1)C(=O)O'), 'Gly': smiles('NCC(=O)O')}
    a, b = monomers(mol), monomers(mol, library=library)
    assert [r.atoms for r in a.residues] == [r.atoms for r in b.residues]
    # Tyr is the N-terminus and Leu the C-terminus: neither is in the library
    assert [b.residues[i].symbol for i in b.main] == [None, 'Gly', 'Gly', 'Phe', None]


def test_frozen():
    p = monomers(molecule('oxytocin'))
    with pytest.raises(FrozenInstanceError):
        p.main = ()
    with pytest.raises(FrozenInstanceError):
        p.residues[0].kind = 'other'
