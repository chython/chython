# -*- coding: utf-8 -*-
#
#  Copyright 2025 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from chython import smiles, from_rdkit_molecule
from pytest import mark
from rdkit import Chem
from rdkit.Chem import AllChem


data = [
    'CCO',
    'C/C=C/C',
    'C[C@H](O)CC',
    'C\C=C/O[C@@H]1OC[C@@H](Oc2ccccc2)[C@@H](O)[C@H]1O\C=C\C',
    '[nH]1cccc1',
    'C\C=C\C=C',
    'C[C@@H](O)[C@H](O)[C@H](C)O'
]

@mark.parametrize('source', data)
def test_to_rdkit(source):
    mol = smiles(source)
    rd_mol = mol.to_rdkit(keep_mapping=False)
    rd_mol_mapping = mol.to_rdkit(keep_mapping=True)

    assert format(smiles(Chem.MolToSmiles(rd_mol)), 'h') == format(mol, 'h')
    assert format(smiles(Chem.MolToSmiles(rd_mol_mapping)), 'm') == format(mol, 'm')


@mark.parametrize('source', data)
def test_from_rdkit(source):
    assert format(from_rdkit_molecule(Chem.MolFromSmiles(source)), 'h') == format(smiles(source), 'h')


def test_coordinates():
    rd_mol = smiles('CCO').to_rdkit(keep_mapping=False)

    AllChem.Compute2DCoords(rd_mol)
    mol = from_rdkit_molecule(rd_mol)
    assert any(a.x for _, a in mol.atoms())

    rd_mol_h = Chem.AddHs(rd_mol)
    AllChem.EmbedMolecule(rd_mol_h)
    rd_mol_nh = Chem.RemoveHs(rd_mol_h)

    mol = from_rdkit_molecule(rd_mol_nh)
    assert hasattr(mol, '_conformers')
    assert isinstance(mol._conformers, list)
    assert len(mol._conformers) == 1
    assert isinstance(mol._conformers[0], dict)
    assert len(mol._conformers[0]) == 3
    assert all(tuple(x) for x in mol._conformers[0].values())
    assert all(len(x) == 3 for x in mol._conformers[0].values())
    assert all(isinstance(x, float) for x in mol._conformers[0].values() for x in x)
    assert any(x for x in mol._conformers[0].values() for x in x)
