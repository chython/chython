# -*- coding: utf-8 -*-
#
#  Copyright 2025, 2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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
    # xy of a 3d conformer is a projection, not a layout. it must not be imported
    assert not any(a.x or a.y for _, a in from_rdkit_molecule(rd_mol_nh).atoms())

    mol = from_rdkit_molecule(rd_mol_nh)
    assert mol._conformers is not None
    assert isinstance(mol._conformers, list)
    assert len(mol._conformers) == 1
    assert isinstance(mol._conformers[0], dict)
    assert len(mol._conformers[0]) == 3
    assert all(tuple(x) for x in mol._conformers[0].values())
    assert all(len(x) == 3 for x in mol._conformers[0].values())
    assert all(isinstance(x, float) for x in mol._conformers[0].values() for x in x)
    assert any(x for x in mol._conformers[0].values() for x in x)


def test_no_layout_no_conformer():
    """A molecule without a 2d layout must not export a degenerate all-zero conformer.

    RDKit dumps it into CXSMILES and derives bond directions from it, losing stereo.
    """
    mol = smiles('C/C=C/CO')
    rd_mol = mol.to_rdkit(keep_mapping=False)
    assert rd_mol.GetNumConformers() == 0
    assert '|' not in Chem.MolToCXSmiles(rd_mol)
    # cis-trans survives the mol block round trip only without the degenerate conformer
    reparsed = Chem.MolFromMolBlock(Chem.MolToMolBlock(rd_mol))
    assert Chem.MolToSmiles(reparsed) == Chem.MolToSmiles(rd_mol)

    mol.clean2d()
    rd_mol = mol.to_rdkit(keep_mapping=False)
    assert rd_mol.GetNumConformers() == 1
    assert not rd_mol.GetConformer(0).Is3D()

    assert smiles('CCO').to_rdkit(keep_coordinates=True).GetNumConformers() == 1
    assert mol.to_rdkit(keep_coordinates=False).GetNumConformers() == 0


def test_conformers_invariant():
    """Conformers either match the atom set exactly or are None. Never desynced.

    Anything that changes the atom set or the connectivity flushes them; only
    implicify_hydrogens keeps geometry, since hiding an explicit H does not move
    the heavy atoms.
    """
    def build():
        rd = Chem.AddHs(Chem.MolFromSmiles('CC(C)(C)OC(=O)NCC(=O)O'))
        AllChem.EmbedMolecule(rd, randomSeed=42)
        mol = from_rdkit_molecule(Chem.RemoveHs(rd))
        assert mol._conformers[0].keys() == set(mol)
        return mol

    # copy and remap preserve the geometry, keys stay in sync
    mol = build()
    copy = mol.copy()
    assert copy._conformers[0] == mol._conformers[0]
    assert copy._conformers[0] is not mol._conformers[0]
    mol.remap({next(iter(mol)): max(mol) + 10})
    assert mol._conformers[0].keys() == set(mol)

    # explicify_hydrogens adds atoms with unknown coordinates -> flush
    mol = build()
    assert mol.explicify_hydrogens()
    assert mol._conformers is None

    # implicify_hydrogens is the one exception: geometry kept, removed H pruned
    mol = smiles('[H]CCO')
    mol.generate_conformers(limit=1)
    assert mol.implicify_hydrogens()
    assert mol._conformers[0].keys() == set(mol)

    # any edit of the atom set or the connectivity flushes
    for edit in (lambda m: m.delete_atom(next(iter(m))),
                 lambda m: m.add_atom('O'),
                 lambda m: m.delete_bond(*next((n, k) for n, k, _ in m.bonds())),
                 lambda m: m.add_bond(next(iter(m)), m.add_atom('O'), 1),
                 lambda m: m.remove_protection(),
                 lambda m: m.union(smiles('CC'), remap=True, copy=False)):
        mol = build()
        edit(mol)
        assert mol._conformers is None, edit

    # substructure / split drop geometry: bonds are cut
    mol = build()
    assert mol.substructure(list(mol)[:3])._conformers is None
    assert mol.union(smiles('CC'), remap=True)._conformers is None
    mol = build()
    assert all(m._conformers is None for m in mol.split())

    # a rolled back transaction restores the geometry it flushed
    mol = build()
    saved = mol._conformers[0].copy()
    try:
        with mol as m:
            m.delete_atom(next(iter(m)))
            assert m._conformers is None
            raise ValueError
    except ValueError:
        pass
    assert mol._conformers[0] == saved
