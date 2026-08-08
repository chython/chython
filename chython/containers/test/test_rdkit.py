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


enhanced = [
    # source, expected extended_stereo per stereocenter in chython atom order
    ('C[C@H](O)[C@@H](N)[C@H](F)C |&1:1,o2:3|', [1, -2, None]),
    ('C[C@H](O)[C@@H](N)C |&1:1,&1:3|', [1, 1]),
    ('C[C@H](O)[C@@H](N)C |&3:1,o2:3|', [3, -2]),
    ('C[C@H](O)CC |o1:1|', [-1]),
    ('C[C@H](O)[C@@H](N)C', [None, None]),
]


@mark.parametrize('source, expected', enhanced)
def test_enhanced_stereo_roundtrip(source, expected):
    """AND/OR groups survive the round trip, group ids included."""
    mol = smiles(source)
    assert [a.extended_stereo for _, a in mol.atoms() if a.stereo is not None] == expected

    back = from_rdkit_molecule(mol.to_rdkit(keep_mapping=False))
    assert format(back, 'h') == format(mol, 'h')
    assert [a.extended_stereo for _, a in back.atoms() if a.stereo is not None] == expected


def test_enhanced_stereo_export():
    mol = smiles('C[C@H](O)[C@@H](N)[C@H](F)C |&1:1,o2:3|')

    # absolute group is implicit by default and explicit on demand, mirroring the MDL V3000 writers
    assert Chem.MolToCXSmiles(mol.to_rdkit(keep_mapping=False)).endswith('|o2:3,&1:5|')
    assert Chem.MolToCXSmiles(mol.to_rdkit(keep_mapping=False, absolute=True)).endswith('|a:1,o2:3,&1:5|')

    collection = [x[7:] for x in Chem.MolToV3KMolBlock(mol.to_rdkit(keep_mapping=False, absolute=True)).splitlines()
                  if 'MDLV30/STE' in x]
    assert collection == ['MDLV30/STEABS ATOMS=(1 6)', 'MDLV30/STERAC1 ATOMS=(1 2)', 'MDLV30/STEREL2 ATOMS=(1 4)']

    # a molecule without extended stereo gets no groups unless asked for
    assert not smiles('C[C@H](O)CC').to_rdkit(keep_mapping=False).GetStereoGroups()


def test_enhanced_stereo_absolute_dropped_on_import():
    """chython marks a center absolute by the absence of a group, so ABS carries no information."""
    mol = from_rdkit_molecule(Chem.MolFromSmiles('C[C@H](O)[C@@H](N)C |a:1,&1:3|'))
    assert [a.extended_stereo for _, a in mol.atoms() if a.stereo is not None] == [None, 1]


def test_enhanced_stereo_only_tetrahedrons():
    """A group on a non-tetrahedral or non-stereogenic atom must not reach rdkit.

    rdkit happily writes a collection entry for an atom without a chiral tag, which would
    produce a mol block claiming stereo the molecule does not have.
    """
    mol = smiles('CC(O)CC')
    mol.atom(2)._extended_stereo = 1  # no stereo label -> masked by the extended_stereo property
    assert not mol.to_rdkit(keep_mapping=False).GetStereoGroups()

    # allene centers carry extended stereo too, but allenes are not exported at all
    allene = smiles('OC(F)=[C@]=C(O)F |&1:3|')
    assert allene.atom(4).extended_stereo == 1  # not masked, the center is labeled
    assert not allene.to_rdkit(keep_mapping=False).GetStereoGroups()


def test_enhanced_stereo_atom_in_two_groups():
    """rdkit allows an atom in an AND and an OR group at once, the signed chython int does not."""
    mol = from_rdkit_molecule(Chem.MolFromSmiles('C[C@H](O)[C@@H](N)C |&1:1,o1:1|'))
    assert [a.extended_stereo for _, a in mol.atoms() if a.stereo is not None] == [1, None]


@mark.parametrize('source', ['C/C=C/C.O', 'C/C=C/C.C', 'CC/C=C\\CC.[Na+].[Cl-]',
                             'C/C=C/C.C/C=C\\C', 'C/C=C/C=C/C.CC(=O)O'])
def test_cis_trans_multifragment(source):
    """rdkit derives smiles bond directions from the stereo atoms only for single fragment mols.

    Without an explicit SetDoubleBondNeighborDirections every salt and solvate silently
    lost its cis-trans marks on smiles export.
    """
    assert Chem.MolToSmiles(smiles(source).to_rdkit(keep_mapping=False)) == Chem.CanonSmiles(source)


def test_cis_trans_at_boron():
    """The dative bond direction fix reversed every bond, breaking SetStereoAtoms on B and At."""
    assert Chem.MolToSmiles(smiles('C/B=C/C').to_rdkit(keep_mapping=False)) == 'C/B=C/C'


@mark.parametrize('order', [('Fe', 'N'), ('N', 'Fe')])
def test_dative_bond_direction(order):
    """A dative bond always points from the donor to the acceptor, whatever the atom order."""
    from chython import MoleculeContainer

    mol = MoleculeContainer()
    n, m = (mol.add_atom(x) for x in order)
    mol.add_bond(n, m, 8)
    rd_mol = mol.to_rdkit(keep_mapping=False)

    bond = rd_mol.GetBondWithIdx(0)
    assert (bond.GetBeginAtom().GetSymbol(), bond.GetEndAtom().GetSymbol()) == ('N', 'Fe')


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
