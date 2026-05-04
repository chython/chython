# -*- coding: utf-8 -*-
from chython import smiles, MoleculeContainer, ReactionContainer, unpack


def test_molecule_pack_unpack():
    mol = smiles('c1ccccc1')
    data = mol.pack()
    mol2 = MoleculeContainer.unpack(data)
    assert str(mol) == str(mol2)


def test_molecule_pack_uncompressed():
    mol = smiles('CCO')
    data = mol.pack(compressed=False)
    mol2 = unpack(data, compressed=False)
    assert str(mol) == str(mol2)


def test_molecule_bytes():
    mol = smiles('[Cu+2]')
    data = bytes(mol)
    mol2 = unpack(data)
    assert str(mol) == str(mol2)


def test_pack_preserves_coordinates():
    mol = smiles('CCO')
    mol.clean2d()
    data = mol.pack()
    mol2 = MoleculeContainer.unpack(data)
    for n in mol:
        assert abs(mol.atom(n).x - mol2.atom(n).x) < 0.01
        assert abs(mol.atom(n).y - mol2.atom(n).y) < 0.01


def test_pack_preserves_charge():
    mol = smiles('[NH4+]')
    data = mol.pack()
    mol2 = unpack(data)
    assert str(mol) == str(mol2)


def test_pack_preserves_isotope():
    mol = smiles('[2H]C([2H])([2H])[2H]')
    data = mol.pack()
    mol2 = unpack(data)
    assert str(mol) == str(mol2)


def test_reaction_pack_unpack():
    rxn = smiles('[CH3:1][OH:2]>>[CH3:1][NH2:3]')
    data = rxn.pack()
    rxn2 = ReactionContainer.unpack(data)
    assert str(rxn) == str(rxn2)


def test_pack_len():
    mol = smiles('CCCCC')
    data = mol.pack(compressed=False)
    assert MoleculeContainer.pack_len(data, compressed=False) == 5
