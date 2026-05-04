# -*- coding: utf-8 -*-
from chython import smiles


def test_reaction_parse():
    rxn = smiles('[CH3:1][OH:2]>>[CH3:1][NH2:3]')
    assert len(rxn.reactants) == 1
    assert len(rxn.products) == 1


def test_reaction_components():
    rxn = smiles('CCO.CC(=O)O>>CCOC(=O)C.O')
    assert len(rxn.reactants) == 2
    assert len(rxn.products) == 2


def test_reaction_equality():
    rxn1 = smiles('[CH3:1][OH:2]>>[CH3:1][NH2:3]')
    rxn2 = smiles('[CH3:1][OH:2]>>[CH3:1][NH2:3]')
    assert rxn1 == rxn2


def test_reaction_canonicalize():
    rxn = smiles('CCO.CC(=O)O>>CCOC(=O)C.O')
    rxn.canonicalize()
    assert str(rxn)


def test_reaction_molecules():
    rxn = smiles('[CH3:1][OH:2]>>[CH3:1][NH2:3]')
    mols = list(rxn.molecules())
    assert len(mols) == 2


def test_reaction_pickle():
    from pickle import loads, dumps
    rxn = smiles('[CH3:1][OH:2]>>[CH3:1][NH2:3]')
    data = dumps(rxn)
    rxn2 = loads(data)
    assert str(rxn) == str(rxn2)
