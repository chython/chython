# -*- coding: utf-8 -*-
from chython import smiles


def test_canonicalize_nitro():
    mol = smiles('c1ccccc1N(=O)=O')
    mol.canonicalize()
    s = str(mol)
    assert '[N+]' in s and '[O-]' in s


def test_canonicalize_benzene():
    mol = smiles('C1=CC=CC=C1')
    mol.canonicalize()
    assert str(mol) == 'c1ccccc1'


def test_kekule_thiele():
    mol = smiles('c1ccccc1')
    mol.kekule()
    s = str(mol)
    assert 'c' not in s
    assert '=' in s
    mol.thiele()
    assert str(mol) == 'c1ccccc1'


def test_standardize():
    mol = smiles('c1ccccc1N(=O)=O')
    mol.standardize()
    s = str(mol)
    assert '[N+]' in s


def test_neutralize():
    mol = smiles('[NH3+]CC(=O)[O-]')
    mol.neutralize()
    s = str(mol)
    assert '+' not in s or '-' not in s


def test_check_valence_clean():
    mol = smiles('CCO')
    errors = mol.check_valence()
    assert errors == []


def test_canonicalize_idempotent():
    mol = smiles('OCC')
    mol.canonicalize()
    s1 = str(mol)
    mol.canonicalize()
    s2 = str(mol)
    assert s1 == s2
