# -*- coding: utf-8 -*-
from chython import xyz


def test_xyz_water():
    mol = xyz((('O', 0., 0., 0.), ('H', 1., 0., 0.), ('H', 0., 1., 0.)))
    assert len(mol) == 3
    assert 'O' in {a.atomic_symbol for _, a in mol.atoms()}


def test_xyz_methane():
    mol = xyz((
        ('C', 0., 0., 0.),
        ('H', 1.09, 0., 0.),
        ('H', -0.36, 1.03, 0.),
        ('H', -0.36, -0.51, 0.89),
        ('H', -0.36, -0.51, -0.89),
    ))
    assert len(mol) == 5
    assert mol.atoms_count == 5
