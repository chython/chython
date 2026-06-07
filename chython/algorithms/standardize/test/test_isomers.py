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
from chython import smiles
from chython.algorithms.standardize._isomers import rules as isomer_rules
from pytest import mark


def _prepare(raw):
    mol = smiles(raw)
    mol.kekule()
    mol.thiele()
    return mol


def _standardized(raw):
    mol = _prepare(raw)
    mol.standardize_isomers()
    return mol


isomer_data = [
    # fixed charge rules
    ('N1C=CN2[NH+]=CC=C12', 'N1C=CC2=[NH+]C=CN12'),
    ('N(C)1C=CN2[N+](C)=CC=C12', 'N(C)1C=CC2=[N+](C)C=CN12'),
    ('N1C=C[N+]2=C1C=CN2', 'N1C=CC2=[NH+]C=CN12'),
    ('N1C=C2C=CN[N+]2=C1', 'N1C=CC2=C[NH+]=CN12'),
    ('N1C=C2NC=C[N+]2=C1', 'N1C=CN2C=[NH+]C=C12'),
    ('N1C2=[NH+]C=CC2=CC=C1', 'N1C=CC2=CC=C[NH+]=C12'),
    ('N1C=C2C(=CC=[NH+]2)C=C1', 'N1C=CC2=CC=[NH+]C=C12'),
    ('C1=CC=2C(=[NH+]1)C=CNC=2', 'N1C=CC2=C[NH+]=CC=C12'),
    ('C=1C=[NH+]C=2C=1NC=CC=2', 'N1C=CC2=[NH+]C=CC=C12'),
    ('C1=2C=[NH+]C=C1NC=CC=2', 'N1C=C2C=CC=[NH+]C2=C1'),
    ('C1=2C=CNC=C1C=[NH+]C=2', 'N1C=C2C=C[NH+]=CC2=C1'),

    # fixed charge rules with D3([A]) charge holders
    ('C12=CC=[NH+]N1C=CN2C', 'C1=C[N+](C)=C2C=CNN12'),
    ('C=1N2[N+](C)=CC=C2NC=1', 'C=1N(C)N2C(=[NH+]C=C2)C=1'),
    ('N12[N+](=CC=C1N(C=C2)C)C', 'C1=CC2=[N+](C)C=CN2N1C'),
    ('C=1N(C=2C=CN[N+]=2C=1)C', 'C1=C[N+](C)=C2C=CNN12'),
    ('C1=CC=2NC=C[N+]=2N1C', 'C=1N(C)N2C(=[NH+]C=C2)C=1'),
    ('CN1C=2C=CN(C)[N+]=2C=C1', 'C1=CC2=[N+](C)C=CN2N1C'),
    ('C=1C2=CN(C)C=[N+]2NC=1', '[N+]1(C)=CN2NC=CC2=C1'),
    ('[N+]12=CNC=C1C=CN2C', 'C1=C2C=CN(N2C=[NH+]1)C'),
    ('N1(C=C2[N+](N(C)C=C2)=C1)C', 'C1=CC=2N(N1C)C=[N+](C)C=2'),
    ('C=12[N+](=CN(C=1)C)C=CN2', 'C=1N2C(NC=1)=C[N+](C)=C2'),
    ('C=1N(C)C=2[N+](=CNC=2)C=1', 'C1=CN2C(=C[NH+]=C2)N1C'),
    ('N1(C=[N+]2C(N(C)C=C2)=C1)C', 'C1=2N(C=CN1C)C=[N+](C)C=2'),
    ('C=1N(C)C2=[NH+]C=CC2=CC=1', '[N+]=1(C)C=CC=C2C=CNC=12'),
    ('[N+]=1(C)C=CC2=CC=CNC=12', 'C=1N(C)C2=[NH+]C=CC=C2C=1'),
    ('C=1N(C)C2=[N+](C)C=CC2=CC=1', 'C=1N(C)C2=[N+](C)C=CC=C2C=1'),
    ('C1=CC=2C(=CN1C)[NH+]=CC=2', 'C=1C2=C(C=C[N+]=1C)C=CN2'),
    ('C=1C2=CC=[N+](C)C2=CNC=1', 'N1(C)C=2C=[NH+]C=CC=2C=C1'),
    ('C1=CC2=CC=[N+](C)C2=CN1C', 'N1(C)C=2C=[N+](C=CC=2C=C1)C'),
    ('C1=CN(C=C2C=C[NH+]=C12)C', 'C=1C2=C(C=CN2)C=[N+](C=1)C'),
    ('C=1[N+](C)=C2C=CNC=C2C=1', 'N1(C)C=2C=C[NH+]=CC=2C=C1'),
    ('C=1[N+](C)=C2C=CN(C=C2C=1)C', 'C=12C=CN(C)C=1C=C[N+](=C2)C'),
    ('C1=CC=C2C(=CC=[NH+]2)N1C', 'C=12C=CNC=1C=CC=[N+]2C'),
    ('C1=2[N+](C)=CC=C1NC=CC=2', 'N1(C)C2=C([NH+]=CC=C2)C=C1'),
    ('C1=CC=C2C(=CC=[N+]2C)N1C', 'C=12N(C)C=CC=1[N+](=CC=C2)C'),
    ('C1=C2N(C=CC=C2C=[NH+]1)C', 'C=1C2=CNC=C2[N+](=CC=1)C'),
    ('[N+]1(C)=CC2=CC=CNC2=C1', 'C1=2C=CC=[NH+]C1=CN(C=2)C'),
    ('C=1N(C=2C(=CC=1)C=[N+](C)C=2)C', 'C1=2C=CC=[N+](C)C1=CN(C=2)C'),
    ('[NH+]1=CC=2C(=C1)C=CN(C=2)C', 'C=1C2=CNC=C2C=C[N+]=1C'),
    ('C=1C2=C[N+](C)=CC2=CNC=1', 'C12=CN(C)C=C1C=[NH+]C=C2'),
    ('C1=C2C(=CN(C=C2)C)C=[N+]1C', 'C12=CN(C)C=C1C=[N+](C)C=C2'),

    # Morgan charge rules
    ('N1C=CC2=CC=[NH+]N12', 'N1C=CC2=CC=[NH+]N12'),
    ('N1C=CC2=[N+]1NC=C2', 'N1C=CC2=CC=[NH+]N12'),
    ('N1C=CN2C=C[NH+]=C12', 'N1C=CN2C=C[NH+]=C12'),
    ('N1C=C[N+]2=C1NC=C2', 'N1C=CN2C=C[NH+]=C12'),
    ('C=1N(C)C=[N+](CC)C=1', 'C=1N(C=[N+](C)C=1)CC'),
    ('C=1N(C=[N+](C)C=1)CC', 'C=1N(C=[N+](C)C=1)CC'),
    ('C=1N(C)C=[NH+]C=1', '[N+]1(=CNC=C1)C'),
    ('C=1N([N+](C)=CC=1)CC', 'C1=CC=[N+](CC)N1C'),
    ('C1=CC=[N+](CC)N1C', 'C1=CC=[N+](CC)N1C'),
    ('C1=CC=[N+](C)N1C', 'C1=CC=[N+](C)N1C'),

    # Morgan charge rules with D3([A]) charge holders
    ('C1=CC2=CC=[NH+]N2N1C', '[N+]=1(C)N2NC=CC2=CC=1'),
    ('[N+]=1(C)N2NC=CC2=CC=1', '[N+]=1(C)N2NC=CC2=CC=1'),
    ('C1=CC2=CC=[N+](C)N2N1C', 'C1=CC2=CC=[N+](C)N2N1C'),
    ('C1=CC2=[N+](N1)N(C)C=C2', '[N+]=1(C)N2NC=CC2=CC=1'),
    ('CN1C=CC=2C=CN(C)[N+]1=2', 'C1=CC2=CC=[N+](C)N2N1C'),
    ('C1=CN2C(=[NH+]C=C2)N1C', 'C[N+]1=C2NC=CN2C=C1'),
    ('C[N+]1=C2NC=CN2C=C1', 'C[N+]1=C2NC=CN2C=C1'),
    ('C=1N(C)C=2N(C=1)C=C[N+]=2C', 'C=1N(C)C=2N(C=1)C=C[N+]=2C'),
    ('CN1C2=[N+](C=C1)C=CN2', 'C[N+]1=C2NC=CN2C=C1'),
    ('C1=CN(C)C=2N(C)C=C[N+]1=2', 'C=1N(C)C=2N(C=1)C=C[N+]=2C'),
    ('N1C=[NH+]C=C1', 'N1C=[NH+]C=C1'),
    ('N1[NH+]=CC=C1', 'C=1C=[NH+]NC=1'),
    ('N(C)1[NH+]=CC=C1', 'C1=CC=[NH+]N1C'),
    ('N1[N+](C)=CC=C1', 'C1=CC=[NH+]N1C'),

    # ferrocene
    ('[CH-]1C=CC=C1.[Fe+2].[CH-]1C=CC=C1', '[CH-]1C=CC=C1.[Fe+2].[CH-]1C=CC=C1'),
    ('[CH-]1C=CC=C1.[Fe+2].C1=C[CH-]C=C1', '[CH-]1C=CC=C1.[Fe+2].[CH-]1C=CC=C1'),

    # fixed tautomer rules
    ('N1C=NC=N1', 'N1C=NN=C1'),
    ('CC1=NC=NN1', 'N1=CNC(C)=N1'),
    ('CC1=NN=CN1', 'N1=CNC(C)=N1'),
    ('N1N=CN=N1', 'N1C=NN=N1'),
    ('CC1=NNN=N1', 'C=1(C)NN=NN=1'),
    ('CC1=NN=NN1', 'C=1(C)NN=NN=1'),

    # Morgan tautomer rules
    ('CC1=NNC=C1', 'N1C=CC(C)=N1'),
    ('CC1=CC=NN1', 'N1C=CC(C)=N1'),
    ('CC1=CN=CN1', 'C=1N=CNC=1C'),
    ('CC1=CNC=N1', 'C=1N=CNC=1C'),
    ('CC1=CN=NN1', 'N1N=NC(C)=C1'),
    ('CC1=CNN=N1', 'N1N=NC(C)=C1'),
    ('CC1=NNN=C1', 'N1N=NC(C)=C1'),
    ('CC1=NC2=NNC(C)=C2N1', 'C1(C)=NC=2C(=NNC=2N1)C'),
    ('CC1=NC2=C(N1)C(C)=NN2', 'C1(C)=NC=2C(=NNC=2N1)C'),
    ('CC1=NC2=C(C)NN=C2N1', 'C1(C)=NC=2C(=NNC=2N1)C'),

    # amidine, guanidine, etc.
    ('COC(=N)NC', 'COC(N)=NC'),
    ('COC(N)=NC', 'COC(N)=NC'),
    ('CCN=C(N)NC', 'CCNC(=NC)N'),
    ('CCNC(=N)NC', 'CCNC(=NC)N'),
    ('CCNC(N)=NC', 'CCNC(=NC)N'),
    ('CNC(N)=NC(=N)NC', 'CNC(=N)NC(=N)NC'),
    ('CNC(=N)NC(=N)NC', 'CNC(=N)NC(=N)NC'),
    ('CNC(N)=NC(N)=NC', 'CNC(=N)NC(=N)NC'),
    ('CCN=CNC=NC', 'CCN=CN=CNC'),
]


@mark.parametrize('raw,result', isomer_data)
def test_standardize_isomers(raw, result):
    assert _standardized(raw) == _prepare(result), f'{raw} > {_standardized(raw)} != {result}'


changed_data = [
    ('N(C)1C=CN2[N+](C)=CC=C12', True),
    ('N1C=C2C=CN[N+]2=C1', True),
    ('N(C)1C=C2C=CN[N+]2=C1', True),
    ('N1C=C2C=CN(C)[N+]2=C1', True),
    ('N(C)1C=C2C=CN(C)[N+]2=C1', True),
    ('C=1N(C)C=[N+](CC)C=1', True),
    ('C=1N(C=[N+](C)C=1)CC', False),
    ('C=1N([N+](C)=CC=1)CC', True),
    ('C1=CC=[N+](CC)N1C', False),
    ('C1=CC=[N+](C)N1C', False),
    ('[N+]1(=CNC=C1)C', False),
    ('C[N+]1=CC=NC=C1', False),
]


@mark.parametrize('raw,changed', changed_data)
def test_standardize_isomers_changed(raw, changed):
    mol = _prepare(raw)
    assert bool(mol.standardize_isomers()) is changed


idempotent_data = [
    raw for raw, _ in isomer_data
] + [
    raw for raw, _ in changed_data
]


@mark.parametrize('raw', idempotent_data)
def test_standardize_isomers_idempotent(raw):
    mol = _prepare(raw)
    mol.standardize_isomers()
    first = str(mol)
    mol.standardize_isomers()
    assert str(mol) == first


def test_all_isomer_rules_have_representative():
    mols = [(raw, _prepare(raw)) for raw in idempotent_data]
    missing = []
    for index, (query, rule_type) in enumerate(isomer_rules, 1):
        if not any(next(query.get_mapping(mol, automorphism_filter=False), None) for _, mol in mols):
            missing.append((index, rule_type, str(query)))

    assert not missing
