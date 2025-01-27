# -*- coding: utf-8 -*-
#
#  Copyright 2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from itertools import zip_longest


def test_keto_enol_2h_pyrrole():
    """
    2H‐pyrrole. [N:1]=1[C:2][C:3]=,:[C:4][C:5]=1
    """
    for t, v in zip(['C1C=CC=N1', 'C1N=CC2=C1C=CC=C2'], ['C=1C=CNC=1', 'N1C=C2C=CC=CC2=C1']):
        s = smiles(t)
        t = set(s.enumerate_tautomers())
        v = smiles(v)
        s.thiele()
        v.thiele()
        assert len(t) == 2, ' '.join(str(x) for x in t)
        assert t == {s, v}, f'{", ".join(str(x) for x in t)} != {s}, {v}'


def test_acid_protonated_nitrogen():
    """
    Ammonia, H-imine,guanidine,amidine. [H][N+]=,-,:
    """
    for t, v in zip_longest(['C=[NH2+].[Cl-]', 'C=[NH+]C.[Cl-]', 'C[NH3+].[Cl-]', 'C[NH2+]C.[Cl-]', 'C[NH1+](C)C.[Cl-]',
                             '[NH4+].[Cl-]', 'C=[N+](C)C.[Cl-]', 'C[N+](C)(C)C.[Cl-]'],
                            ['C=N.Cl', 'C=NC.Cl', 'CN.Cl', 'CNC.Cl', 'CN(C)C.Cl', 'N.Cl']):
        s = smiles(t)
        t = set(s.enumerate_tautomers())
        if v:
            assert len(t) == 2, ' '.join(str(x) for x in t)
            assert t == {s, smiles(v)}
        else:
            assert len(t) == 1


def test_base_nitrogen():
    for t, v in zip_longest(['N1C=CN=N1.Cl', 'NC(N)=N.Cl', 'CN(C)C(=NO)N(C)C.Cl', 'CN(C)C(=NC)N(C)C.Cl',
                             'CN(C)C(=NN)N(C)C.Cl', 'COC(N)=N.Cl', 'CSC(N)=N.Cl', 'COC(OC)=N.Cl', 'COC(C)=N.Cl',
                             'CNN.Cl', 'CN.Cl',
                             'N.Cl', 'N1C=CC=C1.Cl'],
                            [('N1N=NC=C1.Cl', 'N1=CC=NN1.Cl', '[NH+]=1NC=CN=1.[Cl-]', 'N1=[NH+]C=CN1.[Cl-]',
                              'N1=CC=[NH+]N1.[Cl-]'),
                             ('NC(N)=N.Cl', 'NC(N)=[NH2+].[Cl-]'),
                             ('CN(C)C(=NO)N(C)C.Cl', 'CN(C)C(N(C)C)=[NH+]O.[Cl-]'),
                             ('CN(C)C(=NC)N(C)C.Cl', 'CN(C)C(N(C)C)=[NH+]C.[Cl-]'),
                             ('CN(C)C(=NN)N(C)C.Cl', 'CN(C)C(N(C)C)=[NH+]N.[Cl-]'),
                             ('Cl.NC(=N)OC', '[NH2+]=C(N)OC.[Cl-]'), ('Cl.NC(=N)SC', '[NH2+]=C(N)SC.[Cl-]'),
                             ('COC(OC)=N.Cl', 'COC(OC)=[NH2+].[Cl-]'),
                             ('COC(C)=N.Cl', 'COC(C)=[NH2+].[Cl-]'),
                             ('CNN.Cl', 'CN[NH3+].[Cl-]', 'C[NH2+]N.[Cl-]'),
                             ('CN.Cl', 'C[NH3+].[Cl-]')]):
        s = smiles(t)
        t = set(s.enumerate_tautomers(zwitter=True))
        if v:
            assert len(t) == len(v), ' '.join(str(x) for x in t)
            vs = set()
            for x in v:
                x = smiles(x)
                x.thiele(fix_tautomers=False)
                vs.add(x)
            assert t == vs, ' '.join(str(x) for x in t) + ' != ' + ' '.join(str(x) for x in vs)
        else:
            assert len(t) == 1, ' '.join(str(x) for x in t)
