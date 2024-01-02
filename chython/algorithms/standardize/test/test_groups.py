# -*- coding: utf-8 -*-
#
#  Copyright 2023 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from pytest import mark


data = [
        ('CP(C)(C)C', 'C[P+](C)(C)C'),
        ('CB1(C)[H]B(C)(C)[H]1', 'CB1(C)~[H]B(C)(C)~[H]1'),
        ('[O]N(C)[NH]', '[O-][N+](C)=N'), ('[O]N(C)[CH2]', '[O-][N+](C)=C'),
        ('[O]S(C)(C)[O]', 'O=S(C)(C)=O'), ('[O]S(C)(C)[S]', 'O=S(C)(C)=S'),
        ('BN(C)=C', 'B~N(C)=C'),
        ('B=N(C)(C)C', 'B~N(C)(C)C'),
        ('BS(C)C', 'B~S(C)C'), ('BO(C)C', 'B~O(C)C'),
        ('[B-]=[N+](C)C', 'BN(C)C'), ('C[B-]=[N+]C', 'CBNC'), ('[B-]=[N+]', 'BN'),
        ('[O-][B+3]([O-])([O-])[O-]', 'O[B-](O)(O)O'),
        ('[O-]B(O)(O)O', 'O[B-](O)(O)O'),
        ('OB(O)(O)O', 'O[B-](O)(O)O'),
        ('CN(C)(C)C', 'C[N+](C)(C)C'),
        ('C=N(=O)O', 'C[N+](=O)[O-]'),
        ('C=N(=O)C', 'C=[N+]([O-])C'), ('O=N(=O)C', 'O=[N+]([O-])C'), ('N=N(=O)C', 'N=[N+]([O-])C'),
        ('C=[N+]([O-])O', 'C[N+](=O)[O-]'),
        ('CN(=O)=N(=O)C', 'C[N+]([O-])=[N+]([O-])C'),
        ('CN(=O)=N(=N)C', 'C[N+]([O-])=[N+]([NH-])C'), ('CN(=O)=N(=NC)C', 'C[N+]([O-])=[N+]([N-]C)C'),
        ('C=N(=N)C', 'C=[N+]([N-])C'), ('C=N(=NC)C', 'C=[N+]([N-]C)C'), ('N=N(=N)C', 'N=[N+]([N-])C'),
        ('[N-][N+](=O)C', 'N=[N+]([O-])C'), ('C[N-][N+](=O)C', 'CN=[N+]([O-])C'),
        ('[O-]N(=O)=O', '[O-][N+](=O)[O-]'),
        ('CN(:O):O', 'C[N+](=O)[O-]'),
        ('O=[N-]=O', '[O-]N=O'),
        ('O=N#N', 'O=[N+]=[N-]'), ('C=N#N', 'C=[N+]=[N-]'), ('N=N#N', 'N=[N+]=[N-]'),
        ('[O-][N+]#N', 'O=[N+]=[N-]'), ('C[CH-][N+]#N', 'CC=[N+]=[N-]'), ('[NH-][N+]#N', 'N=[N+]=[N-]'),
        ('C[N+]#N=[N-]', 'CN=[N+]=[N-]'),
        ('CN=N=N', 'CN=[N+]=[N-]'),
        ('CNN#N', 'CN=[N+]=[N-]'),
        ('[N-]#N=NC', '[N-]=[N+]=NC'),
        ('[N-]=N#N', '[N-]=[N+]=[N-]'),
        ('CC#N=N', 'CC=[N+]=[N-]'),
        ('CC#N=NC', 'CC#[N+][N-]C'), ('CC#N=O', 'CC#[N+][O-]'),
        ('NN#C', '[NH-][N+]#C'), ('ON#CC', '[O-][N+]#CC'), ('SN#CC', '[S-][N+]#CC'),
        ('CNN#C', 'C[N-][N+]#C'),
        ('N[N+]#[C-]', '[NH-][N+]#C'), ('O[N+]#[C-]', '[O-][N+]#C'), ('S[N+]#[C-]', '[S-][N+]#C'),
        ('CN[N+]#[C-]', 'C[N-][N+]#C'),
        ('CN#C', 'C[N+]#[C-]'),
        ('C[C-]=[N+]=N', 'CC=[N+]=[N-]'),
        ('CN(C)(C)=O', 'C[N+](C)(C)[O-]'), ('CN(C)(C)=NC', 'C[N+](C)(C)[N-]C'),
        ('C[N](C)=O |^1:1|', 'CN(C)[O] |^1:3|'),
        ('C=N(C)[O] |^1:3|', 'C=[N+](C)[O-]'),
        ('CN(C)#N', 'C[N+](C)=[N-]'),
        ('CN=[N+]', 'C[N+]#N'),
        ('[NH2+][C-]=O', 'N=C=O'), ('C[NH+][C-]=O', 'CN=C=O'),
        ('N#CO', 'N=C=O'),
        ('N#C[O-]', '[N-]=C=O'),
        ('CC(C)(C)[N+][O-]', 'CC(C)(C)N=O'),
        ('CN=O', 'C=NO'),
        ('NC=[O+]C', '[NH2+]=COC'), ('CNC=[O+]C', 'C[NH+]=COC'), ('CN(C)C=[O+]C', 'C[N+](C)=COC'),
        ('N=CO', 'NC=O'), ('N=CS', 'NC=S'),
        ('O=C1NC=CC=C1', 'OC1=NC=CC=C1'), ('OC1=NC=CC=C1', 'OC1=NC=CC=C1'), ('N=C1NC=CC=C1', 'NC1=NC=CC=C1'),
        ('CN=C1NC=CC=C1', 'CNC1=NC=CC=C1'),
        ('O=C1C=CNC=C1', 'OC1=CC=NC=C1'), ('OC1=CC=NC=C1', 'OC1=CC=NC=C1'),
        ('C=C(O)O', 'CC(=O)O'), ('C=C(O)N', 'CC(=O)N'),
        ('OC=C', 'O=CC'), ('OC(C)=C', 'O=C(C)C'),
        ('O=C1N=CC=CC1', 'OC=1N=CC=CC=1'), ('OC=1N=CC=CC=1', 'OC=1N=CC=CC=1'), ('N=C1N=CC=CC1', 'NC=1N=CC=CC=1'),
        ('CN=C1N=CC=CC1', 'CNC=1N=CC=CC=1'),
        ('[O-][P+](C)(C)C', 'O=P(C)(C)C'),
        ('[CH2+][P-](C)(C)C', 'C=P(C)(C)C'),
        ('FP(F)(F)(F)(F)F', 'F[P-](F)(F)(F)(F)F'),
        ('C[P-](C)(=O)=O', 'CP(C)(=O)[O-]'),
        ('CP(C)(=O)[S-]', 'CP(C)(=S)[O-]'),
        ('CP(C)(=O)S', 'CP(C)(=S)O'),
        ('CS(=O)(=O)[S-]', 'CS(=O)(=S)[O-]'),
        ('CS(=O)(=O)S', 'CS(=O)(=S)O'),
        ('C[S+](C)[O-]', 'CS(C)=O'),
        ('O=[S+][O-]', 'O=S=O'), ('O=[S+](C)(C)[O-]', 'O=S(C)(C)=O'),
        ('O=[S+2]([O-])[O-]', 'O=S(=O)=O'),
        ('C[S+2](C)([O-])[O-]', 'CS(C)(=O)=O'),
        ('C[S-](C)[CH2+]', 'CS(C)=C'),
        ('O=[S-](C)=O', 'O=S(C)[O-]'),
        ('O=[S-](C)(=O)=O', 'O=S(C)(=O)[O-]'),
        ('O=[S-](=O)[S-]', 'S=S([O-])[O-]'),
        ('CS(C)(=O)O', 'CS(C)(=O)=O'),
        ('CS(C)(=O)N', 'CS(C)(=O)=N'),
        ('N=S(C)O', 'NS(C)=O'), ('N=S(C)(C)(C)O', 'NS(C)(C)(C)=O'),
        ('C#CO', 'C=C=O'),
        ('C#CNC', 'C=C=NC'),
        ('C=O |^1:0|', '[C-]#[O+]'),
        ('[O]O[O] |^1:0,2|', 'O=[O+][O-]'),
        ('[CH2+]N(C)C', 'C=[N+](C)C'),
        ('[CH2+]N(C)O', 'C=[N+](C)O'),
        ('[CH2+]=NC', 'C#[N+]C'),
        ('O[Cl+][O-]', 'OCl=O'),
        ('O[Cl+2]([O-])[O-]', 'OCl(=O)=O'),
        ('O[Cl+3]([O-])([O-])[O-]', 'OCl(=O)(=O)=O'),
        ('[Cl-]=O', 'Cl[O-]'),
        ('OS(=N)(=N)O', 'O=S(N)(N)=O'), ('OS(=N)(=N)C', 'O=S(N)(=N)C')
]


@mark.parametrize('raw,result', data)
def test_group(raw, result):
    tmp = smiles(raw)
    tmp.standardize()
    assert tmp == smiles(result), f'{raw} > {tmp} != {result}'
