# -*- coding: utf-8 -*-
#
#  Copyright 2023-2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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
        ('O=C1NC=CC=C1', 'OC1=NC=CC=C1'), ('OC1=NC=CC=C1', 'OC1=NC=CC=C1'),
        ('N=C1NC=CC=C1', 'NC1=NC=CC=C1'),
        ('CN=C1NC=CC=C1', 'CNC1=NC=CC=C1'),
        ('O=C1C=CNC=C1', 'OC1=CC=NC=C1'), ('OC1=CC=NC=C1', 'OC1=CC=NC=C1'),
        ('C=C(O)O', 'CC(=O)O'), ('C=C(O)N', 'CC(=O)N'),
        ('OC=C', 'O=CC'), ('OC(C)=C', 'O=C(C)C'),
        ('O=C1N=CC=CC1', 'OC=1N=CC=CC=1'), ('OC=1N=CC=CC=1', 'OC=1N=CC=CC=1'),
        ('N=C1N=CC=CC1', 'NC=1N=CC=CC=1'),
        ('CN1N=CNC1=O', 'CN1N=CN=C1O'), ('CN1N=CN=C1O', 'CN1N=CN=C1O'),
        ('S=C1NC=CN1', 'SC1=NC=CN1'), ('SC1=NC=CN1', 'SC1=NC=CN1'),
        ('S=C1NCCN1', 'S=C1NCCN1'), ('SC1=NCCN1', 'S=C1NCCN1'),
        ('CN=C1NC=CN1', 'CNC1=NC=CN1'), ('CNC1=NC=CN1', 'CNC1=NC=CN1'),
        ('CN=C1NCCN1', 'CNC1=NCCN1'), ('CNC1=NCCN1', 'CNC1=NCCN1'),
        ('S=C1NNC=C1', 'SC1=NNC=C1'), ('SC1=NNC=C1', 'SC1=NNC=C1'),
        ('S=C1NNCC1', 'S=C1NNCC1'), ('SC1=NNCC1', 'S=C1NNCC1'),
        ('CN=C1NNC=C1', 'CNC1=NNC=C1'), ('CNC1=NNC=C1', 'CNC1=NNC=C1'),
        ('CN=C1NNCC1', 'CNC1=NNCC1'), ('CNC1=NNCC1', 'CNC1=NNCC1'),
        ('N=C1NC=CO1', 'NC1=NC=CO1'), ('NC1=NC=CO1', 'NC1=NC=CO1'),
        ('OC1=NC=CO1', 'OC1=NC=CO1'), ('O=C1NC=CO1', 'OC1=NC=CO1'),
        ('O=C1NCCO1', 'O=C1NCCO1'), ('OC1=NCCO1', 'O=C1NCCO1'),
        ('CN=C1NOC=C1', 'CNC1=NOC=C1'),
        ('CN=C1N=CC=CC1', 'CNC=1N=CC=CC=1'),
        ('OC1=CC=NN1', 'OC1=CC=NN1'), ('OC1=CN=CO1', 'OC1=CN=CO1'),
        ('OC1=NN=CO1', 'OC1=NN=CO1'), ('OC1=NC=CS1', 'OC1=NC=CS1'),
        ('OC1=NOC=C1', 'OC1=NOC=C1'), ('CN1N=CC=C1O', 'CN1N=CC=C1O'),
        ('OC1=CC=CN1', 'OC1=CC=CN1'), ('O=C1CC=CN1', 'OC1=CC=CN1'),
        ('SC1=CC=CN1', 'SC1=CC=CN1'), ('S=C1CC=CN1', 'SC1=CC=CN1'),
        ('OC1=CC=CS1', 'OC1=CC=CS1'), ('O=C1CC=CS1', 'OC1=CC=CS1'),
        ('OC1=CSC=C1', 'OC1=CSC=C1'), ('O=C1CSC=C1', 'OC1=CSC=C1'),
        ('OC1=CC=CO1', 'OC1=CC=CO1'), ('O=C1CC=CO1', 'OC1=CC=CO1'),
        ('OC1=COC=C1', 'OC1=COC=C1'), ('O=C1COC=C1', 'OC1=COC=C1'),
        # ('OC1=NC(O)=NC=C1', 'OC1=NC(O)=NC=C1'),
        # ('NC1=CN=C(O)N=C1O', 'NC1=CN=C(O)N=C1O'),
        # ('OC1=NC2=CC=CC=C2C(O)=N1', 'OC1=NC2=CC=CC=C2C(O)=N1'),
        ('OC1=CC2=C(N1)C=NC=C2', 'OC1=CC2=C(N1)C=NC=C2'),
        ('OC1=CN2C=CC=CC2=C1', 'OC1=CN2C=CC=CC2=C1'),
        ('OC1=CN2C=CN=CC2=N1', 'OC1=CN2C=CN=CC2=N1'),
        ('OC1=NC2=NN=CN2C=C1', 'OC1=NC2=NN=CN2C=C1'),
        ('OC1=CC=CC2=CC=CN12', 'OC1=CC=CC2=CC=CN12'),
        ('O=C1CN2C=CC=C2C=C1', 'OC1=CN2C=CC=C2C=C1'),
        ('OC1=CN2C=CC=C2C=C1', 'OC1=CN2C=CC=C2C=C1'),
        ('OC1=NC=CN2C=CC=C12', 'OC1=NC=CN2C=CC=C12'),
        ('OC1=NC=NN2C=CN=C12', 'OC1=NC=NN2C=CN=C12'),
        ('SC1=NC2=CC=CC=C2O1', 'SC1=NC2=CC=CC=C2O1'),
        ('SC1=NC2=C(O1)C=CC=C2', 'SC1=NC2=C(O1)C=CC=C2'),
        ('CN1C(O)=NC2=CC=CC=C12', 'CN1C(O)=NC2=CC=CC=C12'),
        ('CN1C(O)=NC2=C(N)N=CN=C12', 'CN1C(O)=NC2=C(N)N=CN=C12'),
        ('CN1C(O)=NC2=C(N)N=C(N)N=C12', 'CN1C(O)=NC2=C(N)N=C(N)N=C12'),
        ('OC1=CSC2=CC=CN12', 'OC1=CSC2=CC=CN12'),
        ('O=C1CSC2=CC=CN12', 'OC1=CSC2=CC=CN12'),
        ('OC1=COC2=CC=CN12', 'OC1=COC2=CC=CN12'),
        ('O=C1COC2=CC=CN12', 'OC1=COC2=CC=CN12'),
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


tautomer_fix_data = [
    ('O=C1NC(=O)NC(=O)N1', 'N1=C(N=C(N=C1O)O)O'),  # 1,3,5-triketone (cyanuric acid)
    ('O=C1NC(=O)C=CN1', 'OC1=NC=CC(O)=N1'),  # pyridin-1,3-dione (uracil)
    ('O=C1C(=O)NC=CN1', 'C=1(C(=NC=CN=1)O)O'),  # pyridin-1,2-dione
    ('O=C1C=CC(=O)NN1', 'C1=CC(O)=NN=C1O'),  # pyridin-1,4-dione
    ('O=C1NC(=O)CC1', 'C=1C=C(O)NC=1O'),  # pyrrole-2,5-dione (succinimide)
    ('O=C1NCC(=O)C1', 'C=1(O)C=C(O)NC=1'),  # pyrrole-2,4-dione
    ('O=C1CC=CC2=CC=CN12', 'C=1C=C2C=CC=C(O)N2C=1'),  # fused r5+r6 indolizinone
    ('O=C1NC=CN2C=CC=C12', 'C=1C=C2C(O)=NC=CN2C=1'),  # fused r5+r6 imidazo[1,2-a]pyridinone
    ('O=C1CN2C=CC=C2C=C1', 'C=1N2C(=CC=C2)C=CC=1O'),  # fused r5+r6 pyridinone bridge N
    ('O=C1CN2C=CC=C2C(=O)C1', 'C=1N2C=CC=C2C(O)=CC=1O'),  # fused r5+r6 diketone N between
    ('O=C1N2C=CC=C2CC(=O)C1', 'C=1(O)N2C=CC=C2C=C(O)C=1'),  # fused r5+r6 diketone N direct
    ('O=C1N2C=CC=C2C(=O)CC1', 'C=1C=C(O)N2C(=CC=C2)C=1O'),  # fused r5+r6 diketone both adj N
    ('O=C1CSC2=CC=CN12', 'C1=CC=C2SC=C(O)N12'),  # fused r5+r5 thiazolone
    ('O=C1COC2=CC=CN12', 'C1=CC=C2OC=C(O)N12'),  # fused r5+r5 oxazolone
    ('O=C1CNN2C=CC=C12', 'C1=CC=C2C(O)=CNN12'),  # fused r5+r5 N-C bridge
    ('O=C1CC2=CC=CN2N1', 'C1=CC=C2C=C(O)NN12'),  # fused r5+r5 C-N bridge
    ('CN1N=[N+](C)CC1=O', 'C1(O)=C[N+](C)=NN1C'),  # charged triazolium azolone
    ('C[N+]1=NOC(=O)C1', 'C[N+]1=NOC(O)=C1'),  # charged oxadiazolium
    ('CN1CC(=O)[N+](C)=N1', 'C[N+]1=NN(C=C1O)C'),  # charged azolone alt N+
    ('C[N+]1=CC=CC(=O)N1', '[N+]=1(C)N=C(C=CC=1)O'),  # charged pyridinone
]


@mark.parametrize('raw,result', tautomer_fix_data)
def test_tautomer_fix(raw, result):
    tmp = smiles(raw)
    tmp.standardize()
    assert tmp == smiles(result), f'{raw} > {tmp} != {result}'


def test_pyrylium_fix():
    # O in ring with C=[N+] -> [O+] in ring with C-N
    tmp = smiles('O1C(=[NH+]C)C=CC=C1')
    tmp.standardize()
    expected = smiles('C=1[O+]=C(NC)C=CC=1')
    assert tmp == expected
