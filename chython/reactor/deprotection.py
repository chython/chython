# -*- coding: utf-8 -*-
#
#  Copyright 2022-2025 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from .. import smarts, MoleculeContainer
from .transformer import Transformer

"""
Predefined transformers for most common protection groups cleavage.
"""

_hydroxyl_thiocarbamate = (  # NaIO4 or H2O2/NaOH
    ('[O;D2;x0:1]-;!@[C;x3;z2](=[S;D1])[N;D3;x0]([C;D1])[C;D1]', '[A:1]',  # rule
     'CC(C)OC(=S)N(C)C', 'CC(C)O'),  # test
)

_hydroxyl_fmoc = (  # Et3N pKa ~ 10
    ('[O;D2;x0:1]-;!@[C;z2;x3](=O)[O;D2;x0]-[C;D2;x1;z1][C;D3;z1;x0;r5]1[C;a;r6]:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D3]:2-[C;a;r6]:3:[C;D2]:[C;D2]:[C;D2]:[C;D2]:C1:3',
     '[A:1]',
     'CC(C)OC(=O)OCC1C2=CC=CC=C2C2=C1C=CC=C2', 'CC(C)O'),
)

_hydroxyl_troc = (  # [Zn]
    ('[O;D2;x0:1]-;!@[C;z2;x3](=O)[O;D2;x0]-[C;D2][C;D4;x3]([Cl;D1])([Cl;D1])[Cl;D1]', '[A:1]',
     'CC(C)OC(=O)OCC(Cl)(Cl)Cl', 'CC(C)O'),
)

_hydroxyl_teoc = (  # [F-]
    ('[O;D2;x0:1]-;!@[C;z2;x3](=O)[O;D2;x0]-[C;D2;z1;x1][C;D2;x1;z1][Si;D4;z1;x0]([C;D1])([C;D1])[C;D1]', '[A:1]',
     'CC(C)OC(=O)OCC[Si](C)(C)C', 'CC(C)O'),
)

_hydroxyl_alloc = (  # [Pd] + NuH
    ('[O;D2;x0:1]-;!@[C;z2;x3](=O)[O;D2;x0]-[C;D2;z1;x1][C;D2;x0;z2]=[C;D1]', '[A:1]',
     'CC(C)OC(=O)OCC=C', 'CC(C)O'),
)

_hydroxyl_allyl = (  # basic or Metal isomerization + hydrolysis
    ('[O;D2;x0:1]-;!@[C;D2;z1;x1][C;D2;x0;z2]=[C;D1]', '[A:1]',
     'CC(C)OCC=C', 'CC(C)O'),
)

_hydroxyl_tms = (  # [F-] ion substitution
    ('[O;D2;x1:1]-;!@[Si;D4;z1;x1]([C;D1])([C;D1])[C;D1]', '[A:1]',
     'CC(C)O[Si](C)(C)C', 'CC(C)O', 'CC(C)O[SiH](C)C', 'CC(C)O[Si](C)(C)OC', 'CC(C)O[Si](C)(C)CC'),
)

_hydroxyl_tes = (
    ('[O;D2;x1:1]-;!@[Si;D4;z1;x1]([C;D2;x1;z1][C;D1])([C;D2;x1;z1][C;D1])[C;D2;x1;z1][C;D1]', '[A:1]',
     'CC(C)O[Si](CC)(CC)CC', 'CC(C)O', 'CC(C)O[SiH](C)C', 'CC(C)O[Si](C)(C)OC'),
)

_hydroxyl_tbs = (  # TBS / TBDMS
    ('[O;D2;x1:1]-;!@[Si;D4;z1;x1]([C;D1])([C;D1])[C;D4;x1;z1]([C;D1])([C;D1])[C;D1]', '[A:1]',
     'CC(C)O[Si](C)(C)C(C)(C)C', 'CC(C)O', 'CC(C)O[SiH](C)C', 'CC(C)O[Si](C)(C)OC'),
)

_hydroxyl_tips = (
    ('[O;D2;x1:1]-;!@[Si;D4;z1;x1]([C;D3;z1;x1]([C;D1])[C;D1])([C;D3;z1;x1]([C;D1])[C;D1])[C;D3;z1;x1]([C;D1])[C;D1]', '[A:1]',
     'CC(C)O[Si](C(C)C)(C(C)C)C(C)C', 'CC(C)O', 'CC(C)O[SiH](C)C', 'CC(C)O[Si](C)(C)OC'),
)

_hydroxyl_tbdps = (
    ('[O;D2;x1:1]-;!@[Si;D4;z1;x1]([C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)([C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)[C;D4;x1;z1]([C;D1])([C;D1])[C;D1]', '[A:1]',
     'CC(C)O[Si](c1ccccc1)(c1ccccc1)C(C)(C)C', 'CC(C)O', 'CC(C)O[SiH](C)C', 'CC(C)O[Si](c1ccc(C)cc1)(c1ccccc1)C(C)(C)C'),
)

_hydroxyl_benzyl = (  # [H], ...
    ('[O;D2;x0:1]-;!@[C;D2;z1;x1]-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1', '[A:1]',
     'CC(C)OCc1ccccc1', 'CC(C)O', 'CC(C)OCc1cccc(C)c1', 'CC(C)OC(C)c1ccccc1'),  # test + decoys
)

_hydroxyl_o_nitrobenzyl = (  # UV-light
    ('[O;D2;x0:1]-;!@[C;D2;z1;x1]-[C;a;r6]:1:[C;D3;x1]([N+](=O)[O-]):[C;D2]:[C;D2]:[C;D2]:[C;D2]:1', '[A:1]',
     'CC(C)OCc1c(N(=O)=O)cccc1', 'CC(C)O', 'CC(C)OC(OC)c1c(N(=O)=O)cccc1'),
)

_hydroxyl_methoxy_benzyl = (  # PMB or MPM
    ('[O;D2;x0:1]-;!@[C;D2;z1;x1]-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D2]:1', '[A:1]',
     'CC(C)OCc1ccc(OC)cc1', 'CC(C)O', 'CC(C)OCc1ccc(OCC)cc1', 'CC(C)OCc1cc(OC)ccc1'),
)

_hydroxyl_dimethoxybenzyl = (
    ('[O;D2;x0:1]-;!@[C;D2;z1;x1]-[C;a;r6]:1:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D2]:1', '[A:1]',
     'CC(C)OCc1c(OC)cc(OC)cc1', 'CC(C)O'),
)

_hydroxyl_naphthyl = (  # Nap
    ('[O;D2;x0:1]-;!@[C;D2;z1;x1]-[C;a;r6]:1:[C;D2]:[C;D3]:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D3]:2:[C;D2]:[C;D2]:1', '[A:1]',
     'CC(C)OCC1=CC2=C(C=CC=C2)C=C1', 'CC(C)O', 'CC(C)OCc1ccccc1'),
)

_hydroxyl_bom = (  # like Bn
    ('[O;D2;x0:1]-;!@[C;D2;x2;z1][O;D2;x0][C;D2;z1;x1]-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1', '[A:1]',
     'CC(C)OCOCc1ccccc1', 'CC(C)O'),
)

_hydroxyl_piv = (
    ('[O;D2;x0:1]-;!@[C;z2;x2](=O)-[C;D4;x0;z1]([C;D1])([C;D1])[C;D1]', '[A:1]',
     'CC(C)OC(=O)C(C)(C)C', 'CC(C)O'),
)

_hydroxyl_methoxy_benzoate = (  # pMeO-Bz
    ('[O;D2;x0:1]-;!@[C;z2;x2](=O)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D2]:1', '[A:1]',
     'COC1=CC=C(C=C1)C(=O)OC(C)C', 'CC(C)O', 'C1=CC=C(C=C1)C(=O)OC(C)C'),
)

_hydroxyl_benzoate = (  # Bz
    ('[O;D2;x0:1]-;!@[C;z2;x2](=O)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1', '[A:1]',
     'C1=CC=C(C=C1)C(=O)OC(C)C', 'CC(C)O', 'COC1=CC=C(C=C1)C(=O)OC(C)C'),
)

_hydroxyl_acyl = (  # Ac
    ('[O;D2;x0:1]-;!@[C;z2;x2](=O)-[C;D1]', '[A:1]',
     'CC(C)OC(=O)C', 'CC(C)O', 'CC(C)OC(=O)CC'),
)

_hydroxyl_tfa = (
    ('[O;D2;x0:1]-;!@[C;z2;x2](=O)-[C;D4;z1;x3](F)(F)F', '[A:1]',
     'CC(C)OC(=O)C(F)(F)F', 'CC(C)O'),
)

_hydroxyl_mom = (
   ('[O;D2;x0:1]-;!@[C;D2;x2;z1][O;D2;x0][C;D1]', '[A:1]',
    'CC(C)OCOC', 'CC(C)O', 'CC(C)OC(C)OC'),
)

_hydroxyl_mem = (
    ('[O;D2;x0:1]-;!@[C;D2;x2;z1][O;D2;x0][C;D2;z1;x1][C;D2;z1;x1][O;D2;x0][C;D1]', '[A:1]',
     'COCCOCOC(C)C', 'CC(C)O'),
)

_hydroxyl_thp = (
    ('[O;D2;x0:1]-;!@[C;D3;x2;z1;r6]1[O;D2][C;D2][C;D2][C;D2][C;D2]1', '[A:1]',
     'CC(C)OC1CCCCO1', 'CC(C)O'),
)

_hydroxyl_ee = (
    ('[O;D2;x0:1]-;!@[C;D3;x2;z1]([O;D2;x0][C;D2;x1;z1][C;D1])[C;D1]', '[A:1]',
     'CC(C)OC(C)OCC', 'CC(C)O', 'CC(C)OC(CC)OCC'),
)

_hydroxyl_mop = (
    ('[O;D2;x0:1]-;!@[C;D4;x2;z1]([O;D2;x0][C;D1])([C;D1])[C;D1]', '[A:1]',
     'CC(C)OC(C)(C)OC', 'CC(C)O'),
)

_hydroxyl_sem = (
    ('[O;D2;x0:1]-;!@[C;D2;x2;z1][O;D2;x0][C;D2;z1;x1][C;D2;z1;x1][Si;D4;z1;x0]([C;D1])([C;D1])[C;D1]', '[A:1]',
     'CC(C)OCOCC[Si](C)(C)C', 'CC(C)O'),
)

_hydroxyl_tritil = (
    ('[O;D2;x0:1]-;!@[C;D4;z1;x1](-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)(-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1', '[A:1]',
     'CC(C)OC(c1ccccc1)(c1ccccc1)c1ccccc1', 'CC(C)O', 'COc1ccc(cc1)C(OC(C)C)(c1ccccc1)c1ccc(OC)cc1'),
)

_hydroxyl_dimetoxy_tritil = (
    ('[O;D2;x0:1]-;!@[C;D4;z1;x1](-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D2]:1)(-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D2]:1)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1', '[A:1]',
     'COc1ccc(cc1)C(OC(C)C)(c1ccccc1)c1ccc(OC)cc1', 'CC(C)O', 'CC(C)OC(c1ccccc1)(c1ccccc1)c1ccccc1'),
)

_hydroxyl_chloro_tritil = (
    ('[O;D2;x0:1]-;!@[C;D4;z1;x1](-[C;a;r6]:1:[C;D3;x1]([Cl;D1]):[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)(-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1', '[A:1]',
     'CC(C)OC(c1c(Cl)cccc1)(c1ccccc1)c1ccccc1', 'CC(C)O'),
)

_hydroxyl_mmt = (
    ('[O;D2;x0:1]-;!@[C;D4;z1;x1](-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D2]:1)(-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1', '[A:1]',
     'c1ccc(cc1)C(OC(C)C)(c1ccccc1)c1ccc(OC)cc1', 'CC(C)O', 'CC(C)OC(c1ccccc1)(c1ccccc1)c1ccccc1'),
)

_hydroxyl_tbu = (
    ('[O;D2;x0:1]-;!@[C;D4;x1;z1]([C;D1])([C;D1])[C;D1]', '[A:1]',
     'CC(C)OC(C)(C)C', 'CC(C)O'),
)

_hydroxyl_methyl = (
    ('[O;D2;x0:1]-;!@[C;D1]', '[A:1]',
     'CC(C)OC', 'CC(C)O', 'CC(C)OCC'),
)

_hydroxyl_mpe = (
    ('[O;D2;x0:1]-;!@[C;D4;x1;z1]([C;D1])([C;D2;x0;z1][C;D1])[C;D2;x0;z1][C;D1]', '[A:1]',
     'CC(C)OC(CC)(CC)C', 'CC(C)O'),
)

_hydroxyl_trifluoroethyl = (
    ('[O;D2;x0:1]-;!@[C;D2;x1;z1][C;D4;x3;z1](F)(F)F', '[A:1]',
     'CC(C)OCC(F)(F)F', 'CC(C)O'),
)


_hydroxyl_dmab = (
    ('[O;D2;x0:1]-;!@[C;D2;x1;z1]-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D3;x1](:[C;D2]:[C;D2]:1)-[N;D2;x0;z1]-[C;z2;x1;D3]([C;D2;x0;z1][C;D3;x0;z1]([C;D1])[C;D1])=[C;D3;r6;x0;z2]1[C;x1;z2;D3](=O)[C;D2][C;D4;x0;z1]([C;D1])([C;D1])[C;D2][C;D3;x1;z2]1=O', '[A:1]',
     'CC(C)CC(NC1=CC=C(COC(C)=O)C=C1)=C1C(=O)CC(C)(C)CC1=O', 'CC(O)=O'),
    ('[O;D2;x0:1]-;!@[C;D2;x1;z1]-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D3;x1](:[C;D2]:[C;D2]:1)-[N;D2;x0;z2]=[C;x1;D3]([C;D2;x0;z1][C;D3;x0;z1]([C;D1])[C;D1])-[C;D3;r6;x0;z1]1[C;x1;z2;D3](=O)[C;D2][C;D4;x0;z1]([C;D1])([C;D1])[C;D2][C;D3;x1;z2]1=O', '[A:1]',
     'CC(C)CC(=NC1=CC=C(COC(C)=O)C=C1)C1C(=O)CC(C)(C)CC1=O', 'CC(O)=O'),
)

_diol_acetone = (
    # 1,2
    ('[O;D2;x0:1]1-;@[C;D4;x2;z1]([C;D1])([C;D1])-[O;D2;x0:2][C;M]!#[C;M]1', '[A:1].[A:2]',
     'CC1COC(C)(C)O1', 'CC(O)CO'),
    # 1,3
    ('[O;D2;x0:1]1-;@[C;D4;x2;z1]([C;D1])([C;D1])-[O;D2;x0:2][C;M][C;M]!#[C;M]1', '[A:1].[A:2]',
     'CC1CCOC(C)(C)O1', 'CC(O)CCO'),
)

_hydroxyl_amine_acetone = (
    ('[O;D2;x0:1]1-;@[C;D4;x2;z1]([C;D1])([C;D1])-[N;z1:2][C;M]!#[C;M]1', '[A:1].[A:2]',
     'CC1CN(C(C)=O)C(C)(C)O1', 'CC(O)CNC(C)=O'),
)

_diol_formalin = (
    # 1,2
    ('[O;D2;x0:1]1-;@[C;D2;x2;z1]-[O;D2;x0:2][C;M]!#[C;M]1', '[A:1].[A:2]',
     'CC1COCO1', 'CC(O)CO'),
    # 1,3
    ('[O;D2;x0:1]1-;@[C;D2;x2;z1]-[O;D2;x0:2][C;M][C;M]!#[C;M]1', '[A:1].[A:2]',
     'CC1CCOCO1', 'CC(O)CCO'),
)

_diol_cyclopentanone = (
    # 1,2
    ('[O;D2;x0:1]1-;@[C;D4;x2;z1]2([C;D2;x0;z1][C;D2;x0;z1][C;D2;x0;z1][C;D2;x0;z1]2)-[O;D2;x0:2][C;M]!#[C;M]1', '[A:1].[A:2]',
     'CC1COC2(CCCC2)O1', 'CC(O)CO'),
    # 1,3
    ('[O;D2;x0:1]1-;@[C;D4;x2;z1]2([C;D2;x0;z1][C;D2;x0;z1][C;D2;x0;z1][C;D2;x0;z1]2)-[O;D2;x0:2][C;M][C;M]!#[C;M]1', '[A:1].[A:2]',
     'CC1CCOC2(CCCC2)O1', 'CC(O)CCO'),
)

_diol_cyclohexanone = (
    # 1,2
    ('[O;D2;x0:1]1-;@[C;D4;x2;z1]2([C;D2;x0;z1][C;D2;x0;z1][C;D2;x0;z1][C;D2;x0;z1][C;D2;x0;z1]2)-[O;D2;x0:2][C;M]!#[C;M]1', '[A:1].[A:2]',
     'CC1COC2(CCCCC2)O1', 'CC(O)CO'),
    # 1,3
    ('[O;D2;x0:1]1-;@[C;D4;x2;z1]2([C;D2;x0;z1][C;D2;x0;z1][C;D2;x0;z1][C;D2;x0;z1][C;D2;x0;z1]2)-[O;D2;x0:2][C;M][C;M]!#[C;M]1', '[A:1].[A:2]',
     'CC1CCOC2(CCCCC2)O1', 'CC(O)CCO'),
)

_diol_diacetal = (
    # 1,2
    ('[O;D2;x0:1]1-;@[C;D4;x2;z1]([O;D2;x0][C;D1])([C;D1])[C;D4;x2;z1]([O;D2;x0][C;D1])([C;D1])-[O;D2;x0:2][C;M]!#[C;M]1', '[A:1].[A:2]',
     'COC1(C)OCC(C)OC1(C)OC', 'CC(O)CO'),
    # 1,3
    ('[O;D2;x0:1]1-;@[C;D4;x2;z1]([O;D2;x0][C;D1])([C;D1])[C;D4;x2;z1]([O;D2;x0][C;D1])([C;D1])-[O;D2;x0:2][C;M][C;M]!#[C;M]1', '[A:1].[A:2]',
     'COC1(C)OCCC(C)OC1(C)OC', 'CC(O)CCO'),
)

_diol_benzylidene = (
    # 1,2
    ('[O;D2;x0:1]1-;@[C;D3;x2;z1]([C;a;r6]:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:2)-[O;D2;x0:2][C;M]!#[C;M]1', '[A:1].[A:2]',
     'CC1COC(O1)c1ccccc1', 'CC(O)CO'),
    # 1,3
    ('[O;D2;x0:1]1-;@[C;D3;x2;z1]([C;a;r6]:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:2)-[O;D2;x0:2][C;M][C;M]!#[C;M]1', '[A:1].[A:2]',
     'CC1CCOC(O1)c1ccccc1', 'CC(O)CCO'),
)

_carbonyl_dithiolane = (  # MeI - S methylation + hydrolysis
    ('[C;D3,D4;z1;x2:1]1[S;D2:3][C;D2][C;D2][S;D2]1', '[A:1]=O'),
)

_carbonyl_dithiane = (  # MeI - S methylation + hydrolysis
    ('[C;D3,D4;z1;x2:1]1[S;D2:3][C;D2][C;D2][C;D2][S;D2]1', '[A:1]=O'),
)

_carbonyl_dimethylsulfide = (  # MeI - S methylation + hydrolysis
    ('[C;D3,D4;z1;x2:1]([S;D2:3][C;D1])[S;D2][C;D1]', '[A:1]=O'),
)

_carbonyl_dioxolane = (
    ('[C;D3,D4;z1;x2:1]1[O;D2:3][C;D2][C;D2][O;D2]1', '[A:1]=O'),
)

_carbonyl_dioxane = (
    ('[C;D3,D4;z1;x2:1]1[O;D2:3][C;D2][C;D2][C;D2][O;D2]1', '[A:1]=O'),
)

_carbonyl_dimethoxy = (
    ('[C;D3,D4;z1;x2:1]([O;D2:3][C;D1])[O;D2][C;D1]', '[A:1]=O'),
)

_carboxyl_trioxabicyclooctane = (  # [H+]. Note! a second step of basic hydrolysis is required.
    ('[C;D4;x3:1]12[O;D2:4][C;D2;x1;z1][C;D4;x0;z1]([C;D1])([C;D2;x1;z1][O;D2]1)[C;D2;x1;z1][O;D2]2', '[A:1](=O)O',
     'CC(C)C12OCC(C)(CO1)CO2', 'CC(C)C(O)=O', 'CC(C)C12OCC(CC)(CO1)CO2', 'CC(C)C12OC(C)C(C)(CO1)CO2'),
)

_amine_methylcarbamate = (  # PrSLi or [OH-]
    # Ar-NH2
    ('[C;a:1][N;D2:2]-C(=O)O[C;D1]', '[A:1][A:2]',
     'c1ccccc1NC(=O)OC', 'c1ccccc1N', 'c1ccccc1NC(=O)OCC', 'c1cccn1NC(=O)OC'),
    # Alk-NH2
    ('[C;D2,D3,D4;z1;x1:1][N;D2:2]-C(=O)O[C;D1]', '[A:1][A:2]',
     'CC(C)NC(=O)OC', 'CC(C)N', 'CC(C)NC(=O)OCC', 'CONC(=O)OC', 'C=NC(=O)OC'),
    # Alk-NH-Ar
    ('[C;a:1][N:2]([C;z1;x1:3])-C(=O)O[C;D1]', '[A:1][A:2][A:3]',
     'c1ccccc1N(C(C)C)C(=O)OC', 'c1ccccc1NC(C)C'),
    # Alk2NH
    ('[C;D2,D3,D4;z1;x1:1][N:2]([C;z1;x1:3])-C(=O)O[C;D1]', '[A:1][A:2][A:3]',
     'CC(C)N(C(C)C)C(=O)OC', 'CC(C)NC(C)C'),
)

_amine_teoc = (  # [F-]
    # Ar-NH2
    ('[C;a:1][N;D2:2]-C(=O)O[C;D2][C;D2][Si]([C;D1])([C;D1])[C;D1]', '[A:1][A:2]',
     'c1ccccc1NC(=O)OCC[Si](C)(C)C', 'c1ccccc1N'),
    # Alk-NH2
    ('[C;D2,D3,D4;z1;x1:1][N;D2:2]-C(=O)O[C;D2][C;D2][Si]([C;D1])([C;D1])[C;D1]', '[A:1][A:2]'),
    ('[C;a:1][N:2]([C;z1;x1:3])-C(=O)O[C;D2][C;D2][Si]([C;D1])([C;D1])[C;D1]', '[A:1][A:2][A:3]'),  # Alk-NH-Ar
    ('[C;D2,D3,D4;z1;x1:1][N:2]([C;z1;x1:3])-C(=O)O[C;D2][C;D2][Si]([C;D1])([C;D1])[C;D1]', '[A:1][A:2][A:3]'),  # Alk2NH
)

_amine_troc = (  # [Zn]
    # Ar-NH2
    ('[C;a:1][N;D2:2]-C(=O)O[C;D2]C(Cl)(Cl)Cl', '[A:1][A:2]',
     'c1ccccc1NC(=O)OCC(Cl)(Cl)Cl', 'c1ccccc1N', 'c1ccccc1NC(=O)OC(C)C(Cl)(Cl)Cl'),
    # Alk-NH2
    ('[C;D2,D3,D4;z1;x1:1][N;D2:2]-C(=O)O[C;D2]C(Cl)(Cl)Cl', '[A:1][A:2]'),
    # Alk-NH-Ar
    ('[C;a:1][N:2]([C;z1;x1:3])-C(=O)O[C;D2]C(Cl)(Cl)Cl', '[A:1][A:2][A:3]'),
    # Alk2NH
    ('[C;D2,D3,D4;z1;x1:1][N:2]([C;z1;x1:3])-C(=O)O[C;D2]C(Cl)(Cl)Cl', '[A:1][A:2][A:3]'),
)

_amine_alloc = (  # [Pd]
    # Ar-NH2
    ('[C;a:1][N;D2:2]-C(=O)O[C;D2][C;D2]=[C;D1]', '[A:1][A:2]',
     'c1ccccc1NC(=O)OCC=C', 'c1ccccc1N', 'c1ccccc1NC(=O)OCC=CC'),
    # Alk-NH2
    ('[C;D2,D3,D4;z1;x1:1][N;D2:2]-C(=O)O[C;D2][C;D2]=[C;D1]', '[A:1][A:2]'),
    # Alk-NH-Ar
    ('[C;a:1][N:2]([C;z1;x1:3])-C(=O)O[C;D2][C;D2]=[C;D1]', '[A:1][A:2][A:3]'),
    # Alk2NH
    ('[C;D2,D3,D4;z1;x1:1][N:2]([C;z1;x1:3])-C(=O)O[C;D2][C;D2]=[C;D1]', '[A:1][A:2][A:3]'),
)

_amine_cbz = (  # [Pd] or Na/NH3
    # Ar-NH2
    ('[C;a:1][N;D2:2]-C(=O)O[C;D2]C:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1', '[A:1][A:2]',
     'c1ccccc1NC(=O)OCc2ccccc2', 'c1ccccc1N', 'c1ccccc1NC(=O)OC(C)c2ccccc2'),
    # Alk-NH2
    ('[C;D2,D3,D4;z1;x1:1][N;D2:2]-C(=O)O[C;D2]C:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1', '[A:1][A:2]'),
    # Alk-NH-Ar
    ('[C;a:1][N:2]([C;z1;x1:3])-C(=O)O[C;D2]C:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1', '[A:1][A:2][A:3]'),
    # Alk2NH
    ('[C;D2,D3,D4;z1;x1:1][N:2]([C;z1;x1:3])-C(=O)O[C;D2]C:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1', '[A:1][A:2][A:3]'),
)

_amine_chloro_cbz = (
    # Ar-NH2
    ('[C;a;M][N;D2:1]-C(=O)O[C;D2]C:1:C([Cl;D1]):[C;D2]:[C;D2]:[C;D2]:[C;D2]:1', '[A:1]',
     'Clc1ccccc1COC(=O)Nc1ccccc1', 'c1ccccc1N', 'c1ccccc1NC(=O)OC(C)c2ccccc2'),
    # Alk-NH2
    ('[C;D2,D3,D4;z1;x1;M][N;D2:1]-C(=O)O[C;D2]C:1:C([Cl;D1]):[C;D2]:[C;D2]:[C;D2]:[C;D2]:1', '[A:1]'),
    # Alk-NH-Ar
    ('[C;a;M][N:1]([C;z1;x1;M])-C(=O)O[C;D2]C:1:C([Cl;D1]):[C;D2]:[C;D2]:[C;D2]:[C;D2]:1', '[A:1]'),
    # Alk2NH
    ('[C;D2,D3,D4;z1;x1;M][N:1]([C;z1;x1;M])-C(=O)O[C;D2]C:1:C([Cl;D1]):[C;D2]:[C;D2]:[C;D2]:[C;D2]:1', '[A:1]'),
)

_amine_nosyl = (  # NS. With SH-CH2-CH2-OH
    # Ar-NH2
    ('[C;a:1][N;D2:2]-S(=O)(=O)C:1:C([N+](=O)[O-]):[C;D2]:[C;D2]:[C;D2]:[C;D2]:1', '[A:1][A:2]',
     '[O-][N+](=O)c1ccccc1S(=O)(=O)Nc1ccccc1', 'c1ccccc1N'),
    # Alk-NH2
    ('[C;D2,D3,D4;z1;x1:1][N;D2:2]-S(=O)(=O)C:1:C([N+](=O)[O-]):[C;D2]:[C;D2]:[C;D2]:[C;D2]:1', '[A:1][A:2]'),
    # Alk-NH-Ar
    ('[C;a:1][N:2]([C;z1;x1:3])-S(=O)(=O)C:1:C([N+](=O)[O-]):[C;D2]:[C;D2]:[C;D2]:[C;D2]:1', '[A:1][A:2][A:3]'),
    # Alk2NH
    ('[C;D2,D3,D4;z1;x1:1][N:2]([C;z1;x1:3])-S(=O)(=O)C:1:C([N+](=O)[O-]):[C;D2]:[C;D2]:[C;D2]:[C;D2]:1',
     '[A:1][A:2][A:3]'),
)

_amine_boc = (
    # Ar-NH2
    ('[C;a;M][N;D2;x0;z1:1]-[C;x3;z2](=O)[O;x0;z1][C;D4;x1]([C;D1])([C;D1])[C;D1]', '[A:1]',
     'c1ccccc1NC(=O)OC(C)(C)C', 'c1ccccc1N'),
    # Alk-NH2
    ('[C;D2,D3,D4;z1;x1;M][N;D2;x0;z1:1]-[C;x3;z2](=O)[O;x0;z1][C;D4;x1]([C;D1])([C;D1])[C;D1]', '[A:1]'),
    # Alk-NH-Ar
    ('[C;a;M][N;D3;x0;z1:1]([C;z1;x1;M])-[C;x3;z2](=O)[O;x0;z1][C;D4;x1]([C;D1])([C;D1])[C;D1]', '[A:1]'),
    # Alk2NH
    ('[C;D2,D3,D4;z1;x1;M][N;D3;x0;z1:1]([C;z1;x1;M])-[C;x3;z2](=O)[O;x0;z1][C;D4;x1]([C;D1])([C;D1])[C;D1]', '[A:1]'),
    # Alk[NH]-O-Alk
    ('[C;D2,D3,D4;z1;x1;M][N;D3;x1;z1:1]([O;D2;x1;z1;M][C;z1;x1;M])-[C;x3;z2](=O)[O;x0;z1][C;D4;x1]([C;D1])([C;D1])[C;D1]', '[A:1]'),
    # Alk-[NH]-COC
    ('[C;D2,D3,D4;z1;x1;r5;M]1-;@[N;D3;x0;z1:1](-;@[C;z1;x2;M]-;@[O,S;D2;x0;z1;M][C;z1;x1;M]1)-[C;x3;z2](=O)[O;x0;z1][C;D4;x1]([C;D1])([C;D1])[C;D1]', '[A:1]'),
    ('[C;D2,D3,D4;z1;x1;r6,r7,r8,r9;M]-;@[N;D3;x0;z1:1](-;@[C;z1;x2;M]-;@[O,S;D2;x0;z1;M][C;z1;x1;M])-[C;x3;z2](=O)[O;x0;z1][C;D4;x1]([C;D1])([C;D1])[C;D1]', '[A:1]'),
    # amino-pyrrolidine
    ('[C;D3;z1;x2;r5;M](-;@[N;M])[N;D2;x0;z1:1]-[C;x3;z2](=O)[O;x0;z1][C;D4;x1]([C;D1])([C;D1])[C;D1]', '[A:1]',
     'CC(C)(C)OC(=O)NC1CCCN1', 'NC1CCCN1'),
)

_amine_tfa = (
    # Ar-NH2
    ('[C;a:1][N;D2:2]-C(=O)C(F)(F)F', '[A:1][A:2]'),
    # Alk-NH2
    ('[C;D2,D3,D4;z1;x1:1][N;D2:2]-C(=O)C(F)(F)F', '[A:1][A:2]'),
    # Alk-NH-Ar
    ('[C;a:1][N:2]([C;z1;x1:3])-C(=O)C(F)(F)F', '[A:1][A:2][A:3]'),
    # Alk2NH
    ('[C;D2,D3,D4;z1;x1:1][N:2]([C;z1;x1:3])-C(=O)C(F)(F)F', '[A:1][A:2][A:3]'),
)

_amine_fmoc = (
    # Ar-NH2
    ('[C;a:1][N;D2:2]-C(=O)O[C;D2][C;D3]1C:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:C:2-C:3:[C;D2]:[C;D2]:[C;D2]:[C;D2]:C1:3',
     '[A:1][A:2]',
     'O=C(Nc1ccccc1)OCC1c2ccccc2-c2ccccc12', 'c1ccccc1N'),
    # Alk-NH2
    ('[C;D2,D3,D4;z1;x1:1][N;D2:2]-C(=O)O[C;D2][C;D3]1C:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:C:2-C:3:[C;D2]:[C;D2]:[C;D2]:[C;D2]:C1:3',  # noqa
     '[A:1][A:2]'),
    # Alk-NH-Ar
    ('[C;a:1][N:2]([C;z1;x1:3])-C(=O)O[C;D2][C;D3]1C:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:C:2-C:3:[C;D2]:[C;D2]:[C;D2]:[C;D2]:C1:3',  # noqa
     '[A:1][A:2][A:3]'),
    # Alk2NH
    ('[C;D2,D3,D4;z1;x1:1][N:2]([C;z1;x1:3])-C(=O)O[C;D2][C;D3]1C:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:C:2-C:3:[C;D2]:[C;D2]:[C;D2]:[C;D2]:C1:3',  # noqa
     '[A:1][A:2][A:3]'),
)

_amine_pbf_pmc_mtr = (
    # Ar-NH2
    ('[C;a;M][N;D2:1]-S(=O)(=O)-C:1:C([C;D1]):C([C;D1]):C(O[C;x1;z1]):C:C([C;D1]):1',
     '[A:1]'),
    # Alk-NH2
    ('[C;D2,D3,D4;z1;x1;M][N;D2:1]-S(=O)(=O)-C:1:C([C;D1]):C([C;D1]):C(O[C;x1;z1]):C:C([C;D1]):1',
    '[A:1]'),
    # Alk-NH-Ar
    ('[C;a;M][N:1]([C;z1;x1;M])-S(=O)(=O)-C:1:C([C;D1]):C([C;D1]):C(O[C;x1;z1]):C:C([C;D1]):1',
    '[A:1]'),
    # Alk2NH
    ('[C;D2,D3,D4;z1;x1;M][N:1]([C;z1;x1;M])-S(=O)(=O)-C:1:C([C;D1]):C([C;D1]):C(O[C;x1;z1]):C:C([C;D1]):1',
    '[A:1]'),
    # Guanidine
    ('[N;M]=[C;M]([N;M])[N;D2:1]-S(=O)(=O)-C:1:C([C;D1]):C([C;D1]):C(O[C;x1;z1]):C:C([C;D1]):1',
     '[A:1]'),
    ('[N;M][C;M]([N;M])=[N:1]-S(=O)(=O)-C:1:C([C;D1]):C([C;D1]):C(O[C;x1;z1]):C:C([C;D1]):1',
     '[A:1]'),
)

_amine_dde = (
    # Ar-NH2
    ('[C;a;M][N;D2:1]-C([C;D1])=C1C(=O)[C;D2]C([C;D1])([C;D1])[C;D2]C1=O', '[A:1]'),
    # Alk-NH2
    ('[C;D2,D3,D4;z1;x1;M][N;D2:1]-C([C;D1])=C1C(=O)[C;D2]C([C;D1])([C;D1])[C;D2]C1=O', '[A:1]'),
    # Alk-NH-Ar
    ('[C;a;M][N:1]([C;z1;x1;M])-C([C;D1])=C1C(=O)[C;D2]C([C;D1])([C;D1])[C;D2]C1=O', '[A:1]'),
    # Alk2NH
    ('[C;D2,D3,D4;z1;x1;M][N:1]([C;z1;x1;M])-C([C;D1])=C1C(=O)[C;D2]C([C;D1])([C;D1])[C;D2]C1=O', '[A:1]'),
)

_amine_ivdde = (
    # Ar-NH2
    ('[C;a;M][N;D2:1]-C([C;D2][C;D3]([C;D1])[C;D1])=C1C(=O)[C;D2]C([C;D1])([C;D1])[C;D2]C1=O', '[A:1]'),
    # Alk-NH2
    ('[C;D2,D3,D4;z1;x1;M][N;D2:1]-C([C;D2][C;D3]([C;D1])[C;D1])=C1C(=O)[C;D2]C([C;D1])([C;D1])[C;D2]C1=O', '[A:1]'),
    # Alk-NH-Ar
    ('[C;a;M][N:1]([C;z1;x1;M])-C([C;D2][C;D3]([C;D1])[C;D1])=C1C(=O)[C;D2]C([C;D1])([C;D1])[C;D2]C1=O', '[A:1]'),
    # Alk2NH
    ('[C;D2,D3,D4;z1;x1;M][N:1]([C;z1;x1;M])-C([C;D2][C;D3]([C;D1])[C;D1])=C1C(=O)[C;D2]C([C;D1])([C;D1])[C;D2]C1=O',
     '[A:1]'),
)

_amine_mtt = (
    # Ar-NH2
    ('[C;a;M][N;D2:1]-C(C:1:[C;D2]:[C;D2]:C([C;D1]):[C;D2]:[C;D2]:1)(C:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:2)C:3:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:3',
     '[A:1]'),
    # Alk-NH2
    ('[C;D2,D3,D4;z1;x1;M][N;D2:1]-C(C:1:[C;D2]:[C;D2]:C([C;D1]):[C;D2]:[C;D2]:1)(C:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:2)C:3:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:3',
     '[A:1]'),
    # Alk-NH-Ar
    ('[C;a;M][N:1]([C;z1;x1;M])-C(C:1:[C;D2]:[C;D2]:C([C;D1]):[C;D2]:[C;D2]:1)(C:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:2)C:3:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:3',
     '[A:1]'),
    # Alk2NH
    ('[C;D2,D3,D4;z1;x1;M][N:1]([C;z1;x1;M])-C(C:1:[C;D2]:[C;D2]:C([C;D1]):[C;D2]:[C;D2]:1)(C:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:2)C:3:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:3',
     '[A:1]'),
)

_amine_bhoc = (
    # Ar-NH2
    ('[C;a;M][N;D2:1]-C(=O)O[C;D3](C:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)C:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:2',
     '[A:1]'),
    # Alk-NH2
    ('[C;D2,D3,D4;z1;x1;M][N;D2:1]-C(=O)O[C;D3](C:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)C:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:2',
     '[A:1]'),
    # Alk-NH-Ar
    ('[C;a;M][N:1]([C;z1;x1;M])-C(=O)O[C;D3](C:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)C:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:2',
     '[A:1]'),
    # Alk2NH
    ('[C;D2,D3,D4;z1;x1;M][N:1]([C;z1;x1;M])-C(=O)O[C;D3](C:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)C:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:2',
     '[A:1]'),
    # Guanidine
    ('[N;M]=[C;M]([N;M])[N;D2:1]-C(=O)O[C;D3](C:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)C:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:2',
     '[A:1]'),
    ('[N;M][C;M]([N;M])=[N:1]-C(=O)O[C;D3](C:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)C:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:2',
     '[A:1]'),
)

_amide_tritil = (
    ('[C;D3;x2;M](=[O;M])[N;D2:1]-C(C:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)(C:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:2)C:3:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:3',  # noqa
     '[A:1]'),
)

_amide_boc = (
    ('[N;D2,D3;x0;z1:1]([C;D3;x2;z2;M]=[O;M])-[C;x3;z2](=O)[O;x0;z1][C;D4;x1]([C;D1])([C;D1])[C;D1]', '[A:1]',
     'CC(=O)NC(=O)OC(C)(C)C', 'CC(N)=O'),
)

_thiol_tritil = (
    ('[S;D2;x0;z1:1]-;!@[C;D4;z1;x1](-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)(-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1', '[A:1]',
     'CC(C)SC(c1ccccc1)(c1ccccc1)c1ccccc1', 'CC(C)S'),
)

_thiol_mmt = (
    ('[S;D2;x0;z1:1]-;!@[C;D4;z1;x1](-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D2]:1)(-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1', '[A:1]',
     'c1ccc(cc1)C(SC(C)C)(c1ccccc1)c1ccc(OC)cc1', 'CC(C)S'),
)

_thiol_dimetoxy_tritil = (
    ('[S;D2;x0;z1:1]-;!@[C;D4;z1;x1](-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D2]:1)(-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D2]:1)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1', '[A:1]',
     'COc1ccc(cc1)C(SC(C)C)(c1ccccc1)c1ccc(OC)cc1', 'CC(C)S'),
)

_thiol_chloro_tritil = (
    ('[S;D2;x0;z1:1]-;!@[C;D4;z1;x1](-[C;a;r6]:1:[C;D3;x1]([Cl;D1]):[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)(-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1', '[A:1]',
     'CC(C)SC(c1c(Cl)cccc1)(c1ccccc1)c1ccccc1', 'CC(C)S'),
)

_thiol_benzyl = (
    ('[S;D2;x0;z1:1]-;!@[C;D2;z1;x1]-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1', '[A:1]',
     'CC(C)SCc1ccccc1', 'CC(C)S'),
)

_thiol_tbu = (
    ('[S;D2;x0;z1:1]-;!@[C;D4;x1;z1]([C;D1])([C;D1])[C;D1]', '[A:1]',
     'CC(C)SC(C)(C)C', 'CC(C)S'),
)

_thiol_stbu = (
    ('[S;D2;x1;z1:1]-;!@[S;D2;z1;x1]-[C;D4;x1;z1]([C;D1])([C;D1])[C;D1]', '[A:1]',
     'CC(C)SSC(C)(C)C', 'CC(C)S'),
)

_thiol_strimethoxyphenyl = (
    ('[S;D2;x1;z1:1]-;!@[S;D2;z1;x1]-[C;a;r6]:1:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D3;x1]:1[O;D2;x0][C;D1]', '[A:1]',
     'COc1cc(OC)c(SSC(C)C)c(OC)c1', 'CC(C)S'),
)

_thiol_amine_dimethoxybenzyl = (
    ('[S;D2;r5;x0;z1:1]1[C;M][C;M][N;z1:2]-[C;D3;x2;z1]1-[C;a;r6]:1:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D2]:1', '[A:1].[A:2]',
     'COC1=CC=C(C2NC(C)CS2)C(OC)=C1', 'NC(C)CS'),
)

_pyrrole_boc = (
    ('[N;a;r5:1]-C(=O)OC([C;D1])([C;D1])[C;D1]', '[A:1]'),
)

_pyrrole_chloro_tritil = (
    ('[N;a;r5:1]-C(C:1:C([Cl;D1]):[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)(C:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:2)C:3:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:3',
     '[A:1]'),
)


#################
# Magic Factory #
#################

_groups = [k[1:] for k, v in globals().items() if k.startswith('_') and isinstance(v, tuple) and v]
__all__ = ['apply_all'] + _groups
_cache = {}


def _prepare_reactor(rules, name):
    rxn = [Transformer(smarts(r), smarts(p)) for r, p, *_ in rules]

    def w(molecule: MoleculeContainer, /) -> MoleculeContainer:
        """
        Remove protective groups from the given molecule if applicable.
        """
        for r in rxn:
            while True:
                try:
                    molecule = next(r(molecule))
                except StopIteration:
                    break
        return molecule

    w.__module__ = __name__
    w.__qualname__ = w.__name__ = name
    return w


def apply_all(molecule: MoleculeContainer, /) -> MoleculeContainer:
    """
    Remove all found protective groups from the given molecule.
    """
    for name in _groups:
        molecule = __getattr__(name)(molecule)
    return molecule


def __getattr__(name):
    try:
        return _cache[name]
    except KeyError:
        if name in _groups:
            _cache[name] = t = _prepare_reactor(globals()[f'_{name}'], name)
            return t
        raise AttributeError


def __dir__():
    return __all__
