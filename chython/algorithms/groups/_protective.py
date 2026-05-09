# -*- coding: utf-8 -*-
#
#  Copyright 2022-2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from lazy_object_proxy import Proxy


def _rules():
    from ... import smarts

    rules = {}

    rules['hydroxyl_thiocarbamate'] = (
        smarts('[O;D2:1]-;!@[C;x3;z2](=[S;D1])[N;D3;x0]([C;D1])[C;D1]'),
        [1],  # atoms to keep
        [],  # atoms to add (atom_number, atom_type, bond_type)
        'CC(C)OC(=S)N(C)C',  # protected
        'CC(C)O',  # cleaved
        []  # optional decoys
    )

    rules['hydroxyl_fmoc'] = (
        smarts('[O;D2:1]-;!@[C;z2;x3](=O)[O;D2;x0]-[C;D2;x1;z1][C;D3;z1;x0;r5]1[C;a;r6]:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D3]:2-[C;a;r6]:3:[C;D2]:[C;D2]:[C;D2]:[C;D2]:C1:3'),
        [1],
        [],
        'CC(C)OC(=O)OCC1C2=CC=CC=C2C2=C1C=CC=C2',
        'CC(C)O',
        []
    )

    rules['hydroxyl_troc'] = (
        smarts('[O;D2:1]-;!@[C;z2;x3](=O)[O;D2;x0]-[C;D2][C;D4;x3]([Cl;D1])([Cl;D1])[Cl;D1]'),
        [1],
        [],
        'CC(C)OC(=O)OCC(Cl)(Cl)Cl',
        'CC(C)O',
        []
    )

    rules['hydroxyl_teoc'] = (
        smarts('[O;D2:1]-;!@[C;z2;x3](=O)[O;D2;x0]-[C;D2;z1;x1][C;D2;x1;z1][Si;D4;z1;x0]([C;D1])([C;D1])[C;D1]'),
        [1],
        [],
        'CC(C)OC(=O)OCC[Si](C)(C)C',
        'CC(C)O',
        []
    )

    rules['hydroxyl_alloc'] = (
        smarts('[O;D2:1]-;!@[C;z2;x3](=O)[O;D2;x0]-[C;D2;z1;x1][C;D2;x0;z2]=[C;D1]'),
        [1],
        [],
        'CC(C)OC(=O)OCC=C',
        'CC(C)O',
        []
    )

    rules['hydroxyl_tms'] = (
        smarts('[O;D2:1]-;!@[Si;D4;z1;x1]([C;D1])([C;D1])[C;D1]'),
        [1],
        [],
        'CC(C)O[Si](C)(C)C',
        'CC(C)O',
        ['CC(C)O[SiH](C)C', 'CC(C)O[Si](C)(C)OC', 'CC(C)O[Si](C)(C)CC']
    )

    rules['hydroxyl_tes'] = (
        smarts('[O;D2:1]-;!@[Si;D4;z1;x1]([C;D2;x1;z1][C;D1])([C;D2;x1;z1][C;D1])[C;D2;x1;z1][C;D1]'),
        [1],
        [],
        'CC(C)O[Si](CC)(CC)CC',
        'CC(C)O',
        ['CC(C)O[SiH](C)C', 'CC(C)O[Si](C)(C)OC']
    )

    rules['hydroxyl_tbs'] = (
        smarts('[O;D2:1]-;!@[Si;D4;z1;x1]([C;D1])([C;D1])[C;D4;x1;z1]([C;D1])([C;D1])[C;D1]'),
        [1],
        [],
        'CC(C)O[Si](C)(C)C(C)(C)C',
        'CC(C)O',
        ['CC(C)O[SiH](C)C', 'CC(C)O[Si](C)(C)OC']
    )

    rules['hydroxyl_tips'] = (
        smarts('[O;D2:1]-;!@[Si;D4;z1;x1]([C;D3;z1;x1]([C;D1])[C;D1])([C;D3;z1;x1]([C;D1])[C;D1])[C;D3;z1;x1]([C;D1])[C;D1]'),
        [1],
        [],
        'CC(C)O[Si](C(C)C)(C(C)C)C(C)C',
        'CC(C)O',
        ['CC(C)O[SiH](C)C', 'CC(C)O[Si](C)(C)OC']
    )

    rules['hydroxyl_tbdps'] = (
        smarts('[O;D2:1]-;!@[Si;D4;z1;x1]([C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)([C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)[C;D4;x1;z1]([C;D1])([C;D1])[C;D1]'),
        [1],
        [],
        'CC(C)O[Si](c1ccccc1)(c1ccccc1)C(C)(C)C',
        'CC(C)O',
        ['CC(C)O[SiH](C)C', 'CC(C)O[Si](c1ccc(C)cc1)(c1ccccc1)C(C)(C)C']
    )

    rules['hydroxyl_o_nitrobenzyl'] = (
        smarts('[O;D2:1]-;!@[C;D2;z1;x1]-[C;a;r6]:1:[C;D3;x1]([N+](=O)[O-]):[C;D2]:[C;D2]:[C;D2]:[C;D2]:1'),
        [1],
        [],
        'CC(C)OCc1c(N(=O)=O)cccc1',
        'CC(C)O',
        ['CC(C)OC(OC)c1c(N(=O)=O)cccc1']
    )

    rules['hydroxyl_methoxy_benzyl'] = (
        smarts('[O;D2:1]-;!@[C;D2;z1;x1]-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D2]:1'),
        [1],
        [],
        'CC(C)OCc1ccc(OC)cc1',
        'CC(C)O',
        ['CC(C)OCc1ccc(OCC)cc1', 'CC(C)OCc1cc(OC)ccc1']
    )

    rules['hydroxyl_dimethoxybenzyl'] = (
        smarts('[O;D2:1]-;!@[C;D2;z1;x1]-[C;a;r6]:1:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D2]:1'),
        [1],
        [],
        'CC(C)OCc1c(OC)cc(OC)cc1',
        'CC(C)O',
        []
    )

    rules['hydroxyl_naphthyl'] = (
        smarts('[O;D2:1]-;!@[C;D2;z1;x1]-[C;a;r6]:1:[C;D2]:[C;D3]:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D3]:2:[C;D2]:[C;D2]:1'),
        [1],
        [],
        'CC(C)OCC1=CC2=C(C=CC=C2)C=C1',
        'CC(C)O',
        ['CC(C)OCc1ccccc1']
    )

    rules['hydroxyl_bom'] = (
        smarts('[O;D2:1]-;!@[C;D2;x2;z1][O;D2;x0][C;D2;z1;x1]-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1'),
        [1],
        [],
        'CC(C)OCOCc1ccccc1',
        'CC(C)O',
        []
    )

    rules['hydroxyl_piv'] = (
        smarts('[O;D2:1]-;!@[C;z2;x2](=O)-[C;D4;x0;z1]([C;D1])([C;D1])[C;D1]'),
        [1],
        [],
        'CC(C)OC(=O)C(C)(C)C',
        'CC(C)O',
        []
    )

    rules['hydroxyl_methoxy_benzoate'] = (
        smarts('[O;D2:1]-;!@[C;z2;x2](=O)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D2]:1'),
        [1],
        [],
        'COC1=CC=C(C=C1)C(=O)OC(C)C',
        'CC(C)O',
        ['C1=CC=C(C=C1)C(=O)OC(C)C']
    )

    rules['hydroxyl_benzoate'] = (
        smarts('[O;D2:1]-;!@[C;z2;x2](=O)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1'),
        [1],
        [],
        'C1=CC=C(C=C1)C(=O)OC(C)C',
        'CC(C)O',
        ['COC1=CC=C(C=C1)C(=O)OC(C)C']
    )

    rules['hydroxyl_tfa'] = (
        smarts('[O;D2:1]-;!@[C;z2;x2](=O)-[C;D4;z1;x3](F)(F)F'),
        [1],
        [],
        'CC(C)OC(=O)C(F)(F)F',
        'CC(C)O',
        []
    )

    rules['hydroxyl_mom'] = (
        smarts('[O;D2:1]-;!@[C;D2;x2;z1][O;D2;x0][C;D1]'),
        [1],
        [],
        'CC(C)OCOC',
        'CC(C)O',
        ['CC(C)OC(C)OC']
    )

    rules['hydroxyl_mem'] = (
        smarts('[O;D2:1]-;!@[C;D2;x2;z1][O;D2;x0][C;D2;z1;x1][C;D2;z1;x1][O;D2;x0][C;D1]'),
        [1],
        [],
        'COCCOCOC(C)C',
        'CC(C)O',
        []
    )

    rules['hydroxyl_thp'] = (
        smarts('[O;D2:1]-;!@[C;D3;x2;z1;r6]1[O;D2][C;D2][C;D2][C;D2][C;D2]1'),
        [1],
        [],
        'CC(C)OC1CCCCO1',
        'CC(C)O',
        []
    )

    rules['hydroxyl_ee'] = (
        smarts('[O;D2:1]-;!@[C;D3;x2;z1]([O;D2;x0][C;D2;x1;z1][C;D1])[C;D1]'),
        [1],
        [],
        'CC(C)OC(C)OCC',
        'CC(C)O',
        ['CC(C)OC(CC)OCC']
    )

    rules['hydroxyl_mop'] = (
        smarts('[O;D2:1]-;!@[C;D4;x2;z1]([O;D2;x0][C;D1])([C;D1])[C;D1]'),
        [1],
        [],
        'CC(C)OC(C)(C)OC',
        'CC(C)O',
        []
    )

    rules['hydroxyl_sem'] = (
        smarts('[O;D2:1]-;!@[C;D2;x2;z1][O;D2;x0][C;D2;z1;x1][C;D2;z1;x1][Si;D4;z1;x0]([C;D1])([C;D1])[C;D1]'),
        [1],
        [],
        'CC(C)OCOCC[Si](C)(C)C',
        'CC(C)O',
        []
    )

    rules['hydroxyl_tritil'] = (
        smarts('[O;D2:1]-;!@[C;D4;z1;x1](-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)(-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1'),
        [1],
        [],
        'CC(C)OC(c1ccccc1)(c1ccccc1)c1ccccc1',
        'CC(C)O',
        ['COc1ccc(cc1)C(OC(C)C)(c1ccccc1)c1ccc(OC)cc1']
    )

    rules['hydroxyl_dimetoxy_tritil'] = (
        smarts('[O;D2:1]-;!@[C;D4;z1;x1](-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D2]:1)(-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D2]:1)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1'),
        [1],
        [],
        'COc1ccc(cc1)C(OC(C)C)(c1ccccc1)c1ccc(OC)cc1',
        'CC(C)O',
        ['CC(C)OC(c1ccccc1)(c1ccccc1)c1ccccc1']
    )

    rules['hydroxyl_chloro_tritil'] = (
        smarts('[O;D2:1]-;!@[C;D4;z1;x1](-[C;a;r6]:1:[C;D3;x1]([Cl;D1]):[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)(-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1'),
        [1],
        [],
        'CC(C)OC(c1c(Cl)cccc1)(c1ccccc1)c1ccccc1',
        'CC(C)O',
        []
    )

    rules['hydroxyl_mmt'] = (
        smarts('[O;D2:1]-;!@[C;D4;z1;x1](-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D2]:1)(-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1'),
        [1],
        [],
        'c1ccc(cc1)C(OC(C)C)(c1ccccc1)c1ccc(OC)cc1',
        'CC(C)O',
        ['CC(C)OC(c1ccccc1)(c1ccccc1)c1ccccc1']
    )

    rules['hydroxyl_mpe'] = (
        smarts('[O;D2:1]-;!@[C;D4;x1;z1]([C;D1])([C;D2;x0;z1][C;D1])[C;D2;x0;z1][C;D1]'),
        [1],
        [],
        'CC(C)OC(CC)(CC)C',
        'CC(C)O',
        []
    )

    rules['hydroxyl_trifluoroethyl'] = (
        smarts('[O;D2:1]-;!@[C;D2;x1;z1][C;D4;x3;z1](F)(F)F'),
        [1],
        [],
        'CC(C)OCC(F)(F)F',
        'CC(C)O',
        []
    )

    rules['hydroxyl_dmab_enamine'] = (
        smarts('[O;D2:1]-;!@[C;D2;x1;z1]-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D3;x1](:[C;D2]:[C;D2]:1)-[N;D2;x0;z1]-[C;z2;x1;D3]([C;D2;x0;z1][C;D3;x0;z1]([C;D1])[C;D1])=[C;D3;r6;x0;z2]1[C;x1;z2;D3](=O)[C;D2][C;D4;x0;z1]([C;D1])([C;D1])[C;D2][C;D3;x1;z2]1=O'),
        [1],
        [],
        'CC(C)CC(NC1=CC=C(COC(C)=O)C=C1)=C1C(=O)CC(C)(C)CC1=O',
        'CC(O)=O',
        []
    )

    rules['hydroxyl_dmab_imine'] = (
        smarts('[O;D2:1]-;!@[C;D2;x1;z1]-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D3;x1](:[C;D2]:[C;D2]:1)-[N;D2;x0;z2]=[C;x1;D3]([C;D2;x0;z1][C;D3;x0;z1]([C;D1])[C;D1])-[C;D3;r6;x0;z1]1[C;x1;z2;D3](=O)[C;D2][C;D4;x0;z1]([C;D1])([C;D1])[C;D2][C;D3;x1;z2]1=O'),
        [1],
        [],
        'CC(C)CC(=NC1=CC=C(COC(C)=O)C=C1)C1C(=O)CC(C)(C)CC1=O',
        'CC(O)=O',
        []
    )

    rules['diol_12_acetone'] = (
        smarts('[O;D2;x0;r5:1]1-;@[C;D4;x2;z1]([C;D1])([C;D1])-[O;D2;x0:2][C:3]!#[C:4]1'),
        [1, 2, 3, 4],
        [],
        'CC1COC(C)(C)O1',
        'CC(O)CO',
        []
    )

    rules['diol_13_acetone'] = (
        smarts('[O;D2;x0;r6:1]1-;@[C;D4;x2;z1]([C;D1])([C;D1])-[O;D2;x0:2][C:3][C:4]!#[C:5]1'),
        [1, 2, 3, 4, 5],
        [],
        'CC1CCOC(C)(C)O1',
        'CC(O)CCO',
        []
    )

    rules['hydroxyl_amine_acetone'] = (
        smarts('[O;D2;x0;r5:1]1-;@[C;D4;x2;z1]([C;D1])([C;D1])-[N;z1:2][C:3]!#[C:4]1'),
        [1, 2, 3, 4],
        [],
        'CC1CN(C(C)=O)C(C)(C)O1',
        'CC(O)CNC(C)=O',
        []
    )

    rules['diol_12_formalin'] = (
        smarts('[O;D2;x0;r5:1]1-;@[C;D2;x2;z1]-[O;D2;x0:2][C:3]!#[C:4]1'),
        [1, 2, 3, 4],
        [],
        'CC1COCO1',
        'CC(O)CO',
        []
    )

    rules['diol_13_formalin'] = (
        smarts('[O;D2;x0;r6:1]1-;@[C;D2;x2;z1]-[O;D2;x0:2][C:3][C:4]!#[C:5]1'),
        [1, 2, 3, 4, 5],
        [],
        'CC1CCOCO1',
        'CC(O)CCO',
        []
    )

    rules['diol_12_cyclopentanone'] = (
        smarts('[O;D2;x0;r5:1]1-;@[C;D4;x2;z1]2([C;D2;x0;z1][C;D2;x0;z1][C;D2;x0;z1][C;D2;x0;z1]2)-[O;D2;x0:2][C:3]!#[C:4]1'),
        [1, 2, 3, 4],
        [],
        'CC1COC2(CCCC2)O1',
        'CC(O)CO',
        []
    )

    rules['diol_13_cyclopentanone'] = (
        smarts('[O;D2;x0;r6:1]1-;@[C;D4;x2;z1]2([C;D2;x0;z1][C;D2;x0;z1][C;D2;x0;z1][C;D2;x0;z1]2)-[O;D2;x0:2][C:3][C:4]!#[C:5]1'),
        [1, 2, 3, 4, 5],
        [],
        'CC1CCOC2(CCCC2)O1',
        'CC(O)CCO',
        []
    )

    rules['diol_12_cyclohexanone'] = (
        smarts('[O;D2;x0;r5:1]1-;@[C;D4;x2;z1]2([C;D2;x0;z1][C;D2;x0;z1][C;D2;x0;z1][C;D2;x0;z1][C;D2;x0;z1]2)-[O;D2;x0:2][C:3]!#[C:4]1'),
        [1, 2, 3, 4],
        [],
        'CC1COC2(CCCCC2)O1',
        'CC(O)CO',
        []
    )

    rules['diol_13_cyclohexanone'] = (
        smarts('[O;D2;x0;r6:1]1-;@[C;D4;x2;z1]2([C;D2;x0;z1][C;D2;x0;z1][C;D2;x0;z1][C;D2;x0;z1][C;D2;x0;z1]2)-[O;D2;x0:2][C:3][C:4]!#[C:5]1'),
        [1, 2, 3, 4, 5],
        [],
        'CC1CCOC2(CCCCC2)O1',
        'CC(O)CCO',
        []
    )

    rules['diol_12_diacetal'] = (
        smarts('[O;D2;x0;r6:1]1-;@[C;D4;x2;z1]([O;D2;x0][C;D1])([C;D1])[C;D4;x2;z1]([O;D2;x0][C;D1])([C;D1])-[O;D2;x0:2][C:3]!#[C:4]1'),
        [1, 2, 3, 4],
        [],
        'COC1(C)OCC(C)OC1(C)OC',
        'CC(O)CO',
        []
    )

    rules['diol_13_diacetal'] = (
        smarts('[O;D2;x0;r7:1]1-;@[C;D4;x2;z1]([O;D2;x0][C;D1])([C;D1])[C;D4;x2;z1]([O;D2;x0][C;D1])([C;D1])-[O;D2;x0:2][C:3][C:4]!#[C:5]1'),
        [1, 2, 3, 4, 5],
        [],
        'COC1(C)OCCC(C)OC1(C)OC',
        'CC(O)CCO',
        []
    )

    rules['diol_12_benzylidene'] = (
        smarts('[O;D2;x0;r5:1]1-;@[C;D3;x2;z1]([C;a;r6]:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:2)-[O;D2;x0:2][C:3]!#[C:4]1'),
        [1, 2, 3, 4],
        [],
        'CC1COC(O1)c1ccccc1',
        'CC(O)CO',
        []
    )

    rules['diol_13_benzylidene'] = (
        smarts('[O;D2;x0;r6:1]1-;@[C;D3;x2;z1]([C;a;r6]:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:2)-[O;D2;x0:2][C:3][C:4]!#[C:5]1'),
        [1, 2, 3, 4, 5],
        [],
        'CC1CCOC(O1)c1ccccc1',
        'CC(O)CCO',
        []
    )

    rules['carbonyl_dithiolane'] = (
        smarts('[C;D3,D4;z1;x2;r5:1]1[S;D2;x0;z1][C;D2;x1;z1][C;D2;x1;z1][S;D2;x0;z1]1'),
        [1],
        [(1, 'O', 2)],
        'CC1SCCS1',
        'CC=O',
        []
    )

    rules['carbonyl_dithiane'] = (
        smarts('[C;D3,D4;z1;x2;r6:1]1[S;D2;x0;z1][C;D2;x1;z1][C;D2;x0;z1][C;D2;x1;z1][S;D2;x0;z1]1'),
        [1],
        [(1, 'O', 2)],
        'CC1SCCCS1',
        'CC=O',
        []
    )

    rules['carbonyl_dimethylsulfide'] = (
        smarts('[C;D3,D4;z1;x2:1](-;!@[S;D2][C;D1])-;!@[S;D2][C;D1]'),
        [1],
        [(1, 'O', 2)],
        'CSC(CC)SC',
        'CCC=O',
        []
    )

    rules['carbonyl_dioxolane'] = (
        smarts('[C;D3,D4;z1;x2;r5:1]1[O;D2;x0][C;D2;x1;z1][C;D2;x1;z1][O;D2;x0]1'),
        [1],
        [(1, 'O', 2)],
        'CC1OCCO1',
        'CC=O',
        []
    )

    rules['carbonyl_dioxane'] = (
        smarts('[C;D3,D4;z1;x2;r6:1]1[O;D2;x0][C;D2;x1;z1][C;D2;x0;z1][C;D2;x1;z1][O;D2;x0]1'),
        [1],
        [(1, 'O', 2)],
        'CC1OCCCO1',
        'CC=O',
        []
    )

    rules['carbonyl_dimethoxy'] = (
        smarts('[C;D3,D4;z1;x2:1](-;!@[O;D2;x0][C;D1])-;!@[O;D2;x0][C;D1]'),
        [1],
        [(1, 'O', 2)],
        'COC(C)OC',
        'CC=O',
        []
    )

    rules['carboxyl_trioxabicyclooctane'] = (
        smarts('[C;D4;x3;r6:1]12-;@[O;D2][C;D2;x1;z1][C;D4;x0;z1]([C;D1])([C;D2;x1;z1][O;D2]1)[C;D2;x1;z1][O;D2]2'),
        [1],
        [(1, 'O', 2), (1, 'O', 1)],
        'CC(C)C12OCC(C)(CO1)CO2',
        'CC(C)C(O)=O',
        ['CC(C)C12OCC(CC)(CO1)CO2', 'CC(C)C12OC(C)C(C)(CO1)CO2']
    )

    rules['amine_methylcarbamate'] = (
        smarts('[N;D2,D3:1]-;!@[C;z2;x3](=O)[O;D2;x0][C;D1]'),
        [1],
        [],
        'c1ccccc1NC(=O)OC',
        'c1ccccc1N',
        ['c1ccccc1NC(=O)OCC']
    )

    rules['amine_ethylcarbamate'] = (
        smarts('[N;D2,D3:1]-;!@[C;z2;x3](=O)[O;D2;x0][C;D2;x1;z1][C;D1]'),
        [1],
        [],
        'c1ccccc1NC(=O)OCC',
        'c1ccccc1N',
        []
    )

    rules['amine_alloc'] = (
        smarts('[N;D2,D3:1]-;!@[C;z2;x3](=O)[O;D2;x0]-[C;D2;z1;x1][C;D2;x0;z2]=[C;D1]'),
        [1],
        [],
        'c1ccccc1NC(=O)OCC=C',
        'c1ccccc1N',
        ['c1ccccc1NC(=O)OCC=CC']
    )

    rules['amine_teoc'] = (
        smarts('[N;D2,D3:1]-;!@[C;z2;x3](=O)[O;D2;x0]-[C;D2;z1;x1][C;D2;x1;z1][Si;D4;z1;x0]([C;D1])([C;D1])[C;D1]'),
        [1],
        [],
        'c1ccccc1NC(=O)OCC[Si](C)(C)C',
        'c1ccccc1N',
        []
    )

    rules['amine_sem'] = (
        smarts('[N;D2,D3:1]-;!@[C;D2;x2;z1][O;D2;x0]-[C;D2;z1;x1][C;D2;x1;z1][Si;D4;z1;x0]([C;D1])([C;D1])[C;D1]'),
        [1],
        [],
        'CN(C)COCC[Si](C)(C)C',
        'CNC',
        []
    )

    rules['amine_troc'] = (
        smarts('[N;D2,D3:1]-;!@[C;z2;x3](=O)[O;D2;x0]-[C;D2][C;D4;x3]([Cl;D1])([Cl;D1])[Cl;D1]'),
        [1],
        [],
        'c1ccccc1NC(=O)OCC(Cl)(Cl)Cl',
        'c1ccccc1N',
        ['c1ccccc1NC(=O)OC(C)C(Cl)(Cl)Cl']
    )

    rules['amine_cbz'] = (
        smarts('[N;D2,D3:1]-;!@[C;z2;x3](=O)-[O;D2;x0][C;D2;x1;z1][C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1'),
        [1],
        [],
        'c1ccccc1NC(=O)OCc2ccccc2',
        'c1ccccc1N',
        ['c1ccccc1NC(=O)OC(C)c2ccccc2']
    )

    rules['amine_chloro_cbz'] = (
        smarts('[N;D2,D3:1]-;!@[C;z2;x3](=O)-[O;D2;x0][C;D2;x1;z1][C;a;r6]:1:[C;D3;x1]([Cl;D1]):[C;D2]:[C;D2]:[C;D2]:[C;D2]:1'),
        [1],
        [],
        'Clc1ccccc1COC(=O)Nc1ccccc1',
        'c1ccccc1N',
        ['c1ccccc1NC(=O)OC(C)c2ccccc2']
    )

    rules['amine_phenylsulfonyl'] = (
        smarts('[N;D2,D3:1]-;!@[S;D4;x3](=O)(=O)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1'),
        [1],
        [],
        'c1ccccc1S(=O)(=O)Nc1ccccc1',
        'c1ccccc1N',
        []
    )

    rules['amine_tosyl'] = (
        smarts('[N;D2,D3:1]-;!@[S;D4;x3](=O)(=O)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D3;x0]([C;D1]):[C;D2]:[C;D2]:1'),
        [1],
        [],
        'Cc1ccc(cc1)S(=O)(=O)Nc1ccccc1',
        'c1ccccc1N',
        []
    )

    rules['amine_nosyl'] = (
        smarts('[N;D2,D3:1]-;!@[S;D4;x3](=O)(=O)-[C;a;r6]:1:[C;D3;x1]([N+](=O)[O-]):[C;D2]:[C;D2]:[C;D2]:[C;D2]:1'),
        [1],
        [],
        '[O-][N+](=O)c1ccccc1S(=O)(=O)Nc1ccccc1',
        'c1ccccc1N',
        []
    )

    rules['amine_boc'] = (
        smarts('[N;D2,D3:1]-;!@[C;z2;x3](=O)-[O;D2;x0]-[C;D4;x1;z1]([C;D1])([C;D1])[C;D1]'),
        [1],
        [],
        'c1ccccc1NC(=O)OC(C)(C)C',
        'c1ccccc1N',
        []
    )

    rules['amine_tfa'] = (
        smarts('[N;D2,D3:1]-;!@[C;z2;x2](=O)-[C;D4;z1;x3](F)(F)F'),
        [1],
        [],
        'CNC(=O)C(F)(F)F',
        'CN',
        []
    )

    rules['amine_fmoc'] = (
        smarts('[N;D2,D3:1]-;!@[C;z2;x3](=O)[O;D2;x0]-[C;D2;x1;z1][C;D3;z1;x0;r5]1[C;a;r6]:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D3]:2-[C;a;r6]:3:[C;D2]:[C;D2]:[C;D2]:[C;D2]:C1:3'),
        [1],
        [],
        'O=C(Nc1ccccc1)OCC1c2ccccc2-c2ccccc12',
        'c1ccccc1N',
        []
    )

    rules['amine_pbf'] = (
        smarts('[N;D2,D3:1]-;!@[S;D4;x3](=O)(=O)-[C;a;r6]:1:[C;D3;x0]([C;D1]):[C;D3;x0]([C;D1]):[C;D3;x1]:2-[O;D2;x0;r5][C;D4;x1]([C;D1])([C;D1])[C;D2;x0;z1][C;D3]:2:[C;D3;x0]([C;D1]):1'),
        [1],
        [],
        'CN(C)S(=O)(=O)c1c(C)c2CC(C)(C)Oc2c(C)c1C',
        'CNC',
        []
    )

    rules['amine_mtr'] = (
        smarts('[N;D2,D3:1]-;!@[S;D4;x3](=O)(=O)-[C;a;r6]:1:[C;D3;x0]([C;D1]):[C;D3;x0]([C;D1]):[C;D3;x1](-;!@[O;D2;x0][C;D1]):[C;D2]:[C;D3;x0]([C;D1]):1'),
        [1],
        [],
        'COC1=C(C)C(C)=C(C(C)=C1)S(=O)(=O)N(C)C',
        'CNC',
        []
    )

    rules['amine_dde_enamine'] = (
        smarts('[N;D2,D3:1]-;!@[C;z2;x1;D3]([C;D1])=[C;D3;r6;x0;z2]1[C;x1;z2;D3](=O)[C;D2][C;D4;x0;z1]([C;D1])([C;D1])[C;D2][C;D3;x1;z2]1=O'),
        [1],
        [],
        'CC(C)NC(C)=C1C(=O)CC(C)(C)CC1=O',
        'CC(C)N',
        []
    )

    rules['amine_dde_imine'] = (
        smarts('[N;D2:1]=;!@[C;x1;D3]([C;D1])-[C;D3;r6;x0;z1]1[C;x1;z2;D3](=O)[C;D2][C;D4;x0;z1]([C;D1])([C;D1])[C;D2][C;D3;x1;z2]1=O'),
        [1],
        [],
        'CC(C)N=C(C)C1C(=O)CC(C)(C)CC1=O',
        'CC(C)N',
        []
    )

    rules['amine_ivdde_enamine'] = (
        smarts('[N;D2,D3:1]-;!@[C;z2;x1;D3]([C;D2;x0;z1][C;D3;x0;z1]([C;D1])[C;D1])=[C;D3;r6;x0;z2]1[C;x1;z2;D3](=O)[C;D2][C;D4;x0;z1]([C;D1])([C;D1])[C;D2][C;D3;x1;z2]1=O'),
        [1],
        [],
        'CC(C)CC(NC(C)C)=C1C(=O)CC(C)(C)CC1=O',
        'CC(C)N',
        []
    )

    rules['amine_ivdde_imine'] = (
        smarts('[N;D2:1]=;!@[C;x1;D3]([C;D2;x0;z1][C;D3;x0;z1]([C;D1])[C;D1])-[C;D3;r6;x0;z1]1[C;x1;z2;D3](=O)[C;D2][C;D4;x0;z1]([C;D1])([C;D1])[C;D2][C;D3;x1;z2]1=O'),
        [1],
        [],
        'CC(C)CC(=NC(C)C)C1C(=O)CC(C)(C)CC1=O',
        'CC(C)N',
        []
    )

    rules['amine_benzyl'] = (
        smarts('[N;D2,D3:1]-;!@[C;D2;z1;x1]-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1'),
        [1],
        [],
        'CC(C)NCc1ccccc1',
        'CC(C)N',
        []
    )

    rules['amine_methoxy_benzyl'] = (
        smarts('[N;D2,D3:1]-;!@[C;D2;z1;x1]-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D2]:1'),
        [1],
        [],
        'CC(C)NCc1ccc(OC)cc1',
        'CC(C)N',
        []
    )

    rules['amine_dimethoxybenzyl'] = (
        smarts('[N;D2,D3:1]-;!@[C;D2;z1;x1]-[C;a;r6]:1:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D2]:1'),
        [1],
        [],
        'CC(C)NCc1c(OC)cc(OC)cc1',
        'CC(C)N',
        []
    )

    rules['amine_mtt'] = (
        smarts('[N;D2,D3:1]-;!@[C;D4;z1;x1](-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D3;x0]([C;D1]):[C;D2]:[C;D2]:1)(-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1'),
        [1],
        [],
        'CC(C)NC(c1ccccc1)(c1ccccc1)c1ccc(C)cc1',
        'CC(C)N',
        []
    )

    rules['amine_bhoc'] = (
        smarts('[N;D2,D3:1]-;!@[C;z2;x3](=O)[O;D2;x0]-[C;D3;z1;x1](-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1'),
        [1],
        [],
        'CNC(=O)OC(c1ccccc1)c1ccccc1',
        'CN',
        []
    )

    rules['amine_tritil'] = (
        smarts('[N;D2,D3:1]-;!@[C;D4;z1;x1](-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)(-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1'),
        [1],
        [],
        'CC(C)NC(c1ccccc1)(c1ccccc1)c1ccccc1',
        'CC(C)N',
        []
    )

    rules['amine_chloro_tritil'] = (
        smarts('[N;D2,D3:1]-;!@[C;D4;z1;x1](-[C;a;r6]:1:[C;D3;x1]([Cl;D1]):[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)(-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1'),
        [1],
        [],
        'CC(C)NC(c1c(Cl)cccc1)(c1ccccc1)c1ccccc1',
        'CC(C)N',
        []
    )

    rules['thiol_tritil'] = (
        smarts('[S;D2;x0;z1:1]-;!@[C;D4;z1;x1](-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)(-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1'),
        [1],
        [],
        'CC(C)SC(c1ccccc1)(c1ccccc1)c1ccccc1',
        'CC(C)S',
        []
    )

    rules['thiol_mmt'] = (
        smarts('[S;D2;x0;z1:1]-;!@[C;D4;z1;x1](-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D2]:1)(-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1'),
        [1],
        [],
        'c1ccc(cc1)C(SC(C)C)(c1ccccc1)c1ccc(OC)cc1',
        'CC(C)S',
        []
    )

    rules['thiol_dimetoxy_tritil'] = (
        smarts('[S;D2;x0;z1:1]-;!@[C;D4;z1;x1](-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D2]:1)(-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D2]:1)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1'),
        [1],
        [],
        'COc1ccc(cc1)C(SC(C)C)(c1ccccc1)c1ccc(OC)cc1',
        'CC(C)S',
        []
    )

    rules['thiol_chloro_tritil'] = (
        smarts('[S;D2;x0;z1:1]-;!@[C;D4;z1;x1](-[C;a;r6]:1:[C;D3;x1]([Cl;D1]):[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)(-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1'),
        [1],
        [],
        'CC(C)SC(c1c(Cl)cccc1)(c1ccccc1)c1ccccc1',
        'CC(C)S',
        []
    )

    rules['thiol_benzyl'] = (
        smarts('[S;D2;x0;z1:1]-;!@[C;D2;z1;x1]-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1'),
        [1],
        [],
        'CC(C)SCc1ccccc1',
        'CC(C)S',
        []
    )

    rules['thiol_tbu'] = (
        smarts('[S;D2;x0;z1:1]-;!@[C;D4;x1;z1]([C;D1])([C;D1])[C;D1]'),
        [1],
        [],
        'CC(C)SC(C)(C)C',
        'CC(C)S',
        []
    )

    rules['thiol_stbu'] = (
        smarts('[S;D2;x1;z1:1]-;!@[S;D2;z1;x1]-[C;D4;x1;z1]([C;D1])([C;D1])[C;D1]'),
        [1],
        [],
        'CC(C)SSC(C)(C)C',
        'CC(C)S',
        []
    )

    rules['thiol_strimethoxyphenyl'] = (
        smarts('[S;D2;x1;z1:1]-;!@[S;D2;z1;x1]-[C;a;r6]:1:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D3;x1]:1[O;D2;x0][C;D1]'),
        [1],
        [],
        'COc1cc(OC)c(SSC(C)C)c(OC)c1',
        'CC(C)S',
        []
    )

    rules['thiol_amine_dimethoxybenzyl'] = (
        smarts('[S;D2;r5;x0;z1:1]1[C:3][C:4][N;z1:2]-[C;D3;x2;z1]1-[C;a;r6]:1:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D3;x1]([O;D2;x0][C;D1]):[C;D2]:[C;D2]:1'),
        [1, 2, 3, 4],
        [],
        'COC1=CC=C(C2NC(C)CS2)C(OC)=C1',
        'NC(C)CS',
        []
    )

    # General patterns that are subsets of more specific ones above.
    # Must come last to prevent false matches due to atom overlap.

    rules['hydroxyl_tbu'] = (
        smarts('[O;D2:1]-;!@[C;D4;x1;z1]([C;D1])([C;D1])[C;D1]'),
        [1],
        [],
        'CC(C)OC(C)(C)C',
        'CC(C)O',
        []
    )

    rules['hydroxyl_allyl'] = (
        smarts('[O;D2:1]-;!@[C;D2;z1;x1][C;D2;x0;z2]=[C;D1]'),
        [1],
        [],
        'CC(C)OCC=C',
        'CC(C)O',
        []
    )

    rules['hydroxyl_benzyl'] = (
        smarts('[O;D2:1]-;!@[C;D2;z1;x1]-[C;a;r6]:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1'),
        [1],
        [],
        'CC(C)OCc1ccccc1',
        'CC(C)O',
        ['CC(C)OCc1cccc(C)c1', 'CC(C)OC(C)c1ccccc1']
    )

    rules['hydroxyl_acyl'] = (
        smarts('[O;D2:1]-;!@[C;z2;x2](=O)-[C;D1]'),
        [1],
        [],
        'CC(C)OC(=O)C',
        'CC(C)O',
        ['CC(C)OC(=O)CC']
    )

    rules['hydroxyl_methyl'] = (
        smarts('[O;D2:1]-;!@[C;D1]'),
        [1],
        [],
        'CC(C)OC',
        'CC(C)O',
        ['CC(C)OCC']
    )

    return rules


rules = Proxy(_rules)


__all__ = ['rules']
