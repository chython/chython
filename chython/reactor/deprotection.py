# -*- coding: utf-8 -*-
#
#  Copyright 2022-2024 Ramil Nugmanov <nougmanoff@protonmail.com>
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

_alcohol_thiocarbamate = (  # NaIO4 or H2O2/NaOH
    ('[C;D2,D3,D4;z1;x1:1][O:2]C(=[S;D1])N([C;D1])[C;D1]', '[A:1][A:2]',  # rule
     'CC(C)OC(=S)N(C)C', 'CC(C)O'),  # test
)

_alcohol_fmoc = (  # Et3N pKa ~ 10
    ('[C;D2,D3,D4;z1;x1:1][O:2]C(=O)O[C;D2][C;D3]1C:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:C:2-C:3:[C;D2]:[C;D2]:[C;D2]:[C;D2]:C1:3',  # noqa
     '[A:1][A:2]',
     'CC(C)OC(=O)OCC1C2=CC=CC=C2C2=C1C=CC=C2', 'CC(C)O'),
)

_alcohol_troc = (  # [Zn]
    ('[C;D2,D3,D4;z1;x1:1][O:2]C(=O)O[C;D2]C([Cl;D1])([Cl;D1])[Cl;D1]', '[A:1][A:2]',
     'CC(C)OC(=O)OCC(Cl)(Cl)Cl', 'CC(C)O'),
)

_alcohol_teoc = (  # [F-]
    ('[C;D2,D3,D4;z1;x1:1][O:2]C(=O)O[C;D2][C;D2][Si]([C;D1])([C;D1])[C;D1]', '[A:1][A:2]',
     'CC(C)OC(=O)OCC[Si](C)(C)C', 'CC(C)O'),
)

_alcohol_alloc = (  # [Pd] + NuH
    ('[C;D2,D3,D4;z1;x1:1][O:2]C(=O)O[C;D2][C;D2]=[C;D1]', '[A:1][A:2]',
     'CC(C)OC(=O)OCC=C', 'CC(C)O'),
)

_alcohol_allyl = (  # basic or Metal isomerization + hydrolysis
    ('[C;D2,D3,D4;z1;x1:1][O:2][C;D2][C;D2]=[C;D1]', '[A:1][A:2]',
     'CC(C)OCC=C', 'CC(C)O'),
)

_alcohol_silyl = (  # TMS TES TBS TBDMS TIPS TBDPS: [F-] ion substitution
    ('[C;D2,D3,D4;z1;x1:1][O:2][Si;D4;z1;x1]', '[A:1][A:2]',
     'CC(C)O[Si](C)(C)CC', 'CC(C)O', 'CC(C)O[SiH](C)C', 'CC(C)O[Si](C)(C)OC'),
)

_alcohol_benzyl = (  # [H], ...
    ('[C;D2,D3,D4;z1;x1:1][O:2][C;D2]C:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1', '[A:1][A:2]',
     'CC(C)OCc1ccccc1', 'CC(C)O', 'CC(C)OCc1cccc(C)c1', 'CC(C)OC(C)c1ccccc1'),  # test + decoys
)

_alcohol_o_nitrobenzyl = (  # UV-light
    ('[C;D2,D3,D4;z1;x1:1][O:2][C;D2,D3;z1;x1]C:1:C([N+](=O)[O-]):[C;D2]:[C;D2]:[C;D2]:[C;D2]:1', '[A:1][A:2]',
     'CC(C)OCc1c(N(=O)=O)cccc1', 'CC(C)O', 'CC(C)OC(OC)c1c(N(=O)=O)cccc1'),
)

_alcohol_methoxy_benzyl = (  # PMB or MPM
    ('[C;D2,D3,D4;z1;x1:1][O:2][C;D2]C:1:[C;D2]:[C;D2]:C(O[C;D1]):[C;D2]:[C;D2]:1', '[A:1][A:2]',
     'CC(C)OCc1ccc(OC)cc1', 'CC(C)O', 'CC(C)OCc1ccc(OCC)cc1', 'CC(C)OCc1cc(OC)ccc1'),
)

_alcohol_bom = (  # like Bn
    ('[C;D2,D3,D4;z1;x1:1][O:2][C;D2]O[C;D2]C:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1', '[A:1][A:2]',
     'CC(C)OCOCc1ccccc1', 'CC(C)O'),
)

_alcohol_piv = (
    ('[C;D2,D3,D4;z1;x1:1][O:2]C(=O)C([C;D1])([C;D1])[C;D1]', '[A:1][A:2]'),  # Piv
)

_alcohol_methoxy_benzoate = (
    ('[C;D2,D3,D4;z1;x1:1][O:2]C(=O)C:1:[C;D2]:[C;D2]:C(O[C;D1]):[C;D2]:[C;D2]:1', '[A:1][A:2]'),  # pMeO-Bz
)

_alcohol_benzoate = (
    ('[C;D2,D3,D4;z1;x1:1][O:2]C(=O)C:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1', '[A:1][A:2]'),  # Bz
)

_alcohol_acyl = (
    ('[C;D2,D3,D4;z1;x1:1][O:2]C(=O)[C;D1]', '[A:1][A:2]'),  # Ac
)

_alcohol_tfa = (
    ('[C;D2,D3,D4;z1;x1:1][O:2]C(=O)C(F)(F)F', '[A:1][A:2]'),  # TFA
)

_alcohol_mom = (
    ('[C;D2,D3,D4;z1;x1:1][O:2][C;D2]O[C;D1]', '[A:1][A:2]'),  # MOM
)

_alcohol_thp = (
    ('[C;D2,D3,D4;z1;x1:1][O:2][C;D3]1[O][C;D2][C;D2][C;D2][C;D2]1', '[A:1][A:2]'),  # THP
)

_alcohol_tritil = (
    ('[C;D2,D3,D4;z1;x1:1][O:2]C(C:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)(C:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:2)C:3:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:3',  # noqa
     '[A:1][A:2]'),
)

_alcohol_tbu = (
    ('[C;D2,D3,D4;z1;x1;M][O:1]-C([C;D1])([C;D1])[C;D1]', '[A:1]'),
)

_alcohol_amide_acetone = (  # N1CCOC1
    ('[O:1]1[C;x1;z1;M][C;x1;z1;M][N:2]([C;M]=[O;M])-C1([C;D1])[C;D1]',
     '[A:1].[A:2]'),
)

_diol12_acetone = (
    ('[C;D3,D4;z1;x1:1]1[O:2]C([C;D1])([C;D1])[O:3][C;z1;x1:4]1', '[A:3][A:4][A:1][A:2]'),  # ketone
)

_diol12_formalin = (
    ('[C;D3,D4;z1;x1:1]1[O:2][C;D2][O:3][C;z1;x1:4]1', '[A:3][A:4][A:1][A:2]'),  # formaldehyde
)

_diol12_cyclopentanone = (
    ('[C;D3,D4;z1;x1:1]1[O:2]C2([C;D2][C;D2][C;D2][C;D2]2)[O:3][C;z1;x1:4]1', '[A:3][A:4][A:1][A:2]'),
)

_diol12_cyclohexanone = (
    ('[C;D3,D4;z1;x1:1]1[O:2]C2([C;D2][C;D2][C;D2][C;D2][C;D2]2)[O:3][C;z1;x1:4]1', '[A:3][A:4][A:1][A:2]'),
)

_diol12_diacetal = (
    ('[C;D3,D4;z1;x1:1]1[O:2]-C([C;D1])(O[C;D1])-C([C;D1])(O[C;D1])-[O:3][C;z1;x1:4]1',
     '[A:3][A:4][A:1][A:2]'),
)

_diol13_formalin = (
    ('[C:5]1[C;D3,D4;z1;x1:1][O:2][C;D2][O:3][C;z1;x1:4]1', '[A:3][A:4][A:5][A:1][A:2]'),  # formaldehyde
)

_diol13_acetone = (
    ('[C:5]1[C;D3,D4;z1;x1:1][O:2]C([C;D1])([C;D1])[O:3][C;z1;x1:4]1', '[A:3][A:4][A:5][A:1][A:2]'),
)

_diol12_benzylidene = (
    ('[C;D3,D4;z1;x1:1]1[O:2][C;D3](C:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:2)[O:3][C;z1;x1:4]1',
     '[A:3][A:4][A:1][A:2]',
     'CC1COC(O1)c1ccccc1', 'CC(O)CO'),
)

_diol13_cyclopentanone = (
    ('[C:5]1[C;D3,D4;z1;x1:1][O:2]C2([C;D2][C;D2][C;D2][C;D2]2)[O:3][C;z1;x1:4]1',
     '[A:3][A:4][A:5][A:1][A:2]'),
)

_diol13_cyclohexanone = (
    ('[C:5]1[C;D3,D4;z1;x1:1][O:2]C2([C;D2][C;D2][C;D2][C;D2][C;D2]2)[O:3][C;z1;x1:4]1',
     '[A:3][A:4][A:5][A:1][A:2]'),
)

_diol13_diacetal = (
    ('[C:5]1[C;D3,D4;z1;x1:1][O:2]-C([C;D1])(O[C;D1])-C([C;D1])(O[C;D1])-[O:3][C;z1;x1:4]1',
     '[A:3][A:4][A:5][A:1][A:2]'),
)

_diol13_benzylidene = (
    ('[C:5]1[C;D3,D4;z1;x1:1][O:2][C;D3](C:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:2)[O:3][C;z1;x1:4]1',
     '[A:3][A:4][A:5][A:1][A:2]',
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

_carboxyl_tbu = (
    ('[C;D3;x2:1](=[O:2])[O:3]-C([C;D1])([C;D1])[C;D1]', '[A:1](=[A:2])[A:3]'),
)

_carboxyl_mpe = (
    ('[C;D3;x2;M](=[O;M])[O:1]-C([C;D1])([C;D2][C;D1])[C;D2][C;D1]', '[A:1]'),
)

_carboxyl_methyl = (
    ('[C;D3;x2:1](=[O:2])-[O:4]-[C;D1]', '[A:1](=[A:2])O'),  # Me
)

_carboxyl_trifluoroethyl = (
    ('[C;D3;x2:1](=[O:2])-[O:4]-[C;D2]C(F)(F)F', '[A:1](=[A:2])O'),  # CF3-CH2-
)

_carboxyl_trioxabicyclooctane = (  # [H+]. Note! second step of basic hydrolysis required.
    ('[C;D4;x3:1]12[O:4][C;D2]C([C;D1])([C;D2]O1)[C;D2]O2', '[A:1](=O)O',
     'CC(C)C12OCC(C)(CO1)CO2', 'CC(C)C(O)=O', 'CC(C)C12OCC(CC)(CO1)CO2', 'CC(C)C12OC(C)C(C)(CO1)CO2'),
)

_carboxyl_allyl = (  # [Pd] + NuH
    ('[C;D3;x2:1](=[O:2])[O:3]-[C;D2][C;D2]=[C;D1]', '[A:1](=[A:2])[A:3]'),
)

_carboxyl_benzyl = (  # [H] or Li/NH3
    ('[C;D3;x2:1](=[O:2])[O:3]-[C;D2]C:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1', '[A:1](=[A:2])[A:3]'),
)

_carboxyl_fm = (
    ('[C;x2;M](=[O;M])[O:1][C;D2][C;D3]1C:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:C:2-C:3:[C;D2]:[C;D2]:[C;D2]:[C;D2]:C1:3',
     '[A:1]'),
)

_carboxyl_dmab = (
    ('[C;x2;M](=[O;M])[O:1]-[C;D2]-C:1:[C;D2]:[C;D2]:C(:[C;D2]:[C;D2]:1)-[N;D2]-C([C;D2][C;D3]([C;D1])[C;D1])=C1C(=O)[C;D2]C([C;D1])([C;D1])[C;D2]C1=O',
     '[A:1]'),
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
    ('[C;D2,D3,D4;z1;x1;M][S;D2:1]-C(C:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)(C:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:2)C:3:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:3',  # noqa
     '[A:1]'),
)

_thiol_mmt = (
    ('[C;D2,D3,D4;z1;x1;M][S;D2:1]-C(C:1:[C;D2]:[C;D2]:C(-O[C;D1]):[C;D2]:[C;D2]:1)(C:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:2)C:3:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:3',  # noqa
     '[A:1]'),
)

_thiol_benzyl = (
    ('[C;D2,D3,D4;z1;x1;M][S;D2:1]-[C;D2]C:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1', '[A:1]'),
)

_thiol_tbu = (
    ('[C;D2,D3,D4;z1;x1;M][S;D2:1]-C([C;D1])([C;D1])[C;D1]', '[A:1]'),
)

_thiol_strimethoxyphenyl = (
    ('[C;D2,D3,D4;z1;x1;M][S;D2:1]-[S;D2]-C:1:C(O[C;D1]):[C;D2]:C(O[C;D1]):[C;D2]:C:1-O[C;D1]', '[A:1]'),
)

_thiol_stbu = (
    ('[C;D2,D3,D4;z1;x1;M][S;D2:1]-[S;D2]-C([C;D1])([C;D1])[C;D1]', '[A:1]'),
)

_thiol_amide_dimethoxybenzyl = (  # N1CCSC1
    ('[S;D2:1]1[C;x1;z1;M][C;x1;z1;M][N:2]([C;M]=[O;M])-[C;D3]1-C:2:C(O[C;D1]):[C;D2]:C(O[C;D1]):[C;D2]:[C;D2]:2',
     '[A:1].[A:2]'),
)

_pyrrole_boc = (
    ('[N;a;r5:1]-C(=O)OC([C;D1])([C;D1])[C;D1]', '[A:1]'),
)

_pyrrole_chloro_tritil = (
    ('[N;a;r5:1]-C(C:1:C([Cl;D1]):[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)(C:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:2)C:3:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:3',
     '[A:1]'),
)

_phenol_tbu = (
    ('[C;a;M][O:1]-C([C;D1])([C;D1])[C;D1]', '[A:1]'),
)

_phenol_hydroxymethyl_acetone = (
    ('[C;M]1:[C;M][O:1]C([C;D1])([C;D1])[O:2][C;z1;x1;M]1', '[A:1].[A:2]'),
)

_phosphate_benzyl = (
    ('[O;M][P;M](=[O;M])([O;M])[O:1]-[C;D2]C:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1', '[A:1]'),
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
