# -*- coding: utf-8 -*-
#
#  Copyright 2022 Ramil Nugmanov <nougmanoff@protonmail.com>
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

###########################
# orthogonal deprotection #
###########################

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
     'CC(C)OCc1ccccc1', 'CC(C)O', 'CC(C)OCc1cccc(C)c1', 'CC(C)OC(C)c1ccccc1'),
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

_diol12_benzylidene = (
    ('[C;D3,D4;z1;x1:1]1[O:2][C;D3](C:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:2)[O:3][C;z1;x1:4]1',
     '[A:3][A:4][A:1][A:2]',
     'CC1COC(O1)c1ccccc1', 'CC(O)CO'),
)

_diol13_benzylidene = (
    ('[C:5]1[C;D3,D4;z1;x1:1][O:2][C;D3](C:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:2)[O:3][C;z1;x1:4]1',
     '[A:3][A:4][A:5][A:1][A:2]',
     'CC1CCOC(O1)c1ccccc1', 'CC(O)CCO'),
)

_carbonyl_dithioacetal = (  # MeI - S methylation + hydrolysis
    ('[C;D3,D4;z1;x2:1]1[S;D2:3][C;D2][C;D2][S;D2]1', '[A:1]=O'),
    ('[C;D3,D4;z1;x2:1]1[S;D2:3][C;D2][C;D2][C;D2][S;D2]1', '[A:1]=O'),
    ('[C;D3,D4;z1;x2:1]([S;D2:3][C;D1])[S;D2][C;D1]', '[A:1]=O'),
)

_carboxyl_allyl = (  # [Pd] + NuH
    ('[C;D3;x2:1](=[O:2])[O:3]-[C;D2][C;D2]=[C;D1]', '[A:1](=[A:2])[A:3]'),
)

_carboxyl_benzyl = (  # [H] or Li/NH3
    ('[C;D3;x2:1](=[O:2])[O:3]-[C;D2]C:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1', '[A:1](=[A:2])[A:3]'),
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
    ('[C;D2,D3,D4;z1;x1:1][N;D2:2]-C(=O)O[C;z1;x1]', '[A:1][A:2]'),
    ('[C;a:1][N:2]([C;z1;x1:3])-C(=O)O[C;z1;x1]', '[A:1][A:2][A:3]'),  # Alk-NH-Ar
    ('[C;D2,D3,D4;z1;x1:1][N:2]([C;z1;x1:3])-C(=O)O[C;z1;x1]', '[A:1][A:2][A:3]'),  # Alk2NH
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

#######################
# acidic deprotection #
#######################

_alcohol_acetal = (
    ('[C;D2,D3,D4;z1;x1:1][O:2][C;D2]O[C;D1]', '[A:1][A:2]'),  # MOM
    ('[C;D2,D3,D4;z1;x1:1][O:2][C;D3]1[O][C;D2][C;D2][C;D2][C;D2]1', '[A:1][A:2]'),  # THP
    ('[C;D2,D3,D4;z1;x1:1][O:2][C;D2]O[C;D2]C:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1', '[A:1][A:2]'),  # BOM
)

_alcohol_tritil = (  # Ph3C-
    ('[C;D2,D3,D4;z1;x1:1][O:2]C(C:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1)(C:2:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:2)C:3:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:3',  # noqa
     '[A:1][A:2]'),
)

_diol12 = (  # acetals based
    ('[C;D3,D4;z1;x1:1]1[O:2][C;D2][O:3][C;z1;x1:4]1', '[A:3][A:4][A:1][A:2]'),  # formaldehyde
    ('[C;D3,D4;z1;x1:1]1[O:2]C([C;D1])([C;D1])[O:3][C;z1;x1:4]1', '[A:3][A:4][A:1][A:2]'),  # ketone
    # cyclopentanone
    ('[C;D3,D4;z1;x1:1]1[O:2]C2([C;D2][C;D2][C;D2][C;D2]2)[O:3][C;z1;x1:4]1', '[A:3][A:4][A:1][A:2]'),
    # cyclohexanone
    ('[C;D3,D4;z1;x1:1]1[O:2]C2([C;D2][C;D2][C;D2][C;D2][C;D2]2)[O:3][C;z1;x1:4]1', '[A:3][A:4][A:1][A:2]'),
    # diacetal
    ('[C;D3,D4;z1;x1:1]1[O:2]-C([C;D1])(O[C;D1])-C([C;D1])(O[C;D1])-[O:3][C;z1;x1:4]1',
     '[A:3][A:4][A:1][A:2]'),
)

_diol13 = (  # acetals based
    ('[C:5]1[C;D3,D4;z1;x1:1][O:2][C;D2][O:3][C;z1;x1:4]1', '[A:3][A:4][A:5][A:1][A:2]'),  # formaldehyde
    # ketone
    ('[C:5]1[C;D3,D4;z1;x1:1][O:2]C([C;D1])([C;D1])[O:3][C;z1;x1:4]1', '[A:3][A:4][A:5][A:1][A:2]'),
    # cyclopentanone
    ('[C:5]1[C;D3,D4;z1;x1:1][O:2]C2([C;D2][C;D2][C;D2][C;D2]2)[O:3][C;z1;x1:4]1',
     '[A:3][A:4][A:5][A:1][A:2]'),
    # cyclohexanone
    ('[C:5]1[C;D3,D4;z1;x1:1][O:2]C2([C;D2][C;D2][C;D2][C;D2][C;D2]2)[O:3][C;z1;x1:4]1',
     '[A:3][A:4][A:5][A:1][A:2]'),
    # diacetal
    ('[C:5]1[C;D3,D4;z1;x1:1][O:2]-C([C;D1])(O[C;D1])-C([C;D1])(O[C;D1])-[O:3][C;z1;x1:4]1',
     '[A:3][A:4][A:5][A:1][A:2]'),
)

_carbonyl = (  # acetals
    ('[C;D3,D4;z1;x2:1]1[O,S;D2:3][C;D2][C;D2][O,S;D2]1', '[A:1]=O'),  # dioxolane
    ('[C;D3,D4;z1;x2:1]1[O,S;D2:3][C;D2][C;D2][C;D2][O,S;D2]1', '[A:1]=O'),  # dioxane
    ('[C;D3,D4;z1;x2:1]([O,S;D2:3][C;D1])[O,S;D2][C;D1]', '[A:1]=O'),  # dimethoxy
)

_carboxyl_tbu = (
    ('[C;D3;x2:1](=[O:2])[O:3]-C([C;D1])([C;D1])[C;D1]', '[A:1](=[A:2])[A:3]'),
)

_carboxyl_trioxabicyclooctane = (  # Note! second step of basic hydrolysis required.
    ('[C;D4;x3:1]12[O:4][C;D2]C([C;D1])([C;D2]O1)[C;D2]O2', '[A:1](=O)O',
     'CC(C)C12OCC(C)(CO1)CO2', 'CC(C)C(O)=O', 'CC(C)C12OCC(CC)(CO1)CO2', 'CC(C)C12OC(C)C(C)(CO1)CO2'),
)

_amine_boc = (
    # Ar-NH2
    ('[C;a:1][N;D2:2]-C(=O)OC([C;D1])([C;D1])[C;D1]', '[A:1][A:2]',
     'c1ccccc1NC(=O)OC(C)(C)C', 'c1ccccc1N'),
    # Alk-NH2
    ('[C;D2,D3,D4;z1;x1:1][N;D2:2]-C(=O)OC([C;D1])([C;D1])[C;D1]', '[A:1][A:2]'),
    # Alk-NH-Ar
    ('[C;a:1][N:2]([C;z1;x1:3])-C(=O)OC([C;D1])([C;D1])[C;D1]', '[A:1][A:2][A:3]'),
    # Alk2NH
    ('[C;D2,D3,D4;z1;x1:1][N:2]([C;z1;x1:3])-C(=O)OC([C;D1])([C;D1])[C;D1]', '[A:1][A:2][A:3]'),
)

######################
# basic deprotection #
######################

_alcohol_ester = (
    ('[C;D2,D3,D4;z1;x1:1][O:2]C(=O)C([C;D1])([C;D1])[C;D1]', '[A:1][A:2]'),  # Piv
    ('[C;D2,D3,D4;z1;x1:1][O:2]C(=O)C:1:[C;D2]:[C;D2]:C(O[C;D1]):[C;D2]:[C;D2]:1', '[A:1][A:2]'),  # pMeO-Bz
    ('[C;D2,D3,D4;z1;x1:1][O:2]C(=O)C:1:[C;D2]:[C;D2]:[C;D2]:[C;D2]:[C;D2]:1', '[A:1][A:2]'),  # Bz
    ('[C;D2,D3,D4;z1;x1:1][O:2]C(=O)[C;D1]', '[A:1][A:2]'),  # Ac
    ('[C;D2,D3,D4;z1;x1:1][O:2]C(=O)C(F)(F)F', '[A:1][A:2]'),  # TFA
)

_carboxyl_basic = (
    ('[C;D3;x2:1](=[O:2])-[O:4]-[C;D1]', '[A:1](=[A:2])O'),  # Me
    ('[C;D3;x2:1](=[O:2])-[O:4]-[C;D2]C(F)(F)F', '[A:1](=[A:2])O'),  # CF3-CH2-
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

#################
# Magic Factory #
#################

__all__ = [k[1:] for k, v in globals().items() if k.startswith('_') and isinstance(v, tuple) and v]

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


def __getattr__(name):
    try:
        return _cache[name]
    except KeyError:
        if name in __all__:
            _cache[name] = t = _prepare_reactor(globals()[f'_{name}'], name)
            return t
        raise AttributeError


def __dir__():
    return __all__
