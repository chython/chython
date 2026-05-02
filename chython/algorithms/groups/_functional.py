# -*- coding: utf-8 -*-
#
#  Copyright 2024-2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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

    # carbohydrides
    rules['terminal_alkene'] = smarts('[C;z2;x0;D1:1]=[C;z2;x0;D2,D3:2]')
    rules['alkene'] = smarts('[C;z2;x0;D2,D3:1]=[C;z2;x0;D2,D3:2]')
    rules['terminal_alkyne'] = smarts('[C;z3;x0;D1:1]#[C;x0;D2:2]')
    rules['alkyne'] = smarts('[C;z3;x0;D2:1]#[C;x0;D2:2]')

    # halides
    rules['aryl_fluoride'] = smarts('[F;D1:100]-[C;a:1]')
    rules['aryl_chloride'] = smarts('[Cl;D1:100]-[C;a:1]')
    rules['aryl_bromide'] = smarts('[Br;D1:100]-[C;a:1]')
    rules['aryl_iodide'] = smarts('[I;D1:100]-[C;a:1]')

    rules['alkyl_fluoride'] = smarts('[F;D1:100][C;z1;x1:1]')
    rules['alkyl_chloride'] = smarts('[Cl;D1:100][C;z1;x1:1]')
    rules['alkyl_bromide'] = smarts('[Br;D1:100][C;z1;x1:1]')
    rules['alkyl_iodide'] = smarts('[I;D1:100][C;z1;x1:1]')

    rules['alkenyl_fluoride'] = smarts('[F;D1:100][C;z2;x1:1]=[C:2]')
    rules['alkenyl_chloride'] = smarts('[Cl;D1:100][C;z2;x1:1]=[C:2]')
    rules['alkenyl_bromide'] = smarts('[Br;D1:100][C;z2;x1:1]=[C:2]')
    rules['alkenyl_iodide'] = smarts('[I;D1:100][C;z2;x1:1]=[C:2]')

    rules['alkynyl_fluoride'] = smarts('[F;D1:100][C;z3;x1:1]')
    rules['alkynyl_chloride'] = smarts('[Cl;D1:100][C;z3;x1:1]')
    rules['alkynyl_bromide'] = smarts('[Br;D1:100][C;z3;x1:1]')
    rules['alkynyl_iodide'] = smarts('[I;D1:100][C;z3;x1:1]')

    # pseudohalides
    rules['aryl_triflate'] = smarts('[S;D4](=O)(=O)(-[O:100]-;!@[C;a:1])-[C;D4](F)(F)F')
    rules['aryl_mesylate'] = smarts('[S;D4](=O)(=O)(-[O:100]-;!@[C;a:1])-[C;D1]')
    rules['aryl_tosylate'] = smarts('[S;D4](=O)(=O)(-[O:100]-;!@[C;a:1])-[C;a]:1:[C;D2]:[C;D2]:[C](-[C;D1]):[C;D2]:[C;D2]:1')

    rules['alkyl_triflate'] = smarts('[S;D4](=O)(=O)(-[O:100]-;!@[C;z1;x1:1])-[C;D4](F)(F)F')
    rules['alkyl_mesylate'] = smarts('[S;D4](=O)(=O)(-[O:100]-;!@[C;z1;x1:1])-[C;D1]')
    rules['alkyl_tosylate'] = smarts('[S;D4](=O)(=O)(-[O:100]-;!@[C;z1;x1:1])-[C;a]:1:[C;D2]:[C;D2]:[C](-[C;D1]):[C;D2]:[C;D2]:1')

    # boronic acids and esters
    rules['aryl_boronic_acid'] = smarts('[B;D3;z1;x2:100](-[O;D1])(-[O;D1])-;!@[C;a:1]')
    rules['aryl_boronic_ester'] = smarts('[B;D3;z1;x2:100](-[O;D2;x1])(-[O;D2;x1])-;!@[C;a:1]')

    rules['alkyl_boronic_acid'] = smarts('[B;D3;z1;x2:100](-[O;D1])(-[O;D1])-;!@[C;z1;x1:1]')
    rules['alkyl_boronic_ester'] = smarts('[B;D3;z1;x2:100](-[O;D2;x1])(-[O;D2;x1])-;!@[C;z1;x1:1]')

    rules['alkenyl_boronic_acid'] = smarts('[B;D3;z1;x2:100](-[O;D1])(-[O;D1])-;!@[C;z2;x1:1]=[C:2]')
    rules['alkenyl_boronic_ester'] = smarts('[B;D3;z1;x2:100](-[O;D2;x1])(-[O;D2;x1])-;!@[C;z2;x1:1]=[C:2]')

    rules['alkynyl_boronic_acid'] = smarts('[B;D3;z1;x2:100](-[O;D1])(-[O;D1])-;!@[C;z3;x1:1]')
    rules['alkynyl_boronic_ester'] = smarts('[B;D3;z1;x2:100](-[O;D2;x1])(-[O;D2;x1])-;!@[C;z3;x1:1]')

    # molander salts (trifluoroborates)
    rules['aryl_molander_salt'] = smarts('[B;D4;z1;x3;-:100](F)(F)(F)-;!@[C;a:1]')
    rules['alkyl_molander_salt'] = smarts('[B;D4;z1;x3;-:100](F)(F)(F)-;!@[C;z1:1]')
    rules['alkenyl_molander_salt'] = smarts('[B;D4;z1;x3;-:100](F)(F)(F)-;!@[C;z2:1]=[C:2]')
    rules['alkynyl_molander_salt'] = smarts('[B;D4;z1;x3;-:100](F)(F)(F)-;!@[C;z3;x1:1]')

    # alcohols and phenols
    rules['primary_alcohol'] = smarts('[O;D1;z1;x0:1][C;D2;x1;z1:2]')
    rules['secondary_alcohol'] = smarts('[O;D1;z1;x0:1][C;D3;x1;z1:2]')
    rules['tertiary_alcohol'] = smarts('[O;D1;z1;x0:1][C;D4;x1;z1:2]')
    rules['phenol'] = smarts('[O;D1;z1;x0:1]-[C;a:2]')

    # aldehydes and ketones
    rules['aldehyde'] = smarts('[O;z2;x0:2]=[C;D2;x1;z2:1]')
    rules['ketone'] = smarts('[O;z2;x0:2]=[C;D3;x1;z2:1]')

    # acids
    rules['alkyl_carboxylic_acid'] = smarts('[O;D1;z1;x0:100][C;z2;x2;D3:1](=[O:2])[C;z1:3]')
    rules['aryl_carboxylic_acid'] = smarts('[O;D1;z1;x0:100][C;z2;x2;D3:1](=[O:2])[C;a:3]')
    rules['carboxylic_acid'] = smarts('[O;D1;z1;x0:100][C;z2;x2;D3:1]=[O:2]')
    rules['acyl_chloride'] = smarts('[Cl:100][C;z2;x2;D3:1]=[O:2]')
    rules['acyl_fluoride'] = smarts('[F:100][C;z2;x2;D3:1]=[O:2]')

    # amines
    rules['primary_amine'] = smarts('[N;D1;z1;x0:1][C;z1:2]')
    rules['primary_aniline'] = smarts('[N;D1;z1;x0:1][C;a:2]')
    rules['secondary_amine'] = smarts('[N;D2;z1;x0:1]([C;z1:2])[C;z1:3]')
    rules['secondary_aniline'] = smarts('[N;D2;z1;x0:1]([C;a:2])[C;z1:3]')
    rules['biaryl_aniline'] = smarts('[N;D2;z1;x0:1]([C;a:2])[C;a:3]')

    # esters and amides
    rules['ester'] = smarts('[O;z2;x0:2]=[C;D3;x2;z2:1]-[O;D2;x0:100]')
    rules['primary_amide'] = smarts('[N;D1;z1;x0:1][C;z2;x2;D3:2]=[O:3]')
    rules['secondary_amide'] = smarts('[N;D2;z1;x0:1][C;z2;x2;D3:2]=[O:3]')

    # sulfonyl
    rules['sulfonyl_chloride'] = smarts('[S;x3;D4:1](=[O:2])(=[O:3])[Cl;D1:100]')
    rules['sulfonyl_fluoride'] = smarts('[S;x3;D4:1](=[O:2])(=[O:3])[F;D1:100]')
    rules['sulfonamide'] = smarts('[S;x3;D4:1](=[O:2])(=[O:3])-[N;z1:100]')
    rules['sulfonyl_anhydride'] = smarts('[S;x3;D4:1](=[O:2])(=[O:3])-[O:100]-[S;x3;D4](=O)(=O)')

    # nitrogen functional groups
    rules['nitrile'] = smarts('[N;D1;z3;x0:2]#[C;D2;x1:1]')
    rules['azide'] = smarts('[N;x1;D2:1]=[N+:2]=[N-:3]')
    rules['isocyanate'] = smarts('[N;z2;x0;D2:1]=[C:2]=[O:3]')
    rules['isocyano'] = smarts('[N;D2;x0;+:1]#[C;-:2]')
    rules['guanidine'] = smarts('[N;z1;x0:1][C;!R:2]([N;z1;x0:3])=[N;x0:4]')
    rules['nitro'] = smarts('[N;D3;x2;+:1]([O;-:2])=[O:3]')

    # sulfur
    rules['thiol'] = smarts('[S;x0;D1;z1:1][C;z1:2]')

    # pyrrole. for tautomerism handling H not in template.
    rules['pyrrole'] = smarts('[N;h1;D2;a;r5:1]')
    rules['pyrazole'] = smarts('[N;h1;D2;a;r5:1]:[N;h0;D2;r5:2]')
    rules['imidazole'] = smarts('[N;h1;D2;a;r5:1]:[A:2]:[N;h0;D2;r5:3]')

    return rules


rules = Proxy(_rules)


__all__ = ['rules']
