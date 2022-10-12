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
from lazy_object_proxy import Proxy


def _rules():
    from ... import smiles, smarts
    raw_rules = []

    # phenol etherification
    r = smiles('[OH:5][cH:1]:[cH2:2].[CH3:4][OH:3]>>[CH3:4][O:3][cH:1]:[cH2:2]')
    m = smarts('[O;D1:5].[O;D1:3]')  # reactant restriction.
    raw_rules.append((r, {3: 5}, m))

    # esterification
    r = smiles('[OH:2][CH:1]=[O:3].[CH3:4][OH:5]>>[CH3:4][O:2][CH:1]=[O:3]')
    m = smarts('[O;D1:2].[O;D1:5]')
    raw_rules.append((r, {2: 5}, m))

    # acid reduction
    r = smiles('[CH3:2][C:1]([OH:3])=[O:4]>>[CH3:2][CH2:1][OH:3]')
    m = smarts('[O;D1:3]')
    raw_rules.append((r, {3: 4, 4: 3}, m))

    # ester reduction
    r = smiles('[CH3:5][O:3][C:1]([CH3:2])=[O:4]>>[CH3:2][CH2:1][OH:3]')
    raw_rules.append((r, {3: 4}, None))

    # ozonolysis
    r = smiles('[O-:4][O+:1]=[O:6].[CH3:3][CH:2]=[CH2:5]>>[CH3:3][CH:2]=[O:1]')
    raw_rules.append((r, {1: 4}, None))

    # peroxide oxirane
    r = smiles('[CH3:6][O:1][OH:5].[CH3:4][CH:2]=[CH2:3]>>[CH3:4][CH:2]1[CH2:3][O:1]1')
    m = smarts('[O;D1:5]')
    raw_rules.append((r, {1: 5}, m))

    # claisen rearrangement
    r = smiles('[cH2:4]:[cH:3]:[cH:2][O:1][CH2:5][CH:6]=[CH2:7]>>[OH:1][cH:2]:[c:3](:[cH2:4])[CH2:5][CH:6]=[CH2:7]')
    raw_rules.append((r, {5: 7, 7: 5}, None))

    # metathesis
    r = smiles('[CH3:1][CH:2]=[CH2:3].[CH3:5][CH:4]=[CH2:6]>>[CH3:5][CH:3]=[CH:2][CH3:1]')
    raw_rules.append((r, {3: 4, 4: 3}, None))

    rules = []
    for r, f, m in raw_rules:
        c = ~r
        rules.append((c, str(c.substructure(c.center_atoms)), m, f))
    return rules


rules = Proxy(_rules)


__all__ = ['rules']
