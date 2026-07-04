# -*- coding: utf-8 -*-
#
#  Copyright 2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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


def _make_reactor(rxn_name, fg_name, output_fg, product_smarts):
    from ._functional import rules
    from ...reactor import Reactor
    from ... import smarts

    q = rules[fg_name]
    product = smarts(product_smarts)
    return rxn_name, fg_name, output_fg, Reactor((q,), (product,),
                                                 delete_atoms=True, one_shot=True, ignore_pyrrole_hydrogen=True)


def _rules():
    rules = []

    # aldehyde → primary_alcohol (NaBH4, LiAlH4)
    rules.append(_make_reactor('aldehyde_to_alcohol', 'aldehyde', 'primary_alcohol', '[A:2]-[A:1]'))

    # ketone → secondary_alcohol (NaBH4, LiAlH4)
    rules.append(_make_reactor('ketone_to_alcohol', 'ketone', 'secondary_alcohol', '[A:2]-[A:1]'))

    # carboxylic_acid → primary_alcohol (LiAlH4)
    rules.append(_make_reactor('acid_to_alcohol', 'carboxylic_acid', 'primary_alcohol', '[A:2]-[A:1]'))

    # ester → alcohol (LiAlH4). OR' group is lost.
    rules.append(_make_reactor('ester_to_alcohol', 'ester', 'primary_alcohol', '[A:2]-[A:1]'))

    # primary_amide → amine (LiAlH4, BH3)
    rules.append(_make_reactor('amide_to_amine', 'primary_amide', 'primary_amine', '[A:1]-[A:2]'))

    # secondary_amide → amine (LiAlH4, BH3)
    rules.append(_make_reactor('amide_to_amine', 'secondary_amide', 'secondary_amine', '[A:1]-[A:2]'))

    # nitrile → primary_amine (LiAlH4, H2/cat)
    rules.append(_make_reactor('nitrile_to_amine', 'nitrile', 'primary_amine', '[A:2]-[A:1]'))

    # nitro → primary_amine (H2/Pd, SnCl2, Fe/HCl)
    rules.append(_make_reactor('nitro_to_amine', 'nitro', 'primary_aniline', '[A:1]'))

    # azide → primary_amine (PPh3/H2O, H2/Pd)
    rules.append(_make_reactor('azide_to_amine', 'azide', 'primary_amine', '[A:1]'))

    # alkene hydrogenation: C=C → C-C (H2/Pd, H2/Pt)
    rules.append(_make_reactor('alkene_hydrogenation', 'alkene', None, '[A:1]-[A:2]'))
    rules.append(_make_reactor('alkene_hydrogenation', 'terminal_alkene', None, '[A:1]-[A:2]'))
    rules.append(_make_reactor('alkene_hydrogenation', 'enamine', None, '[A:1]-[A:2]-[A:3]'))
    rules.append(_make_reactor('alkene_hydrogenation', 'enol_ether', None, '[A:1]-[A:2]-[A:3]'))

    # deoxygenation: alcohol → alkane (Barton-McCombie, Appel+reduction)
    rules.append(_make_reactor('deoxygenation', 'primary_alcohol', None, '[A:2]'))
    rules.append(_make_reactor('deoxygenation', 'secondary_alcohol', None, '[A:2]'))

    # carbonyl → primary amine (reductive amination with NH3)
    rules.append(_make_reactor('carbonyl_to_amine', 'aldehyde', 'primary_amine', '[A:1]-[N:20]'))
    rules.append(_make_reactor('carbonyl_to_amine', 'ketone', 'primary_amine', '[A:1]-[N:20]'))

    # carbonyl → methylene/methyl (Wolff-Kishner, Clemmensen). C=O fully deoxygenated to CH2/CH3.
    rules.append(_make_reactor('ketone_to_methylene', 'ketone', None, '[A:1]'))
    rules.append(_make_reactor('aldehyde_to_methyl', 'aldehyde', None, '[A:1]'))

    # reductive dehalogenation: Ar-X → Ar-H (H2/Pd, n-Bu3SnH, Zn/AcOH)
    rules.append(_make_reactor('debromination', 'aryl_bromide', None, '[A:1]'))
    rules.append(_make_reactor('dechlorination', 'aryl_chloride', None, '[A:1]'))
    rules.append(_make_reactor('deiodination', 'aryl_iodide', None, '[A:1]'))

    return rules


rules = Proxy(_rules)


__all__ = ['rules']
