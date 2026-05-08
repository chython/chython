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


def _make_reactor(rxn_name, fg_name, product_smarts):
    from ._functional import rules
    from ...reactor import Reactor
    from ... import smarts

    q = rules[fg_name]
    product = smarts(product_smarts)
    return rxn_name, fg_name, Reactor((q,), (product,),
                                       delete_atoms=True, one_shot=True, fix_broken_pyrroles=True)


def _rules():
    rules = []

    # 1,3-diketone → isoxazole (+ hydroxylamine, implicit)
    rules.append(_make_reactor('isoxazole', '1_3_diketone', '[A:1]:1:[O:20]:[N:21]:[A:3]:[A:5]:1'))

    # 1,4-diketone → pyridazine (+ hydrazine, implicit)
    rules.append(_make_reactor('pyridazine', '1_4_diketone', '[A:1]:1:[N:20]:[N:21]:[A:3]:[A:6]:[A:5]:1'))

    # Appel: alcohol → alkyl halide (CCl4/CBr4 + PPh3)
    rules.append(_make_reactor('appel', 'primary_alcohol', '[Br:20]-[A:2]'))
    rules.append(_make_reactor('appel', 'secondary_alcohol', '[Br:20]-[A:2]'))

    # Miyaura borylation: ArX → ArBpin (Pd, B2pin2)
    rules.append(_make_reactor('borylation', 'aryl_bromide', '[A:1]-[B:20](-[O:21])-[O:22]'))
    rules.append(_make_reactor('borylation', 'aryl_iodide', '[A:1]-[B:20](-[O:21])-[O:22]'))
    rules.append(_make_reactor('borylation', 'aryl_chloride', '[A:1]-[B:20](-[O:21])-[O:22]'))

    # Nitrile hydrolysis: R-C≡N → R-C(=O)-NH2 (partial, to amide)
    rules.append(_make_reactor('nitrile_hydrolysis', 'nitrile', '[A:2]-[A:1]=[O:20]'))

    # Electrophilic aromatic nitration: Ar-H → Ar-NO2
    rules.append(_make_reactor('nitration', 'arene_ch', '[A:1]-[N;+:20](=[O:21])[O;-:22]'))

    # Electrophilic aromatic halogenation: Ar-H → Ar-X
    rules.append(_make_reactor('bromination', 'arene_ch', '[A:1]-[Br:20]'))
    rules.append(_make_reactor('chlorination', 'arene_ch', '[A:1]-[Cl:20]'))
    rules.append(_make_reactor('iodination', 'arene_ch', '[A:1]-[I:20]'))

    return rules


rules = Proxy(_rules)


__all__ = ['rules']
