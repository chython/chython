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

    return rules


rules = Proxy(_rules)


__all__ = ['rules']
