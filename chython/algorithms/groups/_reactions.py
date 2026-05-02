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


def _make_reactor(rxn_name, reactants, product_smarts):
    from ._functional import rules
    from ...reactor import Reactor
    from ... import smarts

    queries = []
    fg_names = []
    for fg_name, remap in reactants:
        q = rules[fg_name]
        if remap:
            q = q.copy()
            q.remap(remap)
        queries.append(q)
        fg_names.append(fg_name)

    product = smarts(product_smarts)
    return rxn_name, fg_names, Reactor(tuple(queries), (product,),
                                       delete_atoms=True, one_shot=True, fix_broken_pyrroles=True)


def _rules():
    rules = []

    # --- amidation: RCOX + R'NH2 -> RC(=O)NHR' ---
    # acid atoms: 100=LG, 1=C(=O), 2=O
    # primary_amine/primary_aniline atoms: 1=N, 2=C
    # secondary_amine/secondary_aniline atoms: 1=N, 2=C, 3=C
    for acyl in ('carboxylic_acid', 'acyl_chloride'):
        for amine in ('primary_amine', 'primary_aniline'):
            rules.append(_make_reactor('amidation',
                                       [(acyl, None),
                                        (amine, {1: 3, 2: 4})],
                                       '[A:1](=[A:2])-[A:3]-[A:4]'))
        for amine in ('secondary_amine', 'secondary_aniline'):
            rules.append(_make_reactor('amidation',
                                       [(acyl, None),
                                        (amine, {1: 3, 2: 4, 3: 5})],
                                       '[A:1](=[A:2])-[A:3](-[A:4])-[A:5]'))

    # --- suzuki: ArX + ArB(OH)2 -> Ar-Ar ---
    # aryl_halide atoms: 100=X(LG), 1=C(ar)
    # aryl_boronic_acid atoms: 100=B(LG), 1=C(ar)
    for halide in ('aryl_chloride', 'aryl_bromide', 'aryl_iodide'):
        for boron in ('aryl_boronic_acid', 'aryl_boronic_ester', 'aryl_molander_salt'):
            rules.append(_make_reactor('suzuki',
                                       [(halide, {100: 200}),
                                        (boron, {1: 2})],
                                       '[A:1]-[A:2]'))

    # --- buchwald-hartwig: ArX + R'NH2 -> ArNHR' ---
    # aryl_halide atoms: 100=X(LG), 1=C(ar)
    # primary_amine/primary_aniline atoms: 1=N, 2=C
    # secondary_amine/secondary_aniline atoms: 1=N, 2=C, 3=C
    for halide in ('aryl_chloride', 'aryl_bromide', 'aryl_iodide'):
        for amine in ('primary_amine', 'primary_aniline'):
            rules.append(_make_reactor('buchwald_hartwig',
                                       [(halide, None),
                                        (amine, {1: 3, 2: 4})],
                                       '[A:1]-[A:3]-[A:4]'))
        for amine in ('secondary_amine', 'secondary_aniline'):
            rules.append(_make_reactor('buchwald_hartwig',
                                       [(halide, None),
                                        (amine, {1: 3, 2: 4, 3: 5})],
                                       '[A:1]-[A:3](-[A:4])-[A:5]'))

    # --- ugi 3CR: RCHO + R'NH2 + R''NC -> R'NH-CH(R)-C(=O)NHR'' ---
    # aldehyde atoms: 1=C, 2=O
    # amine atoms: 1=N, 2=C
    # isocyano atoms: 1=N(+), 2=C(-)
    rules.append(_make_reactor('ugi_3cr',
                               [('aldehyde', None),
                                ('primary_amine', {1: 3, 2: 4}),
                                ('isocyano', {1: 5, 2: 6})],
                               '[A:1](-[A:3](-[A:4]))-[A:6](=[O:20])-[A:5]'))

    return rules


rules = Proxy(_rules)


__all__ = ['rules']
