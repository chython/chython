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

    # amidation: RCOX + R'NH2 -> RC(=O)NHR'
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

    # suzuki
    # aryl: ArX + ArB(OH)2 -> Ar-Ar
    # alkenyl: ArX + alkenyl_boronic -> Ar-CH=CH-R
    # alkyl: ArX + alkyl_boronic -> Ar-R
    for halide in ('aryl_chloride', 'aryl_bromide', 'aryl_iodide'):
        for boron in ('aryl_boronic_acid', 'aryl_boronic_ester', 'aryl_molander_salt'):
            rules.append(_make_reactor('suzuki',
                                       [(halide, {100: 200}),
                                        (boron, {1: 2})],
                                       '[A:1]-[A:2]'))

        for boron in ('alkenyl_boronic_acid', 'alkenyl_boronic_ester', 'alkenyl_molander_salt'):
            rules.append(_make_reactor('suzuki',
                                       [(halide, {100: 200}),
                                        (boron, {1: 2, 2: 3})],
                                       '[A:1]-[A:2]=[A:3]'))

        for boron in ('alkyl_boronic_acid', 'alkyl_boronic_ester', 'alkyl_molander_salt'):
            rules.append(_make_reactor('suzuki',
                                       [(halide, {100: 200}),
                                        (boron, {1: 2})],
                                       '[A:1]-[A:2]'))

    # buchwald-hartwig: ArX + R'NH2 -> ArNHR'
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

    # ugi 3CR: RCHO + R'NH2 + R''NC -> R'NH-CH(R)-C(=O)NHR''
    rules.append(_make_reactor('ugi_3cr',
                               [('aldehyde', None),
                                ('primary_amine', {1: 3, 2: 4}),
                                ('isocyano', {1: 5, 2: 6})],
                               '[A:1](-[A:3](-[A:4]))-[A:6](=[O:20])-[A:5]'))

    # mitsunobu ether: ROH + ArOH -> R-O-Ar; ester: ROH + RCOOH -> R-OC(=O)R'
    for alcohol in ('primary_alcohol', 'secondary_alcohol'):
        rules.append(_make_reactor('mitsunobu',
                                   [(alcohol, None),
                                    ('phenol', {1: 3, 2: 4})],
                                   '[A:2]-[A:3]-[A:4]'))
        rules.append(_make_reactor('mitsunobu',
                                   [(alcohol, None),
                                    ('carboxylic_acid', {1: 3, 2: 4})],
                                   '[A:2]-[A:100]-[A:3](=[A:4])'))

    # deoxygenative coupling: ROH + ArX -> R-Ar
    for alcohol in ('primary_alcohol', 'secondary_alcohol'):
        for halide in ('aryl_chloride', 'aryl_bromide', 'aryl_iodide'):
            rules.append(_make_reactor('deoxygenative_coupling',
                                       [(alcohol, None),
                                        (halide, {1: 3})],
                                       '[A:2]-[A:3]'))

    # decarboxylative coupling: R-COOH + ArX -> R-Ar
    for acid in ('alkyl_carboxylic_acid', 'aryl_carboxylic_acid'):
        for halide in ('aryl_chloride', 'aryl_bromide', 'aryl_iodide'):
            rules.append(_make_reactor('decarboxylative_coupling',
                                       [(acid, None),
                                        (halide, {100: 200, 1: 4})],
                                       '[A:3]-[A:4]'))

    # XEC (cross-electrophile coupling): ArX + R'X -> Ar-R'
    for halide in ('aryl_chloride', 'aryl_bromide', 'aryl_iodide'):
        for alkyl_halide in ('alkyl_chloride', 'alkyl_bromide', 'alkyl_iodide'):
            rules.append(_make_reactor('xec',
                                       [(halide, None),
                                        (alkyl_halide, {100: 200, 1: 2})],
                                       '[A:1]-[A:2]'))

    # reductive amination: RCHO/R2CO + amine -> amine
    for carbonyl in ('aldehyde', 'ketone'):
        for amine in ('primary_amine', 'primary_aniline'):
            rules.append(_make_reactor('reductive_amination',
                                       [(carbonyl, None),
                                        (amine, {1: 3, 2: 4})],
                                       '[A:1]-[A:3]-[A:4]'))
        for amine in ('secondary_amine', 'secondary_aniline'):
            rules.append(_make_reactor('reductive_amination',
                                       [(carbonyl, None),
                                        (amine, {1: 3, 2: 4, 3: 5})],
                                       '[A:1]-[A:3](-[A:4])-[A:5]'))

    # ullmann phenol: ArX + ArOH -> Ar-O-Ar
    for halide in ('aryl_chloride', 'aryl_bromide', 'aryl_iodide'):
        rules.append(_make_reactor('ullmann_phenol',
                                   [(halide, None),
                                    ('phenol', {1: 3, 2: 4})],
                                   '[A:1]-[A:3]-[A:4]'))

    # ullmann pyrrole: ArX + pyrrole-NH -> Ar-N(pyrrole)
    for halide in ('aryl_chloride', 'aryl_bromide', 'aryl_iodide'):
        rules.append(_make_reactor('ullmann_pyrrole',
                                   [(halide, None),
                                    ('pyrrole', {1: 3})],
                                   '[A:1]-[A:3]'))
        rules.append(_make_reactor('ullmann_pyrrole',
                                   [(halide, None),
                                    ('pyrazole', {1: 3, 2: 4})],
                                   '[A:1]-[A:4]:[A:3]'))
        rules.append(_make_reactor('ullmann_pyrrole',
                                   [(halide, None),
                                    ('imidazole', {1: 3, 2: 4, 3: 5})],
                                   '[A:1]-[A:5]:[A:4]:[A:3]'))

    # chan-lam: ArB(OH)2 + amine/phenol -> Ar-N/Ar-O
    for boron in ('aryl_boronic_acid', 'aryl_boronic_ester'):
        for amine in ('primary_amine', 'primary_aniline'):
            rules.append(_make_reactor('chan_lam',
                                       [(boron, None),
                                        (amine, {1: 3, 2: 4})],
                                       '[A:1]-[A:3]-[A:4]'))
        for amine in ('secondary_amine', 'secondary_aniline'):
            rules.append(_make_reactor('chan_lam',
                                       [(boron, None),
                                        (amine, {1: 3, 2: 4, 3: 5})],
                                       '[A:1]-[A:3](-[A:4])-[A:5]'))
        rules.append(_make_reactor('chan_lam',
                                   [(boron, None),
                                    ('phenol', {1: 3, 2: 4})],
                                   '[A:1]-[A:3]-[A:4]'))

    # sulfonylation: RSO2X + alcohol/phenol -> sulfonate ester
    for sulfonyl in ('sulfonyl_chloride', 'sulfonyl_fluoride'):
        for alcohol in ('primary_alcohol', 'secondary_alcohol'):
            rules.append(_make_reactor('sulfonylation',
                                       [(sulfonyl, None),
                                        (alcohol, {1: 4, 2: 5})],
                                       '[A:1](=[A:2])(=[A:3])-[A:4]-[A:5]'))
        rules.append(_make_reactor('sulfonylation',
                                   [(sulfonyl, None),
                                    ('phenol', {1: 4, 2: 5})],
                                   '[A:1](=[A:2])(=[A:3])-[A:4]-[A:5]'))

    # sulfonamide: RSO2X + amine -> sulfonamide
    for sulfonyl in ('sulfonyl_chloride', 'sulfonyl_fluoride'):
        for amine in ('primary_amine', 'primary_aniline'):
            rules.append(_make_reactor('sulfonamide_formation',
                                       [(sulfonyl, None),
                                        (amine, {1: 4, 2: 5})],
                                       '[A:1](=[A:2])(=[A:3])-[A:4]-[A:5]'))
        for amine in ('secondary_amine', 'secondary_aniline'):
            rules.append(_make_reactor('sulfonamide_formation',
                                       [(sulfonyl, None),
                                        (amine, {1: 4, 2: 5, 3: 6})],
                                       '[A:1](=[A:2])(=[A:3])-[A:4](-[A:5])-[A:6]'))

    # aminolysis: ester + amine -> amide
    for amine in ('primary_amine', 'primary_aniline'):
        rules.append(_make_reactor('aminolysis',
                                   [('ester', None),
                                    (amine, {1: 3, 2: 4})],
                                   '[A:1](=[A:2])-[A:3]-[A:4]'))
    for amine in ('secondary_amine', 'secondary_aniline'):
        rules.append(_make_reactor('aminolysis',
                                   [('ester', None),
                                    (amine, {1: 3, 2: 4, 3: 5})],
                                   '[A:1](=[A:2])-[A:3](-[A:4])-[A:5]'))

    # grignard: RX + RCHO/R2CO -> alcohol
    for alkyl_halide in ('alkyl_chloride', 'alkyl_bromide', 'alkyl_iodide'):
        for carbonyl in ('aldehyde', 'ketone'):
            rules.append(_make_reactor('grignard',
                                       [(alkyl_halide, None),
                                        (carbonyl, {1: 2, 2: 3})],
                                       '[A:1]-[A:2]-[A:3]'))

    # sonogashira: ArX + terminal alkyne -> Ar-C≡C-R
    for halide in ('aryl_chloride', 'aryl_bromide', 'aryl_iodide'):
        rules.append(_make_reactor('sonogashira',
                                   [(halide, None),
                                    ('terminal_alkyne', {1: 2, 2: 3})],
                                   '[A:1]-[A:2]#[A:3]'))

    # urea synthesis: R-N=C=O + amine -> R-NH-C(=O)-NH-R'
    for amine in ('primary_amine', 'primary_aniline'):
        rules.append(_make_reactor('urea_synthesis',
                                   [('isocyanate', None),
                                    (amine, {1: 4, 2: 5})],
                                   '[A:1]-[A:2](=[A:3])-[A:4]-[A:5]'))
    for amine in ('secondary_amine', 'secondary_aniline'):
        rules.append(_make_reactor('urea_synthesis',
                                   [('isocyanate', None),
                                    (amine, {1: 4, 2: 5, 3: 6})],
                                   '[A:1]-[A:2](=[A:3])-[A:4](-[A:5])-[A:6]'))

    # CuAAC: azide + terminal alkyne -> 1,2,3-triazole
    rules.append(_make_reactor('cuaac',
                               [('azide', None),
                                ('terminal_alkyne', {1: 4, 2: 5})],
                               '[A:1]:1:[A:2]:[A:3]:[A:5]:[A:4]:1'))

    # SNAr: ArF + amine -> Ar-NR2 (nucleophilic aromatic substitution)
    for amine in ('primary_amine', 'primary_aniline'):
        rules.append(_make_reactor('snar',
                                   [('aryl_fluoride', None),
                                    (amine, {1: 3, 2: 4})],
                                   '[A:1]-[A:3]-[A:4]'))
    for amine in ('secondary_amine', 'secondary_aniline'):
        rules.append(_make_reactor('snar',
                                   [('aryl_fluoride', None),
                                    (amine, {1: 3, 2: 4, 3: 5})],
                                   '[A:1]-[A:3](-[A:4])-[A:5]'))

    # SNAr: ArF + phenol -> Ar-O-Ar
    rules.append(_make_reactor('snar',
                               [('aryl_fluoride', None),
                                ('phenol', {1: 3, 2: 4})],
                               '[A:1]-[A:3]-[A:4]'))

    # SNAr: ArF + thiol -> Ar-S-R
    rules.append(_make_reactor('snar',
                               [('aryl_fluoride', None),
                                ('thiol', {1: 3, 2: 4})],
                               '[A:1]-[A:3]-[A:4]'))

    # Williamson ether:
    # alkyl_halide + phenol -> R-O-Ar
    # alkyl_halide + primary_alcohol -> R-O-R'
    for halide in ('alkyl_chloride', 'alkyl_bromide', 'alkyl_iodide'):
        rules.append(_make_reactor('williamson',
                                   [(halide, None),
                                    ('phenol', {1: 3, 2: 4})],
                                   '[A:1]-[A:3]-[A:4]'))

        for alcohol in ('primary_alcohol', 'secondary_alcohol'):
            rules.append(_make_reactor('williamson',
                                       [(halide, None),
                                        (alcohol, {1: 3, 2: 4})],
                                       '[A:1]-[A:3]-[A:4]'))

    # acylation (ester): acyl_chloride + alcohol -> ester
    for alcohol in ('primary_alcohol', 'secondary_alcohol'):
        rules.append(_make_reactor('acylation',
                                   [('acyl_chloride', None),
                                    (alcohol, {1: 3, 2: 4})],
                                   '[A:1](=[A:2])-[A:3]-[A:4]'))
    rules.append(_make_reactor('acylation',
                               [('acyl_chloride', None),
                                ('phenol', {1: 3, 2: 4})],
                               '[A:1](=[A:2])-[A:3]-[A:4]'))

    # thioether: thiol + alkyl_halide -> R-S-R'
    for halide in ('alkyl_chloride', 'alkyl_bromide', 'alkyl_iodide'):
        rules.append(_make_reactor('thioether',
                                   [('thiol', None),
                                    (halide, {100: 200, 1: 3})],
                                   '[A:1](-[A:2])-[A:3]'))

    # carbamate: isocyanate + alcohol -> R-NH-C(=O)-OR'
    for alcohol in ('primary_alcohol', 'secondary_alcohol'):
        rules.append(_make_reactor('carbamate',
                                   [('isocyanate', None),
                                    (alcohol, {1: 4, 2: 5})],
                                   '[A:1]-[A:2](=[A:3])-[A:4]-[A:5]'))
    rules.append(_make_reactor('carbamate',
                               [('isocyanate', None),
                                ('phenol', {1: 4, 2: 5})],
                               '[A:1]-[A:2](=[A:3])-[A:4]-[A:5]'))

    # tetrazole: azide + nitrile -> 2,5-disubstituted tetrazole (ring formation)
    rules.append(_make_reactor('tetrazole',
                               [('azide', None),
                                ('nitrile', {1: 4, 2: 5})],
                               '[A:1]:1:[A:2]:[A:3]:[A:4]:[A:5]:1'))

    return rules


rules = Proxy(_rules)


__all__ = ['rules']
