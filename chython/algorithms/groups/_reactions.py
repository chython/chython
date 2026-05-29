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

    # carbamoylation: ROC(=O)X + amine -> ROC(=O)NR'R''  (carbamate from chloroformate)
    for formate in ('chloroformate', 'fluoroformate'):
        for amine in ('primary_amine', 'primary_aniline'):
            rules.append(_make_reactor('carbamoylation',
                                       [(formate, None),
                                        (amine, {1: 4, 2: 5})],
                                       '[A:3]-[A:1](=[A:2])-[A:4]-[A:5]'))
        for amine in ('secondary_amine', 'secondary_aniline'):
            rules.append(_make_reactor('carbamoylation',
                                       [(formate, None),
                                        (amine, {1: 4, 2: 5, 3: 6})],
                                       '[A:3]-[A:1](=[A:2])-[A:4](-[A:5])-[A:6]'))
        rules.append(_make_reactor('carbamoylation',
                                   [(formate, None),
                                    ('pyrrole', {1: 4})],
                                   '[A:3]-[A:1](=[A:2])-[A:4]'))
        rules.append(_make_reactor('carbamoylation',
                                   [(formate, None),
                                    ('imidazole', {1: 4, 2: 5, 3: 6})],
                                   '[A:3]-[A:1](=[A:2])-[A:6]:[A:5]:[A:4]'))

    # urea from carbamoyl chloride: R2NC(=O)X + amine -> R2NC(=O)NR'R''
    for carbamoyl in ('carbamoyl_chloride', 'carbamoyl_fluoride'):
        for amine in ('primary_amine', 'primary_aniline'):
            rules.append(_make_reactor('urea_from_carbamoyl',
                                       [(carbamoyl, None),
                                        (amine, {1: 4, 2: 5})],
                                       '[A:3]-[A:1](=[A:2])-[A:4]-[A:5]'))
        for amine in ('secondary_amine', 'secondary_aniline'):
            rules.append(_make_reactor('urea_from_carbamoyl',
                                       [(carbamoyl, None),
                                        (amine, {1: 4, 2: 5, 3: 6})],
                                       '[A:3]-[A:1](=[A:2])-[A:4](-[A:5])-[A:6]'))

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
        for amine in ('primary_amine', 'primary_aniline', 'primary_amidine_amine'):
            rules.append(_make_reactor('buchwald_hartwig',
                                       [(halide, None),
                                        (amine, {1: 3, 2: 4})],
                                       '[A:1]-[A:3]-[A:4]'))
        for amine in ('secondary_amine', 'secondary_aniline'):
            rules.append(_make_reactor('buchwald_hartwig',
                                       [(halide, None),
                                        (amine, {1: 3, 2: 4, 3: 5})],
                                       '[A:1]-[A:3](-[A:4])-[A:5]'))

    # buchwald-hartwig: lactam halides + amines
    for _x in ('fluoride', 'chloride', 'bromide', 'iodide'):
        for _n in ('1', '2', '3', '4', '5', '6'):
            lactam = f'lactam_{_n}_{_x}'
            for amine in ('primary_amine', 'primary_aniline'):
                rules.append(_make_reactor('buchwald_hartwig',
                                           [(lactam, None),
                                            (amine, {1: 3, 2: 4})],
                                           '[A:1]-[A:3]-[A:4]'))
            for amine in ('secondary_amine', 'secondary_aniline'):
                rules.append(_make_reactor('buchwald_hartwig',
                                           [(lactam, None),
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

    # ullmann alcohol: ArX + ROH -> Ar-O-R (SNAr/Cu-mediated with aliphatic alcohols)
    for halide in ('aryl_chloride', 'aryl_bromide', 'aryl_iodide', 'aryl_fluoride'):
        for alcohol in ('primary_alcohol', 'secondary_alcohol'):
            rules.append(_make_reactor('snar',
                                       [(halide, None),
                                        (alcohol, {1: 3, 2: 4})],
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

    # ullmann pyrrole: lactam halides + N-heterocycles
    for _x in ('fluoride', 'chloride', 'bromide', 'iodide'):
        for _n in ('1', '2', '3', '4', '5', '6'):
            lactam = f'lactam_{_n}_{_x}'
            rules.append(_make_reactor('ullmann_pyrrole',
                                       [(lactam, None),
                                        ('pyrrole', {1: 3})],
                                       '[A:1]-[A:3]'))
            rules.append(_make_reactor('ullmann_pyrrole',
                                       [(lactam, None),
                                        ('pyrazole', {1: 3, 2: 4})],
                                       '[A:1]-[A:4]:[A:3]'))
            rules.append(_make_reactor('ullmann_pyrrole',
                                       [(lactam, None),
                                        ('imidazole', {1: 3, 2: 4, 3: 5})],
                                       '[A:1]-[A:5]:[A:4]:[A:3]'))

    # ullmann pyridol: ArX + hydroxypyridine -> Ar-N-C=O (CN coupling on pyridol tautomer)
    for halide in ('aryl_fluoride', 'aryl_chloride', 'aryl_bromide', 'aryl_iodide'):
        rules.append(_make_reactor('ullmann_pyrrole',
                                   [(halide, None),
                                    ('pyridol', {1: 3, 2: 4, 3: 5})],
                                   '[A:1]-[A:3]-[A:4]=[A:5]'))

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

    # sulfonamide from amide N-H: RSO2X + R'C(=O)NHR" -> R'C(=O)N(SO2R)R"
    for sulfonyl in ('sulfonyl_chloride', 'sulfonyl_fluoride'):
        for amide in ('primary_amide', 'secondary_amide'):
            rules.append(_make_reactor('sulfonamide_formation',
                                       [(sulfonyl, None),
                                        (amide, {1: 4, 2: 5, 3: 6})],
                                       '[A:1](=[A:2])(=[A:3])-[A:4]-[A:5]=[A:6]'))

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

    # grignard: RMgX + RCHO/R2CO -> alcohol
    for grignard in ('alkyl_grignard', 'aryl_grignard'):
        for carbonyl in ('aldehyde', 'ketone'):
            rules.append(_make_reactor('grignard',
                                       [(grignard, None),
                                        (carbonyl, {1: 3, 2: 4})],
                                       '[A:1]-[A:3]-[A:4]'))
    for carbonyl in ('aldehyde', 'ketone'):
        rules.append(_make_reactor('grignard',
                                   [('alkenyl_grignard', None),
                                    (carbonyl, {1: 3, 2: 4})],
                                   '[A:1]=[A:2]-[A:3]-[A:4]'))
    # grignard from halide surrogates
    for alkyl_halide in ('alkyl_chloride', 'alkyl_bromide', 'alkyl_iodide'):
        for carbonyl in ('aldehyde', 'ketone'):
            rules.append(_make_reactor('grignard',
                                       [(alkyl_halide, None),
                                        (carbonyl, {1: 2, 2: 3})],
                                       '[A:1]-[A:2]-[A:3]'))
    for aryl_halide in ('aryl_chloride', 'aryl_bromide', 'aryl_iodide'):
        for carbonyl in ('aldehyde', 'ketone'):
            rules.append(_make_reactor('grignard',
                                       [(aryl_halide, None),
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

        for alcohol in ('primary_alcohol', 'secondary_alcohol', 'tertiary_alcohol'):
            rules.append(_make_reactor('williamson',
                                       [(halide, None),
                                        (alcohol, {1: 3, 2: 4})],
                                       '[A:1]-[A:3]-[A:4]'))

    # Williamson ether with pseudohalides:
    for pseudohalide in ('alkyl_triflate', 'alkyl_mesylate', 'alkyl_tosylate'):
        rules.append(_make_reactor('williamson',
                                   [(pseudohalide, None),
                                    ('phenol', {1: 3, 2: 4})],
                                   '[A:1]-[A:3]-[A:4]'))
        for alcohol in ('primary_alcohol', 'secondary_alcohol', 'tertiary_alcohol'):
            rules.append(_make_reactor('williamson',
                                       [(pseudohalide, None),
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

    # wittig: phosphonium_ylide + aldehyde/ketone -> alkene
    for carbonyl in ('aldehyde', 'ketone'):
        rules.append(_make_reactor('wittig',
                                   [('phosphonium_ylide', None),
                                    (carbonyl, {1: 2, 2: 3})],
                                   '[A:1]=[A:2]'))

    # hwe: phosphonate + aldehyde/ketone -> alkene (E-selective)
    for carbonyl in ('aldehyde', 'ketone'):
        rules.append(_make_reactor('hwe',
                                   [('phosphonate', None),
                                    (carbonyl, {1: 2, 2: 3})],
                                   '[A:1]=[A:2]'))

    # weinreb: weinreb_amide + grignard -> ketone
    for grignard in ('alkyl_grignard', 'aryl_grignard'):
        rules.append(_make_reactor('weinreb',
                                   [('weinreb_amide', None),
                                    (grignard, {100: 200, 101: 201, 1: 3})],
                                   '[A:2]=[A:1]-[A:3]'))
    # weinreb from ester surrogate: ester + grignard -> ketone
    for grignard in ('alkyl_grignard', 'aryl_grignard'):
        rules.append(_make_reactor('weinreb',
                                   [('ester', None),
                                    (grignard, {100: 200, 101: 201, 1: 3})],
                                   '[A:2]=[A:1]-[A:3]'))
    # weinreb from acyl_chloride surrogate: acyl_chloride + grignard -> ketone
    for grignard in ('alkyl_grignard', 'aryl_grignard'):
        rules.append(_make_reactor('weinreb',
                                   [('acyl_chloride', None),
                                    (grignard, {100: 200, 101: 201, 1: 3})],
                                   '[A:2]=[A:1]-[A:3]'))

    # friedel-crafts acylation: acyl_chloride + arene -> Ar-COR
    rules.append(_make_reactor('friedel_crafts',
                               [('acyl_chloride', None),
                                ('arene_ch', {1: 3})],
                               '[A:1](=[A:2])-[A:3]'))

    # heck: ArX + alkene -> Ar-CH=CH-R
    for halide in ('aryl_chloride', 'aryl_bromide', 'aryl_iodide'):
        for alkene in ('terminal_alkene', 'alkene'):
            rules.append(_make_reactor('heck',
                                       [(halide, None),
                                        (alkene, {1: 2, 2: 3})],
                                       '[A:1]-[A:2]=[A:3]'))

    # cross metathesis: terminal_alkene + terminal_alkene -> internal alkene
    rules.append(_make_reactor('cross_metathesis',
                               [('terminal_alkene', None),
                                ('terminal_alkene', {1: 3, 2: 4})],
                               '[A:2]=[A:4]'))

    # kumada: ArX + RMgX -> Ar-R
    for halide in ('aryl_chloride', 'aryl_bromide', 'aryl_iodide'):
        for grignard in ('alkyl_grignard', 'aryl_grignard'):
            rules.append(_make_reactor('kumada',
                                       [(halide, None),
                                        (grignard, {100: 200, 101: 201, 1: 2})],
                                       '[A:1]-[A:2]'))
        rules.append(_make_reactor('kumada',
                                   [(halide, None),
                                    ('alkenyl_grignard', {100: 200, 101: 201, 1: 2, 2: 3})],
                                   '[A:1]-[A:2]=[A:3]'))

    # negishi: ArX + RZnX -> Ar-R
    for halide in ('aryl_chloride', 'aryl_bromide', 'aryl_iodide'):
        for zinc in ('alkyl_zinc', 'aryl_zinc'):
            rules.append(_make_reactor('negishi',
                                       [(halide, None),
                                        (zinc, {100: 200, 101: 201, 1: 2})],
                                       '[A:1]-[A:2]'))
        rules.append(_make_reactor('negishi',
                                   [(halide, None),
                                    ('alkenyl_zinc', {100: 200, 101: 201, 1: 2, 2: 3})],
                                   '[A:1]-[A:2]=[A:3]'))

    # stille: ArX + R-SnR3 -> Ar-R
    for halide in ('aryl_chloride', 'aryl_bromide', 'aryl_iodide'):
        rules.append(_make_reactor('stille',
                                   [(halide, None),
                                    ('aryl_stannane', {100: 200, 1: 2})],
                                   '[A:1]-[A:2]'))
        rules.append(_make_reactor('stille',
                                   [(halide, None),
                                    ('alkenyl_stannane', {100: 200, 1: 2, 2: 3})],
                                   '[A:1]-[A:2]=[A:3]'))
        rules.append(_make_reactor('stille',
                                   [(halide, None),
                                    ('alkyl_stannane', {100: 200, 1: 2})],
                                   '[A:1]-[A:2]'))

    # hiyama: ArX + R-SiR3 -> Ar-R
    for halide in ('aryl_chloride', 'aryl_bromide', 'aryl_iodide'):
        rules.append(_make_reactor('hiyama',
                                   [(halide, None),
                                    ('aryl_silane', {100: 200, 1: 2})],
                                   '[A:1]-[A:2]'))
        rules.append(_make_reactor('hiyama',
                                   [(halide, None),
                                    ('alkenyl_silane', {100: 200, 1: 2, 2: 3})],
                                   '[A:1]-[A:2]=[A:3]'))

    # hantzsch thiazole: alpha_haloketone + thioamide -> thiazole
    rules.append(_make_reactor('hantzsch_thiazole',
                               [('alpha_haloketone', None),
                                ('thioamide', {1: 4, 2: 5, 3: 6})],
                               '[A:5]:1:[A:4]:[A:6]:[A:1]:[A:3]:1'))

    # knorr pyrazole: 1,3-diketone + alkyl_hydrazine -> pyrazole
    rules.append(_make_reactor('knorr_pyrazole',
                               [('1_3_diketone', None),
                                ('alkyl_hydrazine', {1: 6, 2: 7, 3: 8})],
                               '[A:6](:1:[A:7]:[A:1]:[A:5]:[A:3]:1)-[A:8]'))

    # paal_knorr pyrrole: 1,4-diketone + primary_amine -> pyrrole
    for amine in ('primary_amine', 'primary_aniline'):
        rules.append(_make_reactor('paal_knorr',
                                   [('1_4_diketone', None),
                                    (amine, {1: 7, 2: 8})],
                                   '[A:7](:1:[A:1]:[A:5]:[A:6]:[A:3]:1)-[A:8]'))

    # fischer indole: aryl_hydrazine + alpha_ketone -> indole
    rules.append(_make_reactor('fischer_indole',
                               [('aryl_hydrazine', None),
                                ('alpha_ketone', {1: 5, 2: 6, 3: 7})],
                               '[A:1]:1:[A:3]:[A:4]:[A:7]:[A:5]:1'))

    # benzimidazole: o_diaminoarene + aldehyde/acid -> benzimidazole
    for carbonyl in ('aldehyde', 'carboxylic_acid'):
        rules.append(_make_reactor('benzimidazole',
                                   [('o_diaminoarene', None),
                                    (carbonyl, {1: 5, 2: 6})],
                                   '[A:1]:1:[A:3]:[A:4]:[A:2]:[A:5]:1'))

    # benzoxazole: o_aminophenol + aldehyde/acid -> benzoxazole
    for carbonyl in ('aldehyde', 'carboxylic_acid'):
        rules.append(_make_reactor('benzoxazole',
                                   [('o_aminophenol', None),
                                    (carbonyl, {1: 5, 2: 6})],
                                   '[A:1]:1:[A:3]:[A:4]:[A:2]:[A:5]:1'))

    # benzothiazole: o_aminothiophenol + aldehyde/acid -> benzothiazole
    for carbonyl in ('aldehyde', 'carboxylic_acid'):
        rules.append(_make_reactor('benzothiazole',
                                   [('o_aminothiophenol', None),
                                    (carbonyl, {1: 5, 2: 6})],
                                   '[A:1]:1:[A:3]:[A:4]:[A:2]:[A:5]:1'))

    # quinoxaline: o_diaminoarene + 1,2-diketone -> quinoxaline
    rules.append(_make_reactor('quinoxaline',
                               [('o_diaminoarene', None),
                                ('1_2_diketone', {1: 5, 2: 6, 3: 7, 4: 8})],
                               '[A:1]:1:[A:3]:[A:4]:[A:2]:[A:7]:[A:5]:1'))

    # friedlander quinoline: o_aminobenzaldehyde + alpha_ketone -> quinoline
    rules.append(_make_reactor('friedlander',
                               [('o_aminobenzaldehyde', None),
                                ('alpha_ketone', {1: 7, 2: 8, 3: 9})],
                               '[A:1]:1:[A:3]:[A:4]:[A:5]:[A:9]:[A:7]:1'))

    # pictet_spengler: beta_arylethylamine + aldehyde -> THIQ
    rules.append(_make_reactor('pictet_spengler',
                               [('beta_arylethylamine', None),
                                ('aldehyde', {1: 6, 2: 7})],
                               '[A:1]1-[A:2]-[A:3]-[A:4]:[A:5]-[A:6]-1'))

    # 1,2,4-oxadiazole: amidoxime + acyl source -> 1,2,4-oxadiazole
    for acyl in ('acyl_chloride', 'acyl_fluoride', 'carboxylic_acid'):
        rules.append(_make_reactor('oxadiazole_124',
                                   [('amidoxime', None),
                                    (acyl, {1: 5, 2: 6, 100: 200})],
                                   '[A:4]:1:[A:1]:[A:3]:[A:5]:[A:2]:1'))

    # pyrimidine: 1,3-diketone + amidine -> pyrimidine
    rules.append(_make_reactor('pyrimidine',
                               [('1_3_diketone', None),
                                ('amidine', {1: 6, 2: 7, 3: 8})],
                               '[A:7]:1:[A:6]:[A:8]:[A:1]:[A:5]:[A:3]:1'))

    # van_leusen oxazole: alpha_isocyano + aldehyde -> oxazole
    rules.append(_make_reactor('van_leusen_oxazole',
                               [('tosyl_isocyanide', None),
                                ('aldehyde', {1: 4, 2: 5})],
                               '[A:1]:1:[A:2]:[A:3]:[A:5]:[A:4]:1'))

    # imidazo[1,2-a]pyridine: aminopyridine + alpha_haloketone
    rules.append(_make_reactor('imidazopyridine',
                               [('aminopyridine', None),
                                ('alpha_haloketone', {1: 4, 2: 5, 3: 6, 100: 200})],
                               '[A:3]:1:[A:6]:[A:4]:[A:1]:[A:2]:1'))
    # imidazo[1,2-a]pyridine: aminopyridine + alpha_haloester (ester O masked, carbonyl O deleted)
    rules.append(_make_reactor('imidazopyridine',
                               [('aminopyridine', None),
                                ('alpha_haloester', {1: 4, 2: 5, 3: 6, 100: 200})],
                               '[A:3]:1:[A:6]:[A:4]:[A:1]:[A:2]:1'))

    # biginelli: aldehyde + beta_ketoester + urea/thiourea -> DHPM (3-component)
    for urea_fg in ('urea', 'thiourea'):
        rules.append(_make_reactor('biginelli',
                                   [('aldehyde', None),
                                    ('beta_ketoester', {1: 3, 2: 4, 3: 5, 4: 6, 5: 7, 100: 200}),
                                    (urea_fg, {1: 8, 2: 9, 3: 10, 4: 11})],
                                   '[A:8]1-[A:9](=[A:10])-[A:11]-[A:3]=[A:7]-[A:1]-1'))

    # niementowski quinazoline: anthranilic_acid + primary_amide -> 4-oxoquinazoline
    rules.append(_make_reactor('niementowski',
                               [('anthranilic_acid', None),
                                ('primary_amide', {1: 7, 2: 8, 3: 9})],
                               '[A:1]:1:[A:3]:[A:4]:[A:5](-[A:6]):[A:7]:[A:8]:1'))

    # ugi 4CR: RCHO + R'NH2 + R"COOH + R"'NC -> R"C(=O)N(R')CH(R)C(=O)NHR"'
    rules.append(_make_reactor('ugi_4cr',
                               [('aldehyde', None),
                                ('primary_amine', {1: 3, 2: 4}),
                                ('carboxylic_acid', {100: 200, 1: 7, 2: 8}),
                                ('isocyano', {1: 5, 2: 6})],
                               '[A:7](=[A:8])-[A:3](-[A:4])-[A:1]-[A:6](=[O:20])-[A:5]'))

    # liebeskind-srogl: R-C(=O)-SR' + ArB(OH)2 -> R-C(=O)-Ar (ketone synthesis)
    for boron in ('aryl_boronic_acid', 'aryl_boronic_ester'):
        rules.append(_make_reactor('liebeskind_srogl',
                                   [('thioester', None),
                                    (boron, {100: 200, 101: 201, 102: 202, 1: 3})],
                                   '[A:1](=[A:2])-[A:3]'))

    # knoevenagel: active_methylene + aldehyde -> alkene (condensation)
    rules.append(_make_reactor('knoevenagel',
                               [('active_methylene', None),
                                ('aldehyde', {1: 4, 2: 5})],
                               '[A:4]=[A:1](-[A:2])-[A:3]'))

    # aldol: alpha_ketone + aldehyde -> beta-hydroxy carbonyl
    rules.append(_make_reactor('aldol',
                               [('alpha_ketone', None),
                                ('aldehyde', {1: 4, 2: 5})],
                               '[A:1](=[A:2])-[A:3]-[A:4]-[A:5]'))

    # larock indole: o_haloaniline + alkyne -> indole
    rules.append(_make_reactor('larock_indole',
                               [('o_haloaniline', None),
                                ('alkyne', {1: 4, 2: 5})],
                               '[A:1]:1:[A:2]:[A:3]:[A:5]:[A:4]:1'))

    # doebner-miller: aniline_ortho_ch + enal -> quinoline
    rules.append(_make_reactor('doebner_miller',
                               [('aniline_ortho_ch', None),
                                ('enal', {1: 4, 2: 5, 3: 6, 4: 7})],
                               '[A:1]:1:[A:2]:[A:3]:[A:7]:[A:6]:[A:4]:1'))

    # hantzsch pyridine: aldehyde + 2x beta_ketoester -> pyridine (NH3 implicit)
    rules.append(_make_reactor('hantzsch_pyridine',
                               [('aldehyde', {1: 11, 2: 12}),
                                ('beta_ketoester', None),
                                ('beta_ketoester', {1: 6, 2: 7, 3: 8, 4: 9, 5: 10, 100: 200})],
                               '[N:20]:1:[A:1]:[A:5](-[A:3](=[A:4])-[A:100]):[A:11]:[A:10](-[A:8](=[A:9])-[A:200]):[A:6]:1'))

    # nitrile_grignard: R-C≡N + RMgX → ketone (via imine hydrolysis)
    for grignard in ('alkyl_grignard', 'aryl_grignard'):
        rules.append(_make_reactor('nitrile_grignard',
                                   [('nitrile', None),
                                    (grignard, {100: 200, 101: 201, 1: 3})],
                                   '[A:1](=[O:20])-[A:3]'))

    # N-alkylation: alkyl_halide/triflate/tosylate + pyrrole/pyrazole/imidazole/pyridol
    for halide in ('alkyl_chloride', 'alkyl_bromide', 'alkyl_iodide', 'alkyl_triflate', 'alkyl_tosylate'):
        rules.append(_make_reactor('n_alkylation',
                                   [(halide, None),
                                    ('pyrrole', {1: 3})],
                                   '[A:1]-[A:3]'))
        rules.append(_make_reactor('n_alkylation',
                                   [(halide, None),
                                    ('pyrazole', {1: 3, 2: 4})],
                                   '[A:1]-[A:4]:[A:3]'))
        rules.append(_make_reactor('n_alkylation',
                                   [(halide, None),
                                    ('imidazole', {1: 3, 2: 4, 3: 5})],
                                   '[A:1]-[A:5]:[A:4]:[A:3]'))
        rules.append(_make_reactor('n_alkylation',
                                   [(halide, None),
                                    ('pyridol', {1: 3, 2: 4, 3: 5})],
                                   '[A:1]-[A:3]-[A:4]=[A:5]'))

    # N-alkylation: alkyl_halide/triflate + primary/secondary amine
    for halide in ('alkyl_chloride', 'alkyl_bromide', 'alkyl_iodide', 'alkyl_triflate',
                   'boronate_alkyl_chloride', 'boronate_alkyl_bromide', 'boronate_alkyl_iodide'):
        for amine in ('primary_amine', 'primary_aniline'):
            rules.append(_make_reactor('n_alkylation',
                                       [(halide, None),
                                        (amine, {1: 3, 2: 4})],
                                       '[A:1]-[A:3]-[A:4]'))
        for amine in ('secondary_amine', 'secondary_aniline'):
            rules.append(_make_reactor('n_alkylation',
                                       [(halide, None),
                                        (amine, {1: 3, 2: 4, 3: 5})],
                                       '[A:1]-[A:3](-[A:4])-[A:5]'))

    # urea from 2 amines (CDI/phosgene implicit)
    for amine1 in ('primary_amine', 'primary_aniline'):
        for amine2 in ('primary_amine', 'primary_aniline'):
            rules.append(_make_reactor('urea_from_amines',
                                       [(amine1, None),
                                        (amine2, {1: 3, 2: 4})],
                                       '[A:1](-[A:2])-[C:20](=[O:21])-[A:3]-[A:4]'))

    # oxazoline: amino_alcohol + carboxylic_acid → 2-oxazoline ring
    rules.append(_make_reactor('oxazoline',
                               [('amino_alcohol', None),
                                ('carboxylic_acid', {100: 200, 1: 5, 2: 6})],
                               '[A:5]1=[A:1]-[A:2]-[A:3]-[A:4]-1'))

    # oxime O-alkylation: oxime + alkyl_halide → oxime ether
    for halide in ('alkyl_chloride', 'alkyl_bromide', 'alkyl_iodide'):
        rules.append(_make_reactor('oxime_alkylation',
                                   [('oxime', None),
                                    (halide, {100: 200, 1: 4})],
                                   '[A:1](-[A:4])-[A:2]=[A:3]'))

    # oxime ether formation: O-alkylhydroxylamine + aldehyde → R-O-N=CH-R'
    rules.append(_make_reactor('oxime_ether',
                               [('O_alkylhydroxylamine', None),
                                ('aldehyde', {1: 4, 2: 5})],
                               '[A:3]-[A:2]-[A:1]=[A:4]'))

    # aldol with ketone (extension): alpha_ketone + ketone
    rules.append(_make_reactor('aldol',
                               [('alpha_ketone', None),
                                ('ketone', {1: 4, 2: 5})],
                               '[A:1](=[A:2])-[A:3]-[A:4]-[A:5]'))

    # hydrazone formation: aldehyde/ketone + aryl_hydrazine -> C=N-NH-Ar
    for carbonyl in ('aldehyde', 'ketone'):
        rules.append(_make_reactor('hydrazone',
                                   [(carbonyl, None),
                                    ('aryl_hydrazine', {1: 3, 2: 4, 3: 5, 4: 6})],
                                   '[A:1]=[A:4]-[A:3]-[A:5]:[A:6]'))

    return rules


rules = Proxy(_rules)


__all__ = ['rules']
