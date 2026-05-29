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
                                                  delete_atoms=True, one_shot=True, fix_broken_pyrroles=True)


def _rules():
    rules = []

    # 1,3-diketone → isoxazole (+ hydroxylamine, implicit)
    rules.append(_make_reactor('isoxazole', '1_3_diketone', 'isoxazole', '[A:1]:1:[O:20]:[N:21]:[A:3]:[A:5]:1'))

    # 1,4-diketone → pyridazine (+ hydrazine, implicit)
    rules.append(_make_reactor('pyridazine', '1_4_diketone', 'pyridazine', '[A:1]:1:[N:20]:[N:21]:[A:3]:[A:6]:[A:5]:1'))

    # Appel: alcohol → alkyl halide (CCl4/CBr4 + PPh3)
    rules.append(_make_reactor('appel', 'primary_alcohol', 'alkyl_bromide', '[Br:20]-[A:2]'))
    rules.append(_make_reactor('appel', 'secondary_alcohol', 'alkyl_bromide', '[Br:20]-[A:2]'))
    rules.append(_make_reactor('appel_chloride', 'primary_alcohol', 'alkyl_chloride', '[Cl:20]-[A:2]'))
    rules.append(_make_reactor('appel_chloride', 'secondary_alcohol', 'alkyl_chloride', '[Cl:20]-[A:2]'))

    # Miyaura borylation: ArX → ArB(OH)2
    rules.append(_make_reactor('borylation_acid', 'aryl_bromide', 'aryl_boronic_acid', '[A:1]-[B:20](-[O:21])-[O:22]'))
    rules.append(_make_reactor('borylation_acid', 'aryl_iodide', 'aryl_boronic_acid', '[A:1]-[B:20](-[O:21])-[O:22]'))
    rules.append(_make_reactor('borylation_acid', 'aryl_chloride', 'aryl_boronic_acid', '[A:1]-[B:20](-[O:21])-[O:22]'))

    # Miyaura borylation: ArX → ArBpin (Pd, B2pin2)
    rules.append(_make_reactor('borylation_ester', 'aryl_bromide', 'aryl_boronic_ester',
                               '[A:1]-[B:20]1-[O:21]-[C:22]([C:23])([C:24])-[C:25]([C:26])([C:27])-[O:28]-1'))
    rules.append(_make_reactor('borylation_ester', 'aryl_iodide', 'aryl_boronic_ester',
                               '[A:1]-[B:20]1-[O:21]-[C:22]([C:23])([C:24])-[C:25]([C:26])([C:27])-[O:28]-1'))
    rules.append(_make_reactor('borylation_ester', 'aryl_chloride', 'aryl_boronic_ester',
                               '[A:1]-[B:20]1-[O:21]-[C:22]([C:23])([C:24])-[C:25]([C:26])([C:27])-[O:28]-1'))

    # Nitrile hydrolysis: R-C≡N → R-C(=O)-NH2 (partial, to amide)
    rules.append(_make_reactor('nitrile_hydrolysis', 'nitrile', 'primary_amide', '[A:2]-[A:1]=[O:20]'))

    # Electrophilic aromatic nitration: Ar-H → Ar-NO2
    rules.append(_make_reactor('nitration', 'arene_ch', 'nitro', '[A:1]-[N;+:20](=[O:21])[O;-:22]'))

    # Electrophilic aromatic halogenation: Ar-H → Ar-X
    rules.append(_make_reactor('bromination', 'arene_ch', 'aryl_bromide', '[A:1]-[Br:20]'))
    rules.append(_make_reactor('chlorination', 'arene_ch', 'aryl_chloride', '[A:1]-[Cl:20]'))
    rules.append(_make_reactor('iodination', 'arene_ch', 'aryl_iodide', '[A:1]-[I:20]'))

    # benzylic halogenation: Ar-CH2R → Ar-CHBrR (NBS, Br2/hv)
    rules.append(_make_reactor('benzylic_bromination', 'benzylic_ch', 'alkyl_bromide', '[A:1](-[Br:20])-[A:2]'))

    # alpha-halogenation: alpha_ketone → alpha-bromoketone (NBS/Br2)
    rules.append(_make_reactor('alpha_halogenation', 'alpha_ketone', 'alpha_haloketone', '[A:1](=[A:2])-[A:3](-[Br:20])'))

    # alkyl borylation: R-X → R-B(OH)2
    rules.append(_make_reactor('borylation_acid', 'alkyl_bromide', 'alkyl_boronic_acid', '[A:1]-[B:20](-[O:21])-[O:22]'))
    rules.append(_make_reactor('borylation_acid', 'alkyl_iodide', 'alkyl_boronic_acid', '[A:1]-[B:20](-[O:21])-[O:22]'))
    rules.append(_make_reactor('borylation_acid', 'alkyl_chloride', 'alkyl_boronic_acid', '[A:1]-[B:20](-[O:21])-[O:22]'))

    # alkyl borylation: R-X → R-Bpin
    rules.append(_make_reactor('borylation_ester', 'alkyl_bromide', 'alkyl_boronic_ester',
                               '[A:1]-[B:20]1-[O:21]-[C:22]([C:23])([C:24])-[C:25]([C:26])([C:27])-[O:28]-1'))
    rules.append(_make_reactor('borylation_ester', 'alkyl_iodide', 'alkyl_boronic_ester',
                               '[A:1]-[B:20]1-[O:21]-[C:22]([C:23])([C:24])-[C:25]([C:26])([C:27])-[O:28]-1'))
    rules.append(_make_reactor('borylation_ester', 'alkyl_chloride', 'alkyl_boronic_ester',
                               '[A:1]-[B:20]1-[O:21]-[C:22]([C:23])([C:24])-[C:25]([C:26])([C:27])-[O:28]-1'))

    # ester to amide: R-C(=O)-OR' → R-C(=O)-NH2 (NH3)
    rules.append(_make_reactor('ester_to_amide', 'ester', 'primary_amide', '[A:1](=[A:2])-[N:20]'))

    # ester to hydroxamic acid: R-C(=O)-OR' → R-C(=O)-NH-OH (NH2OH)
    rules.append(_make_reactor('ester_to_hydroxamic_acid', 'ester', 'hydroxamic_acid', '[A:1](=[A:2])-[N:20]-[O:21]'))

    # amide hydrolysis: R-C(=O)-NHR → R-C(=O)-OH
    rules.append(_make_reactor('amide_hydrolysis', 'primary_amide', 'carboxylic_acid', '[A:2](=[A:3])-[O:20]'))
    rules.append(_make_reactor('amide_hydrolysis', 'secondary_amide', 'carboxylic_acid', '[A:2](=[A:3])-[O:20]'))

    # cyanation: R-X → R-CN (NaCN/CuCN)
    rules.append(_make_reactor('cyanation', 'alkyl_bromide', 'nitrile', '[A:1]-[C:20]#[N:21]'))
    rules.append(_make_reactor('cyanation', 'alkyl_iodide', 'nitrile', '[A:1]-[C:20]#[N:21]'))
    rules.append(_make_reactor('cyanation', 'alkyl_chloride', 'nitrile', '[A:1]-[C:20]#[N:21]'))
    rules.append(_make_reactor('cyanation', 'aryl_bromide', 'nitrile', '[A:1]-[C:20]#[N:21]'))
    rules.append(_make_reactor('cyanation', 'aryl_iodide', 'nitrile', '[A:1]-[C:20]#[N:21]'))
    rules.append(_make_reactor('cyanation', 'aryl_chloride', 'nitrile', '[A:1]-[C:20]#[N:21]'))

    # dehydration: R-C(=O)-NH2 → R-C≡N (P2O5/SOCl2)
    rules.append(_make_reactor('dehydration', 'primary_amide', 'nitrile', '[A:2]#[A:1]'))

    # triflation: phenol/alcohol → triflate (Tf2O)
    rules.append(_make_reactor('triflation', 'phenol', 'aryl_triflate', '[A:2]-[A:1]-[S:20](=[O:21])(=[O:22])-[C:23]([F:24])([F:25])[F:26]'))
    rules.append(_make_reactor('triflation', 'primary_alcohol', 'alkyl_triflate', '[A:2]-[A:1]-[S:20](=[O:21])(=[O:22])-[C:23]([F:24])([F:25])[F:26]'))
    rules.append(_make_reactor('triflation', 'secondary_alcohol', 'alkyl_triflate', '[A:2]-[A:1]-[S:20](=[O:21])(=[O:22])-[C:23]([F:24])([F:25])[F:26]'))

    # nitrile full hydrolysis: R-C≡N → R-COOH
    rules.append(_make_reactor('nitrile_to_acid', 'nitrile', 'carboxylic_acid', '[A:1](=[O:20])-[O:21]'))

    # tertiary alcohol dehydration: R3C-OH → R2C=CR (loss of water)
    rules.append(_make_reactor('alcohol_dehydration', 'tertiary_alcohol_with_alpha_h', 'alkene', '[A:2]=[A:3]'))

    # oximation: ketone/aldehyde → oxime (+ NH2OH implicit)
    rules.append(_make_reactor('oximation', 'ketone', 'oxime', '[O:20]-[N:21]=[A:1]'))
    rules.append(_make_reactor('oximation', 'aldehyde', 'oxime', '[O:20]-[N:21]=[A:1]'))

    # acid chloride formation: RCOOH → RCOCl (SOCl2, oxalyl chloride)
    rules.append(_make_reactor('acid_chlorination', 'carboxylic_acid', 'acyl_chloride', '[Cl:20]-[A:1]=[A:2]'))

    return rules


rules = Proxy(_rules)


__all__ = ['rules']
