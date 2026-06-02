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
    rules.append(_make_reactor('borylation_acid', 'aryl_triflate', 'aryl_boronic_acid', '[A:1]-[B:20](-[O:21])-[O:22]'))
    rules.append(_make_reactor('borylation_ester', 'aryl_triflate', 'aryl_boronic_ester',
                               '[A:1]-[B:20]1-[O:21]-[C:22]([C:23])([C:24])-[C:25]([C:26])([C:27])-[O:28]-1'))

    # Nitrile hydrolysis: R-C≡N → R-C(=O)-NH2 (partial, to amide)
    rules.append(_make_reactor('nitrile_hydrolysis', 'nitrile', 'primary_amide', '[A:2]-[A:1]=[O:20]'))

    # Electrophilic aromatic nitration: Ar-H → Ar-NO2
    rules.append(_make_reactor('nitration', 'arene_ch', 'nitro', '[A:1]-[N;+:20](=[O:21])[O;-:22]'))

    # Electrophilic aromatic halogenation: Ar-H → Ar-X
    rules.append(_make_reactor('bromination', 'arene_ch', 'aryl_bromide', '[A:1]-[Br:20]'))
    rules.append(_make_reactor('chlorination', 'arene_ch', 'aryl_chloride', '[A:1]-[Cl:20]'))
    rules.append(_make_reactor('iodination', 'arene_ch', 'aryl_iodide', '[A:1]-[I:20]'))
    rules.append(_make_reactor('fluorination', 'arene_ch', 'aryl_fluoride', '[A:1]-[F:20]'))

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

    # mesylation: ROH → ROMs (MsCl, Et3N)
    rules.append(_make_reactor('mesylation', 'primary_alcohol', 'alkyl_mesylate', '[A:2]-[A:1]-[S:20](=[O:21])(=[O:22])-[C:23]'))
    rules.append(_make_reactor('mesylation', 'secondary_alcohol', 'alkyl_mesylate', '[A:2]-[A:1]-[S:20](=[O:21])(=[O:22])-[C:23]'))

    # tosylation: ROH → ROTs (TsCl, pyridine)
    rules.append(_make_reactor('tosylation', 'primary_alcohol', 'alkyl_tosylate',
                               '[A:2]-[A:1]-[S:20](=[O:21])(=[O:22])-[C:23]:1:[C:24]:[C:25]:[C:26](-[C:27]):[C:28]:[C:29]:1'))
    rules.append(_make_reactor('tosylation', 'secondary_alcohol', 'alkyl_tosylate',
                               '[A:2]-[A:1]-[S:20](=[O:21])(=[O:22])-[C:23]:1:[C:24]:[C:25]:[C:26](-[C:27]):[C:28]:[C:29]:1'))

    # halide hydrolysis: ArX → ArOH (NaOH/H2O, Cu cat, high T)
    rules.append(_make_reactor('halide_hydrolysis', 'aryl_chloride', 'phenol', '[A:1]-[O:20]'))
    rules.append(_make_reactor('halide_hydrolysis', 'aryl_bromide', 'phenol', '[A:1]-[O:20]'))
    rules.append(_make_reactor('halide_hydrolysis', 'aryl_iodide', 'phenol', '[A:1]-[O:20]'))
    rules.append(_make_reactor('halide_hydrolysis', 'aryl_fluoride', 'phenol', '[A:1]-[O:20]'))

    # alkyl halide hydrolysis: RX → ROH (SN2, NaOH/H2O)
    rules.append(_make_reactor('halide_hydrolysis', 'alkyl_chloride', 'primary_alcohol', '[A:1]-[O:20]'))
    rules.append(_make_reactor('halide_hydrolysis', 'alkyl_bromide', 'primary_alcohol', '[A:1]-[O:20]'))
    rules.append(_make_reactor('halide_hydrolysis', 'alkyl_iodide', 'primary_alcohol', '[A:1]-[O:20]'))

    # DAST fluorination: ROH → RF (DAST, DeoxoFluor)
    rules.append(_make_reactor('fluorination', 'primary_alcohol', 'alkyl_fluoride', '[F:20]-[A:2]'))
    rules.append(_make_reactor('fluorination', 'secondary_alcohol', 'alkyl_fluoride', '[F:20]-[A:2]'))
    rules.append(_make_reactor('fluorination', 'tertiary_alcohol', 'alkyl_fluoride', '[F:20]-[A:2]'))

    # Appel iodination: ROH → RI (I2/PPh3)
    rules.append(_make_reactor('appel_iodide', 'primary_alcohol', 'alkyl_iodide', '[I:20]-[A:2]'))
    rules.append(_make_reactor('appel_iodide', 'secondary_alcohol', 'alkyl_iodide', '[I:20]-[A:2]'))

    # Hydroxy to chloro: ROH → RCl (SOCl2, PCl5 — broader than Appel)
    # phenol → ArCl (PCl5, POCl3)
    rules.append(_make_reactor('hydroxy_to_chloro', 'phenol', 'aryl_chloride', '[Cl:20]-[A:2]'))
    # pyridol → chloropyridine (POCl3) — covers pyrimidinol, pyridazinol, etc.
    rules.append(_make_reactor('hydroxy_to_chloro', 'pyridol', 'aryl_chloride', '[A:1]:[A:2]-[Cl:20]'))
    # tertiary_alcohol → RCl
    rules.append(_make_reactor('hydroxy_to_chloro', 'tertiary_alcohol', 'alkyl_chloride', '[Cl:20]-[A:2]'))

    # acid fluoride formation: RCOOH → RCOF (cyanuric fluoride, DAST)
    rules.append(_make_reactor('acid_fluorination', 'carboxylic_acid', 'acyl_fluoride', '[F:20]-[A:1]=[A:2]'))

    # Sandmeyer: ArNH2 → ArX (NaNO2/HX, CuX)
    rules.append(_make_reactor('sandmeyer_bromination', 'primary_aniline', 'aryl_bromide', '[Br:20]-[A:2]'))
    rules.append(_make_reactor('sandmeyer_iodination', 'primary_aniline', 'aryl_iodide', '[I:20]-[A:2]'))
    rules.append(_make_reactor('sandmeyer_chlorination', 'primary_aniline', 'aryl_chloride', '[Cl:20]-[A:2]'))

    # Oxo to difluoro: C=O → CF2 (DAST, DeoxoFluor)
    rules.append(_make_reactor('oxo_to_difluoro', 'ketone', None, '[A:1](-[F:20])-[F:21]'))
    rules.append(_make_reactor('oxo_to_difluoro', 'aldehyde', None, '[A:1](-[F:20])-[F:21]'))

    # Oxo to thioxo: C=O → C=S (Lawesson's reagent, P4S10)
    rules.append(_make_reactor('oxo_to_thioxo', 'primary_amide', 'thioamide', '[S:20]=[A:2]-[A:1]'))
    rules.append(_make_reactor('oxo_to_thioxo', 'secondary_amide', None, '[S:20]=[A:2]-[A:1]'))
    rules.append(_make_reactor('oxo_to_thioxo', 'ketone', None, '[S:20]=[A:1]'))

    # Amino to isocyanate: ArNH2/RNH2 → ArNCO (triphosgene/COCl2)
    # primary_aniline: [N:1][C;a:2] → need N=C=O, keep Ar bond
    rules.append(_make_reactor('amino_to_isocyanate', 'primary_aniline', 'isocyanate', '[A:2]-[A:1]=[C:20]=[O:21]'))
    rules.append(_make_reactor('amino_to_isocyanate', 'primary_amine', 'isocyanate', '[A:2]-[A:1]=[C:20]=[O:21]'))

    # Amino to isothiocyanate: ArNH2/RNH2 → ArNCS (CS2/DCC, thiophosgene)
    rules.append(_make_reactor('amino_to_isothiocyanate', 'primary_aniline', None, '[A:2]-[A:1]=[C:20]=[S:21]'))
    rules.append(_make_reactor('amino_to_isothiocyanate', 'primary_amine', None, '[A:2]-[A:1]=[C:20]=[S:21]'))

    # Amino to azide: ArNH2 → ArN3 (diazotization + NaN3)
    rules.append(_make_reactor('amino_to_azide', 'primary_aniline', 'azide', '[A:2]-[N:20]=[N;+:21]=[N;-:22]'))
    rules.append(_make_reactor('amino_to_azide', 'primary_amine', 'azide', '[A:2]-[N:20]=[N;+:21]=[N;-:22]'))

    # Amino to hydroxy: ArNH2 → ArOH (diazotization + H2O)
    rules.append(_make_reactor('amino_to_hydroxy', 'primary_aniline', 'phenol', '[A:2]-[O:20]'))
    rules.append(_make_reactor('amino_to_hydroxy', 'primary_amine', 'primary_alcohol', '[A:2]-[O:20]'))

    # Bpin hydrolysis: ArBpin → ArB(OH)2 (NaIO4, HCl, diethanolamine)
    rules.append(_make_reactor('bpin_hydrolysis', 'aryl_boronic_ester', 'aryl_boronic_acid', '[A:1]-[B:20](-[O:21])-[O:22]'))
    rules.append(_make_reactor('bpin_hydrolysis', 'alkyl_boronic_ester', 'alkyl_boronic_acid', '[A:1]-[B:20](-[O:21])-[O:22]'))

    # Ester hydrolysis: RCOOR' → RCOOH (LiOH, NaOH, HCl)
    rules.append(_make_reactor('ester_hydrolysis', 'ester', 'carboxylic_acid', '[A:1](=[A:2])-[O:20]'))

    # Carboxy to amide: RCOOH → RCONH2 (EDC/HOBt+NH3, or CDI+NH3)
    rules.append(_make_reactor('acid_to_amide', 'carboxylic_acid', 'primary_amide', '[N:20]-[A:1]=[A:2]'))

    # Acyl chloride to acid: RCOCl → RCOOH (H2O)
    rules.append(_make_reactor('acyl_chloride_hydrolysis', 'acyl_chloride', 'carboxylic_acid', '[O:20]-[A:1]=[A:2]'))
    rules.append(_make_reactor('acyl_chloride_hydrolysis', 'acyl_fluoride', 'carboxylic_acid', '[O:20]-[A:1]=[A:2]'))

    # Acyl chloride to amide: RCOCl → RCONH2 (NH3)
    rules.append(_make_reactor('acyl_chloride_to_amide', 'acyl_chloride', 'primary_amide', '[N:20]-[A:1]=[A:2]'))
    rules.append(_make_reactor('acyl_chloride_to_amide', 'acyl_fluoride', 'primary_amide', '[N:20]-[A:1]=[A:2]'))

    # Halide to methoxy: ArX/RX → ArOMe/ROMe (NaOMe, SNAr)
    rules.append(_make_reactor('halide_to_methoxy', 'aryl_fluoride', None, '[A:1]-[O:20]-[C:21]'))
    rules.append(_make_reactor('halide_to_methoxy', 'aryl_chloride', None, '[A:1]-[O:20]-[C:21]'))
    rules.append(_make_reactor('halide_to_methoxy', 'aryl_bromide', None, '[A:1]-[O:20]-[C:21]'))
    rules.append(_make_reactor('halide_to_methoxy', 'alkyl_chloride', None, '[A:1]-[O:20]-[C:21]'))
    rules.append(_make_reactor('halide_to_methoxy', 'alkyl_bromide', None, '[A:1]-[O:20]-[C:21]'))
    rules.append(_make_reactor('halide_to_methoxy', 'alkyl_iodide', None, '[A:1]-[O:20]-[C:21]'))

    # Halide to azide: RX → RN3 (NaN3)
    rules.append(_make_reactor('azidation', 'alkyl_chloride', 'azide', '[A:1]-[N:20]=[N;+:21]=[N;-:22]'))
    rules.append(_make_reactor('azidation', 'alkyl_bromide', 'azide', '[A:1]-[N:20]=[N;+:21]=[N;-:22]'))
    rules.append(_make_reactor('azidation', 'alkyl_iodide', 'azide', '[A:1]-[N:20]=[N;+:21]=[N;-:22]'))
    rules.append(_make_reactor('azidation', 'aryl_chloride', 'azide', '[A:1]-[N:20]=[N;+:21]=[N;-:22]'))
    rules.append(_make_reactor('azidation', 'aryl_bromide', 'azide', '[A:1]-[N:20]=[N;+:21]=[N;-:22]'))
    rules.append(_make_reactor('azidation', 'aryl_fluoride', 'azide', '[A:1]-[N:20]=[N;+:21]=[N;-:22]'))

    # Hydroxy to azide: ROH → RN3 (DPPA/Mitsunobu or Ms→N3)
    rules.append(_make_reactor('azidation', 'primary_alcohol', 'azide', '[A:2]-[N:20]=[N;+:21]=[N;-:22]'))
    rules.append(_make_reactor('azidation', 'secondary_alcohol', 'azide', '[A:2]-[N:20]=[N;+:21]=[N;-:22]'))

    # Halide exchange (Finkelstein): RCl/RBr → RI (NaI/acetone)
    rules.append(_make_reactor('finkelstein', 'alkyl_chloride', 'alkyl_iodide', '[A:1]-[I:20]'))
    rules.append(_make_reactor('finkelstein', 'alkyl_bromide', 'alkyl_iodide', '[A:1]-[I:20]'))
    rules.append(_make_reactor('finkelstein', 'aryl_chloride', 'aryl_iodide', '[A:1]-[I:20]'))
    rules.append(_make_reactor('finkelstein', 'aryl_bromide', 'aryl_iodide', '[A:1]-[I:20]'))

    # Halex fluorination: ArCl/ArBr → ArF (KF, CsF, spray-dried KF)
    rules.append(_make_reactor('halex', 'aryl_chloride', 'aryl_fluoride', '[A:1]-[F:20]'))
    rules.append(_make_reactor('halex', 'aryl_bromide', 'aryl_fluoride', '[A:1]-[F:20]'))

    # Halide to chloride: ArBr/ArI → ArCl; RI → RCl
    rules.append(_make_reactor('halide_to_chloride', 'aryl_bromide', 'aryl_chloride', '[A:1]-[Cl:20]'))
    rules.append(_make_reactor('halide_to_chloride', 'aryl_iodide', 'aryl_chloride', '[A:1]-[Cl:20]'))
    rules.append(_make_reactor('halide_to_chloride', 'alkyl_iodide', 'alkyl_chloride', '[A:1]-[Cl:20]'))

    # Halide to bromide: ArCl → ArBr (not common but exists)
    rules.append(_make_reactor('halide_to_bromide', 'aryl_chloride', 'aryl_bromide', '[A:1]-[Br:20]'))

    # Halide to thiol: RX → RSH (NaSH, thiourea hydrolysis)
    rules.append(_make_reactor('halide_to_thiol', 'aryl_chloride', None, '[A:1]-[S:20]'))
    rules.append(_make_reactor('halide_to_thiol', 'aryl_bromide', None, '[A:1]-[S:20]'))
    rules.append(_make_reactor('halide_to_thiol', 'alkyl_bromide', 'thiol', '[A:1]-[S:20]'))
    rules.append(_make_reactor('halide_to_thiol', 'alkyl_chloride', 'thiol', '[A:1]-[S:20]'))

    # Halide to thioether: RX → RSMe (NaSMe)
    rules.append(_make_reactor('halide_to_thioether', 'aryl_chloride', 'thioether', '[A:1]-[S:20]-[C:21]'))
    rules.append(_make_reactor('halide_to_thioether', 'aryl_bromide', 'thioether', '[A:1]-[S:20]-[C:21]'))
    rules.append(_make_reactor('halide_to_thioether', 'alkyl_bromide', 'thioether', '[A:1]-[S:20]-[C:21]'))
    rules.append(_make_reactor('halide_to_thioether', 'alkyl_chloride', 'thioether', '[A:1]-[S:20]-[C:21]'))

    return rules


rules = Proxy(_rules)


__all__ = ['rules']
