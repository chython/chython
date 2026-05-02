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
import pytest
from chython import smiles


# (test_id, reaction_name, reaction_smiles, expected_count)
# Reaction SMILES are parsed, canonicalized, then reactants are fed to @ operator.
# expected_count: number of results from @ (2 for tautomeric azoles).
_two_component = [
    ('amidation_acyl_chloride_primary_amine',
     'amidation', 'ClC(=O)C.NCC>>CCNC(=O)C', 1),
    ('suzuki_aryl_bromide_boronic_acid',
     'suzuki', 'Brc1ccc(cc1)C.OB(O)c1ccccc1>>c1cc(ccc1-c1ccccc1)C', 1),
    ('suzuki_alkenyl_boronic_acid',
     'suzuki', 'Brc1ccc(cc1)C.OB(O)C=CC>>c1cc(ccc1C=CC)C', 1),
    ('suzuki_alkyl_boronic_acid',
     'suzuki', 'Brc1ccc(cc1)C.OB(O)CCC>>c1cc(ccc1CCC)C', 1),
    ('buchwald_hartwig_aryl_bromide_primary_amine',
     'buchwald_hartwig', 'Brc1ccc(cc1)C.NCC>>c1cc(ccc1C)NCC', 1),
    ('mitsunobu_ether_primary_alcohol_phenol',
     'mitsunobu', 'OCC.Oc1ccc(cc1)C>>c1cc(ccc1C)OCC', 1),
    ('mitsunobu_ester_primary_alcohol_acid',
     'mitsunobu', 'OCC.OC(=O)c1ccccc1>>c1ccccc1C(OCC)=O', 1),
    ('deoxygenative_coupling_alcohol_aryl_bromide',
     'deoxygenative_coupling', 'OCC.Brc1ccc(cc1)C>>c1cc(ccc1C)CC', 1),
    ('decarboxylative_coupling_acid_aryl_bromide',
     'decarboxylative_coupling', 'OC(=O)CC.Brc1ccc(cc1)C>>c1cc(ccc1C)CC', 1),
    ('xec_aryl_bromide_alkyl_bromide',
     'xec', 'Brc1ccc(cc1)C.BrCCC>>c1cc(ccc1CCC)C', 1),
    ('reductive_amination_aldehyde_primary_amine',
     'reductive_amination', 'O=CC.NCc1ccccc1>>c1cc(ccc1)CNCC', 1),
    ('reductive_amination_ketone_secondary_amine',
     'reductive_amination', 'O=C(C)C.N(C)CC>>CC(C)N(C)CC', 1),
    ('ullmann_phenol_aryl_bromide_cresol',
     'ullmann_phenol', 'Brc1ccc(cc1)C.Oc1ccc(cc1)CC>>c1cc(ccc1C)Oc1ccc(cc1)CC', 1),
    ('ullmann_pyrrole_aryl_bromide_pyrrole',
     'ullmann_pyrrole', 'Brc1ccc(cc1)C.[nH]1cccc1>>c1ccn(c1)-c1ccc(cc1)C', 1),
    ('ullmann_pyrrole_aryl_bromide_pyrazole',
     'ullmann_pyrrole', 'Brc1ccc(cc1)C.[nH]1nccc1>>c1cc(ccc1-n1nccc1)C', 2),
    ('ullmann_pyrrole_aryl_bromide_imidazole',
     'ullmann_pyrrole', 'Brc1ccc(cc1)C.[nH]1cncc1>>c1cc(ccc1-n1ccnc1)C', 2),
    ('chan_lam_boronic_acid_primary_amine',
     'chan_lam', 'OB(O)c1ccc(cc1)C.NCC>>c1cc(ccc1C)NCC', 1),
    ('chan_lam_boronic_acid_phenol',
     'chan_lam', 'OB(O)c1ccc(cc1)C.Oc1ccc(cc1)CC>>c1cc(ccc1C)Oc1ccc(cc1)CC', 1),
    ('sulfonylation_sulfonyl_chloride_alcohol',
     'sulfonylation', 'ClS(=O)(=O)c1ccccc1.OCCC>>c1ccccc1S(=O)(=O)OCCC', 1),
    ('sulfonamide_formation_sulfonyl_chloride_amine',
     'sulfonamide_formation', 'ClS(=O)(=O)c1ccccc1.NCC>>c1ccccc1S(=O)(=O)NCC', 1),
    ('aminolysis_ester_primary_amine',
     'aminolysis', 'COC(=O)c1ccc(cc1)C.NCC>>c1cc(ccc1C(NCC)=O)C', 1),
    ('grignard_alkyl_bromide_aldehyde',
     'grignard', 'BrCC.O=Cc1ccc(cc1)C>>c1cc(ccc1C)C(CC)O', 1),
    ('grignard_alkyl_bromide_ketone',
     'grignard', 'BrCC.O=C(C)c1ccccc1>>c1cc(ccc1)C(O)(C)CC', 1),
    ('sonogashira_aryl_bromide_terminal_alkyne',
     'sonogashira', 'Brc1ccc(cc1)C.C#CC>>c1cc(ccc1C)C#CC', 1),
    ('urea_synthesis_isocyanate_primary_amine',
     'urea_synthesis', 'O=C=Nc1ccccc1.NCC>>c1ccccc1NC(NCC)=O', 1),
    ('cuaac_azide_terminal_alkyne',
     'cuaac', 'CC#C.CN=[N+]=[N-]>>Cc1cn(C)nn1', 1),
    ('snar_aryl_fluoride_primary_amine',
     'snar', 'Fc1ccc(cc1)C(F)(F)F.NCC>>c1cc(ccc1C(F)(F)F)NCC', 1),
    ('snar_aryl_fluoride_phenol',
     'snar', 'Fc1ccc(cc1)C(F)(F)F.Oc1ccc(cc1)C>>c1cc(ccc1C(F)(F)F)Oc1ccc(cc1)C', 1),
    ('snar_aryl_fluoride_thiol',
     'snar', 'Fc1ccc(cc1)C(F)(F)F.SCC>>c1cc(ccc1C(F)(F)F)SCC', 1),
    ('williamson_alkyl_bromide_phenol',
     'williamson', 'BrCC.Oc1ccc(cc1)C>>c1cc(ccc1C)OCC', 1),
    ('williamson_alkyl_bromide_alcohol',
     'williamson', 'BrCCC.OCC>>CCCOCC', 1),
    ('acylation_acyl_chloride_alcohol',
     'acylation', 'ClC(=O)c1ccccc1.OCC>>c1ccccc1C(OCC)=O', 1),
    ('thioether_thiol_alkyl_bromide',
     'thioether', 'SCC.BrCCC>>CCSCCC', 1),
    ('carbamate_isocyanate_alcohol',
     'carbamate', 'O=C=Nc1ccccc1.OCC>>c1ccccc1NC(OCC)=O', 1),
    ('tetrazole_azide_nitrile',
     'tetrazole', 'CN=[N+]=[N-].N#Cc1ccccc1>>c1ccc(cc1)-c1nnn(C)n1', 1),
    ('kumada_aryl_bromide_alkyl_grignard',
     'kumada', 'Brc1ccc(cc1)C.C[Mg]Br>>c1cc(ccc1C)C', 1),
    ('kumada_aryl_bromide_aryl_grignard',
     'kumada', 'Brc1ccc(cc1)C.c1ccccc1[Mg]Br>>c1cc(ccc1-c1ccccc1)C', 1),
    ('negishi_aryl_bromide_alkyl_zinc',
     'negishi', 'Brc1ccc(cc1)C.CC[Zn]Cl>>c1cc(ccc1C)CC', 1),
    ('stille_aryl_bromide_aryl_stannane',
     'stille', 'Brc1ccc(cc1)C.c1ccccc1[Sn](C)(C)C>>c1cc(ccc1-c1ccccc1)C', 1),
    ('stille_aryl_bromide_alkenyl_stannane',
     'stille', 'Brc1ccc(cc1)C.C=C[Sn](C)(C)C>>c1cc(ccc1C=C)C', 1),
    ('hiyama_aryl_bromide_aryl_silane',
     'hiyama', 'Brc1ccc(cc1)C.c1ccccc1[Si](C)(C)C>>c1cc(ccc1-c1ccccc1)C', 1),
    ('grignard_explicit_alkyl_grignard_aldehyde',
     'grignard', 'CC[Mg]Br.O=Cc1ccccc1>>c1ccccc1C(CC)O', 1),
    ('heck_aryl_bromide_terminal_alkene',
     'heck', 'Brc1ccc(cc1)C.C=CC>>c1cc(ccc1C=CC)C', 1),
    ('cross_metathesis_two_terminal_alkenes',
     'cross_metathesis', 'C=CC.C=CCCC>>CC=CCCC', 1),
    ('wittig_ylide_aldehyde',
     'wittig', 'c1ccccc1P(c1ccccc1)(c1ccccc1)=C.c1cc(ccc1C=O)C>>c1cc(ccc1C)C=C', 1),
    ('hwe_phosphonate_aldehyde',
     'hwe', 'CCP(=O)(OC)OC.c1cc(ccc1C=O)C>>c1cc(ccc1C=CC)C', 1),
    ('weinreb_explicit_grignard',
     'weinreb', 'CC(=O)N(C)OC.CC[Mg]Br>>CCC(=O)C', 1),
    ('weinreb_ester_surrogate_grignard',
     'weinreb', 'COC(=O)c1ccccc1.CC[Mg]Br>>c1ccccc1C(=O)CC', 1),
    ('weinreb_acyl_chloride_surrogate_grignard',
     'weinreb', 'c1ccccc1C(=O)Cl.CC[Mg]Br>>c1ccccc1C(=O)CC', 1),
    ('friedel_crafts_acyl_chloride_benzene',
     'friedel_crafts', 'CC(=O)Cl.c1ccccc1>>c1ccccc1C(=O)C', 1),
]

_three_component = [
    ('ugi_3cr_aldehyde_amine_isocyano',
     'ugi_3cr', 'O=CC.NCC.[C-]#[N+]C>>CCNC(C(=O)NC)C', 1),
]

_oxidations = [
    ('alcohol_to_aldehyde', 'alcohol_to_aldehyde', 'OCC>>O=CC', 1),
    ('alcohol_to_ketone', 'alcohol_to_ketone', 'OC(C)C>>O=C(C)C', 1),
    ('aldehyde_to_acid', 'aldehyde_to_acid', 'O=CC>>OC(=O)C', 1),
    ('dihydroxylation', 'dihydroxylation', 'C=Cc1ccccc1>>OC(CO)c1ccccc1', 1),
]

_reductions = [
    ('aldehyde_to_alcohol', 'aldehyde_to_alcohol', 'O=Cc1ccccc1>>OCc1ccccc1', 1),
    ('ketone_to_alcohol', 'ketone_to_alcohol', 'O=C(C)c1ccccc1>>OC(C)c1ccccc1', 1),
    ('acid_to_alcohol', 'acid_to_alcohol', 'OC(=O)c1ccccc1>>OCc1ccccc1', 1),
    ('ester_to_alcohol', 'ester_to_alcohol', 'COC(=O)c1ccccc1>>OCc1ccccc1', 1),
    ('amide_to_amine_primary', 'amide_to_amine', 'NC(=O)c1ccccc1>>NCc1ccccc1', 1),
    ('amide_to_amine_secondary', 'amide_to_amine', 'CNC(=O)c1ccccc1>>CNCc1ccccc1', 1),
    ('nitrile_to_amine', 'nitrile_to_amine', 'N#Cc1ccccc1>>NCc1ccccc1', 1),
    ('nitro_to_amine', 'nitro_to_amine', '[O-][N+](=O)c1ccccc1>>Nc1ccccc1', 1),
    ('azide_to_amine', 'azide_to_amine', 'CN=[N+]=[N-]>>CN', 1),
    ('deoxygenation_primary', 'deoxygenation', 'OCC1CCCCC1>>CC1CCCCC1', 1),
    ('deoxygenation_secondary', 'deoxygenation', 'OC(C)c1ccccc1>>CCc1ccccc1', 1),
]


@pytest.mark.parametrize(
    'test_id,rxn_name,rxn_smi,expected_count',
    _two_component,
    ids=[x[0] for x in _two_component]
)
def test_two_component(test_id, rxn_name, rxn_smi, expected_count):
    r = smiles(rxn_smi)
    r.canonicalize()
    r1, r2 = r.reactants

    results = list(r1 @ r2)

    # only one reaction name should fire
    names = set(n for n, _ in results)
    assert names == {rxn_name}, f'expected {rxn_name}, got {sorted(names)}'

    # expected number of results
    assert len(results) == expected_count, f'expected {expected_count}, got {len(results)}'

    # all products match expected
    for _, rxn in results:
        assert rxn == r, f'expected {r}, got {rxn}'


@pytest.mark.parametrize(
    'test_id,rxn_name,rxn_smi,expected_count',
    _three_component,
    ids=[x[0] for x in _three_component]
)
def test_three_component(test_id, rxn_name, rxn_smi, expected_count):
    r = smiles(rxn_smi)
    r.canonicalize()
    r1, r2, r3 = r.reactants

    results = list(r1 @ [r2, r3])

    names = set(n for n, _ in results)
    assert names == {rxn_name}, f'expected {rxn_name}, got {sorted(names)}'

    assert len(results) == expected_count, f'expected {expected_count}, got {len(results)}'

    for _, rxn in results:
        assert rxn == r, f'expected {r}, got {rxn}'


@pytest.mark.parametrize(
    'test_id,rxn_name,rxn_smi,expected_count',
    _oxidations,
    ids=[x[0] for x in _oxidations]
)
def test_oxidation(test_id, rxn_name, rxn_smi, expected_count):
    r = smiles(rxn_smi)
    r.canonicalize()
    mol = r.reactants[0]

    results = list(mol.oxidize())

    names = set(n for n, _ in results)
    assert names == {rxn_name}, f'expected {rxn_name}, got {sorted(names)}'

    assert len(results) == expected_count, f'expected {expected_count}, got {len(results)}'

    for _, rxn in results:
        assert rxn == r, f'expected {r}, got {rxn}'


@pytest.mark.parametrize(
    'test_id,rxn_name,rxn_smi,expected_count',
    _reductions,
    ids=[x[0] for x in _reductions]
)
def test_reduction(test_id, rxn_name, rxn_smi, expected_count):
    r = smiles(rxn_smi)
    r.canonicalize()
    mol = r.reactants[0]

    results = list(mol.reduce())

    names = set(n for n, _ in results)
    assert names == {rxn_name}, f'expected {rxn_name}, got {sorted(names)}'

    assert len(results) == expected_count, f'expected {expected_count}, got {len(results)}'

    for _, rxn in results:
        assert rxn == r, f'expected {r}, got {rxn}'
