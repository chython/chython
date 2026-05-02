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
     'suzuki', 'Brc1ccc(cc1)F.OB(O)c1ccccc1>>c1cc(ccc1-c1ccccc1)F', 1),
    ('buchwald_hartwig_aryl_bromide_primary_amine',
     'buchwald_hartwig', 'Brc1ccc(cc1)F.NCC>>c1cc(ccc1F)NCC', 1),
    ('mitsunobu_ether_primary_alcohol_phenol',
     'mitsunobu', 'OCC.Oc1ccc(cc1)F>>c1cc(ccc1F)OCC', 1),
    ('mitsunobu_ester_primary_alcohol_acid',
     'mitsunobu', 'OCC.OC(=O)c1ccccc1>>c1ccccc1C(OCC)=O', 1),
    ('deoxygenative_coupling_alcohol_aryl_bromide',
     'deoxygenative_coupling', 'OCC.Brc1ccc(cc1)F>>c1cc(ccc1F)CC', 1),
    ('decarboxylative_coupling_acid_aryl_bromide',
     'decarboxylative_coupling', 'OC(=O)CC.Brc1ccc(cc1)F>>c1cc(ccc1F)CC', 1),
    ('xec_aryl_bromide_alkyl_bromide',
     'xec', 'Brc1ccc(cc1)F.BrCCC>>c1cc(ccc1CCC)F', 1),
    ('reductive_amination_aldehyde_primary_amine',
     'reductive_amination', 'O=CC.NCc1ccccc1>>c1cc(ccc1)CNCC', 1),
    ('reductive_amination_ketone_secondary_amine',
     'reductive_amination', 'O=C(C)C.N(C)CC>>CC(C)N(C)CC', 1),
    ('ullmann_phenol_aryl_bromide_cresol',
     'ullmann_phenol', 'Brc1ccc(cc1)F.Oc1ccc(cc1)C>>c1cc(ccc1F)Oc1ccc(cc1)C', 1),
    ('ullmann_pyrrole_aryl_bromide_pyrrole',
     'ullmann_pyrrole', 'Brc1ccc(cc1)F.[nH]1cccc1>>c1ccn(c1)-c1ccc(cc1)F', 1),
    ('ullmann_pyrrole_aryl_bromide_pyrazole',
     'ullmann_pyrrole', 'Brc1ccc(cc1)F.[nH]1nccc1>>c1cc(ccc1-n1nccc1)F', 2),
    ('ullmann_pyrrole_aryl_bromide_imidazole',
     'ullmann_pyrrole', 'Brc1ccc(cc1)F.[nH]1cncc1>>c1cc(ccc1-n1ccnc1)F', 2),
    ('chan_lam_boronic_acid_primary_amine',
     'chan_lam', 'OB(O)c1ccc(cc1)F.NCC>>c1cc(ccc1F)NCC', 1),
    ('chan_lam_boronic_acid_phenol',
     'chan_lam', 'OB(O)c1ccc(cc1)F.Oc1ccc(cc1)C>>c1cc(ccc1F)Oc1ccc(cc1)C', 1),
    ('sulfonylation_sulfonyl_chloride_alcohol',
     'sulfonylation', 'ClS(=O)(=O)c1ccccc1.OCCC>>c1ccccc1S(=O)(=O)OCCC', 1),
    ('sulfonamide_formation_sulfonyl_chloride_amine',
     'sulfonamide_formation', 'ClS(=O)(=O)c1ccccc1.NCC>>c1ccccc1S(=O)(=O)NCC', 1),
    ('aminolysis_ester_primary_amine',
     'aminolysis', 'COC(=O)c1ccc(cc1)F.NCC>>c1cc(ccc1C(NCC)=O)F', 1),
    ('grignard_alkyl_bromide_aldehyde',
     'grignard', 'BrCC.O=Cc1ccc(cc1)F>>c1cc(ccc1F)C(CC)O', 1),
    ('grignard_alkyl_bromide_ketone',
     'grignard', 'BrCC.O=C(C)c1ccccc1>>c1cc(ccc1)C(O)(C)CC', 1),
    ('sonogashira_aryl_bromide_terminal_alkyne',
     'sonogashira', 'Brc1ccc(cc1)F.C#CC>>c1cc(ccc1F)C#CC', 1),
    ('urea_synthesis_isocyanate_primary_amine',
     'urea_synthesis', 'O=C=Nc1ccccc1.NCC>>c1ccccc1NC(NCC)=O', 1),
    ('cuaac_azide_terminal_alkyne',
     'cuaac', 'CC#C.CN=[N+]=[N-]>>Cc1cn(C)nn1', 1),
]

_three_component = [
    ('ugi_3cr_aldehyde_amine_isocyano',
     'ugi_3cr', 'O=CC.NCC.[C-]#[N+]C>>CCNC(C(=O)NC)C', 1),
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
