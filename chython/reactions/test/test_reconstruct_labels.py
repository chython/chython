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
"""A `react:` label names one family of conditions.

Two records that need different conditions get different labels: Pd/Cu coupling of an unactivated aryl
electrophile, SNAr of an activated Cl/Br/I/sulfonate, SNAr of a fluoride, the carbonyl leaving group of an
amide bond, the hybridization of a Suzuki boron, the boron source of a borylation and the origin of a
Grignard reagent.
"""
from collections import Counter

import pytest

from .._tables import reaction_rules
from ...core import read_smiles as smiles


LABELS = [
    # an alkene meeting an aryl electrophile: hydroboration only where nothing conjugates it
    ('Brc1ccccc1.C=CCOc1ccccc1>>c1ccc(cc1)CCCOc1ccccc1', 'hydroboration_suzuki'),
    ('Brc1ccccc1.C=CC(=O)OC>>COC(=O)CCc1ccccc1', 'reductive_heck'),
    ('Brc1ccccc1.C=CC#N>>N#CCCc1ccccc1', 'reductive_heck'),
    # Pd/Cu: the unactivated aryl electrophile
    ('Brc1ccccc1.CCN>>CCNc1ccccc1', 'buchwald_hartwig'),
    ('Brc1cccnc1.CCN>>CCNc1cccnc1', 'buchwald_hartwig'),                    # meta to the ring N
    ('FC(F)(F)S(=O)(=O)Oc1ccccc1.CCN>>CCNc1ccccc1', 'buchwald_hartwig'),
    ('Ic1ccccc1.Oc1ccccc1>>c1ccc(Oc2ccccc2)cc1', 'ullmann_phenol'),
    ('Brc1ccccc1.c1cn[nH]c1>>c1cnn(c1)-c1ccccc1', 'ullmann_pyrrole'),
    ('Brc1ccccc1.CCO>>CCOc1ccccc1', 'ullmann_alcohol'),
    # SNAr: Cl/Br/I/sulfonate ortho or para to a ring N or to NO2
    ('Clc1ncccn1.CCN>>CCNc1ncccn1', 'snar_activated'),
    ('Clc1ncnc2ccccc12.CCN>>CCNc1ncnc2ccccc12', 'snar_activated'),
    ('Clc1ccncc1.CCN>>CCNc1ccncc1', 'snar_activated'),
    ('Clc1ccc(cc1)[N+](=O)[O-].CCN>>CCNc1ccc(cc1)[N+](=O)[O-]', 'snar_activated'),
    ('Clc1ccccc1[N+](=O)[O-].CCNCC>>CCN(CC)c1ccccc1[N+](=O)[O-]', 'snar_activated'),
    ('Clc1ncccn1.Oc1ccccc1>>c1ccc(Oc2ncccn2)cc1', 'snar_activated'),
    ('Clc1ncccn1.c1cn[nH]c1>>c1cnn(c1)-c1ncccn1', 'snar_activated'),
    ('Clc1ncccn1.CCO>>CCOc1ncccn1', 'snar_activated'),
    ('FC(F)(F)S(=O)(=O)Oc1ccccn1.CCN>>CCNc1ccccn1', 'snar_activated'),
    # SNAr: the fluoride, activated or not
    ('Fc1ccc(cc1)[N+](=O)[O-].CCN>>CCNc1ccc(cc1)[N+](=O)[O-]', 'snar'),
    ('Fc1ccccn1.CCN>>CCNc1ccccn1', 'snar'),
    ('Fc1ccccc1.c1cn[nH]c1>>c1cnn(c1)-c1ccccc1', 'snar'),
    # Suzuki: sp2 boron, sp3 boron, sp3 trifluoroborate
    ('Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1', 'suzuki'),
    ('Brc1ccccc1.CC1(C)OB(C=C)OC1(C)C>>C=Cc1ccccc1', 'suzuki'),
    ('Brc1ccccc1.[K+].F[B-](F)(F)c1ccccc1>>c1ccc(-c2ccccc2)cc1', 'suzuki'),
    ('Brc1ccccc1.OB(O)C1CC1>>c1ccc(cc1)C1CC1', 'suzuki_alkyl'),
    ('Brc1ccccc1.CCCCB1OC(C)(C)C(C)(C)O1>>CCCCc1ccccc1', 'suzuki_alkyl'),
    ('Brc1ccccc1.CB1OB(C)OB(C)O1>>Cc1ccccc1', 'suzuki_alkyl'),
    ('Brc1ccccc1.[K+].F[B-](F)(F)CCCC>>CCCCc1ccccc1', 'suzuki_alkyl_bf3k'),
    # borylation: the boron source is an agent, and the product form names the conditions
    ('Brc1ccccc1.CC1(C)OB(B2OC(C)(C)C(C)(C)O2)OC1(C)C>>CC1(C)OB(c2ccccc2)OC1(C)C', 'miyaura_borylation'),
    ('Brc1ccccc1>>CC1(C)OB(c2ccccc2)OC1(C)C', 'miyaura_borylation'),
    ('Brc1ccccc1.CC1(C)COB(B2OCC(C)(C)CO2)OC1>>CC1(C)COB(c2ccccc2)OC1', 'miyaura_borylation'),
    ('Brc1ccccc1.CCC1(CC)OB(B2OC(CC)(CC)C(CC)(CC)O2)OC1(CC)CC>>CCC1(CC)OB(c2ccccc2)OC1(CC)CC', 'miyaura_borylation'),
    ('Brc1ccccc1.CC(C)OB(OC(C)C)OC(C)C>>OB(O)c1ccccc1', 'lithiation_borylation'),
    ('Brc1ccccc1.CC1(C)OB(B2OC(C)(C)C(C)(C)O2)OC1(C)C>>OB(O)c1ccccc1', 'lithiation_borylation'),
    # amide bond: by what leaves the carbonyl
    ('CCC(=O)O.CCN>>CCNC(=O)CC', 'amidation'),
    ('CCC(=O)[O-].[Na+].CCN>>CCNC(=O)CC', 'amidation'),
    ('CCC(=O)Cl.CCN>>CCNC(=O)CC', 'acyl_halide_amidation'),
    ('CCC(=O)Br.CCN>>CCNC(=O)CC', 'acyl_halide_amidation'),
    ('CCC(=O)OC(=O)CC.CCN>>CCNC(=O)CC', 'anhydride_amidation'),
    ('O=C1CCC(=O)O1.CCN>>CCNC(=O)CCC(=O)O', 'anhydride_amidation'),
    ('CCC(=O)ON1C(=O)CCC1=O.CCN>>CCNC(=O)CC', 'activated_ester_amidation'),
    ('CCC(=O)Oc1c(F)c(F)c(F)c(F)c1F.CCN>>CCNC(=O)CC', 'activated_ester_amidation'),
    ('CCC(=O)Oc1ccc(cc1)[N+](=O)[O-].CCN>>CCNC(=O)CC', 'activated_ester_amidation'),
    ('CCC(=O)OC.CCN>>CCNC(=O)CC', 'aminolysis'),
    # organometallic addition: preformed, or made in the pot from the halide
    ('Br[Mg]c1ccccc1.O=Cc1ccccc1>>OC(c1ccccc1)c1ccccc1', 'grignard'),
    ('Brc1ccccc1.O=Cc1ccccc1>>OC(c1ccccc1)c1ccccc1', 'grignard_from_halide'),
    ('CCBr.CC(C)=O>>CCC(C)(C)O', 'grignard_from_halide'),
]


@pytest.mark.parametrize('record, label', LABELS, ids=[label for _, label in LABELS])
def test_a_record_gets_the_label_of_its_conditions(record, label):
    assert smiles(record).reconstruct_mapping() == ('react:' + label,)


@pytest.mark.parametrize('record', [
    'Brc1ccccc1.CC1(C)OB(B2OC(C)(C)C(C)(C)O2)OC1(C)C>>CC1(C)OB(c2ccccc2)OC1(C)C',
    'Brc1ccccc1.CC(C)OB(OC(C)C)OC(C)C>>OB(O)c1ccccc1',
])
def test_the_boron_source_is_an_agent(record):
    """No atom of the diboron or the borate is mapped: the row creates the boron it installs."""
    rxn = smiles(record)
    assert rxn.reconstruct_mapping()
    boron = rxn.reactants[1]
    assert all(not boron.map_number_of(n) for n in boron.atom_numbers)


METAL_CATALYSED = {'buchwald_hartwig', 'ullmann_phenol', 'ullmann_pyrrole', 'ullmann_alcohol'}


def test_an_aryl_fluoride_never_carries_a_metal_label():
    rows = [rule for family in reaction_rules().values() for rule in family]
    metal = [rule.id for rule in rows if 'aryl_fluoride' in rule.groups and rule.name in METAL_CATALYSED]
    assert not metal


def test_an_activated_electrophile_row_is_snar_and_an_unactivated_one_is_not():
    rows = [rule for family in reaction_rules().values() for rule in family]
    for rule in rows:
        if any(group.startswith('activated_') and group.split('_', 2)[1] in ('aza', 'azole', 'nitro')
               for group in rule.groups):
            assert rule.name == 'snar_activated', rule.id
        elif {'aryl_halide', 'aryl_sulfonate'} & set(rule.groups):
            assert rule.name != 'snar', rule.id


# An aryl halide is SNAr-activated beside or para to a six-ring aza N, beside or para to NO2, or between
# two ring heteroatoms of a five-ring, one of them an aza N.  Every ring system class, against ethylamine.
HETEROARYL = [
    # 6
    ('Clc1ccccn1', 'snar_activated'), ('Clc1ccncc1', 'snar_activated'), ('Brc1cccnc1', 'buchwald_hartwig'),
    ('Clc1ncccn1', 'snar_activated'), ('Clc1ccncn1', 'snar_activated'), ('Brc1cncnc1', 'buchwald_hartwig'),
    ('Clc1cnccn1', 'snar_activated'), ('Clc1cccnn1', 'snar_activated'), ('Clc1ncncn1', 'snar_activated'),
    # 5
    ('Clc1nccs1', 'snar_activated'), ('Clc1ncco1', 'snar_activated'), ('Cn1ccnc1Br', 'snar_activated'),
    ('Cn1nnnc1Cl', 'snar_activated'), ('Cn1cnc(Br)n1', 'snar_activated'), ('Clc1nncs1', 'snar_activated'),
    ('Clc1ncns1', 'snar_activated'),
    ('Brc1cnn(C)c1', 'buchwald_hartwig'), ('Brc1cccs1', 'buchwald_hartwig'), ('Brc1ccsc1', 'buchwald_hartwig'),
    ('Brc1cn(C)cn1', 'buchwald_hartwig'), ('Brc1ccn(C)n1', 'buchwald_hartwig'),
    # 5+6
    ('Clc1nc2ccccc2s1', 'snar_activated'), ('Clc1nc2ccccc2o1', 'snar_activated'),
    ('Cn1c(Cl)nc2ccccc21', 'snar_activated'), ('Clc1ncnc2[nH]cnc12', 'snar_activated'),
    ('Cn1c(Br)nc2cncnc21', 'snar_activated'), ('Clc1ncnc2[nH]ccc12', 'snar_activated'),
    ('Clc1ccc2nccn2n1', 'snar_activated'), ('Clc1ncnc2ccsc12', 'snar_activated'),
    ('Brc1ccc2[nH]ccc2c1', 'buchwald_hartwig'), ('Cn1nc(Br)c2ccccc21', 'buchwald_hartwig'),
    # 6+6
    ('Clc1ccc2ccccc2n1', 'snar_activated'), ('Clc1ccnc2ccccc12', 'snar_activated'),
    ('Clc1nccc2ccccc12', 'snar_activated'), ('Clc1cnc2ccccc2n1', 'snar_activated'),
    ('Clc1ncnc2ccccc12', 'snar_activated'), ('Clc1ccnc2cccnc12', 'snar_activated'),
    ('Brc1cnc2ccccc2c1', 'buchwald_hartwig'), ('Brc1ccc2ncccc2c1', 'buchwald_hartwig'),
]


@pytest.mark.parametrize('halide, label', HETEROARYL, ids=[h for h, _ in HETEROARYL])
def test_every_heteroaromatic_ring_class_is_classified(halide, label):
    for x in ('Cl', 'Br'):
        if halide.startswith(x):
            product = 'CCN' + halide[len(x):]
            break
    else:
        product = halide.replace('Cl', 'NCC').replace('Br', 'NCC')
    assert smiles(f'{halide}.CCN>>{product}').reconstruct_mapping() == ('react:' + label,)


# A protecting group a corpus row installs from a recorded reagent is that row's reaction: the same
# conditions as a substituted twin no protective row names.  Only where no row reads the reagent, or no
# reagent is recorded, is the record a protection.
INSTALLED = [
    ('BrCc1ccccc1.C1COCCN1>>C(c1ccccc1)N1CCOCC1', 'react:n_alkylation'),
    ('ClCc1ccccc1.C1COCCN1>>C(c1ccccc1)N1CCOCC1', 'react:n_alkylation'),
    ('O=Cc1ccccc1.C1CCNCC1>>C(c1ccccc1)N1CCCCC1', 'react:reductive_amination'),
    ('BrCc1ccccc1.Nc1ccccc1>>C(Nc1ccccc1)c1ccccc1', 'react:n_alkylation'),
    ('O=Cc1ccccc1.Nc1ccccc1>>C(Nc1ccccc1)c1ccccc1', 'react:reductive_amination'),
    ('BrCc1ccccc1.c1cn[nH]c1>>C(n1cccn1)c1ccccc1', 'react:azole_n_alkylation'),
    ('BrCc1ccccc1.Oc1ccccc1>>C(Oc1ccccc1)c1ccccc1', 'react:phenol_alkylation'),
    ('ClCc1ccccc1.Oc1ccccc1>>C(Oc1ccccc1)c1ccccc1', 'react:phenol_alkylation'),
    ('BrCc1ccccc1.OCC1CCCCC1>>C(OCC1CCCCC1)c1ccccc1', 'react:williamson'),
    ('BrCc1ccccc1.OC(=O)c1ccncc1>>O=C(OCc1ccccc1)c1ccncc1', 'react:acid_alkylation'),
    ('OC(=O)c1ccccc1.C1COCCN1>>O=C(N1CCOCC1)c1ccccc1', 'react:amidation'),
    ('ClC(=O)c1ccccc1.C1COCCN1>>O=C(N1CCOCC1)c1ccccc1', 'react:acyl_halide_amidation'),
    ('OC(=O)c1ccccc1.Nc1ccccc1>>O=C(Nc1ccccc1)c1ccccc1', 'react:amidation'),
    ('ClC(=O)c1ccccc1.Nc1ccccc1>>O=C(Nc1ccccc1)c1ccccc1', 'react:acyl_halide_amidation'),
    ('ClC(=O)c1ccccc1.OCC1CCCCC1>>O=C(OCC1CCCCC1)c1ccccc1', 'react:acylation'),
    ('CI.Oc1ccccc1>>COc1ccccc1', 'react:phenol_alkylation'),
    ('CI.OCC1CCCCC1>>COCC1CCCCC1', 'react:williamson'),
    ('CI.OC(=O)c1ccncc1>>COC(=O)c1ccncc1', 'react:acid_alkylation'),
    ('CO.OC(=O)c1ccncc1>>COC(=O)c1ccncc1', 'react:esterification'),
    ('CCI.Oc1ccccc1>>CCOc1ccccc1', 'react:phenol_alkylation'),
    ('CCO.OC(=O)c1ccncc1>>CCOC(=O)c1ccncc1', 'react:esterification'),
    ('CC(=O)O.C1COCCN1>>CC(=O)N1CCOCC1', 'react:amidation'),
    ('CC(=O)Cl.C1COCCN1>>CC(=O)N1CCOCC1', 'react:acyl_halide_amidation'),
    ('CC(=O)OC(C)=O.C1COCCN1>>CC(=O)N1CCOCC1', 'react:anhydride_amidation'),
    ('CC(=O)O.Nc1ccccc1>>CC(=O)Nc1ccccc1', 'react:amidation'),
    ('CC(=O)Cl.Nc1ccccc1>>CC(=O)Nc1ccccc1', 'react:acyl_halide_amidation'),
    ('CS(=O)(=O)Cl.OCC1CCCCC1>>CS(=O)(=O)OCC1CCCCC1', 'react:sulfonylation'),
    ('CS(=O)(=O)Cl.C1COCCN1>>CS(=O)(=O)N1CCOCC1', 'react:sulfonamide_formation'),
    ('Cc1ccc(S(=O)(=O)Cl)cc1.C1COCCN1>>Cc1ccc(S(=O)(=O)N2CCOCC2)cc1', 'react:sulfonamide_formation'),
    ('O=C(Cl)OCc1ccccc1.C1COCCN1>>O=C(OCc1ccccc1)N1CCOCC1', 'react:carbamoylation'),
    ('COc1ccc(CCl)cc1.OCC1CCCCC1>>COc1ccc(COCC2CCCCC2)cc1', 'react:williamson'),
    ('O=S(=O)(OS(=O)(=O)C(F)(F)F)C(F)(F)F.Oc1ccccc1>>O=S(=O)(Oc1ccccc1)C(F)(F)F', 'react:triflation'),
    ('O=C1CCCCC1.Clc1ccc(N(S(=O)(=O)C(F)(F)F)S(=O)(=O)C(F)(F)F)nc1>>O=S(=O)(OC1=CCCCC1)C(F)(F)F',
     'react:enol_triflation'),
    # an acidic N-H or O-H, a glycosyl donor and a trityl cation are their own conditions
    ('BrCc1ccccc1.O=C1CCCN1>>O=C1CCCN1Cc1ccccc1', 'react:amide_n_alkylation'),
    ('CI.CS(=O)(=O)NCc1ccccc1>>CN(Cc1ccccc1)S(C)(=O)=O', 'react:amide_n_alkylation'),
    ('CCOC(=O)C(F)(F)Br.Oc1ccccc1>>CCOC(=O)C(F)(F)Oc1ccccc1', 'react:difluoroalkylation'),
    ('FC(F)OS(=O)(=O)C(F)(F)F.Oc1ccccc1>>FC(F)Oc1ccccc1', 'react:difluoromethylation'),
    ('CS(=O)(=O)F.C1COCCN1>>CS(=O)(=O)N1CCOCC1', 'react:sulfonyl_fluoride_sulfonamide_formation'),
    ('CS(=O)(=O)OS(C)(=O)=O.OCC1CCCCC1>>CS(=O)(=O)OCC1CCCCC1', 'react:anhydride_sulfonylation'),
    ('C[C@@H]1CCC[C@H](Br)O1.OCc1ccccc1>>c1cc(ccc1)CO[C@@H]1CCC[C@@H](C)O1', 'react:glycosylation'),
    ('ClC(c1ccccc1)(c1ccccc1)c1ccccc1.OCC1CCCCC1>>C(OCC1CCCCC1)(c1ccccc1)(c1ccccc1)c1ccccc1', 'react:tritylation'),
    ('COc1ccc(cc1)C(Cl)(c1ccccc1)c1ccc(OC)cc1.OCC1CCCCC1>>COc1ccc(cc1)C(OCC1CCCCC1)(c1ccccc1)c1ccc(OC)cc1',
     'react:tritylation'),
    ('ClC(c1ccccc1)(c1ccccc1)c1ccccc1.c1nc[nH]c1>>c1ncn(c1)C(c1ccccc1)(c1ccccc1)c1ccccc1', 'react:tritylation'),
    # no row reads the reagent: the protection answers
    ('CC(C)(C)OC(=O)OC(=O)OC(C)(C)C.C1COCCN1>>CC(C)(C)OC(=O)N1CCOCC1', 'protect:amine_boc'),
    ('CC(C)(C)[Si](C)(C)Cl.OCC1CCCCC1>>CC(C)(C)[Si](C)(C)OCC1CCCCC1', 'protect:hydroxyl_tbs'),
    ('C[Si](C)(C)CCOCCl.c1cn[nH]c1>>C[Si](C)(C)CCOCn1cccn1', 'protect:azole_sem'),
    ('C1=COCCC1.OCC1CCCCC1>>C(OC1CCCCO1)C1CCCCC1', 'protect:hydroxyl_thp'),
    ('COS(=O)(=O)OC.Oc1ccccc1>>COc1ccccc1', 'protect:phenol_methyl'),
    # no reagent recorded: the protection answers
    ('Oc1ccccc1>>O=S(=O)(Oc1ccccc1)C(F)(F)F', 'protect:phenol_triflyl'),          # the row's agent, unrecorded
    ('O=C1CCCCC1>>O=S(=O)(OC1=CCCCC1)C(F)(F)F', 'protect:enol_triflyl'),
    ('C1COCCN1>>C(c1ccccc1)N1CCOCC1', 'protect:amine_benzyl'),
    ('Oc1ccccc1>>COc1ccccc1', 'protect:phenol_methyl'),
    ('C1COCCN1>>CC(=O)N1CCOCC1', 'protect:amine_acyl'),
    ('OCC1CCCCC1>>C(OCC1CCCCC1)(c1ccccc1)(c1ccccc1)c1ccccc1', 'protect:hydroxyl_tritil'),
    # no protecting group: the controls
    ('CI.C1CCNCC1>>CN1CCCCC1', 'react:n_alkylation'),
    ('CI.c1cn[nH]c1>>Cn1cccn1', 'react:azole_n_alkylation'),
    ('C=O.C1CCNCC1>>CN1CCCCC1', 'react:reductive_amination'),
    ('Brc1ccccc1>>CC1(C)OB(c2ccccc2)OC1(C)C', 'react:miyaura_borylation'),       # an agent no group names
]


@pytest.mark.parametrize('record, label', INSTALLED, ids=[record for record, _ in INSTALLED])
def test_a_group_installed_from_a_recorded_reagent_is_the_reaction(record, label):
    assert smiles(record).reconstruct_mapping() == (label,)


def _centre(record):
    """The label, and the reaction centre the mapping implies: the changed bonds, and the before/after state
    of their atoms and of the first shell around them, read from `transition_view`."""
    rxn = smiles(record)
    label = rxn.reconstruct_mapping()
    view = rxn.transition_view()
    elements = view.elements.tolist()
    states = list(zip(elements, view.h_before.tolist(), view.n_before.tolist(), view.h_after.tolist(),
                      view.n_after.tolist()))
    changed = [(a, b, x, y) for (a, b), x, y in zip(view.bonds.tolist(), view.bond_before.tolist(),
                                                    view.bond_after.tolist()) if x != y]
    core = {i for a, b, _, _ in changed for i in (a, b)}
    shell = {j for i in core for j in range(len(elements)) if view.distances[i, j] == 1} - core
    return label, (Counter((tuple(sorted((elements[a], elements[b]))), x, y) for a, b, x, y in changed),
                   Counter(states[i] for i in core), Counter(states[i] for i in shell))


@pytest.mark.parametrize('record, twin', [
    ('BrCc1ccccc1.C1COCCN1>>C(c1ccccc1)N1CCOCC1', 'BrCc1ccc(F)cc1.C1COCCN1>>Fc1ccc(CN2CCOCC2)cc1'),
    ('O=Cc1ccccc1.C1CCNCC1>>C(c1ccccc1)N1CCCCC1', 'O=Cc1ccc(F)cc1.C1CCNCC1>>Fc1ccc(CN2CCCCC2)cc1'),
    ('OC(=O)c1ccccc1.NCc1ccccc1>>O=C(NCc1ccccc1)c1ccccc1',
     'OC(=O)c1ccc(F)cc1.NCc1ccccc1>>O=C(NCc1ccccc1)c1ccc(F)cc1'),
    ('ClC(=O)c1ccccc1.NCc1ccccc1>>O=C(NCc1ccccc1)c1ccccc1',
     'ClC(=O)c1ccc(F)cc1.NCc1ccccc1>>O=C(NCc1ccccc1)c1ccc(F)cc1'),
    ('BrCc1ccccc1.Oc1ccccc1>>C(Oc1ccccc1)c1ccccc1', 'Oc1ccccc1.BrCc1ccc(F)cc1>>Fc1ccc(COc2ccccc2)cc1'),
], ids=['benzyl_bromide', 'benzaldehyde', 'benzoic_acid', 'benzoyl_chloride', 'benzyl_ether'])
def test_a_group_installed_by_a_row_maps_as_its_substituted_twin(record, twin):
    label, centre = _centre(record)
    assert label[0].startswith('react:')
    assert _centre(twin) == (label, centre)


# Every electrophile class alkylates either nitrogen of an unsymmetrical azole, under one label: the `pyrrole`
# row takes the N drawn with the H, the `pyrazole` or `imidazole` row the free N.
AZOLE_ELECTROPHILES = [
    ('BrCc2ccccc2', 'Cc2ccccc2', 'azole_n_alkylation'),
    ('COS(C)(=O)=O', 'C', 'azole_n_alkylation'),
    ('ClCOC(=O)C(C)(C)C', 'COC(=O)C(C)(C)C', 'azole_n_alkylation'),
    ('CCOP(=O)(COS(=O)(=O)c2ccc(C)cc2)OCC', 'CP(=O)(OCC)OCC', 'azole_n_alkylation'),
    ('ClC[Si](C)(C)C', 'C[Si](C)(C)C', 'azole_n_alkylation'),
    ('ClCSC', 'CSC', 'azole_n_alkylation'),
    ('ClCOc2ccccc2', 'COc2ccccc2', 'azole_n_alkylation'),
    ('FC(F)OS(=O)(=O)C(F)(F)F', 'C(F)F', 'difluoromethylation'),
    ('CCOC(=O)C(F)(F)Br', 'C(F)(F)C(=O)OCC', 'difluoroalkylation'),
    ('ClC(c2ccccc2)(c2ccccc2)c2ccccc2', 'C(c2ccccc2)(c2ccccc2)c2ccccc2', 'tritylation'),
]


@pytest.mark.parametrize('reagent, group, label', AZOLE_ELECTROPHILES,
                         ids=[x[2] + ':' + x[1] for x in AZOLE_ELECTROPHILES])
@pytest.mark.parametrize('azole, products', [
    ('CC1=NNC=C1', ('Cc1ccn({})n1', 'Cc1ccnn1{}')),           # 3-methylpyrazole: the 1,3- and the 1,5-isomer
    ('Cc1c[nH]cn1', ('Cc1cn({})cn1', 'Cc1cncn1{}')),         # 4-methylimidazole: the 1,4- and the 1,5-isomer
], ids=['pyrazole', 'imidazole'])
def test_either_azole_nitrogen_is_one_label(reagent, group, label, azole, products):
    # both input orders: a template's working union renumbers every input after the first
    for product in products:
        for inputs in (f'{azole}.{reagent}', f'{reagent}.{azole}'):
            assert smiles(f'{inputs}>>{product.format(group)}').reconstruct_mapping() == ('react:' + label,), inputs


def test_react_reaches_either_nitrogen_of_a_mobile_hydrogen():
    """One row, every form: `react()` builds both N-benzyl regioisomers of 3-methylpyrazole and both
    N-aryl products of an N-methylguanidine, whichever form the input was stored in."""
    def built(text, partner, reaction):
        molecule, other = smiles(text), smiles(partner)
        molecule.canonicalize()
        other.canonicalize()
        out = set()
        for outcome in molecule.react(other, reaction=reaction):
            product = outcome.reaction.products[0]
            product.canonicalize()
            out.add(product.canonical_bytes)
        return out

    def canonical(*texts):
        out = set()
        for text in texts:
            molecule = smiles(text)
            molecule.canonicalize()
            out.add(molecule.canonical_bytes)
        return out

    assert built('CC1=NNC=C1', 'BrCc1ccccc1', 'azole_n_alkylation') == canonical(
        'Cc1ccn(Cc2ccccc2)n1', 'Cc1ccnn1Cc1ccccc1')
    assert canonical('CNC(=N)Nc1ccccc1', 'CN=C(N)Nc1ccccc1') <= built('NC(=N)NC', 'Brc1ccccc1', 'buchwald_hartwig')
