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
import pytest

from .._tables import reaction_rules
from ...core import read_smiles as smiles


LABELS = [
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
    # borylation: the diboron with Pd, the borate after lithiation
    ('Brc1ccccc1.CC1(C)OB(B2OC(C)(C)C(C)(C)O2)OC1(C)C>>CC1(C)OB(c2ccccc2)OC1(C)C', 'miyaura_borylation'),
    ('Brc1ccccc1.CC(C)OB(OC(C)C)OC(C)C>>OB(O)c1ccccc1', 'lithiation_borylation'),
    ('Brc1ccccc1.CC(C)OB1OC(C)(C)C(C)(C)O1>>CC1(C)OB(c2ccccc2)OC1(C)C', 'lithiation_borylation'),
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
    'Brc1ccccc1.CC1(C)OB(B2OC(C)(C)C(C)(C)O2)OC1(C)C>>OB(O)c1ccccc1',
    'Brc1ccccc1>>OB(O)c1ccccc1',
])
def test_no_borylation_row_supplies_its_own_boron(record):
    """A diboron in the record is not read as a borate quench, and a record naming no boron source is not
    read as either: the boron source is a slot of every borylation row."""
    assert smiles(record).reconstruct_mapping() == ()


METAL_CATALYSED = {'buchwald_hartwig', 'ullmann_phenol', 'ullmann_pyrrole', 'ullmann_alcohol'}


def test_an_aryl_fluoride_never_carries_a_metal_label():
    rows = [rule for family in reaction_rules().values() for rule in family]
    metal = [rule.id for rule in rows if 'aryl_fluoride' in rule.groups and rule.name in METAL_CATALYSED]
    assert not metal


def test_an_activated_electrophile_row_is_snar_and_an_unactivated_one_is_not():
    rows = [rule for family in reaction_rules().values() for rule in family]
    for rule in rows:
        if any(group.startswith('activated_') and group.split('_', 2)[1] in ('aza', 'nitro')
               for group in rule.groups):
            assert rule.name == 'snar_activated', rule.id
        elif {'aryl_halide', 'aryl_sulfonate'} & set(rule.groups):
            assert rule.name != 'snar', rule.id
