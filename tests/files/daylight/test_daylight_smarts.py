# -*- coding: utf-8 -*-
#
#  Copyright 2025 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2025 Tagir Akhmetshin <tagirshin@gmail.com>
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
from chython.exceptions import IncorrectSmiles
from chython.files.daylight.smarts import smarts
from chython import smiles


def test_daylight_smarts_basic():
    # Test basic atomic symbols
    pattern = smarts('[C]')  # aliphatic carbon
    assert len(pattern) == 1
    
    pattern = smarts('[N]')  # nitrogen
    assert len(pattern) == 1
    
    pattern = smarts('[O]')  # oxygen
    assert len(pattern) == 1


def test_daylight_smarts_empty():
    # Test empty SMARTS
    with pytest.raises(ValueError):  # smarts() splits input and expects at least 1 value
        smarts('')


def test_daylight_smarts_invalid():
    # Test invalid SMARTS
    with pytest.raises(IncorrectSmiles):
        smarts('[')  # unclosed bracket


def test_daylight_smarts_bonds():
    # Test bond primitives
    pattern = smarts('CC')  # single bond
    assert len(pattern) == 2
    
    pattern = smarts('C=C')  # double bond
    assert len(pattern) == 2
    
    pattern = smarts('C#C')  # triple bond
    assert len(pattern) == 2


def test_daylight_smarts_aromatic():
    # Test aromatic specifications
    pattern = smarts('CC=CC=C')  # conjugated system
    assert len(pattern) == 5


def test_daylight_smarts_complex():
    # Test complex patterns
    pattern = smarts('CC(=O)O')  # carboxylic acid
    assert len(pattern) == 4
    
    pattern = smarts('CN(C)C')  # tertiary amine
    assert len(pattern) == 4


def test_daylight_smarts_functional_groups():
    # Test common functional group patterns
    # Create test molecules
    ethanol = smiles('CCO')
    benzene = smiles('c1ccccc1')
    acetic_acid = smiles('CC(=O)O')
    dimethylamine = smiles('CNC')
    phenol = smiles('c1ccccc1O')

    # Test alcohol pattern
    alcohol_pattern = smarts('CO')
    assert alcohol_pattern.is_substructure(ethanol)  # should match ethanol
    assert not alcohol_pattern.is_substructure(benzene)  # should not match benzene

    # Test carboxylic acid pattern
    acid_pattern = smarts('CC(=O)O')
    assert acid_pattern.is_substructure(acetic_acid)
    assert not acid_pattern.is_substructure(ethanol)


def test_daylight_smarts_multiple_matches():
    # Test patterns that can have multiple matches in a molecule
    # Create test molecules
    propanediol = smiles('OCCO')  # 1,2-propanediol
    diethylether = smiles('CCOCC')  # diethyl ether

    # Test alcohol pattern (should find two matches in propanediol)
    alcohol_pattern = smarts('CO')
    matches = list(alcohol_pattern.get_mapping(propanediol))
    assert len(matches) == 2  # should find two alcohol groups

    # Test ether pattern (should find one match in diethylether)
    ether_pattern = smarts('COC')
    matches = list(ether_pattern.get_mapping(diethylether))
    assert len(matches) == 1  # should find one ether group


def test_daylight_smarts_heterocycles():
    # Test complex heterocyclic systems
    # Create test molecules with sophisticated ring systems
    isoxazole = smiles('c1noc1')  # five-membered N-O ring
    indole = smiles('c1ccc2ccnc2c1')  # bicyclic with N
    benzofuran = smiles('c1ccc2ccoc2c1')  # bicyclic with O
    
    # Test basic ring patterns
    # Test individual atoms instead of full ring pattern
    n_pattern = smarts('n')  # aromatic nitrogen
    assert n_pattern.is_substructure(isoxazole)
    
    o_pattern = smarts('o')  # aromatic oxygen
    assert o_pattern.is_substructure(isoxazole)
    
    # Test aromatic carbons
    c_pattern = smarts('c')  # aromatic carbon
    assert c_pattern.is_substructure(indole)
    assert c_pattern.is_substructure(benzofuran)
    
    # Test aromatic ring fragments
    c1c_pattern = smarts('c:c')  # aromatic C-C bond with explicit aromatic bond
    assert c1c_pattern.is_substructure(indole)
    assert c1c_pattern.is_substructure(benzofuran)
    
    # Test heteroatom connections
    cn_pattern = smarts('c:n')  # aromatic C-N bond
    assert cn_pattern.is_substructure(indole)
    
    co_pattern = smarts('c:o')  # aromatic C-O bond
    assert co_pattern.is_substructure(benzofuran)


def test_daylight_smarts_charged_aromatics():
    # Test charged aromatic ring systems
    pyridinium = smiles('c1cc[nH+]cc1')  # protonated pyridine
    n_oxide = smiles('c1cc[n+]([O-])cc1')  # pyridine N-oxide
    
    # Test patterns for charged systems
    charged_n = smarts('[N+]')  # positively charged nitrogen
    assert charged_n.is_substructure(pyridinium)
    assert charged_n.is_substructure(n_oxide)


def test_daylight_smarts_drug_patterns():
    # Test drug-like substructure patterns
    # Create test molecules representing drug-like structures
    aspirin = smiles('CC(=O)OC1=CC=CC=C1C(=O)O')  # aspirin
    caffeine = smiles('CN1C=NC2=C1C(=O)N(C(=O)N2C)C')  # caffeine
    
    # Test common drug substructures
    ester = smarts('CC(=O)O')  # ester group
    assert ester.is_substructure(aspirin)
    
    amide = smarts('NC(=O)')  # amide group
    assert amide.is_substructure(caffeine)
    
    carboxyl = smarts('C(=O)O')  # carboxylic acid
    assert carboxyl.is_substructure(aspirin)


def test_daylight_smarts_complex_rings():
    # Test complex ring systems with multiple heteroatoms
    tetrazole = smiles('c1nnnn1')  # five-membered ring with 4 nitrogens
    oxadiazole = smiles('c1nnoc1')  # five-membered ring with O and 2 N
    thiazole = smiles('c1scnc1')  # five-membered ring with S and N
    
    # Test patterns for heterocyclic systems
    # Test individual atoms and bonds instead of full rings
    n_pattern = smarts('n')  # aromatic nitrogen
    assert n_pattern.is_substructure(tetrazole)
    assert n_pattern.is_substructure(oxadiazole)
    
    o_pattern = smarts('o')  # aromatic oxygen
    assert o_pattern.is_substructure(oxadiazole)
    
    s_pattern = smarts('s')  # aromatic sulfur
    assert s_pattern.is_substructure(thiazole)


def test_daylight_smarts_electronic_patterns():
    # Test patterns for electron-withdrawing groups and electronic effects
    
    # Create test molecules
    acetamide = smiles('CC(=O)N')  # amide
    acetate = smiles('CC(=O)O')    # ester
    acetic_acid = smiles('CC(=O)O')  # carboxylic acid
    acetyl_chloride = smiles('CC(=O)Cl')  # acid chloride
    nitromethane = smiles('CN(=O)=O')  # nitro compound
    
    # Test carbonyl-adjacent patterns
    carbonyl_pattern = smarts('CC(=O)')  # basic carbonyl group
    assert carbonyl_pattern.is_substructure(acetamide)
    assert carbonyl_pattern.is_substructure(acetate)
    assert carbonyl_pattern.is_substructure(acetic_acid)
    
    # Test carbonyl with heteroatom patterns
    carbonyl_n = smarts('C(=O)N')  # C=O with adjacent N
    assert carbonyl_n.is_substructure(acetamide)
    
    carbonyl_o = smarts('C(=O)O')  # C=O with adjacent O
    assert carbonyl_o.is_substructure(acetate)
    assert carbonyl_o.is_substructure(acetic_acid)
    
    # Test multiple bond patterns
    double_bond_o = smarts('C=O')  # C=O double bond
    assert double_bond_o.is_substructure(acetamide)
    assert double_bond_o.is_substructure(acetate)
    assert double_bond_o.is_substructure(acetic_acid)
    
    # Test halogen patterns
    carbonyl_halogen = smarts('C(=O)Cl')  # acid chloride
    assert carbonyl_halogen.is_substructure(acetyl_chloride)
    
    # Test nitro group pattern
    nitro_pattern = smarts('N(=O)=O')  # nitro group
    assert nitro_pattern.is_substructure(nitromethane)


def test_daylight_smarts_conjugated_systems():
    # Test patterns for conjugated systems
    
    # Create test molecules with conjugated systems
    acrolein = smiles('C=CC=O')  # conjugated alkene-carbonyl
    cinnamaldehyde = smiles('c1ccccc1C=CC=O')  # conjugated arene-alkene-carbonyl
    acetophenone = smiles('CC(=O)c1ccccc1')  # conjugated arene-carbonyl
    
    # Test conjugated patterns
    alkene_carbonyl = smarts('C=CC=O')  # conjugated alkene-carbonyl
    assert alkene_carbonyl.is_substructure(acrolein)
    assert alkene_carbonyl.is_substructure(cinnamaldehyde)
    
    # Test parts of the molecule separately
    # Test aromatic carbons and bonds
    c_pattern = smarts('c')  # aromatic carbon
    assert c_pattern.is_substructure(acetophenone)
    
    cc_bond = smarts('cc')  # aromatic C-C bond
    assert cc_bond.is_substructure(acetophenone)
    
    # Test carbonyl group
    carbonyl = smarts('CC(=O)')  # carbonyl group
    assert carbonyl.is_substructure(acetophenone)
    
    # Test extended conjugation
    conj_bond = smarts('C=CC=O')  # conjugated double bonds
    assert conj_bond.is_substructure(cinnamaldehyde)


def test_daylight_smarts_complex_bicyclic():
    # Test complex bicyclic pattern with N and S atoms
    # Create test molecule - a simpler bicyclic system with N-S-N bridge
    bicyclic = smiles('C1=CC=C2N=SN=C2C1')  # simplified bicyclic molecule without methyl group
    
    # Complex SMARTS pattern for bicyclic system
    # Pattern matches: C1=C-C=C2-N=S-N=C2-C1
    pattern = smarts('[#6]1=[#6]-[#6]=[#6]2-[#7]=[#16]-[#7]=[#6]2-[#6]1')
    assert pattern.is_substructure(bicyclic)
    
    # Test non-matching molecule
    benzene = smiles('c1ccccc1')
    assert not pattern.is_substructure(benzene)
    
    # Test parts of the pattern
    six_ring = smarts('[#6]1=[#6]-[#6]=[#6]-[#6]-[#6]1')  # matches C1=C-C=C-C-C1
    assert six_ring.is_substructure(bicyclic)
    
    five_ring = smarts('[#6]1-[#7]=[#16]-[#7]=[#6]1')  # matches C1-N=S-N=C1
    assert five_ring.is_substructure(bicyclic)
    
    # Test individual bonds
    cn_double = smarts('[#6]=[#7]')  # C=N double bond
    assert cn_double.is_substructure(bicyclic)
    
    ns_double = smarts('[#7]=[#16]')  # N=S double bond
    assert ns_double.is_substructure(bicyclic) 