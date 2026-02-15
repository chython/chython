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


# ===== New tests for 5 SMARTS features =====


def test_any_atom_wildcard():
    """Test [*] and bare * parse and match any atom."""
    # [*] in brackets
    pattern_bracket = smarts('[*]')
    assert len(pattern_bracket) == 1

    # bare * outside brackets
    pattern_bare = smarts('*')
    assert len(pattern_bare) == 1

    # Match any atom in ethanol
    ethanol = smiles('CCO')
    matches_bracket = list(pattern_bracket.get_mapping(ethanol))
    assert len(matches_bracket) == 3  # C, C, O

    matches_bare = list(pattern_bare.get_mapping(ethanol))
    assert len(matches_bare) == 3

    # * in context: *-O pattern
    star_o = smarts('*O')
    assert star_o.is_substructure(ethanol)


def test_any_aromatic_atom():
    """Test [a] matches any aromatic atom, doesn't match aliphatic."""
    pattern = smarts('[a]')
    assert len(pattern) == 1

    # Should match aromatic C in benzene
    benzene = smiles('c1ccccc1')
    matches = list(pattern.get_mapping(benzene))
    assert len(matches) == 6  # all 6 aromatic carbons

    # Should not match aliphatic C in ethane
    ethane = smiles('CC')
    matches = list(pattern.get_mapping(ethane))
    assert len(matches) == 0

    # Should match aromatic N in pyridine
    pyridine = smiles('c1ccncc1')
    matches = list(pattern.get_mapping(pyridine))
    assert len(matches) == 6  # 5 aromatic C + 1 aromatic N


def test_total_connectivity():
    """Test [CX3] matches sp2 carbon (3 total connections), [CX4] matches sp3."""
    # X4 carbon: 4 total connections (neighbors + implicit H)
    cx4 = smarts('[C;X4]')
    ethanol = smiles('CCO')
    # First C has 1 neighbor (C) + 3 implicit H = 4 total
    # Second C has 2 neighbors (C, O) + 2 implicit H = 4 total
    matches = list(cx4.get_mapping(ethanol))
    assert len(matches) == 2

    # X3 carbon: 3 total connections
    cx3 = smarts('[C;X3]')
    # Formaldehyde: C=O, C has 1 neighbor (O) + ... actually let's use acetic acid
    # In acetaldehyde CH3CHO: the carbonyl C has 2 neighbors (C, O via double bond) + 1 implicit H = 3
    acetaldehyde = smiles('CC=O')
    matches = list(cx3.get_mapping(acetaldehyde))
    assert len(matches) == 1  # only the carbonyl C

    # X4 should not match the carbonyl C
    cx4_matches = list(cx4.get_mapping(acetaldehyde))
    assert len(cx4_matches) == 1  # only the methyl C


def test_ring_membership():
    """Test [R], [R0], [R1], [R2] ring membership count."""
    # [R] matches any atom in a ring
    r_pattern = smarts('[C;R]')
    cyclohexane = smiles('C1CCCCC1')
    matches = list(r_pattern.get_mapping(cyclohexane))
    assert len(matches) == 6  # all 6 carbons in ring

    # [R] should not match non-ring atom
    ethane = smiles('CC')
    matches = list(r_pattern.get_mapping(ethane))
    assert len(matches) == 0

    # [R0] matches non-ring atom
    r0_pattern = smarts('[C;R0]')
    matches = list(r0_pattern.get_mapping(ethane))
    assert len(matches) == 2  # both carbons not in ring

    matches = list(r0_pattern.get_mapping(cyclohexane))
    assert len(matches) == 0  # all in ring

    # [R2] matches bridgehead atoms (in 2 rings) - naphthalene
    r2_pattern = smarts('[#6;R2]')
    naphthalene = smiles('c1ccc2ccccc2c1')
    matches = list(r2_pattern.get_mapping(naphthalene))
    assert len(matches) == 2  # two bridgehead carbons

    # [R1] matches atoms in exactly 1 ring
    r1_pattern = smarts('[#6;R1]')
    matches_r1 = list(r1_pattern.get_mapping(naphthalene))
    assert len(matches_r1) == 8  # 8 non-bridgehead carbons


def test_negated_charge():
    """Test [N;!+] matches neutral N, doesn't match protonated N."""
    # [N;!+] - nitrogen that is not positively charged
    n_not_pos = smarts('[N;!+]')

    # Should match neutral nitrogen
    methylamine = smiles('CN')
    assert n_not_pos.is_substructure(methylamine)

    # Should NOT match protonated nitrogen
    ammonium = smiles('C[NH3+]')
    matches = list(n_not_pos.get_mapping(ammonium))
    assert len(matches) == 0

    # [O;!-] - oxygen that is not negatively charged
    o_not_neg = smarts('[O;!-]')

    # Should match neutral oxygen
    ethanol = smiles('CCO')
    assert o_not_neg.is_substructure(ethanol)

    # Should NOT match deprotonated oxygen
    alkoxide = smiles('CC[O-]')
    matches = list(o_not_neg.get_mapping(alkoxide))
    assert len(matches) == 0

    # [N;!+;!-] - both negated: must be neutral (charge == 0)
    n_neutral = smarts('[N;!+;!-]')
    assert n_neutral.is_substructure(methylamine)
    assert not n_neutral.is_substructure(ammonium)


def test_combined_features():
    """Test combining multiple new features in one pattern."""
    # [N;X3;!+;D1] - nitrogen with 3 total connections, not positive, 1 heavy neighbor
    pattern = smarts('[N;X3;!+;D1]')
    methylamine = smiles('CN')
    # N in methylamine: 1 neighbor (C) + 2 implicit H = 3 total, not charged, D1
    assert pattern.is_substructure(methylamine)

    # Ammonium should NOT match (!+)
    ammonium = smiles('C[NH3+]')
    assert not pattern.is_substructure(ammonium)


# ===== Recursive SMARTS $() tests =====


def test_recursive_smarts_basic():
    """Test [N;$(NC=O)] matches amide N but not amine N."""
    pattern = smarts('[N;$(NC=O)]')
    assert len(pattern) == 1

    # Amide N: should match
    acetamide = smiles('CC(=O)N')
    assert pattern.is_substructure(acetamide)

    # Amine N: should NOT match
    methylamine = smiles('CN')
    assert not pattern.is_substructure(methylamine)


def test_recursive_smarts_negated():
    """Test [N;!$(NC=O)] matches amine N but not amide N."""
    pattern = smarts('[N;!$(NC=O)]')

    # Amine N: should match
    methylamine = smiles('CN')
    assert pattern.is_substructure(methylamine)

    # Amide N: should NOT match
    acetamide = smiles('CC(=O)N')
    assert not pattern.is_substructure(acetamide)


def test_recursive_smarts_multiple():
    """Test multiple recursive constraints AND-ed together."""
    # N connected to C AND not double-bonded to anything
    pattern = smarts('[N;$(NC);!$(N=*)]')

    # Methylamine: N-C, no N=X → match
    methylamine = smiles('CN')
    assert pattern.is_substructure(methylamine)

    # Imine: N=C → should NOT match (excluded by !$(N=*))
    imine = smiles('C=NC')
    assert not pattern.is_substructure(imine)


def test_recursive_smarts_with_primitives():
    """Test mixing $() with traditional SMARTS primitives."""
    # N with 3 total connections, not positive, connected to C
    pattern = smarts('[N;X3;!+;$(NC)]')

    # Primary amine: N-C, X3 (1 neighbor + 2H), neutral → match
    methylamine = smiles('CN')
    assert pattern.is_substructure(methylamine)

    # Ammonium: N+, even though X4 → should NOT match (!+)
    ammonium = smiles('C[NH3+]')
    assert not pattern.is_substructure(ammonium)


def test_recursive_smarts_nested_brackets():
    """Test recursive SMARTS with brackets inside $()."""
    # N connected to [C;X4] (sp3 carbon)
    pattern = smarts('[N;$(N[C;X4])]')

    # Methylamine: N-CH3 (X4 carbon) → match
    methylamine = smiles('CN')
    assert pattern.is_substructure(methylamine)

    # N connected to carbonyl C (X3) → should NOT match
    acetamide = smiles('CC(=O)N')
    # The amide N is connected to carbonyl C (X3), not X4
    # But also connected to H, so check carefully
    matches = list(pattern.get_mapping(acetamide))
    assert len(matches) == 0


def test_recursive_smarts_pure_recursive():
    """Test patterns that are purely recursive with no element specification."""
    # [$(NC=O)] should default to any element but constrained by recursive
    pattern = smarts('[$(NC=O)]')
    acetamide = smiles('CC(=O)N')
    assert pattern.is_substructure(acetamide)


def test_recursive_smarts_amino_acid():
    """Test recursive SMARTS for alpha-amino acid pattern."""
    # [N;D1;$(N[C;X4]C(=O)[O;D1])] — primary amine on sp3 C adjacent to carboxylic acid
    pattern = smarts('[N;D1;$(N[C;X4]C(=O)[O;D1])]')

    # Glycine (simplest amino acid): H2N-CH2-COOH
    glycine = smiles('NCC(=O)O')
    assert pattern.is_substructure(glycine)

    # Simple amine (no alpha-amino acid pattern) → should NOT match
    methylamine = smiles('CN')
    assert not pattern.is_substructure(methylamine)


def test_recursive_smarts_halogen():
    """Test recursive SMARTS for sp3-halogen pattern."""
    # [Cl;$([Cl][C;X4])] — chlorine on sp3 carbon only
    pattern = smarts('[Cl;$([Cl][C;X4])]')

    # Chloromethane: Cl-CH3 (sp3 C) → match
    chloromethane = smiles('CCl')
    assert pattern.is_substructure(chloromethane)

    # Acyl chloride: Cl-C(=O) (sp2 C, X3) → should NOT match
    acyl_chloride = smiles('CC(=O)Cl')
    assert not pattern.is_substructure(acyl_chloride)


def test_recursive_smarts_element_constraint():
    """Test element constraints in primitive positions (e.g., [*;O,S,P,N])."""
    # [*;O,S,P,N] should match only O, S, P, N — not C
    pattern = smarts('[*;O,S,P,N]')
    ethanol = smiles('CCO')
    matches = list(pattern.get_mapping(ethanol))
    assert len(matches) == 1  # only the O

    methane = smiles('C')
    assert not pattern.is_substructure(methane)

    # [*;!#6] — not carbon
    not_carbon = smarts('[*;!#6]')
    matches = list(not_carbon.get_mapping(ethanol))
    assert len(matches) == 1  # only the O

    # Semantic: [N;!$(N[*]=[*;O,S,P,N])] — N not adjacent to heteroatom double bond
    pattern2 = smarts('[N;!$(N[*]=[*;O,S,P,N])]')

    # Vinyl amine: N-C=C → C not in {O,S,P,N} → should match
    vinyl_amine = smiles('NC=C')
    assert pattern2.is_substructure(vinyl_amine)

    # Amide: N-C=O → O in {O,S,P,N} → should NOT match
    acetamide = smiles('NC(=O)C')
    assert not pattern2.is_substructure(acetamide)


def test_reaction_smarts_basic():
    """Parse reaction SMARTS and verify structure."""
    from chython.containers import ReactionContainer

    rxn = smarts('[C;D2:1]-[Br;D1:2]>>[C;D2:1]-[O;D1:3]')
    assert isinstance(rxn, ReactionContainer)
    assert len(rxn.reactants) == 1
    assert len(rxn.products) == 1


def test_reaction_smarts_multi_component():
    """Parse reaction SMARTS with multiple reactants/products."""
    from chython.containers import ReactionContainer

    rxn = smarts('[C;D3:1]-[C;D3:2](=[O;D1:3])-[Cl;D1:4]>>[C;D3:1]-[C;D3:2](=[O;D1:3])-[N;D2:5].[Cl;D0:4]')
    assert isinstance(rxn, ReactionContainer)
    assert len(rxn.reactants) == 1
    assert len(rxn.products) == 2


def test_reaction_smarts_roundtrip():
    """Parse reaction SMARTS, format back, parse again — strings must match."""
    rxn1 = smarts('[C;D2:1]-[Br;D1:2]>>[C;D2:1]-[O;D1:3]')
    s1 = str(rxn1)
    rxn2 = smarts(s1)
    s2 = str(rxn2)
    assert s1 == s2


def test_reaction_smarts_equality():
    """Two parses of the same reaction SMARTS are equal."""
    rxn1 = smarts('[C;D2:1]-[Br;D1:2]>>[C;D2:1]-[O;D1:3]')
    rxn2 = smarts('[C;D2:1]-[Br;D1:2]>>[C;D2:1]-[O;D1:3]')
    assert rxn1 == rxn2


def test_reaction_smarts_hashing():
    """Identical reaction SMARTS produce same hash; can be used in sets."""
    rxn1 = smarts('[C;D2:1]-[Br;D1:2]>>[C;D2:1]-[O;D1:3]')
    rxn2 = smarts('[C;D2:1]-[Br;D1:2]>>[C;D2:1]-[O;D1:3]')
    assert hash(rxn1) == hash(rxn2)
    assert len({rxn1, rxn2}) == 1


def test_reaction_smarts_inequality():
    """Different reaction SMARTS are not equal."""
    rxn1 = smarts('[C;D2:1]-[Br;D1:2]>>[C;D2:1]-[O;D1:3]')
    rxn2 = smarts('[C;D2:1]-[Cl;D1:2]>>[C;D2:1]-[O;D1:3]')
    assert rxn1 != rxn2


def test_reaction_smarts_to_reactor():
    """Convert reaction SMARTS to Reactor and apply it."""
    from chython.reactor import Reactor

    rxn = smarts('[C;D2:1]-[Br;D1:2]>>[C;D2:1]-[O;D1:3]')
    reactor = rxn.to_reactor()
    assert isinstance(reactor, Reactor)


def test_reaction_smarts_compose():
    """compose() on QueryContainer reaction returns QueryCGRContainer."""
    from chython.containers import QueryCGRContainer

    rxn = smarts('[C;D2:1]-[Br;D1:2]>>[C;D2:1]-[O;D1:3]')
    cgr = rxn.compose()
    assert isinstance(cgr, QueryCGRContainer)


def test_molecule_smarts_regression():
    """Molecule SMARTS (no >>) still returns QueryContainer."""
    from chython.containers import QueryContainer

    q = smarts('[C;D2:1]-[Br;D1:2]')
    assert isinstance(q, QueryContainer)
