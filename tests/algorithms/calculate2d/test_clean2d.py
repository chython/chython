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
from chython import smiles, ReactionContainer
import pytest


def test_molecule_clean2d():
    # Test basic molecule coordinate calculation
    mol = smiles('CCO')
    mol.clean2d()
    
    # Check that coordinates exist for all atoms
    assert all(n in mol._plane for n in mol)
    
    # Check that coordinates are reasonable (finite numbers)
    for x, y in mol._plane.values():
        assert isinstance(x, float)
        assert isinstance(y, float)
        assert -100 < x < 100  # reasonable bounds
        assert -100 < y < 100


def test_complex_molecule_clean2d():
    # Test with more complex molecule
    mol = smiles('c1ccccc1CC(=O)O')
    mol.clean2d()
    
    # Check coordinates exist
    assert all(n in mol._plane for n in mol)
    
    # Verify ring atoms are roughly coplanar
    ring_atoms = [n for n in mol if len(mol._bonds[n]) == 2]
    if ring_atoms:
        coords = [mol._plane[n] for n in ring_atoms]
        # Calculate variance in y coordinates - should be small for planar ring
        y_coords = [y for _, y in coords]
        y_mean = sum(y_coords) / len(y_coords)
        y_variance = sum((y - y_mean) ** 2 for y in y_coords) / len(y_coords)
        assert y_variance < 1.0  # reasonable threshold for planarity


def test_disconnected_components():
    # Test molecule with multiple disconnected components
    mol = smiles('CCO.c1ccccc1')
    mol.clean2d()
    
    # Check all atoms have coordinates
    assert all(n in mol._plane for n in mol)
    
    # Components should be separated in space
    components = list(mol.connected_components)
    assert len(components) == 2
    
    # Get bounding boxes for each component
    def get_bounds(atoms):
        xs = [mol._plane[n][0] for n in atoms]
        ys = [mol._plane[n][1] for n in atoms]
        return min(xs), max(xs), min(ys), max(ys)
    
    bounds1 = get_bounds(components[0])
    bounds2 = get_bounds(components[1])
    
    # Check components don't overlap in x-direction
    assert bounds1[1] < bounds2[0] or bounds2[1] < bounds1[0]


def test_reaction_clean2d():
    # Create a simple reaction
    reactant = smiles('CCO')
    product = smiles('CC=O')
    reaction = ReactionContainer([reactant], [product])
    
    # Clean coordinates
    reaction.clean2d()
    
    # Check that all molecules have coordinates
    for molecule in reaction.molecules():
        assert all(n in molecule._plane for n in molecule)
    
    # Check that reactants are positioned before products
    reactant_max_x = max(x for mol in reaction.reactants 
                        for x, _ in mol._plane.values())
    product_min_x = min(x for mol in reaction.products 
                       for x, _ in mol._plane.values())
    assert reactant_max_x < product_min_x
    
    # Check arrow exists and is positioned between reactants and products
    assert hasattr(reaction, '_arrow')
    arrow_start, arrow_end = reaction._arrow
    assert reactant_max_x < arrow_start < arrow_end < product_min_x


def test_reaction_with_reagents():
    # Create reaction with reagents
    reactant = smiles('CCO')
    reagent = smiles('Cl')
    product = smiles('CCCl')
    reaction = ReactionContainer([reactant], [reagent], [product])
    
    reaction.clean2d()
    
    # Check all molecules have coordinates
    for molecule in reaction.molecules():
        assert all(n in molecule._plane for n in molecule)
    
    # Check reagents are positioned above the arrow
    reagent_coords = [(x, y) for mol in reaction.reagents 
                     for x, y in mol._plane.values()]
    assert all(y > 0 for _, y in reagent_coords)  # reagents should be above
    
    # Verify arrow position
    arrow_start, arrow_end = reaction._arrow
    assert arrow_start < arrow_end
    
    # Check relative positioning
    reactant_max_x = max(x for mol in reaction.reactants 
                        for x, _ in mol._plane.values())
    product_min_x = min(x for mol in reaction.products 
                       for x, _ in mol._plane.values())
    assert reactant_max_x < arrow_start < arrow_end < product_min_x


def test_fix_positions():
    # Test just the position fixing functionality
    reaction = ReactionContainer([smiles('CCO')], [smiles('CC=O')])
    
    # Clean individual molecules first
    for mol in reaction.molecules():
        mol.clean2d()
    
    # Then fix positions
    reaction.fix_positions()
    
    # Check arrow exists
    assert hasattr(reaction, '_arrow')
    
    # Check molecules are properly spaced
    reactant_coords = [(x, y) for mol in reaction.reactants 
                      for x, y in mol._plane.values()]
    product_coords = [(x, y) for mol in reaction.products 
                     for x, y in mol._plane.values()]
    
    # Verify no overlap between reactants and products
    reactant_max_x = max(x for x, _ in reactant_coords)
    product_min_x = min(x for x, _ in product_coords)
    assert reactant_max_x < product_min_x 