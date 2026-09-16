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
from chython.core import read_smiles
from chython.core._core import _cip_subtree_probe


def _key(smiles, branch, anchor):
    """The canonical unfolding of one branch, comparable ACROSS molecules.

    An interned id is only meaningful inside its own context -- it is a counter -- so two molecules are
    compared by key and one molecule's shared subtrees by id.
    """
    return _cip_subtree_probe(read_smiles(smiles), ((branch, anchor),), keys=True)[0]


def test_two_identical_ligands_intern_to_one_id():
    """The two methyls of isopropanol are one subtree, and that identity IS the rule-1 tie proof."""
    mol = read_smiles('CC(O)C')
    ids = _cip_subtree_probe(mol, ((1, 2), (4, 2)))
    assert ids[0] == ids[1]


def test_two_different_ligands_intern_apart():
    mol = read_smiles('CC(O)CC')
    ids = _cip_subtree_probe(mol, ((1, 2), (4, 2)))
    assert ids[0] != ids[1]


def test_a_ring_closure_becomes_a_duplicate_and_terminates():
    """The walk is finite in a ring, and the closure leaves a mark: the ring branch is not the open
    chain of the same length, because the duplicate at the closure is a child the chain does not have."""
    assert _key('C1CCCCC1', 2, 1) != _key('CCCCCC', 2, 1)


def test_a_double_bond_gives_the_far_atom_a_duplicate_of_its_partner():
    """Rule 1a sees C=O as C(O,O,dup) -- so an aldehyde branch is not an alcohol branch."""
    assert _key('CC=O', 2, 1) != _key('CCO', 2, 1)


def test_implicit_hydrogens_are_children():
    """They are a count in the arena and atoms nowhere, so only expansion clause 3 puts them in the
    digraph.  A methyl and a methanide differ in nothing else a rule-1/2 word can see -- charge is
    outside CIP's first three rules -- so this pair fails exactly when the clause is missing."""
    assert _key('CC', 2, 1) != _key('C[CH2-]', 2, 1)


def test_the_two_ring_paths_of_a_symmetric_ring_intern_to_one_id():
    mol = read_smiles('C1CCC(O)CC1')
    ids = _cip_subtree_probe(mol, ((2, 1), (7, 1)))
    assert ids[0] == ids[1]


def test_one_ring_branch_has_one_key_however_deep_it_hangs():
    """Rule 1b is stored as a BACK count -- spheres up to the duplicated atom -- and not as an absolute
    root distance, so a ring closure does not reprice its branch once per sphere.  The cyclohexyl below
    is the same node at sphere 1 and at sphere 2, which is what lets the centres of a steroid share
    their subtrees rather than each paying for their own copy."""
    mol = read_smiles('C1CCCCC1CCC(C)C')
    shallow = _cip_subtree_probe(mol, ((6, 7),), keys=True)[0]
    deep = _cip_subtree_probe(mol, ((7, 8),), keys=True)[0]
    assert shallow in deep[2], 'the cyclohexyl one sphere further out is a different node'


def test_the_anchor_is_on_the_path_so_a_ring_closes_onto_a_duplicate():
    """Entering a ring at one atom and walking round it must end at a DUPLICATE of the atom it started
    from.  Left off the path, the anchor would be expanded a second time and the branch would carry a
    whole extra lap of real atoms."""
    keys = _cip_subtree_probe(read_smiles('C1CCCCC1'), ((2, 1),), keys=True)
    # follow the highest-ranked child down: five real ring carbons, then the closing duplicate
    depth = 0
    node = keys[0]
    while node[2]:
        node = max(node[2])
        depth += 1
    assert depth == 5, depth
