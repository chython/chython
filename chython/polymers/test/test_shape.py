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
"""The shape key: a candidate unit's lookup key read off the molecule itself, with no fragment built."""

from chython.core import smiles

from .._library import get_monomer_library
from .._shape import Shape, shape_of_core


def key_of(mol, atoms, boundary):
    """The key of `atoms` of `mol`, with `boundary` the cut bonds (inside atom, outside atom)."""
    mol.thiele()
    shape = Shape(mol)
    return shape.join([shape.piece(atoms)], internal=(), boundary=boundary)


def test_a_residue_in_a_peptide_has_its_monomer_core_key():
    # Ala-Gly-Gly: the middle glycine is N6 C7 C8 O9, cut from C4 and N10
    m = smiles('C[C@H](N)C(=O)NCC(=O)NCC(=O)O')
    lib = get_monomer_library()
    key = key_of(m, (6, 7, 8, 9), ((6, 4), (8, 10)))
    assert ('G', frozenset((1, 2))) in [(x.symbol, u) for x, u in lib.by_shape(key)]


def test_the_key_ignores_stereo_and_atom_order():
    a = key_of(smiles('CN[C@@H](C)C(C)=O'), (2, 3, 4, 5, 7), ((2, 1), (5, 6)))
    b = key_of(smiles('CC(=O)[C@H](C)NC'), (2, 3, 4, 5, 6), ((2, 1), (6, 7)))
    assert a == b


def test_a_group_key_is_the_union_of_its_pieces():
    m = smiles('C[C@H](N)C(=O)NCC(=O)NCC(=O)O')
    m.thiele()
    shape = Shape(m)
    whole = shape.join([shape.piece((6, 7, 8, 9))], internal=(), boundary=((6, 4), (8, 10)))
    split = shape.join([shape.piece((6, 7)), shape.piece((8, 9))], internal=((7, 8),),
                       boundary=((6, 4), (8, 10)))
    assert whole == split


def test_core_and_molecule_keys_agree():
    lib = get_monomer_library()
    core, names = lib.labelled_core(lib.get('PEPTIDE', 'P'), frozenset((1, 2)))
    m = smiles('CC(=O)N1CCC[C@H]1C(=O)NC')
    assert shape_of_core(core) == key_of(m, (4, 5, 6, 7, 8, 9, 10), ((4, 2), (9, 11)))
