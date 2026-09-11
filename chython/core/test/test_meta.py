# -*- coding: utf-8 -*-
#
#  Copyright 2026 Ramil Nugmanov <nougmanoff@protonmail.com>
#
#  This file is part of chython.
#
#  chython is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
"""`meta` is ONE implementation held by both containers, and it is a plain `dict`.

The wrapper record types existed only because the molecule had no `meta`; a second mapping class would
have been the same duplication wearing a core-shaped hat.
"""
from pickle import dumps, loads
from pytest import raises
from chython.core import ReactionContainer, read_smiles as smiles


def test_meta_is_created_on_first_access():
    mol = smiles('CCO')
    assert mol.meta == {}
    mol.meta['boiling_point'] = '78.37'
    assert mol.meta == {'boiling_point': '78.37'}


def test_meta_is_the_same_type_on_both_containers():
    """A plain dict, and literally the same behaviour -- that is what `zero duplication` means."""
    mol = smiles('CCO')
    rxn = ReactionContainer([mol], [smiles('CC=O')])
    assert type(mol.meta) is dict and type(rxn.meta) is dict


def test_meta_has_no_setter():
    """`mol.meta = {}` would let a caller swap the identity a live reference points at."""
    with raises(AttributeError):
        smiles('CCO').meta = {'a': '1'}


def test_copy_carries_meta_and_leaves_the_log_behind():
    mol = smiles('CCO')
    mol.meta['k'] = 'v'
    copy = mol.copy()
    assert copy.meta == {'k': 'v'}
    copy.meta['k'] = 'w'
    assert mol.meta == {'k': 'v'}, 'shallow copy, not a shared dict'


def test_substructure_starts_with_no_meta():
    """A part of a molecule is not the record the metadata described."""
    mol = smiles('CCO')
    mol.meta['k'] = 'v'
    assert mol.substructure([mol.number_of(0), mol.number_of(1)]).meta == {}


def test_pickle_carries_meta():
    """`to_bytes` is the arena and the arena has no field for it, so `__reduce__` carries it beside."""
    mol = smiles('CCO')
    mol.meta['k'] = 'v'
    assert loads(dumps(mol)).meta == {'k': 'v'}


def test_meta_is_not_identity():
    a, b = smiles('CCO'), smiles('CCO')
    a.meta['k'] = 'v'
    assert a == b and hash(a) == hash(b)
