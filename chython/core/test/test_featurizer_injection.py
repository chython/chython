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
"""Core-side injection tests: the slot, the hook and the error path.

These three tests belong in the core test suite because they are facts about what the CORE
owns: the ten names on MoleculeContainer, the `_set_featurizer_fns` export, and the error that
fires when a key is absent.  If they lived only under `chython/chemistry/test/`, then
`pytest chython/core/` would prove nothing about the core's ownership of the ten names.

DUPLICATION OF `PROPERTIES` AND `METHODS` IS DELIBERATE.  A shared helper would have to live
under one layer and be imported from the other, which is the import this split exists to prevent.
The chemistry-side file carries the same two tuples with the same comment.
"""
import pytest
from chython.core import MoleculeContainer
from chython.core import _core


PROPERTIES = ('rotatable_bonds_count', 'hydrogen_bond_donors_count',
              'hydrogen_bond_acceptors_count', 'tpsa', 'crippen_logp', 'crippen_mr', 'qed')
METHODS = ('maccs_keys', 'maccs_bit_set', 'pharmacophore_invariants')


@pytest.mark.parametrize('name', PROPERTIES + METHODS)
def test_name_exists_on_the_sealed_container(name):
    assert hasattr(MoleculeContainer, name), name


def test_the_hook_is_exported():
    assert callable(_core._set_featurizer_fns)


def test_an_unregistered_family_names_its_package_in_the_error():
    # the hook is called by `chython.chemistry`; ask for a key nobody registers
    with pytest.raises(ImportError, match='chython.chemistry'):
        _core._featurizer_fn_for_test('no_such_featurizer')
