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
"""An R survives the lossless form and is refused by the lossy one.  Nothing writes it unreadably."""
from pytest import raises
from chython.core import MoleculeContainer


def _r_and_carbon(index=7):
    mol = MoleculeContainer()
    with mol.edit() as e:
        r = e.add_atom('R')
        c = e.add_atom('C')
        e.add_bond(r, c, 1)
    if index:
        with mol.edit() as e:
            e.set_r_index(r, index)
    return mol, r


def test_to_bytes_round_trips_an_r_with_its_index():
    mol, r = _r_and_carbon()
    back = MoleculeContainer.unpack(mol.to_bytes())
    assert [back.atom(sid).atomic_symbol for sid in back] == ['R7', 'C']
    assert back.canonical_bytes == mol.canonical_bytes


def test_pach_refuses_a_molecule_holding_an_r():
    mol, r = _r_and_carbon()
    with raises(ValueError, match='pach'):
        mol.pack(version=2)


def test_the_pach_refusal_names_the_lossless_alternative():
    mol, r = _r_and_carbon()
    with raises(ValueError) as caught:
        mol.pack(version=2)
    assert 'to_bytes' in str(caught.value)


def test_dropping_every_field_does_not_waive_the_r_refusal():
    # `drop` waives losing a FIELD; losing an atom is not a field.
    mol, r = _r_and_carbon()
    with raises(ValueError):
        mol.pack(drop='*', version=2)


def test_an_unindexed_r_is_refused_too():
    mol, r = _r_and_carbon(index=0)
    with raises(ValueError):
        mol.pack(version=2)


def test_a_molecule_with_no_r_still_packs():
    mol = MoleculeContainer()
    with mol.edit() as e:
        e.add_bond(e.add_atom('C'), e.add_atom('O'), 1)
    assert MoleculeContainer.unpack(mol.pack(version=2)).brutto == mol.brutto
