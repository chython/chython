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
"""An R is added and indexed through the ordinary edit session, like any other atom kind."""
from pytest import raises
from chython.core import MoleculeContainer


def _benzene_with_r():
    mol = MoleculeContainer()
    with mol.edit() as e:
        ring = [e.add_atom('C') for _ in range(6)]
        for i in range(6):
            e.add_bond(ring[i], ring[(i + 1) % 6], 4)
        r = e.add_atom('R')
        e.add_bond(ring[0], r, 1)
    return mol, ring, r


def test_add_atom_by_symbol():
    mol, ring, r = _benzene_with_r()
    assert mol.atom(r).is_r
    assert mol.atom(r).atomic_symbol == 'R'
    assert mol.atom_count == 7


def test_add_atom_with_an_index_in_the_symbol():
    mol = MoleculeContainer()
    with mol.edit() as e:
        r = e.add_atom('R2')
    assert mol.atom(r).r_index == 2
    assert mol.atom(r).atomic_symbol == 'R2'


def test_add_atom_by_number_zero():
    mol = MoleculeContainer()
    with mol.edit() as e:
        r = e.add_atom(0)
    assert mol.atom(r).is_r


def test_set_r_index():
    mol, ring, r = _benzene_with_r()
    with mol.edit() as e:
        e.set_r_index(r, 7)
    assert mol.atom(r).r_index == 7
    assert mol.atom(r).atomic_symbol == 'R7'


def test_set_r_index_on_an_element_is_refused():
    # Raised on replay, at scope exit -- the element is only final there.
    mol, ring, r = _benzene_with_r()
    with raises(ValueError, match='not an R'):
        with mol.edit() as e:
            e.set_r_index(ring[0], 1)


def test_r_index_above_the_domain_is_refused():
    from chython.core import R_INDEX_MAX

    mol = MoleculeContainer()
    with raises(ValueError, match=str(R_INDEX_MAX)):
        with mol.edit() as e:
            e.add_atom(f'R{R_INDEX_MAX + 1}')


def test_r_index_survives_bytes_round_trip():
    # The R index lives in atom_t.reserved bits 4-11, which enter to_bytes() directly.
    # Two different indices in one molecule: a single index cannot pass by always returning the same value.
    mol = MoleculeContainer()
    with mol.edit() as e:
        r3 = e.add_atom('R3')
        r7 = e.add_atom('R7')
    mol2 = MoleculeContainer.from_bytes(mol.to_bytes())
    assert mol2.atom(r3).r_index == 3
    assert mol2.atom(r7).r_index == 7
    assert mol2.atom(r3).atomic_symbol == 'R3'
    assert mol2.atom(r7).atomic_symbol == 'R7'
