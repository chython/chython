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
"""The R index lives in `atom_t.reserved` bits 4-11 and survives a bytes round trip.

Byte surgery rather than a public setter: these tests pin the wire format, so they must fail if the
field moves even when every accessor agrees with itself.
"""
import struct
from pytest import raises
from chython.core import MoleculeContainer


_SEG_ATOMS = 0
_ATOM_RECORD = 24
_ATOM_RESERVED = 20


def _atoms_offset(data):
    return struct.unpack_from('<I', data, 24 + 8 * _SEG_ATOMS)[0]


def _reserved(data, index):
    return struct.unpack_from('<I', data, _atoms_offset(data) + _ATOM_RECORD * index + _ATOM_RESERVED)[0]


def _set_reserved(data, index, value):
    out = bytearray(data)
    struct.pack_into('<I', out, _atoms_offset(data) + _ATOM_RECORD * index + _ATOM_RESERVED, value)
    return bytes(out)


def test_element_zero_with_an_index_loads():
    from chython.core import R_INDEX_MAX

    mol = MoleculeContainer()
    with mol.edit() as e:
        e.add_atom('C')
    data = mol.to_bytes()
    # element 0, R index R_INDEX_MAX: the highest index the domain allows.
    out = bytearray(_set_reserved(data, 0, R_INDEX_MAX << 4))
    struct.pack_into('<B', out, _atoms_offset(data) + _ATOM_RECORD * 0, 0)
    back = MoleculeContainer.from_bytes(bytes(out))
    assert back.atom_count == 1
    assert _reserved(back.to_bytes(), 0) == R_INDEX_MAX << 4


def test_reserved_bit_12_is_still_undefined():
    mol = MoleculeContainer()
    with mol.edit() as e:
        e.add_atom('C')
    data = _set_reserved(mol.to_bytes(), 0, 0x00001000)
    with raises(ValueError, match='reserved bits'):
        MoleculeContainer.from_bytes(data)


def test_r_index_without_element_zero_is_refused():
    mol = MoleculeContainer()
    with mol.edit() as e:
        e.add_atom('C')
    data = _set_reserved(mol.to_bytes(), 0, 1 << 4)
    with raises(ValueError, match='R index'):
        MoleculeContainer.from_bytes(data)


def test_element_symbols_index_zero_is_r():
    from chython.core._core import element_symbols

    symbols = element_symbols()
    assert symbols[0] == 'R'
    assert symbols[1] == 'H'
    assert len(symbols) == 119


def test_atomic_symbol_carries_the_index():
    mol = MoleculeContainer()
    with mol.edit() as e:
        c = e.add_atom('C')
    data = mol.to_bytes()
    out = bytearray(_set_reserved(data, 0, 3 << 4))
    struct.pack_into('<B', out, _atoms_offset(data) + _ATOM_RECORD * 0, 0)
    atom = MoleculeContainer.from_bytes(bytes(out)).atom(c)
    assert atom.is_r
    assert atom.r_index == 3
    assert atom.atomic_symbol == 'R3'


def test_plain_r_has_no_index():
    mol = MoleculeContainer()
    with mol.edit() as e:
        c = e.add_atom('C')
    data = mol.to_bytes()
    out = bytearray(data)
    struct.pack_into('<B', out, _atoms_offset(data) + _ATOM_RECORD * 0, 0)
    atom = MoleculeContainer.from_bytes(bytes(out)).atom(c)
    assert atom.is_r
    assert atom.r_index == 0
    assert atom.atomic_symbol == 'R'


def test_carbon_is_not_an_r():
    mol = MoleculeContainer()
    with mol.edit() as e:
        c = e.add_atom('C')
    assert not mol.atom(c).is_r
    assert mol.atom(c).r_index == 0


def test_the_domain_is_two_decimal_digits():
    from chython.core import R_INDEX_MAX

    # Two digits, so an index fits a three-character CTfile symbol column and a two-character label.
    assert R_INDEX_MAX == 99


def test_the_field_is_wider_than_the_domain():
    # The field holds eight bits and the domain uses 100 of its 256 values.  The relation is the fact
    # worth asserting: the wire format did not shrink, the set of legal values did.
    from chython.core import R_INDEX_MAX

    assert R_INDEX_MAX < 256


def test_the_highest_legal_index_round_trips():
    from chython.core import MoleculeContainer, R_INDEX_MAX

    mol = MoleculeContainer()
    with mol.edit() as e:
        n = e.add_atom('R')
        mol.set_r_index(n, R_INDEX_MAX)
    back = MoleculeContainer.from_bytes(mol.to_bytes())
    assert back.atom(n).r_index == R_INDEX_MAX
    assert back.atom(n).atomic_symbol == 'R99'


def test_an_index_past_the_domain_is_refused_by_the_setter():
    from chython.core import MoleculeContainer, R_INDEX_MAX

    mol = MoleculeContainer()
    with mol.edit() as e:
        n = e.add_atom('R')
    with raises(ValueError, match=str(R_INDEX_MAX)):
        mol.set_r_index(n, R_INDEX_MAX + 1)


def test_an_index_past_the_domain_is_refused_by_the_symbol_parser():
    from chython.core import MoleculeContainer, R_INDEX_MAX

    mol = MoleculeContainer()
    with mol.edit() as e:
        with raises(ValueError, match=str(R_INDEX_MAX)):
            e.add_atom(f'R{R_INDEX_MAX + 1}')


def test_the_symbol_conversion_bounds_the_index_it_accepts():
    """`R500` is not the spelling of an element, so the conversion refuses it rather than reading 0.

    The 0-99 domain belongs to the field, so the answer to "which element is this symbol" honours it
    too.  `Rb` is the guard that the `R` prefix test does not capture an element.
    """
    from chython.core import MoleculeContainer, R_INDEX_MAX

    mol = MoleculeContainer()
    with mol.edit() as e:
        with raises(ValueError, match='unknown element symbol'):
            e.add_atom('R500')
    mol = MoleculeContainer()
    with mol.edit() as e:
        top = e.add_atom(f'R{R_INDEX_MAX}')
        bare = e.add_atom('R')
        rubidium = e.add_atom('Rb')
    assert (mol.atom(top).element, mol.atom(top).r_index) == (0, R_INDEX_MAX)
    assert (mol.atom(bare).element, mol.atom(bare).r_index) == (0, 0)
    assert mol.atom(rubidium).element == 37


def test_a_stored_value_past_the_domain_is_refused_on_load():
    """A record whose field holds 100-255 is rejected rather than read as an index.

    The field is eight bits wide and the domain is 0-99, so those values are representable and
    illegal -- exactly the class `structure_from_bytes` exists to catch.
    """
    mol = MoleculeContainer()
    with mol.edit() as e:
        e.add_atom('C')
    data = mol.to_bytes()
    out = bytearray(_set_reserved(data, 0, 100 << 4))
    struct.pack_into('<B', out, _atoms_offset(data) + _ATOM_RECORD * 0, 0)
    with raises(ValueError, match='R index'):
        MoleculeContainer.from_bytes(bytes(out))
