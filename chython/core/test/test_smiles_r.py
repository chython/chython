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
"""`[R]` and `[R<n>]` read as the marker; only inside brackets, and never as a query."""
from pytest import raises
from chython.core import read_smarts, read_smiles as smiles
from chython.core._core import IncorrectSmiles, MoleculeContainer


def test_bracket_r_reads_as_the_marker():
    mol = smiles('[R]c1ccccc1')
    r = next(a for a in mol.atoms() if a.is_r)
    assert r.atomic_symbol == 'R'
    assert r.r_index == 0
    assert mol.atom_count == 7


def test_bracket_r_with_an_index():
    mol = smiles('[R12]c1ccccc1')
    r = next(a for a in mol.atoms() if a.is_r)
    assert r.r_index == 12
    assert r.atomic_symbol == 'R12'


def test_two_indices_in_one_molecule():
    mol = smiles('[R1]CCC[R2]')
    assert sorted(a.r_index for a in mol.atoms() if a.is_r) == [1, 2]


def test_rubidium_is_still_rubidium():
    for spelling in ('[Rb]', '[Ru]', '[Rh]', '[Rn]', '[Re]', '[Ra]', '[Rf]', '[Rg]'):
        mol = smiles(spelling)
        atom = next(iter(mol.atoms()))
        assert not atom.is_r
        assert atom.atomic_symbol == spelling[1:-1]


def test_an_index_past_the_domain_is_refused():
    # The domain is 0..R_INDEX_MAX, and a record naming a higher index is a syntax error rather than a
    # silently truncated marker.  100 is the first illegal value.
    with raises(IncorrectSmiles, match='99'):
        smiles('[R100]C')
    with raises(IncorrectSmiles, match='99'):
        smiles('[R1234]C')


def test_bare_r_outside_brackets_is_refused():
    with raises(IncorrectSmiles, match=r'\[R\]'):
        smiles('Rc1ccccc1')


def test_star_is_the_marker_bare_and_bracketed():
    # `*` is what a stored record spells an attachment point with, so it reads as the marker at index
    # 0 and not as a wildcard query.  Both spellings, and the bracket's other fields keep working.
    for spelling in ('*c1ccccc1', '[*]c1ccccc1'):
        mol = smiles(spelling)
        r = next(a for a in mol.atoms() if a.is_r)
        assert r.atomic_symbol == 'R'
        assert r.r_index == 0
        assert mol.atom_count == 7
    mol = smiles('[13*+:5]C')
    r = next(a for a in mol.atoms() if a.is_r)
    assert (r.isotope, r.charge, r.map_number) == (13, 1, 5)


def test_r_in_smarts_is_still_ring_count():
    # `[R]` in a query means "in at least one ring".  There is no query spelling for the marker.
    query = read_smarts('[R]')
    assert any(query.get_mapping(smiles('c1ccccc1')))
    assert not any(query.get_mapping(smiles('CCC')))


def test_r_carries_a_charge_and_a_map_number():
    # Input is garbage by default: an odd marker is stored, not refused.
    mol = smiles('[R1+:5]C')
    r = next(a for a in mol.atoms() if a.is_r)
    assert r.r_index == 1
    assert r.charge == 1
    assert r.map_number == 5


def test_the_written_form_names_the_marker():
    out = str(smiles('[R]c1ccccc1'))
    assert '[R]' in out
    assert '\x00' not in out


def test_the_written_form_carries_the_index():
    assert '[R12]' in str(smiles('[R12]c1ccccc1'))


def test_round_trip_plain_r():
    once = str(smiles('[R]c1ccccc1'))
    assert str(smiles(once)) == once


def test_index_survives_the_round_trip():
    mol = smiles(str(smiles('[R1]CCC[R2]')))
    assert sorted(a.r_index for a in mol.atoms() if a.is_r) == [1, 2]


def test_charge_and_index_together():
    assert str(smiles('[R1+]C')) in ('C[R1+]', '[R1+]C')


def test_the_marker_is_never_lowercased():
    # `[r]` is the ring-count query primitive, so a lowercase marker would name something else.  An
    # order-4 bond onto a marker is what a MOL file can produce, and it must not change the spelling.
    mol = MoleculeContainer()
    with mol.edit() as e:
        ring = [e.add_atom('C') for _ in range(6)]
        for i in range(6):
            e.add_bond(ring[i], ring[(i + 1) % 6], 4)
        r = e.add_atom('R')
        e.add_bond(ring[0], r, 4)
    assert '[R]' in str(mol)
    assert '[r]' not in str(mol)


def test_an_aromatic_bond_onto_a_marker_survives_the_round_trip():
    # The atom token is uppercase, so the bond token may no longer be suppressed: `:` says order 4.
    mol = MoleculeContainer()
    with mol.edit() as e:
        ring = [e.add_atom('C') for _ in range(6)]
        for i in range(6):
            e.add_bond(ring[i], ring[(i + 1) % 6], 4)
        e.add_bond(ring[0], e.add_atom('R'), 4)
    text = str(mol)
    assert ':[R]' in text or '[R]:' in text, text
    assert sorted(b.order for b in smiles(text).bonds()) == [4] * 7, text


def test_a_single_bond_onto_a_marker_stays_unmarked():
    # The counterpart: order 1 is the empty token, and an R must not acquire a `-` it does not need.
    mol = MoleculeContainer()
    with mol.edit() as e:
        ring = [e.add_atom('C') for _ in range(6)]
        for i in range(6):
            e.add_bond(ring[i], ring[(i + 1) % 6], 4)
        e.add_bond(ring[0], e.add_atom('R'), 1)
    assert sorted(b.order for b in smiles(str(mol)).bonds()) == [1] + [4] * 6


def test_the_written_symbol_table_names_the_marker():
    from chython.core._core import smw_symbol_table
    assert smw_symbol_table()[0] == 'R'
