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
"""`set_element` -- the field the mutation surface could not write.

The one thing every test here is really about is the STABLE ID.  `set_element` keeps it: deleting the
atom and adding it back with its bonds is the same change spelled so that the atom comes back with a
NEW id, leaving every mapping the caller was holding pointed at nothing.  A reaction patcher pairs a
product atom with the reactant atom it came from BY that id, so a mutator that allocates a new one
is not usable there at all.
"""
import pytest

from chython.core import Atom, MoleculeContainer
from chython.core._core import JOURNAL_OPS, read_smiles as smiles


def test_the_stable_id_survives_the_element_change():
    m = smiles('CC(=O)O')
    before = list(m.atom_numbers)
    with m.edit():
        m.set_element(4, 'N')
    assert list(m.atom_numbers) == before
    assert m.element_of(4) == 7


def test_bonds_and_charge_are_untouched():
    m = smiles('C[O-]')
    with m.edit():
        m.set_element(2, 'S')
    assert m.element_of(2) == 16
    assert m.charge_of(2) == -1
    assert m.order_of(1, 2) == 1
    assert m.degree_of(2) == 1


def test_the_hydrogen_count_is_left_exactly_as_it_was():
    """AND THIS IS THE DESIGN, not an omission.

    A count derived for the old element is almost certainly wrong for the new one, but this layer
    does not derive counts -- `chython.chemistry.calc_implicit` does, and `core` cannot import it.
    Writing `H_UNKNOWN` here would throw away a count the caller is often about to write correctly
    itself; deriving one is what the core refuses everywhere else.  So the count is the caller's,
    and `set_element`'s docstring says so out loud.
    """
    m = smiles('CO')
    assert m.implicit_h_of(2) == 1
    with m.edit():
        m.set_element(2, 'N')
    assert m.implicit_h_of(2) == 1        # nitrogen's would be 2; nobody was asked to derive it


def test_an_unknown_count_stays_unknown():
    m = MoleculeContainer()
    with m.edit():
        n = m.add_atom(6)                 # implicit_h omitted -> H_UNKNOWN
    assert m.implicit_h_of(n) is None
    with m.edit():
        m.set_element(n, 'N')
    assert m.implicit_h_of(n) is None
    assert m.unknown_h_count == 1


def test_element_takes_a_number_or_a_symbol_exactly_as_add_atom_does():
    m = MoleculeContainer()
    with m.edit():
        a = m.add_atom(6)
        b = m.add_atom(6)
    with m.edit():
        m.set_element(a, 7)
        m.set_element(b, 'P')
    assert [m.element_of(i) for i in (a, b)] == [7, 15]
    # and the same refusal for an `Atom`, because `add_atom` does not take one either
    with pytest.raises(NotImplementedError):
        with m.edit():
            m.set_element(a, Atom(16))


def test_a_bad_element_is_refused_and_writes_nothing():
    m = smiles('C')
    with pytest.raises(ValueError):
        with m.edit():
            m.set_element(1, 'Xx')
    assert m.element_of(1) == 6
    with pytest.raises(KeyError):
        with m.edit():
            m.set_element(99, 'N')


def test_one_journal_record_named_like_every_other_mutator():
    m = MoleculeContainer()
    with m.edit():
        n = m.add_atom(6)
    with m.edit():
        m.set_element(n, 'N')
        assert m.journal_length == 1
        assert m.journal_record(0) == (JOURNAL_OPS['set_element'], n, 0, 7, 0, 0)


def test_a_parity_survives_because_its_frame_does():
    """A parity is a statement about a frame of NEIGHBOURS, and an element change leaves the frame
    alone -- same atoms, same CSR order, same directions.  So the sign is kept, and the way to see
    that it is the same configuration and not a coincidence is to build the target directly."""
    m = smiles('C[C@H](N)O')
    assert m.parity_of(2) == 2
    with m.edit():
        m.set_element(4, 'S')
    assert m.parity_of(2) == 2
    assert m.canonical_bytes == smiles('C[C@H](N)S').canonical_bytes


def test_a_stored_cip_descriptor_is_dropped_and_an_isotope_edit_does_not_drop_it():
    """THE ASYMMETRY IS THE POINT.  Both fields feed a CIP ranking, but atomic number is Rule 1 --
    the primary criterion, ahead of mass -- so an atom whose element changed is not the atom the
    input made its assertion about.  An isotope edit only moves Rule 2 and the stored descriptor
    survives it, as it does today."""
    m = smiles('C[C@H](N)O')
    with m.edit():
        m.set_atom_cip(2, 'S')
    assert m.atom_cip_of(2) == 'S'
    with m.edit():
        m.set_isotope(1, 13)
    assert m.atom_cip_of(2) == 'S'
    with m.edit():
        m.set_element(4, 'S')
    assert m.atom_cip_of(2) is None


def test_derived_data_is_rebuilt_around_the_new_element():
    """`heteroatoms_of` counts neighbours that are not carbon, so it moves when an element does --
    which is the cheapest visible proof that the apply rebuilt the derived segments rather than
    patching one field in place."""
    m = smiles('CCC')
    assert m.heteroatoms_of(2) == 0
    with m.edit():
        m.set_element(1, 'O')
    assert m.heteroatoms_of(2) == 1
    # C2H8O and not C2H6O: the methyl's three hydrogens are still stored on what is now an oxygen,
    # because nothing here derives a count.  The formula is the loudest place that shows, and it is
    # the reason a patcher recomputes hydrogens itself.
    assert m.brutto_formula == 'C2H8O'
