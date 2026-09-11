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
"""`transition_view` on a molecule: the before == after case, which fixes the layout."""
from pytest import raises

from chython.core import H_UNKNOWN, MoleculeContainer, TensorEncoding, read_smiles


def test_a_molecule_reports_the_same_state_on_both_sides():
    """A molecule is a reaction that does nothing, so one vocabulary serves both containers."""
    view = read_smiles('CCO').transition_view()
    assert view.elements.tolist() == [6, 6, 8]
    assert view.h_before.tolist() == [3, 2, 1]
    assert view.h_after.tolist() == [3, 2, 1]
    assert view.n_before.tolist() == [1, 2, 1]
    assert view.n_after.tolist() == [1, 2, 1]
    assert view.distances.tolist() == [[0, 1, 2], [1, 0, 1], [2, 1, 0]]


def test_the_bond_list_is_pairs_of_union_indices_with_i_less_than_j():
    view = read_smiles('CCO').transition_view()
    assert view.bonds.tolist() == [[0, 1], [1, 2]]
    assert view.bonds.shape == (2, 2)
    assert view.bonds.shape[0] == 2
    assert view.bond_before.tolist() == [1, 1]
    assert view.bond_after.tolist() == [1, 1]


def test_a_double_and_an_aromatic_bond_report_their_orders_on_both_sides():
    view = read_smiles('CC=O').transition_view()
    assert view.bond_before.tolist() == [1, 2]
    assert view.bond_after.tolist() == [1, 2]
    view = read_smiles('c1ccccc1').transition_view()
    assert set(view.bond_before.tolist()) == {4}, 'aromatic order is 4 on both sides'


def test_map_numbers_come_out_as_the_record_states_them():
    view = read_smiles('[CH4:7]').transition_view()
    assert view.map_numbers.tolist() == [7]
    assert read_smiles('C').transition_view().map_numbers.tolist() == [0], 'unmapped is 0'


def test_a_molecule_reports_nothing_unmapped_and_no_collision():
    """The two fields exist on both containers so a consumer needs one branch, not two."""
    view = read_smiles('CCO').transition_view()
    assert view.unmapped == {'reactants': 0, 'products': 0}
    assert view.collisions == {'reactants': (), 'products': ()}


def test_every_array_is_contiguous_int32():
    view = read_smiles('CCO').transition_view()
    for name in ('elements', 'h_before', 'n_before', 'h_after', 'n_after', 'map_numbers',
                 'distances', 'bonds', 'bond_before', 'bond_after'):
        array = getattr(view, name)
        assert array.dtype.name == 'int32', name
        assert array.flags['C_CONTIGUOUS'], name


def test_an_unstated_hydrogen_count_reports_unknown_h_on_both_sides():
    mol = MoleculeContainer()
    with mol.edit():
        n = mol.add_atom('N', implicit_h=H_UNKNOWN)
        for _ in range(4):
            mol.add_bond(n, mol.add_atom('F', implicit_h=0), 1)
    assert mol.transition_view().h_before.tolist()[0] == 0
    assert mol.transition_view().h_after.tolist()[0] == 0
    sentinel = mol.transition_view(TensorEncoding(unknown_h=15))
    assert sentinel.h_before.tolist()[0] == 15
    assert sentinel.h_after.tolist()[0] == 15


def test_the_encoding_applies_the_same_way_it_does_to_the_state_view():
    """There is no `hydrogens` column here: the four per-side columns replace it."""
    enc = TensorEncoding(element_shift=2, hydrogen_shift=1, neighbor_shift=3, distance_shift=4,
                         disconnected=1)
    view = read_smiles('CO.[Na+]').transition_view(enc)
    assert not hasattr(view, 'hydrogens')
    assert view.elements.tolist() == [8, 10, 13]
    assert view.h_before.tolist() == [4, 2, 1]
    assert view.n_before.tolist() == [4, 4, 3]
    assert view.distances[0].tolist() == [4, 5, 1], 'the sodium is disconnected, verbatim'


def test_a_bond_order_is_never_shifted():
    """A bond order is a chemical value with its own small domain; a shift would collide with 0."""
    view = read_smiles('CC=O').transition_view(TensorEncoding(element_shift=2, distance_shift=2))
    assert view.bond_before.tolist() == [1, 2]


def test_a_token_key_reads_both_sides():
    enc = TensorEncoding(vocabulary={(6, 3, 1, 3, 1): 10, (8, 1, 1, 1, 1): 12}, unknown=999)
    assert read_smiles('CO').transition_view(enc).tokens.tolist() == [10, 12]
    assert read_smiles('CO').transition_view().tokens is None


def test_width_pads_the_atom_columns_and_leaves_the_bond_list_alone():
    """A bond list is ragged by nature; padding it would need a sentinel row nothing asked for."""
    view = read_smiles('CCO').transition_view(TensorEncoding(width=5, pad=0, pad_diagonal=1))
    assert view.elements.shape == (5,)
    assert view.distances.shape == (5, 5)
    assert view.bonds.shape == (2, 2), 'two bonds, unpadded'
    assert view.distances[4].tolist() == [0, 0, 0, 0, 1]


def test_a_molecule_larger_than_width_raises_and_names_its_atom_count():
    with raises(ValueError, match='10 atoms'):
        read_smiles('C' * 10).transition_view(TensorEncoding(width=6))


def test_a_bond_free_and_an_empty_molecule_give_empty_bond_arrays_of_the_right_rank():
    view = read_smiles('[He].[He]').transition_view()
    assert view.bonds.shape == (0, 2)
    assert view.bond_before.shape == (0,)
    view = MoleculeContainer().transition_view()
    assert view.elements.shape == (0,)
    assert view.bonds.shape == (0, 2)
    assert view.distances.shape == (0, 0)


def test_transition_view_agrees_with_state_view_where_they_overlap():
    """A molecule's transition view must agree cell-for-cell with its state view.

    `state_view` is differentially verified against an external implementation, so agreement
    is the strongest available check that the transition kernel is correct.  The H_UNKNOWN
    cases exercise the `unknown_h` branch, which has a structurally different path from the
    counting branch; a disagreement there would not surface from SMILES-only inputs.
    """
    from chython.core import mol_state_view

    def _check(mol, enc, label):
        sv = mol_state_view(mol, enc)
        tv = mol.transition_view(enc)
        assert sv.elements.tolist() == tv.elements.tolist(), label
        assert sv.hydrogens.tolist() == tv.h_before.tolist(), label
        assert sv.hydrogens.tolist() == tv.h_after.tolist(), label
        assert sv.neighbors.tolist() == tv.n_before.tolist(), label
        assert sv.neighbors.tolist() == tv.n_after.tolist(), label
        assert sv.distances.tolist() == tv.distances.tolist(), label

    enc = TensorEncoding(element_shift=1, hydrogen_shift=1, neighbor_shift=1, distance_shift=1,
                         disconnected=-1, max_distance=5)
    for smi in ('CCO', 'c1ccccc1', 'CO.[Na+]', 'C(=O)O', '[CH4:3]'):
        _check(read_smiles(smi), enc, smi)

    # H_UNKNOWN: unreachable from SMILES; exercise the unknown_h branch under two encodings.
    mol = MoleculeContainer()
    with mol.edit():
        n = mol.add_atom('N', implicit_h=H_UNKNOWN)
        for _ in range(4):
            mol.add_bond(n, mol.add_atom('F', implicit_h=0), 1)
    _check(mol, TensorEncoding(unknown_h=0), 'H_UNKNOWN/unknown_h=0')
    _check(mol, TensorEncoding(unknown_h=5), 'H_UNKNOWN/unknown_h=5')


def test_padded_diagonal_is_pad_when_pad_diagonal_is_off():
    """pad_diagonal=0 means off: padded rows carry `pad`, not 0, on the diagonal.

    Fails if `_ml_fill_transition_arrays` writes 0 unconditionally instead of only when
    `pad_diagonal` is non-zero — the same condition `mol_state_view` uses.
    """
    enc = TensorEncoding(width=5, pad=99, pad_diagonal=0)
    view = read_smiles('CCO').transition_view(enc)
    for i in range(3, 5):
        assert view.distances[i, i] == 99, f'padded diagonal row {i} should be pad=99'
