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
"""`state_view`: the per-atom columns, the distance block, and what padding puts where."""
from pytest import raises

from chython.core import H_UNKNOWN, MoleculeContainer, TensorEncoding, read_smiles


def test_the_columns_of_ethanol_are_its_physical_numbers():
    """The default encoding moves nothing, so this is a readout and not a convention."""
    view = read_smiles('CCO').state_view()
    assert view.elements.tolist() == [6, 6, 8]
    assert view.hydrogens.tolist() == [3, 2, 1]
    assert view.neighbors.tolist() == [1, 2, 1], 'heavy degree, hydrogens excluded'
    assert view.distances.tolist() == [[0, 1, 2], [1, 0, 1], [2, 1, 0]]
    assert view.tokens is None, 'no vocabulary, no token column'


def test_every_column_is_contiguous_int32():
    """A framework wraps these without a copy, which a non-contiguous or int64 array breaks."""
    view = read_smiles('CCO').state_view()
    for name in ('elements', 'hydrogens', 'neighbors', 'distances'):
        array = getattr(view, name)
        assert array.dtype.name == 'int32', name
        assert array.flags['C_CONTIGUOUS'], name


def test_hydrogens_and_neighbors_stay_separate_columns():
    """A single degree+h number cannot say which half was unstated."""
    view = read_smiles('CC(C)C').state_view()
    assert view.neighbors.tolist() == [1, 3, 1, 1]
    assert view.hydrogens.tolist() == [3, 1, 3, 3]


def test_an_unstated_hydrogen_count_reports_unknown_h_and_not_the_sentinel():
    mol = MoleculeContainer()
    with mol.edit():
        n = mol.add_atom('N', implicit_h=H_UNKNOWN)
        for _ in range(4):
            mol.add_bond(n, mol.add_atom('F', implicit_h=0), 1)
    assert mol.state_view().hydrogens.tolist()[0] == 0, 'default unknown_h'
    assert mol.state_view(TensorEncoding(unknown_h=15)).hydrogens.tolist()[0] == 15
    assert mol.state_view(TensorEncoding(unknown_h=7)).hydrogens.tolist()[0] == 7


def test_an_r_marker_reports_element_zero():
    """R is element 0 and matches nothing; it reads as carbon for a neighbour, never for itself."""
    view = read_smiles('C[*]').state_view()
    assert view.elements.tolist() == [6, 0]


def test_a_salt_reports_disconnected_across_its_components():
    view = read_smiles('[Na+].[Cl-]').state_view()
    assert view.distances.tolist() == [[0, -1], [-1, 0]]
    view = read_smiles('[Na+].[Cl-]').state_view(TensorEncoding(disconnected=1))
    assert view.distances.tolist() == [[0, 1], [1, 0]]


def test_the_cross_component_value_is_verbatim_and_never_shifted():
    """A shifted sentinel collides with a real distance, which no consumer can detect."""
    enc = TensorEncoding(distance_shift=2, disconnected=1)
    view = read_smiles('[Na+].[Cl-]').state_view(enc)
    assert view.distances.tolist() == [[2, 1], [1, 2]], 'diagonal shifted, sentinel not'


def test_the_clamp_applies_before_the_shift():
    enc = TensorEncoding(max_distance=2, distance_shift=2)
    view = read_smiles('CCCCC').state_view(enc)
    assert view.distances[0].tolist() == [2, 3, 4, 4, 4], 'clamped to 2, then +2'


def test_the_neighbor_clamp_is_the_heavy_degree_alone():
    """chytorch clamps degree+hydrogens; this column is degree, so the consumer sums then clamps."""
    mol = MoleculeContainer()
    with mol.edit():
        s = mol.add_atom('S', implicit_h=0)
        for _ in range(6):
            mol.add_bond(s, mol.add_atom('F', implicit_h=0), 1)
    assert mol.state_view().neighbors.tolist()[0] == 6
    assert mol.state_view(TensorEncoding(max_neighbors=4)).neighbors.tolist()[0] == 4


def test_each_shift_moves_only_its_own_column():
    enc = TensorEncoding(element_shift=2, hydrogen_shift=1, neighbor_shift=3, distance_shift=4)
    view = read_smiles('CO').state_view(enc)
    assert view.elements.tolist() == [8, 10]
    assert view.hydrogens.tolist() == [4, 2]
    assert view.neighbors.tolist() == [4, 4]
    assert view.distances.tolist() == [[4, 5], [5, 4]]


def test_width_pads_to_a_stackable_shape():
    view = read_smiles('CCO').state_view(TensorEncoding(width=6))
    assert view.elements.shape == (6,)
    assert view.distances.shape == (6, 6)
    assert view.elements.tolist() == [6, 6, 8, 0, 0, 0]
    assert view.hydrogens.tolist() == [3, 2, 1, 0, 0, 0]
    assert view.neighbors.tolist() == [1, 2, 1, 0, 0, 0]


def test_the_pad_value_reaches_every_padded_cell_including_the_distance_block():
    view = read_smiles('CC').state_view(TensorEncoding(width=4, pad=7))
    assert view.elements.tolist() == [6, 6, 7, 7]
    assert view.distances.tolist() == [[0, 1, 7, 7], [1, 0, 7, 7], [7, 7, 7, 7], [7, 7, 7, 7]]


def test_pad_diagonal_leaves_one_unmasked_cell_on_a_padded_row():
    """A fully padded row through a softmax is NaN; one non-masked cell is the fix."""
    view = read_smiles('CC').state_view(TensorEncoding(width=4, pad=0, pad_diagonal=1))
    assert view.distances.tolist() == [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]


def test_a_structure_larger_than_width_raises_and_names_its_atom_count():
    """Truncating to fit is a wrong training example nothing downstream can detect."""
    with raises(ValueError, match='10 atoms'):
        read_smiles('C' * 10).state_view(TensorEncoding(width=6))


def test_a_single_atom_and_a_bond_free_pair_are_records_and_not_errors():
    view = read_smiles('C').state_view()
    assert view.elements.tolist() == [6]
    assert view.distances.tolist() == [[0]]
    view = read_smiles('[He].[He]').state_view()
    assert view.distances.tolist() == [[0, -1], [-1, 0]]


def test_an_empty_molecule_gives_empty_arrays_of_the_right_rank():
    view = MoleculeContainer().state_view()
    assert view.elements.shape == (0,)
    assert view.distances.shape == (0, 0)


def test_an_empty_molecule_still_fills_a_width():
    view = MoleculeContainer().state_view(TensorEncoding(width=3, pad=5))
    assert view.elements.tolist() == [5, 5, 5]
    assert view.distances.shape == (3, 3)


def test_the_view_reads_the_arena_and_not_the_python_atom_objects():
    """An edit session's reseal must be visible; a cached array would not be."""
    mol = read_smiles('CCO')
    first = mol.state_view().elements.tolist()
    with mol.edit() as e:
        e.delete_atom(next(iter(mol.atoms())).n)
    assert mol.state_view().elements.tolist() != first


def test_the_token_column_appears_only_with_a_vocabulary():
    """A key is (element, h_before, n_before, h_after, n_after); a molecule is before == after."""
    vocabulary = {(6, 3, 1, 3, 1): 10, (6, 2, 2, 2, 2): 11, (8, 1, 1, 1, 1): 12}
    view = read_smiles('CCO').state_view(TensorEncoding(vocabulary=vocabulary))
    assert view.tokens.tolist() == [10, 11, 12]
    assert view.tokens.dtype.name == 'int32'


def test_a_missing_key_falls_to_the_unknown_token():
    enc = TensorEncoding(vocabulary={(6, 3, 1, 3, 1): 10}, unknown=999)
    view = read_smiles('CCO').state_view(enc)
    assert view.tokens.tolist() == [10, 999, 999]


def test_the_unknown_token_defaults_to_minus_one_and_not_to_zero():
    """0 is PAD upstream, so a 0 default is a plausible wrong value nobody sees."""
    view = read_smiles('C').state_view(TensorEncoding(vocabulary={(7, 0, 0, 0, 0): 1}))
    assert view.tokens.tolist() == [-1]


def test_a_token_is_never_shifted():
    """The vocabulary owns its id space; a shift would move ids into each other."""
    enc = TensorEncoding(vocabulary={(6, 4, 0, 4, 0): 10}, element_shift=2, hydrogen_shift=2,
                         neighbor_shift=2, distance_shift=2)
    assert read_smiles('C').state_view(enc).tokens.tolist() == [10]


def test_the_key_is_built_before_the_clamps_so_a_clamped_atom_still_finds_its_token():
    """Keying post-clamp yields UNK for every atom in a structure that trips a clamp."""
    mol = MoleculeContainer()
    with mol.edit():
        s = mol.add_atom('S', implicit_h=0)
        for _ in range(6):
            mol.add_bond(s, mol.add_atom('F', implicit_h=0), 1)
    enc = TensorEncoding(vocabulary={(16, 0, 6, 0, 6): 10}, max_neighbors=4, unknown=999)
    assert mol.state_view(enc).tokens.tolist()[0] == 10


def test_the_key_carries_unknown_h_and_not_the_sentinel():
    """The key sees what the hydrogen column reports, which is the whole point of `unknown_h`."""
    mol = MoleculeContainer()
    with mol.edit():
        n = mol.add_atom('N', implicit_h=H_UNKNOWN)
        for _ in range(4):
            mol.add_bond(n, mol.add_atom('F', implicit_h=0), 1)
    zero = TensorEncoding(vocabulary={(7, 0, 4, 0, 4): 10}, unknown=999)
    sentinel = TensorEncoding(vocabulary={(7, 15, 4, 15, 4): 11}, unknown_h=15, unknown=999)
    assert mol.state_view(zero).tokens.tolist()[0] == 10
    assert mol.state_view(sentinel).tokens.tolist()[0] == 11


def test_a_padded_slot_gets_the_pad_value_and_not_the_unknown_token():
    """PAD and UNK are different upstream ids and a padded slot is not a vocabulary miss."""
    enc = TensorEncoding(vocabulary={(6, 4, 0, 4, 0): 10}, width=3, pad=0, unknown=999)
    assert read_smiles('C').state_view(enc).tokens.tolist() == [10, 0, 0]


def test_an_r_marker_keys_on_element_zero():
    enc = TensorEncoding(vocabulary={(0, 0, 1, 0, 1): 10}, unknown=999)
    assert read_smiles('C[*]').state_view(enc).tokens.tolist()[1] == 10


def test_a_lone_r_marker_packs_to_the_empty_slot_key_and_reads_unknown():
    """Element 0, h 0, degree 0 packs to ML_VOCAB_KEY(0,0,0,0,0) = 0 = ML_KEY_EMPTY.

    The entry cannot be stored, so the atom always gets the unknown token regardless of the
    vocabulary supplied.
    """
    enc = TensorEncoding(vocabulary={(6, 4, 0, 4, 0): 10}, unknown=777)
    assert read_smiles('*').state_view(enc).tokens.tolist() == [777]


def test_the_table_capacity_is_at_least_twice_the_entry_count():
    """At load factor 1 every slot is occupied; a probe for a miss loops forever.

    8 entries is the smallest count that fills the mutated table (capacity drops from 16 to 8,
    load factor 1.0).  CH4 is a known miss from this vocabulary.
    """
    vocabulary = {(z, 0, 0, 0, 0): z for z in range(1, 9)}
    enc = TensorEncoding(vocabulary=vocabulary, unknown=999)
    # CH4: element=6, h=4, degree=0 → key (6, 4, 0, 4, 0) is not in the vocabulary
    assert read_smiles('C').state_view(enc).tokens.tolist() == [999]


def test_a_three_hundred_entry_vocabulary_resolves_every_key_it_holds():
    """A probe table degenerates silently when it is too full; 323 entries is the shipped size."""
    mols = [read_smiles(line) for line in ('CCO', 'c1ccccc1', 'CC(=O)Oc1ccccc1C(=O)O', 'NC(=O)N')]
    vocabulary = {}
    for mol in mols:
        view = mol.state_view()
        for z, h, n in zip(view.elements, view.hydrogens, view.neighbors):
            vocabulary.setdefault((int(z), int(h), int(n), int(h), int(n)), len(vocabulary) + 1)
    filler = ((z, h, n, h, n) for z in range(1, 100) for h in range(5) for n in range(5))
    for key in filler:
        if len(vocabulary) >= 320:
            break
        vocabulary.setdefault(key, len(vocabulary) + 1)
    enc = TensorEncoding(vocabulary=vocabulary, unknown=-1)
    for mol in mols:
        assert -1 not in mol.state_view(enc).tokens.tolist(), 'a key present in the table missed'
