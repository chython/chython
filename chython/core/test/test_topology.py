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
"""The topology surface: `split`, the augmented environment, and the two matrices."""
from importlib.util import find_spec
from pytest import mark, raises

from chython.core import MoleculeContainer, read_smiles


# THE TWO MATRICES ARE THE NUMPY HALF OF THIS FILE AND `split`/`augmented_*` ARE NOT, which is why the
# marker is per-test rather than on the module: numpy is an optional dependency (`chython[ml]`) and a
# minimal install must still be held to `split`'s components and the augmented environment's meaning of
# radius.  `test_the_topology_surface_agrees_with_chython_two` is marked because it compares the
# adjacency matrix among other things -- V2's answer is fetched in one subprocess call, so the witness
# is all-or-nothing rather than split into a matrix half and a shells half.
#
# `find_spec` rather than `importorskip`, so collection does not import numpy at all -- the same
# reasoning as `interop/test/conftest.py` gives for the optional toolkits.
needs_numpy = mark.skipif(find_spec('numpy') is None,
                          reason='numpy is not installed; the matrices are arrays')


# --- split ---------------------------------------------------------------------------------------

def test_split_of_two_components_returns_both_in_first_seen_order():
    # ethanol and sodium chloride written as one record
    mol = read_smiles('CCO.[Na+].[Cl-]')
    parts = mol.split()
    assert len(parts) == 3
    assert [p.atom_count for p in parts] == [3, 1, 1]
    assert [sorted(p.atom_numbers) for p in parts] == [[1, 2, 3], [4], [5]]


def test_split_of_a_salt_gives_the_ion_pair():
    # sodium acetate: the carboxylate and the counter-ion
    mol = read_smiles('CC(=O)[O-].[Na+]')
    anion, cation = mol.split()
    assert str(anion) == 'C(C)([O-])=O'
    assert str(cation) == '[Na+]'
    assert anion.union(cation) == mol


def test_split_of_one_component_returns_a_list_of_one_molecule():
    mol = read_smiles('c1ccccc1')
    parts = mol.split()
    assert isinstance(parts, list)
    assert len(parts) == 1
    assert isinstance(parts[0], MoleculeContainer)
    assert parts[0] is not mol
    assert parts[0] == mol


def test_split_preserves_atom_numbers():
    mol = read_smiles('CCO.[Na+]')
    mol.remap({1: 10, 2: 20, 3: 30, 4: 40})
    parts = mol.split()
    assert [sorted(p.atom_numbers) for p in parts] == [[10, 20, 30], [40]]
    assert parts[0].element_of(30) == 8


def test_split_preserves_stereo():
    # (S)-butan-2-ol with two spectator ions; the split cuts no bond, so no parity frame moves
    mol = read_smiles('C[C@H](O)CC.[Na+].[Cl-]')
    before = mol.parity_of(2)
    assert before  # the fixture must actually carry a configuration
    part = mol.split()[0]
    assert part.parity_of(2) == before
    assert str(part) == 'C(C)[C@@H](O)C'


def test_split_of_a_single_atom_is_one_component():
    assert len(read_smiles('[Na+]').split()) == 1


def test_split_of_an_empty_molecule_is_an_empty_list():
    assert MoleculeContainer().split() == []


# --- augmented_substructure ----------------------------------------------------------------------

def test_augmented_substructure_radius_zero_is_the_seed_itself():
    mol = read_smiles('CC(=O)OCC')      # ethyl acetate; atoms 1..6 in written order
    sub = mol.augmented_substructure([2], 0)
    assert sorted(sub.atom_numbers) == [2]


def test_augmented_substructure_radius_one_is_the_seed_plus_its_neighbours():
    mol = read_smiles('CC(=O)OCC')
    # atom 2 is the carbonyl carbon, bonded to 1 (methyl), 3 (=O) and 4 (ester O)
    assert sorted(mol.neighbors_of(2)) == [1, 3, 4]
    sub = mol.augmented_substructure([2], 1)
    assert sorted(sub.atom_numbers) == [1, 2, 3, 4]
    assert sub.bond_count == 3


def test_augmented_substructure_defaults_to_radius_one():
    mol = read_smiles('CC(=O)OCC')
    assert sorted(mol.augmented_substructure([2]).atom_numbers) == [1, 2, 3, 4]


def test_augmented_substructure_keeps_bonds_inside_the_selection_only():
    mol = read_smiles('CC(=O)OCC')
    sub = mol.augmented_substructure([2], 1)
    # 4-5 leaves the selection and is gone; 2-4 is inside and stays
    assert sub.order_of(2, 4) == 1
    assert 5 not in sub


def test_augmented_substructure_seed_of_two_atoms_unions_both_environments():
    mol = read_smiles('CC(=O)OCC')
    sub = mol.augmented_substructure([1, 6], 1)
    assert sorted(sub.atom_numbers) == [1, 2, 5, 6]
    # the two halves are separate: the cut left no path between them
    assert sub.connected_components_count == 2


def test_augmented_substructure_saturates_at_the_whole_component():
    mol = read_smiles('CCO')
    assert sorted(mol.augmented_substructure([1], 99).atom_numbers) == [1, 2, 3]


def test_augmented_substructure_does_not_leave_its_component():
    mol = read_smiles('CCO.[Na+]')
    assert sorted(mol.augmented_substructure([1], 99).atom_numbers) == [1, 2, 3]


def test_augmented_substructures_yields_every_shell_including_the_seed():
    mol = read_smiles('CC(=O)OCC')
    levels = mol.augmented_substructures([2], 2)
    assert [sorted(s.atom_numbers) for s in levels] == [[2], [1, 2, 3, 4], [1, 2, 3, 4, 5]]


def test_augmented_substructures_stops_growing_when_the_component_is_covered():
    mol = read_smiles('CCO')
    levels = mol.augmented_substructures([1], 99)
    assert [sorted(s.atom_numbers) for s in levels] == [[1], [1, 2], [1, 2, 3]]


def test_augmented_substructure_returns_a_molecule_not_a_view():
    mol = read_smiles('CCO')
    sub = mol.augmented_substructure([1], 0)
    assert isinstance(sub, MoleculeContainer)
    assert not sub.shares_arena_with(mol)


def test_augmented_substructure_rejects_an_unknown_atom():
    with raises(KeyError):
        read_smiles('CCO').augmented_substructure([9])


def test_augmented_substructure_rejects_an_empty_seed():
    with raises(ValueError):
        read_smiles('CCO').augmented_substructure([])


def test_augmented_substructure_rejects_a_negative_radius():
    with raises(ValueError):
        read_smiles('CCO').augmented_substructure([1], -1)


# --- adjacency_matrix ----------------------------------------------------------------------------

@needs_numpy
def test_adjacency_matrix_is_ones_and_symmetric():
    mol = read_smiles('CC=O')
    adj = mol.adjacency_matrix()
    assert adj.shape == (3, 3)
    assert adj.tolist() == [[0, 1, 0], [1, 0, 1], [0, 1, 0]]
    assert (adj == adj.T).all()


@needs_numpy
def test_adjacency_matrix_set_bonds_writes_the_order():
    mol = read_smiles('CC=O')
    adj = mol.adjacency_matrix(True)
    assert adj.tolist() == [[0, 1, 0], [1, 0, 2], [0, 2, 0]]
    assert (adj == adj.T).all()


@needs_numpy
def test_adjacency_matrix_set_bonds_is_also_a_keyword():
    mol = read_smiles('CC=O')
    assert (mol.adjacency_matrix(set_bonds=True) == mol.adjacency_matrix(True)).all()


@needs_numpy
def test_adjacency_matrix_rows_follow_the_molecules_own_atom_order():
    mol = read_smiles('CCO')
    mol.remap({1: 7, 2: 8, 3: 9})
    assert mol.atom_numbers == [7, 8, 9]
    assert mol.adjacency_matrix().tolist() == [[0, 1, 0], [1, 0, 1], [0, 1, 0]]


@needs_numpy
def test_adjacency_matrix_of_an_aromatic_ring_carries_order_four():
    mol = read_smiles('c1ccccc1')
    adj = mol.adjacency_matrix(True)
    assert {int(x) for x in adj.flatten()} == {0, 4}


@needs_numpy
def test_adjacency_matrix_of_a_single_atom_is_a_zero():
    assert read_smiles('[Na+]').adjacency_matrix().tolist() == [[0]]


@needs_numpy
def test_adjacency_matrix_of_an_empty_molecule_is_empty():
    assert MoleculeContainer().adjacency_matrix().shape == (0, 0)


# --- distance_matrix -----------------------------------------------------------------------------

@needs_numpy
def test_distance_matrix_of_a_chain_is_exact():
    mol = read_smiles('CCCC')      # butane: a path of four
    d = mol.distance_matrix()
    assert d.shape == (4, 4)
    assert d.tolist() == [[0, 1, 2, 3],
                          [1, 0, 1, 2],
                          [2, 1, 0, 1],
                          [3, 2, 1, 0]]


@needs_numpy
def test_distance_matrix_of_a_ring_takes_the_short_way_round():
    mol = read_smiles('C1CCCCC1')  # cyclohexane: the far pair is 3, never 5
    d = mol.distance_matrix()
    assert d.tolist() == [[0, 1, 2, 3, 2, 1],
                          [1, 0, 1, 2, 3, 2],
                          [2, 1, 0, 1, 2, 3],
                          [3, 2, 1, 0, 1, 2],
                          [2, 3, 2, 1, 0, 1],
                          [1, 2, 3, 2, 1, 0]]


@needs_numpy
def test_distance_matrix_diagonal_is_zero():
    d = read_smiles('c1ccccc1O').distance_matrix()
    assert [d[i, i] for i in range(7)] == [0] * 7


@needs_numpy
def test_distance_matrix_marks_a_disconnected_pair_with_minus_one():
    mol = read_smiles('CCO.[Na+]')
    d = mol.distance_matrix()
    assert d.dtype.name == 'int32'
    assert d.tolist() == [[0, 1, 2, -1],
                          [1, 0, 1, -1],
                          [2, 1, 0, -1],
                          [-1, -1, -1, 0]]


@needs_numpy
def test_distance_matrix_is_symmetric():
    d = read_smiles('CC(C)C(=O)Nc1ccccc1').distance_matrix()
    assert (d == d.T).all()


@needs_numpy
def test_distance_matrix_of_a_single_atom_is_a_zero():
    assert read_smiles('[Na+]').distance_matrix().tolist() == [[0]]


@needs_numpy
def test_distance_matrix_of_an_empty_molecule_is_empty():
    assert MoleculeContainer().distance_matrix().shape == (0, 0)


@needs_numpy
def test_distance_matrix_rows_follow_the_molecules_own_atom_order():
    mol = read_smiles('CCO')
    mol.remap({1: 7, 2: 8, 3: 9})
    assert mol.distance_matrix().tolist() == [[0, 1, 2], [1, 0, 1], [2, 1, 0]]


@needs_numpy
def test_distance_matrix_counts_a_dative_bond_as_a_bond():
    # order 8 is excluded from ring perception but is still a path for a walk
    mol = MoleculeContainer()
    with mol.edit():
        n = mol.add_atom('N')
        b = mol.add_atom('B')
        f = mol.add_atom('F')
        mol.add_bond(n, b, 8)
        mol.add_bond(b, f, 1)
    assert mol.distance_matrix().tolist() == [[0, 1, 2], [1, 0, 1], [2, 1, 0]]


@needs_numpy
def test_the_shifted_distance_matrix_is_the_encoding_chytorch_documents():
    # chytorch's `graph_distances` adds 2 and clamps: 1 means "different components",
    # 2 "an atom with itself", 3 "neighbours".  -1 is chosen so that shift lands them exactly.
    mol = read_smiles('CCO.[Na+]')
    d = mol.distance_matrix() + 2
    assert d[0, 3] == 1
    assert d[0, 0] == 2
    assert d[0, 1] == 3


# --- the chython 2 witness -----------------------------------------------------------------------
#
# Not an oracle: every answer above is stated outright.  What this catches is the thing a stated
# answer cannot -- that `split`'s components and `augmented_substructure`'s meaning of "radius" agree
# with chython 2 on molecules nobody thought to write a case for.  Reached through `oracle`, an
# installed chython 2 in another interpreter, so this file imports no chython 2.
#
# `distance_matrix` IS ABSENT HERE, because chython 2 has no counterpart to compare against.

COMPOUNDS = ('CCO', 'CC(=O)OCC', 'c1ccccc1O', 'CC(C)C(=O)Nc1ccccc1', 'CCO.[Na+]',
             'C1CCCCC1', 'OCC1OC(O)C(O)C(O)C1O', 'CN1C=NC2=C1C(=O)N(C)C(=O)N2C',
             'CC(=O)Nc1ccc(O)cc1', 'C[N+](C)(C)C.[Cl-]', 'N', 'OS(=O)(=O)O',
             'C1CC2CCC1CC2', 'c1ccc2ccccc2c1')

SEEDS = ((1,), (1, 2), (2,))
RADII = (0, 1, 2, 99)

# `_augmented_substructure` and not `augmented_substructure`: chython 2's private method returns the
# atom-number LEVELS, which is what the shells have to be compared as.  Comparing the molecules it
# builds would compare two SMILES writers instead.
V2_VALUES = """
from chython import smiles

out = []
for smi in _payload:
    mol = smiles(smi)
    levels = []
    for seed in ((1,), (1, 2), (2,)):
        for deep in (0, 1, 2, 99):
            try:
                levels.append([sorted(x) for x in mol._augmented_substructure(list(seed), deep)])
            except ValueError:                 # a seed atom this molecule does not have
                levels.append(None)
    out.append({'levels': levels,
                'components': [sorted(x) for x in mol.connected_components],
                'adjacency': mol.adjacency_matrix().tolist(),
                'orders': mol.adjacency_matrix(True).tolist()})
_emit(out)
"""


@needs_numpy
def test_distance_matrix_of_a_two_hundred_atom_chain_is_exact_end_to_end():
    """Large enough that a truncated half-edge copy leaves a reachable atom unreachable."""
    mol = read_smiles('C' * 200)
    dist = mol.distance_matrix()
    assert dist.shape == (200, 200)
    assert dist[0, 199] == 199
    assert dist[0, 100] == 100
    assert int(dist.max()) == 199
    assert int(dist.min()) == 0


@needs_numpy
def test_distance_matrix_of_three_components_marks_every_cross_pair():
    """Three components, so a component index off by one shows as a real distance across the gap."""
    mol = read_smiles('CC.OO.NN')
    dist = mol.distance_matrix()
    assert dist.shape == (6, 6)
    for i in range(6):
        for j in range(6):
            same = i // 2 == j // 2
            assert (dist[i, j] >= 0) == same, (i, j, dist[i, j])


@needs_numpy
def test_the_topology_surface_agrees_with_chython_two():
    from .oracle import ask

    answers = ask(V2_VALUES, list(COMPOUNDS))
    assert len(answers) == len(COMPOUNDS)
    for smi, old in zip(COMPOUNDS, answers):
        mol = read_smiles(smi)
        assert [sorted(p.atom_numbers) for p in mol.split()] == old['components'], smi
        assert mol.adjacency_matrix().tolist() == old['adjacency'], smi
        assert mol.adjacency_matrix(True).tolist() == old['orders'], smi

        i = 0
        for seed in SEEDS:
            for deep in RADII:
                want = old['levels'][i]
                i += 1
                if want is None:
                    # chython 2 refuses an unknown seed atom with ValueError and this refuses with
                    # KeyError; both refuse, and which exception is not what this witness is for
                    with raises(KeyError):
                        mol.augmented_substructures(list(seed), deep)
                    continue
                got = [sorted(s.atom_numbers) for s in mol.augmented_substructures(list(seed), deep)]
                assert got == want, f'{smi} seed={seed} deep={deep}'
