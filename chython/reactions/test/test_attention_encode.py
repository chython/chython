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
"""The four tensors, cell for cell, against the algorithm chython 2 built them with.

WHY A SECOND IMPLEMENTATION IS RIGHT HERE and wrong in the library.  The encoding is the frozen
artefact's input specification, so a reference stating it a second way is a specification and not a
duplicate: the two differ in HOW they get the distances -- Floyd-Warshall over a dense adjacency matrix
here, per-source BFS in C inside `state_view` -- so agreeing on every cell is real evidence rather than
two copies of one mistake.  `chython/core/test/chytorch_oracle.py` makes the same argument.

chython 2's version called `scipy.sparse.csgraph.shortest_path(method='FW')`.  The call is written out
below instead, so the reference costs the tree no dependency; the algorithm is the same one.

NEEDS NO MODEL.  Nothing here loads the weights.
"""
from pytest import importorskip, mark

from ...core import ReactionContainer, read_smiles as smiles


numpy = importorskip('numpy', reason='the encoder answers numpy arrays -- `chython[mapping]`')

MAX_DISTANCE = 10
MAX_NEIGHBORS = 14

#: Reactions covering what the encoder has to get right: a salt inside one container (an unreachable
#: pair), a multi-molecule side, an aromatic ring, a charged atom, an isotope, a single atom, and a
#: record whose product side is bigger than its reactant side.
CORPUS = [(['CC(=O)O', 'CCN'], ['CC(=O)NCC', 'O']),
          (['CC(=O)O.[Na+].[Cl-]'], ['CC(=O)[O-].[Na+]', 'Cl']),
          (['c1ccccc1', 'ClCl'], ['c1ccccc1Cl', 'Cl']),
          (['[13CH4]', 'Cl'], ['[13CH3]Cl']),
          (['O'], ['O']),
          (['CCO'], ['CC=O', 'O', 'O']),
          (['[NH4+].[OH-]'], ['N', 'O'])]


def _reference_molecule(molecule):
    """chython 2's `_encode_molecule`, with its scipy call written out.

    `atoms`: 0 padding, 1 mol_cls, else atomic number + 2.  `neighbors`: 0 at mol_cls, else the heavy
    degree plus the implicit hydrogen count, capped at 14, plus 2.  `distances`: 1 between two atoms with
    no path, else the bond count clamped at 10, plus 2; the mol_cls row and column are 1 throughout.
    """
    n_atoms = len(molecule)
    size = n_atoms + 1

    atoms = numpy.zeros(size, dtype='int32')
    neighbors = numpy.zeros(size, dtype='int32')
    atoms[0] = 1
    for i, n in enumerate(molecule.atom_numbers, 1):
        atoms[i] = molecule.atom(n).element + 2
        count = molecule.degree_of(n) + (molecule.atom(n).implicit_h or 0)
        neighbors[i] = min(count, MAX_NEIGHBORS) + 2

    distance = _floyd_warshall(molecule.adjacency_matrix().astype(bool))
    # -1 BEFORE THE SHIFT, so an unreachable pair lands on 1 -- the value below every real distance,
    # which starts at 2 for an atom to itself.
    numpy.nan_to_num(distance, copy=False, posinf=-1.)
    numpy.clip(distance, None, MAX_DISTANCE, out=distance)
    distance = (distance + 2).astype('int32')

    distances = numpy.ones((size, size), dtype='int32')
    distances[1:, 1:] = distance
    return atoms, neighbors, distances


def _floyd_warshall(adjacency):
    """All-pairs shortest paths over an unweighted undirected graph.  `inf` where there is no path."""
    n = adjacency.shape[0]
    out = numpy.full((n, n), numpy.inf)
    numpy.fill_diagonal(out, 0.)
    out[adjacency] = 1.
    for k in range(n):
        numpy.minimum(out, out[:, k, None] + out[None, k, :], out=out)
    return out


def _reference_reaction(reactants, products):
    """chython 2's `_encode_reaction`: one token row over the record, distances block diagonal in it."""
    atoms = [numpy.zeros(1, dtype='int32')]     # rxn_cls, sharing the padding value
    neighbors = [numpy.zeros(1, dtype='int32')]
    roles = [1]                                 # rxn_cls
    blocks = []

    for molecules, role in ((reactants, 2), (products, 3)):
        for molecule in molecules:
            a, n, d = _reference_molecule(molecule)
            atoms.append(a)
            neighbors.append(n)
            blocks.append(d)
            roles.append(0)                     # mol_cls
            roles.extend([role] * len(molecule))

    atoms = numpy.concatenate(atoms)
    neighbors = numpy.concatenate(neighbors)
    roles = numpy.array(roles, dtype='int32')

    total = len(roles)
    distances = numpy.zeros((total, total), dtype='int32')
    distances[0, 0] = 1                         # the rxn_cls self-loop
    position = 1
    for block in blocks:
        end = position + block.shape[0]
        distances[position:end, position:end] = block
        position = end
    return atoms, neighbors, distances, roles


def _reaction(record):
    left, right = record
    return ReactionContainer([smiles(s) for s in left], [smiles(s) for s in right])


@mark.parametrize('record', CORPUS, ids=lambda r: '.'.join(r[0]) + '>>' + '.'.join(r[1]))
def test_every_tensor_agrees_with_the_reference(record):
    from ..attention._encode import encode_reaction

    rxn = _reaction(record)
    got = encode_reaction(rxn.reactants, rxn.products)
    want = _reference_reaction(rxn.reactants, rxn.products)
    for name, mine, theirs in zip(('atoms', 'neighbors', 'distances', 'roles'), got[:4], want):
        assert mine.shape == theirs.shape, name
        assert numpy.array_equal(mine, theirs), f'{name} disagrees at ' \
                                                f'{numpy.argwhere(mine != theirs).tolist()}'


def test_the_token_row_is_the_documented_layout():
    from ..attention._encode import encode_reaction

    rxn = ReactionContainer([smiles('CC')], [smiles('O')])
    encoded = encode_reaction(rxn.reactants, rxn.products)
    #        rxn_cls  mol_cls  C  C   mol_cls  O
    assert encoded.roles.tolist() == [1, 0, 2, 2, 0, 3]
    assert encoded.atoms.tolist() == [0, 1, 8, 8, 1, 10]        # 6 + 2 carbon, 8 + 2 oxygen
    assert encoded.neighbors.tolist() == [0, 0, 6, 6, 0, 4]     # CH3 is 1 + 3 + 2, water 0 + 2 + 2
    assert encoded.distances[0].tolist() == [1, 0, 0, 0, 0, 0]  # the rxn_cls self-loop and padding
    assert encoded.distances[1].tolist() == [0, 1, 1, 1, 0, 0]  # mol_cls, distance 1 to its own atoms
    assert encoded.distances[2].tolist() == [0, 1, 2, 3, 0, 0]  # C: 0 + 2 to itself, 1 + 2 to the other


def test_an_unreachable_pair_is_one_and_not_a_shifted_distance():
    from ..attention._encode import encode_reaction

    rxn = ReactionContainer([smiles('[Na+].[Cl-]')], [smiles('[Na+].[Cl-]')])
    block = encode_reaction(rxn.reactants, rxn.products).distances[1:4, 1:4]
    assert block.tolist() == [[1, 1, 1],        # mol_cls reaches both ions
                              [1, 2, 1],        # Na to itself is 0 + 2, to Cl there is no path
                              [1, 1, 2]]


def test_a_distance_past_the_clamp_saturates():
    from ..attention._encode import encode_reaction

    chain = smiles('C' * 20)
    rxn = ReactionContainer([chain], [smiles('C')])
    distances = encode_reaction(rxn.reactants, rxn.products).distances
    assert distances[1:21, 1:21].max() == MAX_DISTANCE + 2


def test_an_unstated_hydrogen_count_counts_as_zero_and_not_as_the_sentinel():
    from ...core import H_UNKNOWN, MoleculeContainer
    from ..attention._encode import encode_reaction

    # H_UNKNOWN is 15 and the column's domain is 0..14, so the sentinel must not reach the tensor: an
    # unstated count contributes nothing, which is what `TensorEncoding.unknown_h` defaults to.
    mol = MoleculeContainer()
    with mol.edit() as e:
        n = e.add_atom('N', implicit_h=H_UNKNOWN)
        c = e.add_atom('C', implicit_h=2)
        e.add_bond(n, c, 1)
    assert mol.atom(n).implicit_h is None

    neighbors = encode_reaction([mol], [smiles('C')]).neighbors
    assert neighbors[2] == 1 + 0 + 2                     # the nitrogen: one heavy neighbour, no count
    assert neighbors[3] == 1 + 2 + 2                     # the carbon, which stated two


def test_an_r_marker_encodes_as_two_and_not_as_padding():
    from ..attention._encode import encode_reaction

    # Element 0 shifted by 2.  The weights never saw the token; the record is still encoded, because a
    # reader that refuses an R marker refuses a Markush record it was handed.
    rxn = ReactionContainer([smiles('*CC')], [smiles('C')])
    marker = smiles('*CC')
    assert marker.atom(marker.atom_numbers[0]).element == 0
    assert encode_reaction(rxn.reactants, rxn.products).atoms[2] == 2


def test_the_equality_mask_forbids_a_cross_element_correspondence():
    from ..attention._encode import encode_reaction

    rxn = ReactionContainer([smiles('CO')], [smiles('CO')])
    equal = encode_reaction(rxn.reactants, rxn.products).equal_atoms
    assert equal.shape == (2, 2)
    assert equal.diagonal().all()               # C to C and O to O
    assert not equal[0, 1] and not equal[1, 0]  # C to O never


def test_agents_are_absent_from_the_tensors():
    from ..attention._encode import encode_reaction

    without = encode_reaction([smiles('CC')], [smiles('CC')])
    rxn = ReactionContainer([smiles('CC')], [smiles('CC')], [smiles('[Pd]')])
    withal = encode_reaction(rxn.reactants, rxn.products)
    assert numpy.array_equal(without.atoms, withal.atoms)
    assert numpy.array_equal(without.roles, withal.roles)
