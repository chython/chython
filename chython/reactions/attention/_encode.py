# -*- coding: utf-8 -*-
#
#  Copyright 2026 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2024 Philippe Gantzer <p.gantzer@icredd.hokudai.ac.jp>
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
"""The reaction as the four tensors the weights take, built on `state_view`.

EVERY NUMBER HERE IS THE ARTEFACT'S CONTRACT AND NOT A DESIGN CHOICE.  The shifts, the clamps, the
sentinel for an unreachable pair and the four role codes are what `chython-rxnmap` was trained with;
changing one silently returns a plausible wrong attention matrix rather than an error.

Three of the four tensors are `state_view` columns under one `TensorEncoding`.  `neighbors` is the
exception and stays a sum here, because the clamp is on degree-plus-hydrogens and the ML views refuse a
combined column: one number cannot say which half of it the record left unstated.

    tensor      value
    atoms       0 padding and rxn_cls, 1 mol_cls, else atomic number + 2
    neighbors   0 at both cls tokens, else min(heavy degree + implicit H, 14) + 2
    distances   0 padding, 1 unreachable or cross-component, else bond count + 2, clamped at 10 first
    roles       0 mol_cls, 1 rxn_cls, 2 reactant atom, 3 product atom

The token layout is one row of tokens over the whole record:

    [rxn_cls] [mol_cls, atoms...] x each reactant  [mol_cls, atoms...] x each product

`distances` is block diagonal over that row -- one block per molecule, `1` everywhere between blocks,
so a reactant atom and a product atom attend to each other through the cross-component value rather
than through a path that does not exist.
"""
from typing import NamedTuple

from numpy import array, bool_, concatenate, empty, int64, ix_, minimum, ones, zeros

from ._session import get_session
from ...core import TensorEncoding


#: The distance clamp, applied before the shift.
MAX_DISTANCE = 10
#: The clamp on heavy degree plus implicit hydrogens.  ALSO the hypervalence guard's limit, on heavy
#: degree ALONE -- the two tests are not the same one, and a record passing the guard can still clamp.
MAX_NEIGHBORS = 14

ROLE_MOL_CLS = 0
ROLE_RXN_CLS = 1
ROLE_REACTANT = 2
ROLE_PRODUCT = 3

#: Built once and read field by field on every call: the shifts are the artefact's, and a fresh
#: encoding per reaction would be a kwargs dict against a budget measured in microseconds.
_ENCODING = TensorEncoding(element_shift=2, distance_shift=2, disconnected=1,
                           max_distance=MAX_DISTANCE)


class Encoded(NamedTuple):
    """The model's four inputs, plus what the assignment needs to read its output.

    The masks and the equality matrix are built in the same pass as the tensors ON PURPOSE.  They are
    all statements about one token layout, and two functions describing that layout separately is one
    place for it to drift.
    """
    atoms: object
    neighbors: object
    distances: object
    roles: object
    #: `ix_` selectors picking the product-by-reactant and reactant-by-product blocks of the `[seq, seq]`
    #: attention matrix.
    p2r: tuple
    r2p: tuple
    #: `[product atoms, reactant atoms]` bool: a correspondence between two different elements is not a
    #: correspondence, whatever the attention says.
    equal_atoms: object


def encode_molecule(molecule):
    """`(atoms, neighbors, distances)` for one molecule, its `mol_cls` token first.

    `mol_cls` carries atoms 1, neighbors 0 and distance 1 to every atom of the molecule, which is how
    the graph gets a per-molecule summary token to attend through.
    """
    view = molecule.state_view(_ENCODING)
    size = view.elements.shape[0] + 1

    atoms = empty(size, dtype=view.elements.dtype)
    atoms[0] = 1
    atoms[1:] = view.elements

    neighbors = empty(size, dtype=view.neighbors.dtype)
    neighbors[0] = 0
    neighbors[1:] = minimum(view.neighbors + view.hydrogens, MAX_NEIGHBORS) + 2

    distances = ones((size, size), dtype=view.distances.dtype)
    distances[1:, 1:] = view.distances
    return atoms, neighbors, distances


def encode_reaction(reactants, products) -> Encoded:
    """The whole record as one token row.

    AGENTS ARE NOT PASSED IN, and the caller numbering them afterwards is not an oversight: the weights
    were trained on reactants and products only.  A record may carry a mapped catalyst; the model has
    never been shown one.
    """
    atoms = [zeros(1, dtype='int32')]  # rxn_cls, which shares the padding value
    neighbors = [zeros(1, dtype='int32')]
    roles = [ROLE_RXN_CLS]
    blocks = []
    r_mask = [False]
    p_mask = [False]

    for molecules, role, mask, other in ((reactants, ROLE_REACTANT, r_mask, p_mask),
                                         (products, ROLE_PRODUCT, p_mask, r_mask)):
        for molecule in molecules:
            a, n, d = encode_molecule(molecule)
            atoms.append(a)
            neighbors.append(n)
            blocks.append(d)
            roles.append(ROLE_MOL_CLS)
            roles.extend([role] * (a.shape[0] - 1))
            mask.append(False)  # this side's mol_cls is not one of its atoms
            mask.extend([True] * (a.shape[0] - 1))
            other.extend([False] * a.shape[0])

    atoms = concatenate(atoms)
    neighbors = concatenate(neighbors)
    roles = array(roles, dtype='int32')

    total = atoms.shape[0]
    distances = zeros((total, total), dtype='int32')
    distances[0, 0] = 1  # the rxn_cls self-loop; every other cell touching it stays padding
    position = 1
    for block in blocks:
        end = position + block.shape[0]
        distances[position:end, position:end] = block
        position = end

    r_mask = array(r_mask, dtype=bool_)
    p_mask = array(p_mask, dtype=bool_)
    r_elements = atoms[r_mask]
    p_elements = atoms[p_mask]
    return Encoded(atoms, neighbors, distances, roles,
                   ix_(p_mask, r_mask), ix_(r_mask, p_mask),
                   p_elements[:, None] == r_elements)


def run_model(encoded: Encoded, threads: int):
    """The `[product atoms, reactant atoms]` attention, symmetrized across the arrow and element masked.

    Symmetrized because the two blocks are two readings of one correspondence: `p2r` is how much each
    product atom attends to each reactant atom and `r2p` the reverse, and their sum is the only score
    that does not depend on which side the question was asked from.
    """
    attention = get_session(threads).run(None, {
        'atoms': encoded.atoms[None].astype(int64),
        'neighbors': encoded.neighbors[None].astype(int64),
        'distances': encoded.distances[None].astype(int64),
        'roles': encoded.roles[None].astype(int64)})[0]
    return (attention[encoded.p2r] + attention[encoded.r2p].T) * encoded.equal_atoms


__all__ = ['MAX_DISTANCE', 'MAX_NEIGHBORS', 'ROLE_MOL_CLS', 'ROLE_PRODUCT', 'ROLE_REACTANT',
           'ROLE_RXN_CLS', 'Encoded', 'encode_molecule', 'encode_reaction', 'run_model']
