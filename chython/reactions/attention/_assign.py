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
"""One product atom to one reactant atom, greedily, from the attention matrix.

NOT A SOLVER, and the difference is the point: this walks outward from the strongest correspondence in
the matrix, taking the best remaining cell in the neighbourhood of what it has already placed.  A
maximum-weight matching would score higher on the matrix and worse on the chemistry, because attention
between two atoms is evidence about their neighbourhoods and not an independent cost.

numpy and nothing else here.  The model, the encoding and the container never enter, which is what makes
the walk testable on a matrix written out by hand.
"""
from numpy import argmax, bool_, full, isclose, ix_, mean, nonzero, ones, unravel_index, zeros


def side_adjacency(molecules):
    """The block-diagonal `[n, n]` bool adjacency over one side of the reaction.

    Block diagonal because a side is a set of molecules and the walk's neighbourhood must not cross
    between two of them.  Row order is each molecule's `atom_numbers`, concatenated -- the order
    `encode_reaction` lays the tokens out in, so a row index means the same atom in both.
    """
    total = sum(len(m) for m in molecules)
    out = zeros((total, total), dtype=bool_)
    position = 0
    for molecule in molecules:
        end = position + len(molecule)
        out[position:end, position:end] = molecule.adjacency_matrix().astype(bool_)
        position = end
    return out


def greedy_mapping(attention, r_adj, p_adj, multiplier):
    """`(assignment, score)`: per product atom the reactant atom index it takes, or -1, and the mean.

    `attention` IS CONSUMED -- accepted rows and columns are zeroed in place and neighbourhoods are
    scaled, so a caller wanting it afterwards copies it first.

    The walk:

    1. the strongest cell anywhere starts it;
    2. afterwards the search is restricted to the frontier -- product atoms bonded to one already
       placed -- and falls back to the whole matrix when the frontier empties, which is how a second
       disconnected product component gets started;
    3. an accepted cell scales its own neighbourhood by `multiplier`, biasing the next choice towards a
       correspondence that keeps a bond intact;
    4. the accepted row and column are zeroed, so each atom is used once;
    5. a maximum of zero ends it, and every product atom still unplaced stays unplaced.

    The score is read off a COPY TAKEN BEFORE ANY SCALING, so a cell's contribution is what the model
    said about it and not what the walk did to its neighbourhood afterwards.
    """
    products, reactants = attention.shape
    assignment = full(products, -1, dtype='int64')
    if not products or not reactants:
        return assignment, 0.

    raw = attention.copy()
    frontier = zeros(products, dtype=bool_)
    unplaced = ones(products, dtype=bool_)
    score = []

    for step in range(products):
        if step and frontier.any():
            rows = nonzero(frontier)[0]
            i, j = unravel_index(argmax(attention[frontier]), (rows.shape[0], reactants))
            i = rows[i]
        else:
            i, j = unravel_index(argmax(attention), attention.shape)

        if isclose(attention[i, j], 0.):
            break

        score.append(raw[i, j])
        assignment[i] = j
        attention[ix_(p_adj[i], r_adj[j])] *= multiplier
        attention[i] = 0
        attention[:, j] = 0
        unplaced[i] = False
        frontier[i] = False
        frontier[p_adj[i] & unplaced] = True

    return assignment, float(mean(score)) if score else 0.


__all__ = ['greedy_mapping', 'side_adjacency']
