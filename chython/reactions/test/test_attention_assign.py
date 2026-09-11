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
"""The greedy walk on matrices written out by hand.

WHY THE MATRICES ARE HAND-WRITTEN.  `greedy_mapping` takes numpy and nothing else, so every rule it
follows -- the frontier, the neighbourhood multiplier, the zero stop -- is stateable as four numbers and
an expected assignment.  A test that ran the model to produce the matrix would prove the pair works and
say nothing about which of the two decided the outcome.

NEEDS NO MODEL.
"""
from pytest import importorskip

from ...core import read_smiles as smiles


numpy = importorskip('numpy', reason='the assignment is numpy -- `chython[mapping]`')


def _chain(n):
    """Path adjacency over `n` atoms: 0-1-2-...  The multiplier's neighbourhood needs a real one."""
    out = numpy.zeros((n, n), dtype=bool)
    for i in range(n - 1):
        out[i, i + 1] = out[i + 1, i] = True
    return out


def _isolated(n):
    """No bonds, so the frontier never fills and every step searches the whole matrix."""
    return numpy.zeros((n, n), dtype=bool)


def test_a_clean_diagonal_maps_straight_through():
    from ..attention._assign import greedy_mapping

    attention = numpy.array([[.9, .1, .1],
                             [.1, .8, .1],
                             [.1, .1, .7]])
    assignment, score = greedy_mapping(attention, r_adj=_isolated(3), p_adj=_isolated(3), multiplier=1.75)
    assert assignment.tolist() == [0, 1, 2]
    assert score == numpy.mean([.9, .8, .7])


def test_the_score_is_read_before_the_multiplier_touches_the_cell():
    from ..attention._assign import greedy_mapping

    # Both product atoms are bonded, so accepting (0, 0) scales the whole 2x2 by 1.75 -- including the
    # cell taken next.  The score must report .4 and not .7.
    attention = numpy.array([[.9, .2],
                             [.2, .4]])
    _, score = greedy_mapping(attention, r_adj=_chain(2), p_adj=_chain(2), multiplier=1.75)
    assert score == numpy.mean([.9, .4])


def test_the_multiplier_changes_which_cell_wins():
    from ..attention._assign import greedy_mapping

    # (0, 0) is the strongest cell and starts the walk.  Product 1 is bonded to product 0 and reactant 1
    # to reactant 0, so the neighbourhood cell (1, 1) is scaled; (1, 2) is not.
    def matrix():
        return numpy.array([[.9, .1, .1],
                            [.1, .30, .40]])

    without = greedy_mapping(matrix(), r_adj=_chain(3), p_adj=_chain(2), multiplier=1.)[0]
    assert without.tolist() == [0, 2], 'unscaled, the raw maximum .40 wins'

    scaled = greedy_mapping(matrix(), r_adj=_chain(3), p_adj=_chain(2), multiplier=1.75)[0]
    assert scaled.tolist() == [0, 1], '.30 x 1.75 = .525 beats .40, keeping the bond intact'


def test_a_zero_row_leaves_its_product_atom_unplaced():
    from ..attention._assign import greedy_mapping

    attention = numpy.array([[.9, .1],
                             [0., 0.]])
    assignment, score = greedy_mapping(attention, r_adj=_isolated(2), p_adj=_isolated(2), multiplier=1.75)
    assert assignment.tolist() == [0, -1]
    assert score == .9, 'a step that never happened contributes nothing to the mean'


def test_an_all_zero_matrix_places_nothing():
    from ..attention._assign import greedy_mapping

    assignment, score = greedy_mapping(numpy.zeros((3, 2)), r_adj=_isolated(2), p_adj=_isolated(3), multiplier=1.75)
    assert assignment.tolist() == [-1, -1, -1]
    assert score == 0.


def test_more_product_atoms_than_reactant_atoms_exhausts_the_columns():
    from ..attention._assign import greedy_mapping

    attention = numpy.array([[.9], [.5], [.3]])
    assignment, _ = greedy_mapping(attention, r_adj=_isolated(1), p_adj=_isolated(3), multiplier=1.75)
    assert assignment.tolist() == [0, -1, -1], 'one column, so two product atoms have nothing left'


def test_an_empty_side_returns_an_empty_assignment_and_not_a_numpy_error():
    from ..attention._assign import greedy_mapping

    assignment, score = greedy_mapping(numpy.zeros((0, 3)), r_adj=_isolated(3), p_adj=_isolated(0), multiplier=1.75)
    assert assignment.tolist() == []
    assert score == 0.
    assignment, score = greedy_mapping(numpy.zeros((3, 0)), r_adj=_isolated(0), p_adj=_isolated(3), multiplier=1.75)
    assert assignment.tolist() == [-1, -1, -1]
    assert score == 0.


def test_the_frontier_restarts_on_a_second_component():
    from ..attention._assign import greedy_mapping

    # Products 0-1 are bonded and 2 stands alone.  After both of the first pair are placed the frontier
    # is empty; the walk must fall back to the whole matrix rather than stop.
    p_adj = numpy.zeros((3, 3), dtype=bool)
    p_adj[0, 1] = p_adj[1, 0] = True
    attention = numpy.array([[.9, .1, .1],
                             [.1, .8, .1],
                             [.1, .1, .2]])
    assignment, _ = greedy_mapping(attention, r_adj=_isolated(3), p_adj=p_adj, multiplier=1.75)
    assert assignment.tolist() == [0, 1, 2]


def test_each_atom_is_used_once():
    from ..attention._assign import greedy_mapping

    # One reactant column dominates every row; only the first product atom may take it.
    attention = numpy.array([[.9, .1],
                             [.8, .2],
                             [.7, .3]])
    assignment, _ = greedy_mapping(attention, r_adj=_isolated(2), p_adj=_isolated(3), multiplier=1.75)
    assert sorted(assignment.tolist()) == [-1, 0, 1]
    assert assignment[0] == 0, 'the strongest cell anywhere starts the walk'


def test_the_side_adjacency_is_block_diagonal_over_the_molecules():
    from ..attention._assign import side_adjacency

    out = side_adjacency([smiles('CC'), smiles('CC')])
    assert out.shape == (4, 4)
    assert out.tolist() == [[False, True, False, False],
                            [True, False, False, False],
                            [False, False, False, True],
                            [False, False, True, False]]


def test_the_side_adjacency_row_order_is_atom_numbers_concatenated():
    from ..attention._assign import side_adjacency

    # Row `i` must mean the same atom as token `i` of that side; the two orders agreeing is what makes
    # the assignment's indices readable back onto atoms.
    propanol = smiles('CCO')
    out = side_adjacency([propanol])
    assert numpy.array_equal(out, propanol.adjacency_matrix().astype(bool))


def test_an_empty_side_has_an_empty_adjacency():
    from ..attention._assign import side_adjacency

    assert side_adjacency([]).shape == (0, 0)
