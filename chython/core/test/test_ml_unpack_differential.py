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
"""`state_view` against chytorch's compiled `_unpack`, cell for cell.

FLOYD-WARSHALL AGAINST PER-SOURCE BFS.  Two unrelated shortest-path implementations agreeing on every
pair of every record is the strongest available check that the arrays are faithful; a same-algorithm
comparison would agree on a shared mistake.

THE NEIGHBOUR COLUMN IS SUMMED HERE AND NOT CLAMPED THERE.  `_unpack` clamps `degree + hydrogens` as
one number; `state_view.neighbors` is the heavy degree alone, by design, because a single number cannot
say which half was unstated.  So the encoding leaves `max_neighbors` unset and this file does the sum,
which is what a consumer does.
"""
from pytest import mark

from chython.core import TensorEncoding, read_smiles

from . import chytorch_oracle
from .chytorch_oracle import requires_oracle


MAX_NEIGHBORS = 14
MAX_DISTANCE = 10

# public compounds, each reaching something: a chain longer than max_distance, two components, a
# heteroaromatic ring, a charged pair, a metal, a single atom, and a branch-heavy centre.
SMILES = [
    'CCO',                                          # ethanol
    'C' * 30,                                       # triacontane: past max_distance
    'c1ccccc1',                                     # benzene
    'c1ccncc1',                                     # pyridine
    'CC(=O)Oc1ccccc1C(=O)O',                        # aspirin
    'CN1C=NC2=C1C(=O)N(C)C(=O)N2C',                 # caffeine
    'CC(C)(C)c1ccccc1',                             # tert-butylbenzene: degree 4 centre
    '[Na+].[Cl-]',                                  # two components
    'C(=O)([O-])[O-].[Ca+2]',                       # calcium carbonate: charges plus a metal
    'O',                                            # water
    '[Fe]',                                         # one atom, no bond
    'NC(=O)N',                                      # urea
    'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O',       # glucose
    'CC(N)C(=O)NC(C)C(=O)NC(C)C(=O)O',              # a trialanine peptide
]


def _chython_columns(mol):
    """The three columns under an encoding reproducing `_unpack`'s conventions."""
    enc = TensorEncoding(element_shift=2, neighbor_shift=2, distance_shift=2,
                         disconnected=1, unknown_h=0, max_distance=MAX_DISTANCE)
    view = mol.state_view(enc)
    neighbors = [min(int(n) - 2 + int(h), MAX_NEIGHBORS) + 2
                 for n, h in zip(view.neighbors, view.hydrogens)]
    return view.elements.tolist(), neighbors, view.distances.tolist()


@requires_oracle
@mark.parametrize('line', SMILES)
def test_the_three_columns_agree_cell_for_cell(line):
    mol = read_smiles(line)
    record = mol.pack(compressed=False, version=2)
    theirs, = chytorch_oracle.unpack([record], MAX_NEIGHBORS, MAX_DISTANCE)
    elements, neighbors, distances = _chython_columns(mol)

    assert elements == theirs['atoms']
    assert neighbors == theirs['neighbors']
    assert distances == theirs['distances']


@requires_oracle
def test_the_atom_order_of_a_pach_v2_record_is_the_container_order():
    """The premise every other assertion in this file rests on, established rather than assumed.

    A cell-for-cell agreement would also hold under a shared permutation, so the columns are checked
    against a molecule whose every atom has a distinct element.
    """
    mol = read_smiles('BCNOFP')
    record = mol.pack(compressed=False, version=2)
    theirs, = chytorch_oracle.unpack([record], MAX_NEIGHBORS, MAX_DISTANCE)
    assert theirs['atoms'] == [a.element + 2 for a in mol.atoms()]


@requires_oracle
def test_the_whole_corpus_crosses_in_one_child():
    """One subprocess per record is the harness's cost; check the batch path agrees too."""
    mols = [read_smiles(line) for line in SMILES]
    theirs = chytorch_oracle.unpack([m.pack(compressed=False, version=2) for m in mols],
                                    MAX_NEIGHBORS, MAX_DISTANCE)
    assert len(theirs) == len(mols)
    for mol, other in zip(mols, theirs):
        elements, neighbors, distances = _chython_columns(mol)
        assert (elements, neighbors, distances) == (other['atoms'], other['neighbors'],
                                                    other['distances'])


@requires_oracle
def test_the_harness_can_disagree():
    """A negative control: a differential that cannot fail is not measuring anything."""
    mol = read_smiles('CCO')
    theirs, = chytorch_oracle.unpack([mol.pack(compressed=False, version=2)],
                                     MAX_NEIGHBORS, MAX_DISTANCE)
    elements, _, _ = _chython_columns(read_smiles('CCN'))
    assert elements != theirs['atoms'], 'the channel reports agreement for two different molecules'


def test_the_oracle_child_has_chytorch_and_no_chython():
    """Skips when chytorch is absent; fails when the isolation the differential rests on is gone."""
    chytorch_oracle.require()
    info = chytorch_oracle.verify()
    assert not info['chython']
    assert not info['chytorch']
