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
"""The monomer graph: the one intermediate between HELM, BILN and a molecule.

A node is one monomer instance; an edge joins an R group of one node to an R group of another.  In a
chain, neighbours are R2 of the earlier to R1 of the later, and a branch -- a nucleotide's base -- is R3 of
its backbone monomer to R1 of the branch.  A unit that names no one monomer (a BLOB, a repeat range, an
ambiguity list, ``?``) is still a node, listed in `opaque` with the text it was written as.
"""

from typing import NamedTuple


__all__ = ['Node', 'Edge', 'PolymerGraph', 'POLYMER_TYPES']


#: The HELM polymer types.
POLYMER_TYPES = ('PEPTIDE', 'RNA', 'DNA', 'CHEM', 'BLOB')


class Node(NamedTuple):
    chain: str             # 'PEPTIDE1'
    position: int          # 1-based, in chain order
    polymer_type: str
    symbol: str            # a library symbol, an inline SMILES, or the opaque text
    inline: bool           # `symbol` is an inline SMILES monomer


class Edge(NamedTuple):
    a: int
    ra: int
    b: int
    rb: int


class PolymerGraph:
    __slots__ = ('nodes', 'edges', 'opaque', 'log')

    def __init__(self):
        self.nodes = []
        self.edges = []
        self.opaque = []   # (node index, original text)
        self.log = []

    def __repr__(self):
        return f'PolymerGraph({len(self.nodes)} nodes, {len(self.edges)} edges, {len(self.opaque)} opaque)'
