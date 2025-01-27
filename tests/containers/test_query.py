# -*- coding: utf-8 -*-
#
#  Copyright 2025 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2025 Tagir Akhmetshin <tagirshin@gmail.com>
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
import pytest
from chython.containers.query import QueryContainer
from chython.containers.bonds import QueryBond
from chython.periodictable import QueryElement


def test_query_container_basic():
    # Test basic container creation and atom addition
    qc = QueryContainer()
    
    # Add atoms with different input types
    n1 = qc.add_atom('C')  # from symbol
    assert isinstance(qc._atoms[n1], QueryElement)
    
    n2 = qc.add_atom(7)  # from atomic number (N)
    assert isinstance(qc._atoms[n2], QueryElement)
    
    n3 = qc.add_atom(QueryElement.from_symbol('O')())  # from QueryElement
    assert isinstance(qc._atoms[n3], QueryElement)


def test_query_container_neighbors():
    # Test neighbors validation and storage
    qc = QueryContainer()
    
    # Test valid neighbors
    n1 = qc.add_atom('C', neighbors=2)  # single value
    assert qc._neighbors[n1] == (2,)
    
    n2 = qc.add_atom('C', neighbors=[1, 2, 3])  # list of values
    assert qc._neighbors[n2] == (1, 2, 3)
    
    # Test invalid neighbors
    with pytest.raises(ValueError):
        qc.add_atom('C', neighbors=-1)  # negative value
    
    with pytest.raises(ValueError):
        qc.add_atom('C', neighbors=15)  # value too large
    
    with pytest.raises(ValueError):
        qc.add_atom('C', neighbors=[1, 1])  # duplicate values


def test_query_container_hybridization():
    # Test hybridization validation and storage
    qc = QueryContainer()
    
    # Test valid hybridization
    n1 = qc.add_atom('C', hybridization=1)  # sp3
    assert qc._hybridizations[n1] == (1,)
    
    n2 = qc.add_atom('C', hybridization=[1, 2])  # sp3 and sp2
    assert qc._hybridizations[n2] == (1, 2)
    
    # Test invalid hybridization
    with pytest.raises(ValueError):
        qc.add_atom('C', hybridization=0)  # invalid value
    
    with pytest.raises(ValueError):
        qc.add_atom('C', hybridization=5)  # invalid value


def test_query_container_rings():
    # Test ring size validation and storage
    qc = QueryContainer()
    
    # Test valid ring sizes
    n1 = qc.add_atom('C', rings_sizes=3)  # 3-membered ring
    assert qc._rings_sizes[n1] == (3,)
    
    n2 = qc.add_atom('C', rings_sizes=[5, 6])  # 5 and 6-membered rings
    assert qc._rings_sizes[n2] == (5, 6)
    
    # Test invalid ring sizes
    with pytest.raises(ValueError):
        qc.add_atom('C', rings_sizes=2)  # too small
    
    with pytest.raises(ValueError):
        qc.add_atom('C', rings_sizes=[5, 5])  # duplicate values


def test_query_container_bonds():
    # Test bond addition and validation
    qc = QueryContainer()
    n1 = qc.add_atom('C')
    n2 = qc.add_atom('C')
    
    # Add bond with different input types
    qc.add_bond(n1, n2, 1)  # from int (single bond)
    assert isinstance(qc._bonds[n1][n2], QueryBond)
    
    qc = QueryContainer()
    n1 = qc.add_atom('C')
    n2 = qc.add_atom('C')
    qc.add_bond(n1, n2, (1, 2))  # from tuple (single or double bond)
    assert isinstance(qc._bonds[n1][n2], QueryBond)


def test_query_container_copy():
    # Test container copying
    qc = QueryContainer()
    n1 = qc.add_atom('C', neighbors=2, hybridization=1)
    n2 = qc.add_atom('N', rings_sizes=6)
    qc.add_bond(n1, n2, 1)
    
    # Make a copy
    copy = qc.copy()
    
    # Verify all attributes are copied
    assert copy._neighbors == qc._neighbors
    assert copy._hybridizations == qc._hybridizations
    assert copy._rings_sizes == qc._rings_sizes
    assert len(copy._bonds) == len(qc._bonds)


def test_query_container_union():
    # Test container union
    qc1 = QueryContainer()
    n1 = qc1.add_atom('C', neighbors=2)
    n2 = qc1.add_atom('O')
    qc1.add_bond(n1, n2, 1)
    
    qc2 = QueryContainer()
    n3 = qc2.add_atom('N', rings_sizes=5)
    n4 = qc2.add_atom('C')
    qc2.add_bond(n3, n4, 2)
    
    # Create union with remapping to avoid collisions
    union = qc1.union(qc2, remap=True)
    
    # Verify union properties
    assert len(union._atoms) == 4  # total number of atoms
    assert len(union._bonds) == 4  # each bond is stored twice (bidirectional)
    assert sum(len(bonds) for bonds in union._bonds.values()) == 4  # total number of bond entries
    assert len(union._neighbors) == len(qc1._neighbors) + len(qc2._neighbors)
    assert len(union._rings_sizes) == len(qc1._rings_sizes) + len(qc2._rings_sizes)


def test_query_container_enumerate():
    # Test query enumeration
    qc = QueryContainer()
    n1 = qc.add_atom('C')
    # Add N and O separately to test enumeration
    n2 = qc.add_atom('N')
    n3 = qc.add_atom('O')
    qc.add_bond(n1, n2, (1, 2))  # single or double bond
    qc.add_bond(n1, n3, 1)  # single bond
    
    # Enumerate all possible combinations
    queries = list(qc.enumerate_queries())
    assert len(queries) >= 2  # at least 2 combinations due to bond types
    
    # Test with mark enumeration
    queries = list(qc.enumerate_queries(enumerate_marks=True))
    assert len(queries) >= 2  # should include mark combinations 