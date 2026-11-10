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
"""BILN text into the monomer graph (Fox et al., J. Chem. Inf. Model. 2022, 62, 3942).

Chains split on ``.``, monomers on ``-``; ``name(bond,Rn)`` marks one end of a bond, and every bond id
must occur exactly twice.  A name is bare up to ``(``, or ``[...]``.  Every chain is a PEPTIDE chain:
BILN has no other polymer type.
"""

from ._errors import BilnError
from ._graph import Edge, Node, PolymerGraph


__all__ = ['parse_biln']


def parse_biln(text):
    """The monomer graph of a BILN string.  Raises :class:`BilnError` for text that does not parse."""
    if not isinstance(text, str) or not text.strip():
        raise BilnError('an empty string is not BILN')
    graph = PolymerGraph()
    ends = {}              # bond id -> [(node, R)]
    for number, chain in enumerate(_split(text.strip(), '.'), 1):
        name = f'PEPTIDE{number}'
        previous = None
        for position, token in enumerate(_split(chain, '-'), 1):
            symbol, bonds = _monomer(token)
            index = len(graph.nodes)
            graph.nodes.append(Node(name, position, 'PEPTIDE', symbol, False))
            if previous is not None:
                graph.edges.append(Edge(previous, 2, index, 1))
            previous = index
            for bond, label in bonds:
                ends.setdefault(bond, []).append((index, label))
    for bond, pair in sorted(ends.items()):
        if len(pair) != 2:
            raise BilnError(f'bond {bond} has {len(pair)} end(s); a BILN bond id occurs exactly twice')
        (a, ra), (b, rb) = pair
        graph.edges.append(Edge(a, ra, b, rb))
    return graph


def _split(text, separator):
    out, depth, start = [], 0, 0
    for i, char in enumerate(text):
        if char in '[(':
            depth += 1
        elif char in '])':
            depth -= 1
            if depth < 0:
                raise BilnError(f'unbalanced {char!r} at {i} in {text!r:.80}')
        elif char == separator and not depth:
            out.append(text[start:i])
            start = i + 1
    if depth:
        raise BilnError(f'unclosed bracket in {text!r:.80}')
    out.append(text[start:])
    if any(not part for part in out):
        raise BilnError(f'an empty monomer or chain in {text!r:.80}')
    return out


def _monomer(token):
    """``(name, [(bond id, R number)])``."""
    if token.startswith('['):
        end = token.find(']')
        name, rest = token[1:end], token[end + 1:]
    else:
        cut = token.find('(')
        name, rest = (token, '') if cut < 0 else (token[:cut], token[cut:])
    if not name:
        raise BilnError(f'{token!r:.40} names no monomer')
    bonds = []
    while rest:
        if not rest.startswith('(') or ')' not in rest:
            raise BilnError(f'{token!r:.40}: {rest!r} is not (bond,R)')
        inner, rest = rest[1:rest.index(')')], rest[rest.index(')') + 1:]
        fields = inner.split(',')
        if len(fields) != 2 or not fields[0].strip().isdigit() or not fields[1].strip().isdigit():
            raise BilnError(f'{token!r:.40}: ({inner}) is not (bond,R) with two integers')
        bonds.append((int(fields[0]), int(fields[1])))
    return name, bonds
