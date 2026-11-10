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
"""The monomer graph as HELM and as BILN.

The graph arrives in canonical node order (:func:`~._decompose.decompose`), and every choice below is
made from that order, so the strings are canonical by construction:

| Construct | Written |
|---|---|
| chain | an R2 -> R1 path within one polymer type, from the node nothing links into |
| cyclic backbone | opened at its least node; the closing bond is a connection |
| RNA/DNA base | a branch ``R(A)``: R1 of a node bonded to nothing else, on R3 of a backbone node |
| CHEM monomer | its own chain |
| any other bond | a connection, lower end first, the list in order of its lower end |

Chains are numbered per type in order of their first node.  HELM single-letter symbols are bare, longer
ones and inline SMILES in brackets; the output is ``CHAINS$CONNECTIONS$$$``.  BILN refuses an inline
monomer or a non-PEPTIDE chain: it has neither.
"""

from ._errors import BilnError


__all__ = ['write_helm', 'write_biln']


def write_helm(graph):
    chains, connections = _layout(graph)
    texts = []
    for name, polymer_type, units in chains:
        body = '.'.join(''.join(_helm_symbol(graph.nodes[i]) if not branch else f'({_helm_symbol(graph.nodes[i])})'
                                for i, branch in unit) for unit in units)
        texts.append(f'{name}{{{body}}}')
    links = '|'.join(f'{ca},{cb},{pa}:R{ra}-{pb}:R{rb}' for (ca, pa, ra), (cb, pb, rb) in connections)
    return '|'.join(texts) + '$' + links + '$$$'


def write_biln(graph):
    for node in graph.nodes:
        if node.inline:
            raise BilnError(f'BILN has no inline SMILES monomer: {node.symbol!r:.60} cannot be written')
        if node.polymer_type != 'PEPTIDE':
            raise BilnError(f'BILN has no {node.polymer_type} chain: {node.symbol!r:.40} cannot be written')
    chains, connections = _layout(graph)
    ends = {}
    for k, ((ca, pa, ra), (cb, pb, rb)) in enumerate(connections, 1):
        ends.setdefault((ca, pa), []).append(f'({k},{ra})')
        ends.setdefault((cb, pb), []).append(f'({k},{rb})')
    texts = []
    for name, _, units in chains:
        tokens = []
        for unit in units:
            for i, _ in unit:
                position = len(tokens) + 1
                tokens.append(_biln_symbol(graph.nodes[i].symbol) + ''.join(ends.get((name, position), ())))
        texts.append('-'.join(tokens))
    return '.'.join(texts)


def _layout(graph):
    """``(chains, connections)``: ``[(name, type, units)]`` with a unit ``[(node, is_branch)]``, and
    connections ``[((chain, position, R), (chain, position, R))]`` lower end first, sorted."""
    nodes = graph.nodes
    following, preceding, branch, linked = {}, {}, {}, set()
    degree = [0] * len(nodes)
    for e in graph.edges:
        degree[e.a] += 1
        degree[e.b] += 1
    for k, e in enumerate(graph.edges):
        if {e.ra, e.rb} == {1, 2} and nodes[e.a].polymer_type == nodes[e.b].polymer_type != 'CHEM':
            source, target = (e.a, e.b) if e.ra == 2 else (e.b, e.a)
            if source not in following and target not in preceding:
                following[source], preceding[target] = target, source
                linked.add(k)
    for k, e in enumerate(graph.edges):
        if k in linked:
            continue
        for host, rh, leaf, rl in ((e.a, e.ra, e.b, e.rb), (e.b, e.rb, e.a, e.ra)):
            if (rh, rl) == (3, 1) and degree[leaf] == 1 and nodes[host].polymer_type in ('RNA', 'DNA') \
                    and nodes[leaf].polymer_type == nodes[host].polymer_type and host not in branch:
                branch[host] = leaf
                linked.add(k)
                break
    leaves = set(branch.values())

    paths, seen = [], set()
    for start in range(len(nodes)):
        if start in seen or start in leaves or start in preceding:
            continue
        path = [start]
        seen.add(start)
        while path[-1] in following and following[path[-1]] not in seen:
            path.append(following[path[-1]])
            seen.add(path[-1])
        paths.append(path)
    closing = []
    for start in range(len(nodes)):
        if start in seen or start in leaves:
            continue
        path = [start]
        seen.add(start)
        while following[path[-1]] not in seen:
            path.append(following[path[-1]])
            seen.add(path[-1])
        closing.append((path[-1], path[0]))
        paths.append(path)
    paths.sort(key=lambda p: p[0])

    chains, where, counters = [], {}, {}
    for path in paths:
        polymer_type = nodes[path[0]].polymer_type
        counters[polymer_type] = counters.get(polymer_type, 0) + 1
        name = f'{polymer_type}{counters[polymer_type]}'
        units, position = [], 0
        for i in path:
            nucleic = polymer_type in ('RNA', 'DNA')
            if not units or not nucleic or i in branch:
                units.append([])
            position += 1
            where[i] = (name, position, len(chains))
            units[-1].append((i, False))
            if i in branch:
                position += 1
                where[branch[i]] = (name, position, len(chains))
                units[-1].append((branch[i], True))
        chains.append((name, polymer_type, units))

    cyclic = {frozenset(c) for c in closing}
    connections = []
    for k, e in enumerate(graph.edges):
        if k in linked and not (frozenset((e.a, e.b)) in cyclic and _closes(e, closing)):
            continue
        ends = sorted(((where[e.a][2], where[e.a][1], e.ra, where[e.a][0]),
                       (where[e.b][2], where[e.b][1], e.rb, where[e.b][0])))
        connections.append(tuple((name, position, label) for _, position, label, name in ends))
    connections.sort(key=lambda c: tuple((_index_of(chains, x[0]), x[1], x[2]) for x in c))
    return chains, connections


def _closes(edge, closing):
    return any({edge.a, edge.b} == {a, b} for a, b in closing)


def _index_of(chains, name):
    return next(i for i, chain in enumerate(chains) if chain[0] == name)


def _helm_symbol(node):
    return node.symbol if len(node.symbol) == 1 and not node.inline else f'[{node.symbol}]'


def _biln_symbol(symbol):
    return f'[{symbol}]' if any(c in symbol for c in '-().[],') else symbol
