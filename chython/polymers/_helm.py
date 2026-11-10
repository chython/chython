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
"""HELM text into the monomer graph.

``CHAINS$CONNECTIONS$GROUPS$ANNOTATIONS$VERSION``; chains split on ``|``, units on ``.``, a unit being
backbone monomers with ``(...)`` branches -- ``R(A)P`` is sugar, base, phosphate.  ``[...]`` nests, and
holds a multi-letter symbol or an inline SMILES (one with a ``*`` attachment).

| Written | Read as |
|---|---|
| ``A'3'`` | three units |
| ``A'2-5'``, ``(A,G)``, ``(A+G)``, ``?``, ``*``, a BLOB chain | one opaque node, ``polymers:opaque-unit`` |
| a connection at ``pair`` | no bond, ``polymers:hydrogen-pair-not-built`` |
| a connection at ``?`` | no bond, ``polymers:attachment-unknown`` |
| a ``"..."`` annotation or an annotation section | not kept, ``polymers:annotation-not-kept`` |
| a group section | not kept, ``polymers:group-not-kept`` |
"""

from re import compile as re_compile

from ._errors import HelmError
from ._graph import Edge, Node, PolymerGraph, POLYMER_TYPES
from ..core import LogRecord, LOST


__all__ = ['parse_helm']


_CHAIN = re_compile(r'^(PEPTIDE|RNA|DNA|CHEM|BLOB)([1-9][0-9]*)$')
_ATTACH = re_compile(r'^([0-9]+|\?):(R[0-9]+|pair|\?)$')
_OPENING = {'[': ']', '{': '}', '(': ')'}


def parse_helm(text):
    """The monomer graph of a HELM string.  Raises :class:`HelmError` for text that does not parse."""
    if not isinstance(text, str) or not text.strip():
        raise HelmError('an empty string is not HELM')
    sections = _split(text.strip(), '$')
    graph = PolymerGraph()
    chains = {}            # chain name -> node indices in order
    for part in _split(sections[0], '|'):
        _chain(part, graph, chains)
    if not chains:
        raise HelmError('HELM names no polymer')
    if len(sections) > 1 and sections[1]:
        for part in _split(sections[1], '|'):
            _connection(part, graph, chains)
    if len(sections) > 2 and sections[2]:
        graph.log.append(LogRecord('polymers:group-not-kept', (),
                                   f'unsupported: polymer group section {sections[2]!r:.60} not kept', LOST))
    annotation = sections[3] if len(sections) > 3 else ''
    if len(sections) == 4 and annotation.upper().startswith('V'):
        annotation = ''    # a four-section HELM 1 string ending in its version
    if annotation:
        graph.log.append(LogRecord('polymers:annotation-not-kept', (),
                                   f'unsupported: annotation section {annotation!r:.60} not kept', LOST))
    return graph


def _split(text, separator):
    """`text` split on `separator` outside brackets, braces, parentheses and double quotes."""
    out, depth, quoted, start = [], [], False, 0
    for i, char in enumerate(text):
        if quoted:
            quoted = char != '"'
        elif char == '"':
            quoted = True
        elif char in _OPENING:
            depth.append(_OPENING[char])
        elif depth and char == depth[-1]:
            depth.pop()
        elif char in ')]}':
            raise HelmError(f'unbalanced {char!r} at {i} in {text!r:.80}')
        elif char == separator and not depth:
            out.append(text[start:i])
            start = i + 1
    if depth or quoted:
        raise HelmError(f'unclosed {"quote" if quoted else depth[-1]!r} in {text!r:.80}')
    out.append(text[start:])
    return out


def _chain(part, graph, chains):
    head, brace, rest = part.partition('{')
    match = _CHAIN.match(head.strip())
    if not brace or not match or '}' not in rest:
        raise HelmError(f'{part!r:.60} is not a polymer: TYPE<n>{{...}} with TYPE one of {", ".join(POLYMER_TYPES)}')
    body, _, tail = rest.rpartition('}')
    name, polymer_type = head.strip(), match.group(1)
    if name in chains:
        raise HelmError(f'polymer {name} is defined twice')
    if tail.strip():
        graph.log.append(LogRecord('polymers:annotation-not-kept', (),
                                   f'unsupported: annotation {tail.strip()!r:.60} on {name} not kept', LOST))
    indices = chains[name] = []
    if polymer_type == 'BLOB':
        _add(graph, indices, name, polymer_type, body, opaque=True)
        return
    previous = None
    for unit in _split(body, '.'):
        if not unit:
            raise HelmError(f'{name}: an empty unit in {body!r:.60}')
        for first, last in _unit(unit, graph, indices, name, polymer_type):
            if previous is not None:
                graph.edges.append(Edge(previous, 2, first, 1))
            previous = last


def _unit(unit, graph, indices, name, polymer_type):
    """Add one unit's nodes; yield ``(first backbone, last backbone)`` per copy of it."""
    body, count = unit, 1
    if unit.endswith("'") and unit.count("'") >= 2:
        body, _, repeat = unit[:-1].rpartition("'")
        if repeat.isdigit() and int(repeat) > 0:
            count = int(repeat)
        else:
            yield _opaque(graph, indices, name, polymer_type, unit)
            return
    tokens = _tokens(body, name)
    if tokens is None:
        yield _opaque(graph, indices, name, polymer_type, unit)
        return
    for _ in range(count):
        backbone = None
        first = None
        for symbol, branch in tokens:
            index = _add(graph, indices, name, polymer_type, symbol)
            if branch:
                if backbone is None:
                    raise HelmError(f'{name}: a branch with no monomer before it in {unit!r:.60}')
                graph.edges.append(Edge(backbone, 3, index, 1))
                continue
            if backbone is not None:
                graph.edges.append(Edge(backbone, 2, index, 1))
            backbone = index
            if first is None:
                first = index
        yield first, backbone


def _tokens(body, name):
    """``[(symbol, is_branch)]`` of a unit, or ``None`` when it names no one monomer."""
    out, i = [], 0
    while i < len(body):
        char = body[i]
        if char == '[':
            end = _close(body, i, '[', ']', name)
            out.append((body[i + 1:end], False))
            i = end + 1
        elif char == '(':
            end = _close(body, i, '(', ')', name)
            inner = body[i + 1:end]
            if not inner or any(c in inner for c in ',+') or inner in '?*_':
                return None
            if inner.startswith('['):
                if not inner.endswith(']'):
                    return None
                inner = inner[1:-1]
            out.append((inner, True))
            i = end + 1
        elif char == '"':
            end = body.find('"', i + 1)
            i = len(body) if end < 0 else end + 1
        elif char in '?*_':
            return None
        else:
            out.append((char, False))
            i += 1
    return out or None


def _close(text, start, opening, closing, name):
    depth = 0
    for i in range(start, len(text)):
        if text[i] == opening:
            depth += 1
        elif text[i] == closing:
            depth -= 1
            if not depth:
                return i
    raise HelmError(f'{name}: unclosed {opening!r} in {text!r:.60}')


def _add(graph, indices, name, polymer_type, symbol, opaque=False):
    index = len(graph.nodes)
    graph.nodes.append(Node(name, len(indices) + 1, polymer_type, symbol, not opaque and '*' in symbol))
    indices.append(index)
    if opaque:
        graph.opaque.append((index, symbol))
        graph.log.append(LogRecord('polymers:opaque-unit', (),
                                   f'unsupported: {name} unit {symbol!r:.60} names no one monomer; kept as a '
                                   f'marker', LOST))
    return index


def _opaque(graph, indices, name, polymer_type, text):
    index = _add(graph, indices, name, polymer_type, text, opaque=True)
    return index, index


def _connection(part, graph, chains):
    text = part.strip()
    if '"' in text:
        text, _, note = text.partition('"')
        graph.log.append(LogRecord('polymers:annotation-not-kept', (),
                                   f'unsupported: annotation "{note.rstrip(chr(34))!s:.60}" on a connection not '
                                   f'kept', LOST))
        text = text.strip()
    fields = text.split(',')
    if len(fields) != 3 or '-' not in fields[2]:
        raise HelmError(f'{part!r:.60} is not a connection: A,B,i:Rn-j:Rm')
    left, _, right = fields[2].partition('-')
    ends = []
    for chain, attach in ((fields[0].strip(), left.strip()), (fields[1].strip(), right.strip())):
        if chain not in chains:
            raise HelmError(f'connection {part!r:.60} names {chain}, which no polymer defines')
        match = _ATTACH.match(attach)
        if not match:
            raise HelmError(f'connection {part!r:.60}: {attach!r} is not position:R<n>')
        position, label = match.groups()
        if position != '?' and not 1 <= int(position) <= len(chains[chain]):
            raise HelmError(f'connection {part!r:.60}: {chain} has no position {position}')
        ends.append((chain, position, label))
    if any(label == 'pair' for _, _, label in ends):
        graph.log.append(LogRecord('polymers:hydrogen-pair-not-built', (),
                                   f'unsupported: hydrogen pair {text!r:.60} is not a bond and is not built', LOST))
        return
    if any('?' in (position, label) for _, position, label in ends):
        graph.log.append(LogRecord('polymers:attachment-unknown', (),
                                   f'unsupported: connection {text!r:.60} names no definite attachment; not built',
                                   LOST))
        return
    (ca, pa, la), (cb, pb, lb) = ends
    graph.edges.append(Edge(chains[ca][int(pa) - 1], int(la[1:]), chains[cb][int(pb) - 1], int(lb[1:])))
