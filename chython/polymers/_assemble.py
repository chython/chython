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
"""The monomer graph into a molecule, through ONE SMILES parse.

Each node writes its monomer's fragment (:meth:`MonomerLibrary.fragment`) with every used R left as an
open ring closure; the fragments join with ``.`` and each edge is one closure id shared by its two ends.
Ids are 10..99, so they are handed out like registers: an id is busy from its first occurrence in the
joined text to its second, then free again -- what bounds the count is the bonds open at once, not the
bonds in the molecule.  The SMILES parser places every bond and keeps every configuration as written; an
unused R kept its cap, and the caps' hydrogens go implicit after the parse.

| Node | Writes |
|---|---|
| library monomer | its fragment |
| inline SMILES | cut like a monomer, ``[*:n]`` its leaving atoms; an unused one capped ``[H]`` (logged) |
| unknown symbol, opaque unit | the marker ``*`` with one closure per edge, the symbol as its alias (logged) |
"""

from re import compile as re_compile

from ._errors import HelmError
from ._library import Monomer
from ..core import LogRecord, LOST, REPAIRED, smiles


__all__ = ['assemble']


_INLINE_R = re_compile(r'\[\*:([0-9]+)\]')
_TOKEN = re_compile(r'\{R([0-9]+)\}')
_IDS = range(10, 100)


def assemble(graph, library, log):
    """The molecule of `graph`, its monomers taken from `library`."""
    log.extend(graph.log)
    opaque = {index for index, _ in graph.opaque}
    incident = [{} for _ in graph.nodes]          # node -> {R number: edge index}
    for k, edge in enumerate(graph.edges):
        for node, label in ((edge.a, edge.ra), (edge.b, edge.rb)):
            if label in incident[node]:
                raise HelmError(f'{_name(graph, node)}: R{label} carries two bonds')
            incident[node][label] = k

    templates, reserved, markers = [], set(), []
    for index, node in enumerate(graph.nodes):
        used = incident[index]
        if index in opaque or (not node.inline and library.get(node.polymer_type, node.symbol) is None):
            if index not in opaque:
                log.append(LogRecord('polymers:monomer-unknown', (),
                                     f'{_name(graph, index)}: monomer {node.symbol!r:.40} is not in the library; '
                                     f'written as a marker carrying its symbol', LOST))
            templates.append('*' + ''.join(f'{{R{label}}}' for label in sorted(used)))
            markers.append(node.symbol)
        elif node.inline:
            fragment = library.fragment(_inline(graph, index, node, used, log), frozenset(used))
            reserved.update(fragment.closures)
            templates.append(fragment.text)
        else:
            monomer = library.get(node.polymer_type, node.symbol)
            missing = [label for label in used if label not in monomer.labels]
            if missing:
                raise HelmError(f'{_name(graph, index)}: monomer {node.symbol} has no '
                                f'{", ".join(f"R{x}" for x in missing)}')
            fragment = library.fragment(monomer, frozenset(used))
            reserved.update(fragment.closures)
            templates.append(fragment.text)

    free = [i for i in _IDS if i not in reserved]
    assigned = {}                                   # edge index -> closure id
    texts = []
    for index, template in enumerate(templates):
        closing = []

        def token(match, index=index, closing=closing):
            edge = incident[index][int(match.group(1))]
            number = assigned.get(edge)
            if number is None:
                if not free:
                    raise HelmError('more than 90 bonds are open at one point of the chain')
                number = assigned[edge] = free.pop(0)
            else:
                closing.append(number)
            return f'%{number}'

        texts.append(_TOKEN.sub(token, template))
        free = sorted(free + closing)

    mol = smiles('.'.join(texts))
    if markers:
        atoms = [n for n in mol.atom_numbers if mol.atom(n).is_r]
        aliases = dict(mol.aliases)
        aliases.update(zip(atoms, markers))
        mol.set_aliases(aliases)
    mol.implicify_hydrogens()
    return mol


def _inline(graph, index, node, used, log):
    """An inline SMILES as a :class:`Monomer` whose ``[*:n]`` atoms are the leaving atoms.

    An R the graph leaves unused is capped with hydrogen, logged; the SMILES is the cache key, so one
    inline monomer written twice is cut once.
    """
    text = node.symbol
    labels = sorted({int(x) for x in _INLINE_R.findall(text)})
    missing = [label for label in used if label not in labels]
    if missing:
        raise HelmError(f'{_name(graph, index)}: inline monomer {text!r:.40} has no '
                        f'{", ".join(f"R{x}" for x in missing)}')
    unused = [label for label in labels if label not in used]
    if unused:
        log.append(LogRecord('polymers:inline-cap', (),
                             f'{_name(graph, index)}: inline monomer {text!r:.40}: unused '
                             f'{", ".join(f"R{x}" for x in unused)} capped with hydrogen', REPAIRED))
        text = _INLINE_R.sub(lambda m: m.group(0) if int(m.group(1)) in used else '[H]', text)
    return Monomer(text, node.polymer_type, '', text, tuple(sorted(used)), '')


def _name(graph, index):
    node = graph.nodes[index]
    return f'{node.chain} position {node.position}'
