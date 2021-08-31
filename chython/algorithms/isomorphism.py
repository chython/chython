# -*- coding: utf-8 -*-
#
#  Copyright 2018-2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from abc import abstractmethod
from collections import defaultdict, deque
from functools import cached_property
from itertools import permutations
from typing import Any, Dict, Iterator, Set, TYPE_CHECKING, Union
from .._functions import lazy_product
from ..periodictable.element import Element, Query


if TYPE_CHECKING:
    from chython.containers.graph import Graph


frequency = {1: 10,  # H
             6: 9,  # C
             8: 8,  # O
             7: 7,  # N
             15: 6, 16: 6,  # P, S
             9: 5,  # F
             17: 4, 35: 4,  # Cl, Br
             53: 3,  # I
             5: 2, 14: 2,  # B, Si
             11: 1, 12: 1, 19: 1, 20: 1}  # Na, Mg, K, Ca


def atom_frequency(x):
    return frequency.get(x.atomic_number, 0)


class Isomorphism:
    __slots__ = ()

    def __lt__(self, other):
        if len(self) >= len(other):
            return False
        return self.is_substructure(other)

    def __le__(self, other):
        return self.is_substructure(other)

    def __gt__(self, other):
        if len(self) <= len(other):
            return False
        return other.is_substructure(self)

    def __ge__(self, other):
        return other.is_substructure(self)

    def __contains__(self: 'Graph', other: Union[Element, Query, str]):
        """
        Atom in Structure test.
        """
        if isinstance(other, str):
            return any(other == x.atomic_symbol for x in self._atoms.values())
        return any(other == x for x in self._atoms.values())

    def is_substructure(self, other, /, *, optimize: bool = True) -> bool:
        """
        Test self is substructure of other

        :param optimize: Morgan weights based automorphism preventing.
        """
        try:
            next(self.get_mapping(other, optimize=optimize, automorphism_filter=False))
        except StopIteration:
            return False
        return True

    def is_equal(self, other, /, *, optimize: bool = True) -> bool:
        """
        Test self is same structure as other

        :param optimize: Morgan weights based automorphism preventing.
        """
        if len(self) != len(other):
            return False
        try:
            next(self.get_mapping(other, optimize=optimize, automorphism_filter=False))
        except StopIteration:
            return False
        return True

    @abstractmethod
    def get_mapping(self, other, /, *, automorphism_filter: bool = True,
                    optimize: bool = True, fallback: bool = False) -> Iterator[Dict[int, int]]:
        """
        Get self to other substructure mapping generator.

        :param automorphism_filter: Skip matches to same atoms.
        :param optimize: Morgan weights based automorphism preventing.
        :param fallback: Try without optimization then nothing matched.
        """
        if optimize:
            g = self.__components_mapping(other, other.atoms_order, automorphism_filter)
            m = next(g, None)
            if m is not None:
                yield m
                yield from g
                return
            elif not fallback:
                return
        yield from self.__components_mapping(other, {n: i for i, n in enumerate(other)}, automorphism_filter)

    def __components_mapping(self, other, o_order, automorphism_filter):
        components, closures = self._compiled_query
        o_atoms = other._atoms
        o_bonds = other._bonds

        seen = set()
        if len(components) == 1:
            for other_component in range(other.connected_components_count):
                candidate = self._isomorphism_candidates(other, 0, other_component)
                if not candidate:
                    continue
                for mapping in _get_mapping(components[0], closures, o_atoms, o_bonds, candidate, o_order):
                    if automorphism_filter:
                        atoms = frozenset(mapping.values())
                        if atoms in seen:
                            continue
                        seen.add(atoms)
                    yield mapping
        else:
            for cs in permutations(range(other.connected_components_count), len(components)):
                mappers = []
                for self_component, component, other_component in zip(range(len(components)), components, cs):
                    candidate = self._isomorphism_candidates(other, self_component, other_component)
                    if not candidate:
                        break
                    mappers.append(_get_mapping(component, closures, o_atoms, o_bonds, candidate, o_order))
                else:
                    for match in lazy_product(*mappers):
                        mapping = match[0].copy()
                        for m in match[1:]:
                            mapping.update(m)
                        if automorphism_filter:
                            atoms = frozenset(mapping.values())
                            if atoms in seen:
                                continue
                            seen.add(atoms)
                        yield mapping

    def _isomorphism_candidates(self, other, self_component: int, other_component: int) -> Set[int]:
        """
        Filter of atoms with possible isomorphism.
        By default do nothing.
        """
        return set(other.connected_components[other_component])

    @cached_property
    def _compiled_query(self: 'Graph'):
        components, closures = _compile_query(self._atoms, self._bonds,
                                              {n: atom_frequency(a) for n, a in self._atoms.items()})
        if self.connected_components_count > 1:
            order = {x: n for n, c in enumerate(self.connected_components) for x in c}
            components.sort(key=lambda x: order[x[0][0]])
        return components, closures

    def is_automorphic(self):
        """
        Test for automorphism symmetry of graph.
        """
        try:
            next(self.get_automorphism_mapping())
        except StopIteration:
            return False
        return True

    def get_automorphism_mapping(self: 'Graph') -> Iterator[Dict[int, int]]:
        """
        Iterator of all possible automorphism mappings.
        """
        return _get_automorphism_mapping(self.atoms_order, self._bonds,
                                         {n: atom_frequency(a) for n, a in self._atoms.items()})


def _get_automorphism_mapping(atoms: Dict[int, int], bonds: Dict[int, Dict[int, Any]],
                              atoms_frequencies: Dict[int, int]) -> Iterator[Dict[int, int]]:
    if len(atoms) == len(set(atoms.values())):
        return  # all atoms unique

    components, closures = _compile_query(atoms, bonds, atoms_frequencies)
    groups = {x: n for n, x in enumerate(atoms)}
    mappers = [_get_mapping(order, closures, atoms, bonds, {x for x, *_ in order}, groups)
               for order in components]
    if len(mappers) == 1:
        for mapping in mappers[0]:
            if any(k != v for k, v in mapping.items()):
                yield mapping
    for match in lazy_product(*mappers):
        mapping = match[0].copy()
        for m in match[1:]:
            mapping.update(m)
        if any(k != v for k, v in mapping.items()):
            yield mapping


def _get_mapping(linear_query, query_closures, o_atoms, o_bonds, scope, groups):
    size = len(linear_query) - 1
    order_depth = {v[0]: k for k, v in enumerate(linear_query)}

    stack = deque()
    path = []
    mapping = {}
    reversed_mapping = {}

    s_n, _, s_atom, _ = linear_query[0]
    for n, o_atom in o_atoms.items():
        if n in scope and s_atom == o_atom:
            stack.append((n, 0))

    while stack:
        n, depth = stack.pop()
        current = linear_query[depth][0]
        if depth == size:
            yield {current: n, **mapping}
        else:
            if len(path) != depth:
                for x in path[depth:]:
                    del mapping[reversed_mapping.pop(x)]
                path = path[:depth]

            path.append(n)
            mapping[current] = n
            reversed_mapping[n] = current

            depth += 1
            s_n, back, s_atom, s_bond = linear_query[depth]
            if back != current:
                n = path[order_depth[back]]

            uniq = set()
            for o_n, o_bond in o_bonds[n].items():
                if o_n in scope and o_n not in reversed_mapping and s_bond == o_bond and groups[o_n] not in uniq:
                    uniq.add(groups[o_n])
                    if s_atom == o_atoms[o_n]:
                        # check closures equality
                        o_closures = o_bonds[o_n].keys() & reversed_mapping.keys()
                        o_closures.discard(n)
                        if o_closures == {mapping[m] for m, _ in query_closures[s_n]}:
                            obon = o_bonds[o_n]
                            if all(bond == obon[mapping[m]] for m, bond in query_closures[s_n]):
                                stack.append((o_n, depth))


def _compile_query(atoms, bonds, atoms_frequencies):
    closures = defaultdict(list)
    components = []
    seen = set()
    while len(seen) < len(atoms):
        start = min(atoms.keys() - seen, key=lambda x: atoms_frequencies[x])
        seen.add(start)
        stack = [(n, start, atoms[n], bond) for n, bond in sorted(bonds[start].items(), reverse=True,
                                                                  key=lambda x: atoms_frequencies[x[0]])]
        order = [(start, None, atoms[start], None)]
        components.append(order)

        while stack:
            front, back, *_ = atom = stack.pop()
            if front not in seen:
                order.append(atom)
                for n, bond in sorted(bonds[front].items(), reverse=True, key=lambda x: atoms_frequencies[x[0]]):
                    if n != back:
                        if n not in seen:
                            stack.append((n, front, atoms[n], bond))
                        else:
                            closures[front].append((n, bond))
                seen.add(front)
    return components, closures


__all__ = ['Isomorphism']
