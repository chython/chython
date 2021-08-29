# -*- coding: utf-8 -*-
#
#  Copyright 2020, 2021 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2020 Ravil Mukhametgaleev <sonic-mc@mail.ru>
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
from collections import defaultdict, deque
from functools import cached_property
from numpy import uint, zeros
from typing import List, Tuple, Dict, Set, Any, Union, Iterable, TYPE_CHECKING
from ...containers.bonds import Bond, QueryBond, DynamicBond


if TYPE_CHECKING:
    from chython.containers.graph import Graph


class GraphComponents:
    __slots__ = ()

    @cached_property
    def connected_components(self: 'Graph') -> Tuple[Tuple[int, ...], ...]:
        """
        Isolated components of single graph. E.g. salts as ion pair.
        """
        if not self._atoms:
            return ()
        return tuple(tuple(x) for x in self._connected_components(self._bonds))

    @staticmethod
    def _connected_components(bonds: Dict[int, Union[Set[int], Dict[int, Any]]]) -> List[Set[int]]:
        atoms = set(bonds)
        components = []
        while atoms:
            start = atoms.pop()
            seen = {start}
            queue = deque([start])
            while queue:
                current = queue.popleft()
                for i in bonds[current]:
                    if i not in seen:
                        queue.append(i)
                        seen.add(i)
            components.append(seen)
            atoms.difference_update(seen)
        return components

    @property
    def connected_components_count(self) -> int:
        """
        Number of components in graph
        """
        return len(self.connected_components)

    @cached_property
    def skin_atoms(self: 'Graph') -> Tuple[int, ...]:
        """
        Atoms of skin graph.
        """
        return tuple(self._skin_graph(self._bonds))

    @cached_property
    def skin_graph(self: 'Graph') -> Dict[int, Set[int]]:
        """
        Graph without terminal atoms. Only rings and linkers
        """
        return self._skin_graph(self._bonds)

    @staticmethod
    def _skin_graph(bonds: Dict[int, Union[Set[int], Dict[int, Any]]]) -> Dict[int, Set[int]]:
        """
        Graph without terminal nodes. Only rings and linkers
        """
        bonds = {n: set(ms) for n, ms in bonds.items() if ms}
        while True:  # skip not-cycle chains
            try:
                n = next(n for n, ms in bonds.items() if len(ms) <= 1)
            except StopIteration:
                break
            for m in bonds.pop(n):
                bonds[m].discard(n)
        return bonds

    def adjacency_matrix(self: 'Graph', set_bonds=False):
        """
        Adjacency matrix of Graph.

        :param set_bonds: if True set bond orders instead of 1.
        """
        adj = zeros((len(self), len(self)), dtype=uint)
        mapping = {n: x for x, n in enumerate(self._atoms)}
        if set_bonds:
            for n, ms in self._bonds.items():
                n = mapping[n]
                for m, b in ms.items():
                    adj[n, mapping[m]] = int(b)
        else:
            for n, ms in self._bonds.items():
                n = mapping[n]
                for m, b in ms.items():
                    adj[n, mapping[m]] = 1
        return adj

    @cached_property
    def ring_atoms(self):
        """
        Atoms in rings
        """
        bonds = self._skin_graph(self.not_special_connectivity)
        if not bonds:
            return frozenset()

        in_rings = set()
        atoms = set(bonds)
        while atoms:
            stack = deque([(atoms.pop(), 0, 0)])
            path = []
            seen = set()
            while stack:
                c, p, d = stack.pop()
                if len(path) > d:
                    path = path[:d]
                if c in in_rings:
                    continue
                path.append(c)
                seen.add(c)

                d += 1
                for n in bonds[c]:
                    if n == p:
                        continue
                    elif n in seen:
                        in_rings.update(path[path.index(n):])
                    else:
                        stack.append((n, c, d))

            atoms.difference_update(seen)
        return in_rings

    @cached_property
    def rings_count(self) -> int:
        """
        SSSR rings count. Ignored rings with special bonds.
        """
        bonds = self.not_special_connectivity
        return sum(len(x) for x in bonds.values()) // 2 - len(bonds) + len(self._connected_components(bonds))

    @cached_property
    def not_special_bonds(self: 'Graph') -> Dict[int, Dict[int, Union[Bond, QueryBond, DynamicBond]]]:
        """
        Bonds without special.
        """
        bonds = {}
        for n, ms in self._bonds.items():
            bonds[n] = ngb = {}
            for m, b in ms.items():
                if b != 8:
                    ngb[m] = b
        return bonds

    @cached_property
    def not_special_connectivity(self: 'Graph') -> Dict[int, Set[int]]:
        """
        Graph connectivity without special bonds.
        """
        bonds = {}
        for n, ms in self._bonds.items():
            bonds[n] = ngb = set()
            for m, b in ms.items():
                if b != 8:
                    ngb.add(m)
        return bonds

    @cached_property
    def atoms_rings(self: 'Graph') -> Dict[int, Tuple[Tuple[int, ...]]]:
        """
        Dict of atoms rings which contains it.
        """
        rings = defaultdict(list)
        for r in self.sssr:
            for n in r:
                rings[n].append(r)
        return {n: tuple(rs) for n, rs in rings.items()}

    @cached_property
    def atoms_rings_sizes(self) -> Dict[int, Tuple[int, ...]]:
        """
        Sizes of rings containing atom.
        """
        return {n: tuple(len(r) for r in rs) for n, rs in self.atoms_rings.items()}

    @cached_property
    def int_adjacency(self: 'Graph') -> Dict[int, Dict[int, int]]:
        """
        Adjacency with integer-coded bonds.
        """
        return {n: {m: int(b) for m, b in mb.items()} for n, mb in self._bonds.items()}

    def _augmented_substructure(self: 'Graph', atoms: Iterable[int], deep: int):
        atoms = set(atoms)
        bonds = self._bonds
        if atoms - self._atoms.keys():
            raise ValueError('invalid atom numbers')
        nodes = [atoms]
        for i in range(deep):
            n = {y for x in nodes[-1] for y in bonds[x]} | nodes[-1]
            if n in nodes:
                break
            nodes.append(n)
        return nodes


__all__ = ['GraphComponents']
