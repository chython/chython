# -*- coding: utf-8 -*-
#
#  Copyright 2017-2025 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from typing import Any, Dict, List, Set, Tuple, Union, TYPE_CHECKING
from ._rings import sssr


if TYPE_CHECKING:
    from chython.containers import MoleculeContainer


class Rings:
    __slots__ = ()

    @cached_property
    def sssr(self) -> List[Tuple[int, ...]]:
        """
        Smallest Set of Smallest Rings. Special bonds ignored.

        Based on idea of PID matrices from:
        Lee, C. J., Kang, Y.-M., Cho, K.-H., & No, K. T. (2009).
        A robust method for searching the smallest set of smallest rings with a path-included distance matrix.
        Proceedings of the National Academy of Sciences of the United States of America, 106(41), 17355–17358.
        https://doi.org/10.1073/pnas.0813040106

        :return rings atoms numbers
        """
        if self.rings_count:
            return sssr(_skin_graph(self.not_special_connectivity), self.rings_count)
        return []

    @cached_property
    def atoms_rings(self) -> Dict[int, List[Tuple[int, ...]]]:
        """
        A dictionary with atom numbers as keys and a list of tuples (representing SSSR rings) as values.
        """
        rings = defaultdict(list)
        for r in self.sssr:
            for n in r:
                rings[n].append(r)
        return dict(rings)

    @cached_property
    def atoms_rings_sizes(self) -> Dict[int, Set[int]]:
        """
        Sizes of SSSR rings containing atom.
        """
        return {n: {len(r) for r in rs} for n, rs in self.atoms_rings.items()}

    @cached_property
    def rings_count(self) -> int:
        """
        SSSR rings count. Ignored rings with special bonds.
        """
        bonds = self.not_special_connectivity
        return sum(len(x) for x in bonds.values()) // 2 - len(bonds) + len(_connected_components(bonds))

    @cached_property
    def not_special_connectivity(self: 'MoleculeContainer') -> Dict[int, Set[int]]:
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
    def connected_components(self: 'MoleculeContainer') -> List[Set[int]]:
        """
        Isolated components of single graph. E.g. salts as ion pair.
        """
        return _connected_components(self._bonds)

    @property
    def connected_components_count(self) -> int:
        """
        Number of components in graph
        """
        return len(self.connected_components)

    @cached_property
    def skin_graph(self: 'MoleculeContainer') -> Dict[int, Set[int]]:
        """
        Graph without terminal atoms. Only rings and linkers
        """
        return _skin_graph(self._bonds)

    @cached_property
    def rings_graph(self: 'MoleculeContainer'):
        """
        Graph of rings. Linkers are not included. Special bonds are considered.
        """
        bonds = {n: ms.copy() for n, ms in self.skin_graph.items()}
        if not bonds:
            return {}

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
        for n in bonds.keys() - in_rings:
            for m in bonds.pop(n):
                bonds[m].discard(n)
        return bonds


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


__all__ = ['Rings']
