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
from collections import defaultdict
from itertools import product, combinations, islice
from typing import Dict, Set, Iterator, Tuple
from ..containers import molecule


class MCS:
    __slots__ = ()

    def get_mcs_mapping(self, other: 'molecule.MoleculeContainer', /, *, limit=10000) -> Iterator[Dict[int, int]]:
        """
        Find maximum common substructure. Based on clique searching in product graph.

        :param limit: limit tested cliques
        """
        if not isinstance(other, molecule.MoleculeContainer):
            raise TypeError('MoleculeContainer expected')

        core_product, full_product = self.__get_product(other)
        if not core_product:
            return

        # search maximum bonded substructures
        hits = []
        max_atoms = 0
        max_bonds = 0
        for mapping in islice(_clique(full_product), limit):
            if len(mapping) < max_atoms:
                continue
            # search bonds count
            bonds = 0
            seen = set()
            for n in mapping:
                seen.add(n)
                for m in core_product[n]:
                    if m not in seen and m in mapping:
                        bonds += 1
            if bonds > max_bonds:
                max_bonds = bonds
                max_atoms = len(mapping) - 1  # -1 is ad-hoc
                hits = [mapping]
            elif bonds == max_bonds:
                hits.append(mapping)

        # search maximal components in substructures
        hits2 = []
        max_component = 0
        for mapping in hits:
            # search components
            components = []
            atoms = mapping.copy()
            while atoms:
                n = atoms.pop()
                seen = {n}
                queue = [n]
                component = []
                while queue:
                    n = queue.pop(0)
                    component.append(n)
                    for m in core_product[n]:
                        if m not in seen and m in mapping:
                            queue.append(m)
                            seen.add(m)

                components.append(component)
                atoms.difference_update(component)

            # get max component
            component = max(len(x) for x in components)
            if component > max_component:
                max_component = component
                hits2 = [mapping]
            elif component == max_component:
                hits2.append(mapping)
        yield from (dict(x) for x in hits2)

    def __get_product(self: 'molecule.MoleculeContainer', other: 'molecule.MoleculeContainer'):
        bonds = self._bonds
        o_bonds = other._bonds

        s_equal = defaultdict(list)  # equal self atoms
        for n, atom in self.atoms():
            s_equal[atom].append(n)
        p_equal = defaultdict(list)  # equal other atoms
        for n, atom in other.atoms():
            p_equal[atom].append(n)

        full_product = {}
        core_product = {}
        equal_atoms = {}
        for atom, ns in s_equal.items():
            ms = p_equal[atom]
            if ms:
                for nm in product(ns, ms):
                    full_product[nm] = set()
                    core_product[nm] = set()
                for n in ns:
                    equal_atoms[n] = ms  # memory save

        seen = set()
        for n, o_ns in equal_atoms.items():
            seen.add(n)
            for m, b in bonds[n].items():
                if m in equal_atoms and m not in seen:
                    o_ms = equal_atoms[m]
                    for o_n in o_ns:
                        node1 = (n, o_n)
                        fms = full_product[node1]
                        cms = core_product[node1]
                        for o_m, o_b in o_bonds[o_n].items():
                            if o_m in o_ms and b == o_b:
                                node2 = (m, o_m)
                                full_product[node2].add(node1)
                                core_product[node2].add(node1)
                                fms.add(node2)
                                cms.add(node2)

        atoms = core_product
        while atoms:
            new_atoms = set()
            for n in atoms:
                core = core_product[n]
                for nm1, nm2 in combinations(full_product[n], 2):
                    n1, m1 = nm1
                    n2, m2 = nm2
                    if n1 == n2 or m1 == m2:
                        continue
                    if nm1 in full_product[nm2]:
                        continue
                    if nm1 not in core and nm2 not in core:
                        continue

                    full_product[nm1].add(nm2)
                    full_product[nm2].add(nm1)
                    new_atoms.add(nm1)
                    new_atoms.add(nm2)
            atoms = new_atoms

        return core_product, full_product


def _clique(graph) -> Iterator[Set[Tuple[int, int]]]:
    """
    clique search

    adopted from networkx algorithms.clique.find_cliques
    """
    subgraph = {x for x, y in graph.items() if y}  # skip isolated nodes
    if not subgraph:
        return  # empty or fully disconnected
    elif len(subgraph) == 2:  # dimer
        yield set(subgraph)
        return

    stack = []
    clique_atoms = [None]
    candidates = subgraph.copy()
    roots = candidates - graph[max(subgraph, key=lambda x: len(graph[x]))]

    while True:
        if roots:
            root = roots.pop()
            candidates.remove(root)
            clique_atoms[-1] = root
            neighbors = graph[root]
            neighbors_subgraph = subgraph & neighbors
            if not neighbors_subgraph:
                yield set(clique_atoms)
            else:
                neighbors_candidates = candidates & neighbors
                if neighbors_candidates:
                    stack.append((subgraph, candidates, roots))
                    clique_atoms.append(None)
                    subgraph = neighbors_subgraph
                    candidates = neighbors_candidates
                    roots = candidates - graph[max(subgraph, key=lambda x: len(candidates & graph[x]))]
        elif not stack:
            return
        else:
            clique_atoms.pop()
            subgraph, candidates, roots = stack.pop()


__all__ = ['MCS']
