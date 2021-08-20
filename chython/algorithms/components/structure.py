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
from CachedMethods import cached_property
from collections import defaultdict
from typing import Tuple


_heteroatoms = {5, 6, 7, 8, 14, 15, 16, 17, 33, 34, 35, 52, 53}


class StructureComponents:
    __slots__ = ()

    @cached_property
    def aromatic_rings(self) -> Tuple[Tuple[int, ...], ...]:
        """
        Aromatic rings atoms numbers
        """
        bonds = self._bonds
        return tuple(ring for ring in self.sssr if bonds[ring[0]][ring[-1]] == 4
                     and all(bonds[n][m] == 4 for n, m in zip(ring, ring[1:])))

    @cached_property
    def cumulenes(self) -> Tuple[Tuple[int, ...], ...]:
        """
        Alkenes, allenes and cumulenes atoms numbers
        """
        return self._cumulenes()

    @cached_property
    def tetrahedrons(self) -> Tuple[int, ...]:
        """
        Carbon sp3 atoms numbers
        """
        atoms = self._atoms
        bonds = self._bonds
        charges = self._charges
        radicals = self._radicals

        tetra = []
        for n, atom in atoms.items():
            if atom.atomic_number == 6 and not charges[n] and not radicals[n]:
                env = bonds[n]
                if all(x == 1 for x in env.values()):
                    if sum(int(x) for x in env.values()) > 4:
                        continue
                    tetra.append(n)
        return tetra

    def _cumulenes(self, heteroatoms=False):
        atoms = self._atoms
        bonds = self._bonds

        adj = defaultdict(set)  # double bonds adjacency matrix
        if heteroatoms:
            for n, atom in atoms.items():
                if atom.atomic_number in _heteroatoms:
                    adj_n = adj[n].add
                    for m, bond in bonds[n].items():
                        if bond == 2 and atoms[m].atomic_number in _heteroatoms:
                            adj_n(m)
        else:
            for n, atom in atoms.items():
                if atom.atomic_number == 6:
                    adj_n = adj[n].add
                    for m, bond in bonds[n].items():
                        if bond == 2 and atoms[m].atomic_number == 6:
                            adj_n(m)
        if not adj:
            return ()

        terminals = [x for x, y in adj.items() if len(y) == 1]
        cumulenes = []
        while terminals:
            n = terminals.pop(0)
            m = adj[n].pop()
            path = [n, m]
            while m not in terminals:
                adj_m = adj[m]
                if len(adj_m) > 2:  # not cumulene. SO3 etc.
                    cumulenes.extend(zip(path, path[1:]))  # keep single double bonds.
                    break
                adj_m.discard(n)
                n, m = m, adj_m.pop()
                path.append(m)
            else:
                terminals.remove(m)
                adj[m].pop()
                cumulenes.append(tuple(path))
        return cumulenes


__all__ = ['StructureComponents']
