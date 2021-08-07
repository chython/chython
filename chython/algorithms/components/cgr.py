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


class CGRComponents:
    __slots__ = ()

    @cached_property
    def centers_list(self) -> Tuple[Tuple[int, ...], ...]:
        """ Get a list of lists of atoms of reaction centers
        """
        radicals = self._radicals
        p_charges = self._p_charges
        p_radicals = self._p_radicals
        center = set()
        adj = defaultdict(set)

        for n, c in self._charges.items():
            if c != p_charges[n] or radicals[n] != p_radicals[n]:
                center.add(n)

        for n, m_bond in self._bonds.items():
            for m, bond in m_bond.items():
                if bond.order != bond.p_order:
                    adj[n].add(m)
        center.update(adj)

        # changes in condensed aromatic rings.
        if self.aromatic_rings:
            adj_update = defaultdict(set)
            for r in self.aromatic_rings:
                if not center.isdisjoint(r):
                    n = r[0]
                    m = r[-1]
                    if n in adj and m in adj[n]:
                        for n, m in zip(r, r[1:]):
                            if m not in adj[n]:
                                adj_update[n].add(m)
                                adj_update[m].add(n)
                    elif any(n in adj and m in adj[n] for n, m in zip(r, r[1:])):
                        adj_update[n].add(m)
                        adj_update[m].add(n)
                        for n, m in zip(r, r[1:]):
                            if m not in adj[n]:
                                adj_update[n].add(m)
                                adj_update[m].add(n)
            for n, ms in adj_update.items():
                adj[n].update(ms)

        out = []
        while center:
            n = center.pop()
            if n in adj:
                seen = set()
                nextlevel = {n}
                while nextlevel:
                    thislevel = nextlevel
                    nextlevel = set()
                    for v in thislevel:
                        if v not in seen:
                            seen.add(v)
                            nextlevel.update(adj[v])
                out.append(tuple(seen))
                center.difference_update(seen)
            else:
                out.append((n,))
        return tuple(out)

    @cached_property
    def center_atoms(self) -> Tuple[int, ...]:
        """ Get list of atoms of reaction center (atoms with dynamic: bonds, charges, radicals).
        """
        radicals = self._radicals
        p_charges = self._p_charges
        p_radicals = self._p_radicals

        center = set()
        for n, c in self._charges.items():
            if c != p_charges[n] or radicals[n] != p_radicals[n]:
                center.add(n)

        for n, m_bond in self._bonds.items():
            if any(bond.order != bond.p_order for bond in m_bond.values()):
                center.add(n)

        return tuple(center)

    @cached_property
    def center_bonds(self) -> Tuple[Tuple[int, int], ...]:
        """ Get list of bonds of reaction center (bonds with dynamic orders).
        """
        return tuple((n, m) for n, m, bond in self.bonds() if bond.order != bond.p_order)

    @cached_property
    def aromatic_rings(self) -> Tuple[Tuple[int, ...], ...]:
        """
        existed or formed aromatic rings atoms numbers
        """
        adj = self._bonds
        return tuple(ring for ring in self.sssr if
                     adj[ring[0]][ring[-1]].order == 4 and all(adj[n][m].order == 4 for n, m in zip(ring, ring[1:])) or
                     adj[ring[0]][ring[-1]].p_order == 4 and all(adj[n][m].p_order == 4 for n, m in zip(ring, ring[1:]))
                     )


__all__ = ['CGRComponents']
