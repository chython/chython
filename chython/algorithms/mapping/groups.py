# -*- coding: utf-8 -*-
#
#  Copyright 2022 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from itertools import chain, repeat
from typing import List, Tuple, TYPE_CHECKING, Union
from ._groups import *


if TYPE_CHECKING:
    from chython import ReactionContainer


class GroupsFix:
    __slots__ = ()

    def fix_groups_mapping(self: 'ReactionContainer', *, logging: bool = False) -> \
            Union[bool, List[Tuple[str, Tuple[int, ...]]]]:
        """
        Fix atom-to-atom mapping of some functional groups. Return True if found AAM errors.
        """
        if not self:
            if logging:
                return []
            return False

        log = []
        seen = set()
        remap = {}
        pamer = {}
        r_groups = set()
        p_groups = set()
        r_subs = set()
        p_subs = set()

        # find xonyl groups. any charged-neutral combinations
        for pattern in xonyl_groups:
            for m, g in chain(zip(self.reactants, repeat(r_groups)), zip(self.products, repeat(p_groups))):
                atoms = m._atoms
                for mapping in pattern.get_mapping(m, automorphism_filter=False):
                    n1, n2, n3 = mapping[1], mapping[2], mapping[3]
                    if (t := atoms[n2].atomic_number) == atoms[n3].atomic_number:
                        g.add((n1, n2, n3, atoms[n1].atomic_number, t))

        for pattern, _map in substituents_groups:
            for m, g in chain(zip(self.reactants, repeat(r_subs)), zip(self.products, repeat(p_subs))):
                atoms = m._atoms
                for mapping in pattern.get_mapping(m, automorphism_filter=False):
                    g.add((n := mapping[1], atoms[n].atomic_number,
                           tuple((n := mapping[x], y - 2, atoms[n].atomic_number) for x, y in _map), m))

        r_groups = list(r_groups)
        p_groups = list(p_groups)

        # find pairs
        if r_groups and p_groups:
            for n1, n2, n3, x1, x2 in r_groups:
                if n1 in seen:  # already remapped
                    continue
                for i, (m1, m2, m3, y1, y2) in enumerate(p_groups):
                    if m1 not in seen and n1 == m1 and x1 == y1 and x2 == y2:  # found pair
                        if n2 == m3 and n3 == m2:  # found switch
                            remap[m2] = m3
                            remap[m3] = m2
                            seen.add(n1)
                        break
                else:
                    continue
                del p_groups[i]

        if not p_groups:  # optimize
            r_subs.clear()

        # hydrolysis, etc.
        for (n1, x1, _map, m), g, r in chain(zip(r_subs, repeat(p_groups), repeat(remap)),
                                             zip(p_subs, repeat(r_groups), repeat(pamer))):
            if n1 in seen:
                continue
            for i, (m1, *m23, y1, y2) in enumerate(g):
                if m1 not in seen and n1 == m1 and x1 == y1:  # found center
                    if len(_map) == 1:  # acids substitutions.
                        ni, _, xi = _map[0]
                        # second neighbor should be disconnected from central atom.
                        if xi == y2 and m23[0] == ni and (m23[1] not in m._atoms or m23[1] not in m._bonds[n1]):
                            m2, m3 = m23
                            r[m2] = m3
                            r[m3] = m2
                            seen.add(n1)
                            break
                    elif all(xi == y2 and m23[mi] == ni for ni, mi, xi in _map):
                        m2, m3 = m23
                        r[m2] = m3
                        r[m3] = m2
                        seen.add(n1)
                        break
            else:
                continue
            del g[i]

        if remap:
            seen = set(remap)
            for m in self.products:
                if not seen.isdisjoint(m):
                    m.remap(remap)
            log.append(('products groups remapped', tuple(remap)))
        if pamer:
            seen = set(pamer)
            for m in self.reactants:
                if not seen.isdisjoint(m):
                    m.remap(pamer)
            log.append(('reactants groups remapped', tuple(pamer)))

        if log:
            self.flush_cache()
            if logging:
                return log
            return True
        elif logging:
            return []
        return False


__all__ = ['GroupsFix']
