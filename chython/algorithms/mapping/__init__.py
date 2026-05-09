# -*- coding: utf-8 -*-
#
#  Copyright 2022-2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from collections import ChainMap
from itertools import count , chain, repeat
from .attention import Attention
from .reconstruct import Reconstruct
from ._groups import xonyl_groups, substituents_groups
from ._reactions import rules


class Mapping(Reconstruct, Attention):
    __slots__ = ()

    def reset_mapping(self) -> bool:
        """
        Reset atom-to-atom mapping by remapping atoms to unique numbers.
        """
        r = [n for m in chain(self.reactants, self.reagents) for n in m._atoms]
        p = [n for m in self.products for n in m._atoms]
        c = count(1)
        if len(r) != len(set(r)):
            for m in chain(self.reactants, self.reagents):
                m.remap({n: next(c) for n in m._atoms})
        if len(p) != len(set(p)):
            for m in self.products:
                m.remap({n: next(c) for n in m._atoms})
        if next(c) != 1:
            self.flush_cache()
            return True
        return False

    def fix_mapping(self, *, logging: bool = False) -> bool | list:
        """
        Fix mapping of plipped functional groups and common mechanism errors.
        """
        if not self:
            if logging:
                return []
            return False
        log = self.__fix_groups()
        log.extend(self.__fix_cgr())
        if logging:
            return log
        return bool(log)

    def __fix_cgr(self):
        cgr = ~self
        if not cgr.center_atoms:
            return []

        log = []
        free_number = count(max(cgr) + 1)
        components = [(cgr.substructure(c),
                       cgr.augmented_substructure(c, 2),  # deep DEPENDS on rules!
                       c)
                      for c in cgr.substructure(cgr.center_atoms).connected_components]

        r_atoms = ChainMap(*(x._atoms for x in self.reactants))
        for c, ac, cs in components:
            for rule_num, (query, signature, restrict, fix) in enumerate(rules):
                if str(c) == signature:
                    for mapping in query.get_mapping(ac, automorphism_filter=False):
                        if not cs.issubset(mapping.values()):
                            continue
                        if restrict is not None and any(a != r_atoms.get(mapping[n]) for n, a in restrict.atoms()):
                            continue
                        mapping = {mapping[n]: next(free_number) if m is None else mapping[m] for n, m in fix.items()}
                        for m in self.products:
                            m.remap(mapping)
                        log.append((rule_num, signature, tuple(mapping.values())))
                        break
                    else:
                        continue
                    break  # component remapped!

        if log:
            self.flush_cache()
        return log

    def __fix_groups(self):
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
        return log


__all__ = ['Mapping']
