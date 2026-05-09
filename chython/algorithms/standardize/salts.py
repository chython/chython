# -*- coding: utf-8 -*-
#
#  Copyright 2021-2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from ._salts import acids, rules


# atomic number constants
H = 1
N = 7
s_metals = {3, 4, 11, 12, 19, 20, 37, 38, 55, 56, 87, 88}
s_metals_ammonia = {3, 4, 7, 11, 12, 19, 20, 37, 38, 55, 56, 87, 88}


class Salts:
    __slots__ = ()

    def remove_metals(self, *, skip_elements: list[int] = None, logging=False) -> bool | list:
        """
        Remove disconnected S-metals and ammonia.

        :param logging: return deleted atoms list.
        :param skip_elements: skip elements from removing.
        """
        atoms = self._atoms
        bonds = self._bonds
        filters = s_metals_ammonia.difference(skip_elements) if skip_elements else s_metals_ammonia

        metals = []
        for n, a in atoms.items():
            if not bonds[n] and a.atomic_number in filters:
                metals.append(n)

        if 0 < len(metals) < len(self):
            for n in metals:
                del atoms[n]
                del bonds[n]

            self.flush_cache(keep_sssr=True)
            if logging:
                return metals
            return True
        elif logging:
            return []
        return False

    def remove_acids(self, *, logging=False) -> bool | list:
        """
        Remove common acids from organic bases salts.
        Works only for neutral pairs like HA+B. Use `neutralize` before.

        :param logging: return deleted atoms list.
        """
        if self.connected_components_count > 1:
            log = []
            for c in self.connected_components:
                if self.substructure(c, recalculate_hydrogens=False) in acids:
                    log.extend(c)
            if 0 < len(log) < len(self):  # prevent singularity
                atoms = self._atoms
                bonds = self._bonds

                for n in log:
                    del atoms[n]
                    del bonds[n]

                self.flush_cache()
                if logging:
                    return log
                return True
        if logging:
            return []
        return False

    def split_metal_salts(self, *, skip_elements: list[int] = None, logging=False) -> bool | list:
        """
        Split connected S-metal salts to cation/anion pairs.

        :param logging: return deleted bonds list.
        :param skip_elements: skip elements from splitting.
        """
        atoms = self._atoms
        bonds = self._bonds
        filters = s_metals.difference(skip_elements) if skip_elements else s_metals

        metals = [n for n, a in atoms.items() if a.atomic_number in filters]
        if metals:
            acceptors = set()
            log = []
            for q in rules:
                for mapping in q.get_mapping(self, automorphism_filter=False):
                    acceptors.add(mapping[1])

            for n in metals:
                for m in acceptors & bonds[n].keys():
                    if atoms[n].charge == 4:  # prevent overcharging
                        break
                    del bonds[n][m]
                    del bonds[m][n]
                    atoms[n]._charge += 1
                    atoms[m]._charge -= 1
                    log.append((n, m))
            if log:
                self.flush_cache()
                self.fix_stereo()
                if logging:
                    return log
                return True
        if logging:
            return []
        return False


__all__ = ['Salts']
