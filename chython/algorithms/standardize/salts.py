# -*- coding: utf-8 -*-
#
#  Copyright 2021-2024 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from typing import TYPE_CHECKING, List, Tuple, Union
from ._salts import acids, rules
from ...periodictable import GroupI, GroupII


if TYPE_CHECKING:
    from chython import MoleculeContainer


# atomic number constants
H = 1
N = 7


class Salts:
    __slots__ = ()

    def remove_metals(self: 'MoleculeContainer', *, logging=False) -> Union[bool, List]:
        """
        Remove disconnected S-metals and ammonia.

        :param logging: return deleted atoms list.
        """
        atoms = self._atoms
        bonds = self._bonds

        metals = []
        for n, a in atoms.items():
            if not bonds[n] and (a == N or isinstance(a, (GroupI, GroupII)) and a != H):
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

    def remove_acids(self: 'MoleculeContainer', *, logging=False) -> Union[bool, List[int]]:
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

    def split_metal_salts(self: 'MoleculeContainer', *, logging=False) -> Union[bool, List[Tuple[int, int]]]:
        """
        Split connected S-metal salts to cation/anion pairs.

        :param logging: return deleted bonds list.
        """
        atoms = self._atoms
        bonds = self._bonds

        metals = [n for n, a in atoms.items() if isinstance(a, (GroupI, GroupII)) and a != H]
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
