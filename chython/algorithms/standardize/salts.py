# -*- coding: utf-8 -*-
#
#  Copyright 2021, 2022 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from ._salts import rules


if TYPE_CHECKING:
    from chython import MoleculeContainer


class Salts:
    __slots__ = ()

    def remove_metals(self: 'MoleculeContainer', *, logging=False) -> Union[bool, List]:
        """
        Remove disconnected S-metals.

        :param logging: return deleted atoms list.
        """
        bonds = self._bonds

        metals = []
        for n, a in self._atoms.items():
            if a.atomic_symbol not in {3, 4, 11, 12, 19, 20, 37, 38, 55, 56} and not bonds[n]:
                metals.append(n)

        if metals:
            for n in metals:
                self.delete_atom(n)
            if logging:
                return metals
            return True
        elif logging:
            return []
        return False

    def split_metal_salts(self: 'MoleculeContainer', *, logging=False) -> Union[bool, List[Tuple[int, int]]]:
        """
        Split connected S-metal/lanthanides/actinides salts to cation/anion pairs.

        :param logging: return deleted bonds list.
        """
        bonds = self._bonds
        charges = self._charges

        metals = [n for n, a in self._atoms.items() if a.atomic_number in
                  {3, 4, 11, 12, 19, 20, 37, 38, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 87, 88,
                   89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102}]
        if metals:
            acceptors = set()
            log = []
            for q in rules:
                for mapping in q.get_mapping(self, automorphism_filter=False):
                    acceptors.add(mapping[1])

            for n in metals:
                for m in acceptors & bonds[n].keys():
                    del bonds[n][m]
                    del bonds[m][n]
                    charges[n] += 1
                    charges[m] -= 1
                    log.append((n, m))
            if log:
                self.flush_cache()
            if logging:
                return log
            return bool(log)
        if logging:
            return []
        return False


__all__ = ['Salts']
