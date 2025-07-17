# -*- coding: utf-8 -*-
#
#  Copyright 2021-2025 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from typing import List, TYPE_CHECKING, Union
from ...exceptions import ValenceError


if TYPE_CHECKING:
    from chython import MoleculeContainer


# atomic number constants
B = 5
C = 6
N = 7
O = 8
Si = 14
P = 15
S = 16
As = 33
Se = 34
Te = 52


class Resonance:
    __slots__ = ()

    def fix_resonance(self: Union['MoleculeContainer', 'Resonance'], *, logging=False,
                      _fix_stereo=True) -> Union[bool, List[int]]:
        """
        Transform biradical or dipole resonance structures into neutral form. Return True if structure form changed.

        :param logging: return list of changed atoms.
        """
        atoms = self._atoms
        bonds = self._bonds
        entries, exits, rads, constrains, nitrogen_cat, nitrogen_ani, sulfur_cat = self.__entries()
        hs = set()
        while len(rads) > 1:
            n = rads.pop()
            for path in self.__find_delocalize_path(n, rads, constrains, True):
                atoms[n]._is_radical = False
                hs.add(n)
                for n, m, b in path:
                    hs.add(m)
                    bonds[n][m]._order = b
                atoms[m]._is_radical = False  # noqa
                rads.discard(m)
                break  # path found
            # path not found. atom n keep as is
        while entries and exits:
            n = entries.pop()
            for path in self.__find_delocalize_path(n, exits, constrains, False):
                l, m, b = path[-1]
                if n in nitrogen_cat and m in nitrogen_ani:
                    continue

                if m in sulfur_cat:  # prevent X-[S+]=X >> X=S=X
                    if b != 1:
                        continue
                    atoms[m]._charge -= 1
                else:  # check cations end valence.
                    atoms[m]._charge -= 1  # reduce atom change and check valence
                    try:
                        atoms[m].valence_rules(sum(int(y) for x, y in bonds[m].items() if x != l) + b)
                    except ValenceError:
                        atoms[m]._charge += 1  # roll back
                        continue

                # succeed!
                atoms[n]._charge += 1
                hs.add(n)
                for n, m, b in path:
                    hs.add(m)
                    bonds[n][m]._order = b
                exits.discard(m)
                break  # path from negative atom to positive atom found.
            # path not found. keep negative atom n as is
        if hs:
            for n in hs:
                self.calc_implicit(n)
            self.flush_cache(keep_sssr=True, keep_components=True)
            if _fix_stereo:
                self.fix_stereo()
            if logging:
                return list(hs)
            return True
        elif logging:
            return []
        return False

    def __find_delocalize_path(self: 'MoleculeContainer', start, finish, constrains, odd_only):
        bonds = self._bonds
        stack = [(start, n, 0, b.order + 1) for n, b in bonds[start].items() if n in constrains and b.order < 3]
        path = []
        seen = {start}
        while stack:
            last, current, depth, order = stack.pop()
            if len(path) > depth:
                seen.difference_update(x for _, x, _ in path[depth:])
                path = path[:depth]

            path.append((last, current, order))

            if current in finish:
                if odd_only:  # radicals!
                    if len(path) % 2:
                        yield path
                    else:  # invalid path
                        continue
                elif depth:  # one bonded ignored. we search double bond transfer! A=A-A >> A-A=A.
                    yield path

            depth += 1
            seen.add(current)
            diff = -1 if depth % 2 else 1
            stack.extend((current, n, depth, bo) for n, b in bonds[current].items()
                         if n not in seen and n in constrains and 1 <= (bo := b.order + diff) <= 3)

    def __entries(self: 'MoleculeContainer'):
        atoms = self._atoms
        bonds = self._bonds
        errors = {n for n, a in self.atoms() if a.implicit_hydrogens is None}

        transfer = set()
        entries = set()
        exits = set()
        rads = set()
        nitrogen_cat = set()
        nitrogen_ani = set()
        sulfur_cat = set()
        for n, a in self.atoms():
            if a not in (B, C, N, O, Si, P, S, As, Se, Te):
                # filter non-organic set, halogens and aromatics
                continue
            elif a.is_radical:
                rads.add(n)
            elif a.charge == -1:
                if a == B and a.total_hydrogens:  # skip [BHx-]
                    continue
                elif (lb := len(bonds[n])) == 4 and a == B:  # skip [B-]X4
                    continue
                elif lb == 6 and a == P:  # skip [P-]X6
                    continue
                if n in errors:  # only valid anions accepted
                    continue
                entries.add(n)
            elif a.charge == 1:
                lb = len(bonds[n])
                if a == N:
                    if lb == 4:  # skip ammonia
                        continue
                    elif lb == 2 and a.hybridization == 3:  # skip Azide
                        (n1, b1), (n2, b2) = bonds[n].items()
                        an1 = atoms[n1]
                        an2 = atoms[n2]
                        if b1 == b2 == 2 and (an1.charge == -1 and an1 == N or an2.charge == -1 and an2 == N):
                            continue
                    elif lb == 3 and a.hybridization == 2:  # X=[N+](-X)-X - prevent N-N migration
                        nitrogen_ani.add(n)
                elif a == P and lb == 4:  # skip [P+]R4
                    continue
                elif a == S:
                    if lb == 2 and a.hybridization == 2:  # ad-hoc for X-[S+]=X
                        sulfur_cat.add(n)
                    elif lb == 3 and a.hybridization == 1:  # ad-hoc for X-[S+](-X)-X
                        continue
                exits.add(n)
            transfer.add(n)

        if exits or entries:  # try to move cation to nitrogen. saturation fixup.
            for n, a in self.atoms():
                if a == N and not a.charge:
                    if a.hybridization == 1 and a.neighbors <= 3:  # any amine - potential e-donor
                        entries.add(n)
                        nitrogen_cat.add(n)
                    elif a.hybridization == 3 and a.neighbors == 1:  # N#X-[X-] >> [N-]=X=X
                        exits.add(n)
                        nitrogen_ani.add(n)
        return entries, exits, rads, transfer, nitrogen_cat, nitrogen_ani, sulfur_cat


__all__ = ['Resonance']
