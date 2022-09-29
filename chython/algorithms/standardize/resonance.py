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
from typing import List, TYPE_CHECKING, Union
from ...exceptions import ValenceError


if TYPE_CHECKING:
    from chython import MoleculeContainer


class Resonance:
    __slots__ = ()

    def fix_resonance(self: Union['MoleculeContainer', 'Resonance'], *, logging=False,
                      _fix_stereo=True) -> Union[bool, List[int]]:
        """
        Transform biradical or dipole resonance structures into neutral form. Return True if structure form changed.

        :param logging: return list of changed atoms.
        """
        atoms = self._atoms
        charges = self._charges
        radicals = self._radicals
        bonds = self._bonds
        calc_implicit = self._calc_implicit
        entries, exits, rads, constrains, nitrogen_cat, nitrogen_ani, sulfur_cat = self.__entries()
        hs = set()
        while len(rads) > 1:
            n = rads.pop()
            for path in self.__find_delocalize_path(n, rads, constrains, True):
                radicals[n] = False
                hs.add(n)
                for n, m, b in path:
                    hs.add(m)
                    bonds[n][m]._Bond__order = b  # noqa
                radicals[m] = False  # noqa
                rads.discard(m)
                break  # path found
            # path not found. atom n keep as is
        while entries and exits:
            n = entries.pop()
            for path in self.__find_delocalize_path(n, exits, constrains, False):
                l, m, b = path[-1]
                if n in nitrogen_cat and m in nitrogen_ani:
                    continue

                c_m = charges[m] - 1
                if m in sulfur_cat:  # prevent X-[S+]=X >> X=S=X
                    if b != 1:
                        continue
                else:  # check cations end valence.
                    try:
                        atoms[m].valence_rules(c_m, radicals[m], sum(int(y) for x, y in bonds[m].items() if x != l) + b)
                    except ValenceError:
                        continue

                charges[n] += 1
                hs.add(n)
                for n, m, b in path:
                    hs.add(m)
                    bonds[n][m]._Bond__order = b  # noqa
                charges[m] = c_m
                exits.discard(m)
                break  # path from negative atom to positive atom found.
            # path not found. keep negative atom n as is
        if hs:
            for n in hs:
                calc_implicit(n)
            self.flush_cache()
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
        hybridization = self.hybridization
        neighbors = self.neighbors
        charges = self._charges
        radicals = self._radicals
        bonds = self._bonds
        atoms = self._atoms
        errors = {n for n, h in self._hydrogens.items() if h is None}

        transfer = set()
        entries = set()
        exits = set()
        rads = set()
        nitrogen_cat = set()
        nitrogen_ani = set()
        sulfur_cat = set()
        for n, a in atoms.items():
            if a.atomic_number not in {5, 6, 7, 8, 14, 15, 16, 33, 34, 52}:
                # filter non-organic set, halogens and aromatics
                continue
            elif radicals[n]:
                rads.add(n)
            elif charges[n] == -1:
                if (lb := len(bonds[n])) == 4 and a.atomic_number == 5:  # skip boron
                    continue
                elif lb == 6 and a.atomic_number == 15:  # skip [P-]X6
                    continue
                if n in errors:  # only valid anions accepted
                    continue
                entries.add(n)
            elif charges[n] == 1:
                lb = len(bonds[n])
                if a.atomic_number == 7:
                    if lb == 4:  # skip ammonia
                        continue
                    elif lb == 2 and hybridization(n) == 3:  # skip Azide
                        (n1, b1), (n2, b2) = bonds[n].items()
                        if b1.order == b2.order == 2 and (charges[n1] == -1 and atoms[n1].atomic_number == 7 or
                                                          charges[n2] == -1 and atoms[n2].atomic_number == 7):
                            continue
                    elif lb == 3 and hybridization(n) == 2:  # X=[N+](-X)-X - prevent N-N migration
                        nitrogen_ani.add(n)
                elif a.atomic_number == 15 and lb == 4:  # skip [P+]R4
                    continue
                elif a.atomic_number == 16:
                    if lb == 2 and hybridization(n) == 2:  # ad-hoc for X-[S+]=X
                        sulfur_cat.add(n)
                    elif lb == 3 and hybridization(n) == 1:  # ad-hoc for X-[S+](-X)-X
                        continue
                exits.add(n)
            transfer.add(n)

        if exits or entries:  # try to move cation to nitrogen. saturation fixup.
            for n, a in self._atoms.items():
                if a.atomic_number == 7 and not charges[n]:
                    if hybridization(n) == 1 and neighbors(n) <= 3:  # any amine - potential e-donor
                        entries.add(n)
                        nitrogen_cat.add(n)
                    elif hybridization(n) == 3 and neighbors(n) == 1:  # N#X-[X-] >> [N-]=X=X
                        exits.add(n)
                        nitrogen_ani.add(n)
        return entries, exits, rads, transfer, nitrogen_cat, nitrogen_ani, sulfur_cat


__all__ = ['Resonance']
