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
from collections import defaultdict
from typing import TYPE_CHECKING
from ._rules import freak_rules
from .._rings import sssr
from ..rings import _connected_components


if TYPE_CHECKING:
    from chython import MoleculeContainer


# atomic number constants
B = 5
C = 6
N = 7
O = 8
P = 15
S = 16
Se = 34


class Thiele:
    __slots__ = ()

    def thiele(self: 'MoleculeContainer', *, fix_tautomers=True) -> bool:
        """
        Convert structure to aromatic form (Huckel rule ignored). Return True if found any kekule ring.
        Also marks atoms as aromatic.

        :param fix_tautomers: try to fix condensed rings with pyrroles.
            N1C=CC2=NC=CC2=C1>>N1C=CC2=CN=CC=C12
        """
        atoms = self._atoms
        bonds = self._bonds
        nsc = self.not_special_connectivity

        rings = defaultdict(set)  # aromatic? skeleton. include quinones
        tetracycles = []
        pyrroles = set()
        acceptors = set()
        donors = []
        freaks = []
        for ring in self.sssr:
            lr = len(ring)
            if not 3 < lr < 8:  # skip 3-membered and big rings
                continue
            # only B C N O P S with 2-3 neighbors. detects this: C1=CC=CP12=CC=CC=C2
            if any(atoms[n] not in (C, N, O, S, B, P) or len(nsc[n]) > 3 for n in ring):
                continue
            sp2 = sum(atoms[n].hybridization == 2 for n in ring)
            if sp2 == lr:  # benzene like
                if lr == 4:  # two bonds condensed aromatic rings
                    tetracycles.append(ring)
                else:
                    if fix_tautomers and lr % 2:  # find potential pyrroles
                        acceptors.update(n for n in ring if (a := atoms[n]) == N and not a.charge)
                    n, *_, m = ring
                    rings[n].add(m)
                    rings[m].add(n)
                    for n, m in zip(ring, ring[1:]):
                        rings[n].add(m)
                        rings[m].add(n)
            elif 4 < lr == sp2 + 1:  # pyrroles, furanes, etc
                try:
                    n = next(n for n in ring if atoms[n].hybridization == 1)
                except StopIteration:  # exotic, just skip
                    continue
                if (a := atoms[n]).charge == -1:
                    if a != C or lr != 5:  # skip any but ferrocene
                        continue
                elif a.charge:  # skip any charged
                    continue
                elif lr == 7:  # skip electron-rich 7-membered rings
                    if a != 5:  # not B?
                        continue
                # below lr == 5 or 6 only
                elif a in (O, S, Se):
                    if len(bonds[n]) != 2:  # like CS1(C)C=CC=C1
                        continue
                elif a == N:
                    if (b := len(bonds[n])) > 3:  # extra check for invalid N(IV)
                        continue
                    elif fix_tautomers and lr == 6 and b == 2:
                        donors.append(n)
                elif a in (B, P):
                    if len(bonds[n]) > 3:
                        continue
                else:  # only B, [C-], N, O, P, S, Se
                    continue

                pyrroles.add(n)
                n, *_, m = ring
                rings[n].add(m)
                rings[m].add(n)
                for n, m in zip(ring, ring[1:]):
                    rings[n].add(m)
                    rings[m].add(n)
            # like N1C=Cn2cccc12, S1C=Cn2cccc12
            elif lr == 5 and sp2 == 3:
                freaks.append(ring)
        if not rings:
            return False

        # check out-of-ring double bonds
        double_bonded = {n for n in rings if any(m not in rings[n] and b == 2 for m, b in bonds[n].items())}

        # fix_tautomers
        if fix_tautomers and acceptors and donors:
            for start in donors:
                stack = [(start, n, 0, 2) for n in rings[start] if n not in double_bonded]
                path = []
                seen = {start}
                while stack:
                    last, current, depth, order = stack.pop()
                    if len(path) > depth:
                        seen.difference_update(x for _, x, _ in path[depth:])
                        path = path[:depth]
                    path.append((last, current, order))
                    if current in acceptors:  # we found
                        if order == 1:
                            acceptors.discard(current)
                            pyrroles.discard(start)
                            pyrroles.add(current)
                            atoms[current]._implicit_hydrogens = 1
                            atoms[start]._implicit_hydrogens = 0
                            break
                        else:
                            continue

                    depth += 1
                    seen.add(current)
                    new_order = 1 if order == 2 else 2
                    stack.extend((current, n, depth, new_order) for n in rings[current] if
                                 n not in seen and n not in double_bonded and bonds[current][n] == order)
                else:  # path not found
                    continue
                for n, m, o in path:
                    bonds[n][m]._order = o
                if not acceptors:
                    break
            self.flush_cache(keep_sssr=True, keep_components=True)
            self.calc_labels()

        if double_bonded:  # delete quinones
            for n in double_bonded:
                for m in rings.pop(n):
                    rings[m].discard(n)

            for n in [n for n, ms in rings.items() if not ms]:  # imide leads to isolated atoms
                del rings[n]
            if not rings:
                return False
            while True:
                try:
                    n = next(n for n, ms in rings.items() if len(ms) == 1)
                except StopIteration:
                    break
                m = rings.pop(n).pop()
                if n in pyrroles:
                    rings[m].discard(n)
                else:
                    pm = rings.pop(m)
                    pm.discard(n)
                    for x in pm:
                        rings[x].discard(m)
            if not rings:
                return False

        n_sssr = sum(len(x) for x in rings.values()) // 2 - len(rings) + len(_connected_components(rings))
        if not n_sssr:
            return False
        rings = sssr(rings, n_sssr)  # search rings again

        seen = set()
        for ring in rings:
            seen.update(ring)

        # reset bonds to single
        for ring in tetracycles:
            if seen.issuperset(ring):
                n, *_, m = ring
                bonds[n][m]._order = 1
                for n, m in zip(ring, ring[1:]):
                    bonds[n][m]._order = 1

        for ring in rings:
            n, *_, m = ring
            bonds[n][m]._order = 4
            for n, m in zip(ring, ring[1:]):
                bonds[n][m]._order = 4

        self.flush_cache(keep_sssr=True, keep_components=True)
        self.calc_labels()
        for ring in freaks:  # aromatize rule based
            for q in freak_rules:
                if next(q.get_mapping(self, searching_scope=ring, automorphism_filter=False), None):
                    n, *_, m = ring
                    bonds[n][m]._order = 4
                    for n, m in zip(ring, ring[1:]):
                        bonds[n][m]._order = 4
                    break
        if freaks:
            self.flush_cache(keep_sssr=True, keep_components=True)  # flush again
            self.calc_labels()
        self.fix_stereo()  # check if any stereo centers vanished.
        return True


__all__ = ['Thiele']
