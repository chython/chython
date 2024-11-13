# -*- coding: utf-8 -*-
#
#  Copyright 2022-2024 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from collections import deque, defaultdict
from itertools import product
from typing import TYPE_CHECKING
from ..aromatics.kekule import _kekule_component
from ...exceptions import InvalidAromaticRing


if TYPE_CHECKING:
    from chython import MoleculeContainer


# atomic number constants
B = 5
C = 6
N = 7
P = 15


class HeteroArenes:
    __slots__ = ()

    def _enumerate_hetero_arene_tautomers(self: 'MoleculeContainer'):
        atoms = self._atoms
        bonds = self._bonds

        rings = defaultdict(list)  # aromatic skeleton
        for n, m_bond in bonds.items():
            for m, bond in m_bond.items():
                if bond == 4:
                    rings[n].append(m)
        if not rings:
            return

        acceptors = set()
        donors = set()
        single_bonded = set()
        for n, ms in rings.items():
            a = atoms[n]
            if len(ms) == 2:
                if a in (B, N, P):
                    if not a.charge and not a.is_radical:
                        # only neutral B, N, P
                        if a.implicit_hydrogens:  # pyrrole
                            donors.add(n)
                        elif len(bonds[n]) == 2:  # pyridine
                            acceptors.add(n)
                        else:
                            single_bonded.add(n)
                elif a.charge == -1 and a == C:  # ferrocene
                    single_bonded.add(n)
            elif len(ms) == 3 and a in (B, N, P) and not a.charge and not a.is_radical:
                single_bonded.add(n)
        if not donors or not acceptors:
            return

        atoms = set(rings)
        components = []
        while atoms:
            start = atoms.pop()
            component = {start: rings[start]}
            queue = deque([start])
            while queue:
                current = queue.popleft()
                for n in rings[current]:
                    if n not in component:
                        queue.append(n)
                        component[n] = rings[n]

            atoms.difference_update(component)
            if donors.isdisjoint(component) or acceptors.isdisjoint(component):
                continue
            components.append(component)

        if not components:
            return
        for component in components:
            for d, a in product(component.keys() & donors, component.keys() & acceptors):
                sb = component.keys() & single_bonded
                sb.add(a)  # now pyrrole
                try:
                    next(_kekule_component(component, sb, (), 0))
                except InvalidAromaticRing:
                    continue
                mol = self.copy(keep_sssr=True, keep_components=True)
                mol._atoms[d]._implicit_hydrogens = 0
                mol._atoms[a]._implicit_hydrogens = 1
                yield mol


__all__ = ['HeteroArenes']
