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
from collections import deque, defaultdict
from itertools import product
from typing import TYPE_CHECKING
from ..aromatics.kekule import _kekule_component
from ...exceptions import InvalidAromaticRing


if TYPE_CHECKING:
    from chython import MoleculeContainer


class HeteroArenes:
    __slots__ = ()

    def _enumerate_hetero_arene_tautomers(self: 'MoleculeContainer'):
        atoms = self._atoms
        bonds = self._bonds
        hydrogens = self._hydrogens
        charges = self._charges
        radicals = self._radicals

        rings = defaultdict(list)  # aromatic skeleton
        for n, m_bond in bonds.items():
            for m, bond in m_bond.items():
                if bond.order == 4:
                    rings[n].append(m)
        if not rings:
            return

        acceptors = set()
        donors = set()
        single_bonded = set()
        for n, ms in rings.items():
            if len(ms) == 2:
                if atoms[n].atomic_number in (5, 7, 15):
                    if not charges[n] and not radicals[n]:
                        # only neutral B, N, P
                        if hydrogens[n]:  # pyrrole
                            donors.add(n)
                        elif len(bonds[n]) == 2:  # pyridine
                            acceptors.add(n)
                        else:
                            single_bonded.add(n)
                elif charges[n] == -1 and atoms[n].atomic_number == 6:  # ferrocene
                    single_bonded.add(n)
            elif len(ms) == 3 and atoms[n].atomic_number in (5, 7, 15) and not charges[n] and not radicals[n]:
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
                mol = self.copy()
                mol._hydrogens[d] = 0
                mol._hydrogens[a] = 1
                yield mol


__all__ = ['HeteroArenes']
