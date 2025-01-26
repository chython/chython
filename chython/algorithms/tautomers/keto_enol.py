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
from collections import defaultdict
from functools import cached_property
from itertools import chain, repeat
from typing import TYPE_CHECKING, Union
from ._keto_enol import *


if TYPE_CHECKING:
    from chython import MoleculeContainer


# atomic number constants
C = 6


class KetoEnol:
    __slots__ = ()

    def _enumerate_keto_enol_tautomers(self: Union['MoleculeContainer', 'KetoEnol'], partial=False):
        for fix, ket in self.__enumerate_bonds(partial):
            if ket:
                a = fix[-1][1]
                d = fix[0][0]
            else:
                a = fix[0][0]
                d = fix[-1][1]

            mol = self.copy(keep_sssr=True, keep_components=True)
            m_bonds = mol._bonds
            for n, m, b in fix:
                m_bonds[n][m]._order = b

            a = mol._atoms[a]
            d = mol._atoms[d]
            a._implicit_hydrogens += 1
            d._implicit_hydrogens -= 1
            a._hybridization -= 1  # -C=X>=C-X or -C=C=X>=C-C=X
            d._hybridization += 1
            yield mol, ket

    @cached_property
    def _sugar_groups(self):
        ek = []
        for mapping in sugar_group.get_mapping(self, automorphism_filter=False):
            e, k = mapping[1], mapping[2]
            ek.append((e, k))
        return ek

    def __enumerate_bonds(self: 'MoleculeContainer', partial):
        atoms = self._atoms
        bonds = self._bonds
        rings = self.atoms_rings_sizes

        # search neutral oxygen and nitrogen
        donors = defaultdict(set)
        acceptors = defaultdict(set)
        for q in enol_rules:
            for mapping in q.get_mapping(self, automorphism_filter=False):
                donors[mapping[1]].add(mapping[2])
        for q in keto_rules:
            for mapping in q.get_mapping(self, automorphism_filter=False):
                acceptors[mapping[1]].add(mapping[2])

        for (atom, dirs), hydrogen, anti in chain(zip(donors.items(), repeat(True), repeat(acceptors)),
                                                  zip(acceptors.items(), repeat(False), repeat(donors))):
            path = []
            seen = {atom}
            stack = [(atom, n, 2 if hydrogen else 1, 0) for n in dirs]
            while stack:
                last, current, bond, depth = stack.pop()

                if partial and path and not len(path) % 2 and \
                        (hydrogen or  # enol > ketone
                         atoms[(x := path[-1][1])].implicit_hydrogens and
                         (x not in rings or all(x > 7 for x in rings[x]))):  # ketone>
                    # return partial hops. ignore allenes in small rings.
                    yield path, hydrogen
                if len(path) > depth:  # fork found
                    if not partial and not len(path) % 2 and (hydrogen or atoms[path[-1][1]].implicit_hydrogens):
                        # end of path found. return it and start new one.
                        yield path, hydrogen
                    seen.difference_update(x for _, x, _ in path[depth:])
                    path = path[:depth]

                path.append((last, current, bond))

                # adding neighbors
                depth += 1
                seen.add(current)
                if bond == 2:
                    next_bond = 1
                else:
                    next_bond = 2

                for n, b in bonds[current].items():
                    if n == last:
                        continue
                    elif n in seen:  # aromatic ring destruction. pyridine double bonds shift
                        continue
                    elif n in anti:  # enol-ketone switch
                        if current in anti[n]:  # keton or enol bond
                            if hydrogen:
                                cp = path.copy()
                                cp.append((current, n, 1))  # double to single in keton end
                                yield cp, True
                            else:
                                cp = path.copy()
                                cp.append((current, n, 2))  # single to double in enol end
                                yield cp, False
                    elif b == bond and (a := atoms[n]) == C:  # classic keto-enol route
                        if a.hybridization == 2:  # grow up
                            stack.append((current, n, next_bond, depth))
                        elif hydrogen:
                            if a.hybridization == 3:  # OC=CC=C=C case
                                cp = path.copy()
                                cp.append((current, n, 1))
                                yield cp, True  # ketone found
                        elif a.hybridization == 1 and a.implicit_hydrogens:  # ketone >> enol
                            cp = path.copy()
                            cp.append((current, n, 2))
                            yield cp, False

            if path and not len(path) % 2 and \
                    (hydrogen or  # enol > ketone
                     atoms[(x := path[-1][1])].implicit_hydrogens and
                     (x not in rings or all(x > 7 for x in rings[x]))):
                yield path, hydrogen


__all__ = ['KetoEnol']
