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
from ..tautomers._base import stripped_rules
from ...periodictable.element.query import _inorganic


if TYPE_CHECKING:
    from chython import MoleculeContainer


class Salts:
    __slots__ = ()

    def neutralize_metal_salts(self: 'MoleculeContainer', *, logging=False) -> Union[bool, Tuple[List[int], List[int]]]:
        """
        Convert metal salts to mixture of metal base and neutral anion form. Works only for stripped salts.

        Example: [K+].CC(=O)[O-] >> [K+].[OH-].CC(=O)O
        Note: do '.neutralize()' procedure before for preventing ambiguous results.

        :param logging: return changed atoms list and hydroxides list
        """
        charges = self._charges
        bonds = self._bonds
        hydrogens = self._hydrogens

        metals = []
        for n, a in self._atoms.items():
            if a.atomic_symbol not in _inorganic and charges[n] > 0 and not bonds[n]:
                metals.append(n)

        if metals:
            acceptors = set()
            for q in stripped_rules[:-1]:  # except halogenides and hydroxy group.
                for mapping in q.get_mapping(self, automorphism_filter=False):
                    acceptors.add(mapping[1])

            # for imbalanced structures neutralize only part.
            acceptors = list(acceptors)[:sum(charges[n] for n in metals)]
            hydroxides = []
            for n in acceptors:
                hydrogens[n] += 1
                charges[n] += 1
                hydroxides.append(self.add_atom('O', charge=-1))

            if logging:
                return acceptors, hydroxides
            return bool(acceptors)
        elif logging:
            return [], []
        return False

    def remove_metals(self: 'MoleculeContainer', *, logging=False) -> Union[bool, List]:
        """
        Remove disconnected metals.

        :param logging: return deleted atoms list.
        """
        bonds = self._bonds

        metals = []
        for n, a in self._atoms.items():
            if a.atomic_symbol not in _inorganic and not bonds[n]:
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


__all__ = ['Salts']
