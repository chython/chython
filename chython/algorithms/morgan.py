# -*- coding: utf-8 -*-
#
#  Copyright 2017-2024 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from functools import cached_property
from itertools import groupby
from logging import getLogger
from operator import itemgetter
from typing import Dict, TYPE_CHECKING


logger = getLogger('chython.morgan')


if TYPE_CHECKING:
    from chython.containers import MoleculeContainer


class Morgan:
    __slots__ = ()

    @cached_property
    def atoms_order(self: 'MoleculeContainer') -> Dict[int, int]:
        """
        Morgan like algorithm for graph nodes ordering

        :return: dict of atom-order pairs
        """
        if not self:  # for empty containers
            return {}
        elif len(self) == 1:  # optimize single atom containers
            return dict.fromkeys(self, 1)
        return _morgan({n: hash(a) for n, a in self.atoms()}, self.int_adjacency)

    @cached_property
    def int_adjacency(self: 'MoleculeContainer') -> Dict[int, Dict[int, int]]:
        """
        Adjacency with integer-coded bonds.
        """
        return {n: {m: hash(b) for m, b in mb.items()} for n, mb in self._bonds.items()}


def _morgan(atoms: Dict[int, int], bonds: Dict[int, Dict[int, int]]) -> Dict[int, int]:
    tries = len(atoms) - 1
    numb = len(set(atoms.values()))
    stab = old_numb = 0

    for _ in range(tries):
        atoms = {n: hash((atoms[n], *(x for x in sorted((atoms[m], b) for m, b in ms.items()) for x in x)))
                 for n, ms in bonds.items()}
        old_numb, numb = numb, len(set(atoms.values()))
        if numb == len(atoms):  # each atom now unique
            break
        elif numb == old_numb:  # not changed. molecules like benzene
            if stab == 3:
                break
            stab += 1
        elif stab:  # changed unique atoms number. reset stability check.
            stab = 0
    else:
        if numb < old_numb:
            logger.warning('number of attempts exceeded. uniqueness has decreased.')

    return {n: i for i, (_, g) in enumerate(groupby(sorted(atoms.items(), key=itemgetter(1)), key=itemgetter(1)),
                                            start=1) for n, _ in g}


__all__ = ['Morgan']
