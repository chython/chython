# -*- coding: utf-8 -*-
#
#  Copyright 2014-2022 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2019 Adelia Fatykhova <adelik21979@gmail.com>
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
from typing import Union
from .base import BaseReactor
from ..containers import QueryContainer, MoleculeContainer


class Transformer(BaseReactor):
    """
    Editor for molecules.
    generates modified molecules from input molecule using template.
    Transformer calling returns generator of all possible replacements.
    """
    def __init__(self, pattern: QueryContainer, replacement: Union[MoleculeContainer, QueryContainer],
                 delete_atoms: bool = True, automorphism_filter: bool = True, fix_aromatic_rings: bool = True,
                 fix_tautomers: bool = True):
        """
        :param pattern: Search pattern.
        :param replacement: Resulted structure.
        :param delete_atoms: If True atoms exists in reactants but not exists in products will be removed.
        :param fix_aromatic_rings: Proceed kekule and thiele on products.
        :param fix_tautomers: See `thiele()` docs.
        :param automorphism_filter: Skip matches to same atoms.
        """
        if not isinstance(pattern, QueryContainer) or not isinstance(replacement, (MoleculeContainer, QueryContainer)):
            raise TypeError('invalid params')

        self.pattern = pattern
        self.replacement = replacement
        self.__automorphism_filter = automorphism_filter
        super().__init__({n for n, h in pattern._masked.items() if not h}, replacement, delete_atoms,
                         fix_aromatic_rings, fix_tautomers)

    def __call__(self, structure: MoleculeContainer):
        if not isinstance(structure, MoleculeContainer):
            raise TypeError('only Molecules possible')

        for mapping in self.pattern.get_mapping(structure, automorphism_filter=self.__automorphism_filter):
            yield from self._patcher(structure, mapping)


__all__ = ['Transformer']
