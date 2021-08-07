# -*- coding: utf-8 -*-
#
#  Copyright 2014-2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
                 delete_atoms: bool = True):
        """
        :param pattern: search pattern
        :param replacement: atoms and bonds replacement
        :param delete_atoms: if True atoms exists in pattern but not exists in replacement will be removed
        """
        if not isinstance(pattern, QueryContainer) or not isinstance(replacement, (MoleculeContainer, QueryContainer)):
            raise TypeError('invalid params')

        self.__pattern = pattern
        super().__init__(pattern, replacement, delete_atoms)

    def __call__(self, structure: MoleculeContainer, automorphism_filter: bool = True):
        if not isinstance(structure, MoleculeContainer):
            raise TypeError('only Molecules possible')

        for mapping in self.__pattern.get_mapping(structure, automorphism_filter=automorphism_filter):
            yield self._patcher(structure, mapping)

    def __getstate__(self):
        return {'pattern': self.__pattern, **super().__getstate__()}

    def __setstate__(self, state):
        self.__pattern = state['pattern']
        super().__setstate__(state)


__all__ = ['Transformer']
