# -*- coding: utf-8 -*-
#
#  Copyright 2014-2025 Ramil Nugmanov <nougmanoff@protonmail.com>
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
                 fix_broken_pyrroles: bool = False, fix_tautomers: bool = True, copy_metadata: bool = False):
        """
        :param pattern: Search pattern.
        :param replacement: Resulted structure.
        :param delete_atoms: If True atoms exists in reactants but not exists in products will be removed.
        :param fix_aromatic_rings: Proceed kekule and thiele on products.
        :param fix_tautomers: See `thiele()` docs.
        :param automorphism_filter: Skip matches to same atoms.
        :param copy_metadata: Copy metadata from structure to transformed.
        :param fix_broken_pyrroles: fix invalid rings like Cn1cc[nH]c1.
        """
        if not isinstance(pattern, QueryContainer) or not isinstance(replacement, (MoleculeContainer, QueryContainer)):
            raise TypeError('invalid params')

        self._pattern = pattern
        self._automorphism_filter = automorphism_filter
        self._copy_metadata = copy_metadata
        super().__init__(pattern, replacement, delete_atoms, fix_aromatic_rings, fix_tautomers, fix_broken_pyrroles)

    def __call__(self, structure: MoleculeContainer):
        if not isinstance(structure, MoleculeContainer):
            raise TypeError('only Molecules possible')

        for mapping in self._pattern.get_mapping(structure, automorphism_filter=self._automorphism_filter):
            transformed = self._patcher(structure, mapping)
            if self._copy_metadata:
                transformed.meta.update(structure.meta)
            yield transformed


__all__ = ['Transformer']
