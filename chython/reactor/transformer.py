# -*- coding: utf-8 -*-
#
#  Copyright 2014-2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from logging import getLogger
from typing import Union
from .base import BaseReactor
from ..containers import QueryContainer, MoleculeContainer


logger = getLogger('chython.reactor')


class Transformer(BaseReactor):
    """
    Editor for molecules.
    generates modified molecules from input molecule using template.
    Transformer calling returns generator of all possible replacements.
    """
    def __init__(self, pattern: QueryContainer, replacement: Union[MoleculeContainer, QueryContainer],
                 delete_atoms: bool = True, automorphism_filter: bool = True, canonicalize: bool = True,
                 fix_tautomers: bool = True, ignore_pyrrole_hydrogen: bool = False,
                 copy_metadata: bool = False):
        """
        :param pattern: Search pattern.
        :param replacement: Resulted structure.
        :param delete_atoms: If True atoms exists in reactants but not exists in products will be removed.
        :param automorphism_filter: Skip matches to same atoms.
        :param canonicalize: Run full canonicalization on products (kekule, standardize, thiele, etc.).
        :param fix_tautomers: Canonicalize tautomer forms. Passed to canonicalize().
        :param ignore_pyrrole_hydrogen: Fix invalid rings like Cn1cc[nH]c1. Passed to canonicalize().
        :param copy_metadata: Copy metadata from structure to transformed.
        """
        if not isinstance(pattern, QueryContainer) or not isinstance(replacement, (MoleculeContainer, QueryContainer)):
            raise TypeError('invalid params')

        self._pattern = pattern
        self._automorphism_filter = automorphism_filter
        self._copy_metadata = copy_metadata
        super().__init__(pattern, replacement, delete_atoms, canonicalize, fix_tautomers, ignore_pyrrole_hydrogen)

    def __call__(self, structure: MoleculeContainer):
        if not isinstance(structure, MoleculeContainer):
            raise TypeError('only Molecules possible')

        for mapping in self._pattern.get_mapping(structure, automorphism_filter=self._automorphism_filter):
            try:
                transformed = self._patcher(structure, mapping)
            except Exception:
                logger.info('invalid product structure, skipping')
                continue
            if self._copy_metadata:
                transformed.meta.update(structure.meta)
            yield transformed


__all__ = ['Transformer']
