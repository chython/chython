# -*- coding: utf-8 -*-
#
#  Copyright 2025 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from typing import Literal, TYPE_CHECKING


if TYPE_CHECKING:
    from chython import MoleculeContainer


class Conformers:
    __slots__ = ()

    def generate_conformers(self: 'MoleculeContainer', limit: int = 10, *, optimize: bool = False,
                            engine: Literal['rdkit'] = None, **kwargs) -> int:
        """
        Generate conformers for the molecule ignoring implicit hydrogens. Set them manually to have a full 3D structure.
        By default, the RDKit engine is used. Can be changed globally with the `chython.conformer_engine` parameter.

        :param limit: maximum number of conformers to generate
        :param optimize: optimize conformers using MMFF94 force field (only for RDKit engine)
        :param engine: override globally set engine
        :param kwargs: additional arguments for the engine

        :return: number of generated conformers
        """
        if engine is None:
            from chython import conformer_engine as engine
        if engine == 'rdkit':
            from rdkit.Chem.AllChem import EmbedMultipleConfs, MMFFOptimizeMolecule

            copy = self.copy()
            copy.explicify_hydrogens()
            rmol = copy.to_rdkit(keep_mapping=False, keep_hydrogens=False)
            ids = EmbedMultipleConfs(rmol, numConfs=limit, **kwargs)
            if optimize:
                for i in ids:
                    MMFFOptimizeMolecule(rmol, confId=i)

            conformers = [
                {n: tuple(v) for n, v in zip(self, conf.GetPositions())}
                for conf in rmol.GetConformers() if conf.Is3D()
            ]
        else: raise ValueError(f'Invalid conformer generation engine: {engine}')
        if conformers:
            self._conformers = conformers
        return len(conformers)


__all__ = ['Conformers']
