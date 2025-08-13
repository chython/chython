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
from io import StringIO
from typing import Literal, TYPE_CHECKING


if TYPE_CHECKING:
    from chython import MoleculeContainer


class Conformers:
    __slots__ = ()

    def generate_conformers(self: 'MoleculeContainer', limit: int = 10, *, optimize: bool = False,
                            engine: Literal['rdkit', 'cdpkit'] = None, **kwargs) -> int:
        """
        Generate conformers for the molecule ignoring implicit hydrogens. Set them manually to have a full 3D structure.
        By default, the RDKit engine is used. Can be changed globally with the `chython.conformer_engine` parameter.

        Two conformer generation engines are supported:
        - 'rdkit': Uses RDKit's ETKDG algorithm for conformer generation
        - 'cdpkit': Uses CDPKit's ConformerGenerator for conformer generation (https://pubs.acs.org/doi/10.1021/acs.jcim.3c00563)

        :param limit: maximum number of conformers to generate
        :param optimize: optimize conformers using MMFF94 force field (only for RDKit engine)
        :param engine: override globally set engine ('rdkit' or 'cdpkit')
        :param kwargs: additional arguments for the engine:
            - timeout: timeout for the engine in seconds (only for CDPKit engine, default: 60)
            - min_rmsd: minimum RMSD between generated conformers (only for CDPKit engine, default: .5)
            - energy_window: energy window for the engine (only for CDPKit engine, default: 20)

            Check EmbedMultipleConfs API for RDKit engine
        :return: number of generated conformers
        """
        if engine is None:
            from chython import conformer_engine as engine

        copy = self.copy()
        copy.explicify_hydrogens()

        if engine == 'rdkit':
            from rdkit.Chem.AllChem import EmbedMultipleConfs, MMFFOptimizeMolecule

            rmol = copy.to_rdkit(keep_mapping=False, keep_hydrogens=False)
            ids = EmbedMultipleConfs(rmol, numConfs=limit, **kwargs)
            if optimize:
                for i in ids:
                    MMFFOptimizeMolecule(rmol, confId=i)

            conformers = [
                {n: tuple(v) for n, v in zip(self, conf.GetPositions())}
                for conf in rmol.GetConformers() if conf.Is3D()
            ]
        elif engine == 'cdpkit':
            from CDPL import Base, Chem, ConfGen
            from chython import SDFWrite, SDFRead  # to prevent circular imports
            from chython.files.mdl import parse_mol_v2000

            # the easiest way is just to provide intermediate SDF
            f = StringIO()
            SDFWrite(f, mapping=False).write(copy)
            cmol = Chem.BasicMolecule()
            if not Chem.SDFMoleculeReader(Base.StringIOStream(f.getvalue())).read(cmol):
                return 0

            ConfGen.prepareForConformerGeneration(cmol)
            gen = ConfGen.ConformerGenerator()
            gen.settings.timeout = kwargs.get('timeout', 60) * 1000
            gen.settings.minRMSD = kwargs.get('min_rmsd', .5)
            gen.settings.energyWindow = kwargs.get('energy_window', 20.)
            gen.settings.maxNumOutputConformers = limit
            if gen.generate(cmol) != ConfGen.ReturnCode.SUCCESS:
                return 0

            gen.setConformers(cmol)
            c = gen.getNumConformers()
            f = Base.StringIOStream(mode='w')
            Chem.SDFMolecularGraphWriter(f).write(cmol)
            s = SDFRead(StringIO(f.getvalue()))

            conformers = [
                {
                    n: (a['x'], a['y'], a['z'])
                    for n, a in zip(self, parse_mol_v2000(s._read_mol(current=False))['atoms'])
                }
                for _ in range(c)
            ]
        else: raise ValueError(f'Invalid conformer generation engine: {engine}')
        if conformers:
            self._conformers = conformers
        return len(conformers)


__all__ = ['Conformers']
