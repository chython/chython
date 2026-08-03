# -*- coding: utf-8 -*-
#
#  Copyright 2025, 2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from typing import Literal


class Conformers:
    __slots__ = ()

    def generate_conformers(self, limit: int = 10, *, optimize: bool = False,
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
            from CDPL import Chem, ConfGen

            # build the CDPKit molecule directly from the graph. an SDF proxy loses
            # chirality: chython has no 2D layout, so the wedge bonds it writes are
            # ambiguous and CDPKit generates arbitrary handedness. setting stereo
            # descriptors explicitly makes conformers stereochemically correct.
            cmol = Chem.BasicMolecule()
            pos = {}  # chython atom number -> CDPKit atom index (0-based), atom order preserved
            for i, (n, a) in enumerate(copy.atoms()):
                ca = cmol.addAtom()
                Chem.setType(ca, a.atomic_number)
                Chem.setFormalCharge(ca, a.charge)
                Chem.setImplicitHydrogenCount(ca, 0)  # hydrogens are explicit
                if a.isotope:
                    Chem.setIsotope(ca, a.isotope)
                if a.is_radical:
                    Chem.setRadicalType(ca, Chem.RadicalType.DOUBLET)
                    Chem.setUnpairedElectronCount(ca, 1)
                pos[n] = i

            bond_of = {}  # (n, m) and (m, n) -> CDPKit bond, for cis/trans reference lookup
            for n, m, b in copy.bonds():
                cb = cmol.addBond(pos[n], pos[m])
                Chem.setOrder(cb, 1 if b.order == 4 else b.order)
                bond_of[n, m] = bond_of[m, n] = cb
                if b.order == 4:
                    Chem.setAromaticityFlag(cb, True)
            Chem.calcBasicProperties(cmol, False)

            for n, a in copy.atoms():
                if a.stereo is None or n not in copy.stereogenic_tetrahedrons:
                    continue  # allenes are not supported
                env = list(copy._bonds[n])
                s = copy._translate_tetrahedron_sign(n, env)
                refs = [cmol.getAtom(pos[x]) for x in env]
                if len(refs) < 4:  # implicit hydrogen is encoded as the central atom itself
                    refs.append(cmol.getAtom(pos[n]))
                cfg = Chem.AtomConfiguration.S if s else Chem.AtomConfiguration.R
                Chem.setStereoDescriptor(cmol.getAtom(pos[n]), Chem.StereoDescriptor(cfg, *refs))

            for n, m, b in copy.bonds():
                if b.stereo is None:
                    continue
                nm = copy._stereo_cis_trans_centers.get(n)
                if nm is None or n not in nm or m not in nm:
                    continue
                n1, m1, *_ = copy.stereogenic_cis_trans[nm]
                cfg = Chem.BondConfiguration.CIS if b.stereo else Chem.BondConfiguration.TRANS  # True -> cis/Z
                Chem.setStereoDescriptor(bond_of[n, m], Chem.StereoDescriptor(
                    cfg, cmol.getAtom(pos[n1]), cmol.getAtom(pos[n]), cmol.getAtom(pos[m]), cmol.getAtom(pos[m1])))

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
            # read coordinates directly per conformer; zip(self, ...) keeps only heavy atoms
            atom_of = {n: cmol.getAtom(pos[n]) for n in self}
            conformers = [
                {n: tuple(Chem.getConformer3DCoordinates(atom_of[n], i)) for n in self}
                for i in range(c)
            ]
        else: raise ValueError(f'Invalid conformer generation engine: {engine}')
        if conformers:
            self._conformers = conformers
        return len(conformers)


__all__ = ['Conformers']
