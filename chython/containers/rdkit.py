# -*- coding: utf-8 -*-
#
#  Copyright 2019-2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from typing import TYPE_CHECKING, Type
from ..periodictable import Element


if TYPE_CHECKING:
    from chython import MoleculeContainer


class RDkit:
    __slots__ = ()

    @classmethod
    def from_rdkit(cls: Type['MoleculeContainer'], data):
        """
        RDKit molecule object to MoleculeContainer converter

        AND/OR stereo groups are imported as extended stereo. Absolute groups are dropped:
        chython treats any set stereo label outside an AND/OR group as absolute already.
        """
        from rdkit.Chem import BondStereo, BondType, ChiralType, StereoGroupType

        _rdkit_bond_map = {BondType.SINGLE: 1, BondType.DOUBLE: 2, BondType.TRIPLE: 3, BondType.AROMATIC: 4,
                           BondType.ZERO: 8, BondType.UNSPECIFIED: 8, BondType.DATIVE: 8}
        _chiral_cw = ChiralType.CHI_TETRAHEDRAL_CW
        _chiral_ccw = ChiralType.CHI_TETRAHEDRAL_CCW
        _trans = BondStereo.STEREOE
        _cis = BondStereo.STEREOZ
        _and = StereoGroupType.STEREO_AND
        _or = StereoGroupType.STEREO_OR

        mol = cls()

        mapping = {}
        tetrahedron_stereo = []
        for ra in data.GetAtoms():
            e = Element.from_symbol(ra.GetSymbol())
            a = e(ra.GetIsotope() or None, charge=ra.GetFormalCharge(), is_radical=bool(ra.GetNumRadicalElectrons()),
                  parsed_mapping=ra.GetAtomMapNum(), implicit_hydrogens=ra.GetNumExplicitHs() + ra.GetNumImplicitHs())
            mapping[ra.GetIdx()] = mol.add_atom(a, _skip_calculation=True)
            s = ra.GetChiralTag()
            if s in (_chiral_cw, _chiral_ccw):
                tetrahedron_stereo.append((ra.GetIdx(), [x.GetIdx() for x in ra.GetNeighbors()], s == _chiral_ccw))

        cis_trans_stereo = []
        for b in data.GetBonds():
            n, m = mapping[b.GetBeginAtomIdx()], mapping[b.GetEndAtomIdx()]
            mol.add_bond(n, m, _rdkit_bond_map[b.GetBondType()], _skip_calculation=True)
            s = b.GetStereo()
            if s in (_cis, _trans):
                nn, nm = b.GetStereoAtoms()
                cis_trans_stereo.append((n, m, mapping[nn], mapping[nm], s == _cis))

        if cs := data.GetConformers():
            conformers = []
            plane = None
            for c in cs:
                if c.Is3D():
                    conformers.append({n: tuple(v) for n, v in enumerate(c.GetPositions(), 1)})
                elif plane is None:
                    plane = c  # first 2d conformer is the layout
            if plane is not None:
                # the xy of a 3d conformer is a meaningless projection, not a layout
                for (_, atom), (x, y, _) in zip(mol.atoms(), plane.GetPositions()):
                    atom.xy = (x, y)
            if conformers:
                mol._conformers = conformers

        # move stereo labels as is
        for n, env, s in tetrahedron_stereo:
            n = mapping[n]
            try:
                mol.atom(n)._stereo = mol._translate_tetrahedron_sign(n, [mapping[x] for x in env], s)
            except KeyError:
                pass
        for n, m, nn, nm, s in cis_trans_stereo:
            try:
                mol.bond(n, m)._stereo = mol._translate_cis_trans_sign(n, m, nn, nm, s)
            except KeyError:
                pass

        # AND groups are positive, OR groups negative. group ids are taken as is
        for sg in data.GetStereoGroups():
            if (gt := sg.GetGroupType()) is _and:
                sign = 1
            elif gt is _or:
                sign = -1
            else:  # absolute group. chython has no explicit mark for it
                continue
            # bond-level extended stereo has no place in the chython model, atoms only
            es = sign * (sg.GetReadId() or 1)
            for ra in sg.GetAtoms():
                a = mol.atom(mapping[ra.GetIdx()])
                if a._extended_stereo is None:  # rdkit allows an atom in AND and OR at once. first wins
                    a._extended_stereo = es

        mol.fix_structure(recalculate_hydrogens=False)
        if tetrahedron_stereo or cis_trans_stereo:
            mol.fix_stereo()
        return mol

    def to_rdkit(self: 'MoleculeContainer', *, keep_mapping=True, keep_hydrogens=True, keep_coordinates=None,
                 absolute=False):
        """
        Convert into RDKit molecule object

        Aromatic, tetrahedral, cis-trans and extended (AND/OR) stereo are preserved.
        Allene/cumulene stereo is not exported: RDKit cannot round-trip it (the 3D
        embedder ignores allene chirality), mirroring `to_openbabel` and `to_indigo`.

        :param keep_mapping: set atom numbers
        :param keep_hydrogens: set implicit hydrogens
        :param keep_coordinates: export the 2D layout as an RDKit conformer.
         `None` (default) exports it only when a layout exists, i.e. some atom has nonzero coordinates.
        :param absolute: add an ABS stereo group for stereocenters without extended stereo groups.
            By default chython treats any set chirality flag as absolute if not in an AND/OR group,
            but some tools require the ABS group to be explicitly present. Mirrors the `absolute`
            flag of the MDL V3000 writers.
        """
        from rdkit.Chem import (AssignStereochemistry, Atom, BondStereo, BondType, ChiralType, Conformer,
                                CreateStereoGroup, RWMol, SanitizeMol, SetDoubleBondNeighborDirections,
                                StereoGroupType)

        _bond_map = {1: BondType.SINGLE, 2: BondType.DOUBLE, 3: BondType.TRIPLE,
                     4: BondType.AROMATIC, 8: BondType.DATIVE}

        _chiral_cw = ChiralType.CHI_TETRAHEDRAL_CW
        _chiral_ccw = ChiralType.CHI_TETRAHEDRAL_CCW
        _trans = BondStereo.STEREOE
        _cis = BondStereo.STEREOZ

        mol = RWMol()
        mapping = {}

        for n, a in self.atoms():
            ra = Atom(a.atomic_number)
            if keep_hydrogens and a.implicit_hydrogens is not None:
                ra.SetNumExplicitHs(a.implicit_hydrogens)
            if keep_mapping:
                ra.SetAtomMapNum(n)
            if a.charge:
                ra.SetFormalCharge(a.charge)
            if a.isotope:
                ra.SetIsotope(a.isotope)
            if a.is_radical:
                ra.SetNumRadicalElectrons(1)
            mapping[n] = mol.AddAtom(ra)

        inverted = {v: k for k, v in mapping.items()}

        for n, m, b in self.bonds():
            if b.order == 8 and self.atom(n).atomic_symbol not in _inorganic:
                n, m = m, n  # dative bond points from the donor to the acceptor
            mol.AddBond(mapping[n], mapping[m], _bond_map[b.order])

        for n, a in self.atoms():
            if a.stereo is None:
                continue
            if n not in self.stereogenic_tetrahedrons:
                continue  # allenes are not supported
            ra = mol.GetAtomWithIdx(mapping[n])
            env = [inverted[x.GetIdx()] for x in ra.GetNeighbors()]
            s = self._translate_tetrahedron_sign(n, env)
            ra.SetChiralTag(_chiral_ccw if s else _chiral_cw)

        for n, m, b in self.bonds():
            if b.stereo is None:
                continue
            # check for simple cis-trans
            nm = self._stereo_cis_trans_centers.get(n)
            if nm is None or n not in nm or m not in nm:
                continue

            n1, m1, *_ = self.stereogenic_cis_trans[nm]
            rb = mol.GetBondBetweenAtoms(mapping[n], mapping[m])
            rb.SetStereoAtoms(mapping[n1], mapping[m1])
            rb.SetStereo(_cis if b.stereo else _trans)

        # extended stereo groups. only tetrahedrons, matching the exported chirality tags
        rac = defaultdict(list)  # AND groups: positive extended_stereo
        rel = defaultdict(list)  # OR groups: negative extended_stereo
        ast = []  # absolute stereo
        for n, a in self.atoms():
            if a.stereo is None or n not in self.stereogenic_tetrahedrons:
                continue
            if (es := a.extended_stereo) is None:
                ast.append(mapping[n])
            elif es > 0:
                rac[es].append(mapping[n])
            else:
                rel[-es].append(mapping[n])

        groups = []
        if absolute and ast:
            groups.append((StereoGroupType.STEREO_ABSOLUTE, 0, ast))
        groups.extend((StereoGroupType.STEREO_AND, gid, rac[gid]) for gid in sorted(rac))
        groups.extend((StereoGroupType.STEREO_OR, gid, rel[gid]) for gid in sorted(rel))
        if groups:
            sgs = []
            for gt, gid, al in groups:
                sg = CreateStereoGroup(gt, mol, al, [], gid)
                sg.SetWriteId(gid)  # otherwise rdkit renumbers the groups from one on write
                sgs.append(sg)
            mol.SetStereoGroups(sgs)

        if keep_coordinates is None:
            keep_coordinates = any(a.x or a.y for _, a in self.atoms())
        if keep_coordinates:
            conf = Conformer()
            for n, a in self.atoms():
                conf.SetAtomPosition(mapping[n], (a.x, a.y, 0))
            conf.Set3D(False)
            mol.AddConformer(conf, assignId=True)

        for c in self._conformers or ():
            conf = Conformer()
            for n, xyz in c.items():
                conf.SetAtomPosition(mapping[n], xyz)
            mol.AddConformer(conf, assignId=True)

        SanitizeMol(mol)
        AssignStereochemistry(mol, flagPossibleStereoCenters=True, force=True)
        # rdkit derives the smiles bond directions from the stereo atoms only for single fragment
        # molecules. without this any salt or solvate loses its cis-trans marks on smiles export
        SetDoubleBondNeighborDirections(mol)
        return mol


_inorganic = {'He', 'Ne', 'Ar', 'Kr', 'Xe', 'F', 'Cl', 'Br', 'I', 'C', 'N', 'O',
              'H', 'Si', 'P', 'S', 'Se', 'Ge', 'As', 'Sb', 'Te'}


__all__ = ['RDkit']
