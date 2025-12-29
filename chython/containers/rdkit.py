# -*- coding: utf-8 -*-
#
#  Copyright 2019-2025 Ramil Nugmanov <nougmanoff@protonmail.com>
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
        """
        from rdkit.Chem import BondStereo, BondType, ChiralType

        _rdkit_bond_map = {BondType.SINGLE: 1, BondType.DOUBLE: 2, BondType.TRIPLE: 3, BondType.AROMATIC: 4,
                           BondType.ZERO: 8, BondType.UNSPECIFIED: 8, BondType.DATIVE: 8}
        _chiral_cw = ChiralType.CHI_TETRAHEDRAL_CW
        _chiral_ccw = ChiralType.CHI_TETRAHEDRAL_CCW
        _trans = BondStereo.STEREOE
        _cis = BondStereo.STEREOZ

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
            # set coordinates from the first rdkit conformer. usually it's 2d layout
            for (_, atom), (x, y, _) in zip(mol.atoms(), cs[0].GetPositions()):
                atom.xy = (x, y)

            conformers = []
            for c in cs:
                if c.Is3D():
                    conformers.append({n: tuple(v) for n, v in enumerate(c.GetPositions(), 1)})
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

        mol.fix_structure(recalculate_hydrogens=False)
        if tetrahedron_stereo or cis_trans_stereo:
            mol.fix_stereo()
        return mol

    def to_rdkit(self: 'MoleculeContainer', *, keep_mapping=True, keep_hydrogens=True):
        """
        Convert into RDKit molecule object

        :param keep_mapping: set atom numbers
        :param keep_hydrogens: set implicit hydrogens
        """
        from rdkit.Chem import (AssignStereochemistry, Atom, BondStereo, BondType, ChiralType,
                                Conformer, RWMol, SanitizeMol)

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
            if self.atom(n).atomic_symbol not in _inorganic:
                n, m = m, n  # fix direction of dative bond
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

        conf = Conformer()
        for n, a in self.atoms():
            conf.SetAtomPosition(mapping[n], (a.x, a.y, 0))
        conf.Set3D(False)
        mol.AddConformer(conf, assignId=True)

        if hasattr(self, '_conformers'):
            for c in self._conformers:
                conf = Conformer()
                for n, xyz in c.items():
                    conf.SetAtomPosition(mapping[n], xyz)
                mol.AddConformer(conf, assignId=True)

        SanitizeMol(mol)
        AssignStereochemistry(mol, flagPossibleStereoCenters=True, force=True)
        return mol


_inorganic = {'He', 'Ne', 'Ar', 'Kr', 'Xe', 'F', 'Cl', 'Br', 'I', 'C', 'N', 'O',
              'H', 'Si', 'P', 'S', 'Se', 'Ge', 'As', 'Sb', 'Te'}


__all__ = ['RDkit']
