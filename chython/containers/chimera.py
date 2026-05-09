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
from .._java import get_cdk


class Chimera:
    __slots__ = ()

    def to_cdk(self):
        """
        Convert molecule to CDK Molecule object.

        Due to translation through SMILES string, atom order is not preserved.
        Use `self.smiles_atoms_order` to map atoms back.
        """
        cdk = get_cdk()
        parser = cdk.smiles.SmilesParser(cdk.DefaultChemObjectBuilder.getInstance())
        return parser.parseSmiles(str(self))

    def to_openbabel(self):
        """
        Convert molecule to OpenBabel OBMol object.

        Due to translation through SMILES string, atom order is not preserved.
        Use `self.smiles_atoms_order` to map atoms back.
        """
        from openbabel import openbabel

        conv = openbabel.OBConversion()
        conv.SetInFormat('smi')
        mol = openbabel.OBMol()
        assert conv.ReadString(mol, str(self)), 'OpenBabel failed to parse smiles'
        return mol

    def to_indigo(self):
        """
        Convert molecule to Indigo molecule object.

        Due to translation through SMILES string, atom order is not preserved.
        Use `self.smiles_atoms_order` to map atoms back.
        """
        from indigo import Indigo

        return Indigo().loadMolecule(str(self))


__all__ = ['Chimera']
