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
from CachedMethods import class_cached_property


class Chimera:
    __slots__ = ()

    def to_cdk(self):
        """
        Convert molecule to CDK Molecule object.

        Due to translation through SMILES string, atom order is not preserved.
        Use `self.smiles_atoms_order` to map atoms back.
        """
        parser = self._cdk_engine.smiles.SmilesParser(self._cdk_engine.DefaultChemObjectBuilder.getInstance())
        return parser.parseSmiles(str(self))

    def to_openbabel(self):
        """
        Convert molecule to OpenBabel OBMol object.

        Due to translation through SMILES string, atom order is not preserved.
        Use `self.smiles_atoms_order` to map atoms back.
        """
        from openbabel import openbabel

        mol = openbabel.OBMol()
        assert self._obparser(mol, str(self)), 'OpenBabel failed to parse smiles'
        return mol

    def to_indigo(self):
        """
        Convert molecule to Indigo molecule object.

        Due to translation through SMILES string, atom order is not preserved.
        Use `self.smiles_atoms_order` to map atoms back.
        """
        return self._indigo_engine.loadMolecule(str(self))

    @class_cached_property
    def _cdk_engine(self):
        try:
            from jpype import isJVMStarted, startJVM, JPackage

            if not isJVMStarted():
                from chython import class_paths

                startJVM('--enable-native-access=ALL-UNNAMED', classpath=class_paths)

            return JPackage('org').openscience.cdk
        except (ImportError, AttributeError):
            raise ImportError('Java/JPype/CDK.jar is not installed or broken. make sure CDK_PATH env variable is set')

    @class_cached_property
    def _indigo_engine(self):
        from indigo import Indigo

        return Indigo()

    @class_cached_property
    def _obparser(self):
        from openbabel import openbabel

        obparser = openbabel.OBConversion()
        obparser.SetInFormat('smi')
        return obparser.ReadString

    @class_cached_property
    def _obgen2d(self):
        from openbabel import openbabel

        return openbabel.OBOp.FindType('gen2D').Do


__all__ = ['Chimera']
