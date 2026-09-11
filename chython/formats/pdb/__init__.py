# -*- coding: utf-8 -*-
#
#  Copyright 2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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
"""The PDB family: PDBx/mmCIF and legacy PDB, plus the pass that turns a record into a molecule.

Both readers return :class:`PDBRecord` rather than a :class:`~chython.core.MoleculeContainer`: the
record holds the z coordinate and the residue annotation the container has nowhere to put.  Every bond
they build is one the file states (``_chem_comp_bond``, ``_struct_conn``, ``CONECT``, ``SSBOND``,
``LINK``); no interatomic distance is ever computed here.  :func:`build_molecule` is the separately
invoked opt-in that consults ``chemistry/tables/residues.tsv`` -- a lookup on a stated component id.
"""
from ._builder import build_molecule
from ._legacy import pdb, read_pdb
from ._mmcif import mmcif, read_mmcif
from ._records import PDBAtom, PDBBond, PDBRecord
from ._star import INAPPLICABLE, UNKNOWN, StarBlock, StarLoop, is_null, parse_star


__all__ = ['PDBAtom', 'PDBBond', 'PDBRecord', 'build_molecule', 'mmcif', 'read_mmcif', 'pdb',
           'read_pdb',
           'INAPPLICABLE', 'UNKNOWN', 'StarBlock', 'StarLoop', 'is_null', 'parse_star']
