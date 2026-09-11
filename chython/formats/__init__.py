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
"""File formats: MDL V2000/V3000 with SDF and RDfile framing, Tripos MOL2, XYZ, PDBx/mmCIF, legacy
PDB and the XML family.  Strings -- SMILES, SMARTS, InChI, names -- are not files and live elsewhere.
A reader returns what the file said and normalises nothing; a writer neither repairs nor mutates its
argument; an unparsable record is stored with its error (``SDFRead.failed``) rather than raised past.
Unmodelled constructs round-trip, stored on the container so a molecule-only consumer cannot lose them.
"""

from .ctfile import (ERDFWrite, ESDFWrite, FIELDDISP_TAIL, FailedRecord, RDFRead, RDFWrite, SDFRead,
                     SDFWrite, SGroup, SGroupStore, add_data_sgroup, data_sgroups, mol, rxn)
from .mol2 import Mol2ParseError, mol2, mol2_mol, read_mol2
from .xml import (ForbiddenXml, MalformedXml, UnsupportedXml, XmlError, cml, mrv, read_cml, read_mrv,
                  read_xml, write_cml, write_mrv)
from .xyz import XYZAtom, XYZFrame, xyz, xyz_conformers
from ..core._core import _set_sgroup_fns


# The one registration this package makes, and it goes the same way `chemistry`'s do: the core owns
# S-group storage and the method names, this package owns what a CTfile `DAT` record means.  See
# `_set_sgroup_fns`.
_set_sgroup_fns(add_data_sgroup=add_data_sgroup, data_sgroups=data_sgroups)

# `xyz` and the PDB-family entry points return their own intermediates, and that is the interface: XYZ
# states no bonds, both state a z coordinate the container cannot yet hold, and a PDB record carries a
# residue annotation too.  Turning one into chemistry is an explicit later call: `build_molecule` in the
# subpackage that read the record -- `xyz.build_molecule` or `pdb.build_molecule` -- then
# `chython.chemistry.perceive_bonds` where the format stated no bond, then `chython.chemistry.saturate`.
# One facade name cannot serve two record types, so `chython.build_molecule` is the PDB record's and the
# XYZ builder is reached by its full path.

# `pdb` and `xml` are kept out of `__all__`: the facade star-imports these names, so exporting either
# would make `chython.formats.pdb` mean the function -- hiding the STAR/CIF tokeniser, reachable only
# there -- and `chython.xml` read like `xml.etree`.  The facade imports the four PDB entry points by full
# path, pinned by `chython/test/test_facade_names.py`.  `xyz` may shadow its module, `xyz.py`
# re-exporting all of its own names.

# An intermediate is exported only when it is the *only* thing a format returns: hence `XYZFrame` and
# `PDBRecord`, but not `.ctfile`'s `Ctab` nor `.xml`'s `Record`, since those readers hand back molecules.
# `read_xml` sniffs the XML dialect; `mol()` deliberately does not, one entry point per question.

# `SGroup` and `SGroupStore` are on the facade because `add_data_sgroup` hands one back and a caller
# reads it: the alternative is `mol.sgroups`, which is the arena's dicts and the core's alphabet.
__all__ = ['SDFRead', 'SDFWrite', 'ESDFWrite', 'RDFRead', 'RDFWrite', 'ERDFWrite',
           'FailedRecord', 'mol', 'rxn',
           'SGroup', 'SGroupStore', 'add_data_sgroup', 'data_sgroups', 'FIELDDISP_TAIL',
           # `mol2` shadows its own module inside this package, as `xyz` already does: the facade is the
           # name a caller wants, and `from .mol2 import ...` here is by full path anyway.
           'mol2', 'read_mol2', 'mol2_mol', 'Mol2ParseError',
           'cml', 'mrv', 'read_cml', 'write_cml', 'read_mrv', 'write_mrv', 'read_xml',
           'XmlError', 'MalformedXml', 'UnsupportedXml', 'ForbiddenXml',
           'xyz', 'xyz_conformers', 'XYZFrame', 'XYZAtom']
