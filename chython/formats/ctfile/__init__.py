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
"""CTfile (MDL) reading and writing: V2000/V3000 CTAB, SDF, RXN and RDfile.

:mod:`._v2000` and :mod:`._v3000` are the only modules that know a column offset or a ``M  V30``
keyword; :mod:`._ctab` is the shared intermediate every path goes through.  Wedge sign conventions
live below, in :mod:`chython.core.wedge`.  Input is stored and logged, never rejected for being
chemically wrong; writing stays inside the specification.
"""

from ._ctab import Ctab, CtabAtom, CtabBond
from ._errors import CtfileError, MalformedCtfile, UnsupportedCtfile
from ._facade import mol, needs_v3000, rxn
from ._rdf import (ERDFWrite, RDF_HEADER, RDFRead, RDFWrite, parse_rdf_fields, parse_rdf_record,
                   split_rdf_records)
from ._rxn import (RXN_HEADER_LINES, emit_rxn, parse_rxn, parse_rxn_record, sniff_rxn_version,
                   split_rxn_v2000)
from ._sdf import (UNPARSED_KEY, V2000_STAMP, V3000_STAMP, emit_record, parse_record,
                   split_records)
from ._sgroup import FIELDDISP_TAIL, SGroup, SGroupStore, add_data_sgroup, data_sgroups
from ._stream import ESDFWrite, FailedRecord, SDFRead, SDFWrite
from ._v2000 import emit_v2000, parse_v2000
from ._v3000 import emit_v3000, parse_v3000


__all__ = ['SDFRead', 'SDFWrite', 'ESDFWrite', 'FailedRecord',
           'mol', 'rxn', 'needs_v3000',
           'Ctab', 'CtabAtom', 'CtabBond',
           'CtfileError', 'MalformedCtfile', 'UnsupportedCtfile',
           'SGroup', 'SGroupStore', 'add_data_sgroup', 'data_sgroups', 'FIELDDISP_TAIL',
           'V2000_STAMP', 'V3000_STAMP', 'UNPARSED_KEY',
           'parse_record', 'emit_record', 'split_records',
           'parse_v2000', 'emit_v2000', 'parse_v3000', 'emit_v3000',
           'emit_rxn', 'parse_rxn', 'parse_rxn_record', 'sniff_rxn_version', 'split_rxn_v2000',
           'RXN_HEADER_LINES',
           'RDF_HEADER', 'split_rdf_records', 'parse_rdf_fields', 'parse_rdf_record',
           'RDFRead', 'RDFWrite', 'ERDFWrite']
