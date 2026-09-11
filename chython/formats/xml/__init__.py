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
"""XML chemical formats: one hardened tokenizer, one table-driven engine, two dialects.

:mod:`._tree` parses untrusted bytes, :mod:`._dialect` drives the vocabulary tables, :mod:`._cml` and
:mod:`._mrv` are the dialects and :mod:`._errors` the failure kinds.  A dialect is chosen by XML
namespace, never by vendor or idiom -- Marvin also writes CML-namespaced files.  Every record lands in
:class:`~chython.formats.ctfile._ctab.Ctab`, the intermediate every MDL path goes through."""

from ._cml import CML, CML_NS, parse_cml, read_cml, write_cml, write_cml_element
from ._cml import record_from_molecule as cml_record_from_molecule
from ._dialect import (Dialect, Field, NotModelled, Record, Tags, apply_fields, dialect, dialects,
                       parse_xml_document, read_document, read_molecule, read_xml, register, sniff,
                       write_molecule)
from ._errors import ForbiddenXml, MalformedXml, UnsupportedXml, XmlError
from ._facade import cml, mrv
# `record_from_molecule` is named per dialect on purpose: MRV's `hydrogenCount` is the implicit count
# and CML's is the total, so a bare name here would silently be whichever import line came second.
from ._mrv import MRV, MRV_NS, parse_mrv, read_mrv, write_mrv, write_mrv_element
from ._mrv import record_from_molecule as mrv_record_from_molecule
from ._tree import MAX_DEPTH, available_engines, parse_xml


__all__ = ['cml', 'mrv',
           'parse_cml', 'read_cml', 'write_cml', 'write_cml_element', 'cml_record_from_molecule',
           'CML', 'CML_NS',
           'parse_mrv', 'read_mrv', 'write_mrv', 'write_mrv_element', 'mrv_record_from_molecule',
           'MRV', 'MRV_NS',
           'read_xml', 'parse_xml_document',
           'parse_xml', 'available_engines', 'MAX_DEPTH',
           'Dialect', 'Field', 'Tags', 'Record', 'NotModelled',
           'register', 'dialect', 'dialects', 'sniff',
           'apply_fields', 'read_document', 'read_molecule', 'write_molecule',
           'XmlError', 'MalformedXml', 'UnsupportedXml', 'ForbiddenXml']
