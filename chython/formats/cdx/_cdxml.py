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
"""CDXML text to the object tree.

Parsing goes through :func:`chython.formats.xml._tree.parse_xml`, so the entity policy and the depth cap
are the XML family's.  A CDXML document always declares the external ``CDXML.dtd``; that declaration is
dropped without a log line, being the format's signature rather than anything the document lost.
"""

from ._errors import MalformedCdx
from ._tree import CdxObject
from ..xml._errors import XmlError
from ..xml._tree import local, parse_xml
from ...core import LogRecord, REPAIRED


__all__ = ['parse_cdxml', 'serialize_cdxml']

#: The declaration every CDXML document opens with; no reader fetches it.
DOCTYPE = '<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd" >'


def parse_cdxml(source, log):
    """The tree of a CDXML document.  `source` is text, bytes, a path or an open file."""
    own = []
    try:
        root = parse_xml(source, log=own)
    except XmlError as e:
        raise MalformedCdx(f'not a well-formed CDXML document: {e}') from e
    log.extend(r for r in own if r.rule != 'xml:doctype-dropped')
    if local(root.tag) != 'CDXML':
        raise MalformedCdx(f'root element is <{local(root.tag)}>, a CDXML document has <CDXML>')
    return _convert(root, log)


def _convert(element, log):
    kind = local(element.tag)
    attributes = dict(element.attrib)
    raw_id = attributes.pop('id', None)
    oid = 0
    if raw_id is not None:
        try:
            oid = int(raw_id)
        except ValueError:
            log.append(LogRecord('cdx:bad-id', (),
                                 f'<{kind}> id {raw_id!r:.20} is not an integer; read as unassigned', REPAIRED))
    node = CdxObject(kind, oid, attributes)
    if kind == 't':
        # A text object's string is its <s> style runs in order; the runs' font and size are layout.
        node.text = ''.join(''.join(s.itertext()) for s in element if local(s.tag) == 's')
        node.children = [_convert(child, log) for child in element if local(child.tag) != 's']
        return node
    node.children = [_convert(child, log) for child in element]
    return node


def serialize_cdxml(root, log):
    """CDXML text of a tree, one element per line."""
    lines = ['<?xml version="1.0" encoding="UTF-8" ?>', DOCTYPE]
    _emit(root, lines, 0)
    return '\n'.join(lines) + '\n'


def _emit(node, lines, depth):
    pad = ' ' * depth
    head = node.kind
    if node.id:
        head += f' id="{node.id}"'
    for name, value in node.attributes.items():
        head += f' {name}="{_escape(value)}"'
    if node.kind == 't' and node.text is not None:
        lines.append(f'{pad}<{head}><s>{_escape(node.text)}</s></t>')
        return
    if not node.children:
        lines.append(f'{pad}<{head}/>')
        return
    lines.append(f'{pad}<{head}>')
    for child in node.children:
        _emit(child, lines, depth + 1)
    lines.append(f'{pad}</{node.kind}>')


def _escape(value):
    return (str(value).replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
            .replace('"', '&quot;'))
