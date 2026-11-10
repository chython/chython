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
"""The one CDX object tree, which both spellings tokenize into and serialize from.

Every attribute value is a string **in its CDXML spelling** -- ``Order="1.5"``, ``p="72 144"``,
``Display="WedgeBegin"`` -- whichever spelling the document arrived in, so one builder reads both.
"""


__all__ = ['CdxObject', 'MAX_DEPTH']


#: How deeply objects may nest.  The XML tokenizer's own cap, applied to the binary stream as well.
MAX_DEPTH = 200


class CdxObject:
    """One CDX object: an element of a CDXML document or a tagged object of a CDX stream.

    `kind` is the CDXML element name (``'n'``, ``'b'``, ``'fragment'``); `id` is the object id, 0 when
    the document assigned none; `attributes` maps CDXML attribute names to CDXML-spelt values; `text`
    is a text object's string, its style runs joined.  `raw` holds binary properties the CDXML table
    does not name, under their constant name, as hex -- so a log line can name what was not read.
    """
    __slots__ = ('kind', 'id', 'attributes', 'children', 'text', 'raw')

    def __init__(self, kind, id=0, attributes=None, children=None, text=None):
        self.kind = kind
        self.id = id
        self.attributes = {} if attributes is None else attributes
        self.children = [] if children is None else children
        self.text = text
        self.raw = {}

    def walk(self):
        """This object and every descendant, depth first, in document order."""
        stack = [self]
        while stack:
            node = stack.pop()
            yield node
            stack.extend(reversed(node.children))

    def __repr__(self):
        return f'CdxObject({self.kind}, id={self.id}, {len(self.attributes)} attrs, {len(self.children)} children)'
