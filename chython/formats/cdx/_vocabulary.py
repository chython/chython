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
"""The CDX vocabulary: object tags, property tags and enum spellings, from ``vocabulary.tsv``.

Read lazily on first use and never at import.  One table serves both spellings: a binary tag resolves
to the CDXML element or attribute name, so the binary tokenizer hands the builder the same tree the
CDXML tokenizer does.
"""

from pathlib import Path
from typing import NamedTuple


__all__ = ['Property', 'Vocabulary', 'vocabulary']


class Property(NamedTuple):
    """One CDX property: its tag, the constant name, the CDXML attribute and the CDX data type.

    `xml` is empty for a property the CDXML table does not name; the binary tokenizer then files the
    value under `name`, so a log line can still say what was dropped.  `enum` is the enum type whose
    spellings `Vocabulary.enums` holds, empty for a property that is not enumerated.
    """
    tag: int
    name: str
    xml: str
    type: str
    enum: str


class Vocabulary(NamedTuple):
    objects: dict          # tag -> CDXML element name
    elements: dict         # CDXML element name -> tag
    properties: dict       # tag -> Property
    attributes: dict       # CDXML attribute name -> Property
    enums: dict            # enum type -> {value: spelling}
    spellings: dict        # enum type -> {spelling: value}


_VOCABULARY_CACHE = {}


def vocabulary() -> Vocabulary:
    """The compiled table, built once."""
    if _VOCABULARY_CACHE:
        return _VOCABULARY_CACHE['table']
    objects, elements, properties, attributes, enums, spellings = {}, {}, {}, {}, {}, {}
    with (Path(__file__).parent / 'vocabulary.tsv').open(encoding='utf-8') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            row = line.rstrip('\n').split('\t')
            kind = row[0]
            if kind == 'object':
                tag = int(row[1], 16)
                objects[tag] = row[3]
                elements.setdefault(row[3], tag)
            elif kind == 'property':
                row += [''] * (6 - len(row))
                p = Property(int(row[1], 16), row[2], row[3], row[4], row[5])
                properties[p.tag] = p
                if p.xml:
                    attributes.setdefault(p.xml, p)
            elif kind == 'enum':
                value, etype, spelling = int(row[1]), row[2], row[3]
                enums.setdefault(etype, {}).setdefault(value, spelling)
                spellings.setdefault(etype, {}).setdefault(spelling, value)
    table = Vocabulary(objects, elements, properties, attributes, enums, spellings)
    _VOCABULARY_CACHE['table'] = table
    return table
