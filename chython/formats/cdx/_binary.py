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
"""CDX bytes to the object tree, every known property decoded into its CDXML spelling.

The stream, little-endian throughout:

| Item | Bytes |
|---|---|
| header | ``VjCD0100``, 4 reserved, 10 zero |
| object | tag (bit 15 set), INT32 id, then properties and objects, then ``00 00`` |
| property | tag (bit 15 clear), UINT16 length -- ``FFFF`` means a UINT32 length follows -- then data |

The reference library reads a 28-byte "header" that swallows the Document object's own tag and id;
this reader starts the object stream at byte 22, which is the same bytes.  Bit 14 set marks a
user-defined tag, skipped.  Coordinate order follows the reference library and real files, which for 3D
differ from the spec text:

| Type | CDX bytes | CDXML |
|---|---|---|
| CDXPoint2D | y, x | ``x y`` |
| CDXPoint3D | x, y, z (the spec says z, y, x) | ``x y z`` |
| CDXRectangle | top, left, bottom, right | ``left top right bottom`` |
"""

from struct import pack, unpack_from

from ._errors import MalformedCdx
from ._tree import CdxObject, MAX_DEPTH
from ._vocabulary import vocabulary
from ...core import LogRecord, LOST, REPAIRED


__all__ = ['parse_cdx', 'serialize_cdx', 'MAGIC', 'HEADER']


MAGIC = b'VjCD0100'
#: The bytes before the Document object: magic, 4 reserved (``04 03 02 01``), 10 zero.
HEADER = MAGIC + b'\x04\x03\x02\x01' + bytes(10)
#: One CDXCoordinate unit is 1/65536 point; CDXML states points.
SCALE = 65536.0

_SIGNED = {'INT8': True, 'INT16': True, 'INT32': True, 'UINT8': False, 'UINT16': False, 'UINT32': False}


def parse_cdx(data, log):
    """The tree of a CDX document."""
    if not isinstance(data, (bytes, bytearray, memoryview)):
        raise TypeError('a CDX document is bytes; for CDXML text use cdxml()')
    data = bytes(data)
    if not data.startswith(MAGIC):
        raise MalformedCdx('not a CDX document: the first 8 bytes are not VjCD0100')
    if len(data) < len(HEADER) + 6:
        raise MalformedCdx(f'CDX document of {len(data)} bytes holds no Document object')
    reader = _Reader(data, log)
    pos = len(HEADER)
    (tag,) = unpack_from('<H', data, pos)
    if tag != 0x8000:
        log.append(LogRecord('cdx:no-document-object', (),
                             f'byte {pos}: object tag 0x{tag:04X} where the Document object 0x8000 belongs; '
                             f'read as the document', REPAIRED))
    root, pos = reader.object(pos, 0)
    if root is None:
        raise MalformedCdx('CDX document holds no readable Document object')
    root.kind = 'CDXML'
    if reader.user_tags:
        log.append(LogRecord('cdx:user-defined-tags', (),
                             f'{reader.user_tags} user-defined property or object tag(s) skipped: they are '
                             f'private to the program that wrote them', LOST))
    return root


class _Reader:
    __slots__ = ('data', 'log', 'table', 'user_tags')

    def __init__(self, data, log):
        self.data = data
        self.log = log
        self.table = vocabulary()
        self.user_tags = 0

    def object(self, pos, depth):
        """Read the object whose tag is at `pos`.  Returns ``(object or None, position after it)``."""
        data = self.data
        end = len(data)
        if pos + 6 > end:
            self.log.append(LogRecord('cdx:truncated', (), f'byte {pos}: stream ends inside an object header', LOST))
            return None, end
        tag, oid = unpack_from('<Hi', data, pos)
        pos += 6
        if depth >= MAX_DEPTH:
            raise MalformedCdx(f'objects nest deeper than {MAX_DEPTH}')
        user = bool(tag & 0x4000)
        kind = self.table.objects.get(tag)
        if kind is None and not user:
            kind = f'0x{tag:04X}'
        node = CdxObject(kind or '', oid)
        while True:
            if pos + 2 > end:
                self.log.append(LogRecord('cdx:truncated', (),
                                          f'stream ends inside object 0x{tag:04X} id {oid}; read so far kept', LOST))
                return (None if user else node), end
            (child,) = unpack_from('<H', data, pos)
            if child == 0:
                pos += 2
                break
            if child & 0x8000:
                inner, pos = self.object(pos, depth + 1)
                if inner is not None:
                    node.children.append(inner)
                continue
            pos = self.property(pos, node)
        if user:
            self.user_tags += 1
            return None, pos
        return node, pos

    def property(self, pos, node):
        data = self.data
        end = len(data)
        if pos + 4 > end:
            self.log.append(LogRecord('cdx:truncated', (), f'byte {pos}: stream ends inside a property header', LOST))
            return end
        tag, length = unpack_from('<HH', data, pos)
        pos += 4
        if length == 0xFFFF:
            if pos + 4 > end:
                return end
            (length,) = unpack_from('<I', data, pos)
            pos += 4
        value = data[pos:pos + length]
        if len(value) < length:
            self.log.append(LogRecord('cdx:truncated', (),
                                      f'property 0x{tag:04X} claims {length} bytes, {len(value)} remain', LOST))
        pos += length
        if tag & 0x4000:
            self.user_tags += 1
            return pos
        prop = self.table.properties.get(tag)
        if prop is None:
            node.raw[f'0x{tag:04X}'] = value.hex()
            return pos
        if prop.name == 'kCDXProp_Text':
            node.text = _string(value)
            return pos
        if not prop.xml:
            node.raw[prop.name] = value.hex()
            return pos
        try:
            spelt = self.spell(prop, value)
        except (ValueError, IndexError) as e:
            self.log.append(LogRecord('cdx:bad-property', (),
                                      f'{prop.xml} on object id {node.id}: {e}; dropped', LOST))
            return pos
        if spelt is None:
            node.raw[prop.name] = value.hex()
        else:
            node.attributes[prop.xml] = spelt
        return pos

    def spell(self, prop, value):
        """`value` in the CDXML spelling of `prop`, or ``None`` for a type kept raw."""
        kind = prop.type
        n = len(value)
        if prop.name == 'kCDXProp_Bond_Order':
            return _bond_order(_integer(value, True))
        if prop.name == 'kCDXProp_Atom_Charge' and n == 4:
            # ChemDraw 8.0 wrote the charge as a 4-byte Fixed; the spec's reading rule.
            return str(unpack_from('<i', value)[0] >> 24)
        if prop.enum:
            number = _integer(value, _SIGNED.get(kind, True))
            names = self.table.enums.get(prop.enum, {})
            return names.get(number, str(number))
        if kind in _SIGNED:
            return str(_integer(value, _SIGNED[kind]))
        if kind == 'int':
            return str(_integer(value, True))
        if kind == 'bool':
            return 'yes' if not n or value[0] else 'no'
        if kind == 'string':
            return _text(value)
        if kind == 'FLOAT64':
            return _number(unpack_from('<d', value)[0])
        if kind == 'CDXCoordinate':
            return _number(_integer(value, True) / SCALE)
        if kind == 'CDXPoint2D':
            y, x = unpack_from('<ii', value)
            return f'{_number(x / SCALE)} {_number(y / SCALE)}'
        if kind == 'CDXPoint3D':
            x, y, z = unpack_from('<iii', value)
            return f'{_number(x / SCALE)} {_number(y / SCALE)} {_number(z / SCALE)}'
        if kind == 'CDXRectangle':
            top, left, bottom, right = unpack_from('<iiii', value)
            return ' '.join(_number(v / SCALE) for v in (left, top, right, bottom))
        if kind == 'CDXBoolean':
            return 'yes' if n and value[0] else 'no'
        if kind == 'CDXBooleanImplied':
            return 'yes' if not n or value[0] else 'no'
        if kind == 'CDXObjectID':
            return str(_integer(value, True))
        if kind == 'CDXObjectIDArray':
            return ' '.join(str(v) for v in unpack_from(f'<{n // 4}i', value))
        if kind == 'CDXObjectIDArrayWithCounts':
            (count,) = unpack_from('<H', value)
            return ' '.join(str(v) for v in unpack_from(f'<{count}i', value, 2))
        if kind == 'CDXString':
            return _string(value)
        if kind == 'CDXElementList':
            (count,) = unpack_from('<h', value)
            elements = unpack_from(f'<{abs(count)}h', value, 2)
            return ('NOT ' if count < 0 else '') + ' '.join(str(e) for e in elements)
        return None


def _integer(value, signed):
    if not value:
        return 0
    return int.from_bytes(value[:8], 'little', signed=signed)


def _number(v):
    """A coordinate or float as CDXML writes it: up to four decimals, no trailing zeros."""
    text = f'{v:.4f}'.rstrip('0').rstrip('.')
    return '0' if text in ('', '-0') else text


def _string(value):
    """A CDXString: a UINT16 style-run count, 10 bytes per run, then the text.

    The reference library also writes some CDXString properties as bare text with no run count (a
    `chemicalproperty` Name, for one), so a count whose runs do not fit the value reads the value as text.
    """
    if len(value) < 2:
        return _text(value)
    (runs,) = unpack_from('<H', value)
    if 2 + 10 * runs > len(value):
        return _text(value)
    return _text(value[2 + 10 * runs:])


def _text(value):
    """Bare text: UTF-8 where it decodes, Windows-1252 otherwise -- older ChemDraw wrote the latter."""
    try:
        return value.decode('utf-8')
    except UnicodeDecodeError:
        return value.decode('cp1252', errors='replace')


#: CDXBondOrder bits in their CDXML spelling, low bit first.
_ORDER_BITS = ((0x0001, '1'), (0x0002, '2'), (0x0004, '3'), (0x0008, '4'), (0x0010, '5'), (0x0020, '6'),
               (0x0040, '0.5'), (0x0080, '1.5'), (0x0100, '2.5'), (0x0200, '3.5'), (0x0400, '4.5'),
               (0x0800, '5.5'), (0x1000, 'dative'), (0x2000, 'ionic'), (0x4000, 'hydrogen'),
               (0x8000, 'threecenter'))


def _bond_order(value):
    if value == -1 or value == 0xFFFF:
        return 'any'
    value &= 0xFFFF
    return ' '.join(name for bit, name in _ORDER_BITS if value & bit) or '1'


_WIDTH = {'INT8': 'b', 'UINT8': 'B', 'INT16': 'h', 'UINT16': 'H', 'INT32': 'i', 'UINT32': 'I', 'int': 'i'}
#: CDXML elements that are properties in CDX, or that CDX folds into their parent: written by neither
#: tag nor property here.  A colour table changes no structure, and a `<s>` run is the text it holds.
_NOT_OBJECTS = frozenset(('colortable', 'color', 'fonttable', 'font', 's'))


def serialize_cdx(root, log):
    """CDX bytes of a tree.  The inverse of :func:`parse_cdx` for every type it decodes."""
    table = vocabulary()
    out = bytearray(HEADER)
    skipped = {}
    _put_object(root, out, table, skipped, log)
    out += b'\x00\x00'
    for (kind, name), count in sorted(skipped.items()):
        log.append(LogRecord('cdx:attribute-unwritable', (),
                             f'{kind} attribute {name} on {count} object(s) has no CDX encoding here; not '
                             f'written', LOST))
    return bytes(out)


def _put_object(node, out, table, skipped, log):
    tag = 0x8000 if node.kind == 'CDXML' else table.elements.get(node.kind)
    if tag is None or tag < 0x8000:
        if node.kind not in _NOT_OBJECTS:
            skipped[(node.kind, '<element>')] = skipped.get((node.kind, '<element>'), 0) + 1
        return
    out += pack('<Hi', tag, node.id)
    for name, value in node.attributes.items():
        prop = table.attributes.get(name)
        data = None if prop is None else _encode(prop, value, table)
        if data is None:
            skipped[(node.kind, name)] = skipped.get((node.kind, name), 0) + 1
            continue
        _put_property(prop.tag, data, out)
    if node.text is not None:
        _put_property(0x0700, pack('<H', 0) + node.text.encode('utf-8'), out)
    for child in node.children:
        _put_object(child, out, table, skipped, log)
    out += b'\x00\x00'


def _put_property(tag, data, out):
    if len(data) < 0xFFFF:
        out += pack('<HH', tag, len(data))
    else:
        out += pack('<HHI', tag, 0xFFFF, len(data))
    out += data


def _encode(prop, value, table):
    """`value`, a CDXML spelling, as the bytes of `prop`; ``None`` for a type this module does not write."""
    kind = prop.type
    try:
        if prop.name == 'kCDXProp_Bond_Order':
            if value == 'any':
                return pack('<h', -1)
            bits = dict((name, bit) for bit, name in _ORDER_BITS)
            return pack('<H', sum(bits[t] for t in value.split()))
        if prop.enum:
            number = table.spellings.get(prop.enum, {}).get(value)
            if number is None:
                number = int(value)
            return pack('<' + _WIDTH.get(kind, 'h'), number)
        if kind in _WIDTH:
            return pack('<' + _WIDTH[kind], int(value))
        if kind == 'FLOAT64':
            return pack('<d', float(value))
        if kind == 'CDXCoordinate':
            return pack('<i', round(float(value) * SCALE))
        if kind == 'CDXPoint2D':
            x, y = (float(t) for t in value.split()[:2])
            return pack('<ii', round(y * SCALE), round(x * SCALE))
        if kind == 'CDXPoint3D':
            x, y, z = (float(t) for t in value.split()[:3])
            return pack('<iii', round(x * SCALE), round(y * SCALE), round(z * SCALE))
        if kind == 'CDXRectangle':
            left, top, right, bottom = (float(t) for t in value.split()[:4])
            return pack('<iiii', *(round(v * SCALE) for v in (top, left, bottom, right)))
        if kind == 'CDXBoolean':
            return pack('<b', value == 'yes')
        if kind in ('CDXBooleanImplied', 'bool'):
            return b'' if value == 'yes' else pack('<b', 0)
        if kind == 'CDXObjectID':
            return pack('<i', int(value))
        if kind == 'CDXObjectIDArray':
            ids = [int(t) for t in value.split()]
            return pack(f'<{len(ids)}i', *ids)
        if kind == 'CDXString':
            return pack('<H', 0) + value.encode('utf-8')
        if kind == 'string':
            return value.encode('utf-8')
    except (ValueError, KeyError, OverflowError):
        return None
    return None
