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
"""The one place this library turns bytes into an element tree, and the entity policy that makes that
safe.  :func:`parse_xml` applies both rules itself, whichever backend parses, so an optional
dependency cannot change what a reader accepts: (1) no entity declarations -- a ``<!DOCTYPE`` with an
internal subset is refused, an external-only one dropped and logged, since nothing here fetches a URL;
(2) a depth cap, nesting being what turns a small document into deep recursion downstream."""

from os import PathLike
from xml.etree.ElementTree import Element, ParseError

from ._errors import ForbiddenXml, MalformedXml
from ...core import LogRecord, LOST


__all__ = ['parse_xml', 'available_engines', 'local', 'namespace_of', 'text_of', 'MAX_DEPTH',
           'ENGINES']


#: How deeply an element tree may nest before :func:`parse_xml` refuses it.  Chemical XML is shallow:
#: CML's deepest path is document/molecule/bondArray/bond/bondStereo, five, and MRV's is about eight.
MAX_DEPTH = 200

#: Backend names :func:`parse_xml` understands for its ``engine`` keyword, most preferred first.
#: ``None`` means "the first of these that imports", which is what every caller but the tests wants.
ENGINES = ('defusedxml', 'stdlib')

# name -> callable(data) -> Element.  Populated on first use, never at import, so an optional
# dependency is not imported by importing this package.
_BACKEND_CACHE = {}


def _backends():
    """``{name: parse function}`` for every backend available here, built once."""
    if _BACKEND_CACHE:
        return _BACKEND_CACHE
    try:
        from defusedxml.ElementTree import fromstring as defused
    except ImportError:
        pass
    else:
        def _defused(data, allow_dtd):
            # The three guards duplicate rule 1 as defence in depth, and they must track `allow_dtd`
            # rather than being pinned on: with all three off the two backends agree byte for byte,
            # including on an external DOCTYPE, since `forbid_external` guards entity *references* and
            # neither backend ever fetches a SYSTEM identifier.
            return defused(data, forbid_dtd=not allow_dtd, forbid_entities=not allow_dtd,
                           forbid_external=not allow_dtd)
        _BACKEND_CACHE['defusedxml'] = _defused

    from xml.etree.ElementTree import fromstring as stdlib
    # Safe *under rule 1* and only under it: with no DTD there is no entity for `fromstring` to
    # expand, and ElementTree never installs an external-entity handler.
    _BACKEND_CACHE['stdlib'] = lambda data, allow_dtd: stdlib(data)
    return _BACKEND_CACHE


def available_engines():
    """The backend names usable in this interpreter, most preferred first.

    ``'stdlib'`` is always in it; ``'defusedxml'`` when the optional dependency is installed.
    """
    have = _backends()
    return tuple(name for name in ENGINES if name in have)


# rule 1 -- the prolog scan

_DOCTYPE = '<!DOCTYPE'


def _doctype_span(text, start):
    """``(has_internal_subset, end)`` for the DOCTYPE declaration beginning at `start`.

    `end` is one past the closing ``>``, meaningful only when there is no internal subset.  An
    unterminated declaration counts as having one, which is the safe way for the ambiguity to fall.
    """
    i = start + len(_DOCTYPE)
    quote = ''
    n = len(text)
    while i < n:
        c = text[i]
        if quote:
            if c == quote:
                quote = ''
        elif c in '"\'':
            quote = c
        elif c == '[':
            return True, -1
        elif c == '>':
            return False, i + 1
        i += 1
    return True, -1


def _find_doctype(text):
    """``(start, has_internal_subset, end)`` of the prolog's DOCTYPE, or ``None``.

    Scans only the prolog and stops at the root element's ``<``: that is the only place a DOCTYPE may
    legally appear, so a ``<!DOCTYPE`` in element content or CDATA is text.  Comments and processing
    instructions are stepped over, since either may precede the declaration and contain a ``<``.
    """
    i = 0
    n = len(text)
    while i < n:
        j = text.find('<', i)
        if j < 0:
            return None
        if text.startswith('<!--', j):
            k = text.find('-->', j + 4)
            if k < 0:
                return None  # unterminated comment; the backend will report it
            i = k + 3
        elif text.startswith('<?', j):
            k = text.find('?>', j + 2)
            if k < 0:
                return None
            i = k + 2
        elif text.startswith(_DOCTYPE, j):
            internal, end = _doctype_span(text, j)
            return j, internal, end
        else:
            return None  # the root element, or junk the backend will refuse
    return None


def _apply_entity_policy(data, allow_dtd, log):
    """`data` with any external-only DOCTYPE blanked out, or :class:`ForbiddenXml`.

    Blanked rather than cut so every byte offset a backend reports afterwards still points at the same
    place in the caller's file.
    """
    if isinstance(data, bytes):
        # latin-1 is byte-for-byte, so an index into this string is an index into `data`.  Every
        # encoding an XML document may declare is ASCII-compatible in the prolog except UTF-16, and
        # `parse_xml` has already decoded that case to `str`.
        text = data.decode('latin-1')
        blank = b' '
    else:
        text = data
        blank = ' '
    found = _find_doctype(text)
    if found is None:
        return data
    start, internal, end = found
    if internal:
        if not allow_dtd:
            raise ForbiddenXml(
                'the document declares an internal DTD subset, which is the only place an XML '
                'entity can be declared; refused rather than expanded. Pass allow_dtd=True to read '
                'a file you trust')
        return data
    if not allow_dtd:
        log.append(LogRecord('xml:doctype-dropped', (),
                             'record: DOCTYPE declaration dropped; an external DTD is never fetched, so any '
                             'entity it declares would be undefined anyway', LOST))
        return data[:start] + blank * (end - start) + data[end:]
    return data


# rule 2 -- the depth limit

def _check_depth(root, limit):
    """Refuse a tree nesting deeper than `limit`.  Iterative: recursion here would be the bug."""
    stack = [(root, 1)]
    while stack:
        node, depth = stack.pop()
        if depth > limit:
            raise ForbiddenXml(f'element <{local(node.tag)}> nests {depth} deep, past the limit of '
                               f'{limit}; refused. Raise max_depth to read a document you trust')
        for child in node:
            stack.append((child, depth + 1))


# the entry point

def _decode_utf16(data):
    """`data`, decoded if its first bytes are UTF-16 of either endianness, with a BOM or without.

    The one encoding that must be handled before the bytes reach a backend: every other encoding an XML
    document may declare is ASCII-compatible in its prolog, while a UTF-16 prolog is invisible to a byte
    scan, so rule 1 would not run.  A BOM-less UTF-16 document is not conforming but expat accepts one,
    so the mark is not the only spelling to test for.  The test is the byte-pattern table of **XML 1.0
    Appendix F**: a conforming document begins with ``<`` or a mark, so ``3C 00`` / ``00 3C`` in the
    first two bytes is UTF-16 and nothing else.  UCS-4 shares those bytes and is ruled out on the other
    two -- expat refuses it either way, so it is left to say so rather than mis-decoded here.
    """
    if not isinstance(data, bytes):
        return data
    head = data[:4]
    if head in (b'\x3c\x00\x00\x00', b'\x00\x00\x00\x3c',    # UCS-4, little and big endian
                b'\x00\x00\x3c\x00', b'\x00\x3c\x00\x00',    # UCS-4, the two unusual octet orders
                b'\xff\xfe\x00\x00', b'\x00\x00\xfe\xff'):   # UCS-4 with a mark
        return data
    if head[:2] in (b'\xff\xfe', b'\xfe\xff'):
        return data.decode('utf-16')
    if head[:2] == b'\x3c\x00':
        return data.decode('utf-16-le')
    if head[:2] == b'\x00\x3c':
        return data.decode('utf-16-be')
    return data


def _read_source(source):
    """`source` as ``str`` or ``bytes``, from a string, a path or an open file.

    A ``str`` or ``bytes`` with no ``<`` in it is treated as a path, not a document: callers pass both
    spellings, and this discrimination cannot misfire, since every XML document holds at least one
    ``<`` for its root element.
    """
    if isinstance(source, (str, bytes, bytearray)):
        data = bytes(source) if isinstance(source, bytearray) else source
        probe = b'<' if isinstance(data, bytes) else '<'
        if probe not in data:
            with open(source, 'rb') as f:
                return _decode_utf16(f.read())
    elif isinstance(source, PathLike) or hasattr(source, '__fspath__'):
        with open(source, 'rb') as f:
            data = f.read()
    elif hasattr(source, 'read'):
        data = source.read()
    else:
        raise TypeError(f'cannot read XML from {type(source).__name__}')
    return _decode_utf16(data)


def parse_xml(source, *, log=None, engine=None, max_depth=MAX_DEPTH, allow_dtd=False):
    """Parse `source` into an ``Element``, applying the entity and depth policy.  Returns the root.

    `source` is a ``str``, ``bytes``, a path, or an open file in either mode.  `engine` names the
    backend -- see :data:`ENGINES` -- and exists so a test can pin one; ``None`` picks the first
    available, the two behaving the same by construction.

    Raises :class:`~._errors.ForbiddenXml` for a document the policy declines to expand and
    :class:`~._errors.MalformedXml` for one that is not well-formed.  Never returns ``None``.
    """
    out = [] if log is None else log
    data = _apply_entity_policy(_read_source(source), allow_dtd, out)

    have = _backends()
    if engine is None:
        name = next((e for e in ENGINES if e in have), None)
    elif engine in have:
        name = engine
    else:
        raise ValueError(f'XML engine {engine!r} is not available here; have '
                         f'{available_engines()}')
    if name is None:  # pragma: no cover - 'stdlib' is unconditional, so this cannot fire
        raise RuntimeError('no XML backend available')

    try:
        root = have[name](data, allow_dtd)
    except ParseError as e:
        raise MalformedXml(f'not well-formed XML: {e}') from e
    except ForbiddenXml:
        raise
    except Exception as e:
        # `defusedxml` raises its own `EntitiesForbidden` / `DTDForbidden`, which are not `ParseError`
        # and must not surface as themselves, or the exception type would depend on which optional
        # package is installed.  Anything else a backend raises on bad bytes lands here for that reason.
        forbidden = type(e).__name__ in ('DTDForbidden', 'EntitiesForbidden',
                                         'ExternalReferenceForbidden')
        cls = ForbiddenXml if forbidden else MalformedXml
        raise cls(f'{"refused" if forbidden else "not well-formed"} by the {name} backend: '
                  f'{type(e).__name__}: {e}') from e
    if root is None:  # pragma: no cover - defensive; no backend documents this
        raise MalformedXml('the document has no root element')
    _check_depth(root, max_depth)
    return root


# namespace helpers -- every dialect needs these and none of them should own a copy

def local(tag):
    """The local name of an ``Element`` tag, with any ``{namespace}`` prefix removed.

    Every dialect matches on local names: chemical XML in the wild routinely declares the right
    vocabulary under the wrong namespace, so a reader keyed on the URI refuses files it understands.
    The namespace is read by :func:`namespace_of` and picks the *dialect*, not the element.
    """
    if isinstance(tag, str) and tag.startswith('{'):
        return tag[tag.index('}') + 1:]
    return tag


def namespace_of(tag):
    """The namespace URI of an ``Element`` tag, or ``''`` when it has none."""
    if isinstance(tag, str) and tag.startswith('{'):
        return tag[1:tag.index('}')]
    return ''


def text_of(node):
    """`node`'s text content, stripped, or ``''``.

    ``Element.text`` is ``None`` for ``<bondStereo/>`` and ``'\\n  '`` for a pretty-printed
    ``<bondStereo>\\n  W\\n</bondStereo>``, so every dialect reading element content needs this.
    """
    return (node.text or '').strip()
