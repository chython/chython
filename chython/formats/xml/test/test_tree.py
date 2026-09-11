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
"""The tokenizer and its safety policy, tested below every dialect and under both backends.

Two claims: an entity cannot expand because it cannot be declared (a prolog scan, not a parser
budget), and the guarantee does not depend on ``defusedxml``.  The rest is the ``str``/``bytes``/path/
file-handle question, answered by one rule -- no ``<``, no document.
"""

from io import BytesIO, StringIO
from pathlib import Path

from pytest import mark, raises

from .._errors import ForbiddenXml, MalformedXml
from .._tree import (ENGINES, MAX_DEPTH, available_engines, local, namespace_of, parse_xml, text_of)


#: A document with no DTD, no namespace and one element.  Small on purpose: tests below vary exactly
#: one thing about it, so a second interesting feature would make the assertion ambiguous.
PLAIN = '<molecule id="m1"/>'


# rule 1, the entities

def test_the_entity_bomb_does_not_expand(data, engine):
    """Nine nested entities, 700 bytes on disk, a gigabyte expanded.

    Refused for having an internal DTD subset at all -- the only place XML lets an entity be declared
    -- so nothing begins to expand and there is no memory ceiling or timeout to tune.
    """
    with raises(ForbiddenXml, match='internal DTD subset'):
        parse_xml(data('cml_entity_bomb.cml'), engine=engine)


def test_the_refusal_is_ours_and_not_the_backend_s(data):
    """Both backends refuse with the same message, so the refusal came from :mod:`.._tree`.

    The stdlib parser would expand this file on its own.  Asserted as one statement over
    ``available_engines()``, not parameterized, so that a disagreement is what fails.
    """
    seen = {}
    for name in available_engines():
        with raises(ForbiddenXml) as info:
            parse_xml(data('cml_entity_bomb.cml'), engine=name)
        seen[name] = str(info.value)
    assert 'stdlib' in seen, available_engines()
    assert len(set(seen.values())) == 1, seen  # the same message, so the same code path refused


def test_a_bare_entity_declaration_is_refused_even_with_nothing_to_expand(engine):
    """An unreferenced, harmless entity declaration is refused too: the policy is "no entity
    declarations", not "no dangerous ones", so no judgement call needs auditing."""
    with raises(ForbiddenXml, match='internal DTD subset'):
        parse_xml('<!DOCTYPE m [<!ENTITY x "y">]>\n<molecule id="m1"/>', engine=engine)


def test_an_unterminated_doctype_counts_as_having_a_subset(engine):
    """A ``<!DOCTYPE`` with no closing ``>`` could hide a subset, so the ambiguity falls to refusal
    rather than to handing the parser an unblanked declaration."""
    with raises(ForbiddenXml):
        parse_xml('<!DOCTYPE molecule SYSTEM "cml.dtd"', engine=engine)


def test_a_quoted_bracket_in_a_doctype_is_not_a_subset(engine):
    """``SYSTEM "a[b].dtd"``: the scan tracks quotes, so a ``[`` inside the public identifier is not
    an internal subset and the document is read."""
    log = []
    root = parse_xml('<!DOCTYPE molecule SYSTEM "a[b].dtd">\n<molecule id="m1"/>', log=log,
                     engine=engine)
    assert local(root.tag) == 'molecule'
    assert len(log) == 1 and 'DOCTYPE' in log[0], log


def test_an_external_doctype_is_read_and_logged(data, engine):
    """An external DOCTYPE -- the shape a real CML 1 file has -- is read and logged.

    Nothing here fetches a URL, so an external DTD is unreachable and any entity it declares would be
    undefined; refusing it would refuse documents that parse perfectly.
    """
    log = []
    root = parse_xml(data('cml_external_dtd.cml'), log=log, engine=engine)
    assert local(root.tag) == 'molecule'
    assert [str(x) for x in log] == ['record: DOCTYPE declaration dropped; an external DTD is never fetched, so any '
                                     'entity it declares would be undefined anyway']


def test_the_dropped_doctype_preserves_byte_offsets(engine):
    """The declaration is blanked in place, not cut out, so a parse error further down still names
    the line and column it has in the caller's file."""
    text = '<!DOCTYPE molecule SYSTEM "cml.dtd">\n<molecule id="m1">\n<atomArray>\n</molecule>'
    with raises(MalformedXml) as info:
        parse_xml(text, engine=engine)
    assert 'line 4' in str(info.value), str(info.value)


def test_allow_dtd_widens_rule_one_and_not_rule_two(engine):
    """``allow_dtd=True`` admits the internal subset and does not raise the depth limit.

    Two keywords on purpose: trusting a file's provenance says nothing about how deeply it nests.
    """
    text = '<!DOCTYPE m [<!ENTITY x "y">]>\n<molecule id="m1"/>'
    assert local(parse_xml(text, engine=engine, allow_dtd=True).tag) == 'molecule'
    deep = '<a>' * 12 + '</a>' * 12
    with raises(ForbiddenXml, match='nests'):
        parse_xml(deep, engine=engine, allow_dtd=True, max_depth=10)


def test_allow_dtd_means_the_same_thing_under_both_backends():
    """With ``allow_dtd=True`` both backends admit the subset and expand the entity.

    The backend guards must track the keyword, or a caller's contract depends on which optional
    package is installed.
    """
    text = '<!DOCTYPE m [<!ENTITY x "chython">]>\n<molecule id="&x;"/>'
    seen = {name: parse_xml(text, engine=name, allow_dtd=True).attrib for name in available_engines()}
    assert 'stdlib' in seen, available_engines()
    assert list(seen.values()) == [{'id': 'chython'}] * len(seen), seen


def test_a_doctype_inside_element_content_is_text(engine):
    """The scan stops at the root element's ``<``, so a ``<!DOCTYPE`` in a CDATA section is content:
    the prolog is the only place the declaration is legal."""
    root = parse_xml('<molecule id="m1"><name><![CDATA[<!DOCTYPE x [<!ENTITY e "v">]>]]></name>'
                     '</molecule>', engine=engine)
    assert '<!DOCTYPE' in text_of(root[0])


def test_a_comment_before_the_doctype_is_stepped_over(engine):
    """A comment may contain a ``<``, so a scan stopping at the first one would fail *open* and hand
    the declaration behind it straight to the parser."""
    with raises(ForbiddenXml, match='internal DTD subset'):
        parse_xml('<!-- a <b> in a comment -->\n<!DOCTYPE m [<!ENTITY x "y">]>\n<molecule/>')


def test_the_xml_declaration_is_stepped_over(engine):
    """``<?xml version="1.0"?>`` precedes the DOCTYPE in every real file -- same fail-open risk as a
    comment before it."""
    with raises(ForbiddenXml, match='internal DTD subset'):
        parse_xml('<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE m [<!ENTITY x "y">]><molecule/>',
                  engine=engine)


# rule 2, the depth

def test_nesting_past_the_limit_is_refused(engine):
    """Rule 1 makes the tree linear in input size, so only nesting can turn a small document into
    deep recursion downstream.  The limit is a keyword and the message says so."""
    with raises(ForbiddenXml, match='nests 11 deep, past the limit of 10'):
        parse_xml('<a>' * 11 + '</a>' * 11, engine=engine, max_depth=10)


def test_the_limit_is_inclusive_at_the_boundary(engine):
    """Exactly ``max_depth`` passes, one more fails: an off-by-one here is invisible at any other
    nesting count."""
    assert parse_xml('<a>' * 10 + '</a>' * 10, engine=engine, max_depth=10) is not None
    with raises(ForbiddenXml):
        parse_xml('<a>' * 11 + '</a>' * 11, engine=engine, max_depth=10)


def test_the_default_limit_has_three_orders_of_magnitude_of_headroom():
    """Chemical XML is shallow -- CML's deepest path is five levels, MRV's about eight -- so the
    default has headroom, and lowering it towards real files takes a deliberate edit."""
    assert MAX_DEPTH == 200


def test_the_depth_check_does_not_recurse(engine):
    """The walk is iterative: a recursive check crashes on a document nesting deeper than the
    interpreter's recursion limit, before it can refuse anything.  Nested past the usual 1000."""
    with raises(ForbiddenXml, match='nests'):
        parse_xml('<a>' * 1500 + '</a>' * 1500, engine=engine)


# what is a document, what a path

def test_a_string_with_no_angle_bracket_is_a_path(data):
    """The rule: callers pass both content and filenames as strings, and the only discrimination
    that cannot misfire is the ``<`` every XML document must hold for its root element."""
    root = parse_xml(str(data('cml_external_dtd.cml')), log=[])
    assert local(root.tag) == 'molecule'


def test_a_string_with_an_angle_bracket_is_a_document():
    """A document is never opened as a file, so passing content cannot raise ``FileNotFoundError``
    naming three kilobytes of XML."""
    assert local(parse_xml(PLAIN).tag) == 'molecule'


def test_a_missing_path_raises_the_file_error_and_not_a_parse_error():
    """A nonexistent filename fails as a filename, not as ``MalformedXml`` in a document never
    read."""
    with raises(FileNotFoundError):
        parse_xml('no_such_file.cml')


def test_bytes_are_accepted_both_ways(data):
    """``bytes`` is treated like ``str``: content if it holds a ``<``, a path if not.  A socket read
    and ``os.fsencode`` of a filename are both ``bytes``."""
    assert local(parse_xml(PLAIN.encode()).tag) == 'molecule'
    assert local(parse_xml(str(data('cml_external_dtd.cml')).encode(), log=[]).tag) == 'molecule'


def test_a_path_object_is_read_without_the_angle_bracket_rule(data):
    """A ``Path`` dispatches on its type and is never inspected for a ``<``."""
    assert isinstance(data('cml_external_dtd.cml'), Path)
    assert local(parse_xml(data('cml_external_dtd.cml'), log=[]).tag) == 'molecule'


def test_an_open_file_is_read_in_either_mode():
    """Text and binary handles both, so a caller who already opened the file need not know which
    mode this package prefers."""
    assert local(parse_xml(StringIO(PLAIN)).tag) == 'molecule'
    assert local(parse_xml(BytesIO(PLAIN.encode())).tag) == 'molecule'


def test_a_utf16_byte_order_mark_is_decoded_before_the_backend_sees_it():
    """UTF-16 is the one encoding decoded here rather than left to the declaration.

    Every other encoding is ASCII-compatible in its prolog; a UTF-16 prolog is invisible to a byte
    scan, so the DOCTYPE check would miss an internal subset entirely.
    """
    doc = ('<?xml version="1.0" encoding="UTF-16"?>' + PLAIN)
    assert local(parse_xml(doc.encode('utf-16')).tag) == 'molecule'
    with raises(ForbiddenXml, match='internal DTD subset'):
        parse_xml('<!DOCTYPE m [<!ENTITY x "y">]><molecule/>'.encode('utf-16'))


@mark.parametrize('encoding', ['utf-16-le', 'utf-16-be'])
def test_utf16_with_no_byte_order_mark_is_decoded_too(encoding, engine):
    """BOM-less UTF-16 is non-conforming and expat auto-detects it anyway, so rule 1 must too.

    Detected by the byte-pattern table of XML 1.0 Appendix F: a document begins with ``<`` or a mark,
    so ``3C 00`` or ``00 3C`` is UTF-16 of a known endianness.
    """
    doc = '<?xml version="1.0" encoding="UTF-16"?>' + PLAIN
    assert local(parse_xml(doc.encode(encoding), engine=engine).tag) == 'molecule'
    with raises(ForbiddenXml, match='internal DTD subset'):
        parse_xml('<!DOCTYPE m [<!ENTITY x "y">]><molecule/>'.encode(encoding), engine=engine)


def test_the_backends_agree_about_a_bom_less_utf16_entity_declaration():
    """Which optional package is installed must not decide what a reader accepts.

    Asserted as an agreement over ``available_engines()``, like the entity bomb above, so that a
    disagreement is what fails.
    """
    doc = ('<?xml version="1.0" encoding="UTF-16"?>'
           '<!DOCTYPE cml [<!ENTITY a "AAAAAAAAAA"><!ENTITY b "&a;&a;&a;">]>'
           '<cml><molecule title="&b;"/></cml>').encode('utf-16-le')
    seen = {}
    for name in available_engines():
        with raises(ForbiddenXml) as info:
            parse_xml(doc, engine=name)
        seen[name] = str(info.value)
    assert 'stdlib' in seen, available_engines()
    assert len(set(seen.values())) == 1, seen


@mark.parametrize('encoding', ['utf-32-le', 'utf-32-be'])
def test_ucs4_is_left_to_the_backend_to_refuse(encoding, engine):
    """UCS-4 opens with the same two bytes as UTF-16 and is ruled out on the other two: decoding it
    as UTF-16 would blame the file for an unsupported encoding."""
    doc = ('<?xml version="1.0" encoding="UTF-32"?>' + PLAIN).encode(encoding)
    with raises(MalformedXml, match='not well-formed'):
        parse_xml(doc, engine=engine)


def test_an_unreadable_type_is_a_type_error():
    """An integer is a caller mistake, so it reads as ``TypeError`` rather than as a bad document."""
    with raises(TypeError, match='cannot read XML from int'):
        parse_xml(42)


# malformation, and backends

def test_a_malformed_document_is_malformed_and_not_forbidden(engine):
    """Two distinct exceptions on purpose: a corpus sweep counts "this file is broken" apart from
    "we declined to expand this file"."""
    with raises(MalformedXml, match='not well-formed'):
        parse_xml('<molecule><atomArray></molecule>', engine=engine)


def test_the_empty_document_is_malformed(engine):
    """No root element.  Both backends must arrive as ``MalformedXml``, so a caller's ``except``
    clause does not depend on which package is installed."""
    with raises(MalformedXml):
        parse_xml('<', engine=engine)


def test_a_backend_specific_exception_never_escapes_as_itself(data):
    """``defusedxml``'s ``DTDForbidden``/``EntitiesForbidden`` are not ``ParseError``, and are
    translated: the exception type a caller catches must not depend on an optional dependency."""
    for name in available_engines():
        with raises(ForbiddenXml) as info:
            parse_xml(data('cml_entity_bomb.cml'), engine=name)
        assert type(info.value) is ForbiddenXml, type(info.value)


def test_an_unavailable_engine_is_named_in_the_error():
    """An uninstalled backend is a caller error, not a silent fall back -- which would let a test
    pinning one backend pass while measuring another."""
    with raises(ValueError, match='is not available here'):
        parse_xml(PLAIN, engine='lxml')


def test_the_stdlib_backend_is_always_available():
    """``'stdlib'`` is always in ``available_engines()``: there is no hard dependency on
    ``defusedxml`` and the package must be testable without it."""
    assert 'stdlib' in available_engines()
    assert set(available_engines()) <= set(ENGINES)
    assert available_engines() == tuple(e for e in ENGINES if e in available_engines())  # preference


def test_the_backends_are_not_imported_at_import_time():
    """The backend cache is built on first use: importing the optional dependency at import time
    would make "no hard dependency" unmeasurable."""
    from .. import _tree
    assert isinstance(_tree._BACKEND_CACHE, dict)
    assert _tree._backends() is _tree._BACKEND_CACHE  # built once, returned by identity


# namespace helpers

@mark.parametrize('tag,name,ns', [('{http://www.xml-cml.org/schema}molecule', 'molecule',
                                   'http://www.xml-cml.org/schema'),
                                  ('molecule', 'molecule', ''),
                                  ('{}molecule', 'molecule', '')])
def test_local_and_namespace_split_a_tag(tag, name, ns):
    """``ElementTree`` produces the ``{}`` form for ``xmlns=""``, and it must read as "no namespace"
    rather than as an empty URI: the difference decides whether a dialect claims the document."""
    assert local(tag) == name
    assert namespace_of(tag) == ns


def test_text_of_normalizes_the_three_shapes_of_element_content():
    """``Element.text`` is ``None`` for ``<bondStereo/>``, ``'\\n  W\\n'`` pretty-printed and ``'W'``
    compact; the ``None`` case is the one a dialect's own copy forgets."""
    root = parse_xml('<b><x/><y>\n  W\n</y><z>W</z></b>')
    assert [text_of(c) for c in root] == ['', 'W', 'W']
