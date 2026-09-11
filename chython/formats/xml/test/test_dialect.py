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
"""The engine, tested through a dialect that is not CML.

``TOY`` below is a whole dialect in forty lines -- no subclassing, atoms spelled ``<a e="C">`` inside
``<atoms>`` -- so a CML assumption growing into the engine fails over ``TOY`` and not over CML.  Also
measured: the acceptance rule, apply a construct or name it in the log, which the engine enforces.
"""

from xml.etree.ElementTree import tostring

from pytest import raises

from chython.core import write_smiles

from .._dialect import (Dialect, Field, NotModelled, Record, Tags, apply_fields, dialect, dialects,
                        molecule_nodes, parse_xml_document, read_document, read_molecule, read_xml,
                        register, sniff, write_molecule)
from .._errors import ForbiddenXml, MalformedXml, UnsupportedXml
from .._tree import local, parse_xml, text_of
from ...ctfile._hydrogens import StatedChannels


TOY_NS = 'urn:chython-test:toy'


def _element(text):
    """``elementType`` in a dialect that spells it ``e``, refusing a pseudo-atom the way CML does."""
    value = text.strip()
    if not value:
        raise ValueError('empty')
    if value == '*':
        raise UnsupportedXml('a pseudo-atom has nowhere to go in a molecule')
    return value


def _order(text):
    """A bond order, with one legal value the container cannot hold -- ``NotModelled``'s whole point."""
    value = text.strip()
    if value == 'half':
        raise NotModelled('a half bond is not modelled')
    return int(value)


def _toy_bond_child(node, record, position, log):
    """``<flip/>`` on a bond, standing in for CML's ``<bondStereo>``: a child element a table cannot
    express, its meaning being in its content rather than in an attribute."""
    if local(node.tag) != 'flip':
        return False
    record.bond_extras[position]['flip'] = text_of(node) or 'yes'
    return True


def _toy_finish(record, log):
    """Turn what the walk parked into what the build needs, once per molecule.

    In CML it is the atom parities, which name atoms the walker has not reached.  Here it does nothing
    else, so the contract -- after the whole molecule, before the build -- is what gets measured.
    """
    record.extras['finished'] = len(record.ctab.atoms)


TOY = Dialect(
    name='toy',
    # A made-up spelling on purpose: the advice half of an unknown-count log line comes from here and
    # not from the shared derivation, which is MDL's.  `valence=None` is the case CML is a real one of.
    channels=StatedChannels(count='an `hs` attribute', valence=None),
    namespaces=frozenset((TOY_NS,)),
    ns=TOY_NS,
    tags=Tags(molecule='mol', atom_array='atoms', atom='a', bond_array='bonds', bond='b',
              roots=frozenset(('toy', 'mol'))),
    atom_id='n',
    bond_refs=('ends', 'atomRefs2'),
    molecule_fields=(Field('label', 'title', str, lambda v: v or None),),
    atom_fields=(Field('e', 'element', _element, str),
                 Field('q', 'charge', int, lambda v: str(v) if v else None),
                 Field('depth', 'z', float, None)),
    bond_fields=(Field('o', 'order', _order, str),),
    atom_ignored=frozenset(('pretty',)),
    molecule_children_ignored=frozenset(('meta',)),
    bond_child=_toy_bond_child,
    finish=_toy_finish)

#: One molecule in `TOY`, exercising every slot the dialect declares.  Deliberately not valid CML.
TOY_DOC = f'''<toy xmlns="{TOY_NS}">
  <mol label="toy ethanolate">
    <meta>ignored by declaration</meta>
    <atoms>
      <a n="x1" e="C" pretty="yes"/>
      <a n="x2" e="O" q="-1"/>
    </atoms>
    <bonds>
      <b n="k1" ends="x1 x2" o="1"><flip>W</flip></b>
    </bonds>
  </mol>
</toy>'''


def _read_toy(text, log=None):
    out = [] if log is None else log
    return read_document(parse_xml(text), TOY, out)


# the design test: a dialect is a table

def test_a_dialect_that_is_not_cml_reads_with_no_engine_change():
    """``TOY`` shares not one element or attribute name with CML -- ``<mol>``, ``<atoms>``, ``<a>``,
    ``n``, ``ends`` -- and is a ``Dialect`` value plus two functions nothing in the engine knows of."""
    record, = _read_toy(TOY_DOC)
    assert record.ctab.title == 'toy ethanolate'
    assert [a.element for a in record.ctab.atoms] == ['C', 'O']
    assert [a.charge for a in record.ctab.atoms] == [0, -1]
    assert record.ids == ['x1', 'x2']
    assert [(b.a, b.b, b.order) for b in record.ctab.bonds] == [(0, 1, 1)]


def test_the_toy_dialect_declares_no_handler_it_does_not_need():
    """Every handler slot defaults to ``None``, so a dialect declares only what it has.

    ``TOY`` has a bond child and a ``finish`` and nothing else, and it reads.  Twelve required
    callables would make the next dialect a copy of CML's.
    """
    assert TOY.atom_child is None and TOY.molecule_child is None
    assert TOY.atom_array_hook is None and TOY.bond_array_hook is None
    assert TOY.emit_atom is None and TOY.emit_bond is None and TOY.emit_molecule is None
    assert _read_toy(TOY_DOC)


def test_finish_runs_once_after_the_walk_and_before_the_build():
    """``finish`` resolves a field that depends on the finished graph, so it must see every atom.
    ``TOY``'s counts them: a hook running per-element, or before the arrays, records fewer."""
    record, = _read_toy(TOY_DOC)
    assert record.extras['finished'] == 2


def test_a_bond_child_handler_claims_its_element_and_the_engine_stays_silent():
    """A handler returning ``True`` means claimed, and a claimed element must not also be logged.  Read
    *and* reported unsupported is invisible without the log assertion."""
    log = []
    record, = _read_toy(TOY_DOC, log)
    assert record.bond_extras[0] == {'flip': 'W'}
    assert log == [], log


# the acceptance rule

def test_an_attribute_with_no_row_is_logged_unsupported():
    """``spin="3"`` is a real thing a file may say and ``TOY`` has no row for it, so it is named --
    prefixed ``unsupported: ``, because the file is fine and we are the limitation."""
    log = []
    _read_toy(TOY_DOC.replace('e="O" q="-1"', 'e="O" q="-1" spin="3"'), log)
    assert [str(x) for x in log] == ['unsupported: atom x2: attribute spin=\'3\' is not modelled'], log


def test_an_ignored_attribute_is_silent_and_that_is_a_claim():
    """``pretty`` is in ``atom_ignored``, so it produces no line -- a reader honouring it would build
    the same molecule.  Widening that set until the log goes quiet is what makes the prefix meaningless."""
    log = []
    _read_toy(TOY_DOC, log)
    assert log == [], log
    assert 'pretty' in TOY.atom_ignored


def test_a_structural_attribute_is_silent_without_being_declared_ignored():
    """The identity attribute and the bond endpoints are consumed by the engine, so they need no row
    and no ignore entry -- otherwise every atom and bond in the file would cost two log lines."""
    log = []
    _read_toy(TOY_DOC, log)
    assert not any('ends' in x or ' n=' in x for x in log), log


def test_a_malformed_value_costs_the_attribute_and_not_the_atom():
    """``ValueError`` from a decode is a plain line, no prefix: the file was broken and we read it
    anyway.  The atom survives with its default charge -- a nonsense charge is not a lost atom."""
    log = []
    record, = _read_toy(TOY_DOC.replace('q="-1"', 'q="minus one"'), log)
    assert len(record.ctab.atoms) == 2
    assert record.ctab.atoms[1].charge == 0
    assert len(log) == 1 and str(log[0]).startswith('atom x2: attribute q='), log
    assert not str(log[0]).startswith('unsupported'), log


def test_not_modelled_is_the_unsupported_half_of_the_same_branch():
    """Why ``NotModelled`` is a class: ``o="half"`` is legal in this format with no representation in a
    molecule, ``o="wrong"`` is garbage.  Both come out of one ``decode`` and cost the attribute, and only
    the prefix tells a corpus sweep "our reader is short" from "your files are broken"."""
    unsupported, malformed = [], []
    _read_toy(TOY_DOC.replace('o="1"', 'o="half"'), unsupported)
    _read_toy(TOY_DOC.replace('o="1"', 'o="wrong"'), malformed)
    assert len(unsupported) == 1 and str(unsupported[0]).startswith('unsupported: bond k1:'), unsupported
    assert len(malformed) == 1 and str(malformed[0]).startswith('bond k1:'), malformed
    assert not str(malformed[0]).startswith('unsupported'), malformed


def test_an_xml_error_from_a_decode_refuses_the_whole_record():
    """The one answer that must not be a log line: an element naming a pseudo-atom leaves the atom
    nowhere to go, and defaulting it to carbon is the silent wrong answer, so it is re-raised through
    ``apply_fields`` rather than caught with the ``ValueError``s."""
    with raises(UnsupportedXml, match='pseudo-atom'):
        _read_toy(TOY_DOC.replace('e="C"', 'e="*"'))


def test_a_child_element_no_handler_claims_is_logged_unsupported():
    """The element half of the rule, in all three positions, and the element name is in the line -- a
    log saying only "unsupported child" is one nobody can act on."""
    log = []
    _read_toy(TOY_DOC.replace('<atoms>', '<electron/><atoms>'), log)
    assert [str(x) for x in log] == ['unsupported: record: <electron> in <mol> is not modelled'], log

    log = []
    _read_toy(TOY_DOC.replace('e="C" pretty="yes"/>', 'e="C"><shell/></a>'), log)
    assert [str(x) for x in log] == ['unsupported: atom x1: <shell> is not modelled'], log

    log = []
    _read_toy(TOY_DOC.replace('<atoms>', '<atoms><ghost/>'), log)
    assert [str(x) for x in log] == ['unsupported: atom: <ghost> in <atoms> is not modelled'], log


def test_an_ignored_child_element_is_silent():
    """``<meta>`` is in ``molecule_children_ignored`` and ``TOY_DOC`` contains one, so every other
    test here would carry a spurious line if the set were not honoured."""
    assert 'ignored by declaration' in TOY_DOC
    log = []
    _read_toy(TOY_DOC, log)
    assert log == [], log


def test_every_log_line_starts_with_a_token_from_the_closed_set():
    """After any ``unsupported: `` comes one of a closed set of location words.  ``molecule``,
    ``atomArray`` and ``bondArray`` are not in it, which is why the engine says ``record`` and puts the
    element name inside the message."""
    log = []
    _read_toy(TOY_DOC.replace('<atoms>', '<electron/><atoms>').replace('o="1"', 'o="half"'), log)
    assert log
    for line in log:
        s = str(line)
        body = s[len('unsupported: '):] if s.startswith('unsupported: ') else s
        assert body.split(':')[0].split()[0] in ('atom', 'bond', 'record', 'sgroup', 'coordinates',
                                                 'stereo'), line


# the engine's two attributes

def test_the_engine_consumes_exactly_two_attributes_on_its_own():
    """Identity and endpoints, since it has to resolve one against the other; everything else is a
    row.  A third attribute in the engine is a third thing every future dialect must spell its way."""
    assert TOY.atom_id == 'n'
    assert TOY.bond_refs[0] == 'ends'
    record, = _read_toy(TOY_DOC)
    assert record.index_of == {'x1': 0, 'x2': 1}


def test_bond_refs_is_a_preference_list():
    """Several spellings of the endpoints, most preferred first -- CML 2 writes ``atomRefs2`` and some
    writers emit ``atomRefs``.  ``TOY`` declares two, and the second works."""
    record, = _read_toy(TOY_DOC.replace('ends="x1 x2"', 'atomRefs2="x1 x2"'))
    assert [(b.a, b.b) for b in record.ctab.bonds] == [(0, 1)]


def test_a_bond_naming_an_unknown_atom_is_dropped_and_named():
    """A bond to an absent atom cannot be stored, and losing one silently is how a reader produces a
    plausible wrong molecule.  The rest of the record survives."""
    log = []
    record, = _read_toy(TOY_DOC.replace('ends="x1 x2"', 'ends="x1 x99"'), log)
    assert len(record.ctab.atoms) == 2 and not record.ctab.bonds
    assert len(log) == 1 and 'unknown atom' in log[0], log


def test_a_bond_naming_one_atom_or_three_is_dropped_and_named():
    """``atomRefs2`` with the wrong cardinality: there is nothing to guess about which end was meant,
    so the line says what it saw."""
    for refs in ('x1', 'x1 x2 x1'):
        log = []
        record, = _read_toy(TOY_DOC.replace('ends="x1 x2"', f'ends="{refs}"'), log)
        assert not record.ctab.bonds
        assert len(log) == 1 and 'does not name two atoms' in log[0], log


def test_a_bond_with_no_endpoints_at_all_is_dropped_and_named():
    log = []
    record, = _read_toy(TOY_DOC.replace('ends="x1 x2" ', ''), log)
    assert not record.ctab.bonds
    assert len(log) == 1 and 'no ends, dropped' in log[0], log


def test_a_duplicate_atom_id_lets_the_first_declaration_win():
    """Letting the later atom claim the name would silently *move* every bond written before it.

    The duplicated atom is kept -- the file states it -- and is simply unreachable by name.
    """
    doc = f'''<mol xmlns="{TOY_NS}"><atoms>
      <a n="x1" e="C"/><a n="x2" e="O"/><a n="x1" e="N"/>
    </atoms><bonds><b ends="x1 x2" o="1"/></bonds></mol>'''
    record, = _read_toy(doc, [])
    assert [a.element for a in record.ctab.atoms] == ['C', 'O', 'N']
    assert record.index_of == {'x1': 0, 'x2': 1}
    assert [(b.a, b.b) for b in record.ctab.bonds] == [(0, 1)]


def test_a_row_whose_slot_is_not_on_the_intermediate_spills_into_extras():
    """How a dialect states a field the neutral ``Ctab`` has no business growing.

    ``depth`` does have a ``CtabAtom`` slot, so the spill is measured on the *molecule* row instead.
    The intermediate is ``__slots__``-ed, so an unknown slot lands in the parallel dict.
    """
    log = []
    record, = _read_toy(TOY_DOC.replace('label=', 'nick='), log)
    # `nick` has no row at all, so it is reported rather than spilled: the spill is for a row that
    # exists and names a slot the target lacks.
    assert any('nick' in x for x in log), log
    spilling = TOY._replace(molecule_fields=(Field('label', 'nickname', str),))
    record = read_molecule(molecule_nodes(parse_xml(TOY_DOC), TOY)[0], spilling, [])
    assert record.extras == {'nickname': 'toy ethanolate', 'finished': 2}


def test_a_spill_with_nowhere_to_go_raises():
    """No ``spill`` plus a row naming an absent slot is a dialect bug and must surface as one: swallowed,
    a mistyped slot name becomes a field that silently never arrives."""
    class Target:
        __slots__ = ('kept',)
    row = Field('v', 'missing', str)
    with raises(AttributeError):
        apply_fields(parse_xml('<x v="1"/>'), {'v': row}, Target(), [], 'record')


def test_a_write_only_row_is_not_read_and_a_read_only_row_is_not_written():
    """Both cases are real -- an extension honoured on the way in need not be written, and a value
    emitted for a consumer need not be read back.  A row with no ``decode`` reports as a missing row."""
    log = []
    write_only = TOY._replace(atom_fields=(Field('e', 'element', _element, str),
                                           Field('q', 'charge', None, str)))
    record = read_molecule(molecule_nodes(parse_xml(TOY_DOC), write_only)[0], write_only, log)
    assert record.ctab.atoms[1].charge == 0
    assert any('attribute q=' in x and str(x).startswith('unsupported') for x in log), log
    # `depth` is read-only in TOY: it decodes and has no encode, so it never appears in output.
    assert TOY.atom_fields[2].encode is None
    written = tostring(write_molecule(_read_toy(TOY_DOC)[0], TOY, []), encoding='unicode')
    assert 'depth' not in written


# molecules in a tree

def test_only_the_outermost_molecule_is_read():
    """A nested molecule is CML's assembly, and the atoms of the parts are not also atoms of the whole
    -- reading both doubles every atom.  The nesting is still reported, by the walker that meets it."""
    log = []
    nested = TOY_DOC.replace('</atoms>', '</atoms><mol label="part"><atoms>'
                                         '<a n="y1" e="N"/></atoms></mol>')
    records = _read_toy(nested, log)
    assert len(records) == 1
    assert [a.element for a in records[0].ctab.atoms] == ['C', 'O']
    assert any('<mol> in <mol>' in x for x in log), log


def test_molecules_come_back_in_document_order():
    """Order is load-bearing for a caller matching a multi-molecule file against a list of names, so it
    is document order, breadth-first from the root."""
    two = TOY_DOC.replace('</mol>', '</mol><mol label="second"><atoms><a n="z1" e="N"/></atoms></mol>')
    assert [r.ctab.title for r in _read_toy(two)] == ['toy ethanolate', 'second']


def test_a_document_whose_root_is_the_molecule_needs_no_wrapper():
    """A bare ``<mol>`` at the root is a single-record file, so ``molecule_nodes`` matches the root
    itself and not only its children."""
    single = f'<mol xmlns="{TOY_NS}" label="alone"><atoms><a n="q1" e="C"/></atoms></mol>'
    record, = _read_toy(single)
    assert record.ctab.title == 'alone' and len(record.ctab.atoms) == 1


# the registry

def test_the_registry_is_lazy_and_holds_cml():
    """Tables load on first use, never at import.  ``dialects()`` is sorted, so listing is stable."""
    assert 'cml' in dialects()
    assert dialects() == tuple(sorted(dialects()))
    assert dialect('cml').name == 'cml'


def test_an_unknown_dialect_name_lists_what_there_is():
    """The message names the alternatives, since the likely cause is a typo or a dialect that has not
    landed.  Both shipped dialects must appear: a partial list sends the reader to the wrong place."""
    with raises(ValueError, match="no XML dialect named 'cdxml'; have"):
        dialect('cdxml')
    try:
        dialect('cdxml')
    except ValueError as e:
        assert "'cml'" in str(e) and "'mrv'" in str(e), str(e)


def test_registering_a_dialect_outside_this_package_needs_no_fork():
    """``register`` is public, so a house format or vendor variant of CML is a table rather than a patch
    to this module.  Unregistered afterwards, so suite order cannot matter."""
    from .. import _dialect
    register(TOY)
    try:
        assert 'toy' in dialects()
        assert dialect('toy') is TOY
    finally:
        _dialect._DIALECT_CACHE.pop('toy', None)
    assert 'toy' not in dialects()


def test_sniff_picks_a_dialect_by_namespace():
    """The namespace identifies a vocabulary, so it is matched first.  ``TOY``'s document declares its
    own URI and no CML element, so nothing but registration decides this."""
    from .. import _dialect
    register(TOY)
    try:
        log = []
        assert sniff(parse_xml(TOY_DOC), log) is TOY
        assert log == [], log
    finally:
        _dialect._DIALECT_CACHE.pop('toy', None)


def test_sniff_falls_back_on_the_root_name_and_says_so():
    """A document naming a vocabulary nobody claims is read by root name, with a line saying so.

    Synthetic rather than a fixture: what is under test is the *absence* of a dialect for a declared
    namespace, and a real file stops testing that the day its vocabulary lands.  ``urn:`` never can.
    """
    log = []
    dial = sniff(parse_xml('<cml xmlns="urn:example:unclaimed"><molecule/></cml>'), log)
    assert dial.name == 'cml'
    assert len(log) == 1 and str(log[0]).startswith('unsupported: record: no dialect claims '), log
    assert 'read as cml' in log[0]


def test_sniff_gives_a_marvin_document_to_the_mrv_dialect_without_a_line(data):
    """``test/implicit.mrv`` declares ``http://www.chemaxon.com``, which the MRV dialect claims, so the
    file never reaches the root-name fallback and earns no line.  The absence is checked against a log
    the sibling above proves capable of holding one."""
    log = []
    dial = sniff(parse_xml(data('implicit.mrv')), log)
    assert dial.name == 'mrv'
    assert log == [], log


def test_the_fallback_names_the_namespace_it_did_not_recognise():
    """"Unrecognised namespace" cannot be acted on: the URI is what tells a reader of the log which
    dialect to ask for."""
    log = []
    sniff(parse_xml('<cml xmlns="urn:not-a-format"><molecule/></cml>'), log)
    assert 'urn:not-a-format' in log[0], log


def test_a_document_with_no_namespace_at_all_is_read_without_a_line():
    """Plenty of CML in the wild declares no namespace, and reading one loses nothing.

    The fallback's line reports an *unclaimed* namespace -- a URI the document named and no dialect
    answered to.  A document naming no URI is not that: its root is a name this dialect claims and
    anything unmodelled inside it gets its own line.  The shape is a hand-written or pre-schema file;
    :func:`~.._cml.write_cml` declares the namespace, so our own output takes the other branch.
    """
    log = []
    assert sniff(parse_xml('<cml><molecule/></cml>'), log).name == 'cml'
    assert log == [], log


def test_a_namespace_that_was_declared_and_unclaimed_still_gets_its_line():
    """The pair of the test above: same root name, same molecules, one URI of difference -- and whether
    the document named a vocabulary is exactly what the line reports."""
    log = []
    assert sniff(parse_xml('<cml xmlns="urn:not-a-format"><molecule/></cml>'), log).name == 'cml'
    assert len(log) == 1 and str(log[0]).startswith('unsupported: record: no dialect claims '), log
    assert 'urn:not-a-format' in log[0] and 'read as cml' in log[0], log


def test_a_root_no_dialect_reads_is_refused():
    """The one refusal here: a ``<svg>`` document is a different kind of document, not a chemical file
    in a vocabulary we lack, and reading it as CML gives zero molecules and no explanation."""
    with raises(MalformedXml, match='<svg> is not a document root'):
        sniff(parse_xml('<svg><molecule/></svg>'), [])


# the dialect-agnostic entry point

def test_the_sniffing_reader_routes_to_a_registered_dialect():
    """``read_xml`` gives the registry a production caller -- every named reader states its dialect, so
    without it ``sniff`` and ``register`` are reachable only from a test.  ``TOY`` is the subject because
    a document only this test's registration can read proves the routing is by the table."""
    from .. import _dialect
    register(TOY)
    try:
        log = []
        records = parse_xml_document(TOY_DOC, log=log)
        assert [a.element for a in records[0].ctab.atoms] == ['C', 'O']
        assert log == [], log
    finally:
        _dialect._DIALECT_CACHE.pop('toy', None)


def test_the_sniffing_reader_builds_and_the_parsing_one_does_not():
    """The ``parse_``/``read_`` split, as every named dialect offers it.

    ``parse_xml_document`` stops at the record, so a caller sees the file's own atom ids and title before
    committing to a build; a reader that only ever built would make those identifiers unreachable.
    """
    doc = '<cml><molecule title="t"><atomArray><atom id="a1" elementType="C"/></atomArray></molecule></cml>'
    records = parse_xml_document(doc, log=[])
    assert records[0].ctab.title == 't' and records[0].ids == ['a1']
    molecules = read_xml(doc, log=[])
    assert len(molecules) == 1 and write_smiles(molecules[0]) == 'C'


def test_the_sniffing_reader_applies_the_entity_policy():
    """The hardening is not bypassable by choosing this entry point instead of a named one.

    It goes through ``parse_xml`` like everything else here; asserted rather than assumed, since a new
    entry point forgetting the policy is the failure one hardened tokenizer exists to make impossible.
    """
    bomb = ('<?xml version="1.0"?><!DOCTYPE cml [<!ENTITY a "aaaaaaaaaa">'
            '<!ENTITY b "&a;&a;&a;&a;&a;&a;&a;&a;&a;&a;">]><cml><molecule title="&b;"/></cml>')
    with raises(ForbiddenXml):
        read_xml(bomb, log=[])


def test_the_sniffing_reader_passes_its_keywords_to_the_tokenizer():
    """``engine``, ``max_depth`` and ``allow_dtd`` behave here as on the named readers.

    Checked with ``max_depth``, the one whose effect is visible without a DTD: a limit of two cannot
    reach an ``<atom>`` three levels down, and the refusal comes from the tokenizer.
    """
    doc = '<cml><molecule><atomArray><atom id="a1" elementType="C"/></atomArray></molecule></cml>'
    with raises(ForbiddenXml, match='past the limit of 2'):
        read_xml(doc, log=[], max_depth=2)
    assert read_xml(doc, log=[], max_depth=10)


# writing

def test_the_writer_is_driven_by_the_same_table_as_the_reader():
    """The classic round-trip defect is a value read one way and written another by a second body of
    code that agrees everywhere except one row.  A codec per row leaves nowhere for the two to diverge,
    measured on a dialect with no writer of its own."""
    record, = _read_toy(TOY_DOC)
    out = tostring(write_molecule(record, TOY._replace(ns=''), []), encoding='unicode')
    assert out.startswith('<mol label="toy ethanolate">')
    assert '<a n="x1" e="C" />' in out
    assert '<a n="x2" e="O" q="-1" />' in out
    assert '<b n="b1" ends="x1 x2" o="1" />' in out


def test_attributes_are_written_in_declaration_order():
    """Two runs over one molecule must be byte-identical.  Attribute order carries no meaning in XML,
    which is why it has to come from the table and not from whatever dict the table was built from."""
    record, = _read_toy(TOY_DOC)
    plain = TOY._replace(ns='')
    first = tostring(write_molecule(record, plain, []), encoding='unicode')
    assert first == tostring(write_molecule(record, plain, []), encoding='unicode')
    assert first.index(' e="O"') < first.index(' q="-1"')  # e is declared before q


def test_an_encode_returning_none_omits_the_attribute():
    """How a default stays unwritten: the charge encoder returns ``None`` for zero, so a neutral atom
    carries no ``q`` -- a file full of ``q="0"`` states something the writer did not mean."""
    record, = _read_toy(TOY_DOC)
    out = tostring(write_molecule(record, TOY._replace(ns=''), []), encoding='unicode')
    assert out.count('q=') == 1


def test_no_empty_bond_array_is_written_for_a_single_atom():
    """Real files omit it, and an empty array is a statement a reader then has to decide about."""
    single = f'<mol xmlns="{TOY_NS}"><atoms><a n="q1" e="C"/></atoms></mol>'
    record, = _read_toy(single)
    out = tostring(write_molecule(record, TOY._replace(ns=''), []), encoding='unicode')
    assert 'bonds' not in out


def test_the_writer_qualifies_with_the_dialect_s_own_namespace():
    """``namespaces`` is a set because real files declare several historical URIs for one vocabulary;
    ``ns`` is one value because we write exactly one."""
    record, = _read_toy(TOY_DOC)
    out = tostring(write_molecule(record, TOY, []), encoding='unicode')
    assert f'{{{TOY_NS}}}mol' in out or TOY_NS in out


def test_a_record_written_into_a_parent_becomes_a_subelement():
    """A document is a root with molecules under it, so the writer takes the parent rather than making
    the caller re-parent an element built standalone."""
    from xml.etree.ElementTree import Element
    root = Element('toy')
    record, = _read_toy(TOY_DOC)
    node = write_molecule(record, TOY._replace(ns=''), [], root)
    assert list(root) == [node]


# the Record type

def test_a_record_numbers_its_atoms_when_the_file_does_not():
    """A file may omit the identity attribute -- CML 1 files do -- so the record invents ``a1``,
    ``a2``, following the file's own convention so what we write back looks like what we read."""
    record = Record()
    assert record.add_atom(object()) == 0
    assert record.add_atom(object(), 'named') == 1
    assert record.ids == ['a1', 'named']
    assert len(record) == 2
    assert 'Record(2 atoms, 0 bonds' in repr(record)


# the dialect's stated channels

#: A bonded atom of an element the valence collection has no row for.  Xenon because the point is the
#: advice, not the chemistry: nothing derives a count, so the message naming both channels is reached.
_XE_CML = ('<molecule><atomArray>'
           '<atom id="a1" elementType="Xe" x2="0" y2="0"/>'
           '<atom id="a2" elementType="C" x2="1" y2="0"/>'
           '</atomArray><bondArray><bond atomRefs2="a1 a2" order="1"/></bondArray></molecule>')

_XE_MOLFILE = ['xenon', '', '',
               '  2  1  0  0  0  0  0  0  0  0999 V2000',
               '    0.0000    0.0000    0.0000 Xe  0  0  0  0  0  0  0  0  0  0  0  0',
               '    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
               '  1  2  1  0  0  0  0',
               'M  END']


def _unknown_count_line(log):
    """The one line that ends in advice, or a failure naming what the log said instead."""
    lines = [str(x) for x in log if 'count not known' in x]
    assert len(lines) == 1, log
    return lines[0]


def test_the_advice_names_the_channel_the_reading_dialect_actually_HAS():
    """A dialect's own spelling reaches the log, and MDL's does not reach an XML reader.

    The derivation is shared across both CTAB versions, MRV and CML, and every "count not known" line
    ends in advice -- so shared advice is wrong for two of the four callers: an `MRV_IMPLICIT_H` data
    S-group is how a *molfile* carries a stated count, and a Marvin document cannot hold one.
    """
    from ...ctfile import parse_v2000
    from .._cml import read_cml
    from .._mrv import read_mrv

    log = []
    read_mrv(_XE_CML, log=log)
    line = _unknown_count_line(log)
    assert '`hydrogenCount` attribute' in line and '`mrvValence` attribute' in line, line
    assert 'MRV_IMPLICIT_H' not in line and 'VAL=' not in line, line

    log = []
    read_cml(_XE_CML, log=log)
    line = _unknown_count_line(log)
    assert '`hydrogenCount` attribute' in line, line
    assert 'MRV_IMPLICIT_H' not in line and 'VAL=' not in line, line

    # The control that makes the two above a routing test and not a deletion: the molfile reader still
    # gets MDL's two spellings.
    log = []
    ctab = parse_v2000(_XE_MOLFILE, log)
    _mol, _store, build_log = ctab.build()
    line = _unknown_count_line(build_log)
    assert 'MRV_IMPLICIT_H data S-group' in line and 'VAL=' in line, line


def test_a_dialect_with_no_valence_channel_drops_the_clause_rather_than_respelling_it():
    """CML has `hydrogenCount` and no total-valence field at all, and the sentence shortens.

    MRV has both channels, so its advice offers a choice; CML's second option does not exist, and
    naming a CML attribute that would hold a valence would be inventing a field.
    `StatedChannels.valence is None` is the third answer, and this is it in the log.
    """
    from .._cml import read_cml
    from .._mrv import read_mrv

    log = []
    read_mrv(_XE_CML, log=log)
    mrv_line = _unknown_count_line(log)
    log = []
    read_cml(_XE_CML, log=log)
    cml_line = _unknown_count_line(log)

    assert 'neither a hydrogen count nor a usable valence' in mrv_line, mrv_line
    assert ' or ' in mrv_line.rsplit('State it with', 1)[1], mrv_line
    # The CML half: the file states no *count* -- not "no count and no valence", which would report
    # the absence of a field the document has no way to have.
    assert 'the file states no hydrogen count' in cml_line, cml_line
    assert ' or ' not in cml_line.rsplit('State it with', 1)[1], cml_line
    assert 'valence' not in cml_line.rsplit('State it with', 1)[1], cml_line


def test_every_registered_dialect_declares_its_own_channels():
    """The field is required rather than defaulted, and nothing may be left at the default.

    A `Dialect` declaring nothing would inherit `MOLFILE_CHANNELS`, right for the CTAB versions and
    wrong for every XML vocabulary, so forgetting is one misleading sentence in one log line.

    Not checked, deliberately: that two dialects spell their channels *differently*.  CML and MRV both
    say `hydrogenCount` and are both right, so a uniqueness assertion would forbid the true answer.
    """
    from ...ctfile._hydrogens import MOLFILE_CHANNELS

    names = dialects()
    assert len(names) >= 2, names
    for name in names:
        dial = dialect(name)
        assert dial.channels is not None, name
        assert dial.channels != MOLFILE_CHANNELS, (
            f'{name} left the molfile default, so its advice names an MDL data S-group')
        assert 'MRV_IMPLICIT_H' not in dial.channels.count, name


def test_the_channels_a_dialect_declares_are_what_the_record_carries():
    """`read_molecule` is the one place this is handed over, so a dialect cannot miss it.

    Asserted through `TOY`: a mechanism tested only on its two present users has not been shown to work
    for the third.
    """
    record, = _read_toy(TOY_DOC)
    assert record.ctab.channels == TOY.channels
    assert record.ctab.channels.count == 'an `hs` attribute'
    assert record.ctab.channels.valence is None
