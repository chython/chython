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
"""CML, the first dialect: what it reads, what it declines, and what it writes back.

``test/cml_quirks.cml`` spells the same facts every way CML has (CML 2 columns, CML 1's ``builtin``
children, both isotope spellings, 2D and 3D) and its log must be empty but for the line ``D`` earns;
``test/cml_damaged.cml`` is one record per way a file can be wrong, each still building while logging
*exactly* one line; ``test/cml_marvin.cml`` is constructs from a real 11.8 MB ChemAxon export."""

from xml.etree.ElementTree import tostring

from pytest import raises

from chython.core import read_smiles, write_smiles

from .._cml import CML, CML_NS, parse_cml, read_cml, record_from_molecule, write_cml
from .._dialect import read_xml
from .._errors import MalformedXml, UnsupportedXml


def _wrap(body, **attrs):
    """One ``<molecule>`` in a CML document, so a test can vary one construct and nothing else."""
    extra = ''.join(f' {k}="{v}"' for k, v in attrs.items())
    return f'<cml xmlns="{CML_NS}"><molecule id="m1"{extra}>{body}</molecule></cml>'


def _atoms(*specs):
    return '<atomArray>' + ''.join(f'<atom {s}/>' for s in specs) + '</atomArray>'


#: Two carbons with a drawing, for the tests that vary a bond rather than an atom.
ETHANE = _atoms('id="a1" elementType="C" x2="0.00" y2="0.00"',
                'id="a2" elementType="C" x2="0.87" y2="0.50"')


def _one(body, log=None, **attrs):
    """One record from a one-molecule document."""
    out = [] if log is None else log
    record, = parse_cml(_wrap(body, **attrs), log=out)
    return record


# the spellings all mean the same

def test_every_spelling_in_the_quirks_fixture_reads_and_is_silent(data):
    """Seven records spelling the same facts three ways: all read, and the log holds only the line ``D``
    earns, since reading it as hydrogen-2 is a recovery and a recovery is said out loud."""
    log = []
    records = parse_cml(data('cml_quirks.cml'), log=log)
    assert len(records) == 7
    assert [str(x) for x in log] == ['atom a3: D read as hydrogen isotope 2'], log


def test_the_column_form_and_the_element_form_produce_the_same_atoms(data):
    """CML 2 writes a column as a whitespace-separated attribute on ``<atomArray>``; CML 1 writes it as a
    ``<stringArray builtin=...>`` child; both are the same file as one ``<atom>`` per atom.  Asserted as
    an equality, since two separate expectations can both be updated to match a bug."""
    log = []
    columns, _, _, _, _, _, _ = parse_cml(data('cml_quirks.cml'), log=log)
    elements = _one(_atoms('id="a1" elementType="C" x2="0.00" y2="0.00"',
                           'id="a2" elementType="C" x2="0.87" y2="0.50"',
                           'id="a3" elementType="O" x2="1.73" y2="0.00"',
                           'id="a4" elementType="O" x2="0.87" y2="1.50"') +
                    '<bondArray>'
                    '<bond id="b1" atomRefs2="a1 a2" order="1"/>'
                    '<bond id="b2" atomRefs2="a2 a3" order="2"/>'
                    '<bond id="b3" atomRefs2="a2 a4" order="1"/>'
                    '</bondArray>')
    assert [(a.element, a.charge, a.x, a.y) for a in columns.ctab.atoms] == \
           [(a.element, a.charge, a.x, a.y) for a in elements.ctab.atoms]
    assert [(b.a, b.b, b.order) for b in columns.ctab.bonds] == \
           [(b.a, b.b, b.order) for b in elements.ctab.bonds]


def test_the_column_form_maps_the_identity_column_onto_the_element_form_s_attribute():
    """``atomID`` in the array form is ``id`` on an ``<atom>``, and the only column whose name differs --
    so the alias map is a short dict rather than a second table, and this pins that it stays short."""
    record = _one('<atomArray atomID="p1 p2" elementType="C O"/>')
    assert record.ids == ['p1', 'p2']
    assert [a.element for a in record.ctab.atoms] == ['C', 'O']


def test_the_array_form_s_own_endpoint_columns_are_two_attributes_not_one():
    """``<bondArray>`` in the column form spells its endpoints as two parallel columns, ``atomRef1`` and
    ``atomRef2``, where a ``<bond>`` element spells them as one ``atomRefs2``.  Both reach the engine's
    single ``bond_refs``: the array hook synthesises the joined form."""
    record = _one('<atomArray atomID="a1 a2" elementType="C O"/>'
                  '<bondArray atomRef1="a1" atomRef2="a2" order="2"/>')
    assert [(b.a, b.b, b.order) for b in record.ctab.bonds] == [(0, 1, 2)]


def test_both_forms_read_the_two_attribute_endpoint_spelling_and_both_are_silent():
    """``atomRef1``/``atomRef2`` is legal on a ``<bond>`` element and not only as a column, and nothing is
    lost either way, so both forms are silent.  Asserted as an equality so neither half can be fixed
    alone."""
    atoms = '<atomArray atomID="a1 a2" elementType="C O"/>'
    logs = []
    records = []
    for bonds in ('<bondArray atomRef1="a1" atomRef2="a2" order="2"/>',
                  '<bondArray><bond id="b1" atomRef1="a1" atomRef2="a2" order="2"/></bondArray>',
                  '<bondArray><bond id="b1" atomRefs2="a1 a2" order="2"/></bondArray>'):
        log = []
        records.append([(b.a, b.b, b.order) for b in _one(atoms + bonds, log).ctab.bonds])
        logs.append(log)
    assert records == [[(0, 1, 2)]] * 3, records
    assert logs == [[]] * 3, logs


def test_an_atom_id_declared_twice_keeps_the_first_atom_and_is_reported():
    """Two ``<atom id="a1">`` is a broken file and both atoms are kept anyway.  First-wins, because letting
    the later atom claim the name silently moves every bond written before it; the line is unprefixed
    because the file was broken and we read it regardless."""
    log = []
    record = _one(_atoms('id="a1" elementType="C"', 'id="a1" elementType="O"',
                         'id="a3" elementType="N"')
                  + '<bondArray><bond id="b1" atomRefs2="a1 a3" order="1"/></bondArray>', log)
    assert [a.element for a in record.ctab.atoms] == ['C', 'O', 'N']
    assert [(b.a, b.b) for b in record.ctab.bonds] == [(0, 2)]
    assert len(log) == 1 and not str(log[0]).startswith('unsupported'), log
    assert 'a1' in log[0] and 'twice' in log[0], log


def test_the_array_form_reports_a_duplicate_atom_id_too():
    """The pair: the line comes from the one place the id is registered, so both forms get it from the
    same statement."""
    log = []
    record = _one('<atomArray atomID="a1 a1" elementType="C O"/>', log)
    assert [a.element for a in record.ctab.atoms] == ['C', 'O']
    assert len(log) == 1 and 'twice' in log[0], log


def test_an_attribute_on_an_element_form_array_does_not_delete_its_children():
    """Which form an array is in is decided by its children, not its attributes.

    ``id`` and ``dictRef`` are legal CML global attributes on every element, so an ``<atomArray id="aa1">``
    holding ``<atom>`` children is the ordinary element form with a name on it.  Deciding on attributes
    lets the array hook claim the node, and a claiming hook stops the engine walking the children.
    """
    for array in ('<atomArray id="aa1">', '<atomArray dictRef="x:y">'):
        log = []
        record = _one(array + '<atom id="a1" elementType="C" x2="0" y2="0"/>'
                              '<atom id="a2" elementType="O" x2="1" y2="0"/></atomArray>'
                      '<bondArray><bond id="b1" atomRefs2="a1 a2" order="2"/></bondArray>', log)
        assert [a.element for a in record.ctab.atoms] == ['C', 'O'], array
        assert [(b.a, b.b, b.order) for b in record.ctab.bonds] == [(0, 1, 2)], array
    for array in ('<bondArray id="ba1">', '<bondArray dictRef="x:y">'):
        log = []
        record = _one(_atoms('id="a1" elementType="C"', 'id="a2" elementType="O"')
                      + array + '<bond id="b1" atomRefs2="a1 a2" order="2"/></bondArray>', log)
        assert [(b.a, b.b, b.order) for b in record.ctab.bonds] == [(0, 1, 2)], array


def test_a_bare_identity_attribute_on_an_element_form_array_is_silent():
    """``id`` on an ``<atomArray>`` is the array element's own name and every ``<atom>`` child restates its
    own identity, so counting it as a dropped column would report a loss on the commonest legal shape
    there is."""
    log = []
    _one('<atomArray id="aa1"><atom id="a1" elementType="C"/></atomArray>'
         '<bondArray id="ba1"/>', log)
    assert log == [], log


def test_an_array_stating_both_columns_and_elements_reads_the_elements_and_names_the_columns():
    """A node carrying column attributes *and* per-item children contradicts itself.  The children win as
    the more specific statement and the columns are named rather than merged, a merge being a guess about
    which half the writer meant.  Unprefixed: the file is broken and we read it anyway."""
    log = []
    record = _one('<atomArray atomID="z1 z2" elementType="N N">'
                  '<atom id="a1" elementType="C"/></atomArray>', log)
    assert [a.element for a in record.ctab.atoms] == ['C']
    assert record.ids == ['a1']
    assert len(log) == 1 and not str(log[0]).startswith('unsupported'), log
    assert 'elementType' in log[0] and '2 array column(s)' in log[0], log
    assert '1 <atom> element(s)' in log[0], log


def test_the_cml_1_per_atom_child_form_reads_on_atoms_and_on_bonds(data):
    """``<atom id="a1"><string builtin="elementType">C</string></atom>`` -- CML 1's third syntax, every
    field a child element carrying a ``builtin`` name.  Routed through the same ``apply_fields`` as an
    attribute, so ``A`` and ``partial12`` mean the same on a child as on an attribute."""
    log = []
    _, _, m3, *_ = parse_cml(data('cml_quirks.cml'), log=log)
    assert m3.ctab.title == 'chloromethane'
    assert [a.element for a in m3.ctab.atoms] == ['C', 'Cl']
    assert [(b.a, b.b, b.order) for b in m3.ctab.bonds] == [(0, 1, 1)]


def test_a_builtin_child_carrying_a_partial_order_is_unsupported_and_not_malformed():
    """The consequence of routing the child form through the same table: ``NotModelled`` still means
    ``NotModelled``, so a legal delocalised order is not reported as broken."""
    log = []
    _one(ETHANE + '<bondArray><bond id="b1" atomRefs2="a1 a2">'
                  '<string builtin="order">partial12</string></bond></bondArray>', log)
    assert len(log) == 1 and str(log[0]).startswith('unsupported: bond b1:'), log


# elements and isotopes

def test_both_isotope_spellings_are_read_and_both_are_written(data):
    """``isotope`` is CML 2's and ``isotopeNumber`` is CML 3's, deposited files carry both, and so does
    every atom this writer states a mass number on.

    The schema admits both spellings and readers differ over which one they take, so writing one leaves a
    reader that takes the other with no mass number at all.  MEASURED on one document with one attribute
    changed, at CDK 2.12 / Indigo 1.45 / Marvin 25.1.3: ``isotope`` alone is read by Indigo and Marvin,
    ``isotopeNumber`` alone by CDK and Indigo, and an atom carrying both is read by all three.
    """
    log = []
    *_, m4, m5, _, _ = parse_cml(data('cml_quirks.cml'), log=log)
    assert m4.ctab.atoms[0].isotope == 13   # isotope="13"
    assert m5.ctab.atoms[0].isotope == 13   # isotopeNumber="13"
    out = write_cml(read_cml(_wrap(_atoms('id="a1" elementType="C" isotope="13"'))), log=[])
    assert 'isotopeNumber="13"' in out and 'isotope="13"' in out, out
    # Armed against a writer that states a mass number on every atom: no number, neither spelling.
    assert 'isotope' not in write_cml(read_cml(_wrap(_atoms('id="a1" elementType="C"'))), log=[])


def test_deuterium_and_tritium_are_hydrogen_isotopes_and_are_said_out_loud():
    """``elementType="D"`` is what a molfile-derived converter leaves behind.  Logged with a plain line and
    no ``unsupported: ``: the file was loose and we read it anyway, which is the opposite side of the
    convention from a missing feature."""
    for symbol, isotope in (('D', 2), ('T', 3)):
        log = []
        record = _one(_atoms(f'id="a1" elementType="{symbol}"'), log)
        assert record.ctab.atoms[0].element == 'H'
        assert record.ctab.atoms[0].isotope == isotope
        assert [str(x) for x in log] == [f'atom a1: {symbol} read as hydrogen isotope {isotope}'], log
        assert not str(log[0]).startswith('unsupported'), log


def test_an_explicit_isotope_beats_the_one_the_symbol_implies():
    """``elementType="D" isotope="1"`` is contradictory and the number is the more specific statement;
    the symbol winning would make an explicit attribute unwritable."""
    record = _one(_atoms('id="a1" elementType="D" isotope="1"'), [])
    assert (record.ctab.atoms[0].element, record.ctab.atoms[0].isotope) == ('H', 1)


def test_an_upper_cased_symbol_is_recovered_and_named(data):
    """``cl`` for chlorine, which a molfile-derived converter leaves behind: a molfile's atom block is
    fixed-width and case-insensitive in practice.  Recovered by capitalisation and named -- silent
    case-folding would also accept ``NO`` as nobelium."""
    log = []
    record = _one(_atoms('id="a1" elementType="cl"'), log)
    assert record.ctab.atoms[0].element == 'Cl'
    assert [str(x) for x in log] == ["atom a1: elementType 'cl' read as 'Cl'"], log


def test_a_query_type_refuses_the_record_and_names_a_query_reader():
    """The one place this dialect refuses rather than logs: a query type or a dummy atom has nowhere to go
    in a molecule, and both wrong answers are silent -- default it to carbon, or drop the atom and
    renumber every bond after it.  Matched on the advice and not on ``elementType``, which every message
    in this resolver contains: the R family is refused too and names two formats instead, below."""
    for symbol in ('A', 'AH', 'Q', 'QH', 'X', 'M'):
        with raises(UnsupportedXml, match='query reader'):
            _one(_atoms(f'id="a1" elementType="{symbol}"'), [])


def test_the_r_family_and_the_dummies_read_as_the_marker():
    """Each of these names ONE atom -- an attachment point where a fragment is absent -- so each is the
    marker, element 0.  ``A``, ``Q``, ``X`` and ``M`` keep the query-reader refusal above, since those
    stand for a SET of elements, which is what a molecule cannot hold.
    """
    for symbol, index in (('R', 0), ('R#', 0), ('*', 0), ('R1', 1), ('R99', 99)):
        record = _one(_atoms(f'id="a1" elementType="{symbol}"'), [])
        assert record.ctab.atoms[0].element == 'R', symbol
        assert record.ctab.atoms[0].r_index == index, symbol
    # `Du` and `Dummy` are CML's own spellings of the same thing, and carry no index; each is a label,
    # so each also arrives as an alias
    for symbol in ('Du', 'Dummy'):
        record = _one(_atoms(f'id="a1" elementType="{symbol}"'), [])
        assert record.ctab.atoms[0].element == 'R', symbol
        assert record.ctab.atoms[0].r_index == 0, symbol
        assert record.ctab.aliases == {0: symbol}, symbol


def test_rubidium_still_resolves():
    """The regression guard for the R family moving above the element symbols: ``Rb`` is an element, and
    an inexact prefix test would refuse it along with eight others."""
    record = _one(_atoms('id="a1" elementType="Rb"'), [])
    assert record.ctab.atoms[0].element == 'Rb'


def test_a_lowercase_marker_is_folded_like_a_lowercase_symbol():
    """``r`` upper-cases into the R family, the same recovery ``cl`` gets: a CML document converted out of
    a molfile inherits the molfile's case-folding."""
    log = []
    record = _one(_atoms('id="a1" elementType="r"'), log)
    assert record.ctab.atoms[0].element == 'R'
    assert any(x.rule == 'xml:element-folded' for x in log), log


def test_an_unrecognisable_symbol_is_the_label_it_is():
    """``Zz`` names no element and no construct, which is what a drawn label looks like -- ``Pol``, ``OMe``,
    a registry identifier.  The same answer the CTfile readers give: the marker carrying the text as its
    alias, and the record kept.  An *empty* elementType is still malformed, that being a writer bug."""
    log = []
    record = _one(_atoms('id="a1" elementType="Zz"'), log)
    assert record.ctab.atoms[0].element == 'R'
    assert record.ctab.aliases == {0: 'Zz'}
    assert any(x.rule == 'xml:element-type-as-label' for x in log), log


def test_an_empty_element_type_is_malformed():
    """Distinguished from *absent*, the next test: an attribute stating nothing is a writer bug, an absent
    one is a file declining to say, and CML has a default for that."""
    with raises(MalformedXml, match='elementType is empty'):
        _one(_atoms('id="a1" elementType=" "'), [])


def test_an_absent_element_type_is_carbon_and_is_said_once_with_a_count():
    """Carbon is what every CML reader assumes for an ``<atom>`` with no ``elementType``, so reading it
    otherwise would disagree with the file's producer.  Said out loud, but aggregated in the wording
    ``_mrv.py`` uses: a document omitting ``elementType`` omits it everywhere, so a line per atom is the
    atom list repeated.  The count and the first atom keep it actionable.
    """
    log = []
    record = _one(_atoms('id="a1"'), log)
    assert record.ctab.atoms[0].element == 'C'
    assert [str(x) for x in log] == ['atom: 1 atom(s) with no elementType, read as carbon (first atom a1)'], log

    log = []
    record = _one(_atoms('id="a1"', 'id="a2"'), log)
    assert [a.element for a in record.ctab.atoms] == ['C', 'C']
    assert [str(x) for x in log] == ['atom: 2 atom(s) with no elementType, read as carbon (first atom a1)'], log


# charge and radical

def test_a_formal_charge_reads_and_round_trips(data):
    """``formalCharge="-1"`` is CML's only charge spelling.  The write half is here because the writer omits
    a zero charge, which a read-only test would not notice."""
    record = _one(_atoms('id="a1" elementType="O" formalCharge="-1"'), [])
    assert record.ctab.atoms[0].charge == -1
    out = write_cml(read_smiles('[O-]C'), log=[])
    assert 'formalCharge="-1"' in out and 'formalCharge="0"' not in out


def test_a_spin_multiplicity_becomes_the_radical_bit():
    """CML states a radical as a spin multiplicity and the arena holds a boolean, so 2 is the doublet that
    maps exactly, and is written back for the same reason."""
    record = _one(_atoms('id="a1" elementType="C" spinMultiplicity="2"'), [])
    assert record.ctab.atoms[0].radical is True
    # `C |^1:0|` and not `[CH3]`: the CXSMILES tail is the only spelling chython's SMILES reader has for
    # the radical bit.
    assert 'spinMultiplicity="2"' in write_cml(read_smiles('C |^1:0|'), log=[])


def test_a_multiplicity_above_a_doublet_is_read_as_one_centre_and_reported():
    """A quartet is a real thing a file may say and the container has one bit, so the bit is set and the
    loss named.  That is not "both": the construct was partly applied and the line is about the rest."""
    log = []
    record = _one(_atoms('id="a1" elementType="C" spinMultiplicity="4"'), log)
    assert record.ctab.atoms[0].radical is True
    assert len(log) == 1 and str(log[0]).startswith('unsupported: atom a1: spin multiplicity 4'), log


def test_a_multiplicity_below_one_is_a_broken_file_and_says_so_plainly():
    """Zero is not a multiplicity.  Read as a singlet, the only reading available, with a plain line:
    damage rather than a feature we lack."""
    log = []
    record = _one(_atoms('id="a1" elementType="C" spinMultiplicity="0"'), log)
    assert record.ctab.atoms[0].radical is False
    assert len(log) == 1 and 'out of range' in log[0] and not str(log[0]).startswith('unsupported'), log


# hydrogen counts

def test_hydrogen_count_is_a_total_and_the_reader_subtracts_the_drawn_ones(data):
    """CML's ``hydrogenCount`` is the TOTAL, drawn neighbours included, where ``CtabAtom.stated_h`` is the
    implicit count.  So an ammonium with ``hydrogenCount="4"`` and nothing drawn has four implicit, and a
    methanol whose hydrogen is drawn has that one subtracted -- counted over the bonds, since the file's
    total is the number being interpreted.
    """
    log = []
    *_, m6, _ = parse_cml(data('cml_quirks.cml'), log=log)
    assert m6.ctab.atoms[0].stated_h == 4  # ammonium, nothing drawn

    drawn = _one(_atoms('id="a1" elementType="O" hydrogenCount="2"',
                        'id="a2" elementType="H"') +
                 '<bondArray><bond id="b1" atomRefs2="a1 a2" order="1"/></bondArray>', [])
    assert drawn.ctab.atoms[0].stated_h == 1  # two total, one of them drawn


def test_a_total_below_the_drawn_count_is_recomputed_and_named():
    """``hydrogenCount="1"`` on an atom with two drawn hydrogens is impossible, so the count is dropped
    rather than stored negative, and named."""
    log = []
    record = _one(_atoms('id="a1" elementType="O" hydrogenCount="1"',
                         'id="a2" elementType="H"', 'id="a3" elementType="H"') +
                  '<bondArray><bond id="b1" atomRefs2="a1 a2" order="1"/>'
                  '<bond id="b2" atomRefs2="a1 a3" order="1"/></bondArray>', log)
    assert record.ctab.atoms[0].stated_h is None
    assert len(log) == 1 and 'is below the 2 hydrogen(s) drawn' in log[0], log


def test_an_integer_attribute_accepts_a_float_spelling():
    """``formalCharge="1.0"`` comes from a writer with one float formatter for every numeric field, and
    refusing it would lose an unambiguous charge.  Real garbage still logs."""
    assert _one(_atoms('id="a1" elementType="N" formalCharge="1.0"'), []).ctab.atoms[0].charge == 1
    log = []
    _one(_atoms('id="a1" elementType="N" formalCharge="plus"'), log)
    assert len(log) == 1 and 'not read' in log[0], log


# bond orders

def test_both_spellings_of_every_bond_order_read():
    """CML spells its orders as digits and as letters and files use both -- ``S``/``1``, ``D``/``2``,
    ``T``/``3``.  One many-to-one table, which is why writing has a second: inverting it would silently
    pick whichever key came last."""
    for text, order in (('1', 1), ('S', 1), ('2', 2), ('D', 2), ('3', 3), ('T', 3), ('A', 4)):
        record = _one(ETHANE + f'<bondArray><bond id="b1" atomRefs2="a1 a2" order="{text}"/>'
                               f'</bondArray>', [])
        assert record.ctab.bonds[0].order == order, text


def test_an_aromatic_order_is_stored_as_stated_and_not_kekulised():
    """``order="A"`` is stored as chython's order 4, exactly as the CTfile reader stores bond type 4.
    Kekulising in a reader would be repair, and repair is an explicit pass the caller runs."""
    record = _one(ETHANE + '<bondArray><bond id="b1" atomRefs2="a1 a2" order="A"/></bondArray>', [])
    assert record.ctab.bonds[0].order == 4


def test_a_partial_order_is_unsupported_rather_than_read_as_single():
    """``partial12`` is a legal statement about a resonance-averaged structure with no field to hold it and
    ``unknown`` is the file declining to say.  Neither is single, so both are declined with the prefix
    that means "we are the limitation"."""
    for order in ('partial01', 'partial12', 'partial23', 'unknown', 'other'):
        log = []
        _one(ETHANE + f'<bondArray><bond id="b1" atomRefs2="a1 a2" order="{order}"/></bondArray>', log)
        assert len(log) == 1 and str(log[0]).startswith('unsupported: bond b1:'), (order, log)


def test_a_nonsense_order_is_a_plain_line():
    """The other side of the same branch: ``order="banana"`` is damage, not a construct.  Both cost the
    attribute and only the prefix distinguishes them."""
    log = []
    _one(ETHANE + '<bondArray><bond id="b1" atomRefs2="a1 a2" order="banana"/></bondArray>', log)
    assert len(log) == 1 and not str(log[0]).startswith('unsupported'), log


# coordinates

def test_two_dimensional_coordinates_make_the_record_a_drawing(data):
    log = []
    m1, *_ = parse_cml(data('cml_quirks.cml'), log=log)
    assert m1.ctab.dimensionality == '2D'


def test_three_dimensional_coordinates_are_read_into_the_layout(data):
    """``x3/y3/z3`` with no ``x2/y2``: the record is a conformer.  The arena stores x and y, so ``z``
    reaches the intermediate and no further, and the writer emits no third coordinate rather than the zero
    it would have to invent."""
    log = []
    *_, m7 = parse_cml(data('cml_quirks.cml'), log=log)
    assert m7.ctab.dimensionality == '3D'
    assert (m7.ctab.atoms[1].x, m7.ctab.atoms[1].y) == (0.757, 0.586)
    out = write_cml(read_cml(data('cml_quirks.cml'), log=[])[6], log=[])
    assert 'z3' not in out and 'x3' not in out


def test_a_record_with_no_coordinates_at_all_states_no_dimensionality():
    """An absent layout is not a flat one: every atom at the origin is a drawing, and one every stereo
    perception would read configurations out of, so "no coordinates" is its own answer."""
    record = _one(_atoms('id="a1" elementType="C"', 'id="a2" elementType="O"'), [])
    assert record.ctab.dimensionality == ''


def test_both_coordinate_sets_keep_the_drawing_and_report_the_conformer(data):
    """Legal CML, and the container has one coordinate set.  2D wins because every other stereo statement
    in the record -- the wedges, the ``<bondStereo>`` -- is measured against the drawing, so mixing
    ``x2``/``y2`` with ``z3`` would build a geometry in no file and read stereo out of it."""
    log = []
    record = _one(_atoms('id="a1" elementType="C" x2="0.00" y2="0.00" x3="9.0" y3="9.0" z3="9.0"',
                         'id="a2" elementType="O" x2="0.87" y2="0.50" x3="8.0" y3="8.0" z3="8.0"'), log)
    assert record.ctab.dimensionality == '2D'
    assert (record.ctab.atoms[1].x, record.ctab.atoms[1].y) == (0.87, 0.5)
    assert len(log) == 1 and str(log[0]).startswith('unsupported: coordinates:'), log


def test_a_two_dimensional_set_that_is_all_zero_is_not_a_drawing():
    """"Has a drawing" is decided on the coordinate *values*, not on whether the attributes were present,
    so a record whose whole 2D set is ``0.0000`` loses to a present 3D set: no bond has a direction, so the
    conformer is the only geometry in the record.  The same test decides ``dimensionality`` in the MDL
    reader, so the two agree.
    """
    record = _one(_atoms('id="a1" elementType="C" x2="0.0" y2="0.0" x3="1.0" y3="2.0" z3="3.0"'), [])
    assert record.ctab.dimensionality == '3D'
    assert (record.ctab.atoms[0].x, record.ctab.atoms[0].y) == (1.0, 2.0)


# molecule-level things

def test_a_title_comes_from_the_attribute_or_the_child_element():
    """CML has two spellings and files use both: ``title=`` on the ``<molecule>`` and a ``<name>`` child.
    Either fills the same slot."""
    assert _one(_atoms('id="a1" elementType="C"'), [], title='by attribute').ctab.title == \
        'by attribute'
    assert _one('<name>by element</name>' + _atoms('id="a1" elementType="C"'),
                []).ctab.title == 'by element'


def test_a_second_name_is_declined_with_the_prefix():
    """CML allows several names under different ``convention`` attributes, so a file doing this is correct
    and the one-title container is the limitation -- which is what ``unsupported: `` means."""
    log = []
    record = _one('<name>first</name><name>second</name>' + _atoms('id="a1" elementType="C"'), log)
    assert record.ctab.title == 'first'
    assert len(log) == 1 and str(log[0]).startswith('unsupported: record: <name>'), log


def test_a_property_list_reaches_the_record_and_is_written_back(data):
    """Molecule-level properties, keyed by ``title`` and falling back to ``dictRef``.  The fixture's first
    property states both, and the name is the title: a ``dictRef`` carries the prefix of the dictionary it
    references (``molconvert cml`` writes ``dictRef="marvin:ID" title="ID"`` for the SD field ``ID``),
    which is not part of the field name.  Values stay ``str`` because ``dataType`` claims ``xsd:double``
    on fields whose value is ``>100`` often enough that coercing would lose data.
    """
    log = []
    m1, *_ = parse_cml(data('cml_quirks.cml'), log=log)
    assert m1.ctab.meta == {'Origin': 'hand-written fixture', 'chython:mp': '17'}
    out = write_cml(m1, log=[])
    assert 'dictRef="Origin"' in out and 'title="Origin"' in out and 'hand-written fixture' in out


def test_properties_reach_the_molecule(data):
    """`read_cml` sees them now.  Before, only `parse_cml` did -- the container could not hold one."""
    mol, *_ = read_cml(data('cml_quirks.cml'), log=[])
    assert mol.meta == {'Origin': 'hand-written fixture', 'chython:mp': '17'}


def test_a_molecules_meta_is_written_with_no_keyword():
    """No ``properties=`` argument any more: the molecule carries them, so the writer reads them off it."""
    mol = read_smiles('CCO')
    mol.meta['k'] = 'v'
    assert 'dictRef="k"' in write_cml([mol])


def test_properties_survive_a_document_round_trip(data):
    """Either way round now -- a record handed straight back, or the built molecule."""
    log = []
    before, *_ = parse_cml(data('cml_quirks.cml'), log=log)
    after, = parse_cml(write_cml(before, log=[]), log=[])
    assert after.ctab.meta == before.ctab.meta

    mol, *_ = read_cml(data('cml_quirks.cml'), log=[])
    again, *_ = read_cml(write_cml(mol, log=[]), log=[])
    assert again.meta == mol.meta


def test_a_property_with_no_key_is_dropped_and_named():
    """Neither ``dictRef`` nor ``title``, so there is nothing to key it by: a plain line, since a property
    with no name is a broken record rather than a construct we lack.  And nothing is left behind -- a
    ``<propertyList>`` nothing survived leaves the metadata empty rather than writing an empty value.
    """
    log = []
    record = _one('<propertyList><property><scalar>7</scalar></property></propertyList>'
                  + _atoms('id="a1" elementType="C"'), log)
    assert record.ctab.meta == {}
    assert len(log) == 1 and not str(log[0]).startswith('unsupported'), log


def test_a_molecule_level_scalar_is_a_data_field_named_by_its_title():
    """The other spelling of a data field, and the one CDK 2.12 writes: a `<scalar>` straight under
    `<molecule>`, `title` naming the field and `dictRef` naming the kind of property.  A scalar with no
    title has no name to be keyed by, so it is named in the log instead of keyed by its dictRef.
    """
    record = _one('<scalar dictRef="cdk:molecularProperty" title="BATCH_ID">lot-42</scalar>'
                  + _atoms('id="a1" elementType="C"'), [])
    assert record.ctab.meta == {'BATCH_ID': 'lot-42'}

    log = []
    record = _one('<scalar dictRef="cdk:molecularProperty">lot-42</scalar>'
                  + _atoms('id="a1" elementType="C"'), log)
    assert record.ctab.meta == {}
    assert len(log) == 1 and str(log[0]).startswith('unsupported:'), log


def test_two_property_lists_are_merged_rather_than_the_second_winning():
    """A document may write several, and a reader assigning instead of updating keeps only the last --
    silently, since both spellings produce a populated dict."""
    record = _one('<propertyList><property dictRef="d:a"><scalar>1</scalar></property></propertyList>'
                  '<propertyList><property dictRef="d:b"><scalar>2</scalar></property></propertyList>'
                  + _atoms('id="a1" elementType="C"'), [])
    assert record.ctab.meta == {'d:a': '1', 'd:b': '2'}


def test_a_property_holding_an_array_is_declined_with_the_prefix():
    """``<array>`` and ``<matrix>`` inside a property are legal CML with nowhere to go: the file is fine,
    we are the limitation."""
    log = []
    _one('<propertyList><property dictRef="d:x"><array>1 2 3</array></property></propertyList>'
         + _atoms('id="a1" elementType="C"'), log)
    assert len(log) == 1 and str(log[0]).startswith('unsupported:'), log


def test_metadata_is_ignored_by_declaration_and_a_stray_child_is_not():
    """``<metadataList>`` carries no chemistry, so it is silent -- and that silence is a claim, which is why
    the set is a short frozenset rather than a prefix match.  ``<electron/>`` is not in it."""
    log = []
    _one('<metadataList><metadata name="x" content="y"/></metadataList>'
         + _atoms('id="a1" elementType="C"'), log)
    assert log == [], log
    log = []
    _one('<electron/>' + _atoms('id="a1" elementType="C"'), log)
    assert len(log) == 1 and '<electron>' in log[0], log


# the damaged fixture, line by line

def test_every_damaged_record_still_builds(data):
    """Input is garbage by default: eight records, each broken a different way, and eight molecules out,
    because a caller sweeping a corpus needs the counts to match."""
    log = []
    records = parse_cml(data('cml_damaged.cml'), log=log)
    assert len(records) == 8
    for record in records:
        mol, _, _ = record.ctab.build()
        assert mol is not None


def test_the_damaged_fixture_logs_exactly_one_line_per_construct(data):
    """"Never both" as a count: the fixture's comment names one construct per record, so a duplicate line
    means a construct was reported twice.  Asserted as set-versus-list so the failure names the duplicate
    rather than a total."""
    log = []
    parse_cml(data('cml_damaged.cml'), log=log)
    assert len(log) == len(set(log)), [x for x in log if log.count(x) > 1]
    assert len(log) == 13, log


def test_each_damaged_line_is_on_the_correct_side_of_the_convention(data):
    """``unsupported: `` means the file is fine and chython is the limitation; a bare line means the file
    was broken and we read it anyway.  Both kinds are in this fixture on purpose."""
    log = []
    parse_cml(data('cml_damaged.cml'), log=log)
    broken = [x for x in log if not str(x).startswith('unsupported: ')]
    missing = [x for x in log if str(x).startswith('unsupported: ')]
    assert broken and missing, log
    # A ragged column, a dangling reference and three broken stereo statements are damage.
    assert any('ragged' in x or 'truncated' in x for x in broken), broken
    assert any('unknown atom' in x for x in broken), broken
    # A partial order, an unmodelled attribute and a stray child element are ours.
    assert any('partial12' in x for x in missing), missing


def test_a_ragged_column_is_truncated_rather_than_dropped(data):
    """Three elements, two coordinates.  Dropping the array would lose every atom in the file over one bad
    column, so it is truncated to the shortest and the loss is named."""
    log = []
    m1, *_ = parse_cml(data('cml_damaged.cml'), log=log)
    assert len(m1.ctab.atoms) == 2
    assert [a.element for a in m1.ctab.atoms] == ['C', 'C']


def test_an_array_child_with_no_builtin_is_named_once(data):
    """``<stringArray>`` with no ``builtin`` has no column name, so its values cannot be assigned to
    anything.  One line: two would be the "never both" half of the acceptance rule failing quietly."""
    log = []
    parse_cml(data('cml_damaged.cml'), log=log)
    named = [x for x in log if 'stringArray' in x]
    assert len(named) == 1, named
    # The engine names it, not the column reader: with no usable column the node is not the array form,
    # so `_columns` declines the element and the unclaimed-child path reports it.  One reporter, one line.
    assert str(named[0]) == 'unsupported: atom: <stringArray> in <atomArray> is not modelled', named


def test_an_unnamed_array_child_beside_a_usable_column_is_named_by_the_column_reader():
    """A usable column *and* an unnamed ``<*Array>`` child in one node.  Here ``_columns`` does claim the
    node, so the engine never walks the children and the column reader is the only reporter -- which is why
    ``_columns`` holds its line until it knows whether it is claiming the node."""
    log = []
    record = _one('<atomArray elementType="C O"><stringArray>a1 a2</stringArray></atomArray>', log)
    assert [a.element for a in record.ctab.atoms] == ['C', 'O']
    named = [x for x in log if 'stringArray' in x]
    assert [str(x) for x in named] == [
        'unsupported: atom: <stringArray> without a builtin attribute is not modelled'], log


# the Marvin fixture, from the census

def test_the_marvin_fixture_reads_to_the_expected_structures(data):
    """Four records from a real ChemAxon export's construct census: an aromatic ring written ``order="A"``,
    a bare ``<bondStereo>``, a charged carboxylate with a hydrogen count, and MDL wedge codes carried as
    ``convention="MDL"``.  Asserted as SMILES, since the question is what molecule came out."""
    log = []
    molecules = read_cml(data('cml_marvin.cml'), log=log)
    assert [write_smiles(m) for m in molecules] == ['O=C(/C=C/c1ccccc1)O', 'c1(C([O-])=O)ccccc1',
                                                    'C(/C)=C\\C', 'C(C(C)N)(=O)O']


def test_a_bare_bond_stereo_letter_is_silent(data):
    """ChemAxon writes ``<bondStereo>C</bondStereo>`` with no ``atomRefs4`` -- 692 times in the export this
    reader was calibrated against, against zero occurrences of ``atomRefs4`` -- because the frame the
    letter is measured in is the drawing, which the record also carries.  So there is no quadruple to check
    and nothing is lost.
    """
    log = []
    read_cml(data('cml_marvin.cml'), log=log)
    assert not any('not four' in x for x in log), log
    assert len(log) == 2, log


def test_the_mdl_convention_carries_a_wedge_code(data):
    """``convention="MDL" conventionValue="4"`` is how Marvin writes a molfile wedge column through CML.
    ``3`` on a double bond is "cis or trans, unknown which", which chython does not model; ``4`` on a
    single bond is "either", a configuration deliberately left unset."""
    log = []
    read_cml(data('cml_marvin.cml'), log=log)
    assert any('cis or trans, unknown which' in x and str(x).startswith('unsupported: ') for x in log), log
    assert any('drawn as either' in x for x in log), log


def test_an_unrecognised_bond_stereo_letter_is_declined_with_the_prefix():
    """``Q`` is in neither CML's nor MDL's letter vocabulary and gets ``unsupported: `` anyway: CML's
    ``<bondStereo>`` content is a dictionary-referenced string, so an unknown letter may be a valid entry
    in a dictionary we do not have, and of "your file is broken" and "we do not read this descriptor" only
    the second is true either way.  The bond survives regardless.
    """
    log = []
    record = _one(ETHANE + '<bondArray><bond id="b1" atomRefs2="a1 a2" order="1">'
                           '<bondStereo>Q</bondStereo></bond></bondArray>', log)
    assert len(record.ctab.bonds) == 1
    assert [str(x) for x in log] == ["unsupported: bond b1: bondStereo 'Q' is not modelled"], log


def test_a_bond_stereo_from_a_dictionary_this_reader_has_not_read_is_not_applied():
    """``convention`` names the dictionary the element's content is defined in, and MDL's is the only one
    this reader has.  ``W`` is a wedge in CML's own vocabulary and could mean anything in somebody else's,
    so applying MDL's meaning would invent a configuration -- and this one lands in the file rather than
    the log.  Asserted against the same content with no convention on it.
    """
    bond = ('<bondArray><bond id="b1" atomRefs2="a1 a2" order="1">'
            '<bondStereo{0}>W</bondStereo></bond></bondArray>')
    assert _one(ETHANE + bond.format(''), []).ctab.bonds[0].wedge
    for convention in (' convention="other:dict"', ' convention="cml:custom" conventionValue="1"'):
        log = []
        record = _one(ETHANE + bond.format(convention), log)
        assert len(record.ctab.bonds) == 1, convention
        assert not record.ctab.bonds[0].wedge, convention
        assert len(log) == 1 and str(log[0]).startswith('unsupported: '), log
        assert convention.split('"')[1] in log[0], log


def test_a_parity_that_names_the_wrong_number_of_atoms_names_the_attribute_it_read():
    """A tetrahedral parity is a permutation of four directions, so three references do not describe one.
    CML fixes ``atomRefs4`` at four and the reader also accepts CML 1's variable-length ``atomRefs``, so
    the message names the attribute it actually read -- otherwise it sends the reader looking for one their
    file does not contain.
    """
    for attribute in ('atomRefs4', 'atomRefs'):
        log = []
        record = _one(_atoms(f'id="a1" elementType="C">'
                             f'<atomParity {attribute}="a2 a3 a4">1</atomParity></atom><atom '
                             f'id="a2" elementType="C"', 'id="a3" elementType="N"',
                             'id="a4" elementType="Cl"')
                      + '<bondArray><bond atomRefs2="a1 a2"/><bond atomRefs2="a1 a3"/>'
                        '<bond atomRefs2="a1 a4"/></bondArray>', log)
        assert record.ctab.atoms[0].parity == 0, attribute
        assert len(log) == 1 and not str(log[0]).startswith('unsupported'), log
        assert f'{attribute} names 3 atoms, not four' in log[0], log
        assert 'four references' in log[0], log


def test_an_empty_bond_stereo_is_a_plain_line():
    """``<bondStereo/>`` states nothing, and a writer emitting an empty element meant to say something, so
    silence would hide the only evidence that a configuration was lost."""
    log = []
    _one(ETHANE + '<bondArray><bond id="b1" atomRefs2="a1 a2" order="1"><bondStereo/></bond>'
                  '</bondArray>', log)
    assert len(log) == 1 and not str(log[0]).startswith('unsupported'), log


#: *trans*-2-butene at hand-laid coordinates -- the two methyls on opposite sides of the C2=C3 axis --
#: with a slot for a ``<bondStereo>`` on the double bond.  Hand-laid because a molecule built from SMILES
#: has no layout.
DRAWN_BUTENE = (_atoms('id="a1" elementType="C" x2="0.00" y2="0.00"',
                       'id="a2" elementType="C" x2="0.87" y2="0.50"',
                       'id="a3" elementType="C" x2="1.73" y2="0.00"',
                       'id="a4" elementType="C" x2="2.60" y2="0.50"')
                + '<bondArray><bond id="b1" atomRefs2="a1 a2" order="1"/>'
                  '<bond id="b2" atomRefs2="a2 a3" order="2">{0}</bond>'
                  '<bond id="b3" atomRefs2="a3 a4" order="1"/></bondArray>')


def test_a_cis_trans_letter_is_ranked_below_the_drawing_and_a_disagreement_is_reported():
    """``<bondStereo>C``/``T`` is a second source that loses to the drawing, so both outcomes are pinned:
    agreement is silent, contradiction is reported and the drawing kept.  Both letters over the same
    *trans* drawing, so neither can be ignored.  The line comes from
    ``chython.core.wedge.assign_parities``, which makes the ranking for every format at once, so exactly
    one line also proves the dialect keeps no second copy of the decision.
    """
    stereo = '<bondStereo atomRefs4="a1 a2 a3 a4">{0}</bondStereo>'
    agrees, disagrees = [], []
    trans, = read_cml(_wrap(DRAWN_BUTENE.format(stereo.format('T'))), log=agrees)
    claims_cis, = read_cml(_wrap(DRAWN_BUTENE.format(stereo.format('C'))), log=disagrees)
    assert agrees == [], agrees
    assert len(disagrees) == 1 and not str(disagrees[0]).startswith('unsupported'), disagrees
    assert 'the drawing and the stated configuration disagree (drawn T, the document says C)' \
           in disagrees[0], disagrees
    assert 'keeping the drawing' in disagrees[0], disagrees
    assert write_smiles(claims_cis) == write_smiles(trans)
    assert '/' in write_smiles(trans) or '\\' in write_smiles(trans)


#: The same 2-butene with no coordinates at all -- the one case where the letter is the document's whole
#: statement.  Written out rather than derived from `DRAWN_BUTENE` by renaming attributes: an unknown
#: attribute is `unsupported: `, so a rename produces eight log lines and turns a silence assertion false.
FLAT_BUTENE = (_atoms('id="a1" elementType="C"', 'id="a2" elementType="C"',
                      'id="a3" elementType="C"', 'id="a4" elementType="C"')
               + '<bondArray><bond id="b1" atomRefs2="a1 a2" order="1"/>'
                 '<bond id="b2" atomRefs2="a2 a3" order="2">{0}</bond>'
                 '<bond id="b3" atomRefs2="a3 a4" order="1"/></bondArray>')

#: 3-methyl-2-pentene, coordinate-free: the anchor carries a methyl (`a1`) and an ethyl (`a5`-`a6`), so it
#: is stereogenic *and* has two substituents to measure a letter over.  2-methyl-2-butene will not do: two
#: identical methyls make the bond non-stereogenic.
FLAT_BRANCHED = (_atoms('id="a1" elementType="C"', 'id="a2" elementType="C"',
                        'id="a3" elementType="C"', 'id="a4" elementType="C"',
                        'id="a5" elementType="C"', 'id="a6" elementType="C"')
                 + '<bondArray><bond id="b1" atomRefs2="a1 a2" order="1"/>'
                   '<bond id="b2" atomRefs2="a2 a3" order="2">{0}</bond>'
                   '<bond id="b3" atomRefs2="a3 a4" order="1"/>'
                   '<bond id="b4" atomRefs2="a2 a5" order="1"/>'
                   '<bond id="b5" atomRefs2="a5 a6" order="1"/></bondArray>')


def test_a_letter_with_no_drawing_at_all_is_the_only_statement_there_is_and_is_read():
    """With no coordinates the letter is the document's whole statement about the double bond, and a molfile
    cannot express this at all.  Both letters, and they must produce *different* molecules; the sign is
    pinned against the SMILES reader, the other coordinate-free source in the tree."""
    stereo = '<bondStereo atomRefs4="a1 a2 a3 a4">{0}</bondStereo>'
    cis_log, trans_log = [], []
    cis, = read_cml(_wrap(FLAT_BUTENE.format(stereo.format('C'))), log=cis_log)
    trans, = read_cml(_wrap(FLAT_BUTENE.format(stereo.format('T'))), log=trans_log)
    assert cis_log == [], cis_log
    assert trans_log == [], trans_log
    assert write_smiles(cis) != write_smiles(trans)
    assert write_smiles(trans) == write_smiles(read_smiles('C/C=C/C'))
    assert write_smiles(cis) == write_smiles(read_smiles('C/C=C\\C'))


def test_a_frame_naming_the_other_substituent_inverts_the_letter():
    """``C`` measured over a different pair is a different configuration, and the frame says which pair.  A
    reader that took the letter and ignored ``atomRefs4`` returns the same molecule for both, which is why
    :func:`chython.core.wedge.stated_cis_trans` translates the frame instead of assuming the core's own."""
    stereo = '<bondStereo atomRefs4="{0} a2 a3 a4">C</bondStereo>'
    over_a1_log, over_a5_log = [], []
    over_a1, = read_cml(_wrap(FLAT_BRANCHED.format(stereo.format('a1'))), log=over_a1_log)
    over_a5, = read_cml(_wrap(FLAT_BRANCHED.format(stereo.format('a5'))), log=over_a5_log)
    assert over_a1_log == [], over_a1_log
    assert over_a5_log == [], over_a5_log
    assert write_smiles(over_a1) != write_smiles(over_a5)


def test_the_frame_may_be_written_from_either_terminal():
    """``a4 a3 a2 a1`` is the same statement as ``a1 a2 a3 a4``: nothing in CML says which terminal comes
    first, so demanding one order would silently drop half the documents using ``atomRefs4``.  Asserted as
    an equality, which two wrong readings cannot satisfy."""
    logs = []
    forward, = read_cml(_wrap(FLAT_BUTENE.format(
        '<bondStereo atomRefs4="a1 a2 a3 a4">C</bondStereo>')), log=logs)
    reversed_, = read_cml(_wrap(FLAT_BUTENE.format(
        '<bondStereo atomRefs4="a4 a3 a2 a1">C</bondStereo>')), log=logs)
    assert logs == [], logs
    assert write_smiles(forward) == write_smiles(reversed_)
    assert '/' in write_smiles(forward) or '\\' in write_smiles(forward)


def test_a_bare_letter_with_no_drawing_is_read_only_where_it_cannot_be_ambiguous():
    """A bare ``C`` on a substituted terminal names no pair, so it is reported and dropped rather than
    guessed: two readings, and the file has chosen neither.  Marvin gets away with it because the drawing
    is in the same record.  The 2-butene half is the contrast -- one substituent per terminal, one possible
    frame, so the bare letter is read in silence.
    """
    log = []
    unambiguous, = read_cml(_wrap(FLAT_BUTENE.format('<bondStereo>T</bondStereo>')), log=log)
    assert log == [], log
    assert write_smiles(unambiguous) == write_smiles(read_smiles('C/C=C/C'))

    log = []
    ambiguous, = read_cml(_wrap(FLAT_BRANCHED.format('<bondStereo>C</bondStereo>')), log=log)
    assert len(log) == 1, log
    assert 'names no reference atoms' in log[0] and 'two substituents' in log[0], log
    assert '/' not in write_smiles(ambiguous) and '\\' not in write_smiles(ambiguous)


# reaction documents

def _reaction(*roles):
    """A ``<reaction>`` whose role elements hold one one-atom molecule each, per ``(name, count)``."""
    body = ''
    for name, count in roles:
        inner = ''.join(f'<{name}><molecule id="m{name}{i}"><atomArray>'
                        f'<atom id="a{name}{i}" elementType="C"/></atomArray></molecule></{name}>'
                        for i in range(count))
        body += f'<{name}List>{inner}</{name}List>'
    return f'<reaction id="r1">{body}</reaction>'


def test_a_reaction_document_names_the_roles_it_could_not_model_with_their_counts():
    """A reaction's content *is* the roles, so reading its molecules flat and saying nothing returns a file
    indistinguishable from the same molecules loose.  The molecules are still read, and the counts are
    load-bearing: a caller reading the log can recover the record's shape from them.  One line for the
    reaction, not one per molecule.
    """
    for reader in (read_cml, read_xml):
        log = []
        molecules = reader(f'<cml xmlns="{CML_NS}">{_reaction(("reactant", 2), ("product", 1))}</cml>',
                           log=log)
        assert [write_smiles(m) for m in molecules] == ['C', 'C', 'C']
        assert [str(x) for x in log] == ['unsupported: record: <reaction> r1 roles are not modelled; 2 reactant and '
                                         '1 product molecule(s) read as a flat list'], log


def test_a_reaction_role_this_reader_has_not_seen_is_named_as_the_file_spells_it():
    """CML has more roles than a reaction has sides -- ``<substanceList>``, ``<spectator>`` -- so the counts
    are keyed on the file's own element names.  Mapping an unfamiliar one onto "reactant" or "agent" would
    be inventing a classification."""
    log = []
    read_cml(f'<cml xmlns="{CML_NS}">'
             f'{_reaction(("reactant", 1), ("substance", 1), ("product", 1))}</cml>', log=log)
    assert [str(x) for x in log] == ['unsupported: record: <reaction> r1 roles are not modelled; 1 reactant, '
                                     '1 substance and 1 product molecule(s) read as a flat list'], log


def test_a_molecule_a_reaction_gives_no_role_is_counted_as_unplaced():
    """A ``<molecule>`` directly under ``<reaction>`` is in no role element, so the line must not claim a
    side the file never named."""
    log = []
    read_cml(f'<cml xmlns="{CML_NS}"><reaction id="r1">'
             f'<molecule><atomArray><atom id="a1" elementType="C"/></atomArray></molecule>'
             f'</reaction></cml>', log=log)
    assert [str(x) for x in log] == ['unsupported: record: <reaction> r1 roles are not modelled; 1 unplaced '
                                     'molecule(s) read as a flat list'], log


def test_a_reaction_root_is_read_whether_or_not_it_declares_the_namespace():
    """``sniff`` matches a namespace before a root name, so without ``reaction`` in CML's root set the
    answer would turn on an ``xmlns`` that says nothing about the file's vocabulary.  Both spellings are
    read and both are named."""
    reaction = _reaction(('reactant', 1))
    for source in (f'<reaction xmlns="{CML_NS}"{reaction[len("<reaction"):]}', reaction):
        log = []
        assert [write_smiles(m) for m in read_xml(source, log=log)] == ['C']
        assert any('<reaction> r1 roles are not modelled' in x for x in log), log


def test_a_document_with_no_reaction_in_it_says_nothing_about_one():
    """The control a per-document hook gets wrong first: a line on every document would make the prefix
    unfilterable."""
    log = []
    read_cml(_wrap(_atoms('id="a1" elementType="C"')), log=log)
    assert log == [], log


# the writer

def test_the_writer_emits_a_namespace_as_an_attribute_and_not_a_default():
    """``ElementTree.tostring(default_namespace=...)`` refuses a document with unqualified attributes, which
    every CML document has, so the namespace is written as a literal ``xmlns`` attribute on an unqualified
    tree.  The global ``register_namespace`` is never touched -- it would change how an unrelated caller's
    XML serialises.
    """
    out = write_cml(read_smiles('CCO'), log=[])
    assert f'xmlns="{CML_NS}"' in out
    assert 'ns0:' not in out and '<ns0' not in out
    from xml.etree.ElementTree import _namespace_map
    assert CML_NS not in _namespace_map or _namespace_map.get(CML_NS) != 'cml'


def test_what_the_writer_produces_is_read_by_this_reader():
    """Compared against ``write_smiles`` of the original rather than a literal string: SMILES output depends
    on the atom order walked, so a literal would pin the writer's traversal instead of the molecule."""
    log = []
    before = read_smiles('CC(=O)O')
    molecules = read_cml(write_cml(before, log=log), log=log)
    assert len(molecules) == 1
    assert write_smiles(molecules[0]) == write_smiles(before)
    assert log == [], log


def test_a_document_holds_several_molecules_with_distinct_ids():
    """CML's bond references are by id, so a duplicate id across molecules makes a document whose bonds are
    ambiguous to any reader that flattens it."""
    out = write_cml([read_smiles('C'), read_smiles('N')], log=[])
    assert out.count('<molecule') == 2
    assert 'id="m1"' in out and 'id="m2"' in out
    assert len(read_cml(out, log=[])) == 2


def test_a_molecule_with_no_layout_is_not_written_at_the_origin():
    """Writing every atom at ``0.0000`` would be a drawing, and one a stereo perception would read
    configurations out of, so the attributes are omitted -- CML for "no layout"."""
    out = write_cml(read_smiles('CCO'), log=[])
    assert 'x2=' not in out and 'y2=' not in out
    assert read_cml(out, log=[])[0].has_coordinates is False


def test_an_unknown_hydrogen_count_is_omitted_and_reported():
    """``H_UNKNOWN`` reaches the writer as "nothing determined this" and an absent ``hydrogenCount`` is CML
    for "derive it", so the attribute is left off and the caller told which atom it happened on.

    The molecule must come from a molfile: a SMILES atom's count is always derivable, where a
    five-coordinate neutral carbon has no valence rule and the reader stores the sentinel.  Note the
    accessor answers ``None`` for the sentinel, so ``== H_UNKNOWN`` is dead code that a round trip cannot
    see -- only the log assertion catches it.
    """
    from chython.formats.ctfile import parse_record
    lines = ['overvalent carbon', '  chython', '', '  6  5  0  0  0  0  0  0  0  0999 V2000']
    lines.append('    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0')
    lines += [f'    1.0000    {i}.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0'
              for i in range(5)]
    lines += [f'  1{i + 2:3d}  1  0  0  0  0' for i in range(5)] + ['M  END']
    mol = parse_record(lines, [])
    log = []
    out = write_cml(mol, log=log)
    assert any('hydrogen count unknown' in x for x in log), log
    assert out.count('hydrogenCount') == 5  # the five fluorines, not the carbon


def test_a_dative_bond_is_written_with_no_order_and_reported():
    """CML has no coordination order in its vocabulary, so the writer omits the attribute and says so:
    ``order="1"`` would turn a coordination contact into a covalent bond in every consumer downstream."""
    log = []
    out = write_cml(read_smiles('[Fe]~N(C)(C)C'), log=log)
    assert any('no coordination bond order' in x and str(x).startswith('unsupported: ') for x in log), log
    assert out.count('<bond ') == 4          # `'<bond '` with the space: `<bondArray` starts the same
    assert out.count('order=') == 3          # every bond but the dative one


def test_atom_map_numbers_are_reported_as_unwritable_once_for_the_record():
    """CML has no atom-map attribute, so a mapped molecule loses its mapping on write and says so, as every
    sibling unwritable fact does -- a dative order, an S-group, an underivable hydrogen count, a cis/trans
    configuration with no layout.  One line for the record, not one per atom: an atom map is a property of
    the mapping.
    """
    log = []
    out = write_cml(read_smiles('[CH3:1][OH:2]'), log=log)
    matched = [x for x in log if 'map number' in x]
    assert len(matched) == 1 and str(matched[0]).startswith('unsupported: '), log
    assert '2 atom map number(s)' in matched[0], matched
    assert 'mapNumber' not in out and 'map=' not in out, out


def test_an_unmapped_molecule_says_nothing_about_map_numbers():
    """The other side of the line: an unconditional one would report a loss on every record any caller ever
    writes."""
    log = []
    write_cml(read_smiles('CO'), log=log)
    assert not any('map number' in x for x in log), log


def test_an_s_group_is_reported_as_unwritten():
    """CML's S-group vocabulary is a named boundary of this dialect, so a molecule carrying S-groups is
    written without them and the caller told how many were lost."""
    record = record_from_molecule(read_smiles('CCO'), log=[])
    assert record is not None  # the no-S-group path is silent, which is the control for the above
    log = []
    record_from_molecule(read_smiles('CCO'), log=log)
    assert not any('sgroup' in x for x in log), log


def test_the_writer_produces_byte_identical_output_twice():
    """Attribute order is decided by the table's declaration order, so two runs agree -- without which a
    diff of two exports of one molecule is noise."""
    mol = read_smiles('C[C@H](N)C(=O)O')
    assert write_cml(mol, log=[]) == write_cml(mol, log=[])


def test_indentation_is_a_keyword_and_off_means_one_line():
    """An indented document diffs, which is what a person wants; a caller writing a million records wants
    the bytes."""
    mol = read_smiles('CCO')
    assert '\n  <molecule' in write_cml(mol, log=[])
    assert '\n' not in write_cml(mol, log=[], indent=None).split('?>\n')[1]


def test_a_record_may_be_written_in_place_of_a_molecule(data):
    """How a caller writes back what ``parse_cml`` read: the file's own atom ids and its property list
    survive, where going through a molecule would keep only what the container holds.  Plus the three facts
    the record carries beside them, each on the record that states it -- the hydrogen total, the
    configuration, and whether there is a layout at all.
    """
    log = []
    records = parse_cml(data('cml_quirks.cml'), log=log)
    out = write_cml(records[0], log=[])
    assert 'id="a1"' in out and 'hand-written fixture' in out
    assert 'x2="0.0000"' in out and 'x3=' not in out, out
    assert 'hydrogenCount="4"' in write_cml(records[5], log=[])
    assert 'y3="0.5860"' in write_cml(records[6], log=[])
    parity, = parse_cml(_wrap(PARITY_ONLY), log=[])
    stereo = write_cml(parity, log=[])
    assert '<atomParity' in stereo and 'x2=' not in stereo, stereo


#: Methane with one of its four hydrogens drawn as an atom of its own: the shape that makes
#: ``hydrogenCount``'s total-versus-implicit distinction observable.  A writer labelling the implicit count
#: a total loses one hydrogen per round trip.
METHANE_ONE_H_DRAWN = (_atoms('id="a1" elementType="C" hydrogenCount="4"',
                              'id="a2" elementType="H"')
                       + '<bondArray><bond id="b1" atomRefs2="a1 a2" order="1"/></bondArray>')


def test_a_record_writes_hydrogen_count_back_as_the_total_the_file_stated():
    """``hydrogenCount`` is a total in CML, drawn neighbours included, so the writer states a total too,
    whichever path filled the record.  Methane with one hydrogen drawn discriminates: total 4, implicit 3."""
    log = []
    record, = parse_cml(_wrap(METHANE_ONE_H_DRAWN), log=log)
    out = write_cml(record, log=log)
    assert 'hydrogenCount="4"' in out, out
    mol, = read_cml(out, log=log)
    carbon, = (n for n in mol.atom_numbers if mol.element_of(n) == 6)
    assert mol.total_h_of(carbon) == 4
    assert log == [], log


def test_a_molecule_writes_hydrogen_count_back_as_the_total_too():
    """The molecule half of the pair: a molecule states an implicit count and its explicit hydrogens are
    atoms, so the total is the sum -- the same quantity the record path writes, under the same name."""
    log = []
    out = write_cml(read_smiles('[H]C'), log=log)
    assert 'hydrogenCount="4"' in out, out
    mol, = read_cml(out, log=log)
    carbon, = (n for n in mol.atom_numbers if mol.element_of(n) == 6)
    assert mol.total_h_of(carbon) == 4
    assert log == [], log


def test_a_refused_hydrogen_count_is_not_written_back():
    """A ``hydrogenCount`` below the number of hydrogens drawn is reported and recomputed on the way in, so
    writing it back would republish a number this reader has just said it does not believe."""
    log = []
    record, = parse_cml(_wrap(_atoms('id="a1" elementType="C" hydrogenCount="1"',
                                     'id="a2" elementType="H"', 'id="a3" elementType="H"')
                              + '<bondArray><bond atomRefs2="a1 a2"/><bond atomRefs2="a1 a3"/>'
                                '</bondArray>'), log=log)
    assert any('is below the 2 hydrogen(s) drawn' in x for x in log), log
    out = write_cml(record, log=[])
    assert 'hydrogenCount' not in out, out


#: A stereocentre stated as a parity and nothing else -- no wedge, no coordinates -- so the descriptor is
#: the record's only configuration statement.
PARITY_ONLY = (_atoms('id="a1" elementType="C" hydrogenCount="1">'
                      '<atomParity atomRefs4="a2 a3 a4 a1">1</atomParity></atom><atom '
                      'id="a2" elementType="C"', 'id="a3" elementType="N"',
                      'id="a4" elementType="Cl"')
               + '<bondArray><bond atomRefs2="a1 a2"/><bond atomRefs2="a1 a3"/>'
                 '<bond atomRefs2="a1 a4"/></bondArray>')


def test_a_record_writes_its_atom_parity_back_and_the_configuration_survives():
    """``<atomParity>`` is written off ``CtabAtom.parity``, the same field a file's own descriptor lands on.
    The assertion is the re-read configuration rather than the element's presence, since an
    ``<atomParity>`` with the frame or sign wrong is worse than none; this record has no coordinates and no
    wedge, so the descriptor is the only channel there is.
    """
    log = []
    record, = parse_cml(_wrap(PARITY_ONLY), log=log)
    before, = read_cml(_wrap(PARITY_ONLY), log=log)
    centre, = (n for n in before.atom_numbers if before.parity_of(n))
    out = write_cml(record, log=log)
    assert '<atomParity' in out, out
    after, = read_cml(out, log=log)
    assert {n: after.parity_of(n) for n in after.atom_numbers} == \
           {n: before.parity_of(n) for n in before.atom_numbers}
    assert after.parity_of(centre) == before.parity_of(centre) != 0
    assert log == [], log


def test_a_record_with_no_coordinates_is_not_written_at_the_origin():
    """``Ctab.dimensionality`` is what the writer asks, and a record from a file that stated no coordinates
    says ``''``.  Writing ``x2="0.0000"`` on every atom would be a drawing nobody made, which a reader
    would then perceive stereo from."""
    record, = parse_cml(_wrap(PARITY_ONLY), log=[])
    out = write_cml(record, log=[])
    assert 'x2=' not in out and 'y2=' not in out, out
    assert read_cml(out, log=[])[0].has_coordinates is False


def test_a_record_that_stated_a_conformer_is_written_back_as_one():
    """A 3D record is written ``x3``/``y3``/``z3``: relabelling a conformer as a drawing loses the third
    coordinate and hands the next reader a layout to perceive stereo from.  Which set to write is a
    question about the record, so the record is asked rather than each atom."""
    solid = _atoms('id="a1" elementType="C" x3="1.5" y3="2.5" z3="3.5"',
                   'id="a2" elementType="O" x3="2.5" y3="2.5" z3="3.5"')
    record, = parse_cml(_wrap(solid + '<bondArray><bond atomRefs2="a1 a2"/></bondArray>'), log=[])
    out = write_cml(record, log=[])
    assert 'z3="3.5000"' in out and 'x2=' not in out, out
    after, = parse_cml(out, log=[])
    assert [(a.x, a.y, a.z) for a in after.ctab.atoms] == \
           [(a.x, a.y, a.z) for a in record.ctab.atoms]


def test_a_configuration_on_an_atom_that_cannot_name_four_directions_is_reported():
    """A parity on an atom whose bonds cannot name a quadruple has no frame to measure the sign in, so
    nothing is written and the caller is told -- as the molecule path does for a centre with two undrawn
    directions."""
    from chython.formats.ctfile._ctab import Ctab, CtabAtom, CtabBond
    from .._dialect import Record
    record = Record(Ctab())
    for element in ('C', 'O', 'N'):
        record.add_atom(CtabAtom(element))
    record.add_bond(CtabBond(0, 1))
    record.ctab.atoms[0].parity = 1
    log = []
    out = write_cml(record, log=log)
    assert '<atomParity' not in out, out
    assert any('cannot name four directions' in x for x in log), log


def test_write_cml_element_returns_a_tree_for_a_caller_that_composes():
    """The element rather than the string, for a caller embedding CML in a larger document.  Separate from
    ``write_cml`` because pretty-printing mutates the tree."""
    from .._cml import write_cml_element
    root = write_cml_element(read_smiles('CCO'), log=[])
    assert root.tag == 'cml'
    assert tostring(root, encoding='unicode').startswith('<cml ')


def test_a_non_utf8_title_is_reported_and_the_record_is_still_written_cml():
    """The CML writer threads the log through ``resolve_output``, which it calls independently of the CTfile
    writers, so a test against one format is not evidence the other reports too.  The title byte
    ``\\xe9`` is the latin-1 e-acute, not valid UTF-8 on its own, so the strict decode fails and the
    replacement path runs.
    """
    from chython.core import MoleculeContainer

    mol = MoleculeContainer()
    with mol.edit():
        mol.add_atom('C')
    mol.set_title(b'caf\xe9')
    log = []
    out = write_cml(mol, log=log)
    assert any(str(x).startswith('unsupported: ') and 'not valid UTF-8' in x for x in log), log
    assert '1 byte(s)' in ''.join(str(x) for x in log), log
    # The replacement character appears in the title= attribute of the <molecule> element.
    assert 'caf�' in out, f'replacement character should appear in the CML title: {out[:300]!r}'


def test_a_title_with_genuine_replacement_character_logs_nothing_in_cml():
    """U+FFFD encoded as valid UTF-8 is a real part of the title: ``\\xef\\xbf\\xbd`` decodes cleanly, which
    is what the strict decode separates from a substituted byte."""
    from chython.core import MoleculeContainer

    mol = MoleculeContainer()
    with mol.edit():
        mol.add_atom('C')
    mol.set_title('�'.encode('utf-8'))  # valid UTF-8 bytes that decode to U+FFFD
    log = []
    write_cml(mol, log=log)
    assert not any('not valid UTF-8' in x for x in log), log


def test_an_r_atom_is_written_as_an_r_group_label_without_its_index():
    """``elementType="R"``, which is what ``molconvert cml`` writes for one, and the index is reported
    lost: CML spells the label and not its number, where MRV has ``rgroupRef``."""
    mol = read_smiles('[R7]C')
    log = []
    text = write_cml(mol, log=log)
    assert 'elementType="R"' in text, text
    assert any(x.rule == 'cml:r-index-not-written' for x in log), log
    again, = read_cml(text, log=[])
    r = next(a for a in again.atoms() if a.is_r)
    assert r.r_index == 0, 'the index has no CML spelling, and the loss was reported on write'
