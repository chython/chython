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
"""The MRV dialect: what it reads, what it refuses, and what it says it lost.

A silence has to be armed -- an assertion on an empty log passes against a reader that read nothing -- so
every one here sits beside one on the same construct's content.  And a round trip proves only that the
two halves agree, so every codec is also checked against MRV's literal attribute text.
"""

from pytest import raises

from .._dialect import dialect, read_xml, sniff
from .._errors import MalformedXml, UnsupportedXml
from .._mrv import MRV_NS, parse_mrv, read_mrv, record_from_molecule, write_mrv
from .._tree import parse_xml
from ...ctfile import mol
from ...ctfile._ctab import WEDGE_DOWN, WEDGE_UP
from ....core import read_smiles, write_smiles
from ....core.wedge import wedges_for_write


def _wrap(body, molecule_attributes=''):
    """`body` inside the ``<cml><MDocument><MChemicalStruct><molecule>`` wrapper every MRV file has."""
    return (f'<cml xmlns="{MRV_NS}"><MDocument><MChemicalStruct>'
            f'<molecule molID="m1" {molecule_attributes}>{body}</molecule>'
            f'</MChemicalStruct></MDocument></cml>')


def _wrap_document(body):
    """`body` where a ``<molecule>`` would sit, for the constructs that are not one."""
    return (f'<cml xmlns="{MRV_NS}"><MDocument><MChemicalStruct>{body}'
            f'</MChemicalStruct></MDocument></cml>')


def _atoms(*atoms):
    return '<atomArray>' + ''.join(atoms) + '</atomArray>'


def _unsupported(log):
    return [x for x in log if str(x).startswith('unsupported: ')]


# the containers and the wrapper

def test_the_document_wrapper_needs_no_handler():
    """``<MDocument>`` and ``<MChemicalStruct>`` are pure containers and must not acquire code.

    The engine reads through whatever lies between the root and a ``<molecule>``, and *free* means no log
    line either: a container stating nothing about its molecules loses nothing when a flat list replaces
    it.  Armed by the reaction test below, where a document element that does state something is reported.
    """
    log = []
    mol, = read_mrv(_wrap(_atoms('<atom id="a1" elementType="O"/>')), log=log)
    assert mol.atom_count == 1
    assert log == [], log


def test_the_root_is_called_cml_and_the_namespace_is_what_routes_it():
    """The same root name as CML, so only the namespace can decide.

    Asserted on which dialect was chosen and not on an empty log: both documents are clean, so both read
    silently whichever dialect took them.  Both halves are needed -- either alone is satisfied by a reader
    that ignores namespaces and always answers the same dialect.  The second document is a Marvin idiom in
    the CML namespace, so a dialect routed on "this looks like Marvin" would steal a real file from CML.
    """
    body = _atoms('<atom id="a1" elementType="O"/>')
    assert dialect('mrv').tags.roots >= {'cml'}
    log = []
    assert sniff(parse_xml(_wrap(body)), log).name == 'mrv'
    read_xml(_wrap(body), log=log)
    assert log == [], log

    cml_document = (f'<cml xmlns="http://www.xml-cml.org/schema"><molecule id="m1">{body}</molecule>'
                    f'</cml>')
    log = []
    assert sniff(parse_xml(cml_document), log).name == 'cml'
    read_xml(cml_document, log=log)
    assert log == [], log


# the atom vocabulary

def test_the_atom_attributes_are_read_from_their_own_spellings():
    """One assertion per row of the atom table, against MRV's literal attribute text.

    The five that carry chemistry, each spelled as a Marvin file spells it.  Every one has a row, so a
    line about any of them means the row is not reached.
    """
    log = []
    mol, = read_mrv(_wrap(_atoms(
        '<atom id="a1" elementType="N" formalCharge="1" hydrogenCount="3" mrvMap="7"/>',
        '<atom id="a2" elementType="C" isotope="13" radical="monovalent" hydrogenCount="3"/>',
    ) + '<bondArray><bond id="b1" atomRefs2="a1 a2" order="1"/></bondArray>'), log=log)
    n, c = mol.atom_numbers
    assert mol.charge_of(n) == 1
    assert mol.implicit_h_of(n) == 3
    assert mol.map_number_of(n) == 7
    assert mol.isotope_of(c) == 13
    assert mol.radical_of(c)
    assert log == [], log


def test_a_radical_of_more_than_one_electron_is_read_as_one_and_reported():
    """``divalent`` is two unpaired electrons and a molecule holds one radical bit, so the bit is set and
    the count is named.  ``monovalent`` is the armed control: same bit, no line."""
    log = []
    mol, = read_mrv(_wrap(_atoms('<atom id="a1" elementType="C" radical="divalent"/>')), log=log)
    assert mol.radical_of(next(iter(mol.atom_numbers)))
    assert [x for x in _unsupported(log) if 'unpaired electrons' in x], log

    log = []
    read_mrv(_wrap(_atoms('<atom id="a1" elementType="C" radical="monovalent"/>')), log=log)
    assert not [x for x in log if 'unpaired electrons' in x], log


def test_a_radical_name_the_format_does_not_have_costs_the_attribute_and_not_the_atom():
    """An unknown ``radical`` value is a malformed attribute: the atom stays, unprefixed line.

    ``0`` and ``-`` are the two nulls a column writes for an atom that states nothing, so neither may be
    reported and neither may make a radical.
    """
    log = []
    mol, = read_mrv(_wrap(_atoms('<atom id="a1" elementType="C" radical="pentavalent"/>')), log=log)
    assert mol.atom_count == 1
    assert not mol.radical_of(next(iter(mol.atom_numbers)))
    assert log and not _unsupported(log), log

    for null in ('0', '-'):
        log = []
        mol, = read_mrv(_wrap(_atoms(f'<atom id="a1" elementType="C" radical="{null}"/>')), log=log)
        assert not mol.radical_of(next(iter(mol.atom_numbers)))
        assert log == [], log


def test_the_whole_r_family_reads_as_the_marker():
    """``R``, ``R#``, ``R<n>`` and ``*`` all name one atom -- a fragment's attachment point -- so each is
    the marker, element 0.  ``elementType="R" rgroupRef="1"`` is Marvin's own spelling: it is what
    ``molconvert mrv`` writes for an ``M  RGP`` atom, so this is the round trip with Marvin and not a
    convention of ours.  MRV has no set-valued query type, which is why nothing here refuses.
    """
    for symbol, index in (('R', 0), ('R#', 0), ('*', 0), ('R1', 1), ('R99', 99)):
        mol, = read_mrv(_wrap(_atoms(f'<atom id="a1" elementType="{symbol}"/>')), log=[])
        atom = next(iter(mol.atoms()))
        assert atom.is_r, symbol
        assert atom.r_index == index, symbol


def test_the_index_travels_in_rgroup_ref():
    """The attribute Marvin puts it in, and the one the writer puts it back in."""
    mol, = read_mrv(_wrap(_atoms('<atom id="a1" elementType="R" rgroupRef="7"/>')), log=[])
    assert next(iter(mol.atoms())).r_index == 7
    text = write_mrv(mol)
    assert 'elementType="R"' in text and 'rgroupRef="7"' in text, text
    again, = read_mrv(text, log=[])
    assert next(iter(again.atoms())).r_index == 7


def test_an_index_past_the_domain_leaves_the_marker_unindexed():
    log = []
    mol, = read_mrv(_wrap(_atoms('<atom id="a1" elementType="R100"/>')), log=log)
    assert next(iter(mol.atoms())).is_r
    assert next(iter(mol.atoms())).r_index == 0
    assert any(x.rule == 'xml:r-index-too-wide' for x in log), log


def test_rubidium_still_resolves():
    """The regression guard for the R family moving above the element symbols: ``Rb`` is an element, and
    an inexact prefix test would refuse it along with eight others."""
    mol, = read_mrv(_wrap(_atoms('<atom id="a1" elementType="Rb"/>')), log=[])
    assert mol.atom(next(iter(mol.atom_numbers))).atomic_symbol == 'Rb'


def test_a_lowercase_marker_is_folded_like_a_lowercase_symbol():
    """``r1`` is the R family upper-cased, the same recovery ``cl`` gets and for the same reason -- a
    document converted out of a molfile inherits the molfile's case-folding -- and it is reported."""
    log = []
    mol, = read_mrv(_wrap(_atoms('<atom id="a1" elementType="r1"/>')), log=log)
    atom = next(iter(mol.atoms()))
    assert atom.is_r and atom.r_index == 1
    assert any(x.rule == 'xml:element-folded' for x in log), log


def test_the_two_recoveries_an_element_type_gets_are_reached_from_this_dialect():
    """``D``/``T`` as hydrogen isotopes, and an upper-cased symbol folded back to its element.

    The fixture is ``CL`` because an XML document converted out of a molfile inherits the molfile's
    upper-cased symbols.  Both are recovered *and* reported: the atom that comes out is not the atom the
    text names.  Per dialect rather than on the shared resolver, since what it proves is the routing.
    """
    log = []
    mol, = read_mrv(_wrap(_atoms('<atom id="a1" elementType="D"/>',
                                 '<atom id="a2" elementType="T"/>',
                                 '<atom id="a3" elementType="CL"/>')), log=log)
    d, t, cl = mol.atom_numbers
    assert (mol.element_of(d), mol.isotope_of(d)) == (1, 2)
    assert (mol.element_of(t), mol.isotope_of(t)) == (1, 3)
    assert (mol.element_of(cl), mol.isotope_of(cl)) == (17, 0)
    assert [str(x) for x in log] == ['atom a1: D read as hydrogen isotope 2',
                                     'atom a2: T read as hydrogen isotope 3',
                                     "atom a3: elementType 'CL' read as 'Cl'"], log


def test_no_element_type_is_read_as_carbon_and_said_once_for_the_record():
    """MRV requires an ``elementType``; an ``<atom>`` without one is broken and carbon is all there is.

    Said once -- a count and the first atom -- rather than a line per atom, which on a file whose whole
    array lost the column is one line per atom.  The oxygen is the armed control: the count is 2, not 3.
    """
    log = []
    mol, = read_mrv(_wrap(_atoms('<atom id="a1"/>', '<atom id="a2" elementType="O"/>',
                                 '<atom id="a3"/>')), log=log)
    assert [mol.element_of(s) for s in mol.atom_numbers] == [6, 8, 6]
    assert [str(x) for x in log] == ['atom: 2 atom(s) with no elementType, read as carbon (first atom a1)'], log


def test_a_malformed_number_costs_the_attribute_and_not_the_atom_or_the_bond():
    """A negative ``mrvValence`` and a ``convention`` naming neither of MRV's two, one per side.

    Both are codec refusals the engine turns into an unprefixed line naming the attribute and its value.
    The atom and bond survive, and the bond keeps the ``order`` stated beside the bad ``convention``.
    """
    log = []
    record, = parse_mrv(_wrap(_atoms('<atom id="a1" elementType="C" mrvValence="-2"/>')), log=log)
    assert record.ctab.atoms[0].valence is None
    assert [x for x in log if 'mrvValence' in x and 'negative' in x] and not _unsupported(log), log

    log = []
    record, = parse_mrv(_wrap(_atoms('<atom id="a1" elementType="C"/>', '<atom id="a2" elementType="C"/>')
                              + '<bondArray><bond id="b1" atomRefs2="a1 a2" order="1" '
                                'convention="cxn:nonsense"/></bondArray>'), log=log)
    assert record.ctab.bonds[0].order == 1
    assert [x for x in log if 'cxn:coord' in x] and not _unsupported(log), log


def test_an_unmodelled_atom_attribute_is_named_by_the_engine():
    """The acceptance rule, which is the engine's property and not this dialect's diligence.

    ``mrvQueryProps`` is a real MRV construct with no row here, so it earns a prefixed line -- and the
    atom is still read.  Reported by decision rather than for want of a home: MRV does not read into a
    ``QueryContainer``, so there is nowhere for a query primitive to land.
    """
    log = []
    mol, = read_mrv(_wrap(_atoms(
        '<atom id="a1" elementType="C" mrvQueryProps="A"/>')), log=log)
    assert mol.atom_count == 1
    assert 'mrvQueryProps' in ' '.join(str(x) for x in _unsupported(log)), log


def test_a_modelled_atom_attribute_is_not_named_by_the_engine():
    """The other side of the rule above, and the half that decays silently without a test: a row that is
    *read* must stop being reported, or the log tells a caller a construct was lost that was not.  The
    atom is read either way, so nothing but the log distinguishes the two states."""
    log = []
    mol, = read_mrv(_wrap(_atoms(
        '<atom id="a1" elementType="C" mrvStereoGroup="and1"/>')), log=log)
    assert mol.atom_count == 1
    assert log == [], log


def test_the_gui_selection_flag_is_ignored_rather_than_reported():
    """``isSelected`` is the Marvin editor's selection state: the same molecule either way, which is the
    claim an ``*_ignored`` entry makes.  The sibling above proves the log can speak about an attribute."""
    log = []
    mol, = read_mrv(_wrap(_atoms('<atom id="a1" elementType="C" isSelected="true"/>')), log=log)
    assert mol.atom_count == 1
    assert log == [], log


# mrvAlias.  `MoleculeContainer.aliases` documents itself as "V2000 `A  <n>` lines, MRV mrvAlias", so
# these tests compare the two paths against each other.

#: A two-atom molfile with an ``A  <n>`` alias, for the parity test below.  The alias text is on the
#: line *after* the header, which is the whole shape of the V2000 construct.
_V2000_WITH_ALIAS = '\n'.join((
    'ethanol-ish',
    '  chython',
    '',
    '  2  1  0  0  0  0  0  0  0  0999 V2000',
    '    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0',
    '    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
    '  1  2  1  0  0  0  0',
    'A    1',
    'OMe',
    'M  END',
))


def test_an_atom_alias_is_read_onto_the_molecule_the_v2000_way():
    """``mrvAlias`` reaches ``mol.aliases`` keyed by the atom's stable id, and earns no log line.

    Both halves: an empty log is satisfied by an attribute silently dropped, and the alias being present
    says nothing about whether the reader also reported it as a loss it did not take.
    """
    log = []
    mol, = read_mrv(_wrap(_atoms('<atom id="a1" elementType="O" mrvAlias="OMe"/>',
                                 '<atom id="a2" elementType="C"/>')), log=log)
    sid = next(iter(mol.atom_numbers))
    assert mol.aliases == {sid: b'OMe'}
    assert log == [], log


def test_the_alias_is_the_same_value_the_v2000_reader_produces_for_the_same_label():
    """One label, two formats, one answer, measured against the V2000 path.

    ``bytes``, which codec, and the key being a stable id rather than a file index are all decisions
    ``_v2000._parse_properties`` and ``SGroupStore`` already made; a second convention invented here
    would be invisible to every assertion that only looked at MRV.
    """
    log = []
    from_mrv, = read_mrv(_wrap(_atoms('<atom id="a1" elementType="O" mrvAlias="OMe"/>',
                                      '<atom id="a2" elementType="C"/>')), log=log)
    from_v2000 = mol(_V2000_WITH_ALIAS)
    assert list(from_v2000.aliases.values()) == [b'OMe'], from_v2000.aliases
    assert dict(from_mrv.aliases) == dict(from_v2000.aliases)
    assert log == [], log


def test_an_empty_alias_is_stored_rather_than_dropped_because_v2000_stores_one():
    """``mrvAlias=""`` is a stated empty label, and the V2000 path stores its blank line the same way.

    Dropping a falsy alias is a codec deciding a statement was silence, which is the distinction
    :func:`~.._dialect.encode_stated` keeps.  Round-tripped too, or the storage is a dead end.
    """
    log = []
    mol, = read_mrv(_wrap(_atoms('<atom id="a1" elementType="O" mrvAlias=""/>')), log=log)
    assert mol.aliases == {next(iter(mol.atom_numbers)): b''}
    assert log == [], log
    assert 'mrvAlias=""' in write_mrv(mol, indent=None)


def test_an_alias_survives_a_write_and_is_not_counted_as_an_s_group():
    """The write half, and the second assertion is the point: an alias reaching the S-group loss line
    would be a writer reporting a loss it did not take."""
    log = []
    mol = read_smiles('CO')
    sid = sorted(mol.atom_numbers)[-1]
    mol.set_aliases({sid: 'OMe'})
    text = write_mrv(mol, log=log, indent=None)
    assert 'mrvAlias="OMe"' in text, text
    assert log == [], log


def test_an_alias_round_trips_through_the_document_and_back_onto_its_own_atom():
    """Write, read, and the label is on the same atom -- not merely somewhere in the molecule.

    An assertion on the *set* of alias texts passes for any permutation, and a permutation is the defect
    a position-keyed intermediate can introduce, so the atom is named by its element.
    """
    mol = read_smiles('CO')
    # `element_of` answers an atomic number; 8 is oxygen.
    oxygen = next(n for n in mol.atom_numbers if mol.element_of(n) == 8)
    mol.set_aliases({oxygen: 'OMe'})

    log = []
    again, = read_mrv(write_mrv(mol, indent=None), log=log)
    assert log == [], log
    labelled, = again.aliases
    assert again.element_of(labelled) == 8
    assert again.aliases[labelled] == b'OMe'


def test_the_zero_placeholder_is_no_alias_in_either_form():
    """``mrvAlias="0"`` is MRV's "no alias here", the placeholder ``mrvStereoGroup="0 and1 0"`` uses.

    The column form needs one -- every cell of a column is filled -- and MEASURED at Marvin 25.1.3, the
    element form reads it the same way: ``<atom mrvAlias="0"/>`` comes back out of ``molconvert -g mrv``
    with no ``mrvAlias`` at all, byte-identical to the atom that stated none.  So the placeholder is the
    attribute's and not the column's, and one codec answers for both forms.

    Bounded by the two neighbouring statements, which are not the placeholder: ``mrvAlias=""`` is a stated
    empty label (see above), and ``zero``/``.`` -- ChemAxon's escape spellings for the string ``0`` and for
    empty -- come back out of the same ``molconvert`` verbatim, so they are read as the text they are.
    """
    log = []
    mol, = read_mrv(_wrap('<atomArray atomID="a1 a2 a3" elementType="C C O" mrvAlias="Me 0 0" '
                          'x2="0 -0.7 -1.4" y2="0 -0.4 0"/>'), log=log)
    labelled, = mol.aliases
    assert mol.aliases[labelled] == b'Me'
    assert mol.element_of(labelled) == 6
    assert log == [], log

    log = []
    mol, = read_mrv(_wrap(_atoms('<atom id="a1" elementType="O" mrvAlias="0"/>',
                                 '<atom id="a2" elementType="C" mrvAlias="zero"/>',
                                 '<atom id="a3" elementType="C" mrvAlias="."/>')), log=log)
    assert sorted(mol.aliases.values()) == [b'.', b'zero'], mol.aliases
    assert log == [], log


def test_an_alias_whose_text_is_the_placeholder_cannot_be_written_and_says_so():
    """The one label this dialect cannot state, reported rather than written in silence.

    ``0`` is the placeholder in both forms, so an atom whose alias text *is* ``0`` has no MRV spelling --
    Marvin 25.1.3 passes ``zero`` through verbatim rather than reading it as the string ``0``, so that
    escape is no channel either.  The attribute is still written, the document then carrying the text for a
    reader that takes it literally, and the line is what keeps it from being a silent loss.
    """
    log = []
    mol = read_smiles('CO')
    mol.set_aliases({sorted(mol.atom_numbers)[-1]: '0'})
    text = write_mrv(mol, log=log, indent=None)
    assert 'mrvAlias="0"' in text, text
    assert [x for x in _unsupported(log) if 'mrvAlias' in x], log


def test_a_record_read_from_a_document_writes_its_aliases_back():
    """The ``Record`` path, which a caller uses to keep the file's own atom ids.

    ``_finish`` leaves the decoded alias in the atom's spill *and* fills ``Ctab.aliases``, so a record
    that never becomes a molecule still writes it back; one that popped the spill loses it exactly here.
    """
    log = []
    record, = parse_mrv(_wrap(_atoms('<atom id="a7" elementType="O" mrvAlias="OMe"/>')), log=log)
    assert record.ctab.aliases == {0: 'OMe'}
    text = write_mrv(record, log=log, indent=None)
    assert 'id="a7"' in text and 'mrvAlias="OMe"' in text, text
    assert log == [], log


# hydrogenCount

def test_hydrogen_count_is_the_implicit_count_and_not_the_total(data):
    """MRV's ``hydrogenCount`` is the count of hydrogens **not** drawn, unlike CML's total.

    ``test/mrv_hydrogens.mrv`` is the discriminating case: a nitrogen carrying ``hydrogenCount="1"``
    *and* one drawn ``<atom elementType="H">``.  Read as implicit, which is MRV's meaning, the total is
    two -- methylamine; read as a total, which is CML's, there is no implicit hydrogen and it is an
    aminyl.  Hence the assertion on the total, stated as a number.
    """
    log = []
    mol, = read_mrv(data('mrv_hydrogens.mrv'), log=log)
    nitrogen, = [s for s in mol.atom_numbers if mol.element_of(s) == 7]
    assert mol.implicit_h_of(nitrogen) == 1
    assert mol.total_h_of(nitrogen) == 2
    # The drawn hydrogen is an atom of the graph, as in the file, so the comparison keeps it -- and it
    # is a structure comparison, never a SMILES one.
    assert mol == read_smiles('[H]NC')
    assert log == [], log


def test_a_hydrogen_count_outside_what_an_atom_can_hold_costs_the_attribute():
    """A malformed count is a fact about the file, so the atom is kept and the line is unprefixed.

    Both ends, plus a legal count as the control -- a reader dropping every ``hydrogenCount`` would pass
    the first half alone.
    """
    for value in ('-1', '99'):
        log = []
        mol, = read_mrv(_wrap(_atoms(f'<atom id="a1" elementType="C" hydrogenCount="{value}"/>')),
                        log=log)
        assert mol.atom_count == 1
        assert log and not _unsupported(log), log
    log = []
    mol, = read_mrv(_wrap(_atoms('<atom id="a1" elementType="C" hydrogenCount="4"/>')), log=log)
    assert mol.implicit_h_of(next(iter(mol.atom_numbers))) == 4
    assert log == [], log


def test_a_stated_valence_reaches_the_one_hydrogen_derivation():
    """``mrvValence`` is a stated total valence, the third rank of the shared derivation.

    Read onto ``CtabAtom.valence``, the field V2000's ``vvv`` and V3000's ``VAL=`` land on, so the count
    is the CTfile reader's answer and not a second one.  ``_hydrogens.py`` ranks a stated valence *below*
    the derivation, so an atom with bonds drawn is not discriminating -- an atom with **no** bonds drawn
    is the one place the statement wins, so a bare carbon reads 3 with the attribute and 4 without it.
    """
    log = []
    mol, = read_mrv(_wrap(_atoms('<atom id="a1" elementType="C" mrvValence="3"/>')), log=log)
    assert mol.implicit_h_of(next(iter(mol.atom_numbers))) == 3
    assert [x for x in log if 'stated valence 3' in x and 'no bonds drawn' in x], log
    assert not _unsupported(log), log

    log = []
    mol, = read_mrv(_wrap(_atoms('<atom id="a1" elementType="C"/>')), log=log)
    assert mol.implicit_h_of(next(iter(mol.atom_numbers))) == 4
    assert log == [], log

    record, = parse_mrv(_wrap(_atoms('<atom id="a1" elementType="C" mrvValence="3"/>')), log=[])
    assert record.ctab.atoms[0].valence == 3


# the bond vocabulary

def test_the_four_bond_orders_are_read_and_aromatic_is_stored_as_stated():
    """``1``, ``2``, ``3`` and ``A``, and ``A`` is stored as aromatic rather than kekulised -- a reader
    that kekulised here would be repairing on a read path."""
    log = []
    record, = parse_mrv(_wrap(_atoms(*(f'<atom id="a{i}" elementType="C"/>' for i in range(1, 6)))
                              + '<bondArray>'
                                '<bond id="b1" atomRefs2="a1 a2" order="1"/>'
                                '<bond id="b2" atomRefs2="a2 a3" order="2"/>'
                                '<bond id="b3" atomRefs2="a3 a4" order="3"/>'
                                '<bond id="b4" atomRefs2="a4 a5" order="A"/>'
                                '</bondArray>'), log=log)
    assert [b.order for b in record.ctab.bonds] == [1, 2, 3, 4]
    assert log == [], log


def test_a_bond_order_the_format_does_not_have_costs_the_order_and_not_the_bond():
    """``order="9"`` is malformed, so the bond survives as a single and the line is unprefixed."""
    log = []
    record, = parse_mrv(_wrap(_atoms('<atom id="a1" elementType="C"/>',
                                     '<atom id="a2" elementType="C"/>')
                              + '<bondArray><bond id="b1" atomRefs2="a1 a2" order="9"/></bondArray>'),
                        log=log)
    assert len(record.ctab.bonds) == 1
    assert log and not _unsupported(log), log


def test_a_coordination_bond_is_the_convention_attribute_and_not_an_order():
    """MRV writes a dative bond as ``convention="cxn:coord"`` with **no** ``order``.

    The convention outranks an order that is present, so it is parked and resolved rather than applied
    where it is read -- the attributes arrive in whatever sequence the file wrote them, and a row
    assigning straight onto the order would give this bond order 1 half the time.
    """
    for bond in ('<bond id="b1" atomRefs2="a1 a2" convention="cxn:coord"/>',
                 '<bond id="b1" atomRefs2="a1 a2" convention="cxn:coord" order="1"/>',
                 '<bond id="b1" atomRefs2="a1 a2" order="1" convention="cxn:coord"/>'):
        log = []
        record, = parse_mrv(_wrap(_atoms('<atom id="a1" elementType="N"/>',
                                         '<atom id="a2" elementType="B"/>')
                                  + f'<bondArray>{bond}</bondArray>'), log=log)
        assert [b.order for b in record.ctab.bonds] == [8], bond
        assert log == [], log


def test_a_hydrogen_bond_is_dropped_rather_than_read_as_a_single_bond():
    """``convention="cxn:hydrogen"`` has no order in a molecule, and single is the wrong guess.

    Kept as the single its missing ``order`` defaults to, it joins two molecules the file drew apart -- a
    different compound, not a damaged one.  The coordination case above is the control.
    """
    log = []
    record, = parse_mrv(_wrap(_atoms('<atom id="a1" elementType="O"/>',
                                     '<atom id="a2" elementType="O"/>')
                              + '<bondArray><bond id="b1" atomRefs2="a1 a2" '
                                'convention="cxn:hydrogen"/></bondArray>'), log=log)
    assert len(record.ctab.atoms) == 2 and record.ctab.bonds == []
    assert len(record.bond_extras) == 0
    assert [x for x in _unsupported(log) if 'hydrogen bond' in x], log


def test_a_query_bond_type_is_named_rather_than_guessed_at():
    """``queryType="SD"`` is a query feature with no row, so the engine names it: single-or-double is not
    a bond order and picking either invents chemistry.  The bond is still read as a single."""
    log = []
    record, = parse_mrv(_wrap(_atoms('<atom id="a1" elementType="C"/>',
                                     '<atom id="a2" elementType="C"/>')
                              + '<bondArray><bond id="b1" atomRefs2="a1 a2" '
                                'queryType="SD"/></bondArray>'), log=log)
    assert len(record.ctab.bonds) == 1
    assert [x for x in _unsupported(log) if 'queryType' in x], log


# the column form of <atomArray>

def test_the_column_form_of_an_atom_array_is_read():
    """One attribute per column, ``atomID`` naming the atoms, which real Marvin files write.

    The failure this guards is silent: with no hook the engine walks children a column form does not
    have, so the molecule comes back with no atoms, every bond references one that was never declared,
    and the log holds nothing but dropped bonds.  Compared as structures against the element form.
    """
    log = []
    columns, = read_mrv(_wrap('<atomArray atomID="a1 a2 a3" elementType="C C O" '
                              'formalCharge="0 0 -1" x2="0.0 1.3 2.6" y2="0.0 0.75 0.0"/>'
                              '<bondArray>'
                              '<bond id="b1" atomRefs2="a1 a2" order="1"/>'
                              '<bond id="b2" atomRefs2="a2 a3" order="1"/>'
                              '</bondArray>'), log=log)
    assert columns.atom_count == 3
    assert log == [], log
    elements, = read_mrv(_wrap(_atoms('<atom id="a1" elementType="C" x2="0.0" y2="0.0"/>',
                                      '<atom id="a2" elementType="C" x2="1.3" y2="0.75"/>',
                                      '<atom id="a3" elementType="O" formalCharge="-1" x2="2.6" '
                                      'y2="0.0"/>')
                               + '<bondArray>'
                                 '<bond id="b1" atomRefs2="a1 a2" order="1"/>'
                                 '<bond id="b2" atomRefs2="a2 a3" order="1"/>'
                                 '</bondArray>'), log=[])
    assert columns == elements


def test_a_column_of_nulls_states_nothing_and_is_not_reported():
    """``-`` is the column form's null, and ``0`` is the null the ``radical`` column uses.

    Such a cell states nothing, so it must neither set a field nor earn a line.  Armed by the third atom,
    which states a real value in the same columns.
    """
    log = []
    mol, = read_mrv(_wrap('<atomArray atomID="a1 a2 a3" elementType="C C C" '
                          'hydrogenCount="- - 1" mrvValence="- - 4" radical="0 0 monovalent" '
                          'lonePair="0 0 0"/>'), log=log)
    first, second, third = mol.atom_numbers
    assert mol.implicit_h_of(first) == 4 and mol.implicit_h_of(second) == 4
    assert not mol.radical_of(first)
    assert mol.implicit_h_of(third) == 1 and mol.radical_of(third)
    assert log == [], log


def test_a_zero_in_a_column_with_a_row_is_a_stated_zero():
    """The other side of the null rule: ``hydrogenCount="0"`` is an absence of hydrogens, not a silence,
    so only a column with **no** row treats ``0`` as nothing said.  A carbon stating zero hydrogens is a
    carbene, and reading it as unstated would quietly give it four."""
    log = []
    mol, = read_mrv(_wrap('<atomArray atomID="a1" elementType="C" hydrogenCount="0"/>'), log=log)
    assert mol.implicit_h_of(next(iter(mol.atom_numbers))) == 0
    assert log == [], log


def test_columns_beside_atom_elements_are_the_document_contradicting_itself():
    """The children win, being the more specific statement, and the dropped columns are named.

    The form is decided by the children and never by the attributes: the vocabulary a global XML
    attribute can come from is open, so a skip list would let the next such name claim an element-form
    array and lose every atom in it.
    """
    log = []
    mol, = read_mrv(_wrap('<atomArray elementType="N N"><atom id="a1" elementType="O"/></atomArray>'),
                    log=log)
    assert mol.atom_count == 1
    assert [x for x in log if 'column' in x and 'dropped' in x], log


def test_ragged_columns_are_truncated_rather_than_dropped():
    """Losing every atom in a molecule over one short column is the opposite of the input posture."""
    log = []
    mol, = read_mrv(_wrap('<atomArray atomID="a1 a2 a3" elementType="C C"/>'), log=log)
    assert mol.atom_count == 2
    assert [x for x in log if 'different lengths' in x], log


def test_an_array_wide_convention_is_named_and_a_title_is_not():
    """The acceptance rule on the *array* element: an attribute is applied or it is named.

    ``title`` is a label, so a reader honouring it builds the same atoms and it is silent.  ``convention``
    names the dictionary the columns beside it are defined in, and a foreign dictionary can redefine every
    one of them -- the columns are read under MRV's own meanings regardless.  Both arrays take the
    attribute; the ``title``-only document is the armed control.
    """
    log = []
    mol, = read_mrv(_wrap('<atomArray title="tt" convention="cxn:whatever" atomID="a1 a2 a3" '
                          'elementType="C C O" formalCharge="0 0 -1" x2="0 1 2" y2="0 0 0"/>'), log=log)
    assert mol.atom_count == 3
    assert [str(x) for x in log] == ["unsupported: atom: <atomArray> convention='cxn:whatever' names a dictionary this "
                                     "reader has not read; the columns are read under MRV's own meanings"], log

    log = []
    record, = parse_mrv(_wrap(_atoms('<atom id="a1" elementType="C"/>', '<atom id="a2" elementType="C"/>')
                              + '<bondArray convention="foo" title="t"/>'), log=log)
    assert len(record.ctab.atoms) == 2
    assert [str(x) for x in log] == [
        "unsupported: bond: <bondArray> convention='foo' names a dictionary this reader has "
        "not read; the columns are read under MRV's own meanings"], log

    log = []
    mol, = read_mrv(_wrap('<atomArray title="tt" atomID="a1 a2" elementType="C O"/>'), log=log)
    assert mol.atom_count == 2
    assert log == [], log


def test_the_column_form_of_a_bond_array_is_named_rather_than_guessed_at():
    """No source describes one, so there is no vocabulary to write a hook against.

    Inventing the endpoint columns from the atom array's shape is the guess this table avoids, and the
    bonds are then lost, so the line stops them being lost silently.  An empty ``<bondArray/>`` is the
    control: no columns, so nothing stated and nothing said.
    """
    log = []
    record, = parse_mrv(_wrap(_atoms('<atom id="a1" elementType="C"/>',
                                     '<atom id="a2" elementType="C"/>')
                              + '<bondArray atomRefs2="a1 a2" order="1"/>'), log=log)
    assert record.ctab.bonds == []
    assert [x for x in _unsupported(log) if 'column form' in x], log

    log = []
    record, = parse_mrv(_wrap(_atoms('<atom id="a1" elementType="C"/>') + '<bondArray/>'), log=log)
    assert record.ctab.bonds == []
    assert log == [], log


# <bondStereo>

def _wedge_record(child, log):
    body = _atoms('<atom id="a1" elementType="C" x2="0.0" y2="0.0"/>',
                  '<atom id="a2" elementType="F" x2="1.3" y2="0.75"/>')
    record, = parse_mrv(_wrap(body + f'<bondArray><bond id="b1" atomRefs2="a1 a2" order="1">{child}'
                                     f'</bond></bondArray>'), log=log)
    return record


def test_the_three_spellings_of_a_wedge_all_reach_one_field():
    """A bare letter, a ``dictRef`` and the MDL dictionary reference are one construct three ways.

    All three are MRV's own, each a closed set in both sources, and all three land on ``CtabBond.wedge``.
    The MDL numbers decode through the CTfile reader's own table rather than a copy of it.
    """
    for child, wedge in (('<bondStereo>W</bondStereo>', WEDGE_UP),
                         ('<bondStereo>H</bondStereo>', WEDGE_DOWN),
                         ('<bondStereo>w</bondStereo>', WEDGE_UP),
                         ('<bondStereo dictRef="cml:W"/>', WEDGE_UP),
                         ('<bondStereo dictRef="cml:H"/>', WEDGE_DOWN),
                         ('<bondStereo convention="MDL" conventionValue="1"/>', WEDGE_UP),
                         ('<bondStereo convention="MDL" conventionValue="6"/>', WEDGE_DOWN)):
        log = []
        record = _wedge_record(child, log)
        assert record.ctab.bonds[0].wedge == wedge, child
        assert log == [], (child, log)


def test_a_cis_trans_bond_stereo_is_read_as_a_configuration_and_not_as_a_wedge():
    """``C`` and ``T`` state a double-bond configuration, which is a different field from a wedge.

    Both of MRV's spellings, folded into one letter before either is looked at.  The frame is **empty** --
    no source describes an ``atomRefs4`` on an MRV ``<bondStereo>`` -- and storing it that way rather than
    guessing is what lets :func:`chython.core.wedge.stated_cis_trans` refuse the letter where it would be
    ambiguous.  In a coordinate-free record the letter is the document's only statement.
    """
    for child in ('<bondStereo>C</bondStereo>', '<bondStereo dictRef="cml:C"/>'):
        log = []
        record = _wedge_record(child, log)
        assert record.ctab.bonds[0].wedge == 0
        assert record.ctab.bonds[0].configuration == ('C', ()), child
        assert log == [], log


def test_a_bare_cis_trans_letter_with_no_drawing_carries_the_configuration():
    """The case the letter is worth reading for, end to end through this dialect.

    A coordinate-free MRV record states its configuration in the letter and nowhere else, and MRV writes
    it bare -- so the trip closes only if a bare letter is read where it cannot be ambiguous, one
    substituent per terminal, which 2-butene is.  Pinned against the SMILES reader, the tree's other
    coordinate-free source and so the only independent check of the sign.
    """
    body = _atoms(*(f'<atom id="a{i}" elementType="C"/>' for i in range(1, 5)))
    doc = _wrap(body + '<bondArray><bond id="b1" atomRefs2="a1 a2" order="1"/>'
                       '<bond id="b2" atomRefs2="a2 a3" order="2">'
                       '<bondStereo>{0}</bondStereo></bond>'
                       '<bond id="b3" atomRefs2="a3 a4" order="1"/></bondArray>')
    cis_log, trans_log = [], []
    cis, = read_mrv(doc.format('C'), log=cis_log)
    trans, = read_mrv(doc.format('T'), log=trans_log)
    assert cis_log == [], cis_log
    assert trans_log == [], trans_log
    assert write_smiles(trans) == write_smiles(read_smiles('C/C=C/C'))
    assert write_smiles(cis) == write_smiles(read_smiles('C/C=C\\C'))


def test_a_bond_stereo_in_a_dictionary_this_reader_has_not_read_applies_nothing():
    """A ``convention`` names the dictionary its value is defined in, and MDL's is the only one read.

    ``W`` could mean anything in somebody else's, so the check runs before the MDL branch -- a
    ``conventionValue`` under a foreign dictionary is not decoded as a CTfile code either.
    """
    log = []
    record = _wedge_record('<bondStereo convention="ACME" conventionValue="1"/>', log)
    assert record.ctab.bonds[0].wedge == 0
    assert [x for x in _unsupported(log) if 'has not read' in x], log


def test_a_malformed_bond_stereo_value_costs_the_descriptor_and_not_the_bond():
    """An unnumbered ``conventionValue``, a number outside MDL's set, and an empty element.

    All three are the file being broken rather than a feature we lack, so all three are unprefixed and
    the bond survives.  ``3`` is the one value in MDL's own set with no field here, so it is the armed
    opposite: prefixed, because the construct is real.
    """
    for child in ('<bondStereo convention="MDL" conventionValue="x"/>',
                  '<bondStereo convention="MDL" conventionValue="7"/>',
                  '<bondStereo></bondStereo>'):
        log = []
        record = _wedge_record(child, log)
        assert len(record.ctab.bonds) == 1 and record.ctab.bonds[0].wedge == 0
        assert log and not _unsupported(log), (child, log)

    log = []
    record = _wedge_record('<bondStereo convention="MDL" conventionValue="3"/>', log)
    assert [x for x in _unsupported(log) if 'unknown which' in x], log


def test_a_bond_stereo_letter_outside_mrv_s_own_set_applies_nothing():
    """The bare-letter branch's else: ``W``, ``H``, ``C`` and ``T`` are the whole vocabulary, so anything
    else applies nothing and the bond survives.

    The line is **bare**, and the prefix's absence is the point: ``unsupported: `` promises the file was
    right and this reader is the limitation, which a caller can act on by converting the file or waiting.
    ``<bondStereo>`` text is a closed set, so a fifth letter is the file being broken and there is no
    feature to wait for.  The MDL ``3`` case two tests up keeps the prefix, being a real construct.
    """
    log = []
    record = _wedge_record('<bondStereo>Z</bondStereo>', log)
    assert len(record.ctab.bonds) == 1 and record.ctab.bonds[0].wedge == 0
    assert [str(x) for x in log] == ["bond b1: bondStereo 'Z' is not one of W, H, C or T, dropped"], log


def test_a_bond_child_no_handler_claims_is_named():
    """The acceptance rule one level in: a child element of ``<bond>`` that is not a ``<bondStereo>``."""
    log = []
    _wedge_record('<mrvQueryProps queryString="x"/>', log)
    assert [x for x in _unsupported(log) if 'mrvQueryProps' in x], log


# the document level, and coordinates

def test_a_reaction_s_roles_are_counted_in_the_log(data):
    """The molecules are read flat and the roles are reported, which is CML's treatment generalised.

    The walker reads straight through a ``<reaction>``, so the record comes back indistinguishable from
    the same molecules loose -- and the counts are what let a caller holding only the log recover the
    record's shape, which is why the line names the roles rather than saying a reaction was here.
    """
    log = []
    molecules = read_mrv(data('mrv_reaction.mrv'), log=log)
    assert len(molecules) == 3
    roles, = [x for x in _unsupported(log) if '<reaction>' in x]
    assert '1 reactant' in roles and '1 agent' in roles and '1 product' in roles, roles

    # A `<reaction>` holding no molecule still says so, and says something other than a count of
    # nothing: an empty role list is a shape a caller holding only the log cannot recover.
    log = []
    assert read_mrv(_wrap_document('<reaction id="r1"><reactantList/></reaction>'), log=log) == []
    assert [str(x) for x in log] == ['unsupported: record: <reaction> r1 roles are not modelled; it holds no '
                                     'molecules'], log


def test_the_furniture_marvin_draws_around_a_structure_is_named_once(data):
    """An arrow and a reaction sign are real constructs with nothing here to hold them.

    Once per document with a count, not once per element: a line per text box says the same thing as many
    times as the drawing is elaborate.  The container test above is the arm.
    """
    log = []
    read_mrv(data('mrv_reaction.mrv'), log=log)
    furniture, = [x for x in _unsupported(log) if 'furniture' in x]
    assert '<arrow>' in furniture and '<MReactionSign>' in furniture, furniture


def test_a_conformer_is_read_from_the_three_dimensional_attributes():
    """``x3``/``y3``/``z3`` is a conformer, and a record carrying only that is three-dimensional."""
    log = []
    record, = parse_mrv(_wrap(_atoms('<atom id="a1" elementType="C" x3="1.0" y3="2.0" z3="3.0"/>')),
                        log=log)
    assert record.ctab.dimensionality == '3D'
    assert (record.ctab.atoms[0].x, record.ctab.atoms[0].y, record.ctab.atoms[0].z) == (1.0, 2.0, 3.0)
    assert log == [], log


def test_a_drawing_and_a_conformer_together_keep_the_drawing_and_say_so():
    """Both sets is a drawing *and* a conformer, and every stereo statement is measured against the
    drawing -- mixing them gives a geometry that exists in no file and reads stereo out of it."""
    log = []
    record, = parse_mrv(_wrap(_atoms('<atom id="a1" elementType="C" x2="1.0" y2="2.0" '
                                     'x3="4.0" y3="5.0" z3="6.0"/>')), log=log)
    assert record.ctab.dimensionality == '2D'
    assert (record.ctab.atoms[0].x, record.ctab.atoms[0].y) == (1.0, 2.0)
    assert [x for x in _unsupported(log) if 'both 2D and 3D' in x], log


# the write side

def test_the_wrapper_is_the_one_marvin_writes_and_claims_no_version():
    """``<cml><MDocument><MChemicalStruct>``, a ``molID`` per molecule, and no provenance claim.

    Literal text rather than a read-back: this reader walks through containers on purpose, so it reads a
    document with none of them just as happily.  ``version=`` is checked *absent* -- real files put
    ``version="ChemAxon file format v18.11.0, generated by v19.7.0"`` there, naming the writing program,
    and a downstream reader is entitled to act on that string.
    """
    text = write_mrv([read_smiles('O'), read_smiles('N')])
    assert f'<cml xmlns="{MRV_NS}">' in text
    assert '<MDocument>' in text and '<MChemicalStruct>' in text
    assert 'molID="m1"' in text and 'molID="m2"' in text
    # The XML declaration's own `version` is the document's, so the root element is checked and not the
    # whole text.
    assert text.startswith('<?xml version="1.0" encoding="UTF-8"?>\n<cml ')
    root_tag = text.split('>', 2)[1]
    assert 'version=' not in root_tag, root_tag


def test_every_atom_quantity_goes_out_in_mrv_s_own_spelling():
    """The exact attribute names, against a molecule stating all five.

    ``mrvMap`` and not ``id``, ``radical`` as a name and not a count.  A round trip passes with any
    spelling as long as both halves agree on it, so the strings are pinned here.
    """
    text = write_mrv(read_smiles('[13CH3][NH3+]'), indent=None)
    assert 'elementType="C" isotope="13" hydrogenCount="3"' in text
    assert 'elementType="N" formalCharge="1" hydrogenCount="3"' in text

    text = write_mrv(read_smiles('[CH3:5][O-:2]'), indent=None)
    assert 'mrvMap="5"' in text and 'formalCharge="-1" mrvMap="2"' in text

    # `|^1:0|` and not `[CH3]`: chython does not infer a radical from a short valence, so the flag must
    # be stated for the attribute to have anything to write.
    text = write_mrv(read_smiles('[CH3]C |^1:0|'), indent=None)
    assert 'elementType="C" radical="monovalent" hydrogenCount="3"' in text
    assert text.count('radical=') == 1  # and not on the methyl carbon next to it


def test_a_stated_zero_and_an_unknown_hydrogen_count_are_written_differently():
    """``hydrogenCount="0"`` is a statement, so an unknown count must not borrow it.

    Two halves of one decision: the anion states 0 and it goes out, and the aromatic nitrogen whose class
    only the ring decides gets no attribute plus a line.  Defaulting the unknown to 0 passes the first.
    """
    assert 'hydrogenCount="0"' in write_mrv(read_smiles('C[O-]'), indent=None)

    log = []
    # An aromatic pyrrole-or-pyridine nitrogen: the one atom class this tree leaves `H_UNKNOWN`.
    mol, = read_mrv(_wrap(_atoms(*(f'<atom id="a{i}" elementType="{e}"/>'
                                   for i, e in enumerate('CCCCN', 1)))
                          + '<bondArray>'
                          + ''.join(f'<bond id="b{i}" atomRefs2="a{a} a{b}" order="A"/>'
                                    for i, (a, b) in enumerate(((1, 2), (2, 3), (3, 4), (4, 5), (5, 1)), 1))
                          + '</bondArray>'))
    text = write_mrv(mol, log=log, indent=None)
    assert text.count('hydrogenCount=') == 4  # the four carbons, and not the nitrogen
    assert [x for x in log if 'hydrogen count unknown' in x], log


def test_the_bond_orders_and_the_coordination_convention_go_out_as_read():
    """``order="A"`` for aromatic, and a dative bond as ``convention`` with no ``order`` beside it.

    Aromatic is a letter and not a number here, the one order whose spelling cannot be guessed from the
    internal value.  The coordination bond writes a *different attribute* than it reads, so the absence of
    an ``order`` is asserted too -- a bond carrying both reads back by the precedence rule and hides it.
    """
    mol = read_smiles('c1ccccc1')
    mol.thiele()
    text = write_mrv(mol, indent=None)
    assert text.count('order="A"') == 6 and 'order="4"' not in text

    text = write_mrv(read_smiles('[Fe]~N(C)(C)C'), indent=None)
    assert '<bond id="b1" atomRefs2="a1 a2" convention="cxn:coord" />' in text
    assert text.count('order=') == 3  # the three N-C bonds only


def test_a_wedge_goes_out_as_the_bare_letter_this_reader_reads_back():
    """``<bondStereo>W</bondStereo>``, of the three spellings the read path accepts: the bare letters are
    MRV's own, and ``dictRef`` and the MDL reference are accommodations this package has not measured."""
    mol = read_smiles('[C@H](N)(O)C')
    mol.clean2d()
    text = write_mrv(mol, indent=None)
    assert '<bondStereo>W</bondStereo>' in text or '<bondStereo>H</bondStereo>' in text
    assert 'dictRef' not in text and 'convention=' not in text


def test_all_three_wedges_a_bond_can_carry_go_out_in_their_own_spelling():
    """Up, down and *either*, three branches the ``or`` in the test above cannot separate.

    A molecule read from a drawing keeps the wedges it was drawn with -- ``wedges_for_write`` returns
    stored ones untouched -- so reading each spelling back and writing it out pins the emitter branch by
    branch.  *Either* has no letter in MRV's vocabulary, so it goes out as the MDL dictionary reference,
    the one place this writer emits a ``convention`` for stereo, and its loss is asserted on the read log.
    """
    def drawn(stereo, log):
        body = _atoms('<atom id="a1" elementType="C" x2="0.0" y2="0.0"/>',
                      '<atom id="a2" elementType="N" x2="1.3" y2="0.75"/>',
                      '<atom id="a3" elementType="O" x2="-1.3" y2="0.75"/>',
                      '<atom id="a4" elementType="F" x2="0.0" y2="-1.5"/>')
        bonds = (f'<bondArray><bond id="b1" atomRefs2="a1 a2" order="1">{stereo}</bond>'
                 f'<bond id="b2" atomRefs2="a1 a3" order="1"/>'
                 f'<bond id="b3" atomRefs2="a1 a4" order="1"/></bondArray>')
        mol, = read_mrv(_wrap(body + bonds), log=log)
        return mol

    log = []
    assert '<bondStereo>W</bondStereo>' in write_mrv(drawn('<bondStereo>W</bondStereo>', log), log=log,
                                                     indent=None)
    assert '<bondStereo>H</bondStereo>' in write_mrv(drawn('<bondStereo>H</bondStereo>', log), log=log,
                                                     indent=None)
    assert log == [], log

    read_log = []
    mol = drawn('<bondStereo convention="MDL" conventionValue="4"/>', read_log)
    assert [x for x in read_log if 'drawn as either' in x], read_log
    log = []
    text = write_mrv(mol, log=log, indent=None)
    assert '<bondStereo convention="MDL" conventionValue="4" />' in text, text
    assert log == [], log


def test_a_wedge_is_written_with_its_narrow_end_first():
    """MDL's rule, and MRV inherits it: the wedge starts at the atom ``atomRefs2`` names first.

    Every other bond goes out low id first, so a wedged bond whose narrow end is the *higher* atom is
    written against that order deliberately -- emitting the bond's own order states the wedge from the
    wrong end, a different configuration and not a different spelling.  The ``narrow > wide`` assertion is
    on the fixture: it failing means the layout or the chooser moved and this test covers nothing.
    """
    mol = read_smiles('C[C@H](Br)CC')
    mol.clean2d()
    (narrow, wide, _), = wedges_for_write(mol, [])[0]
    assert narrow > wide, (narrow, wide)

    text = write_mrv(mol, indent=None)
    assert f'atomRefs2="a{narrow} a{wide}" order="1"><bondStereo>' in text, text
    assert 'atomRefs2="a2 a3"' in text  # an unwedged bond of the same atom, still low id first


def test_an_s_group_a_molecule_carries_is_written_as_a_nested_molecule():
    """The S-group a molecule carries survives a write, as MRV's nested ``<molecule>``, and is not
    reported as lost.

    The honesty claim is held by ``_emit_sgroups``' per-record reporting rather than by this test: a type
    with no MRV role, or a field this dialect cannot spell, earns its own ``unsupported:`` line naming
    *that* record.  The control is the same molecule without the record, which must produce no nested
    ``<molecule>`` at all -- so ``molID="m2"`` is about the S-group and not about every molecule written.
    """
    log = []
    mol = read_smiles('CCO')
    mol.set_sgroups([{'type': b'DAT', 'atoms': (next(iter(mol.atom_numbers)),), 'name': b'FIELD'}])
    text = write_mrv(mol, log=log, indent=None)
    assert log == [], log
    assert ('<molecule molID="m2" id="sg1" role="DataSgroup" atomRefs="a1" fieldName="FIELD" '
            'x="0.0000" y="0.0000" />') in text, text

    log = []
    text = write_mrv(read_smiles('CCO'), log=log, indent=None)
    assert log == [], log
    assert 'molID="m2"' not in text, text


def test_a_data_s_group_is_written_with_the_label_anchor_marvin_reads_it_by():
    """``x``/``y`` on ``<molecule role="DataSgroup">``: the anchor of the drawn label.

    Marvin 25.1.3 needs the pair -- ``molconvert -g mrv`` on a document without it exits 0 and returns
    ``<MDocument></MDocument>``, so the whole record depends on two attributes -- and Marvin's own writer
    states ``x="0.0000" y="0.0000"`` for a group whose source file carried no ``FIELDDISP``.  So the pair
    is always written: the group's own anchor where CTfile stated one, and that zero pair where it did not.

    The zero is not an invented drawing the way ``x2``/``y2`` on an atom would be: it anchors a label at
    the frame's origin and states nothing about where an atom is.
    """
    log = []
    mol = read_smiles('CCO')
    first = next(iter(mol.atom_numbers))
    mol.set_sgroups([{'type': b'DAT', 'atoms': (first,), 'name': b'FIELD',
                      'disp': (1.25, -2.5), 'disp_tail': b'    DA    ALL  1       5'}])
    text = write_mrv(mol, log=log, indent=None)
    assert 'x="1.2500" y="-2.5000"' in text, text
    # The anchor is written; the styling columns after it are what MRV has no spelling for.
    assert [x for x in log if 'FIELDDISP styling' in x], log
    assert not [x for x in log if 'FIELDDISP has no MRV spelling' in x], log

    # ... and the anchor comes back off the same attributes, so an SDF anchor survives the trip out.
    record, = parse_mrv(text, log=[])
    sgroup, = record.ctab.sgroups
    assert sgroup.disp == (1.25, -2.5, ''), sgroup.disp


def test_a_molecule_with_no_layout_gets_no_coordinates_and_its_stereo_is_reported():
    """No invented drawing, and the loss is reported for the configuration that has no other channel.

    ``x2="0.0000" y2="0.0000"`` on every atom is not a neutral default -- it is a drawing, and this reader
    would read stereo out of it.  The tetrahedral configuration's only channel here *is* the geometry, so
    it says so through the same wedge chooser the MDL emitters use.  A double bond is not symmetric with
    it: the letter is a channel the coordinates are not, per the test below.
    """
    log = []
    text = write_mrv(read_smiles('[C@H](N)(O)C'), log=log, indent=None)
    assert 'x2=' not in text and 'y2=' not in text and 'bondStereo' not in text
    assert [x for x in log if 'configured stereocentre(s) but no coordinates' in x], log


def test_a_double_bond_configuration_survives_a_layoutless_write_as_a_bare_letter():
    """The write half of the bare letter, bounded by the same ambiguity test the read half applies.

    With no coordinates the letter is the document's only possible statement and the parity is on the
    container, so writing it is neither repair nor invention.  Bounded, hence the two halves: a bare letter
    names no reference atoms, so on a terminal carrying a second substituent it does not say which pair is
    cis.  ``stated_cis_trans`` refuses to read one and ``cis_trans_letter(framed=False)`` refuses to write
    one, from one condition in :mod:`chython.core.wedge`.  3-methyl-2-pentene is that case.
    """
    log = []
    text = write_mrv(read_smiles('C/C=C/C'), log=log, indent=None)
    assert 'x2=' not in text, 'no drawing was invented to carry it'
    assert '<bondStereo>T</bondStereo>' in text, text
    assert 'atomRefs4' not in text, 'Marvin writes the letter bare and so does this writer'
    assert [x for x in _unsupported(log) if 'double-bond configuration' in x] == [], log

    log = []
    text = write_mrv(read_smiles('C/C=C(\\C)CC'), log=log, indent=None)
    assert 'bondStereo' not in text, text
    assert [x for x in _unsupported(log) if 'double-bond configuration' in x], log


def test_the_bare_letter_this_writer_states_is_the_one_it_reads_back():
    """The round trip the two halves exist for, over a molecule with no drawing.

    Both ways round: checking only trans passes on a writer that emits ``T`` unconditionally.
    """
    for smiles in ('C/C=C/C', 'C/C=C\\C'):
        mol = read_smiles(smiles)
        log = []
        text = write_mrv(mol, log=log, indent=None)
        assert log == [], log
        back, = read_mrv(text)
        assert back == mol, smiles


def test_a_conformer_is_written_back_in_the_three_dimensional_attributes():
    """A record's ``dimensionality`` decides the set, so a 3D record does not come back flattened.

    The container holds x and y only, so this is reachable only by writing a *record*; a writer asking
    the molecule instead loses the third coordinate silently.
    """
    record, = parse_mrv(_wrap(_atoms('<atom id="a1" elementType="C" x3="1.0" y3="2.0" z3="-3.5"/>')))
    text = write_mrv(record, indent=None)
    assert 'x3="1.0000" y3="2.0000" z3="-3.5000"' in text
    assert 'x2=' not in text


def test_a_record_is_written_back_with_the_file_s_own_atom_ids_and_statements(data):
    """Writing what ``parse_mrv`` read reproduces the file's statements, not this tree's derivations.

    The fixture's drawn hydrogen states no ``hydrogenCount`` and the written record states none, where
    writing the *molecule* would state the derived count.  The ids are the visible half of the same thing.
    """
    record, = parse_mrv(data('mrv_hydrogens.mrv'))
    text = write_mrv(record, indent=None)
    assert 'title="methylamine"' in text
    assert '<atom id="a3" elementType="H" x2="1.2990" y2="2.2500" />' in text
    assert text.count('hydrogenCount=') == 2  # the carbon's 3 and the nitrogen's 1, and not the H's


def test_the_round_trip_returns_the_same_structure():
    """Every quantity at once, compared as structures: two SMILES strings differ for reasons that are
    not chemistry, so a string comparison fails on those and passes on a lost isotope."""
    for smiles in ('CCO', '[13CH3][NH3+]', 'C[O-].[Na+]', 'CC(=O)Nc1ccccc1', '[Fe]~N(C)(C)C',
                   'C#CC(F)(Cl)Br', '[CH3]C |^1:0|'):
        mol = read_smiles(smiles)
        mol.thiele()
        log = []
        back, = read_mrv(write_mrv(mol, log=log), log=log)
        assert back == mol, (smiles, log)
        assert not _unsupported(log), (smiles, log)


def test_a_document_this_writer_produced_routes_back_to_this_dialect():
    """The namespace the writer declares is the one ``sniff`` matches, with no naming line -- a writer
    emitting a namespace nothing claims round-trips perfectly and is unreadable to everyone else."""
    log = []
    dial = sniff(parse_xml(write_mrv(read_smiles('CCO'))), log)
    assert dial.name == 'mrv'
    assert log == [], log


def test_a_title_is_written_where_the_reader_looks_for_it():
    """``title`` on ``<molecule>``, and an empty one is no attribute rather than ``title=""``."""
    assert 'title="ethanol"' in write_mrv(read_smiles('CCO'), title='ethanol', indent=None)
    assert 'title=' not in write_mrv(read_smiles('CCO'), title='', indent=None)


def test_the_package_exports_this_dialect_under_a_name_that_names_it():
    """``mrv_record_from_molecule``, and no bare ``record_from_molecule`` for either dialect.

    The two builders disagree -- MRV states the implicit hydrogen count and CML the total -- so a bare
    name would be whichever import line came second, and the shadow has no symptom at the import.
    """
    from ... import xml

    assert xml.MRV_NS == MRV_NS and xml.MRV.name == 'mrv'
    for name in ('parse_mrv', 'read_mrv', 'write_mrv', 'write_mrv_element', 'mrv_record_from_molecule',
                 'cml_record_from_molecule'):
        assert name in xml.__all__ and hasattr(xml, name), name
    # The attribute and not only `__all__`: a bare name bound on the module is reachable by an explicit
    # import whatever `__all__` says.
    assert 'record_from_molecule' not in xml.__all__ and not hasattr(xml, 'record_from_molecule')
    assert xml.mrv_record_from_molecule is record_from_molecule


def test_a_non_utf8_title_is_reported_and_the_record_is_still_written_mrv():
    """The same loss CML takes, for the same reason: XML admits no lone surrogate."""
    from ....core import MoleculeContainer

    mol = MoleculeContainer()
    with mol.edit() as e:
        e.add_atom('C')
    mol.set_title(b'caf\xe9')
    log = []
    out = write_mrv(mol, log=log)
    assert any(str(x).startswith('unsupported: ') and 'not valid UTF-8' in x for x in log), log
    assert 'caf�' in out


def test_data_fields_go_out_as_a_property_list_and_come_back_on_meta():
    """The channel Marvin's own ``sdf``-to-``mrv`` conversion writes.

    Both name attributes go out because Marvin 25.1.3 writes both; a reader taking either gets the name.
    """
    molecule = read_smiles('CCO')
    molecule.meta['BATCH_ID'] = 'lot-42'
    text = write_mrv(molecule, indent=None)
    assert '<propertyList><property dictRef="BATCH_ID" title="BATCH_ID">' in text
    # Before the arrays, where Marvin puts it.
    assert text.index('<propertyList>') < text.index('<atomArray')
    assert read_mrv(text)[0].meta == {'BATCH_ID': 'lot-42'}


def test_a_value_with_a_space_in_it_goes_out_inside_a_cdata_section():
    """Measured against Marvin 25.1.3: a plain-text ``<scalar>`` value comes back as its first
    whitespace-delimited token, and the same value in a CDATA section comes back whole.  A section and
    escaped text are one document to a conforming parser, so this costs nothing.
    """
    molecule = read_smiles('C')
    molecule.meta['NOTES'] = 'first line\nsecond line'
    molecule.meta['MARKUP'] = 'a & b <c>'
    text = write_mrv(molecule, indent=None)
    assert '<scalar><![CDATA[first line\nsecond line]]></scalar>' in text
    # Verbatim inside the section: escaping it there would put `&amp;` in the value.
    assert '<scalar><![CDATA[a & b <c>]]></scalar>' in text
    assert read_mrv(text)[0].meta == {'NOTES': 'first line\nsecond line', 'MARKUP': 'a & b <c>'}


def test_a_value_holding_the_end_of_a_cdata_section_is_written_as_two_sections():
    """``]]>`` has no spelling inside a section, and two sections are one value to every parser."""
    molecule = read_smiles('C')
    molecule.meta['K'] = 'before ]]> after'
    text = write_mrv(molecule, indent=None)
    assert '<scalar><![CDATA[before ]]]]><![CDATA[> after]]></scalar>' in text
    assert read_mrv(text)[0].meta == {'K': 'before ]]> after'}


def test_an_r_atom_round_trips_through_this_dialect():
    """Marvin's spelling both ways, index included, so a molecule read as SMILES writes as MRV."""
    mol = read_smiles('[R7]C')
    text = write_mrv(mol)
    assert 'elementType="R"' in text and 'rgroupRef="7"' in text, text
    again, = read_mrv(text, log=[])
    r = next(a for a in again.atoms() if a.is_r)
    assert r.r_index == 7
    assert [a.atomic_symbol for a in again.atoms()] == [a.atomic_symbol for a in mol.atoms()]
