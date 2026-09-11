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
"""Every record a CML or MRV reader makes about a molecule it returns is on that molecule's `log`.

One document holds many ``<molecule>`` elements and the reader's `log=` list is one flat sequence for
all of them, so this is where per-record scoping is easiest to get wrong: the walk writes to the
record's own ``Ctab.log``, which is what the build folds onto the molecule, and the caller's list gets a
copy.  A document-level line -- an entity policy, a dialect fallback, a ``<reaction>`` whose roles this
layer does not model -- names no molecule and stays on the caller's list alone.
"""

from chython.formats.xml import cml, mrv, parse_cml, read_cml, read_mrv, read_xml


#: Three molecules, and only the middle one is damaged: an unmodelled child element (a parse-phase
#: line) and a bond written twice (a build-phase line), so both halves are checked for scope at once.
_CML = '''<cml>
 <molecule id="m1"><atomArray><atom id="a1" elementType="C" x2="0" y2="0"/></atomArray></molecule>
 <molecule id="m2">
  <atomArray>
   <atom id="a1" elementType="C" x2="0" y2="0"/>
   <atom id="a2" elementType="O" x2="1" y2="0"><bogus/></atom>
  </atomArray>
  <bondArray>
   <bond atomRefs2="a1 a2" order="1"/>
   <bond atomRefs2="a2 a1" order="1"/>
  </bondArray>
 </molecule>
 <molecule id="m3"><atomArray><atom id="a1" elementType="C" x2="0" y2="0"/></atomArray></molecule>
</cml>'''

_MRV = '''<MDocument><MChemicalStruct>
 <molecule molID="m1"><atomArray>
  <atom id="a1" elementType="C" x2="0" y2="0" bogusAttr="7"/>
 </atomArray></molecule>
</MChemicalStruct></MDocument>'''


def test_a_repair_is_on_the_molecule_with_nothing_passed_in():
    molecules = read_cml(_CML)
    assert [x.rule for x in molecules[1].log] == ['xml:element-not-modelled', 'ctab:duplicate-bond']
    assert {x.stage for x in molecules[1].log} == {'read'}


def test_only_the_damaged_molecule_of_a_document_carries_the_records():
    """The parse phase and the build phase are both scoped; a flat fold would put 2 on all three."""
    molecules = read_cml(_CML)
    assert [len(x.log) for x in molecules] == [0, 2, 0], [list(x.log) for x in molecules]


def test_the_callers_list_gets_the_documents_lines_once():
    log = []
    molecules = read_cml(_CML, log=log)
    assert [str(x) for x in log] == [str(x) for x in molecules[1].log]


def test_the_facade_and_the_dialect_free_reader_agree():
    for molecules in (cml(_CML), read_xml(_CML)):
        assert [len(x.log) for x in molecules] == [0, 2, 0]


def test_mrv_records_on_the_molecule_too():
    molecule, = read_mrv(_MRV)
    assert [x.rule for x in molecule.log] == ['xml:attribute-not-modelled']
    assert molecule.log[0].stage == 'read'
    assert [x.rule for x in mrv(_MRV)[0].log] == ['xml:attribute-not-modelled']


def test_a_document_level_line_names_no_molecule(data):
    """MRV states a reaction's roles and this layer models none of it -- a fact about the document.

    It must not be copied onto the three molecules the document does yield: a line every molecule
    carries answers no question about any of them.
    """
    log = []
    molecules = read_mrv(data('mrv_reaction.mrv'), log=log)
    assert [x.rule for x in log] == ['mrv:reaction-roles-not-modelled', 'mrv:document-furniture']
    assert molecules and not any(x.log for x in molecules)


def test_a_parsed_record_carries_the_walks_lines_to_its_build():
    """`parse_cml` stops before the build, so the record's own log is where the walk left them."""
    records = parse_cml(_CML)
    assert [len(x.ctab.log) for x in records] == [0, 1, 0]
    molecule, _, log = records[1].ctab.build()
    assert [x.rule for x in molecule.log] == ['xml:element-not-modelled', 'ctab:duplicate-bond']
    assert [str(x) for x in log] == [str(x) for x in molecule.log]


def test_building_one_record_twice_gives_each_molecule_the_same_records():
    """`Ctab.build` copies `ctab.log` and never appends to it, so a second build is not cumulative."""
    record = parse_cml(_CML)[1]
    first, _, _ = record.ctab.build()
    second, _, _ = record.ctab.build()
    assert [x.rule for x in first.log] == [x.rule for x in second.log] == \
           ['xml:element-not-modelled', 'ctab:duplicate-bond']
