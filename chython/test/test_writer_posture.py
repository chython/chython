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
"""A WRITER WRITES THE RECORD and says what it could not carry.  Tree-wide, one rule, every writer.

A refusal is answer-boundary behaviour and a serialiser is not an answer boundary: a loop writing forty
thousand records must not be stopped by a field one of them holds.  Three classes, decided by what a
legal record would say and never by how much is lost:

| Class                                                 | Behaviour                                |
| ----------------------------------------------------- | ---------------------------------------- |
| **No field for the data**, a legal record is a subset | write, log one sentence naming the field |
| **No legal record exists**                            | raise, naming the fix                    |
| **Degraded geometry**                                 | normalize, then write                    |

Class 2 is narrow by construction -- *no legal record exists*, never *this field would be lossy*.  A
capacity a format can dodge by choosing a version is class 1 and the writer's choice to make, which is
why `mol()` answers 1200 atoms with a V3000 record and refuses only the pinned `version=2000`.

This file is the tree-wide half: one class-1 case per writer and the class-2 boundary beside it, in one
place, so a writer added later is measured against the same rule.  The per-format detail lives beside
each writer -- `formats/ctfile/test/test_v2000.py` for the V2000 collection line, `core/test/test_pach.py`
and `test_pach3.py` for pach's field-by-field table.

NOT WRITERS UNDER THIS RULE, stated rather than left to be inferred: an identity string (`smiles`,
`inchi`, `inchikey`) is a projection of the container and not a record of it -- being lossy is the point,
and a line per record would say nothing a caller did not ask for -- and `to_bytes()` is the arena
verbatim, which loses nothing but `meta`.
"""
from io import StringIO
from pytest import mark, raises
from chython import (ERDFWrite, ESDFWrite, FIELDDISP_TAIL, RDFWrite, SDFWrite, STEREO_AND,
                     add_data_sgroup, cml, inchi, mol, mrv, pach, read_reaction_smiles, rxn, smiles)
# `MalformedCtfile` is the CTfile package's own error and is not a facade name: the class-2 refusals
# below are the format's, so the test names the format's exception rather than a base class that would
# also catch a `ValueError` from anywhere else in the call.
from chython.formats.ctfile import MalformedCtfile


def _rich():
    """2-amino-2-phenylacetic acid carrying every field some writer below has no room for.

    Aromatic bonds, a stated parity, an AND collection over it, a stored CIP descriptor, a map number,
    a data S-group with FIELDDISP styling, metadata and a plane.
    """
    m = smiles('c1ccccc1[C@H](N)C(=O)O')
    m.meta['ID'] = 'phenylglycine'
    centre = [a.n for a in m.atoms() if a.parity][0]
    first = m.atom_numbers[0]
    with m.edit():
        m.set_stereo_group(centre, STEREO_AND, 1)
        m.set_atom_cip(centre, 'R')
        m.set_map_number(first, 7)
    add_data_sgroup(m, 'FOO', 'bar', atoms=[first], disp=(1., 2., FIELDDISP_TAIL))
    m.clean2d()
    return m


def _stream(writer, record):
    buf = StringIO()
    with writer(buf) as f:
        log = f.write(record)
    return buf.getvalue(), log


# ----------------------------------------------------------------------------------------------
# nothing refuses it
# ----------------------------------------------------------------------------------------------

@mark.parametrize('name,write', [
    ('mol', lambda m: mol(m)),
    ('mol v2000', lambda m: mol(m, version=2000)),
    ('mol v3000', lambda m: mol(m, version=3000)),
    ('SDFWrite', lambda m: _stream(SDFWrite, m)[0]),
    ('ESDFWrite', lambda m: _stream(ESDFWrite, m)[0]),
    ('mrv', lambda m: mrv(m)),
    ('cml', lambda m: cml(m)),
    ('pach', lambda m: pach(m)),
    ('pach v2', lambda m: pach(m, version=2)),
    ('pach v4', lambda m: pach(m, version=4)),
    ('pack', lambda m: m.pack()),
    ('to_bytes', lambda m: m.to_bytes()),
    ('depict', lambda m: m.depict()),
])
def test_no_writer_refuses_a_record_it_cannot_fully_carry(name, write):
    """The molecule holds something every one of these has no field for, and every one writes it."""
    assert write(_rich()), name


@mark.parametrize('name,write', [
    ('rxn', lambda r: rxn(r)),
    ('RDFWrite', lambda r: _stream(RDFWrite, r)[0]),
    ('ERDFWrite', lambda r: _stream(ERDFWrite, r)[0]),
    ('pach', lambda r: pach(r)),
    ('pack', lambda r: r.pack()),
])
def test_a_reaction_is_written_by_the_same_rule(name, write):
    assert write(_reaction()), name


def _reaction():
    """A reduction whose metadata no structure record has a field for."""
    reaction = read_reaction_smiles('CC=O.[H][H]>>CCO')
    reaction.meta['note'] = 'trailing '        # whitespace the RDfile reader eats on the way back
    return reaction


# ----------------------------------------------------------------------------------------------
# and each one says what it left out
# ----------------------------------------------------------------------------------------------

def test_v2000_names_the_collection_it_has_no_syntax_for():
    log = []
    assert mol(_rich(), version=2000, log=log)
    assert log[0].rule == 'v2000:enhanced-stereo-not-written', log
    assert 'AND/OR' in log[0].message


def test_v3000_carries_the_collection_and_says_nothing_about_it():
    """The other half of the pair: a loss line is a fact about the format and not about the writer, so
    the version with a `COLLECTION` block reports no collection.  The one line it does write is about
    the metadata framing -- `mol()` writes the fields and no `$$$$` to walk them by -- and V2000 writes
    that same line, which is what makes the collection line the difference between the two."""
    log = []
    assert mol(_rich(), version=3000, log=log)
    assert [x.rule for x in log] == ['ctfile:data-fields-after-end'], log


def test_cml_names_the_sgroup_and_the_map_numbers():
    log = []
    assert cml(_rich(), log=log)
    assert {'cml:sgroup-not-written', 'cml:map-numbers-not-written'} <= {x.rule for x in log}, log


def test_mrv_names_the_fielddisp_styling_it_has_no_attribute_for():
    """MRV holds the S-group, the map numbers and the anchor, so the styling columns after the anchor
    are the one thing left to report."""
    log = []
    assert mrv(_rich(), log=log)
    assert [x.rule for x in log] == ['mrv:sgroup-field-not-written'], log
    assert 'FIELDDISP styling' in log[0].message


def test_pach_names_each_field_on_the_molecules_own_log():
    """`pach` has no `log=`: it is reached through a container method, and `molecule.log` is the one
    destination.  The stage is the writer's name, so `by_stage('pach')` answers for the record."""
    m = _rich()
    assert pach(m)
    assert sorted(x.rule for x in m.log) == ['pach:cip-lost', 'pach:meta-lost', 'pach:sgroups-lost']
    assert {x.stage for x in m.log} == {'pach'}
    # version 2 is narrower, and every field it cannot hold is named in the same way
    m = _rich()
    assert pach(m, version=2)
    assert sorted(x.rule for x in m.log) == ['pach:cip-lost', 'pach:map-number-lost', 'pach:meta-lost',
                                             'pach:sgroups-lost', 'pach:stereo-groups-lost']


def test_the_rdfile_writer_names_what_the_reader_will_eat():
    """Whitespace at a value's edge and a `$`-led continuation are both written and both reported: the
    record is legal and a reader will read it back changed, which is the caller's to know."""
    reaction = _reaction()
    reaction.meta['multi'] = 'first line\n$DTYPE not a keyword'
    text, log = _stream(RDFWrite, reaction)
    assert text
    assert {'rdf:trailing-whitespace', 'rdf:dollar-in-value'} <= {x.rule for x in log}, log


# ----------------------------------------------------------------------------------------------
# class 2: the pinned container, and nothing else
# ----------------------------------------------------------------------------------------------

def _dodecane_times_a_hundred():
    """1200 atoms, past the V2000 3-character count field and nothing else's limit."""
    return smiles('C' * 1200)


def test_a_capacity_the_writer_can_dodge_is_not_a_refusal():
    big = _dodecane_times_a_hundred()
    assert big.atom_count == 1200
    assert '999 V3000' in mol(big).splitlines()[3], 'the writer chooses the version that fits'
    assert _stream(ESDFWrite, big)[0]


def test_a_pinned_version_that_cannot_hold_the_record_is_the_refusal_that_stays():
    """A caller naming a version has named a container, and a V2000 stamp over a CTAB the format cannot
    spell is the false record class 2 exists to prevent.  Both spellings of the pin refuse."""
    big = _dodecane_times_a_hundred()
    with raises(MalformedCtfile, match='will not fit the V2000'):
        mol(big, version=2000)
    with raises(MalformedCtfile, match='will not fit the V2000'):
        _stream(SDFWrite, big)


def test_a_coordinate_past_the_v2000_column_is_the_same_boundary():
    m = smiles('CC')
    with m.edit():
        for i, sid in enumerate(m):
            m.set_xy(sid, 100000. * (i + 1), 0.)
    assert '999 V3000' in mol(m).splitlines()[3]
    with raises(MalformedCtfile, match='10-character column'):
        mol(m, version=2000)


# ----------------------------------------------------------------------------------------------
# the all-or-nothing caller
# ----------------------------------------------------------------------------------------------

def test_strict_is_the_all_or_nothing_switch_and_pach_is_where_it_is():
    """One writer offers the refusal as an option, because one writer is a store rather than a document:
    a caller keeping pach records as the structure's only copy may want a field it cannot hold to stop
    the write.  The text writers have no such switch -- the record they produce is readable either way.
    """
    with raises(ValueError, match='cip descriptor'):
        pach(_rich(), strict=True)
    for writer in (mol, mrv, cml):
        with raises(TypeError, match='strict'):
            writer(_rich(), strict=True)


def test_a_field_named_at_the_door_is_a_waiver_and_logs_nothing():
    m = _rich()
    assert pach(m, drop='*')
    assert m.log == [], 'the caller already said it; saying it back is noise'


def test_an_explicit_version_4_is_a_waiver_of_the_drawing():
    """`4` is the coordinate-free layout asked for by name, so the plane it does not carry is not
    reported either -- the same rule as `drop=`, stated by a version instead of by a field name."""
    m = _rich()
    assert m.has_coordinates
    assert pach(m, version=4)
    assert 'pach:coordinates-lost' not in [x.rule for x in m.log]


# ----------------------------------------------------------------------------------------------
# the declared boundary
# ----------------------------------------------------------------------------------------------

def test_an_identity_string_is_not_a_writer():
    """`inchi` answers what a compound IS, so a molecule it cannot name is a refusal and not a loss
    line: the R marker stands for a substituent nobody stated, and an identifier that quietly read it
    as something else would be an identity for a different compound.  The file writers are unaffected
    by the same molecule -- that is the whole distinction.
    """
    marked = smiles('[R]C(=O)O')
    with raises(ValueError, match='R marker'):
        inchi(marked)
    assert mol(marked)
    assert mrv(marked)
    assert pach(marked)
