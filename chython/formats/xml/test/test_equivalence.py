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
"""The same molecule read through CML and through MDL gives the same answer.

``test/cml_stereo.mol`` and ``test/cml_stereo.cml`` are one drawing written twice -- alanine, same atom
order and coordinates -- with the wedge as MDL bond-block stereo ``1`` in one and
``<bondStereo>W</bondStereo>`` in the other.  The second half is a round trip over every written field.
"""

from pytest import mark

from chython.core import read_smiles, write_smiles
from chython.formats.ctfile import parse_record

from .._cml import parse_cml, read_cml, write_cml
from .._dialect import read_xml
from .._mrv import read_mrv


def _mdl(path):
    """One molfile as a molecule plus its log, through the CTfile reader."""
    log = []
    with open(path, encoding='utf8') as f:
        molecule = parse_record(f.read().splitlines(), log)
    return molecule, log


def _cml(path):
    """The same, through this package."""
    log = []
    molecules = read_cml(path, log=log)
    assert len(molecules) == 1
    return molecules[0], log


def _parities(mol):
    """``{stable id: parity}`` for every configured centre.  The comparable form of "what stereo"."""
    return {sid: mol.parity_of(sid) for sid in mol.atom_numbers if mol.parity_of(sid)}


# the equivalence bar

def test_the_wedge_and_the_bond_stereo_give_the_same_configuration(data):
    """Same SMILES, same parity map, both logs empty.

    The empty logs matter as much as the parities: a reader reaching the right answer by way of a
    repair would be agreeing for a different reason.
    """
    mdl, mdl_log = _mdl(data('cml_stereo.mol'))
    cml, cml_log = _cml(data('cml_stereo.cml'))
    assert write_smiles(mdl) == write_smiles(cml)
    assert _parities(mdl) == _parities(cml)
    assert mdl_log == [] and cml_log == [], (mdl_log, cml_log)


def test_the_configuration_is_actually_there(data):
    """The control for the test above: two readers perceiving *nothing* would agree perfectly, so the
    answer has to be non-empty -- one tetrahedral centre, on the atom the wedge is drawn from."""
    cml, _ = _cml(data('cml_stereo.cml'))
    assert _parities(cml) == {1: 2}
    assert '@' in write_smiles(cml)


def test_the_two_files_describe_the_same_drawing(data):
    """The two fixtures share atom order and coordinates -- asserted, not trusted, since an edit
    making it false would surface as a stereo bug."""
    mdl, _ = _mdl(data('cml_stereo.mol'))
    cml, _ = _cml(data('cml_stereo.cml'))
    assert [mdl.element_of(s) for s in mdl.atom_numbers] == [cml.element_of(s) for s in cml.atom_numbers]
    for a, b in zip(mdl.atom_numbers, cml.atom_numbers):
        assert mdl.xy_of(a) == cml.xy_of(b), (a, b)


def test_the_parity_channel_agrees_with_the_wedge_channel(data):
    """An ``<atomParity>`` gives the configuration the wedge drew, which pins the sign calibration.

    The value ``-1`` in the frame ``atomRefs4="a2 a3 a4 a1"`` is measured, not chosen: it is what the
    writer emits for this molecule, and ``1`` is its mirror.  Hence the literal frame.
    """
    drawn, _ = _cml(data('cml_stereo.cml'))
    with open(data('cml_stereo.cml'), encoding='utf8') as f:
        text = f.read()
    stated = text.replace('<bondStereo>W</bondStereo>', '').replace(
        '<atom id="a1" elementType="C" x2="0.0000" y2="0.0000"/>',
        '<atom id="a1" elementType="C" x2="0.0000" y2="0.0000">'
        '<atomParity atomRefs4="a2 a3 a4 a1">-1</atomParity></atom>')
    log = []
    molecules = read_cml(stated, log=log)
    assert _parities(molecules[0]) == _parities(drawn), log


def test_the_negated_parity_gives_the_other_configuration(data):
    """The discriminator: the two signs in one frame must give opposite centres, or the test above
    would pass against a reader that ignored the value and read the drawing."""
    with open(data('cml_stereo.cml'), encoding='utf8') as f:
        text = f.read()
    base = text.replace('<bondStereo>W</bondStereo>', '')
    out = {}
    for value in ('1', '-1'):
        doc = base.replace(
            '<atom id="a1" elementType="C" x2="0.0000" y2="0.0000"/>',
            f'<atom id="a1" elementType="C" x2="0.0000" y2="0.0000">'
            f'<atomParity atomRefs4="a2 a3 a4 a1">{value}</atomParity></atom>')
        out[value] = _parities(read_cml(doc, log=[])[0])
    assert out['1'] and out['-1']
    assert out['1'] != out['-1'], out


def test_the_marvin_fixture_agrees_with_the_mdl_wedge_codes(data):
    """The MDL wedge vocabulary reaching CML through the ``convention="MDL"`` escape hatch.

    ``conventionValue="4"`` is the molfile's "either" code, ``"3"`` its "cis or trans, unknown which";
    both decode through the CTfile package's own ``WEDGE_FROM_V2000`` rather than a copy of it.
    """
    from .._cml import WEDGE_FROM_V2000
    log = []
    read_cml(data('cml_marvin.cml'), log=log)
    assert set(WEDGE_FROM_V2000) == {0, 1, 4, 6}
    assert any('drawn as either' in x for x in log), log


# the round trip

#: One molecule per field the writer has to carry.  Named rather than generated: each entry is here for
#: the channel it exercises, and a random corpus would exercise the same three every time.
ROUND_TRIP = [
    ('a plain chain', 'CCO'),
    ('a charge', 'CC(=O)[O-]'),
    ('a cation', 'C[N+](C)(C)C'),
    ('a radical', 'C |^1:0|'),
    ('an isotope', '[13CH4]'),
    ('an isotope and a charge together', '[15NH4+]'),
    ('an aromatic ring', 'c1ccccc1O'),
    ('a fused aromatic system', 'c1ccc2ccccc2c1'),
    ('a tetrahedral centre', 'C[C@H](N)C(=O)O'),
    ('a hydrogen count worth stating', '[nH]1cccc1'),
    ('a metal salt', '[Na+].CC(=O)[O-]'),
]

#: What the round trip is allowed to say, per entry, and nothing else.  An allow-list, so a line that
#: is not an accounted-for boundary is a loss nobody has explained.
EXPECTED_LOG = {
    # The writer's line only: the configuration goes out as an `<atomParity>` and comes back from one, so
    # the read half is silent.  Information rather than loss, hence the exact text pinned.
    'a tetrahedral centre': ['1 configured stereocentre(s) but no coordinates; no wedges written'],
}


@mark.parametrize('why,smiles', ROUND_TRIP, ids=[x[0].replace(' ', '_') for x in ROUND_TRIP])
def test_a_molecule_survives_being_written_and_read(why, smiles):
    """``read(write(mol))`` is the same molecule, for each field in turn.

    Compared by ``write_smiles`` of both sides rather than against a literal, which would also pin the
    SMILES writer's traversal.
    """
    before = read_smiles(smiles)
    log = []
    after, = read_cml(write_cml(before, log=log), log=log)
    assert write_smiles(after) == write_smiles(before), (why, log)


@mark.parametrize('why,smiles', ROUND_TRIP, ids=[x[0].replace(' ', '_') for x in ROUND_TRIP])
def test_a_round_trip_says_only_what_it_is_allowed_to(why, smiles):
    """The trip says nothing on the way that ``EXPECTED_LOG`` does not account for.

    Compared exactly rather than filtered by substring: a construct is either applied or named, so an
    unaccounted line is the failure.  Separate from the test above so a failure names which claim broke.
    """
    log = []
    read_cml(write_cml(read_smiles(smiles), log=log), log=log)
    assert [str(x) for x in log] == EXPECTED_LOG.get(why, []), (why, log)


def test_a_double_bond_configuration_with_no_layout_survives_in_the_letter_and_says_nothing():
    """The one round trip CML can make and a molfile cannot: a coordinate-free cis/trans descriptor.

    With no coordinates the drawing states nothing, so ``<bondStereo>C``/``T`` with an ``atomRefs4`` is
    the only channel.  Asserted as a pair -- the configuration survives and the log is empty -- since
    either half alone is satisfied by the wrong code.
    """
    before = read_smiles('C/C=C/C')
    assert '/' in write_smiles(before) or '\\' in write_smiles(before)
    log = []
    document = write_cml(before, log=log)
    assert 'x2=' not in document      # armed: with a layout the letter would not be the only channel
    after, = read_cml(document, log=log)
    assert write_smiles(after) == write_smiles(before)
    assert log == [], log


#: *trans*-2-butene at hand-laid coordinates, the two methyls on opposite sides of the C2=C3 axis.  Not
#: from a SMILES: the point is a molecule that has a layout, and one built from SMILES has none.
DRAWN_BUTENE = """<cml xmlns="http://www.xml-cml.org/schema">
  <molecule id="m1">
    <atomArray>
      <atom id="a1" elementType="C" x2="0.00" y2="0.00"/>
      <atom id="a2" elementType="C" x2="0.87" y2="0.50"/>
      <atom id="a3" elementType="C" x2="1.73" y2="0.00"/>
      <atom id="a4" elementType="C" x2="2.60" y2="0.50"/>
    </atomArray>
    <bondArray>
      <bond id="b1" atomRefs2="a1 a2" order="1"/>
      <bond id="b2" atomRefs2="a2 a3" order="2"/>
      <bond id="b3" atomRefs2="a3 a4" order="1"/>
    </bondArray>
  </molecule>
</cml>"""


def test_a_drawn_double_bond_configuration_round_trips_through_the_coordinates():
    """The discriminator for the test above: the loss is about the layout, not about double bonds.

    With coordinates the descriptor goes out through both channels, as a Marvin file does, and they
    cannot disagree -- the letter is derived from the parity, which was read from the drawing.
    """
    log = []
    before, = read_cml(DRAWN_BUTENE, log=log)
    assert log == [], log
    assert '/' in write_smiles(before) or '\\' in write_smiles(before)
    out = write_cml(before, log=log)
    assert log == [], log
    assert '<bondStereo atomRefs4="a1 a2 a3 a4">T</bondStereo>' in out
    assert 'x2=' in out
    after, = read_cml(out, log=log)
    assert write_smiles(after) == write_smiles(before)
    assert log == [], log


def test_a_three_dimensional_structure_survives_as_a_drawing(data):
    """A 3D record is read into the layout -- the arena stores x and y -- so a conformer written back
    out is a projection.  The molecule survives exactly, the third coordinate does not."""
    log = []
    water = read_cml(data('cml_quirks.cml'), log=log)[6]
    again, = read_cml(write_cml(water, log=[]), log=[])
    assert write_smiles(again) == write_smiles(water)
    assert [again.xy_of(s) for s in again.atom_numbers] == [water.xy_of(s) for s in water.atom_numbers]


def test_a_stereocentre_survives_the_round_trip_by_both_channels():
    """With no coordinates a wedge has no drawing to sit in, so the parity is the only channel.  A
    drawn molecule gets both, and they agree: the two order translations are the same permutation
    applied twice, and cancel."""
    mol = read_smiles('C[C@H](N)C(=O)O')
    out = write_cml(mol, log=[])
    assert '<atomParity' in out and 'x2=' not in out
    again, = read_cml(out, log=[])
    assert write_smiles(again) == write_smiles(mol)
    assert _parities(again) == _parities(mol)


def test_a_drawn_stereocentre_is_written_with_a_wedge_and_a_parity(data):
    """Both channels on a drawn molecule: a CML file stating only one loses whichever consumer reads
    the other."""
    drawn, _ = _cml(data('cml_stereo.cml'))
    out = write_cml(drawn, log=[])
    assert '<bondStereo>' in out and '<atomParity' in out
    again, = read_cml(out, log=[])
    assert _parities(again) == _parities(drawn)


def test_the_round_trip_is_stable_under_repetition():
    """Twice through equals once through: a writer that shifted an atom order or a wedge's narrow end
    would converge on the second pass, which one trip cannot see."""
    mol = read_smiles('C[C@H](N)C(=O)O')
    once, = read_cml(write_cml(mol, log=[]), log=[])
    twice, = read_cml(write_cml(once, log=[]), log=[])
    assert write_cml(once, log=[]) == write_cml(twice, log=[])


def _shape(record):
    """Everything a record states, as one comparable value.

    Over every field the intermediate has, including the hydrogen count parked in the spill rather than
    in a ``CtabAtom`` slot: a subset comparison passes while the trip loses data.
    """
    return (record.ids, record.ctab.dimensionality,
            [(a.element, a.charge, a.isotope, a.radical, a.x, a.y, a.z, a.parity, a.stated_h,
              spill.get('hydrogen_total'))
             for a, spill in zip(record.ctab.atoms, record.atom_extras)],
            [(b.a, b.b, b.order, b.wedge) for b in record.ctab.bonds])


def test_a_record_round_trips_with_everything_it_states(data):
    """Through ``parse_cml`` and back, including the file's own atom identifiers.

    A molecule cannot carry those ids -- the arena has its own -- so this is the only path that
    preserves what the file called its atoms.  Every record in the fixture goes through.
    """
    log = []
    for before in parse_cml(data('cml_quirks.cml'), log=log):
        after, = parse_cml(write_cml(before, log=[]), log=[])
        assert _shape(after) == _shape(before)


def test_the_mrv_fixture_at_least_tokenizes(data):
    """``test/implicit.mrv`` reads through :func:`read_xml`, which asks the file what it is rather
    than being told -- the weakest claim there is: whatever routes it, the file reads."""
    log = []
    molecules = read_xml(data('implicit.mrv'), log=log)
    assert molecules
    for mol in molecules:
        assert write_smiles(mol)


def test_the_named_reader_and_the_sniffing_reader_agree_on_the_mrv_fixture(data):
    """Routing is a routing decision, not a different reader: sniffed, named or read as CML the
    molecules are identical, since MRV spells these atoms and bonds as CML does.  Compared as
    structures -- a SMILES string is not an identity in this tree."""
    named = read_mrv(data('implicit.mrv'), log=[])
    sniffed = read_xml(data('implicit.mrv'), log=[])
    as_cml = read_cml(data('implicit.mrv'), log=[])
    assert sniffed == named
    assert sniffed == as_cml


def test_a_document_in_an_unclaimed_namespace_says_which_dialect_read_it():
    """The fallback names the dialect that read a document in a namespace nobody claims.

    Prefixed ``unsupported: `` rather than as damage: the document is good XML in a vocabulary this
    library has not written.  Synthetic on purpose -- a real fixture stops being unclaimed the day its
    dialect lands.
    """
    log = []
    molecules = read_xml('<cml xmlns="urn:example:unclaimed"><molecule><atomArray>'
                         '<atom id="a1" elementType="C"/></atomArray></molecule></cml>', log=log)
    assert molecules
    matched = [x for x in log if 'no dialect claims' in x]
    assert len(matched) == 1, log
    assert str(matched[0]).startswith('unsupported: '), matched
    assert 'urn:example:unclaimed' in matched[0] and 'read as cml' in matched[0], matched
