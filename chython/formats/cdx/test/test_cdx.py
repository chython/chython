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
"""CDX and CDXML: the vocabulary, both tokenizers, the reader's rules and the writer, on hand-drawn
public compounds.  Documents are built by :func:`document` from node and bond tuples, so each test
states only the attributes it is about."""

from struct import pack

from pytest import mark, raises

from chython.core import STEREO_ABS, STEREO_AND, STEREO_OR, smiles

from .. import CdxError, MalformedCdx, UnsupportedCdx, cdx, cdxml
from .._binary import HEADER, parse_cdx, serialize_cdx
from .._cdxml import parse_cdxml
from .._vocabulary import vocabulary
from ...ctfile import mol


def document(nodes, bonds=(), extra=''):
    """CDXML text.  `nodes` is ``(id, x, y, attributes)``, `bonds` is ``(id, begin, end, attributes)``;
    attributes are a ready-made string such as ``'Element="8"'``."""
    body = ''.join(f'<n id="{i}" p="{x} {y}" {a}/>' for i, x, y, a in nodes)
    body += ''.join(f'<b id="{i}" B="{b}" E="{e}" {a}/>' for i, b, e, a in bonds)
    return ('<?xml version="1.0" encoding="UTF-8" ?>\n'
            '<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd" >\n'
            f'<CDXML id="1"><page id="2"><fragment id="3">{body}</fragment>{extra}</page></CDXML>')


def rules(log):
    return {r.rule for r in log}


# bromochlorofluoromethane, F wedged up from the centre; CDXML points are V3000 units x30, y negated
CHFClBr_V3000 = '''
  chython

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0 0 0 0
M  V30 2 F 0 1.5 0 0
M  V30 3 Cl 1.3 -0.75 0 0
M  V30 4 Br -1.3 -0.75 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2 CFG=1
M  V30 2 1 1 3
M  V30 3 1 1 4
M  V30 END BOND
M  V30 END CTAB
M  END
'''


def chfclbr(sign=1, display='Display="WedgeBegin"', reverse=False):
    nodes = [(4, 0, 0, ''), (5, 0, -45 * sign, 'Element="9"'), (6, 39, 22.5 * sign, 'Element="17"'),
             (7, -39, 22.5 * sign, 'Element="35"')]
    wedge = (8, 5, 4, display) if reverse else (8, 4, 5, display)
    return document(nodes, [wedge, (9, 4, 6, ''), (10, 4, 7, '')])


# ---------------------------------------------------------------------------------------------------
# vocabulary
# ---------------------------------------------------------------------------------------------------

def test_vocabulary_takes_tags_from_the_header_enum():
    table = vocabulary()
    assert table.objects[0x8000] == 'CDXML'
    assert table.objects[0x8004] == 'n' and table.objects[0x8005] == 'b'
    assert table.objects[0x8021] == 'arrow'           # the docs table says 0x8027; real files say 0x8021
    assert table.objects[0x8029] == 'rlogic'
    assert table.properties[0x0600].xml == 'Order'
    assert table.properties[0x0200].type == 'CDXPoint2D'
    assert table.properties[0x1400].xml == 'RLogicGroup'  # header comment says 0x1401; the enum says 0x1400
    assert table.enums['CDXBondDisplay'][6] == 'WedgeBegin'


# ---------------------------------------------------------------------------------------------------
# reading
# ---------------------------------------------------------------------------------------------------

def test_a_drawn_molecule_reads_with_derived_hydrogens():
    m = cdxml(document([(4, 0, 0, ''), (5, 30, -15, ''), (6, 60, 0, 'Element="8"')], [(7, 4, 5, ''), (8, 5, 6, '')]))
    assert m == smiles('CCO')
    assert m.implicit_h_of(3) == 1


def test_coordinates_are_points_over_thirty_with_y_flipped():
    m = cdxml(document([(4, 30, 60, ''), (5, 60, 60, '')], [(6, 4, 5, '')]))
    assert m.xy_of(1) == (1.0, -2.0)


def test_a_wedge_reads_as_the_same_molfile_reads():
    assert cdxml(chfclbr()) == mol(CHFClBr_V3000)


def test_the_mirror_image_drawing_is_the_other_enantiomer():
    assert cdxml(chfclbr(sign=-1)) != cdxml(chfclbr())
    assert cdxml(chfclbr(sign=-1)).parity_of(1)


def test_wedge_end_puts_the_narrow_end_at_e():
    assert cdxml(chfclbr(display='Display="WedgeEnd"', reverse=True)) == cdxml(chfclbr())


def test_a_hashed_wedge_is_the_opposite_configuration():
    assert cdxml(chfclbr(display='Display="WedgedHashBegin"')) == cdxml(chfclbr(sign=-1))


def test_a_nickname_with_a_fragment_is_expanded_and_keeps_its_label():
    inner = ('<fragment id="20"><n id="21" p="60 0"/><n id="22" p="90 0" Element="8"/>'
             '<n id="23" p="60 -30" Element="8"/><n id="24" NodeType="ExternalConnectionPoint"/>'
             '<b id="25" B="21" E="22"/><b id="26" B="21" E="23" Order="2"/><b id="27" B="24" E="21"/>'
             '</fragment><t id="28"><s>COOH</s></t>')
    text = document([(4, 0, 0, ''), (5, 30, 0, '')], [(7, 4, 5, ''), (8, 5, 6, '')])
    text = text.replace('</fragment></page>', f'<n id="6" p="60 0" NodeType="Nickname">{inner}</n>'
                                              '</fragment></page>')
    m = cdxml(text)
    assert m == smiles('CCC(=O)O')
    (group,) = m.sgroups
    assert group['type'] == b'SUP' and (b'LABEL', b'COOH') in group['fields']
    assert sorted(group['atoms']) == [3, 4, 5]


def test_an_r_atom_keeps_its_group_number():
    m = cdxml(document([(4, 0, 0, ''), (5, 30, 0, 'NodeType="GenericNickname" GenericNickname="R2"')],
                       [(6, 4, 5, '')]))
    assert m.atom(2).is_r and m.atom(2).r_index == 2


def test_a_polymer_bead_is_the_marker_aliased_pol():
    m = cdxml(document([(4, 0, 0, ''), (5, 30, 0, 'NodeType="ExternalConnectionPoint" '
                                                  'ExternalConnectionType="PolymerBead"')], [(6, 4, 5, '')]))
    assert m.atom(2).is_r and m.aliases[2] in ('Pol', b'Pol')


@mark.parametrize('attributes, reason', [
    ('NodeType="GenericNickname" GenericNickname="A"', 'any heavy atom'),
    ('NodeType="GenericNickname" GenericNickname="X"', 'any halogen'),
    ('NodeType="ElementList" ElementList="7 8"', 'atom list'),
])
def test_a_query_atom_is_refused_by_name(attributes, reason):
    with raises(UnsupportedCdx, match=reason):
        cdxml(document([(4, 0, 0, ''), (5, 30, 0, attributes)], [(6, 4, 5, '')]))


@mark.parametrize('order', ['1 2', '1 1.5', 'any'])
def test_a_query_bond_is_refused(order):
    with raises(UnsupportedCdx, match='query bond'):
        cdxml(document([(4, 0, 0, ''), (5, 30, 0, '')], [(6, 4, 5, f'Order="{order}"')]))


def test_bond_orders():
    benzene = [(i, 30 * i, 0, '') for i in range(4, 10)]
    ring = [(10 + i, 4 + i, 4 + (i + 1) % 6, 'Order="1.5"') for i in range(6)]
    assert cdxml(document(benzene, ring)).aromatic_bond_count == 6
    m = cdxml(document([(4, 0, 0, 'Element="7"'), (5, 30, 0, 'Element="5"')], [(6, 4, 5, 'Order="dative"')]))
    assert m.bond(1, 2).order == 8
    log = []
    m = cdxml(document([(4, 0, 0, 'Element="8"'), (5, 30, 0, 'Element="1"')], [(6, 4, 5, 'Order="hydrogen"')]),
              log=log)
    assert m.bond(1, 2).order == 8 and 'cdx:bond-order-lost' in rules(log)


def test_charge_isotope_radical_and_stated_hydrogens():
    m = cdxml(document([(4, 0, 0, 'Element="7" Charge="1" NumHydrogens="3"'), (5, 30, 0, 'Isotope="13"'),
                        (6, 60, 0, 'Radical="Doublet"')], [(7, 4, 5, ''), (8, 5, 6, '')]))
    assert m.atom(1).charge == 1 and m.implicit_h_of(1) == 3
    assert m.atom(2).isotope == 13
    assert m.atom(3).is_radical


@mark.parametrize('spelling, kind', [('Absolute', STEREO_ABS), ('Or', STEREO_OR), ('And', STEREO_AND)])
def test_enhanced_stereo_is_read_from_the_node(spelling, kind):
    text = chfclbr().replace('<n id="4" p="0 0" />',
                             f'<n id="4" p="0 0" EnhancedStereoType="{spelling}" EnhancedStereoGroupNum="1"/>')
    m = cdxml(text)
    assert m.stereo_group_of(1) == (kind, 0 if kind == STEREO_ABS else 1)


def test_drawing_objects_are_counted_not_read():
    log = []
    cdxml(document([(4, 0, 0, '')], extra='<t id="9" p="0 0"><s>note</s></t><graphic id="10"/>'), log=log)
    assert 'cdx:drawing-objects' in rules(log)


def test_a_document_without_atoms_is_refused():
    with raises(MalformedCdx):
        cdxml(document([]))


def test_not_cdxml_is_refused():
    with raises(MalformedCdx):
        cdxml('<cml/>')


# ---------------------------------------------------------------------------------------------------
# binary
# ---------------------------------------------------------------------------------------------------

def test_cdx_and_cdxml_are_one_tree():
    text = chfclbr()
    a = parse_cdxml(text, [])
    b = parse_cdx(serialize_cdx(a, []), [])
    ia = {(o.kind, o.id): o.attributes for o in a.walk() if o.id}
    ib = {(o.kind, o.id): o.attributes for o in b.walk() if o.id}
    assert ia.keys() == ib.keys()
    for key, attributes in ia.items():
        for name, value in attributes.items():
            assert [float(t) if t[-1].isdigit() else t for t in ib[key][name].split()] == \
                   [float(t) if t[-1].isdigit() else t for t in value.split()], (key, name)
    assert cdx(serialize_cdx(a, [])) == cdxml(text)


def test_a_binary_point_is_y_then_x():
    raw = HEADER + pack('<Hi', 0x8000, 1) + pack('<Hi', 0x8001, 2) + pack('<Hi', 0x8003, 3)
    raw += pack('<Hi', 0x8004, 4) + pack('<HHii', 0x0200, 8, 60 * 65536, 30 * 65536) + b'\0\0'
    raw += b'\0\0' * 3
    assert cdx(raw).xy_of(1) == (1.0, -2.0)


def test_a_four_byte_charge_is_shifted_as_chemdraw_8_wrote_it():
    raw = HEADER + pack('<Hi', 0x8000, 1) + pack('<Hi', 0x8001, 2) + pack('<Hi', 0x8003, 3)
    raw += pack('<Hi', 0x8004, 4) + pack('<HHh', 0x0402, 2, 7) + pack('<HHi', 0x0421, 4, 1 << 24) + b'\0\0'
    raw += b'\0\0' * 3
    assert cdx(raw).atom(1).charge == 1


def test_bad_magic_is_refused():
    with raises(MalformedCdx):
        cdx(b'not a cdx document at all')


def test_a_truncated_stream_keeps_what_it_read():
    data = cdx(cdxml(chfclbr()))
    log = []
    m = cdx(data[:-8], log=log)
    assert len(m) == 4 and 'cdx:truncated' in rules(log)


def test_cdx_takes_bytes_and_cdxml_takes_text():
    with raises(TypeError):
        cdx('VjCD0100')
    assert issubclass(UnsupportedCdx, CdxError) and issubclass(MalformedCdx, ValueError)


# ---------------------------------------------------------------------------------------------------
# writing
# ---------------------------------------------------------------------------------------------------

@mark.parametrize('writer', [cdx, cdxml])
def test_a_molecule_round_trips(writer):
    m = mol(CHFClBr_V3000)
    assert writer(writer(m)) == m


@mark.parametrize('writer', [cdx, cdxml])
def test_markers_round_trip(writer):
    m = cdxml(document([(4, 0, 0, ''), (5, 30, 0, 'NodeType="GenericNickname" GenericNickname="R3"'),
                        (6, 60, 0, 'NodeType="ExternalConnectionPoint" ExternalConnectionType="PolymerBead"')],
                       [(7, 4, 5, ''), (8, 4, 6, '')]))
    back = writer(writer(m))
    assert back.atom(2).r_index == 3 and back.aliases[3] in ('Pol', b'Pol')


@mark.parametrize('kind', [STEREO_OR, STEREO_AND])
def test_an_enhanced_stereo_group_round_trips(kind):
    spelling = {STEREO_OR: 'Or', STEREO_AND: 'And'}[kind]
    m = cdxml(chfclbr().replace('<n id="4" p="0 0" />',
                                f'<n id="4" p="0 0" EnhancedStereoType="{spelling}" EnhancedStereoGroupNum="2"/>'))
    assert cdxml(cdxml(m)).stereo_group_of(1) == (kind, 2)


def test_writing_without_coordinates_is_logged():
    log = []
    text = cdxml(smiles('CCO'), log=log)
    assert 'cdx:no-coordinates' in rules(log)
    assert cdxml(text) == smiles('CCO')


def test_heteroatoms_state_their_hydrogens():
    text = cdxml(smiles('CN'))
    assert 'NumHydrogens="2"' in text and text.count('NumHydrogens') == 1


# ---------------------------------------------------------------------------------------------------
# reactions
# ---------------------------------------------------------------------------------------------------

def esterification(step_extra='', page_extra=''):
    """Acetic acid + methanol -> methyl acetate, HCl above the arrow, the carbonyl carbon and the
    methanol oxygen mapped."""
    acid = ('<fragment id="10"><n id="11" p="0 0"/><n id="12" p="30 0"/><n id="13" p="45 -26" Element="8"/>'
            '<n id="14" p="45 26" Element="8"/><b id="15" B="11" E="12"/><b id="16" B="12" E="13" Order="2"/>'
            '<b id="17" B="12" E="14"/></fragment>')
    methanol = ('<fragment id="20"><n id="21" p="90 0"/><n id="22" p="120 0" Element="8"/>'
                '<b id="23" B="21" E="22"/></fragment>')
    hcl = '<fragment id="30"><n id="31" p="180 -40" Element="17"/></fragment>'
    ester = ('<fragment id="40"><n id="41" p="300 0"/><n id="42" p="330 0"/><n id="43" p="345 -26" Element="8"/>'
             '<n id="44" p="345 26" Element="8"/><n id="45" p="375 26"/><b id="46" B="41" E="42"/>'
             '<b id="47" B="42" E="43" Order="2"/><b id="48" B="42" E="44"/><b id="49" B="44" E="45"/></fragment>')
    text = '<t id="50" p="170 -10"><s>reflux</s></t>'
    step = (f'<scheme id="60"><step id="61" ReactionStepReactants="10 20" ReactionStepProducts="40" '
            f'ReactionStepObjectsAboveArrow="30 50" ReactionStepAtomMap="12 42 22 44" {step_extra}/></scheme>')
    return ('<?xml version="1.0" encoding="UTF-8" ?>\n<CDXML id="1"><page id="2">'
            f'{acid}{methanol}{hcl}{ester}{text}{page_extra}{step}</page></CDXML>')


def test_a_document_with_a_step_reads_as_a_reaction():
    r = cdxml(esterification())
    assert [str(m) for m in r.reactants] == [str(smiles('CC(=O)O')), str(smiles('CO'))]
    assert [str(m) for m in r.products] == [str(smiles('CC(=O)OC'))]
    assert [m.atom(1).atomic_symbol for m in r.agents] == ['Cl']
    assert r.meta == {'above arrow': 'reflux'}


def test_the_atom_map_pairs_become_map_numbers():
    r = cdxml(esterification())
    acid, methanol = r.reactants
    (ester,) = r.products
    assert acid.atom(2).map_number == ester.atom(2).map_number == 1
    assert methanol.atom(2).map_number == ester.atom(4).map_number == 2
    assert not ester.atom(5).map_number


def test_a_structure_no_side_names_is_logged():
    log = []
    cdxml(esterification(page_extra='<fragment id="70"><n id="71" p="0 200"/></fragment>'), log=log)
    assert 'cdx:structure-outside-reaction' in rules(log)


def test_step_selects_a_reaction_and_out_of_range_raises():
    text = esterification().replace('</scheme>', '<step id="62" ReactionStepReactants="40" '
                                                 'ReactionStepProducts="10"/></scheme>')
    assert [str(m) for m in cdxml(text, step=1).reactants] == [str(smiles('CC(=O)OC'))]
    assert str(cdxml(text, step=-1)) == str(cdxml(text, step=1))
    assert str(cdxml(text)) == str(cdxml(text, step=0))
    with raises(IndexError, match='2 reaction step'):
        cdxml(text, step=2)
    with raises(IndexError, match='no reaction step'):
        cdxml(document([(4, 0, 0, '')]), step=0)


def test_a_second_step_is_logged():
    log = []
    text = esterification().replace('</scheme>', '<step id="62" ReactionStepReactants="40"/></scheme>')
    cdxml(text, log=log)
    assert 'cdx:steps-dropped' in rules(log)


@mark.parametrize('writer', [cdx, cdxml])
def test_a_reaction_round_trips_with_its_map_and_conditions(writer):
    r = cdxml(esterification())
    back = writer(writer(r))
    assert str(back) == str(r) and back.meta == r.meta
    assert back.reactants[1].atom(2).map_number == back.products[0].atom(4).map_number


@mark.parametrize('writer', [cdx, cdxml])
def test_a_salt_participant_stays_one_molecule(writer):
    from chython.core.reaction import ReactionContainer
    r = ReactionContainer([smiles('C[NH3+].[Cl-]')], [smiles('CN')], [])
    back = writer(writer(r))
    assert len(back.reactants) == 1 and len(back.reactants[0].connected_components) == 2


def test_a_dative_bond_to_a_metal_starts_at_the_ligand():
    m = cdxml(document([(4, 0, 0, 'Element="78"'), (5, 30, 0, 'Element="7"')], [(6, 4, 5, 'Order="dative"')]))
    tree = parse_cdxml(cdxml(m), [])
    (bond,) = [o for o in tree.walk() if o.kind == 'b']
    nodes = {o.id: o.attributes.get('Element') for o in tree.walk() if o.kind == 'n'}
    assert nodes[int(bond.attributes['B'])] == '7' and nodes[int(bond.attributes['E'])] == '78'
