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
"""PDBx/mmCIF reader tests, and the STAR tokeniser under it.

``[mutant: ...]`` names the implementation line whose alteration makes the assertion fail.  Fixture
geometry is deliberately unreal: ``mmcif_disulfide.cif`` puts two bonded atoms 3.8 A apart, so a
reader that grew a distance cutoff fails here."""
from pathlib import Path

import pytest

from ..pdb import INAPPLICABLE, UNKNOWN, PDBRecord, is_null, mmcif, parse_star, read_mmcif
from ..pdb._mmcif import _Damage, _charge, _number, _order


_DATA = Path(__file__).resolve().parent.parent.parent.parent / 'test'


def _read(name):
    log = []
    return mmcif((_DATA / name).read_text(encoding='utf-8'), log=log), log


def _blocks(text):
    log = []
    return list(parse_star(text, log=log)), log


def _reported(damage):
    """The lines a ``_Damage`` gives a caller -- the only way its counts are observable."""
    log = []
    damage.report(log)
    return log


# the STAR grammar underneath: `_star.py` knows no chemistry -- no tag name, no category, no element


def test_the_two_nulls_are_distinct_and_neither_is_a_value():
    """``.`` means the item does not apply; ``?`` means its value is not known.

    Both are falsy so that ``value or default`` reads naturally; identity is how a caller that cares
    asks.  [mutant: `_value` returning None for either]
    """
    blocks, log = _blocks('''
data_NULLS
_a.applies   .
_a.unknown   ?
_a.stated    0
_a.quoted    '.'
''')
    block = blocks[0]
    assert block.get('_a.applies') is INAPPLICABLE
    assert block.get('_a.unknown') is UNKNOWN
    assert is_null(block.get('_a.applies')) and is_null(block.get('_a.unknown'))
    assert INAPPLICABLE is not UNKNOWN
    assert not block.get('_a.applies') and not block.get('_a.unknown')
    # `0` is a stated value that happens to be falsy, and a quoted `.` is the string, not the null
    assert block.get('_a.stated') == '0' and not is_null(block.get('_a.stated'))
    assert block.get('_a.quoted') == '.' and not is_null(block.get('_a.quoted'))
    assert not log


def test_a_quote_closes_only_before_whitespace():
    """``O5'`` is a real atom name and ``can't`` is a real word; CIF needs no escape for either.

    A quote closes the value only when the next character is whitespace or end of line, which is what
    lets an apostrophe sit inside a single-quoted value.  [mutant: the `line[j + 1] in ' \\t'` guard]
    """
    blocks, log = _blocks('''
data_QUOTES
_a.word    'can't stop'
_a.bare    O5'
_a.plain   "O5' and more"
_a.two     'first' 'second'
''')
    block = blocks[0]
    assert block.get('_a.word') == "can't stop"
    # a token that merely *contains* a quote needs none: only a leading quote opens a value
    assert block.get('_a.bare') == "O5'"
    assert block.get('_a.plain') == "O5' and more"
    # and the guard does not swallow a following value: two quoted tokens on one line are two tokens
    assert block.get('_a.two') == 'first'
    assert [str(x) for x in log] == ["record: value 'second' on line 6 belongs to no item"]


def test_a_semicolon_text_field_is_only_a_text_field_in_column_one():
    """A ``;`` in column 1 opens a multi-line value; a ``;`` anywhere else is an ordinary character.

    A keyword written inside such a field is data, not syntax.  [mutant: the `line[:1] == ';'` test
    widened to `';' in line`]
    """
    blocks, log = _blocks('''
data_TEXT
_a.long
;first line
loop_ is data here, not a keyword
;
_a.short   a;b
_a.after   stated
''')
    block = blocks[0]
    assert block.get('_a.long') == 'first line\nloop_ is data here, not a keyword'
    assert block.get('_a.short') == 'a;b'
    assert block.get('_a.after') == 'stated'
    assert not block.loops
    assert not log


def test_a_comment_runs_to_the_end_of_its_line_and_not_inside_a_value():
    """``#`` starts a comment, except where it is part of a quoted or multi-line value.

    [mutant: the `ch == '#'` break]
    """
    blocks, log = _blocks('''
# a whole-line comment
data_HASH
_a.one   value    # trailing comment
_a.two   'has # inside'
''')
    block = blocks[0]
    assert block.get('_a.one') == 'value'
    assert block.get('_a.two') == 'has # inside'
    assert not log


def test_cif_two_announces_itself_and_is_read_with_the_older_grammar():
    """The CIF 2.0 magic promises constructs this grammar does not have, so it is named.

    Reading on is right: every construct the two versions share parses identically.  [mutant: the
    magic-comment branch]
    """
    blocks, log = _blocks('#\\#CIF_2.0\ndata_TWO\n_a.b   c\n')
    assert blocks[0].get('_a.b') == 'c'
    assert any(str(m).startswith('unsupported: CIF 2.0 syntax') for m in log), log


def test_a_dictionary_save_frame_is_skipped_rather_than_merged():
    """``save_`` frames hold dictionary definitions, and merging them into a block corrupts it.

    A data file has none, so this is a dictionary handed to the reader by mistake.  [mutant: the
    `save_` branch]
    """
    blocks, log = _blocks('''
data_DICT
_a.real   kept
save__frame_definition
_item.name    '_not.real'
_item.units   angstroms
save_
_a.after   also_kept
''')
    block = blocks[0]
    assert block.get('_a.real') == 'kept' and block.get('_a.after') == 'also_kept'
    assert block.get('_item.name') is None
    assert any(str(m).startswith('unsupported: STAR save frame') for m in log), log


def test_a_line_that_kept_its_carriage_return_does_not_put_one_in_a_value():
    """CRLF and a bare CR both appear in files written on other platforms.

    This lexer splits on space and tab only, so an unstripped ``\\r`` becomes part of the last bare
    token on every line.  The path that reaches it is an iterable of lines that kept their endings.
    [mutant: the `rstrip('\\r')` in `_tokens`]
    """
    text = 'data_CR\r\n_a.b   value\r\n_a.c   other\r\n'
    for source in (text, text.split('\n')):
        blocks = list(parse_star(source))
        assert blocks[0].name == 'CR'
        assert blocks[0].get('_a.b') == 'value'
        assert blocks[0].get('_a.c') == 'other'


def test_multiple_data_blocks_come_back_separately():
    """One CIF file can hold many blocks.  [mutant: `yield block` on `data_`]"""
    blocks, log = _blocks('data_ONE\n_a.b  1\ndata_TWO\n_a.b  2\n')
    assert [b.name for b in blocks] == ['ONE', 'TWO']
    assert [b.get('_a.b') for b in blocks] == ['1', '2']
    assert not log


# the shape of the answer


def test_a_record_is_not_a_container():
    """The reader stops at a neutral record: three coordinates and every annotation, no chemistry.

    ``MoleculeContainer`` has nowhere to put z, a residue name or an alt_loc, so building one here
    would be lossy.  [mutant: `_split_models` returning containers]
    """
    records, _ = _read('mmcif_ligand_water.cif')
    assert len(records) == 1
    assert isinstance(records[0], PDBRecord)
    assert [(a.x, a.y, a.z) for a in records[0].atoms][3] == (8.0, 8.0, 8.0)
    assert records[0].atoms[3].residue_name == 'HOH'


def test_annotations_survive_in_both_numberings():
    """label_* and auth_* are two parallel identifier sets and both are stored.

    The label asym id is a *string*, not a single character, and it is not the auth chain: this
    fixture has label ``A``/auth ``B`` for the ligand, label ``C``/auth ``W`` for the water.
    [mutant: `auth_asym_id` read into `chain`]
    """
    records, _ = _read('mmcif_ligand_water.cif')
    ligand, water = records[0].atoms[0], records[0].atoms[3]
    assert (ligand.chain, ligand.auth_chain, ligand.auth_seq) == ('A', 'B', 501)
    assert (water.chain, water.auth_chain, water.auth_seq) == ('C', 'W', 601)


def test_entity_type_classifies_water_and_ligand():
    """A residue is identifiable as water, polymer or ligand from `_entity.type`, not a name list.

    [mutant: `_entity_types` returning `{}`]
    """
    records, _ = _read('mmcif_ligand_water.cif')
    assert [a.is_water for a in records[0].atoms] == [False, False, False, True]
    assert records[0].atoms[0].is_ligand
    assert not records[0].atoms[0].is_polymer


# the element


def test_element_comes_from_type_symbol_and_never_from_the_atom_name():
    """A haem pyrrole nitrogen is named ``NA`` and is nitrogen.

    ``type_symbol`` is the only source of the element; the atom name is not one, since reading it
    turns four nitrogens per haem into sodium.  [mutant: `type_symbol` falling back to
    `label_atom_id`]
    """
    records = mmcif('''
data_HEM
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
HETATM 1 N  NA HEM A 0.000 0.000 0.000
HETATM 2 N  NB HEM A 1.000 0.000 0.000
HETATM 3 N  NC HEM A 2.000 0.000 0.000
HETATM 4 N  ND HEM A 3.000 0.000 0.000
HETATM 5 FE FE HEM A 1.500 1.500 0.000
''')
    assert [a.element for a in records[0].atoms] == ['N', 'N', 'N', 'N', 'Fe']
    assert [a.atom_name for a in records[0].atoms] == ['NA', 'NB', 'NC', 'ND', 'FE']


def test_an_atom_with_no_element_is_stored_without_one_and_logged():
    """``?`` in the element column is unknown, not carbon.

    The atom keeps its coordinates and annotations and loses only the fact the file does not state.
    [mutant: the `atom: no element stated` branch]
    """
    records, log = _read('mmcif_damaged.cif')
    assert records[0].atoms[1].element is None
    assert records[0].atoms[1].atom_name == 'C2'
    assert any(str(m).startswith('atom: no element stated') for m in log), log


def test_pdbx_formal_charge_reads_both_spellings():
    """The dictionary spells it ``-2``; writers also spell it ``2-``.  [mutant: the `_charge` swap]"""
    damage = _Damage()
    assert _charge('-2', 'x', damage) == -2
    assert _charge('2-', 'x', damage) == -2
    assert _charge('2+', 'x', damage) == 2
    assert _charge('+2', 'x', damage) == 2
    assert not _reported(damage)
    assert _charge('two', 'x', damage) == 0
    assert str(_reported(damage)[0]).startswith('atom: formal charge')


def test_a_number_carrying_an_esd_is_read_and_the_uncertainty_named():
    """``1.234(5)`` is a value with an estimated standard deviation, common in small-molecule CIF.

    The value is read and the esd is named.  [mutant: the esd branch of `_number`]
    """
    damage = _Damage()
    assert _number('1.234(5)', 'x coordinate', 'atom 1', damage) == 1.234
    assert 'carries an uncertainty' in _reported(damage)[0]


# bonds the file states


def test_chem_comp_bond_bonds_every_instance_of_its_component():
    """One ``_chem_comp_bond`` block bonds the ligand, and its double bond stays double.

    Counting bonds alone is not enough: a reader reading every order as single gives the same count.
    [mutant: `_order` returning 1 for DOUB]
    """
    records, _ = _read('mmcif_disulfide.cif')
    record = records[0]
    assert len(record.bonds) == 11
    doubles = [b for b in record.bonds if b.order == 2]
    assert len(doubles) == 2
    assert {record.atoms[b.a].atom_name for b in doubles} == {'C'}
    assert {record.atoms[b.b].atom_name for b in doubles} == {'O'}
    assert all(b.stated_order for b in record.bonds)


def test_a_stated_bond_outranks_any_distance():
    """The disulfide is built because ``_struct_conn`` states it, and for no other reason.

    CB and SG of the second monomer are 3.8 A apart, past every cutoff any toolkit uses, and that
    bond is built too.  [mutant: any distance test added to `_add_struct_conn`]
    """
    records, _ = _read('mmcif_disulfide.cif')
    record = records[0]
    disulfide = [b for b in record.bonds if b.conn_type == 'disulf']
    assert len(disulfide) == 1
    assert {record.atoms[disulfide[0].a].atom_name,
            record.atoms[disulfide[0].b].atom_name} == {'SG'}
    assert record.atoms[disulfide[0].a].residue_seq != record.atoms[disulfide[0].b].residue_seq
    far = next(b for b in record.bonds
               if {record.atoms[b.a].atom_name, record.atoms[b.b].atom_name} == {'CB', 'SG'}
               and record.atoms[b.a].residue_seq == 2)
    one, other = record.atoms[far.a], record.atoms[far.b]
    assert (one.x - other.x) ** 2 + (one.y - other.y) ** 2 + (one.z - other.z) ** 2 > 9.


def test_a_file_stating_no_connectivity_yields_its_atoms_and_says_so():
    """Zero bonds is an answer.  Silence about zero bonds is the defect.

    The record comes back with its atoms and the log names the unbonded count, which is what tells a
    caller the *file* stated no connectivity.  [mutant: the `bond: the file states no connectivity`
    branch]
    """
    records, log = _read('mmcif_no_bonds.cif')
    assert len(records) == 1
    assert len(records[0].atoms) == 3 and not records[0].bonds
    assert any('bond: the file states no connectivity for this record; 3 atom(s) are unbonded' in m
               for m in log), log


def test_a_partly_bonded_record_names_the_atoms_no_bond_touches():
    """The water in this fixture is bonded to nothing, and one line says exactly that.

    Separate from the no-connectivity case: a file with some bonds and an unbonded atom is where an
    unreported gap hides.  [mutant: the `touched by no stated bond` branch]
    """
    records, log = _read('mmcif_ligand_water.cif')
    assert records[0].unbonded_count() == 1
    assert any('bond: 1 atom(s) are touched by no stated bond' in m for m in log), log


def test_alternate_conformers_never_bond_to_each_other():
    """``CB`` of conformer A bonds ``OG`` of conformer A, and never ``OG`` of conformer B.

    ``_chem_comp_bond`` names ``CB-OG`` once, the model holds two of each, and to a name match the
    cross pairs are as plausible as the right ones.  [mutant: `_alt_compatible` returning True]
    """
    records, _ = _read('mmcif_altloc.cif')
    record = records[0]
    pairs = {(record.atoms[b.a].atom_name, record.atoms[b.a].alt_loc,
              record.atoms[b.b].atom_name, record.atoms[b.b].alt_loc) for b in record.bonds}
    assert ('CB', 'A', 'OG', 'A') in pairs
    assert ('CB', 'B', 'OG', 'B') in pairs
    assert ('CB', 'A', 'OG', 'B') not in pairs
    assert ('CB', 'B', 'OG', 'A') not in pairs
    # a conformer-free atom bonds both conformers: the rule is "compatible", not "equal"
    assert sum(1 for b in record.bonds if record.atoms[b.a].atom_name == 'CA'
               and record.atoms[b.b].atom_name == 'CB') == 2


def test_the_aromatic_flag_is_applied_when_it_is_the_only_order_stated():
    """``pdbx_aromatic_flag Y`` with no ``value_order`` is an aromatic bond, not a single one.

    Where both are stated the Kekule ``value_order`` is more information and wins; this row has only
    the flag.  [mutant: the `order, stated = 4, True` branch]
    """
    records = mmcif('''
data_AROM
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.pdbx_aromatic_flag
LIG C1 C2 Y
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
HETATM 1 C C1 LIG A 0.000 0.000 0.000
HETATM 2 C C2 LIG A 1.390 0.000 0.000
''')
    assert [b.order for b in records[0].bonds] == [4]
    assert records[0].bonds[0].stated_order


def test_an_order_chython_cannot_hold_is_stored_single_and_named():
    """``QUAD`` has no chython order, so the bond is a single one and the approximation is a log line.

    ``stated_order`` carries the same fact on the bond, which is how a caller tells this single bond
    from a stated one.  [mutant: the `_UNHELD_ORDERS` branch of `_order`]
    """
    records, log = _read('mmcif_damaged.cif')
    approximated = [b for b in records[0].bonds if not b.stated_order]
    assert len(approximated) == 1 and approximated[0].order == 1
    assert any('quadruple bond order' in m and str(m).startswith('unsupported: ') for m in log), log


def test_a_partner_under_a_non_identity_symmetry_operator_gets_no_bond():
    """The atom the row names is a symmetry image, and images are not among the deposited coordinates.

    Bonding to the copy in the asymmetric unit would join two atoms that are not neighbours, so the
    construct is named instead.  [mutant: the `_IDENTITY_SYMMETRY` comparison]
    """
    records, log = _read('mmcif_damaged.cif')
    assert not any(b.conn_type == 'metalc' for b in records[0].bonds)
    assert any('symmetry operator 2_555' in m and str(m).startswith('unsupported: ') for m in log), log


def test_metal_coordination_is_a_dative_bond_with_no_direction_claimed():
    """``metalc`` is order 8, and the log says the file states no donor direction.

    ``_struct_conn`` names two partners and does not say which donates, so the reader stores them in
    the order the file names them.  Choosing a direction needs a periodic table and belongs to the
    container-build pass.  [mutant: `order, stated = 8, True` for metalc, or the `metal` counter's
    log line]
    """
    log = []
    records = mmcif('''
data_ZN
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_auth_seq_id
_struct_conn.ptnr1_auth_asym_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_auth_seq_id
_struct_conn.ptnr2_auth_asym_id
_struct_conn.ptnr2_label_atom_id
m1 metalc A CYS 10 A SG B ZN 300 B ZN
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM   1 S  SG CYS A A 10  3.750 -0.900 1.300
HETATM 2 ZN ZN ZN  B B 300 5.500 -1.200 2.000
''', log=log)
    bond = records[0].bonds[0]
    assert bond.order == 8 and bond.conn_type == 'metalc'
    assert (records[0].atoms[bond.a].atom_name, records[0].atoms[bond.b].atom_name) == ('SG', 'ZN')
    assert any('states no donor direction' in m for m in log), log


def test_a_hydrogen_bond_is_not_a_bond():
    """``hydrog`` is a stated contact chython has no order for, so it is named and not built.

    Building it as a single bond would put a hydrogen bond in the connectivity table.  [mutant:
    `_NON_COVALENT_CONN` membership test]
    """
    records, log = _read('mmcif_damaged.cif')
    assert not any(b.conn_type == 'hydrog' for b in records[0].bonds)
    assert any('hydrogen bond' in m and str(m).startswith('unsupported: ') for m in log), log


def test_a_struct_conn_partner_that_is_not_in_the_model_is_named():
    """A ``_struct_conn`` row naming an absent atom loses its bond, and the row is named.

    [mutant: the `names an atom that is not in this model` branch]
    """
    records, log = _read('mmcif_damaged.cif')
    assert any(str(m).startswith('bond: _struct_conn row 1 names an atom') for m in log), log


def test_a_struct_conn_row_with_no_sequence_number_says_that_and_not_something_else():
    """Two ways a partner comes back empty, and the log tells them apart.

    A sequence number that misses points outside the model; a row stating none names no atom to search
    for, and the atoms it meant may well be present.  [mutant: the `searched` flag in `_partner`]
    """
    log = []
    records = mmcif('''
data_NOSEQ
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_atom_id
c1 covale A LIG C1 A LIG C2
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
HETATM 1 C C1 LIG A 0.000 0.000 0.000
HETATM 2 C C2 LIG A 1.500 0.000 0.000
''', log=log)
    assert not records[0].bonds
    assert any('states neither a label nor an auth sequence number for partner 1' in m for m in log), \
        log


def test_the_polymer_backbone_linkage_is_a_named_boundary():
    """mmCIF states the peptide bond as a *sequence*, not as a pair of atoms, so it is not built.

    Building it needs the per-component attachment points -- the leaving atom of one monomer and the
    entry atom of the next -- a component-dictionary fact this reader does not have.  Each monomer
    comes back bonded within itself and the gap is one log line.  [mutant:
    `_report_polymer_linkage`]
    """
    records, log = _read('mmcif_disulfide.cif')
    record = records[0]
    assert not any({record.atoms[b.a].atom_name, record.atoms[b.b].atom_name} == {'C', 'N'}
                   for b in record.bonds)
    assert any(str(m).startswith('unsupported: the file states its polymer linkage as a sequence')
               and '2 polymer monomer(s)' in m for m in log), log


# models and blocks


def test_each_model_becomes_its_own_record():
    """An NMR ensemble is N records, and neither model is the other's conformer.

    Merging them puts every atom in twice at two positions and bonds them to each other.  [mutant:
    `_split_models` keying on nothing]
    """
    records, log = _read('mmcif_two_models.cif')
    assert [r.model for r in records] == [1, 2]
    assert [len(r.atoms) for r in records] == [2, 2]
    assert records[0].atoms[0].x != records[1].atoms[0].x
    assert sum(1 for m in log if str(m).startswith('record: the block states 2 models')) == 2


def test_a_block_level_message_reaches_the_caller_once_and_every_record_always():
    """A record has to be readable on its own; the caller's flat log must not repeat a block-level line.

    The two halves pull opposite ways: the per-model line appears once per record in the caller's log,
    while the block-level line appears once there but in *both* records' own logs.  [mutant:
    `log.extend(block_log)` inside the record loop]
    """
    records, log = _read('mmcif_two_models.cif')
    assert sum(1 for m in log if str(m).startswith('record: the block states 2 models')) == 2
    for record in records:
        assert any(str(m).startswith('record: the block states 2 models') for m in record.log)
    records, log = _read('mmcif_ligand_water.cif')
    assert sum(1 for m in log if str(m).startswith('unsupported: 1 mmCIF category')) == 1
    assert any(str(m).startswith('unsupported: 1 mmCIF category') for m in records[0].log)


def test_an_unread_category_is_named_with_its_own_name():
    """``_cell`` is in the file and not in the reader, so it is named.

    The list names what was in the file rather than what the reader happens to support.  [mutant:
    `_report_categories`]
    """
    _, log = _read('mmcif_ligand_water.cif')
    assert any(str(m) == 'unsupported: 1 mmCIF category(ies) not modelled: _cell' for m in log), log


def test_an_unread_atom_site_item_is_named_too():
    """One level finer: an ``_atom_site`` item the reader does not read is named individually.

    ``_atom_site`` is a category the reader claims to read, so an unread item inside it is not covered
    by the category line.  [mutant: `_report_atom_site_items`]
    """
    _, log = _read('mmcif_damaged.cif')
    assert any(str(m) == 'unsupported: 1 _atom_site item(s) not modelled: _atom_site.cartn_x_esd'
               for m in log), log


def test_a_single_atom_written_as_items_rather_than_a_loop_is_read():
    """STAR lets a one-row category be scalar items, and small files use it.

    A reader that only understands ``loop_`` reads no atoms at all here.  [mutant: the scalar branch
    of `_atom_site_rows`]
    """
    records = mmcif('''
data_ONE
_atom_site.group_PDB     HETATM
_atom_site.id            1
_atom_site.type_symbol   O
_atom_site.label_atom_id O
_atom_site.label_comp_id HOH
_atom_site.label_asym_id A
_atom_site.Cartn_x       1.000
_atom_site.Cartn_y       2.000
_atom_site.Cartn_z       3.000
''')
    assert len(records) == 1 and len(records[0].atoms) == 1
    assert (records[0].atoms[0].element, records[0].atoms[0].z) == ('O', 3.0)


def test_reading_from_a_path_matches_reading_from_text():
    """`read_mmcif` takes a path and yields lazily; `mmcif` takes the text.  [mutant: `_iter_lines`]"""
    log = []
    records = list(read_mmcif(_DATA / 'mmcif_ligand_water.cif', log=log))
    text_records, text_log = _read('mmcif_ligand_water.cif')
    assert len(records) == len(text_records)
    assert [a.atom_name for a in records[0].atoms] == [a.atom_name for a in text_records[0].atoms]
    assert log == text_log


# damage, none of it fatal


def test_every_malformation_in_the_damaged_fixture_is_read_and_logged():
    """The input posture in one file: nothing is refused, everything is named.

    Asserted as prefixes rather than a count, so a new log line does not break it.  [mutant: any
    single damage branch]
    """
    records, log = _read('mmcif_damaged.cif')
    assert len(records) == 1
    record = records[0]
    assert len(record.atoms) == 4 and len(record.bonds) == 2
    # the atom whose x coordinate is the word 'middle' keeps y, z and every annotation
    assert (record.atoms[2].x, record.atoms[2].y, record.atoms[2].z) == (None, 1.2, 0.0)
    assert record.atoms[2].occupancy == -0.5
    # a 5-character component id is not damage: wwPDB issues them and no mmCIF field has a width
    assert record.atoms[0].residue_name == 'LIGXY'
    for expected in ('record: duplicate item _entry.id',
                     'record: item _struct.title on line 14 has no value',
                     'record: loop_ over _atom_site states 15 column(s)',
                     'atom: no element stated',
                     'atom: occupancy -0.5',
                     "atom: x coordinate 'middle'",
                     'unsupported: 1 _atom_site item(s) not modelled',
                     'unsupported: quadruple bond order',
                     'bond: _struct_conn row 1 names an atom',
                     'unsupported: 1 _struct_conn row(s) state a hydrogen bond',
                     'unsupported: 1 _struct_conn row(s) join an atom under symmetry operator',
                     'bond: 1 atom(s) are touched by no stated bond'):
        assert any(str(m).startswith(expected) for m in log), (expected, log)


def test_a_duplicated_item_keeps_the_first_value():
    """Two ``_entry.id`` values: the first wins and the duplicate is named.  [mutant: the
    `duplicate item` branch]
    """
    records, log = _read('mmcif_damaged.cif')
    assert records[0].entry_id == 'TESTDAMAGE'
    assert any(str(m).startswith('record: duplicate item _entry.id on line 13') for m in log), log


def test_a_short_final_loop_row_is_padded_and_the_row_named():
    """A ``loop_`` row short of its header is padded with unknown, not dropped.

    The iron here has no occupancy and no auth numbering, and is still an atom with coordinates.
    [mutant: the padding branch]
    """
    records, log = _read('mmcif_damaged.cif')
    iron = records[0].atoms[3]
    assert (iron.element, iron.x) == ('Fe', 5.0)
    assert iron.occupancy is None and iron.auth_seq is None
    assert any('its last row holds 11, padded to width with unknown' in m for m in log), log


def test_an_unterminated_text_field_runs_to_the_end_of_the_file_and_says_so():
    """A lost closing ``;`` swallows the rest of the file, which is what the grammar says happens.

    There is nothing to resynchronise on, so the reader keeps the atom that came before and names the
    line the field opened on.  [mutant: the unterminated-field branch of `_tokens`]
    """
    records, log = _read('mmcif_unterminated_text.cif')
    assert len(records) == 1 and len(records[0].atoms) == 1
    assert 'HETATM 2' in records[0].title
    assert any(str(m).startswith('record: multi-line text field opened on line 25 is never closed')
               for m in log), log


@pytest.mark.parametrize('value,order,stated,named', [
    ('SING', 1, True, False), ('DOUB', 2, True, False), ('TRIP', 3, True, False),
    ('AROM', 4, True, False), ('sing', 1, True, False),
    ('QUAD', 1, False, True), ('POLY', 1, False, True), ('DELO', 1, False, True),
    ('PI', 1, False, True), ('NONSENSE', 1, False, True),
    (None, 1, False, False)])
def test_value_order_spellings(value, order, stated, named):
    """The four orders chython holds, the four it does not, an unknown token, and nothing at all.

    ``stated`` distinguishes a single bond the file stated from one the reader fell back to.  ``None``
    is the one case with no line of its own: a field stating nothing is counted by the caller and
    reported once per table.  [mutant: `_ORDERS`, `_UNHELD_ORDERS`]
    """
    damage = _Damage()
    assert _order(value, damage, 'a row') == (order, stated)
    assert bool(_reported(damage)) is named


def test_a_bond_whose_order_the_file_omits_is_counted_once_for_the_table():
    """The aggregate line `_order` leaves to its caller, on a row that states no order at all.

    [mutant: the `unstated` counter in `_component_bonds`]
    """
    log = []
    records = mmcif('''
data_NOORDER
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
LIG C1 C2
LIG C2 C3
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
HETATM 1 C C1 LIG A 0.000 0.000 0.000
HETATM 2 C C2 LIG A 1.500 0.000 0.000
HETATM 3 C C3 LIG A 3.000 0.000 0.000
''', log=log)
    assert [b.order for b in records[0].bonds] == [1, 1]
    assert not any(b.stated_order for b in records[0].bonds)
    assert any(str(m) == 'bond: 2 _chem_comp_bond row(s) state no bond order at all; each is stored as a '
               'single bond' for m in log), log


# every connection row, answered


def test_a_row_naming_a_conformer_the_model_does_not_have_says_so():
    """A third way a partner comes back empty: the atom is in the model, the conformer is not.

    Neither of the other two reasons is true about it, so a caller told "names an atom that is not in
    this model" would go looking for an atom that is there.  Asserted as the exact line, since the
    failure mode is a reason of ``None`` interpolated into the message.  [mutant: the
    narrowed-to-nothing branch of `_partner`]
    """
    log = []
    records = mmcif('''
data_ALTGONE
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.pdbx_ptnr1_label_alt_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id
c1 covale A LIG 1 C1 B A LIG 1 C2 .
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
HETATM 1 C C1 A LIG A 1 0.000 0.000 0.000
HETATM 2 C C2 . LIG A 1 1.500 0.000 0.000
''', log=log)
    assert not records[0].bonds
    assert any(str(m) == "bond: _struct_conn row 1 names alternate conformer 'B' of partner 1, which this "
               'model does not contain; no bond built' for m in log), log
    assert not any('None' in m for m in log), log


@pytest.mark.parametrize('kind', ['covale', 'covale_base', 'covale_phosphate', 'covale_sugar',
                                  'modres'])
def test_every_covalent_connection_type_reads_as_a_covalent_bond(kind):
    """Five ordinary ``conn_type_id`` values; four are how a nucleic acid and a modified residue
    spell a covalent link.

    The row's own spelling stays on the bond, because a caller filtering for a modified-residue link
    needs to see ``modres`` and not ``covale``.  [mutant: `_COVALENT_CONN` membership]
    """
    log = []
    records = mmcif(f'''
data_COVALENT
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id
c1 {kind} A LIG 1 C1 A LIG 2 C2
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
HETATM 1 C C1 LIG A 1 0.000 0.000 0.000
HETATM 2 C C2 LIG A 2 1.500 0.000 0.000
''', log=log)
    bond, = records[0].bonds
    assert (bond.order, bond.stated_order, bond.source) == (1, False, 'struct_conn')
    assert bond.conn_type == kind
    assert not any('not a connection type this reader knows' in m for m in log), log
    # Counted with every other order-less covalent row, which is the line the arithmetic below checks.
    assert any(str(m) == 'bond: 1 _struct_conn covalent link(s) state no bond order; each is stored as a '
               'single bond' for m in log), log


def test_a_row_joining_an_atom_to_itself_is_named_and_counted_as_no_bond():
    """A row whose two partners resolve to one atom builds nothing, so nothing may be reported built.

    Counting the order before discarding the pair makes the aggregate line claim a bond the record
    does not hold.  Both halves are asserted: the line the discard earns and the absence of the one it
    must not.  A broken file read anyway carries no ``unsupported:``.  [mutant: the order counted
    before the discard]
    """
    log = []
    records = mmcif('''
data_SELF
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id
c1 covale A LIG 1 C1 A LIG 1 C1
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
HETATM 1 C C1 LIG A 1 0.000 0.000 0.000
''', log=log)
    assert not records[0].bonds
    assert any(str(m) == 'bond: _struct_conn row 1 joins an atom to itself; no bond built' for m in log), log
    assert not any('covalent link(s) state no bond order' in m for m in log), log


# the log a large file gives back


def test_damage_of_one_kind_is_one_line_however_many_rows_carry_it():
    """A column a writer left out is missing from every row, and 6000 rows are one broken writer.

    The line is the first message with the rest counted onto it: the count is what a caller acts on
    and the first row is what points into the file.  [mutant: `_Damage.report`]
    """
    log = []
    records = mmcif('''
data_NOELEMENT
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
HETATM 1 C1 LIG A 0.000 0.000 0.000
HETATM 2 C2 LIG A 1.500 0.000 0.000
HETATM 3 C3 LIG A 3.000 0.000 0.000
''', log=log)
    assert [a.element for a in records[0].atoms] == [None, None, None]
    element = [m for m in log if 'no element stated' in m]
    assert [str(x) for x in element] == [
        'atom: no element stated on atom 1; element not stored (and 2 more row(s))'], log


def test_a_file_with_no_data_block_is_told_that_once():
    """Text that is not CIF, and the empty string, both get one answer and it is structural.

    "No ``data_`` block anywhere" is a property of the parse and checkable; how CIF-like the text
    looked is a guess.  [mutant: the `saw_data_block` test, the orphan counter]
    """
    for text in ('', 'This file is prose and not CIF at all.\nIt goes on for a second line.'):
        log = []
        assert mmcif(text, log=log) == []
        assert any(str(m) == 'record: no data_ block is stated anywhere in the text, so there is no CIF '
                   'here; no block is read' for m in log), (text, log)
        orphans = [m for m in log if 'belongs to no item' in m]
        assert len(orphans) <= 1, log                   # one line for the lot, or none for the empty
    assert orphans and str(orphans[0]).endswith('(and 15 more value(s))'), orphans


def test_the_unread_author_identifiers_are_named_like_any_other_unread_item():
    """Two of the four ``auth_*`` items are stored per atom and two are not, so two are reported.

    ``_ATOM_SITE_READ`` is a drift detector: an item listed in it and read by nothing points the
    detector the wrong way.  [mutant: `auth_comp_id`/`auth_atom_id` back inside `_ATOM_SITE_READ`]
    """
    records, log = _read('mmcif_ligand_water.cif')
    assert any(str(m) == 'unsupported: 2 _atom_site item(s) not modelled: _atom_site.auth_atom_id, '
               '_atom_site.auth_comp_id' for m in log), log
    # the two that *are* stored, so the line above is a statement about which four
    assert records[0].atoms[0].auth_chain is not None and records[0].atoms[0].auth_seq is not None


# what the file states, and only that


def test_a_polymer_is_recognised_by_its_sequence_numbers_when_the_file_states_no_entity():
    """``_entity`` is not mandatory, so ``label_seq_id`` is the fallback.

    ``label_seq_id`` is dictionary-defined for a polymer entity and null for everything else, so a
    non-null one is the file stating the same fact in its other place.  [mutant:
    `_is_polymer_monomer` reading `entity_type` only]
    """
    log = []
    records = mmcif('''
data_TWOGLY
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 N N  GLY A 1 0.000 0.000 0.000
ATOM 2 C CA GLY A 1 1.450 0.000 0.000
ATOM 3 N N  GLY A 2 3.000 0.000 0.000
ATOM 4 C CA GLY A 2 4.450 0.000 0.000
''', log=log)
    assert [a.entity_type for a in records[0].atoms] == [None] * 4
    assert any(str(m).startswith('unsupported: the file states its polymer linkage as a sequence')
               and '2 polymer monomer(s)' in m for m in log), log


def test_a_water_is_not_a_polymer_monomer_when_the_file_states_no_entity():
    """The pair of the test above: a null ``label_seq_id`` is the file saying "not a polymer".

    A fallback counting every residue reports a polymer linkage for a box of waters.  [mutant:
    `_is_polymer_monomer` returning True unconditionally]
    """
    log = []
    records = mmcif('''
data_WATERS
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
HETATM 1 O O HOH A . 1 0.000 0.000 0.000
HETATM 2 O O HOH A . 2 3.000 0.000 0.000
''', log=log)
    assert len(records[0].atoms) == 2
    assert not any('polymer linkage' in m for m in log), log


def test_an_insertion_code_is_read_and_keeps_two_residues_apart():
    """``pdbx_PDB_ins_code`` distinguishes 100 from 100A, and it is part of a residue's key.

    Not reading it merges the two residues an insertion code exists to separate, and then applies one
    component's bond table across both.  [mutant: the ``ins_code`` read]
    """
    log = []
    records = mmcif('''
data_INSCODE
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 N N GLY A 100 . 0.000 0.000 0.000
ATOM 2 N N GLY A 100 A 3.000 0.000 0.000
''', log=log)
    first, second = records[0].atoms
    assert (first.ins_code, second.ins_code) == (None, 'A')
    assert first.residue_key != second.residue_key


def test_a_distance_decides_nothing_in_either_direction():
    """Two atoms half an Angstrom apart with no bond stated, and two 500 Angstroms apart with one.

    Every other fixture states the bonds its geometry implies, so only this one fails the moment a
    coordinate reaches a bond decision -- as a cutoff that adds a bond or a check that drops one.
    There is no covalent radius and no distance function in this package.  [mutant: any distance test
    inside the bond builders]
    """
    log = []
    records = mmcif('''
data_GEOMETRY
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
LIG C1 C3 SING
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
HETATM 1 C C1 LIG A 0.000 0.000 0.000
HETATM 2 C C2 LIG A 0.500 0.000 0.000
HETATM 3 C C3 LIG A 500.000 0.000 0.000
''', log=log)
    assert {bond.key for bond in records[0].bonds} == {(0, 2)}
    assert any(str(m).startswith('bond: 1 atom(s) are touched by no stated bond') for m in log), log
