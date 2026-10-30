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
from chython.chemistry import audit_stereo
from chython.chemistry._stereo_audit import _PASS_THROUGH, _apply_fix, _audit
from chython.chemistry._tables import read_table, stereo_audit_codes
from chython.core import MoleculeContainer, read_smiles


def canonical(m):
    """`m` (a molecule or a SMILES) canonicalized on a copy: the meaning the audit compares, by the explicit
    rewrite, since `==` and `str()` follow the stored labels."""
    m = (read_smiles(m) if isinstance(m, str) else m).copy()
    m.canonicalize_stereo_groups()
    return m
from pathlib import Path
from pytest import mark
from re import DOTALL, MULTILINE, findall, search


def codes(text):
    return sorted(f.code for f in audit_stereo(read_smiles(text)).findings)


def test_the_table_names_every_spec_code():
    assert set(stereo_audit_codes()) == {
        'unknown-unit', 'unlabelled-sign', 'dead-unconfigured', 'dead-not-stereogenic', 'dead-waits-for',
        'dead-restating-group', 'and-single-member-racemate', 'and-split', 'and-absorbed',
        'bridgehead-unrealisable', 'fused-strained', 'bridgehead-implied', 'bridgehead-redundant-label',
        'bridgehead-half-labelled', 'bridgehead-merged-groups', 'meso-absolute', 'identity-approximate',
        'truncated'}


def test_empty_and_achiral_molecules_are_clean():
    assert audit_stereo(MoleculeContainer()).clean and not audit_stereo(MoleculeContainer()).findings
    assert audit_stereo(read_smiles('CCO')).findings == ()


def test_an_unsigned_stereocentre_is_unknown():
    report = audit_stereo(read_smiles('CC(O)CC'))                  # butan-2-ol, CID 6568
    [f] = report.findings
    assert f.code == 'stereo-audit:unknown-unit' and f.atoms == (2,) and report.clean
    assert sorted(x.args for x in f.fixes) == [(1,), (2,)]
    assert all(x.action == 'set-sign' and not x.preserves_meaning for x in f.fixes)
    assert {read_smiles(x.result) for x in f.fixes} == {read_smiles('C[C@H](O)CC'), read_smiles('C[C@@H](O)CC')}


GLUCOSE = 'C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O'        # D-glucopyranose, CID 5793: anomeric C unsigned
BETA_GLUCOSE = 'C([C@@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O)O'  # beta-D-glucopyranose, CID 64689


def test_glucose_leaves_its_anomeric_centre_unknown():
    [f] = audit_stereo(read_smiles(GLUCOSE)).findings
    assert f.code == 'stereo-audit:unknown-unit' and f.atoms == (6,)
    assert audit_stereo(read_smiles(BETA_GLUCOSE)).findings == ()


def test_a_label_without_a_sign_is_dead_unconfigured():
    [f] = [f for f in audit_stereo(read_smiles('CC(O)CC |&1:1|')).findings if f.severity == 'warning']
    assert f.code == 'stereo-audit:dead-unconfigured' and f.group == (3, 1)
    [fix] = f.fixes
    assert fix.action == 'drop-label' and fix.preserves_meaning
    assert canonical(fix.result) == canonical('CC(O)CC |&1:1|')


def test_an_unsigned_double_bond_is_unknown():
    # but-2-ene, CID 12220: cis/trans axis, owner pair reported in atoms
    [f] = audit_stereo(read_smiles('CC=CC')).findings
    assert f.code == 'stereo-audit:unknown-unit' and isinstance(f.atoms[0], tuple)
    assert len(f.fixes) == 2
    results = [read_smiles(x.result) for x in f.fixes]
    assert results[0] != results[1]


def test_an_unsigned_allene_is_unknown():
    # penta-2,3-diene: axial stereocentre; owner pair in atoms, no crash on set_parity
    [f] = audit_stereo(read_smiles('CC=C=CC')).findings
    assert f.code == 'stereo-audit:unknown-unit' and isinstance(f.atoms[0], tuple)
    assert len(f.fixes) == 2


def test_a_sign_on_a_non_stereogenic_unit_is_dead():
    [f] = audit_stereo(read_smiles('C[C@H](C)CC |&1:1|')).findings
    assert f.code == 'stereo-audit:dead-not-stereogenic'
    [fix] = f.fixes
    assert fix.action == 'drop-label' and fix.preserves_meaning


def test_equal_readings_are_one_fix():
    [f] = [f for f in audit_stereo(read_smiles('CC(O)CC |&1:1|')).findings if f.code == 'stereo-audit:unknown-unit']
    assert len(f.fixes) == 1


def test_messages_name_atoms_as_numbers():
    report = audit_stereo(read_smiles('C[C@H]1CC2(C1)C[C@H](C)C2'))   # CID 21304360
    assert report.findings[0].message == 'sign at 2 counts only once 4 is signed'
    [f] = audit_stereo(read_smiles('CC=CC')).findings                  # but-2-ene, CID 12220
    assert f.message == 'unit at 2-3 has no sign'


def test_the_audit_is_a_read():
    m = read_smiles('CC(O)CC |&1:1|')
    before = (str(m), len(m.log))
    audit_stereo(m)
    assert (str(m), len(m.log)) == before

    # str() on a bridgehead molecule writes stereo:bridgehead-unrealisable to its log;
    # the audit must not add any record to the input
    bh = read_smiles('C[C@@]12CC[C@H](C1(C)C)CC2(O)C')
    bh_log_before = len(bh.log)
    audit_stereo(bh)
    assert len(bh.log) == bh_log_before


def test_read_time_stereo_records_pass_through():
    m = read_smiles('C[C@H](O)CC |&1:1,&1:9|')                    # names atom 9 of 5
    [f] = [f for f in audit_stereo(m).findings if f.code == 'smiles:stereo-group-bad-index']
    assert f.severity == 'warning' and f.fixes == ()
    m = read_smiles('C[Pt@SP1](F)(Cl)Br')
    [f] = [f for f in audit_stereo(m).findings if f.code == 'smiles:chiral-class-unsupported']
    assert f.severity == 'warning' and not audit_stereo(m).clean


def test_a_non_stereo_read_record_does_not_pass_through():
    m = read_smiles('CCO |$a;b;c;d;e$|')                          # labels for atoms 4 and 5 of 3
    assert [r.rule for r in m.log] == ['smiles:cx-label-bad-index'] * 2
    report = audit_stereo(m)
    assert report.findings == () and report.clean


def test_the_pass_through_list_is_the_read_rows_of_the_docs_index():
    text = Path(__file__).parents[3].joinpath('docs', 'stereo.rst').read_text(encoding='utf-8')
    section = search(r'\nIndex of stereo log rules\n-+\n(.*?)\n\nFive families', text, DOTALL).group(1)
    rows = findall(r'^``([a-z0-9:-]+)``\s+[A-Z]+(?: or [A-Z]+)?\s+(\S+)\s', section, MULTILINE)
    assert len(rows) > 100
    assert _PASS_THROUGH == {rule for rule, stage in rows if stage == 'read'}
    assert not [rule for rule in _PASS_THROUGH if rule.startswith('stereo:')]   # the audit computes these


IBUPROFEN = 'CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O'                   # CID 3672


def test_racemic_ibuprofen_is_a_stated_racemate():
    report = audit_stereo(read_smiles(IBUPROFEN + ' |&1:10|'))
    assert [f.code for f in report.findings] == ['stereo-audit:and-single-member-racemate'] and report.clean


def test_meso_under_one_and_group_reads_as_absolute():
    meso = '[C@@H]([C@@H](C(=O)O)O)(C(=O)O)O'                   # meso-tartaric acid, CID 447315
    [f] = audit_stereo(read_smiles(meso + ' |&1:0,1|')).findings
    assert f.code == 'stereo-audit:meso-absolute' and f.group == (3, 1) and f.severity == 'info'
    assert [(x.action, x.args, x.preserves_meaning) for x in f.fixes] == [('set-group', ('ABS',), True)]


def test_two_and_groups_in_one_component_are_a_split():
    text = 'C[C@H](O)[C@@H](C)CC |&1:1,&2:3|'                    # 3-methylpentan-2-ol, CID 11261
    [f] = [f for f in audit_stereo(read_smiles(text)).findings if f.code == 'stereo-audit:and-split']
    assert f.severity == 'warning'
    merge = [x for x in f.fixes if x.action == 'merge-groups']
    assert len(merge) == 1 and not merge[0].preserves_meaning
    assert len(f.fixes) == 1                                     # keep is the absence of a fix


def test_one_and_group_per_component_is_not_a_split():
    text = IBUPROFEN + '.' + IBUPROFEN + ' |&1:10,&2:25|'
    assert 'stereo-audit:and-split' not in codes(text)


def test_one_and_group_over_two_components_is_not_a_racemate():
    text = IBUPROFEN + '.' + IBUPROFEN + ' |&1:10,25|'
    assert 'stereo-audit:and-single-member-racemate' not in codes(text)


def test_an_and_group_another_group_phases_is_absorbed():
    text = 'C[C@H]1CC[C@@H](C)CC1 |&1:1,&2:4|'                  # 1,4-dimethylcyclohexane, CID 11523
    report = audit_stereo(read_smiles(text))
    [f] = [f for f in report.findings if f.code == 'stereo-audit:and-absorbed']
    assert f.atoms == (2,) and f.group == (3, 1) and f.because == (5,)
    [fix] = f.fixes
    assert fix.action == 'set-group' and fix.args == ('ABS',) and fix.preserves_meaning
    assert 'stereo-audit:and-split' not in [x.code for x in report.findings]


def test_cis_dimethylcyclohexane_under_one_and_group_reads_as_absolute():
    m = read_smiles('C[C@H]1CC[C@H](C)CC1 |&1:1,4|')             # 1,4-dimethylcyclohexane, CID 11523
    [f] = audit_stereo(m).findings
    assert f.code == 'stereo-audit:meso-absolute' and f.atoms == (2, 5)


def test_an_unlabelled_sign_beside_a_group_is_flagged():
    [f] = [f for f in audit_stereo(read_smiles('C[C@H](O)[C@@H](C)CC |&1:1|')).findings
           if f.code == 'stereo-audit:unlabelled-sign']
    assert f.atoms == (4,)
    assert {x.args[0] for x in f.fixes} == {'ABS', 'AND', 'OR'}


def test_an_unlabelled_sign_alone_is_not_flagged():
    assert 'stereo-audit:unlabelled-sign' not in codes('C[C@H](O)[C@@H](C)CC')
    assert 'stereo-audit:unlabelled-sign' not in codes('C[C@H](O)CC.C[C@H](O)CC |&1:6|')


def test_an_axis_member_is_named_by_its_owner_pair():
    # BINOL, CID 11762 (racemic record) -- the axis as an AND member
    m = read_smiles('Oc1ccc2ccccc2c1-c1c(O)ccc2ccccc12')
    axis = [u for u in m.stereogenic_units() if isinstance(u['owners'], tuple)]
    assert axis
    [f] = [f for f in audit_stereo(m).findings if f.code == 'stereo-audit:unknown-unit']
    assert f.atoms == (axis[0]['owners'],)
    for fix in f.fixes:
        assert read_smiles(fix.result) != m


def build(smiles, signs, groups):
    """`smiles` with `signs` {element: parity} and `groups` {element: (kind, id)} set through the edit session."""
    m = read_smiles(smiles)
    anchor = {u['owners']: u['anchor'] for u in m.stereo_units()}
    with m.edit() as e:
        for x, p in signs.items():
            e.set_parity(anchor[x], p)
        for x, (kind, gid) in groups.items():
            e.set_stereo_group(x, kind, gid)
    return m


DEAD = ('stereo-audit:dead-unconfigured', 'stereo-audit:dead-not-stereogenic', 'stereo-audit:dead-restating-group')


def test_a_respelled_label_is_not_dead():
    # Fecht acid, CID 240836: the canonical form states ABS 4 and OR 6 on the spiro atom; neither is dead -- an
    # ABS is never muted and the OR group restating the molecule reads as ABS
    m = build('OC(=O)C1CC2(C1)CC(C2)C(=O)O', {4: 1, 6: 1, 9: 1}, {4: (1, 0), 6: (2, 1), 9: (2, 2)})
    assert m.dead_stereo_groups() == {} and str(canonical(m)).endswith(' |o1:5|')
    assert not [f for f in audit_stereo(m).findings if f.code in DEAD]


def test_every_dead_finding_drops_without_moving_the_meaning():
    from itertools import product
    from chython.core.test.test_stereo_group_meaning import CASES, Reading, labellings

    for _, smiles, elements in CASES:
        r = Reading(smiles, elements)
        for labels in labellings(len(elements)):
            for parities in product((1, 2), repeat=len(elements)):
                for f in audit_stereo(r.build(list(parities), labels)).findings:
                    if f.code in DEAD:
                        assert f.fixes and all(x.preserves_meaning for x in f.fixes), (smiles, labels, f)


def test_a_dead_label_in_another_component_is_not_a_relative_neighbour():
    meso = '[C@@H]([C@@H](C(=O)O)O)(C(=O)O)O'                   # meso-tartaric acid, CID 447315
    assert 'stereo-audit:unlabelled-sign' not in codes(meso + '.C[C@H](O)CC |&1:0,1,13|')


def test_an_axis_in_an_and_group_is_named_by_its_owner_pair():
    # BINOL, CID 11762: the signed axis alone under AND 1
    m = build('Oc1ccc2ccccc2c1-c1c(O)ccc2ccccc12', {(11, 12): 1}, {(11, 12): (3, 1)})
    [f] = audit_stereo(m).findings
    assert f.code == 'stereo-audit:and-single-member-racemate' and f.atoms == ((11, 12),) and f.group == (3, 1)


ALLENOL = 'CC(O)C=C=CC'                                          # hexa-3,4-dien-2-ol, CID 542268


def test_an_axis_and_a_centre_in_two_and_groups_are_a_split():
    m = build(ALLENOL, {2: 1, (4, 6): 1}, {2: (3, 1), (4, 6): (3, 2)})
    [f] = [f for f in audit_stereo(m).findings if f.code == 'stereo-audit:and-split']
    assert f.atoms == (2, (4, 6)) and f.group == ((3, 1), (3, 2))
    [fix] = f.fixes
    assert fix.args == ('AND', 1)
    assert [set(ms) for ms in read_smiles(fix.result).stereo_groups().values()] == [{5, (1, 3)}]


def test_an_unlabelled_axis_takes_the_next_free_ids():
    m = build(ALLENOL, {2: 1, (4, 6): 1}, {2: (3, 1)})
    [f] = [f for f in audit_stereo(m).findings if f.code == 'stereo-audit:unlabelled-sign']
    assert f.atoms == ((4, 6),)
    assert [x.args for x in f.fixes] == [('ABS',), ('AND', 2), ('OR', 1)]
    for fix, kind in zip(f.fixes, (1, 3, 2)):
        assert fix.action == 'set-group' and fix.atoms == ((4, 6),)
        back = read_smiles(fix.result)                           # the written axis is renumbered (1, 3)
        assert back == canonical(_apply_fix(m, fix.action, fix.atoms, fix.args))
        assert any(k == kind and (1, 3) in ms for (k, _), ms in back.stereo_groups().items())


SPIRO = 'C[C@H]1CC2(C1)C[C@H](C)C2'                              # 2,6-dimethylspiro[3.3]heptane, CID 21304360
CAMPHOR = 'C[C@@]12CC[C@@H](C1(C)C)CC2=O'                        # (1R,4R)-camphor, CID 159055


def test_spiro_terminals_wait_for_the_centre():
    found = [f for f in audit_stereo(read_smiles(SPIRO)).findings if f.code == 'stereo-audit:dead-waits-for']
    assert sorted(f.atoms for f in found) == [(2,), (7,)]
    assert all(f.because == (4,) for f in found)
    for f in found:
        assert sorted(x.args for x in f.fixes) == [(1,), (2,)]
        assert all(x.action == 'set-sign' and x.atoms == (4,) for x in f.fixes)
        assert all(read_smiles(x.result) != read_smiles(SPIRO) for x in f.fixes)


def test_waiting_terminals_replace_their_dead_label():
    report = audit_stereo(read_smiles(SPIRO + ' |a:1,6|'))
    assert sorted(f.atoms for f in report.findings if f.code == 'stereo-audit:dead-waits-for') == [(2,), (7,)]
    assert not [f for f in report.findings if f.code in DEAD]


def test_a_signed_spiro_centre_is_clean():
    assert codes('C[C@H]1C[C@]2(C1)C[C@H](C)C2') == []


def test_mixed_and_spiro_terminals_are_a_racemate_not_waiting():
    # the designed spelling of the racemate is the centre under one AND: C[C@]2 |&1:3|
    assert 'stereo-audit:dead-waits-for' not in codes(SPIRO + ' |&1:1,&2:6|')


def test_an_impossible_bridgehead_pair_is_an_error():
    report = audit_stereo(read_smiles('C[C@@]12CC[C@H](C1(C)C)CC2=O'))
    [f] = [f for f in report.findings if f.code == 'stereo-audit:bridgehead-unrealisable']
    assert f.severity == 'error' and not report.clean and sorted(f.atoms) == [2, 5]
    assert len(f.fixes) == 2 and all(read_smiles(x.result) == read_smiles(CAMPHOR) or
                                     read_smiles(x.result) == read_smiles(CAMPHOR.replace('@@', '@'))
                                     for x in f.fixes)


NORBORNANE_BAD = 'C1C[C@@H]2CC[C@@H]1C2'                        # norbornane, CID 9233, signs no bicycle realises
TROPINONE_BAD = 'CN1[C@@H]2CC[C@@H]1CC(=O)C2'                    # tropinone, CID 446337, C5 inverted


def test_an_unrealisable_symmetric_bicycle_is_an_error():
    for text in (NORBORNANE_BAD, TROPINONE_BAD):
        m = read_smiles(text)
        str(m)                                                   # a cached write logs nothing on a copy
        report = audit_stereo(m)
        [f] = report.findings
        assert f.code == 'stereo-audit:bridgehead-unrealisable' and f.atoms == (3, 6) and not report.clean
        [fix] = f.fixes                                          # flipping either member gives one compound
        fixed = _apply_fix(m, fix.action, fix.atoms, fix.args)
        assert fix.action == 'set-sign' and audit_stereo(fixed).clean


CAMPHOR_BAD = 'C[C@@]12CC[C@H](C1(C)C)CC2=O'                     # camphor, CID 2537, C4 inverted


def test_two_unrealisable_bicycles_are_each_fixed_alone():
    for text in (CAMPHOR_BAD + '.' + CAMPHOR_BAD, NORBORNANE_BAD + '.' + CAMPHOR_BAD):
        m = read_smiles(text)
        found = [f for f in audit_stereo(m).findings if f.code == 'stereo-audit:bridgehead-unrealisable']
        assert len(found) == 2 and all(f.fixes for f in found), text
        fixed = m
        for f in found:                                          # one fix per set reaches a clean report
            fix = f.fixes[0]
            fixed = _apply_fix(fixed, fix.action, fix.atoms, fix.args)
        assert audit_stereo(fixed).clean, text


def test_a_set_no_single_flip_realises_is_cleared():
    # adamantane, CID 9238, bridgeheads signed (1, 1, 1, 2)
    m = build('C1C2CC3CC1CC(C2)C3', {2: 1, 4: 1, 6: 1, 8: 2}, {})
    [f] = audit_stereo(m).findings
    assert f.code == 'stereo-audit:bridgehead-unrealisable' and f.atoms == (2, 4, 6, 8)
    [fix] = f.fixes
    assert fix.action == 'clear-sign' and fix.atoms == (2, 4, 6, 8)
    fixed = _apply_fix(m, fix.action, fix.atoms, fix.args)
    assert audit_stereo(fixed).clean and str(canonical(fixed)) == fix.result


def test_a_realisable_or_achiral_bridgehead_pair_is_not_an_error():
    assert codes(CAMPHOR) == []
    assert 'stereo-audit:bridgehead-unrealisable' not in codes('C1C[C@@H]2CC[C@H]1C2')   # norbornane, CID 9233


def test_a_half_signed_bridgehead_pair_is_implied():
    [f] = [f for f in audit_stereo(read_smiles('C[C@@]12CCC(C1(C)C)CC2=O')).findings
           if f.code == 'stereo-audit:bridgehead-implied']
    assert f.atoms == (5,) and f.because == (2,) and f.fixes == ()


def test_an_implied_bridgehead_is_not_an_unlabelled_sign():
    # the set is one element labelled on every member, so the unlabelled partner is reported, and it is no
    # unlabelled sign
    assert codes('C[C@@]12CCC(C1(C)C)CC2=O |&1:1|') == ['stereo-audit:bridgehead-half-labelled',
                                                         'stereo-audit:bridgehead-implied']
    assert 'stereo-audit:unlabelled-sign' not in codes(CAMPHOR + ' |&1:1|')


@mark.parametrize('tail, code, atom, source, action', [
    ('|&1:1,o1:4|', 'redundant-label', 5, 2, 'set-group'),      # AND outranks OR
    ('|o1:1,&1:4|', 'redundant-label', 2, 5, 'set-group'),
    ('|a:1,&1:4|', 'redundant-label', 2, 5, 'set-group'),       # AND outranks ABS
    ('|&1:1|', 'half-labelled', 5, 2, 'set-group'),             # one element, labelled on every member
    ('|&1:1,&2:4|', 'merged-groups', 5, 2, 'merge-groups'),     # the cage inverts with both groups
])
def test_a_second_label_on_a_coupled_set_is_reported_with_a_meaning_preserving_fix(tail, code, atom, source,
                                                                                    action):
    [f] = [f for f in audit_stereo(read_smiles(f'{CAMPHOR} {tail}')).findings
           if f.code.startswith('stereo-audit:bridgehead-') and f.code != 'stereo-audit:bridgehead-implied']
    assert f.code == f'stereo-audit:bridgehead-{code}' and f.atoms == (atom,) and f.because == (source,)
    assert f.fixes[0].action == action and all(fix.preserves_meaning for fix in f.fixes)
    assert not [g for g in audit_stereo(read_smiles(f.fixes[0].result)).findings
                if g.code.startswith('stereo-audit:bridgehead-') and g.code != 'stereo-audit:bridgehead-implied']


@mark.parametrize('tail', ['|&1:1,o1:4|', '|o1:1,&1:4|', '|a:1,&1:4|'])
def test_the_free_bridgehead_is_the_strongest_label_whatever_the_numbering(tail):
    m = read_smiles(f'{CAMPHOR} {tail}')
    [f] = [f for f in audit_stereo(m).findings if f.code == 'stereo-audit:bridgehead-redundant-label']
    m.remap({2: 5, 5: 2})
    [g] = [g for g in audit_stereo(m).findings if g.code == 'stereo-audit:bridgehead-redundant-label']
    assert g.atoms == ({2: 5, 5: 2}[f.atoms[0]],) and g.fixes[0].result == f.fixes[0].result


def test_a_half_signed_set_takes_no_group_from_its_unsigned_member():
    # the set's group comes from its grouped signed member; C4 is unsigned, so the set writes no group
    m = read_smiles('C[C@@]12CCC(C1(C)C)CC2=O |&1:4|')
    assert m.stereo_groups() and m.dead_stereo_groups() == {}
    assert codes('C[C@@]12CCC(C1(C)C)CC2=O |&1:4|') == [
        'stereo-audit:bridgehead-implied', 'stereo-audit:bridgehead-redundant-label']        # no unlabelled-sign


def test_a_stated_centre_is_not_unknown():
    # the mixed-AND racemate states the spiro centre: signing it keeps the meaning
    assert 'stereo-audit:unknown-unit' not in codes(SPIRO + ' |&1:1,&2:6|')


def test_an_unrealisable_pair_is_named_only_by_its_bridgehead_finding():
    report = audit_stereo(read_smiles('C[C@@]12CC[C@H](C1(C)C)CC2=O |&1:4|'))
    assert [f.code for f in report.findings] == ['stereo-audit:bridgehead-unrealisable']


def test_a_waited_for_centre_is_named_by_the_waits_for_finding_only():
    assert codes(SPIRO) == ['stereo-audit:dead-waits-for'] * 2
    [f] = audit_stereo(read_smiles('CC1CC[C@H](C)CC1')).findings     # 1,4-dimethylcyclohexane, CID 11523
    assert f.code == 'stereo-audit:dead-waits-for' and f.because == (2,) and len(f.fixes) == 2


TABLE = [
    'CC(O)CC', 'CC(O)CC |&1:1|', 'C[C@H](C)CC |&1:1|', IBUPROFEN + ' |&1:10|',
    '[C@@H]([C@@H](C(=O)O)O)(C(=O)O)O |&1:0,1|', 'C[C@H](O)[C@@H](C)CC |&1:1,&2:3|',
    'C[C@H](O)[C@@H](C)CC |&1:1|', SPIRO, 'C[C@@]12CC[C@H](C1(C)C)CC2=O', 'C[C@@]12CCC(C1(C)C)CC2=O',
    CAMPHOR + ' |&1:1,&2:4|', SPIRO + ' |&1:1,&2:6|', 'C[C@@]12CC[C@H](C1(C)C)CC2=O |&1:4|',
    'C[C@H]1CC[C@@H](C)CC1 |&1:1,&2:4|', 'CC12CCC(C1(C)C)CC2=O', NORBORNANE_BAD, TROPINONE_BAD,
    'C[C@@]12CCC(C1(C)C)CC2=O |&1:4|', CAMPHOR_BAD + '.' + CAMPHOR_BAD, NORBORNANE_BAD + '.' + CAMPHOR_BAD,
]


def _computed(report):
    return [f for f in report.findings if f.code.startswith('stereo-audit:')]


@mark.parametrize('text', TABLE)
def test_every_fix_removes_its_finding(text):
    m = read_smiles(text)
    for f in _computed(audit_stereo(m)):
        for fix in f.fixes:
            fixed = _apply_fix(m, fix.action, fix.atoms, fix.args)          # same atom numbers as `m`
            after = audit_stereo(fixed)
            assert (f.code, f.atoms) not in {(g.code, g.atoms) for g in _computed(after)}, (text, f, fix)
            assert str(canonical(fixed)) == fix.result
            if fix.preserves_meaning:
                assert canonical(fixed) == canonical(m)


# A write drops a sign or fills an implied one and keeps every label: the finding before, the ones it may become
# after, and whether their atoms stay.
AFTER_A_WRITE = {'stereo-audit:bridgehead-unrealisable': (('stereo-audit:unknown-unit', True),
                                                          ('stereo-audit:dead-unconfigured', True)),
                 'stereo-audit:dead-not-stereogenic': (('stereo-audit:dead-unconfigured', False),),
                 'stereo-audit:bridgehead-redundant-label': (('stereo-audit:bridgehead-half-labelled', False),),
                 'stereo-audit:and-absorbed': (('stereo-audit:and-single-member-racemate', False),)}


@mark.parametrize('text', TABLE)
def test_the_audit_survives_a_round_trip(text):
    m = read_smiles(text)
    again = read_smiles(str(m))
    to_m = again.isomorphism(m)
    before = _computed(audit_stereo(m))
    kept = {(f.code, f.severity) for f in before}
    for g in _computed(audit_stereo(again)):
        if (g.code, g.severity) in kept:
            continue
        atoms = {to_m[x] for x in g.atoms if isinstance(x, int)}
        assert any(code == g.code and (not stays or atoms <= set(f.atoms))
                   for f in before for code, stays in AFTER_A_WRITE.get(f.code, ())), (text, g)


def test_the_docs_table_names_every_code_and_severity():
    text = Path(__file__).parents[3].joinpath('docs', 'stereo.rst').read_text(encoding='utf-8')
    section = search(r'\nAuditing labels\n-+\n(.*?)\n[^\n]+\n-+\n', text, DOTALL).group(1)
    rows = dict(findall(r'^``stereo-audit:([a-z-]+)``\s+(error|warning|info)\s', section, MULTILINE))
    assert rows == {code: row.severity for code, row in stereo_audit_codes().items()}


def test_an_unsigned_bridgehead_pair_is_one_unknown_unit():
    m = read_smiles('CC12CCC(C1(C)C)CC2=O')                       # camphor, CID 2537
    [f] = audit_stereo(m).findings
    assert f.code == 'stereo-audit:unknown-unit' and f.atoms == (2,) and f.because == (5,)
    assert [x.args for x in f.fixes] == [(1,), (2,)]
    for x in f.fixes:
        assert codes(x.result) == []                            # the written pair is realisable


#: cis-2,6-dimethylmorpholine, N-phenyl: a configuration that is its own mirror image
MORPHOLINE = 'C[C@H]1CN(C[C@@H](C)O1)c1ccccc1'


@mark.parametrize('tail', ['|&1:1,5|', '|o1:1,5|'])
def test_a_group_on_a_meso_configuration_reads_as_absolute(tail):
    """The group separates the drawn diastereomer from the other, as ABS does; set-group ABS keeps the meaning."""
    [f] = [f for f in audit_stereo(read_smiles(f'{MORPHOLINE} {tail}')).findings
           if f.code == 'stereo-audit:meso-absolute']
    assert f.atoms == (2, 6) and f.fixes[0].action == 'set-group' and f.fixes[0].preserves_meaning


def test_abs_on_a_meso_configuration_is_clean_and_is_not_the_unlabelled_drawing():
    assert codes(f'{MORPHOLINE} |a:1,5|') == []
    assert read_smiles(f'{MORPHOLINE} |a:1,5|') != read_smiles(MORPHOLINE)


def test_a_meso_group_beside_another_group_reads_as_absolute():
    """cis-1,3-Dimethylcyclobutane in AND beside an OR butan-2-ol: the pair is its own mirror image, so the label
    states the cis isomer as ABS does -- not a phase another group fixes -- and is kept as ABS."""
    m = read_smiles('C[C@H]1C[C@@H](C)C1.C[C@H](O)CC |&1:1,3,o1:7|')
    assert str(canonical(m)) == 'C1[C@@H](C)C[C@H]1C.C[C@H](O)CC |a:1,4,o1:7|'
    assert [f.code for f in audit_stereo(m).findings] == ['stereo-audit:meso-absolute']


def test_a_copy_is_audited_as_the_original_is():
    """Reader findings come from `molecule.log`, which `copy()` carries."""
    mol = read_smiles('CC(O)CC {A1=q}')
    codes = [f.code for f in audit_stereo(mol).findings]
    assert 'smiles:cip-undetermined' in codes
    assert [f.code for f in audit_stereo(mol.copy()).findings] == codes


def _oligolactate(n):
    """HO-[CH(CH3)-C(=O)-O]n-H, L throughout, the first two centres in `&1`: every other sign is unlabelled
    beside a live relative group, one `unlabelled-sign` each, three fixes apiece -- the audit's widest fan
    of copies of one constitution."""
    return read_smiles('O' + '[C@@H](C)C(=O)O' * n + ' |&1:1,6|')


@mark.parametrize('text', ['C[C@H]1C[C@@H](C)C1.C[C@H](O)CC |&1:1,3,o1:7|', 'C[C@@]12CC[C@H](C1(C)C)CC2=O',
                           'NC(=O)C1[C@H](C(N)=O)[C@H]1C(N)=O |&1:4,8|', 'C[C@H](O)[C@H](O)[C@@H](C)O |&1:1,o1:3|'])
def test_the_identity_memo_changes_no_finding(text):
    """`audit_stereo` runs in an `identity_memo()` scope; without it, every finding and fix is the same."""
    for m in (read_smiles(text), _oligolactate(8)):
        assert audit_stereo(m).findings == _audit(m).findings


def test_reported_fixes_are_built_on_first_read_and_equal_the_eager_ones(monkeypatch):
    """4-aminopentan-2-ol with C2 in `&1` and C4 signed alone: `unlabelled-sign` only reports its fixes, so
    they wait until read; the snapshot keeps them those of the audited molecule, and they pickle as a tuple."""
    import pickle
    from chython.chemistry._stereo_audit import Fixes, _Context
    text = 'C[C@H](O)C[C@H](N)C |&1:1|'
    mol = read_smiles(text)
    [finding] = audit_stereo(mol).findings
    assert finding.code == 'stereo-audit:unlabelled-sign' and isinstance(finding.fixes, Fixes)
    assert finding.fixes._fixes is None, 'nothing read yet, nothing built'
    with mol.edit() as e:
        e.set_parity(5, 0)
    lazy = tuple(finding.fixes)
    monkeypatch.setattr(_Context, 'later', _Context.fix)
    [eager] = audit_stereo(read_smiles(text)).findings
    assert lazy == tuple(eager.fixes) and finding == eager and hash(finding.fixes) == hash(lazy)
    assert len(lazy) == 3 and pickle.loads(pickle.dumps(finding.fixes)) == lazy


@mark.parametrize('trans, cis, code', [
    ('[C@H]12C[C@@H]1CCC2', '[C@H]12C[C@H]1CCC2', 'stereo-audit:bridgehead-unrealisable'),      # [3.1.0]
    ('[C@H]12CC[C@@H]1CC2', '[C@H]12CC[C@H]1CC2', 'stereo-audit:bridgehead-unrealisable'),      # [2.2.0]
    ('[C@H]12C[C@@H]1CCCC2', '[C@H]12C[C@H]1CCCC2', 'stereo-audit:fused-strained'),            # [4.1.0]
    ('[C@H]12CC[C@@H]1CCC2', '[C@H]12CC[C@H]1CCC2', 'stereo-audit:fused-strained'),            # [3.2.0]
    ('[C@H]12C[C@@H]1CCCCC2', '[C@H]12C[C@H]1CCCCC2', None),                                  # [5.1.0]
    ('[C@H]12CC[C@@H]1CCCC2', '[C@H]12CC[C@H]1CCCC2', None),                                  # [4.2.0]
])
def test_a_trans_fused_small_bicycle_is_unrealisable_or_strained(trans, cis, code):
    """Trans fusion by size, MMFF over an RDKit embedding against the cis isomer: [3.1.0] and [2.2.0] do not
    embed; [4.1.0] +25 kcal/mol and [3.2.0] +20 are made but strained; [5.1.0] +9 and [4.2.0] +4 are ordinary.
    Cis is never reported."""
    found = {f.code for f in audit_stereo(read_smiles(trans)).findings}
    assert (code in found) if code else not found & {'stereo-audit:bridgehead-unrealisable',
                                                       'stereo-audit:fused-strained'}
    assert not {f.code for f in audit_stereo(read_smiles(cis)).findings}


def test_a_strained_fusion_keeps_its_signs_and_offers_the_cis_fusion():
    """`fused-strained` is a warning: both signs stay live, and either fix flips one fusion atom to cis."""
    mol = read_smiles('[C@H]12C[C@@H]1CCCC2')
    [finding] = audit_stereo(mol).findings
    assert finding.severity == 'warning' and finding.atoms == (1, 3)
    assert all(mol.stereo_is_live(a) for a in (1, 3)) and '@' in str(mol)
    for fix in finding.fixes:
        assert not audit_stereo(read_smiles(fix.result)).findings


@mark.parametrize('name,text,count', [
    ('4-methylcyclohexylideneacetic acid', 'CC1CCC(=CC(=O)O)CC1', 2),
    ('spiro[3.3]heptane-2,6-dicarboxylic acid', 'OC(=O)C1CC2(C1)CC(C2)C(=O)O', 3),
    ('cis/trans-1,4-dimethylcyclohexane', 'CC1CCC(C)CC1', 2),
])
def test_units_stating_a_configuration_only_together_are_one_unknown_unit(name, text, count):
    """An exocyclic C=C and its ring centre, a spiro atom and the ring centres beside it: one sign alone
    restates the molecule, so the set is one finding naming every member."""
    findings = [f for f in audit_stereo(read_smiles(text)).findings if f.code == 'stereo-audit:unknown-unit']
    assert len(findings) == 1 and len(findings[0].atoms) == count, (name, findings)


def test_a_partly_signed_axial_set_waits_for_its_missing_member():
    """C4 signed and the C=C open: the sign states nothing until the double bond is signed."""
    [f] = audit_stereo(read_smiles('C[C@H]1CCC(=CC(=O)O)CC1')).findings
    assert f.code == 'stereo-audit:dead-waits-for' and f.because == ((5, 6),)
