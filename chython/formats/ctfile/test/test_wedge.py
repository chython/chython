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
"""Drawing to configuration, and back.  The corpus oracle -- chython 2 reading the same bytes, bridged by
V2's own ``_translate_tetrahedron_sign(n, env)`` -- runs first: a round trip and a hand-written geometry
check both pass with every parity flipped, so only a second stack catches a global sign inversion.
"""

from pytest import fixture

from ....core import WEDGE_NONE, MoleculeContainer, read_smiles
from .._sdf import sniff_version, split_records
from .._v2000 import emit_v2000, parse_v2000
from .._v3000 import V3000_STAMP, emit_v3000, parse_v3000
from ....core.wedge import (SU_ALLENE, SU_ATROPISOMER, SU_CIS_TRANS, SU_TETRA, _near_collinear,
                            allene_parity, assign_parities, cis_trans_for_write, cis_trans_letter,
                            cis_trans_parity, signed_volume, stated_cis_trans, stated_parity,
                            tetrahedral_parity, wedges_for_write)


def _build(record):
    parse = parse_v3000 if sniff_version(record, []) == V3000_STAMP else parse_v2000
    return parse(record, []).build()


def _tetra(mol):
    return [u for u in mol.stereo_units() if u['kind'] == SU_TETRA and u['stereogenic']]


# --- the oracle

def test_every_configuration_agrees_with_chython2_or_the_record_contradicts_itself(corpus,
                                                                                   v2_molecules):
    """Two independent stacks, same bytes, same configurations.  Each disagreement must be a centre
    whose wedge drawing and atom parity field contradict each other, with a log line naming it -- the
    drawing wins there, which is a decision and not a tolerance.
    """
    agree = disagree = 0
    for name, records in corpus.items():
        v2 = v2_molecules.get(name)
        if v2 is None or len(v2) != len(records):
            continue
        for n, (record, m2) in enumerate(zip(records, v2)):
            try:
                mol, _, log = _build(record)
            except Exception:
                continue
            sids = list(mol.atom_numbers)
            nums = list(m2)
            if len(sids) != len(nums):
                continue
            to_num = dict(zip(sids, nums))
            for unit in _tetra(mol):
                anchor = unit['anchor']
                ours = mol.parity_of(anchor)
                num = to_num[anchor]
                if m2.atom(num).stereo is None or not ours:
                    continue
                env = tuple(to_num[r] for r in unit['refs'] if r is not None)
                # `None` is V2 declining the frame: it cannot express this one, so it is no oracle
                # for this centre.
                sign = m2.tetrahedron_sign(num, env)
                if sign is None:
                    continue
                theirs = 2 if sign else 1
                if theirs == ours:
                    agree += 1
                    continue
                disagree += 1
                assert any(f'atom {anchor}:' in x and 'disagree' in x for x in log), (
                    f'{name} record {n} atom {anchor}: parity {ours} against chython 2\'s {theirs}, '
                    f'and nothing in the log explains it')
    assert agree > 900, f'only {agree} configurations compared; the oracle stopped covering'
    assert disagree < agree / 100, f'{disagree} unexplained-by-contradiction disagreements'


def test_no_centre_that_chython2_configures_is_left_unset(corpus, v2_molecules):
    """The half of the oracle that catches dropped stereo rather than inverted stereo.  One gap is
    admitted: 6 ring-constrained centres the core calls non-stereogenic -- a configuration stored there is
    one the core will not translate.
    """
    missing = []
    not_stereogenic = 0
    for name, records in corpus.items():
        v2 = v2_molecules.get(name)
        if v2 is None or len(v2) != len(records):
            continue
        for n, (record, m2) in enumerate(zip(records, v2)):
            try:
                mol, _, _ = _build(record)
            except Exception:
                continue
            sids = list(mol.atom_numbers)
            nums = list(m2)
            if len(sids) != len(nums):
                continue
            units = {u['anchor']: u for u in mol.stereo_units()}
            for anchor, num in zip(sids, nums):
                if m2.atom(num).stereo is None or mol.parity_of(anchor):
                    continue
                unit = units.get(anchor)
                if unit is None or not unit['stereogenic']:
                    not_stereogenic += 1
                else:
                    missing.append(f'{name} record {n} atom {num}')
    assert not missing, f'{len(missing)} centres configured by chython 2 and unset here: {missing[:5]}'
    assert not_stereogenic == 6, (f'{not_stereogenic} centres the core calls non-stereogenic and '
                                  f'chython 2 configures; the perception gap moved')


def test_the_corpus_actually_exercises_every_kind_this_module_reads(corpus):
    """Coverage guard on the two oracle tests above: without it they pass on a corpus with no stereo.
    The two axial counts are exact.  The allene count is the number
    ``test_every_allene_configuration_agrees_with_chython2`` compares; one of the 11 (``stereo.sdf``
    record 140) is a five-carbon cumulene, so the axial reading is exercised past three carbons.  The
    atropisomer count has no oracle -- chython 2 has no axial perception of a biaryl -- so it is a guard on
    the corpus keeping the 6 axes that ``atropisomer_parity`` reads, of the 7 it perceives.
    """
    tetra = bond = allene = atropisomer = 0
    for records in corpus.values():
        for record in records:
            try:
                mol, _, _ = _build(record)
            except Exception:
                continue
            for unit in mol.stereo_units():
                if not (unit['stereogenic'] and mol.parity_of(unit['anchor'])):
                    continue
                if unit['kind'] == SU_TETRA:
                    tetra += 1
                elif unit['kind'] == SU_CIS_TRANS:
                    bond += 1
                elif unit['kind'] == SU_ALLENE:
                    allene += 1
                elif unit['kind'] == SU_ATROPISOMER:
                    atropisomer += 1
    assert tetra > 1000, tetra
    assert bond > 100, bond
    assert allene == 11, allene
    assert atropisomer == 6, atropisomer


def test_the_stereo_unit_kind_constants_match_what_the_core_emits():
    """``SU_*`` are mirrored from Cython ``DEF`` constants, which do not exist at runtime, so a
    renumbering in the core would silently turn every tetrahedron here into an allene."""
    mol = MoleculeContainer()
    with mol.edit():
        c = mol.add_atom('C')
        for element in ('F', 'Cl', 'Br', 'I'):
            mol.add_bond(c, mol.add_atom(element), 1)
    kinds = {u['kind'] for u in mol.stereo_units()}
    assert kinds == {SU_TETRA}, f'a four-substituent carbon is kind {kinds}, not SU_TETRA'
    assert (SU_TETRA, SU_CIS_TRANS, SU_ALLENE, SU_ATROPISOMER) == (0, 1, 2, 3)


# --- the geometry

def test_signed_volume_is_positive_for_a_clockwise_triple_seen_from_the_first_point():
    """The sign convention everything else is expressed in, stated once on numbers small enough to
    check by hand."""
    assert signed_volume((0, 0, 1), (1, 0, 0), (0, 1, 0), (-1, -1, 0)) < 0
    assert signed_volume((0, 0, 1), (0, 1, 0), (1, 0, 0), (-1, -1, 0)) > 0


def test_swapping_two_frame_directions_inverts_the_parity():
    """The defining property of a parity, and the reason every frame in this package is explicit."""
    v = [(0, 0, 1), (1, 0, 0), (0, 1, 0), (-1, -1, 0)]
    a = signed_volume(*v)
    v[1], v[2] = v[2], v[1]
    assert a * signed_volume(*v) < 0


def test_a_flat_drawing_with_no_wedge_states_no_configuration():
    """A flat drawing of a stereocentre states no configuration; inventing one from the coordinates
    alone would fabricate an enantiomer."""
    mol, _, log = _build(_TETRA_FLAT)
    assert [u for u in mol.stereo_units() if u['stereogenic']]
    assert all(not mol.parity_of(u['anchor']) for u in mol.stereo_units())
    assert not any('either' in x or 'degenerate' in x for x in log), \
        'a flat drawing is ordinary; it needs no explanation'


def test_a_wavy_either_bond_leaves_the_centre_unset_and_says_so():
    """Bond stereo 4 means "up or down, unknown which", which is a statement -- and it is not the same
    statement as an unmarked bond, so it gets a log line where a flat drawing does not."""
    lines = list(_TETRA_FLAT)
    lines[8] = '  1  2  1  4  0  0  0'
    mol, _, log = _build(lines)
    assert all(not mol.parity_of(u['anchor']) for u in mol.stereo_units())
    assert any('either' in x for x in log), log


def test_an_up_wedge_and_a_down_wedge_on_the_same_bond_give_opposite_parities():
    """The minimum guarantee: the drawing is being read at all, and its two directions differ."""
    up = list(_TETRA_FLAT)
    up[8] = '  1  2  1  1  0  0  0'
    down = list(_TETRA_FLAT)
    down[8] = '  1  2  1  6  0  0  0'
    a, _, _ = _build(up)
    b, _, _ = _build(down)
    pa = [a.parity_of(u['anchor']) for u in _tetra(a)]
    pb = [b.parity_of(u['anchor']) for u in _tetra(b)]
    assert pa and pa != [0]
    assert pa == [3 - x for x in pb], f'{pa} against {pb}'


def test_which_direction_gets_which_label_is_the_one_chython2_gives(v2_reader):
    """The sign pinned on one small fixture, expected value fetched from chython 2 rather than written
    down: ``volume < 0 == anticlockwise == SMILES @ == parity 2``.
    """
    lines = list(_TETRA_FLAT)
    lines[8] = '  1  2  1  1  0  0  0'
    mol, _, _ = _build(lines)
    m2 = v2_reader('\n'.join(lines))
    unit, = _tetra(mol)
    # our frame, expressed in V2's numbering -- the two stacks number this fixture's atoms alike
    env = tuple(r for r in unit['refs'] if r is not None)
    assert m2.atom(unit['anchor']).stereo is not None, 'chython 2 read no configuration here'
    sign = m2.tetrahedron_sign(unit['anchor'], env)
    assert sign is not None, 'chython 2 cannot express this frame, so it pins nothing'
    assert mol.parity_of(unit['anchor']) == (2 if sign else 1)


def test_a_wedge_on_a_non_stereogenic_centre_is_reported_and_not_stored():
    """Storing the parity would claim a configuration the constitution cannot distinguish."""
    lines = list(_TETRA_FLAT)
    lines[7] = '    1.0000    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0'
    lines[8] = '  1  2  1  1  0  0  0'
    mol, _, log = _build(lines)
    assert all(not mol.parity_of(u['anchor']) for u in mol.stereo_units())
    assert any('non-stereogenic' in x for x in log), log


def test_a_degenerate_layout_is_reported_rather_than_resolved():
    """A zero determinant is no configuration to read, and the message must say the file needs a 2D
    clean.  Fixture: four *distinct* substituents (duplicating one makes the centre non-stereogenic and
    changes the code path) with the three unwedged neighbours drawn collinear -- two substituents at the
    *same point* is not this case, since the wedge still lifts one out of the plane.
    """
    mol, _, log = _build(_TETRA_DEGENERATE)
    assert _tetra(mol), 'the fixture must keep a stereogenic centre'
    assert all(not mol.parity_of(u['anchor']) for u in mol.stereo_units())
    assert any('degenerate' in x for x in log), log


# --- the wedge drawn from the wrong end

def test_a_wedge_whose_point_is_at_the_wrong_end_is_re_anchored_and_read():
    """The whole claim, and it is an equality against the same drawing done right: moving the point to
    the end that can hold a configuration reads the configuration the drawing shows, and the repair is
    logged rather than silent.
    """
    mol, _, log = _build(_TETRA_WEDGE_AT_WIDE_END)
    lines = list(_TETRA_WEDGE_AT_WIDE_END)
    lines[9] = '  2  1  1  6  0  0  0'                      # the same bond, point at the centre
    control, _, control_log = _build(lines)
    unit, = _tetra(mol)
    assert mol.parity_of(unit['anchor']) == control.parity_of(unit['anchor']) != 0
    assert format(mol, 'A') == format(control, 'A') == '[C@@H](CBr)(C)N'
    assert list(mol.wedges()) == list(control.wedges()), 'the drawing itself is repaired, not just read'
    assert any('re-anchored' in x for x in log), log
    assert not any('re-anchored' in x for x in control_log), 'nothing to repair in the control'


def test_a_re_anchored_wedge_is_what_gets_written_back():
    """The repair is on the record, so the bond line the writer emits names the centre first -- which is
    what makes the round trip give the configuration back rather than losing it a second time."""
    mol, _, _ = _build(_TETRA_WEDGE_AT_WIDE_END)
    unit, = _tetra(mol)
    lines, _ = emit_v2000(mol)
    again, _, _ = _build(lines)
    assert again.parity_of(unit['anchor']) == mol.parity_of(unit['anchor'])
    assert any(line.startswith('  2  1  1  6') for line in lines), lines


def test_a_wedge_between_two_centres_is_left_where_the_file_drew_it():
    """The narrow end holds a configuration of its own, so the drawing means what it says and there is
    nothing to repair.  Which end of a wedge a reader keys on is the convention; guessing past it on a
    bond where both ends are stereogenic would move the statement to the other centre.
    """
    mol, _, log = _build(_DIBROMOBUTANE)
    assert list(mol.wedges()) == [(2, 3, 1)]
    assert mol.parity_of(2) and not mol.parity_of(3)
    assert not any('re-anchored' in x for x in log), log


def test_a_wedge_pointing_at_a_centre_that_already_has_one_is_left_alone():
    """Two wedges at one centre is how a drawing contradicts itself, and this pass must not be the thing
    that creates one."""
    lines = list(_TETRA_WEDGE_AT_WIDE_END)
    lines[11] = '  2  4  1  1  0  0  0'                     # the centre states itself already
    mol, _, log = _build(lines)
    assert (1, 2, 2) in list(mol.wedges()), 'the misanchored wedge stays as drawn'
    assert not any('re-anchored' in x for x in log), log


# --- the atom parity field

class _Unit(dict):
    """The two keys :func:`stated_parity` reads, so the frame arithmetic can be tested on its own."""
    def __init__(self, anchor, refs):
        super().__init__(anchor=anchor, refs=refs, kind=SU_TETRA, stereogenic=True)


def test_a_field_that_is_not_a_configuration_yields_nothing():
    """Atom parity 0 is "unstated" and 3 is "either or racemic"; neither is a configuration."""
    positions = {1: 0, 2: 1, 3: 2}
    for field in (0, 3, -1, 9):
        assert stated_parity(_Unit(1, (2, 3, None, None)), field, positions) == 0


def test_the_field_is_taken_verbatim_when_the_frames_coincide():
    """The measured mapping: odd (1) is core parity 1, even (2) is core parity 2, no inversion.  See
    ``stated_parity``."""
    positions = {10: 0, 20: 1, 30: 2, 40: 3}
    unit = _Unit(10, (20, 30, 40, None))
    assert stated_parity(unit, 1, positions) == 1
    assert stated_parity(unit, 2, positions) == 2


def test_a_core_frame_out_of_file_order_inverts_the_field_by_the_permutation():
    """The field is stated against ascending atom-block position while the core reports its own order,
    so one transposition between them flips the sign.  The two coincide often enough that assuming it
    survives casual testing."""
    positions = {10: 0, 20: 1, 30: 2, 40: 3}
    swapped = _Unit(10, (30, 20, 40, None))  # one transposition from ascending
    assert stated_parity(swapped, 1, positions) == 2
    assert stated_parity(swapped, 2, positions) == 1
    rotated = _Unit(10, (30, 40, 20, None))  # a 3-cycle: even, so no flip
    assert stated_parity(rotated, 1, positions) == 1


def test_the_undrawn_direction_ranks_above_every_atom_in_the_field_s_frame():
    """MDL counts an implicit hydrogen as the highest-numbered neighbour, so the hole sorts last.
    Placing it first is a cyclic rotation of four -- odd -- and inverts every three-substituent
    centre."""
    positions = {10: 0, 20: 1, 30: 2, 40: 3}
    assert stated_parity(_Unit(10, (20, 30, 40, None)), 1, positions) == 1
    assert stated_parity(_Unit(10, (None, 20, 30, 40)), 1, positions) == 2


def test_two_undrawn_directions_make_the_frame_ambiguous_and_are_reported():
    log = []
    assert stated_parity(_Unit(10, (20, 30, None, None)), 1, {10: 0, 20: 1, 30: 2}, log) == 0
    assert any('ambiguous' in x for x in log), log


def test_the_field_agrees_with_the_drawing_on_the_whole_corpus(corpus):
    """The oracle for :func:`stated_parity`'s frame.  Every centre stating a configuration twice is a
    paired measurement of two independent encodings: 1046 of 1053 agree for this reading and 7 for the
    specification's own "highest-numbered neighbour first", the two being exact opposites.  A frame or
    sign error inverts the rate rather than nudging it, so the bound can be tight.
    """
    both = agree = 0
    for records in corpus.values():
        for record in records:
            parse = parse_v3000 if sniff_version(record, []) == V3000_STAMP else parse_v2000
            try:
                ctab = parse(record, [])
                mol, _, _ = ctab.build()
            except Exception:
                continue
            positions = {sid: i for i, sid in enumerate(mol.atom_numbers)}
            for unit in _tetra(mol):
                anchor = unit['anchor']
                said = stated_parity(unit, ctab.atoms[positions[anchor]].parity, positions)
                drawn = mol.parity_of(anchor)
                if said and drawn:
                    both += 1
                    agree += said == drawn
    assert both > 1000, f'only {both} doubly-stated centres; the corpus stopped covering this'
    assert agree > both * 0.99, f'{agree} of {both} agree; the field is being read in a wrong frame'


def test_the_field_recovers_nothing_the_drawing_did_not_already_say(corpus):
    """A pinned negative measurement: the fixed-point search reaches every corpus centre from the
    drawing, so reading the field changes no configuration here.  What the field is kept for is the case
    the corpus has none of -- a record with no layout, covered by hand above.
    """
    for records in corpus.values():
        for record in records:
            parse = parse_v3000 if sniff_version(record, []) == V3000_STAMP else parse_v2000
            try:
                ctab = parse(record, [])
                full, _, _ = ctab.build()
                bare_ctab = parse(record, [])
                bare, _, _ = bare_ctab.build(ignore_stereo=True)
            except Exception:
                continue
            sids = list(bare.atom_numbers)
            z = {sid: a.z for sid, a in zip(sids, bare_ctab.atoms)} \
                if any(a.z for a in bare_ctab.atoms) else None
            assign_parities(bare, z, [])  # no `stated`: the drawing alone
            for sid in sids:
                assert bool(bare.parity_of(sid)) == bool(full.parity_of(sid)), \
                    f'atom {sid} of {ctab.title!r} is configured by the field alone'


def test_where_the_drawing_and_the_field_disagree_the_drawing_wins_and_it_is_logged(corpus):
    """The stored parity is the drawn one and the contradiction is named per atom; a caller has no other
    way to learn that the record disagrees with itself."""
    seen = 0
    for records in corpus.values():
        for record in records:
            parse = parse_v3000 if sniff_version(record, []) == V3000_STAMP else parse_v2000
            try:
                ctab = parse(record, [])
                mol, _, log = ctab.build()
            except Exception:
                continue
            positions = {sid: i for i, sid in enumerate(mol.atom_numbers)}
            for unit in _tetra(mol):
                anchor = unit['anchor']
                field = ctab.atoms[positions[anchor]].parity
                if field not in (1, 2):
                    continue
                said = stated_parity(unit, field, positions)
                drawn = mol.parity_of(anchor)
                if not said or said == drawn:
                    continue
                seen += 1
                assert mol.parity_of(anchor) == drawn, 'the field overwrote the drawing'
                assert any(f'atom {anchor}:' in x and 'disagree' in x for x in log), log
    assert seen, 'no record in the corpus contradicts itself; this test proved nothing'


def test_a_record_with_no_layout_still_yields_its_stated_configurations():
    """Every coordinate at the origin means no wedge can be read, so the parity column is the record's
    only statement about configuration and returning early on "no coordinates" loses all of it."""
    lines = [x for x in _TETRA_FLAT]
    for i in range(4, 8):
        lines[i] = '    0.0000    0.0000    0.0000' + lines[i][30:]
    lines[4] = lines[4][:39] + '  1' + lines[4][42:]
    mol, _, log = _build(lines)
    assert not mol.has_coordinates
    assert [mol.parity_of(u['anchor']) for u in _tetra(mol)] == [1]
    assert not log, f'nothing about this record needs reporting: {log}'


def test_a_wedge_without_coordinates_is_reported_before_the_field_is_read():
    """A wedge with no layout cannot be honoured, but it is evidence the writer meant to state a
    configuration, so the caller is told the stereo came from elsewhere."""
    lines = list(_TETRA_FLAT)
    for i in range(4, 8):
        lines[i] = '    0.0000    0.0000    0.0000' + lines[i][30:]
    lines[4] = lines[4][:39] + '  1' + lines[4][42:]
    lines[8] = '  1  2  1  1  0  0  0'
    mol, _, log = _build(lines)
    assert [mol.parity_of(u['anchor']) for u in _tetra(mol)] == [1]
    assert any('no coordinates' in x for x in log), log


def test_the_v3000_atom_cfg_is_the_same_channel_and_states_the_same_configuration():
    """The V3000 atom line's ``CFG=`` fills the same ``CtabAtom.parity`` slot as V2000's ``sss``, and
    every other test here states the field in V2000 columns.  Asserted as an equality between the two
    versions rather than as a bare ``== 1``: a sign convention applied to one version only is invisible
    to two independent single-version assertions.
    """
    two = [x for x in _TETRA_FLAT]
    for i in range(4, 8):
        two[i] = '    0.0000    0.0000    0.0000' + two[i][30:]
    two[4] = two[4][:39] + '  1' + two[4][42:]
    flat, _, _ = _build(two)
    assert not flat.has_coordinates

    three, _, log = _build(_TETRA_FLAT_V3000_CFG)
    assert not three.has_coordinates, 'the V3000 fixture has to be layout-free too'
    assert not any(three.wedges()), 'a wedge in the fixture would make the field unnecessary'
    assert [three.element_of(s) for s in three.atom_numbers] == [6, 9, 17, 35]

    said, = [three.parity_of(u['anchor']) for u in _tetra(three)]
    assert said == 1, f'the atom CFG was not read: {log}'
    assert said == [flat.parity_of(u['anchor']) for u in _tetra(flat)][0], \
        'the two versions of one field disagree about the same molecule'
    assert not log, f'nothing about this record needs reporting: {log}'


def test_a_v3000_atom_cfg_of_three_is_not_a_configuration():
    """``CFG=3`` is the spec's "either", the same statement as V2000's ``sss`` 3, and neither is a
    configuration.  Pinned separately because the reader's ``_int`` hands any integer through.
    """
    lines = [x.replace('CFG=1', 'CFG=3') for x in _TETRA_FLAT_V3000_CFG]
    assert any('CFG=3' in x for x in lines), 'the substitution has to have landed'
    mol, _, log = _build(lines)
    assert [mol.parity_of(u['anchor']) for u in _tetra(mol)] == [0]
    assert not log, f'an unstated configuration is not a complaint: {log}'


# --- resolving one centre at a time

def test_a_chain_of_centres_resolves_over_several_passes(corpus):
    """Configuring a neighbour is what makes two branches differ, so a chain resolves from the outside
    in and a single-pass reader drops its middle.  ``test/stereo.sdf`` record 198 is a five-centre
    chain, addressed by index so a change to the file fails here rather than testing another molecule.
    """
    records = corpus['stereo.sdf']
    assert len(records) == 300, 'stereo.sdf changed; the record index below is no longer meaningful'
    mol, _, log = _build(records[198])
    parities = [mol.parity_of(u['anchor']) for u in _tetra(mol)]
    assert len(parities) == 5 and all(parities), parities
    assert not any('non-stereogenic' in x for x in log), \
        'a centre that resolves on a later pass must not be reported as non-stereogenic'


def test_the_search_terminates_on_a_molecule_where_nothing_resolves():
    """The loop is bounded by the number of stereo units, so a record where no centre ever resolves
    costs one pass and not an infinite number."""
    mol, _, log = _build(_TETRA_FLAT)
    assert not any('did not settle' in x for x in log), log


def test_log_lines_are_not_repeated_once_per_pass():
    """Logging is suppressed during the search and the reasons collected afterwards against the settled
    molecule; otherwise an undetermined centre is explained once per iteration."""
    lines = list(_TETRA_FLAT)
    lines[8] = '  1  2  1  4  0  0  0'
    mol, _, log = _build(lines)
    assert len(log) == len(set(log)), log


# --- double bonds

def test_a_cis_double_bond_reads_differently_from_a_trans_one():
    """A CTfile has no field for double-bond geometry -- bond stereo 3 says only "cis or trans,
    unknown which" -- so the coordinates are the entire statement."""
    cis, _, _ = _build(_BUTENE_CIS)
    trans, _, _ = _build(_BUTENE_TRANS)
    a = [cis.parity_of(u['anchor']) for u in cis.stereo_units() if u['kind'] == SU_CIS_TRANS]
    b = [trans.parity_of(u['anchor']) for u in trans.stereo_units() if u['kind'] == SU_CIS_TRANS]
    assert a and all(a) and b and all(b), (a, b)
    assert a == [3 - x for x in b], (a, b)


def test_the_same_side_pair_gets_the_label_chython2_gives_it(oracle_session):
    """The bond sign, pinned against chython 2's *SMILES* stack: its SDF reader derives no double-bond
    geometry from coordinates, so there is nothing there to compare against.  Both stacks are asked the
    same frame-free question -- one named terminal substituent on each end.
    """
    ours = {}
    for name, lines in (('cis', _BUTENE_CIS), ('trans', _BUTENE_TRANS)):
        mol, _, _ = _build(lines)
        unit, = [u for u in mol.stereo_units() if u['kind'] == SU_CIS_TRANS and u['stereogenic']]
        assert cis_trans_parity(mol, unit) == mol.parity_of(unit['anchor']), \
            'the stored parity is not what the geometry function returns'
        ours[name] = mol.parity_of(unit['anchor'])

    # `C/C=C\C` is the same-side pair; atoms 1 and 4 are the terminal substituents in both notations
    read = dict(oracle_session.read_smiles([r'C/C=C\C', 'C/C=C/C']))
    for name, smi in (('cis', r'C/C=C\C'), ('trans', 'C/C=C/C')):
        sign = read[smi].cis_trans_sign(2, 3, 1, 4)
        assert sign is not None, f'{name}: chython 2 states no geometry for {smi}'
        theirs = 2 if sign else 1
        assert ours[name] == theirs, f'{name}: {ours[name]} against chython 2\'s {theirs}'
    assert ours['cis'] == 2, 'the same-side pair is parity 2; see cis_trans_parity'


def test_a_collinear_double_bond_layout_is_reported_rather_than_resolved():
    lines = list(_BUTENE_CIS)
    lines[7] = '    3.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0'
    mol, _, log = _build(lines)
    assert any('collinear' in x for x in log), log


# --- allenes

def _allene(mol):
    return [u for u in mol.stereo_units() if u['kind'] == SU_ALLENE and u['stereogenic']]


def test_an_allene_wedge_is_read_and_the_two_wedge_directions_disagree():
    """``_ALLENE_WEDGE`` is penta-2,3-diene drawn along x -- C2 at the origin, centre C3 at (1, 0), C4
    at (2, 0) -- with an UP wedge from the terminal C2 to its methyl C1, giving the core frame
    ``refs == (1, None, 5, None)``.  A flat drawing puts the wedged terminal's substituent plane
    perpendicular to the paper, so its in-plane offsets carry nothing and only the wedge sign does:
    collapsing it onto its own position leaves ``p0 = (0, 0, +1)``, ``p1 = (0, 0, -1)`` and the
    determinant ``-2 * (2.87 * -0.5 - 0.5 * 2.13) = +5``, a positive volume being parity 1.
    """
    up, _, log = _build(_ALLENE_WEDGE)
    unit, = _allene(up)
    assert unit['refs'] == (1, None, 5, None), unit
    assert up.parity_of(unit['anchor']) == 1
    assert not any('allene' in x for x in log), log

    down = list(_ALLENE_WEDGE)
    down[9] = '  2  1  1  6  0  0  0'
    hashed, _, _ = _build(down)
    unit, = _allene(hashed)
    assert hashed.parity_of(unit['anchor']) == 2


def test_a_mirrored_allene_drawing_inverts_the_parity():
    """The enantiomer drawn as one -- layout reflected in x, wedge kept -- rather than made by flipping
    the wedge: a reader that read only the wedge would answer the same for both."""
    mol, _, _ = _build(_ALLENE_WEDGE)
    unit, = _allene(mol)
    mirror, _, _ = _build(_ALLENE_MIRROR)
    other, = _allene(mirror)
    assert mol.parity_of(unit['anchor']) == 3 - mirror.parity_of(other['anchor'])
    assert mirror.parity_of(other['anchor']) == 2


def test_a_tetrasubstituted_allene_reads_from_a_pair_of_wedges_on_one_terminal():
    """The shape every corpus allene has: both directions of one terminal named and both wedged, one up
    one down, which is the drawing convention for "this terminal's plane is perpendicular to the paper".
    The two wedges are redundant -- the pair collapses to the same two z values a single wedge implies.
    1,3-dibromo-1,3-difluoroallene, the compound ``_inchi.pxi`` pins the core's axial sign on:
    ``-2 * (2.87 * -0.5 - 0.5 * 2.87) = +5.74``, so parity 1.
    """
    mol, _, log = _build(_ALLENE_TETRA)
    unit, = _allene(mol)
    assert unit['refs'] == (1, 3, 6, 7), unit
    assert unit['unnamed_mask'] == 0, 'this fixture has no implicit hydrogen; the holes are untested'
    assert mol.parity_of(unit['anchor']) == 1
    assert not log, log

    one, = [i for i, x in enumerate(_ALLENE_TETRA) if x.startswith('  2  3  1  6')]
    single = list(_ALLENE_TETRA)
    single[one] = '  2  3  1  0  0  0  0'
    lone, _, _ = _build(single)
    unit, = _allene(lone)
    assert lone.parity_of(unit['anchor']) == 1, 'dropping the redundant hash changed the answer'


def test_the_allene_sign_is_the_one_chython2_gives(v2_reader):
    """The axial sign pinned against chython 2, which reads allene wedges too.  The frame bridge is
    V2's ``_translate_allene_sign(c, nn, nm)``, which takes one substituent from each terminal:
    ``refs[0]`` and ``refs[2]`` are exactly that and are never holes.  True is parity 2, as in the
    tetrahedral and cis/trans cases.
    """
    for lines in (_ALLENE_WEDGE, _ALLENE_MIRROR, _ALLENE_TETRA):
        mol, _, _ = _build(lines)
        unit, = _allene(mol)
        anchor = unit['anchor']
        m2 = v2_reader('\n'.join(lines))
        assert m2.atom(anchor).stereo is not None, 'chython 2 read no configuration here'
        sign, = m2.translate([[2, anchor, unit['refs'][0], unit['refs'][2]]])
        assert sign is not None, 'chython 2 cannot express this frame, so it pins nothing'
        assert mol.parity_of(anchor) == (2 if sign else 1)


def test_every_allene_configuration_agrees_with_chython2(corpus, v2_molecules):
    """The corpus oracle for the axial sign.  All 11 configured corpus allenes are compared and the
    assertion is exact: 11 is small enough that a frame error on the two carrying a hole would hide
    inside any percentage.
    """
    agree = 0
    wrong = []
    for name, records in corpus.items():
        v2 = v2_molecules.get(name)
        if v2 is None or len(v2) != len(records):
            continue
        for n, (record, m2) in enumerate(zip(records, v2)):
            try:
                mol, _, _ = _build(record)
            except Exception:
                continue
            sids = list(mol.atom_numbers)
            nums = list(m2)
            if len(sids) != len(nums):
                continue
            to_num = dict(zip(sids, nums))
            for unit in _allene(mol):
                anchor = unit['anchor']
                ours = mol.parity_of(anchor)
                assert ours, f'{name} record {n} atom {anchor}: a drawn allene left unread'
                sign, = m2.translate([[SU_ALLENE, to_num[anchor], to_num[unit['refs'][0]],
                                       to_num[unit['refs'][2]]]])
                if sign is None:
                    continue
                if ours == (2 if sign else 1):
                    agree += 1
                else:
                    wrong.append(f'{name} record {n} atom {anchor}: {ours} against {2 if sign else 1}')
    assert not wrong, wrong
    assert agree == 11, f'{agree} allenes compared; the oracle stopped covering'


def test_a_five_carbon_cumulene_is_axial_and_read_the_same_way():
    """``SU_ALLENE`` is emitted for any odd chain and anchored on the chain's *centre*, which for five
    carbons is adjacent to neither terminal, so the reader must walk the chain: a neighbour lookup on
    the anchor works on every allene and on no longer cumulene.  Hepta-2,3,4,5-tetraene drawn like
    ``_ALLENE_WEDGE`` with two more chain atoms -- ``-2 * (4.87 * -0.5 - 0.5 * 4.13) = +9``, parity 1.
    """
    mol, _, log = _build(_CUMULENE5)
    unit, = _allene(mol)
    anchor = unit['anchor']
    assert sorted(mol.neighbors_of(anchor)) == [3, 5], 'the anchor is not the chain centre'
    assert unit['refs'] == (1, None, 7, None), unit
    assert all(r not in mol.neighbors_of(anchor) for r in unit['refs'] if r is not None), \
        'this fixture must keep the substituents off the anchor, or the walk is untested'
    assert mol.parity_of(anchor) == 1
    assert not log, log

    down = list(_CUMULENE5)
    down[11] = '  2  1  1  6  0  0  0'
    hashed, _, _ = _build(down)
    other, = _allene(hashed)
    assert hashed.parity_of(other['anchor']) == 2


def test_a_three_dimensional_allene_states_its_configuration_without_a_wedge():
    """Given a non-flat `z` map the coordinates state the configuration outright and the wedges are
    ignored.  Called directly rather than through ``assign_parities`` because the 3D path belongs to the
    XYZ and V3000-with-z callers; a flat allene lifted at one substituent is the molecule the up wedge
    draws, so the two branches must return the same parity or one has its frame reversed.
    """
    mol, _, _ = _build(_ALLENE_FLAT)
    unit, = _allene(mol)
    assert allene_parity(mol, unit) == 0, 'flat and unwedged states nothing'
    assert allene_parity(mol, unit, {1: 1.0}) == 1
    assert allene_parity(mol, unit, {1: -1.0}) == 2

    wedged, _, _ = _build(_ALLENE_WEDGE)
    other, = _allene(wedged)
    assert allene_parity(wedged, other) == allene_parity(mol, unit, {1: 1.0})


def test_a_non_stereogenic_allene_gets_no_parity_and_no_log_noise():
    """Two shapes of "not stereogenic", both drawn with a wedge.  4-methylpenta-2,3-diene duplicates one
    of a terminal's directions, so the core emits no axial unit at all; 3-ethylhexa-3,4-diene emits the
    unit and the automorphism search unmarks it.  Neither may log: an allene's anchor is the chain centre
    while a wedge's narrow end is a terminal, so "wedge on a non-stereogenic centre" cannot fire here.
    """
    mol, _, log = _build(_ALLENE_SYMMETRIC)
    assert not [u for u in mol.stereo_units() if u['kind'] == SU_ALLENE], \
        'the duplicated direction should keep this out of the unit list entirely'
    assert all(not mol.parity_of(s) for s in mol.atom_numbers)
    assert not log, log

    mol, _, log = _build(_ALLENE_NOT_STEREOGENIC)
    unit, = [u for u in mol.stereo_units() if u['kind'] == SU_ALLENE]
    assert not unit['stereogenic'], 'this fixture must keep the unit and lose the mark'
    assert not mol.parity_of(unit['anchor'])
    assert not log, log


def test_a_wedge_that_does_not_determine_the_allene_twist_is_reported():
    """Three drawings that state something without stating a configuration: a wavy bond ("up or down,
    unknown which"); a wedge on each terminal, two statements about a frame that needs one; and two
    wedges the same way on one terminal, which puts both its substituents toward the viewer.
    """
    either = list(_ALLENE_WEDGE)
    either[9] = '  2  1  1  4  0  0  0'
    mol, _, log = _build(either)
    unit, = _allene(mol)
    assert not mol.parity_of(unit['anchor'])
    assert any('either' in x for x in log), log

    both = list(_ALLENE_WEDGE)
    both[12] = '  4  5  1  1  0  0  0'
    mol, _, log = _build(both)
    unit, = _allene(mol)
    assert not mol.parity_of(unit['anchor'])
    assert any('both terminals' in x for x in log), log

    same = list(_ALLENE_TETRA)
    one, = [i for i, x in enumerate(_ALLENE_TETRA) if x.startswith('  2  3  1  6')]
    same[one] = '  2  3  1  1  0  0  0'
    mol, _, log = _build(same)
    unit, = _allene(mol)
    assert not mol.parity_of(unit['anchor'])
    assert any('the same way' in x for x in log), log


def test_a_flat_allene_drawing_states_no_configuration_and_needs_no_explanation():
    """A CTfile has no parity field for an axial unit, so a flat allene is simply an unspecified one."""
    mol, _, log = _build(_ALLENE_FLAT)
    assert _allene(mol)
    assert all(not mol.parity_of(s) for s in mol.atom_numbers)
    assert not log, log


def test_a_degenerate_allene_layout_is_reported_rather_than_resolved():
    """The far terminal's substituent drawn along the axis leaves the paper-plane half of the frame with
    no width, so the determinant is zero however far the wedge lifts the other terminal."""
    mol, _, log = _build(_ALLENE_DEGENERATE)
    assert _allene(mol), 'the fixture must keep a stereogenic allene'
    assert all(not mol.parity_of(s) for s in mol.atom_numbers)
    assert any('degenerate' in x for x in log), log


def test_an_allene_with_no_coordinates_is_left_unset_and_says_so():
    """With every atom at the origin there is no picture, and unlike a tetrahedral centre an allene has
    no second channel: no CTfile field states an axial configuration."""
    lines = list(_ALLENE_WEDGE)
    for i in range(4, 9):
        lines[i] = '    0.0000    0.0000    0.0000' + lines[i][30:]
    mol, _, log = _build(lines)
    assert not mol.has_coordinates
    assert all(not mol.parity_of(s) for s in mol.atom_numbers)
    assert any('no coordinates' in x for x in log), log


def test_an_allene_survives_a_molfile_round_trip():
    """Two round trips exercising different halves of the writer: the first re-emits the wedges the file
    was drawn with, the second must *choose* one, the parity having been set onto a flat layout -- the
    case of any molecule that did not come from a CTfile.
    """
    mol, _, _ = _build(_ALLENE_WEDGE)
    unit, = _allene(mol)
    target = mol.parity_of(unit['anchor'])
    lines, _ = emit_v2000(mol)
    again, _, _ = _build(lines)
    other, = _allene(again)
    assert again.parity_of(other['anchor']) == target

    for parity in (1, 2):
        bare, _, _ = _build(_ALLENE_FLAT)
        unit, = _allene(bare)
        with bare.edit():
            bare.set_parity(unit['anchor'], parity)
        wedges, log = wedges_for_write(bare)
        assert wedges, log
        lines, _ = emit_v2000(bare)
        back, _, _ = _build(lines)
        other, = _allene(back)
        assert back.parity_of(other['anchor']) == parity


def test_a_chosen_allene_wedge_lands_on_a_terminal_and_not_on_the_centre():
    """An axial unit is anchored on the chain's *centre* while the wedge that states it sits on a bond
    from a *terminal*: the tetrahedral rule would aim at the anchor, whose every bond is double and can
    carry no wedge.  Checked on the five-carbon cumulene too, where the terminal is not adjacent to it.
    """
    flat_cumulene = _CUMULENE5[:11] + ['  2  1  1  0  0  0  0'] + _CUMULENE5[12:]
    for lines, expected in ((_ALLENE_FLAT, {2, 4}), (flat_cumulene, {2, 6})):
        mol, _, _ = _build(lines)
        assert not any(mol.wedges()), 'wedges_for_write returns existing wedges; nothing is chosen'
        unit, = _allene(mol)
        with mol.edit():
            mol.set_parity(unit['anchor'], 1)
        (narrow, wide, wedge), = wedges_for_write(mol)[0]
        assert narrow in expected, f'the narrow end {narrow} is not a chain terminal'
        assert mol.order_of(narrow, wide) == 1, 'the allene writer picks a single bond, and one exists here'
        assert wedge in (1, 2)


# --- atropisomers

def _atropo(mol):
    return [u for u in mol.stereo_units() if u['kind'] == SU_ATROPISOMER and u['stereogenic']]


def test_an_atropisomer_wedge_is_read_and_the_two_wedge_directions_disagree():
    """``_ATROPO_WEDGE`` is 2-chloro-2'-fluorobiphenyl as two unit hexagons along x, pivots at (1, 0) and
    (2, 0), with an UP wedge from pivot 1 to its ortho carbon 2.  The frame is the core's ``refs``, each
    pivot's two ring directions in slot order, so ``(2, 6, 9, 13)``, and the wedged pivot collapses onto
    its own position exactly as an allene terminal does: ``p0 = (1, 0, +1)``, ``p1 = (1, 0, -1)``, the far
    pair at their drawn ``(2.5, +-0.866, 0)``, and the determinant ``-2 * (1.5 * -0.866 - 0.866 * 1.5) =
    +5.196``, a positive volume being parity 1.
    """
    up, _, log = _build(_ATROPO_WEDGE)
    unit, = _atropo(up)
    assert unit['refs'] == (2, 6, 9, 13), unit
    assert not unit['unnamed_mask'], 'a pivot has no hydrogen, so no ref of this kind is ever a hole'
    assert up.parity_of(unit['anchor']) == 1
    assert not log, log

    down = list(_ATROPO_WEDGE)
    down[18] = '  1  2  1  6  0  0  0'
    hashed, _, _ = _build(down)
    unit, = _atropo(hashed)
    assert hashed.parity_of(unit['anchor']) == 2


def test_a_mirrored_atropisomer_drawing_inverts_the_parity():
    """The enantiomer drawn as one -- layout reflected in x, wedge kept -- rather than made by flipping the
    wedge: a reader that read only the wedge would answer the same for both."""
    mol, _, _ = _build(_ATROPO_WEDGE)
    unit, = _atropo(mol)
    mirror, _, _ = _build(_ATROPO_MIRROR)
    other, = _atropo(mirror)
    assert mol.parity_of(unit['anchor']) == 3 - mirror.parity_of(other['anchor'])
    assert mirror.parity_of(other['anchor']) == 2


def test_either_ring_bond_of_a_pivot_carries_the_statement():
    """A drawing puts the wedge on whichever of the pivot's two ring bonds suits it, and the two state
    opposite configurations.  Here the other one is C1=C6, a *double* bond in this Kekule form: the wedge
    is read there because it is where a package drew it, the bond order not entering the frame.  Wedged
    both ways at once the two agree, and that is one statement rather than a contradiction.
    """
    other = list(_ATROPO_WEDGE)
    other[18] = '  1  2  1  0  0  0  0'
    other[23] = '  1  6  2  1  0  0  0'
    mol, _, log = _build(other)
    unit, = _atropo(mol)
    assert mol.order_of(1, 6) == 2
    assert mol.parity_of(unit['anchor']) == 2
    assert not log, log

    both = list(_ATROPO_WEDGE)
    both[23] = '  1  6  2  6  0  0  0'
    mol, _, log = _build(both)
    unit, = _atropo(mol)
    assert mol.parity_of(unit['anchor']) == 1
    assert not log, log


def test_a_wedge_at_the_far_pivot_is_read_too_and_the_layout_decides_which_way():
    """Either pivot may carry the axial statement, so the reader collapses whichever end is wedged.  That
    the frame is geometric and not a wedge lookup is what the second half shows: which of the two answers
    an UP wedge gives depends on which side of the axis the *other* pivot's first ring direction is drawn
    on, so reflecting ring B alone brings the two placements into agreement.
    """
    far = list(_ATROPO_WEDGE)
    far[18] = '  1  2  1  0  0  0  0'
    far[26] = '  8  9  1  1  0  0  0'
    mol, _, log = _build(far)
    unit, = _atropo(mol)
    assert mol.parity_of(unit['anchor']) == 2
    assert not log, log

    flipped_far = _ATROPO_RING_B_FLIPPED[:18] + ['  1  2  1  0  0  0  0'] + _ATROPO_RING_B_FLIPPED[19:26] \
        + ['  8  9  1  1  0  0  0'] + _ATROPO_RING_B_FLIPPED[27:]
    near, _, _ = _build(_ATROPO_RING_B_FLIPPED)
    far, _, _ = _build(flipped_far)
    assert near.parity_of(_atropo(near)[0]['anchor']) == far.parity_of(_atropo(far)[0]['anchor']) == 2


def test_a_flat_biaryl_states_no_configuration_and_needs_no_explanation():
    """No CTfile field states an axial configuration, so a biaryl drawn without a wedge on the axis is
    simply an unspecified one -- which all but 88 of the 4368 axes in a 119,534-molfile sample of a
    production corpus are, and a log line each would be 4280 lines of noise."""
    mol, _, log = _build(_ATROPO_FLAT)
    assert _atropo(mol)
    assert all(not mol.parity_of(s) for s in mol.atom_numbers)
    assert not log, log


def test_a_contradictory_atropisomer_drawing_is_reported_and_left_unset():
    """Three ways a drawing says two things at once, each named by its own line: a wedge at each pivot,
    both ring bonds of one pivot lifted the same way, and a bond drawn as either."""
    both_pivots = list(_ATROPO_WEDGE)
    both_pivots[26] = '  8  9  1  1  0  0  0'
    mol, _, log = _build(both_pivots)
    unit, = _atropo(mol)
    assert not mol.parity_of(unit['anchor'])
    assert any('both pivots' in x for x in log), log

    same = list(_ATROPO_WEDGE)
    same[23] = '  1  6  2  1  0  0  0'
    mol, _, log = _build(same)
    unit, = _atropo(mol)
    assert not mol.parity_of(unit['anchor'])
    assert any('the same way' in x for x in log), log

    wavy = list(_ATROPO_WEDGE)
    wavy[18] = '  1  2  1  4  0  0  0'
    mol, _, log = _build(wavy)
    unit, = _atropo(mol)
    assert not mol.parity_of(unit['anchor'])
    assert any('either' in x for x in log), log


def test_a_corpus_atropisomer_is_read_from_its_drawing(corpus):
    """``stereo.sdf`` record 54 is 2,2'-dibromobiphenyl-6,6'-dicarboxylic acid, addressed by index so a
    change to the file fails here.  Its axis is one of the 6 the corpus coverage guard counts."""
    records = corpus['stereo.sdf']
    assert len(records) == 300, 'stereo.sdf changed; the record index below is no longer meaningful'
    mol, _, log = _build(records[54])
    unit, = _atropo(mol)
    assert mol.parity_of(unit['anchor']) == 1
    assert not log, log


def test_an_atropisomer_survives_a_molfile_round_trip():
    """Two round trips exercising different halves of the writer, as for the allene: the first re-emits the
    wedge the file was drawn with, the second must *choose* one, the parity having been set onto a flat
    layout.  A chosen wedge lands on a bond from a pivot, the anchor being a pivot itself here.
    """
    mol, _, _ = _build(_ATROPO_WEDGE)
    unit, = _atropo(mol)
    target = mol.parity_of(unit['anchor'])
    assert wedges_for_write(mol)[0] == [(1, 2, 1)], 'the drawn wedge must be returned untouched'
    lines, _ = emit_v2000(mol)
    again, _, _ = _build(lines)
    other, = _atropo(again)
    assert again.parity_of(other['anchor']) == target

    for parity in (1, 2):
        bare, _, _ = _build(_ATROPO_FLAT)
        unit, = _atropo(bare)
        with bare.edit():
            bare.set_parity(unit['anchor'], parity)
        (narrow, wide, wedge), = wedges_for_write(bare)[0]
        assert (narrow, wide) in ((1, 2), (1, 6), (8, 9), (8, 13)), 'not a pivot ring bond'
        assert wedge in (1, 2)
        lines, _ = emit_v2000(bare)
        back, _, _ = _build(lines)
        other, = _atropo(back)
        assert back.parity_of(other['anchor']) == parity


# --- choosing a wedge

def test_the_wedges_a_file_was_drawn_with_are_returned_untouched(corpus):
    """A round trip must not redraw: re-derived wedges land on different bonds, which is a different
    picture of the same molecule."""
    checked = 0
    for records in corpus.values():
        for record in records:
            try:
                mol, _, _ = _build(record)
            except Exception:
                continue
            existing = list(mol.wedges())
            if not existing:
                continue
            wedges, _ = wedges_for_write(mol)
            assert wedges == existing
            checked += 1
    assert checked > 50, checked


def test_a_chosen_wedge_reads_back_as_the_parity_it_was_chosen_for():
    """The writer does not compute the sign; it asks the reader which wedge would read back as the
    stored parity, so a round trip cannot invert one."""
    for wedge_line in ('  1  2  1  1  0  0  0', '  1  2  1  6  0  0  0'):
        lines = list(_TETRA_FLAT)
        lines[8] = wedge_line
        mol, _, _ = _build(lines)
        target = {u['anchor']: mol.parity_of(u['anchor']) for u in _tetra(mol)}
        assert any(target.values())

        # a copy with the same geometry and no wedges, so the writer has to choose
        bare, _, _ = _build(_TETRA_FLAT)
        with bare.edit():
            for anchor, parity in target.items():
                bare.set_parity(anchor, parity)
        wedges, log = wedges_for_write(bare)
        assert wedges, log
        with bare.edit():
            for narrow, wide, w in wedges:
                bare.set_wedge(narrow, wide, w)
        for unit in _tetra(bare):
            assert tetrahedral_parity(bare, unit) == target[unit['anchor']]


def test_no_wedges_are_written_for_a_molecule_with_no_coordinates():
    """Without a picture a wedge says nothing: emitting one produces a file whose stereo depends on
    whatever layout a later tool invents."""
    mol = MoleculeContainer()
    with mol.edit():
        c = mol.add_atom('C')
        for element in ('F', 'Cl', 'Br', 'I'):
            mol.add_bond(c, mol.add_atom(element), 1)
    with mol.edit():
        mol.set_parity(next(iter(mol.atom_numbers)), 1)
    wedges, log = wedges_for_write(mol)
    assert wedges == []
    assert any('no coordinates' in x for x in log), log


def test_an_unconfigured_molecule_needs_no_wedges_and_says_nothing():
    mol, _, _ = _build(_TETRA_FLAT)
    assert wedges_for_write(mol) == ([], [])


# --- a double-bond configuration no coordinate format can carry

def _cis_trans_lines(log):
    return [x for x in log if 'a double-bond configuration' in x]


def _configured_double_bonds(mol):
    return sum(1 for u in mol.stereo_units()
               if u['kind'] == SU_CIS_TRANS and mol.parity_of(u['anchor']))


def test_a_double_bond_configuration_with_no_coordinates_is_reported_by_both_writers():
    """A CTAB states double-bond geometry only through the coordinates, so a molecule with no layout --
    ``C/C=C/C`` from SMILES carries the configuration with every coordinate unset -- has nowhere to put
    it.  Prefixed ``unsupported: `` because the input is fine and we are the limitation; the caller can
    act on it by laying the molecule out first.
    """
    mol = read_smiles('C/C=C/C')
    assert _configured_double_bonds(mol) == 1, 'the fixture states no configuration'
    assert not mol.has_coordinates, 'and the point is that it states one without a drawing'

    for emit in (emit_v2000, emit_v3000):
        log = []
        emit(mol, log=log)
        lines = _cis_trans_lines(log)
        assert len(lines) == 1, log
        assert str(lines[0]).startswith('unsupported: stereo: atom '), lines
        assert str(lines[0]).endswith('a double-bond configuration on a molecule with no coordinates '
                                      'is not written'), lines


def test_every_writer_that_calls_the_chooser_reports_the_same_lines():
    """Four writers may state this loss and must state it identically, since a caller screening for
    ``unsupported:`` sees one stream whatever dialect produced it.  The set is over *wordings* and not
    over writers: CML and MRV have ``<bondStereo>C``/``T``, a channel the coordinates are not, and say so
    per anchor through ``cis_trans_stated``, so asserting all four logs equal would pin the CTAB
    limitation onto them.  The XML dialects are imported here because the shared chooser is under test.
    """
    from ...xml._cml import record_from_molecule as cml_record
    from ...xml._mrv import record_from_molecule as mrv_record

    mol = read_smiles('C/C=C/C=C/C')
    assert _configured_double_bonds(mol) == 2, 'the fixture does not carry two configurations'

    logs = {}
    for name, call in (('v2000', lambda log: emit_v2000(mol, log=log)),
                       ('v3000', lambda log: emit_v3000(mol, log=log)),
                       ('cml', lambda log: cml_record(mol, log=log)),
                       ('mrv', lambda log: mrv_record(mol, log=log))):
        log = []
        call(log)
        logs[name] = sorted(_cis_trans_lines(log))

    # One wording: every line, stripped of the atom it names -- the fourth field of
    # `unsupported: stereo: atom N: ...` -- must be the same sentence.
    assert len({str(x).split(': ', 3)[-1] for lines in logs.values() for x in lines}) == 1, logs
    # The two CTAB versions, whose only channel is the coordinates, report both configurations.
    assert logs['v2000'] == logs['v3000'], logs
    assert len(logs['v2000']) == 2, logs['v2000']
    # Both descriptors are writable as a bare letter here -- every terminal carries one substituent --
    # so neither XML dialect loses anything.
    assert logs['cml'] == [] and logs['mrv'] == [], logs


def test_a_dialect_with_a_letter_reports_only_the_anchors_the_letter_cannot_carry():
    """The per-anchor half of the rule above, which a per-format flag could not express.  MRV's
    ``<bondStereo>`` names no reference atoms, so its letter states a configuration only where a terminal
    carries no second substituent; CML's takes an ``atomRefs4`` and always can.  The fixture holds one of
    each, so MRV must report exactly one loss and CML none.
    """
    from ...xml._cml import record_from_molecule as cml_record
    from ...xml._mrv import record_from_molecule as mrv_record

    mol = read_smiles('C/C=C/C=C(\\C)CC')
    assert _configured_double_bonds(mol) == 2, 'the fixture does not carry two configurations'

    log = []
    mrv_record(mol, log=log)
    lines = _cis_trans_lines(log)
    assert len(lines) == 1, log

    log = []
    cml_record(mol, log=log)
    assert _cis_trans_lines(log) == [], log


def test_the_report_is_one_line_per_unit_and_names_the_anchor():
    """One line per configuration lost, naming the anchor, since an aggregate count cannot be acted
    on."""
    mol = read_smiles('C/C=C/C=C/C')
    log = []
    emit_v2000(mol, log=log)
    lines = _cis_trans_lines(log)
    assert len(lines) == 2, lines
    anchors = {mol.parity_of(u['anchor']) and u['anchor'] for u in mol.stereo_units()
               if u['kind'] == SU_CIS_TRANS}
    assert {int(str(x).split('atom ')[1].split(':')[0]) for x in lines} == anchors, (lines, anchors)


def test_a_stereogenic_double_bond_that_states_nothing_is_not_a_loss():
    """``CC=CC`` has a stereogenic double bond and no configuration on it, so nothing is dropped.
    Counting stereogenic units rather than configured ones reports a loss on every unspecified alkene."""
    mol = read_smiles('CC=CC')
    assert any(u['kind'] == SU_CIS_TRANS and u['stereogenic'] for u in mol.stereo_units()), \
        'the fixture has no stereogenic double bond, so it cannot arm this'
    assert _configured_double_bonds(mol) == 0
    for emit in (emit_v2000, emit_v3000):
        log = []
        emit(mol, log=log)
        assert not _cis_trans_lines(log), log


def test_a_drawn_double_bond_configuration_is_written_and_not_reported(corpus):
    """A record that arrived with a layout expresses its double bonds in the coordinates, which are
    re-emitted, so nothing is lost and nothing is said.  Taken from the corpus by searching rather than
    by index, a hand-picked record number stopping testing what it was chosen for when the file changes.
    """
    found = 0
    for name, records in corpus.items():
        for record in records:
            try:
                mol, _, _ = _build(record)
            except Exception:
                continue
            if not mol.has_coordinates or not _configured_double_bonds(mol):
                continue
            found += 1
            log = []
            written = emit_v2000(mol, log=log)
            assert not _cis_trans_lines(log), (name, log)
            back, _, _ = _build(written[0] if isinstance(written, tuple) else written)
            assert _configured_double_bonds(back) == _configured_double_bonds(mol), \
                f'{name}: the drawn configuration did not survive the round trip'
            if found == 3:
                return
    assert found, 'no corpus record carries a drawn double-bond configuration, so this proves nothing'


def test_the_tetrahedral_report_is_unchanged_and_the_two_do_not_borrow_each_others_count():
    """A coordinate-free tetrahedral molecule gets its own line with its own count and no double-bond
    line.  Widening ``configured`` to include ``SU_CIS_TRANS`` is the tempting wrong fix: that list
    drives wedge selection and a cis/trans unit wants no wedge."""
    mol = MoleculeContainer()
    with mol.edit():
        c = mol.add_atom('C')
        for element in ('F', 'Cl', 'Br', 'I'):
            mol.add_bond(c, mol.add_atom(element), 1)
    with mol.edit():
        mol.set_parity(next(iter(mol.atom_numbers)), 1)
    wedges, log = wedges_for_write(mol)
    assert wedges == []
    assert any('1 configured stereocentre(s) but no coordinates' in x for x in log), log
    assert not _cis_trans_lines(log), log


# --- the non-geometric channel: a letter instead of a drawing

def _cis_trans_unit(mol):
    unit, = (u for u in mol.stereo_units() if u['kind'] == SU_CIS_TRANS and mol.parity_of(u['anchor']))
    return unit


def test_the_letter_and_the_parity_are_one_mapping_in_both_directions():
    """``C`` is parity 2 and ``T`` is parity 1, asserted where the convention lives rather than in a
    dialect.  The letter/parity pair must be the identity; a global inversion of the convention passes
    this, so the anchoring statement is the SMILES -- a third stack's spelling of the same fact,
    ``C/C=C/C`` being trans.
    """
    for smiles, letter in (('C/C=C/C', 'T'), ('C/C=C\\C', 'C')):
        mol = read_smiles(smiles)
        unit = _cis_trans_unit(mol)
        got, frame = cis_trans_letter(mol, unit)
        assert got == letter, smiles
        assert stated_cis_trans(mol, unit, got) == mol.parity_of(unit['anchor']), smiles
        # And with the frame stated, which is the channel CML uses: same answer either way.
        assert stated_cis_trans(mol, unit, got, refs=frame) == mol.parity_of(unit['anchor']), smiles


def test_a_frame_naming_the_other_substituent_reads_as_the_opposite_configuration():
    """``C`` about one pair of substituents is ``T`` about the pair that swaps one end, so a reader that
    accepted an ``atomRefs4`` and ignored it would be right half the time.  3-methyl-2-pentene's far
    terminal carries a methyl and an ethyl, so there are two frames to tell apart."""
    mol = read_smiles('C/C=C(\\C)CC')
    unit = _cis_trans_unit(mol)
    letter, frame = cis_trans_letter(mol, unit, framed=True)
    near, anchor, partner, far = frame
    other, = (n for n in mol.neighbors_of(partner) if n not in (anchor, far))

    parity = stated_cis_trans(mol, unit, letter, refs=frame)
    assert parity == mol.parity_of(anchor)
    assert stated_cis_trans(mol, unit, letter, refs=(near, anchor, partner, other)) == 3 - parity


def test_a_bare_letter_is_refused_in_both_directions_where_it_would_be_ambiguous():
    """A bare letter names no reference atoms, so on a terminal carrying two substituents it does not say
    which pair is cis.  :func:`stated_cis_trans` will not read one and :func:`cis_trans_letter` with
    ``framed=False`` will not write one -- one test called twice, so a dialect with nowhere to put a frame
    cannot state more than it can take back.  ``framed=True`` is unaffected.
    """
    mol = read_smiles('C/C=C(\\C)CC')
    unit = _cis_trans_unit(mol)
    letter, _ = cis_trans_letter(mol, unit, framed=True)

    assert cis_trans_letter(mol, unit, framed=False) is None
    log = []
    assert stated_cis_trans(mol, unit, letter, log=log) == 0
    assert len(log) == 1 and 'not stated' in log[0], log


def test_the_frame_may_be_written_from_either_terminal():
    """``a b c d`` and ``d c b a`` are the same statement -- the quadruple is an ordered path through the
    bond -- and a document is free to write either."""
    mol = read_smiles('C/C=C/C')
    unit = _cis_trans_unit(mol)
    letter, frame = cis_trans_letter(mol, unit)
    parity = mol.parity_of(unit['anchor'])
    assert stated_cis_trans(mol, unit, letter, refs=frame) == parity
    assert stated_cis_trans(mol, unit, letter, refs=tuple(reversed(frame))) == parity


def test_a_frame_that_is_not_four_atoms_or_names_a_stranger_is_not_read():
    """A malformed quadruple is reported and dropped rather than guessed at, leaving the drawing, if
    there is one, as the only statement."""
    mol = read_smiles('C/C=C/C')
    unit = _cis_trans_unit(mol)
    letter, frame = cis_trans_letter(mol, unit)
    for refs in (frame[:3], frame + frame, (frame[0], frame[1], frame[2], 999)):
        log = []
        assert stated_cis_trans(mol, unit, letter, refs=refs, log=log) == 0, refs
        assert len(log) == 1, (refs, log)


def test_a_letter_is_offered_for_every_configured_double_bond_and_no_others():
    """What :func:`cis_trans_for_write` omits is what no dialect can state: an unconfigured stereogenic
    double bond, and -- with ``framed`` false -- a bond whose configuration a bare letter cannot carry."""
    mol = read_smiles('C/C=C/C=C(\\C)CC.CC=CC')
    configured = {u['anchor'] for u in mol.stereo_units()
                  if u['kind'] == SU_CIS_TRANS and mol.parity_of(u['anchor'])}
    assert len(configured) == 2, 'the fixture does not carry two configurations'
    assert {a for a, _, _ in cis_trans_for_write(mol, framed=True)} == configured
    bare = {a for a, _, _ in cis_trans_for_write(mol, framed=False)}
    assert len(bare) == 1 and bare < configured


# --- wedges for a given layout

def _undrawn(lines, wedge_line, index=8):
    """The fixture's configuration, on a copy of the fixture that carries no wedge.

    Two reads rather than one: the writer only chooses when the molecule is undrawn, and reading it
    drawn is how the expected parity is obtained without writing one down by hand.
    """
    drawn = list(lines)
    drawn[index] = wedge_line
    mol, _, _ = _build(drawn)
    target = {u['anchor']: mol.parity_of(u['anchor']) for u in _tetra(mol)}
    assert any(target.values()), 'the fixture as drawn states no configuration'

    bare, _, _ = _build(lines)
    with bare.edit():
        for anchor, parity in target.items():
            bare.set_parity(anchor, parity)
    return bare, target


def _plane_of(mol, f):
    """``{stable id: f(x, y)}`` -- a layout derived from the molecule's own, as a renderer's would be."""
    return {n: f(*mol.xy_of(n)) for n in mol.atom_numbers}


def test_a_supplied_plane_and_not_the_stored_coordinates_decides_the_wedge():
    """Which way a wedge points is a property of the drawing, so a renderer laying a molecule out into a
    temporary must be able to ask for that layout's wedges.  Reflecting the plane inverts its handedness,
    so the same bond carries the opposite wedge -- a relative check, which is what a function ignoring
    its `plane` argument could not pass."""
    mol, _ = _undrawn(_TETRA_FLAT, '  1  2  1  1  0  0  0')
    stored, _ = wedges_for_write(mol)
    mirrored, _ = wedges_for_write(mol, plane=_plane_of(mol, lambda x, y: (x, -y)))
    assert stored and mirrored
    # same bond -- the reflection does not change which bond reads most cleanly ...
    assert [(a, b) for a, b, _ in stored] == [(a, b) for a, b, _ in mirrored]
    # ... and the opposite code, which is the whole of the statement
    assert [w for _, _, w in stored] != [w for _, _, w in mirrored]


def test_wedges_chosen_for_a_plane_read_back_as_the_parity_against_that_plane():
    """The round trip run against the supplied layout: a wedge set self-consistent against coordinates
    nobody is drawing states the wrong configuration in the picture that gets drawn."""
    mol, target = _undrawn(_TETRA_FLAT, '  1  2  1  1  0  0  0')
    # a transposition, which is a reflection and so genuinely a different drawing
    plane = _plane_of(mol, lambda x, y: (y, x))
    wedges, log = wedges_for_write(mol, plane=plane)
    assert wedges, log

    with mol.edit():
        for n, (x, y) in plane.items():
            mol.set_xy(n, x, y)
        for narrow, wide, w in wedges:
            mol.set_wedge(narrow, wide, w)
    for unit in _tetra(mol):
        assert tetrahedral_parity(mol, unit) == target[unit['anchor']]


def test_a_plane_makes_wedges_choosable_for_a_molecule_that_has_no_coordinates():
    """With no stored coordinates every sign is 0 and the configuration is dropped, so a renderer that
    has just computed a layout must be able to hand it over."""
    mol = MoleculeContainer()
    with mol.edit():
        c = mol.add_atom('C')
        ids = [mol.add_atom(e) for e in ('F', 'Cl', 'Br', 'I')]
        for i in ids:
            mol.add_bond(c, i, 1)
    with mol.edit():
        mol.set_parity(c, 1)
    assert not mol.has_coordinates

    wedges, log = wedges_for_write(mol)
    assert wedges == [] and any('no coordinates' in x for x in log), log

    plane = dict(zip([c] + ids, [(0.0, 0.0), (0.0, 1.0), (0.87, 0.5), (0.87, -0.5), (0.0, -1.0)]))
    wedges, log = wedges_for_write(mol, plane=plane)
    assert wedges, log
    # The molecule is not written to: a chooser that stored the plane on the way past would leave the
    # caller's molecule claiming coordinates it never had.  `_Planar` is a proxy for this reason.
    assert not mol.has_coordinates


# --- which bond gets the wedge

@fixture(scope='module')
def drawings(root):
    """``{title: record lines}`` for ``test/wedge_stereo.sdf``.

    Public textbook stereochemistry laid out in 2d by an outside tool and committed, so the layout these
    tests judge is fixed -- the metrics below are pinned numbers and a layout that moved between runs
    would make them meaningless.  The fused and bridged polycyclics are the point: a steroid ring
    junction has every heavy neighbour in a ring, so it tests the fallback from the acyclic preference.
    """
    out = {}
    with (root / 'test' / 'wedge_stereo.sdf').open(encoding='utf8', errors='replace') as f:
        for record in split_records(f):
            out[(record[0] or '').strip()] = record
    assert len(out) > 80, len(out)
    return out


def _undrawn_molecule(record):
    """The record's molecule with its configurations kept and its wedges cleared -- `wedges_for_write`
    returns stored wedges untouched, so the outside tool's have to go before it will choose anything.
    """
    mol, _, _ = _build(record)
    existing = list(mol.wedges())
    if existing:
        with mol.edit():
            for narrow, wide, _ in existing:
                mol.set_wedge(narrow, wide, WEDGE_NONE)
    return mol


def _configured(mol):
    return [u for u in mol.stereo_units()
            if u['kind'] == SU_TETRA and mol.parity_of(u['anchor'])]


def _adjacent_pairs(wedges):
    """Unordered pairs of chosen wedges that share an atom.  Counted over atoms and not over bonds: a
    wedge whose wide end is another's narrow end, and two wedges pointing at one atom, read alike.
    """
    return [(i, j) for i in range(len(wedges)) for j in range(i + 1, len(wedges))
            if set(wedges[i][:2]) & set(wedges[j][:2])]


def _near_collinear_wedges(mol, wedges):
    """Chosen wedges whose bond has another within 30 degrees of it at the narrow end: a wedge is
    attributed to a bond by lying along it, so two bonds closing at the stereocentre leave the triangle
    claiming both, stating opposite configurations.  Measured through the module's own predicate.
    """
    return [(n, w) for n, w, _ in wedges if _near_collinear(mol, n, w)]


def _quality(records):
    """The objective metric over a corpus.  ``round_trip_failures`` is the gate and the rest are the
    goal, in that order: a prettier wedge set encoding the wrong configuration is a regression.
    """
    m = dict(molecules=0, centres=0, wedges=0, ring_wedges=0, adjacent_pairs=0,
             near_collinear=0, unencoded=0, round_trip_failures=0)
    for record in records:
        mol = _undrawn_molecule(record)
        configured = _configured(mol)
        if not configured:
            continue
        target = {u['anchor']: mol.parity_of(u['anchor']) for u in configured}
        m['molecules'] += 1
        m['centres'] += len(configured)

        wedges, _ = wedges_for_write(mol)
        m['wedges'] += len(wedges)
        m['ring_wedges'] += sum(1 for a, b, _ in wedges if mol.bond_in_ring(a, b))
        m['adjacent_pairs'] += len(_adjacent_pairs(wedges))
        m['near_collinear'] += len(_near_collinear_wedges(mol, wedges))
        m['unencoded'] += len(target) - len({a for a, _, _ in wedges})

        # the gate: draw what was chosen, then read it back with the reader
        with mol.edit():
            for narrow, wide, w in wedges:
                mol.set_wedge(narrow, wide, w)
        if any(tetrahedral_parity(mol, u) != target[u['anchor']]
               for u in mol.stereo_units()
               if u['kind'] == SU_TETRA and u['anchor'] in target):
            m['round_trip_failures'] += 1
    return m


def _wedge_at(mol, anchor):
    """The wide end and code of the wedge chosen for `anchor`, or None."""
    wedges, _ = wedges_for_write(mol)
    for narrow, wide, w in wedges:
        if narrow == anchor:
            return wide, w
    return None


def test_an_acyclic_bond_is_preferred_to_a_ring_bond(drawings):
    """A ring's drawn bonds are what the viewer reads as its plane, so a wedge on one asks for a ring
    atom to be lifted out of a plane the same picture insists is flat -- twice over on a fused bond.
    Menthol's isopropyl-bearing ring carbon is the case: no neighbour of it is terminal, so the degree
    tie-break below cannot separate them and only the ring test can.
    """
    mol = _undrawn_molecule(drawings['menthol'])
    for unit in _configured(mol):
        anchor = unit['anchor']
        acyclic = [r for r in unit['refs'] if r is not None and not mol.bond_in_ring(anchor, r)]
        if not acyclic:
            continue
        wide, _ = _wedge_at(mol, anchor)
        assert not mol.bond_in_ring(anchor, wide), \
            f'centre {anchor} took the ring bond to {wide} with {acyclic} available'


def test_no_centre_takes_a_ring_bond_while_an_acyclic_one_is_free(drawings):
    """The same rule as a sweep: the single-molecule test above passes on an implementation that happens
    to order menthol's neighbours luckily."""
    offenders = []
    for title, record in sorted(drawings.items()):
        mol = _undrawn_molecule(record)
        wedges, _ = wedges_for_write(mol)
        chosen = {a: b for a, b, _ in wedges}
        for unit in _configured(mol):
            anchor = unit['anchor']
            if anchor not in chosen or not mol.bond_in_ring(anchor, chosen[anchor]):
                continue
            if any(r is not None and not mol.bond_in_ring(anchor, r) for r in unit['refs']):
                offenders.append((title, anchor))
    assert not offenders, f'{len(offenders)} avoidable ring wedges: {offenders[:10]}'


def test_a_terminal_neighbour_is_preferred_to_a_branching_one(drawings):
    """The tie-break among acyclic bonds, decided by what lies past the wide end: nothing is drawn
    beyond a leaf, so the lift is unambiguous, while past a branch point the whole substituent is
    implicitly lifted too.  So degree 1 first.
    """
    mol = _undrawn_molecule(drawings['menthol'])
    for unit in _configured(mol):
        anchor = unit['anchor']
        terminal = [r for r in unit['refs'] if r is not None and mol.degree_of(r) == 1
                    and not mol.bond_in_ring(anchor, r)]
        if not terminal:
            continue
        wide, _ = _wedge_at(mol, anchor)
        assert wide in terminal, f'centre {anchor} chose {wide} over the terminal {terminal}'


def test_an_explicit_hydrogen_is_the_preferred_wide_end():
    """Among bonds equally good on every other count the hydrogen wins, which is the conventional
    drawing.  An *implicit* hydrogen is no candidate at all: it appears in the core's frame as a hole and
    a CTfile wedge is a pair of atom-block entries, so there is nothing to draw the wide end to.
    """
    mol, _ = _undrawn(_TETRA_WITH_H, '  1  3  1  1  0  0  0')
    anchor = next(u['anchor'] for u in _tetra(mol))
    wide, _ = _wedge_at(mol, anchor)
    assert mol.element_of(wide) == 1, f'chose atom {wide}, element {mol.element_of(wide)}'


def test_two_wedges_never_share_an_atom_when_they_need_not(drawings):
    """A shared atom is simultaneously "lifted" by one wedge and the base of another, and a reader cannot
    take both claims at once.  Norbornanol is the case: its carbinol carbon neighbours a bridgehead, so
    the naive choice aims the bridgehead's wedge at an atom that already carries one.
    """
    mol = _undrawn_molecule(drawings['norbornanol'])
    wedges, log = wedges_for_write(mol)
    assert not _adjacent_pairs(wedges), f'{[(wedges[i], wedges[j]) for i, j in _adjacent_pairs(wedges)]}\n{log}'


def test_adjacent_wedges_fall_across_the_corpus(drawings):
    """The sweep behind the test above.  Not zero: a centre whose every drawable bond leads to another
    centre's wedge takes one anyway rather than lose its configuration, and one pair is bought on purpose
    to remove an unattributable mark -- see `_EXPECTED_QUALITY`."""
    m = _quality(list(drawings.values()))
    assert m['adjacent_pairs'] <= _EXPECTED_QUALITY['adjacent_pairs'], m


# --- a wedge has to be attributable, and drawable

def test_a_multiple_bond_is_never_asked_to_carry_a_wedge():
    """``sss`` in the bond block is the wedge -- 1 up, 6 down, 4 either -- only on a *single* bond; on a
    double bond the same column is cis/trans (0 "use the coordinates", 3 "either"), so a 1 written there
    is an out-of-domain value and the configuration is lost.  A sulfoxide's S=O is terminal, acyclic and
    often the longest bond at the sulfur, so it wins every term of `_draw_cost`.  Both stable-id orders
    are tested: with the methyl lower the honest bond wins by accident.
    """
    for smiles_ish, expect_multiple_first in (('O=[S@](C)CC', True), ('C[S@](=O)CC', False)):
        mol = read_smiles(smiles_ish)
        ids = sorted(mol)
        # every bond exactly one unit long, so no length tie-break can decide this
        h = 3 ** .5 / 2
        with mol.edit():
            for n, (x, y) in zip(ids, [(0., 1.), (0., 0.), (-h, -.5), (h, -.5), (2 * h, 0.)]):
                mol.set_xy(n, x, y)
        sulfur = [u['anchor'] for u in _configured(mol) if mol.element_of(u['anchor']) == 16]
        assert len(sulfur) == 1, sulfur
        anchor = sulfur[0]
        double = [r for r in mol.neighbors_of(anchor) if mol.order_of(anchor, r) == 2]
        assert double, f'{smiles_ish} has no multiple bond at the sulfur to be tempted by'
        # the double bond is a terminal, acyclic, non-centre neighbour, so it is tempting
        assert mol.degree_of(double[0]) == 1 and not mol.bond_in_ring(anchor, double[0])
        assert (min(ids) == double[0]) is expect_multiple_first, 'the id order under test moved'

        wedges, _ = wedges_for_write(mol)
        chosen = [(n, w) for n, w, _ in wedges if n == anchor]
        assert chosen, f'{smiles_ish}: the sulfur lost its configuration entirely'
        for n, w in chosen:
            assert mol.order_of(n, w) == 1, \
                f'{smiles_ish}: wedge on a bond of order {mol.order_of(n, w)}'


def test_the_double_bond_guard_holds_when_the_S_O_is_the_LONGEST_bond():
    """Dimethyl sulfoxide with S=O drawn twice as long as either methyl, so the double bond wins the
    length term outright -- the term `_draw_cost` weighs most heavily, and the one the equal-length
    companion above cannot exercise.  Stable ids: 1=C, 2=S, 3=O, 4=C.
    """
    mol = read_smiles('C[S@@](=O)C')
    mol.set_xy(2, 0., 0.)       # S at the origin
    mol.set_xy(1, -1., .6)      # C(1), about 1.17 away
    mol.set_xy(4, 1., .6)       # C(4), the same
    mol.set_xy(3, 0., -2.)      # O(3), twice as far -- and it is the double-bond end
    assert mol.order_of(2, 3) == 2, 'premise: the longest bond at the sulfur is the double one'
    wedges, _ = wedges_for_write(mol)
    assert wedges, 'the sulfur lost its configuration entirely'
    for narrow, wide, _ in wedges:
        assert mol.order_of(narrow, wide) == 1, \
            f'wedge placed on bond {narrow}-{wide} which has order {mol.order_of(narrow, wide)}'


def test_a_bond_drawn_on_top_of_a_sibling_is_not_chosen_to_carry_the_wedge():
    """A wedge is attributed to a bond by lying along it, so a sibling drawn close to it leaves the
    reader unable to say which the triangle belongs to.  Built so the angle is the only difference: the
    near-collinear pair are the *longer* bonds with the *lower* stable ids, so every other term of the
    cost prefers the wrong bond.  Read from SMILES rather than assembled, because an atom built in code
    has an unknown implicit hydrogen count and a centre with an unknown direction is not yet a unit --
    `set_parity` stores the byte anyway, which is what makes that mistake quiet.
    """
    mol = read_smiles('[C@H](Cl)(Br)F')
    c, near, also, clean = sorted(mol)          # C, Cl, Br, F -- Cl and Br 8 degrees apart, both longer
    with mol.edit():
        for i, (x, y) in zip((c, near, also, clean),
                             [(0., 0.), (1.386, .195), (1.4, 0.), (-.5, -.87)]):
            mol.set_xy(i, x, y)
    assert _configured(mol), 'the fixture is not a perceived stereocentre'
    assert _near_collinear(mol, c, near) and _near_collinear(mol, c, also), 'fixture is not collinear'
    assert not _near_collinear(mol, c, clean)
    assert clean > near and clean > also, 'the id tie-break no longer prefers the wrong bond'

    wedges, _ = wedges_for_write(mol)
    assert len(wedges) == 1, wedges
    narrow, wide, _ = wedges[0]
    assert (narrow, wide) == (c, clean), \
        f'took {wide}, which is drawn on top of a sibling, over the one that is not'


def test_no_wedge_on_the_corpus_is_drawn_on_top_of_a_sibling(drawings):
    """The sweep, and why separation is a *tier* rather than only a cost term.  Morphine's centre 20 has
    three ring bonds, two within 11 degrees of a sibling and the third pointing at an atom centre 19's
    wedge owns; the adjacency preference lives in the pass and outranks anything the cost says, so only
    exhausting the well-separated bonds first takes this to zero.
    """
    offenders = []
    for title, record in sorted(drawings.items()):
        mol = _undrawn_molecule(record)
        wedges, _ = wedges_for_write(mol)
        for narrow, wide in _near_collinear_wedges(mol, wedges):
            offenders.append((title, narrow, wide))
    assert not offenders, f'{len(offenders)} unattributable wedges: {offenders}'


def test_the_collinearity_threshold_sits_in_a_gap_in_the_corpus_and_is_not_a_tuned_knob(drawings):
    """Why 30 degrees.  Nearest-sibling separation over the corpus is bimodal -- a mode at 110-120
    degrees, a few degenerate directions below 20, nothing between 20 and 30 -- so every cut in that band
    classifies this corpus identically.  A future corpus filling the gap fails here.
    """
    from math import atan2, degrees

    band = []
    for record in drawings.values():
        mol = _undrawn_molecule(record)
        for unit in _configured(mol):
            anchor = unit['anchor']
            neighbours = list(mol.neighbors_of(anchor))
            ax, ay = mol.xy_of(anchor)
            angles = {r: degrees(atan2(mol.xy_of(r)[1] - ay, mol.xy_of(r)[0] - ax))
                      for r in neighbours}
            for r in neighbours:
                seps = []
                for other in neighbours:
                    if other == r:
                        continue
                    d = abs(angles[r] - angles[other]) % 360.
                    seps.append(min(d, 360. - d))
                if seps and 20. <= min(seps) < 30.:
                    band.append((anchor, r, min(seps)))
    assert not band, f'the 20-30 degree band is no longer empty: {band}'


def test_every_shared_atom_that_remains_is_one_the_rules_could_not_avoid(drawings):
    """A centre whose only acyclic bond leads to an atom another wedge already claims must either share
    that atom or move onto a ring bond, and sharing is the lesser defect: a ring wedge is a false claim
    about the ring's plane, a shared atom only a crowded drawing.  So for each surviving pair at least
    one member has no untouched acyclic bond left -- lactose, whose two anomeric carbons have the one
    glycosidic oxygen between them, is the clean example.
    """
    unforced = []
    for title, record in sorted(drawings.items()):
        mol = _undrawn_molecule(record)
        wedges, _ = wedges_for_write(mol)
        refs = {u['anchor']: u['refs'] for u in _configured(mol)}
        for i, j in _adjacent_pairs(wedges):
            forced = False
            for k in (i, j):
                anchor = wedges[k][0]
                # atoms another wedge already claims: a wedge to one of them shares an atom as well
                claimed = {x for n, w in enumerate(wedges) if n != k for x in w[:2]}
                if not [r for r in refs[anchor] if r is not None
                        and not mol.bond_in_ring(anchor, r) and r not in claimed]:
                    forced = True
            if not forced:
                unforced.append((title, wedges[i], wedges[j]))
    assert not unforced, f'{len(unforced)} avoidable shared atoms: {unforced}'


def test_a_bridgehead_takes_a_ring_wedge_rather_than_losing_its_configuration(drawings):
    """A norbornane bridgehead's three heavy neighbours are all ring bonds and its fourth direction is an
    implicit hydrogen, which cannot be a wide end.  So the preference gives way, a ring bond is used and
    the fact is logged; refusing would write a molfile stating less stereochemistry than the molecule has.
    """
    mol = _undrawn_molecule(drawings['norbornanol'])
    bridgeheads = [u['anchor'] for u in _configured(mol)
                   if all(r is None or mol.bond_in_ring(u['anchor'], r) for r in u['refs'])]
    assert len(bridgeheads) == 2, bridgeheads

    wedges, log = wedges_for_write(mol)
    chosen = {a: b for a, b, _ in wedges}
    for anchor in bridgeheads:
        assert anchor in chosen, f'bridgehead {anchor} lost its configuration'
        assert mol.bond_in_ring(anchor, chosen[anchor])
    assert any('ring bond' in x for x in log), log
    # and the configuration really is expressed, not merely drawn
    assert _quality([drawings['norbornanol']])['round_trip_failures'] == 0


def test_the_choice_is_deterministic_across_repeated_calls(drawings):
    """No set iteration, no dict-order dependence, no float comparison that could go either way: the same
    layout must give the same picture every time, or a round trip becomes a diff."""
    for title in sorted(drawings)[:20]:
        mol = _undrawn_molecule(drawings[title])
        first, _ = wedges_for_write(mol)
        for _ in range(3):
            again, _ = wedges_for_write(mol)
            assert again == first, title
        # and independent of the arena a second read builds
        other = _undrawn_molecule(drawings[title])
        assert wedges_for_write(other)[0] == first, title


#: The measured quality of the write path over `drawings`, as a ratchet.  `ring_wedges` cannot go below
#: `_RING_WEDGE_FLOOR`: that many centres have no drawable bond that is not a ring bond, so the number is
#: a property of the molecules and not of the algorithm.  One of the 9 `adjacent_pairs` is bought
#: deliberately -- morphine's centre 20 can offer only a collinear bond or a shared atom, and ambiguity is
#: the worse failure -- which is what takes `near_collinear` to 0.
_EXPECTED_QUALITY = dict(molecules=92, centres=282, wedges=282, ring_wedges=83,
                         adjacent_pairs=9, near_collinear=0, unencoded=0, round_trip_failures=0)

#: Centres whose every heavy neighbour is a ring neighbour -- bridgeheads and ring-fusion carbons --
#: counted from the corpus rather than asserted, below.
_RING_WEDGE_FLOOR = 83


def test_the_corpus_quality_does_not_regress(drawings):
    """The before/after table as a test.  Every column is a ratchet in the direction that means
    better, and the two correctness columns are exact."""
    m = _quality(list(drawings.values()))
    assert m['round_trip_failures'] == 0, m
    assert m['unencoded'] == 0, m
    assert m['molecules'] == _EXPECTED_QUALITY['molecules'], m
    assert m['centres'] == _EXPECTED_QUALITY['centres'], m
    assert m['wedges'] == _EXPECTED_QUALITY['wedges'], m
    assert m['ring_wedges'] <= _EXPECTED_QUALITY['ring_wedges'], m
    assert m['adjacent_pairs'] <= _EXPECTED_QUALITY['adjacent_pairs'], m
    assert m['near_collinear'] <= _EXPECTED_QUALITY['near_collinear'], m


def test_the_ring_wedge_floor_is_a_property_of_the_molecules(drawings):
    """A centre whose every heavy neighbour is joined by a ring bond has nothing else to offer, its
    remaining direction being an implicit hydrogen with no atom to draw.  Counting those gives the floor
    the metric is measured against, so "83 ring wedges" does not read as a failure.
    """
    floor = 0
    for record in drawings.values():
        mol = _undrawn_molecule(record)
        for unit in _configured(mol):
            if all(r is None or mol.bond_in_ring(unit['anchor'], r) for r in unit['refs']):
                floor += 1
    assert floor == _RING_WEDGE_FLOOR, floor
    assert _quality(list(drawings.values()))['ring_wedges'] == floor, \
        'every ring wedge should be a forced one'


def test_every_configured_centre_on_the_corpus_survives_a_written_round_trip(drawings, tmp_path):
    """The gate, through the real writer and reader.  `_quality` asks the same parity function the chooser
    used, so on its own it cannot catch an emitter that writes the wedge down wrongly -- a swapped narrow
    and wide end, say; this goes out through `emit_v2000` and back through `parse_v2000`.
    """
    for title, record in sorted(drawings.items()):
        mol = _undrawn_molecule(record)
        configured = _configured(mol)
        if not configured:
            continue
        target = {u['anchor']: mol.parity_of(u['anchor']) for u in configured}

        lines, _ = emit_v2000(mol)
        back, _, _ = _build(lines)
        # the writer emits atoms in `atom_numbers` order, so position i of one is position i of the other
        forward = dict(zip(mol.atom_numbers, back.atom_numbers))
        got = {a: back.parity_of(forward[a]) for a in target}
        assert got == {a: p for a, p in target.items()}, \
            f'{title}: wrote {target}, read back {got}'


# --- the fixtures

_TETRA_WITH_H = """fluorochlorobromomethane with the hydrogen drawn
  test
comment
  5  4  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.9500    0.3100    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.5900   -0.8100    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -0.9500    0.3100    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
M  END""".split('\n')

_TETRA_FLAT = """bromochlorofluoromethane
  test
comment
  4  3  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.8700   -0.5000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -0.8700   -0.5000    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
M  END""".split('\n')


#: ``_TETRA_FLAT``'s molecule in V3000, every atom at the origin and the configuration stated only in the
#: atom line's ``CFG=``.  Written out rather than derived from ``_TETRA_FLAT``: the two versions state the
#: field in unrelated syntaxes, three fixed columns against a keyword in a free-format line.
_TETRA_FLAT_V3000_CFG = """bromochlorofluoromethane
  test
comment
  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0 0.0 0.0 0 CFG=1
M  V30 2 F 0.0 0.0 0.0 0
M  V30 3 Cl 0.0 0.0 0.0 0
M  V30 4 Br 0.0 0.0 0.0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 1 3
M  V30 3 1 1 4
M  V30 END BOND
M  V30 END CTAB
M  END""".split('\n')


_TETRA_DEGENERATE = """three neighbours on one line
  test
comment
  5  4  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0000   -1.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.0000    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
    1.0000   -1.0000    0.0000 I   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  1  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
M  END""".split('\n')


#: 1-bromopropan-2-amine, drawn with its one wedge running from the amine nitrogen INTO the stereocentre.
#: The point is at atom 1, which has nothing to configure; the centre is atom 2.  Line 9 is the wedged
#: bond and line 11 is the centre's bond to the CH2Br arm, both edited by the tests above.
_TETRA_WEDGE_AT_WIDE_END = """1-bromopropan-2-amine
  test
comment
  5  4  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.8700   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8700   -1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7400    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6100   -0.5000    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  6  0  0  0
  2  3  1  0  0  0  0
  2  4  1  0  0  0  0
  4  5  1  0  0  0  0
M  END""".split('\n')


#: 2,3-dibromobutane: one wedge on the bond between its two stereocentres, so both ends of it can hold a
#: configuration and the convention -- the narrow end -- is the only thing that says which one does.
_DIBROMOBUTANE = """2,3-dibromobutane
  test
comment
  6  5  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8700    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7400    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6100    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8700    1.5000    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
    1.7400   -1.0000    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  1  0  0  0
  3  4  1  0  0  0  0
  2  5  1  0  0  0  0
  3  6  1  0  0  0  0
M  END""".split('\n')


_BUTENE_CIS = """but-2-ene
  test
comment
  4  3  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.8000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5000    0.8000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  1  4  1  0  0  0  0
M  END""".split('\n')


# The same molecule with atom 3 reflected across the double bond, so the two methyls are on opposite
# sides.
_BUTENE_TRANS = _BUTENE_CIS[:6] + ['    1.5000   -0.8000    0.0000 C   0  0  0  0  0  0  0  0  0  0'
                                   '  0  0'] + _BUTENE_CIS[7:]


# Penta-2,3-diene, the textbook axially chiral allene: the C=C=C axis along x and the wedge on the bond
# from the near *terminal* to its methyl, the anchor being the central carbon with no single bond to wedge.
_ALLENE_WEDGE = """penta-2,3-diene
  test
comment
  5  4  0  0  0  0            999 V2000
   -0.8700    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8700    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  1  0  0  0
  2  3  2  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
M  END""".split('\n')


# The other enantiomer: both methyls reflected across the axis, wedge unchanged.  Written out in full
# rather than sliced from the fixture above, since the point is that the coordinates decide.
_ALLENE_MIRROR = """(mirror image of penta-2,3-diene)
  test
comment
  5  4  0  0  0  0            999 V2000
   -0.8700   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8700   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  1  0  0  0
  2  3  2  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
M  END""".split('\n')


# The same drawing with no wedge: a flat allene, which is what most real files contain.
_ALLENE_FLAT = _ALLENE_WEDGE[:9] + ['  2  1  1  0  0  0  0'] + _ALLENE_WEDGE[10:]


# 1,3-dibromo-1,3-difluoroallene -- the compound the core's axial sign convention is documented against
# in `core/_inchi.pxi`.  Every direction is a heavy atom, so the frame has no holes, and both bonds of the
# near terminal are wedged up and down, which is how drawing packages state "perpendicular to the paper".
_ALLENE_TETRA = """1,3-dibromo-1,3-difluoropropa-1,2-diene
  test
comment
  7  6  0  0  0  0            999 V2000
   -0.8700    0.5000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8700   -0.5000    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8700    0.5000    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
    2.8700   -0.5000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  1  0  0  0
  2  3  1  6  0  0  0
  2  4  2  0  0  0  0
  4  5  2  0  0  0  0
  5  6  1  0  0  0  0
  5  7  1  0  0  0  0
M  END""".split('\n')


# Hepta-2,3,4,5-tetraene: five cumulated carbons, still axially chiral, and the reason this module
# cannot find a terminal by looking at the anchor's neighbours.
_CUMULENE5 = """hepta-2,3,4,5-tetraene
  test
comment
  7  6  0  0  0  0            999 V2000
   -0.8700    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.8700    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  1  0  0  0
  2  3  2  0  0  0  0
  3  4  2  0  0  0  0
  4  5  2  0  0  0  0
  5  6  2  0  0  0  0
  6  7  1  0  0  0  0
M  END""".split('\n')


# The far terminal's methyl drawn straight along the axis: the half of the frame that is supposed to
# lie in the paper has no width, and no wedge on the other terminal can supply it.
_ALLENE_DEGENERATE = _ALLENE_WEDGE[:8] + \
    ['    3.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0'] + _ALLENE_WEDGE[9:]


# 4-methylpenta-2,3-diene: the near terminal carries two methyls, so it has one direction twice and
# the core emits no axial unit at all.
_ALLENE_SYMMETRIC = """4-methylpenta-2,3-diene
  test
comment
  6  5  0  0  0  0            999 V2000
   -0.8700    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8700   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8700    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  1  0  0  0
  2  3  1  0  0  0  0
  2  4  2  0  0  0  0
  4  5  2  0  0  0  0
  5  6  1  0  0  0  0
M  END""".split('\n')


# 3-ethylhexa-3,4-diene: two ethyls on one terminal.  Both directions are named, so the unit IS
# emitted, and it is the automorphism search that unmarks it -- the other half of "not stereogenic".
_ALLENE_NOT_STEREOGENIC = """3-ethylhexa-3,4-diene
  test
comment
  8  7  0  0  0  0            999 V2000
   -1.7400    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8700    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8700   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7400   -1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8700    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  3  2  1  1  0  0  0
  2  1  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  1  0  0  0  0
  3  6  2  0  0  0  0
  6  7  2  0  0  0  0
  7  8  1  0  0  0  0
M  END""".split('\n')


# 2-chloro-2'-fluorobiphenyl, the smallest biaryl whose two ortho pairs are distinguishable: two regular
# hexagons of unit edge along x, pivots at (1, 0) and (2, 0), both ortho substituents drawn UP.  The wedge
# is on a RING bond, pivot 1 to its ortho 2, which is where a drawing states an axial configuration; the
# Kekule form puts a single bond there, so the code sits on a bond MDL defines it for.
_ATROPO_WEDGE = """2-chloro-2'-fluorobiphenyl
  test
comment
 14 15  0  0  0  0            999 V2000
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5000    0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5000    0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5000   -0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5000   -0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    1.7320    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5000    0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5000    0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5000   -0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5000   -0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    1.7320    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  1  0  0  0
  2  3  2  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
  5  6  1  0  0  0  0
  6  1  2  0  0  0  0
  2  7  1  0  0  0  0
  1  8  1  0  0  0  0
  8  9  1  0  0  0  0
  9 10  2  0  0  0  0
 10 11  1  0  0  0  0
 11 12  2  0  0  0  0
 12 13  1  0  0  0  0
 13  8  2  0  0  0  0
  9 14  1  0  0  0  0
M  END""".split('\n')


# The other enantiomer, drawn as one: the whole layout reflected in x, the wedge unchanged.  Written out
# in full rather than sliced from the fixture above, since the point is that the coordinates decide.
_ATROPO_MIRROR = """(mirror image of 2-chloro-2'-fluorobiphenyl)
  test
comment
 14 15  0  0  0  0            999 V2000
   -1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5000    0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5000    0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5000   -0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5000   -0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0000    1.7320    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -2.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5000    0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5000    0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5000   -0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5000   -0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0000    1.7320    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  1  0  0  0
  2  3  2  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
  5  6  1  0  0  0  0
  6  1  2  0  0  0  0
  2  7  1  0  0  0  0
  1  8  1  0  0  0  0
  8  9  1  0  0  0  0
  9 10  2  0  0  0  0
 10 11  1  0  0  0  0
 11 12  2  0  0  0  0
 12 13  1  0  0  0  0
 13  8  2  0  0  0  0
  9 14  1  0  0  0  0
M  END""".split('\n')


# The same molecule with ring B drawn the other way up: atoms 9 to 14 reflected in y, so the far pivot's
# first ring direction moves to the other side of the axis.  Same constitution, same wedge, opposite
# configuration -- which is what `test_which_side_of_the_axis_the_far_pair_is_drawn_on_decides_the_sign`
# is about.
_ATROPO_RING_B_FLIPPED = _ATROPO_WEDGE[:12] + [
    '    2.5000   -0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
    '    3.5000   -0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
    '    4.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
    '    3.5000    0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
    '    2.5000    0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
    '    2.0000   -1.7320    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0'] + _ATROPO_WEDGE[18:]


# The same drawing with no wedge anywhere: a flat biaryl, which is what all but 88 of the 4368 axes in a
# 119,534-molfile sample of a production corpus are.
_ATROPO_FLAT = _ATROPO_WEDGE[:18] + ['  1  2  1  0  0  0  0'] + _ATROPO_WEDGE[19:]


def test_the_allene_fixtures_are_what_the_tests_above_assume():
    """An allene fixture encodes a hand-computed determinant, so a moved coordinate does not fail: it
    silently changes the right answer."""
    mol, _, _ = _build(_ALLENE_WEDGE)
    assert [mol.element_of(s) for s in mol.atom_numbers] == [6] * 5
    assert [(x, y) for x, y in (mol.xy_of(s) for s in mol.atom_numbers)] == \
        [(-0.87, 0.5), (0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (2.87, 0.5)]
    assert list(mol.wedges()) == [(2, 1, 1)], 'the wedge must be UP, narrow end on the terminal'

    # the mirror differs from the fixture above in the two substituent y values and nothing else
    assert [i for i, (a, b) in enumerate(zip(_ALLENE_WEDGE, _ALLENE_MIRROR)) if a != b] == [0, 4, 8]
    mirror, _, _ = _build(_ALLENE_MIRROR)
    assert list(mirror.wedges()) == [(2, 1, 1)]

    flat, _, _ = _build(_ALLENE_FLAT)
    assert not any(flat.wedges())
    assert [i for i, (a, b) in enumerate(zip(_ALLENE_WEDGE, _ALLENE_FLAT)) if a != b] == [9]

    degenerate, _, _ = _build(_ALLENE_DEGENERATE)
    assert degenerate.xy_of(5) == (3.0, 0.0), 'the far methyl must be ON the axis'
    assert [i for i, (a, b) in enumerate(zip(_ALLENE_WEDGE, _ALLENE_DEGENERATE)) if a != b] == [8]

    tetra, _, _ = _build(_ALLENE_TETRA)
    assert [tetra.element_of(s) for s in tetra.atom_numbers] == [9, 6, 35, 6, 6, 35, 9]
    assert sorted(tetra.wedges()) == [(2, 1, 1), (2, 3, 2)]

    cumulene, _, _ = _build(_CUMULENE5)
    assert sorted(b.order for b in cumulene.bonds()) == [1, 1, 2, 2, 2, 2]


def test_the_atropisomer_fixtures_are_what_the_tests_above_assume():
    """The four biaryls are one constitution in four layouts, and the tests above rewrite bond lines by
    index: a shifted line would wedge a different bond and still read a parity."""
    for lines in (_ATROPO_WEDGE, _ATROPO_MIRROR, _ATROPO_RING_B_FLIPPED, _ATROPO_FLAT):
        mol, _, _ = _build(lines)
        assert format(mol) == 'C1=C(C(C2=CC=CC=C2Cl)=CC=C1)F'
        assert len(_atropo(mol)) == 1, 'the fixture must keep a stereogenic axis'
    assert _ATROPO_WEDGE[18] == '  1  2  1  1  0  0  0'      # pivot 1 to its ortho, the drawn wedge
    assert _ATROPO_WEDGE[23] == '  6  1  2  0  0  0  0'      # pivot 1's other ring bond, and a double one
    assert _ATROPO_WEDGE[26] == '  8  9  1  0  0  0  0'      # pivot 8 to its ortho
    assert _ATROPO_WEDGE[33] == 'M  END' and len(_ATROPO_WEDGE) == 34

    # the mirror is a reflection and nothing else: same bond block, every x negated
    assert _ATROPO_WEDGE[18:] == _ATROPO_MIRROR[18:]
    mol, _, _ = _build(_ATROPO_WEDGE)
    mirror, _, _ = _build(_ATROPO_MIRROR)
    assert all(mirror.xy_of(s) == (-x, y) for s in mol.atom_numbers for x, y in [mol.xy_of(s)])

    # ring B alone is reflected, in y, so the axis and ring A are where they were
    flipped, _, _ = _build(_ATROPO_RING_B_FLIPPED)
    assert [i for i, (a, b) in enumerate(zip(_ATROPO_WEDGE, _ATROPO_RING_B_FLIPPED)) if a != b] \
        == [12, 13, 15, 16, 17], 'index 14 is atom 11, para to the axis and at y = 0, so a reflection ' \
                                 'leaves its line as it was'
    assert all(flipped.xy_of(s) == mol.xy_of(s) for s in (1, 2, 3, 4, 5, 6, 7, 8))
    assert all(flipped.xy_of(s) == (x, -y) for s in (9, 10, 11, 12, 13, 14) for x, y in [mol.xy_of(s)])

    flat, _, _ = _build(_ATROPO_FLAT)
    assert not any(flat.wedges())
    assert [i for i, (a, b) in enumerate(zip(_ATROPO_WEDGE, _ATROPO_FLAT)) if a != b] == [18]


def test_the_fixtures_are_what_the_tests_above_assume():
    """A wrong column in a fixture makes several tests above pass for the wrong reason: the atom line is 31
    characters before its flag columns, and an off-by-one there moves a wedge to a different bond."""
    mol, _, _ = _build(_TETRA_FLAT)
    assert [mol.element_of(s) for s in mol.atom_numbers] == [6, 9, 17, 35]
    assert _TETRA_FLAT[8].startswith('  1  2  1')
    # the V3000 sibling must be the same molecule, or the cross-version equality above compares two
    # different questions.  Read through the version sniffer, so a wrong stamp fails here.
    three, _, _ = _build(_TETRA_FLAT_V3000_CFG)
    assert sniff_version(_TETRA_FLAT_V3000_CFG, []) == V3000_STAMP
    assert [three.element_of(s) for s in three.atom_numbers] == [6, 9, 17, 35]
    assert sorted((min(b.n, b.m), max(b.n, b.m), b.order) for b in three.bonds()) \
        == sorted((min(b.n, b.m), max(b.n, b.m), b.order) for b in mol.bonds())
    for lines in (_BUTENE_CIS, _BUTENE_TRANS):
        mol, _, _ = _build(lines)
        assert sorted(b.order for b in mol.bonds()) == [1, 1, 2]
        assert len(lines) == 12 and lines[-1] == 'M  END'
    # the reflection is the only difference, and it is on the atom the docstring says it is
    assert [i for i, (a, b) in enumerate(zip(_BUTENE_CIS, _BUTENE_TRANS)) if a != b] == [6]


def test_wedge_in_file_order_puts_the_narrow_end_first():
    """CTfile writes the wedge's point at the *first* atom, so the bond may be written reversed.  Shared
    by both writers, and it cannot read the arena instead: `wedges_for_write` may have chosen these
    wedges, so they are not in the arena to read.
    """
    from chython.core.wedge import wedge_in_file_order

    wedge_of = {(7, 3): 1}
    assert wedge_in_file_order(wedge_of, 7, 3) == (7, 3, 1)
    assert wedge_in_file_order(wedge_of, 3, 7) == (7, 3, 1), 'the bond is written narrow end first'
    assert wedge_in_file_order(wedge_of, 4, 5) == (4, 5, None), 'no wedge: the pair is left alone'
