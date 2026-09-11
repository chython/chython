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
"""Atom labels: what gets written, and how big the ink is.

Two jobs in one module because they cannot disagree: the label decides what is drawn, the box says where
bonds must stop.  Deriving the box from the font size instead of the composed text trims a wide `NH2+` as
if it were as narrow as `N`, and the bond runs through the glyphs.
"""
from math import hypot
from pytest import approx, raises
from chython import smiles
from chython.core import H_UNKNOWN
from chython.depict._config import cpk
from chython.depict.label import CPK, _DEGENERATE_SPAN, element_colour, is_labelled, labels
from chython.depict.scene import to_hex
from chython.depict.style import DepictStyle


def _laid_out(text):
    mol = smiles(text)
    mol.clean2d()
    return mol, mol.coordinates()


def test_is_labelled_table():
    """direct table: carbon flags and non-carbon are covered case-by-case"""
    mol = smiles('CCO')
    mol.clean2d()
    carbon = next(a for a in mol.atoms() if a.atomic_symbol == 'C' and a.degree > 0)
    oxygen = next(a for a in mol.atoms() if a.atomic_symbol == 'O')
    default = DepictStyle()
    with_carbon = DepictStyle().tuned(**{'atom.carbon': True})
    without_radicals = DepictStyle().tuned(**{'atom.radicals': False})

    assert is_labelled(oxygen, default) is True,    'non-carbon is always written'
    assert is_labelled(carbon, default) is False,   'skeletal carbon is a bare vertex'
    assert is_labelled(carbon, with_carbon) is True, 'carbon=True forces a symbol'

    # A radical carbon is written only when the style asks for radicals
    mol2 = smiles('CC |^1:1|')
    radical_c = next(a for a in mol2.atoms() if a.is_radical)
    assert is_labelled(radical_c, default) is True
    assert is_labelled(radical_c, without_radicals) is False, \
        'radicals=False should hide the radical from the labelling decision'


def test_every_atom_gets_an_entry_even_when_it_is_not_labelled():
    """the caller indexes this by atom and must never have to check for a missing key

    `bonds.py` looks up both ends of every bond; a dict with holes puts a `.get(n)` and a `None` branch at
    every trimming site, which is the shape that loses a trim.
    """
    mol, plane = _laid_out('CCO')
    result = labels(mol, plane, DepictStyle())
    assert set(result) == set(mol)
    assert sum(label.text is not None for label in result.values()) == 1, 'the O, and neither carbon'


def test_a_skeletal_drawing_labels_no_plain_carbon():
    mol, plane = _laid_out('CCCC')
    assert all(label.text is None for label in labels(mol, plane, DepictStyle()).values())


def test_carbon_is_labelled_when_the_style_says_so():
    mol, plane = _laid_out('CCCC')
    style = DepictStyle().tuned(**{'atom.carbon': True})
    result = labels(mol, plane, style)
    assert all(label.text is not None for label in result.values())
    assert result[next(iter(mol))].text.runs[0].text == 'C'


def test_a_charged_carbon_is_always_labelled():
    """a skeletal drawing hides plain carbons; a carbanion is not plain, and hiding it hides the charge"""
    mol, plane = _laid_out('C[CH2-]')
    labelled = [sid for sid, label in labels(mol, plane, DepictStyle()).items() if label.text]
    assert len(labelled) == 1
    assert mol.atom(labelled[0]).charge == -1


def test_a_radical_carbon_is_always_labelled():
    mol = smiles('CC |^1:1|')
    mol.clean2d()
    labelled = [sid for sid, label in labels(mol, mol.coordinates(), DepictStyle()).items() if label.text]
    assert len(labelled) == 1


def test_an_isotopic_carbon_is_always_labelled():
    mol, plane = _laid_out('[13CH4]')
    assert labels(mol, plane, DepictStyle())[next(iter(mol))].text is not None


def test_a_lone_atom_is_labelled_whatever_it_is():
    """methane drawn skeletally is an empty picture, which is not a drawing of methane"""
    mol, plane = _laid_out('C')
    assert labels(mol, plane, DepictStyle())[next(iter(mol))].text is not None


def test_hydrogen_count_is_written_as_a_subscript_run():
    # smilesdrawer places the heteroatom to the left of carbon in C-X molecules, so the carbon
    # neighbour is to the right and the hydrogen flips left: the label reads H first.
    mol, plane = _laid_out('CO')
    label = [lb for lb in labels(mol, plane, DepictStyle()).values() if lb.text][0]
    texts = [run.text for run in label.text.runs]
    assert set(texts) == {'O', 'H'}, f'expected O and H runs, got {texts}'
    assert len(label.text.runs) == 2, 'OH has one hydrogen and writes no count'

    mol, plane = _laid_out('CN')
    label = [lb for lb in labels(mol, plane, DepictStyle()).values() if lb.text][0]
    texts = [run.text for run in label.text.runs]
    assert set(texts) == {'N', 'H', '2'}, f'expected N, H, 2 runs, got {texts}'
    subscript = [run for run in label.text.runs if run.text == '2'][0]
    assert subscript.dy < 0., 'the count is a subscript, and the scene is y-up'


def test_one_hydrogen_writes_no_count():
    """H1 is not a thing chemists write"""
    mol, plane = _laid_out('CO')
    label = [lb for lb in labels(mol, plane, DepictStyle()).values() if lb.text][0]
    assert '1' not in ''.join(run.text for run in label.text.runs)


def test_an_unknown_hydrogen_count_writes_nothing_rather_than_zero():
    """H_UNKNOWN means "not derivable", and rendering it as H0 asserts a fact the input refused"""
    mol = smiles('CO')
    sid = [a.n for a in mol.atoms() if a.atomic_symbol == 'O'][0]
    with mol.edit():
        mol.set_hydrogens(sid, H_UNKNOWN)
    mol.clean2d()
    assert mol.atom(sid).implicit_h is None, 'expected an unknown count to set up the test'
    label = labels(mol, mol.coordinates(), DepictStyle())[sid]
    assert ''.join(run.text for run in label.text.runs) == 'O'


def test_an_unknown_hydrogen_count_can_be_marked_when_asked():
    mol = smiles('CO')
    sid = [a.n for a in mol.atoms() if a.atomic_symbol == 'O'][0]
    with mol.edit():
        mol.set_hydrogens(sid, H_UNKNOWN)
    mol.clean2d()
    style = DepictStyle().tuned(**{'atom.unknown_h_marks': True})
    label = labels(mol, mol.coordinates(), style)[sid]
    assert '?' in ''.join(run.text for run in label.text.runs)


def test_a_charge_is_a_superscript_and_uses_the_typographic_minus():
    """a hyphen is not a minus sign, and at 8 pt in a printed figure the difference is visible"""
    mol, plane = _laid_out('C[O-]')
    label = [lb for lb in labels(mol, plane, DepictStyle()).values() if lb.text][0]
    charge = [run for run in label.text.runs if run.dy > 0.]
    assert charge and charge[0].text == '−'
    assert '-' not in ''.join(run.text for run in label.text.runs)


def test_a_double_charge_writes_the_magnitude_before_the_sign():
    mol, plane = _laid_out('[Ca+2]')
    label = labels(mol, plane, DepictStyle())[next(iter(mol))]
    assert ''.join(run.text for run in label.text.runs).endswith('2+')


def test_an_isotope_is_a_leading_superscript():
    mol, plane = _laid_out('[13CH4]')
    label = labels(mol, plane, DepictStyle())[next(iter(mol))]
    assert label.text.runs[0].text == '13'
    assert label.text.runs[0].dy > 0.
    assert label.text.runs[1].text == 'C'


def test_a_radical_is_not_a_text_run():
    """a dot glyph's size and position depend on the font; a drawn disc does not

    The dot is geometry, returned by `figure.py` from the label's anchor and the style's radius.
    """
    mol = smiles('CC |^1:1|')
    mol.clean2d()
    label = [lb for lb in labels(mol, mol.coordinates(), DepictStyle()).values() if lb.text][0]
    all_text = ''.join(run.text for run in label.text.runs)
    assert '•' not in all_text
    assert '·' not in all_text
    assert 'C' in all_text, 'the radical carbon must still have its symbol in the runs'


def test_the_box_is_measured_and_grows_with_the_label():
    """the defect this module exists to fix: a wide label needs a wide clearance"""
    mol, plane = _laid_out('CO')
    narrow = [lb for lb in labels(mol, plane, DepictStyle()).values() if lb.text][0]
    mol2, plane2 = _laid_out('CS(=O)(=O)N')
    wide = max((lb for lb in labels(mol2, plane2, DepictStyle()).values() if lb.text),
               key=lambda lb: lb.box.width)
    assert wide.box.width > narrow.box.width


def test_the_box_is_padded_by_exactly_the_style_pad():
    mol, plane = _laid_out('CO')
    tight = [lb for lb in labels(mol, plane, DepictStyle().tuned(**{'label.pad': 0.})).values() if lb.text][0]
    padded = [lb for lb in labels(mol, plane, DepictStyle().tuned(**{'label.pad': .1})).values() if lb.text][0]
    assert padded.box.width == approx(tight.box.width + .2)
    assert padded.box.height == approx(tight.box.height + .2)


def test_the_box_straddles_the_atom_point():
    """a label is centred ON the atom, so a bond arriving from any direction is trimmed symmetrically"""
    mol, plane = _laid_out('CO')
    sid = [a.n for a in mol.atoms() if a.atomic_symbol == 'O'][0]
    style = DepictStyle()
    label = labels(mol, plane, style)[sid]
    x, y = plane[sid]
    assert label.box.min_x < x < label.box.max_x
    assert label.box.min_y < y < label.box.max_y
    # Vertically the label is centred near the atom point (within one-tenth of a font height)
    mid_y = (label.box.min_y + label.box.max_y) / 2
    assert abs(mid_y - y) < style.label.size / 10


def test_an_unlabelled_atom_has_a_degenerate_box_at_its_point():
    """so `bonds.py` trims against it with the same code path and no `if label.text is None`"""
    mol, plane = _laid_out('CCCC')
    sid = next(iter(mol))
    label = labels(mol, plane, DepictStyle())[sid]
    assert tuple(label.box) == approx((plane[sid][0], plane[sid][1]) * 2)


def test_box_includes_hydrogen_and_charge_runs():
    """the bond trimming box covers the full composite, not just the element symbol"""
    mol, plane = _laid_out('CO')
    o_sid = [a.n for a in mol.atoms() if a.atomic_symbol == 'O'][0]
    with_h = labels(mol, plane, DepictStyle())[o_sid]
    no_h = labels(mol, plane, DepictStyle().tuned(**{'atom.hydrogens': False}))[o_sid]
    assert with_h.box.width > no_h.box.width, 'H run should widen the box'

    mol2, plane2 = _laid_out('C[O-]')
    o2_sid = [a.n for a in mol2.atoms() if a.atomic_symbol == 'O'][0]
    charged = labels(mol2, plane2, DepictStyle())[o2_sid]
    plain = labels(mol2, plane2, DepictStyle().tuned(**{'atom.charges': False}))[o2_sid]
    assert charged.box.width > plain.box.width, 'superscript charge run should widen the box'


def test_the_hydrogen_count_goes_on_the_side_with_more_room():
    """NH2 on the left of a chain must read H2N, or the H sits on top of the bond"""
    mol, plane = _laid_out('NCCCCCC')
    sid = [a.n for a in mol.atoms() if a.atomic_symbol == 'N'][0]
    neighbour = mol.atom(next(iter(mol.neighbors_of(sid))))
    label = labels(mol, plane, DepictStyle())[sid]
    written = ''.join(run.text for run in label.text.runs)
    if neighbour.x > mol.atom(sid).x:
        assert written.startswith('H'), f'the neighbour is to the right, so the H goes left: {written}'
    else:
        assert written.startswith('N'), written


def test_flip_direction_is_determined_by_neighbour_position():
    """hand-built planes: no clean2d, no layout dependency"""
    mol = smiles('NC')
    n_sid = [a.n for a in mol.atoms() if a.atomic_symbol == 'N'][0]
    c_sid = [a.n for a in mol.atoms() if a.atomic_symbol == 'C'][0]

    # C strictly to the right → every neighbour right → flip → H precedes N
    plane_right = {n_sid: (0.0, 0.0), c_sid: (1.0, 0.0)}
    label = labels(mol, plane_right, DepictStyle())[n_sid]
    runs = [r.text for r in label.text.runs]
    assert runs[0] == 'H', f'with C to the right, H should lead; got {runs}'

    # C strictly to the left → no flip → N precedes H
    plane_left = {n_sid: (0.0, 0.0), c_sid: (-1.0, 0.0)}
    label = labels(mol, plane_left, DepictStyle())[n_sid]
    runs = [r.text for r in label.text.runs]
    assert runs[0] == 'N', f'with C to the left, N should lead; got {runs}'


def test_isotope_stays_beside_the_symbol_when_the_label_is_flipped():
    """when H groups move left, 13 still reads immediately before C, not 13 H2 C"""
    mol = smiles('[13CH3]CC')
    c13_sid = [a.n for a in mol.atoms() if a.isotope][0]
    other_sids = [a.n for a in mol.atoms() if not a.isotope]
    # All non-C13 atoms strictly to the right → flip=True
    plane = {c13_sid: (0.0, 0.0)}
    for i, sid in enumerate(other_sids, 1):
        plane[sid] = (float(i), 0.0)
    label = labels(mol, plane, DepictStyle())[c13_sid]
    runs = [r.text for r in label.text.runs]
    isotope_idx = next(i for i, r in enumerate(label.text.runs) if r.text == '13')
    symbol_idx = next(i for i, r in enumerate(label.text.runs) if r.text == 'C')
    assert symbol_idx == isotope_idx + 1, (
        f'isotope must be immediately left of symbol; got order {runs}'
    )


def test_the_label_text_has_correct_fill_and_anchor():
    """fill comes from element_colour, Text.anchor mode is middle, Label.anchor is the atom point"""
    mol, plane = _laid_out('CO')
    style = DepictStyle()
    o_sid = [a.n for a in mol.atoms() if a.atomic_symbol == 'O'][0]
    label = labels(mol, plane, style)[o_sid]
    assert label.text.anchor == 'middle', f'expected middle, got {label.text.anchor!r}'
    assert label.text.fill == to_hex(cpk[7]), 'oxygen fill must come from the CPK table'
    # Label.anchor is the atom point — figure.py centres radical dots, halos and ribbons on it
    assert label.anchor == approx(plane[o_sid])


def test_element_colour_comes_from_the_table_and_carbon_from_the_style():
    mol, plane = _laid_out('CO')
    style = DepictStyle()
    oxygen = [a for a in mol.atoms() if a.atomic_symbol == 'O'][0]
    carbon = [a for a in mol.atoms() if a.atomic_symbol == 'C'][0]
    assert element_colour(oxygen, style) == to_hex(CPK[7])
    assert element_colour(carbon, style) == style.atom.carbon_colour


def test_colouring_can_be_turned_off_entirely():
    mol, plane = _laid_out('CO')
    style = DepictStyle().tuned(**{'atom.colour_by_element': False,
                                   'atom.default_colour': '#112233'})
    oxygen = [a for a in mol.atoms() if a.atomic_symbol == 'O'][0]
    assert element_colour(oxygen, style) == '#112233'


def test_the_cpk_table_is_the_one_in_config_not_a_copy():
    assert CPK is cpk
    assert len(CPK) == 118
    assert all(c.startswith('#') and len(c) == 7 for c in CPK)


def test_a_map_number_is_written_by_default_and_can_be_turned_off():
    """on by default, because a mapped structure whose mapping is invisible is a misleading picture

    The off branch TUNES the flag off rather than reading a bare `DepictStyle()`: the assertion has to be
    about the flag, not about whichever way the default happens to point.
    """
    mol = smiles('[CH3:3]CCC')
    mol.clean2d()
    sid = [a.n for a in mol.atoms() if a.map_number == 3][0]

    default = labels(mol, mol.coordinates(), DepictStyle())[sid]
    assert any(run.text == '3' for ann in default.annotations for run in ann.runs), \
        'the default style must write the map numbers the molecule carries'
    assert not any(':' in run.text for ann in default.annotations for run in ann.runs), \
        'the number is drawn bare: the colon is SMILES punctuation, not part of a mapping'

    off = labels(mol, mol.coordinates(), DepictStyle().tuned(**{'atom.map_numbers': False}))[sid]
    assert not any(run.text == '3' for ann in off.annotations for run in ann.runs), \
        'atom.map_numbers=False must withhold it'


def test_an_unmapped_atom_gets_no_map_annotation_under_the_default():
    """"numbers only where mapped": a `map_number` of 0 is nothing to write, not a `0`"""
    mol = smiles('[CH3:3]CCC')
    mol.clean2d()
    bare = [a.n for a in mol.atoms() if not a.map_number]
    assert len(bare) == 3, 'premise: three unmapped carbons'
    result = labels(mol, mol.coordinates(), DepictStyle())
    for sid in bare:
        assert result[sid].annotations == (), f'atom {sid} is unmapped and must carry no annotation'


def test_a_stored_cip_descriptor_is_written_only_when_asked():
    mol = smiles('C[C@H](N)O')
    mol.clean2d()
    sid = [a.n for a in mol.atoms() if a.parity][0]
    with mol.edit():
        mol.set_atom_cip(sid, 'R')

    off = labels(mol, mol.coordinates(), DepictStyle())[sid]
    assert not any('R' in run.text
                   for ann in off.annotations for run in ann.runs), \
        'CIP annotation must not appear without atom.stereo_labels=True'

    on = labels(mol, mol.coordinates(),
                DepictStyle().tuned(**{'atom.stereo_labels': True}))[sid]
    assert any('R' in run.text
               for ann in on.annotations for run in ann.runs), \
        'CIP annotation must appear when atom.stereo_labels=True'


def test_cip_annotation_is_not_in_the_bond_trim_box():
    """annotations are excluded from Label.box: a descriptor beside a bond is acceptable"""
    mol = smiles('C[C@H](N)O')
    mol.clean2d()
    sid = [a.n for a in mol.atoms() if a.parity][0]
    with mol.edit():
        mol.set_atom_cip(sid, 'R')
    style = DepictStyle().tuned(**{'atom.stereo_labels': True})
    lbl_with = labels(mol, mol.coordinates(), style)[sid]
    lbl_without = labels(mol, mol.coordinates(), DepictStyle())[sid]
    assert lbl_with.box.width == approx(lbl_without.box.width)
    assert lbl_with.box.height == approx(lbl_without.box.height)


def test_annotation_fill_placement_and_size():
    """CIP and map annotations: fill, side and size, on a hand-built plane so the flip is deterministic"""
    mol = smiles('[NH:3]CC')
    n_sid = [a.n for a in mol.atoms() if a.atomic_symbol == 'N'][0]
    c1_sid = [a.n for a in mol.atoms() if a.atomic_symbol == 'C'][0]
    c2_sid = [a.n for a in mol.atoms() if a.atomic_symbol == 'C'][-1]
    with mol.edit():
        mol.set_atom_cip(n_sid, 'R')

    # N at origin, both carbons to the right → flip=True → H goes left, annotation follows it left
    plane = {n_sid: (0.0, 0.0), c1_sid: (1.0, 0.0), c2_sid: (2.0, 0.0)}
    style = DepictStyle().tuned(**{'atom.stereo_labels': True, 'atom.map_numbers': True})
    lbl = labels(mol, plane, style)[n_sid]

    assert len(lbl.annotations) == 2, f'expected CIP and map annotations, got {len(lbl.annotations)}'
    cip_ann, map_ann = lbl.annotations[0], lbl.annotations[1]

    assert cip_ann.fill == style.atom.default_colour, \
        f'CIP fill should be default_colour, got {cip_ann.fill!r}'
    assert map_ann.fill == style.label.map_colour, \
        f'map fill should be map_colour, got {map_ann.fill!r}'

    # Placement: the N's one bond goes right, so the column is turned 120 degrees off it -- to the LOWER
    # of the two turns, which is down and left -- rather than set straight back along the bond, where the
    # hydrogens are written and where a mark reads as the chain continuing.  dx is negative, so 'end'.
    assert cip_ann.anchor == 'end', f'expected end anchor, got {cip_ann.anchor!r}'
    # The INK stops at the box, not the anchor: a glyph's right side bearing is not clearance.  A diagonal
    # slide stops at the FIRST axis that separates, and for a turn this steep that is y -- so each row
    # hangs under the label rather than reaching around it, and each stops at its own ink.
    inks = [ann.bounds for ann in lbl.annotations]
    assert all(ink.max_y == approx(lbl.box.min_y) for ink in inks), \
        f'each row\'s ink should stop at box.min_y={lbl.box.min_y}, got {[tuple(i) for i in inks]}'
    assert all(ink.max_x < lbl.box.max_x for ink in inks), 'the turn is to the left, so is the column'

    cip_size = style.label.size * style.label.stereo_scale
    map_size = style.label.size * style.label.map_scale
    assert cip_ann.runs[0].size == approx(cip_size), \
        f'CIP run size {cip_ann.runs[0].size} should be size*stereo_scale={cip_size}'
    assert map_ann.runs[0].size == approx(map_size), \
        f'map run size {map_ann.runs[0].size} should be size*map_scale={map_size}'


def test_the_stereo_row_and_the_map_row_do_not_overlap():
    """THE regression: both annotations were placed at the same x, anchor and y=baseline

    Measured ink boxes and not just "the y values differ": two rows .001 apart also differ, and the claim
    is that a reader can tell the descriptor from the map number.  The hand-built plane fixes the flip so
    the side is not the variable under test.

    WHICH row ends up outside is the sector's business and not this test's: the map number is the row that
    keeps its own distance and the descriptor is stacked beyond it, so a column sliding along a diagonal
    ends staggered rather than stacked.  What is asserted is the separation, on the measured ink.
    """
    mol = smiles('[NH:3]CC')
    n_sid = [a.n for a in mol.atoms() if a.atomic_symbol == 'N'][0]
    c1, c2 = [a.n for a in mol.atoms() if a.atomic_symbol == 'C']
    with mol.edit():
        mol.set_atom_cip(n_sid, 'R')
    plane = {n_sid: (0., 0.), c1: (1., 0.), c2: (2., 0.)}
    style = DepictStyle().tuned(**{'atom.stereo_labels': True, 'atom.map_numbers': True})
    stereo, mapping = labels(mol, plane, style)[n_sid].annotations

    assert (stereo.x, stereo.y) != (mapping.x, mapping.y), 'the two rows are set at one point again'
    assert _clearance(stereo.bounds, mapping.bounds) > .09, \
        f'the two are not a reader apart: stereo {tuple(stereo.bounds)}, map {tuple(mapping.bounds)}'


def test_each_annotation_row_is_fixed_and_not_negotiated():
    """a map number must not move because the atom also carries a descriptor

    Were the rows chosen from what is present, one compound would be drawn two ways -- and the overlap
    test above passes on exactly that implementation, which is why this one is beside it.
    """
    mol = smiles('[NH:3]CC')
    n_sid = [a.n for a in mol.atoms() if a.atomic_symbol == 'N'][0]
    c1, c2 = [a.n for a in mol.atoms() if a.atomic_symbol == 'C']
    with mol.edit():
        mol.set_atom_cip(n_sid, 'R')
    plane = {n_sid: (0., 0.), c1: (1., 0.), c2: (2., 0.)}
    both = DepictStyle().tuned(**{'atom.stereo_labels': True, 'atom.map_numbers': True})

    with_stereo = labels(mol, plane, both)[n_sid].annotations[1]
    alone, = labels(mol, plane, both.tuned(**{'atom.stereo_labels': False}))[n_sid].annotations
    assert alone.y == approx(with_stereo.y), 'the map row moved when the descriptor was withheld'
    assert alone.x == approx(with_stereo.x)


def test_a_descriptor_beside_the_number_does_not_push_the_number_out():
    """THE enhanced-stereo half of "stay beside the atom": a wide `(R)&1` must not carry the number with it

    The rows are slid one at a time and the number goes down first, so its distance is its own ink's demand
    and nothing else's.  A single slide for the whole column -- measured from the union of its ink, which
    `(R)` is three times the width of -- put the number further out on exactly the centres that carry a
    descriptor, drawing one series two ways.

    The sector here points straight up, along the very axis the rise and the drop stack the rows on, which
    is the only case where the two rows interact at all: the descriptor has to go outside the number.  It
    takes TWO bonds to aim a sector straight up -- a terminal atom's column is turned 120 degrees off its
    one bond and so always leaves at an angle -- so the plane is the pair below, the usual chain vertex.
    """
    mol = smiles('C[NH:7]C')
    n_sid = [a.n for a in mol.atoms() if a.atomic_symbol == 'N'][0]
    c1, c2 = [a.n for a in mol.atoms() if a.atomic_symbol == 'C']
    with mol.edit():
        mol.set_atom_cip(n_sid, 'R')
    # both neighbours below at 120 degrees to each other, so the widest sector's bisector is straight up
    plane = {n_sid: (0., 0.), c1: (-.5, -.866), c2: (.5, -.866)}
    both = DepictStyle().tuned(**{'atom.stereo_labels': True, 'atom.map_numbers': True})

    stereo, mapping = labels(mol, plane, both)[n_sid].annotations
    alone, = labels(mol, plane, both.tuned(**{'atom.stereo_labels': False}))[n_sid].annotations
    assert (alone.x, alone.y) == approx((mapping.x, mapping.y)), \
        'the descriptor moved the map number'
    assert stereo.bounds.min_y >= mapping.bounds.max_y, \
        'the descriptor must stack OUTSIDE the number, which is the row that keeps its place'
    assert stereo.bounds.min_y - mapping.bounds.max_y == approx(both.label.pad), \
        'and one pad clear of it: tight ink boxes set (R) against 7 with nothing between them'


def test_the_two_rows_of_one_atom_never_overlap_however_the_sector_points():
    """the price of sliding the rows separately, paid over a corpus and not on one hand-built plane

    Every sector is tried on some atom of a fused or bridged skeleton, including the ones pointing along the
    axis the rows are stacked on, where independent slides are what would stack them.
    """
    style = DepictStyle().tuned(**{'atom.stereo_labels': True, 'atom.map_numbers': True})
    pairs = 0
    for smi in CROWDED_STEREO:
        mol = _all_mapped(smi)
        for sid, lbl in labels(mol, mol.coordinates(), style).items():
            if len(lbl.annotations) < 2:
                continue
            pairs += 1
            stereo, mapping = lbl.annotations
            assert _clearance(stereo.bounds, mapping.bounds) > style.label.pad - 1e-9, \
                f'{smi} atom {sid}: the descriptor and the number are not one pad apart'
    assert pairs > 20, f'premise: the corpus really carries two-row columns, got {pairs}'


def _group_marks(smi):
    """The stereo-row text of every atom that got one, keyed by atom, at the default style."""
    mol = smiles(smi)
    mol.clean2d()
    out = {}
    for sid, lbl in labels(mol, mol.coordinates(), DepictStyle()).items():
        for ann in lbl.annotations:
            text = ''.join(run.text for run in ann.runs)
            if not text.isdigit():      # the map row is the bare number; this wants the stereo row
                out[sid] = text
    return mol, out


def test_an_and_group_is_written_as_an_ampersand_and_its_number():
    """the CXSMILES vocabulary, so the picture and the `|&1:...|` the file carried agree"""
    mol, marks = _group_marks('C[C@H](N)[C@H](O)C |&1:1,3|')
    assert mol.stereo_groups(), 'premise: the notation set an AND group'
    assert sorted(marks.values()) == ['&1', '&1'], marks


def test_an_or_group_is_written_as_an_o_and_its_number():
    """two collections at once, so neither mark can be a constant"""
    mol, marks = _group_marks('C[C@H](N)[C@H](O)C |o1:1,&2:3|')
    assert sorted(marks.values()) == ['&2', 'o1'], marks


def test_an_abs_centre_is_written_as_a_bare_a():
    """ABS carries group 0, so a numbered mark would read `a0` and mean nothing"""
    mol, marks = _group_marks('C[C@H](N)O |a:1|')
    assert list(marks.values()) == ['a'], marks


def test_an_atom_in_no_stereo_group_gets_no_mark():
    """the control: a plain stereocentre is UNSPECIFIED, not ABS, so nothing is written

    Without this, an `a` on every centre a container merely happens to know about scores as correct.
    """
    mol, marks = _group_marks('C[C@H](N)O')
    assert not mol.has_stereo_groups, 'premise: a bare @ sets no collection'
    assert marks == {}, marks


def test_the_cip_descriptor_and_the_group_read_as_one_line():
    """one `Text`, two runs: `(R)&1` is one statement about one centre, and the row holds one

    Also the shape assertion -- a second `Text` here is what would put three annotations into two
    corners -- and the italic split, because a group id is a label and not a descriptor.
    """
    mol = smiles('C[C@H](N)[C@H](O)C |&1:1,3|')
    mol.clean2d()
    sid = mol.stereo_groups()[(3, 1)][0]
    with mol.edit():
        mol.set_atom_cip(sid, 'R')
    style = DepictStyle().tuned(**{'atom.stereo_labels': True})
    row, = labels(mol, mol.coordinates(), style)[sid].annotations
    assert ''.join(run.text for run in row.runs) == '(R)&1'
    assert [run.style for run in row.runs] == ['italic', 'normal']
    assert len({run.size for run in row.runs}) == 1, 'one row, one size'


def test_an_annotation_on_a_BARE_VERTEX_is_set_beside_it_and_not_on_it():
    """`Label`'s docstring promises "beside", and a degenerate box's corner IS the vertex

    An unlabelled atom's box collapses to its own point, so a column that stopped at the box edge would sit
    on the converging bond lines; the gap spent is `label.pad`, the clearance a labelled atom's box already
    carries, so no new style field is invented.  BOTH BRANCHES, because an assertion on the bare atoms
    alone also passes when the pad is spent on every atom and every descriptor moves.

    Stated on the measured INK and in whichever direction the annotation went: the column is set into the
    freest sector around the atom, so the axis it clears on is the layout's business, not this test's.
    """
    mol = smiles('[CH3:1][CH2:2][OH:3]')
    mol.clean2d()
    plane = mol.coordinates()
    style = DepictStyle().tuned(**{'atom.map_numbers': True})
    pad = style.label.pad
    assert pad > 0., 'premise: the style has a pad to spend'
    result = labels(mol, plane, style)

    bare = sorted(sid for sid, lbl in result.items() if lbl.text is None)
    assert len(bare) == 2, 'premise: both carbons draw as bare vertices'
    for sid in bare:
        lbl = result[sid]
        ann, = lbl.annotations
        assert lbl.box.width == 0. and lbl.box.height == 0., \
            'premise: a bare vertex box is degenerate at the atom point'
        # `label.pad` around the point IS the box a bare vertex's annotation is set outside of, so the
        # statement is that the ink is outside that box and against it.  Not "exactly touching": the
        # column slides along its sector's direction and stops at the first axis that separates, and a
        # diagonal slide can carry the other axis a little past as well -- bounded here by one more pad.
        gap = _clearance(ann.bounds, lbl.box.inflate(pad))
        assert 0. <= gap <= pad, \
            f'atom {sid}: the annotation must be one label pad clear of the vertex, not on it, got {gap}'

    labelled = [sid for sid, lbl in result.items() if lbl.text is not None]
    assert len(labelled) == 1, 'premise: the hydroxyl oxygen is the one written atom'
    lbl = result[labelled[0]]
    ann, = lbl.annotations
    assert _clearance(ann.bounds, lbl.box) == approx(0.), \
        'a written label\'s own box already holds the pad; spending it twice is a second gap'


def _clearance(box, obstacle):
    """The gap between two boxes on the axis that separates them; negative when they overlap.

    The LARGER of the two axes' gaps, because separation on one axis is separation: a column set to the
    left of a label clears it in x by the gap and in y not at all, and the gap is the answer wanted.
    """
    return max(obstacle.min_x - box.max_x, box.min_x - obstacle.max_x,
               obstacle.min_y - box.max_y, box.min_y - obstacle.max_y)


# Fused, spiro and crowded: every one of these has an atom with a bond in every sector, which is where a
# placement that is allowed to look further out starts walking.
CROWDED = ('c1ccc2ccccc2c1', 'c1ccc2cc3ccccc3cc2c1', 'C1CC2(CC1)CCC2',
           'CC(C)(C)c1ccc(cc1)C(=O)Nc1ccc2ncccc2c1', 'CN1C=NC2=C1C(=O)N(C)C(=O)N2C',
           'OCC1OC(O)C(O)C(O)C1O', 'CC(C)Cc1ccc(cc1)C(C)C(=O)O',
           'C1=CC2=C(C=C1)C(=O)c1ccccc1C2=O')


def _all_mapped(smi):
    """Every atom numbered, so the crowded positions are covered and not just the convenient ones."""
    mol = smiles(smi)
    ids = list(mol)                      # reads inside an open edit raise, so the ids come first
    with mol.edit() as e:
        for i, sid in enumerate(ids, 1):
            e.set_map_number(sid, i)
    mol.clean2d()
    return mol


# The same crowding, with enhanced-stereo collections on it: sugars, a terpenoid, a steroid and two
# bridged bicyclics, all public compounds, so `&N`/`oN`/`a` marks are really drawn beside the numbers.
CROWDED_STEREO = ('OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O |a:2,4,6,8,9|',
                  'OC[C@H]1O[C@@H](O)[C@H](O)[C@@H]1O |&1:2,4,6,7|',
                  'C[C@H]1CC[C@H](C(C)C)CC1 |&1:1,4|',
                  'C[C@]12CC[C@H]3[C@@H](CC[C@@H]4CC(=O)CC[C@]34C)[C@@H]1CCC2 |a:2,6,7,10,15,20|',
                  'C[C@H](O)[C@@H](C)[C@H](C)O |o1:1,3,5|',
                  'O[C@H]1CC[C@@H](O)CC1 |&1:1,4|',
                  'C1[C@H]2CC[C@@H]1CC2 |a:1,4|')


def test_AN_ANNOTATION_NEVER_WALKS_AWAY_FROM_ITS_ATOM():
    """the sector search chooses a SIDE, never a distance -- one crowded sector is not grounds to look out

    The regression this pins is a placement that answered a bond in every sector by stepping the number
    further out until it found room.  It found room; the numbers it moved were up to a whole bond length
    from the atom they name, which is a worse defect than the line they were moved off -- a number nobody
    can attribute is not legible.  `figure.py` answers the crowded atom with a knock-out plate instead.

    Stated twice, because either half alone is weak: the gap is at most the pad the vertex test already
    measures, AND the ink's centre is nearer its own atom than any other atom in the drawing.
    """
    style = DepictStyle().tuned(**{'atom.map_numbers': True})
    pad = style.label.pad
    numbers = 0
    for smi in CROWDED:
        mol = _all_mapped(smi)
        plane = mol.coordinates()
        for sid, lbl in labels(mol, plane, style).items():
            ann, = lbl.annotations
            numbers += 1
            # A bare vertex's box is its own point, so the pad is the box its annotation is set outside of;
            # a written label's box already carries the pad -- the same two branches as the vertex test.
            gap = _clearance(ann.bounds, lbl.box.inflate(pad) if lbl.text is None else lbl.box)
            assert gap <= pad, f'{smi} atom {sid}: the number is {gap} out, more than one pad'

            box = ann.bounds
            centre = ((box.min_x + box.max_x) / 2., (box.min_y + box.max_y) / 2.)
            own = hypot(plane[sid][0] - centre[0], plane[sid][1] - centre[1])
            for other, (x, y) in plane.items():
                if other != sid:
                    assert hypot(x - centre[0], y - centre[1]) >= own, \
                        f'{smi}: atom {sid}\'s number sits nearer atom {other}'
    assert numbers > 100, f'premise: the corpus is worth measuring, got {numbers} numbers'


def test_an_unknown_hydrogen_count_can_be_marked_ON_A_CARBON():
    """the `?` marker's own element: `H_UNKNOWN` arrives on carbon more than anywhere else

    The marker lives in `_compose`, which runs only for a LABELLED atom, so `is_labelled` must write a
    carbon for an underivable count as it does for a charge or an isotope -- otherwise a CH3 vertex is
    drawn where the container says "not derivable".  Three assertions, since the marker is a conditional
    in two places: the marked carbon, a neighbour whose count IS known, and the option off.
    """
    mol = smiles('CC')
    unknown, known = sorted(a.n for a in mol.atoms())
    with mol.edit():
        mol.set_hydrogens(unknown, H_UNKNOWN)
    mol.clean2d()
    plane = mol.coordinates()
    marks = DepictStyle().tuned(**{'atom.unknown_h_marks': True})
    assert mol.atom(unknown).implicit_h is None, 'premise: the count is not derivable'
    assert mol.atom(known).implicit_h is not None, 'premise: the other one is'

    assert is_labelled(mol.atom(unknown), marks) is True
    assert is_labelled(mol.atom(known), marks) is False, \
        'a carbon whose count IS known is still a bare vertex'
    assert is_labelled(mol.atom(unknown), DepictStyle()) is False, \
        'and without the option there is no mark to write, so no symbol either'

    result = labels(mol, plane, marks)
    assert result[unknown].text is not None, 'the marked carbon must get a written label'
    assert '?' in ''.join(run.text for run in result[unknown].text.runs)
    assert result[known].text is None, 'and nothing else in the picture changes'
    assert labels(mol, plane, DepictStyle())[unknown].text is None


def test_labels_refuses_a_molecule_with_no_layout():
    """not "draws it at the origin": coordinates() on an un-laid-out molecule returns {}, so
    the guard fires on the first missing key rather than the degenerate-span check"""
    mol = smiles('CCO')
    with raises(ValueError, match='layout'):
        labels(mol, mol.coordinates(), DepictStyle())


def test_labels_refuses_a_plane_missing_an_atom():
    """caller contract: every atom in mol must have a key in plane"""
    mol = smiles('CCO')
    mol.clean2d()
    plane = mol.coordinates()
    sid = next(iter(mol))
    bad_plane = {k: v for k, v in plane.items() if k != sid}
    with raises(ValueError, match=str(sid)):
        labels(mol, bad_plane, DepictStyle())


def test_labels_refuses_a_degenerate_plane():
    """all atoms at one point is not a layout"""
    mol = smiles('CCO')
    bad_plane = {sid: (0.0, 0.0) for sid in mol}
    with raises(ValueError, match='degenerate'):
        labels(mol, bad_plane, DepictStyle())


def test_the_degenerate_span_is_the_ARENA_S_threshold_and_carries_a_name():
    """`labels()` must refuse exactly the planes the arena itself calls "no layout"

    The literal is a fact about the store and not a rendering preference -- the same threshold
    `has_layout` applies in the arena's fixed-point check -- which is what exempts it from "no number in a
    drawing module that is not from `DepictStyle`".  THE EXPECTED VALUE COMES FROM THE CORE: probe planes
    built out of `_DEGENERATE_SPAN` would move with the literal and pass for any value of it.
    """
    mol = smiles('CC')
    a, b = sorted(mol)
    style = DepictStyle()
    assert isinstance(_DEGENERATE_SPAN, float), 'the threshold is a named module constant'
    for span in (0., 0.001, 0.005, 0.009, 0.011, 0.05, 1.):
        plane = {a: (0., 0.), b: (span, 0.)}
        with mol.edit():
            mol.set_xy(a, 0., 0.)
            mol.set_xy(b, span, 0.)
        drawable = mol.has_layout
        try:
            assert set(labels(mol, plane, style)) == {a, b}
            refused = False
        except ValueError as exc:
            assert 'degenerate' in str(exc)
            refused = True
        assert refused is not drawable, \
            f'span {span}: the arena says has_layout={drawable} and labels() ' \
            f'{"refused" if refused else "accepted"} it -- one threshold, two answers'
