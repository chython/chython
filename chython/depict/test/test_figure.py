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
"""Assembly: a laid-out molecule becomes a Scene, and `mol.depict()` returns SVG.

Painter's order is asserted here and nowhere else, being a property of the assembly rather than of any one
module: highlights and fields UNDER the structure, bonds next, labels last.  Two rulings are pinned here
too -- `depict()` stores nothing but logs the layout it computed, and a query is not depictable.
"""
from re import findall
from xml.etree.ElementTree import fromstring
from pytest import approx, raises
from chython import ReactionContainer, smarts, smiles
from chython.depict.bonds import bond_paths
from chython.depict.figure import molecule_scene, reaction_scene
from chython.depict.label import element_colour, labels
from chython.depict.layout import molecule as _molecule_layout, reaction as _reaction_layout
from chython.depict.scene import Group, Path, Scene, Text
from chython.depict.style import DepictStyle, set_depict_style
from chython.depict.wedge import wedge_paths


def test_a_molecule_becomes_a_scene():
    mol = smiles('CCO')
    mol.clean2d()
    scene = molecule_scene(mol)
    assert isinstance(scene, Scene)
    assert scene.children


def test_labels_are_drawn_after_bonds():
    """painter's order: a label under a bond is unreadable, and no amount of trimming fixes it"""
    mol = smiles('CCO')
    mol.clean2d()
    flat = _flatten(molecule_scene(mol))
    last_path = max(i for i, node in enumerate(flat) if isinstance(node, Path))
    first_text = min(i for i, node in enumerate(flat) if isinstance(node, Text))
    assert first_text > last_path


def test_labels_are_drawn_after_bonds_across_a_whole_reaction():
    """painter's order is a property of the FIGURE, not of each molecule in it

    Assembled per molecule, the second reactant's bonds land on top of the first one's labels; nothing
    overlaps in a left-to-right arrangement, but an overlay under the structure would.
    """
    rxn = ReactionContainer([smiles('CCO'), smiles('CC(=O)O')], [smiles('CCOC(C)=O')])
    rxn.clean2d()
    flat = _flatten(reaction_scene(rxn))
    last_path = max(i for i, node in enumerate(flat) if isinstance(node, Path))
    first_text = min(i for i, node in enumerate(flat) if isinstance(node, Text))
    assert first_text > last_path


def test_depict_returns_an_svg_document_that_parses():
    mol = smiles('c1ccccc1C(=O)O')          # benzoic acid
    mol.clean2d()
    fromstring(mol.depict())


def test_depict_draws_a_molecule_with_no_layout_and_stores_nothing():
    """a caller who wants a picture of a parsed SMILES should get one, not an exception

    And without the molecule changing underneath: `clean2d()` is the explicit call for a layout the caller
    means to KEEP.  The log line says a layout was computed for this drawing alone.
    """
    mol = smiles('CCO')
    assert not mol.has_layout
    log = []
    svg = mol.depict(log=log)
    assert svg.startswith('<svg')
    assert not mol.has_layout, 'drawing is not a mutation: a temporary layout is not stored'
    assert not mol.coordinates()
    assert [record.rule for record in log] == ['depict:layout']
    assert 'layout' in log[0].message
    assert log[0].atoms == tuple(mol)


def test_a_stored_layout_is_drawn_with_no_log_line_at_all():
    """the negative half of the test above: an unconditional `log.append` satisfies it otherwise"""
    mol = smiles('CCO')
    mol.clean2d()
    log = []
    mol.depict(log=log)
    assert log == []


def test_depict_can_draw_a_plane_it_was_given_without_storing_it():
    mol = smiles('CCO')
    plane = mol.layout2d()
    assert not mol.has_layout
    log = []
    mol.depict(plane=plane, log=log)
    assert not mol.has_layout, 'a plane passed in is not a plane stored'
    assert log == [], 'a caller who supplied the plane is told nothing about a layout'
    # a shifted copy differs from the layout the engine would recompute, so a dropped `plane` argument
    # makes the children diverge
    shifted = {sid: (x + 10., y) for sid, (x, y) in plane.items()}
    assert mol.scene(plane=shifted).children == molecule_scene(mol, plane=shifted).children


def test_the_scene_is_not_cached_on_the_molecule():
    """two styles, two pictures, in either order -- a cached scene would return the first twice"""
    mol = smiles('CCO')
    mol.clean2d()
    thin = mol.depict(style=DepictStyle().tuned(**{'bond.width': .03}))
    thick = mol.depict(style=DepictStyle().tuned(**{'bond.width': .09}))
    assert thin != thick
    assert mol.depict(style=DepictStyle().tuned(**{'bond.width': .03})) == thin


def test_repr_svg_uses_the_process_default_style():
    mol = smiles('CCO')
    mol.clean2d()
    original = DepictStyle()
    try:
        set_depict_style(DepictStyle().tuned(**{'bond.width': .123}))
        assert 'stroke-width="0.123"' in mol._repr_svg_()
    finally:
        set_depict_style(original)


def test_a_reaction_repr_svg_uses_the_process_default_style():
    rxn = ReactionContainer([smiles('CCO')], [smiles('CC=O')])
    rxn.clean2d()
    original = DepictStyle()
    try:
        set_depict_style(DepictStyle().tuned(**{'bond.width': .123}))
        assert 'stroke-width="0.123"' in rxn._repr_svg_()
    finally:
        set_depict_style(original)


def test_a_query_is_not_depictable_and_says_so():
    """a SMARTS has no label rule, so drawing one would put guessed notation on a printed page

    `[C,N,O]` has no atomic symbol and `[!R;D2]` has no element at all.  Deferred deliberately, and the
    REFUSAL is what is pinned so the gap stays visible.
    """
    query = smarts('[N;D1;z1;x0:1][C;z1:2]')
    with raises(AttributeError):
        query.depict()
    with raises(AttributeError):
        query.clean2d()


def test_a_disconnected_molecule_draws_both_components():
    mol = smiles('CCO.CCN')
    mol.clean2d()
    svg = mol.depict()
    assert svg.count('<text') >= 2


def test_a_reaction_becomes_one_scene_with_an_arrow():
    rxn = ReactionContainer([smiles('CCO')], [smiles('CC=O')])
    rxn.clean2d()
    scene = reaction_scene(rxn)
    flat = _flatten(scene)
    assert any(isinstance(node, Path) and node.fill for node in flat), 'the arrow head is filled'


def test_the_arrow_spans_the_gap_the_layout_measured():
    """the arrow is drawn where `layout2d` put it, at the length the style states, and it points RIGHT

    Every number here is one a wrong arrow still looks like an arrow with: a head on the left end draws the
    retrosynthesis, a shaft from the origin runs back through the last reactant, and a head measured from
    `x1` is a triangle the length of the whole gap.
    """
    rxn = ReactionContainer([smiles('CCO')], [smiles('CC=O')])
    rxn.clean2d()
    _, arrow, _ = rxn.layout2d()
    x1, x2, y = arrow
    style = DepictStyle()
    # the arrow is the LAST thing appended under the type, so its two paths close the path list
    shaft, head = [node for node in _flatten(reaction_scene(rxn, style=style))
                   if isinstance(node, Path)][-2:]

    assert shaft.stroke and not shaft.fill, 'the shaft is stroked'
    assert shaft.subpaths[0][0][1] == approx(x1), 'the shaft starts where the layout said the gap does'
    assert shaft.subpaths[0][-1][1] == approx(x2 - style.reaction.head_length)

    tip = head.subpaths[0][0]                      # the M of the triangle is its point: ('M', x, y)
    assert head.fill and not head.stroke, 'the head is filled'
    assert tip[2] == approx(y)
    assert tip[1] == approx(x2), 'the head sits at the far end of the span the layout returned'
    xs = [segment[1] for segment in head.subpaths[0] if segment[0] != 'Z']
    assert x2 - min(xs) == approx(style.reaction.head_length), 'the head is a head, not the whole gap'
    # all four ReactionStyle fields are read by _arrow_paths
    ys = [segment[2] for segment in head.subpaths[0] if segment[0] != 'Z']
    assert max(ys) - min(ys) == approx(style.reaction.head_width)
    assert shaft.width == approx(style.reaction.arrow_width)
    assert shaft.stroke == style.reaction.colour
    assert head.fill == style.reaction.colour
    assert head.subpaths[0][-1][0] == 'Z', 'the head is a closed polygon: open leaves it unfilled'


def test_the_arrow_colour_is_taken_from_the_style():
    """shaft stroke and head fill both read `ReactionStyle.colour`, not a hardcoded value

    The default reaction colour is '#000000', which a hardcoded black satisfies, so the style is tuned.
    """
    rxn = ReactionContainer([smiles('CCO')], [smiles('CC=O')])
    rxn.clean2d()
    style = DepictStyle().tuned(**{'reaction.colour': '#0000ff'})
    shaft, head = [node for node in _flatten(reaction_scene(rxn, style=style))
                   if isinstance(node, Path)][-2:]
    assert shaft.stroke == style.reaction.colour
    assert head.fill == style.reaction.colour


def test_a_reaction_draws_a_plus_between_two_reactants():
    rxn = ReactionContainer([smiles('CCO'), smiles('CC(=O)O')], [smiles('CCOC(C)=O')])
    rxn.clean2d()
    style = DepictStyle()
    _, _, signs = rxn.layout2d()
    texts = [node for node in _flatten(reaction_scene(rxn, style=style)) if isinstance(node, Text)]
    assert any('+' in ''.join(run.text for run in node.runs) for node in texts)
    # all four _sign_texts parameters are wired through the style
    plus_nodes = [n for n in texts if any('+' in r.text for r in n.runs)]
    assert len(plus_nodes) == len(signs)
    x, y = signs[0]
    plus = plus_nodes[0]
    assert plus.x == approx(x)
    assert plus.y == approx(y - style.label.baseline_shift * style.reaction.sign_size)
    assert plus.anchor == 'middle'
    assert plus.runs[0].size == approx(style.reaction.sign_size)


def test_one_reactant_and_one_product_draw_no_plus():
    """the premise of the test above: a `+` per GAP, not one per figure"""
    rxn = ReactionContainer([smiles('CCO')], [smiles('CC=O')])
    rxn.clean2d()
    texts = [node for node in _flatten(reaction_scene(rxn)) if isinstance(node, Text)]
    assert not any('+' in ''.join(run.text for run in node.runs) for node in texts)


def test_an_empty_reaction_side_does_not_crash():
    rxn = ReactionContainer([smiles('CCO')], [])
    rxn.clean2d()
    fromstring(rxn.depict())


def test_a_dropped_bond_reaches_the_callers_log():
    """the log threads all the way out: a picture that silently lost a bond is worse than an error

    Measured, not assumed: at `label.size 2.5` the five glyphs of sulfuric acid swallow all four of its
    bonds and `bonds.py` files a `depict:crowded` record for each.
    """
    mol = smiles('OS(=O)(=O)O')
    mol.clean2d()
    log = []
    mol.depict(style=DepictStyle().tuned(**{'label.size': 2.5, 'label.pad': .5}), log=log)
    assert [record.rule for record in log] == ['depict:crowded'] * 4


def test_a_stereo_bond_the_wedge_chooser_could_not_honour_reaches_the_callers_log():
    """the OTHER log this module threads: `wedge_paths(..., log=log)`, not `bond_paths`

    Two independent `log=` arguments, so one can be dropped without the crowded-bond test noticing.
    (2E,4Z)-hexadiene's two configured double bonds cannot both be drawn as declared from these
    coordinates, so `wedge.py` files `depict:crossed-contradiction`, a rule id `bonds.py` never emits.
    """
    mol = smiles('C/C=C/C=C\\C')
    mol.clean2d()
    log = []
    mol.depict(log=log)
    assert [record.rule for record in log] == ['depict:crossed-contradiction']


def test_a_reactions_log_reaches_the_caller_too():
    """the reaction entry points thread `log` as well, and the member that reported it is identifiable

    `rxn.depict()` and `rxn.scene()` are separate bodies with their own argument lists, either of which
    could drop `log=`.  The `+` and the arrow file nothing, so what arrives is the members' own records.
    """
    rxn = ReactionContainer([smiles('C/C=C/C=C\\C')], [smiles('CCCCCC')])
    rxn.clean2d()
    for draw in (rxn.depict, rxn.scene):
        log = []
        draw(log=log)
        assert [record.rule for record in log] == ['depict:crossed-contradiction'], draw.__name__


def test_a_roomy_reaction_reports_nothing():
    """the negative half: the reaction furniture itself files nothing"""
    rxn = ReactionContainer([smiles('CCO'), smiles('CC(=O)O')], [smiles('CCOC(C)=O')])
    rxn.clean2d()
    log = []
    rxn.depict(log=log)
    assert log == []


def test_a_roomy_molecule_reports_nothing():
    """and the same molecule at a sane size logs nothing, or the assertion above tests nothing"""
    mol = smiles('OS(=O)(=O)O')
    mol.clean2d()
    log = []
    mol.depict(log=log)
    assert log == []


def test_a_wedge_bond_is_not_drawn_twice():
    """`wedge_paths` returns the keys it claimed and the assembly MUST pass them on as `skip=`

    Otherwise the wedge is drawn over a plain line along the same axis: at print weight the line
    escapes past the triangle's edges and the figure has a wedge with whiskers.
    """
    mol = smiles('C[C@H](N)O')
    mol.clean2d()
    plane = mol.coordinates()
    style = DepictStyle()
    boxes = labels(mol, plane, style)
    wedges, claimed = wedge_paths(mol, plane, boxes, style)
    assert claimed, 'this molecule draws no wedge, so the test has lost its subject'
    with_skip = bond_paths(mol, plane, boxes, style, skip=claimed)
    without = bond_paths(mol, plane, boxes, style)
    assert len(without) > len(with_skip), 'the skip makes no difference here; pick another molecule'

    drawn = [node for node in _flatten(molecule_scene(mol)) if isinstance(node, Path)]
    assert len(drawn) == len(wedges) + len(with_skip)


def test_a_radical_gets_a_dot_on_its_own_atom():
    """the dot is geometry from the label's box, not a glyph, and it is clear of the ink

    `Label.anchor` is the atom's point, so the dot sits on the anchor's x; the label's cap is centred on
    that point, so a dot AT the anchor is inside the symbol.  `AtomStyle.radical_gap` sets the clearance.
    """
    mol = smiles('CC |^1:1|')
    mol.clean2d()
    style = DepictStyle()
    plane = mol.coordinates()
    boxes = labels(mol, plane, style)
    radical = next(atom.n for atom in mol.atoms() if atom.is_radical)
    label = boxes[radical]

    flat = _flatten(molecule_scene(mol))
    dots = [node for node in flat
            if isinstance(node, Path) and node.fill and len(node.subpaths[0]) == 6]
    assert len(dots) == 1
    start = dots[0].subpaths[0][0]                  # the M of `circle()`, at (cx + r, cy)
    assert start[1] - style.atom.radical_radius == approx(label.anchor[0])
    assert start[2] == approx(label.box.max_y + style.atom.radical_gap)
    # the dot is the last path in painter's order -- over the bonds, under the type
    non_dot_path_idxs = [i for i, n in enumerate(flat) if isinstance(n, Path) and n is not dots[0]]
    assert non_dot_path_idxs, 'the molecule draws no bonds, so the ordering assertion has no subject'
    assert flat.index(dots[0]) > max(non_dot_path_idxs)


def test_a_non_radical_gets_no_dot():
    """the half that makes the test above mean something"""
    mol = smiles('CC')
    mol.clean2d()
    dots = [node for node in _flatten(molecule_scene(mol))
            if isinstance(node, Path) and node.fill]
    assert not dots


def test_a_style_that_withholds_radicals_draws_no_dot():
    """`AtomStyle.radicals` off is a picture that does not CLAIM the radical, so it must not mark it

    The same radical as above, so the style is the only thing that changed.  The mark is a chemical
    statement, so drawing it against the setting is worse than ignoring an ordinary preference.
    """
    mol = smiles('CC |^1:1|')
    mol.clean2d()
    style = DepictStyle().tuned(**{'atom.radicals': False})
    dots = [node for node in _flatten(molecule_scene(mol, style=style))
            if isinstance(node, Path) and node.fill]
    assert not dots


def test_annotations_reach_the_scene():
    """`Label.annotations` is a field, and a field nothing reads is a defect

    Map numbers, because nothing in chython 3 COMPUTES a CIP descriptor (`set_atom_cip` is storage only).
    Methanol maps to two atoms and only the oxygen is labelled, so one symbol and two annotations; its
    spelling (`OH` or `HO`) belongs to `test_label.py`.
    """
    mol = smiles('[CH3:1][OH:2]')
    mol.clean2d()
    style = DepictStyle().tuned(**{'atom.map_numbers': True})
    texts = [node for node in _flatten(molecule_scene(mol, style=style)) if isinstance(node, Text)]
    written = [''.join(run.text for run in node.runs) for node in texts]
    assert len(written) == 3, written
    assert sorted(t for t in written if t.isdigit()) == ['1', '2']


def test_a_labels_symbol_precedes_its_annotations():
    """within each label's contribution to the type layer, the symbol comes before its annotations

    Swapping the two appends in `_molecule_nodes` leaves the presence-and-count test above green, but a
    symbol drawn after its annotation is invisible under the descriptor beside it.
    """
    mol = smiles('[NH2:1][OH:2]')
    mol.clean2d()
    style = DepictStyle().tuned(**{'atom.map_numbers': True})
    flat = _flatten(molecule_scene(mol, style=style))
    texts = [node for node in flat if isinstance(node, Text)]
    # a map annotation is the bare number, so a digit-only text is what tells the two layers apart
    symbols = [n for n in texts if not ''.join(r.text for r in n.runs).isdigit()]
    annots = [n for n in texts if ''.join(r.text for r in n.runs).isdigit()]
    assert len(symbols) == 2 and len(annots) == 2, texts
    for sym, ann in zip(symbols, annots):
        assert flat.index(sym) < flat.index(ann), 'symbol must precede its own annotation'


def test_a_mapped_molecule_shows_its_mapping_with_no_style_argument():
    """the user-facing claim: `mol.depict()` on a mapped structure shows the mapping"""
    mol = smiles('[CH3:1][OH:2]')
    mol.clean2d()
    texts = [node for node in _flatten(molecule_scene(mol)) if isinstance(node, Text)]
    written = [''.join(run.text for run in node.runs) for node in texts]
    assert sorted(t for t in written if t.isdigit()) == ['1', '2'], written


def test_an_annotation_the_style_withholds_is_not_drawn():
    """the control for the test above, which would pass on a figure that always drew the mapping"""
    mol = smiles('[CH3:1][OH:2]')
    mol.clean2d()
    style = DepictStyle().tuned(**{'atom.map_numbers': False})
    texts = [node for node in _flatten(molecule_scene(mol, style=style)) if isinstance(node, Text)]
    written = [''.join(run.text for run in node.runs) for node in texts]
    assert len(written) == 1, written                # the hydroxyl symbol alone, and no `1`/`2`
    assert not [t for t in written if t.isdigit()]


def test_a_racemic_centre_reaches_the_scene():
    """`Label.annotations` is flushed by `_molecule_nodes` at figure.py:201, so the mark must arrive"""
    mol = smiles('C[C@H](N)[C@H](O)C |&1:1,3|')
    mol.clean2d()
    texts = [node for node in _flatten(molecule_scene(mol)) if isinstance(node, Text)]
    written = [''.join(run.text for run in node.runs) for node in texts]
    assert sorted(t for t in written if t.startswith('&')) == ['&1', '&1'], written


def test_an_unmapped_molecule_stays_clean_with_no_style_argument():
    """and the other half of the ruling: nothing is drawn where nothing is mapped"""
    mol = smiles('CCO')
    mol.clean2d()
    texts = [node for node in _flatten(molecule_scene(mol)) if isinstance(node, Text)]
    assert not [n for n in texts if ''.join(r.text for r in n.runs).isdigit()]


def _mapped(text):
    """The molecule with 1..N on its atoms and a layout: what a mapped record looks like to the figure.

    Read before the edit session opens -- a container with pending edits refuses to be read -- and
    `map_number` is set rather than the ids renumbered, the two being separate fields.
    """
    mol = smiles(text)
    ids = list(mol)
    with mol.edit() as edit:
        for number, sid in enumerate(ids, 1):
            edit.set_map_number(sid, number)
    mol.clean2d()
    return mol


def _plates(flat):
    """The knock-out plates: filled, unstroked paths.  The fixtures below carry no wedge and no radical,
    which are the other two filled shapes a molecule can produce."""
    return [node for node in flat if isinstance(node, Path) and node.fill and node.stroke is None]


def test_EVERY_MAP_NUMBER_GETS_A_PLATE_UNDER_IT():
    """naphthalene, mapped: ten numbers, ten plates, and the LINE gives way at each of them

    Not the number's position -- `label.py` will not walk a number away from the atom it names, since half
    a bond out it stops naming it -- so what makes a number readable is a plate in the background colour
    between it and the drawing.  Unconditional: a plate over nothing is invisible, and a digit in the
    corridor between a ring's perimeter and its inner line touches neither line and is unreadable anyway.
    White because `page.background` is None: a knock-out on an unpainted page has to assume the white the
    figure will be put on.
    """
    flat = _flatten(molecule_scene(_mapped('c1ccc2c(c1)cccc2')))
    plates = _plates(flat)
    numbers = [n for n in flat if isinstance(n, Text) and ''.join(r.text for r in n.runs).isdigit()]
    assert len(numbers) == 10, numbers
    assert len(plates) == len(numbers), 'one plate per annotation, whatever the drawing does around it'
    assert all(plate.fill == '#ffffff' for plate in plates), [p.fill for p in plates]
    # A plate that does not cover the glyph it was drawn for is ink for nothing, so each one is matched
    # to a number whose ink it contains.
    for plate in plates:
        box = plate.bounds
        assert any(box.min_x <= ink.min_x and box.max_x >= ink.max_x
                   and box.min_y <= ink.min_y and box.max_y >= ink.max_y
                   for ink in (number.bounds for number in numbers)), 'a plate covers a whole number'


def test_A_STEREO_MARK_IS_PLATED_LIKE_A_NUMBER():
    """the other annotation row, on the fixtures whose marks land in the worst place there is

    A stereo statement is `(R)`, `&1`, `o1` or `a` set above the atom point, in the same column a map
    number is set below it -- so it is the same kind of mark over the same line work and gets the same
    knock-out.  A sugar has a mark on five adjacent centres and a wedge at each of them, so its stereo row
    is the crowded case; a plate that reached only the numbers would leave every `a` on a hash line.

    Asserted per annotation and by count: every stereo mark's ink is inside some plate, AND there are as
    many plates as there are annotations of both kinds, which is what "unconditional" means.
    """
    style = DepictStyle().tuned(**{'atom.map_numbers': True, 'atom.stereo_labels': True})
    for smi in ('OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O |a:2,4,6,8,9|',
                'C[C@H](N)[C@H](O)C |&1:1,3|', 'C[C@H](O)[C@@H](C)[C@H](C)O |o1:1,3,5|',
                'C1[C@H]2CC[C@@H]1CC2 |a:1,4|'):
        mol = _mapped(smi)
        annotations = [text for label in labels(mol, mol.coordinates(), style).values()
                       for text in label.annotations]
        marks = [text for text in annotations if not ''.join(r.text for r in text.runs).isdigit()]
        assert marks, f'{smi}: premise -- this fixture is drawn with stereo marks on it'
        plates = _plates(_flatten(molecule_scene(mol, style=style)))
        # A wedge is filled and unstroked too, so the count is a lower bound on this fixture; the
        # containment below is what pins each mark to a plate of its own.
        assert len(plates) >= len(annotations), f'{smi}: {len(plates)} plates for {len(annotations)} rows'
        for mark in marks:
            ink = mark.bounds
            assert any(plate.bounds.min_x <= ink.min_x and plate.bounds.max_x >= ink.max_x
                       and plate.bounds.min_y <= ink.min_y and plate.bounds.max_y >= ink.max_y
                       for plate in plates), \
                f'{smi}: {"".join(r.text for r in mark.runs)!r} is drawn on the structure with no plate'


def test_a_plate_is_drawn_over_the_structure_and_under_every_glyph():
    """the plate's place in painter's order, which is the whole of what makes it work

    Over the bonds or it knocks nothing out; under EVERY glyph and not just its own number, or a plate
    that happens to reach beneath a neighbouring symbol would erase the chemistry instead of the line.
    """
    flat = _flatten(molecule_scene(_mapped('CC(=O)Nc1ccccc1')))
    plates = [i for i, node in enumerate(flat)
              if isinstance(node, Path) and node.fill and node.stroke is None]
    strokes = [i for i, node in enumerate(flat) if isinstance(node, Path) and node.stroke is not None]
    texts = [i for i, node in enumerate(flat) if isinstance(node, Text)]
    assert plates and strokes and texts
    assert min(plates) > max(strokes), 'a plate under the bonds knocks nothing out'
    assert max(plates) < min(texts), 'a glyph must not be painted over by any plate'


def test_an_uncrowded_map_number_IS_PLATED_TOO():
    """ethanol: three numbers with a clear side each, and three plates all the same

    The plate is not measured against the drawing.  Over the page it knocks out nothing and is invisible,
    which is the whole reason it can be drawn always -- and the case a per-number test got wrong is the
    number a line comes NEAR rather than crosses, which is most of them in a ring.
    """
    plates = _plates(_flatten(molecule_scene(smiles('[CH3:1][CH2:2][OH:3]'))))
    assert len(plates) == 3, plates


def test_the_plate_can_be_withheld():
    """the one control on a layer that is otherwise always drawn: a figure with no knock-outs in it"""
    style = DepictStyle().tuned(**{'label.annotation_plate': 'none'})
    flat = _flatten(molecule_scene(_mapped('c1ccc2c(c1)cccc2'), style=style))
    assert not _plates(flat)
    assert [n for n in flat if isinstance(n, Text) and ''.join(r.text for r in n.runs).isdigit()], \
        'the numbers are still drawn; only the knock-out under them is gone'


def test_the_plate_takes_the_pages_background_colour():
    """a knock-out is the colour of what it knocks out, so a painted page changes it"""
    style = DepictStyle().tuned(**{'page.background': '#ffeecc'})
    plates = _plates(_flatten(molecule_scene(_mapped('c1ccc2c(c1)cccc2'), style=style)))
    assert plates and all(plate.fill == '#ffeecc' for plate in plates), [p.fill for p in plates]


def test_the_plate_colour_can_be_stated_outright():
    """for the figure whose background the page does not know -- a slide, a printed panel"""
    style = DepictStyle().tuned(**{'page.background': '#ffeecc',
                                   'label.annotation_plate_colour': '#112233'})
    plates = _plates(_flatten(molecule_scene(_mapped('c1ccc2c(c1)cccc2'), style=style)))
    assert plates and all(plate.fill == '#112233' for plate in plates), [p.fill for p in plates]


def test_the_plate_has_two_shapes_and_the_rounded_one_is_the_default():
    """`'rounded'` is a stadium -- lines and corner arcs -- and `'ellipse'` is four arcs and nothing else

    The default is the tighter of the two: an ellipse has to pass through the corners of the number's box
    to keep the digits inside it, so it covers half again as much of the drawing.
    """
    mol = _mapped('c1ccc2c(c1)cccc2')
    rounded = _plates(_flatten(molecule_scene(mol)))
    style = DepictStyle().tuned(**{'label.annotation_plate': 'ellipse'})
    elliptic = _plates(_flatten(molecule_scene(mol, style=style)))
    assert rounded and len(elliptic) == len(rounded)
    verbs = {seg[0] for plate in rounded for sub in plate.subpaths for seg in sub}
    assert verbs == {'M', 'L', 'C', 'Z'}, verbs
    verbs = {seg[0] for plate in elliptic for sub in plate.subpaths for seg in sub}
    assert verbs == {'M', 'C', 'Z'}, verbs
    for round_plate, ellipse_plate in zip(rounded, elliptic):
        assert ellipse_plate.bounds.width > round_plate.bounds.width


def test_a_radical_dot_uses_the_elements_cpk_colour():
    """the dot inherits the atom's CPK colour, not a hardcoded black

    CPK carbon is black, so a carbon radical cannot tell `element_colour` from '#000000'.  Oxygen is red.
    """
    mol = smiles('CO |^1:1|')      # methanol with the radical on oxygen (index 1)
    mol.clean2d()
    style = DepictStyle()
    radical_atom = next(atom for atom in mol.atoms() if atom.is_radical)
    dots = [node for node in _flatten(molecule_scene(mol, style=style))
            if isinstance(node, Path) and node.fill and len(node.subpaths[0]) == 6]
    assert len(dots) == 1
    assert dots[0].fill == element_colour(radical_atom, style)


def test_the_output_is_deterministic():
    mol = smiles('CC(=O)OC1=CC=CC=C1C(=O)O')     # aspirin
    mol.clean2d()
    assert mol.depict() == mol.depict()
    assert not findall(r'[0-9a-f]{8}-[0-9a-f]{4}-', mol.depict())


def test_an_acs_figure_is_the_stated_width():
    mol = smiles('CC(=O)OC1=CC=CC=C1C(=O)O')
    mol.clean2d()
    assert 'width="83mm"' in mol.depict(style=DepictStyle.preset('acs'))


def test_an_acs_reaction_figure_is_the_stated_width():
    """the reaction side's `style=` is wired to `to_svg()`, not just to `reaction_scene()`

    The mirror of the molecule-side test: dropping `style=` from the reaction's `to_svg()` call otherwise
    falls back to the process default and the figure comes out at the wrong physical size.
    """
    rxn = ReactionContainer([smiles('CCO')], [smiles('CC=O')])
    rxn.clean2d()
    assert 'width="83mm"' in rxn.depict(style=DepictStyle.preset('acs'))


def test_the_reaction_scene_is_not_cached():
    """two styles, two pictures, in either order -- the reaction mirror of the molecule-side test

    The difference is in the member drawings, so a reaction drawing them at the process default returns
    the same SVG for both widths.
    """
    rxn = ReactionContainer([smiles('CCO')], [smiles('CC=O')])
    rxn.clean2d()
    thin = rxn.depict(style=DepictStyle().tuned(**{'bond.width': .03}))
    thick = rxn.depict(style=DepictStyle().tuned(**{'bond.width': .09}))
    assert thin != thick
    assert rxn.depict(style=DepictStyle().tuned(**{'bond.width': .03})) == thin


def test_a_partial_drawing_registration_is_refused_at_the_hook():
    """the drawing four are ONE registration too, and the guard is the second one, not a widened first

    Both halves matter: part of the drawing group must refuse, and the layout group alone must NOT start
    refusing because the drawing group went unmentioned -- one all-or-nothing check across all nine
    arguments passes the first three assertions and fails the last.  Only refusals are exercised, since a
    call that got past the guard would replace this process's real registration with a stub.
    """
    from chython.core._core import _set_depict_fns

    def stub(*args, **kwargs):
        raise AssertionError('the stub must never be reachable: the call above has to refuse first')

    for offered in ('depict', 'scene', 'reaction_depict', 'reaction_scene'):
        with raises(ValueError, match='ONE registration'):
            _set_depict_fns(**{offered: stub})

    # all four 3-of-4 combinations: each omitted name must appear in the error message
    drawing_names = ('depict', 'scene', 'reaction_depict', 'reaction_scene')
    for omitted in drawing_names:
        with raises(ValueError, match=omitted):
            _set_depict_fns(**{n: stub for n in drawing_names if n != omitted})

    # the groups stay independent: the layout five on their own are a complete registration
    _set_depict_fns(clean2d=_molecule_layout.clean2d, layout2d=_molecule_layout.layout2d,
                    rescale2d=_molecule_layout.rescale2d,
                    reaction_clean2d=_reaction_layout.clean2d,
                    reaction_layout2d=_reaction_layout.layout2d)
    smiles('CCO').clean2d()
    assert smiles('CCO').depict().startswith('<svg'), 'the drawing group was collateral damage'


def test_the_container_methods_are_the_module_functions():
    """`mol.scene()` is `molecule_scene`, registered -- not a second assembly beside it"""
    mol = smiles('CCO')
    mol.clean2d()
    assert mol.scene().children == molecule_scene(mol).children
    assert mol.depict() == molecule_scene(mol).to_svg()

    rxn = ReactionContainer([smiles('CCO')], [smiles('CC=O')])
    rxn.clean2d()
    assert rxn.scene().children == reaction_scene(rxn).children


def test_depict_accepts_overlays_and_puts_them_under_the_structure():
    from chython.depict.overlay import Highlight

    mol = smiles('Oc1ccccc1')
    mol.clean2d()
    plain = mol.depict()
    with_highlight = mol.depict(overlays=[Highlight(atoms=[1])])
    assert len(with_highlight) > len(plain)
    # measured: a plain scene's children[0] is a stroked bond Path, a highlighted one's is the Group
    assert isinstance(mol.scene(overlays=[Highlight(atoms=[1])]).children[0], Group)
    assert isinstance(mol.scene().children[0], Path)


def test_a_field_overlay_produces_a_document_that_parses():
    from chython.depict.overlay import AtomField

    mol = smiles('Oc1ccccc1')
    mol.clean2d()
    values = {a.n: (-.45 if a.atomic_symbol == 'O' else .05) for a in mol.atoms()}
    fromstring(mol.depict(overlays=[AtomField(values)]))


def test_bond_scaling_changes_the_bond_widths_in_the_output():
    """BondScale does not draw; it must reach bond_paths, and this is the test that it does"""
    from chython.depict.overlay import BondScale

    mol = smiles('Oc1ccccc1')
    mol.clean2d()
    # two values, so the range maps end to end: measured, {(2, 3): .02, (3, 4): .12}.  ONE value gives
    # .07, the midpoint of a degenerate domain, and asserts nothing about the range
    scaled = mol.depict(overlays=[BondScale({(2, 3): 1., (3, 4): 2.}, width_range=(.02, .12))])
    assert 'stroke-width="0.12"' in scaled
    assert 'stroke-width="0.02"' in scaled
    assert 'stroke-width="0.04"' in mol.depict(), 'the default width, so the two above are the override'


def test_overlays_are_not_stored_on_the_molecule():
    from chython.depict.overlay import Highlight

    mol = smiles('Oc1ccccc1')
    mol.clean2d()
    highlighted = mol.depict(overlays=[Highlight(atoms=[1])])
    assert mol.depict() != highlighted
    assert mol._repr_svg_() == mol.depict()


def test_a_reaction_takes_overlays_indexed_by_position_in_molecules():
    from chython.depict.overlay import Highlight

    rxn = ReactionContainer([smiles('CCO')], [smiles('CC=O')])
    rxn.clean2d()
    reactant = list(rxn.molecules())[0]
    svg = rxn.depict(overlays={0: [Highlight(atoms=[a.n for a in reactant.atoms()][:1])]})
    fromstring(svg)
    assert len(svg) > len(rxn.depict())


def test_the_index_is_a_position_in_molecules_so_agents_are_reachable():
    """reactants -> agents -> products: an agent is index len(reactants), not something unaddressable"""
    from chython.depict.overlay import Highlight

    rxn = ReactionContainer([smiles('CCO')], [smiles('CC=O')], [smiles('O')])
    rxn.clean2d()
    agent = list(rxn.molecules())[1]
    assert agent is rxn.agents[0]
    svg = rxn.depict(overlays={1: [Highlight(atoms=[a.n for a in agent.atoms()])]})
    fromstring(svg)
    assert len(svg) > len(rxn.depict())


def test_two_identical_reactants_get_separate_overlays():
    """the reason the key is an index: as molecule objects these two would hash to one entry"""
    from chython.depict.overlay import Highlight

    rxn = ReactionContainer([smiles('CCO'), smiles('CCO')], [smiles('CCOCC')])
    rxn.clean2d()
    first = [a.n for a in list(rxn.molecules())[0].atoms()][:1]
    both = rxn.depict(overlays={0: [Highlight(atoms=first)], 1: [Highlight(atoms=first)]})
    one = rxn.depict(overlays={0: [Highlight(atoms=first)]})
    fromstring(both)
    assert len(both) > len(one) > len(rxn.depict())


def test_an_index_outside_the_reaction_is_refused():
    from chython.depict.overlay import Highlight

    rxn = ReactionContainer([smiles('CCO')], [smiles('CC=O')])
    rxn.clean2d()
    with raises(IndexError, match='2 molecules'):
        rxn.depict(overlays={5: [Highlight(atoms=[1])]})


def test_a_molecule_used_as_a_key_is_refused_by_type_not_silently_ignored():
    """the natural wrong guess -- say so, and say what to write instead"""
    from chython.depict.overlay import Highlight

    rxn = ReactionContainer([smiles('CCO')], [smiles('CC=O')])
    rxn.clean2d()
    with raises(TypeError, match='index into molecules'):
        rxn.depict(overlays={list(rxn.molecules())[0]: [Highlight(atoms=[1])]})


def test_a_field_brings_its_colorbar_and_a_highlight_does_not():
    from chython.depict.overlay import AtomField, Highlight

    mol = smiles('c1ccccc1O')
    mol.clean2d()
    charge = {a.n: -.4 if a.atomic_symbol == 'O' else .05 for a in mol.atoms()}
    with_bar = mol.depict(overlays=[AtomField(charge)])
    no_bar = mol.depict(overlays=[Highlight(atoms=[a.n for a in mol.atoms()][:2])])
    assert with_bar.count('<text') > no_bar.count('<text'), 'the ticks are the extra text'


def test_the_legend_can_be_switched_off_and_the_structure_is_unchanged_by_it():
    """the bar is placed OUTSIDE the content box, so turning it off must not move one atom"""
    from chython.depict.overlay import AtomField

    mol = smiles('c1ccccc1O')
    mol.clean2d()
    charge = {a.n: -.4 if a.atomic_symbol == 'O' else .05 for a in mol.atoms()}
    style = DepictStyle.preset('acs')
    with_bar = mol.scene(style=style, overlays=[AtomField(charge)])
    without = mol.scene(style=style.tuned(**{'page.legend': 'none'}), overlays=[AtomField(charge)])
    # one extra child at the end, and every structure node keeps its position
    assert len(with_bar.children) == len(without.children) + 1
    assert [n.bounds for n in with_bar.children[:-1]] == [n.bounds for n in without.children]


def _legend_swatches(scene):
    """The legend's swatch fills, low to high along the strip.

    The bar is one Group appended last, so it is always `children[-1]`; a swatch is a filled, unstroked
    Path of it.
    """
    from chython.depict.scene import Path

    bar = scene.children[-1]
    swatches = [n for n in bar.children if isinstance(n, Path) and n.fill is not None and n.stroke is None]
    return [s.fill for s in sorted(swatches, key=lambda s: (s.bounds.min_y, s.bounds.min_x))]


def _legend_rules(scene):
    """The legend's RULE strokes, low to high: a level with no filled band marks its one value with a line.

    `fill=False` makes every level one of these, matching a picture of stroked contours, which has no
    block of colour standing for an interval of values anywhere in it.
    """
    from chython.depict.scene import Path

    bar = scene.children[-1]
    rules = [n for n in bar.children if isinstance(n, Path) and n.fill is None and n.stroke is not None]
    return [s.stroke for s in sorted(rules, key=lambda s: (s.bounds.min_y, s.bounds.min_x))]


def _field_ink(scene):
    """Every colour the field group actually put on the page, fills and strokes alike."""
    from chython.depict.scene import Path

    bar = scene.children[-1]
    out = set()
    for node in scene.children:
        if node is bar or not isinstance(node, Group):
            continue
        for child in node.children:
            if isinstance(child, Path):
                out.add(child.fill or child.stroke)
    return out


def test_the_bar_advertises_exactly_the_bands_the_field_drew():
    """the bar is built from `bands_of` -- the levels this figure put ink on -- never from a recount

    With three explicit levels, falling through to `_default_levels`' 5 draws five swatches over three
    bands, and only the middle one is a colour the picture contains.
    """
    from chython.depict.overlay import AtomField

    mol = smiles('c1ccccc1O')  # phenol
    mol.clean2d()
    charge = {a.n: -.4 if a.atomic_symbol == 'O' else .05 for a in mol.atoms()}
    scene = mol.scene(overlays=[AtomField(charge, levels=[-.35, -.25, -.15], fill=False)])
    # `fill=False`, so the bar is three rules at their values, not three blocks of colour that would each
    # claim an interval nothing was filled over
    assert _legend_swatches(scene) == []
    assert _legend_rules(scene) == ['#4a61ce', '#6b88e3', '#8eabed']
    assert set(_legend_rules(scene)) == _field_ink(scene), 'and not one swatch more or fewer'


def test_the_default_nine_levels_draw_nine_bands_and_the_bar_labels_those_nine():
    """the default case: levels are of the SAMPLED RANGE, so `levels=9` means nine bands and nine swatches

    Taking them from the colormap's symmetrized domain (±0.6) instead leaves four of the nine reaching
    nothing, since this field's own range is −0.4388…+0.0340: five bands drawn, nine swatches advertised.
    """
    from chython.depict.overlay import _BAND_OUTLINE_DARKEN, AtomField, _darker

    mol = smiles('c1ccccc1O')  # phenol
    mol.clean2d()
    charge = {a.n: -.4 if a.atomic_symbol == 'O' else .05 for a in mol.atoms()}
    scene = mol.scene(overlays=[AtomField(charge)])
    swatches = _legend_swatches(scene)
    assert len(swatches) == 9, swatches
    # equality both ways, which `<=` alone would not give: every colour on the page is a swatch or that
    # swatch's own band outline, and every swatch is on the page
    assert _field_ink(scene) == set(swatches) | {_darker(c, _BAND_OUTLINE_DARKEN) for c in swatches}


def test_a_halo_only_figure_still_gets_a_bar_because_it_has_no_bands_to_match():
    """`_default_levels`' one honest case: a halo encodes its value as a RADIUS, so there are no bands to
    line the swatches up with and any count is as correct as any other."""
    from chython.depict.overlay import AtomHalo

    mol = smiles('c1ccccc1O')  # phenol
    mol.clean2d()
    charge = {a.n: -.4 if a.atomic_symbol == 'O' else .05 for a in mol.atoms()}
    scene = mol.scene(overlays=[AtomHalo(charge)])
    assert len(_legend_swatches(scene)) == 5


def test_the_colorbar_never_overlaps_the_structure():
    from chython.depict.colorbar import colorbar, legend_side, place_colorbar
    from chython.depict.overlay import AtomField
    from chython.depict.scene import Box

    mol = smiles('CC(=O)Oc1ccccc1C(=O)O')  # aspirin
    mol.clean2d()
    charge = {a.n: (-.5 if a.atomic_symbol == 'O' else .1) for a in mol.atoms()}
    overlays = (AtomField(charge),)
    scene = mol.scene(overlays=overlays)
    style = DepictStyle()

    def _is_legend(n):
        # the legend is one Group and z-order puts it last in the scene
        return isinstance(n, Group) and n is scene.children[-1]

    content = Box.of([n.bounds for n in scene.children if not _is_legend(n)])
    legend = Box.of([n.bounds for n in scene.children if _is_legend(n)])
    assert legend.min_x >= content.max_x or legend.max_y <= content.min_y


def _flatten(scene):
    out = []

    def walk(node):
        if isinstance(node, Group):
            for child in node.children:
                walk(child)
        else:
            out.append(node)

    for child in scene.children:
        walk(child)
    return out
