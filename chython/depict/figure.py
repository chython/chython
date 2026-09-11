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
"""Assembly: labels, bonds, wedges and a reaction's furniture composed into one `Scene`.

Painter's order, global over a whole figure and not per molecule: overlays, bonds, wedges (which claim
their bond keys so `bonds.py` skips them), radical dots, type.  Drawing stores nothing -- a molecule
with no layout gets a temporary one, reported through `log`; `clean2d()` keeps a layout.  No caching.
"""

from math import sqrt

from ..core import LogRecord
from .bonds import bond_paths, has_ink
from .colorbar import colorbar, legend_side, place_colorbar
from .field import contour_levels
from .label import element_colour, labels
from .layout import molecule as _molecule_layout, reaction as _reaction_layout
from .overlay import BondScale, bands_of, render_overlays, scale_of, tiled_swatches
from .scene import (Box, Path, Scene, Text, TextRun, WHITE, circle, ellipse, polyline, rounded_box)
from .style import DepictStyle, get_depict_style
from .wedge import wedge_paths


__all__ = ['molecule_depict', 'molecule_scene', 'reaction_depict', 'reaction_scene']


# What an ellipse's radii are scaled by to pass through the corners of the box they were measured from.
SQRT2 = sqrt(2.)


def molecule_scene(mol, *, style: DepictStyle | None = None, plane=None, overlays=(),
                   log=None) -> Scene:
    """One molecule as a `Scene`, in molecule coordinates, y-up.

    :param plane: draw a layout the molecule does not carry.  Given `None`, its own layout is used when
        it has one, else a temporary one is computed, not stored, and reported through `log`.
    :param overlays: a sequence of `Highlight`, `AtomHalo`, `AtomField`, `BondScale` or `ValueLabels`,
        rendered in order.  A colorbar is added when any of them carries a colour scale.
    """
    if style is None:
        style = get_depict_style()
    if plane is None:
        plane = _plane_for(mol, log)
    under, over = _molecule_nodes(mol, plane, style, log, overlays)

    # The colorbar comes after the content box is known: its position is derived from the content bounds.
    legend_nodes = []
    children = under + over
    if overlays:
        content = Box.of(n.bounds for n in children)
        side = legend_side(style, overlays, content)
        if side is not None:
            cmap = next(s for s in (scale_of(o) for o in overlays) if s is not None)
            # The bands THIS figure drew, from the same helper `AtomField.render` draws from -- never a
            # level count re-derived beside them.
            swatches = [s for o in overlays for s in bands_of(o, mol, plane, style)]
            if not swatches:
                # No bands to line up with, so the bar is a plain scale whose swatches tile the domain
                # edge to edge: there are no missing levels to leave a gap for.
                swatches = tiled_swatches(cmap, contour_levels(_default_levels(overlays),
                                                               cmap.vmin, cmap.vmax))
            nodes, box = colorbar(cmap, swatches=swatches, style=style,
                                  side=side,
                                  length=content.height if side == 'right' else content.width)
            legend_nodes = [place_colorbar(nodes, box, content, style, side)]

    return Scene(children + legend_nodes)


def molecule_depict(mol, *, style: DepictStyle | None = None, plane=None, overlays=(),
                    log=None) -> str:
    """The SVG document for one molecule.  `mol.depict()`'s body."""
    if style is None:
        style = get_depict_style()
    return molecule_scene(mol, style=style, plane=plane, overlays=overlays, log=log).to_svg(
        style=style)


def reaction_scene(rxn, *, style: DepictStyle | None = None, overlays=None, log=None) -> Scene:
    """A whole reaction as one `Scene`: its molecules arranged left to right, the arrow, the `+` signs.

    The arrangement comes from `layout/reaction.py`, which shifts plane dicts and never touches the
    molecules, so the members' coordinates do not change.  The arrow and the signs belong to this
    drawing and are recomputed per figure.

    :param overlays: `{index: [overlay, ...]}`, keyed by position in `molecules()`
        (reactants->agents->products).
    """
    if style is None:
        style = get_depict_style()
    planes, arrow, signs = _reaction_layout.layout2d(rxn)
    molecules = list(rxn.molecules())
    per_molecule = _reaction_overlays(molecules, overlays) if overlays else [() for _ in molecules]

    under: list[Path] = []
    over: list[Text] = []
    for mol, plane, mol_overlays in zip(molecules, planes, per_molecule):
        molecule_under, molecule_over = _molecule_nodes(mol, plane, style, log, mol_overlays)
        under.extend(molecule_under)
        over.extend(molecule_over)
    under.extend(_arrow_paths(arrow, style))
    over.extend(_sign_texts(signs, style))
    return Scene(under + over)


def reaction_depict(rxn, *, style: DepictStyle | None = None, overlays=None, log=None) -> str:
    """The SVG document for a whole reaction.  `rxn.depict()`'s body."""
    if style is None:
        style = get_depict_style()
    return reaction_scene(rxn, style=style, overlays=overlays, log=log).to_svg(style=style)


def _default_levels(overlays):
    """Swatch count for a bar with no bands to match -- a halo-only or bond-only figure.

    Any count is as correct as any other there; five reads as an axis without crowding the ticks.
    """
    return 5


def _reaction_overlays(molecules, overlays):
    """Resolve `{index: [overlay, ...]}` against `molecules()`, refusing before anything is drawn.

    The key is a position, not a molecule: a molecule key would call `__hash__` -- the canonical form --
    on every component of every reaction depicted, so it raises `TypeError` rather than silently
    drawing a plain picture.
    """
    per_molecule = [() for _ in molecules]
    for key, group in overlays.items():
        if not isinstance(key, int) or isinstance(key, bool):
            raise TypeError(f'reaction overlays are keyed by index into molecules(), got {key!r}')
        if not -len(molecules) <= key < len(molecules):
            raise IndexError(f'overlay index {key} but this reaction has {len(molecules)} molecules')
        per_molecule[key] = tuple(group)
    return per_molecule


def _plane_for(mol, log) -> dict[int, tuple[float, float]]:
    """The layout to draw: the stored one, or a temporary that is reported and thrown away.

    `layout2d` never stores, so the branch here is only about filing the log line -- which happens only
    when a layout was really computed, since a line on every picture would be noise.
    """
    if mol.has_layout:
        return _molecule_layout.layout2d(mol)
    plane = _molecule_layout.layout2d(mol)
    if log is not None:
        log.append(LogRecord('depict:layout', tuple(mol),
                             'this molecule carried no 2D layout, so one was computed for this drawing '
                             'only and NOT stored; call clean2d() for a layout that is kept'))
    return plane


def _molecule_nodes(mol, plane, style: DepictStyle, log, overlays=()) -> tuple[list, list]:
    """One molecule's nodes, split into what goes under the type and what is the type.

    Two lists rather than one, so a caller composing several molecules keeps painter's order over the
    whole figure -- see the module docstring.  Overlays bracket the structure: field bands, highlight
    regions and halos go under the bonds, value labels over the atom labels.
    """
    boxes = labels(mol, plane, style)

    # Overlays go under the structure: a highlight over a bond hides the chemistry it points at.
    overlay_under, overlay_over = render_overlays(overlays, mol, plane, boxes, style, log=log)

    # BondScale overrides go on the bond's own path, so bond_paths applies them atomically -- two
    # stacked strokes at different widths would show as an outline.
    widths: dict[tuple[int, int], float] = {}
    colours: dict[tuple[int, int], str] = {}
    for o in overlays:
        if isinstance(o, BondScale):
            if o.encode in ('width', 'both'):
                widths.update(o.bond_widths(mol, style))
            if o.encode in ('color', 'both'):
                colours.update(o.bond_colours(mol, style))

    # The stereo bonds claim their keys, and the plain bonds skip them.
    wedges, claimed = wedge_paths(mol, plane, boxes, style, log=log)
    bond_nodes = bond_paths(mol, plane, boxes, style, skip=claimed,
                            widths=widths if widths else None,
                            colours=colours if colours else None,
                            log=log)

    radical_nodes = _radical_dots(mol, boxes, style)

    label_nodes: list[Text] = []
    annotations: list[Text] = []
    for label in boxes.values():
        if label.text is not None:
            label_nodes.append(label.text)
        annotations.extend(label.annotations)

    # Every plate before every glyph, not each plate before its own: a plate that reaches under a
    # neighbouring symbol must not knock that symbol out, and the whole type layer is one painter's step.
    plates = _annotation_plates(annotations, style)

    # The z-order is fixed -- no node carries a z and nothing is sorted.  Wedges come AFTER bonds
    # because a wedge's wide base overlaps the adjacent bonds at the shared vertex and must win that
    # overlap; `skip=claimed` only covers the wedge's own axis.
    under = [*overlay_under, *bond_nodes, *wedges, *radical_nodes]
    over = [*plates, *label_nodes, *annotations, *overlay_over]
    return under, over


def _annotation_plates(annotations, style: DepictStyle) -> list[Path]:
    """A knock-out plate under EVERY annotation, `label.annotation_plate_pad` around its measured ink.

    The number stays beside its atom -- `label.py` will not walk it away from the atom it names -- so at a
    fused branch point, where every sector has a bond in it, something has to give between the number and
    the line.  It is the line: a plate in the background colour, drawn over the structure and under the type.

    UNCONDITIONAL, and not per number measured against the drawn paths.  A plate over nothing is invisible,
    a plate the line merely grazes is the one a reader needs, and "does the ink touch" is not the question a
    reader asks: a digit in the corridor between a ring's perimeter and its inner line touches neither and
    is still unreadable.  `label.annotation_plate = 'none'` withholds the whole layer.
    """
    label = style.label
    if label.annotation_plate == 'none' or not annotations:
        return []
    # The plate is the colour of what is behind it, and on an unpainted page that is the white the
    # figure is going to be put on.  Stated in the style when that guess is wrong.
    colour = label.annotation_plate_colour or style.page.background or WHITE
    out = []
    for text in annotations:
        box = text.bounds.inflate(label.annotation_plate_pad)
        if label.annotation_plate == 'ellipse':
            # An ellipse THROUGH the box's corners, or the digits would hang out of the sides of it.
            plate = ellipse((box.min_x + box.max_x) / 2., (box.min_y + box.max_y) / 2.,
                            box.width / 2. * SQRT2, box.height / 2. * SQRT2)
        else:
            # Half the shorter side: a one-digit plate is a disc, a three-digit one a stadium, and
            # neither has a corner pointing at the line it was drawn to cover.
            plate = rounded_box(box, min(box.width, box.height) / 2.)
        out.append(Path([plate], fill=colour))
    return out


def _radical_dots(mol, boxes, style: DepictStyle) -> list[Path]:
    """A filled dot per radical atom, centred on the atom's x and clear of its label's ink.

    Not at the anchor: `LabelStyle.baseline_shift` centres a label's cap height on the atom point, so a
    dot there would sit inside the glyph.  It goes `AtomStyle.radical_gap` above the measured ink box.
    """
    if not style.atom.radicals:
        return []
    out = []
    for atom in mol.atoms():
        if not atom.is_radical:
            continue
        label = boxes[atom.n]
        x = label.anchor[0]
        top = label.box.max_y if has_ink(label.box) else label.anchor[1]
        out.append(Path([circle(x, top + style.atom.radical_gap, style.atom.radical_radius)],
                        fill=element_colour(atom, style)))
    return out


def _arrow_paths(arrow, style: DepictStyle) -> list[Path]:
    """The reaction arrow: a stroked shaft and a filled head, from the span the layout returned.

    `(x1, x2, y)` is the whole span and the head is inside it, so a head never hangs past `x2` into the
    clearance before the first product.  The head is filled so its size does not track `bond.width`.
    """
    x1, x2, y = arrow
    reaction = style.reaction
    base = x2 - reaction.head_length
    half = reaction.head_width / 2.
    return [Path([polyline([(x1, y), (base, y)])], stroke=reaction.colour,
                 width=reaction.arrow_width, cap='butt', join='miter'),
            Path([polyline([(x2, y), (base, y + half), (base, y - half)], closed=True)],
                 fill=reaction.colour)]


def _sign_texts(signs, style: DepictStyle) -> list[Text]:
    """A `+` per gap between two members of one side, at `ReactionStyle.sign_size`.

    The label family's typographic plus, so it matches the charges beside it, and vertically centred by
    the same `baseline_shift` fraction a label uses.
    """
    reaction = style.reaction
    return [Text([TextRun('+', family=style.label.family, size=reaction.sign_size)],
                 x=x, y=y - style.label.baseline_shift * reaction.sign_size, anchor='middle',
                 fill=reaction.colour)
            for x, y in signs]
