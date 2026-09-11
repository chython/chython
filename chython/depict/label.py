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
"""Atom labels: which atoms get one, what it says, and how much room its ink takes.

The decision and the measurement live together so they cannot disagree: a label reports the box it
actually occupies, measured through `metrics`, and `bonds.py` trims against that.  An unlabelled atom
still gets a `Label` -- `text=None` and a degenerate box at its point -- so bond trimming has no
`None` branch to miss.
"""
from math import atan2, cos, hypot, pi, sin
from typing import NamedTuple

from ..core import STEREO_ABS, STEREO_AND, STEREO_OR
from ._config import R_COLOUR, cpk as CPK
from .bonds import has_ink, segment_hits_box, trim_ink
from .metrics import text_box
from .scene import Box, Text, TextRun, to_hex
from .style import DepictStyle


__all__ = ['CPK', 'Label', 'element_colour', 'is_labelled', 'labels']


# `CPK` is the alias under which this module re-exports `_config.cpk`, the 118-entry palette indexed by
# Z - 1.  Imported, not copied: a second table would drift.

# 0.01 molecule units: the same threshold `has_layout` applies in the arena's fixed-point check, so
# `labels()` refuses exactly the planes the store calls "no layout".  Not a style field, or a caller
# could disagree with the container about whether a molecule has been laid out.
_DEGENERATE_SPAN = 0.01

# The typographic minus, not a hyphen: at 8 pt the two are visibly different lengths and heights.
MINUS = '−'

# CXSMILES' own vocabulary, so the picture and the `|&1:...|` the file carried spell it the same way.
_GROUP_MARKS = {STEREO_AND: '&', STEREO_OR: 'o'}

# Below this the annotation column's direction is vertical enough that a left/right anchor would hang the
# text off to one side of the sector it was placed in; it is centred on the direction instead.
_SIDE_LIMIT = .35

# What a crossed bond and an overlapped label cost when the candidate directions are compared.  A label
# is the heavier of the two because two sets of glyphs on one another are unreadable, while a number over
# a single line still reads -- badly, which is what the weight says.
_BOND_COST = 1.
_LABEL_COST = 2.

# Two atoms on one point: no direction is derivable from that neighbour, so it does not vote.  The same
# threshold `bonds.py` refuses a bond at, since it is the same question.
_MIN_LENGTH = 1e-9

# How far off its one bond a terminal atom's annotation is turned, and the corner a crowded one falls back
# to.  120 degrees is a drawing's own angle -- the place the next substituent would have been drawn -- and
# the bottom right is where a reader of a mapped structure looks first.  `_free_directions` explains both.
_TERMINAL_TURN = 2. * pi / 3.
_BOTTOM_RIGHT = (cos(-pi / 4.), sin(-pi / 4.))

# An annotation is never pushed further out than the atom's own ink to find room.  A number half a bond
# away from its atom has stopped saying WHICH atom it belongs to, which is a worse defect than the line it
# was moved off -- so the sector search chooses the emptiest side at that one distance and the crowded
# atom is answered by `figure.py`'s knock-out plate instead, `label.annotation_plate`.


def is_labelled(atom, style: DepictStyle) -> bool:
    """Does this atom get a written label?

    The skeletal convention: carbons are implied by the vertices and everything else is written.  A
    carbon carrying a fact the vertex cannot show -- a charge, a radical, an isotope, an unknown
    implicit-hydrogen count, no neighbours at all -- is written anyway.

    Annotations (CIP descriptors, map numbers) are placed separately by `labels()` and do not force a
    symbol here, so a chain carbon whose only extra fact is a stored CIP stays a bare vertex.
    """
    if atom.atomic_symbol != 'C':
        return True
    if style.atom.carbon:
        return True
    if atom.charge or (style.atom.radicals and atom.is_radical) or atom.isotope:
        return True
    # an unknown hydrogen count is such a fact, and carbon is where it arrives.  `_compose` writes the
    # `?` marker only for a labelled atom and only under `hydrogens`, so both flags are required here or
    # the symbol would be drawn for a fact still not shown.
    if style.atom.hydrogens and style.atom.unknown_h_marks and atom.implicit_h is None:
        return True
    if not atom.degree:                 # a lone atom: methane skeletally is an empty picture
        return True
    return False


def element_colour(atom, style: DepictStyle) -> str:
    """The fill colour for this atom's label, as a normalized lowercase hex string."""
    if not style.atom.colour_by_element:
        return style.atom.default_colour
    if atom.atomic_symbol == 'C':
        return style.atom.carbon_colour
    if atom.element == 0:
        return to_hex(R_COLOUR)
    return to_hex(CPK[atom.element - 1])


class Label(NamedTuple):
    """One atom's label.  `text` is None when the atom is drawn as a bare vertex.

    `box` is the ink a bond must stop clear of -- measured over ALL runs of the composite, then
    inflated by `LabelStyle.pad` -- and for an unlabelled atom it is a degenerate box at `anchor`.
    `anchor` is the atom's point, which is where a radical dot, a halo or a highlight ribbon centres.

    `annotations` carries the stereo statement and the map number, in that order, as fully positioned
    `Text` objects, set into the emptiest sector around the atom by `_place_annotations` -- BESIDE the
    atom in every case, close enough that which atom they belong to is not in question.  They are not
    included in `box`: trimming bonds for them would open a gap at both ends of every bond in a mapped
    figure, so they are moved out of the bonds' way instead of the bonds out of theirs, and where no side
    is clear `figure.py` puts a plate under them rather than moving them away.
    """
    atom: int
    text: Text | None
    box: Box
    anchor: tuple[float, float]
    annotations: tuple[Text, ...] = ()


def labels(mol, plane, style: DepictStyle) -> dict[int, Label]:
    """One `Label` per atom, in arena order.

    `plane` is passed rather than read off the molecule so a caller can draw a layout it has not stored.
    The guard therefore validates `plane` itself -- every atom present, non-degenerate span -- and not
    `mol.has_layout`, which would refuse that very case.

    :raises ValueError: an atom is missing from `plane`, or every atom sits near one point.
    """
    for sid in mol:
        if sid not in plane:
            raise ValueError('plane has no entry for atom %d: call clean2d() or pass a proper '
                             '2D layout plane with a key for every atom' % sid)
    ids = list(mol)
    if len(ids) > 1:
        xs = [plane[sid][0] for sid in ids]
        ys = [plane[sid][1] for sid in ids]
        if max(xs) - min(xs) < _DEGENERATE_SPAN and max(ys) - min(ys) < _DEGENERATE_SPAN:
            raise ValueError('plane is degenerate (all atoms near one point): call clean2d() first')
    out = {}
    for atom in mol.atoms():
        sid = atom.n
        point = plane[sid]
        baseline = point[1] - style.label.baseline_shift * style.label.size
        if is_labelled(atom, style):
            runs, anchor_mode = _compose(atom, mol, plane, style)
            text = Text(runs, x=point[0], y=baseline, anchor=anchor_mode,
                        fill=element_colour(atom, style))
            box = text_box(text).inflate(style.label.pad)
        else:
            text = None
            box = Box(point[0], point[1], point[0], point[1])
        out[sid] = Label(sid, text, box, point)
    # A second pass, because where an annotation goes is decided against the OTHER labels' ink and the
    # trimmed bonds -- neither of which is known while the first atom is still being measured.
    if style.atom.stereo_labels or style.atom.map_numbers or style.atom.stereo_groups:
        _place_annotations(mol, plane, out, style)
    return out


def _compose(atom, mol, plane, style):
    """The runs of one label, in reading order, and the anchor mode that centres them.

    Order: isotope, symbol, hydrogens, charge -- or hydrogens, isotope, symbol, charge when the label
    reads right-to-left.  The isotope immediately precedes the symbol in both directions, so a flipped
    ¹³C reads ¹³C and not ¹³H₃C.  The direction depends on where the neighbours are: `NH2` at the left
    end of a chain has to read `H2N`, or the H sits on top of the bond.

    The stereo statement and the map number are separate `Text` annotations from `_place_annotations`, not
    runs here, so they keep their own colour, stay out of the bond-trim box, and are free to sit somewhere
    else entirely when the label's own side is occupied.
    """
    label = style.label
    size = label.size
    flip = _reads_right_to_left(atom, mol, plane)
    symbol = TextRun(atom.atomic_symbol, family=label.family, size=size)

    hydrogens = []
    count = atom.implicit_h
    if style.atom.hydrogens:
        if count is None:
            if style.atom.unknown_h_marks:
                hydrogens.append(TextRun('?', family=label.family, size=size * label.superscript_scale,
                                         dy=size * label.superscript_rise))
        elif count:
            hydrogens.append(TextRun('H', family=label.family, size=size))
            if count > 1:               # H1 is not written; nobody writes it
                hydrogens.append(TextRun(str(count), family=label.family,
                                         size=size * label.subscript_scale,
                                         dy=-size * label.subscript_drop))

    head = []
    if style.atom.isotopes and atom.isotope:   # isotope is 0 for "none", not None
        head.append(TextRun(str(atom.isotope), family=label.family,
                            size=size * label.superscript_scale, dy=size * label.superscript_rise))

    tail = []
    if style.atom.charges and atom.charge:
        magnitude = abs(atom.charge)
        sign = '+' if atom.charge > 0 else MINUS
        tail.append(TextRun(('' if magnitude == 1 else str(magnitude)) + sign, family=label.family,
                            size=size * label.superscript_scale, dy=size * label.superscript_rise))

    if flip:
        return tuple(hydrogens + head + [symbol] + tail), 'middle'
    return tuple(head + [symbol] + hydrogens + tail), 'middle'


def _place_annotations(mol, plane, out: dict[int, Label], style: DepictStyle) -> None:
    """Fill in `Label.annotations` for every atom, in place, INTO THE ROOM THERE IS.

    A map number is on every atom of a mapped record, so the bond cannot be trimmed for it the way it is
    for a symbol -- that would open a gap at both ends of every bond in the figure.  The annotation picks
    its SIDE instead: the widest sector between the atom's own bonds, where by construction there is no
    bond ink, with the candidate sectors MEASURED against the trimmed bonds and the other atoms' ink so
    the widest can lose to a clear narrower one.  Number over line was the complaint.

    WHAT IT WILL NOT DO IS WALK AWAY.  Every candidate sits at the same distance -- where the atom's own
    ink ends -- so the choice is which side, never how far.  At a fused branch point every sector has a
    bond in it and the least-bad one is chosen; the number stays where it belongs and `figure.py` draws
    the knock-out plate that makes it readable.

    TWO FIXED ROWS within the chosen sector: the stereo statement `label.annotation_rise` above the
    sector's line, the map number `annotation_drop` below it.  Fixed and not chosen from what is present,
    or one compound would be drawn two ways depending on whether a descriptor happened to be stored -- and
    each row's distance is its own, `_annotation_column`, so an enhanced-stereo `(R)&1` beside the number
    does not push the number out too.  The side can still differ: the descriptor is part of the column the
    sectors are scored with, and a wider column may fit a different one.  Both rows are annotations, so the
    stereo statement is placed by this, is kept as close as the number is, and is plated the same way.

    Neither row is added to `Label.box`: the boxes here are what bonds are kept OFF, not what trims them.
    """
    label = style.label
    segments = _bond_segments(mol, plane, out, style)
    boxes = [lb.box for lb in out.values() if lb.text is not None]
    for atom in mol.atoms():
        rows = []
        stereo = _stereo_runs(atom, style)
        if stereo:
            rows.append((stereo, style.atom.default_colour, label.annotation_rise * label.size,
                         label.size * label.stereo_scale))
        if style.atom.map_numbers and atom.map_number:
            # THE NUMBER AND NOTHING ELSE, and no leading colon: the colon is SMILES punctuation --
            # `[CH3:1]` needs it to separate the label from the atom -- and a drawing has already
            # separated them by putting the number beside the atom, so it reads as a stray mark.
            size = label.size * label.map_scale
            rows.append(([TextRun(str(atom.map_number), family=label.family, size=size)],
                         label.map_colour, -label.annotation_drop * label.size, size))
        if not rows:
            continue

        current = out[atom.n]
        best = None
        for direction in _free_directions(atom, mol, plane):
            column = _annotation_column(rows, current, direction, style)
            cost = _crowding(column, segments, boxes)
            if best is None or cost < best[0]:
                best = (cost, column)
            if not cost:                    # a clear sector, and the widest is tried first
                break
        out[atom.n] = current._replace(annotations=best[1])
        boxes.extend(text_box(text) for text in best[1])   # the next atom's number keeps off this one


def _bond_segments(mol, plane, out: dict[int, Label], style: DepictStyle) -> list[tuple]:
    """The bond axes as DRAWN -- trimmed at whichever ends carry ink, the same call `bond_paths` makes --
    each with the half-width of the shape that will be drawn on it.

    The drawn line and not the full centre-to-centre span, or the room a label's knock-out just opened
    would be scored as occupied and the annotation pushed out of the one place it fits.

    Only the axis is known here: a double bond's second line and a triple's outer pair are `bond_paths`'
    geometry, drawn later off a kekule form this module has no business computing.  So the axis carries
    the furthest THAT ORDER's lines stray from it, and a plain single bond stays the thin thing it is --
    one margin for every bond would be `bond.spacing` wide everywhere and would find no room anywhere.
    """
    bond = style.bond
    half = bond.width / 2.
    strays = {2: bond.spacing + half, 3: bond.triple_spacing + half, 4: bond.spacing + half}
    segments = []
    for b in mol.bonds():
        segment = trim_ink(plane[b.n], plane[b.m], out[b.n].box, out[b.m].box, bond.trim)
        if segment is not None:
            segments.append((segment[0], segment[1], strays.get(b.order, half)))
    return segments


def _free_directions(atom, mol, plane) -> list[tuple[float, float]]:
    """Where this atom's annotations may go, best first, as unit vectors.  One rule per degree.

    ONE BOND: 120 degrees off it, either way, and only then straight away.  Straight away is the bisector
    of the one sector and the widest room there is, but it is also the atom's own bond line continued, and a
    number on that line reads as the chain going on -- and at a terminal atom it is exactly where the
    hydrogens are written.  120 degrees is where a drawing would have put the next substituent, so it is
    where a reader looks for something belonging to this atom.

    TWO BONDS: the bisector of the widest sector, then the other, which is the usual pair of choices at a
    chain vertex -- outside the elbow first, inside it second.

    THREE OR MORE: still widest sector first, and then the bottom right, which is the corner a mapped
    drawing conventionally puts a number in.  At a fused or bridged centre every sector is narrow and every
    one of them holds something; the fixed corner is a candidate rather than the answer, so it wins only by
    scoring better than the sectors, and `figure.py` plates whichever one wins.

    A lone atom gets the one direction a label reads in.
    """
    x, y = plane[atom.n]
    angles = []
    for neighbour in mol.neighbors_of(atom.n):
        nx, ny = plane[neighbour]
        if hypot(nx - x, ny - y) > _MIN_LENGTH:
            angles.append(atan2(ny - y, nx - x))
    if not angles:
        return [(1., 0.)]
    if len(angles) == 1:
        # Sorted so the drawing is not decided by which way the bond happens to point: of the two 120s the
        # lower one is offered first, and the right one of two equally low, which is the corner a reader of
        # mapped structures is used to.
        turns = sorted((angles[0] + _TERMINAL_TURN, angles[0] - _TERMINAL_TURN),
                       key=lambda a: (round(sin(a), 9), -round(cos(a), 9)))
        return [(cos(a), sin(a)) for a in (*turns, angles[0] + pi)]
    angles.sort()
    sectors = []
    for i, start in enumerate(angles):
        span = (angles[(i + 1) % len(angles)] - start) % (2. * pi)
        sectors.append((span, start + span / 2.))
    sectors.sort(key=lambda sector: -sector[0])
    out = [(cos(middle), sin(middle)) for _, middle in sectors]
    if len(angles) > 2:
        out.append(_BOTTOM_RIGHT)
    return out


def _annotation_column(rows, label: Label, direction, style: DepictStyle) -> tuple[Text, ...]:
    """The rows as positioned `Text`s, set on the atom's own point and then slid along `direction` until
    their INK is clear of the atom's -- that far and NO FURTHER, whatever else is in the way.

    Slid rather than offset by a distance, because a distance would have to be the label's reach in this
    direction plus each row's own rise or drop plus the baseline shift, and the three do not add up to
    anything a reader can check.  What the reader can check is the measured statement: the number's ink is
    outside the box the bonds already stop at.  A bare vertex has no ink to be outside of, so it spends
    `label.pad` -- the same gap a labelled atom's inflated box already carries.

    PER ROW, and not by the union of the column's ink.  `(R)&1` is three times the width of `12`, so a
    union slide is the descriptor's slide, and the map number would stand off further on the centres that
    happen to carry a descriptor than on the ones that do not -- one series drawn two ways, which is the
    defect the fixed rows exist to prevent.  So each row is slid by what its OWN ink asks for.

    Where that stacks two rows on each other -- the sector pointing along the axis the rise and the drop
    already separate them on -- the LAST row keeps its place and the earlier ones go outside it.  Last is
    the map number, and it keeps its place because it is the row on every atom of a mapped record: a number
    at its own distance whatever else the centre carries, with the descriptor, which qualifies it and is
    wider anyway, stacked beyond.  Nothing is ever pushed out to look for room, only to stay off ink.
    """
    dx, dy = direction
    x, y = label.anchor
    if dx > _SIDE_LIMIT:
        anchor = 'start'
    elif dx < -_SIDE_LIMIT:
        anchor = 'end'
    else:
        anchor = 'middle'
    column = tuple(Text(runs, x=x, y=y + offset - style.label.baseline_shift * size, anchor=anchor,
                        fill=fill)
                   for runs, fill, offset, size in rows)
    keep_off = label.box if has_ink(label.box) else label.box.inflate(style.label.pad)
    slides = [_clearing_slide(text_box(text), keep_off, direction) for text in column]

    # Innermost LAST-ROW-FIRST: a row only ever has to clear the rows already placed, and the one that must
    # not be moved by its neighbours goes down first.
    placed: dict[int, Text] = {}
    for i in reversed(range(len(column))):
        ink = text_box(column[i])
        distance = slides[i]
        for done in placed.values():
            # Only a row this one would really land on, and `label.pad` past that one: ink boxes are tight,
            # so a slide that stopped at contact would set `a` against `10` with nothing between them.  The
            # rise and the drop separate the rows wherever the sector does not point along them, and there
            # a row that is already clear must not be pushed sideways past a row it never touched.
            moved = ink.translated(dx * distance, dy * distance)
            keep_apart = text_box(done).inflate(style.label.pad)
            if _boxes_overlap(moved, keep_apart):
                distance += _clearing_slide(moved, keep_apart, direction)
        placed[i] = column[i].translated(dx * distance, dy * distance)
    return tuple(placed[i] for i in range(len(column)))


def _clearing_slide(box: Box, obstacle: Box, direction) -> float:
    """How far along `direction` `box` must slide to stop overlapping `obstacle`.  0 if it already does.

    Either axis separating them is enough, so the answer is the SMALLER of the two axes' demands -- the
    x-only case being the one the drawing has always done: a label's annotation set just past its box.

    An overlap is NOT tested for: a box already clear of the obstacle in the direction of travel gets a
    negative demand and is clamped to zero, and one clear only across that direction is still slid, which
    is the "beside" the drawing wants -- an annotation that merely clears the atom's box in y belongs past
    it in x, not under it on the bonds converging at the vertex.  Row against row is the other case and
    asks the caller to check first: two rows the rise and the drop have already separated must not move.
    """
    demands = []
    for delta, near, far, other_near, other_far in ((direction[0], box.min_x, box.max_x,
                                                     obstacle.min_x, obstacle.max_x),
                                                    (direction[1], box.min_y, box.max_y,
                                                     obstacle.min_y, obstacle.max_y)):
        if delta > _MIN_LENGTH:
            demands.append((other_far - near) / delta)
        elif delta < -_MIN_LENGTH:
            demands.append((other_near - far) / delta)
    if not demands:                      # a zero direction, which `_free_directions` never returns
        return 0.
    return max(min(demands), 0.)


def _boxes_overlap(box: Box, other: Box) -> bool:
    """Do two boxes share any area?  Touching counts, which is what the placement wants: ink against ink."""
    return (box.min_x <= other.max_x and other.min_x <= box.max_x
            and box.min_y <= other.max_y and other.min_y <= box.max_y)


def _crowding(column: tuple[Text, ...], segments, boxes) -> float:
    """What this column lands on: bond lines crossed, plus other atoms' and annotations' ink overlapped.

    A count and not a distance, because the answer wanted is "does it read", and one crossing already
    means no.  Own box excluded by construction -- `_annotation_column` starts past it.
    """
    cost = 0.
    for text in column:
        box = text_box(text)
        for p, q, stray in segments:
            if segment_hits_box(p, q, box.inflate(stray)):
                cost += _BOND_COST
        for other in boxes:
            if (other.min_x < box.max_x and other.max_x > box.min_x
                    and other.min_y < box.max_y and other.max_y > box.min_y):
                cost += _LABEL_COST
    return cost


def _stereo_runs(atom, style) -> list:
    """The stereo row's runs: the stored CIP descriptor, then the enhanced-stereo group mark.

    One `Text` and not two, because `(R)&1` is one statement about one centre and reads as one line.  The
    group id is upright: it is a label, not a descriptor.
    """
    label = style.label
    size = label.size * label.stereo_scale
    runs = []
    if style.atom.stereo_labels and atom.cip is not None:
        runs.append(TextRun(f'({atom.cip})', family=label.family, size=size,
                            style='italic' if label.stereo_italic else 'normal'))
    if style.atom.stereo_groups:
        mark = _stereo_group_mark(atom)
        if mark is not None:
            runs.append(TextRun(mark, family=label.family, size=size))
    return runs


def _stereo_group_mark(atom) -> str | None:
    """`&N` for AND, `oN` for OR, `a` for ABS, or None when the atom is in no collection.

    ABS carries group 0, so it draws bare.  No `has_stereo_groups` gate: an unspecified atom answers
    `(STEREO_UNSPECIFIED, 0)` through the same O(1) lookup, and a bare `@` sets no collection at all.
    """
    kind, group = atom.stereo_group
    if kind == STEREO_ABS:
        return 'a'
    mark = _GROUP_MARKS.get(kind)
    return None if mark is None else f'{mark}{group}'


def _reads_right_to_left(atom, mol, plane) -> bool:
    """True when the hydrogens belong on the left of the symbol.

    Only when every neighbour is strictly to the right, leaving the left side clear: one neighbour at or
    to the left already occupies it.  A degree-0 atom reads left to right by convention.
    """
    x = plane[atom.n][0]
    has_neighbours = False
    for neighbour in mol.neighbors_of(atom.n):
        has_neighbours = True
        if plane[neighbour][0] <= x:
            return False   # at least one neighbour is at or left of this atom: do not flip
    return has_neighbours  # True only when every neighbour is strictly to the right
