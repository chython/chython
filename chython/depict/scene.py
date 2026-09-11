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
"""The scene IR -- three primitives every backend can draw, in MOLECULE coordinates, y-up.

`Path`, `Text`, `Group`, and no chemistry: a `Path` does not know it is a bond.  The y-flip and the unit
scale are device facts and happen once, in each backend's emitter.  Geometry is absolute (no transform
stack) and everything is `frozen=True, slots=True`: a `Scene` is shared between backends and caches.
"""
from collections.abc import Sequence
from dataclasses import dataclass, replace
from math import inf


__all__ = ['Box', 'EMPTY_BOX', 'Group', 'Path', 'Scene', 'Text', 'TextRun', 'BLACK', 'WHITE',
           'circle', 'close', 'curve', 'ellipse', 'line', 'move', 'polyline', 'rounded_box', 'rgb',
           'to_hex']


# A cubic Bezier approximates a quarter circle to within .00027 r at this handle length.  Every round
# thing in the package is built from `circle()`, so none of them re-derives it.
KAPPA = 0.5522847498307936

# Two names, and not the start of a CSS colour list: `scene.py` needs a default stroke and a default
# knock-out fill without importing the style tree it sits below.  Everything else is hex from `style.py`.
_NAMED = {'black': '#000000', 'white': '#ffffff'}

BLACK = '#000000'
WHITE = '#ffffff'


def rgb(r: int, g: int, b: int) -> str:
    """`(255, 0, 128)` -> `'#ff0080'`.  Channels are 0-255 and out-of-range is an error, not a clamp."""
    for channel in (r, g, b):
        if not 0 <= channel <= 255:
            raise ValueError(f'colour channel out of range: {(r, g, b)}')
    return '#%02x%02x%02x' % (r, g, b)


def to_hex(colour: str | None) -> str | None:
    """Normalize a colour to lowercase `'#rrggbb'`.  None passes through -- it means "do not paint"."""
    if colour is None:
        return None
    if not isinstance(colour, str):
        raise ValueError(f'colour must be a string or None, got {colour!r}')
    lowered = colour.lower()
    if lowered in _NAMED:
        return _NAMED[lowered]
    if len(lowered) == 7 and lowered[0] == '#':
        try:
            int(lowered[1:], 16)
        except ValueError:
            pass
        else:
            return lowered
    if len(lowered) == 4 and lowered[0] == '#':  # #abc -> #aabbcc
        try:
            int(lowered[1:], 16)
        except ValueError:
            pass
        else:
            return '#' + lowered[1] * 2 + lowered[2] * 2 + lowered[3] * 2
    raise ValueError(f'unparsable colour {colour!r}: expected #rrggbb, #rgb, "black" or "white"')


@dataclass(frozen=True, slots=True)
class Box:
    """An axis-aligned bounding box in molecule coordinates.  y-up, so `max_y` is the top."""
    min_x: float
    min_y: float
    max_x: float
    max_y: float

    @property
    def width(self) -> float:
        return self.max_x - self.min_x

    @property
    def height(self) -> float:
        return self.max_y - self.min_y

    def inflate(self, distance: float) -> 'Box':
        return Box(self.min_x - distance, self.min_y - distance,
                   self.max_x + distance, self.max_y + distance)

    def translated(self, dx: float, dy: float) -> 'Box':
        return Box(self.min_x + dx, self.min_y + dy, self.max_x + dx, self.max_y + dy)

    def contains(self, other: 'Box') -> bool:
        """Does this box cover `other` entirely?  The empty box is contained by anything."""
        if other.min_x > other.max_x:
            return True
        return (self.min_x <= other.min_x and self.min_y <= other.min_y
                and self.max_x >= other.max_x and self.max_y >= other.max_y)

    @classmethod
    def of(cls, boxes) -> 'Box':
        """The union of an iterable of boxes, starting from `EMPTY_BOX`.

        Seeding from the INVERTED empty box rather than `Box(0, 0, 0, 0)` is what keeps the origin out of
        every bounding box.
        """
        result = EMPTY_BOX
        for box in boxes:
            result = result.union(box)
        return result

    def union(self, other: 'Box') -> 'Box':
        if other.min_x > other.max_x:   # the empty box, which must not drag a union to infinity
            return self
        if self.min_x > self.max_x:
            return other
        return Box(min(self.min_x, other.min_x), min(self.min_y, other.min_y),
                   max(self.max_x, other.max_x), max(self.max_y, other.max_y))

    def __iter__(self):
        """so `approx((a, b, c, d))` compares against it, and `min_x, min_y, max_x, max_y = box` works"""
        yield from (self.min_x, self.min_y, self.max_x, self.max_y)

    def __len__(self):
        """required for pytest `approx` sequence comparison alongside `__iter__`"""
        return 4


EMPTY_BOX = Box(inf, inf, -inf, -inf)


Segment = (tuple[str, float, float] | tuple[str, float, float, float, float, float, float] |
           tuple[str])


def move(x: float, y: float) -> Segment:
    return 'M', float(x), float(y)


def line(x: float, y: float) -> Segment:
    return 'L', float(x), float(y)


def curve(x1: float, y1: float, x2: float, y2: float, x: float, y: float) -> Segment:
    """A cubic: two control points then the end point.  The ONLY curve in the IR.

    Quadratics and arcs are absent -- PostScript has neither, and one curve type means the bounds code
    and every emitter have one case.
    """
    return 'C', float(x1), float(y1), float(x2), float(y2), float(x), float(y)


def close() -> Segment:
    return ('Z',)


def polyline(points: Sequence[tuple[float, float]], *, closed: bool = False) -> tuple[Segment, ...]:
    """The whole-subpath shorthand: a move to the first point and a line to each of the rest."""
    if not points:
        raise ValueError('polyline needs at least one point')
    out = [move(*points[0])]
    for x, y in points[1:]:
        out.append(line(x, y))
    if closed:
        out.append(close())
    return tuple(out)


def circle(cx: float, cy: float, r: float) -> tuple[Segment, ...]:
    """A closed circle as four cubics, starting at (cx + r, cy) and going counter-clockwise.

    Counter-clockwise in molecule coordinates (y-up), so an even-odd hole punched in a filled shape keeps
    a consistent winding whichever backend draws it.
    """
    if r <= 0.:
        raise ValueError(f'circle radius must be positive, got {r}')
    return ellipse(cx, cy, r, r)


def ellipse(cx: float, cy: float, rx: float, ry: float) -> tuple[Segment, ...]:
    """A closed ellipse as four cubics, starting at (cx + rx, cy) and going counter-clockwise.

    `circle` is this with one radius, so there is one arc approximation in the package and not two: KAPPA
    is per axis, the quarter arcs being independent in x and y.
    """
    if rx <= 0. or ry <= 0.:
        raise ValueError(f'ellipse radii must be positive, got {rx} and {ry}')
    kx = KAPPA * rx
    ky = KAPPA * ry
    return (move(cx + rx, cy),
            curve(cx + rx, cy + ky, cx + kx, cy + ry, cx, cy + ry),
            curve(cx - kx, cy + ry, cx - rx, cy + ky, cx - rx, cy),
            curve(cx - rx, cy - ky, cx - kx, cy - ry, cx, cy - ry),
            curve(cx + kx, cy - ry, cx + rx, cy - ky, cx + rx, cy),
            close())


def rounded_box(box: 'Box', radius: float) -> tuple[Segment, ...]:
    """`box` with its corners rounded, counter-clockwise from the middle of its right edge.

    `radius` is CLAMPED to half the shorter side rather than refused: the caller is a knock-out plate
    behind a one- or three-digit number, and the shape wanted at the limit is the stadium -- fully round
    ends -- not an exception for having asked for more than the box has room for.
    """
    if radius < 0.:
        raise ValueError(f'corner radius must not be negative, got {radius}')
    if box.min_x > box.max_x or box.min_y > box.max_y:
        raise ValueError('cannot round the corners of an empty box')
    r = min(radius, (box.max_x - box.min_x) / 2., (box.max_y - box.min_y) / 2.)
    if not r:
        return polyline(((box.max_x, box.min_y), (box.max_x, box.max_y), (box.min_x, box.max_y),
                         (box.min_x, box.min_y)), closed=True)
    k = KAPPA * r
    return (move(box.max_x, box.min_y + r),
            line(box.max_x, box.max_y - r),
            curve(box.max_x, box.max_y - r + k, box.max_x - r + k, box.max_y, box.max_x - r, box.max_y),
            line(box.min_x + r, box.max_y),
            curve(box.min_x + r - k, box.max_y, box.min_x, box.max_y - r + k, box.min_x, box.max_y - r),
            line(box.min_x, box.min_y + r),
            curve(box.min_x, box.min_y + r - k, box.min_x + r - k, box.min_y, box.min_x + r, box.min_y),
            line(box.max_x - r, box.min_y),
            curve(box.max_x - r + k, box.min_y, box.max_x, box.min_y + r - k, box.max_x, box.min_y + r),
            close())


_ANCHORS = frozenset(('start', 'middle', 'end'))
_CAPS = frozenset(('butt', 'round', 'square'))
_JOINS = frozenset(('miter', 'round', 'bevel'))


def _freeze_subpaths(subpaths) -> tuple[tuple[Segment, ...], ...]:
    out = []
    for subpath in subpaths:
        segments = tuple(subpath)
        if not segments:
            raise ValueError('empty subpath: a path with nothing in it is a caller bug')
        if segments[0][0] != 'M':
            raise ValueError(f'a subpath must start with M, got {segments[0][0]!r}')
        for segment in segments:
            if segment[0] not in ('M', 'L', 'C', 'Z'):
                raise ValueError(f'unknown segment {segment[0]!r}')
        out.append(segments)
    if not out:
        raise ValueError('a path needs at least one subpath')
    return tuple(out)


@dataclass(frozen=True, slots=True)
class Path:
    """Filled and/or stroked geometry.  `subpaths` is a sequence of segment sequences.

    Several subpaths in ONE path wherever the shapes are one drawing operation: a bond chained through
    unlabelled atoms, both halves of a double bond, a contour band and the hole inside it.  A filled band
    with a hole is only expressible that way.
    """
    subpaths: tuple[tuple[Segment, ...], ...]
    fill: str | None = None
    stroke: str | None = None
    width: float | None = None
    dashes: tuple[float, ...] | None = None
    cap: str = 'butt'
    join: str = 'miter'
    miter_limit: float | None = None
    even_odd: bool = False

    def __init__(self, subpaths, *, fill=None, stroke=None, width=None, dashes=None, cap='butt',
                 join='miter', miter_limit=None, even_odd=False):
        fill = to_hex(fill)
        stroke = to_hex(stroke)
        if fill is None and stroke is None:
            raise ValueError('a path must state a fill or stroke: an invisible path is a caller bug')
        if stroke is not None:
            if width is None:
                raise ValueError('a stroked path must state a width')
            if width <= 0.:
                raise ValueError(f'stroke width must be positive, got {width}')
        if cap not in _CAPS:
            raise ValueError(f'cap must be one of {sorted(_CAPS)}, got {cap!r}')
        if join not in _JOINS:
            raise ValueError(f'join must be one of {sorted(_JOINS)}, got {join!r}')
        if dashes is not None:
            dashes = tuple(float(d) for d in dashes)
            if not dashes or any(d <= 0. for d in dashes):
                raise ValueError(f'dash lengths must all be positive, got {dashes}')
        if miter_limit is not None:
            # Converted and range-checked here because the value is otherwise read only by a backend's
            # number formatter, which would raise far from this call site.  The floor is 1 because the
            # limit IS the ratio of miter length to stroke width, and that ratio cannot be below 1.
            miter_limit = float(miter_limit)
            if miter_limit < 1.:
                raise ValueError(f'miter limit is a ratio of miter length to stroke width and cannot '
                                 f'be below 1, got {miter_limit}')
        object.__setattr__(self, 'subpaths', _freeze_subpaths(subpaths))
        object.__setattr__(self, 'fill', fill)
        object.__setattr__(self, 'stroke', stroke)
        object.__setattr__(self, 'width', None if width is None else float(width))
        object.__setattr__(self, 'dashes', dashes)
        object.__setattr__(self, 'cap', cap)
        object.__setattr__(self, 'join', join)
        object.__setattr__(self, 'miter_limit', miter_limit)   # already float or None; see above
        object.__setattr__(self, 'even_odd', even_odd)

    @property
    def bounds(self) -> Box:
        """The control-point hull, inflated by half the stroke width.

        The hull, not the true curve extent: it is never smaller than the curve, and a slightly generous
        viewBox costs whitespace where a slightly tight one crops the picture.
        """
        min_x = min_y = inf
        max_x = max_y = -inf
        for subpath in self.subpaths:
            for segment in subpath:
                for i in range(1, len(segment), 2):
                    x = segment[i]
                    y = segment[i + 1]
                    if x < min_x:
                        min_x = x
                    if x > max_x:
                        max_x = x
                    if y < min_y:
                        min_y = y
                    if y > max_y:
                        max_y = y
        if min_x > max_x:
            return EMPTY_BOX
        box = Box(min_x, min_y, max_x, max_y)
        if self.stroke is not None:
            box = box.inflate(self.width / 2.)
        return box

    def translated(self, dx: float, dy: float) -> 'Path':
        """The same path, moved.  Every coordinate in a segment shifts; nothing else changes.

        A shift on the NODE rather than a transform attribute on a group, so `bounds` needs no
        composition and no backend has to express one.
        """
        moved = []
        for subpath in self.subpaths:
            segments = []
            for segment in subpath:
                shifted = [segment[0]]
                for i in range(1, len(segment), 2):
                    shifted.append(segment[i] + dx)
                    shifted.append(segment[i + 1] + dy)
                segments.append(tuple(shifted))
            moved.append(tuple(segments))
        return replace(self, subpaths=tuple(moved))


@dataclass(frozen=True, slots=True)
class TextRun:
    """One span of a label: its own family, size, weight and offset from the label's anchor.

    `dy` is in molecule units and y-UP, so a subscript has a NEGATIVE dy; the backends flip it with the
    rest of the scene.  The two offsets compose differently, and the measurer (`metrics.text_box`) and
    every backend must agree: `dy` is ABSOLUTE -- each run's own shift from the LABEL's baseline, so a
    slot's height cannot depend on which other slots the label holds -- while `dx` is CUMULATIVE, extra
    advance where the previous run left the pen.  SVG's `<tspan dy>` is cumulative, so `render/svg.py`
    emits the difference between consecutive runs' `dy`.
    """
    text: str
    family: str = 'helvetica'
    size: float = .4
    weight: str = 'normal'
    style: str = 'normal'
    dx: float = 0.
    dy: float = 0.

    def __post_init__(self):
        if not self.text:
            raise ValueError('an empty text run draws nothing; leave it out')
        if self.size <= 0.:
            raise ValueError(f'text size must be positive, got {self.size}')
        if self.weight not in ('normal', 'bold'):
            raise ValueError(f'weight must be normal or bold, got {self.weight!r}')
        if self.style not in ('normal', 'italic'):
            raise ValueError(f'style must be normal or italic, got {self.style!r}')


@dataclass(frozen=True, slots=True)
class Text:
    """One anchored label, made of runs.  `x`/`y` is the anchor point on the text's BASELINE.

    A whole label -- `CH3`, `NH2+`, `13C` -- is one `Text`, because it is positioned as a unit and every
    backend has its own way of advancing between runs.
    """
    runs: tuple[TextRun, ...]
    x: float = 0.
    y: float = 0.
    anchor: str = 'start'
    fill: str | None = None

    def __init__(self, runs, *, x=0., y=0., anchor='start', fill=None):
        runs = tuple(runs)
        if not runs:
            raise ValueError('a text needs at least one run')
        if anchor not in _ANCHORS:
            raise ValueError(f'anchor must be one of {sorted(_ANCHORS)}, got {anchor!r}')
        object.__setattr__(self, 'runs', runs)
        object.__setattr__(self, 'x', float(x))
        object.__setattr__(self, 'y', float(y))
        object.__setattr__(self, 'anchor', anchor)
        object.__setattr__(self, 'fill', to_hex(fill))

    @property
    def bounds(self) -> Box:
        """The tight ink box, measured through `metrics`.

        `metrics` is imported here and not at module scope: it reads a shipped table off disk on first
        use, and `scene.py` is imported to build a two-line path as often as to draw a molecule.
        """
        from .metrics import text_box

        return text_box(self)

    def translated(self, dx: float, dy: float) -> 'Text':
        return replace(self, x=self.x + dx, y=self.y + dy)


@dataclass(frozen=True, slots=True)
class Group:
    """Children that composite together.  `opacity` applies to the group's flattened result.

    GROUP opacity is the overlap mechanism, and why no polygon boolean union exists in this package:
    overlapping translucent children each composite separately and show their seams, while the same
    children inside one `Group(opacity=...)` are flattened first and read as one translucent shape.

    `clip` is an optional `Path` whose interior is what shows.
    """
    # THE WHOLE UNION IS QUOTED, not just `Group`: this is a dataclass, so the annotation is
    # evaluated when the class is created, and `'Group' | Path` is a `str.__or__` -- a TypeError.
    children: tuple['Group | Path | Text', ...]
    opacity: float | None = None
    clip: Path | None = None

    def __init__(self, children, *, opacity=None, clip=None):
        if opacity is not None and not 0. <= opacity <= 1.:
            raise ValueError(f'opacity must be within [0, 1], got {opacity}')
        object.__setattr__(self, 'children', tuple(children))
        object.__setattr__(self, 'opacity', None if opacity is None else float(opacity))
        object.__setattr__(self, 'clip', clip)

    @property
    def bounds(self) -> Box:
        box = EMPTY_BOX
        for child in self.children:
            box = box.union(child.bounds)
        if self.clip is not None:      # clipped geometry cannot exceed the clip
            box = box.union(EMPTY_BOX) if box.min_x > box.max_x else _intersect(box, self.clip.bounds)
        return box

    def translated(self, dx: float, dy: float) -> 'Group':
        return replace(self, children=tuple(c.translated(dx, dy) for c in self.children),
                       clip=None if self.clip is None else self.clip.translated(dx, dy))


def _intersect(a: Box, b: Box) -> Box:
    box = Box(max(a.min_x, b.min_x), max(a.min_y, b.min_y), min(a.max_x, b.max_x), min(a.max_y, b.max_y))
    return EMPTY_BOX if box.min_x > box.max_x or box.min_y > box.max_y else box


@dataclass(frozen=True, slots=True)
class Scene:
    """A whole picture, and the object a caller holds.  Serialization hangs off it.

    `bounds` may be STATED, and then it wins over the union of the children -- a grid cell, a frame in a
    series, or a figure that must match another one need a fixed frame.
    """
    children: tuple[Group | Path | Text, ...]
    _bounds: Box | None = None

    def __init__(self, children, *, bounds=None):
        object.__setattr__(self, 'children', tuple(children))
        object.__setattr__(self, '_bounds', bounds)

    @property
    def bounds(self) -> Box:
        if self._bounds is not None:
            return self._bounds
        box = EMPTY_BOX
        for child in self.children:
            box = box.union(child.bounds)
        return Box(0., 0., 0., 0.) if box.min_x > box.max_x else box

    def frame(self, margin: float) -> Box:
        """The box to render: the STATED bounds as given, or the computed union inflated by `margin`.

        A caller states `bounds` to FIX the frame, so no margin is added to it; state the box with the
        breathing room already in it.  Lives here so every backend answers identically.
        """
        if self._bounds is not None:
            return self._bounds
        return self.bounds.inflate(margin)

    def to_svg(self, **kwargs) -> str:
        from .render.svg import to_svg

        return to_svg(self, **kwargs)

    def to_svgz(self, **kwargs) -> bytes:
        from .render.svg import to_svgz

        return to_svgz(self, **kwargs)

    def to_pdf(self, **kwargs) -> bytes:
        raise NotImplementedError('no PDF backend before chython 3.1; use to_svg() and convert, or '
                                  'set page sizes in mm and let the SVG carry them')

    def to_eps(self, **kwargs) -> bytes:
        raise NotImplementedError('no EPS backend before chython 3.1; use to_svg() and convert, or '
                                  'set page sizes in mm and let the SVG carry them')

    def _repr_svg_(self) -> str:
        return self.to_svg()
