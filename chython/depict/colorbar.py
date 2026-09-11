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
"""The colour legend: a discrete swatch strip with min / max / zero ticks.

Discrete, and the swatches are passed in rather than recomputed, so the bar cannot advertise a colour
the picture does not contain.  The bar's axis is the colormap's *domain* -- comparable between figures
sharing one scale, and it leaves a gap where the field never reached a level.  Coordinates are local to
the strip's lower-left corner, y up; only `place_colorbar` knows where on the page it goes.
"""
from collections.abc import Sequence

from .label import MINUS
from .overlay import Swatch, scale_of
from .scene import Box, Group, Path, Text, TextRun, polyline


__all__ = ['colorbar', 'legend_side', 'place_colorbar']


def colorbar(cmap, *, swatches: Sequence[Swatch], style, side: str,
             length: float) -> tuple[list, Box]:
    """The swatches at their places along `cmap`'s domain, plus the ticks, in local coordinates.

    `swatches` is a sequence of `overlay.Swatch` -- exactly what the figure drew, from `overlay.bands_of`;
    nothing is re-derived here.  They are sorted ascending by value, because a bar reads low-to-high while
    `bands_of` is in painter's order, and those differ for a diverging field.

    `cmap` supplies the axis as well as the ticks: a swatch spans `normalised(lo) … normalised(hi)` of
    `length`, the same clamped map the fills were coloured through.  A swatch with no interval -- an open
    contour, or every level in `fill=False` mode -- is drawn as a rule at its level.
    """
    swatches = sorted(swatches)
    if not swatches:
        return [], Box(0., 0., 0., 0.)

    breadth = style.page.legend_breadth
    family = style.label.family
    size = style.label.size * style.field.value_scale

    nodes: list = []
    for swatch in swatches:
        lo, hi = cmap.normalised(swatch.lo) * length, cmap.normalised(swatch.hi) * length
        # `hi <= lo` is a band the domain cannot show -- a level past `vmax`, which `normalised` clamps to
        # its neighbour's end.  Drawn as a rule, because a zero-width rectangle is an invisible swatch.
        if swatch.filled and hi > lo:
            if side == 'bottom':
                rect = ((lo, 0.), (hi, 0.), (hi, breadth), (lo, breadth))
            else:  # 'right'
                rect = ((0., lo), (breadth, lo), (breadth, hi), (0., hi))
            nodes.append(Path([polyline(rect, closed=True)], fill=swatch.colour))
        else:
            at = cmap.normalised(swatch.level) * length
            rule = ((at, 0.), (at, breadth)) if side == 'bottom' else ((0., at), (breadth, at))
            nodes.append(Path([polyline(rule)], stroke=swatch.colour,
                              width=style.field.contour.line_width))

    # the '+' in the tick format is kept for a two-sided domain and dropped for a one-sided one
    spans_zero = cmap.vmin < 0. < cmap.vmax
    fmt = style.page.legend_fmt if spans_zero else style.page.legend_fmt.replace('+', '')

    # proportional to the strip breadth, so narrowing the strip narrows the gap: one knob, not two
    tick_gap = breadth * .15

    tick_values = [cmap.vmin, cmap.vmax]
    if spans_zero:
        tick_values.append(0.)

    for tick_val in tick_values:
        formatted = fmt.format(tick_val)
        # the typographic minus, not the ASCII hyphen: at print scale the two differ visibly
        if formatted.startswith('-'):
            formatted = MINUS + formatted[1:]

        t = cmap.normalised(tick_val)

        # Set on the baseline first and then moved by its own MEASURED ink, because a baseline is not a
        # box: glyphs sit above it, so a horizontal bar's tick placed at `-tick_gap` would have its
        # digits inside the strip, and a vertical bar's would read half a cap height above its value.
        if side == 'bottom':
            tick = Text([TextRun(formatted, family=family, size=size)],
                        x=t * length, y=0., anchor='middle')
            tick = tick.translated(0., -tick_gap - tick.bounds.max_y)
        else:  # 'right'
            tick = Text([TextRun(formatted, family=family, size=size)],
                        x=breadth + tick_gap, y=0., anchor='start')
            ink = tick.bounds
            tick = tick.translated(0., t * length - (ink.min_y + ink.max_y) / 2.)

        nodes.append(tick)

    box = Box.of([node.bounds for node in nodes])
    return nodes, box


def place_colorbar(nodes: list, box: Box, content: Box, style, side: str) -> Group:
    """Translate the colorbar so it sits outside `content` on `side`, separated by `style.page.margin`.

    Returns one `Group`, so the legend is a single scene child at a known position (last).
    """
    if side == 'right':
        dx = content.max_x + style.page.margin - box.min_x
        dy = content.min_y + (content.height - box.height) / 2. - box.min_y
    else:  # 'bottom'
        dx = content.min_x + (content.width - box.width) / 2. - box.min_x
        dy = content.min_y - style.page.margin - box.max_y
    return Group([node.translated(dx, dy) for node in nodes])


def legend_side(style, overlays, content: Box) -> str | None:
    """Which side the bar goes on, or None for no bar.

    `'auto'` reads the content's shape, not the page's: a tall molecule leaves free width and a wide one
    free height, so the bar goes where the space already is and the figure does not grow lengthwise.

    :raises ValueError: the overlays carry two different scales, which one bar cannot label.
    """
    if style.page.legend == 'none':
        return None
    scales = [s for s in (scale_of(o) for o in overlays) if s is not None]
    if not scales:
        return None
    first = scales[0]
    for other in scales[1:]:
        if (other.vmin, other.vmax) != (first.vmin, first.vmax) or other.stops != first.stops:
            raise ValueError('these overlays carry two different scales and one colorbar cannot label '
                             'both -- draw them as separate figures, or set page.legend to "none"')
    if style.page.legend != 'auto':
        return style.page.legend
    return 'right' if content.height >= content.width else 'bottom'
