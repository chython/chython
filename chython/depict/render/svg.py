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
"""Scene -> SVG.

The y-flip is per coordinate, not a document `scale(S, -S)`: a negative scale mirrors glyphs.  No masks,
CSS or generated ids -- clip ids are sequential, so the output is deterministic and diffable.  The viewBox
is in molecule units (one bond is 1.0) while `width`/`height` are physical.
"""
from gzip import compress
from math import isfinite
from xml.sax.saxutils import escape, quoteattr

from ..metrics import SVG_FAMILY
from ..scene import Group, Path, Scene, Text
from ..style import DepictStyle, get_depict_style


__all__ = ['format_number', 'to_svg', 'to_svgz']


def format_number(value: float) -> str:
    """A coordinate as the shortest decimal that means it, with no exponent and no negative zero.

    Shared by all backends, so one number is written one way everywhere.  No exponent because not every
    SVG path parser accepts one; no `-0`, so two renders of one picture cannot differ as strings.

    :raises ValueError: `value` is not finite -- `inf`/`nan` are not valid SVG path data.
    """
    if not isfinite(value):
        raise ValueError(f'coordinate is not finite: {value!r}')
    if value == 0.:            # catches -0.0, which formats as '-0'
        return '0'
    text = f'{value:.4f}'.rstrip('0').rstrip('.')
    return '0' if text in ('-0', '') else text


def _path_data(path: Path) -> str:
    """A `Path`'s subpaths as SVG path data, y negated as each number is written."""
    out = []
    for subpath in path.subpaths:
        for segment in subpath:
            command = segment[0]
            if command == 'Z':
                out.append('Z')
            elif command == 'C':
                out.append('C%s %s %s %s %s %s' % (
                    format_number(segment[1]), format_number(-segment[2]),
                    format_number(segment[3]), format_number(-segment[4]),
                    format_number(segment[5]), format_number(-segment[6])))
            else:
                out.append('%s%s %s' % (command, format_number(segment[1]),
                                        format_number(-segment[2])))
    return ''.join(out)


def _paint_attributes(path: Path) -> str:
    """Fill, stroke and the line attributes -- only the ones that differ from SVG's defaults.

    `fill="none"` IS written when there is no fill: SVG's default fill is black, so a stroked outline
    would otherwise come out as a filled blob.
    """
    out = [' fill=%s' % quoteattr(path.fill if path.fill is not None else 'none')]
    if path.even_odd:
        out.append(' fill-rule="evenodd"')
    if path.stroke is not None:
        out.append(' stroke=%s' % quoteattr(path.stroke))
        out.append(' stroke-width="%s"' % format_number(path.width))
        if path.cap != 'butt':
            out.append(' stroke-linecap=%s' % quoteattr(path.cap))
        if path.join != 'miter':
            out.append(' stroke-linejoin=%s' % quoteattr(path.join))
        elif path.miter_limit is not None:
            out.append(' stroke-miterlimit="%s"' % format_number(path.miter_limit))
        if path.dashes is not None:
            out.append(' stroke-dasharray="%s"' % ' '.join(format_number(d) for d in path.dashes))
    return ''.join(out)


def _text_element(text: Text) -> str:
    """One `<text>` with a `<tspan>` per run.

    One element per LABEL, not per run: it is anchored once and the runs advance from each other the way
    the metrics measured them.  `dy` is negated with the rest of the scene and emitted as the DIFFERENCE
    between consecutive runs, because `TextRun.dy` is absolute from the label's baseline while SVG's
    `<tspan dy>` is cumulative; writing each run's own `dy` puts every run after a shifted one at the sum
    of the shifts.  This is the only place that conversion happens.
    """
    out = ['<text x="%s" y="%s"' % (format_number(text.x), format_number(-text.y))]
    if text.anchor != 'start':
        out.append(' text-anchor=%s' % quoteattr(text.anchor))
    if text.fill is not None:
        out.append(' fill=%s' % quoteattr(text.fill))
    out.append('>')
    previous_dy = 0.
    for run in text.runs:
        out.append('<tspan font-family=%s font-size="%s"'
                   % (quoteattr(SVG_FAMILY[run.family]),
                      format_number(run.size)))
        if run.weight != 'normal':
            out.append(' font-weight=%s' % quoteattr(run.weight))
        if run.style != 'normal':
            out.append(' font-style=%s' % quoteattr(run.style))
        if run.dx:
            out.append(' dx="%s"' % format_number(run.dx))
        step = run.dy - previous_dy
        if step:
            out.append(' dy="%s"' % format_number(-step))
        previous_dy = run.dy
        out.append('>%s</tspan>' % escape(run.text))
    out.append('</text>')
    return ''.join(out)


def _render(node, out: list, clips: list):
    """Append `node`'s markup to `out`, collecting any clip paths into `clips`.

    A `Group` with neither opacity nor a clip emits no wrapper: only those two make grouping visible.
    """
    if isinstance(node, Path):
        out.append('<path%s d="%s"/>' % (_paint_attributes(node), _path_data(node)))
    elif isinstance(node, Text):
        out.append(_text_element(node))
    elif isinstance(node, Group):
        attributes = []
        if node.opacity is not None:
            attributes.append(' opacity="%s"' % format_number(node.opacity))
        if node.clip is not None:
            identifier = 'clip%d' % len(clips)
            clips.append((identifier, node.clip))
            attributes.append(' clip-path="url(#%s)"' % identifier)
        if attributes:
            out.append('<g%s>' % ''.join(attributes))
            for child in node.children:
                _render(child, out, clips)
            out.append('</g>')
        else:
            for child in node.children:
                _render(child, out, clips)
    else:
        raise TypeError(f'{type(node).__name__} is not a scene primitive')


def to_svg(scene: Scene, *, style: DepictStyle | None = None, standalone: bool = False) -> str:
    """Serialize a scene.  `standalone` adds the XML declaration, for a file rather than a notebook.

    Deterministic: the same scene and style give byte-identical output, every time.
    """
    if style is None:
        style = get_depict_style()
    page = style.page
    # `scene.frame()` inflates only a COMPUTED bounds; a STATED one is the caller's fixed frame.
    box = scene.frame(page.margin)
    width = box.width
    height = box.height

    # `PageStyle.__post_init__` refuses neither-set, so exactly one of `width_mm`/`scale_mm` is set.
    if page.width_mm is not None:
        # A zero width makes the scale meaningless; 0. is safe, since the SVG spec already disables
        # rendering of a zero-extent viewBox.
        scale = page.width_mm / width if width > 0. else 0.
    else:
        scale = page.scale_mm

    body = []
    clips = []
    if page.background is not None:
        body.append('<rect x="%s" y="%s" width="%s" height="%s" fill=%s/>'
                    % (format_number(box.min_x), format_number(-box.max_y), format_number(width),
                       format_number(height), quoteattr(page.background)))
    for child in scene.children:
        _render(child, body, clips)

    # Clips are emitted BEFORE the body: a forward reference to a `<clipPath>` is legal SVG but not
    # universally handled -- Cairo-based rasterizers among the offenders.
    defs = []
    if clips:
        defs.append('<defs>')
        for identifier, clip in clips:
            clip_rule = ' clip-rule="evenodd"' if clip.even_odd else ''
            defs.append('<clipPath id="%s"><path%s d="%s"/></clipPath>'
                        % (identifier, clip_rule, _path_data(clip)))
        defs.append('</defs>')

    header = ('<svg xmlns="http://www.w3.org/2000/svg" width="%smm" height="%smm" viewBox="%s %s %s %s">'
              % (format_number(width * scale), format_number(height * scale),
                 format_number(box.min_x), format_number(-box.max_y),
                 format_number(width), format_number(height)))
    document = header + ''.join(defs) + ''.join(body) + '</svg>'
    if standalone:
        return '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n' + document + '\n'
    return document


def to_svgz(scene: Scene, *, style: DepictStyle | None = None, standalone: bool = True) -> bytes:
    """The gzipped document, which is what `.svgz` is.

    `mtime=0`, so two renders of one figure are byte-identical -- gzip would otherwise write the
    current time.
    """
    return compress(to_svg(scene, style=style, standalone=standalone).encode('utf-8'), mtime=0)
