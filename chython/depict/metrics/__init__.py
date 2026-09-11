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
"""Glyph metrics, read from the shipped TSVs.

Loaded lazily into a per-family cache on first use, never at import.
"""
from csv import reader
from functools import lru_cache
from importlib.resources import files
from typing import NamedTuple

from ..scene import Box, Text


__all__ = ['FAMILIES', 'Glyph', 'METRICS_DIRECTORY', 'PDF_BASE_FONT', 'advance', 'font_metrics',
           'text_box']


# The one way to locate a shipped table.  A `Traversable`, not a `Path`, because that is what works
# from inside a wheel or a zipimport; exported so tests and the loader resolve the TSVs identically.
METRICS_DIRECTORY = files(__package__)

FAMILIES = frozenset(('helvetica', 'times'))

# The PDF base-14 name for each family.  A base-14 font is NAMED in the PDF and not embedded, which is
# why the family list is these two: one metric table is then exact in SVG, PDF and EPS alike.
PDF_BASE_FONT = {'helvetica': 'Helvetica', 'times': 'Times-Roman'}

# The PostScript name per family and (weight, style).  Consumed by the EPS backend; render/svg.py uses
# SVG_FAMILY instead, since CSS matches family names and not PostScript font names.
PS_NAME = {
    ('helvetica', 'normal', 'normal'): 'Helvetica',
    ('helvetica', 'bold', 'normal'): 'Helvetica-Bold',
    ('helvetica', 'normal', 'italic'): 'Helvetica-Oblique',
    ('helvetica', 'bold', 'italic'): 'Helvetica-BoldOblique',
    ('times', 'normal', 'normal'): 'Times-Roman',
    ('times', 'bold', 'normal'): 'Times-Bold',
    ('times', 'normal', 'italic'): 'Times-Italic',
    ('times', 'bold', 'italic'): 'Times-BoldItalic',
}

# The SVG/CSS `font-family` list for each family: real family names first, generic last.  NOT `PS_NAME`
# -- CSS matches FAMILY names, so a viewer handed `Times-Roman` silently falls through to its generic
# serif, whose widths are not `times.tsv`'s.
SVG_FAMILY = {'helvetica': 'Helvetica,Arial,sans-serif', 'times': '"Times New Roman",Times,serif'}


class Glyph(NamedTuple):
    """One glyph's advance width and ink box, in 1/1000 em -- the AFM's own units, unscaled.

    Unscaled because a glyph is measured at many sizes in one picture; the division belongs at the point
    of use, where the size is known.
    """
    char: str
    name: str
    wx: int
    llx: int
    lly: int
    urx: int
    ury: int


@lru_cache(maxsize=None)
def font_metrics(family: str) -> dict[str, Glyph]:
    """`{character: Glyph}` for one family.  Cached; the TSV is read once per process.

    Keyed by CHARACTER, not by glyph name: a caller measuring `'Cl'` has characters.  Glyphs the AFM
    carries outside the default encoding (`minus`, `bullet`) are keyed by the character the generator
    resolved them to; the ones it could not resolve are absent, so measuring one is a KeyError rather
    than a guessed width.

    :raises ValueError: `family` is not one of `FAMILIES`.
    """
    if family not in FAMILIES:
        raise ValueError(f'unknown font family {family!r}: expected one of {sorted(FAMILIES)}')
    out = {}
    text = METRICS_DIRECTORY.joinpath(f'{family}.tsv').read_text(encoding='utf-8')
    rows = reader((line for line in text.splitlines() if line and not line.startswith('#')),
                  delimiter='\t', quotechar=None)
    next(rows)   # the header
    for char, name, wx, llx, lly, urx, ury in rows:
        if char == '-' and name != 'hyphen':
            continue    # in the font but reachable by glyph name only, so no label can contain it
        out[char] = Glyph(char, name, int(wx), int(llx), int(lly), int(urx), int(ury))
    return out


def advance(text: str, family: str, size: float) -> float:
    """The pen advance for `text` at `size`, in molecule units.

    Sum of the glyph widths, NOT kerned: a kerned advance would disagree with what a viewer renders
    unless the backend also emitted the kerns.
    """
    table = font_metrics(family)
    total = 0
    for char in text:
        try:
            total += table[char].wx
        except KeyError:
            raise KeyError(f'{char!r} has no metrics in {family}: it cannot appear in a label') from None
    return total / 1000. * size


def text_box(text: Text) -> Box:
    """The tight INK box of a whole `Text`, with its anchor applied.

    Ink and not the em box: the label knock-out has to be the size of the drawn symbol, or a symbol sits
    in a hole too big for it and the bonds are trimmed too short.  Padding is one style field, applied
    once, to a measured box.
    """
    pen = 0.
    min_x = min_y = float('inf')
    max_x = max_y = float('-inf')
    for run in text.runs:
        table = font_metrics(run.family)
        pen += run.dx
        scale = run.size / 1000.
        for char in run.text:
            try:
                glyph = table[char]
            except KeyError:
                raise KeyError(f'{char!r} has no metrics in {run.family}: '
                               'it cannot appear in a label') from None
            left = pen + glyph.llx * scale
            right = pen + glyph.urx * scale
            bottom = run.dy + glyph.lly * scale
            top = run.dy + glyph.ury * scale
            if left < min_x:
                min_x = left
            if right > max_x:
                max_x = right
            if bottom < min_y:
                min_y = bottom
            if top > max_y:
                max_y = top
            pen += glyph.wx * scale
    if min_x > max_x:                       # every run empty, which `TextRun` already refuses
        return Box(text.x, text.y, text.x, text.y)
    if text.anchor == 'middle':
        shift = -pen / 2.
    elif text.anchor == 'end':
        shift = -pen
    else:
        shift = 0.
    return Box(text.x + min_x + shift, text.y + min_y, text.x + max_x + shift, text.y + max_y)
