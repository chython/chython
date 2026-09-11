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
"""Every importer that discards a source's coordinates says so.

A source-level check, scoped to one function body: the behavioural version needs three toolkits
installed to produce one log line, and would skip on every machine that matters.
"""

import re
from pathlib import Path

from pytest import mark


_INTEROP = Path(__file__).resolve().parent.parent


def _function_body(source, entry):
    """Return source text from ``entry`` to the next top-level ``def``, or EOF.

    The boundary is a newline plus 'def ' at column 0, so a nested def stays inside the body.
    """
    start = source.index(entry)
    rest = source[start + len(entry):]
    m = re.search(r'\ndef ', rest)
    if m is None:
        return source[start:]
    return source[start: start + len(entry) + m.start() + 1]


# _cdpkit.py is deliberately absent: from_cdpkit raises DirectionNotImplemented unconditionally, so it
# discards nothing.  Whoever builds that direction adds both the log line and this row.
@mark.parametrize('filename,entry', [('_indigo.py', 'def from_indigo'),
                                     ('_cdk.py', 'def _from_cdk'),
                                     ('_openbabel.py', 'def from_openbabel')])
def test_importer_mentions_the_coordinate_drop(filename, entry):
    source = (_INTEROP / filename).read_text(encoding='utf8')
    body = _function_body(source, entry)
    assert 'coordinates are not imported' in body, f'{filename}: {entry} drops coordinates silently'
