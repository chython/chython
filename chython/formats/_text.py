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
"""One refusal, shared by the read-only facades.

Each takes a document as text and reads it; a `bytes` or a `Path` is a caller who wanted the file reader.
Saying so by name beats the `AttributeError` a str-only parser raises three frames down.
"""

__all__ = ['require_text']


def require_text(data, what):
    """`data` if it is a `str`, else `TypeError` naming `what`."""
    if isinstance(data, str):
        return data
    raise TypeError(f'{what}() takes {what.upper()} text as a str, not '
                    f'{type(data).__name__}; decode bytes, or use the file reader for a path')
