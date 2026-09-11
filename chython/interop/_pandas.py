# -*- coding: utf-8 -*-
#
#  Copyright 2021-2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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
"""
Teach pandas that a container is one value and not a list of atoms.

`pandas.io.formats.printing.is_sequence(x)` asks whether `iter(x)` and `len(x)` both work, and for a
`MoleculeContainer` they do, so a column of molecules would render as a column of atom-number lists.
The predicate here is the whole of `is_container`, so a container that later grows `__len__` cannot
regress this silently.
"""

__all__ = ['patch_pandas']

_patched = False


def patch_pandas():
    """
    Render chython containers as single values in pandas output rather than as lists of atoms.

    Idempotent: patching twice would wrap the wrapper, so a second call returns.
    """
    global _patched

    if _patched:
        return
    _patched = True

    from pandas.io.formats import printing
    from pandas.io.formats.printing import is_sequence

    from . import is_container

    def patched(obj):
        # `is_container` first: two `isinstance` calls, where `is_sequence` builds an iterator
        return False if is_container(obj) else is_sequence(obj)

    printing.is_sequence = patched
