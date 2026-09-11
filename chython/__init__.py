# -*- coding: utf-8 -*-
#
#  Copyright 2014-2026 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2014-2019 Timur Madzhidov tmadzhidov@gmail.com features and API discussion
#  Copyright 2014-2019 Alexandre Varnek <varnek@unistra.fr> base idea of CGR approach
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
"""chython's public surface: a facade re-exporting the packages below it.

`__all__` is empty by design.  `smarts` is the short spelling of `read_smarts` and the same function
object.  `smiles` is bidirectional: a string in reads -- a `>` that is not a dative `->` makes it a
reaction SMILES and the result a `ReactionContainer` -- and a container in writes one.  `pach` is the
same door for the wire format, and `unpach`/`unpack` its import half under chython 2's two names.
"""
from sys import modules as _modules
from types import ModuleType as _ModuleType
from .core import *
from .core import read_smarts as smarts
from .depict import (Clean2DEngine, DepictStyle, get_clean2d_engine, get_depict_style,
                     set_clean2d_engine, set_depict_style)
from .formats import *
# By full path, not through `formats`' star: `pdb` is that subpackage's name too and the function
# would shadow it.
from .formats.pdb import PDBAtom, PDBBond, PDBRecord, build_molecule, mmcif, pdb, read_mmcif, read_pdb
# Imported for its registration side effect as much as for its names: it calls `_set_standardize_fn`
# at import time, which is what makes `mol.standardize()` exist.  Not an unused import.
from .chemistry import *
# Likewise: `_set_reactions_fns` at import time is what makes `mol.react()` and `mol @ other` exist.
from .reactions import *
from .interop import iupac, patch_pandas
from .interop.config import _facade_alias as _interop_facade_alias


class _Facade(_ModuleType):
    """Gives `chython` itself a property, so `chython.clean2d_engine` forwards both the read and the
    write to its one home in `depict/_config.py` and the setter validates the name on the spot.
    """
    @property
    def clean2d_engine(self) -> Clean2DEngine:
        return get_clean2d_engine()

    @clean2d_engine.setter
    def clean2d_engine(self, engine: Clean2DEngine):
        set_clean2d_engine(engine)


_modules[__name__].__class__ = _Facade

# `conformer_engine` and `class_paths` live in `chython.interop.config`; aliased here rather than
# copied, or `chython.conformer_engine = 'cdpkit'` would be a silent no-op.  Must run AFTER the
# `__class__` assignment above: `_facade_alias` subclasses whatever class the module currently has, so
# the reverse order would replace the aliasing subclass and turn both names into AttributeErrors.
_interop_facade_alias(__name__, 'conformer_engine', 'class_paths')

__all__ = []
