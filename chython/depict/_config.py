# -*- coding: utf-8 -*-
#
#  Copyright 2018-2026 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2019-2020 Dinar Batyrshin <batyrshin-dinar@mail.ru>
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
"""The CPK palette, and the choice of 2D layout engine.

Owned here rather than by the facade, which only proxies `clean2d_engine`: state read by `depict/` is
`depict/`'s, and nothing may import the facade.  That is what makes `import chython.depict` work on
its own.  Every rendering parameter lives in `style.py` instead, as a field of `DepictStyle`.
"""
from typing import Literal


Clean2DEngine = Literal['rdkit', 'smilesdrawer', 'cdk', 'obabel', 'indigo']

# Which 2D layout backend `clean2d()` uses when the caller does not name one.  `smilesdrawer` is a
# JavaScript bundle on QuickJS: ~1 MB installed, and its behaviour does not move between releases.
clean2d_engine: Clean2DEngine = 'smilesdrawer'

cpk = tuple('''
 #909090                                                                                         #D9FFFF
 #CC80FF #C2FF00                                         #FFB5B5 #101010 #3050F8 #FF0D0D #90E050 #B3E3F5
 #AB5CF2 #8AFF00                                         #BFA6A6 #F0C8A0 #FF8000 #C6C600 #1FF01F #80D1E3
 #8F40D4 #3DFF00 #E6E6E6 #BFC2C7 #A6A6AB #8A99C7 #9C7AC7
                 #E06633 #F090A0 #50D050 #C88033 #7D80B0 #C28F8F #668F8F #BD80E3 #FFA100 #A62929 #5CB8D1
 #702EB0 #00FF00 #94FFFF #94E0E0 #73C2C9 #54B5B5 #3B9E9E
                 #248F8F #0A7D8C #006985 #C0C0C0 #FFD98F #A67573 #668080 #9E63B5 #D47A00 #940094 #429EB0
 #57178F #00C900 #70D4FF
                 #FFFFC7 #D9FFC7 #C7FFC7 #A3FFC7 #8FFFC7 #61FFC7 #45FFC7
                 #30FFC7 #1FFFC7 #00FF9C #00E675 #00D452 #00BF38 #00AB24
                         #4DC2FF #4DA6FF #2194D6 #267DAB
                 #266696 #175487 #D0D0E0 #FFD123 #B8B8D0 #A6544D #575961 #9E4FB5 #AB5C00 #754F45 #428296
 #420066 #007D00 #70ABFA
                 #00BAFF #00A1FF #008FFF #0080FF #006BFF #545CF2 #785CE3
                 #8A4FE3 #A136D4 #B31FD4 #B31FBA #B30DA6 #BD0D87 #C70066
                         #CC0059 #D1004F #D90045 #E00038
                 #E6002E #EB0026 #EB0026 #EB0026 #EB0026 #EB0026 #EB0026 #EB0026 #EB0026 #EB0026 #EB0026
'''.split())

# The marker's colour.  Not a `cpk` row: that table is indexed by Z - 1 and a marker has no Z, and a
# 119th entry would be read as an element by everything that iterates it.  Dark grey, so an R reads as
# an attachment point rather than as a halogen or a metal.
R_COLOUR = '#404040'


def get_clean2d_engine(engine: Clean2DEngine = None) -> Clean2DEngine:
    """Resolve the layout engine for one call: the argument if given, otherwise the module default.

    An accessor rather than a direct read of the global, so that a rebinding through
    `chython.clean2d_engine = 'rdkit'` is visible to every downstream caller.
    """
    return clean2d_engine if engine is None else engine


def set_clean2d_engine(engine: Clean2DEngine):
    """Set the default layout engine.  Validated here so a typo fails at assignment, not at draw."""
    if engine not in ('rdkit', 'smilesdrawer', 'cdk', 'obabel', 'indigo'):
        raise ValueError(f'Invalid clean2d engine: {engine}')
    global clean2d_engine
    clean2d_engine = engine


__all__ = ['cpk', 'R_COLOUR', 'Clean2DEngine', 'clean2d_engine', 'get_clean2d_engine', 'set_clean2d_engine']
