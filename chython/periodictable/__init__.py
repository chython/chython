# -*- coding: utf-8 -*-
#
#  Copyright 2018-2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from abc import ABCMeta
from .element import *
from .groups import *
from .periods import *
from .groupI import *
from .groupII import *
from .groupIII import *
from .groupIV import *
from .groupV import *
from .groupVI import *
from .groupVII import *
from .groupVIII import *
from .groupIX import *
from .groupX import *
from .groupXI import *
from .groupXII import *
from .groupXIII import *
from .groupXIV import *
from .groupXV import *
from .groupXVI import *
from .groupXVII import *
from .groupXVIII import *

modules = {v.__name__: v for k, v in globals().items() if k.startswith('group') and k != 'groups'}
elements = {k: v for k, v in globals().items() if isinstance(v, ABCMeta) and k != 'Element' and issubclass(v, Element)}

__all__ = ['Element', 'DynamicElement', 'QueryElement', 'AnyElement', 'ListElement', 'AnyMetal']
__all__.extend(k for k in globals() if k.startswith('Group'))
__all__.extend(k for k in globals() if k.startswith('Period'))
__all__.extend(elements)


for _class in (DynamicElement, QueryElement):
    for k, v in elements.items():
        name = f'{_class.__name__[:-7]}{k}'
        globals()[name] = cls = type(name, (_class, *v.__mro__[-3:-1]),
                                     {'__module__': v.__module__, '__slots__': (), 'atomic_number': v.atomic_number,
                                      'atomic_radius': v.atomic_radius})
        setattr(modules[v.__module__], name, cls)
        modules[v.__module__].__all__.append(name)
        __all__.append(name)
