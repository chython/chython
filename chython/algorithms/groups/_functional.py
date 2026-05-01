# -*- coding: utf-8 -*-
#
#  Copyright 2024-2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from lazy_object_proxy import Proxy


def _rules():
    from ... import smarts

    rules = {}

    rules['aryl_fluoride'] = smarts('[F;D1:1]-[C;a:2]')
    rules['aryl_chloride'] = smarts('[Cl;D1:1]-[C;a:2]')
    rules['aryl_bromide'] = smarts('[Br;D1:1]-[C;a:2]')
    rules['aryl_iodide'] = smarts('[I;D1:1]-[C;a:2]')
    rules['aryl_triflate'] = smarts('[S;D4](=O)(=O)(-[O:1]-;!@[C;a:2])-[C;D4](F)(F)F')
    rules['aryl_tosylate'] = smarts('[S;D4](=O)(=O)(-[O:1]-;!@[C;a:2])-[C;a]:1:[C;D2]:[C;D2]:[C](-[C;D1]):[C;D2]:[C;D2]:1')
    rules['aryl_mesylate'] = smarts('[S;D4](=O)(=O)(-[O:1]-;!@[C;a:2])-[C;D1]')
    rules['phenol'] = smarts('[O;D1;z1;x0:1]-[C;a:2]')

    return rules


rules = Proxy(_rules)


__all__ = ['rules']
