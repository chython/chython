# -*- coding: utf-8 -*-
#
#  Copyright 2024 Timur Gimadiev <timur.gimadiev@gmail.com>
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
from abc import ABC

from .element.markushi_substituent import MarkushiElement
from .groups import Rgroup, Xgroup


class CoreR(MarkushiElement, Rgroup, ABC):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 0

    @property
    def isotopes_distribution(self):
        return {}

    @property
    def isotopes_masses(self):
        return {}

    @property
    def _common_valences(self):
        return 0,

    @property
    def _valences_exceptions(self):
        return ()

    @property
    def atomic_radius(self):
        return None


class CoreX(MarkushiElement, Xgroup, ABC):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 0

    @property
    def isotopes_distribution(self):
        return {}

    @property
    def isotopes_masses(self):
        return {}

    @property
    def _common_valences(self):
        return 0,

    @property
    def _valences_exceptions(self):
        return ()

    @property
    def atomic_radius(self):
        return None


for num in range(1, 100):  # generate 99 classes of R
    locals()[f'R{num}'] = type(f'R{num}', (CoreR, MarkushiElement,), {"__qualname__": f'MarkushiElement.R{num}'})  #

for num in range(1, 100):  # generate 99 classes of X
    locals()[f'X{num}'] = type(f'X{num}', (CoreX, MarkushiElement,), {"__qualname__": f'MarkushiElement.X{num}'})


__all__ = [*[f'R{num}' for num in range(1, 100)], *[f'R{num}' for num in range(1, 100)]]
