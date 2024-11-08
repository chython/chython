# -*- coding: utf-8 -*-
#
#  Copyright 2019-2024 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2019 Tagir Akhmetshin <tagirshin@gmail.com>
#  Copyright 2019 Tansu Nasyrova <tansu.nasurova@gmail.com>
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
from .base import Element
from .base.groups import GroupIX
from .base.periods import PeriodIV, PeriodV, PeriodVI, PeriodVII


class Co(Element, PeriodIV, GroupIX):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 27

    @property
    def isotopes_distribution(self):
        return {55: 0., 57: 0., 58: 0., 59: 1.0, 60: 0.}

    @property
    def isotopes_masses(self):
        return {55: 54.941999, 57: 56.936291, 58: 57.935753, 59: 58.933200, 60: 59.933817}

    @property
    def _common_valences(self):
        return 0, 2, 3

    @property
    def _valences_exceptions(self):
        return ((2, False, 0, ()),
                (3, False, 0, ()),
                (0, False, 0, ((1, 'H'),)),  # HCo(CO)n
                (2, False, 0, ((1, 'N'),)),  # B12

                (-3, False, 0, ((2, 'O'), (1, 'O'), (1, 'O'), (1, 'O'))),  # [CoO4]3-
                (-2, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),  # [CoF6]2-

                (-1, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'))),  # [CoF3]-
                (-1, False, 0, ((1, 'Cl'), (1, 'Cl'), (1, 'Cl'))),
                (-1, False, 0, ((1, 'O'), (1, 'O'), (1, 'O'))),  # [Co(OH)3]-

                (-2, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),  # [CoF4]2-
                (-2, False, 0, ((1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'))),
                (-2, False, 0, ((1, 'Br'), (1, 'Br'), (1, 'Br'), (1, 'Br'))),
                (-2, False, 0, ((1, 'I'), (1, 'I'), (1, 'I'), (1, 'I'))),
                (-2, False, 0, ((1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'))),  # [Co(OH)4]2-

                (-3, False, 0, ((1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'))),  # [CoCl5]3-

                (-3, False, 0, ((1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'))),  # [Co(OH)6]3-
                (-4, False, 0, ((1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'))))  # [Co(OH)6]4-

    @property
    def atomic_radius(self):
        return 1.52

    @property
    def mdl_isotope(self):
        return 59


class Rh(Element, PeriodV, GroupIX):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 45

    @property
    def isotopes_distribution(self):
        return {103: 1.0, 105: 0.}

    @property
    def isotopes_masses(self):
        return {103: 102.905504, 105: 104.905694}

    @property
    def _common_valences(self):
        return 0, 3, 4

    @property
    def _valences_exceptions(self):
        return ((0, False, 0, ((2, 'O'),)),  # RhO
                (0, False, 0, ((1, 'O'), (1, 'O'))),  # Rh(OH)2
                (0, False, 0, ((2, 'S'),)),
                (0, False, 0, ((1, 'S'), (1, 'S'))),
                (-1, False, 0, ((1, 'Br'), (1, 'Br'), (1, 'Br'), (1, 'Br'))),  # [RhBr4]-
                (-3, False, 0, ((1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'))),  # [RhCl6]3-
                (-3, False, 0, ((1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'))),  # [Rh(NO2)6]3-
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'H'),)),  # HRh(CO)4, HRh(CO)[P(Ph)3]3
                (0, False, 0, ((1, 'Cl'),)))  # Rh2Cl2

    @property
    def atomic_radius(self):
        return 1.73

    @property
    def mdl_isotope(self):
        return 103


class Ir(Element, PeriodVI, GroupIX):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 77

    @property
    def isotopes_distribution(self):
        return {191: 0.373, 192: 0., 193: 0.627}

    @property
    def isotopes_masses(self):
        return {191: 190.960591, 192: 191.962605, 193: 192.962924}

    @property
    def _common_valences(self):
        return 0, 3, 4

    @property
    def _valences_exceptions(self):
        return ((0, False, 0, ((1, 'O'),)),
                (0, False, 0, ((1, 'F'),)),
                (0, False, 0, ((1, 'Cl'),)),
                (0, False, 0, ((1, 'Br'),)),
                (0, False, 0, ((1, 'I'),)),
                (0, False, 0, ((1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((1, 'Br'), (1, 'Br'))),
                (0, False, 0, ((1, 'I'), (1, 'I'))),
                (0, False, 0, ((1, 'S'), (1, 'S'))),
                (0, False, 0, ((2, 'S'),)),
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),
                (-3, False, 0, ((1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'))))

    @property
    def atomic_radius(self):
        return 1.8

    @property
    def mdl_isotope(self):
        return 192


class Mt(Element, PeriodVII, GroupIX):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 109

    @property
    def isotopes_distribution(self):
        return {278: 1.0}

    @property
    def isotopes_masses(self):
        return {278: 278.15481}

    @property
    def _common_valences(self):
        return 0,

    @property
    def _valences_exceptions(self):
        return ()

    @property
    def atomic_radius(self):
        return 1.8  # unknown, taken radius of previous element in group

    @property
    def mdl_isotope(self):
        return 278


__all__ = ['Co', 'Rh', 'Ir', 'Mt']
