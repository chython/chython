# -*- coding: utf-8 -*-
#
#  Copyright 2019-2021 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2019 Tagir Akhmetshin <tagirshin@gmail.com>
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
from .element import Element
from .groups import GroupVIII
from .periods import PeriodIV, PeriodV, PeriodVI, PeriodVII


class Fe(Element, PeriodIV, GroupVIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 26

    @property
    def isotopes_distribution(self):
        return {54: 0.05845, 55: 0., 56: 0.91754, 57: 0.02119, 58: 0.00282, 59: 0.}

    @property
    def isotopes_masses(self):
        return {54: 53.939615, 55: 54.938293, 56: 55.934942, 57: 56.935399, 58: 57.933281, 59: 58.934876}

    @property
    def _common_valences(self):
        return 0, 2, 3

    @property
    def _valences_exceptions(self):
        return (2, False, 0, ()), (3, False, 0, ()),

    @property
    def atomic_radius(self):
        return 1.56


class Ru(Element, PeriodV, GroupVIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 44

    @property
    def isotopes_distribution(self):
        return {96: 0.0554, 98: 0.0187, 99: 0.1276, 100: 0.126, 101: 0.1706, 102: 0.3155, 104: 0.1862, 106: 0.}

    @property
    def isotopes_masses(self):
        return {96: 95.907598, 98: 97.905287, 99: 98.905939, 100: 99.90422, 101: 100.905582, 102: 101.904349,
                104: 103.90543, 106: 105.907329}

    @property
    def _common_valences(self):
        return 0,

    @property
    def _valences_exceptions(self):
        return (0, False, 0, ((2, 'O'), (2, 'O'), (2, 'O'), (2, 'O'))),  # RuO4

    @property
    def atomic_radius(self):
        return 1.78


class Os(Element, PeriodVI, GroupVIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 76

    @property
    def isotopes_distribution(self):
        return {184: 0.0002, 186: 0.0159, 187: 0.0196, 188: 0.1324, 189: 0.1615, 190: 0.2626, 191: 0., 192: 0.4078}

    @property
    def isotopes_masses(self):
        return {184: 183.952491, 186: 185.953838, 187: 186.955748, 188: 187.955836, 189: 188.958145, 190: 189.958445,
                191: 190.960930, 192: 191.961479}

    @property
    def _common_valences(self):
        return 0,

    @property
    def _valences_exceptions(self):
        return ((0, False, 0, ((2, 'O'), (2, 'O'), (2, 'O'), (2, 'O'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (2, 'O'), (1, 'O'), (1, 'O'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (2, 'O'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'O'), (1, 'O'))))

    @property
    def atomic_radius(self):
        return 1.85


class Hs(Element, PeriodVII, GroupVIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 108

    @property
    def isotopes_distribution(self):
        return {240: 1.0}

    @property
    def isotopes_masses(self):
        return {270: 270.134293}

    @property
    def _common_valences(self):
        return 0,

    @property
    def _valences_exceptions(self):
        return ()

    @property
    def atomic_radius(self):
        return 1.85  # unknown, taken radius of previous element in group


__all__ = ['Fe', 'Ru', 'Os', 'Hs']
