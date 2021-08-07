# -*- coding: utf-8 -*-
#
#  Copyright 2019-2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from .groups import GroupI
from .periods import *


class H(Element, PeriodI, GroupI):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 1

    @property
    def isotopes_distribution(self):
        return {1: 0.999885, 2: 0.000115, 3: 0.}

    @property
    def isotopes_masses(self):
        return {1: 1.007825, 2: 2.014102, 3: 3.016049}

    @property
    def _common_valences(self):
        return 1,

    @property
    def _valences_exceptions(self):
        return (1, False, 0, ()), (0, True, 0, ()), (-1, False, 0, ())

    @property
    def atomic_radius(self):
        return 0.53


class Li(Element, PeriodII, GroupI):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 3

    @property
    def isotopes_distribution(self):
        return {6: 0.0759, 7: 0.9241}

    @property
    def isotopes_masses(self):
        return {6: 6.015122, 7: 7.016004}

    @property
    def _common_valences(self):
        return 0, 1

    @property
    def _valences_exceptions(self):
        return (1, False, 0, ()),

    @property
    def atomic_radius(self):
        return 167


class Na(Element, PeriodIII, GroupI):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 11

    @property
    def isotopes_distribution(self):
        return {22: 0., 23: 1.0}

    @property
    def isotopes_masses(self):
        return {22: 21.994437, 23: 22.98977}

    @property
    def _common_valences(self):
        return 0, 1

    @property
    def _valences_exceptions(self):
        return (1, False, 0, ()),

    @property
    def atomic_radius(self):
        return 1.9


class K(Element, PeriodIV, GroupI):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 19

    @property
    def isotopes_distribution(self):
        return {39: 0.932581, 40: 0.000117, 41: 0.067302, 42: 0.}

    @property
    def isotopes_masses(self):
        return {39: 38.963707, 40: 39.963999, 41: 40.961826, 42: 41.962402}

    @property
    def _common_valences(self):
        return 0, 1

    @property
    def _valences_exceptions(self):
        return (1, False, 0, ()),

    @property
    def atomic_radius(self):
        return 2.43


class Rb(Element, PeriodV, GroupI):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 37

    @property
    def isotopes_distribution(self):
        return {82: 0., 85: 0.7217, 87: 0.2783}

    @property
    def isotopes_masses(self):
        return {82: 81.918209, 85: 84.911789, 87: 86.909183}

    @property
    def _common_valences(self):
        return 0, 1

    @property
    def _valences_exceptions(self):
        return (1, False, 0, ()),

    @property
    def atomic_radius(self):
        return 2.65


class Cs(Element, PeriodVI, GroupI):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 55

    @property
    def isotopes_distribution(self):
        return {131: 0., 133: 1.0}

    @property
    def isotopes_masses(self):
        return {131: 130.905464, 133: 132.905447}

    @property
    def _common_valences(self):
        return 0, 1

    @property
    def _valences_exceptions(self):
        return (1, False, 0, ()),

    @property
    def atomic_radius(self):
        return 2.98


class Fr(Element, PeriodVII, GroupI):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 87

    @property
    def isotopes_distribution(self):
        return {223: 1.0}

    @property
    def isotopes_masses(self):
        return {223: 223.019736}

    @property
    def _common_valences(self):
        return 0, 1

    @property
    def _valences_exceptions(self):
        return (1, False, 0, ()),

    @property
    def atomic_radius(self):
        return 2.98  # unknown, taken radius of previous element in group


__all__ = ['H', 'Li', 'Na', 'K', 'Rb', 'Cs', 'Fr']
