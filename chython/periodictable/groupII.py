# -*- coding: utf-8 -*-
#
#  Copyright 2019-2024 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from .base import Element
from .base.groups import GroupII
from .base.periods import PeriodII, PeriodIII, PeriodIV, PeriodV, PeriodVI, PeriodVII


class Be(Element, PeriodII, GroupII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 4

    @property
    def isotopes_distribution(self):
        return {9: 1.0}

    @property
    def isotopes_masses(self):
        return {9: 9.012182}

    @property
    def _common_valences(self):
        return 0, 2

    @property
    def _valences_exceptions(self):
        return (2, False, 0, ()),

    @property
    def atomic_radius(self):
        return 1.12

    @property
    def mdl_isotope(self):
        return 9


class Mg(Element, PeriodIII, GroupII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 12

    @property
    def isotopes_distribution(self):
        return {24: 0.7899, 25: 0.1, 26: 0.1101}

    @property
    def isotopes_masses(self):
        return {24: 23.985042, 25: 24.985837, 26: 25.982593}

    @property
    def _common_valences(self):
        return 0, 2

    @property
    def _valences_exceptions(self):
        return ((2, False, 0, ()),
                (1, False, 0, ((1, 'C'),)),
                (1, False, 0, ((1, 'O'),)),
                (1, False, 0, ((1, 'Br'),)),
                (1, False, 0, ((1, 'Cl'),)))

    @property
    def atomic_radius(self):
        return 1.45

    @property
    def mdl_isotope(self):
        return 24


class Ca(Element, PeriodIV, GroupII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 20

    @property
    def isotopes_distribution(self):
        return {40: 0.96941, 42: 0.00647, 43: 0.00135, 44: 0.02086, 45: 0., 46: 4e-05, 47: 0., 48: 0.00187}

    @property
    def isotopes_masses(self):
        return {40: 39.962591, 42: 41.958618, 43: 42.958767, 44: 43.955481, 45: 44.956186, 46: 45.953693, 47: 46.954541,
                48: 47.952534}

    @property
    def _common_valences(self):
        return 0, 2

    @property
    def _valences_exceptions(self):
        return (2, False, 0, ()),

    @property
    def atomic_radius(self):
        return 1.94

    @property
    def mdl_isotope(self):
        return 40


class Sr(Element, PeriodV, GroupII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 38

    @property
    def isotopes_distribution(self):
        return {84: 0.0056, 85: 0., 86: 0.0986, 87: 0.07, 88: 0.8258, 89: 0.}

    @property
    def isotopes_masses(self):
        return {84: 83.913425, 85: 84.912933, 86: 85.909262, 87: 86.908879, 88: 87.905614, 89: 88.907451}

    @property
    def _common_valences(self):
        return 0, 2

    @property
    def _valences_exceptions(self):
        return (2, False, 0, ()),

    @property
    def atomic_radius(self):
        return 2.19

    @property
    def mdl_isotope(self):
        return 88


class Ba(Element, PeriodVI, GroupII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 56

    @property
    def isotopes_distribution(self):
        return {130: 0.00106, 132: 0.00101, 134: 0.02417, 135: 0.06592, 136: 0.07854, 137: 0.11232, 138: 0.71698}

    @property
    def isotopes_masses(self):
        return {130: 129.90631, 132: 131.905056, 134: 133.904503, 135: 134.905683, 136: 135.90457, 137: 136.905821,
                138: 137.905241}

    @property
    def _common_valences(self):
        return 0, 2

    @property
    def _valences_exceptions(self):
        return (2, False, 0, ()),

    @property
    def atomic_radius(self):
        return 2.53

    @property
    def mdl_isotope(self):
        return 137


class Ra(Element, PeriodVII, GroupII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 88

    @property
    def isotopes_distribution(self):
        return {223: 0., 226: 1.0, 228: 0., 233: 0.}

    @property
    def isotopes_masses(self):
        return {223: 223.018502, 226: 226.025410, 228: 228.031070, 233: 233.048065}

    @property
    def _common_valences(self):
        return 0, 2

    @property
    def _valences_exceptions(self):
        return (2, False, 0, ()),

    @property
    def atomic_radius(self):
        return 2.53  # unknown, taken radius of previous element in group

    @property
    def mdl_isotope(self):
        return 226


__all__ = ['Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra']
