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
from .groups import GroupIV
from .periods import PeriodIV, PeriodV, PeriodVI, PeriodVII


class Ti(Element, PeriodIV, GroupIV):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 22

    @property
    def isotopes_distribution(self):
        return {46: 0.0825, 47: 0.0744, 48: 0.7372, 49: 0.0541, 50: 0.0518}

    @property
    def isotopes_masses(self):
        return {46: 45.95263, 47: 46.951764, 48: 47.947947, 49: 48.947871, 50: 49.944792}

    @property
    def _common_valences(self):
        return 0, 4

    @property
    def _valences_exceptions(self):
        return ((-2, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),  # [TiF6]2-
                (-2, False, 0, ((1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'))),  # [TiCl6]2-
                (-2, False, 0, ((1, 'Br'), (1, 'Br'), (1, 'Br'), (1, 'Br'), (1, 'Br'), (1, 'Br'))),  # [TiBr6]2-
                (-2, False, 0, ((1, 'I'), (1, 'I'), (1, 'I'), (1, 'I'), (1, 'I'), (1, 'I'))),  # [TiI6]2-
                (-2, False, 0, ((1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'))),  # [Ti(SO4)3]2-
                (-2, False, 0, ((1, 'S'), (1, 'S'), (1, 'S'), (1, 'S'), (1, 'S'), (1, 'S'))),  # [Ti(SCN)6]2-
                (-2, False, 0, ((1, 'N'), (1, 'N'), (1, 'N'), (1, 'N'), (1, 'N'), (1, 'N'))),  # [Ti(NCS)6]2-
                (-2, False, 0, ((1, 'C'), (1, 'C'), (1, 'C'), (1, 'C'), (1, 'C'), (1, 'C'))),  # [Ti(CN)6]2-

                (-2, False, 0, ((2, 'O'), (2, 'O'), (2, 'O'))),  # [TiO3]2-
                (-2, False, 0, ((2, 'O'), (1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'))),  # [TiO(SO4)2]2-
                (-2, False, 0, ((2, 'O'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),  # [TiOF4]2-
                (-2, False, 0, ((2, 'O'), (1, 'S'), (1, 'S'), (1, 'S'), (1, 'S'))),  # [TiO(SCN)4]2-
                (-2, False, 0, ((2, 'O'), (1, 'N'), (1, 'N'), (1, 'N'), (1, 'N'))),  # [TiO(NCS)4]2-

                (0, False, 0, ((1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((1, 'Br'), (1, 'Br'))),
                (0, False, 0, ((1, 'I'), (1, 'I'))),
                (0, False, 0, ((1, 'H'), (1, 'H'))),
                (0, False, 0, ((1, 'O'), (1, 'O'))),
                (0, False, 0, ((2, 'O'),)),

                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'Cl'), (1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((1, 'Br'), (1, 'Br'), (1, 'Br'))),
                (0, False, 0, ((1, 'I'), (1, 'I'), (1, 'I'))),
                (0, False, 0, ((1, 'O'), (1, 'O'), (1, 'O'))),
                (0, False, 0, ((2, 'O'), (1, 'O'))),
                (0, False, 0, ((1, 'N'), (1, 'N'), (1, 'N'))),  # TiN
                (0, False, 0, ((2, 'N'), (1, 'N'))),
                (0, False, 0, ((3, 'N'),)))

    @property
    def atomic_radius(self):
        return 1.76


class Zr(Element, PeriodV, GroupIV):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 40

    @property
    def isotopes_distribution(self):
        return {89: 0., 90: 0.5145, 91: 0.1122, 92: 0.1715, 94: 0.1738, 96: 0.028}

    @property
    def isotopes_masses(self):
        return {89: 88.908890, 90: 89.904704, 91: 90.905645, 92: 91.905040, 94: 93.906316, 96: 95.908276}

    @property
    def _common_valences(self):
        return 0, 4

    @property
    def _valences_exceptions(self):
        return ((0, False, 0, ((1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((1, 'Br'), (1, 'Br'))),
                (0, False, 0, ((1, 'I'), (1, 'I'))),
                (0, False, 0, ((1, 'H'), (1, 'H'))),

                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'Cl'), (1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((1, 'Br'), (1, 'Br'), (1, 'Br'))),
                (0, False, 0, ((1, 'I'), (1, 'I'), (1, 'I'))),

                (0, False, 0, ((1, 'N'), (1, 'N'), (1, 'N'))),  # ZrN
                (0, False, 0, ((2, 'N'), (1, 'N'))),
                (0, False, 0, ((3, 'N'),)),
                (0, False, 0, ((1, 'N'), (1, 'N'))),  # Zr2N3
                (0, False, 0, ((2, 'N'),)))

    @property
    def atomic_radius(self):
        return 2.06


class Hf(Element, PeriodVI, GroupIV):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 72

    @property
    def isotopes_distribution(self):
        return {174: 0.0016, 176: 0.0526, 177: 0.186, 178: 0.2728, 179: 0.1362, 180: 0.3508}

    @property
    def isotopes_masses(self):
        return {174: 173.94004, 176: 175.941402, 177: 176.94322, 178: 177.943698, 179: 178.945815, 180: 179.946549}

    @property
    def _common_valences(self):
        return 0, 4

    @property
    def _valences_exceptions(self):
        return ((0, False, 0, ((1, 'Cl'),)),
                (0, False, 0, ((1, 'Br'), (1, 'Br'))),
                (0, False, 0, ((1, 'Cl'), (1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((1, 'Br'), (1, 'Br'), (1, 'Br'))),
                (0, False, 0, ((1, 'I'), (1, 'I'), (1, 'I'))),
                (0, False, 0, ((1, 'N'), (1, 'N'), (1, 'N'))),  # HfN
                (0, False, 0, ((2, 'N'), (1, 'N'))),
                (0, False, 0, ((3, 'N'),)))

    @property
    def atomic_radius(self):
        return 2.08


class Rf(Element, PeriodVII, GroupIV):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 104

    @property
    def isotopes_distribution(self):
        return {267: 1.0}

    @property
    def isotopes_masses(self):
        return {267: 267.12153}

    @property
    def _common_valences(self):
        return 0, 4

    @property
    def _valences_exceptions(self):
        return (4, False, 0, ()),

    @property
    def atomic_radius(self):
        return 2.08  # unknown, taken radius of previous element in group


__all__ = ['Ti', 'Zr', 'Hf', 'Rf']
