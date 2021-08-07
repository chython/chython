# -*- coding: utf-8 -*-
#
#  Copyright 2019-2021 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2019 Tagir Akhmetshin <tagirshin@gmail.com>
#  Copyright 2019 Dayana Bashirova <dayana.bashirova@yandex.ru>
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
from .groups import GroupX
from .periods import PeriodIV, PeriodV, PeriodVI, PeriodVII


class Ni(Element, PeriodIV, GroupX):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 28

    @property
    def isotopes_distribution(self):
        return {58: 0.680769, 60: 0.262231, 61: 0.011399, 62: 0.036345, 63: 0., 64: 0.009256}

    @property
    def isotopes_masses(self):
        return {58: 57.935348, 60: 59.930791, 61: 60.931060, 62: 61.928349, 63: 62.929669, 64: 63.927970}

    @property
    def _common_valences(self):
        return 0, 2

    @property
    def _valences_exceptions(self):
        return ((2, False, 0, ()),
                (0, False, 0, ((2, 'O'), (1, 'O'))))  # Ni2O3

    @property
    def atomic_radius(self):
        return 1.49


class Pd(Element, PeriodV, GroupX):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 46

    @property
    def isotopes_distribution(self):
        return {102: 0.0102, 103: 0., 104: 0.1114, 105: 0.2233, 106: 0.2733, 108: 0.2646, 109: 0., 110: 0.1172}

    @property
    def isotopes_masses(self):
        return {102: 101.905608, 103: 102.906087, 104: 103.904035, 105: 104.905084, 106: 105.903483, 108: 107.903894,
                109: 108.905950, 110: 109.905152}

    @property
    def _common_valences(self):
        return 0, 2

    @property
    def _valences_exceptions(self):
        return ((2, False, 0, ()),
                (-2, False, 0, ((1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'))),  # [Pd(OH)4]2-
                (-2, False, 0, ((1, 'C'), (1, 'C'), (1, 'C'), (1, 'C'))),  # [Pd(CN)4]2-
                (-2, False, 0, ((1, 'N'), (1, 'N'), (1, 'N'), (1, 'N'))),  # [Pd(NCS)4]2-
                (-2, False, 0, ((1, 'S'), (1, 'S'), (1, 'S'), (1, 'S'))),  # [Pd(SCN)4]2-
                (-2, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),  # [PdF4]2-
                (-2, False, 0, ((1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'))))  # [PdCl4]2-

    @property
    def atomic_radius(self):
        return 1.69


class Pt(Element, PeriodVI, GroupX):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 78

    @property
    def isotopes_distribution(self):
        return {190: 0.00014, 192: 0.00782, 194: 0.32967, 195: 0.33832, 196: 0.25242, 198: 0.07163}

    @property
    def isotopes_masses(self):
        return {190: 189.95993, 192: 191.961035, 194: 193.962664, 195: 194.964774, 196: 195.964935, 198: 197.967876}

    @property
    def _common_valences(self):
        return 0, 2

    @property
    def _valences_exceptions(self):
        return ((0, False, 0, ((1, 'N'), (1, 'N'), (1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((1, 'N'), (1, 'N'), (1, 'O'), (1, 'O'))),
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),  # PtF6
                (0, False, 0, ((2, 'O'), (2, 'O'), (2, 'O'))))  # PtO3

    @property
    def atomic_radius(self):
        return 1.77


class Ds(Element, PeriodVII, GroupX):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 110

    @property
    def isotopes_distribution(self):
        return {281: 1.0}

    @property
    def isotopes_masses(self):
        return {281: 281.164516}

    @property
    def _common_valences(self):
        return 0,

    @property
    def _valences_exceptions(self):
        return ()

    @property
    def atomic_radius(self):
        return 1.77  # unknown, taken radius of previous element in group


__all__ = ['Ni', 'Pd', 'Pt', 'Ds']
