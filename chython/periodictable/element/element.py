# -*- coding: utf-8 -*-
#
#  Copyright 2020, 2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from abc import ABC, abstractmethod
from CachedMethods import class_cached_property
from collections import defaultdict
from typing import Dict, List, Optional, Set, Tuple, Type
from .core import Core
from ...exceptions import IsNotConnectedAtom, ValenceError


class Element(Core, ABC):
    __slots__ = ()
    __class_cache__ = {}

    def __init__(self, isotope: Optional[int] = None):
        """
        Element object with specified isotope

        :param isotope: Isotope number of element
        """
        if isinstance(isotope, int):
            if isotope not in self.isotopes_distribution:
                raise ValueError(f'isotope number {isotope} impossible or not stable for {self.atomic_symbol}')
        elif isotope is not None:
            raise TypeError('integer isotope number required')
        super().__init__(isotope)

    @property
    def atomic_symbol(self) -> str:
        return self.__class__.__name__

    @property
    def atomic_mass(self) -> float:
        mass = self.isotopes_masses
        if self.__isotope is None:
            return sum(x * mass[i] for i, x in self.isotopes_distribution.items())
        return mass[self.__isotope]

    @property
    @abstractmethod
    def isotopes_distribution(self) -> Dict[int, float]:
        """
        Isotopes distribution in earth
        """

    @property
    @abstractmethod
    def isotopes_masses(self) -> Dict[int, float]:
        """
        Isotopes masses
        """

    @property
    @abstractmethod
    def atomic_radius(self) -> float:
        """
        Valence radius of atom
        """

    @Core.charge.setter
    def charge(self, charge: int):
        if not isinstance(charge, int):
            raise TypeError('formal charge should be int in range [-4, 4]')
        elif charge > 4 or charge < -4:
            raise ValueError('formal charge should be in range [-4, 4]')
        try:
            g = self._graph()
            g._charges[self._map] = charge
        except AttributeError:
            raise IsNotConnectedAtom
        else:
            g._calc_implicit(self._map)
            g.flush_cache()
            g.fix_stereo()

    @Core.is_radical.setter
    def is_radical(self, is_radical: bool):
        if not isinstance(is_radical, bool):
            raise TypeError('bool expected')
        try:
            g = self._graph()
            g._radicals[self._map] = is_radical
        except AttributeError:
            raise IsNotConnectedAtom
        else:
            g._calc_implicit(self._map)
            g.flush_cache()
            g.fix_stereo()

    @property
    def implicit_hydrogens(self) -> Optional[int]:
        try:
            return self._graph()._hydrogens[self._map]
        except AttributeError:
            raise IsNotConnectedAtom

    @property
    def explicit_hydrogens(self) -> int:
        try:
            return self._graph().explicit_hydrogens(self._map)
        except AttributeError:
            raise IsNotConnectedAtom

    @property
    def total_hydrogens(self) -> int:
        try:
            return self._graph().total_hydrogens(self._map)
        except AttributeError:
            raise IsNotConnectedAtom

    @property
    def heteroatoms(self) -> int:
        try:
            return self._graph().heteroatoms(self._map)
        except AttributeError:
            raise IsNotConnectedAtom

    @property
    def neighbors(self) -> int:
        """
        Neighbors count of atom
        """
        try:
            return self._graph().neighbors(self._map)
        except AttributeError:
            raise IsNotConnectedAtom

    @property
    def hybridization(self):
        """
        1 - if atom has zero or only single bonded neighbors, 2 - if has only one double bonded neighbor and any amount
        of single bonded, 3 - if has one triple bonded and any amount of double and single bonded neighbors or
        two double bonded and any amount of single bonded neighbors, 4 - if atom in aromatic ring.
        """
        try:
            return self._graph().hybridization(self._map)
        except AttributeError:
            raise IsNotConnectedAtom

    @property
    def ring_sizes(self) -> Tuple[int, ...]:
        """
        Atom rings sizes.
        """
        try:
            return self._graph().atoms_rings_sizes[self._map]
        except AttributeError:
            raise IsNotConnectedAtom
        except KeyError:
            return ()

    @property
    def in_ring(self) -> bool:
        """
        Atom in any ring.
        """
        try:
            return self._map in self._graph().ring_atoms
        except AttributeError:
            raise IsNotConnectedAtom

    @classmethod
    def from_symbol(cls, symbol: str) -> Type['Element']:
        """
        get Element class by its symbol
        """
        try:
            element = next(x for x in Element.__subclasses__() if x.__name__ == symbol)
        except StopIteration:
            raise ValueError(f'Element with symbol "{symbol}" not found')
        return element

    @classmethod
    def from_atomic_number(cls, number: int) -> Type['Element']:
        """
        get Element class by its number
        """
        try:
            elements = cls.__class_cache__['elements']
        except KeyError:
            elements = {x.atomic_number.fget(None): x for x in Element.__subclasses__()}
            cls.__class_cache__['elements'] = elements
        try:
            return elements[number]
        except KeyError:
            raise ValueError(f'Element with number "{number}" not found')

    @classmethod
    def from_atom(cls, atom: 'Element') -> 'Element':
        """
        get Element copy
        """
        if not isinstance(atom, Element):
            raise TypeError('Element expected')
        return atom.copy()

    def __eq__(self, other):
        """
        compare attached to molecules elements
        """
        return isinstance(other, Element) and self.atomic_number == other.atomic_number and \
            self.isotope == other.isotope and self.charge == other.charge and self.is_radical == other.is_radical and \
            self.implicit_hydrogens == other.implicit_hydrogens

    def __hash__(self):
        return hash((self.isotope or 0, self.atomic_number, self.charge, self.is_radical, self.implicit_hydrogens or 0))

    def valence_rules(self, charge: int, is_radical: bool, valence: int) -> \
            List[Tuple[Set[Tuple[int, 'Element']], Dict[Tuple[int, 'Element'], int], int]]:
        """
        valence rules for element with specific charge/radical state
        """
        try:
            return self._compiled_valence_rules[(charge, is_radical, valence)]
        except KeyError:
            raise ValenceError

    @property
    @abstractmethod
    def _common_valences(self) -> Tuple[int, ...]:
        """
        common valences of element
        """

    @property
    @abstractmethod
    def _valences_exceptions(self) -> Tuple[Tuple[int, bool, int, Tuple[Tuple[int, str], ...]]]:
        """
        exceptions in charges, radical state, implicit H count, and non H neighbors of element
        examples:
        (-1, False, 1, ()) - anion, not radical, has 1 implicit hydrogen if explicit atom not exists: [OH]- or C[O-]
        (0, True, 1, ()) - neutral, radical, has 1 implicit hydrogen if explicit atom not exists: [OH]* or C[O*]
        number of free electrons calculated as diff of default valence and number of connected atoms (include implicit)
        (0, False, 1, ((1, 'C'),)) - neutral, not radical, has 1 implicit hydrogen and 1 single bonded carbon:
        CO - alcohol or e.g. COC - ether
        (0, False, 0, ((1, 'O'), (2, 'O'))) - can be anion/cation or neutral: HNO2 or [NO2]-.
        state of neighbors atoms don't take into account. order of neighbors atoms don't take into account.
        use both: (0, False, 0, ((2, 'O'),)) and (0, False, 0, ((1, 'O'), (1, 'O'))) if chains possible
        use charge transfer for carbonyles, cyanides etc:
        (-1, False, 0, ((1, 'C'),)) - [M-]-C#[O+]
        user both for cyanates and isocyanates etc complexes:
        (-1, False, 0, ((1, 'O'),)) and (-1, False, 0, ((1, 'N'),))
        """

    @class_cached_property
    def _compiled_charge_radical(self) -> Set[Tuple[int, bool]]:
        """
        exceptions in charges, radical state
        examples:
        (-1, False) - anion, not radical
        (0, True) - neutral radical
        """
        return {(c, r) for c, r, *_ in self._valences_exceptions}

    @class_cached_property
    def _compiled_valence_rules(self) -> \
            Dict[Tuple[int, bool, int], List[Tuple[Set[Tuple[int, int]], Dict[Tuple[int, int], int], int]]]:
        """
        dictionary with key = (charge, is_radical, sum_of_bonds) and
        value = list of possible neighbors and implicit H count
        """
        elements_classes = {x.__name__: x.atomic_number.fget(None) for x in Element.__subclasses__()}

        rules = defaultdict(list)
        if self._common_valences[0] and self.atomic_number != 1:  # atom has implicit hydrogens by default except H.
            prev = -1
            for valence in self._common_valences:
                for h in range(valence - prev):
                    rules[(0, False, valence - h)].append((set(), {}, h))  # any atoms and bonds possible
                prev = valence
        else:
            for valence in self._common_valences:
                rules[(0, False, valence)].append((set(), {}, 0))  # any atoms and bonds possible

        for charge, is_radical, implicit, environment in self._valences_exceptions:
            explicit = sum(x for x, _ in environment)
            explicit_dict = defaultdict(int)
            explicit_set = set()
            for b, e in environment:
                be = (b, elements_classes[e])
                explicit_set.add(be)
                explicit_dict[be] += 1
            explicit_dict = dict(explicit_dict)

            if implicit:
                valence = explicit + implicit

                for h in range(implicit + 1):
                    rules[(charge, is_radical, valence - h)].append((explicit_set, explicit_dict, h))
            else:
                rules[(charge, is_radical, explicit)].append((explicit_set, explicit_dict, 0))
        return dict(rules)


# averaged isotopes. 3.*-compatibility
isotopes = tuple(map(int, '''  1                                                                   4
                               7   9                                          11  12  14  16  19  20
                              23  24                                          27  28  31  32  35  40
                              39  40  45  48  51  52  55  56  59  59  64  65  70  73  75  79  80  84
                              85  88  89  91  93  96  98 101 103 106 108 112 115 119 122 128 127 131
                             133 137 139

                                     140 141 144 145 150 152 157 159 163 165 167 169 173 175

                                         178 181 184 186 190 192 195 197 201 204 207 209 209 210 222 
                             223 226 227

                                     232 231 238 237 244 243 247 247 251 252 257 258 259 260

                                         261 270 269 270 270 278 281 281 285 278 289 289 293 297 294
                          '''.split()))


__all__ = ['Element']
