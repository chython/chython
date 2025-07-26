# -*- coding: utf-8 -*-
#
#  Copyright 2020-2025 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from .vector import Vector
from ...exceptions import ValenceError


class Element(ABC):
    __slots__ = ('_isotope', '_charge', '_is_radical', '_xy', '_implicit_hydrogens',
                 '_explicit_hydrogens', '_stereo', '_parsed_mapping',
                 '_neighbors', '_heteroatoms', '_hybridization', '_ring_sizes', '_in_ring', '_extended_stereo')
    __class_cache__ = {}

    def __init__(self, isotope: Optional[int] = None, *,
                 charge: int = 0, is_radical: bool = False, x: float = 0., y: float = 0.,
                 implicit_hydrogens: Optional[int] = None, stereo: Optional[bool] = None,
                 parsed_mapping: Optional[int] = None, delta_isotope: Optional[int] = None,
                 extended_stereo: Optional[int] = None):
        """
        Element object with specified isotope

        :param isotope: Isotope number of element
        """
        if delta_isotope is not None:
            assert isotope is None, 'isotope absolute value and delta value provided'
            isotope = self.mdl_isotope + delta_isotope

        self.isotope = isotope
        self.charge = charge
        self.is_radical = is_radical
        self._xy = Vector(x, y)

        self._implicit_hydrogens = implicit_hydrogens
        self._stereo = stereo
        self._parsed_mapping = parsed_mapping
        self._extended_stereo = extended_stereo

    def __repr__(self):
        if self.isotope:
            return f'{self.__class__.__name__}({self.isotope})'
        return f'{self.__class__.__name__}()'

    @property
    def atomic_symbol(self) -> str:
        return self.__class__.__name__

    @property
    @abstractmethod
    def atomic_number(self) -> int:
        """
        Element number
        """

    @property
    def isotope(self) -> Optional[int]:
        """
        Isotope number
        """
        return self._isotope

    @isotope.setter
    def isotope(self, value: Optional[int]):
        if isinstance(value, int):
            if value not in self.isotopes_distribution:
                raise ValueError(f'isotope number {value} impossible or not stable for {self.atomic_symbol}')
        elif value is not None:
            raise TypeError('integer isotope number required')
        self._isotope = value

    @property
    def atomic_mass(self) -> float:
        mass = self.isotopes_masses
        if self.isotope is None:
            return sum(x * mass[i] for i, x in self.isotopes_distribution.items())
        return mass[self.isotope]

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

    @property
    @abstractmethod
    def mdl_isotope(self) -> int:
        """
        MDL MOL common isotope
        """

    @property
    def is_forming_single_bonds(self) -> bool:
        """
        Atom can form stable covalent single bonds in molecules
        """
        return False

    @property
    def is_forming_double_bonds(self) -> bool:
        """
        Atom can form stable covalent double bonds in molecules
        """
        return False

    @property
    def charge(self) -> int:
        """
        Charge of atom
        """
        return self._charge

    @charge.setter
    def charge(self, value: int):
        """
        Update charge of atom. Make sure to flush cache and recalculate hydrogens count and stereo.
        Or use context manager on molecule:

            with mol:
                mol.atom(1).charge = 1
        """
        if not isinstance(value, int):
            raise TypeError('formal charge should be int in range [-4, 4]')
        elif value > 4 or value < -4:
            raise ValueError('formal charge should be in range [-4, 4]')
        self._charge = value

    @property
    def is_radical(self) -> bool:
        """
        Radical state of atoms
        """
        return self._is_radical

    @is_radical.setter
    def is_radical(self, value: bool):
        """
        Update radical state of atom. Make sure to flush cache and recalculate hydrogens count and stereo.
        Or use context manager on molecule:

            with mol:
                mol.atom(1).is_radical = True
        """
        if not isinstance(value, bool):
            raise TypeError('bool expected')
        self._is_radical = value

    @property
    def x(self) -> float:
        """
        X coordinate of atom on 2D plane
        """
        return self._xy.x

    @x.setter
    def x(self, value: float):
        self._xy.x = value

    @property
    def y(self) -> float:
        """
        Y coordinate of atom on 2D plane
        """
        return self._xy.y

    @y.setter
    def y(self, value: float):
        self._xy.y = value

    @property
    def xy(self) -> Vector:
        """
        (X, Y) coordinates of atom on 2D plane
        """
        return self._xy

    @xy.setter
    def xy(self, value: Tuple[float, float]):
        self._xy = Vector(*value)

    @property
    def implicit_hydrogens(self) -> Optional[int]:
        return self._implicit_hydrogens

    @property
    def explicit_hydrogens(self) -> int:
        return self._explicit_hydrogens

    @property
    def total_hydrogens(self) -> int:
        if self.implicit_hydrogens is None:
            raise ValenceError
        return self.implicit_hydrogens + self.explicit_hydrogens

    @property
    def stereo(self) -> Optional[bool]:
        """
        Tetrahedron or allene stereo label
        """
        return self._stereo

    @property
    def extended_stereo(self) -> Optional[int]:
        """
        Positive int for an AND group id, negative int for an OR group id
        """
        if self._stereo is None:
            return None
        return self._extended_stereo

    @property
    def heteroatoms(self) -> int:
        return self._heteroatoms

    @property
    def neighbors(self) -> int:
        """
        Neighbors count of atom
        """
        return self._neighbors

    @property
    def hybridization(self):
        """
        1 - if atom has zero or only single bonded neighbors, 2 - if has only one double bonded neighbor and any amount
        of single bonded, 3 - if has one triple bonded and any amount of double and single bonded neighbors or
        two double bonded and any amount of single bonded neighbors, 4 - if atom in aromatic ring.
        """
        return self._hybridization

    @property
    def ring_sizes(self) -> Set[int]:
        """
        Atom rings sizes.
        """
        return self._ring_sizes

    @property
    def in_ring(self) -> bool:
        """
        Atom in any ring.
        """
        return self._in_ring

    def copy(self, full=False, hydrogens=False, stereo=False) -> 'Element':
        """
        Get a copy of the Element object with attribute copy control.
        """
        copy = object.__new__(self.__class__)
        copy._isotope = self.isotope
        copy._charge = self.charge
        copy._is_radical = self.is_radical
        copy._xy = self.xy
        if full:
            copy._implicit_hydrogens = self.implicit_hydrogens
            copy._stereo = self.stereo
            copy._extended_stereo = self.extended_stereo
            copy._explicit_hydrogens = self.explicit_hydrogens
            copy._neighbors = self.neighbors
            copy._heteroatoms = self.heteroatoms
            copy._hybridization = self.hybridization
            copy._ring_sizes = self.ring_sizes.copy()
            copy._in_ring = self.in_ring
        else:
            if hydrogens:
                copy._implicit_hydrogens = self.implicit_hydrogens
            else:
                copy._implicit_hydrogens = None
            if stereo:
                copy._stereo = self.stereo
                copy._extended_stereo = self.extended_stereo
            else:
                copy._stereo = copy._extended_stereo = None
        return copy

    def __copy__(self):
        return self.copy()

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

    def __eq__(self, other):
        """
        compare attached to molecules elements
        """
        if isinstance(other, int):
            return self.atomic_number == other
        elif isinstance(other, str):
            return self.atomic_symbol == other
        return isinstance(other, Element) and self.atomic_number == other.atomic_number and \
            self.isotope == other.isotope and self.charge == other.charge and self.is_radical == other.is_radical

    def __hash__(self):
        return hash((self.isotope or 0, self.atomic_number, self.charge, self.is_radical,
                     self.implicit_hydrogens or 0, self.in_ring))

    def valence_rules(self, valence: int) -> \
            List[Tuple[Set[Tuple[int, 'Element']], Dict[Tuple[int, 'Element'], int], int]]:
        """
        valence rules for element with specific charge/radical state
        """
        try:
            return self._compiled_valence_rules[(self.charge, self.is_radical, valence)]
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
    def _valences_exceptions(self) -> Tuple[Tuple[int, bool, int, Tuple[Tuple[int, str], ...]], ...]:
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
            # only first common valence represents implicit H.
            valence = self._common_valences[0]
            for h in range(valence + 1):
                rules[(0, False, valence - h)].append((set(), {}, h))  # any atoms and bonds possible

            for valence in self._common_valences[1:]:
                rules[(0, False, valence)].append((set(), {}, 0))
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

    @class_cached_property
    def _compiled_saturation_rules(self) -> List[Tuple[int, bool, int, int, Optional[Dict[Tuple[int, int], int]]]]:
        """
        dictionary with key = (charge, is_radical, sum_of_bonds) and
        value = list of possible neighbors
        """
        elements_classes = {x.__name__: x.atomic_number.fget(None) for x in Element.__subclasses__()}

        rules = []
        for valence in self._common_valences:
            rules.append((0, False, valence, 0, None))  # any atoms and bonds possible

        for charge, is_radical, implicit, environment in self._valences_exceptions:
            if not environment:
                rules.append((charge, is_radical, implicit, 0, None))
            else:
                explicit_dict = defaultdict(int)
                explicit = 0
                for b, e in environment:
                    explicit_dict[(b, elements_classes[e])] += 1
                    explicit += b
                rules.append((charge, is_radical, implicit + explicit, implicit, dict(explicit_dict)))
        return rules


__all__ = ['Element']
