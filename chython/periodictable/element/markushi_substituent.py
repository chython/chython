from abc import ABC, abstractmethod
from CachedMethods import class_cached_property
from collections import defaultdict
from typing import Dict, List, Optional, Set, Tuple, Type
from .core import Core
from ...exceptions import IsNotConnectedAtom, ValenceError


class MarkushiElement(Core, ABC):
    __slots__ = ()
    __class_cache__ = {}

    def __init__(self, isotope: Optional[int] = None):
        """
        MarkushiElement object with specified isotope

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


    @classmethod
    def from_symbol(cls, symbol: str) -> Type['MarkushiElement']:
        """
        get MarkushiElement class by its symbol
        """
        try:
            element = next(x for x in MarkushiElement.__subclasses__() if x.__name__ == symbol)
        except StopIteration:
            raise ValueError(f'MarkushiElement with symbol "{symbol}" not found')
        return element

    @classmethod
    def from_atomic_number(cls, number: int) -> Type['MarkushiElement']:
        """
        get Element class by its number
        """
        try:
            elements = cls.__class_cache__['elements']
        except KeyError:
            elements = {x.atomic_number.fget(None): x for x in MarkushiElement.__subclasses__()}
            cls.__class_cache__['elements'] = elements
        try:
            return elements[number]
        except KeyError:
            raise ValueError(f'Element with number "{number}" not found')

    @classmethod
    def from_atom(cls, atom: 'MarkushiElement') -> 'MarkushiElement':
        """
        get Element copy
        """
        if not isinstance(atom, MarkushiElement):
            raise TypeError('Element expected')
        return atom.copy()

    @property
    def implicit_hydrogens(self) -> Optional[int]:
        try:
            return self._graph()._hydrogens[self._n]
        except AttributeError:
            raise IsNotConnectedAtom

    @property
    def explicit_hydrogens(self) -> int:
        try:
            return self._graph().explicit_hydrogens(self._n)
        except AttributeError:
            raise IsNotConnectedAtom

    # def __eq__(self, other):
    #     """
    #     compare attached to molecules elements
    #     """
    #     return isinstance(other, Element) and self.atomic_number == other.atomic_number and \
    #         self.isotope == other.isotope and self.charge == other.charge and self.is_radical == other.is_radical
    #
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
        elements_classes = {x.__name__: x.atomic_number.fget(None) for x in MarkushiElement.__subclasses__()}

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
        elements_classes = {x.__name__: x.atomic_number.fget(None) for x in MarkushiElement.__subclasses__()}

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


__all__ = ['MarkushiElement']