# -*- coding: utf-8 -*-
#
#  Copyright 2019-2021 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2019 Adelia Fatykhova <adelik21979@gmail.com>
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
from collections import deque
from functools import reduce
from itertools import count, permutations, combinations
from logging import info
from operator import or_
from typing import Iterable, List, Iterator, Tuple, Union
from .base import BaseReactor
from .._functions import lazy_product
from ..containers import QueryContainer, MoleculeContainer, ReactionContainer


class Reactor(BaseReactor):
    """
    Reactor for molecules transformations.
    Generates reaction from input molecules using transformation template.

    Reactor calling transforms reactants to products and
    returns generator of reaction transformations with all
    possible reactions.
    """
    def __init__(self, patterns: Tuple[QueryContainer, ...],
                 products: Tuple[Union[MoleculeContainer, QueryContainer], ...], *,
                 delete_atoms: bool = True, one_shot: bool = True,
                 polymerise_limit: int = 10, automorphism_filter: bool = True):
        """
        :param patterns: search patterns for each reactant
        :param products: resulted structures
        :param delete_atoms: if True atoms exists in reactants but
                            not exists in products will be removed
        :param one_shot: do only single reaction center then True, else do all possible combinations of reactions.
        :param polymerise_limit: limit of self reactions. Make sense than one_shot = False.
        """
        if not patterns or not products:
            raise ValueError('empty template')

        if not all(isinstance(x, QueryContainer) for x in patterns):
            raise TypeError('invalid params')
        elif not all(isinstance(x, (QueryContainer, MoleculeContainer)) for x in products):
            raise TypeError('invalid params')
        self.__patterns = patterns
        products_ = reduce(or_, products)

        self.__one_shot = one_shot
        self.__polymerise_limit = polymerise_limit
        self.__products_atoms = tuple(set(m) for m in products)
        self.__automorphism_filter = automorphism_filter
        super().__init__(reduce(or_, patterns), products_, delete_atoms)

    def __call__(self, structures: Iterable[MoleculeContainer]):
        if any(not isinstance(structure, MoleculeContainer) for structure in structures):
            raise TypeError('only list of Molecules possible')

        len_patterns = len(self.__patterns)
        structures = self.__remap(structures)
        s_nums = set(range(len(structures)))
        if self.__one_shot:
            for chosen in permutations(s_nums, len_patterns):
                ignored = [structures[x] for x in s_nums.difference(chosen)]
                chosen = [structures[x] for x in chosen]
                for new in self.__single_stage(chosen, {x for x in ignored for x in x}):
                    r = ReactionContainer([x.copy() for x in structures], new + [x.copy() for x in ignored])
                    if len(new) > 1:  # try to keep salts
                        r.contract_ions()
                    yield r
        else:
            queue = deque(([structures[x] for x in chosen], [structures[x] for x in s_nums.difference(chosen)], 0)
                          for chosen in permutations(s_nums, len_patterns))
            seen = set()
            while queue:
                chosen, ignored, depth = queue.popleft()
                depth += 1
                for new in self.__single_stage(chosen, {x for x in ignored for x in x}):
                    r = ReactionContainer([x.copy() for x in structures], new + [x.copy() for x in ignored])
                    if len(new) > 1:
                        r.contract_ions()  # try to keep salts
                        if str(r) in seen:
                            continue
                        seen.add(str(r))
                        if len(r.products) != len(ignored) + len(self.__products_atoms):
                            info('ambiguous multicomponent structures. skip multistage processing')
                            yield r
                            continue
                    elif str(r) in seen:
                        continue
                    else:
                        seen.add(str(r))

                    if depth < self.__polymerise_limit:
                        prod = r.products
                        if len_patterns == 1:  # simple case. only products or ignored can be transformed.
                            for i in range(len(prod)):
                                queue.append(([prod[i]], [*prod[:i], *prod[i + 1:]], depth))
                        else:  # one of products molecule combined with previously chosen
                            for chp in combinations(chosen, len_patterns - 1):
                                for i in range(len(prod)):
                                    for ch in permutations(self.__remap((prod[i], *chp)), len_patterns):
                                        queue.append((ch, [*prod[:i], *prod[i + 1:]], depth))
                    yield r

    def __single_stage(self, chosen, ignored) -> Iterator[List[MoleculeContainer]]:
        max_ignored_number = max(ignored, default=0)
        united_chosen = reduce(or_, chosen)
        split = len(self.__products_atoms) > 1
        for match in lazy_product(*(x.get_mapping(y, automorphism_filter=self.__automorphism_filter) for x, y in
                                    zip(self.__patterns, chosen))):
            mapping = match[0]
            for m in match[1:]:
                mapping.update(m)
            new = self._patcher(united_chosen, mapping)
            collision = set(new).intersection(ignored)
            if collision:
                new.remap(dict(zip(collision, count(max(max_ignored_number, max(new)) + 1))))

            if split:
                yield [new.substructure(c) for c in new._connected_components(new._bonds)]
            else:
                yield [new]

    @staticmethod
    def __remap(structures) -> List[MoleculeContainer]:
        checked = []
        checked_atoms = set()
        for structure in structures:
            intersection = set(structure).intersection(checked_atoms)
            if intersection:
                mapping = dict(zip(intersection, count(max(max(checked_atoms), max(structure)) + 1)))
                structure = structure.remap(mapping, copy=True)
                info('some atoms in input structures had the same numbers.\n'
                     f'atoms {list(mapping)} were remapped to {list(mapping.values())}')
            checked_atoms.update(structure)
            checked.append(structure)
        return checked

    def __getstate__(self):
        return {'patterns': self.__patterns, 'products_atoms': self.__products_atoms,
                'polymerise_limit': self.__polymerise_limit, 'one_shot': self.__one_shot,
                'automorphism_filter': self.__automorphism_filter, **super().__getstate__()}

    def __setstate__(self, state):
        if 'split' in state:
            raise ValueError('Reactor pickled with incompatible version.')
        self.__patterns = state['patterns']
        self.__one_shot = state['one_shot']
        self.__polymerise_limit = state['polymerise_limit']
        self.__products_atoms = state['products_atoms']
        self.__automorphism_filter = state['automorphism_filter']
        super().__setstate__(state)


__all__ = ['Reactor']
