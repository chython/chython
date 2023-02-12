# -*- coding: utf-8 -*-
#
#  Copyright 2019-2023 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from itertools import count, permutations, combinations, chain
from logging import getLogger, INFO
from operator import or_
from typing import List, Iterator, Tuple, Union
from .base import BaseReactor
from .._functions import lazy_product
from ..containers import QueryContainer, MoleculeContainer, ReactionContainer


logger = getLogger('chython.reactor')
logger.setLevel(INFO)


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
                 delete_atoms: bool = True, one_shot: bool = True, polymerise_limit: int = 10,
                 automorphism_filter: bool = True, fix_aromatic_rings: bool = True, fix_tautomers: bool = True):
        """
        :param patterns: Search patterns for each reactant.
        :param products: Resulted structures.
        :param delete_atoms: If True atoms exists in reactants but not exists in products will be removed.
        :param one_shot: Do only single reaction center then True, else do all possible combinations of reactions.
        :param polymerise_limit: Limit of self reactions. Make sense than one_shot = False.
        :param fix_aromatic_rings: Proceed kekule and thiele on products.
        :param fix_tautomers: See `thiele()` docs.
        :param automorphism_filter: Skip matches to same atoms.
        """
        if not patterns or not products:
            raise ValueError('empty template')

        if not all(isinstance(x, QueryContainer) for x in patterns):
            raise TypeError('invalid params')
        elif not all(isinstance(x, (QueryContainer, MoleculeContainer)) for x in products):
            raise TypeError('invalid params')
        self.patterns = patterns
        self.products = products

        self.__one_shot = one_shot
        self.__polymerise_limit = polymerise_limit
        self.__products_atoms = tuple(set(m) for m in products)
        self.__automorphism_filter = automorphism_filter
        super().__init__({n for x in patterns for n, h in x._masked.items() if not h}, reduce(or_, products),
                         delete_atoms, fix_aromatic_rings, fix_tautomers)

    def __call__(self, *structures: MoleculeContainer):
        if any(not isinstance(structure, MoleculeContainer) for structure in structures):
            raise TypeError('only list of Molecules possible')

        len_patterns = len(self.patterns)
        structures = fix_mapping_overlap(structures)
        s_nums = set(range(len(structures)))
        seen = set()
        black_list = [set() for _ in range(len_patterns)]  # unmatched patterns
        if self.__one_shot:
            for chosen in permutations(s_nums, len_patterns):
                if any(x in b for x, b in zip(chosen, black_list)):  # drop unmatched
                    continue

                chosen_molecules = [structures[x] for x in chosen]
                mapper = self.__get_mapping(chosen_molecules)
                if (n := next(mapper)) is not None:
                    black_list[n].add(chosen[n])
                    continue

                ignored = [structures[x] for x in s_nums.difference(chosen)]
                for new in self.__single_stage(chosen_molecules, {x for x in ignored for x in x}, mapper):
                    # store reacted molecules in same order as matched pattern
                    r = ReactionContainer([x.copy() for x in chosen_molecules] + [x.copy() for x in ignored],
                                          new + [x.copy() for x in ignored])
                    if len(new) > 1:  # try to keep salts
                        r.contract_ions()
                    if str(r) in seen:
                        continue
                    seen.add(str(r))
                    yield r
        else:
            queue = deque(([structures[x] for x in chosen], [structures[x] for x in s_nums.difference(chosen)],
                           lazy_product(*(q.get_mapping(structures[x], automorphism_filter=self.__automorphism_filter)
                                          for q, x in zip(self.patterns, chosen))), 0)
                          for chosen in permutations(s_nums, len_patterns))
            while queue:
                chosen, ignored, mapper, depth = queue.popleft()
                depth += 1
                for new in self.__single_stage(chosen, {x for x in ignored for x in x}, mapper):
                    r = ReactionContainer([x.copy() for x in structures], new + [x.copy() for x in ignored])
                    if len(new) > 1:
                        r.contract_ions()  # try to keep salts
                        if str(r) in seen:
                            continue
                        seen.add(str(r))
                        if len(r.products) != len(ignored) + len(self.__products_atoms):
                            logger.info('ambiguous multicomponent structures. skip multistage processing')
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
                                    for ch in permutations(fix_mapping_overlap((prod[i], *chp)), len_patterns):
                                        queue.append((ch, [*prod[:i], *prod[i + 1:]], depth))
                    yield r

    def __get_mapping(self, chosen):
        mappers = []
        for n, (q, x) in enumerate(zip(self.patterns, chosen)):
            iso = q.get_mapping(x, automorphism_filter=self.__automorphism_filter)
            try:
                m = next(iso)
            except StopIteration:
                yield n
                return
            mappers.append(chain((m,), iso))
        yield
        yield from lazy_product(*mappers)

    def __single_stage(self, chosen, ignored, mapper) -> Iterator[List[MoleculeContainer]]:
        max_ignored_number = max(ignored, default=0)
        united_chosen = reduce(or_, chosen)
        split = len(self.__products_atoms) > 1
        for match in mapper:
            mapping = match[0].copy()
            for m in match[1:]:
                mapping.update(m)
            for new in self._patcher(united_chosen, mapping):
                collision = set(new).intersection(ignored)
                if collision:
                    new.remap(dict(zip(collision, count(max(max_ignored_number, max(new)) + 1))))

                if split:
                    yield new.split()
                else:
                    yield [new]


def fix_mapping_overlap(structures) -> List[MoleculeContainer]:
    checked = []
    checked_atoms = set()
    for structure in structures:
        intersection = set(structure).intersection(checked_atoms)
        if intersection:
            mapping = dict(zip(intersection, count(max(max(checked_atoms), max(structure)) + 1)))
            structure = structure.remap(mapping, copy=True)
            logger.info('some atoms in input structures had the same numbers.\n'
                        f'atoms {list(mapping)} were remapped to {list(mapping.values())}')
        checked_atoms.update(structure)
        checked.append(structure)
    return checked


__all__ = ['Reactor', 'fix_mapping_overlap']
