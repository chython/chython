# -*- coding: utf-8 -*-
#
#  Copyright 2018-2023 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from array import array
from collections import defaultdict, deque
from functools import cached_property, partial
from itertools import permutations
from typing import Any, Collection, Dict, Iterator, Optional, TYPE_CHECKING, Union
from .._functions import lazy_product
from ..periodictable.element import Element, Query, AnyElement, AnyMetal, ListElement


if TYPE_CHECKING:
    from chython.containers.graph import Graph
    from chython.containers import MoleculeContainer, QueryContainer


class Isomorphism:
    __slots__ = ()

    def __lt__(self, other):
        if len(self) >= len(other):
            return False
        return self.is_substructure(other)

    def __le__(self, other):
        return self.is_substructure(other)

    def __gt__(self, other):
        if len(self) <= len(other):
            return False
        return other.is_substructure(self)

    def __ge__(self, other):
        return other.is_substructure(self)

    def __contains__(self: 'Graph', other: Union[Element, Query, str]):
        """
        Atom in Structure test.
        """
        if isinstance(other, str):
            return any(other == x.atomic_symbol for x in self._atoms.values())
        return any(other == x for x in self._atoms.values())

    def is_substructure(self, other, /) -> bool:
        """
        Test self is substructure of other
        """
        try:
            next(self.get_mapping(other, automorphism_filter=False))
        except StopIteration:
            return False
        return True

    def is_equal(self, other, /) -> bool:
        """
        Test self is same structure as other
        """
        if len(self) != len(other):
            return False
        try:
            next(self.get_mapping(other, automorphism_filter=False))
        except StopIteration:
            return False
        return True

    def is_automorphic(self):
        """
        Test for automorphism symmetry of graph.
        """
        try:
            next(self.get_automorphism_mapping())
        except StopIteration:
            return False
        return True

    def get_automorphism_mapping(self: 'Graph') -> Iterator[Dict[int, int]]:
        """
        Iterator of all possible automorphism mappings.
        """
        return _get_automorphism_mapping(self.atoms_order, self._bonds)

    def _get_mapping(self, other, /, *, automorphism_filter=True, searching_scope=None,
                     components=None, get_mapping=None) -> Iterator[Dict[int, int]]:
        if components is None:  # ad-hoc for QueryContainer
            components, closures = self._compiled_query
            get_mapping = partial(_get_mapping, query_closures=closures, o_atoms=other._atoms, o_bonds=other._bonds)

        if searching_scope is not None and not isinstance(searching_scope, set):
            searching_scope = set(searching_scope)

        seen = set()
        if len(components) == 1:
            for candidate in other._connected_components:
                if searching_scope:
                    candidate = searching_scope.intersection(candidate)
                    if not candidate:
                        continue
                for mapping in get_mapping(components[0], scope=candidate):
                    if automorphism_filter:
                        atoms = frozenset(mapping.values())
                        if atoms in seen:
                            continue
                        seen.add(atoms)
                    yield mapping
        else:
            for candidates in permutations(other._connected_components, len(components)):
                mappers = []
                for component, candidate in zip(components, candidates):
                    if searching_scope:
                        candidate = searching_scope.intersection(candidate)
                        if not candidate:
                            break
                    mappers.append(get_mapping(component, scope=candidate))
                else:
                    for match in lazy_product(*mappers):
                        mapping = match[0].copy()
                        for m in match[1:]:
                            mapping.update(m)
                        if automorphism_filter:
                            atoms = frozenset(mapping.values())
                            if atoms in seen:
                                continue
                            seen.add(atoms)
                        yield mapping

    @cached_property
    def _compiled_query(self: 'Graph'):
        components, closures = _compile_query(self._atoms, self._bonds)
        if self.connected_components_count > 1:
            order = {x: n for n, c in enumerate(self.connected_components) for x in c}
            components.sort(key=lambda x: order[x[0][0]])
        return components, closures


class MoleculeIsomorphism(Isomorphism):
    def get_mapping(self, other: 'MoleculeContainer', /, *, automorphism_filter: bool = True,
                    searching_scope: Optional[Collection[int]] = None):
        """
        Get self to other Molecule substructure mapping generator.

        :param other: Molecule
        :param automorphism_filter: Skip matches to the same atoms.
        :param searching_scope: substructure atoms list to localize isomorphism.
        """
        if isinstance(other, MoleculeIsomorphism):
            return self._get_mapping(other, automorphism_filter=automorphism_filter, searching_scope=searching_scope)
        raise TypeError('MoleculeContainer expected')

    @cached_property
    def _cython_compiled_structure(self):
        # long I:
        # bond: single, double, triple, aromatic, special = 5 bit
        # bond in ring: 2 bit
        # atom: H-Ba: 56 bit
        # transfer bit

        # long II:
        # atom La-Mc: 59 bit
        # Lv-Ts-Og: 3 elements packed into 1 bit.
        # hybridizations: 1-4 = 4 bit

        # long III:
        # isotope: not specified, isotope - common_isotope = -8 - +8 = 18 bit
        # is_radical: 2 bit
        # charge: -4 - +4: 9 bit
        # implicit_hydrogens: 0-4 = 5 bit
        # neighbors: 0-14 = 15 bit
        # heteroatoms: 0-14 = 15 bit

        # long IV:
        # ring_sizes: not-in-ring bit, 3-atom ring, 4-...., 65-atom ring
        from ..files._mdl.mol import common_isotopes

        charges = self._charges
        radicals = self._radicals
        hydrogens = self._hydrogens
        neighbors = self.neighbors
        heteroatoms = self.heteroatoms
        rings_sizes = self.atoms_rings_sizes
        hybridization = self.hybridization

        mapping = {}
        numbers = []
        bits1 = []
        bits2 = []
        bits3 = []
        bits4 = []
        for i, (n, a) in enumerate(self._atoms.items()):
            mapping[n] = i
            numbers.append(n)
            v2 = 1 << (hybridization(n) - 1)
            if (an := a.atomic_number) > 56:
                if an > 116:  # Ts, Og
                    an = 116
                v1 = 1  # transfer bit
                v2 |= 1 << (120 - an)
            else:
                v1 = 1 << (57 - an)

            if a.isotope:
                v3 = 1 << (a.isotope - common_isotopes[a.atomic_symbol] + 54)
                if radicals[n]:
                    v3 |= 0x200000000000
                else:
                    v3 |= 0x100000000000
            elif radicals[n]:
                v3 = 0x8000200000000000
            else:
                v3 = 0x8000100000000000

            v3 |= 1 << (charges[n] + 39)
            v3 |= 1 << ((hydrogens[n] or 0) + 30)
            v3 |= 1 << (neighbors(n) + 15)
            v3 |= 1 << heteroatoms(n)

            if n in rings_sizes:
                v4 = 0
                for r in rings_sizes[n]:
                    if r > 65:  # big rings not supported
                        continue
                    v4 |= 1 << (65 - r)
                if not v4:  # only 65+ rings. set as rings-free.
                    v4 = 0x8000000000000000
            else:  # not in rings
                v4 = 0x8000000000000000

            bits1.append(v1)
            bits2.append(v2)
            bits3.append(v3)
            bits4.append(v4)

        o_from = [0] * len(mapping)
        o_to = [0] * len(mapping)
        indices = [0] * self.bonds_count * 2
        bonds = [0] * self.bonds_count * 2
        start = 0
        for n, ms in self._bonds.items():
            i = mapping[n]
            o_from[i] = start
            for j, (m, b) in enumerate(ms.items(), start):
                indices[j] = x = mapping[m]
                v = bits1[x]
                o = b.order
                if o == 1:
                    v |= 0x0800000000000000
                elif o == 4:
                    v |= 0x4000000000000000
                elif o == 2:
                    v |= 0x1000000000000000
                elif o == 3:
                    v |= 0x2000000000000000
                else:
                    v |= 0x8000000000000000
                v |= 0x0400000000000000 if b.in_ring else 0x0200000000000000
                bonds[j] = v
            start += len(ms)
            o_to[i] = start

        return (array('L', numbers), array('Q', bits1), array('Q', bits2), array('Q', bits3), array('Q', bits4),
                array('Q', bonds), array('I', o_from), array('I', o_to), array('I', indices))


class QueryIsomorphism(Isomorphism):
    def get_mapping(self, other: Union['MoleculeContainer', 'QueryContainer'], /, *, automorphism_filter: bool = True,
                    searching_scope: Optional[Collection[int]] = None, _cython=True):
        """
        Get self to other Molecule or Query substructure mapping generator.

        :param other: Molecule or Query
        :param automorphism_filter: Skip matches to the same atoms.
        :param searching_scope: substructure atoms list to localize isomorphism.
        """
        # _cython - by default cython implementation enabled.
        # disable it by overriding method if Query Atoms or Containers logic changed.
        # Lv, Ts and Og in cython optimized mode treated as equal.
        if isinstance(other, QueryIsomorphism):
            return self._get_mapping(other, automorphism_filter=automorphism_filter, searching_scope=searching_scope)
        elif isinstance(other, MoleculeIsomorphism):
            if _cython:
                try:  # windows? ;)
                    from ._isomorphism import get_mapping as _cython_get_mapping
                except ImportError:
                    components = get_mapping = None
                else:
                    components = self._cython_compiled_query  # override to cython data

                    def get_mapping(query, scope):
                        return _cython_get_mapping(*query, *other._cython_compiled_structure,
                                                   array('I', [n in scope for n in other]))
            else:
                components = get_mapping = None
            return self._get_mapping(other, automorphism_filter=automorphism_filter, searching_scope=searching_scope,
                                     components=components, get_mapping=get_mapping)
        raise TypeError('MoleculeContainer or QueryContainer expected')

    @cached_property
    def _cython_compiled_query(self):
        # long I:
        # bond: single, double, triple, aromatic, special = 5 bit
        # bond in ring: 2 bit
        # atom: H-Ba: 56 bit
        # transfer bit

        # long II:
        # atom La-Mc: 59 bit
        # Lv-Ts-Og: 3 elements packed into 1 bit.
        # hybridizations: 1-4 = 4 bit

        # long III:
        # isotope: not specified, isotope - common_isotope = -8 - +8 = 18 bit
        # is_radical: 2 bit
        # charge: -4 - +4: 9 bit
        # implicit_hydrogens: 0-4 = 5 bit
        # neighbors: 0-14 = 15 bit
        # heteroatoms: 0-14 = 15 bit

        # long IV:
        # ring_sizes: not-in-ring bit, 3-atom ring, 4-...., 65-atom ring

        # int V: bonds closures
        # padding: 1 bit
        # bond: single, double, triple, aromatic, special = 5 bit
        # bond in ring: 2 bit
        from ..files._mdl.mol import common_isotopes

        _components, _closures = self._compiled_query
        components = []
        for c in _components:
            mapping = {n: i for i, (n, *_) in enumerate(c)}
            masks1 = []
            masks2 = []
            masks3 = []
            masks4 = []
            for *_, a, b in c:
                if isinstance(a, AnyMetal):  # isotope, radical, charge, hydrogens and heteroatoms states ignored
                    # elements except 1, 2, 5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18, 32, 33, 34, 35, 36, 51, 52, 53, 54
                    v1 = 0x0060707ffc1fff87
                    v2 = 0xfffffffffffffff0
                    v3 = 0xffffffffc0007fff
                    v4 = 0xffffffffffffffff
                else:
                    if isinstance(a, AnyElement):
                        v1 = 0x01ffffffffffffff
                        v2 = 0xfffffffffffffff0
                    else:
                        if isinstance(a, ListElement):
                            v1 = v2 = 0
                            for n in a._numbers:
                                if n > 56:
                                    if n > 116:  # Ts, Og
                                        n = 116
                                    v1 |= 1  # set transfer bit
                                    v2 |= 1 << (120 - n)
                                else:
                                    v1 |= 1 << (57 - n)
                        elif (n := a.atomic_number) > 56:
                            if n > 116:  # Ts, Og
                                n = 116
                            v1 = 1  # transfer bit
                            v2 = 1 << (120 - n)
                        else:
                            v1 = 1 << (57 - n)
                            v2 = 0
                    if a.isotope:
                        v3 = 1 << (a.isotope - common_isotopes[a.atomic_symbol] + 54)
                        if a.is_radical:
                            v3 |= 0x200000000000
                        else:
                            v3 |= 0x100000000000
                    elif a.is_radical:  # any isotope
                        v3 = 0xffffe00000000000
                    else:
                        v3 = 0xffffd00000000000

                    v3 |= 1 << (a.charge + 39)

                    if not a.implicit_hydrogens:
                        v3 |= 0x7c0000000
                    else:
                        for h in a.implicit_hydrogens:
                            v3 |= 1 << (h + 30)

                    if not a.heteroatoms:
                        v3 |= 0x7fff
                    else:
                        for n in a.heteroatoms:
                            v3 |= 1 << n

                    if a.ring_sizes:
                        if a.ring_sizes[0]:
                            v4 = 0
                            for r in a.ring_sizes:
                                if r > 65:  # big rings not supported
                                    continue
                                v4 |= 1 << (65 - r)
                            if not v4:  # only 65+ rings. set as rings-free.
                                v4 = 0x8000000000000000
                        else:  # not in rings
                            v4 = 0x8000000000000000
                    else:  # any rings
                        v4 = 0xffffffffffffffff

                if not a.neighbors:
                    v3 |= 0x3fff8000
                else:
                    for n in a.neighbors:
                        v3 |= 1 << (n + 15)

                if not a.hybridization:
                    v2 |= 0xf
                else:
                    for n in a.hybridization:
                        v2 |= 1 << (n - 1)

                if b is not None:
                    for o in b.order:
                        if o == 1:
                            v1 |= 0x0800000000000000
                        elif o == 4:
                            v1 |= 0x4000000000000000
                        elif o == 2:
                            v1 |= 0x1000000000000000
                        elif o == 3:
                            v1 |= 0x2000000000000000
                        else:
                            v1 |= 0x8000000000000000
                    if b.in_ring is None:
                        v1 |= 0x0600000000000000
                    elif b.in_ring:
                        v1 |= 0x0400000000000000
                    else:
                        v1 |= 0x0200000000000000

                masks1.append(v1)
                masks2.append(v2)
                masks3.append(v3)
                masks4.append(v4)

            closures = [0] * len(c)  # closures amount
            q_from = [0] * len(c)
            q_to = [0] * len(c)
            indices = [0] * sum(len(ms) for n, ms in _closures.items() if n in mapping)
            bonds = indices.copy()

            start = 0
            for n, ms in _closures.items():
                if (i := mapping.get(n)) is not None:
                    closures[i] = len(ms)
                    q_from[i] = start
                    for j, (m, b) in enumerate(ms, start):
                        v = 0x01ffffffffffffff  # atom doesn't matter.
                        for o in b.order:
                            if o == 1:
                                v |= 0x0800000000000000
                            elif o == 4:
                                v |= 0x4000000000000000
                            elif o == 2:
                                v |= 0x1000000000000000
                            elif o == 3:
                                v |= 0x2000000000000000
                            else:
                                v |= 0x8000000000000000
                        if b.in_ring is None:
                            v |= 0x0600000000000000
                        elif b.in_ring:
                            v |= 0x0400000000000000
                        else:
                            v |= 0x0200000000000000
                        bonds[j] = v
                        indices[j] = mapping[m]
                    start += len(ms)
                    q_to[i] = start
            components.append((array('L', [n for n, *_ in c]), array('I', [0] + [mapping[x] for _, x, *_ in c[1:]]),
                               array('Q', masks1), array('Q', masks2), array('Q', masks3), array('Q', masks4),
                               array('I', closures), array('I', q_from), array('I', q_to),
                               array('I', indices), array('Q', bonds)))
        return components


def _get_automorphism_mapping(atoms: Dict[int, int], bonds: Dict[int, Dict[int, Any]]) -> Iterator[Dict[int, int]]:
    if len(atoms) == len(set(atoms.values())):
        return  # all atoms unique

    components, closures = _compile_query(atoms, bonds)
    mappers = [_get_mapping(order, closures, atoms, bonds, {x for x, *_ in order})
               for order in components]
    if len(mappers) == 1:
        for mapping in mappers[0]:
            if any(k != v for k, v in mapping.items()):
                yield mapping
    for match in lazy_product(*mappers):
        mapping = match[0].copy()
        for m in match[1:]:
            mapping.update(m)
        if any(k != v for k, v in mapping.items()):
            yield mapping


def _get_mapping(linear_query, query_closures, o_atoms, o_bonds, scope):
    size = len(linear_query) - 1
    order_depth = {v[0]: k for k, v in enumerate(linear_query)}

    stack = deque()
    path = []
    mapping = {}
    reversed_mapping = {}

    s_n, _, s_atom, _ = linear_query[0]
    for n, o_atom in o_atoms.items():
        if n in scope and s_atom == o_atom:
            stack.append((n, 0))

    while stack:
        n, depth = stack.pop()
        current = linear_query[depth][0]
        if depth == size:
            yield {**mapping, current: n}
        else:
            if len(path) != depth:
                for x in path[depth:]:
                    del mapping[reversed_mapping.pop(x)]
                path = path[:depth]

            path.append(n)
            mapping[current] = n
            reversed_mapping[n] = current

            depth += 1
            s_n, back, s_atom, s_bond = linear_query[depth]
            if back != current:
                n = path[order_depth[back]]

            for o_n, o_bond in o_bonds[n].items():
                if o_n in scope and o_n not in reversed_mapping and s_bond == o_bond:
                    if s_atom == o_atoms[o_n]:
                        # check closures equality
                        o_closures = o_bonds[o_n].keys() & reversed_mapping.keys()
                        o_closures.discard(n)
                        if o_closures == {mapping[m] for m, _ in query_closures[s_n]}:
                            obon = o_bonds[o_n]
                            if all(bond == obon[mapping[m]] for m, bond in query_closures[s_n]):
                                stack.append((o_n, depth))


def _compile_query(atoms, bonds):
    closures = defaultdict(list)
    components = []
    seen = set()
    iter_atoms = iter(atoms)
    while len(seen) < len(atoms):
        start = next(x for x in iter_atoms if x not in seen)
        seen.add(start)
        stack = [(n, start, atoms[n], bond) for n, bond in reversed(bonds[start].items())]
        order = [(start, None, atoms[start], None)]
        components.append(order)

        while stack:
            front, back, *_ = atom = stack.pop()
            if front not in seen:
                order.append(atom)
                for n, bond in reversed(bonds[front].items()):
                    if n != back:
                        if n not in seen:
                            stack.append((n, front, atoms[n], bond))
                        else:
                            closures[front].append((n, bond))
                seen.add(front)
    return components, closures


__all__ = ['MoleculeIsomorphism', 'QueryIsomorphism']
