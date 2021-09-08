# -*- coding: utf-8 -*-
#
#  Copyright 2017-2021 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2019 Timur Gimadiev <timur.gimadiev@gmail.com>
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
from CachedMethods import cached_method
from collections import defaultdict
from functools import cached_property
from hashlib import sha512
from itertools import count, product
from random import random
from typing import Tuple, TYPE_CHECKING, Union


if TYPE_CHECKING:
    from chython import MoleculeContainer, CGRContainer, QueryContainer
    from chython.containers.graph import Graph

charge_str = {-4: '-4', -3: '-3', -2: '-2', -1: '-', 0: '0', 1: '+', 2: '+2', 3: '+3', 4: '+4'}
order_str = {1: '-', 2: '=', 3: '#', 4: ':', 8: '~', None: '.'}
organic_set = {'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I', 'B'}
hybridization_str = {4: 'a', 3: 't', 2: 'd', 1: 's', None: 'n'}
dyn_order_str = {(None, 1): '[.>-]', (None, 2): '[.>=]', (None, 3): '[.>#]', (None, 4): '[.>:]', (None, 8): '[.>~]',
                 (1, None): '[->.]', (1, 1): '', (1, 2): '[->=]', (1, 3): '[->#]', (1, 4): '[->:]', (1, 8): '[->~]',
                 (2, None): '[=>.]', (2, 1): '[=>-]', (2, 2): '=', (2, 3): '[=>#]', (2, 4): '[=>:]', (2, 8): '[=>~]',
                 (3, None): '[#>.]', (3, 1): '[#>-]', (3, 2): '[#>=]', (3, 3): '#', (3, 4): '[#>:]', (3, 8): '[#>~]',
                 (4, None): '[:>.]', (4, 1): '[:>-]', (4, 2): '[:>=]', (4, 3): '[:>#]', (4, 4): ':', (4, 8): '[:>~]',
                 (8, None): '[~>.]', (8, 1): '[~>-]', (8, 2): '[~>=]', (8, 3): '[~>#]', (8, 4): '[~>:]', (8, 8): '~'}

dyn_charge_str = {(i, j): f'{charge_str[i]}>{charge_str[j]}' if i != j else charge_str[i]
                  for i, j in product(range(-4, 5), repeat=2)}
dyn_charge_str[(0, 0)] = ''

dyn_radical_str = {(True, True): '*', (True, False): '*>^', (False, True): '^>*'}


class Smiles(ABC):
    __slots__ = ()

    @cached_method
    def __str__(self):
        return ''.join(self._smiles(self._smiles_order))

    def __format__(self: Union['Graph', 'Smiles'], format_spec):
        """
        Signature generation options.

        :param format_spec: String with keys:
            a - Generate asymmetric closures.
            !s - Disable stereo marks.
            A - Use aromatic bonds instead aromatic atoms.
            m - Set atom mapping.
            r - Generate random-ordered smiles.
            o - Old canonic ordering algorithm.

            Combining possible. Order independent. Another keys ignored.
        """
        if format_spec:
            kwargs = {}
            if 'a' in format_spec:
                kwargs['asymmetric_closures'] = True
            if '!s' in format_spec:
                kwargs['stereo'] = False
            if 'A' in format_spec:
                kwargs['aromatic'] = False
            if 'm' in format_spec:
                kwargs['mapping'] = True
            if 'r' in format_spec:
                kwargs['random'] = True

                def w(_):
                    return random()
            elif 'o' in format_spec:
                kwargs['old_order'] = True
                w = self.atoms_order.get
            else:
                w = self._smiles_order
            return ''.join(self._smiles(w, **kwargs))
        return str(self)

    def __eq__(self, other):
        return isinstance(other, Smiles) and str(self) == str(other)

    @cached_method
    def __hash__(self):
        return hash(str(self))

    @cached_method
    def __bytes__(self):
        return sha512(str(self).encode()).digest()

    @cached_property
    def smiles_atoms_order(self) -> Tuple[int, ...]:
        """
        Atoms order in canonic SMILES.
        """
        smiles, order = self._smiles(self._smiles_order, _return_order=True)
        self.__dict__['__cached_method___str__'] = ''.join(smiles)  # cache smiles also
        return tuple(order)

    def _smiles(self: Union['Graph', 'Smiles'], weights, *, asymmetric_closures=False,
                open_parenthesis='(', close_parenthesis=')', delimiter='.', _return_order=False, **kwargs):
        if not self._atoms:
            return []
        bonds = self._bonds
        atoms_set = set(self._atoms)
        seen = {}
        cycles = count()
        casted_cycles = {}
        string = []
        order = []
        if asymmetric_closures:
            visited_bond = set()

        groups = defaultdict(int)
        for n in atoms_set:
            groups[weights(n)] -= 1

        if kwargs.get('random', False):
            mod_weights_start = mod_weights = weights
        elif kwargs.get('old_order', False):
            def mod_weights_start(x):
                lb = len(bonds[x])
                if lb:
                    return (-groups[weights(x)],  # rare groups
                            -lb,  # more neighbors
                            lb / len({weights(x) for x in bonds[x]}),  # more unique neighbors
                            weights(x))  # smallest weight
                else:
                    return -groups[weights(x)], weights(x)  # rare groups > smallest weight

            def mod_weights(x):
                lb = len(bonds[x])
                return (-groups[weights(x)],  # rare groups
                        -lb,  # more neighbors
                        lb / len({weights(x) for x in bonds[x]}),  # more unique neighbors
                        weights(x),  # smallest weight
                        seen[x])  # BFS nearest to starting
        else:
            def mod_weights_start(x):
                return (groups[weights(x)],  # common groups
                        weights(x))  # smallest weight

            def mod_weights(x):
                return (groups[weights(x)],  # common groups
                        weights(x),  # smallest weight
                        seen[x])  # BFS nearest to starting

        while True:
            start = min(atoms_set, key=mod_weights_start)
            seen[start] = 0
            queue = [(start, 1)]
            while queue:
                n, d = queue.pop(0)
                for m in bonds[n].keys() - seen.keys():
                    queue.append((m, d + 1))
                    seen[m] = d

            # modified NX dfs with cycle detection
            stack = [(start, len(atoms_set), iter(sorted(bonds[start], key=mod_weights)))]
            visited = {start: []}  # predecessors for stereo. atom: (visited[atom], *edges[atom])
            disconnected = set()
            edges = defaultdict(list)
            tokens = defaultdict(list)
            while stack:
                parent, depth_now, children = stack[-1]
                try:
                    child = next(children)
                except StopIteration:
                    stack.pop()
                else:
                    if child not in visited:
                        edges[parent].append(child)
                        visited[child] = [parent]
                        if depth_now > 1:
                            front = bonds[child].keys() - {parent}
                            if front:
                                stack.append((child, depth_now - 1, iter(sorted(front, key=mod_weights))))
                    elif (child, parent) not in disconnected:
                        disconnected.add((parent, child))
                        disconnected.add((child, parent))
                        cycle = next(cycles)
                        tokens[parent].append((child, cycle))
                        tokens[child].append((parent, cycle))

            # flatten directed graph: edges
            stack = [[start, 0, [start]]]
            while True:
                tail, closure, smiles = stack[-1]
                if tail in edges:
                    children = edges[tail]
                    if len(children) > 1:  # has side chain
                        child = children[-1]
                        stack_len = len(stack)
                        stack.append([child, 0, [(tail, child), child]])  # end of current chain
                        for child in children[-2::-1]:  # start side chains
                            stack.append([child, stack_len, ['(', (tail, child), child]])
                    else:  # chain grow
                        child = children[-1]
                        stack[-1][0] = child
                        smiles.append((tail, child))
                        smiles.append(child)
                elif closure:  # end of side chain
                    stack.pop()
                    if smiles[-2] == '(':
                        smiles.pop(-2)
                    else:
                        smiles.append(')')
                    stack[closure - 1][2].extend(smiles)
                elif len(stack) > 2:
                    stack.pop()
                    stack[-1][0] = tail
                    stack[-1][2].extend(smiles)
                elif len(stack) == 2:
                    stack[0][2].extend(smiles)
                    smiles = stack[0][2]
                    break
                else:
                    break

            # prepare new neighbors order for stereo sign calculation
            for token in smiles:
                if token in tokens:
                    tokens[token].sort(key=lambda x: casted_cycles.get(x[1]) or
                                                     casted_cycles.setdefault(x[1], len(casted_cycles) + 1))
                    visited[token].extend(n for n, _ in tokens[token])
                if token in edges:
                    visited[token].extend(edges[token])

            for token in smiles:
                if isinstance(token, int):  # atoms
                    string.append(self._format_atom(token, adjacency=visited, **kwargs))
                    order.append(token)
                    if token in tokens:
                        for m, c in tokens[token]:
                            if asymmetric_closures:
                                if (token, m) not in visited_bond:
                                    string.append(self._format_bond(token, m, adjacency=visited, **kwargs))
                                    visited_bond.add((m, token))
                            else:
                                string.append(self._format_bond(token, m, adjacency=visited, **kwargs))
                            string.append(self._format_closure(casted_cycles[c]))
                elif token == '(':
                    string.append(open_parenthesis)
                elif token == ')':
                    string.append(close_parenthesis)
                else:  # bonds
                    string.append(self._format_bond(*token, adjacency=visited, **kwargs))

            atoms_set.difference_update(visited)
            if atoms_set:
                string.append(delimiter)
            else:
                break
        if _return_order:
            return string, order
        return string

    @staticmethod
    def _format_closure(c):
        return str(c) if c < 10 else f'%{c}'

    @abstractmethod
    def _format_atom(self, n, adjacency, **kwargs):
        ...

    @abstractmethod
    def _format_bond(self, n, m, adjacency, **kwargs):
        ...

    @property
    @abstractmethod
    def _smiles_order(self):
        ...


class MoleculeSmiles(Smiles):
    __slots__ = ()

    @property
    def _smiles_order(self: 'MoleculeContainer'):
        return self._chiral_morgan.__getitem__

    def _format_atom(self: 'MoleculeContainer', n, adjacency, **kwargs):
        atom = self._atoms[n]
        charge = self._charges[n]
        ih = self._hydrogens[n]
        hyb = self.hybridization(n)

        smi = ['',  # [
               str(atom.isotope) if atom.isotope else '',  # isotope
               None,
               '',  # stereo
               '',  # hydrogen
               '',  # charge
               f':{n}' if kwargs.get('mapping', False) else '',  # mapping
               '']  # ]

        if kwargs.get('stereo', True):
            if n in self._atoms_stereo:
                if ih and next(x for x in adjacency) == n:  # first atom in smiles has reversed chiral mark
                    smi[3] = '@@' if self._translate_tetrahedron_sign(n, adjacency[n]) else '@'
                else:
                    smi[3] = '@' if self._translate_tetrahedron_sign(n, adjacency[n]) else '@@'
            elif n in self._allenes_stereo:
                t1, t2 = self._stereo_allenes_terminals[n]
                env = self._stereo_allenes[n]
                n1 = next(x for x in adjacency[t1] if x in env)
                n2 = next(x for x in adjacency[t2] if x in env)
                smi[3] = '@' if self._translate_allene_sign(n, n1, n2) else '@@'
            elif charge:
                smi[5] = charge_str[charge]
        elif charge:
            smi[5] = charge_str[charge]

        if any(smi) or atom.atomic_symbol not in organic_set or self._radicals[n]:
            smi[0] = '['
            smi[-1] = ']'
            if ih == 1:
                smi[4] = 'H'
            elif ih:
                smi[4] = f'H{ih}'
        elif hyb == 4 and ih and atom.atomic_number in (5, 7, 15):  # pyrole
            smi[0] = '['
            smi[-1] = ']'
            if ih == 1:
                smi[4] = 'H'
            else:
                smi[4] = f'H{ih}'

        if kwargs.get('aromatic', True) and hyb == 4:
            smi[2] = atom.atomic_symbol.lower()
        else:
            smi[2] = atom.atomic_symbol
        return ''.join(smi)

    def _format_bond(self: 'MoleculeContainer', n, m, adjacency, **kwargs):
        order = self._bonds[n][m].order
        if order == 4:
            if kwargs.get('aromatic', True):
                return ''
            return ':'
        elif kwargs.get('stereo', True) and order == 1:  # cis-trans /\
            ctt = self._stereo_cis_trans_terminals
            if n in ctt:
                ts = ctt[n]
                if ts in self._cis_trans_stereo:
                    env = self._stereo_cis_trans[ts]
                    if m == next(x for x in adjacency[n] if x in env):  # only first neighbor of double bonded atom
                        if n == next(x for x in adjacency if x in ts):  # Cn(\Rm)(X)=C, Cn(=C)(\Rm)X, C(=C=Cn(\Rm)X)=C
                            return '\\'
                        else:  # C=Cn(Rm)(X) cases
                            n2 = ts[1] if ts[0] == n else ts[0]
                            m2 = next(x for x in adjacency[n2] if x in env)
                            return '\\' if self._translate_cis_trans_sign(n2, n, m2, m) else '/'
            elif m in ctt and self._atoms[n].atomic_number != 1:  # Hn-Cm= case ignored!
                ts = ctt[m]
                if ts in self._cis_trans_stereo:
                    if m == next(x for x in adjacency if x in ts):  # Rn-Cm(X)=C case
                        return '/'  # always start with UP R/C=C-X. RUSSIANS POSITIVE!
                    else:  # second RnCm=1X or R1...=C1(X) case
                        env = self._stereo_cis_trans[ts]
                        n2 = ts[1] if ts[0] == m else ts[0]
                        m2 = next(x for x in adjacency[n2] if x in env)
                        return '/' if self._translate_cis_trans_sign(n2, m, m2, n) else '\\'

            if kwargs.get('aromatic', True) and self.hybridization(n) == self.hybridization(m) == 4:
                return '-'
            return ''
        elif order == 2:
            return '='
        elif order == 3:
            return '#'
        elif order == 1:
            if kwargs.get('aromatic', True) and self.hybridization(n) == self.hybridization(m) == 4:
                return '-'
            return ''
        else:  # order == 8
            return '~'


class CGRSmiles(Smiles):
    __slots__ = ()

    @property
    def _smiles_order(self: 'CGRContainer'):
        return self.atoms_order.__getitem__

    def _format_atom(self: 'CGRContainer', n, **kwargs):
        atom = self._atoms[n]
        charge = self._charges[n]
        is_radical = self._radicals[n]
        p_charge = self._p_charges[n]
        p_is_radical = self._p_radicals[n]
        if atom.isotope:
            smi = [str(atom.isotope), atom.atomic_symbol]
        else:
            smi = [atom.atomic_symbol]

        if charge or p_charge:
            smi.append(dyn_charge_str[(charge, p_charge)])
        if is_radical or p_is_radical:
            smi.append(dyn_radical_str[(is_radical, p_is_radical)])

        if len(smi) != 1 or atom.atomic_symbol not in organic_set:
            smi.insert(0, '[')
            smi.append(']')
        return ''.join(smi)

    def _format_bond(self: 'CGRContainer', n, m, **kwargs):
        bond = self._bonds[n][m]
        return dyn_order_str[(bond.order, bond.p_order)]


class QuerySmiles(Smiles):
    __slots__ = ()

    @property
    def _smiles_order(self: 'QueryContainer'):
        return self.atoms_order.__getitem__

    def _format_atom(self: 'QueryContainer', n, **kwargs):
        atom = self._atoms[n]
        charge = self._charges[n]
        hybridization = self._hybridizations[n]
        neighbors = self._neighbors[n]
        hydrogens = self._hydrogens[n]
        rings = self._rings_sizes[n]
        heteroatoms = self._heteroatoms[n]

        if atom.isotope:
            smi = ['[', str(atom.isotope), atom.atomic_symbol]
        else:
            smi = ['[', atom.atomic_symbol]

        if n in self._atoms_stereo:  # carbon only
            smi.append('@' if self._translate_tetrahedron_sign(n, kwargs['adjacency'][n]) else '@@')

        if neighbors:
            smi.append(';D')
            smi.append(''.join(str(x) for x in neighbors))

        if hydrogens:
            smi.append(';H')
            smi.append(''.join(str(x) for x in hydrogens))

        if rings:
            smi.append(';r')
            smi.append(''.join(str(x) for x in rings))

        if hybridization:
            smi.append(';Z')
            smi.append(''.join(hybridization_str[x] for x in hybridization))

        if heteroatoms:
            smi.append(';W')
            smi.append(''.join(str(x) for x in heteroatoms))

        if self._radicals[n]:
            smi.append(';*')

        if charge:
            smi.append(charge_str[charge])

        smi.append(']')
        return ''.join(smi)

    def _format_bond(self: 'QueryContainer', n, m, **kwargs):
        return ','.join(order_str[x] for x in self._bonds[n][m].order)


__all__ = ['MoleculeSmiles', 'CGRSmiles', 'QuerySmiles']
