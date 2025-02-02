# -*- coding: utf-8 -*-
#
#  Copyright 2017-2024 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from heapq import heappop, heappush
from itertools import product
from random import random
from typing import Callable, Optional, Tuple, TYPE_CHECKING, Union


if TYPE_CHECKING:
    from chython import MoleculeContainer, CGRContainer, QueryContainer
    from chython.containers.graph import Graph

charge_str = {-4: '-4', -3: '-3', -2: '-2', -1: '-', 0: '0', 1: '+', 2: '+2', 3: '+3', 4: '+4'}
order_str = {1: '-', 2: '=', 3: '#', 4: ':', 8: '~', None: '.'}
organic_set = {'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I', 'B'}
hybridization_str = {4: '4', 3: '1', 2: '2', 1: '3', None: 'n'}
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

    @property
    def smiles(self):
        """
        Generate SMILES string of the molecule.
        """
        return str(self)

    @cached_method
    def __str__(self):
        smiles, order = self._smiles(self._smiles_order(), _return_order=True)
        if (cx := self._format_cxsmiles(order)) is not None:
            smiles.append(' ')
            smiles.append(cx)

        self.__dict__['smiles_atoms_order'] = tuple(order)  # cache smiles_atoms_order
        return ''.join(smiles)

    def __format__(self: Union['Graph', 'Smiles'], format_spec, *, _return_order=False):
        """
        Signature generation options.

        :param format_spec: String with keys:
            a - Generate asymmetric closures.
            !s - Disable stereo marks.
            A - Use aromatic bonds instead aromatic atoms.
            m - Set atom mapping.
            r - Generate random-ordered smiles.
            h - Show implicit hydrogens.
            !b - Disable bonds tokens.
            !x - Disable CXSMILES extension.
            !z - Disable charge representation.

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
            if 'h' in format_spec:
                kwargs['hydrogens'] = True
            if '!b' in format_spec:
                kwargs['bonds'] = False
            if '!z' in format_spec:
                kwargs['charges'] = False

            if 'r' in format_spec:
                kwargs['random'] = True

                def w(_):
                    return random()
            else:
                w = self._smiles_order('!s' not in format_spec)

            smiles, order = self._smiles(w, _return_order=True, **kwargs)
            if _return_order:
                return ''.join(smiles), order
            elif '!x' in format_spec or (cx := self._format_cxsmiles(order)) is None:
                return ''.join(smiles)
            else:
                smiles.append(' ')
                smiles.append(cx)
                return ''.join(smiles)

        elif _return_order:
            smiles, order = self._smiles(self._smiles_order(), _return_order=True)
            smiles = ''.join(smiles)
            if (cx := self._format_cxsmiles(order)) is not None:  # cache molecule smiles
                self.__dict__['__cached_method___str__'] = f'{smiles} {cx}'
            else:
                self.__dict__['__cached_method___str__'] = smiles
            self.__dict__['smiles_atoms_order'] = tuple(order)  # cache smiles_atoms_order
            return smiles, order
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
        smiles, order = self._smiles(self._smiles_order(), _return_order=True)
        if (cx := self._format_cxsmiles(order)) is not None:
            smiles.append(' ')
            smiles.append(cx)
        self.__dict__['__cached_method___str__'] = ''.join(smiles)  # cache smiles
        return tuple(order)

    def _smiles(self: Union['Graph', 'Smiles'], weights, *, asymmetric_closures=False,
                open_parenthesis='(', close_parenthesis=')', delimiter='.', _return_order=False, **kwargs):
        if not self._atoms:
            return []
        bonds = self._bonds
        atoms_set = set(self._atoms)
        seen = {}
        cycle = 0
        casted_cycles = {}
        string = []
        order = []
        visited_bond = set()
        heap = list(range(1, 100))

        if kwargs.get('random', False):
            mod_weights_start = mod_weights = weights
        else:
            groups = defaultdict(int)
            for n in atoms_set:
                groups[weights(n)] -= 1

            def mod_weights_start(x):
                return (groups[weights(x)],  # common groups
                        weights(x))  # smallest weight

            def mod_weights(x):
                return (groups[weights(x)],  # common groups
                        weights(x),  # smallest weight
                        seen[x])  # BFS nearest to starting

        while True:
            start = min(atoms_set, key=mod_weights_start)
            if not kwargs.get('random', False):
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
                        cycle += 1
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

            # get order of each atom in ring closures
            rings_order = {token: n for n, token in enumerate(smiles) if token in tokens}
            for token in rings_order:  # prepare closure numbers
                released = []
                # order closure atoms as in smiles string
                for _, c in sorted(tokens[token], key=lambda x: rings_order[x[0]]):
                    if c in casted_cycles:  # release ring closure number
                        released.append(casted_cycles[c])
                    else:
                        casted_cycles[c] = heappop(heap)
                for c in released:  # delayed release to avoid duplicates. e.g. C1..C11..C1 instead of C1..C12..C2
                    heappush(heap, c)

            # prepare new neighbors order for stereo sign calculation
            for token in smiles:
                if token in tokens:
                    tokens[token].sort(key=lambda x: casted_cycles[x[1]])  # order closures
                    visited[token].extend(n for n, _ in tokens[token])
                if token in edges:
                    visited[token].extend(edges[token])

            for token in smiles:
                if isinstance(token, int):  # atoms
                    string.append(self._format_atom(token, visited, **kwargs))
                    order.append(token)
                    if token in tokens:
                        for m, c in tokens[token]:
                            if asymmetric_closures:
                                if (token, m) not in visited_bond:
                                    string.append(self._format_bond(token, m, visited, **kwargs))
                                    visited_bond.add((m, token))
                            else:
                                string.append(self._format_bond(token, m, visited, **kwargs))
                            string.append(self._format_closure(casted_cycles[c]))
                elif token == '(':
                    string.append(open_parenthesis)
                elif token == ')':
                    string.append(close_parenthesis)
                else:  # bonds
                    string.append(self._format_bond(*token, visited, **kwargs))

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

    def _smiles_order(self: 'Graph', stereo=True) -> Callable:
        return self.atoms_order.__getitem__

    def _format_cxsmiles(self, order) -> Optional[str]:
        ...


class MoleculeSmiles(Smiles):
    __slots__ = ()

    def sticky_smiles(self: Union['MoleculeContainer', 'MoleculeSmiles'], left: int, right: int = None, *,
                      remove_left: bool = False, remove_right: bool = False, tries: int = 10):
        """
        Generate smiles with fixed left and optionally right terminal atoms.
        Note: Produce expected results only with acyclic terminal atoms.

        :param remove_left: drop terminal atom and corresponding bond
        :param remove_right: drop terminal atom and corresponding bond
        :param tries: number of attempts to generate smiles
        """
        bonds = self._bonds
        bonds[left]  # noqa. check left atom availability
        if right:
            assert tries > 0, 'tries count should be positive'
            assert len(bonds[right]) == 1, 'right atom should be terminal'
            assert left != right, 'left and right atoms the same'
            assert self.connected_components_count == 1, 'only single component structures supported'

            seen = {right: 0}
            queue = [(right, -10)]
            while queue:
                n, d = queue.pop(0)
                for m in bonds[n].keys() - seen.keys():
                    queue.append((m, d - 10))
                    seen[m] = d
            seen[left] = -1_000_000_000  # prioritize left atom

            for _ in range(tries):
                smiles, order = self._smiles(lambda x: seen[x] + random(), _return_order=True, random=True)
                if order[-1] == right:
                    break
            else:
                raise Exception('generation of smiles failed')
            if remove_left:
                smiles = smiles[2:]
            if remove_right:
                smiles = smiles[:-2]
        else:
            smiles = self._smiles(lambda x: x != left, random=True)
            if remove_left:
                smiles = smiles[2:]
        return ''.join(smiles)

    def _smiles_order(self: 'MoleculeContainer', stereo=True) -> Callable:
        if stereo:
            return self._chiral_morgan.__getitem__
        else:
            return self.atoms_order.__getitem__

    def _format_cxsmiles(self: 'MoleculeContainer', order):
        if self.is_radical:
            radical = self._radicals
            return f'|^1:{",".join(str(n) for n, m in enumerate(order) if radical[m])}|'
        return

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
            elif charge and kwargs.get('charges', True):
                smi[5] = charge_str[charge]
        elif charge and kwargs.get('charges', True):
            smi[5] = charge_str[charge]

        if any(smi) or atom.atomic_symbol not in organic_set or self._radicals[n] or kwargs.get('hydrogens', False):
            smi[0] = '['
            smi[-1] = ']'
            if ih == 1:
                smi[4] = 'H'
            elif ih:
                smi[4] = f'H{ih}'
        elif hyb == 4 and ih and atom.atomic_number in (5, 7, 15):  # pyrrole
            smi[0] = '['
            smi[-1] = ']'
            if ih == 1:
                smi[4] = 'H'
            else:
                smi[4] = f'H{ih}'
        elif not ih and atom.atomic_number in (5, 6, 15, 16) and not self.not_special_connectivity[n]:
            # elemental B, C, P, S
            smi[0] = '['
            smi[-1] = ']'
        elif ih and atom.atomic_number == 15 and hyb != 1:
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
        if atom.atomic_symbol in ('R', 'X'):
            # fix markush representation
            smi[1], smi[2] = smi[2], smi[1]
        return ''.join(smi)

    def _format_bond(self: 'MoleculeContainer', n, m, adjacency, **kwargs):
        if not kwargs.get('bonds', True):
            return ''
        bonds = self._bonds
        order = bonds[n][m].order
        if order == 4:
            if kwargs.get('aromatic', True):
                return ''
            return ':'
        elif order == 1:  # cis-trans /\
            if kwargs.get('aromatic', True) and self.hybridization(n) == self.hybridization(m) == 4:
                return '-'
            if kwargs.get('stereo', True):
                if 'cache' in adjacency:
                    ct_map = adjacency['cache']
                else:
                    ct_map = adjacency['cache'] = self.__ct_map(adjacency)

                if (x := ct_map.get((n, m))) is not None:
                    return '/' if x else '\\'
            return ''
        elif order == 2:
            return '='
        elif order == 3:
            return '#'
        else:  # order == 8
            return '~'

    def __ct_map(self, adjacency):
        ct_map = {}
        cts = self._cis_trans_stereo
        if not cts:
            return ct_map
        ctt = self._stereo_cis_trans_terminals
        sct = self._stereo_cis_trans
        ctc = self._stereo_cis_trans_counterpart

        seen = set()
        for k, vs in adjacency.items():
            seen.add(k)
            if (ts := ctt.get(k)) and ts in cts:
                env = sct[ts]
                for v in vs:
                    if v in env:
                        if (k, v) in ct_map:
                            continue
                        elif x := ct_map.get(k):  # second substituent of C=
                            s = ct_map[(k, x)]
                            ct_map[(k, v)] = not s  # X/C(/R)=, C(\X)(/R)=, C(=C(\X)/R)=C=
                            ct_map[(v, k)] = s
                            if y := ctt.get(v):  # =C(\X)/R=, C(\X)(/R=)=
                                ct_map[v] = k
                                seen.add(y)
                        elif ts in seen:
                            o = ctc[k]
                            on = ct_map[o]
                            s = ct_map[(o, on)]
                            if not self._translate_cis_trans_sign(k, o, v, on):
                                s = not s
                            ct_map[(k, v)] = s
                            ct_map[k] = v
                            ct_map[(v, k)] = not s  # C/R=, R\1...C/1
                            if y := ctt.get(v):
                                ct_map[v] = k
                                seen.add(y)
                        else:  # left entry to double bond
                            if y := ctt.get(v):  # 1,3-diene case
                                ct_map[v] = k
                                seen.add(y)
                            ct_map[(v, k)] = True  # R/C=, C\1=...R/1, C(/R=)=, C(=C(/R=))=C=
                            ct_map[(k, v)] = False  # first DOWN
                            ct_map[k] = v
                seen.add(ts)
        return ct_map


class CGRSmiles(Smiles):
    __slots__ = ()

    def _format_atom(self: 'CGRContainer', n, adjacency, **kwargs):
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

    def _format_bond(self: 'CGRContainer', n, m, adjacency, **kwargs):
        bond = self._bonds[n][m]
        return dyn_order_str[(bond.order, bond.p_order)]


class QuerySmiles(Smiles):
    __slots__ = ()

    def _format_cxsmiles(self: 'QueryContainer', order):
        hybridization = self._hybridizations
        heteroatoms = self._heteroatoms
        masked = self._masked
        radical = self._radicals

        hh = ['atomProp']
        cx = []
        if any(radical.values()):
            cx.append(f'^1:{",".join(str(n) for n, m in enumerate(order) if radical[m])}')

        for n, m in enumerate(order):
            if len(hb := hybridization[m]) > 1 or (hb and hb[0] != 4):
                hh.append(f'{n}.hyb.{"".join(hybridization_str[x] for x in hb)}')
            if ha := heteroatoms[m]:
                hh.append(f'{n}.het.{"".join(str(x) for x in ha)}')
            if masked[m]:
                hh.append(f'{n}.msk.1')
        if len(hh) > 1:
            cx.append(':'.join(hh))
        if cx:
            return f'|{",".join(cx)}|'

    def _format_atom(self: 'QueryContainer', n, adjacency, **kwargs):
        atom = self._atoms[n]
        charge = self._charges[n]
        hybridization = self._hybridizations[n]
        neighbors = self._neighbors[n]
        hydrogens = self._hydrogens[n]
        rings = self._rings_sizes[n]

        if atom.isotope:
            smi = ['[', str(atom.isotope), atom.atomic_symbol]
        else:
            smi = ['[', atom.atomic_symbol]

        if n in self._atoms_stereo:  # mark atom as chiral. it's too difficult to set correct sign
            smi.append(';@?')
        if n in self._allenes_stereo:
            smi.append(';@?')

        if charge:
            smi.append(';')
            smi.append(charge_str[charge])

        if hydrogens:  # h<n> 	implicit-H-count 	<n> implicit hydrogens
            smi.append(';')
            smi.append(','.join(f'h{x}' for x in hydrogens))

        if neighbors:  # D<n> 	degree 	<n> explicit connections
            smi.append(';')
            smi.append(','.join(f'D{x}' for x in neighbors))

        if rings:
            smi.append(';')
            if rings[0]:
                smi.append(','.join(f'r{x}' for x in rings))
            else:
                smi.append('!R')

        if len(hybridization) == 1 and hybridization[0] == 4:  # only aromatic. other marks in cx extension
            smi.append(';a')

        smi.append(']')
        return ''.join(smi)

    def _format_bond(self: 'QueryContainer', n, m, adjacency, **kwargs):
        # bond chirality skipped. too difficult to implement.
        b = self._bonds[n][m]
        s = ','.join(order_str[x] for x in b.order)
        if (c := b.in_ring) is not None:
            s += ';@' if c else ';!@'
        return s


__all__ = ['MoleculeSmiles', 'CGRSmiles', 'QuerySmiles']
