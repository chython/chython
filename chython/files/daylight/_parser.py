# -*- coding: utf-8 -*-
#
#  Copyright 2018-2021 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2021 Aleksandr Sizov <murkyrussian@gmail.com>
#  Copyright 2019 Artem Mukanov <nostro32@mail.ru>
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
from collections import defaultdict
from re import split, compile, fullmatch, findall, search
from typing import Union, List
from .._mdl import Parser, parse_error
from ...containers import MoleculeContainer, CGRContainer, QueryContainer, ReactionContainer
from ...exceptions import IncorrectSmiles


# -,= OR bonds supported
# @ bond not supported
# @;!: and any other complicated combinations not supported


# tokens structure:
# (type: int, value)
# types:
# 0: atom
# 1: bond
# 2: open chain (
# 3: close chain )
# 4: dot bond .
# 5: in bracket raw data []
# 6: closure number
# 7: raw closure number
# 8: aromatic atom
# 9: up down bond
# 10: dynamic bond
# 11: dynamic atom
# 12: dynamic aromatic atom
# 13: query bond
# 14: query atom

replace_dict = {'-': 1, '=': 2, '#': 3, ':': 4, '~': 8, '.': None, '(': 2, ')': 3}
charge_dict = {'+': 1, '+1': 1, '++': 2, '+2': 2, '+3': 3, '+++': 3, '+4': 4, '++++': 4,
               '-': -1, '-1': -1, '--': -2, '-2': -2, '-3': -3, '---': -3, '-4': -4, '----': -4}

atom_re = compile(r'([1-9][0-9]{0,2})?([A-IK-PR-Zacnopsbt][a-ik-pr-vy]?)(@@|@)?(H[1-4]?)?([+-][1-4+-]?)?(:[0-9]{1,4})?')


class DaylightParser(Parser, ABC):
    def __init__(self):
        super().__init__()

    def __data(self):
        file = self._file
        parse = self.parse
        try:
            seekable = file.seekable()
        except AttributeError:
            seekable = False
        pos = file.tell() if seekable else None
        for n, line in enumerate(self.__file):
            try:
                x = parse(line)
            except ValueError:
                yield parse_error(n, pos, self._format_log(), line)
            else:
                yield x
            if seekable:
                pos = file.tell()

    @classmethod
    def create_parser(cls, header=None, ignore_stereo=False, *args, **kwargs):
        """
        Create SMILES parser function configured same as SMILESRead object.
        """
        obj = object.__new__(cls)
        obj._SMILESRead__header = header
        obj._SMILESRead__ignore_stereo = ignore_stereo
        super(DaylightParser, obj).__init__(*args, **kwargs)
        return obj.parse

    def close(self, force=False):
        """
        Close opened file.

        :param force: Force closing of externally opened file or buffer.
        """
        if not self.__is_buffer or force:
            self._file.close()

    def __enter__(self):
        return self

    def __exit__(self, _type, value, traceback):
        self.close()

    def read(self) -> List[Union[MoleculeContainer, CGRContainer, ReactionContainer]]:
        """
        Parse whole file.

        :return: List of parsed molecules or reactions.
        """
        return list(iter(self))

    def __iter__(self):
        return (x for x in self._data if not isinstance(x, parse_error))

    def __next__(self):
        return next(iter(self))

    @abstractmethod
    def parse(self, string: str) -> Union[MoleculeContainer, CGRContainer, QueryContainer, ReactionContainer]:
        ...

    def _parse_tokens(self, tokens):
        strong_cycle = not self._ignore
        t1 = tokens[0][0]
        if t1 == 2:
            if tokens[1][0] not in (0, 8, 11, 12):
                raise IncorrectSmiles('not atom started')
        elif t1 not in (0, 8, 11, 12):
            raise IncorrectSmiles('not atom started')

        atoms = []
        bonds = []
        order = defaultdict(list)
        atoms_types = []
        atom_num = 0
        last_num = 0
        stack = []
        cycles = {}
        used_cycles = set()
        cgr = []
        query = []
        stereo_bonds = defaultdict(dict)
        stereo_atoms = {}
        hydrogens = {}
        previous = None

        for token_type, token in tokens:
            if token_type == 2:  # ((((((
                if previous:
                    if previous[0] != 4:
                        raise IncorrectSmiles('bond before side chain')
                    previous = None
                stack.append(last_num)
            elif token_type == 3:  # ))))))
                if previous:
                    raise IncorrectSmiles('bond before closure')
                try:
                    last_num = stack.pop()
                except IndexError:
                    raise IncorrectSmiles('close chain more than open')
            elif token_type in (1, 4, 9, 10, 13):  # bonds. only keeping for atoms connecting
                if previous:
                    raise IncorrectSmiles('2 bonds in a row')
                elif not atoms:
                    raise IncorrectSmiles('started from bond')
                previous = (token_type, token)
            elif token_type == 6:  # cycle
                if previous and previous[0] == 4:
                    raise IncorrectSmiles('dot-cycle pattern invalid')
                elif token not in cycles:
                    if token in used_cycles:
                        if strong_cycle:
                            raise IncorrectSmiles('reused closure number')
                        else:
                            self._info(f'reused closure number: {token}')
                    else:
                        used_cycles.add(token)
                    cycles[token] = (last_num, previous, len(order[last_num]))
                    order[last_num].append(None)  # Reserve a table
                else:
                    a, ob, ind = cycles[token]
                    if ob:
                        if not previous:
                            bt, b = ob
                            if bt == 9:  # closure open is \/ bonded
                                stereo_bonds[a][last_num] = b
                                bt = b = 1
                            elif strong_cycle:
                                raise IncorrectSmiles('not equal cycle bonds')
                        else:
                            bt, b = previous
                            obt, ob = ob
                            if bt == 9:  # \/ bonds can be unequal
                                if obt == 9:
                                    stereo_bonds[a][last_num] = ob
                                elif ob != 1:
                                    raise IncorrectSmiles('not equal cycle bonds')
                                stereo_bonds[last_num][a] = b
                                bt = b = 1
                            elif obt == 9:
                                if b != 1:
                                    raise IncorrectSmiles('not equal cycle bonds')
                                stereo_bonds[a][last_num] = ob
                            elif b != ob:
                                raise IncorrectSmiles('not equal cycle bonds')
                    elif previous:
                        bt, b = previous
                        if bt == 9:  # stereo \/
                            stereo_bonds[last_num][a] = b
                            bt = b = 1
                        elif strong_cycle:
                            raise IncorrectSmiles('not equal cycle bonds')
                    else:
                        bt = 1
                        b = 4 if atoms_types[last_num] in (8, 12) and atoms_types[a] in (8, 12) else 1

                    if bt == 1:
                        bonds.append((last_num, a, b))
                    elif bt == 13:
                        bonds.append((last_num, a, None))
                        query.append((last_num, a, b))
                    else:  # bt == 10
                        bonds.append((last_num, a, None))
                        cgr.append((last_num, a, b))
                    order[a][ind] = last_num
                    order[last_num].append(a)
                    del cycles[token]
                previous = None
            else:  # atom
                if atoms:
                    if not previous:
                        bt = 1
                        b = 4 if token_type in (8, 12) and atoms_types[last_num] in (8, 12) else 1
                        order[last_num].append(atom_num)
                        order[atom_num].append(last_num)
                    else:
                        bt, b = previous
                        if bt != 4:
                            order[last_num].append(atom_num)
                            order[atom_num].append(last_num)
                    if bt == 1:
                        bonds.append((atom_num, last_num, b))
                    elif bt == 9:
                        bonds.append((atom_num, last_num,
                                      4 if token_type in (8, 12) and atoms_types[last_num] in (8, 12) else 1))
                        stereo_bonds[last_num][atom_num] = b
                        stereo_bonds[atom_num][last_num] = not b
                    elif bt == 13:
                        bonds.append((atom_num, last_num, None))
                        query.append((atom_num, last_num, b))
                    elif bt == 10:
                        bonds.append((atom_num, last_num, None))
                        cgr.append((atom_num, last_num, b))

                if token_type not in (11, 12):
                    stereo = token.pop('stereo')
                    if stereo is not None:
                        stereo_atoms[atom_num] = stereo
                    hydrogen = token.pop('hydrogen')
                    if hydrogen is not None:
                        hydrogens[atom_num] = hydrogen

                atoms.append(token)
                atoms_types.append(token_type)

                last_num = atom_num
                atom_num += 1
                previous = None

        if stack:
            raise IncorrectSmiles('number of ( does not equal to number of )')
        elif cycles:
            raise IncorrectSmiles('cycle is not finished')
        elif previous:
            raise IncorrectSmiles('bond on the end')

        stereo_bonds = {n: ms for n, ms in stereo_bonds.items() if len(ms) == 1 or len(ms) == set(ms.values())}
        mol = {'atoms': atoms, 'bonds': bonds, 'order': order,
               'stereo_bonds': stereo_bonds, 'stereo_atoms': stereo_atoms, 'hydrogens': hydrogens, 'meta': {}}
        if cgr or any(x in (11, 12) for x in atoms_types):
            mol['cgr'] = cgr
        elif query or any(x == 14 for x in atoms_types):
            mol['query'] = query
        return mol


def _raw_tokenize(smiles):
    token_type = token = None
    tokens = []
    for s in smiles:
        if s == '[':  # open complex token
            if token_type == 5:  # two opened [
                raise IncorrectSmiles('[..[')
            elif token:
                tokens.append((token_type, token))
            elif token_type == 7:  # empty closure
                raise IncorrectSmiles('invalid closure')
            token = []
            token_type = 5
        elif s == ']':  # close complex token
            if token_type != 5:
                raise IncorrectSmiles(']..]')
            elif not token:
                raise IncorrectSmiles('empty [] brackets')
            tokens.append((5, ''.join(token)))
            token_type = token = None
        elif token_type == 5:  # grow token with brackets. skip validation
            token.append(s)
        elif s == '(':
            if token:
                tokens.append((token_type, token))
                token = None
            elif token_type == 7:  # empty closure
                raise IncorrectSmiles('invalid closure')
            elif token_type == 2:  # barely opened
                raise IncorrectSmiles('((')
            token_type = 2
            tokens.append((2, None))
        elif s == ')':
            if token:
                tokens.append((token_type, token))
                token = None
            elif token_type == 7:  # empty closure
                raise IncorrectSmiles('invalid closure')
            elif token_type == 2:  # barely opened
                raise IncorrectSmiles('()')
            token_type = 3
            tokens.append((3, None))
        elif s.isnumeric():  # closures
            if token_type == 7:  # % already found. collect number
                if not token and s == '0':
                    raise IncorrectSmiles('number starts with 0')
                token.append(s)
                if len(token) == 2:
                    tokens.append((token_type, token))
                    token_type = token = None
            else:
                if s == '0':
                    raise IncorrectSmiles('number starts with 0')
                elif token:
                    tokens.append((token_type, token))
                    token = None
                token_type = 6
                tokens.append((6, int(s)))
        elif s == '%':
            if token:
                tokens.append((token_type, token))
            elif token_type == 7:
                raise IncorrectSmiles('%%')
            token_type = 7
            token = []
        elif s in '=#:-~':  # bonds found
            if token_type == 13:
                token.append(replace_dict[s])
                tokens.append((13, token))
                token = None
            elif token:
                tokens.append((token_type, token))
                token = None
            else:
                token_type = 1
                tokens.append((1, replace_dict[s]))
        elif s in r'\/':
            if token:
                tokens.append((token_type, token))
                token = None
            token_type = 9
            tokens.append((9, s == '/'))  # Up is true
        elif s == '.':
            if token:
                tokens.append((token_type, token))
                token = None
            token_type = 4
            tokens.append((4, None))
        elif s in 'NOPSFI':  # organic atoms
            if token:
                tokens.append((token_type, token))
                token = None
            token_type = 0
            tokens.append((0, s))
        elif s in 'cnopsb':  # aromatic ring atom
            if token:
                tokens.append((token_type, token))
                token = None
            token_type = 8
            tokens.append((8, s.upper()))
        elif s in 'CB':  # flag possible Cl or Br
            if token:
                tokens.append((token_type, token))
            token_type = 0
            token = s
        elif s == ',':  # query bond separator
            if token_type != 1:
                raise IncorrectSmiles('SMARTS query bond invalid')
            token_type = 13
            token = [tokens.pop(-1)[1]]
        elif token_type == 0:
            if s == 'l':
                if token == 'C':
                    tokens.append((0, 'Cl'))
                    token = None
                else:
                    raise IncorrectSmiles('invalid element Bl')
            elif s == 'r':
                if token == 'B':
                    tokens.append((0, 'Br'))
                    token = None
                else:
                    raise IncorrectSmiles('invalid smiles for Cr')
            else:
                raise IncorrectSmiles('invalid smiles')
        else:
            raise IncorrectSmiles('invalid smiles')

    if token_type == 5:
        raise IncorrectSmiles('atom description has not finished')
    elif token_type == 2:
        raise IncorrectSmiles('not closed')
    elif token:
        tokens.append((token_type, token))  # %closure or C or B
    return [(6, int(''.join(y))) if x == 7 else (x, y) for x, y in tokens]  # composite closures folding


def _atom_parse(token):
    # [isotope]Element[element][@[@]][H[n]][+-charge][:mapping]
    match = fullmatch(atom_re, token)
    if match is None:
        raise IncorrectSmiles(f'atom token invalid {{{token}}}')
    isotope, element, stereo, hydrogen, charge, mapping = match.groups()

    if isotope:
        isotope = int(isotope)

    if stereo:
        stereo = stereo == '@'

    if hydrogen:
        if len(hydrogen) > 1:
            hydrogen = int(hydrogen[1:])
        else:
            hydrogen = 1
    else:
        hydrogen = 0

    if charge:
        try:
            charge = charge_dict[charge]
        except KeyError:
            raise IncorrectSmiles('charge token invalid')
    else:
        charge = 0

    if mapping:
        try:
            mapping = int(mapping[1:])
        except ValueError:
            raise IncorrectSmiles('invalid mapping token')
    else:
        mapping = 0

    if element in ('c', 'n', 'o', 'p', 's', 'as', 'se', 'b', 'te'):
        _type = 8
        element = element.capitalize()
    else:
        _type = 0
    return _type, {'element': element, 'charge': charge, 'isotope': isotope, 'is_radical': False,
                   'mapping': mapping, 'x': 0., 'y': 0., 'z': 0., 'hydrogen': hydrogen, 'stereo': stereo}

