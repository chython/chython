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
from .emol import EMOLRead
from ...exceptions import EmptyMolecule


class ERXNRead:
    def __init__(self, line, ignore=False, log_buffer=None):
        tmp = line[13:].split()
        self.__reactants_count = int(tmp[0])
        self.__products_count = int(tmp[1])
        self.__reagents_count = int(tmp[2]) if len(tmp) == 3 else 0

        self.__reactants = []
        self.__products = []
        self.__reagents = []
        self.__ignore = ignore
        if log_buffer is None:
            log_buffer = []
        self.__log_buffer = log_buffer

    def __call__(self, line):
        if self.__empty_skip:
            if not line.startswith('M  V30 END CTAB'):
                return
            self.__empty_skip = False
        elif self.__in_mol:
            try:
                x = self.__parser(line)
            except EmptyMolecule:
                if not self.__ignore:
                    raise
                self.__empty_skip = True
                self.__in_mol -= 1
                if self.__in_mol:
                    self.__parser = EMOLRead(self.__ignore, self.__log_buffer)
                self.__log_buffer.append('empty molecule ignored')
            else:
                if x:
                    x = self.__parser.getvalue()
                    self.__in_mol -= 1
                    if self.__in_mol:
                        self.__parser = EMOLRead(self.__ignore, self.__log_buffer)
                    if self.__parser_group == 'REACTANT':
                        self.__reactants.append(x)
                    elif self.__parser_group == 'PRODUCT':
                        self.__products.append(x)
                    elif self.__parser_group == 'AGENT':
                        self.__reagents.append(x)
        elif self.__rend:
            raise ValueError('parser closed')
        elif line.startswith('M  V30 END'):
            if self.__parser_group != line[11:].strip():
                raise ValueError('invalid CTAB')
        elif line.startswith('M  V30 BEGIN'):
            x = line[13:].strip()
            if x == 'REACTANT':
                self.__in_mol = self.__reactants_count
            elif x == 'PRODUCT':
                self.__in_mol = self.__products_count
            elif x == 'AGENT':
                self.__in_mol = self.__reagents_count
            else:
                raise ValueError('invalid RXN CTAB')
            self.__parser_group = x
            if self.__in_mol:
                self.__parser = EMOLRead(self.__ignore, self.__log_buffer)
        elif line.startswith('M  END'):
            self.__rend = True
            return True
        else:
            raise ValueError('invalid CTAB')

    def getvalue(self):
        if self.__rend:
            return {'reactants': self.__reactants, 'products': self.__products, 'reagents': self.__reagents, 'meta': {}}
        raise ValueError('reaction not complete')

    __parser_group = __parser = None
    __rend = __empty_skip = False
    __in_mol = 0


__all__ = ['ERXNRead']
