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
from .mol import MOLRead
from ...exceptions import EmptyMolecule


class RXNRead:
    def __init__(self, line, ignore=False, log_buffer=None):
        self.__reactants_count = int(line[:3])
        self.__products_count = int(line[3:6]) + self.__reactants_count
        self.__reagents_count = int(line[6:].rstrip() or 0) + self.__products_count
        self.__molecules = []
        self.__ignore = ignore
        if log_buffer is None:
            log_buffer = []
        self.__log_buffer = log_buffer

    def __call__(self, line):
        if self.__parser:
            if self.__parser(line):
                self.__im = 4
                mol = self.__parser.getvalue()
                if self.__title:
                    mol['title'] = self.__title
                self.__molecules.append(mol)
                self.__parser = None
                if len(self.__molecules) == self.__reagents_count:
                    self.__rend = True
                    return True
        elif self.__empty_skip:
            if not line.startswith('$MOL'):
                return
            self.__empty_skip = False
            self.__im = 3
        elif self.__rend:
            raise ValueError('parser closed')
        elif self.__im == 4:
            if not line.startswith('$MOL'):
                raise ValueError('invalid RXN')
            self.__im = 3
        elif self.__im:
            if self.__im == 3:
                self.__title = line.strip()
            self.__im -= 1
        else:
            try:
                self.__parser = MOLRead(line, self.__log_buffer)
            except EmptyMolecule:
                if not self.__ignore:
                    raise
                self.__log_buffer.append('empty molecule ignored')
                if len(self.__molecules) < self.__reactants_count:
                    self.__reactants_count -= 1
                    self.__products_count -= 1
                    self.__reagents_count -= 1
                elif len(self.__molecules) < self.__products_count:
                    self.__products_count -= 1
                    self.__reagents_count -= 1
                elif len(self.__molecules) < self.__reagents_count:
                    self.__reagents_count -= 1
                    if len(self.__molecules) == self.__reagents_count:  # empty molecule is last in list
                        self.__rend = True
                        return True
                self.__empty_skip = True

    def getvalue(self):
        if self.__rend:
            return {'reactants': self.__molecules[:self.__reactants_count],
                    'products': self.__molecules[self.__reactants_count:self.__products_count],
                    'reagents': self.__molecules[self.__products_count:self.__reagents_count],
                    'meta': {}}
        raise ValueError('reaction not complete')

    __parser = None
    __empty_skip = __rend = False
    __im = 4


__all__ = ['RXNRead']
