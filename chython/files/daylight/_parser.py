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
from fileinput import FileInput
from io import StringIO, TextIOWrapper
from pathlib import Path
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


class DaylightParser(Parser, ABC):
    def __init__(self, file, **kwargs):
        if isinstance(file, str):
            self._file = open(file)
            self.__is_buffer = False
        elif isinstance(file, Path):
            self._file = file.open()
            self.__is_buffer = False
        elif isinstance(file, (TextIOWrapper, StringIO, FileInput)):
            self._file = file
            self.__is_buffer = True
        else:
            raise TypeError('invalid file. TextIOWrapper, StringIO subclasses possible')
        super().__init__(**kwargs)
        self.__file = iter(self._file.readline, '')
        self._data = self.__data()

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
