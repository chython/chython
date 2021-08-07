# -*- coding: utf-8 -*-
#
#  Copyright 2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from io import StringIO, TextIOWrapper
from pathlib import Path


class _MDLWrite:
    def __init__(self, file, *, mapping: bool = True, append: bool = False):
        """
        :param mapping: write atom mapping.
        """
        self._mapping = mapping

        if isinstance(file, str):
            self._file = open(file, 'a' if append else 'w')
            self._is_buffer = False
        elif isinstance(file, Path):
            self._file = file.open('a' if append else 'w')
            self._is_buffer = False
        elif isinstance(file, (TextIOWrapper, StringIO)):
            self._file = file
            self._is_buffer = True
        else:
            raise TypeError('invalid file. '
                            'TextIOWrapper, StringIO, BytesIO, BufferedReader and BufferedIOBase subclasses possible')

    def close(self, force=False):
        """
        close opened file

        :param force: force closing of externally opened file or buffer
        """
        self.write = self.__write_closed

        if not self._is_buffer or force:
            self._file.close()

    def __enter__(self):
        return self

    def __exit__(self, _type, value, traceback):
        self.close()

    @staticmethod
    def __write_closed(_):
        raise ValueError('I/O operation on closed writer')


__all__ = ['_MDLWrite']
