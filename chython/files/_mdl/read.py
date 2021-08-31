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
from base64 import urlsafe_b64encode
from fileinput import FileInput
from io import StringIO, TextIOWrapper
from itertools import islice
from os.path import abspath, join
from pathlib import Path
from pickle import dump, load, UnpicklingError
from sys import platform
from tempfile import gettempdir
from typing import Union, Iterator, List
from .parser import parse_error
from .stereo import MDLStereo
from ...containers import ReactionContainer, MoleculeContainer


class MDLReadMeta(type):
    def __call__(cls, *args, **kwargs):
        if kwargs.get('indexable'):
            _cls = type(cls.__name__, (cls,), {'__len__': lambda x: len(x._shifts) - 1, '__module__': cls.__module__})
            obj = object.__new__(_cls)
        else:
            obj = object.__new__(cls)
        obj.__init__(*args, **kwargs)
        return obj


class MDLRead(MDLStereo, metaclass=MDLReadMeta):
    def __init__(self, file, **kwargs):
        if isinstance(file, str):
            self._file = open(file)
            self._is_buffer = False
        elif isinstance(file, Path):
            self._file = file.open()
            self._is_buffer = False
        elif isinstance(file, (TextIOWrapper, StringIO, FileInput)):
            self._file = file
            self._is_buffer = True
        else:
            raise TypeError('invalid file. TextIOWrapper, StringIO subclasses possible')
        super().__init__(**kwargs)

    def close(self, force=False):
        """
        Close opened file

        :param force: force closing of externally opened file or buffer
        """
        if not self._is_buffer or force:
            self._file.close()

    def __enter__(self):
        return self

    def __exit__(self, _type, value, traceback):
        self.close()

    def _load_cache(self):
        """
        Load existing cache or create new. Working only for UNIX-like systems and local files (not buffers).
        """
        if platform == 'win32' or self._is_buffer:
            return
        try:
            with open(self.__cache_path, 'rb') as f:
                self._shifts = load(f)
        except FileNotFoundError:  # cache not found
            self.reset_index()
        except IsADirectoryError as e:
            raise IsADirectoryError(f'Please delete {self.__cache_path} directory') from e
        except (UnpicklingError, EOFError) as e:  # invalid file. ask user to check it.
            raise UnpicklingError(f'Invalid cache file {self.__cache_path}. Please delete it') from e

    def reset_index(self):
        """
        Create (rewrite) indexation table. Implemented only for object that
        is a real file (the path to the file is specified) because the external grep utility is used.
        """
        if platform != 'win32' and not self._is_buffer:
            self._shifts = self._get_shifts(self._file.name)
            with open(self.__cache_path, 'wb') as f:
                dump(self._shifts, f)
        else:
            raise self._implement_error

    @property
    def __cache_path(self):
        return abspath(join(gettempdir(), 'chython_' + urlsafe_b64encode(abspath(self._file.name).encode()).decode()))

    def read(self) -> List[Union[ReactionContainer, MoleculeContainer]]:
        """
        Parse whole file

        :return: list of parsed molecules
        """
        return list(iter(self))

    def __iter__(self) -> Iterator[Union[ReactionContainer, MoleculeContainer]]:
        return (x for x in self._data if not isinstance(x, parse_error))

    def __next__(self) -> Union[ReactionContainer, MoleculeContainer]:
        return next(iter(self))

    def __getitem__(self, item) -> Union[ReactionContainer, MoleculeContainer, parse_error]:
        """
        Getting the item by index from the original file,
        For slices records with errors skipped.
        For indexed access records with errors returned as error container.
        :return: [Molecule, Reaction]Container or list of [Molecule, Reaction]Containers
        """
        if self._shifts:
            _len = len(self._shifts) - 1
            if isinstance(item, int):
                if item >= _len or item < -_len:
                    raise IndexError('List index out of range')
                if item < 0:
                    item += _len
                self.seek(item)
                return next(self._data)
            elif isinstance(item, slice):
                start, stop, step = item.indices(_len)
                if start == stop:
                    return []
                if step == 1:
                    self.seek(start)
                    records = [x for x in islice(self._data, stop - start) if not isinstance(x, parse_error)]
                else:
                    records = []
                    for index in range(start, stop, step):
                        self.seek(index)
                        record = next(self._data)
                        if not isinstance(record, parse_error):
                            records.append(record)
                return records
            else:
                raise TypeError('Indices must be integers or slices')
        raise self._implement_error

    def read_text(self, item):
        """
        Read record block as text
        """
        if self._shifts:
            if not isinstance(item, int):
                raise TypeError('int required')
            _len = len(self._shifts) - 1
            if item >= _len or item < -_len:
                raise IndexError('List index out of range')
            if item < 0:
                item += _len
            start = self._shifts[item]
            end = self._shifts[item + 1]
            current = self._file.tell()
            self._file.seek(start)
            data = self._file.read(end - start)
            self._file.seek(current)
            return data
        raise self._implement_error

    def _prepare_meta(self, meta):
        new_meta = {}
        for k, v in meta.items():
            if v:
                new_meta[k] = '\n'.join(v)
            else:
                self._info(f'invalid metadata entry: {k}: {v}')
        return new_meta

    _shifts = None
    _implement_error = NotImplementedError('Indexable supported in unix-like o.s. and for files stored on disk')


__all__ = ['MDLRead']
