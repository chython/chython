# -*- coding: utf-8 -*-
#
#  Copyright 2021-2023 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from abc import ABCMeta, abstractmethod
from base64 import urlsafe_b64encode
from fileinput import FileInput
from io import StringIO, TextIOWrapper
from itertools import islice
from os.path import abspath, join
from pathlib import Path
from pickle import load, UnpicklingError
from sys import platform
from tempfile import gettempdir
from typing import Union, Iterator, List, Dict, Optional
from ...containers import ReactionContainer, MoleculeContainer


class MDLReadMeta(ABCMeta):
    def __call__(cls, *args, **kwargs):
        if kwargs.get('indexable'):
            _cls = type(cls.__name__, (cls,), {'__len__': lambda x: len(x._shifts), '__module__': cls.__module__})
            obj = object.__new__(_cls)  # noqa
        else:
            obj = object.__new__(cls)  # noqa
        obj.__init__(*args, **kwargs)
        return obj


class MDLRead(metaclass=MDLReadMeta):
    def __init__(self, file, buffer_size=1000, indexable=False, ignore=True, remap=False, ignore_bad_isotopes=False,
                 ignore_stereo=False, calc_cis_trans=False):
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
            raise TypeError('invalid file. TextIOWrapper, StringIO or FileInput subclasses or path to file expected')
        self._shifts = None
        self._tell = 0
        self._buffer_size = buffer_size
        self._buffer = None
        self._ignore = ignore
        self._remap = remap
        self._ignore_bad_isotopes = ignore_bad_isotopes
        self._ignore_stereo = ignore_stereo
        self._calc_cis_trans = calc_cis_trans
        if indexable:
            self._load_cache()

    def read(self, amount: Optional[int] = None) -> List[Union[ReactionContainer, MoleculeContainer]]:
        """
        Parse whole file

        :param amount: number of records to read
        """
        if amount:
            return list(islice(iter(self), amount))
        return list(iter(self))

    @abstractmethod
    def read_structure(self, *, current: bool = True):
        """
        Read Reaction or Molecule container.

        :param current: return current structure if already parsed, otherwise read next
        """

    @abstractmethod
    def read_metadata(self, *, current: bool = True) -> Dict[str, str]:
        """
        Read metadata block
        """

    def close(self, force: bool = False):
        """
        Close opened file

        :param force: force closing of externally opened file or buffer
        """
        if not self._is_buffer or force:
            self._file.close()

    def tell(self):
        """
        Number of records processed from the original file
        """
        return self._tell

    def seek(self, offset):
        """
        Shift to a given record number
        """
        if self._shifts:
            if 0 <= offset < len(self._shifts):
                self._tell = offset
                self._buffer = None
                self._file.seek(self._shifts[offset])
            else:
                raise IndexError('invalid offset')
        else:
            raise NotImplementedError('Indexable supported in unix-like o.s. and for files stored on disk')

    @abstractmethod
    def reset_index(self):
        """
        Create (rewrite) indexation table. Implemented only for object that
        is a real file (the path to the file is specified) because the external grep utility is used.
        """

    def read_block(self, *, current: bool = True) -> str:
        """
        Read full record block with metadata
        """
        return ''.join(self._read_block(current=current))

    @abstractmethod
    def _read_block(self, *, current: bool = True) -> List[str]:
        """
        Read full record block with metadata
        """

    def _load_cache(self):
        """
        Load existing cache or create new. Working only for UNIX-like systems and local files (not buffers).
        """
        if platform == 'win32' or self._is_buffer:
            return
        try:
            with open(self._cache_path, 'rb') as f:
                self._shifts = load(f)
        except FileNotFoundError:  # cache not found
            self.reset_index()
        except IsADirectoryError as e:
            raise IsADirectoryError(f'Please delete {self._cache_path} directory') from e
        except (UnpicklingError, EOFError) as e:  # invalid file. ask user to check it.
            raise UnpicklingError(f'Invalid cache file {self._cache_path}. Please delete it') from e

    @property
    def _cache_path(self):
        return abspath(join(gettempdir(), 'chython_' + urlsafe_b64encode(abspath(self._file.name).encode()).decode()))

    def __enter__(self):
        return self

    def __exit__(self, _type, value, traceback):
        self.close()

    def __iter__(self) -> Iterator[Union[ReactionContainer, MoleculeContainer]]:
        while True:
            try:
                yield self.read_structure(current=False)
            except ValueError:
                pass
            except EOFError:
                return

    def __next__(self) -> Union[ReactionContainer, MoleculeContainer]:
        return next(iter(self))

    def __getitem__(self, item) -> Union[ReactionContainer, MoleculeContainer,
                                         List[Union[ReactionContainer, MoleculeContainer]]]:
        """
        Getting the item by index from the original file,
        For slices records with errors skipped.
        """
        if self._shifts:
            _len = len(self._shifts)
            if isinstance(item, int):
                if item >= _len or item < -_len:
                    raise IndexError('List index out of range')
                if item < 0:
                    item += _len
                self.seek(item)
                return self.read_structure()
            elif isinstance(item, slice):
                start, stop, step = item.indices(_len)
                if start == stop:
                    return []
                if step == 1:
                    self.seek(start)
                    records = []
                    for _ in range(start, stop):
                        try:
                            records.append(self.read_structure(current=False))
                        except EOFError:
                            break
                        except ValueError:
                            pass
                else:
                    records = []
                    for index in range(start, stop, step):
                        self.seek(index)
                        try:
                            records.append(self.read_structure())
                        except EOFError:
                            break
                        except ValueError:
                            pass
                return records
            else:
                raise TypeError('Indices must be integers or slices')
        raise NotImplementedError('Indexable supported in unix-like o.s. and for files stored on disk')


__all__ = ['MDLRead']
