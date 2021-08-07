# -*- coding: utf-8 -*-
#
#  Copyright 2018-2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from ctypes import c_char, c_double, c_short, c_long, create_string_buffer, POINTER, Structure, cdll, byref
from distutils.util import get_platform
from fileinput import FileInput
from io import StringIO, TextIOWrapper
from os import name
from pathlib import Path
from re import split
from sys import prefix, exec_prefix
from traceback import format_exc
from typing import List, Dict, Union, Iterator
from warnings import warn
from ._mdl import Parser, common_isotopes, parse_error
from ..containers import MoleculeContainer


class INCHIRead(Parser):
    """
    INCHI separated per lines files reader. works similar to opened file object. support `with` context manager.
    on initialization accept opened in text mode file, string path to file,
    pathlib.Path object or another buffered reader object.
    line should be start with INCHI string and
    optionally continues with space/tab separated list of key:value [or key=value] data if header=None.
        example:
            InChI=1S/C2H5/c1-2/h1H2,2H3/q+1 id:123 key=value
    if header=True then first line of file should be space/tab separated list of keys including INCHI column key.
        example:
            ignored_inchi_key key1 key2
            InChI=1S/C2H5/c1-2/h1H2,2H3/q+1 1 2
    also possible to pass list of keys (without inchi_pseudo_key) for mapping space/tab separated list
    of INCHI and values: header=['key1', 'key2'] # order depended
    """
    def __init__(self, file, header=None, ignore_stereo=False, **kwargs):
        """
        :param ignore: Skip some checks of data or try to fix some errors.
        :param remap: Remap atom numbers started from one.
        :param store_log: Store parser log if exists messages to `.meta` by key `ParserLog`.
        :param calc_cis_trans: Calculate cis/trans marks from 2d coordinates.
        :param ignore_stereo: Ignore stereo data.
        """
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

        if header is True:
            self.__header = next(self.__file).split()[1:]
        elif header:
            if not isinstance(header, (list, tuple)) or not all(isinstance(x, str) for x in header):
                raise TypeError('expected list (tuple) of strings')
            self.__header = header
        else:
            self.__header = None

        self.__ignore_stereo = ignore_stereo
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
            x = parse(line)
            if isinstance(x, dict):
                yield parse_error(n, pos, self._format_log(), x)
                if seekable:
                    pos = file.tell()
            else:
                yield x

    @classmethod
    def create_parser(cls, *args, **kwargs):
        """
        Create INCHI parser function configured same as INCHIRead object
        """
        obj = object.__new__(cls)
        obj._INCHIRead__header = None
        obj._INCHIRead__ignore_stereo = False
        super(INCHIRead, obj).__init__(*args, **kwargs)
        return obj.parse

    def close(self, force=False):
        """
        close opened file

        :param force: force closing of externally opened file or buffer
        """
        if not self.__is_buffer or force:
            self._file.close()

    def __enter__(self):
        return self

    def __exit__(self, _type, value, traceback):
        self.close()

    def read(self) -> List[MoleculeContainer]:
        """
        parse whole file

        :return: list of parsed molecules
        """
        return list(iter(self))

    def __iter__(self) -> Iterator[MoleculeContainer]:
        return (x for x in self._data if not isinstance(x, parse_error))

    def __next__(self) -> MoleculeContainer:
        return next(iter(self))

    def parse(self, inchi: str) -> Union[MoleculeContainer, Dict[str, str]]:
        """
        convert INCHI string into MoleculeContainer object. string should be start with INCHI and
        optionally continues with space/tab separated list of key:value [or key=value] data.
        """
        self._flush_log()
        inchi, *data = inchi.split()
        if not inchi:
            self._info('empty inchi string')
            return {}
        elif self.__header is None:
            meta = {}
            for x in data:
                try:
                    k, v = split('[=:]', x, 1)
                    meta[k] = v
                except ValueError:
                    self._info(f'invalid metadata entry: {x}')
        else:
            meta = dict(zip(self.__header, data))

        try:
            record = self.__parse_inchi(inchi)
        except ValueError:
            self._info(f'string: {inchi}\nconsist errors:\n{format_exc()}')
            return meta

        record['meta'] = meta
        try:
            container = self._convert_molecule(record)
        except ValueError:
            self._info(f'record consist errors:\n{format_exc()}')
            return meta
        else:
            if self._store_log:
                log = self._format_log()
                if log:
                    container.meta['ParserLog'] = log
            return container

    @staticmethod
    def __parse_inchi(string):
        structure = INCHIStructure()
        if lib.GetStructFromINCHI(byref(InputINCHI(string)), byref(structure)):
            lib.FreeStructFromINCHI(byref(structure))
            raise ValueError('invalid INCHI')

        atoms, bonds = [], []
        seen = set()
        for n in range(structure.num_atoms):
            seen.add(n)
            atom = structure.atom[n]
            element = atom.elname.decode()

            isotope = atom.isotopic_mass
            if isotope in (0, 10000):
                isotope = None
            elif isotope > 10000:
                isotope = isotope - 10000 + common_isotopes[element]

            atoms.append({'element': element, 'charge': int.from_bytes(atom.charge, byteorder='big', signed=True),
                          'mapping': 0, 'x': atom.x, 'y': atom.y, 'z': atom.z, 'isotope': isotope,
                          'is_radical': bool(int.from_bytes(atom.radical, byteorder='big'))})

            for k in range(atom.num_bonds):
                m = atom.neighbor[k]
                if m in seen:
                    continue
                order = atom.bond_type[k]
                if order:
                    bonds.append((n, m, order))

        lib.FreeStructFromINCHI(byref(structure))
        return {'atoms': atoms, 'bonds': bonds}


class InputINCHI(Structure):
    def __init__(self, string, options=None):
        if options is None:
            options = create_string_buffer(1)
        else:
            options = create_string_buffer(' '.join(f'{opt_flag}{x}' for x in options).encode())
        super().__init__(create_string_buffer(string.encode()), options)

    _fields_ = [('szInChI', POINTER(c_char)),  # InChI ASCII string to be converted to a strucure
                ('szOptions', POINTER(c_char))  # InChI options: space-delimited; each is preceded
                                                # by '/' or '-' depending on OS and compiler
                ]


class Atom(Structure):
    _fields_ = [('x', c_double), ('y', c_double), ('z', c_double),  # atom coordinates
                ('neighbor', c_short * 20),  # adjacency list: ordering numbers of the adjacent atoms, >= 0
                ('bond_type', c_char * 20),  # inchi_BondType
                ('bond_stereo', c_char * 20),  # inchi_BondStereo2D; negative if the sharp end points to opposite atom
                ('elname', c_char * 6),  # zero-terminated chemical element name: "H", "Si", etc.
                ('num_bonds', c_short),  # number of neighbors, bond types and bond stereo in the adjacency list
                ('num_iso_H', c_char * 4),  # implicit hydrogen atoms
                                            # [0]: number of implicit non-isotopic H (exception: num_iso_H[0]=-1 means
                                            # INCHI adds implicit H automatically),
                                            # [1]: number of implicit isotopic 1H (protium),
                                            # [2]: number of implicit 2H (deuterium),
                                            # [3]: number of implicit 3H (tritium)
                ('isotopic_mass', c_short),  # 0 => non-isotopic; isotopic mass or 10000 + mass - average atomic mass
                ('radical', c_char),  # inchi_Radical,
                ('charge', c_char)  # positive or negative; 0 => no charge
                ]


class Stereo0D(Structure):
    _fields_ = [('neighbor', c_short * 4),  # 4 atoms always
                ('central_atom', c_short),  # central tetrahedral atom or a central atom of allene; otherwise NO_ATOM
                ('type', c_char),  # inchi_StereoType0D
                ('parity', c_char)  # inchi_StereoParity0D: may be a combination of two parities:
                                    # ParityOfConnected | (ParityOfDisconnected << 3), see Note above
                ]


class INCHIStructure(Structure):
    _fields_ = [('atom', POINTER(Atom)),  # array of num_atoms elements
                ('stereo0D', POINTER(Stereo0D)),  # array of num_stereo0D 0D stereo elements or NULL
                ('num_atoms', c_short),  # number of atoms in the structure
                ('num_stereo0D', c_short),  # number of 0D stereo elements
                ('szMessage', POINTER(c_char)),  # Error/warning ASCII message
                ('szLog', POINTER(c_char)),  # log-file ASCII string, contains a human-readable list
                                             # of recognized options and possibly an Error/warn message
                ('WarningFlags', (c_long * 2) * 2)
                ]


try:
    from site import getuserbase
except ImportError:
    prefixes = {prefix, exec_prefix}
else:
    user_prefix = getuserbase()
    if user_prefix:
        prefixes = {prefix, exec_prefix, user_prefix}
    else:
        prefixes = {prefix, exec_prefix}

sitepackages = []
for pr in prefixes:
    pr = Path(pr)
    if name == 'posix':
        sitepackages.append(pr / 'local' / 'lib')
    else:
        sitepackages.append(pr)
    sitepackages.append(pr / 'lib')

platform = get_platform()
if platform in ('linux-x86_64', 'win-amd64'):
    if platform == 'win-amd64':
        opt_flag = '/'
        libname = 'libinchi.dll'
    else:
        opt_flag = '-'
        libname = 'libinchi.so'

    for site in sitepackages:
        lib_path = site / libname
        if lib_path.exists():
            try:
                lib = cdll.LoadLibrary(str(lib_path))
            except OSError:
                warn('libinchi loading problem', ImportWarning)
                __all__ = []
                del INCHIRead
                break
            __all__ = ['INCHIRead']
            break
    else:
        warn('broken package installation. libinchi not found', ImportWarning)
        __all__ = []
        del INCHIRead
else:
    warn('unsupported platform for libinchi', ImportWarning)
    __all__ = []
    del INCHIRead
