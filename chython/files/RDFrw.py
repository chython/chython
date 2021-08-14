# -*- coding: utf-8 -*-
#
#  Copyright 2014-2021 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2019 Dinar Batyrshin <batyrshin-dinar@mail.ru>
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
from bisect import bisect_left
from collections import defaultdict
from itertools import chain
from os.path import getsize
from subprocess import check_output
from time import strftime
from traceback import format_exc
from typing import Union
from ._mdl import parse_error, MDLRead, MDLWrite, MOLRead, EMOLRead, RXNRead, ERXNRead, EMDLWrite
from ..containers import ReactionContainer, MoleculeContainer


class RDFRead(MDLRead):
    """
    MDL RDF files reader. works similar to opened file object. support `with` context manager.
    on initialization accept opened in text mode file, string path to file,
    pathlib.Path object or another buffered reader object
    """
    def __init__(self, file, indexable=False, **kwargs):
        """
        :param indexable: if True: supported methods seek, tell, object size and subscription, it only works when
            dealing with a real file (the path to the file is specified) because the external grep utility is used,
            supporting in unix-like OS the object behaves like a normal open file.

            if False: works like generator converting a record into ReactionContainer and returning each object in
            order, records with errors are skipped
        :param ignore: Skip some checks of data or try to fix some errors.
        :param remap: Remap atom numbers started from one.
        :param store_log: Store parser log if exists messages to `.meta` by key `ParserLog`.
        :param calc_cis_trans: Calculate cis/trans marks from 2d coordinates.
        :param ignore_stereo: Ignore stereo data.
        """
        super().__init__(file, **kwargs)
        self.__file = iter(self._file.readline, '')
        self._data = self.__reader()

        if indexable:
            if next(self._data):
                self._load_cache()
        else:
            next(self._data)

    @staticmethod
    def _get_shifts(file):
        shifts = [int(x.split(b':', 1)[0]) for x in
                  check_output(['grep', '-boE', r'^\$[RM]FMT', file]).split()]
        shifts.append(getsize(file))
        return shifts

    def seek(self, offset):
        """
        shifts on a given number of record in the original file
        :param offset: number of record
        """
        if self._shifts:
            if 0 <= offset < len(self._shifts):
                current_pos = self._file.tell()
                new_pos = self._shifts[offset]
                if current_pos != new_pos:
                    if current_pos == self._shifts[-1]:  # reached the end of the file
                        self._file.seek(0)
                        self.__file = iter(self._file.readline, '')
                        self._data = self.__reader()
                        next(self._data)
                        if offset:  # move not to the beginning of the file
                            self._file.seek(new_pos)
                            self._data.send(offset)
                            self.__already_seeked = True
                    else:
                        if not self.__already_seeked:
                            self._file.seek(new_pos)
                            self._data.send(offset)
                            self.__already_seeked = True
                        else:
                            raise BlockingIOError('File already seeked. New seek possible only after reading any data')
            else:
                raise IndexError('invalid offset')
        else:
            raise self._implement_error

    def tell(self):
        """
        :return: number of records processed from the original file
        """
        if self._shifts:
            t = self._file.tell()
            if t == self._shifts[0]:
                return 0
            elif t == self._shifts[-1]:
                return len(self._shifts) - 1
            elif t in self._shifts:
                return bisect_left(self._shifts, t)
            else:
                return bisect_left(self._shifts, t) - 1
        raise self._implement_error

    def __reader(self):
        record = parser = mkey = pos = None
        failed = False
        file = self._file
        try:
            seekable = file.seekable()
        except AttributeError:
            seekable = False

        if next(self.__file).startswith('$RXN'):  # parse RXN file
            is_reaction = True
            ir = 3
            meta = defaultdict(list)
            if seekable:
                pos = 0  # $RXN line starting position
            count = 0
            yield False
        elif next(self.__file).startswith('$DATM'):  # skip header
            ir = 0
            is_reaction = meta = None
            seek = yield True
            if seek is not None:
                yield
                count = seek - 1
                self.__already_seeked = False
            else:
                count = -1
        else:
            raise ValueError('invalid file')

        for line in self.__file:
            if failed and not line.startswith(('$RFMT', '$MFMT')):
                continue
            elif parser:
                try:
                    if parser(line):
                        record = parser.getvalue()
                        parser = None
                except ValueError:
                    parser = None
                    self._info(f'line:\n{line}\nconsist errors:\n{format_exc()}')
                    seek = yield parse_error(count, pos, self._format_log(), {})
                    if seek is not None:
                        yield
                        count = seek - 1
                        self.__already_seeked = False
                    else:
                        failed = True
            elif line.startswith('$RFMT'):
                if record:
                    record['meta'].update(self._prepare_meta(meta))
                    if title:
                        record['title'] = title
                    try:
                        if is_reaction:
                            container = self._convert_reaction(record)
                        else:
                            container = self._convert_molecule(record)
                    except ValueError:
                        self._info(f'record consist errors:\n{format_exc()}')
                        seek = yield parse_error(count, pos, self._format_log(), record['meta'])
                    else:
                        seek = yield container

                    record = None
                    if seek is not None:
                        yield
                        count = seek - 1
                        self.__already_seeked = False
                        continue

                if seekable:
                    pos = file.tell()  # $RXN line starting position
                count += 1
                is_reaction = True
                ir = 4
                failed = False
                mkey = None
                meta = defaultdict(list)
            elif line.startswith('$MFMT'):
                if record:
                    record['meta'].update(self._prepare_meta(meta))
                    if title:
                        record['title'] = title
                    try:
                        if is_reaction:
                            container = self._convert_reaction(record)
                        else:
                            container = self._convert_molecule(record)
                    except ValueError:
                        self._info(f'record consist errors:\n{format_exc()}')
                        seek = yield parse_error(count, pos, self._format_log(), record['meta'])
                    else:
                        seek = yield container

                    record = None
                    if seek is not None:
                        yield
                        count = seek - 1
                        self.__already_seeked = False
                        continue

                if seekable:
                    pos = file.tell()  # MOL block line starting position
                count += 1
                ir = 3
                failed = is_reaction = False
                mkey = None
                meta = defaultdict(list)
            elif record:
                if line.startswith('$DTYPE'):
                    mkey = line[7:].strip()
                    if not mkey:
                        self._info(f'invalid metadata entry: {line}')
                elif mkey:
                    data = line.lstrip("$DATUM").strip()
                    if data:
                        meta[mkey].append(data)
            elif ir:
                if ir == 3:  # parse mol or rxn title
                    title = line.strip()
                ir -= 1
            else:
                try:
                    if is_reaction:
                        if line.startswith('M  V30 COUNTS'):
                            parser = ERXNRead(line, self._ignore, self._log_buffer)
                        else:
                            parser = RXNRead(line, self._ignore, self._log_buffer)
                    else:
                        if 'V2000' in line:
                            parser = MOLRead(line, self._log_buffer)
                        elif 'V3000' in line:
                            parser = EMOLRead(self._ignore, self._log_buffer)
                        else:
                            raise ValueError('invalid MOL entry')
                except ValueError:
                    failed = True
                    self._info(f'line:\n{line}\nconsist errors:\n{format_exc()}')
                    seek = yield parse_error(count, pos, self._format_log(), {})
                    if seek is not None:
                        yield
                        count = seek - 1
                        self.__already_seeked = False
        if record:
            record['meta'].update(self._prepare_meta(meta))
            if title:
                record['title'] = title
            try:
                if is_reaction:
                    container = self._convert_reaction(record)
                else:
                    container = self._convert_molecule(record)
            except ValueError:
                self._info(f'record consist errors:\n{format_exc()}')
                yield parse_error(count, pos, self._format_log(), record['meta'])
            else:
                yield container

    __already_seeked = False


class _RDFWrite:
    def __init__(self, file, *, append: bool = False, mapping: bool = True):
        """
        :param append: append to existing file (True) or rewrite it (False). For buffered writer object append = False
            will write RDF header and append = True will omit the header.
        :param mapping: write atom mapping.
        """
        super().__init__(file, append=append, mapping=mapping)
        if not append or not (self._is_buffer or self._file.tell() != 0):
            self.write = self.__write

    def __write(self, data):
        """
        write single molecule or reaction into file
        """
        del self.write
        self._file.write(strftime('$RDFILE 1\n$DATM    %m/%d/%y %H:%M\n'))
        self.write(data)


class RDFWrite(_RDFWrite, MDLWrite):
    """
    MDL RDF files writer. works similar to opened for writing file object. support `with` context manager.
    on initialization accept opened for writing in text mode file, string path to file,
    pathlib.Path object or another buffered writer object
    """
    def write(self, data: Union[ReactionContainer, MoleculeContainer]):
        if isinstance(data, ReactionContainer):
            ag = f'{len(data.reagents):3d}' if data.reagents else ''
            self._file.write(f'$RFMT\n$RXN\n{data.name}\n\n\n{len(data.reactants):3d}{len(data.products):3d}{ag}\n')
            for m in chain(data.reactants, data.products, data.reagents):
                m = self._convert_molecule(m)
                self._file.write('$MOL\n')
                self._file.write(m)
        else:
            m = self._convert_molecule(data)
            self._file.write('$MFMT\n')
            self._file.write(m)
        for k, v in data.meta.items():
            self._file.write(f'$DTYPE {k}\n$DATUM {v}\n')


class ERDFWrite(_RDFWrite, EMDLWrite):
    """
    MDL V3000 RDF files writer. works similar to opened for writing file object. support `with` context manager.
    on initialization accept opened for writing in text mode file, string path to file,
    pathlib.Path object or another buffered writer object
    """
    def write(self, data: Union[ReactionContainer, MoleculeContainer]):
        if isinstance(data, ReactionContainer):
            ag = f' {len(data.reagents)}' if data.reagents else ''
            self._file.write(f'$RFMT\n$RXN V3000\n{data.name}\n\n\n'
                             f'M  V30 COUNTS {len(data.reactants)} {len(data.products)}{ag}\nM  V30 BEGIN REACTANT\n')
            for m in data.reactants:
                self._file.write(self._convert_molecule(m))
            self._file.write('M  V30 END REACTANT\nM  V30 BEGIN PRODUCT\n')
            for m in data.products:
                self._file.write(self._convert_molecule(m))
            self._file.write('M  V30 END PRODUCT\n')
            if data.reagents:
                self._file.write('M  V30 BEGIN AGENT\n')
                for m in data.reagents:
                    self._file.write(self._convert_molecule(m))
                self._file.write('M  V30 END AGENT\n')
            self._file.write('M  END\n')
        else:
            m = self._convert_molecule(data)
            self._file.write(f'$MFMT\n{data.name}\n\n\n  0  0  0     0  0            999 V3000\n')
            self._file.write(m)
            self._file.write('M  END\n')
        for k, v in data.meta.items():
            self._file.write(f'$DTYPE {k}\n$DATUM {v}\n')


__all__ = ['RDFRead', 'RDFWrite', 'ERDFWrite']
