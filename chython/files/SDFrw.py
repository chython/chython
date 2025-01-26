# -*- coding: utf-8 -*-
#
#  Copyright 2014-2024 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from collections import defaultdict
from io import BytesIO
from pickle import dump
from re import compile, match
from subprocess import check_output
from sys import platform
from typing import Optional, List
from .mdl import MDLRead, MOLWrite, EMOLWrite, parse_mol_v2000, parse_mol_v3000, postprocess_molecule
from ._convert import create_molecule
from ._mapping import postprocess_parsed_molecule
from ..containers import MoleculeContainer
from ..exceptions import BufferOverflow, InvalidMolBlock


meta_pattern = compile(r'^>([^<]+)<([^>]+)>([^><]*)$')


class SDFRead(MDLRead):
    """
    MDL SDF files reader. works similar to opened file object. support `with` context manager.
    on initialization accept opened in text mode file, string path to file,
    pathlib.Path object or another buffered reader object
    """
    escape_map = {'&gt;': '>', '&lt;': '<'}
    molecule_cls = MoleculeContainer

    def __init__(self, file, *, buffer_size=10000, indexable: bool = False, ignore: bool = True, remap: bool = False,
                 calc_cis_trans: bool = False, ignore_stereo: bool = False, ignore_bad_isotopes: bool = False):
        """
        :param buffer_size: readahead size. increase if you have big molecules or metadata records.
        :param indexable: if True: supported methods seek, tell, object size and subscription, it only works when
            dealing with a real file (the path to the file is specified) because the external grep utility is used,
            supporting in unix-like OS the object behaves like a normal open file.

            if False: works like generator converting a record into MoleculeContainer and returning each object in
            order, records with errors are skipped
        :param ignore: Skip some checks of data or try to fix some errors.
        :param remap: Remap atom numbers started from one.
        :param calc_cis_trans: Calculate cis/trans marks from 2d coordinates.
        :param ignore_stereo: Ignore stereo data.
        :param ignore_bad_isotopes: reset invalid isotope mark to non-isotopic.
        """
        super().__init__(file, indexable=indexable, ignore=ignore, remap=remap, ignore_bad_isotopes=ignore_bad_isotopes,
                         ignore_stereo=ignore_stereo, calc_cis_trans=calc_cis_trans, buffer_size=buffer_size)
        self.__m_end = None

    def read_structure(self, *, current=True) -> MoleculeContainer:
        data = self._read_mol(current=current)
        if data[4].startswith('M  V30 BEGIN CTAB'):
            tmp = parse_mol_v3000(data)
        else:
            tmp = parse_mol_v2000(data)

        postprocess_parsed_molecule(tmp, remap=self._remap, ignore=self._ignore)
        mol = create_molecule(tmp, ignore_bad_isotopes=self._ignore_bad_isotopes, _cls=self.molecule_cls)
        if not self._ignore_stereo:
            postprocess_molecule(mol, tmp, calc_cis_trans=self._calc_cis_trans)
        meta = self.read_metadata()
        if meta:
            mol.meta.update(meta)
        return mol

    def read_metadata(self, *, current=True):
        mkey = None
        meta = defaultdict(list)
        for line in self._read_metadata(current=current):
            if matched := match(meta_pattern, line):
                mkey = ' '.join(y for x in matched.groups() if (y := x.strip()))
                for e, s in self.escape_map.items():
                    mkey = mkey.replace(e, s)
            elif mkey:
                if line := line.strip():
                    meta[mkey].append(line)
            else:
                meta['chython_unparsed_metadata'].append(line.strip())
        return {k: '\n'.join(v) for k, v in meta.items()}

    def read_mol(self, *, current: bool = True) -> str:
        """
        Read MOL block without metadata
        """
        return ''.join(self._read_mol(current=current))

    def seek(self, offset):
        super().seek(offset)
        self.__m_end = None

    def reset_index(self):
        if platform != 'win32' and not self._is_buffer:
            shifts = [0]
            for x in BytesIO(check_output(['grep', '-bE', r'\$\$\$\$', self._file.name])):
                pos, line = x.split(b':', 1)
                shifts.append(int(pos) + len(line))
            shifts.pop(-1)
            with open(self._cache_path, 'wb') as f:
                dump(shifts, f)
            self._shifts = shifts
        else:
            raise NotImplementedError('Indexable supported in unix-like o.s. and for files stored on disk')

    def _read_block(self, *, current=True) -> List[str]:
        if current and self._buffer:
            return self._buffer
        self.__m_end = m_end = None
        self._buffer = buffer = []
        buffer_size = self._buffer_size

        for n, line in enumerate(self._file):
            if line.startswith('$$$$'):
                break
            elif n == buffer_size:
                raise BufferOverflow
            buffer.append(line)
            if not m_end and line.startswith('M  END'):
                self.__m_end = m_end = len(buffer)
        if buffer:
            self._tell += 1
        else:
            raise EOFError
        return buffer

    def _read_mol(self, *, current: bool = True) -> List[str]:
        data = self._read_block(current=current)
        if not self.__m_end:
            raise InvalidMolBlock
        return data[:self.__m_end]

    def _read_metadata(self, *, current: bool = True):
        data = self._read_block(current=current)
        if not self.__m_end:
            raise InvalidMolBlock
        return data[self.__m_end:]


class SDFWrite(MOLWrite):
    """
    MDL SDF files writer. works similar to opened for writing file object. support `with` context manager.
    on initialization accept opened for writing in text mode file, string path to file,
    pathlib.Path object or another buffered writer object
    """
    escape_map = {'>': '&gt;', '<': '&lt;'}

    def write(self, data: MoleculeContainer, write3d: Optional[int] = None):
        """
        write single molecule into file

        :param write3d: write conformer coordinates with given index
        """
        self._write_molecule(data, write3d=write3d)

        file = self._file
        for k, v in data.meta.items():
            for e, s in self.escape_map.items():
                k = k.replace(e, s)
            file.write(f'>  <{k}>\n{v}\n\n')
        file.write('$$$$\n')


class ESDFWrite(EMOLWrite):
    """
    MDL V3000 SDF files writer. works similar to opened for writing file object. support `with` context manager.
    on initialization accept opened for writing in text mode file, string path to file,
    pathlib.Path object or another buffered writer object
    """
    escape_map = {'>': '&gt;', '<': '&lt;'}

    def write(self, data: MoleculeContainer, write3d: Optional[int] = None):
        """
        write single molecule into file

        :param write3d: write conformer coordinates with given index
        """
        file = self._file
        file.write(f'{data.name}\n\n\n  0  0  0     0  0            999 V3000\n')
        self._write_molecule(data, write3d)
        file.write('M  END\n')

        for k, v in data.meta.items():
            for e, s in self.escape_map.items():
                k = k.replace(e, s)
            file.write(f'>  <{k}>\n{v}\n\n')
        file.write('$$$$\n')


def mdl_mol(data: str, /, *, ignore=True, calc_cis_trans=False, ignore_stereo=False, remap=False,
            ignore_bad_isotopes=False, _cls=MoleculeContainer) -> MoleculeContainer:
    """
    Parse string with mol file.
    """
    data = data.splitlines()
    if data[4].startswith('M  V30 BEGIN CTAB'):
        tmp = parse_mol_v3000(data)
    else:
        tmp = parse_mol_v2000(data)

    postprocess_parsed_molecule(tmp, remap=remap, ignore=ignore)
    mol = create_molecule(tmp, ignore_bad_isotopes=ignore_bad_isotopes, _cls=_cls)
    if not ignore_stereo:
        postprocess_molecule(mol, tmp, calc_cis_trans=calc_cis_trans)
    return mol


__all__ = ['SDFRead', 'SDFWrite', 'ESDFWrite', 'mdl_mol']
