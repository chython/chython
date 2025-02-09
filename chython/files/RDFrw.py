# -*- coding: utf-8 -*-
#
#  Copyright 2014-2024 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from collections import defaultdict
from io import BytesIO
from itertools import chain
from pickle import dump
from subprocess import check_output
from sys import platform
from time import strftime
from typing import Union, Dict, List
from .mdl import (MDLRead, MOLWrite, EMOLWrite, parse_mol_v2000, parse_mol_v3000, parse_rxn_v2000, parse_rxn_v3000,
                  postprocess_molecule)
from ._convert import create_molecule, create_reaction
from ._mapping import postprocess_parsed_molecule, postprocess_parsed_reaction
from ..containers import ReactionContainer, MoleculeContainer
from ..exceptions import BufferOverflow


class RDFRead(MDLRead):
    """
    MDL RDF files reader. works similar to opened file object. support `with` context manager.
    on initialization accept opened in text mode file, string path to file,
    pathlib.Path object or another buffered reader object
    """
    molecule_cls = MoleculeContainer
    reaction_cls = ReactionContainer

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
        self.__m_start = None

    def read_structure(self, *, current=True) -> Union[ReactionContainer, MoleculeContainer]:
        data = self._read_block(current=current)
        meta = self.read_metadata()
        if data[0].startswith('$RXN'):
            if data[4].startswith('M  V30 COUNTS'):
                tmp = parse_rxn_v3000(data, ignore=self._ignore)
            else:
                tmp = parse_rxn_v2000(data, ignore=self._ignore)

            postprocess_parsed_reaction(tmp, remap=self._remap, ignore=self._ignore)
            rxn = create_reaction(tmp, ignore_bad_isotopes=self._ignore_bad_isotopes, _m_cls=self.molecule_cls,
                                  _r_cls=self.reaction_cls)
            if not self._ignore_stereo:
                for mol, tmp in zip(rxn.molecules(), chain(tmp['reactants'], tmp['reagents'], tmp['products'])):
                    postprocess_molecule(mol, tmp, calc_cis_trans=self._calc_cis_trans)
            if meta:
                rxn.meta.update(meta)
            return rxn
        elif data[4].startswith('M  V30 BEGIN CTAB'):
            tmp = parse_mol_v3000(data)
        else:
            tmp = parse_mol_v2000(data)

        postprocess_parsed_molecule(tmp)
        mol = create_molecule(tmp, ignore_bad_isotopes=self._ignore_bad_isotopes, _cls=self.molecule_cls)
        if not self._ignore_stereo:
            postprocess_molecule(mol, tmp, calc_cis_trans=self._calc_cis_trans)
        if meta:
            mol.meta.update(meta)
        return mol

    def read_metadata(self, *, current=True) -> Dict[str, str]:
        mkey = None
        meta = defaultdict(list)
        for line in self._read_metadata(current=current):
            if line.startswith('$DTYPE'):
                mkey = line[7:].strip()
                if not mkey:
                    meta['chython_unparsed_metadata'].append(line.strip())
            elif mkey:
                data = line.lstrip("$DATUM").strip()
                if data:
                    meta[mkey].append(data)
            else:
                meta['chython_unparsed_metadata'].append(line.strip())
        return {k: '\n'.join(v) for k, v in meta.items()}

    def read_rxn(self, *, current: bool = True) -> str:
        """
        Read rxn block without metadata
        """
        return ''.join(self._read_block(current=current)[:self.__m_start])

    def read_mol(self, n: int, /, *, current: bool = True) -> str:
        """
        Read requested MOL block
        """
        data = self._read_block(current=current)
        if data[0].startswith('$RXN'):
            if data[4].startswith('M  V30 COUNTS'):
                idx_l = [i for i, x in enumerate(data) if x.startswith('M  V30 BEGIN CTAB')]
                idx_r = [i for i, x in enumerate(data, 1) if x.startswith('M  V30 END CTAB')]
                if 0 <= n < len(idx_l):
                    ct = ''.join(data[idx_l[n]: idx_r[n]])
                    return f'\n\n\n  0  0  0     0  0            999 V3000\n{ct}M  END\n'
                else:
                    raise IndexError('molecule number is out of range')
            else:
                idx = [i for i, x in enumerate(data) if x.startswith('$MOL')]
                if 0 <= n < len(idx):
                    idx.append(self.__m_start)
                    return ''.join(data[idx[n] + 1: idx[n + 1]])
                else:
                    raise IndexError('molecule number is out of range')
        elif n:
            raise IndexError('molecule number is out of range')
        return ''.join(data[:self.__m_start])

    def seek(self, offset):
        super().seek(offset)
        self.__m_start = None

    def reset_index(self):
        if platform != 'win32' and not self._is_buffer:
            shifts = []
            for x in BytesIO(check_output(['grep', '-bE', r'^\$[RM]FMT', self._file.name])):
                pos, _ = x.split(b':', 1)
                shifts.append(int(pos))
            shifts[0] = 0  # first record parsing always starts from the beginning
            with open(self._cache_path, 'wb') as f:
                dump(shifts, f)
            self._shifts = shifts
        else:
            raise NotImplementedError('Indexable supported in unix-like o.s. and for files stored on disk')

    def _read_block(self, *, current=True) -> List[str]:
        """
        Read RXN or MOL block with metadata
        """
        if current and self._buffer:
            return self._buffer
        self.__m_start = m_start = None
        self._buffer = buffer = []
        buffer_size = self._buffer_size

        drop = not self._tell  # only first record starts from [RM]FMT search
        for n, line in enumerate(self._file):
            if drop:
                if line.startswith('$RXN'):  # RXN file
                    drop = False
                    buffer.append(line)
                elif line.startswith(('$RFMT', '$MFMT')):  # first occurrence found
                    drop = False
                continue
            elif n == buffer_size:
                raise BufferOverflow
            elif not m_start and line.startswith('$DTYPE'):
                self.__m_start = m_start = len(buffer)
            elif line.startswith(('$RFMT', '$MFMT')):  # next record found
                break
            buffer.append(line)
        if buffer:
            self._tell += 1
        else:
            raise EOFError
        return buffer

    def _read_metadata(self, *, current: bool = True):
        data = self._read_block(current=current)
        if not self.__m_start:
            return []
        return data[self.__m_start:]


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


class RDFWrite(_RDFWrite, MOLWrite):
    """
    MDL RDF files writer. works similar to opened for writing file object. support `with` context manager.
    on initialization accept opened for writing in text mode file, string path to file,
    pathlib.Path object or another buffered writer object
    """
    def write(self, data: Union[ReactionContainer, MoleculeContainer]):
        file = self._file
        if isinstance(data, ReactionContainer):
            file.write(f'$RFMT\n$RXN\n{data.name}\n\n\n{len(data.reactants):3d}{len(data.products):3d}')
            if data.reagents:
                file.write(f'{len(data.reagents):3d}\n')
            else:
                file.write('\n')
            for m in chain(data.reactants, data.products, data.reagents):
                file.write('$MOL\n')
                self._write_molecule(m)
        else:
            file.write('$MFMT\n')
            self._write_molecule(data)
        for k, v in data.meta.items():
            file.write(f'$DTYPE {k}\n$DATUM {v}\n')


class ERDFWrite(_RDFWrite, EMOLWrite):
    """
    MDL V3000 RDF files writer. works similar to opened for writing file object. support `with` context manager.
    on initialization accept opened for writing in text mode file, string path to file,
    pathlib.Path object or another buffered writer object
    """
    def write(self, data: Union[ReactionContainer, MoleculeContainer]):
        file = self._file
        if isinstance(data, ReactionContainer):
            file.write(f'$RFMT\n$RXN V3000\n{data.name}\n\n\nM  V30 COUNTS {len(data.reactants)} {len(data.products)}')
            if data.reagents:
                file.write(f' {len(data.reagents)}\nM  V30 BEGIN REACTANT\n')
            else:
                file.write('\nM  V30 BEGIN REACTANT\n')
            for m in data.reactants:
                self._write_molecule(m)
            file.write('M  V30 END REACTANT\nM  V30 BEGIN PRODUCT\n')
            for m in data.products:
                self._write_molecule(m)
            file.write('M  V30 END PRODUCT\n')
            if data.reagents:
                file.write('M  V30 BEGIN AGENT\n')
                for m in data.reagents:
                    self._write_molecule(m)
                file.write('M  V30 END AGENT\n')
            file.write('M  END\n')
        else:
            file.write(f'$MFMT\n{data.name}\n\n\n  0  0  0     0  0            999 V3000\n')
            self._write_molecule(data)
            file.write('M  END\n')
        for k, v in data.meta.items():
            file.write(f'$DTYPE {k}\n$DATUM {v}\n')


def mdl_rxn(data: str, /, *, ignore=True, calc_cis_trans=False, ignore_stereo=False, remap=False,
            ignore_bad_isotopes=False, _r_cls=ReactionContainer, _m_cls=MoleculeContainer) -> ReactionContainer:
    """
    Parse string with rxn file.
    """
    data = data.splitlines()
    if not data[0].startswith('$RXN'):
        raise ValueError('invalid RXN')
    if data[4].startswith('M  V30 COUNTS'):
        tmp = parse_rxn_v3000(data, ignore=ignore)
    else:
        tmp = parse_rxn_v2000(data, ignore=ignore)

    postprocess_parsed_reaction(tmp, remap=remap, ignore=ignore)
    rxn = create_reaction(tmp, ignore_bad_isotopes=ignore_bad_isotopes, _m_cls=_m_cls, _r_cls=_r_cls)
    if not ignore_stereo:
        for mol, tmp in zip(rxn.molecules(), chain(tmp['reactants'], tmp['reagents'], tmp['products'])):
            postprocess_molecule(mol, tmp, calc_cis_trans=calc_cis_trans)
    return rxn


__all__ = ['RDFRead', 'RDFWrite', 'ERDFWrite', 'mdl_rxn']
