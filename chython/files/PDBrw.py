# -*- coding: utf-8 -*-
#
#  Copyright 2020-2022 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from fileinput import FileInput
from io import StringIO, TextIOWrapper
from itertools import repeat
from pathlib import Path
from traceback import format_exc
from typing import Optional, Iterable
from .XYZrw import XYZ
from ..exceptions import ParseError


one_symbol_names = {'ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN',
                    'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR',  # Amino acids
                    'ASH',  # H-asparagine
                    'HID', 'HIE', 'HIP',  # H-histidines
                    'DA', 'DC', 'DG', 'DT', 'DI',  # Deoxyribonucleotides
                    'A', 'C', 'G', 'U', 'I',  # Ribonucleotides
                    'IOD', 'K'}  # Elements
two_symbol_names = {'CL', 'BR', 'NA', 'CA', 'MG', 'CO', 'MN', 'FE', 'CU', 'ZN'}


class PDBRead(XYZ):
    """PDB files reader. Works similar to opened file object. Support `with` context manager.
    On initialization accept opened in text mode file, string path to file,
    pathlib.Path object or another buffered reader object.

    Supported multiple structures in same file separated by ENDMDL. Supported only ATOM and HETATM parsing.
    END or ENDMDL required in the end.
    """
    def __init__(self, file, ignore=False, element_name_priority=False, parse_as_single=False, atom_name_map=None,
                 charge_map: Optional[Iterable[int]] = None, radical_map: Optional[Iterable[int]] = None, **kwargs):
        """
        :param ignore: Skip some checks of data or try to fix some errors.
        :param store_log: Store parser log if exists messages to `.meta` by key `ParserLog`.
        :param element_name_priority: For ligands use element symbol column value and ignore atom name column.
        :param parse_as_single: Usable if all models in file is the same structure.
            2d graph will be restored only from first model.
        :param atom_name_map: dictionary with atom names replacements. e.g.: {'Ow': 'O'}. Keys should be capitalized.
        :param charge_map: iterable with total charges of each model in file.
        :param radical_map: iterable with total radicals count of each model in file.
        """
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
        self.__file = iter(self._file.readline, '')
        self.__ignore = ignore
        self.__element_name_priority = element_name_priority
        self.__parse_as_single = parse_as_single
        self.__parsed_first = None
        self.__atom_name_map = atom_name_map or {}
        self.__charge_map = charge_map
        self.__radical_map = radical_map
        self._data = self.__reader()

    def __reader(self):
        element_name_priority = self.__element_name_priority
        atom_name_map = self.__atom_name_map
        file = self._file
        try:
            seekable = file.seekable()
        except AttributeError:
            seekable = False
        ignore = self.__ignore
        failkey = False
        charges = repeat(0) if self.__charge_map is None else iter(self.__charge_map)
        radicals = repeat(0) if self.__radical_map is None else iter(self.__radical_map)
        atoms = []

        pos = 0 if seekable else None
        count = 0
        for n, line in enumerate(self.__file):
            if failkey:
                if line.startswith('END'):
                    try:  # drop unused charges/radicals
                        next(charges)
                        next(radicals)
                    except StopIteration:
                        raise ValueError('charge_map or radical_map invalid')
                    failkey = False
                    if seekable:
                        pos = file.tell()
                    count += 1
            elif line.startswith(('ATOM', 'HETATM')):
                # COLUMNS        DATA  TYPE    FIELD        DEFINITION
                # -------------------------------------------------------------------------------------
                #  1 -  6        Record name   "ATOM  " or "HETATM"
                #  7 - 11        Integer       serial       Atom  serial number.
                # 13 - 16        Atom          name         Atom name.
                # 17             Character     altLoc       Alternate location indicator.
                # 18 - 20        Residue name  resName      Residue name.
                # 22             Character     chainID      Chain identifier.
                # 23 - 26        Integer       resSeq       Residue sequence number.
                # 27             AChar         iCode        Code for insertion of residues.
                # 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
                # 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
                # 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
                # 55 - 60        Real(6.2)     occupancy    Occupancy.
                # 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
                # 77 - 78        LString(2)    element      Element symbol, right-justified.
                # 79 - 80        LString(2)    charge       Charge  on the atom.
                charge = line[78:80].strip()
                if charge:
                    try:
                        charge = int(charge)
                    except ValueError:
                        self._info(f'Line [{n}]: consist errors in charge')
                        failkey = True
                        atoms = []
                        yield ParseError(count, pos, self._format_log(), line)
                else:
                    charge = None
                try:
                    x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                except ValueError:
                    self._info(f'Line [{n}]: consist errors in coordinates definition')
                    failkey = True
                    atoms = []
                    yield ParseError(count, pos, self._format_log(), line)
                    continue

                element = line[76:78].strip()
                residue = line[17:20].strip()
                atom_name = line[12:16].strip(' 0123456789-+')
                if residue in one_symbol_names:  # bio-polymers and I
                    atom_name = atom_name[0]
                elif residue in two_symbol_names:
                    atom_name = atom_name[:2].capitalize()
                elif residue == 'MSE':
                    if atom_name.startswith('SE'):
                        atom_name = 'Se'
                    else:
                        atom_name = atom_name[0]
                elif residue == 'CBR':
                    if atom_name.startswith('BR'):
                        atom_name = 'Br'
                    else:
                        atom_name = atom_name[0]
                # ligands    
                elif element_name_priority:
                    atom_name = element
                else:
                    atom_name = atom_name.capitalize()
                    atom_name = atom_name_map.get(atom_name, atom_name)

                if atom_name != element:
                    self._info(f'Atom name and Element symbol not equal: {line[:-1]}')
                    if not ignore:
                        failkey = True
                        atoms = []
                        yield ParseError(count, pos, self._format_log(), line)
                        continue
                atoms.append((atom_name, charge, x, y, z, residue))
            elif line.startswith('END'):  # EOF or end of complex
                if atoms:  # convert collected atoms
                    try:
                        container = self._convert_molecule(atoms, next(charges), next(radicals))
                    except ValueError:
                        self._info(f'Structure consist errors:\n{format_exc()}')
                        yield ParseError(count, pos, self._format_log(), None)
                    except StopIteration:
                        raise ValueError('charge_map or radical_map invalid')
                    else:
                        yield container
                    atoms = []
                else:
                    try:  # drop unused charges/radicals
                        next(charges)
                        next(radicals)
                    except StopIteration:
                        raise ValueError('charge_map or radical_map invalid')
                    self._info(f'Line [{n}]: END or ENDMDL before ATOM or HETATM')
                    yield ParseError(count, pos, self._format_log(), line)
                count += 1
                if seekable:
                    pos = file.tell()
        if atoms:  # ENDMDL or END not found
            self._info('PDB not finished')
            yield ParseError(count, pos, self._format_log(), None)

    def _convert_molecule(self, matrix, charge=0, radical=0):
        if self.__parsed_first is None:
            mol = super()._convert_molecule([(e, c, x, y, z) for e, c, x, y, z, _ in matrix])
            mol.meta['RESIDUE'] = {n: x[-1] for n, x in zip(mol, matrix)}
            if self.__parse_as_single:
                self.__parsed_first = mol.copy()
            return mol
        else:
            if len(self.__parsed_first) != len(matrix):
                raise ValueError('models not equal')
            c = {}
            for (n, a), (e, _, x, y, z, _) in zip(self.__parsed_first.atoms(), matrix):
                if a.atomic_symbol != e:
                    raise ValueError('models or atom order not equal')
                c[n] = (x, y, z)
            mol = self.__parsed_first.copy()
            mol._conformers[0] = c
            return mol


__all__ = ['PDBRead']
