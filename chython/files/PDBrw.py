# -*- coding: utf-8 -*-
#
#  Copyright 2020-2024 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from itertools import islice
from pathlib import Path
from typing import Optional, Sequence, Iterator, List
from .xyz import xyz
from ..containers import MoleculeContainer
from ..exceptions import BufferOverflow


one_symbol_names = {'ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN',
                    'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR',  # Amino acids
                    'ASH',  # H-asparagine
                    'HID', 'HIE', 'HIP',  # H-histidines
                    'DA', 'DC', 'DG', 'DT', 'DI',  # Deoxyribonucleotides
                    'A', 'C', 'G', 'U', 'I',  # Ribonucleotides
                    'IOD', 'K'}  # Elements
two_symbol_names = {'CL', 'BR', 'NA', 'CA', 'MG', 'CO', 'MN', 'FE', 'CU', 'ZN'}


class PDBRead:
    """PDB files reader. Works similar to opened file object. Support `with` context manager.
    On initialization accept opened in text mode file, string path to file,
    pathlib.Path object or another buffered reader object.

    Supported multiple structures in same file separated by ENDMDL. Supported only ATOM and HETATM parsing.
    END or ENDMDL required in the end.
    """
    molecule_cls = MoleculeContainer

    def __init__(self, file, *, buffer_size=10000, ignore: bool = True, element_name_priority: bool = False,
                 parse_as_single: bool = False, atom_name_map=None, charge_map: Optional[Sequence[int]] = None,
                 radical_map: Optional[Sequence[int]] = None, radius_multiplier: float = 1.25):
        """
        :param ignore: Skip some checks of data or try to fix some errors.
        :param element_name_priority: For ligands use element symbol column value and ignore atom name column.
        :param parse_as_single: Usable if all models in file is the same structure.
            2d graph will be restored only from first model.
        :param atom_name_map: dictionary with atom names replacements. e.g.: {'Ow': 'O'}. Keys should be capitalized.
        :param charge_map: iterable with total charges of each model in file.
        :param radical_map: iterable with total radicals count of each model in file.
        """
        if isinstance(file, str):
            self.__file = open(file)
            self.__is_buffer = False
        elif isinstance(file, Path):
            self.__file = file.open()
            self.__is_buffer = False
        elif isinstance(file, (TextIOWrapper, StringIO, FileInput)):
            self.__file = file
            self.__is_buffer = True
        else:
            raise TypeError('invalid file. TextIOWrapper, StringIO subclasses expected')

        self.__radius_multiplier = radius_multiplier
        self.__ignore = ignore
        self.__element_name_priority = element_name_priority
        self.__parse_as_single = parse_as_single
        self.__parsed_first = None
        self.__atom_name_map = atom_name_map or {}
        self.__charge_map = charge_map
        self.__radical_map = radical_map

        self.__buffer = None
        self.__buffer_size = buffer_size
        self.__tell = 0

    def read(self, amount: Optional[int] = None) -> List[MoleculeContainer]:
        """
        Parse whole file

        :param amount: number of records to read
        """
        if amount:
            return list(islice(iter(self), amount))
        return list(iter(self))

    def read_structure(self, *, current: bool = True) -> MoleculeContainer:
        """
        Read Molecule container.

        :param current: return current structure if already parsed, otherwise read next
        """
        data = self._read_block(current=current)
        element_name_priority = self.__element_name_priority
        atom_name_map = self.__atom_name_map

        atoms = []
        res = []
        charges = []
        log = []
        for line in data:
            if line.startswith(('ATOM', 'HETATM')):
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
                    charge = int(charge)
                else:
                    charge = None
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
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
                    log.append(f'Atom name and Element symbol is not equal: {line[:-1]}')
                    if not self.__ignore:
                        raise ValueError('Atom name and Element symbol is not equal')
                atoms.append((atom_name, x, y, z))
                res.append(residue)
                charges.append(charge)

        if not atoms:
            raise ValueError('invalid PDB')

        if self.__charge_map:
            t_charge = self.__charge_map[self.__tell - 1]
            radical = self.__radical_map[self.__tell - 1]
        else:
            t_charge = radical = 0

        if self.__parsed_first is None:
            mol = xyz(atoms, charge=t_charge, radical=radical, radius_multiplier=self.__radius_multiplier,
                      atom_charge=charges, _cls=self.molecule_cls)

            mol.meta['RESIDUE'] = dict(enumerate(res, 1))
            if log:
                mol.meta['chython_parsing_log'] = log
            if self.__parse_as_single:
                self.__parsed_first = mol.copy()
            return mol
        else:
            if len(self.__parsed_first) != len(atoms):
                raise ValueError('models not equal')
            c = {}
            for (n, a), (e, x, y, z) in zip(self.__parsed_first.atoms(), atoms):
                if a.atomic_symbol != e:
                    raise ValueError('models or atom order not equal')
                c[n] = (x, y, z)
            mol = self.__parsed_first.copy()
            mol._conformers[0] = c
            if log:
                if 'chython_parsing_log' in mol.meta:
                    mol.meta['chython_parsing_log'] = mol.meta['chython_parsing_log'] + log
                else:
                    mol.meta['chython_parsing_log'] = log
            return mol

    def close(self, force: bool = False):
        """
        Close opened file

        :param force: force closing of externally opened file or buffer
        """
        if not self.__is_buffer or force:
            self.__file.close()

    def tell(self):
        """
        Number of records processed from the original file
        """
        return self.__tell

    def __enter__(self):
        return self

    def __exit__(self, _type, value, traceback):
        self.close()

    def __iter__(self) -> Iterator[MoleculeContainer]:
        while True:
            try:
                yield self.read_structure(current=False)
            except ValueError:
                pass
            except EOFError:
                return

    def __next__(self) -> MoleculeContainer:
        return next(iter(self))

    def _read_block(self, *, current=True) -> List[str]:
        if current and self.__buffer:
            return self.__buffer
        self.__buffer = None
        buffer_size = self.__buffer_size
        buffer = []

        for n, line in enumerate(self.__file):
            buffer.append(line)
            if line.startswith('END'):
                break
            elif n == buffer_size:
                raise BufferOverflow
        else:
            raise EOFError

        self.__tell += 1
        self.__buffer = buffer
        return buffer


__all__ = ['PDBRead']
