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
from ...exceptions import EmptyMolecule, InvalidCharge, InvalidV2000


_ctf_data = {'R': 'is_radical', 'C': 'charge', 'I': 'isotope'}
_charge_map = {'  0': 0, '  1': 3, '  2': 2, '  3': 1, '  4': 0, '  5': -1, '  6': -2, '  7': -3}


def parse_mol_v2000(data):
    line = data[3]

    atoms_count = int(line[0:3])
    bonds_count = int(line[3:6])
    if not atoms_count:
        raise EmptyMolecule

    log = []
    title = data[0].strip() or None
    atoms = []
    bonds = []
    stereo = []
    dat = {}

    for line in data[4: 4 + atoms_count]:
        try:
            charge = _charge_map[line[36:39]]
        except KeyError:
            raise InvalidCharge
        element = line[31:34].strip()
        isotope = line[34:36]
        delta_isotope = None

        if element in 'AL':
            raise ValueError('queries not supported')
        elif element == 'D':
            element = 'H'
            if isotope != ' 0':
                raise ValueError('isotope on deuterium atom')
            isotope = 2
        elif isotope != ' 0':
            delta_isotope = int(isotope)
            isotope = None
        else:
            isotope = None

        mapping = line[60:63]
        atoms.append({'element': element, 'charge': charge, 'isotope': isotope,
                      'parsed_mapping': int(mapping) if mapping else 0, 'x': float(line[0:10]), 'y': float(line[10:20]),
                      'z': float(line[20:30]), 'delta_isotope': delta_isotope})

    for line in data[4 + atoms_count: 4 + atoms_count + bonds_count]:
        a1, a2 = int(line[0:3]) - 1, int(line[3:6]) - 1
        s = line[9:12]
        if s == '  1':
            stereo.append((a1, a2, 1))
        elif s == '  6':
            stereo.append((a1, a2, -1))
        elif s != '  0':
            log.append(f'unsupported or invalid stereo: {line}')
        b = int(line[6:9])
        if b == 9:  # added ad-hoc for bond type 9
            b = 8
            log.append(f'coordinate bond replaced with special: {line}')
        bonds.append((a1, a2, b))

    for line in data[4 + atoms_count + bonds_count:]:
        if line.startswith('M  END'):
            break
        elif line.startswith('M  ALS'):
            raise ValueError('list of atoms not supported')
        elif line.startswith(('M  ISO', 'M  RAD', 'M  CHG')):
            _type = _ctf_data[line[3]]
            for i in range(int(line[6:9])):
                i8 = i * 8
                atom = int(line[10 + i8:13 + i8])
                if not atom or atom > len(atoms):
                    raise InvalidV2000('invalid atoms number')
                atom = atoms[atom - 1]
                atom[_type] = int(line[14 + i8:17 + i8])

        elif line.startswith('M  STY'):
            for i in range(int(line[6:9])):
                i8 = i * 8
                if (st := line[14 + i8:17 + i8]) == 'DAT':
                    dat[int(line[10 + i8:13 + i8])] = {}
                elif st == 'SUP':
                    dat[int(line[10 + i8:13 + i8])] = {'type': 'MDL_SUP'}
        elif line.startswith('M  SAL'):
            i = int(line[7:10])
            if i in dat:
                dat[i]['atoms'] = tuple(int(line[14 + 4 * i:17 + 4 * i]) - 1 for i in range(int(line[10:13])))
        elif line.startswith('M  SDT'):
            i = int(line[7:10])
            if i in dat:
                dat[i]['type'] = line.split()[-1].lower()
        elif line.startswith('M  SED'):
            i = int(line[7:10])
            if i in dat:
                dat[i]['value'] = line[10:].strip().replace('/', '').lower()
        elif line.startswith('M  SMT'):
            i = int(line[7:10])
            if i in dat:
                dat[i]['value'] = line[10:].strip()
        elif not line.startswith('M  SDD'):
            log.append(f'ignored line: {line}')

    for a in atoms:
        if 'is_radical' in a:  # int to bool
            a['is_radical'] = True
    for x in dat.values():
        try:
            _type = x['type']
            if _type == 'mrv_implicit_h':
                _atoms = x['atoms']
                value = x['value']
                if len(_atoms) != 1 or _atoms[0] == -1 or not value:
                    raise InvalidV2000(f'MRV_IMPLICIT_H spec invalid {x}')
                atoms[_atoms[0]]['implicit_hydrogens'] = int(value[6:])
            else:
                log.append(f'ignored data: {x}')
        except KeyError:
            raise InvalidV2000(f'Invalid SGROUP {x}')

    return {'title': title, 'atoms': atoms, 'bonds': bonds, 'stereo': stereo, 'log': log}


__all__ = ['parse_mol_v2000']
