# -*- coding: utf-8 -*-
#
#  Copyright 2020-2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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

    try:
        atoms_count = int(line[0:3])
    except (ValueError, IndexError):
        raise InvalidV2000(f'V2000 counts line: cannot parse atom count from {line!r:.80}')
    try:
        bonds_count = int(line[3:6])
    except (ValueError, IndexError):
        raise InvalidV2000(f'V2000 counts line: cannot parse bond count from {line!r:.80}')

    if not atoms_count:
        raise EmptyMolecule

    available_lines = len(data) - 4
    if available_lines < atoms_count:
        raise InvalidV2000(f'V2000: expected {atoms_count} atom lines but only {available_lines} lines available')
    if available_lines < atoms_count + bonds_count:
        raise InvalidV2000(f'V2000: expected {bonds_count} bond lines but only {available_lines - atoms_count} '
                           f'lines available after atom block')

    log = []
    title = data[0].strip() or None
    atoms = []
    bonds = []
    stereo = []
    dat = {}

    # parse atom block
    _cm = _charge_map
    for i, line in enumerate(data[4: 4 + atoms_count], 1):
        try:
            charge = _cm[line[36:39]]
        except KeyError:
            raise InvalidCharge(f'V2000 atom line {i}: invalid charge field {line[36:39]!r} in: {line!r:.80}')

        element = line[31:34].strip()
        isotope = line[34:36]

        if element in 'AL':
            raise ValueError(f'V2000 atom line {i}: queries not supported')
        elif element == 'D':
            element = 'H'
            if isotope != ' 0':
                raise ValueError(f'V2000 atom line {i}: isotope on deuterium atom')
            isotope = 2
            delta_isotope = None
        elif isotope != ' 0':
            delta_isotope = int(isotope)
            isotope = None
        else:
            isotope = None
            delta_isotope = None

        try:
            x = float(line[0:10])
            y = float(line[10:20])
            z = float(line[20:30])
        except ValueError:
            raise InvalidV2000(f'V2000 atom line {i}: cannot parse coordinates from {line[0:30]!r}')

        mapping = line[60:63]
        atoms.append({'element': element, 'charge': charge, 'isotope': isotope,
                      'parsed_mapping': int(mapping) if mapping.strip() else 0,
                      'x': x, 'y': y, 'z': z, 'delta_isotope': delta_isotope})

    # parse bond block
    for i, line in enumerate(data[4 + atoms_count: 4 + atoms_count + bonds_count], 1):
        try:
            a1 = int(line[0:3]) - 1
            a2 = int(line[3:6]) - 1
        except ValueError:
            raise InvalidV2000(f'V2000 bond line {i}: cannot parse atom indices from {line[0:6]!r}')
        if a1 < 0 or a1 >= atoms_count or a2 < 0 or a2 >= atoms_count:
            raise InvalidV2000(f'V2000 bond line {i}: atom index out of range '
                               f'(a1={a1 + 1}, a2={a2 + 1}, atoms={atoms_count})')
        s = line[9:12]
        if s == '  1':
            stereo.append((a1, a2, 1))
        elif s == '  6':
            stereo.append((a1, a2, -1))
        elif s != '  0':
            log.append(f'unsupported or invalid stereo: {line}')
        b = int(line[6:9])
        if b == 9:
            b = 8
            log.append(f'coordinate bond replaced with special: {line}')
        bonds.append((a1, a2, b))

    # parse properties block
    for line in data[4 + atoms_count + bonds_count:]:
        if line.startswith('M  END'):
            break
        elif line.startswith('M  ALS'):
            raise ValueError('list of atoms not supported')
        elif line.startswith(('M  ISO', 'M  RAD', 'M  CHG')):
            _type = _ctf_data[line[3]]
            try:
                count = int(line[6:9])
            except ValueError:
                raise InvalidV2000(f'V2000 properties: cannot parse entry count from {line!r:.80}')
            for i in range(count):
                i8 = i * 8
                try:
                    atom = int(line[10 + i8:13 + i8])
                except ValueError:
                    raise InvalidV2000(f'V2000 properties: cannot parse atom number in {line!r:.80}')
                if not atom or atom > atoms_count:
                    raise InvalidV2000(f'V2000 properties: atom number {atom} out of range [1, {atoms_count}] '
                                       f'in: {line!r:.80}')
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
        if 'is_radical' in a:
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
