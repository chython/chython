# -*- coding: utf-8 -*-
#
#  Copyright 2020, 2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from ...exceptions import EmptyMolecule


class EMOLRead:
    def __init__(self, ignore=False, log_buffer=None):
        self.__atoms = []
        self.__bonds = []
        self.__atom_map = {}
        self.__stereo = []
        self.__meta = {}
        self.__hydrogens = {}
        self.__star_points = []
        self.__ignore = ignore
        if log_buffer is None:
            log_buffer = []
        self.__log_buffer = log_buffer

    def getvalue(self):
        if self.__in_mol or self.__in_mol is None:
            raise ValueError('molecule not complete')
        return {'atoms': self.__atoms, 'bonds': self.__bonds, 'stereo': self.__stereo, 'meta': self.__meta,
                'hydrogens': self.__hydrogens}

    def __call__(self, line):
        if self.__in_mol:
            if line.startswith('M  V30 END CTAB'):
                self.__in_mol = False
                return True
            elif self.__atoms_count:
                if line.startswith('M  V30 END'):
                    x = line[11:].strip()
                    cp = self.__parser
                    self.__parser = None

                    if x == 'ATOM':
                        if cp == self.__atom_parser and len(self.__atoms) == self.__atoms_count:
                            return
                    elif x == 'BOND':
                        if cp == self.__bond_parser and len(self.__bonds) == self.__bonds_count:
                            return
                    else:
                        return
                    raise ValueError(f'invalid number of {x} records or invalid CTAB')

                elif self.__parser:
                    collected = self.__record_collector(line)
                    if collected:
                        self.__parser(collected)

                elif line.startswith('M  V30 BEGIN ATOM'):
                    self.__parser = self.__atom_parser
                elif line.startswith('M  V30 BEGIN BOND'):
                    self.__parser = self.__bond_parser
                elif line.startswith('M  V30 BEGIN SGROUP'):
                    self.__parser = self.__sgroup_parser
                elif line.startswith('M  V30 BEGIN COLLECTION'):
                    self.__parser = self.__ignored_block_parser
                else:
                    raise ValueError('invalid CTAB')

            else:  # M  V30 COUNTS line expected
                a, b, *meta = line[13:].split()
                atom_count = int(a)
                if not atom_count:
                    raise EmptyMolecule
                self.__bonds_count = int(b)
                self.__atoms_count = atom_count
                for kv in meta:
                    if '=' in kv:
                        k, v = kv.split('=', 1)
                        if k and v:
                            self.__meta[k] = v

        elif self.__in_mol is not None:
            raise SyntaxError('invalid usage')
        elif not line.startswith('M  V30 BEGIN CTAB'):
            raise ValueError('invalid CTAB')
        else:
            self.__in_mol = True

    def __record_collector(self, line):
        if not line.endswith('-\n'):
            line = line[6:]
            if self.__record:
                line = self.__record + line
                self.__record = None

            line = line.strip()
            collect = []
            tmp = []
            until = None
            for s in line:
                if until:
                    tmp.append(s)
                    if s == until:
                        until = None
                elif s == '(':
                    tmp.append('(')
                    until = ')'
                elif s == '"':
                    tmp.append(s)
                    until = '"'
                elif s == ' ':
                    if tmp:
                        collect.append(''.join(tmp))
                        tmp = []
                else:
                    tmp.append(s)
            if tmp:
                collect.append(''.join(tmp))

            return collect

        line = line[6:-2]
        if not self.__record:
            self.__record = line
        else:
            self.__record += line

    def __bond_parser(self, line):
        _, t, a1, a2, *kvs = line
        if a1 in self.__star_points:
            if a2 in self.__star_points:
                if self.__ignore:
                    self.__log_buffer.append('invalid bond ignored: star-point to star-point')
                    self.__bonds_count -= 1
                    return
                else:
                    raise ValueError('invalid bond: star-point to star-point')
            try:
                star = self.__atom_map[a2]
            except KeyError:
                raise ValueError('invalid atoms number')
            endpoints = None
        elif a2 in self.__star_points:
            try:
                star = self.__atom_map[a1]
            except KeyError:
                raise ValueError('invalid atoms number')
            endpoints = None
        else:
            star = None
            try:
                t = int(t)
                if t in (9, 10):  # added ad-hoc for bond type 9
                    t = 8
                    self.__log_buffer.append('coordinate bond replaced to special')
                self.__bonds.append((self.__atom_map[a1], self.__atom_map[a2], t))
            except KeyError:
                raise ValueError('invalid atoms numbers')

        for kv in kvs:
            k, v = kv.split('=')
            if k == 'CFG':
                if v == '1':
                    self.__stereo.append((self.__atom_map[a1], self.__atom_map[a2], 1))
                elif v == '3':
                    self.__stereo.append((self.__atom_map[a1], self.__atom_map[a2], -1))
                else:
                    self.__log_buffer.append('invalid or unsupported stereo')
            elif k == 'ENDPTS':
                endpoints = v[1:-1].split()
                if len(endpoints) != int(endpoints[0]) + 1:
                    raise ValueError('invalid ENDPTS block')
        if star is not None:
            if endpoints:
                self.__bonds_count += int(endpoints[0]) - 1
                for m in endpoints[1:]:
                    self.__bonds.append((star, self.__atom_map[m], 8))
            elif self.__ignore:
                self.__bonds_count -= 1
                self.__log_buffer.append('Bond ignored. Star atom not allowed as endpoint')
            else:
                raise ValueError('Star atom not allowed')

    def __atom_parser(self, line):
        n, a, x, y, z, m, *kvs = line
        if a.startswith(('[', 'NOT')):
            raise ValueError('list of atoms not supported')
        elif a == '*':
            self.__star_points.append(n)
            self.__atoms_count -= 1
            return

        i = None
        c = 0
        r = False
        for kv in kvs:
            k, v = kv.split('=', 1)
            if k == 'CHG':
                c = int(v)
            elif k == 'MASS':
                i = int(v)
            elif k == 'RAD':
                r = True
        if a == 'D':
            if i:
                raise ValueError('isotope on deuterium atom')
            a = 'H'
            i = 2

        self.__atom_map[n] = len(self.__atoms)
        self.__atoms.append({'element': a, 'isotope': i, 'charge': c, 'is_radical': r,
                             'x': float(x), 'y': float(y), 'z': float(z), 'mapping': int(m)})

    def __sgroup_parser(self, line):
        _, _type, i, *kvs = line
        if _type.startswith('DAT'):
            a = f = d = None
            for kv in kvs:
                k, v = kv.split('=', 1)
                if k == 'ATOMS':
                    a = tuple(self.__atom_map[x] for x in v[1:-1].split()[1:])
                elif k == 'FIELDNAME':
                    f = v.strip('"')
                if k == 'FIELDDATA':
                    d = v.strip('"')
            if a and f and d:
                if f == 'MRV_IMPLICIT_H':
                    self.__hydrogens[a[0]] = int(d[6:])
                else:
                    self.__meta[f'SGROUP DAT {i}'] = (a, f, d)
        elif _type.startswith('SRU'):
            raise ValueError('Polymers not supported')

    def __ignored_block_parser(self, line):
        return

    __record = __atoms_count = __in_mol = __parser = None


__all__ = ['EMOLRead']
