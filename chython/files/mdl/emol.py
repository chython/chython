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
from ...exceptions import EmptyMolecule


def parse_mol_v3000(data, *, _header=True):
    if _header:
        title = data[0].strip() or None
        data = data[4:]
    else:
        title = None

    atom_count, bonds_count, *kvs = data[1][13:].split()
    atom_count = int(atom_count)
    if not atom_count:
        raise EmptyMolecule
    bonds_count = int(bonds_count)

    log = []
    atoms = []
    bonds = []
    stereo = []
    meta = {}
    atom_map = {}
    star_points = []

    for kv in kvs:
        if '=' in kv:
            k, v = kv.split('=', 1)
            if k and v:
                meta[k] = v

    # concatenate line breaks
    tmp = []
    keep = None
    for line in data[3:]:
        if line.endswith('-\n'):
            line = line[7:-2]  # skip `M  V30 ` and `-\n`
            if keep:
                keep += line
            else:
                keep = line.lstrip()
        else:
            line = line[7:]  # skip `M  V30 `
            if keep:
                tmp.append(keep + line.rstrip())
                keep = None
            else:
                tmp.append(line.strip())
    data = tmp

    for line in data[:atom_count]:
        n, a, x, y, z, m, *kvs = split(line)
        if a.startswith(('[', 'NOT')):
            raise ValueError('list of atoms not supported')
        elif a == '*':
            star_points.append(n)
            continue
        elif a == 'R#':
            raise ValueError('R-groups not supported')

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

        atom_map[n] = len(atoms)
        atoms.append({'element': a, 'isotope': i, 'charge': c, 'is_radical': r,
                      'x': float(x), 'y': float(y), 'z': float(z), 'parsed_mapping': int(m)})

    for line in data[2 + atom_count: 2 + atom_count + bonds_count]:
        _, t, a1, a2, *kvs = split(line)
        if a1 in star_points:
            if a2 in star_points:
                log.append('invalid bond ignored: star-point to star-point')
                continue
            try:
                star = atom_map[a2]
            except KeyError:
                raise ValueError('invalid atoms number')
            endpoints = None
        elif a2 in star_points:
            try:
                star = atom_map[a1]
            except KeyError:
                raise ValueError('invalid atoms number')
            endpoints = None
        else:
            star = None
            try:
                t = int(t)
                if t in (9, 10):  # added ad-hoc for bond type 9
                    t = 8
                    log.append('coordinate bond replaced to special')
                bonds.append((atom_map[a1], atom_map[a2], t))
            except KeyError:
                raise ValueError('invalid atoms numbers')

        for kv in kvs:
            k, v = kv.split('=')
            if k == 'CFG':
                if v == '1':
                    stereo.append((atom_map[a1], atom_map[a2], 1))
                elif v == '3':
                    stereo.append((atom_map[a1], atom_map[a2], -1))
                else:
                    log.append('invalid or unsupported stereo')
            elif k == 'ENDPTS':
                endpoints = v[1:-1].split()
                if len(endpoints) != int(endpoints[0]) + 1:
                    raise ValueError('invalid ENDPTS block')
        if star is not None:
            if endpoints:  # noqa
                for m in endpoints[1:]:  # noqa
                    try:
                        bonds.append((star, atom_map[m], 8))
                    except KeyError:
                        raise ValueError('invalid atoms numbers in ENDPTS block')
            else:
                log.append('Bond ignored. Star atom not allowed as endpoint')

    drop = True
    for line in data[3 + atom_count + bonds_count:]:
        if line.startswith('END CTAB'):
            break
        elif drop:
            if line.startswith('BEGIN SGROUP'):
                drop = False
            continue
        elif line.startswith('END SGROUP'):
            break

        _, _type, i, *kvs = split(line)
        if _type.startswith('DAT'):
            a = f = d = None
            for kv in kvs:
                k, v = kv.split('=', 1)
                if k == 'ATOMS':
                    a = tuple(atom_map[x] for x in v[1:-1].split()[1:] if x not in star_points)
                elif k == 'FIELDNAME':
                    f = v.strip('"')
                if k == 'FIELDDATA':
                    d = v.strip('"')
            if a and f and d:
                if f == 'MRV_IMPLICIT_H':
                    atoms[a[0]]['implicit_hydrogens'] = int(d[6:])
                else:
                    log.append(f'ignored SGROUP DAT {i}: {a}\t{f}\t{d}')
        elif _type.startswith('SRU'):
            raise ValueError('Polymers not supported')

    return {'title': title, 'atoms': atoms, 'bonds': bonds, 'stereo': stereo, 'meta': meta, 'log': log}


def split(line):  # todo optimize
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


__all__ = ['parse_mol_v3000']
