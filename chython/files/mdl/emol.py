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
from re import compile
from ...exceptions import EmptyMolecule


_stereo_collection = compile(r'MDLV30/STE(ABS|RAC|REL)(\d*)\s+ATOMS=\((.+)\)')
_v3000_token_re = compile(r'(?:[^\s"(]|"[^"]*"|\([^)]*\))+')


def parse_mol_v3000(data, *, _header=True):
    if _header:
        title = data[0].strip() or None
        data = data[4:]
    else:
        title = None

    try:
        counts_tokens = data[1][13:].split()
        atom_count = int(counts_tokens[0])
    except (IndexError, ValueError):
        raise ValueError(f'V3000: cannot parse counts line: {data[1]!r:.80}')
    if not atom_count:
        raise EmptyMolecule
    try:
        bonds_count = int(counts_tokens[1])
    except (IndexError, ValueError):
        raise ValueError(f'V3000: cannot parse bond count from counts line: {data[1]!r:.80}')

    log = []
    atoms = []
    bonds = []
    stereo = []
    meta = {}
    atom_map = {}
    star_points = set()

    for kv in counts_tokens[2:]:
        if '=' in kv:
            k, v = kv.split('=', 1)
            if k and v:
                meta[k] = v

    # concatenate line continuations using list+join
    tmp = []
    parts = []
    for line in data[3:]:
        if line.endswith('-\n'):
            parts.append(line[7:-2])  # skip `M  V30 ` prefix and `-\n` suffix
        else:
            content = line[7:].rstrip()
            if parts:
                parts.append(content)
                tmp.append(''.join(parts))
                parts = []
            else:
                tmp.append(content.lstrip())
    data = tmp
    if len(data) < atom_count:
        raise ValueError(f'V3000: expected {atom_count} atom lines but only {len(data)} lines available')

    # parse atom block
    _split = _v3000_token_re.findall
    for i, line in enumerate(data[:atom_count], 1):
        try:
            tokens = _split(line)
            n, a, x, y, z, m = tokens[0], tokens[1], tokens[2], tokens[3], tokens[4], tokens[5]
            kvs = tokens[6:]
        except (IndexError, ValueError):
            raise ValueError(f'V3000 atom line {i}: cannot parse: {line!r:.80}')

        if a.startswith(('[', 'NOT')):
            raise ValueError(f'V3000 atom line {i}: list of atoms not supported')
        elif a == '*':
            star_points.add(n)
            continue
        elif a == 'R#':
            raise ValueError(f'V3000 atom line {i}: R-groups not supported')

        isotope = None
        charge = 0
        is_radical = False
        for kv in kvs:
            k, v = kv.split('=', 1)
            if k == 'CHG':
                charge = int(v)
            elif k == 'MASS':
                isotope = int(v)
            elif k == 'RAD':
                is_radical = True
        if a == 'D':
            if isotope:
                raise ValueError(f'V3000 atom line {i}: isotope on deuterium atom')
            a = 'H'
            isotope = 2

        try:
            fx, fy, fz = float(x), float(y), float(z)
        except ValueError:
            raise ValueError(f'V3000 atom line {i}: cannot parse coordinates from {x!r}, {y!r}, {z!r}')

        atom_map[n] = len(atoms)
        atoms.append({'element': a, 'isotope': isotope, 'charge': charge, 'is_radical': is_radical,
                      'x': fx, 'y': fy, 'z': fz, 'parsed_mapping': int(m)})

    # parse bond block
    bond_start = 2 + atom_count
    bond_end = bond_start + bonds_count
    if bond_end > len(data):
        raise ValueError(f'V3000: expected {bonds_count} bond lines but only {len(data) - bond_start} available')
    for i, line in enumerate(data[bond_start:bond_end], 1):
        try:
            tokens = _split(line)
            _, t, a1, a2 = tokens[0], tokens[1], tokens[2], tokens[3]
            kvs = tokens[4:]
        except (IndexError, ValueError):
            raise ValueError(f'V3000 bond line {i}: cannot parse: {line!r:.80}')

        if a1 in star_points:
            if a2 in star_points:
                log.append('invalid bond ignored: star-point to star-point')
                continue
            try:
                star = atom_map[a2]
            except KeyError:
                raise ValueError(f'V3000 bond line {i}: invalid atom number {a2}')
            endpoints = None
        elif a2 in star_points:
            try:
                star = atom_map[a1]
            except KeyError:
                raise ValueError(f'V3000 bond line {i}: invalid atom number {a1}')
            endpoints = None
        else:
            star = None
            try:
                bt = int(t)
                if bt in (9, 10):
                    bt = 8
                    log.append('coordinate bond replaced to special')
                bonds.append((atom_map[a1], atom_map[a2], bt))
            except KeyError:
                raise ValueError(f'V3000 bond line {i}: invalid atom numbers {a1}, {a2}')
            except ValueError:
                raise ValueError(f'V3000 bond line {i}: invalid bond type {t!r}')

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
                    raise ValueError(f'V3000 bond line {i}: invalid ENDPTS block')
        if star is not None:
            if endpoints:  # noqa
                for ep in endpoints[1:]:  # noqa
                    try:
                        bonds.append((star, atom_map[ep], 8))
                    except KeyError:
                        raise ValueError(f'V3000 bond line {i}: invalid atom number {ep} in ENDPTS block')
            else:
                log.append('Bond ignored. Star atom not allowed as endpoint')

    # parse remaining blocks (SGROUP, COLLECTION, etc.)
    in_sgroup = False
    in_collection = False
    remaining_start = 3 + atom_count + bonds_count
    for line in data[remaining_start:]:
        if line.startswith('END CTAB'):
            break
        elif line.startswith('BEGIN SGROUP'):
            in_sgroup = True
        elif line.startswith('END SGROUP'):
            in_sgroup = False
        elif line.startswith('BEGIN COLLECTION'):
            in_collection = True
        elif line.startswith('END COLLECTION'):
            in_collection = False
        elif in_sgroup:
            tokens = _split(line)
            try:
                _, _type, idx = tokens[0], tokens[1], tokens[2]
                kvs = tokens[3:]
            except IndexError:
                log.append(f'V3000 SGROUP: cannot parse line: {line!r:.80}')
                continue
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
                        log.append(f'ignored SGROUP DAT {idx}: {a}\t{f}\t{d}')
            elif _type.startswith('SRU'):
                raise ValueError('Polymers not supported')
        elif in_collection:
            m = _stereo_collection.search(line)
            if m:
                stype, gid, atom_list = m.group(1), m.group(2), m.group(3)
                atom_indices = atom_list.split()
                atom_indices = atom_indices[1:]  # skip count
                if stype == 'RAC':  # AND group = positive
                    sg = int(gid) if gid else 1
                    for a in atom_indices:
                        try:
                            atoms[atom_map[a]]['extended_stereo'] = sg
                        except KeyError:
                            log.append(f'invalid atom in STERAC collection: {a}')
                elif stype == 'REL':  # OR group = negative
                    sg = -(int(gid) if gid else 1)
                    for a in atom_indices:
                        try:
                            atoms[atom_map[a]]['extended_stereo'] = sg
                        except KeyError:
                            log.append(f'invalid atom in STEREL collection: {a}')
                # STEABS is the default (no extended_stereo needed)

    return {'title': title, 'atoms': atoms, 'bonds': bonds, 'stereo': stereo, 'meta': meta, 'log': log}


__all__ = ['parse_mol_v3000']
