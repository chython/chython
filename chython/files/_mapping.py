# -*- coding: utf-8 -*-
#
#  Copyright 2014-2025 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from itertools import count
from ..exceptions import MappingError


def postprocess_parsed_molecule(data, *, remap=False, ignore=True):
    if remap:
        remapped = list(range(1, len(data['atoms']) + 1))
    else:
        length = count(max(x.get('parsed_mapping') or 0 for x in data['atoms']) + 1)
        remapped, used = [], set()
        for n, atom in enumerate(data['atoms']):
            m = atom.get('parsed_mapping')
            if not m:
                remapped.append(next(length))
            elif m in used:
                if not ignore:
                    raise MappingError('mapping in molecules should be unique')
                remapped.append(next(length))
                if data.get('log') is None:
                    data['log'] = []
                data['log'].append(f'mapping in molecule changed from {m} to {remapped[n]}')
            else:
                remapped.append(m)
                used.add(m)
    data['mapping'] = remapped


def postprocess_parsed_reaction(data, *, remap=False, ignore=True):
    maps = {'reactants': [], 'products': [], 'reagents': []}
    for i, tmp in maps.items():
        for molecule in data[i]:
            used = set()
            for atom in molecule['atoms']:
                m = atom.get('parsed_mapping')
                if m:
                    if m in used:
                        if not ignore:
                            raise MappingError('mapping in molecules should be unique')
                        molecule['log'].append(f'non-unique mapping in molecule: {m}')
                    else:
                        used.add(m)
                    tmp.append(m)
                else:
                    tmp.append(0)

    length = count(max(max(maps['products'], default=0), max(maps['reactants'], default=0),
                       max(maps['reagents'], default=0)) + 1)

    # map unmapped atoms.
    for i, tmp in maps.items():
        used = set()
        maps[i] = _remap = []
        for m in tmp:
            if not m:
                _remap.append(next(length))
            elif m in used:
                if not ignore:
                    raise MappingError('mapping in reagents or products or reactants should be unique')
                # force remap non unique atoms in molecules.
                _remap.append(next(length))
                if data.get('log') is None:
                    data['log'] = []
                data['log'].append(f'mapping in {i} changed from {m} to {_remap[-1]}')
            else:
                _remap.append(m)
                used.add(m)

    if maps['reagents']:
        tmp = (set(maps['reactants']) | set(maps['products'])) & set(maps['reagents'])
        if tmp:
            e = f'reagents has map intersection with reactants or products: {tmp}'
            if not ignore:
                raise MappingError(e)
            if data.get('log') is None:
                data['log'] = []
            data['log'].append(e)
            maps['reagents'] = [x if x not in tmp else next(length) for x in maps['reagents']]

    # find breaks in map. e.g. 1,2,5,6. 3,4 - skipped
    if remap:
        lose = sorted(
            set(range(1, next(length))) - set(maps['reactants']) - set(maps['products']) - set(maps['reagents']),
            reverse=True)
        if lose:
            for i, tmp in maps.items():
                if not tmp:
                    continue
                for j in lose:
                    maps[i] = tmp = [x if x < j else x - 1 for x in tmp]

    for i, tmp in maps.items():
        shift = 0
        for j in data[i]:
            atom_len = len(j['atoms'])
            j['mapping'] = tmp[shift: atom_len + shift]
            shift += atom_len


__all__ = ['postprocess_parsed_molecule', 'postprocess_parsed_reaction']
