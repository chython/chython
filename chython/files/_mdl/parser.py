# -*- coding: utf-8 -*-
#
#  Copyright 2014-2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from collections import namedtuple
from itertools import count
from ...containers import MoleculeContainer, ReactionContainer
from ...containers.bonds import Bond
from ...exceptions import AtomNotFound, MappingError
from ...periodictable import Element


parse_error = namedtuple('ParseError', ('number', 'position', 'log', 'meta'))


class Parser:
    """
    Override classes below then inheritance used.
    """
    MoleculeContainer = MoleculeContainer
    ReactionContainer = ReactionContainer

    def __init__(self, remap=False, ignore=False, store_log=False):
        self.__remap = remap
        self._ignore = ignore
        self._store_log = store_log
        self._log_buffer = []

    def _info(self, msg):
        self._log_buffer.append(msg)

    def _flush_log(self):
        self._log_buffer.clear()

    def _format_log(self):
        log = '\n'.join(self._log_buffer)
        self._log_buffer.clear()
        return log

    def _convert_reaction(self, data):
        if not (data['reactants'] or data['products'] or data['reagents']):
            raise ValueError('empty reaction')
        maps = {'reactants': [], 'products': [], 'reagents': []}
        for i, tmp in maps.items():
            for molecule in data[i]:
                used = set()
                for atom in molecule['atoms']:
                    m = atom['mapping']
                    if m:
                        if m in used:
                            if not self._ignore:
                                raise MappingError('mapping in molecules should be unique')
                            self._info(f'non-unique mapping in molecule: {m}')
                        else:
                            used.add(m)
                    tmp.append(m)

        length = count(max(max(maps['products'], default=0), max(maps['reactants'], default=0),
                           max(maps['reagents'], default=0)) + 1)

        # map unmapped atoms.
        for i, tmp in maps.items():
            used = set()
            maps[i] = remap = []
            for m in tmp:
                if not m:
                    remap.append(next(length))
                elif m in used:
                    if not self._ignore:
                        raise MappingError('mapping in reagents or products or reactants should be unique')
                    # force remap non unique atoms in molecules.
                    remap.append(next(length))
                    self._info(f'mapping in {i} changed from {m} to {remap[-1]}')
                else:
                    remap.append(m)
                    used.add(m)

        if maps['reagents']:
            tmp = (set(maps['reactants']) | set(maps['products'])) & set(maps['reagents'])
            if tmp:
                e = f'reagents has map intersection with reactants or products: {tmp}'
                if not self._ignore:
                    raise MappingError(e)
                self._info(e)
                maps['reagents'] = [x if x not in tmp else next(length) for x in maps['reagents']]

        # find breaks in map. e.g. 1,2,5,6. 3,4 - skipped
        if self.__remap:
            lose = sorted(set(range(1, next(length))) - set(maps['reactants']) - set(maps['products']) -
                          set(maps['reagents']), reverse=True)
            if lose:
                for i, tmp in maps.items():
                    if not tmp:
                        continue
                    for j in lose:
                        maps[i] = tmp = [x if x < j else x - 1 for x in tmp]

        rc = {'reactants': [], 'products': [], 'reagents': []}
        for i, tmp in maps.items():
            shift = 0
            for j in data[i]:
                atom_len = len(j['atoms'])
                remapped = {x: y for x, y in enumerate(tmp[shift: atom_len + shift])}
                shift += atom_len
                g = self._create_molecule(j, remapped)
                g.meta.update(j['meta'])
                rc[i].append(g)
        if self._store_log:
            log = self._format_log()
            if log:
                data['meta']['ParserLog'] = log
        return ReactionContainer(meta=data['meta'], name=data.get('title'), **rc)

    def _convert_molecule(self, data):
        if self.__remap:
            remapped = {n: k for n, k in enumerate(range(1, len(data['atoms']) + 1))}
        else:
            length = count(max(x['mapping'] for x in data['atoms']) + 1)
            remapped, used = {}, set()
            for n, atom in enumerate(data['atoms']):
                m = atom['mapping']
                if not m:
                    remapped[n] = next(length)
                elif m in used:
                    if not self._ignore:
                        raise MappingError('mapping in molecules should be unique')
                    remapped[n] = next(length)
                    self._info(f'mapping in molecule changed from {m} to {remapped[n]}')
                else:
                    remapped[n] = m
                    used.add(m)
        mol = self._create_molecule(data, remapped)
        if self._store_log:
            log = self._format_log()
            if log:
                mol.meta['ParserLog'] = log
        return mol

    def _create_molecule(self, data, mapping, *, _skip_calc_implicit=False):
        g = object.__new__(self.MoleculeContainer)
        pm = {}
        atoms = {}
        plane = {}
        charges = {}
        radicals = {}
        bonds = {}
        for n, atom in enumerate(data['atoms']):
            n = mapping[n]
            atoms[n] = Element.from_symbol(atom['element'])(atom['isotope'])
            bonds[n] = {}
            if (charge := atom['charge']) > 4 or charge < -4:
                raise ValueError('formal charge should be in range [-4, 4]')
            charges[n] = charge
            radicals[n] = atom['is_radical']
            plane[n] = (atom['x'], atom['y'])
            pm[n] = atom['mapping']
        for n, m, b in data['bonds']:
            n, m = mapping[n], mapping[m]
            if n == m:
                raise ValueError('atom loops impossible')
            if n not in bonds or m not in bonds:
                raise AtomNotFound('atoms not found')
            if n in bonds[m]:
                raise ValueError('atoms already bonded')
            bonds[n][m] = bonds[m][n] = Bond(b)
        if any(a['z'] for a in data['atoms']):
            conformers = [{mapping[n]: (a['x'], a['y'], a['z']) for n, a in enumerate(data['atoms'])}]
        else:
            conformers = []
        g.__setstate__(
                {'atoms': atoms, 'bonds': bonds, 'meta': data['meta'], 'plane': plane, 'parsed_mapping': pm,
                 'charges': charges, 'radicals': radicals, 'name': data.get('title', ''), 'conformers': conformers,
                 'atoms_stereo': {}, 'allenes_stereo': {}, 'cis_trans_stereo': {}, 'hydrogens': {}})
        if not _skip_calc_implicit:
            for n in atoms:
                g._calc_implicit(n)
        return g


__all__ = ['Parser', 'parse_error']
