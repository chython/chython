# -*- coding: utf-8 -*-
#
#  Copyright 2017-2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from collections import defaultdict, namedtuple
from io import StringIO, BytesIO, TextIOWrapper, BufferedIOBase, BufferedReader
from itertools import count
from lxml.etree import iterparse, QName, tostring
from pathlib import Path
from traceback import format_exc
from typing import Union, List, Iterator
from ._mdl import MDLStereo
from ..containers import MoleculeContainer, ReactionContainer
from ..exceptions import EmptyMolecule


parse_error = namedtuple('MRVParseError', ('number', 'json', 'log', 'meta'))
organic_set = {'B', 'C', 'N', 'O', 'P', 'S', 'Se', 'F', 'Cl', 'Br', 'I'}
bond_map = {8: '1" queryType="Any', 4: 'A', 1: '1', 2: '2', 3: '3',
            'Any': 8, 'any': 8, 'A': 4, 'a': 4, '1': 1, '2': 2, '3': 3}


def xml_dict(parent_element, stop_list=None):
    stop_list = set() if stop_list is None else set(stop_list)
    out = {}
    for x, y in parent_element.items():
        y = y.strip()
        if y:
            x = '@%s' % x.strip()
            out[x] = y

    text = []
    if len(parent_element):
        elements_grouped = defaultdict(list)
        for element in parent_element:
            name = QName(element).localname
            if name in stop_list:
                text.append(tostring(element, encoding=str, with_tail=False))
            else:
                elements_grouped[name].append(element)

            if element.tail:
                t = element.tail.strip()
                if t:
                    text.append(t)

        for element_tag, element_group in elements_grouped.items():
            if len(element_group) == 1:
                out[element_tag] = xml_dict(element_group[0], stop_list)
            else:
                out[element_tag] = [xml_dict(x, stop_list) for x in element_group]

    if parent_element.text:
        t = parent_element.text.strip()
        if t:
            text.insert(0, t)
    if text:
        out['$'] = ''.join(text)

    return out


class MRVRead(MDLStereo):
    """
    ChemAxon MRV files reader. works similar to opened file object. support `with` context manager.
    on initialization accept opened in binary mode file, string path to file,
    pathlib.Path object or another binary buffered reader object
    """
    def __init__(self, file, **kwargs):
        """
        :param ignore: Skip some checks of data or try to fix some errors.
        :param remap: Remap atom numbers started from one.
        :param store_log: Store parser log if exists messages to `.meta` by key `ParserLog`.
        :param calc_cis_trans: Calculate cis/trans marks from 2d coordinates.
        :param ignore_stereo: Ignore stereo data.
        """
        if isinstance(file, str):
            self.__file = open(file, 'rb')
            self.__is_buffer = False
        elif isinstance(file, Path):
            self.__file = file.open('rb')
            self.__is_buffer = False
        elif isinstance(file, (BytesIO, BufferedReader, BufferedIOBase)):
            self.__file = file
            self.__is_buffer = True
        else:
            raise TypeError('invalid file. BytesIO, BufferedReader and BufferedIOBase subclasses possible')
        super().__init__(**kwargs)
        self._data = self.__reader()

    def close(self, force=False):
        """
        close opened file

        :param force: force closing of externally opened file or buffer
        """
        if not self.__is_buffer or force:
            self.__file.close()

    def __enter__(self):
        return self

    def __exit__(self, _type, value, traceback):
        self.close()

    def read(self) -> List[Union[ReactionContainer, MoleculeContainer]]:
        """
        parse whole file

        :return: list of parsed molecules or reactions
        """
        return list(iter(self))

    def __iter__(self) -> Iterator[Union[ReactionContainer, MoleculeContainer]]:
        return (x for x in self._data if not isinstance(x, parse_error))

    def __next__(self) -> Union[ReactionContainer, MoleculeContainer]:
        return next(iter(self))

    def __reader(self) -> Iterator[Union[ReactionContainer, MoleculeContainer, parse_error]]:
        for n, (_, element) in enumerate(iterparse(self.__file, tag='{*}MChemicalStruct')):
            parsed = xml_dict(element)
            element.clear()
            if 'molecule' in parsed and isinstance(parsed['molecule'], dict):
                parsed = parsed['molecule']
                if 'propertyList' in parsed and 'property' in parsed['propertyList']:
                    meta = self.__parse_property(parsed['propertyList']['property'])
                else:
                    meta = {}

                try:
                    record = self.__parse_molecule(parsed)
                except (KeyError, ValueError):
                    self._info(f'record consist errors:\n{format_exc()}')
                    yield parse_error(n, parsed, self._format_log(), meta)
                else:
                    record['meta'].update(meta)
                    try:
                        container = self._convert_molecule(record)
                    except ValueError:
                        self._info(f'record consist errors:\n{format_exc()}')
                        yield parse_error(n, parsed, self._format_log(), meta)
                    else:
                        yield container
            elif 'reaction' in parsed and isinstance(parsed['reaction'], dict):
                parsed = parsed['reaction']
                if 'propertyList' in parsed and 'property' in parsed['propertyList']:
                    meta = self.__parse_property(parsed['propertyList']['property'])
                else:
                    meta = {}

                try:
                    record = self.__parse_reaction(parsed)
                except (KeyError, ValueError):
                    self._info(f'record consist errors:\n{format_exc()}')
                    yield parse_error(n, parsed, self._format_log(), meta)
                else:
                    record['meta'] = meta
                    try:
                        container = self._convert_reaction(record)
                    except ValueError:
                        self._info(f'record consist errors:\n{format_exc()}')
                        yield parse_error(n, parsed, self._format_log(), meta)
                    else:
                        yield container
            else:
                self._info('invalid MDocument')
                yield parse_error(n, parsed, self._format_log(), {})

    def __parse_reaction(self, data):
        reaction = {'reactants': [], 'products': [], 'reagents': []}
        title = data.get('@title')
        if title:
            reaction['title'] = title
        for tag, group in (('reactantList', 'reactants'), ('productList', 'products'), ('agentList', 'reagents')):
            if tag in data and 'molecule' in data[tag]:
                molecule = data[tag]['molecule']
                if isinstance(molecule, dict):
                    molecule = (molecule,)
                for m in molecule:
                    try:
                        reaction[group].append(self.__parse_molecule(m))
                    except EmptyMolecule:
                        if not self._ignore:
                            raise
                        self._info('empty molecule ignored')
        return reaction

    def __parse_property(self, data):
        meta = {}
        if isinstance(data, dict):
            key = data['@title']
            val = data['scalar']['$'].strip()
            if key and val:
                meta[key] = val
            else:
                self._info(f'invalid metadata entry: {data}')
        else:
            for x in data:
                key = x['@title']
                val = x['scalar']['$'].strip()
                if key and val:
                    meta[key] = val
                else:
                    self._info(f'invalid metadata entry: {x}')
        return meta

    def __parse_molecule(self, data):
        atoms, bonds, stereo = [], [], []
        hydrogens = {}
        atom_map = {}
        if 'atom' in data['atomArray']:
            da = data['atomArray']['atom']
            if isinstance(da, dict):
                da = (da,)
            for n, atom in enumerate(da):
                atom_map[atom['@id']] = n
                atoms.append({'element': atom['@elementType'],
                              'isotope': int(atom['@isotope']) if '@isotope' in atom else None,
                              'charge': int(atom.get('@formalCharge', 0)),
                              'is_radical': '@radical' in atom,
                              'mapping': int(atom.get('@mrvMap', 0))})
                if '@z3' in atom:
                    atoms[-1].update(x=float(atom['@x3']), y=float(atom['@y3']), z=float(atom['@z3']))
                else:
                    atoms[-1].update(x=float(atom['@x2']) / 2, y=float(atom['@y2']) / 2, z=0.)
                if '@mrvQueryProps' in atom:
                    raise ValueError('queries unsupported')
                if '@hydrogenCount' in atom:
                    hydrogens[n] = int(atom['@hydrogenCount'])
        else:
            atom = data['atomArray']
            for n, (_id, e) in enumerate(zip(atom['@atomID'].split(), atom['@elementType'].split())):
                atom_map[_id] = n
                atoms.append({'element': e, 'charge': 0, 'mapping': 0, 'isotope': None, 'is_radical': False})
            if '@z3' in atom:
                for a, x, y, z in zip(atoms, atom['@x3'].split(), atom['@y3'].split(), atom['@z3'].split()):
                    a['x'] = float(x)
                    a['y'] = float(y)
                    a['z'] = float(z)
            else:
                for a, x, y in zip(atoms, atom['@x2'].split(), atom['@y2'].split()):
                    a['x'] = float(x) / 2
                    a['y'] = float(y) / 2
                    a['z'] = 0.
            if '@isotope' in atom:
                for a, x in zip(atoms, atom['@isotope'].split()):
                    if x != '0':
                        a['isotope'] = int(x)
            if '@formalCharge' in atom:
                for a, x in zip(atoms, atom['@formalCharge'].split()):
                    if x != '0':
                        a['charge'] = int(x)
            if '@mrvMap' in atom:
                for a, x in zip(atoms, atom['@mrvMap'].split()):
                    if x != '0':
                        a['mapping'] = int(x)
            if '@radical' in atom:
                for a, x in zip(atoms, atom['@radical'].split()):
                    if x != '0':
                        a['is_radical'] = True
            if '@mrvQueryProps' in atom:
                raise ValueError('queries unsupported')
        if not atoms:
            raise EmptyMolecule

        if 'bond' in data['bondArray']:
            db = data['bondArray']['bond']
            if isinstance(db, dict):
                db = (db,)
            for bond in db:
                order = bond_map[bond['@queryType' if '@queryType' in bond else '@order']]
                a1, a2 = bond['@atomRefs2'].split()
                if 'bondStereo' in bond:
                    if '$' in bond['bondStereo']:
                        s = bond['bondStereo']['$']
                        if s == 'H':
                            stereo.append((atom_map[a1], atom_map[a2], -1))
                        elif s == 'W':
                            stereo.append((atom_map[a1], atom_map[a2], 1))
                        else:
                            self._info('invalid or unsupported stereo')
                    else:
                        self._info('incorrect bondStereo tag')
                bonds.append((atom_map[a1], atom_map[a2], order))

        mol = {'atoms': atoms, 'bonds': bonds, 'stereo': stereo, 'meta': {}, 'hydrogens': hydrogens}
        if '@title' in data:
            mol['title'] = data['@title']
        return mol


class MRVWrite:
    """
    ChemAxon MRV files writer. works similar to opened for writing file object. support `with` context manager.
    on initialization accept opened for writing in text mode file, string path to file,
    pathlib.Path object or another buffered writer object
    """
    def __init__(self, file):
        if isinstance(file, str):
            self._file = open(file, 'w')
            self._is_buffer = False
        elif isinstance(file, Path):
            self._file = file.open('w')
            self._is_buffer = False
        elif isinstance(file, (TextIOWrapper, StringIO)):
            self._file = file
            self._is_buffer = True
        else:
            raise TypeError('invalid file. '
                            'TextIOWrapper, StringIO, BytesIO, BufferedReader and BufferedIOBase subclasses possible')
        self.__writable = True

    def close(self, force=False):
        """
        write close tag of MRV file and close opened file

        :param force: force closing of externally opened file or buffer
        """
        if not self.__finalized:
            self._file.write('</cml>\n')
            self.__finalized = True
        if self.__writable:
            self.write = self.__write_closed
            self.__writable = False

        if not self._is_buffer or force:
            self._file.close()

    def __enter__(self):
        return self

    def __exit__(self, _type, value, traceback):
        self.close()

    @staticmethod
    def __write_closed(_):
        raise ValueError('I/O operation on closed writer')

    def write(self, data: Union[ReactionContainer, MoleculeContainer]):
        """
        write single molecule or reaction into file
        """
        self._file.write('<cml>\n')
        self.__write(data)
        self.write = self.__write

    def __write(self, data):
        if isinstance(data, ReactionContainer):
            buffer = ['<MDocument><MChemicalStruct>']
            if not data._arrow:
                data.fix_positions()

            if data.name:
                buffer.append(f'<reaction title="{data.name}">')
            else:
                buffer.append('<reaction>')

            if data.meta:
                buffer.append('<propertyList>')
                for k, v in data.meta.items():
                    if isinstance(v, str):
                        v = f'<![CDATA[{v}]]>'
                    buffer.append(f'<property title="{k}"><scalar>{v}</scalar></property>')
                buffer.append('</propertyList>')
            c = count(1)
            for i, j in ((data.reactants, 'reactantList'), (data.products, 'productList'),
                         (data.reagents, 'agentList')):
                if not i:
                    continue
                buffer.append(f'<{j}>')
                for n, m in zip(c, i):
                    if m.name:
                        buffer.append(f'<molecule title="{m.name}" molID="m{n}">')
                    else:
                        buffer.append(f'<molecule molID="m{n}">')
                    buffer.append(self.__convert_structure(m))
                    buffer.append('</molecule>')
                buffer.append(f'</{j}>')

            buffer.append(f'<arrow type="DEFAULT" x1="{data._arrow[0] * 2:.4f}" y1="0" '
                          f'x2="{data._arrow[1] * 2:.4f}" y2="0"/>')
            buffer.append('</reaction>')
            self._file.writelines(buffer)
        elif not isinstance(data, MoleculeContainer):
            raise TypeError('MoleculeContainer expected')
        else:
            m = self.__convert_structure(data)
            self._file.write('<MDocument><MChemicalStruct>')

            if data.name:
                self._file.write(f'<molecule title="{data.name}">')
            else:
                self._file.write('<molecule>')

            if data.meta:
                self._file.write('<propertyList>')
                for k, v in data.meta.items():
                    if isinstance(v, str):
                        v = f'<![CDATA[{v}]]>'
                    self._file.write(f'<property title="{k}"><scalar>{v}</scalar></property>')
                self._file.write('</propertyList>')
            self._file.write(m)
            self._file.write('</molecule>')
        self._file.write('</MChemicalStruct></MDocument>\n')

    @staticmethod
    def __convert_structure(g):
        gp = g._plane
        gc = g._charges
        gr = g._radicals
        bg = g._bonds
        hg = g._hydrogens
        hb = g.hybridization

        out = ['<atomArray>']
        for n, atom in g._atoms.items():
            x, y = gp[n]
            ih = hg[n]
            out.append(f'<atom id="a{n}" elementType="{atom.atomic_symbol}" '
                       f'x2="{x * 2:.4f}" y2="{y * 2:.4f}" mrvMap="{n}"')
            if gc[n]:
                out.append(f' formalCharge="{gc[n]}"')
            if gr[n]:
                out.append(' radical="monovalent"')
            if atom.isotope:
                out.append(f' isotope="{atom.isotope}"')
            if ih and (atom.atomic_symbol not in organic_set or hb(n) == 4 and atom.atomic_number in (5, 7, 15)):
                out.append(f' hydrogenCount="{ih}"')
            out.append('/>')
        out.append('</atomArray>')

        out.append('<bondArray>')
        wedge = defaultdict(set)
        for n, (i, j, s) in enumerate(g._wedge_map, start=1):
            out.append(f'<bond id="b{n}" atomRefs2="a{i} a{j}" order="{bond_map[bg[i][j].order]}">'
                       f'<bondStereo>{s == 1 and "W" or "H"}</bondStereo></bond>')
            wedge[i].add(j)
            wedge[j].add(i)
        for n, (i, j, bond) in enumerate(g.bonds(), start=len(out)):
            if j not in wedge[i]:
                out.append(f'<bond id="b{n}" atomRefs2="a{i} a{j}" order="{bond_map[bond.order]}"/>')
        out.append('</bondArray>')
        return ''.join(out)

    __finalized = False


__all__ = ['MRVRead', 'MRVWrite']
