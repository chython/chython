# -*- coding: utf-8 -*-
#
#  Copyright 2018-2021 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2021 Aleksandr Sizov <murkyrussian@gmail.com>
#  Copyright 2019 Artem Mukanov <nostro32@mail.ru>
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
from functools import reduce
from io import StringIO, TextIOWrapper
from itertools import chain
from operator import or_
from pathlib import Path
from re import split, compile, fullmatch, findall, search
from ._parser import DaylightParser
from ...containers import ReactionContainer
from ...exceptions import IncorrectSmiles, IsChiral, NotChiral, ValenceError

cx_fragments = compile(r'f:(?:[0-9]+(?:\.[0-9]+)+)(?:,(?:[0-9]+(?:\.[0-9]+)+))*')
cx_radicals = compile(r'\^[1-7]:[0-9]+(?:,[0-9]+)*')
delimiter = compile(r'[=:]')


class SMILESRead(DaylightParser):
    """SMILES separated per lines files reader. Works similar to opened file object. Support `with` context manager.
    On initialization accept opened in text mode file, string path to file,
    pathlib.Path object or another buffered reader object.

    Line should be start with SMILES string and optionally continues with space/tab separated list of
    `key:value` [or `key=value`] data if `header=None`. For example::

        C=C>>CC id:123 key=value

    if `header=True` then first line of file should be space/tab separated list of keys including smiles column key.
    For example::

        ignored_smi_key key1 key2
        CCN 1 2

    Also possible to pass list of keys (without smiles_pseudo_key) for mapping space/tab separated list
    of SMILES and values: `header=['key1', 'key2'] # order depended`.

    For reactions . [dot] in bonds should be used only for molecules separation.
    """
    def __init__(self, file, header=None, ignore_stereo=False, **kwargs):
        """
        :param ignore: Skip some checks of data or try to fix some errors.
        :param remap: Remap atom numbers started from one.
        :param store_log: Store parser log if exists messages to `.meta` by key `ParserLog`.
        :param ignore_stereo: Ignore stereo data.
        """
        if isinstance(file, str):
            self._file = open(file)
            self.__is_buffer = False
        elif isinstance(file, Path):
            self._file = file.open()
            self.__is_buffer = False
        elif isinstance(file, (TextIOWrapper, StringIO, FileInput)):
            self._file = file
            self.__is_buffer = True
        else:
            raise TypeError('invalid file. TextIOWrapper, StringIO subclasses possible')
        super().__init__(**kwargs)
        self.__file = iter(self._file.readline, '')

        if header is True:
            self.__header = next(self.__file).split()[1:]
        elif header:
            if not isinstance(header, (list, tuple)) or not all(isinstance(x, str) for x in header):
                raise TypeError('expected list (tuple) of strings')
            self.__header = header
        else:
            self.__header = None

        self.__ignore_stereo = ignore_stereo
        self._data = self.__data()
    def _create_molecule(self, data, mapping):
        mol = super()._create_molecule(data, mapping)
        hydrogens = mol._hydrogens
        radicals = mol._radicals
        calc_implicit = mol._calc_implicit
        for n, h in data['hydrogens'].items():
            n = mapping[n]
            hc = hydrogens[n]
            if hc is None:  # aromatic rings or valence errors. just store given H count.
                hydrogens[n] = h
            elif hc != h:  # H count mismatch. try radical state of atom.
                if radicals[n]:
                    if self._ignore:
                        hydrogens[n] = h  # set parsed hydrogens count
                        self._info(f'implicit hydrogen count ({h}) mismatch with '
                                   f'calculated ({hc}) on atom {n}. calculated count replaced.')
                    else:
                        raise ValueError(f'implicit hydrogen count ({h}) mismatch with '
                                         f'calculated ({hc}) on atom {n}.')
                else:
                    radicals[n] = True
                    calc_implicit(n)
                    if hydrogens[n] != h:  # radical state also has errors.
                        if self._ignore:
                            radicals[n] = False  # reset radical state
                            hydrogens[n] = h  # set parsed hydrogens count
                            self._info(f'implicit hydrogen count ({h}) mismatch with '
                                       f'calculated ({hc}) on atom {n}. calculated count replaced.')
                        else:
                            raise ValueError(f'implicit hydrogen count ({h}) mismatch with '
                                             f'calculated ({hc}) on atom {n}.')

        if self.__ignore_stereo or not data['stereo_atoms'] and not data['stereo_bonds']:
            return mol

        st = mol._stereo_tetrahedrons
        sa = mol._stereo_allenes
        sat = mol._stereo_allenes_terminals
        ctt = mol._stereo_cis_trans_terminals

        order = {mapping[n]: [mapping[m] for m in ms] for n, ms in data['order'].items()}

        stereo = []
        for i, s in data['stereo_atoms'].items():
            n = mapping[i]
            if not i and hydrogens[n]:  # first atom in smiles has reversed chiral mark
                s = not s

            if n in st:
                stereo.append((mol.add_atom_stereo, n, order[n], s))
            elif n in sa:
                t1, t2 = sat[n]
                env = sa[n]
                n1 = next(x for x in order[t1] if x in env)
                n2 = next(x for x in order[t2] if x in env)
                stereo.append((mol.add_atom_stereo, n, (n1, n2), s))

        stereo_bonds = {mapping[n]: {mapping[m]: s for m, s in ms.items()}
                        for n, ms in data['stereo_bonds'].items()}
        seen = set()
        for n, ns in stereo_bonds.items():
            if n in seen:
                continue
            if n in ctt:
                nm = ctt[n]
                m = nm[1] if nm[0] == n else nm[0]
                if m in stereo_bonds:
                    seen.add(m)
                    n2, s2 = stereo_bonds[m].popitem()
                    n1, s1 = ns.popitem()
                    stereo.append((mol.add_cis_trans_stereo, n, m, n1, n2, s1 == s2))

        while stereo:
            fail_stereo = []
            old_stereo = len(stereo)
            for f, *args in stereo:
                try:
                    f(*args, clean_cache=False)
                except NotChiral:
                    fail_stereo.append((f, *args))
                except IsChiral:
                    pass
                except ValenceError:
                    self._info('structure has errors, stereo data skipped')
                    mol.flush_cache()
                    break
            else:
                stereo = fail_stereo
                if len(stereo) == old_stereo:
                    break
                mol.flush_stereo_cache()
                continue
            break
        return mol

    def parse(self, smiles: str):
        """SMILES string parser."""
        if not smiles:
            raise ValueError('Empty string')

        smi, *data = smiles.split()
        if data and data[0].startswith('|') and data[0].endswith('|'):
            fr = search(cx_fragments, data[0])
            if fr is not None:
                contract = [sorted(int(x) for x in x.split('.')) for x in fr.group()[2:].split(',')]
                if len({x for x in contract for x in x}) < len([x for x in contract for x in x]):
                    self._info(f'collisions in cxsmiles fragments description: {data[0]}')
                    contract = None
                elif any(x[0] < 0 for x in contract):
                    self._info(f'invalid cxsmiles fragments description: {data[0]}')
                    contract = None

                radicals = [int(x) for x in findall(cx_radicals, data[0]) for x in x[3:].split(',')]
                if any(x < 0 for x in radicals):
                    self._info(f'invalid cxsmiles radicals description: {data[0]}')
                    radicals = []
                if len(set(radicals)) != len(radicals):
                    self._info(f'collisions in cxsmiles radicals description: {data[0]}')
                    radicals = []
                data = data[1:]
            else:
                radicals = [int(x) for x in findall(cx_radicals, data[0]) for x in x[3:].split(',')]
                if radicals:
                    if any(x < 0 for x in radicals):
                        self._info(f'invalid cxsmiles radicals description: {data[0]}')
                        radicals = []
                    if len(set(radicals)) != len(radicals):
                        self._info(f'collisions in cxsmiles radicals description: {data[0]}')
                        radicals = []
                    data = data[1:]
                contract = None
        else:
            radicals = []
            contract = None

        if self.__header is None:
            meta = {}
            for x in data:
                try:
                    k, v = split(delimiter, x, 1)
                    meta[k] = v
                except ValueError:
                    self._info(f'invalid metadata entry: {x}')
        else:
            meta = dict(zip(self.__header, data))

        if '>' in smi and (smi[smi.index('>') + 1] in '>([' or smi[smi.index('>') + 1].isalpha()):
            record = {'reactants': [], 'reagents': [], 'products': [], 'meta': meta, 'title': ''}
            try:
                reactants, reagents, products = smi.split('>')
            except ValueError as e:
                raise ValueError('invalid reaction smiles') from e

            for k, d in zip(('reactants', 'products', 'reagents'), (reactants, products, reagents)):
                if d:
                    for x in d.split('.'):
                        if not x:
                            if self._ignore:
                                self._info('two dots in line ignored')
                            else:
                                raise ValueError('invalid reaction smiles. two dots in line')
                        else:
                            record[k].append(self.__parse_tokens(x))

            if radicals:
                atom_map = dict(enumerate(a for m in chain(record['reactants'], record['reagents'], record['products'])
                                          for a in m['atoms']))
                for x in radicals:
                    atom_map[x]['is_radical'] = True

            container = self._convert_reaction(record)
            if contract:
                if max(x for x in contract for x in x) >= len(container):
                    self._info(f'skipped invalid contract data: {contract}')
                lr = len(container.reactants)
                reactants = set(range(lr))
                reagents = set(range(lr, lr + len(container.reagents)))
                products = set(range(lr + len(container.reagents),
                                     lr + len(container.reagents) + len(container.products)))
                new_molecules = [None] * len(container)
                for c in contract:
                    if reactants.issuperset(c):
                        new_molecules[c[0]] = reduce(or_, (container.reactants[x] for x in c))
                        reactants.difference_update(c)
                    elif products.issuperset(c):
                        new_molecules[c[0]] = reduce(or_, (container.products[x - len(container)] for x in c))
                        products.difference_update(c)
                    elif reagents.issuperset(c):
                        new_molecules[c[0]] = reduce(or_, (container.reagents[x - lr] for x in c))
                        reagents.difference_update(c)
                    else:
                        self._info(f'impossible to contract different parts of reaction: {contract}')
                for x in reactants:
                    new_molecules[x] = container.reactants[x]
                for x in products:
                    new_molecules[x] = container.products[x - len(container)]
                for x in reagents:
                    new_molecules[x] = container.reagents[x - lr]

                meta = container.meta
                if self._store_log:
                    if log := self._format_log():
                        if 'ParserLog' in meta:
                            meta['ParserLog'] += '\n' + log
                        else:
                            meta['ParserLog'] = log
                else:
                    self._flush_log()
                return ReactionContainer([x for x in new_molecules[:lr] if x is not None],
                                         [x for x in new_molecules[-len(container.products):] if x is not None],
                                         [x for x in new_molecules[lr: -len(container.products)] if x is not None],
                                         meta=meta)
            return container
        else:
            record = self.__parse_tokens(smi)
            record['meta'].update(meta)
            if 'cgr' in record:  # CGR smiles parser
                return self._convert_cgr(record)
            for x in radicals:
                record['atoms'][x]['is_radical'] = True
            return self._convert_molecule(record)


__all__ = ['SMILESRead']
