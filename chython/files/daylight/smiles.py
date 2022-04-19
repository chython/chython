# -*- coding: utf-8 -*-
#
#  Copyright 2022 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from itertools import chain
from pathlib import Path
from re import compile, findall, search, split
from typing import List, Union, Optional
from .parser import parser
from .tokenize import smiles_tokenize
from .._mdl import Parser
from ...containers import MoleculeContainer, ReactionContainer
from ...exceptions import IsChiral, NotChiral, ValenceError, ParseError


delimiter = compile(r'[=:]')
cx_fragments = compile(r'f:(?:[0-9]+(?:\.[0-9]+)+)(?:,(?:[0-9]+(?:\.[0-9]+)+))*')
cx_radicals = compile(r'\^[1-7]:[0-9]+(?:,[0-9]+)*')

implicit_mismatch_key = 'implicit_mismatch'


class SMILESRead(Parser):
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
    def __init__(self, file, header=None, ignore_stereo=False, keep_implicit=False, ignore_carbon_radicals=False,
                 **kwargs):
        """
        :param ignore: Skip some checks of data or try to fix some errors.
        :param remap: Remap atom numbers started from one.
        :param store_log: Store parser log if exists messages to `.meta` by key `ParserLog`.
        :param ignore_stereo: Ignore stereo data.
        :param keep_implicit: keep given in smiles implicit hydrogen count, otherwise ignore on valence error.
        :param ignore_bad_isotopes: reset invalid isotope mark to non-isotopic.
        :param ignore_carbon_radicals: fill carbon radicals with hydrogen (X[C](X)X case).
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
        self.__keep_implicit = keep_implicit
        self.__ignore_carbon_radicals = ignore_carbon_radicals
        self._data = self.__data()

    def __data(self):
        file = self._file
        parse = self.parse
        try:
            seekable = file.seekable()
        except AttributeError:
            seekable = False
        pos = file.tell() if seekable else None
        for n, line in enumerate(self.__file):
            try:
                x = parse(line)
            except ValueError:  # yield exception, not raise
                yield ParseError(n, pos, self._format_log(), line)
            else:
                yield x
            if seekable:
                pos = file.tell()

    @classmethod
    def create_parser(cls, header=None, ignore_stereo=False, keep_implicit=False, ignore_carbon_radicals=False,
                      **kwargs):
        """
        Create SMILES parser function configured same as SMILESRead object.
        """
        obj = object.__new__(cls)
        obj._SMILESRead__header = header
        obj._SMILESRead__ignore_stereo = ignore_stereo
        obj._SMILESRead__keep_implicit = keep_implicit
        obj._SMILESRead__ignore_carbon_radicals = ignore_carbon_radicals
        super(SMILESRead, obj).__init__(**kwargs)
        return obj.parse

    def close(self, force=False):
        """
        Close opened file.

        :param force: Force closing of externally opened file or buffer.
        """
        if not self.__is_buffer or force:
            self._file.close()

    def __enter__(self):
        return self

    def __exit__(self, _type, value, traceback):
        self.close()

    def read(self) -> List[Union[MoleculeContainer, ReactionContainer]]:
        """
        Parse whole file.

        :return: List of parsed molecules or reactions.
        """
        return list(iter(self))

    def __iter__(self):
        return (x for x in self._data if not isinstance(x, ParseError))

    def __next__(self):
        return next(iter(self))

    def parse(self, smiles: str) -> Union[MoleculeContainer, ReactionContainer]:
        """
        SMILES string parser. String should start with SMILES and
        optionally continues with space/tab separated list of key:value [or key=value] data.
        """
        if not smiles:
            raise ValueError('Empty string')

        contract: Optional[List[List[int]]]  # typing

        smi, *data = smiles.split()
        if data and data[0].startswith('|') and data[0].endswith('|'):
            fr = search(cx_fragments, data[0])
            if fr is not None:
                contract = [sorted(int(x) for x in x.split('.')) for x in fr.group()[2:].split(',')]
                if len({x for x in contract for x in x}) < len([x for x in contract for x in x]):
                    self._info(f'collisions in cxsmiles fragments description: {data[0]}')
                    contract = None

                radicals = [int(x) for x in findall(cx_radicals, data[0]) for x in x[3:].split(',')]
                if radicals and len(set(radicals)) != len(radicals):
                    self._info(f'collisions in cxsmiles radicals description: {data[0]}')
                    radicals = []
                data.pop(0)
            else:
                radicals = [int(x) for x in findall(cx_radicals, data[0]) for x in x[3:].split(',')]
                if radicals:
                    if len(set(radicals)) != len(radicals):
                        self._info(f'collisions in cxsmiles radicals description: {data[0]}')
                        radicals = []
                    data.pop(0)
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

        if '>' in smi:
            record = {'reactants': [], 'reagents': [], 'products': [], 'meta': meta, 'title': ''}
            try:
                reactants, reagents, products = smi.split('>')
            except ValueError as e:
                raise ValueError('invalid reaction smiles') from e

            mol_count = 0
            for k, d in zip(('reactants', 'products', 'reagents'), (reactants, products, reagents)):
                if not d:
                    continue
                for x in d.split('.'):
                    if not x:
                        if self._ignore:
                            self._info('two dots in line ignored')
                        else:
                            raise ValueError('invalid reaction smiles. two dots in line')
                    else:
                        record[k].append(x)
                        mol_count += 1

            if contract:
                if max(x for x in contract for x in x) >= mol_count:
                    self._info(f'skipped invalid contract data: {contract}')
                lr = len(record['reactants'])
                lp = len(record['products'])
                reactants = set(range(lr))
                reagents = set(range(lr, mol_count - lp))
                products = set(range(mol_count - lp, mol_count))
                new_molecules: List[Optional[str]] = [None] * mol_count
                for c in contract:
                    if reactants.issuperset(c):
                        new_molecules[c[0]] = '.'.join(record['reactants'][x] for x in c)
                        reactants.difference_update(c)
                    elif products.issuperset(c):
                        new_molecules[c[0]] = '.'.join(record['products'][x - mol_count] for x in c)
                        products.difference_update(c)
                    elif reagents.issuperset(c):
                        new_molecules[c[0]] = '.'.join(record['reagents'][x - lr] for x in c)
                        reagents.difference_update(c)
                    else:
                        self._info(f'impossible to contract different parts of reaction: {contract}')
                for x in reactants:
                    new_molecules[x] = record['reactants'][x]
                for x in products:
                    new_molecules[x] = record['products'][x - mol_count]
                for x in reagents:
                    new_molecules[x] = record['reagents'][x - lr]

                record['reactants'] = [x for x in new_molecules[:lr] if x is not None]
                record['products'] = [x for x in new_molecules[-lp:] if x is not None]
                record['reagents'] = [x for x in new_molecules[lr: -lp] if x is not None]

            for k in ('reactants', 'products', 'reagents'):
                tmp = []
                for x in record[k]:
                    r, log = parser(smiles_tokenize(x), not self._ignore)
                    tmp.append(r)
                    for l in log:
                        self._info(l)
                record[k] = tmp

            if radicals:
                atom_map = dict(enumerate(a for m in chain(record['reactants'], record['reagents'], record['products'])
                                          for a in m['atoms']))
                for x in radicals:
                    atom_map[x]['is_radical'] = True

            return self._convert_reaction(record)
        else:
            record, log = parser(smiles_tokenize(smi), not self._ignore)
            for l in log:
                self._info(l)
            record['meta'].update(meta)
            for x in radicals:
                record['atoms'][x]['is_radical'] = True
            return self._convert_molecule(record)

    def _create_molecule(self, data, mapping):
        mol = super()._create_molecule(data, mapping)
        atoms = mol._atoms
        bonds = mol._bonds
        charges = mol._charges
        hydrogens = mol._hydrogens
        radicals = mol._radicals
        calc_implicit = mol._calc_implicit
        hyb = mol.hybridization
        keep_implicit = self.__keep_implicit
        radicalized = []

        for n, a in enumerate(data['atoms']):
            h = a['hydrogen']
            if h is None:  # simple atom token
                continue
            # bracket token should always contain implicit hydrogens count.
            n = mapping[n]
            if keep_implicit:  # override any calculated hydrogens count.
                hydrogens[n] = h
            elif (hc := hydrogens[n]) is None:  # atom has invalid valence for now.
                if hyb(n) == 4:  # this is aromatic rings. just store given H count.
                    hydrogens[n] = h
                    if not h and not charges[n] and atoms[n].atomic_number in (5, 6, 7, 15) and \
                            sum(b.order != 8 for b in bonds[n].values()) == 2:
                        # c[c]c - aromatic B,C,N,P radical
                        radicals[n] = True
                        radicalized.append(n)
                elif not radicals[n]:  # CXSMILES radical not set.
                    # SMILES doesn't code radicals. so, let's try to guess.
                    radicals[n] = True
                    calc_implicit(n)
                    if hydrogens[n] != h:  # radical state also has errors.
                        if self._ignore:
                            hydrogens[n] = None  # reset hydrogens
                            radicals[n] = False  # reset radical state
                            if implicit_mismatch_key in mol.meta:
                                mol.meta[implicit_mismatch_key][n] = h
                            else:
                                mol.meta[implicit_mismatch_key] = {n: h}
                            self._info(f'implicit hydrogen count ({h}) mismatch with calculated ({hc}) on atom {n}.')
                        else:
                            raise ValueError(f'implicit hydrogen count ({h}) mismatch with '
                                             f'calculated ({hc}) on atom {n}.')
                    else:
                        radicalized.append(n)
            elif hc != h and not radicals[n]:  # H count mismatch. try radical state of atom.
                radicals[n] = True
                calc_implicit(n)
                if hydrogens[n] != h:  # radical state also has errors.
                    if self._ignore:
                        hydrogens[n] = hc  # reset hydrogens to previous valid state
                        radicals[n] = False  # reset radical state
                        if implicit_mismatch_key in mol.meta:
                            mol.meta[implicit_mismatch_key][n] = h
                        else:
                            mol.meta[implicit_mismatch_key] = {n: h}
                        self._info(f'implicit hydrogen count ({h}) mismatch with calculated ({hc}) on atom {n}.')
                    else:
                        raise ValueError(f'implicit hydrogen count ({h}) mismatch with '
                                         f'calculated ({hc}) on atom {n}.')
                else:
                    radicalized.append(n)

        if self.__ignore_carbon_radicals:
            for n in radicalized:
                if atoms[n].atomic_number == 6:
                    radicals[n] = False
                    if (h := 4 - sum(b.order for b in bonds[n].values() if b.order != 8)) >= 0:
                        hydrogens[n] = h
                        self._info(f'Carbon radical replaced with ({h}) implicit hydrogens')

        if self.__ignore_stereo:
            return mol
        stereo_atoms = [(n, s) for n, a in enumerate(data['atoms']) if (s := a['stereo']) is not None]
        if not stereo_atoms and not data['stereo_bonds']:
            return mol

        st = mol._stereo_tetrahedrons
        sa = mol._stereo_allenes
        sat = mol._stereo_allenes_terminals
        ctt = mol._stereo_cis_trans_terminals

        order = {mapping[n]: [mapping[m] for m in ms] for n, ms in data['order'].items()}

        stereo = []
        for i, s in stereo_atoms:
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


__all__ = ['SMILESRead']
