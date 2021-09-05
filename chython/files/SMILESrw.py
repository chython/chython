# -*- coding: utf-8 -*-
#
#  Copyright 2018-2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from collections import defaultdict
from fileinput import FileInput
from functools import reduce
from itertools import permutations, chain
from io import StringIO, TextIOWrapper
from operator import or_
from pathlib import Path
from re import split, compile, fullmatch, findall, search
from typing import Union, List
from ._mdl import Parser, parse_error
from ..containers import MoleculeContainer, CGRContainer, ReactionContainer
from ..containers.bonds import DynamicBond
from ..exceptions import IncorrectSmiles, IsChiral, NotChiral, ValenceError
from ..periodictable import DynamicElement


# tokens structure:
# (type: int, value)
# types:
# 0: atom
# 1: bond
# 2: open chain (
# 3: close chain )
# 4: dot bond .
# 5: in bracket raw data []
# 6: closure number
# 7: raw closure number
# 8: aromatic atom
# 9: up down bond
# 10: dynamic bond
# 11: dynamic atom
# 12: dynamic aromatic atom

replace_dict = {'-': 1, '=': 2, '#': 3, ':': 4, '~': 8, '.': None, '(': 2, ')': 3}
charge_dict = {'+': 1, '+1': 1, '++': 2, '+2': 2, '+3': 3, '+++': 3, '+4': 4, '++++': 4,
               '-': -1, '-1': -1, '--': -2, '-2': -2, '-3': -3, '---': -3, '-4': -4, '----': -4}
dynamic_bonds = {'.>-': (None, 1), '.>=': (None, 2), '.>#': (None, 3), '.>:': (None, 4), '.>~': (None, 8),
                 '->.': (1, None), '->=': (1, 2), '->#': (1, 3), '->:': (1, 4), '->~': (1, 8),
                 '=>.': (2, None), '=>-': (2, 1), '=>#': (2, 3), '=>:': (2, 4), '=>~': (2, 8),
                 '#>.': (3, None), '#>-': (3, 1), '#>=': (3, 2), '#>:': (3, 4), '#>~': (3, 8),
                 ':>.': (4, None), ':>-': (4, 1), ':>=': (4, 2), ':>#': (4, 3), ':>~': (4, 8),
                 '~>.': (8, None), '~>-': (8, 1), '~>=': (8, 2), '~>#': (8, 3), '~>:': (8, 4)}
dyn_charge_dict = {'-4': -4, '-3': -3, '-2': -2, '--': -2, '-': -1, '0': 0, '+': 1, '++': 2, '+2': 2, '+3': 3, '+4': 4}
tmp = {f'{i}>{j}': (x, y) for (i, x), (j, y) in permutations(dyn_charge_dict.items(), 2) if x != y}
dyn_charge_dict = {k: (v, v) for k, v in dyn_charge_dict.items()}
dyn_charge_dict.update(tmp)
dyn_radical_dict = {'*': (True, True), '*>^': (True, False), '^>*': (False, True)}

atom_re = compile(r'([1-9][0-9]{0,2})?([A-IK-PR-Zacnopsbt][a-ik-pr-vy]?)(@@|@)?(H[1-4]?)?([+-][1-4+-]?)?(:[0-9]{1,4})?')
dyn_atom_re = compile(r'([1-9][0-9]{0,2})?([A-IK-PR-Zacnopsbt][a-ik-pr-vy]?)([+-0][1-4+-]?(>[+-0][1-4+-]?)?)?'
                      r'([*^](>[*^])?)?')
delimiter = compile(r'[=:]')
cx_fragments = compile(r'f:(?:[0-9]+(?:\.[0-9]+)+)(?:,(?:[0-9]+(?:\.[0-9]+)+))*')
cx_radicals = compile(r'\^[1-7]:[0-9]+(?:,[0-9]+)*')


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
            except ValueError:
                yield parse_error(n, pos, self._format_log(), line)
            else:
                yield x
            if seekable:
                pos = file.tell()

    @classmethod
    def create_parser(cls, header=None, ignore_stereo=False, *args, **kwargs):
        """
        Create SMILES parser function configured same as SMILESRead object.
        """
        obj = object.__new__(cls)
        obj._SMILESRead__header = header
        obj._SMILESRead__ignore_stereo = ignore_stereo
        super(SMILESRead, obj).__init__(*args, **kwargs)
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

    def read(self) -> List[Union[MoleculeContainer, CGRContainer, ReactionContainer]]:
        """
        Parse whole file.

        :return: List of parsed molecules or reactions.
        """
        return list(iter(self))

    def __iter__(self):
        return (x for x in self._data if not isinstance(x, parse_error))

    def __next__(self):
        return next(iter(self))

    def parse(self, smiles: str) -> Union[MoleculeContainer, CGRContainer, ReactionContainer]:
        """SMILES string parser."""
        if not smiles:
            raise ValueError('Empty string')
        self._flush_log()

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
                    log = self._format_log()
                    if log:
                        if 'ParserLog' in meta:
                            meta['ParserLog'] += '\n' + log
                        else:
                            meta['ParserLog'] = log
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

    def _convert_cgr(self, data):
        atoms = data['atoms']
        bonds = defaultdict(dict)

        for n, m, value in data['cgr']:
            bonds[n][m] = bonds[m][n] = DynamicBond(*value)

        g = object.__new__(CGRContainer)
        g_atoms = {}
        g_bonds = {}
        plane = {}
        charges = {}
        radicals = {}
        p_charges = {}
        p_radicals = {}
        for n, atom in enumerate(atoms, 1):
            g_atoms[n] = DynamicElement.from_symbol(atom['element'])(atom['isotope'])
            g_bonds[n] = {}
            charges[n] = atom['charge']
            radicals[n] = atom['is_radical']
            if 'p_charge' in atom:
                p_charges[n] = atom['p_charge']
                p_radicals[n] = atom['p_is_radical']
            else:
                p_charges[n] = atom['charge']
                p_radicals[n] = atom['is_radical']
            plane[n] = (0., 0.)
        for n, m, b in data['bonds']:
            if m in bonds[n]:
                if b != 8:
                    raise ValueError('CGR spec invalid')
                b = bonds[n][m]
            else:
                b = DynamicBond(b, b)
            n += 1
            m += 1
            if n in g_bonds[m]:
                raise ValueError('atoms already bonded')
            g_bonds[n][m] = g_bonds[m][n] = b
        g.__setstate__({'atoms': g_atoms, 'bonds': g_bonds, 'plane': plane, 'charges': charges, 'radicals': radicals,
                        'p_charges': p_charges, 'p_radicals': p_radicals, 'conformers': []})
        return g

    @staticmethod
    def _raw_tokenize(smiles):
        token_type = token = None
        tokens = []
        for s in smiles:
            if s == '[':  # open complex token
                if token_type == 5:  # two opened [
                    raise IncorrectSmiles('[...[')
                elif token:
                    tokens.append((token_type, token))
                elif token_type == 7:  # empty closure
                    raise IncorrectSmiles('invalid closure')
                token = []
                token_type = 5
            elif s == ']':  # close complex token
                if token_type != 5:
                    raise IncorrectSmiles(']..]')
                elif not token:
                    raise IncorrectSmiles('empty [] brackets')
                tokens.append((5, ''.join(token)))
                token_type = token = None
            elif s in '()':
                if token_type == 5:
                    raise IncorrectSmiles('brackets in atom/bond token')
                elif token:
                    tokens.append((token_type, token))
                    token = None
                elif token_type == 7:  # empty closure
                    raise IncorrectSmiles('invalid closure')
                elif token_type == 2:  # barely opened
                    raise IncorrectSmiles('(( or ()')
                token_type = replace_dict[s]
                tokens.append((token_type, None))
            elif token_type == 5:  # grow token with brackets. skip validation
                token.append(s)
            elif s.isnumeric():  # closures
                if token_type == 7:  # % already found. collect number
                    if not token and s == '0':
                        raise IncorrectSmiles('number starts with 0')
                    token.append(s)
                    if len(token) == 2:
                        tokens.append((token_type, token))
                        token_type = token = None
                else:
                    if s == '0':
                        raise IncorrectSmiles('number starts with 0')
                    elif token:
                        tokens.append((token_type, token))
                        token = None
                    token_type = 6
                    tokens.append((6, int(s)))
            elif s == '%':
                if token:
                    tokens.append((token_type, token))
                elif token_type == 7:
                    raise IncorrectSmiles('%%')
                token_type = 7
                token = []
            elif s in '=#:-~':  # bonds found
                if token:
                    tokens.append((token_type, token))
                    token = None
                token_type = 1
                tokens.append((1, replace_dict[s]))
            elif s in r'\/':
                if token:
                    tokens.append((token_type, token))
                    token = None
                token_type = 9
                tokens.append((9, s == '/'))  # Up is true
            elif s == '.':
                if token:
                    tokens.append((token_type, token))
                    token = None
                token_type = 4
                tokens.append((4, None))
            elif s in 'NOPSFI':  # organic atoms
                if token:
                    tokens.append((token_type, token))
                    token = None
                token_type = 0
                tokens.append((0, s))
            elif s in 'cnopsb':  # aromatic ring atom
                if token:
                    tokens.append((token_type, token))
                    token = None
                token_type = 8
                tokens.append((8, s.upper()))
            elif s in 'CB':  # flag possible Cl or Br
                if token:
                    tokens.append((token_type, token))
                token_type = 0
                token = s
            elif token_type == 0:
                if s == 'l':
                    if token == 'C':
                        tokens.append((0, 'Cl'))
                        token = None
                    else:
                        raise IncorrectSmiles('invalid element Bl')
                elif s == 'r':
                    if token == 'B':
                        tokens.append((0, 'Br'))
                        token = None
                    else:
                        raise IncorrectSmiles('invalid smiles for Cr')
                else:
                    raise IncorrectSmiles('invalid smiles')
            else:
                raise IncorrectSmiles('invalid smiles')

        if token_type == 5:
            raise IncorrectSmiles('atom description has not finished')
        elif token_type == 2:
            raise IncorrectSmiles('not closed')
        elif token:
            tokens.append((token_type, token))  # %closure or C or B
        return [(6, int(''.join(x[1]))) if x[0] == 7 else x for x in tokens]  # composite closures folding

    @classmethod
    def _fix_tokens(cls, tokens):
        out = []
        for token_type, token in tokens:
            if token_type in (0, 8):  # simple atom
                out.append((token_type, {'element': token, 'charge': 0, 'isotope': None, 'is_radical': False,
                                         'mapping': 0, 'x': 0., 'y': 0., 'z': 0., 'hydrogen': None, 'stereo': None}))
            elif token_type == 5:
                if '>' in token:  # dynamic bond or atom
                    if len(token) == 3:  # bond only possible
                        try:
                            out.append((10, dynamic_bonds[token]))
                        except KeyError:
                            raise IncorrectSmiles(f'invalid dynamic bond token {{{token}}}')
                    else:  # dynamic atom token
                        out.append(cls.__dynatom_parse(token))
                elif '*' in token:  # CGR atom radical mark
                    out.append(cls.__dynatom_parse(token))
                else:  # atom token
                    out.append(cls.__atom_parse(token))
            else:  # as is types: 1, 2, 3, 4, 6, 9
                out.append((token_type, token))
        return out

    @staticmethod
    def __atom_parse(token):
        # [isotope]Element[element][@[@]][H[n]][+-charge][:mapping]
        match = fullmatch(atom_re, token)
        if match is None:
            raise IncorrectSmiles(f'atom token invalid {{{token}}}')
        isotope, element, stereo, hydrogen, charge, mapping = match.groups()

        if isotope:
            isotope = int(isotope)

        if stereo:
            stereo = stereo == '@'

        if hydrogen:
            if len(hydrogen) > 1:
                hydrogen = int(hydrogen[1:])
            else:
                hydrogen = 1
        else:
            hydrogen = 0

        if charge:
            try:
                charge = charge_dict[charge]
            except KeyError:
                raise IncorrectSmiles('charge token invalid')
        else:
            charge = 0

        if mapping:
            try:
                mapping = int(mapping[1:])
            except ValueError:
                raise IncorrectSmiles('invalid mapping token')
        else:
            mapping = 0

        if element in ('c', 'n', 'o', 'p', 's', 'as', 'se', 'b', 'te'):
            _type = 8
            element = element.capitalize()
        else:
            _type = 0
        return _type, {'element': element, 'charge': charge, 'isotope': isotope, 'is_radical': False,
                       'mapping': mapping, 'x': 0., 'y': 0., 'z': 0., 'hydrogen': hydrogen, 'stereo': stereo}

    @staticmethod
    def __dynatom_parse(token):
        # [isotope]Element[element][+-charge[>+-charge]][*^[>*^]]
        match = fullmatch(dyn_atom_re, token)
        if match is None:
            raise IncorrectSmiles(f'atom token invalid {{{token}}}')
        isotope, element, charge, _, is_radical, _ = match.groups()

        if isotope:
            isotope = int(isotope)

        if charge:
            try:
                charge, p_charge = dyn_charge_dict[charge]
            except KeyError:
                raise IncorrectSmiles('charge token invalid')
        else:
            charge = p_charge = 0

        if is_radical:
            try:
                is_radical, p_is_radical = dyn_radical_dict[is_radical]
            except KeyError:
                raise IncorrectSmiles('invalid dynamic radical token')
        else:
            is_radical = p_is_radical = False

        if element in ('c', 'n', 'o', 'p', 's', 'as', 'se', 'b', 'te'):
            _type = 12
            element = element.capitalize()
        else:
            _type = 11
        return _type, {'element': element, 'charge': charge, 'isotope': isotope, 'is_radical': is_radical,
                       'p_charge': p_charge, 'p_is_radical': p_is_radical}

    def __parse_tokens(self, smiles):
        tokens = self._raw_tokenize(smiles)
        tokens = self._fix_tokens(tokens)
        return self._parse_tokens(tokens)

    def _parse_tokens(self, tokens):
        strong_cycle = not self._ignore
        t1 = tokens[0][0]
        if t1 == 2:
            if tokens[1][0] not in (0, 8, 11, 12):
                raise IncorrectSmiles('not atom started')
        elif t1 not in (0, 8, 11, 12):
            raise IncorrectSmiles('not atom started')

        atoms = []
        bonds = []
        order = defaultdict(list)
        atoms_types = []
        atom_num = 0
        last_num = 0
        stack = []
        cycles = {}
        used_cycles = set()
        cgr = []
        stereo_bonds = defaultdict(dict)
        stereo_atoms = {}
        hydrogens = {}
        previous = None

        for token_type, token in tokens:
            if token_type == 2:  # ((((((
                if previous:
                    if previous[0] != 4:
                        raise IncorrectSmiles('bond before side chain')
                    previous = None
                stack.append(last_num)
            elif token_type == 3:  # ))))))
                if previous:
                    raise IncorrectSmiles('bond before closure')
                try:
                    last_num = stack.pop()
                except IndexError:
                    raise IncorrectSmiles('close chain more than open')
            elif token_type in (1, 4, 9, 10):  # bonds. only keeping for atoms connecting
                if previous:
                    raise IncorrectSmiles('2 bonds in a row')
                elif not atoms:
                    raise IncorrectSmiles('started from bond')
                previous = (token_type, token)
            elif token_type == 6:  # cycle
                if previous and previous[0] == 4:
                    raise IncorrectSmiles('dot-cycle pattern invalid')
                elif token not in cycles:
                    if token in used_cycles:
                        if strong_cycle:
                            raise IncorrectSmiles('reused closure number')
                        else:
                            self._info(f'reused closure number: {token}')
                    else:
                        used_cycles.add(token)
                    cycles[token] = (last_num, previous, len(order[last_num]))
                    order[last_num].append(None)  # Reserve a table
                else:
                    a, ob, ind = cycles[token]
                    if ob:
                        if not previous:
                            bt, b = ob
                            if bt == 9:  # closure open is \/ bonded
                                stereo_bonds[a][last_num] = b
                                bt = b = 1
                            elif strong_cycle:
                                raise IncorrectSmiles('not equal cycle bonds')
                        else:
                            bt, b = previous
                            obt, ob = ob
                            if bt == 9:  # \/ bonds can be unequal
                                if obt == 9:
                                    stereo_bonds[a][last_num] = ob
                                elif ob != 1:
                                    raise IncorrectSmiles('not equal cycle bonds')
                                stereo_bonds[last_num][a] = b
                                bt = b = 1
                            elif obt == 9:
                                if b != 1:
                                    raise IncorrectSmiles('not equal cycle bonds')
                                stereo_bonds[a][last_num] = ob
                            elif b != ob:
                                raise IncorrectSmiles('not equal cycle bonds')
                    elif previous:
                        bt, b = previous
                        if bt == 9:  # stereo \/
                            stereo_bonds[last_num][a] = b
                            bt = b = 1
                        elif strong_cycle:
                            raise IncorrectSmiles('not equal cycle bonds')
                    else:
                        bt = 1
                        b = 4 if atoms_types[last_num] in (8, 12) and atoms_types[a] in (8, 12) else 1

                    if bt == 1:
                        bonds.append((last_num, a, b))
                    else:  # bt == 10
                        bonds.append((last_num, a, 8))
                        cgr.append((last_num, a, b))
                    order[a][ind] = last_num
                    order[last_num].append(a)
                    del cycles[token]
                previous = None
            else:  # atom
                if atoms:
                    if not previous:
                        bt = 1
                        b = 4 if token_type in (8, 12) and atoms_types[last_num] in (8, 12) else 1
                        order[last_num].append(atom_num)
                        order[atom_num].append(last_num)
                    else:
                        bt, b = previous
                        if bt != 4:
                            order[last_num].append(atom_num)
                            order[atom_num].append(last_num)
                    if bt == 1:
                        bonds.append((atom_num, last_num, b))
                    elif bt == 9:
                        bonds.append((atom_num, last_num,
                                      4 if token_type in (8, 12) and atoms_types[last_num] in (8, 12) else 1))
                        stereo_bonds[last_num][atom_num] = b
                        stereo_bonds[atom_num][last_num] = not b
                    elif bt == 10:
                        bonds.append((atom_num, last_num, 8))
                        cgr.append((atom_num, last_num, b))

                if token_type not in (11, 12):
                    stereo = token.pop('stereo')
                    if stereo is not None:
                        stereo_atoms[atom_num] = stereo
                    hydrogen = token.pop('hydrogen')
                    if hydrogen is not None:
                        hydrogens[atom_num] = hydrogen

                atoms.append(token)
                atoms_types.append(token_type)

                last_num = atom_num
                atom_num += 1
                previous = None

        if stack:
            raise IncorrectSmiles('number of ( does not equal to number of )')
        elif cycles:
            raise IncorrectSmiles('cycle is not finished')
        elif previous:
            raise IncorrectSmiles('bond on the end')

        stereo_bonds = {n: ms for n, ms in stereo_bonds.items() if len(ms) == 1 or len(ms) == set(ms.values())}
        mol = {'atoms': atoms, 'bonds': bonds, 'order': order,
               'stereo_bonds': stereo_bonds, 'stereo_atoms': stereo_atoms, 'hydrogens': hydrogens, 'meta': {}}
        if cgr or any(x in (11, 12) for x in atoms_types):
            mol['cgr'] = cgr
        return mol


__all__ = ['SMILESRead']
