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
from collections import defaultdict
from itertools import permutations
from re import compile, fullmatch
from ...containers import CGRContainer
from ...containers.bonds import DynamicBond
from ...exceptions import IncorrectSmiles
from ...periodictable import DynamicElement


# tokens structure:
# 10: dynamic bond
# 11: dynamic atom
# 12: dynamic aromatic atom


bonds_dict = {'.>-': (None, 1), '.>=': (None, 2), '.>#': (None, 3), '.>:': (None, 4), '.>~': (None, 8),
              '->.': (1, None), '->=': (1, 2), '->#': (1, 3), '->:': (1, 4), '->~': (1, 8),
              '=>.': (2, None), '=>-': (2, 1), '=>#': (2, 3), '=>:': (2, 4), '=>~': (2, 8),
              '#>.': (3, None), '#>-': (3, 1), '#>=': (3, 2), '#>:': (3, 4), '#>~': (3, 8),
              ':>.': (4, None), ':>-': (4, 1), ':>=': (4, 2), ':>#': (4, 3), ':>~': (4, 8),
              '~>.': (8, None), '~>-': (8, 1), '~>=': (8, 2), '~>#': (8, 3), '~>:': (8, 4)}
atom_re = compile(r'([1-9][0-9]{0,2})?([A-IK-PR-Zacnopsbt][a-ik-pr-vy]?)([+-0][1-4+-]?(>[+-0][1-4+-]?)?)?'
                  r'([*^](>[*^])?)?')
charge_dict = {'-4': -4, '-3': -3, '-2': -2, '--': -2, '-': -1, '0': 0, '+': 1, '++': 2, '+2': 2, '+3': 3, '+4': 4}
tmp = {f'{i}>{j}': (x, y) for (i, x), (j, y) in permutations(charge_dict.items(), 2) if x != y}
charge_dict = {k: (v, v) for k, v in charge_dict.items()}
charge_dict.update(tmp)
radical_dict = {'*': (True, True), '*>^': (True, False), '^>*': (False, True)}


def _atom_parse(token):
    # [isotope]Element[element][+-charge[>+-charge]][*^[>*^]]
    match = fullmatch(atom_re, token)
    if match is None:
        raise IncorrectSmiles(f'atom token invalid {{{token}}}')
    isotope, element, charge, _, is_radical, _ = match.groups()

    if isotope:
        isotope = int(isotope)

    if charge:
        try:
            charge, p_charge = charge_dict[charge]
        except KeyError:
            raise IncorrectSmiles('charge token invalid')
    else:
        charge = p_charge = 0

    if is_radical:
        try:
            is_radical, p_is_radical = radical_dict[is_radical]
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


def _process_tokens(tokens):
    out = []
    for token_type, token in tokens:
        if token_type in (0, 8):  # simple atom
            out.append((token_type, {'element': token, 'charge': 0, 'isotope': None, 'is_radical': False, 'p_charge': 0,
                                     'p_is_radical': False}))
        elif token_type == 5:
            if '>' in token and len(token) == 3:  # dynamic bond
                try:
                    out.append((10, bonds_dict[token]))
                except KeyError:
                    raise IncorrectSmiles(f'invalid dynamic bond token {{{token}}}')
            else:  # dynamic atom token
                out.append(_atom_parse(token))
        else:
            out.append((token_type, token))
    return out


def _convert_cgr(data):
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


__all__ = []
