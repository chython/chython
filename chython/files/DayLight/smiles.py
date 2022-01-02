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
from re import compile, fullmatch
from ...exceptions import IncorrectSmiles


atom_re = compile(r'([1-9][0-9]{0,2})?([A-IK-PR-Zacnopsbt][a-ik-pr-vy]?)(@@|@)?(H[1-4]?)?([+-][1-4+-]?)?(:[0-9]{1,4})?')
charge_dict = {'+': 1, '+1': 1, '++': 2, '+2': 2, '+3': 3, '+++': 3, '+4': 4, '++++': 4,
               '-': -1, '-1': -1, '--': -2, '-2': -2, '-3': -3, '---': -3, '-4': -4, '----': -4}


def _atom_parse(token):
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


def _process_tokens(tokens):
    out = []
    for token_type, token in tokens:
        if token_type in (0, 8):  # simple atom
            out.append((token_type, {'element': token, 'charge': 0, 'isotope': None, 'is_radical': False,
                                     'mapping': 0, 'x': 0., 'y': 0., 'z': 0., 'hydrogen': None, 'stereo': None}))
        elif token_type == 5:
            out.append(_atom_parse(token))
        else:
            out.append((token_type, token))
    return out


__all__ = []
