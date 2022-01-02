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
from ...exceptions import IncorrectSmiles


# -,= OR bonds supported
# !: NOT bonds supported
# ~ ANY bonds supported
# @ bond not supported
# @;!: and any other complicated combinations not supported


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
# 13: query OR bond
# 15: query NOT bond


replace_dict = {'-': 1, '=': 2, '#': 3, ':': 4, '~': 8}
not_dict = {'-': [2, 3, 4], '=': [1, 3, 4], '#': [1, 2, 4], ':': [1, 2, 3]}


def tokenize(smiles):
    token_type = token = None
    tokens = []
    for s in smiles:
        if s == '[':  # open complex token
            if token_type == 5:  # two opened [
                raise IncorrectSmiles('[..[')
            elif token_type in (13, 15):
                raise IncorrectSmiles('SMARTS query bond invalid')
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
        elif token_type == 5:  # grow token with brackets. skip validation
            token.append(s)
        elif s == '(':
            if token_type in (13, 15):
                raise IncorrectSmiles('SMARTS query bond invalid')
            elif token:
                tokens.append((token_type, token))
                token = None
            elif token_type == 7:  # empty closure
                raise IncorrectSmiles('invalid closure')
            elif token_type == 2:  # barely opened
                raise IncorrectSmiles('((')
            token_type = 2
            tokens.append((2, None))
        elif s == ')':
            if token_type in (13, 15):
                raise IncorrectSmiles('SMARTS query bond invalid')
            elif token:
                tokens.append((token_type, token))
                token = None
            elif token_type == 7:  # empty closure
                raise IncorrectSmiles('invalid closure')
            elif token_type == 2:  # barely opened
                raise IncorrectSmiles('()')
            token_type = 3
            tokens.append((3, None))
        elif s.isnumeric():  # closures
            if token_type in (13, 15):
                raise IncorrectSmiles('SMARTS query bond invalid')
            elif token_type == 7:  # % already found. collect number
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
            if token_type in (13, 15):
                raise IncorrectSmiles('SMARTS query bond invalid')
            elif token:
                tokens.append((token_type, token))
            elif token_type == 7:
                raise IncorrectSmiles('%%')
            token_type = 7
            token = []
        elif s in '=#:-~':  # bonds found
            if token_type == 13:
                token.append(replace_dict[s])
                tokens.append((13, token))
                token_type = token = None  # finalize token
            elif token_type == 15:
                token_type = None
                tokens.append((13, not_dict[s]))
            elif token:
                tokens.append((token_type, token))
                token = None
            else:
                token_type = 1
                tokens.append((1, replace_dict[s]))
        elif s in r'\/':
            if token_type in (13, 15):
                raise IncorrectSmiles('SMARTS query bond invalid')
            elif token:
                tokens.append((token_type, token))
                token = None
            token_type = 9
            tokens.append((9, s == '/'))  # Up is true
        elif s == '.':
            if token_type in (13, 15):
                raise IncorrectSmiles('SMARTS query bond invalid')
            elif token:
                tokens.append((token_type, token))
                token = None
            token_type = 4
            tokens.append((4, None))
        elif s in 'NOPSFI':  # organic atoms
            if token_type in (13, 15):
                raise IncorrectSmiles('SMARTS query bond invalid')
            elif token:
                tokens.append((token_type, token))
                token = None
            token_type = 0
            tokens.append((0, s))
        elif s in 'cnopsb':  # aromatic ring atom
            if token_type in (13, 15):
                raise IncorrectSmiles('SMARTS query bond invalid')
            elif token:
                tokens.append((token_type, token))
                token = None
            token_type = 8
            tokens.append((8, s.upper()))
        elif s in 'CB':  # flag possible Cl or Br
            if token_type in (13, 15):
                raise IncorrectSmiles('SMARTS query bond invalid')
            elif token:
                tokens.append((token_type, token))
            token_type = 0
            token = s
        elif s == ',':  # query bond separator
            if token_type != 1:
                raise IncorrectSmiles('SMARTS query bond invalid')
            token_type = 13
            token = [tokens.pop(-1)[1]]
        elif s == '!':  # query not bond
            if token_type not in (0, 2, 3, 6, 8, None):
                raise IncorrectSmiles('SMARTS query bond invalid')
            token_type = 15
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
    return [(6, int(''.join(y))) if x == 7 else (x, y) for x, y in tokens]  # composite closures folding


__all__ = ['tokenize']
