# -*- coding: utf-8 -*-
#
#  Copyright 2022-2025 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from re import compile, match, search
from ...containers.bonds import QueryBond
from ...exceptions import IncorrectSmiles, IncorrectSmarts


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

# 10: query OR bond
# 11: query NOT bond
# 12: in ring bond


iso_re = compile(r'^[0-9]+')
chg_re = compile(r'[+-][1-4+-]?')
mpp_re = compile(r':[1-9][0-9]*$')
str_re = compile(r'@[@?]?')

replace_dict = {'-': 1, '=': 2, '#': 3, ':': 4, '~': 8}
not_dict = {'-': [2, 3, 4], '=': [1, 3, 4], '#': [1, 2, 4], ':': [1, 2, 3]}
atom_re = compile(r'([1-9][0-9]{0,2})?([A-IK-PR-Zacnopsbt][a-ik-pr-vy]?)(@@|@)?(H[1-4]?)?([+-][1-4+-]?)?(:[0-9]{1,4})?')
charge_dict = {'+': 1, '+1': 1, '++': 2, '+2': 2, '+3': 3, '+++': 3, '+4': 4, '++++': 4,
               '-': -1, '-1': -1, '--': -2, '-2': -2, '-3': -3, '---': -3, '-4': -4, '----': -4}


def _tokenize(smiles):
    token_type = token = None
    tokens = []
    for s in smiles:
        if token_type == 12:  # -;, =;, #; or :; found.
            if s == '!':  # not in ring
                if token:  # !! case
                    raise IncorrectSmarts('Invalid ring bond token')
                token = True
            elif s == '@':
                tokens.append((12, QueryBond(tokens.pop(-1)[1], not token)))
                token_type = token = None  # finalize [not]ring bond
            else:
                raise IncorrectSmarts('Invalid ring bond token')
        # [atom block parser]
        elif s == '[':  # open complex token
            if token_type == 5:  # two opened [
                raise IncorrectSmiles('[..[')
            elif token_type in (10, 11):
                raise IncorrectSmarts('Query bond invalid')
            elif token_type == 7:  # empty closure
                raise IncorrectSmiles('invalid closure')
            elif token:
                tokens.append((token_type, token))
            token = []
            token_type = 5
        elif s == ']':  # close complex token
            if token_type != 5:
                raise IncorrectSmiles(']..]')
            elif not token:
                raise IncorrectSmiles('empty [] brackets')
            tokens.append((5, ''.join(token)))
            token = None
            token_type = 0  # mark as atom
        elif token_type == 5:  # grow token with brackets. skip validation
            token.append(s)
        # closure parser
        elif s.isnumeric():  # closures
            if token_type in (10, 11):
                raise IncorrectSmarts('Query bond invalid')
            elif token_type == 2:
                raise IncorrectSmiles('(1 case invalid')
            elif token_type == 7:  # % already found. collect number
                if not token and s == '0':
                    raise IncorrectSmiles('number starts with 0')
                token.append(s)
                if len(token) == 2:
                    tokens.append((6, int(''.join(token))))
                    token = None
                    token_type = 6  # mark finished
            else:
                if s == '0':
                    raise IncorrectSmiles('number starts with 0')
                elif token:
                    tokens.append((token_type, token))
                    token = None
                token_type = 6
                tokens.append((6, int(s)))
        elif token_type == 7:
            raise IncorrectSmiles('expected closure number')
        elif s == '%':
            if token_type in (10, 11):
                raise IncorrectSmarts('Query bond invalid')
            elif token_type == 2:
                raise IncorrectSmiles('(%10 case invalid')
            elif token:
                tokens.append((token_type, token))
            token_type = 7
            token = []
        # bonds parser
        elif s in '=#:-~':  # bonds found
            if token_type == 10:
                token.append(replace_dict[s])
                tokens.append((10, token))
                token_type = token = None  # finalize token
            elif token_type == 11:
                token_type = None
                tokens.append((10, not_dict[s]))
            else:
                if token:
                    tokens.append((token_type, token))
                    token = None
                token_type = 1
                tokens.append((1, replace_dict[s]))
        elif token_type in (10, 11):  # expected bond symbol
            raise IncorrectSmarts('query bond invalid')
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
        elif s == ';':  # ;@ or ;!@ - in ring or not in ring bond
            if token_type is not None and token_type != 1:  # start of smiles or bond, list of bonds and not bond.
                raise IncorrectSmarts('Ring bond token invalid')
            token_type = 12
        elif s == ',':  # query bond separator
            if token_type != 1:
                raise IncorrectSmarts('Query bond invalid')
            token_type = 10
            token = [tokens.pop(-1)[1]]
        elif s == '!':  # query not bond
            if token_type not in (0, 2, 3, 6, 8):  # closures, brackets or atoms expected
                raise IncorrectSmarts('Query bond invalid')
            elif token:
                tokens.append((token_type, token))
                token = None
            token_type = 11
        # brackets
        elif s == '(':
            if token_type == 2:  # barely opened
                raise IncorrectSmiles('((')
            elif token:
                tokens.append((token_type, token))
                token = None
            token_type = 2
            tokens.append((2, None))
        elif s == ')':
            if token_type == 2:  # barely opened
                raise IncorrectSmiles('()')
            elif token:
                tokens.append((token_type, token))
                token = None
            token_type = 3
            tokens.append((3, None))
        # simple atoms
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
    elif token_type == 7:
        if token:
            tokens.append((6, int(token[0])))
        else:
            raise IncorrectSmiles('invalid %closure')
    elif token:
        tokens.append((token_type, token))  # C or B
    return tokens


def _atom_parse(token):
    # [isotope]Element[element][@[@]][H[n]][+-charge][:mapping]
    _match = atom_re.fullmatch(token)
    if _match is None:
        raise IncorrectSmiles(f'atom token invalid {token}')
    isotope, element, stereo, hydrogen, charge, mapping = _match.groups()

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

    if element in ('c', 'n', 'o', 'p', 's', 'as', 'se', 'b', 'te'):
        _type = 8
        element = element.capitalize()
    else:
        _type = 0
    return _type, {'element': element, 'isotope': isotope, 'parsed_mapping': mapping, 'charge': charge,
                   'implicit_hydrogens': hydrogen, 'stereo': stereo}


def _query_parse(token):
    out = {}
    if isotope := match(iso_re, token):
        token = token[isotope.end():]  # remove isotope substring
        out['isotope'] = int(isotope.group())
    if charge := search(chg_re, token):
        token = token[:charge.start()] + token[charge.end():]  # remove charge substring
        out['charge'] = charge_dict[charge.group()]

    if mapping := search(mpp_re, token):
        token = token[:mapping.start()]
        out['parsed_mapping'] = int(mapping.group()[1:])

    if stereo := search(str_re, token):  # drop stereo mark. unsupported
        token = token[:stereo.start()] + token[stereo.end():]
        out['stereo'] = stereo.group() == '@'

    # supported only <;> and <,> logic. <&> and silent <&> not supported!
    primitives = token.split(';')
    if element := primitives[0]:
        element = [int(x[1:]) if x.startswith('#') else x for x in element.split(',')]
        if len(element) == 1:
            element = element[0]
    else:
        raise IncorrectSmarts('Empty element')
    out['element'] = element

    for p in primitives[1:]:  # parse hydrogens (h), neighbors (D), rings_sizes (r or !R), hybridization == 4 (a)
        if not p:
            continue
        elif p == 'a':  # aromatic atom
            out['hybridization'] = 4
        elif p == 'A':  # ignore aliphatic mark. Ad-Hoc for Marwin.
            continue
        elif p == '!R':
            out['ring_sizes'] = 0
        elif p == 'M':
            out['masked'] = True
        else:
            p = p.split(',')
            if len(p) != 1 and len({x[0] for x in p}) > 1:
                raise IncorrectSmarts('Unsupported OR statement')
            elif (t := p[0][0]) not in ('D', 'h', 'r', 'x', 'z'):
                raise IncorrectSmarts('Unsupported SMARTS primitive. Use only D, h, r, !R and a.')
                # z and x private chython marks for hybridization and heteroatoms count
            try:
                p = [int(x[1:]) for x in p]
            except ValueError:
                raise IncorrectSmarts('Unsupported SMARTS primitive')

            if t == 'D':
                out['neighbors'] = p
            elif t == 'h':
                out['implicit_hydrogens'] = p
            elif t == 'r':  # r
                out['ring_sizes'] = p
            elif t == 'x':
                out['heteroatoms'] = p
            else:  # z
                out['hybridization'] = p
    return 0, out


def smiles_tokenize(smi):
    tokens = _tokenize(smi)
    out = []
    for token_type, token in tokens:
        if token_type in (0, 8):  # simple atom
            out.append((token_type, {'element': token}))
        elif token_type == 5:
            out.append(_atom_parse(token))
        elif token_type == 10:
            raise IncorrectSmiles('SMARTS detected')
        else:
            out.append((token_type, token))
    return out


def smarts_tokenize(smi):
    tokens = _tokenize(smi)
    out = []
    for token_type, token in tokens:
        if token_type in (0, 8):  # simple atom
            out.append((0, {'element': token}))
        elif token_type == 5:
            out.append(_query_parse(token))
        else:
            out.append((token_type, token))
    return out


__all__ = ['smiles_tokenize', 'smarts_tokenize']
