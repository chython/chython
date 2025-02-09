# -*- coding: utf-8 -*-
#
#  Copyright 2022, 2023 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from ...exceptions import IncorrectSmiles


def parser(tokens, strong_cycle):
    log = []
    t1 = tokens[0][0]
    if t1 == 2:
        if tokens[1][0] not in (0, 8):
            raise IncorrectSmiles('not atom started')
    elif t1 not in (0, 8):
        raise IncorrectSmiles('not atom started')

    atoms = []
    bonds = []
    order = defaultdict(list)
    atoms_types = []
    atom_num = 0
    last_num = 0
    stack = []
    cycles = {}
    stereo_atoms = {}
    stereo_bonds = defaultdict(dict)
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
        elif token_type in (1, 4, 9, 10, 12):  # bonds. only keeping for atoms connecting
            if previous:
                raise IncorrectSmiles('2 bonds in a row')
            elif not atoms:
                raise IncorrectSmiles('started from bond')
            previous = (token_type, token)
        elif token_type == 6:  # cycle
            if previous and previous[0] == 4:
                raise IncorrectSmiles('dot-cycle pattern invalid')
            elif token not in cycles:
                cycles[token] = (last_num, previous, len(order[last_num]))
                order[last_num].append(None)  # Reserve a table
                previous = None
            else:
                a, ob, ind = cycles[token]
                if ob:
                    if not previous:
                        bt, b = ob
                        if bt == 9:  # closure open is \/ bonded
                            stereo_bonds[a][last_num] = b
                            b = 1
                        elif strong_cycle:
                            raise IncorrectSmiles('not equal cycle bonds')
                        else:
                            log.append('ignored difference in cycle bonds')
                    else:
                        bt, b = previous
                        obt, ob = ob
                        if bt == 9:  # \/ bonds can be unequal
                            if obt == 9:
                                stereo_bonds[a][last_num] = ob
                            elif ob != 1:
                                raise IncorrectSmiles('not equal cycle bonds')
                            stereo_bonds[last_num][a] = b
                            b = 1
                        elif obt == 9:
                            if b != 1:
                                raise IncorrectSmiles('not equal cycle bonds')
                            stereo_bonds[a][last_num] = ob
                        elif b != ob:
                            raise IncorrectSmiles('not equal cycle bonds')
                        previous = None
                elif previous:
                    bt, b = previous
                    if bt == 9:  # stereo \/
                        stereo_bonds[last_num][a] = b
                        b = 1
                    elif strong_cycle:
                        raise IncorrectSmiles('not equal cycle bonds')
                    else:
                        log.append('ignored difference in cycle bonds')
                    previous = None
                else:
                    b = 4 if atoms_types[last_num] == atoms_types[a] == 8 else 1

                bonds.append((last_num, a, b))
                order[a][ind] = last_num
                order[last_num].append(a)
                del cycles[token]
        else:  # atom
            if atoms:
                if not previous:
                    bonds.append((atom_num, last_num, 4 if token_type == atoms_types[last_num] == 8 else 1))
                    order[last_num].append(atom_num)
                    order[atom_num].append(last_num)
                else:
                    bt, b = previous
                    if bt == 9:
                        bonds.append((atom_num, last_num, 4 if token_type == atoms_types[last_num] == 8 else 1))
                        order[last_num].append(atom_num)
                        order[atom_num].append(last_num)

                        stereo_bonds[last_num][atom_num] = b
                        stereo_bonds[atom_num][last_num] = not b
                    elif bt in (1, 10, 12):
                        bonds.append((atom_num, last_num, b))
                        order[last_num].append(atom_num)
                        order[atom_num].append(last_num)
                    # else bt == 4 - skip dot
                    previous = None

            if token.get('stereo') is not None:
                stereo_atoms[atom_num] = token.pop('stereo')
            atoms.append(token)
            atoms_types.append(token_type)
            last_num = atom_num
            atom_num += 1

    if stack:
        raise IncorrectSmiles('number of ( does not equal to number of )')
    elif cycles:
        raise IncorrectSmiles('cycle is not finished')
    elif previous:
        raise IncorrectSmiles('bond on the end')

    return {'atoms': atoms, 'bonds': bonds, 'order': order, 'stereo_atoms': stereo_atoms,
            'stereo_bonds': stereo_bonds, 'log': log}


__all__ = ['parser']
