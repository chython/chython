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
from functools import partial
from itertools import count
from re import compile, findall, search
from .parser import parser
from .tokenize import smarts_tokenize
from ...containers import QueryContainer, QueryBond
from ...periodictable import ListElement, QueryElement


cx_radicals = compile(r'\^[1-7]:[0-9]+(?:,[0-9]+)*')

# AD-HOC for masked atoms in SMARTS
# not thread safe
global_free_masked = count(10 ** 9 + 1)


def smarts(data: str):
    """
    Parse SMARTS string.

    * only D, a, h, r and !R atom primitives supported.
    * bond order list (max 2) and not bond supported.
    * [not]ring bond supported only in combination with explicit bonds, not bonds and bonds orders lists.
    * mapping, charge and isotopes supported.
    * list of elements supported.
    * A - treats as any element. <A> primitive (aliphatic) ignored.
    * M - treats as any metal
    * <&> logic operator unsupported.
    * <;> logic operator is mandatory except (however preferable) for charge, isotope, stereo marks.
    * CXSMARTS radicals supported.
    * masked atom - `chython.Reactor` specific mark for masking reactant atoms from deletion.

    For example::

        [C;r5,r6;a]-;!@[C;h1,h2;z2,z4] |^1:1| - aromatic C member of 5 or 6 atoms ring
        connected with non-ring single bond to aromatic or SP2 radical C with 1 or 2 hydrogens.

    * primitive <xN> - heteroatoms (e.g. x2 - two heteroatoms)
    * primitive <zN> - hybridization (N = 1 - sp3, 2 - sp2, 3 - sp, 4 - aromatic)
    * primitive <M> - masked atom

    Note: atom numbers greater than 10 ** 9 forbidden for usage and reserved for masked atoms numbering.
    In multiprocess mode has potential bugs in reaction enumeration task then used templates prepared from components
    from different processes. For avoiding, prepare templates on single process and then share it.
    """
    if not isinstance(data, str):
        raise TypeError('Must be a SMARTS string')
    smr, *cx = data.split()

    parsed = parser(smarts_tokenize(smr), False)

    if cx and cx[0].startswith('|') and cx[0].endswith('|'):
        for x in findall(cx_radicals, cx[0]):
            for i in x[3:].split(','):
                parsed['atoms'][int(i)]['is_radical'] = True

    for i, s in parsed['stereo_atoms'].items():
        parsed['atoms'][i]['stereo'] = s
    stereo_bonds = parsed['stereo_bonds']

    g = QueryContainer(data)

    mapping = {}
    free = count(max(a.get('parsed_mapping', 0) for a in parsed['atoms']) + 1)
    for i, a in enumerate(parsed['atoms']):
        mapping[i] = n = a.pop('parsed_mapping', 0) or next(global_free_masked if a.get('masked') else free)
        e = a.pop('element')
        if isinstance(e, int):
            e = QueryElement.from_atomic_number(e)
        elif isinstance(e, str):
            e = QueryElement.from_symbol(e)
        else:
            e = partial(ListElement, e)
        g.add_atom(e(**a), n)

    for n, m, b in parsed['bonds']:
        if n in stereo_bonds and m in stereo_bonds:
            if m not in stereo_bonds[n]:  # only simple cis-trans supported, not cumulenes
                _, s1 = stereo_bonds[n].popitem()
                _, s2 = stereo_bonds[m].popitem()
                if isinstance(b, int):
                    b = QueryBond(b, stereo=s1 == s2)
                else:
                    b.stereo = s1 == s2
        g.add_bond(mapping[n], mapping[m], b)
    return g


__all__ = ['smarts']
