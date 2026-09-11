# -*- coding: utf-8 -*-
#
#  Copyright 2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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
"""Read a knowledge table's SMARTS with the core lexer, and number its atoms the tables' way.

Parsing is `read_smarts`; all this adds is the numbering a patch column is keyed by -- an atom's
explicit `:N`, else the lowest unclaimed positive integer in declaration order.  The core's own answer
is a stable id, so the translation lives here beside the tables it serves.  A table written for this
numbering must spell `[M]` as `[M;*;^,!^]` if it means any charge, since `[M]` alone is neutral.
"""
from ..core import IncorrectSmarts, QueryContainer, read_smarts


__all__ = ['compile_smarts', 'SmartsSyntaxError']


# the core lexer's own refusal, under the name this package catches.  An alias, not a wrapper class:
# re-raising would put this module in the traceback of a position reported by the lexer.
SmartsSyntaxError = IncorrectSmarts


def compile_smarts(text: str) -> tuple[QueryContainer, dict[int, int], dict[int, str]]:
    """Compile `text` and return `(query, numbers, wildcards)`, both dicts keyed by atom number.

    `numbers` maps a table atom number to the core stable id.  `wildcards` names the atoms that
    constrain no element, as `{number: 'any' | 'metal'}`; a rule collection reads it to tell shared
    context from the site being repaired.

    Raises `SmartsSyntaxError` (the core's `IncorrectSmarts`) at load time, with the offending
    position in the message: a typo in a knowledge file must fail where the typo is.
    """
    query = read_smarts(text)
    numbers = _number_atoms(query)
    wildcards = query.wildcard_atoms()
    return query, numbers, {number: wildcards[sid] for number, sid in numbers.items()
                            if sid in wildcards}


def _number_atoms(query: QueryContainer) -> dict[int, int]:
    """The tables' numbering: explicit map number, else the lowest unclaimed, in declaration order.

    Declaration order is `range(1, atom_count + 1)`, since ids are allocated from 1 as the lexer reads
    left to right.  Not `query_numbers()`, which is the sealed order: that DFS roots at the rarest
    element, so `[C][O]` seals oxygen first.
    """
    explicit = query.map_numbers()            # {stable id: map number}, non-zero only
    claimed = set(explicit.values())
    if len(claimed) != len(explicit):
        raise SmartsSyntaxError('two atoms carry the same map number')
    numbers: dict[int, int] = {}
    nxt = 1
    for sid in range(1, query.atom_count + 1):
        if sid in explicit:
            numbers[explicit[sid]] = sid
        else:
            while nxt in claimed:
                nxt += 1
            numbers[nxt] = sid
            nxt += 1
    return numbers
