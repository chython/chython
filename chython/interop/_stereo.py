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
"""
The stereo translation every converter needs, written once.

A unit's `refs` is four slots: an atom kind (`SU_TETRA`) packs four directions in one order, a bond
kind packs two direction lists of two slots each, anchor's pair first, and an unnamed slot may sit in
the *middle* of a half rather than only at the tail -- hence the helpers here instead of `refs[2]`.
Nothing here normalises: a configuration the foreign tool has no form for is reported, never absorbed.
"""
from collections.abc import Sequence
from typing import NamedTuple
from ..core import SU_ALLENE, SU_ATROPISOMER, SU_CIS_TRANS, SU_HELICAL, SU_TETRA


#: Every stereo-unit kind the core can emit, with a name fit for a log line.  Keyed off
#: `chython.core`'s own constants, so a kind the core grows cannot fall off a converter's chain
#: unreported.
KIND_NAMES = {
    SU_TETRA: 'tetrahedral',
    SU_CIS_TRANS: 'cis/trans',
    SU_ALLENE: 'allene/cumulene',
    SU_ATROPISOMER: 'atropisomer',
    SU_HELICAL: 'helical',
}


def kind_name(kind: int, /) -> str:
    """A name for a stereo-unit kind, including one the core grew after this table was written."""
    return KIND_NAMES.get(kind, f'stereo kind {kind}')


def first_of_each_pair(refs: Sequence[int | None], /) -> tuple[int | None, int | None]:
    """
    The first named atom of each half of a bond unit's four direction slots.

    A fallback and not `refs[0], refs[2]`, because `_stereo.pxi` allows a bond kind to carry an unnamed
    slot *inside* a half and does not promise the halves ascend across each other.

    Returns `(None, None)` for a half that names nothing; the caller must report that as "no frame"
    rather than pass it on to `translate_stereo`.
    """
    near = refs[0] if refs[0] is not None else refs[1]
    far = refs[2] if refs[2] is not None else refs[3]
    return near, far


def cis_trans_partner(mol, anchor: int, refs: Sequence[int | None], /) -> int | None:
    """
    The other end of the bond a cis/trans-like unit is anchored on.

    Found from the refs through the graph -- the partner is the anchor's neighbour that the *far* half's
    directions hang off -- and never as "the neighbour with bond order 2", which is not part of the
    definition.  Returns `None` when no neighbour carries the far reference, which the caller reports.
    """
    far = refs[2] if refs[2] is not None else refs[3]
    if far is None:
        return None
    for p in mol.neighbors_of(anchor):
        if far in mol.neighbors_of(p):
            return p
    return None


class CisTransFrame(NamedTuple):
    """
    A cis/trans unit resolved into the four things a converter needs to name a configuration.

    `near`/`far` are the reference substituents the descriptor is about, on the anchor's and the
    partner's side.  `order` arranges those four directions so `translate_stereo(anchor, order)` accepts
    them, turning the stored parity into a parity in the tool's own frame.
    """
    near: int
    far: int
    partner: int
    order: tuple[int | None, int | None, int | None, int | None]


def cis_trans_frame(mol, anchor: int, refs: Sequence[int | None], /) -> CisTransFrame | None:
    """
    Resolve a cis/trans-like unit into `(near, far, partner, order)`, or `None` when it cannot be framed.

    One function and not three calls, because the order must be built from the same `near`/`far` the
    caller reports to the tool; recomputing one independently names one pair of atoms while the parity
    describes another.
    """
    near, far = first_of_each_pair(refs)
    if near is None or far is None:
        return None
    partner = cis_trans_partner(mol, anchor, refs)
    if partner is None:
        return None
    order = (near, refs[1] if refs[0] == near else refs[0],
             far, refs[3] if refs[2] == far else refs[2])
    return CisTransFrame(near, far, partner, order)


def set_parity_by_probe(mol, anchor: int, order: Sequence[int | None], want: int, /) -> bool:
    """
    Write the stored parity that reads back as `want` when read in `order`.  True when it was written.

    By probing, not by arithmetic: `set_parity` writes in the unit's own refs order while
    `translate_stereo` reads in the caller's, and restating the permutation-parity rule here is the one
    place a sign error yields a silently mirrored molecule rather than a raise.  So write parity 1, ask
    `translate_stereo` what it reads as, flip if wrong.

    On a frame the core rejects the parity is *cleared* and `False` returned -- half a configuration is
    indistinguishable from a real one, so it is worse than none.
    """
    try:
        mol.set_parity(anchor, 1)
        if mol.translate_stereo(anchor, tuple(order)) != want:
            mol.set_parity(anchor, 2)
    except (KeyError, ValueError):
        try:
            mol.set_parity(anchor, 0)
        except (KeyError, ValueError):  # pragma: no cover -- a bad anchor cannot be cleared either
            pass
        return False
    return True


__all__ = ['KIND_NAMES', 'CisTransFrame', 'kind_name', 'first_of_each_pair', 'cis_trans_partner',
           'cis_trans_frame', 'set_parity_by_probe']
