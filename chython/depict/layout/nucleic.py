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
"""A nucleic-acid strand laid out level, 5' on the left: every nucleotide the same drawing, its base up.

| part | rule |
| --- | --- |
| unit | `UNIT`: the furanose a regular pentagon on its bottom edge C4'-C3', C1' at its top, O4' up-left and C2' \
up-right; O3' on C3''s exterior bisector; the phosphorus a square cross, its O3'-O5' axis along `LINK`, its =O \
up and its other atom down |
| strand | unit `k` the unit drawing shifted `k * PERIOD` along x: its 5'-phosphorus at the height of the last \
unit's, one bond from that unit's O3', so the strand never drifts and no link is long |
| groups | bases, 2' groups, what hangs from the phosphorus' two atoms and the conjugates at either end hung from \
the drawn atoms by the chain layout's group placement, `peptide._place_groups` |
| rest | a part left with two attachments or none, a bridged sugar's ring among them, refuses: None |
"""
from math import asin, cos, pi, radians, sin, tan

from .peptide import L, _place_groups, _ring_systems, _Scene, _Tiles


__all__ = ['nucleic_layout', 'UNIT', 'PERIOD']

LINK = radians(18)      # the O3'-P-O5' axis of each phosphorus, its other two atoms square to it, =O up
SWING = radians(-30)    # O5' to C5'


def _unit():
    """The unit drawing in bonds: the pentagon on C4' (0, 0) and C3' (1, 0), O3' on C3''s exterior bisector, the link
    from it along `LINK` through the phosphorus to O5', O5' to C5' along `SWING`, and C5' to the next unit's C4' the
    turn bringing that C4' level with this one; the phosphorus, O5' and C5' are the next unit's, shifted back."""
    out = {'c4': (0., 0.), 'c3': (1., 0.), 'c2': (1 + cos(radians(72)), sin(radians(72))),
           'c1': (.5, .5 * tan(radians(72))), 'o4': (-cos(radians(72)), sin(radians(72)))}
    o3 = (1 + cos(radians(-54)), sin(radians(-54)))
    p = (o3[0] + cos(LINK), o3[1] + sin(LINK))
    o5 = (p[0] + cos(LINK), p[1] + sin(LINK))
    c5 = (o5[0] + cos(SWING), o5[1] + sin(SWING))
    turn = asin(-c5[1])
    period = c5[0] + cos(turn)
    out.update(o3=o3, p=(p[0] - period, p[1]), o5=(o5[0] - period, o5[1]), c5=(c5[0] - period, c5[1]))
    return out, period


UNIT, PERIOD = _unit()  # unit path role -> its point in bonds; from one unit's C4' to the next's


class _Context:
    """What `_place_groups` reads of a chain layout's context."""
    __slots__ = ('peptide', 'rank', 'rings', 'tiles', 'cis_trans')

    def __init__(self, mol, strand, tile):
        self.peptide = strand
        self.rank = mol.canonical_order()
        self.rings = _ring_systems(mol, strand.members)
        self.tiles = _Tiles(mol, tile)
        self.cis_trans = [u for u in mol.stereo_units() if u['kind'] == 1 and u['parity']]


def nucleic_layout(mol, strand, tile):
    """`{atom: (x, y)}` for every atom of `mol` with the units of `strand`, `nucleotides(mol)`, drawn level, or None.

    :param tile: `tile(sub) -> {atom: (x, y)}` at mean bond `L`, for the groups' ring systems.
    """
    sc = _Scene(mol, frozenset())
    core, rank = {}, mol.canonical_order()
    for k, r in enumerate(strand.residues):
        path = (None,) * (6 - len(r.path)) + r.path     # a 5' end's path is the tail of a unit's
        _, _, _, c4, c3, o3 = path
        o4 = next(y for y in mol.neighbors_of(c4) if y in r.atoms and y not in path)
        c2 = next(y for y in mol.neighbors_of(c3) if y in r.atoms and y not in path)
        c1 = next(y for y in mol.neighbors_of(o4) if y in r.atoms and y != c4)
        dx = k * PERIOD
        for role, x in zip(('p', 'o5', 'c5', 'c4', 'c3', 'o3', 'o4', 'c1', 'c2'), path + (o4, c1, c2)):
            if x is not None:
                u, v = UNIT[role]
                core[x] = (L * (dx + u), L * v)
        tail = [y for y in mol.neighbors_of(o3) if mol.element_of(y) == 15 and y not in core] \
            if k == len(strand.residues) - 1 else []
        if tail:                                        # a 3' phosphorus: where the next unit's would be
            u, v = UNIT['p']
            core[tail[0]] = (L * (dx + PERIOD + u), L * v)
        for p in [x for x in (path[0], *tail) if x is not None]:    # each phosphorus a square cross, =O up
            rest = sorted((y for y in mol.neighbors_of(p) if y not in r.atoms and y not in core),
                          key=lambda y: (-mol.order_of(p, y), rank[y]))
            px, py = core[p]
            bridge = [y for y in rest if len(mol.neighbors_of(y)) > 1 and mol.order_of(p, y) == 1]
            if len(rest) == 3 and len(bridge) == 1:     # a cap's or a conjugate's oxygen: on the axis
                t = LINK if p in tail else LINK + pi
                core[bridge[0]] = (px + L * cos(t), py + L * sin(t))
                rest.remove(bridge[0])
            if len(rest) == 2:
                for y, t in zip(rest, (LINK + pi / 2, LINK - pi / 2)):
                    core[y] = (px + L * cos(t), py + L * sin(t))
    sc.commit(core)
    ctx = _Context(mol, strand, tile)
    order = list(range(len(strand.residues)))
    for _ in range(6):
        _, progress = _place_groups(sc, ctx, order)
        if not progress:
            break
    if len(sc.pos) != len(mol):
        return None
    return dict(sc.pos)
