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
"""Nucleic-acid segmentation: nucleotides split at their phosphorus links, the main strand in 5' to 3' order.

Read-only topology beside `monomers.py`; `chains.perceive` asks it after peptides and before macrocycles.

| Step | Rule |
| --- | --- |
| sugar | a five-membered ring of one oxygen, O4': C4' beside it bears an exocyclic C5' bearing an O5', the ring \
carbon C3' beside C4' bears an exocyclic O3'; a 2' substituent, a base or none at all are the unit's groups; a C5' \
with no O5', a 5'-vinylphosphonate's, opens the strand only |
| link | O3' of one unit and O5' of another bonded to one phosphorus: phosphodiester, phosphorothioate, any P(V) |
| unit path | P, O5', C5', C4', C3', O3', the 5'-phosphorus first, a 5' cap's too; a unit with no 5' phosphorus \
starts at O5', or C5' |
| main strand | the longest run of links, 5' to 3'; fewer than `MIN_UNITS` is no strand |
| groups | every other atom, claimed in rank order from the units without crossing the backbone: bases, 2' groups, \
the phosphorus' oxygens and sulfur, the conjugates at either end |
"""
from .chains import Chain, Residue, _claim


__all__ = ['nucleotides', 'MIN_UNITS']

MIN_UNITS = 4


def nucleotides(mol):
    """`mol`'s main strand as a `Chain` of kind nucleic, or None when it has fewer than `MIN_UNITS` linked units."""
    if len(mol.connected_components) != 1:
        return None
    rank = mol.canonical_order()
    units = [u for r in mol.sssr if len(r) == 5 for u in _sugars(mol, r)]
    if len(units) < MIN_UNITS:
        return None
    by_o5 = {u['o5']: u for u in units if u['o5'] is not None}
    after, before = {}, {}
    for u in units:
        for p in mol.neighbors_of(u['o3']):
            if mol.element_of(p) != 15:
                continue
            for o in mol.neighbors_of(p):
                v = by_o5.get(o)
                if v is not None and v is not u and id(v) not in before and id(u) not in after:
                    after[id(u)], before[id(v)] = (p, v), u
    best = []
    for u in sorted(units, key=lambda u: rank[u['c4']]):
        if id(u) in before:
            continue
        run = [(None, u)]
        while id(run[-1][1]) in after:
            run.append(after[id(run[-1][1])])
        if len(run) > len(best):
            best = run
    if len(best) < MIN_UNITS:
        return None
    residues, backbone = [], []
    head = best[0][1]
    if head['o5'] is not None:                          # a 5' cap on a phosphorus, an inverted abasic's: it opens
        caps = [y for y in mol.neighbors_of(head['o5']) if mol.element_of(y) == 15]
        if caps:
            best[0] = (caps[0], head)
    for p, u in best:
        path = ([p] if p is not None else []) + ([u['o5']] if u['o5'] is not None else []) + \
            [u['c5'], u['c4'], u['c3'], u['o3']]
        backbone += path
        residues.append(Residue(frozenset(path) | u['ring'], 'nucleotide', path[0], path[-1], tuple(path)))
    members = _claim(mol, residues, set(backbone), rank)
    for m, r in zip(members, residues):
        m |= r.atoms
    return Chain(tuple(residues), tuple(frozenset(m) for m in members), tuple(range(len(residues))), (), (), False,
                 tuple(backbone), (), 'nucleic')


def _sugars(mol, ring):
    """The furanose readings of `ring`: `{o4, c4, c5, o5, c3, o3, ring}` for each C4' beside its O4'."""
    oxygens = [x for x in ring if mol.element_of(x) == 8]
    if len(oxygens) != 1 or any(mol.element_of(x) != 6 for x in ring if x != oxygens[0]):
        return []
    o4, out = oxygens[0], []
    for c4 in (y for y in mol.neighbors_of(o4) if y in ring):
        c3 = next(y for y in mol.neighbors_of(c4) if y in ring and y != o4)
        for c5 in mol.neighbors_of(c4):
            if c5 in ring or mol.element_of(c5) != 6:
                continue
            o5 = [y for y in mol.neighbors_of(c5) if y != c4 and mol.element_of(y) == 8]
            o3 = [y for y in mol.neighbors_of(c3) if y not in ring and mol.element_of(y) == 8]
            if len(o5) <= 1 and len(o3) == 1:               # no O5': a 5' head alone, a vinylphosphonate's
                out.append({'o4': o4, 'c4': c4, 'c5': c5, 'o5': o5[0] if o5 else None, 'c3': c3, 'o3': o3[0],
                            'ring': frozenset(ring)})
    return out
