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
"""Approved drawings: the user's verdicts, the evaluation set the layout engine is held to; the engine never reads them.

| part | rule |
| --- | --- |
| library | `LIBRARY`: one row per approved molecule, its id, SMILES and plane, in its SMILES atom order |
| entry | each component of two or more atoms, its plane scaled to a mean bond of `L` and centred |
| key | a component's skeleton, `shapes.skeleton`: elements and bond orders ignored, so a full analog matches |
| held | `same_picture`: the engine draws an approved molecule as approved, a mirror image counted the same; \
`held`: it draws each component of any molecule as an approved drawing of that component's skeleton, so a full \
analog is held to its approved twin |
| writes | `approve` and `withdraw`, from the review page's likes and dislikes |
"""
import json
from math import fsum, hypot
from pathlib import Path

from .lattice import L
from .shapes import skeleton


__all__ = ['LIBRARY', 'approve', 'withdraw', 'read_library', 'picture', 'same_picture', 'held']

LIBRARY = Path(__file__).resolve().parents[3] / 'test' / 'layout_approved.tsv'
SAME = 1e-3             # `L`s two points of one picture may differ by
_INDEX = {}             # (path, mtime) -> {skeleton bytes: [(molecule, component atoms, plane)]}


def read_library(path=LIBRARY):
    """`{id: (smiles, plane)}` of the approved drawings."""
    path = Path(path)
    if not path.exists():
        return {}
    out = {}
    for line in path.read_text().splitlines()[1:]:
        if line:
            i, text, plane = line.split('\t')
            out[i] = (text, json.loads(plane))
    return out


def _write(path, library):
    tmp = Path(f'{path}.tmp')
    tmp.write_text('id\tsmiles\tplane\n' + ''.join(f'{i}\t{t}\t{json.dumps(p)}\n'
                                                  for i, (t, p) in sorted(library.items())))
    tmp.replace(path)


def approve(row_id, text, plane, path=LIBRARY):
    """Keep `plane`, a point list in the atom order of `text`, as the approved drawing of `row_id`."""
    library = read_library(path)
    library[row_id] = (text, [[round(v, 6) for v in p] for p in plane])
    _write(path, library)


def withdraw(row_id, path=LIBRARY):
    """Drop the approval of `row_id`, if any; True when there was one."""
    library = read_library(path)
    if library.pop(row_id, None) is None:
        return False
    _write(path, library)
    return True


def _normal(mol, atoms, plane):
    bonds = [hypot(plane[n][0] - plane[m][0], plane[n][1] - plane[m][1]) for n in atoms
             for m in mol.neighbors_of(n) if m in atoms and n < m]
    if not bonds:
        return None
    k = L / (fsum(bonds) / len(bonds))
    cx = fsum(plane[x][0] for x in atoms) / len(atoms)
    cy = fsum(plane[x][1] for x in atoms) / len(atoms)
    return {x: ((plane[x][0] - cx) * k, (plane[x][1] - cy) * k) for x in atoms}


def _index(path=LIBRARY):
    path = Path(path)
    stamp = (str(path), path.stat().st_mtime_ns if path.exists() else None)
    if stamp not in _INDEX:
        from ...core import smiles
        index = {}
        for text, points in read_library(path).values():
            mol = smiles(text)
            plane = dict(zip(mol, map(tuple, points)))
            for c in mol.connected_components:
                atoms = set(c)
                if len(atoms) < 2:
                    continue
                pts = _normal(mol, atoms, plane)
                if pts is not None:
                    index.setdefault(skeleton(mol, atoms)[0].canonical_bytes, []).append((mol, atoms, pts))
        _INDEX.clear()
        _INDEX[stamp] = index
    return _INDEX[stamp]


def picture(mol, plane):
    """The drawing of each component of two or more atoms as its points with their elements, centred and scaled to
    a mean bond of `L`: equal for two drawings showing the same picture, whatever atom a symmetric twin stands for."""
    out = []
    for c in sorted((set(c) for c in mol.connected_components), key=len):
        pts = _normal(mol, c, plane)
        if pts is not None:
            out.append(sorted((mol.atom(x).atomic_symbol, round(p[0] / L, 3) + 0., round(p[1] / L, 3) + 0.)
                              for x, p in pts.items()))
    return out


MIRRORS = ((1, 1), (-1, 1), (1, -1), (-1, -1))      # a component as drawn, mirrored left-right or top-bottom, or turned


def same_picture(mol, plane, other):
    """True when `plane` and `other` draw every component of `mol` as the same picture, up to `MIRRORS`: each point of
    one, centred and scaled to a mean bond of `L`, within `SAME` of a point of the same element in the other."""
    for c in (set(c) for c in mol.connected_components):
        a, b = _normal(mol, c, plane), _normal(mol, c, other)
        if a is None:
            continue
        if not any(_matched(mol, c, {x: (p[0] * sx, p[1] * sy) for x, p in a.items()}, b) for sx, sy in MIRRORS):
            return False
    return True


def _matched(mol, atoms, a, b):
    free = list(atoms)
    for x in atoms:
        e = mol.atom(x).atomic_symbol
        for i, y in enumerate(free):
            if mol.atom(y).atomic_symbol == e and hypot(a[x][0] - b[y][0], a[x][1] - b[y][1]) <= SAME * L:
                del free[i]
                break
        else:
            return False
    return True


def held(mol, plane, path=LIBRARY):
    """None when some component of two or more atoms of `mol` has no approved skeleton, else True when `plane` draws
    each such component as an approved drawing of its skeleton: the same points within `SAME`, up to `MIRRORS`,
    elements aside."""
    index = _index(path)
    for c in (set(c) for c in mol.connected_components):
        if len(c) < 2:
            continue
        found = index.get(skeleton(mol, c)[0].canonical_bytes)
        if not found:
            return None
        mine = _normal(mol, c, plane)
        if not any(_same([(p[0] * sx, p[1] * sy) for p in mine.values()], list(pts.values()))
                   for _, _, pts in found for sx, sy in MIRRORS):
            return False
    return True


def _same(a, b):
    free = list(b)
    for p in a:
        for i, q in enumerate(free):
            if hypot(p[0] - q[0], p[1] - q[1]) <= SAME * L:
                del free[i]
                break
        else:
            return False
    return True
