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
"""The primitive environment: the honeycomb tree embedding with ring-system templates, then flips, rotations, stretches
and whole-drawing turns, nothing else. Macrocycles no template covers, substituents and orientation are what a policy
learns.

| part | rule |
| --- | --- |
| start | `tree.tree_embedding`, then each templated ring system drawn: fewest broken cis/trans units, then \
defects, then the order of its drawings; built once per molecule component by `Initial`, copied for every rollout |
| templates | each ring system `cascade.system_drawing` draws: a shape, a regular polygon, a cached system drawing \
with its rings evened by `relax`'s ring steps, or a macrocycle of small-ring units; the drawing or its mirror image, \
whichever breaks fewer cis/trans units, then adds fewer defects; any other macrocycle stays on the tree |
| entry | a template hanging from the tree turned, with its subtree, about its entry atom so that the bond to the \
tree leaves along the entry atom's first exit direction (a shape: away from its template neighbours); at a spiro \
atom the two systems back to back instead; parents first |
| exits | substituent out: a polygon template atom's substituents, each with its subtree, on the directions splitting \
evenly its widest gap outside every ring corner, a fusion atom's lone one on the extension of its fused bond; a \
shape's, or an atom with neither, on whichever of `EXITS` directions leaves the subtree most room up to `ROOM`, then \
is nearest the direction away from the template; two of one atom at least `APART`; a substituent-out exit leaving \
less than `CLEAR` gives way to the most room when that has fewer defects |
| boundary | every tree bond outside a template then set to `L`, its subtree moved along it, parents first |
| upright | the whole drawing then turned so the largest `shapes.shape` system stands as its shape is drawn |
| template | a drawn system moves only rigidly: no flip, rotation or stretch is offered on a bond inside it; the \
crossings a bridged template draws are `inherent`, not `defects` |
| `flip` | a tree bond's subtree mirrored across the bond |
| `rotate` | a tree bond's subtree turned rigidly about its parent atom until the angle anchor-parent-child is \
one of `ANGLES`, on either side of the anchor bond; the anchor is the parent's own tree parent, or at the root its \
first other child; an angle already drawn is not offered |
| `stretch` | a tree bond's subtree moved along the bond by `STEP` either way, the bond kept within `L` and `LONGEST` |
| `turn` | the whole drawing turned about its centroid by each of `SCENE` |
| `stop` | ends the drawing; offered once no defect is left, or when no other step is |
| masks | `edits.masked`: no stated cis/trans broken, no linear atom bent, no bond leaving `L` to `LONGEST`; a step \
may pass through a defect, `stop` may not |

A macrocycle no template covers stays on the tree: its ring bonds are tree bonds like any other, and rotating or
stretching one moves its closure, which the score prices.
"""
from math import atan2, cos, pi, radians, sin

from .cascade import TEMPLATE, in_ring_angle, system_drawing
from .edits import STEP, Move, broken, embed, masked, outside_angles
from .env import STOP, Step
from .lattice import L
from .relax import _ct_broken, fair, regularize, seat, wheels
from .rings import ring_blocks
from .score import TERMS, Scope, _parts, add, combine, delta_terms, key, terms
from .tree import Grid, tree_embedding


__all__ = ['KINDS', 'RINGS', 'ANGLES', 'SCENE', 'CAP', 'LIMIT', 'LONGEST', 'TEMPLATE', 'Initial', 'PrimitiveEnv']

KINDS = {'flip': 0, 'rotate': 1, 'stretch': 2, 'turn': 3, 'stop': 4}
RINGS = range(3, 10)    # ring sizes whose regular angle a rotation sets
ANGLES = tuple(radians(180 - 360 / n) for n in RINGS) + (pi,)
SCENE = tuple(radians(a) for a in (30, -30, 90, -90))
CAP = 4                 # steps allowed per atom
LIMIT = 64              # and in all
LONGEST = 6 * L         # a stretched tree bond
EXITS = 24              # directions tried for a template atom's substituent
APART = radians(50)     # least angle between two substituents of one template atom
ROOM = 1.5 * L          # clearance beyond which a direction is as good as any
AWAY = .8 * L           # a shape's: past it, the direction straight out of the shape wins
CLEAR = .3 * L          # clearance below which a substituent-out exit gives way to the most room
STRETCH = .25 * L       # the step a blocked substituent-out exit stretches by, up to `L` more
SMALL = 4               # the most atoms a stretched substituent holds: a nitro, a carboxyl, a sulfonyl


class Initial:
    """A component's start, built once: its scope, tree embedding with templates, and every tree bond's subtree."""
    __slots__ = ('mol', 'scope', 'plane', 'root', 'edges', 'anchor', 'rigid', 'shaped', 'inherent', 'upright', 'terms',
                 'ct')

    def __init__(self, mol, atoms, rank=None):
        rank = mol.canonical_order() if rank is None else rank
        tree = tree_embedding(Scope(mol, atoms, rank, []))
        self.mol = mol
        self.scope = Scope(mol, atoms, rank, tree.closures)
        self.scope.orient = True                        # the drawing's orientation is a step here, not a post-step
        self.plane = dict(tree.plane)
        self.rigid = []                                 # the atoms of each templated ring system
        upright, even, kinds = None, [], []
        for b, rs in ring_blocks(mol):                  # canonical order: each system on the tree drawn so far
            if not b <= atoms:
                continue
            found = system_drawing(mol, b, rs, rank)
            if found is None:
                continue
            kind, drawings = found
            best = None
            for k, v in enumerate(drawings):            # the drawing or its mirror: cis/trans, then defects
                moved = embed(self.scope, self.plane, b, v, tree.root)
                score = (broken(self.scope, self.plane, moved), delta_terms(self.scope, self.plane, moved)[0], k)
                if best is None or score < best[0]:
                    best = (score, moved, v)
            self.plane.update(best[1])
            self.rigid.append(frozenset(b))
            kinds.append(kind)
            if kind == 'shape' and (upright is None or len(b) > len(upright[0])):
                upright = (b, best[2])
            elif kind == 'system':
                even += rs
        for fix in (regularize, wheels, fair, seat):    # a system drawing's rings, as `relax` evens them
            fix(mol, self.plane, even, rank)
        self.shaped = frozenset().union(*(b for b, k in zip(self.rigid, kinds) if k == 'shape'))
        self.inherent = sum(self._internal(self.plane, b) for b in self.rigid)
        self.upright = upright is not None
        self.root = tree.root
        children = {}
        for x, p in tree.parent.items():
            if p is not None:
                children.setdefault(p, []).append(x)
        below = {}

        def under(x):
            out, stack = {x}, [x]
            while stack:
                for y in children.get(stack.pop(), ()):
                    out.add(y)
                    stack.append(y)
            return frozenset(out)
        for x in sorted((x for x in tree.parent if tree.parent[x] is not None), key=rank.__getitem__):
            below[x] = under(x)
        self.edges = [(tree.parent[x], x, below[x]) for x in sorted(below, key=rank.__getitem__)
                      if not any(x in b and tree.parent[x] in b for b in self.rigid)]
        self.anchor = {}                                # (parent, child) -> the atom the child's angle is set against
        for p, x, _ in self.edges:
            a = tree.parent[p]
            if a is None:
                a = min((y for y in children.get(p, ()) if y != x), key=rank.__getitem__, default=None)
            if a is not None:
                self.anchor[(p, x)] = a
        for b, kind in sorted(zip(self.rigid, kinds), key=lambda t: self._depth(tree, min(t[0], key=lambda x:
                                                                                        self._depth(tree, x)))):
            self._enter(tree, children, below, b, kind)
        for b, kind in zip(self.rigid, kinds):
            self._exits(tree, children, below, b, kind)
        depth = {}
        for p, x, _ in self.edges:
            k, y = 0, p
            while y is not None:
                k, y = k + 1, tree.parent[y]
            depth[x] = k
        plane = self.plane
        for p, x, part in sorted(self.edges, key=lambda e: depth[e[1]]):     # a template's boundary bonds back to L
            dx, dy = plane[x][0] - plane[p][0], plane[x][1] - plane[p][1]
            d = (dx * dx + dy * dy) ** .5
            if d and abs(d - L) > 1e-9:
                ux, uy = dx / d * (L - d), dy / d * (L - d)
                for y in part:
                    plane[y] = (plane[y][0] + ux, plane[y][1] + uy)
        if upright is not None:                         # the largest shaped system stands as drawn in its shape
            b, ref = upright
            cx, cy = sum(plane[x][0] for x in b) / len(b), sum(plane[x][1] for x in b) / len(b)
            rx, ry = sum(ref[x][0] for x in b) / len(b), sum(ref[x][1] for x in b) / len(b)
            dot = sum((plane[x][0] - cx) * (ref[x][0] - rx) + (plane[x][1] - cy) * (ref[x][1] - ry) for x in b)
            crs = sum((plane[x][0] - cx) * (ref[x][1] - ry) - (plane[x][1] - cy) * (ref[x][0] - rx) for x in b)
            c, s_ = cos(atan2(crs, dot)), sin(atan2(crs, dot))
            for y in plane:
                dx, dy = plane[y][0] - cx, plane[y][1] - cy
                plane[y] = (cx + c * dx - s_ * dy, cy + s_ * dx + c * dy)
        self.terms = terms(self.scope, self.plane)
        self.ct = _ct_broken(mol, self.plane, self.scope.frames, set(self.plane))

    @staticmethod
    def _depth(tree, x):
        k = 0
        while tree.parent[x] is not None:
            k, x = k + 1, tree.parent[x]
        return k

    def _directions(self, a, block, k):
        """`k` exit directions of the template atom `a`, substituent out: the widest outside gap split evenly; at a
        fusion atom with no such gap and one exit, the extension of the fused bond; None when neither holds."""
        scope, plane = self.scope, self.plane
        dirs = outside_angles(scope, plane, a, k, ring=block)
        if dirs:
            return dirs
        if k == 1:
            fused = [y for y in scope.adj[a] if y in block and sum(a in r and y in r for r in scope.rings) > 1]
            if fused:
                y = min(fused, key=scope.rank.__getitem__)
                return [atan2(plane[a][1] - plane[y][1], plane[a][0] - plane[y][0])]
        return None

    def _outline(self, a, block, around):
        """`around` less `a`'s neighbours drawn inside `block`'s outline, while two of `block` on it remain."""
        from .shapes import _inside
        rank, plane = self.scope.rank, self.plane
        pts = sorted(block, key=rank.__getitem__)
        hidden = {pts[i] for i in _inside([plane[x] for x in pts])}
        outline = [y for y in around if y not in hidden]
        return outline if sum(y in block for y in outline) >= 2 else around

    def _away(self, a, block, outline=False):
        """The direction from `a` away from its neighbours in `block`, summed in rank order, those inside the
        drawing's `_outline` left out when `outline`; None when they cancel."""
        plane, ux, uy = self.plane, 0., 0.
        around = sorted((y for y in self.mol.neighbors_of(a) if y in block), key=self.scope.rank.__getitem__)
        for y in self._outline(a, block, around) if outline else around:
            if y in block:
                dx, dy = plane[y][0] - plane[a][0], plane[y][1] - plane[a][1]
                d = (dx * dx + dy * dy) ** .5 or 1.
                ux, uy = ux - dx / d, uy - dy / d
        return atan2(uy, ux) if abs(ux) + abs(uy) > 1e-9 else None

    def _turn_part(self, part, o, t):
        plane, (ox, oy) = self.plane, self.plane[o]
        c, s = cos(t), sin(t)
        for y in part:
            dx, dy = plane[y][0] - ox, plane[y][1] - oy
            plane[y] = (ox + c * dx - s * dy, oy + s * dx + c * dy)

    def _enter(self, tree, children, below, block, kind):
        """Turn a template hanging from the tree about its entry atom so its bond to the tree leaves along the entry
        atom's first exit direction, or away from the template for a shape."""
        e = min(block, key=lambda x: self._depth(tree, x))
        q = tree.parent[e]
        if q is None:
            return
        plane = self.plane
        part = frozenset().union(*(below[c] for c in children.get(e, ()) if c in block))   # a spiro atom's other
        other = frozenset().union(*(b for b in self.rigid if e in b and b is not block))      # system stays
        if other:                                       # spiro: back to back, this system's outside facing the other
            mine, theirs = self._away(e, block), self._away(e, other)
            if mine is not None and theirs is not None:
                self._turn_part(part, e, theirs + pi - mine)
            return
        outside = sum(1 for y in self.mol.neighbors_of(e) if y not in block)
        dirs = None if kind == 'shape' else self._directions(e, block, outside)
        t = dirs[0] if dirs else self._away(e, block, kind == 'shape')
        if t is not None:
            self._turn_part(part, e, atan2(plane[q][1] - plane[e][1], plane[q][0] - plane[e][0]) - t)

    def _exits(self, tree, children, below, block, kind):
        """Turn each substituent of a template atom, with its subtree, onto an exit: a polygon template's atoms
        take `_directions`, the entry atom keeping the one its tree bond already takes; a shape's atoms, or one with
        no direction, take the most room."""
        rank, plane = self.scope.rank, self.plane
        for a in sorted(block, key=rank.__getitem__):
            held = frozenset().union(*(b for b in self.rigid if a in b))   # a spiro atom's other system is placed
            kids = sorted((x for x in children.get(a, ()) if x not in held), key=rank.__getitem__)
            if not kids:
                continue
            parent = tree.parent[a] if tree.parent[a] is not None and tree.parent[a] not in held else None
            dirs = None if kind == 'shape' else self._directions(a, held, len(kids) + (parent is not None))
            if dirs and parent is not None:
                up = atan2(plane[parent][1] - plane[a][1], plane[parent][0] - plane[a][0])
                dirs.remove(min(dirs, key=lambda t: abs((t - up + pi) % (2 * pi) - pi)))
            if not dirs:
                self._roomiest(tree, below, held, a, kids, kind != 'shape')
                continue
            for x, t in zip(kids, dirs):
                self._turn_part(below[x], a, t - atan2(plane[x][1] - plane[a][1], plane[x][0] - plane[a][0]))
            for x in kids:                              # substituent out, unless the most room has fewer defects
                if self._room(below[x], a, 0.) >= CLEAR or self._clear_by_flip(tree, below, held | block, a, x) or \
                        self._clear_by_stretch(tree, below, held | block, a, x):
                    continue
                fixed = {y: plane[y] for y in below[x]}
                self._roomiest(tree, below, held, a, [x])
                moved = {y: plane[y] for y in below[x]}
                plane.update(fixed)
                if delta_terms(self.scope, plane, moved)[0] < 0:
                    plane.update(moved)

    def _clear_by_flip(self, tree, below, fixed, a, x):
        """True when a flip of a placed acyclic substituent in the way, about one of its tree bonds below a template,
        or, for one in the way above it, of the template on its entry bond, leaves `x`'s subtree on its exit from `a`
        at least `CLEAR` of every other atom and adds no defect; the flip is kept then."""
        plane, scope = self.plane, self.scope
        part = below[x]
        others = Grid()
        for y in plane:
            if y not in part and y != a:
                others.add(y, plane[y])
        blockers = sorted({z for y in part for z in others.atoms(plane[y], CLEAR)
                           if z not in fixed and ((plane[y][0] - plane[z][0]) ** 2 +
                                                  (plane[y][1] - plane[z][1]) ** 2) ** .5 < CLEAR},
                          key=scope.rank.__getitem__)
        for z in blockers:
            v = z
            while tree.parent[v] is not None and v not in fixed:   # up to the bond leaving the template
                u = tree.parent[v]
                if frozenset((u, v)) in scope.ring_bonds or a in below[v] or x in below[v]:
                    break
                flip = _mirrored(plane, below[v], plane[u], plane[v])
                if flip and delta_terms(scope, plane, flip)[0] <= 0:
                    old = {y: plane[y] for y in flip}
                    plane.update(flip)
                    if self._room(part, a, 0.) >= CLEAR:
                        return True
                    plane.update(old)
                v = u
        for e in sorted(fixed, key=scope.rank.__getitem__):    # a blocker above: the template mirrored on its entry
            u = tree.parent[e]
            if u is None or u in fixed or a not in below[e] or not any(z not in below[e] for z in blockers):
                continue
            flip = _mirrored(plane, below[e], plane[u], plane[e])
            if flip and delta_terms(scope, plane, flip)[0] <= 0:
                old = {y: plane[y] for y in flip}
                plane.update(flip)
                if self._room(part, a, 0.) >= CLEAR:
                    return True
                plane.update(old)
        return False

    def _clear_by_stretch(self, tree, below, fixed, a, x):
        """True when a substituent-out exit of two to `SMALL` atoms, `x`'s from `a` or that of a substituent in its
        way, moved out along its exit by up to `L` in `STRETCH` steps, the larger first, leaves `x`'s subtree `CLEAR`
        of every other atom and adds no defect; the first such move is kept."""
        plane, rank = self.plane, self.scope.rank
        part = below[x]
        roots = {(a, x)}
        for z in plane:
            if z in part or z in fixed or z == a or all(
                    ((plane[y][0] - plane[z][0]) ** 2 + (plane[y][1] - plane[z][1]) ** 2) ** .5 >= CLEAR for y in part):
                continue
            v = z                                       # the blocker's exit: its tree bond leaving the template
            while tree.parent[v] is not None and tree.parent[v] not in fixed:
                v = tree.parent[v]
            if tree.parent[v] is not None and a not in below[v]:
                roots.add((tree.parent[v], v))
        for u, v in sorted((e for e in roots if 1 < len(below[e[1]]) <= SMALL),
                           key=lambda e: (-len(below[e[1]]), rank[e[1]])):
            (ux, uy), (vx, vy) = plane[u], plane[v]
            d = ((vx - ux) ** 2 + (vy - uy) ** 2) ** .5
            ex, ey = (vx - ux) / d, (vy - uy) / d
            old = {y: plane[y] for y in below[v]}
            for k in range(1, int(L / STRETCH + 1e-9) + 1):
                moved = {y: (q[0] + ex * k * STRETCH, q[1] + ey * k * STRETCH) for y, q in old.items()}
                if delta_terms(self.scope, plane, moved)[0] > 0:
                    continue
                plane.update(moved)
                if self._room(part, a, 0.) >= CLEAR:
                    return True
                plane.update(old)
        return False

    def _room(self, part, a, turn, others=None):
        """The least distance, up to `ROOM`, from `part` turned by `turn` about `a` to every other atom."""
        plane = self.plane
        if others is None:
            others = Grid()
            for y in plane:
                if y not in part and y != a:
                    others.add(y, plane[y])
        ax, ay = plane[a]
        c, s = cos(turn), sin(turn)
        room = ROOM
        for y in part:
            dx, dy = plane[y][0] - ax, plane[y][1] - ay
            px, py = ax + c * dx - s * dy, ay + s * dx + c * dy
            for z in others.atoms((px, py), room):
                q = plane[z]
                room = min(room, ((px - q[0]) ** 2 + (py - q[1]) ** 2) ** .5)
        return room

    def _roomiest(self, tree, below, block, a, kids, flat=True):
        """Turn each of `kids`, in rank order, to whichever of `EXITS` directions leaves its subtree most room from
        every other atom, up to `ROOM`, then is nearest the direction away from `block` and `a`'s tree parent, a shape's
        neighbours inside its outline left out while two on it remain; two at least `APART`; never, when `flat`,
        `in_ring_angle`; room counted up to `AWAY` for a shape's atom of three or more bonds in it or of a double or
        aromatic one."""
        rank, plane, mol = self.scope.rank, self.plane, self.mol
        held = sorted((y for y in mol.neighbors_of(a) if y in block or y == tree.parent[a]), key=rank.__getitem__)
        ax, ay = plane[a]
        ux = uy = 0.
        away = held if flat else self._outline(a, block, held)     # a shape's atom bisects its outline bonds
        for y in away:
            dx, dy = plane[y][0] - ax, plane[y][1] - ay
            d = (dx * dx + dy * dy) ** .5 or 1.
            ux, uy = ux - dx / d, uy - dy / d
        if abs(ux) + abs(uy) < 1e-9:
            cx, cy = sum(plane[y][0] for y in block) / len(block), sum(plane[y][1] for y in block) / len(block)
            ux, uy = ax - cx, ay - cy
        base = atan2(uy, ux)
        inside = [y for y in held if y in block]        # a shape's bridgehead or flat ring atom: straight out wins
        cap = ROOM if flat or len(inside) == 2 and all(mol.bond(a, y).order == 1 for y in inside) else AWAY
        taken = []
        for x in kids:
            part = below[x]
            others = Grid()                             # the rest; only atoms within `ROOM` change a key
            for y in plane:
                if y not in part and y != a:
                    others.add(y, plane[y])
            best = None
            for i in range(EXITS):
                t = base + 2 * pi * i / EXITS
                if any(abs((t - u + pi) % (2 * pi) - pi) < APART for u in taken):
                    continue
                if flat and in_ring_angle(self.scope, plane, a, t):
                    continue
                room = self._room(part, a, t - atan2(plane[x][1] - ay, plane[x][0] - ax), others)
                key = (-round(min(room, cap) / L, 2), round(abs((t - base + pi) % (2 * pi) - pi), 6))
                if best is None or key < best[0]:
                    best = (key, t)
            if best is not None:
                taken.append(best[1])
                self._turn_part(part, a, best[1] - atan2(plane[x][1] - ay, plane[x][0] - ax))

    def _internal(self, plane, block):
        """The defects among the atoms and bonds of `block` drawn as `plane`."""
        return _parts(self.scope, plane, set(block), [x for x in self.scope.atoms if x in block],
                      [(n, m) for n, m in self.scope.bonds if n in block and m in block])[0]

    def env(self):
        """A fresh environment at the start."""
        return PrimitiveEnv(self)


class PrimitiveEnv:
    """One component from its tree embedding to `stop`, by primitive steps."""
    primitive = True
    __slots__ = ('initial', 'mol', 'scope', 'plane', 'terms', 'ct', 'steps', 'cap', 'done', 'moves', 'blocks',
                 'drawn', 'everywhere', '_legal')

    def __init__(self, initial):
        self.initial = initial
        self.mol, self.scope = initial.mol, initial.scope
        self.plane, self.terms, self.ct = dict(initial.plane), initial.terms, initial.ct
        self.steps, self.cap, self.done, self.moves = 0, min(CAP * len(initial.scope.atoms) + 8, LIMIT), False, []
        self.blocks, self.drawn, self.everywhere = [], [], True      # read by `features`
        self._legal = None

    @property
    def defects(self):
        """The defects of the drawing beyond those its templates draw."""
        return self.terms[0] - self.initial.inherent

    @property
    def score(self):
        return combine(self.terms)

    @property
    def rank(self):
        return self.ct, key(self.score)

    def tie(self, move):
        if move.kind == 'stop':
            return (len(KINDS),)
        rank = self.scope.rank
        return KINDS[move.kind], rank.get(move.pivot, -1), rank.get(move.other, -1), round(move.arg, 9)

    def fork(self):
        f = PrimitiveEnv.__new__(PrimitiveEnv)
        f.initial, f.mol, f.scope = self.initial, self.mol, self.scope
        f.plane, f.terms, f.ct = dict(self.plane), self.terms, self.ct
        f.steps, f.cap, f.done, f.moves = self.steps, self.cap, self.done, list(self.moves)
        f.blocks, f.drawn, f.everywhere, f._legal = [], [], True, None
        return f

    def _moved(self, m):
        plane = self.plane
        if m.kind == 'turn':
            n = len(plane)
            cx, cy = sum(p[0] for p in plane.values()) / n, sum(p[1] for p in plane.values()) / n
            c, s = cos(m.arg), sin(m.arg)
            return {x: (cx + c * (p[0] - cx) - s * (p[1] - cy), cy + s * (p[0] - cx) + c * (p[1] - cy))
                    for x, p in plane.items()}
        (px, py), (qx, qy) = plane[m.pivot], plane[m.other]
        if m.kind == 'flip':
            dx, dy = qx - px, qy - py
            d = dx * dx + dy * dy or 1.
            out = {}
            for x in m.part:
                ux, uy = plane[x][0] - px, plane[x][1] - py
                k = 2 * (ux * dx + uy * dy) / d
                out[x] = (px + k * dx - ux, py + k * dy - uy)
            return out
        if m.kind == 'rotate':
            t = self._turn(m)
            c, s = cos(t), sin(t)
            return {x: (px + c * (plane[x][0] - px) - s * (plane[x][1] - py),
                        py + s * (plane[x][0] - px) + c * (plane[x][1] - py)) for x in m.part}
        d = ((qx - px) ** 2 + (qy - py) ** 2) ** .5 or 1.         # stretch
        dx, dy = (qx - px) / d * m.arg, (qy - py) / d * m.arg
        return {x: (plane[x][0] + dx, plane[x][1] + dy) for x in m.part}

    def _turn(self, m):
        """The turn about `m.pivot` that sets the angle anchor-pivot-other to `m.arg`, positive anticlockwise."""
        plane = self.plane
        a, (px, py), (qx, qy) = plane[self.initial.anchor[(m.pivot, m.other)]], plane[m.pivot], plane[m.other]
        t = atan2(a[1] - py, a[0] - px) + m.arg - atan2(qy - py, qx - px)
        return (t + pi) % (2 * pi) - pi

    def _moves(self):
        out = []
        plane = self.plane
        for p, x, part in self.initial.edges:
            out.append(Move('flip', p, x, 0., part))
            if (p, x) in self.initial.anchor:
                for t in ANGLES:
                    for side in ((t, -t) if t < pi - 1e-9 else (t,)):
                        m = Move('rotate', p, x, side, part)
                        if abs(self._turn(m)) > 1e-9:          # an angle already drawn is no step
                            out.append(m)
            d = ((plane[x][0] - plane[p][0]) ** 2 + (plane[x][1] - plane[p][1]) ** 2) ** .5
            out += [Move('stretch', p, x, t, part) for t in (STEP, -STEP) if L - 1e-9 <= d + t <= LONGEST + 1e-9]
        out += [Move('turn', -1, -1, t, frozenset(plane)) for t in SCENE]
        return out

    def legal(self):
        """The legal steps of the state, `stop` last when offered; computed once per state."""
        if self._legal is not None:
            return self._legal
        if self.done:
            self._legal = []
            return self._legal
        out = []
        if self.steps < self.cap:
            for m in self._moves():
                moved = self._moved(m)
                if m.kind != 'turn' and masked(self.scope, self.plane, moved, LONGEST):
                    continue
                d = delta_terms(self.scope, self.plane, moved)
                ct = self.ct if m.kind == 'turn' else self.ct - _ct_broken(
                    self.mol, self.plane, self.scope.frames, set(moved)) + \
                    _ct_broken(self.mol, {**self.plane, **moved}, self.scope.frames, set(moved))
                out.append(Step(m, moved, (ct, key(add(self.score, combine(d)))), d, ct))
        if out and self.defects > 0:                        # a drawing with a defect is not done
            self._legal = out
        else:
            self._legal = out + [Step(STOP, {}, self.rank, (0,) * len(TERMS), self.ct)]
        return self._legal

    def step(self, s):
        if s.move.kind == 'stop':
            self.done = True
        else:
            self.plane.update(s.moved)
            self.terms = tuple(a + b for a, b in zip(self.terms, s.terms))
            self.ct = s.ct
            self.steps += 1
        self.moves.append(s.move)
        self._legal = None


def _mirrored(plane, part, p, q):
    """`part` reflected across the line through `p` and `q`: `{atom: (x, y)}`, empty for a part on the line."""
    dx, dy = q[0] - p[0], q[1] - p[1]
    d = dx * dx + dy * dy
    out = {}
    for y in part:
        rx, ry = plane[y][0] - p[0], plane[y][1] - p[1]
        t = (rx * dx + ry * dy) / d
        fx, fy = p[0] + t * dx, p[1] + t * dy
        out[y] = (2 * fx - plane[y][0], 2 * fy - plane[y][1])
    return {y: v for y, v in out.items() if abs(v[0] - plane[y][0]) + abs(v[1] - plane[y][1]) > 1e-9}
