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
"""Chains: a backbone of units with groups hanging from them and cross-links between them, the segmentation every
chain layout reads.

| field | meaning |
| --- | --- |
| `residues` | the units: peptide residues, a cycle's atoms and beads |
| `members` | `members[i]` is unit i's atoms plus the groups hanging from it |
| `main` | the units of the main chain, in order |
| `branches` | attached chains of a peptide |
| `crosslinks` | atom pairs drawn as chords; when `head_to_tail` the last is the closure `(last exit, first entry)` |
| `head_to_tail` | the chain closes on itself |
| `backbone` | the main chain's atoms, consecutive ones bonded or a bead's entry and exit |
| `beads` | indices `k` of the backbone steps `backbone[k]`-`backbone[k + 1]` that cross a bead, entry to exit |
| `kind` | `'peptide'` or `'cycle'` |
"""
from dataclasses import dataclass


__all__ = ['Chain', 'Circle', 'Residue', 'cycles', 'perceive', 'MACRO', 'BEAD', 'BEAD_RINGS', 'TEMPLATES']

MACRO = 13              # the smallest ring a cycle is
BEAD = 8                # the largest ring of a bead
BEAD_RINGS = 3          # the most rings a bead holds; more is a supercondensed system the cascade keeps
TEMPLATES = (           # macrocycle skeletons drawn by depict/layout/templates.py, refused here
    'C1=CC2=CC3=CC=C(N3)C=C4C=CC(=N4)C=C5C=CC(=N5)C=C1N2',                      # porphine
    'C1=CC2=CC3=CC=C(N3)C=C4C=C6C(=N4)C7=C5C=CC(=N5)C=C1N2.C6C7',               # phorbine
)


@dataclass(frozen=True, slots=True)
class Residue:
    """One unit.  `path` runs entry to exit through the core; `symbol` comes from a library match."""
    atoms: frozenset
    kind: str
    entry: int | None
    exit: int | None
    path: tuple
    symbol: str | None = None


@dataclass(frozen=True, slots=True)
class Circle:
    """A crown-type macrocycle drawn as a circle of its repeating units: its main cycle's atoms in ring order from
    its lowest-ranked heteroatom, and the atoms of one unit."""
    ring: tuple
    period: int


@dataclass(frozen=True, slots=True)
class Chain:
    """A segmentation: see the module table."""
    residues: tuple
    members: tuple
    main: tuple
    branches: tuple
    crosslinks: tuple
    head_to_tail: bool
    backbone: tuple
    beads: tuple = ()
    kind: str = 'peptide'


def perceive(mol):
    """The chain or circle `mol` is, or None: a peptide first, then `cycles`."""
    from .monomers import monomers                      # monomers builds Chains: imported at call time
    return monomers(mol) or cycles(mol)


def cycles(mol):
    """`mol`'s macrocycle as a `Chain` of kind cycle, a `Circle`, or None.

    | rule | value |
    | --- | --- |
    | main cycle | the shortest SSSR ring of `MACRO` or more atoms, ties to the lower ranks |
    | bead | a system of rings of up to `BEAD` atoms the cycle runs through; more than `BEAD_RINGS` rings refuses |
    | template | the main cycle's block a skeleton of `TEMPLATES` refuses, as does a porphyrinoid: every crossed bead \
a five-membered ring with a nitrogen between its entry and exit |
    | circle | no bead entered and left at two unbonded atoms, and the cycle a crown: a unit of 2 to 4 atoms, one a \
heteroatom, repeated three or more times |
    | opening | a C(=O)-O or C(=O)-N bond of the cycle, else one between a bead and another unit, else any; ranks |
    | units | an atom on the cycle; a crossed hexagon entered meta, its three-atom arc on the backbone; any other \
crossed bead, with its entry and exit |
    | groups | every other atom, claimed in rank order from the units without crossing the backbone |
    | links | a bond between two groups off the backbone, moved to start on it; the closure last |

    A mixture is no cycle: each component is asked on its own.
    """
    from .monomers import _carbonyl, _to_backbone       # monomers imports this module: at call time
    if len(mol.connected_components) != 1:              # callers split a mixture into its components
        return None
    rank = mol.canonical_order()
    rings = [tuple(r) for r in mol.sssr]
    macro = [r for r in rings if len(r) >= MACRO]
    if not macro:
        return None
    main = min(macro, key=lambda r: (len(r), sorted(rank[x] for x in r)))
    if _templated(mol, main, rings):
        return None
    systems = []                                        # (atoms, rings) of the small-ring systems
    for r in rings:
        if len(r) > BEAD:
            continue
        atoms, held = set(r), [r]
        for s in [s for s in systems if s[0] & atoms]:
            systems.remove(s)
            atoms |= s[0]
            held += s[1]
        systems.append((atoms, held))
    system = {x: i for i, (s, _) in enumerate(systems) for x in s}
    ring = _ordered(mol, main, rank)
    visits = _visits(ring, system)
    if visits is None:
        return None
    crossed = [v for v in visits if v[0] is not None and len(v[1]) > 1 and
               v[1][-1] not in mol.neighbors_of(v[1][0])]
    if any(len(systems[v[0]][1]) > BEAD_RINGS for v in crossed):
        return None
    if crossed and all(_pyrrolic(mol, systems[v[0]], v[1]) for v in crossed):
        return None                                     # a porphyrinoid: corrin, chlorin, corrole, porphine
    if not crossed:
        period = _period(mol, ring)
        if period:                                      # a crown: from its lowest-ranked heteroatom
            k = min((i for i in range(len(ring)) if mol.element_of(ring[i]) != 6), key=lambda i: rank[ring[i]])
            ring = ring[k:] + ring[:k]
            if rank[ring[-1]] < rank[ring[1]]:
                ring = ring[:1] + ring[1:][::-1]
            return Circle(tuple(ring), period)
    inside = {frozenset(p) for v in crossed for p in zip(v[1], v[1][1:])}
    bead = {x for v in crossed for x in v[1]}
    n = len(ring)

    def opening(i):
        a, b = ring[i], ring[(i + 1) % n]
        acyl = any(_carbonyl(mol, c) and mol.element_of(o) in (7, 8) for c, o in ((a, b), (b, a)))
        return (not acyl, (a in bead) == (b in bead), sorted((rank[a], rank[b])))
    cut = min((i for i in range(n) if frozenset((ring[i], ring[(i + 1) % n])) not in inside), key=opening)
    ring = ring[cut + 1:] + ring[:cut + 1]
    visits = _visits(ring, system, rotate=False)
    residues, backbone, beads = [], [], []
    for s, path in visits:
        if s is not None and len(path) == 3 and len(systems[s][1]) == 1 and len(systems[s][1][0]) == 6:
            backbone += path                            # a meta hexagon: its arc on the lattice, the ring above
            residues.append(Residue(frozenset(systems[s][0]), 'arc', path[0], path[-1], tuple(path)))
        elif s is not None and len(path) > 1 and path[-1] not in mol.neighbors_of(path[0]):
            beads.append(len(backbone))
            backbone += [path[0], path[-1]]
            residues.append(Residue(frozenset(systems[s][0]), 'bead', path[0], path[-1], tuple(path)))
        else:
            for x in path:
                backbone.append(x)
                residues.append(Residue(frozenset((x,)), 'unit', x, x, (x,)))
    members = _claim(mol, residues, set(backbone), rank)
    group = {x: i for i, m in enumerate(members) for x in m}
    steps = {frozenset(p) for p in zip(backbone, backbone[1:])}
    links = sorted({tuple(sorted((x, y), key=rank.__getitem__)) for x in mol for y in mol.neighbors_of(x)
                    if group[x] != group[y] and frozenset((x, y)) not in steps and
                    frozenset((x, y)) != frozenset((backbone[-1], backbone[0]))},
                   key=lambda t: (rank[t[0]], rank[t[1]]))
    links = [_to_backbone(mol, members, set(backbone), link, rank) for link in links]
    links.append((backbone[-1], backbone[0]))
    return Chain(tuple(residues), tuple(frozenset(m) for m in members), tuple(range(len(residues))), (),
                 tuple(links), True, tuple(backbone), tuple(beads), 'cycle')


def _period(mol, ring):
    """The atoms of a crown's repeating unit, 2 to 4 of them, one a heteroatom, three or more units round `ring`;
    0 when `ring` is no crown."""
    n = len(ring)
    het = [mol.element_of(x) != 6 for x in ring]
    for p in (2, 3, 4):
        if n % p or n // p < 3:
            continue
        if all(het[i] == het[(i + p) % n] for i in range(n)) and sum(het[:p]) == 1:
            return p
    return 0


def _pyrrolic(mol, system, path):
    """True for a bead of one five-membered ring whose entry and exit flank a nitrogen: a porphyrinoid's unit."""
    atoms, held = system
    if len(held) != 1 or len(held[0]) != 5:
        return False
    a, b = path[0], path[-1]
    return any(mol.element_of(x) == 7 and a in mol.neighbors_of(x) and b in mol.neighbors_of(x) for x in atoms)


def _ordered(mol, ring, rank):
    """`ring`'s atoms in ring order from its lowest rank, toward the lower-ranked of that atom's two ring
    neighbours."""
    ring = list(ring)
    k = min(range(len(ring)), key=lambda i: rank[ring[i]])
    ring = ring[k:] + ring[:k]
    if rank[ring[-1]] < rank[ring[1]]:
        ring = ring[:1] + ring[1:][::-1]
    return ring


def _visits(ring, system, rotate=True):
    """`[(system or None, atoms)]`: `ring` split into runs of one small-ring system and single other atoms, the
    runs in ring order; None when a system is visited twice. With `rotate`, a run wrapping the ring's start is
    joined first."""
    if rotate:
        k = 0
        while k < len(ring) and system.get(ring[k]) is not None and system.get(ring[k]) == system.get(ring[k - 1]):
            k += 1
        if k == len(ring):
            return None
        ring = ring[k:] + ring[:k]
    out = []
    for x in ring:
        s = system.get(x)
        if out and s is not None and out[-1][0] == s:
            out[-1][1].append(x)
        else:
            out.append((s, [x]))
    seen = [s for s, _ in out if s is not None]
    return None if len(seen) != len(set(seen)) else out


def _claim(mol, residues, backbone, rank):
    """Each residue's atoms with the atoms it reaches first, in rank order, never through a backbone atom."""
    owner = {x: i for i, r in enumerate(residues) for x in r.atoms}
    queue = sorted(owner, key=rank.__getitem__)
    while queue:
        nxt = []
        for x in queue:
            for y in sorted(mol.neighbors_of(x), key=rank.__getitem__):
                if y not in owner and y not in backbone:
                    owner[y] = owner[x]
                    nxt.append(y)
        queue = nxt
    members = [set() for _ in residues]
    for x, i in owner.items():
        members[i].add(x)
    return members


def _templated(mol, main, rings):
    """True when the block of rings sharing atoms with `main` has the skeleton of one of `TEMPLATES`."""
    from . import smiles
    block, grew = set(main), True
    while grew:
        grew = False
        for r in rings:
            if block.intersection(r) and not block.issuperset(r):
                block.update(r)
                grew = True
    mine = _skeleton(mol, block)
    return any(_skeleton(t, set(t)) == mine for t in map(smiles, TEMPLATES))


def _skeleton(mol, atoms):
    """The canonical bytes of `atoms` and their bonds as an all-carbon, all-single-bond graph."""
    from . import MoleculeContainer
    out = MoleculeContainer()
    with out.edit() as e:
        index = {x: e.add_atom('C') for x in sorted(atoms)}
        for x in atoms:
            for y in mol.neighbors_of(x):
                if y in index and x < y:
                    e.add_bond(index[x], index[y], 1)
    return out.canonical_bytes
