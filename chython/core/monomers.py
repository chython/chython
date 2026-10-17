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
"""Peptide segmentation: residues split at amide bonds, the main chain, branches and cross-links.

Read-only topology, so it sits in `core` beside `wedge.py`: `chemistry` exports it and `depict` lays it
out, and neither may import the other.  Exported as `chython.chemistry.monomers`.

| Step | Rule |
| --- | --- |
| amide cut | single C(=O)-N, N not aromatic, N with another heavy neighbour, bond in no ring of size <= 7 |
| residue core | a 3/4/5-atom path (alpha/beta/gamma) through sp3 carbons from an entry N to an exit C |
| partition | multi-source BFS from the cores, a small ring or a non-single bond claimed whole; cross-links between; \
a link beside an acyclic C=C moves a bond off it |
| chain link | a cut from residue A's exit to residue B's entry, each used once |
| main chain | the chain with most non-`other` residues; fewer than `MIN_RESIDUES` is not a peptide |
| branch | an attached chain of >= `MIN_BRANCH` residues; a shorter one merges into its host residue |
| order | every choice and every walk goes by `canonical_order()` rank, so the result ignores atom numbering |
"""
from collections import deque
from dataclasses import dataclass

from .chains import Chain, Residue


__all__ = ['PeptideBranch', 'Peptide', 'Residue', 'monomers', 'MIN_RESIDUES', 'MIN_BRANCH']

MIN_RESIDUES = 4
MIN_BRANCH = 3
_KINDS = {3: 'alpha', 4: 'beta', 5: 'gamma'}


@dataclass(frozen=True, slots=True)
class PeptideBranch:
    """A chain hung on `host`, an atom of residue `residue`.  `path` is its backbone, host side first."""
    host: int
    residue: int
    chain: tuple
    path: tuple


Peptide = Chain                 # the segmentation of a peptide is a `Chain` of kind peptide


def monomers(mol, library=None) -> Peptide | None:
    """Segment `mol` into residues, or `None` when its main chain has fewer than `MIN_RESIDUES`.

    :param library: `{symbol: MoleculeContainer}` of free amino acids.  A residue whose atoms match an
        entry minus its C-terminal OH gets that symbol; symbols never move a boundary.
    """
    rank = mol.canonical_order()
    small = set()
    for ring in mol.sssr:
        if len(ring) <= 7:
            small.update(frozenset(b) for b in zip(ring, ring[1:] + ring[:1]))
    cuts = sorted(_amide_cuts(mol, small), key=lambda t: (rank[t[0]], rank[t[1]]))
    if len(cuts) < MIN_RESIDUES - 1:
        return None
    residues, crosslinks = _residues(mol, cuts, _units(mol, small), rank)
    residue_of = {x: i for i, r in enumerate(residues) for x in r.atoms}

    nxt, prv, attach = {}, {}, []
    for c, n in cuts:
        a, b = residue_of[c], residue_of[n]
        if residues[a].exit == c and residues[b].entry == n and a not in nxt and b not in prv:
            nxt[a], prv[b] = b, a
        else:
            attach.append((c, n))
    chains = _chains(residues, nxt, prv, rank)
    main, cyclic = chains[0]
    if sum(residues[i].kind != 'other' for i in main) < MIN_RESIDUES:
        return None

    members = [set(r.atoms) for r in residues]
    home = {i: i for i in range(len(residues))}            # residue -> the residue whose group holds it
    chain_of = {i: ci for ci, (ch, _) in enumerate(chains) if ci for i in ch}
    placed, done, used, branches = set(main), set(), set(), []
    todo = deque(main)
    while todo:
        i = todo.popleft()
        for c, n in attach:
            for h, b in ((c, n), (n, c)):
                if residue_of[h] != i or residue_of[b] not in chain_of or chain_of[residue_of[b]] in done:
                    continue
                ci = chain_of[residue_of[b]]
                ch = chains[ci][0]
                done.add(ci)
                used.add((c, n))
                if residue_of[b] == ch[-1]:
                    ch = ch[::-1]
                elif residue_of[b] != ch[0]:
                    crosslinks.append((h, b))
                    continue
                if len(ch) >= MIN_BRANCH:
                    branches.append(PeptideBranch(h, i, tuple(ch),
                                                  tuple(_oriented(mol, residues, members, ch, h, rank))))
                    placed.update(ch)
                    todo.extend(ch)
                else:
                    for j in ch:
                        members[home[i]] |= residues[j].atoms
                        home[j] = home[i]
                    placed.update(ch)
    # an amide between two placed residues that is no chain link, e.g. a side-chain lactam
    for c, n in attach:
        if (c, n) not in used and residue_of[c] in placed and residue_of[n] in placed:
            crosslinks.append((c, n))

    backbone = _oriented(mol, residues, members, main, None, rank)
    if len(main) > 1:
        p0, p1 = residues[main[0]].path, residues[main[1]].path
        if p0 and p1 and p0[0] in mol.neighbors_of(p1[0]):
            backbone = list(p0[::-1]) + backbone[len(p0):]
    crosslinks = [_to_backbone(mol, members, set(backbone), link, rank) for link in crosslinks]
    if cyclic:
        crosslinks.append((residues[main[-1]].exit, residues[main[0]].entry))
    if library:
        residues = _symbols(mol, residues, library)
    return Peptide(tuple(residues), tuple(frozenset(m) for m in members), tuple(main), tuple(branches),
                   tuple(crosslinks), cyclic, tuple(backbone))


def _heavy(mol, n):
    return [x for x in mol.neighbors_of(n) if mol.element_of(x) != 1]


def _carbonyl(mol, c):
    return mol.element_of(c) == 6 and any(mol.element_of(o) == 8 and mol.order_of(c, o) == 2
                                          for o in mol.neighbors_of(c))


def _amide_cuts(mol, small):
    cuts = []
    for b in mol.bonds():
        if b.order != 1 or frozenset((b.n, b.m)) in small:
            continue
        for c, n in ((b.n, b.m), (b.m, b.n)):
            if mol.element_of(n) == 7 and mol.hybridization_of(n) != 4 and _carbonyl(mol, c) \
                    and any(x != c for x in _heavy(mol, n)):
                cuts.append((c, n))
    return cuts


def _neighbors(mol, x, rank):
    """`x`'s neighbours in `rank` order."""
    return sorted(mol.neighbors_of(x), key=rank.__getitem__)


def _bfs_path(mol, a, b, allowed, rank):
    """Shortest a..b path with every interior atom in `allowed`, or None; ties go to the lower rank."""
    prev, queue = {a: None}, deque([a])
    while queue:
        x = queue.popleft()
        if x == b:
            path = [b]
            while prev[path[-1]] is not None:
                path.append(prev[path[-1]])
            return path[::-1]
        for y in _neighbors(mol, x, rank):
            if y not in prev and (y == b or y in allowed):
                prev[y] = x
                queue.append(y)
    return None


def _farthest_path(mol, a, atoms, rank):
    """The path from `a` to the last heavy atom BFS reaches inside `atoms`, stopping at ring atoms."""
    prev, queue, last = {a: None}, deque([a]), a
    while queue:
        x = queue.popleft()
        last = x
        for y in _neighbors(mol, x, rank):
            if y in atoms and y not in prev and not mol.in_ring_of(y) and mol.element_of(y) != 1:
                prev[y] = x
                queue.append(y)
    path = [last]
    while prev[path[-1]] is not None:
        path.append(prev[path[-1]])
    return path[::-1]


def _fragments(mol, cut_set, rank):
    seen, out = set(), []
    for s in sorted(mol, key=rank.__getitem__):
        if s in seen:
            continue
        comp, stack = set(), [s]
        seen.add(s)
        while stack:
            x = stack.pop()
            comp.add(x)
            for y in mol.neighbors_of(x):
                if y not in seen and frozenset((x, y)) not in cut_set:
                    seen.add(y)
                    stack.append(y)
        out.append(comp)
    return out


def _units(mol, small):
    """`{atom: frozenset}`: atoms joined by a bond in a small ring or of order other than 1."""
    parent = {}

    def find(x):
        while parent.get(x, x) != x:
            x = parent[x]
        return x
    for b in mol.bonds():
        if b.order != 1 or frozenset((b.n, b.m)) in small:
            u, v = find(b.n), find(b.m)
            parent.setdefault(v, v)
            parent[u] = v
    groups = {}
    for x in parent:
        groups.setdefault(find(x), set()).add(x)
    return {x: frozenset(g) for g in groups.values() for x in g}


def _residues(mol, cuts, units, rank):
    key = rank.__getitem__
    cut_set = {frozenset(x) for x in cuts}
    c_ends = {c for c, _ in cuts}
    n_ends = {n for _, n in cuts}
    residues, crosslinks = [], []
    for atoms in _fragments(mol, cut_set, rank):
        sp3c = {x for x in atoms if mol.element_of(x) == 6 and mol.hybridization_of(x) == 1}
        entries = sorted((x for x in atoms if x in n_ends), key=key)
        exits = sorted((x for x in atoms if x in c_ends), key=key)
        t_n = [x for x in sorted(atoms, key=key) if mol.element_of(x) == 7 and x not in n_ends and
               mol.hybridization_of(x) == 1 and not any(_carbonyl(mol, y) for y in mol.neighbors_of(x))]
        t_c = [x for x in sorted(atoms, key=key) if x not in c_ends and _carbonyl(mol, x) and
               sum(1 for y in _heavy(mol, x) if y in atoms) == 3]
        candidates = []
        for e in entries + t_n:
            for x in exits + t_c:
                p = _bfs_path(mol, e, x, sp3c, rank)
                if p and 3 <= len(p) <= 5:
                    candidates.append(((e not in n_ends) + (x not in c_ends), len(p), [rank[y] for y in p], p))
        candidates.sort()
        cores, used = [], set()
        for *_, p in candidates:
            if not used.intersection(p):
                cores.append(p)
                used.update(p)
        if not cores:
            ent = entries[0] if entries else None
            ext = exits[0] if exits else None
            if ent is not None and ext is not None:
                path = _bfs_path(mol, ent, ext, atoms, rank)
            elif ent is not None:
                path = _farthest_path(mol, ent, atoms, rank)
            elif ext is not None:
                path = _farthest_path(mol, ext, atoms, rank)[::-1]
            else:
                path = []
            residues.append(Residue(frozenset(atoms), 'other', ent, ext, tuple(path or ())))
            continue
        owner, queue = {}, deque()
        for i, p in enumerate(cores):
            for x in p:
                owner[x] = i
                queue.append(x)
        while queue:
            x = queue.popleft()
            for y in _neighbors(mol, x, rank):
                if y in atoms and y not in owner:
                    # a unit is claimed whole: links land on single bonds
                    for z in sorted(units.get(y, (y,)), key=key):
                        if z in atoms and z not in owner:
                            owner[z] = owner[x]
                            queue.append(z)
        _off_double_bonds(mol, atoms, owner, used, units, rank)
        for i, p in enumerate(cores):
            residues.append(Residue(frozenset(x for x in atoms if owner[x] == i), _KINDS[len(p)], p[0], p[-1],
                                    tuple(p)))
        links = [tuple(sorted((b.n, b.m), key=key)) for b in mol.bonds() if b.n in atoms and b.m in atoms and
                 owner[b.n] != owner[b.m]]
        crosslinks.extend(sorted(links, key=lambda t: (rank[t[0]], rank[t[1]])))
    return residues, crosslinks


def _to_backbone(mol, members, backbone, link, rank):
    """`link` moved to start on the backbone: of its two ends off the backbone, each side's bridge, the atoms of its
    group reached from it without the backbone, holding no atom of a ring of up to 8 and meeting the backbone at one
    atom, the one meeting it at the lower rank goes whole to the group across, and the link then joins that backbone
    atom to the bridge: the side chain hangs from one residue, the chord leaves the other's backbone."""
    a, b = link
    if a in backbone or b in backbone:
        return link
    small = {y for r in mol.sssr if len(r) <= 8 for y in r}
    group = {x: i for i, m in enumerate(members) for x in m}
    sides = []
    for x, o in ((a, b), (b, a)):
        own = members[group[x]]
        seen, stack, ends = {x}, [x], set()
        while stack:
            z = stack.pop()
            for y in _heavy(mol, z):
                if y in backbone:
                    ends.add((y, z))
                elif y in own and y not in seen and frozenset((z, y)) != frozenset(link):
                    seen.add(y)
                    stack.append(y)
        if len(ends) == 1 and not seen & small:
            (end, first), = ends
            sides.append((rank[end], end, first, seen, group[o]))
    if not sides:
        return link
    _, end, first, part, to = min(sides, key=lambda t: t[0])
    members[group[first]] -= part
    members[to] |= part
    return end, first


def _off_double_bonds(mol, atoms, owner, cores, units, rank):
    """Move each link between two owners a bond off an acyclic all-carbon unit, a C=C: the plain chain atom across
    the link, of two heavy neighbours and in no core, goes to the unit's owner, while one is there to move."""
    def plain(x):
        return x not in cores and x not in units and sum(1 for y in _heavy(mol, x) if y in atoms) == 2

    def chain_unit(x):
        return x in units and all(mol.element_of(y) == 6 for y in units[x]) and \
            all(mol.order_of(y, z) != 1 for y in units[x] for z in mol.neighbors_of(y) if z in units[x])
    for _ in range(len(atoms)):
        moved = False
        for b in sorted(mol.bonds(), key=lambda b: sorted((rank[b.n], rank[b.m]))):
            for u, x in ((b.n, b.m), (b.m, b.n)):
                if u in atoms and x in atoms and owner[u] != owner[x] and chain_unit(u) and plain(x):
                    owner[x] = owner[u]
                    moved = True
                    break
            if moved:
                break
        if not moved:
            return


def _chains(residues, nxt, prv, rank):
    """Linear chains, then cycles rotated to start at the residue whose entry ranks lowest; most real residues
    first."""
    seen, chains = set(), []
    for i in range(len(residues)):
        if i in seen or i in prv:
            continue
        chain = [i]
        while chain[-1] in nxt:
            chain.append(nxt[chain[-1]])
        seen.update(chain)
        chains.append((chain, False))
    for i in range(len(residues)):
        if i in seen:
            continue
        chain = [i]
        while nxt[chain[-1]] != i:
            chain.append(nxt[chain[-1]])
        seen.update(chain)
        k = min(range(len(chain)), key=lambda q: rank[residues[chain[q]].entry])
        chains.append((chain[k:] + chain[:k], True))
    chains.sort(key=lambda t: -sum(residues[i].kind != 'other' for i in t[0]))
    return chains


def _oriented(mol, residues, members, chain, prev, rank):
    """The chain's residue paths concatenated so that consecutive atoms are bonded; `prev` precedes it."""
    out = []
    for q, j in enumerate(chain):
        r = residues[j]
        p = list(r.path) or [r.entry if r.entry is not None else r.exit]
        last = out[-1] if out else prev
        if last is not None:
            nb = mol.neighbors_of(last)
            if p[-1] in nb and p[0] not in nb:
                p = p[::-1]
            elif p[0] not in nb:
                # entered through a side chain: run to the end bonded to the next residue
                start = next((x for x in sorted(members[j], key=rank.__getitem__) if x in nb), p[0])
                ahead = residues[chain[q + 1]].atoms if q + 1 < len(chain) else ()
                end = next((x for x in (r.exit, r.entry) if x is not None and x != start and
                            any(y in ahead for y in mol.neighbors_of(x))), None)
                p = (end is not None and _bfs_path(mol, start, end, r.atoms, rank)) or \
                    _farthest_path(mol, start, r.atoms, rank)
        out.extend(p)
    return out


def _symbols(mol, residues, library):
    variants = []
    for symbol, entry in library.items():
        for o in entry:
            if entry.element_of(o) != 8 or len(_heavy(entry, o)) != 1:
                continue
            c = _heavy(entry, o)[0]
            if entry.order_of(c, o) == 1 and _carbonyl(entry, c):
                v = entry.copy()
                with v.edit() as e:
                    e.delete_atom(o)
                variants.append((symbol, v))
    out = []
    for r in residues:
        sub = mol.substructure(sorted(r.atoms))
        heavy = sum(mol.element_of(x) != 1 for x in r.atoms)
        bonds = sum(1 for _ in sub.bonds())
        symbol = next((s for s, v in variants if len(v) == heavy and sum(1 for _ in v.bonds()) == bonds
                       and v.is_substructure(sub)), None)
        out.append(Residue(r.atoms, r.kind, r.entry, r.exit, r.path, symbol))
    return out
