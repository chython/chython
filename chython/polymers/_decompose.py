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
"""A molecule into the monomer graph.

1. **Cut** every single bond whose two ends are a pair of library attachment atoms that could join
   (:attr:`MonomerLibrary.cut_scopes`, one pair set per polymer type); the pieces are what is left.
2. **Units** are connected unions of up to four adjacent pieces; each is written by `detached_smiles`
   with a marker per boundary bond -- so its configuration is as drawn -- and looked up by key.
3. **Choose** the covering of all pieces that leaves the fewest atoms unmatched, then uses the fewest
   units: asparagine is one monomer, not an acid and an amide.  Adjacent unmatched pieces are one inline
   SMILES monomer, logged.
4. A unit or region whose constitution is a monomer but whose stereo is not raises
   :class:`StereoMismatch`.

Every order -- pieces, units, the search, the nodes -- follows ``canonical_order``, so the graph is a
function of the molecule and not of the order its atoms were written in.
"""

from collections import Counter

from ._errors import StereoMismatch
from ._graph import Edge, Node, PolymerGraph
from ._library import markers_for_closures, piece_keys, signature
from ._shape import Shape
from ..core import LogRecord, LOST, smiles


__all__ = ['decompose', 'MAX_MERGE']


#: The most pieces one unit is assembled from.
MAX_MERGE = 4
_PREFERENCE = {'PEPTIDE': 0, 'RNA': 1, 'DNA': 2, 'CHEM': 3}


def decompose(mol, library, log):
    """The monomer graph of `mol`, nodes in canonical order.  Raises :class:`StereoMismatch`.

    Each cut scope of the library (:attr:`MonomerLibrary.cut_scopes`) is tried in turn, PEPTIDE first, each
    matching only its own polymer types and CHEM.  The first covering that leaves no atom unmatched wins;
    else the single-type scope with the least unmatched atoms plus units, an earlier scope winning a tie.
    Atoms plus units, because another type's small linkers cover a molecule by fragmenting it; the
    all-types scope wins only by covering everything, for the same reason.
    """
    work = mol.copy()
    if work.has_stereo_groups:
        log.append(LogRecord('polymers:enhanced-stereo-not-written', (),
                             'unsupported: HELM has no enhanced stereo; the stereo groups are not written', LOST))
        work.clean_stereo_groups()
    work.thiele()
    rank = work.canonical_order()
    shape = Shape(work)
    cap_cuts = _cap_cuts(work, rank, library, shape)
    phosphorus = any(work.atom(n).atomic_symbol == 'P' for n in work.atom_numbers)
    best = None
    for scope, pairs in library.cut_scopes:
        if scope != 'PEPTIDE' and best is not None and not phosphorus:
            break       # every nucleotide carries phosphorus; without it the other scopes only fragment
        attempt = _attempt(work, rank, pairs, library, cap_cuts, _types(scope), shape)
        if not attempt[0][0]:
            best = attempt
            break
        if scope != 'ALL' and (best is None or sum(attempt[0]) < sum(best[0])):
            best = attempt
    _, nodes, pieces, cuts = best
    return _graph(work, nodes, pieces, cuts, rank, library, log)


def _types(scope):
    return None if scope == 'ALL' else frozenset((scope, 'CHEM'))


def _attempt(work, rank, pairs, library, cap_cuts, types, shape):
    """``((unmatched atoms, units), nodes, pieces, cuts)`` of one cut scope; lookups limited to `types`.

    Candidates are screened by :mod:`._shape` key, read off the molecule with nothing written; only the
    units a covering chooses are written out and checked exactly.  A chosen unit that fails is dropped
    and the covering rerun.  Raises :class:`StereoMismatch` when a unit matching a monomer by constitution
    only would cover atoms the final covering leaves unmatched: then a configuration, not the library, is
    what stands between that unit and its monomer.
    """
    cuts = list(cap_cuts)
    small = _small_ring_bonds(work)
    for bond in work.bonds():
        n, m = bond.n, bond.m
        if bond.order == 1 and frozenset((n, m)) not in small and not _alkyl(work, n, m) and \
                (signature(work, n, m), signature(work, m, n)) in pairs:
            pair = (n, m) if rank[n] < rank[m] else (m, n)
            if pair not in cap_cuts:
                cuts.append(pair)
    cuts.sort(key=lambda b: (rank[b[0]], rank[b[1]]))

    pieces = _pieces(work, cuts, rank)
    owner = {atom: i for i, piece in enumerate(pieces) for atom in piece}
    neighbours = [set() for _ in pieces]
    for n, m in cuts:
        if owner[n] != owner[m]:
            neighbours[owner[n]].add(owner[m])
            neighbours[owner[m]].add(owner[n])

    units = _Units(work, pieces, cuts, owner, library, types)
    parts = [shape.piece(p) for p in pieces]
    incident = [[] for _ in pieces]
    for bond in cuts:
        a, b = owner[bond[0]], owner[bond[1]]
        if a != b:
            incident[a].append(bond)
            incident[b].append(bond)
    candidates = set()
    for group in _connected_groups(neighbours):
        internal, boundary = [], []
        for p in group:
            for n, m in incident[p]:
                inside_n, inside_m = owner[n] in group, owner[m] in group
                if inside_n and inside_m:
                    if p == min(owner[n], owner[m]):
                        internal.append((n, m))
                else:
                    boundary.append((n, m) if inside_n else (m, n))
        key = shape.join([parts[p] for p in group], internal, boundary)
        if _scoped_entries(library.by_shape(key), types):
            candidates.add(group)

    verified, loose = {}, {}                        # exact matches; constitution-only matches
    while True:
        chosen, unmatched = _cover(pieces, candidates, neighbours)
        failed = False
        for group in chosen:
            result = units.evaluate(group)
            if result.matches:
                verified[group] = result
                continue
            candidates.discard(group)
            failed = True
            if result.constitution:
                loose[group] = result
        if not failed:
            break
    blamed = [r for g, r in loose.items() if g <= unmatched]
    if blamed:
        raise StereoMismatch('; '.join(f'{r.text} is monomer {_names(r.constitution)} by constitution and '
                                       f'not by stereo' for r in blamed))

    nodes = [(group, verified[group]) for group in chosen]
    for region in _regions(unmatched, neighbours):
        nodes.append((region, units.evaluate(region, inline=True)))
    nodes.sort(key=lambda item: min(rank[a] for p in item[0] for a in pieces[p]))
    return (_missed(pieces, unmatched), len(nodes)), nodes, pieces, cuts


#: The largest ring whose bonds are never a monomer junction; a macrocycle's are.
SMALL_RING = 8


def _small_ring_bonds(work):
    """The bonds of every ring of at most SMALL_RING atoms: inside a monomer, never between two."""
    out = set()
    for ring in work.rings:
        if len(ring) <= SMALL_RING:
            out.update(frozenset(pair) for pair in zip(ring, ring[1:] + ring[:1]))
    return out


def _alkyl(work, n, m):
    """Is `n`-`m` a bond between two saturated carbons?  Never cut, by a pair or a cap: a monomer junction has a
    heteroatom, a carbonyl or an aromatic atom at one end, and a library with alkyl attachment sites (an R3
    on an alpha carbon, an alkyl cap) would otherwise cut every side chain into its methylenes."""
    for atom in (n, m):
        if work.atom(atom).atomic_symbol != 'C':
            return False
        if any(work.order_of(atom, x) != 1 for x in work.neighbors_of(atom)):
            return False
    return True


def _missed(pieces, unmatched):
    return sum(len(pieces[p]) for p in unmatched)


#: The most heavy atoms a terminal cap may hold.
MAX_CAP = 12


def _cap_cuts(work, rank, library, shape):
    """Bonds whose far side is a terminal cap.

    The far side is ring-free, at most MAX_CAP atoms, and has the :mod:`._shape` key of a one-R library
    monomer with a marker where the bond was; and the near end is an attachment the cap can join -- an R1
    cap meets an R2 site, an R2 cap an R1 site, of the cap's polymer type.  Without the second test an
    N-terminal amine reads as an amide cap and an alanine's methyl as a methyl cap.
    """
    out = set()
    formulas = library.cap_formulas
    sites = library.cap_sites
    for bond in work.bonds():
        if bond.order != 1 or bond.in_ring or _alkyl(work, bond.n, bond.m):
            continue
        for near, far in ((bond.n, bond.m), (bond.m, bond.n)):
            if signature(work, near, far) not in sites:
                continue
            side = _side(work, far, near, MAX_CAP)
            if side is None or _formula(work, side, 1) not in formulas:
                continue
            # a screen: a wrong cap cut only leaves one more piece, which a unit merges back
            key = shape.join([shape.piece(tuple(side))], (), ((far, near),))
            end = signature(work, near, far)
            for monomer, used in library.caps_by_shape(key):
                (label,) = used
                partner = 'R2' if label == 1 else 'R1' if label == 2 else 'R3+'
                if end in library.link_ends(monomer.polymer_type, partner):
                    out.add((near, far) if rank[near] < rank[far] else (far, near))
                    break
    return out


def _side(work, start, across, limit):
    """The atoms reachable from `start` without crossing to `across`, or ``None`` past `limit` atoms."""
    seen, stack = {start}, [start]
    while stack:
        atom = stack.pop()
        for n in work.neighbors_of(atom):
            if n not in seen and not (atom == start and n == across):
                seen.add(n)
                if len(seen) > limit:
                    return None
                stack.append(n)
    return seen


def _pieces(work, cuts, rank):
    """The components left when every cut bond is removed, each a sorted tuple, in canonical order."""
    cut = {frozenset(b) for b in cuts}
    seen, out = set(), []
    for start in sorted(work.atom_numbers, key=rank.get):
        if start in seen:
            continue
        stack, piece = [start], []
        seen.add(start)
        while stack:
            atom = stack.pop()
            piece.append(atom)
            for n in work.neighbors_of(atom):
                if n not in seen and frozenset((atom, n)) not in cut:
                    seen.add(n)
                    stack.append(n)
        out.append(tuple(sorted(piece, key=rank.get)))
    return out


def _connected_groups(neighbours):
    """Every connected set of up to MAX_MERGE pieces, each once."""
    out, frontier = set(), [frozenset((i,)) for i in range(len(neighbours))]
    out.update(frontier)
    for _ in range(MAX_MERGE - 1):
        grown = set()
        for group in frontier:
            for p in group:
                for q in neighbours[p]:
                    if q not in group:
                        bigger = group | {q}
                        if bigger not in out:
                            grown.add(bigger)
        out.update(grown)
        frontier = grown
    return sorted(out, key=lambda g: (len(g), sorted(g)))


class _Evaluation:
    __slots__ = ('text', 'markers', 'matches', 'constitution')

    def __init__(self, text, markers, matches, constitution):
        self.text = text                    # SMILES with '*' markers, in boundary-bond order
        self.markers = markers              # boundary bonds (inside atom, outside atom), in marker order
        self.matches = matches              # ((monomer, used), ...) by full key
        self.constitution = constitution    # ((monomer, used), ...) by constitution key only


class _Units:
    """Unit texts and lookups, each union of pieces evaluated once."""
    __slots__ = ('work', 'pieces', 'cuts', 'owner', 'library', 'cache', 'formulas', 'types')

    def __init__(self, work, pieces, cuts, owner, library, types):
        self.work, self.pieces, self.cuts, self.owner, self.library = work, pieces, cuts, owner, library
        self.cache = {}
        self.formulas = library.core_formulas
        self.types = types

    def evaluate(self, group, inline=False):
        result = self.cache.get((group, inline))
        if result is not None:
            return result
        atoms = {a for p in group for a in self.pieces[p]}
        boundary = [(n, m) if n in atoms else (m, n) for n, m in self.cuts if (n in atoms) != (m in atoms)]
        if not inline:
            formula = _formula(self.work, atoms, len(boundary))
            if formula not in self.formulas:
                result = self.cache[(group, inline)] = _Evaluation('', boundary, (), ())
                return result
        detached = self.work.detached_smiles({10 + i: bond for i, bond in enumerate(boundary)})
        text = detached.text
        if not boundary and '.' in text:
            text = _component(text, detached.order, atoms)
        if detached.lost:
            raise StereoMismatch(f'{text}: the configuration at {detached.lost} cannot be written for this piece')
        text, order = markers_for_closures(text)
        markers = [boundary[number - 10] for number in order]
        key, flat = piece_keys(text)
        result = self.cache[(group, inline)] = _Evaluation(text, markers, self._scoped(self.library.by_core(key)),
                                                           self._scoped(self.library.by_constitution(flat)))
        return result

    def _scoped(self, entries):
        return _scoped_entries(entries, self.types)


def _scoped_entries(entries, types):
    return entries if types is None else tuple(e for e in entries if e[0].polymer_type in types)


def _formula(work, atoms, markers):
    counts = Counter(work.atom(a).atomic_symbol for a in atoms)
    counts['R'] += markers
    return tuple(sorted(counts.items()))


def _component(text, order, atoms):
    """The one ``.`` segment of `text` that writes `atoms`."""
    position = 0
    for segment in text.split('.'):
        size = smiles(segment).atom_count
        if order[position] in atoms:
            return segment
        position += size
    return text


def _cover(pieces, matched, neighbours):
    """``(chosen units, unmatched pieces)``: fewest atoms unmatched, then fewest units, first in canonical order.

    Pieces interact only through a unit that holds both, so the search runs per component of that
    relation -- about one residue each -- and a piece no unit holds is unmatched without a search.
    """
    owner = list(range(len(pieces)))

    def find(i):
        while owner[i] != i:
            owner[i] = owner[owner[i]]
            i = owner[i]
        return i

    for group in matched:
        first, *rest = sorted(group)
        for other in rest:
            owner[find(other)] = find(first)
    components = {}
    for group in matched:
        components.setdefault(find(min(group)), []).append(group)
    held = {p for group in matched for p in group}
    chosen, unmatched = [], {p for p in range(len(pieces)) if p not in held}
    for root in sorted(components):
        groups = components[root]
        members = sorted({p for g in groups for p in g})
        best = _search(pieces, members, groups)
        chosen.extend(best[1])
        unmatched |= best[2]
    chosen.sort(key=min)
    return chosen, unmatched


def _search(pieces, members, groups):
    """The best covering of `members` by `groups`, leaving pieces unmatched where none fits."""
    containing = {p: [] for p in members}
    for group in sorted(groups, key=lambda g: (-len(g), sorted(g))):
        containing[min(group)].append(group)
    best = [None]

    def dfs(position, covered, chosen, unmatched, missed):
        if best[0] is not None and (missed, len(chosen)) >= best[0][0]:
            return
        while position < len(members) and members[position] in covered:
            position += 1
        if position == len(members):
            best[0] = ((missed, len(chosen)), list(chosen), set(unmatched))
            return
        free = members[position]
        for group in containing[free]:
            if not (group & covered):
                chosen.append(group)
                dfs(position + 1, covered | group, chosen, unmatched, missed)
                chosen.pop()
        dfs(position + 1, covered | {free}, chosen, unmatched | {free}, missed + len(pieces[free]))

    dfs(0, frozenset(), [], frozenset(), 0)
    return best[0]


def _regions(unmatched, neighbours):
    """Connected groups of unmatched pieces."""
    out, seen = [], set()
    for start in sorted(unmatched):
        if start in seen:
            continue
        region, stack = {start}, [start]
        seen.add(start)
        while stack:
            p = stack.pop()
            for q in neighbours[p]:
                if q in unmatched and q not in seen:
                    seen.add(q)
                    region.add(q)
                    stack.append(q)
        out.append(frozenset(region))
    return out


def _names(entries):
    return ', '.join(sorted({m.symbol for m, _ in entries}))


def _graph(work, nodes, pieces, cuts, rank, library, log):
    graph = PolymerGraph()
    unit_of = {}
    marker_label = {}                               # (node, boundary bond) -> R number
    for index, (group, result) in enumerate(nodes):
        for p in group:
            for atom in pieces[p]:
                unit_of[atom] = index
        if result.matches:
            monomer, used = _pick(result.matches)
            if len({m.symbol for m, _ in result.matches}) > 1:
                log.append(LogRecord('polymers:monomer-ambiguous', (),
                                     f'{result.text} is monomers {_names(result.matches)}; {monomer.symbol} is '
                                     f'written'))
            graph.nodes.append(Node('', 0, monomer.polymer_type, monomer.symbol, False))
            for bond, label in _labels(work, result, monomer, used, rank, library, log):
                marker_label[(index, bond)] = label
        else:
            parts = result.text.split('*')
            symbol = parts[0] + ''.join(f'[*:{k}]{rest}' for k, rest in enumerate(parts[1:], 1))
            for k, bond in enumerate(result.markers, 1):
                marker_label[(index, bond)] = k
            graph.nodes.append(Node('', 0, 'CHEM', symbol, True))
            log.append(LogRecord('polymers:inline-monomer', (),
                                 f'{symbol} matches no library monomer; written as an inline SMILES monomer'))
    for index, node in enumerate(graph.nodes):
        if node.inline and any(graph.nodes[unit_of[o]].polymer_type == 'PEPTIDE'
                               for (i, (_, o)) in marker_label if i == index):
            graph.nodes[index] = node._replace(polymer_type='PEPTIDE')
    for n, m in cuts:
        a, b = unit_of[n], unit_of[m]
        if a != b:
            graph.edges.append(Edge(a, marker_label[(a, (n, m))], b, marker_label[(b, (m, n))]))
    return graph


def _pick(matches):
    """The monomer to write: PEPTIDE before RNA, DNA, CHEM; then the latest library's."""
    ordered = sorted(enumerate(matches), key=lambda x: (_PREFERENCE.get(x[1][0].polymer_type, 9), -x[0]))
    return ordered[0][1]


def _labels(work, result, monomer, used, rank, library, log):
    """``[(boundary bond, R number)]`` for a matched unit.

    The unit is mapped onto the monomer's labelled core.  Markers the core holds in one symmetry class are
    interchangeable, so their R numbers go in canonical order of the atom each attaches to.
    """
    unit = smiles(result.text)
    unit.implicify_hydrogens()
    unit.thiele()
    core, names = library.labelled_core(monomer, used)
    mapping = unit.isomorphism(core)
    classes = core.atoms_order
    unit_markers = [n for n in unit.atom_numbers if unit.atom(n).is_r]
    by_class = {}
    for marker, bond in zip(unit_markers, result.markers):
        target = mapping[marker]
        by_class.setdefault(classes[target], []).append((bond, names[target]))
    out = []
    for members in by_class.values():
        if len(members) > 1:
            log.append(LogRecord('polymers:symmetric-attachment', (),
                                 f'{monomer.symbol}: R{", R".join(str(x) for _, x in members)} are symmetric; '
                                 f'numbered in canonical order of their attachment atoms'))
        labels = sorted(label for _, label in members)
        bonds = sorted((bond for bond, _ in members), key=lambda b: rank[b[0]])
        out.extend(zip(bonds, labels))
    return out
