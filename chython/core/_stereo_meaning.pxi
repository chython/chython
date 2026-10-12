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
# THE MEANING OF A STEREO LABELLING, AND ITS LEAST SPELLING.  `_meaning_canonical` is the one entry.
#
# A labelling gives each unit of the searched components one label:
#
# | code       | label                                                                              |
# | ---------- | ---------------------------------------------------------------------------------- |
# | 0          | cleared: the unit is unconfigured, both configurations are open                    |
# | 1          | unlabelled sign: configured as drawn, which of ABS, OR or AND applies is not known |
# | 2          | ABS                                                                                |
# | 3 + 2g     | AND group g                                                                        |
# | 4 + 2g     | OR group g                                                                         |
#
# MEANING.  A world is one choice of every OR group and cleared unit, as the set of pure molecules over
# the AND phases; a reading's meaning is `(drawn worlds, every world)`, the drawn ones having no OR group
# inverted, since an OR group keeps its sign.  The unlabelled signs of one component are read together as
# ABS, as one OR group or as one AND group, and the meaning of a labelling is the set of its readings'.
# Before reading, an unlabelled or ABS sign whose inversion leaves the ABS reading unchanged states
# nothing and is cleared, to a fixpoint: the C3 sign of pentane-2,3,4-triol between two open ends.
#
# THE PURE TABLE.  Every meaning is read off `T[mask]`, the class of the molecule with unit i at parity
# 1 + bit i and no label, so a search weighs a candidate without building it.
#
# THE LEAST SPELLING of a meaning is the least by `_labelling_measure_codes` -- configured units, AND
# and OR members, groups, explicit labels, explicit labels off a central unit -- then by identity.  A
# unit is central when no other unit of its cluster is further from the rest; a cluster is the units
# whose configurations are coupled in `T`.  So a label sits on the spiro atom of Fecht acid, on the C=C
# of a bi(cyclohexylidene) and on both spiro atoms of a dispiro chain: without it the terminal signs do
# not state a configuration.  The racemate of Fecht acid is its spiro sign alone in an AND group, the
# terminals open.

DEF LABELLING_UNIT_LIMIT = 5            # units searched per block; above it `_relabelled` answers
DEF LABELLING_CANDIDATE_LIMIT = 65536   # labellings weighed per block; above it the same
DEF READING_GROUP = 100                 # group index of a component's unlabelled signs read as a group

cdef dict _LABEL_SHAPES = {}            # unit count -> every labelling up to group ids, as code tuples


cdef bytes _unit_codes(Structure view, dict cache):
    """Per slot: 1 at a unit stereogenic with every parity cleared, 2 at one some parities may make
    stereogenic, 0 elsewhere; None on a truncated perception.  Constitution alone decides it, so it
    is kept in `cache` under `b'M'`, where b'' stands for None."""
    cdef bytes codes = cache.get(b'M')
    cdef uint32_t n = view.header.atom_count, k, count
    cdef uint8_t unnamed
    cdef stereo_unit_t *units
    cdef unsigned char *fill
    cdef unsigned char best
    cdef tuple members
    cdef object i
    cdef Structure v
    if codes is not None:
        return codes if codes else None
    v = structure_clone(view)
    memset(structure_parities(v), 0, <size_t> n)
    structure_invalidate_stereo_units(v)
    refresh_parity_features(v)
    ensure_stereo_units(v)
    if structure_stereo_truncated(v):
        cache[b'M'] = b''
        return None
    fill = <unsigned char *> PyMem_Malloc(<size_t> n + 1)
    if fill is NULL:
        raise MemoryError()
    try:
        memset(fill, 0, <size_t> n)
        count = structure_stereo_unit_count(v)
        units = structure_stereo_units(v)
        for k in range(count):
            unnamed = units[k].spare >> SU_UNNAMED_SHIFT
            if units[k].spare & SU_STEREOGENIC:
                fill[units[k].anchor] = 1
            elif units[k].kind != SU_TETRA or not unnamed & (unnamed - 1):
                fill[units[k].anchor] = 2
        # a coupled set is one choice, named on its lowest member: the others are implied by it
        for members in _coupled_sets(_coupled_bridgeheads(v, cache)).values():
            best = 0
            for i in members:
                if fill[<uint32_t> i] and (not best or fill[<uint32_t> i] < best):
                    best = fill[<uint32_t> i]
                fill[<uint32_t> i] = 0
            fill[<uint32_t> members[0]] = best
        codes = (<char *> fill)[:n]
    finally:
        PyMem_Free(fill)
    cache[b'M'] = codes
    return codes


cdef tuple _reading(list T, uint32_t ps, tuple labels, uint32_t k):
    """`(drawn worlds, every world)` of `labels` at parity mask `ps`, an unlabelled sign read as
    drawn."""
    cdef dict members = {}
    cdef list ands = [0], ors = [0], unknown = [0], more
    cdef uint32_t i, fo, m0, x
    cdef set drawn = set(), every = set()
    cdef object label, w
    for i in range(k):
        label = labels[i]
        if label == 0:
            unknown = unknown + [x | (1u << i) for x in unknown]
            ps &= ~(1u << i)
        elif label >= 3:
            members[label] = members.get(label, 0) | (1u << i)
    for label, x in members.items():
        more = ands if label & 1 else ors
        more.extend([m0 ^ x for m0 in more])
    for i, fo in enumerate(ors):
        for x in unknown:
            m0 = (ps ^ fo) | x
            w = frozenset([T[m0 ^ a] for a in ands])
            every.add(w)
            if i == 0:
                drawn.add(w)
    return frozenset(drawn), frozenset(every)


cdef tuple _inert_cleared(list T, uint32_t ps, tuple labels, uint32_t k):
    """`labels` with every unlabelled or ABS sign cleared whose inversion leaves the ABS reading
    unchanged, to a fixpoint."""
    cdef list now = list(labels), inert
    cdef tuple ref, cur
    cdef uint32_t i
    while True:
        cur = tuple(now)
        ref = _reading(T, ps, cur, k)
        inert = [i for i in range(k) if (now[i] == 1 or now[i] == 2)
                 and _reading(T, ps ^ (1u << i), cur, k) == ref]
        if not inert:
            return cur
        for i in inert:
            now[i] = 0


cdef frozenset _readings(list T, uint32_t ps, tuple labels, uint32_t k, list comp):
    """The meaning of an inert-free `labels`: its reading with each component's unlabelled signs as
    ABS, as one OR group and as one AND group."""
    cdef list comps = sorted({comp[i] for i in range(k) if labels[i] == 1}), now
    cdef set out = set()
    cdef uint32_t i, j, choice, pick, total = 1
    for j in range(len(comps)):
        total *= 3
    for choice in range(total):
        now = list(labels)
        pick = choice
        for j in range(len(comps)):
            for i in range(k):
                if labels[i] == 1 and comp[i] == comps[j]:
                    now[i] = 2 if pick % 3 == 0 else 2 + pick % 3 + 2 * (READING_GROUP + j)
            pick //= 3
        out.add(_reading(T, ps, tuple(now), k))
    return frozenset(out)


cdef frozenset _stated_clusters(list T, uint32_t ps, tuple labels, uint32_t k, list clusters):
    """The clusters (`_clusters`) whose labels state something: clearing every ABS or group label of the
    cluster moves the reading, unlabelled signs read as drawn.

    Unlabelled is unknown, so a cluster stating a label keeps one -- the canonical form may move it inside
    the cluster or respell its kind, never drop it.  On a configuration that is its own mirror image
    `|a:1,5|` reads as no label does, and cleared it is either diastereomer: the label is what separates the
    drawn one from the other, a cis 1,3-cyclobutane beside an OR centre elsewhere as much as alone.  Labels
    whose signs determine nothing -- beside an open end of Fecht acid -- clear without moving it."""
    cdef tuple ref = _reading(T, ps, labels, k)
    cdef list cleared
    cdef set out = set(), seen = set()
    cdef object c
    cdef uint32_t j
    for j in range(k):
        c = clusters[j]
        if labels[j] >= 2 and c not in seen:
            seen.add(c)
            cleared = [0 if labels[i] >= 2 and i in c else labels[i] for i in range(k)]
            if _reading(T, ps, tuple(cleared), k) != ref:
                out.add(c)
    return frozenset(out)


cdef list _label_shapes(uint32_t k):
    """Every labelling of `k` units, groups numbered in first-use order, each group of either kind."""
    cdef list out = _LABEL_SHAPES.get(k)
    cdef list stack, shapes = []
    cdef tuple cur
    cdef uint32_t groups, o, kinds
    if out is not None:
        return out
    stack = [((), 0)]
    while stack:
        cur, groups = stack.pop()
        if len(cur) == k:
            shapes.append((cur, groups))
            continue
        for o in range(3 + groups + 1):
            stack.append((cur + (o,), groups + (o == 3 + groups)))
    out = []
    for cur, groups in shapes:
        for kinds in range(1u << groups):
            out.append(tuple([o if o < 3 else 3 + 2 * (o - 3) + ((kinds >> (o - 3)) & 1) for o in cur]))
    _LABEL_SHAPES[k] = out
    return out


cdef tuple _labelling_measure_codes(tuple labels, list central):
    """`(configured, AND and OR members, groups, explicit labels, explicit labels off a central unit)`."""
    cdef uint32_t i, conf = 0, members = 0, explicit = 0, off = 0
    cdef set groups = set()
    for i in range(len(labels)):
        if labels[i] == 0:
            continue
        conf += 1
        if labels[i] >= 3:
            members += 1
            groups.add(labels[i])
        if labels[i] >= 2:
            explicit += 1
            if not central[i]:
                off += 1
    return conf, members, len(groups), explicit, off


cdef list _clusters(list T, uint32_t k):
    """Per unit, its cluster as the set of its members.  Units i and j are coupled when the effect of
    inverting one depends on the other's configuration, or when inverting both restates a molecule
    inverting one does not."""
    cdef list adj = [set() for _ in range(k)], out = [None] * k, stack
    cdef set members
    cdef frozenset frozen
    cdef uint32_t i, j, m, bi, bj, x, y
    for i in range(k):
        for j in range(k):
            if i == j:
                continue
            bi = 1u << i
            bj = 1u << j
            for m in range(1u << k):
                if ((T[m] == T[m ^ bi]) != (T[m ^ bj] == T[m ^ bi ^ bj])
                        or (T[m ^ bi ^ bj] == T[m] and T[m] != T[m ^ bi])):
                    adj[i].add(j)
                    adj[j].add(i)
                    break
    for i in range(k):
        if out[i] is not None:
            continue
        members = {i}
        stack = [i]
        while stack:
            x = stack.pop()
            for y in adj[x]:
                if y not in members:
                    members.add(y)
                    stack.append(y)
        frozen = frozenset(members)
        for x in members:
            out[x] = frozen
    return out


cdef list _central_units(list T, uint32_t k, list dist):
    """Per unit: is it central in its cluster (`_clusters`)?"""
    cdef list central = [False] * k
    cdef dict ecc
    cdef object members
    cdef uint32_t m
    for members in set(_clusters(T, k)):
        ecc = {x: max([dist[x][y] for y in members]) for x in members}
        m = min(ecc.values())
        for x in members:
            central[x] = ecc[x] == m
    return central


cdef list _unit_distances(Structure view, list anchors):
    """Bond distances between units, the least over their owner atoms.  `view` has its unit table."""
    cdef uint32_t n = view.header.atom_count, k = len(anchors), i, j, a = 0, b = 0, x, e
    cdef atom_t *atoms = view.atoms()
    cdef uint32_t *ptr = csr_ptr(view)
    cdef halfedge_t *edges = csr_edges(view)
    cdef stereo_unit_t *u
    cdef list owners = [], dist = [], front, nxt, row
    cdef dict seen
    cdef int nown
    for i in range(k):
        u = stereo_unit_of(view, <uint32_t> anchors[i])
        nown = stereo_unit_owners(atoms, ptr, edges, u, &a, &b)
        if nown == 1:
            owners.append((a,))
        else:
            owners.append((a, b))
    for i in range(k):
        seen = {x: 0 for x in owners[i]}
        front = list(owners[i])
        while front:
            nxt = []
            for x in front:
                for e in range(ptr[x], ptr[x + 1]):
                    if edges[e].to not in seen:
                        seen[edges[e].to] = seen[x] + 1
                        nxt.append(edges[e].to)
            front = nxt
        row = []
        for j in range(k):
            row.append(min([seen.get(x, n) for x in owners[j]]))
        dist.append(row)
    return dist


cdef Structure _labelled_clone(Structure base, list anchors, uint32_t ps, tuple labels, list free_ids):
    """`base` with unit j at parity 1 + bit j of `ps` under label j; group g numbered by `free_ids`
    per kind in first-use order."""
    cdef Structure c = structure_clone(base)
    _label_into(c, anchors, ps, labels, free_ids)
    structure_invalidate_stereo_units(c)
    refresh_parity_features(c)
    return c


cdef int _label_into(Structure c, list anchors, uint32_t ps, tuple labels, list free_ids) except -1:
    """`_labelled_clone` in place, the stereo units left for the caller to refresh."""
    cdef uint8_t *par = structure_parities(c)
    cdef uint8_t *sg = structure_stereo_groups(c)
    cdef dict number = {}
    cdef uint32_t j, a, kind
    cdef object label
    for j in range(len(anchors)):
        a = anchors[j]
        label = labels[j]
        if label == 0:
            par[a] = 0
            continue
        par[a] = 1 + ((ps >> j) & 1)
        if label == 2:
            sg[a] = sg_pack(SG_KIND_ABS, 0)
        elif label >= 3:
            kind = SG_KIND_AND if label & 1 else SG_KIND_ABS + 1
            if label not in number:
                number[label] = free_ids[kind].pop(0)
            sg[a] = sg_pack(<uint8_t> kind, <uint8_t> number[label])
    return 0


cdef bint _may_differ(Structure structure, list inert=None) except -1:
    """False when the view's answer is the axiom's: no sign carries a label and no stereogenic unit is
    open, a member of `inert` counting as none.  Every reading of unlabelled signs inverts them
    together, which keeps an open unit that is not stereogenic in the drawing not stereogenic, so each
    sign is read at fixed neighbours.  A cheap screen before `_unit_codes` perceives again.

    `inert` holds the coupled bridgeheads whose inversion `_bridgeheads_filled` found to restate the molecule."""
    cdef uint32_t n = structure.header.atom_count, i, j, count
    cdef uint8_t *par = structure_parities(structure)
    cdef uint8_t *sg
    cdef stereo_unit_t *units
    cdef bint signed = False
    if structure_has(structure, SEG_STEREO_GROUPS):
        sg = structure_stereo_groups(structure)
        for i in range(n):
            if sg[i] and par[i]:
                return True
    for i in range(n):
        if par[i]:
            signed = True
            break
    if not signed:
        return False
    ensure_stereo_units(structure)
    count = structure_stereo_unit_count(structure)
    units = structure_stereo_units(structure)
    for j in range(count):
        if par[units[j].anchor] or inert and units[j].anchor in inert:
            continue
        if units[j].spare & SU_STEREOGENIC:
            return True
    return False


cdef bint _anchored_swap(Structure v, Structure view, uint32_t *orbits, dict anchors, list rest,
                         set scope) except -1:
    """Add to `scope` each component of `rest` where an automorphism fixing every anchor and owner
    swaps two neighbours of one of them, so inverts units together: the ring C4 and the exocyclic
    C=C of 4-methylcyclohexylideneacetic acid.  `v` is the parity-cleared clone, `orbits` a scratch
    array of its atom count.  True when the search was truncated."""
    cdef uint32_t n = v.header.atom_count, i, x, e, a = 0, b = 0, flags = 0
    cdef Py_ssize_t classes
    cdef atom_t *atoms = view.atoms()
    cdef uint32_t *ptr = csr_ptr(view)
    cdef halfedge_t *edges = csr_edges(view)
    cdef stereo_unit_t *u
    cdef uint32_t *seed = <uint32_t *> PyMem_Malloc(<size_t> n * sizeof(uint32_t))
    cdef dict fixed = {}
    cdef set seen
    cdef object c
    if seed is NULL:
        raise MemoryError()
    try:
        classes = compute_atoms_order(v, seed, NULL)
        if classes < 0:
            raise MemoryError()
        for c in rest:
            fixed[c] = []
            for i in anchors[c]:
                u = stereo_unit_of(view, i)
                fixed[c].append(i)
                if stereo_unit_owners(atoms, ptr, edges, u, &a, &b) == 1:
                    fixed[c].append(a)
                else:
                    fixed[c] += [a, b]
            for x in fixed[c]:
                seed[x] = <uint32_t> (classes + 1 + x)
        mol_automorphisms(v, seed, orbits, &flags)
        if flags & CANON_BUDGET_EXCEEDED:
            return True
        for c in rest:
            for x in fixed[c]:
                seen = set()
                for e in range(ptr[x], ptr[x + 1]):
                    if orbits[edges[e].to] in seen:
                        scope.add(c)
                        break
                    seen.add(orbits[edges[e].to])
                if c in scope:
                    break
    finally:
        PyMem_Free(seed)
    return False


cdef list _coupled_blocks(set scope, uint8_t *sg, uint8_t *par, list cof, uint32_t n):
    """The components of `scope`, each grown by every component a chain of AND or OR groups joins it
    to, as disjoint sets in ascending order of their least component.  An ABS byte joins nothing."""
    cdef dict joined = {}
    cdef list out = [], stack
    cdef set placed = set(), cl, members
    cdef uint32_t i
    cdef object c, x, y
    if sg is not NULL:
        for i in range(n):
            if sg[i] and par[i] and sg_kind(sg[i]) != SG_KIND_ABS:
                joined.setdefault(sg[i], set()).add(cof[i])
    for c in sorted(scope):
        if c in placed:
            continue
        cl = {c}
        stack = [c]
        while stack:
            x = stack.pop()
            for members in joined.values():
                if x in members:
                    for y in members:
                        if y not in cl:
                            cl.add(y)
                            stack.append(y)
        placed |= cl
        out.append(cl)
    return out


cdef object _block_least(Structure structure, Structure view, Structure base, set cl, list anchor,
                           list cof, list free_ids, dict pcache):
    """The least labelling of one block as `(anchors, parity mask, labels, clone, identity)`, the
    clone being `base` with the block labelled; None when no two of its labellings meet, so the view's
    labels stand; False over a limit, MEANING_ATOM_LIMIT included.  `base` has every block cleared.  A
    truncated identity raises `AutomorphismBudgetExceeded`."""
    cdef uint32_t n = view.header.atom_count, i, j, k = 0, mask, count, ps = 0
    cdef uint8_t *par = structure_parities(structure)
    cdef uint8_t *sg = NULL
    cdef uint8_t *bpar
    cdef stereo_unit_t *units
    cdef Structure v, best = None
    cdef bytes identity, best_id = None
    cdef dict classes = {}, groups = {}, levels = {}
    cdef set grow, seen
    cdef list T = None, target_labels = [], unit_comp, central, fixed_free, clusters
    cdef frozenset stated
    cdef tuple lab, ref, measure, best_lab = None
    cdef frozenset target
    cdef object label
    cdef uint32_t q, bits, tried = 0, fixed, best_q = 0
    if structure_has(structure, SEG_STEREO_GROUPS):
        sg = structure_stereo_groups(structure)

    if n > MEANING_ATOM_LIMIT:
        return False
    # the pure table, the units grown by each one some combination makes stereogenic
    while True:
        k = len(anchor)
        if k > LABELLING_UNIT_LIMIT:
            return False
        T = []
        classes = {}
        grow = set()
        for mask in range(1u << k):
            v = structure_clone(base)
            bpar = structure_parities(v)
            for j in range(k):
                bpar[<uint32_t> anchor[j]] = 1 + ((mask >> j) & 1)
            structure_invalidate_stereo_units(v)
            refresh_parity_features(v)
            ensure_stereo_units(v)
            count = structure_stereo_unit_count(v)
            units = structure_stereo_units(v)
            for j in range(count):
                i = units[j].anchor
                if units[j].spare & SU_STEREOGENIC and cof[i] in cl and i not in anchor:
                    grow.add(i)
            if grow:
                break
            identity = mol_identity_bytes(structure_stereogenic_view(v, pcache))
            T.append(classes.setdefault(identity, len(classes)))
        if not grow:
            break
        anchor = sorted(set(anchor) | grow)
    if len(classes) == 1u << k:
        return None                     # every combination is its own compound: no two labellings meet

    for j in range(k):
        i = anchor[j]
        if not par[i]:
            label = 0
        elif sg is NULL or not sg[i]:
            label = 1
        elif sg_kind(sg[i]) == SG_KIND_ABS:
            label = 2
        else:
            label = (3 if sg_kind(sg[i]) == SG_KIND_AND else 4) + 2 * groups.setdefault(sg[i], len(groups))
        target_labels.append(label)
        if par[i] == 2:
            ps |= 1u << j
    unit_comp = [cof[<uint32_t> i] for i in anchor]
    lab = _inert_cleared(T, ps, tuple(target_labels), k)
    if all([x < 2 for x in lab]):
        # unlabelled signs only: a labelling with fewer signs states a different ABS reading, and one
        # with as many is the same compound up to symmetry, so the axiom's answer is the least
        v = _labelled_clone(base, anchor, ps, lab, [None, None, [], []])
        return anchor, ps, lab, v, mol_identity_bytes(v)
    target = _readings(T, ps, lab, k, unit_comp)
    clusters = _clusters(T, k)
    stated = _stated_clusters(T, ps, lab, k, clusters)
    central = _central_units(T, k, _unit_distances(view, anchor))
    for lab in _label_shapes(k):
        levels.setdefault(_labelling_measure_codes(lab, central), []).append(lab)

    for measure in sorted(levels):
        for lab in levels[measure]:
            # an AND group's two phases state one thing: its first member is held at parity 1
            fixed = 0
            seen = set()
            for j in range(k):
                if lab[j] >= 3 and lab[j] & 1 and lab[j] not in seen:
                    seen.add(lab[j])
                    fixed |= 1u << j
            fixed_free = [j for j in range(k) if lab[j] and not (fixed >> j) & 1]
            for bits in range(1u << len(fixed_free)):
                tried += 1
                if tried > LABELLING_CANDIDATE_LIMIT:
                    return False
                q = 0
                for j in range(len(fixed_free)):
                    if (bits >> j) & 1:
                        q |= 1u << <uint32_t> fixed_free[j]
                ref = _reading(T, q, lab, k)
                if ref not in target:
                    continue
                if any([(lab[j] == 1 or lab[j] == 2) and _reading(T, q ^ (1u << j), lab, k) == ref
                        for j in range(k)]):
                    continue
                if _readings(T, q, lab, k, unit_comp) != target:
                    continue
                if _stated_clusters(T, q, lab, k, clusters) != stated:
                    continue
                v = _labelled_clone(base, anchor, q, lab, [None, None, list(free_ids[2]), list(free_ids[3])])
                identity = mol_identity_bytes(v)
                if best is None or identity < best_id:
                    best = v
                    best_id = identity
                    best_q = q
                    best_lab = lab
        if best is not None:
            return anchor, best_q, best_lab, best, best_id
    return False


cdef list _free_group_ids(Structure s):
    """Per kind, the OR and AND group ids `s` does not hold, ascending."""
    cdef uint32_t i, n = s.header.atom_count
    cdef uint8_t *sg = structure_stereo_groups(s)
    cdef list out = [None, None, [], []]
    cdef set seen
    cdef object kind
    for kind in (SG_KIND_ABS + 1, SG_KIND_AND):
        seen = {sg_group(sg[i]) for i in range(n) if sg[i] and sg_kind(sg[i]) == kind}
        out[kind] = [g for g in range(1, STEREO_GROUP_MAX + 1) if g not in seen]
    return out


cdef int _take_components(Structure out, Structure src, set comps, list cof) except -1:
    """Copy the parities and group bytes of the atoms of `comps` from `src` into `out`, each OR and
    AND group of `src` under an id `out` does not hold.  The caller refreshes the stereo units."""
    cdef uint32_t i, n = out.header.atom_count
    cdef uint8_t *opar = structure_parities(out)
    cdef uint8_t *osg = structure_stereo_groups(out)
    cdef uint8_t *spar = structure_parities(src)
    cdef uint8_t *ssg = structure_stereo_groups(src) if structure_has(src, SEG_STEREO_GROUPS) else NULL
    cdef list free_ids = _free_group_ids(out)
    cdef dict number = {}
    cdef uint8_t b
    for i in range(n):
        if cof[i] not in comps:
            continue
        opar[i] = spar[i]
        b = ssg[i] if ssg is not NULL else 0
        if b and sg_kind(b) != SG_KIND_ABS:
            if b not in number:
                number[b] = sg_pack(sg_kind(b), <uint8_t> free_ids[sg_kind(b)].pop(0))
            b = number[b]
        osg[i] = b
    return 0


cdef object _meaning_canonical(Structure structure, Structure view, dict cache, list approximate=None,
                               list inert=None):
    """The least labelling stating what `structure` states: a Structure, `structure` itself when it is
    already least; None when no component needs the search, so `_relabelled` answers; False when no
    searched block is answered and one is over a limit, so `_relabelled` answers too.

    `view` is `structure_stereogenic_view(structure)`, which the components not searched are read
    from.  The searched ones are read from `structure`: a sign the view clears as not stereogenic in
    the drawing can be stereogenic in a reading of an unlabelled sign.  A component is searched when
    it holds a configured unit and two of its labellings may state one thing: a unit some parities
    make stereogenic, two anchors in one constitutional orbit, or a symmetry fixing the anchors that
    inverts units together.

    BLOCKS ARE SEARCHED APART.  A block is a searched component with every component a chain of
    AND or OR groups joins it to (`_coupled_blocks`).  No choice a meaning ranges over reaches across
    two blocks: an OR choice and an AND phase are per group, a reading of unlabelled signs is per
    component, a cleared unit is one unit, and ABS states each member alone.  A block's pure table,
    taken with every other block cleared, partitions its combinations as the whole table would,
    since a mixture's identity is the multiset of its components' and the rest is held fixed.  So the
    meaning of the molecule is the tuple of its blocks' meanings, the measure adds over blocks, and
    the least labelling is each block's least, ties broken by identity against the same cleared rest.
    The limits apply per block.  One over a limit takes `_relabelled`'s labels, which can spell one
    meaning two ways, and its anchors are appended to `approximate`.

    Cost per block: one identity per pure combination of its units, then the candidates of each
    measure in ascending order until one states the meaning, an identity per candidate that does.
    """
    cdef uint32_t n = view.header.atom_count, i, flags = 0
    cdef bytes codes, out_id
    cdef const unsigned char *code
    cdef uint8_t *par
    cdef uint8_t *sg = NULL
    cdef uint8_t *bpar
    cdef uint8_t *bsg
    cdef uint32_t *comp
    cdef uint32_t *orbits = NULL
    cdef Structure v, base, out, fallback
    cdef dict anchors = {}, pcache = {}
    cdef set configured = set(), weak = set(), open_ = set(), scope, seen, searched = set()
    cdef set kept = set(), over = set()
    cdef list cof, rest, blocks, block_anchors = [], results = [], free_ids
    cdef object c, cl, r, anchor
    if n == 0 or not structure_has(structure, SEG_PARITY):
        return None
    if not _may_differ(structure, inert):
        return None
    codes = _unit_codes(view, cache)
    if codes is None:
        return None
    ensure_stereo_units(view)
    ensure_component_labels(view)
    code = <const unsigned char *> codes
    comp = structure_component_labels(view)
    cof = [comp[i] for i in range(n)]
    par = structure_parities(structure)
    if structure_has(structure, SEG_STEREO_GROUPS):
        sg = structure_stereo_groups(structure)
    for i in range(n):
        if not code[i] or inert and i in inert:
            continue
        anchors.setdefault(cof[i], []).append(i)
        if par[i]:
            configured.add(cof[i])
            if sg is not NULL and sg[i]:
                open_.add(cof[i])
        else:
            open_.add(cof[i])
        if code[i] == 2:
            weak.add(cof[i])
    # unlabelled signs and no unit open: the axiom is single-flip stereogenicity, which the view reads
    configured &= open_
    scope = configured & weak
    rest = [c for c in configured if c not in weak and len(anchors[c]) > 1]
    if rest:
        v = structure_clone(view)
        memset(structure_parities(v), 0, <size_t> n)
        structure_invalidate_stereo_units(v)
        refresh_parity_features(v)
        orbits = <uint32_t *> PyMem_Malloc(<size_t> n * sizeof(uint32_t))
        if orbits is NULL:
            raise MemoryError()
        try:
            mol_automorphisms(v, NULL, orbits, &flags)
            if flags & CANON_BUDGET_EXCEEDED:
                return None
            for c in rest:
                seen = set()
                for i in anchors[c]:
                    if orbits[i] in seen:
                        scope.add(c)
                        break
                    seen.add(orbits[i])
            rest = [c for c in rest if c not in scope]
            if rest and _anchored_swap(v, view, orbits, anchors, rest, scope):
                return None
        finally:
            PyMem_Free(orbits)
    if not scope:
        return None
    blocks = _coupled_blocks(scope, sg, par, cof, n)
    for cl in blocks:
        searched |= cl
        block_anchors.append(sorted([i for c in cl for i in anchors.get(c, ()) if code[i] == 1 or par[i]]))

    # every block unlabelled and its units cleared
    if sg is NULL:
        base = structure_with_stereo_groups(view)
        rebuild_derived(base)
    else:
        base = structure_clone(view)
    bsg = structure_stereo_groups(base)
    bpar = structure_parities(base)
    for i in range(n):
        if cof[i] in searched:
            bsg[i] = 0
    for anchor in block_anchors:
        for i in anchor:
            bpar[i] = 0
    for c in _coupled_sets(_coupled_bridgeheads(view, cache)):
        if cof[<uint32_t> c] in searched:
            bpar[<uint32_t> c] = 0      # an implied member, which `_unit_codes` names no choice on
    structure_invalidate_stereo_units(base)
    refresh_parity_features(base)
    # group ids free per kind: those the searched components held, and every unused one
    free_ids = _free_group_ids(base)

    try:
        for cl, anchor in zip(blocks, block_anchors):
            r = _block_least(structure, view, base, cl, anchor, cof, free_ids, pcache)
            results.append(r)
            if r is None:
                kept |= cl
            elif r is False:
                over |= cl
                if approximate is not None:
                    approximate.extend(anchor)
    except AutomorphismBudgetExceeded:
        return None
    if not [r for r in results if r]:
        return False if over else None
    if len(blocks) == 1:
        out = results[0][3]
        out_id = results[0][4]
    else:
        out = structure_clone(base)
        if kept:
            _take_components(out, view, kept, cof)
        if over:
            fallback = _without_dead_groups(structure, view, None, cache)
            fallback = _relabelled(fallback, cache, True)
            fallback = _without_dead_groups(fallback, structure_stereogenic_view(fallback, cache), None, cache)
            _take_components(out, fallback, over, cof)
        free_ids = _free_group_ids(out)
        for r in results:
            if r:
                _label_into(out, r[0], r[1], r[2], free_ids)
        structure_invalidate_stereo_units(out)
        refresh_parity_features(out)
        try:
            out_id = mol_identity_bytes(out)
        except AutomorphismBudgetExceeded:
            return None
    # a symmetric molecule states one labelling in several atom placements: the stored one is kept when
    # it is one of them, so a string read back is already least in its written order
    try:
        if mol_identity_bytes(structure) == out_id:
            return structure
    except AutomorphismBudgetExceeded:
        pass
    return out
