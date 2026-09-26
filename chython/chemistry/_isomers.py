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
"""Canonical placement of a mobile hydrogen or charge over a conjugated nitrogen system.

The shift `standardize()`'s local rules cannot see -- the two N-H forms of 4-methylimidazole are one
compound, and both kekulise, so there is nothing to repair and only a choice to make the same way
whichever form arrived.  Placements are ranked in a placement-stripped frame (`_ranks`), proved by the
kekuliser rather than scored (`_admissible`), and decided per system (`_choose`).

A charge is mobile on all three shapes and in both signs: which nitrogen of an imidazolium holds the `+`,
or which one of a guanidide holds the `-`, is a fact about the spelling.  Hydrogen and charge are dealt
SEPARATELY, since the two do not always travel together -- 1-methylimidazolium is `Cn1cc[nH+]c1` and
`C[n+]1cc[nH]c1` for the one ion, its substituted nitrogen taking the charge and never the hydrogen.  Net
charge is `neutralize()`'s and no placement here creates or destroys one.

Three shapes carry a mobile hydrogen, and only the first is aromatic:

* an **aromatic** ring system -- pyrazole, imidazole, the purines.  `_sites`.
* a **kekule** ring system that is conjugated but stored with alternating orders, which is how a
  lactam is stored: `thiele()` refuses to aromatise a pyridone, deliberately, so `4-methylpyrimidin-2-one`
  arrives as `CC1=CC=NC(=O)N1` or `CC1=NC(=O)NC=C1` and no aromatic bond is anywhere in it.  Rather
  than a second placement algorithm these systems are *spelled* aromatic on a working copy
  (`_aromatized`), decided by the one above, and kekulised back -- so the same canonical order answers
  both, which two algorithms could not promise.
* an **amidine** or **guanidine**, where the hydrogen moves between the nitrogens of one sp2 carbon
  and the C=N moves with it.  No ring, so no kekuliser: which nitrogen may take the double bond is a
  closed-form question about its valence (`_amidines`).  Units sharing a nitrogen -- a biguanide's
  bridge -- are one placement problem, since that nitrogen may accept one C=N and not two
  (`_amidine_systems`).

What is deliberately *not* here: moving a hydrogen from oxygen or sulfur onto nitrogen.  An enol and
its ketone are not two valid spellings of one drawing, they are a drawing to repair, and the repair is
a `tautomer` row in `standardize_groups.tsv` with a `why` -- which is why this module logs `INFO` and
never `REPAIRED`.
"""
from collections.abc import Iterable, Sequence
from itertools import combinations, product
from math import comb
from ..core import INFO, LOST, LogRecord, MoleculeContainer, recording


__all__ = ['standardize_isomers']


#: Table-qualified, as every rule id must be -- never a bare index.  Three and not one: a consumer
#: filtering the log can tell a hydrogen that moved from a ring whose bond orders were rewritten with
#: it, which is the more invasive edit even though the compound is the same either way.
_RULE = 'isomers:placement'
_RULE_KEKULE = 'isomers:kekule-placement'
_RULE_AMIDINE = 'isomers:amidine'
_RULE_BUDGET = 'isomers:placement-budget'
_RULE_REFUSED = 'isomers:placement-refused'

#: The elements whose ring hydrogen is mobile.  O and S are absent deliberately: the classification
#: table answers must-not for them at every count, so furan's and thiophene's hydrogens never moved.
_MOBILE = frozenset((7, 15, 33))

#: Admitted into a ring system that is about to be spelled aromatic, as the lone-pair donor a furan
#: oxygen is.  Never a site: see the module docstring on enols.
_DONOR = frozenset((8, 16))

#: Per group, not per molecule: that is what makes the answer independent of the other groups.  A
#: five-site group with two protons and one charge is 30 trials, and reaching the cap is a `LOST`
#: record rather than a silent truncation.
_TRIALS_MAX = 512


def _sites(molecule: MoleculeContainer) -> list[int]:
    """Every atom whose hydrogen or charge the ring decides.

    A non-radical N, P or As carrying at most one unit of charge, with exactly two aromatic bonds and at
    most one hydrogen -- nothing in the graph says whether it donates its lone pair.  A non-aromatic bond
    has already decided the question, and an atom holding its hydrogen explicitly is not a site either,
    which `implicify_hydrogens()` converts to this spelling first.

    Cationic as well as anionic: which nitrogen of an imidazolium carries the `+` is a fact about the
    spelling and not about the ion, so it is as mobile as the hydrogen beside it.  A SUBSTITUTED nitrogen
    is a site on the charge alone: 1-methylimidazolium is drawn `Cn1cc[nH+]c1` and `C[n+]1cc[nH]c1` for
    the one ion, its `+` mobile where its hydrogen is not, which is why the two are dealt separately in
    `_choose` and why a third bond does not disqualify.

    Pure reads, and it must stay that way: the container refuses a read once an edit session is open.
    """
    out: list[int] = []
    for n in molecule.atom_numbers:
        if molecule.element_of(n) not in _MOBILE or molecule.radical_of(n):
            continue
        if molecule.charge_of(n) not in (-1, 0, 1):
            continue
        h = molecule.implicit_h_of(n)
        if h is None or h > 1:
            continue                      # an unknown count is not a placement anybody may choose
        orders = [molecule.order_of(n, m) for m in molecule.neighbors_of(n)]
        if orders.count(4) != 2 or any(o not in (1, 4) for o in orders):
            continue
        out.append(n)
    return out


def _groups(molecule: MoleculeContainer, sites: list[int]) -> list[list[int]]:
    """The sites split by aromatic ring SYSTEM, so each system is decided on its own.

    Connected components over aromatic bonds, walked from the sites and through any atom, because two
    nitrogens of one fused system are one problem however many carbons lie between them.
    """
    marked = set(sites)
    seen = set()
    out: list[list[int]] = []
    for start in sites:
        if start in seen:
            continue
        stack = [start]
        seen.add(start)
        members = []
        while stack:
            v = stack.pop()
            if v in marked:
                members.append(v)
            for w in molecule.neighbors_of(v):
                if w not in seen and molecule.order_of(v, w) == 4:
                    seen.add(w)
                    stack.append(w)
        members.sort()
        out.append(members)
    return out


def _ranks(molecule: MoleculeContainer, sites: Iterable[int],
           bonds: Sequence[tuple[int, int]] = ()) -> dict[int, int]:
    """`atoms_order` of a copy with every site stripped -- the spelling-independent frame.

    Ranking the molecule as written cannot decide where a hydrogen goes, since the ranks depend on where
    it already is; stripped, two tautomers of one compound are the same graph and rank identically.  The
    stripped copy need not be kekulizable, nothing being asked of it but its canonical order.

    Per connected component, and that is load-bearing rather than an optimisation: `atoms_order` is a
    total order over the whole molecule, so it must separate two automorphic components, and eight
    identical pyrazoles ranked in one frame got eight different placements.  `substructure` preserves
    stable ids, so the dicts merge; ranks collide across components, harmlessly, since two groups are
    never compared with each other.

    `bonds` is stripped too, to order one.  A ring hydrogen is the whole difference between two annular
    spellings, but an amidine's is not: its C=N moves with it, so the frame must forget the double bond
    as well or the two drawings rank differently and each keeps its own.
    """
    work = molecule.copy()
    with work.edit():
        for u, v in bonds:
            work.set_order(u, v, 1)
        for n in sites:
            work.set_charge(n, 0)
            work.set_hydrogens(n, 0)
    components = work.connected_components
    if len(components) == 1:
        return work.atoms_order
    out: dict[int, int] = {}
    for component in components:
        out.update(work.substructure(component).atoms_order)
    return out


def _budget(molecule: MoleculeContainer, sites: Iterable[int]) -> tuple[int, dict[int, int]]:
    """What the sites arrived holding: how many hydrogens in total, and the multiset of their charges.

    A placement deals both back out over the same sites, which is what conserves the formula and every
    component's charge: a hydrogen and a `+` may move from one site to another and neither may be created.
    Net charge is `neutralize()`'s question and none of this pass's.

    Two budgets and not one multiset of `(hydrogens, charge)` pairs, because the two do not travel
    together: a substituted nitrogen takes a `+` and no hydrogen, so 1-methylimidazolium holds the pairs
    `(0, +1), (1, 0)` drawn one way and `(0, 0), (1, +1)` the other.  Dealing pairs could never bring
    those two drawings together; dealing one hydrogen and one `+` over the same two sites does.  The
    charge multiset rather than the net: a `+` and a `-` on one system is a charge separation, which is
    `fix_resonance()`'s to remove and never this pass's to invent.
    """
    charges: dict[int, int] = {}
    hydrogens = 0
    for n in sites:
        hydrogens += molecule.implicit_h_of(n) or 0
        q = molecule.charge_of(n)
        charges[q] = charges.get(q, 0) + 1
    return hydrogens, charges


def _deal(sites: tuple[int, ...], states: tuple[tuple, ...]):
    """Every distinct way to deal the states over the sites, as `{site: state}`.

    One state per site and every state dealt, so an assignment is a permutation of the multiset.  Both
    paths deal their charges with this, zeros included so that every site is covered; hydrogens are a
    count rather than a multiset and are chosen with `combinations` beside it on the ring path, and follow
    from the double bond on the amidine one.
    """
    if not states:
        yield {}
        return
    (state, count), rest = states[0], states[1:]
    for chosen in combinations(sites, count):
        remaining = tuple(n for n in sites if n not in chosen)
        for tail in _deal(remaining, rest):
            yield {**{n: state for n in chosen}, **tail}


def _admissible(molecule: MoleculeContainer, group: list[int],
                placement: dict[int, tuple[int, int]]) -> int | None:
    """How many bonds `thiele()` writes aromatic under this placement, or `None` when it gives no Kekule
    form for this group.  The kekuliser is the oracle.

    This group's systems and not the whole molecule's: one ring nobody can kekulise (`c1cccc1`) would
    otherwise make every placement of every other ring inadmissible.

    Asking the oracle is not enough on its own, because `kekule()` repairs by design: handed a ring that
    cannot carry the hydrogens it was given, it drops one and logs it, and handed a neutral nitrogen that
    has to be a cation it charges it.  A placement that kekulises only because the kekuliser rewrote it
    is a placement of a different molecule, so the counts are read back and compared with what was asked
    for -- which is also what tells a `+` this pass placed from one the kekuliser added.  1,2,4-triazol-3-
    one taught this: protonating both nitrogens flanking its lone ring carbon leaves that carbon no
    partner for a double bond, and the relaxation hid it by taking the hydrogens away.
    """
    work = molecule.copy()
    with work.edit():
        for n, (hydrogens, charge) in placement.items():
            work.set_charge(n, charge)
            work.set_hydrogens(n, hydrogens)
    members = frozenset(group)
    if any(members.intersection(system) for system in work.kekule().unresolved):
        return None
    if any((work.implicit_h_of(n), work.charge_of(n)) != state for n, state in placement.items()):
        return None
    work.thiele()
    return work.aromatic_bond_count


def _choose(molecule: MoleculeContainer, group: list[int], ranks: dict[int, int]):
    """The canonical placement for one group, `None` when nothing must move, `'budget'` when too big.

    The hydrogens and the charges are read off the group rather than assumed, and dealt back over it.  The
    key is the aromatic bond count `thiele()` gives the placement, most first, then the hydrogens beside a
    C=O or C=S carbon, then the hydrogens on five-membered rings, each most first, then the sorted ranks
    of the sites holding the hydrogens, then of those holding each charge in turn -- a strict total order,
    because `atoms_order` is a permutation and those two sets fix the placement, so no two placements
    share a key and no tie is left for an arbitrary rule to break.

    A lactam N-H sits next to its carbonyl: 6-methylpyrimidin-4(3H)-one, not its 1H form.  An azole N-H
    outranks an azine N-H whenever both forms are equally aromatic, so every 7-azaindole reads 1H
    (`c1cnc2[nH]ccc2c1`), its carboxylic acid too.  Ring size and bond order are graph facts, so both
    terms are as spelling-independent as the ranks after them.

    Aromaticity leads because a fused system spelled aromatic joins rings a hydrogen cannot cross for
    free: pyrido[4,3-d]pyrimidine-2,4-dione `O=C1NC(=O)c2cnccc2N1` with its N1-H moved onto the pyridine
    N kekulises, as a quinoid imine with no pyridine sextet.  Only nitrogen sites deal here and `thiele()`
    leaves a lactam non-aromatic, so the term never weighs an enol against its ketone.

    A hydrogen deals over the sites with two heavy neighbours and a charge over all of them, which is the
    whole difference a substituted nitrogen makes.  That much is not left to the oracle: `kekule()` decides
    bond orders and not valences, and it will happily find a Kekule form for a 1-methylpyrazolium whose
    hydrogen was put on its methylated nitrogen -- four bonds on a neutral nitrogen, kekulised and
    over-valent.
    """
    hydrogens, charges = _budget(molecule, group)
    carriers = [n for n in group if len(tuple(molecule.neighbors_of(n))) == 2]
    order = tuple((q, charges[q]) for q in sorted(charges))
    trials, rest = comb(len(carriers), hydrogens), len(group)
    for _, count in order:
        trials *= comb(rest, count)
        rest -= count
        if trials > _TRIALS_MAX:
            return 'budget'
    if trials == 1:
        return None                       # one way to deal them: there is no distribution to choose
    signs = [q for q in sorted(charges) if q]
    azole = {n for n in carriers if 5 in molecule.ring_sizes_of(n)}
    lactam = {n for n in carriers
              if any(any(molecule.order_of(c, x) == 2 and molecule.element_of(x) in _DONOR
                         for x in molecule.neighbors_of(c)) for c in molecule.neighbors_of(n))}
    current = {n: (molecule.implicit_h_of(n) or 0, molecule.charge_of(n)) for n in group}
    best = None
    for protonated in combinations(carriers, hydrogens):
        for dealt in _deal(tuple(group), order):
            placement = {n: (1 if n in protonated else 0, dealt[n]) for n in group}
            aromatic = _admissible(molecule, group, placement)
            if aromatic is None:
                continue
            key = [[-aromatic, -sum(n in lactam for n in protonated), -sum(n in azole for n in protonated)],
                   sorted(ranks[n] for n in protonated)]
            key.extend(sorted(ranks[n] for n in group if dealt[n] == q) for q in signs)
            if best is None or key < best[0]:
                best = (key, placement)
    if best is None or best[1] == current:
        return None
    return best[1]


# ---------------------------------------------------------------------------------------------------
# The kekule ring systems.  A lactam is stored non-aromatic on purpose, so the sites above cannot see
# it; rather than a second placement algorithm the system is spelled aromatic on a working copy and
# handed to the one above.
# ---------------------------------------------------------------------------------------------------

def _ring_systems(molecule: MoleculeContainer) -> list[tuple[frozenset[int], tuple[tuple[int, int], ...]]]:
    """Fused ring systems as `(atoms, ring bonds)`, one per component of the ring-bond graph.

    Built from `sssr` rather than from `bond_in_ring`, because a system is exactly what shares a ring
    bond: two rings joined at one spiro atom are two systems and must be decided as two.
    """
    adjacency: dict[int, set] = {}
    for ring in molecule.sssr:
        for i, a in enumerate(ring):
            b = ring[i - 1]
            adjacency.setdefault(a, set()).add(b)
            adjacency.setdefault(b, set()).add(a)
    seen = set()
    out = []
    for start in adjacency:
        if start in seen:
            continue
        stack, atoms = [start], set()
        seen.add(start)
        while stack:
            v = stack.pop()
            atoms.add(v)
            for w in adjacency[v]:
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
        out.append((frozenset(atoms),
                    tuple(sorted((a, b) for a in atoms for b in adjacency[a] if a < b))))
    return out


def _conjugated(molecule: MoleculeContainer, atoms: frozenset[int]) -> bool:
    """Is every atom of this system sp2, holding either one double bond or a lone pair?

    The gate that keeps this module out of tautomerism it has no business doing.  One sp3 carbon and the
    whole system is refused, so cyclohexa-2,4-dien-1-one never becomes phenol and barbituric acid never
    loses the hydrogens on its CH2 -- both of which an aromatic spelling would quietly do.  Read on the
    kekule form, where a double bond still means what it says.

    An atom already carrying an aromatic bond is passed over: it belongs to an aromatic system, whose
    own gates are `_sites`, and its valence bookkeeping is the one `kekule()` owns.  That is what lets a
    part-aromatic fused system -- guanine's imidazole beside its pyrimidinone -- be decided as one.
    """
    for n in atoms:
        if molecule.charge_of(n) or molecule.radical_of(n):
            return False
        h = molecule.implicit_h_of(n)
        if h is None:
            return False
        orders = [molecule.order_of(n, m) for m in molecule.neighbors_of(n)]
        if any(o not in (1, 2, 4) for o in orders):
            return False                  # a triple or a dative bond is not this shape
        if 4 in orders:
            continue
        element = molecule.element_of(n)
        doubles = orders.count(2)
        if element == 6:
            if doubles != 1 or len(orders) + h != 3:
                return False              # exactly one double bond, and nothing sp3
        elif element in _MOBILE:
            if doubles > 1 or len(orders) + h != (2 if doubles else 3):
                return False              # `=N-` accepts a hydrogen, `-NH-` and `-N(C)-` donate one
        elif element in _DONOR:
            if doubles or len(orders) != 2 or h:
                return False              # a ring oxygen or sulfur lends its lone pair, as furan does
        else:
            return False
    return True


def _aromatized(molecule: MoleculeContainer):
    """A working copy with every eligible kekule ring system spelled aromatic.

    Returns `(copy, systems)`, `systems` being `(atoms, the bonds rewritten)` per system taken, or
    `(None, ())` when none qualifies.  Systems are taken one at a time and each is verified before it is
    kept: the aromatic spelling has to kekulise back and has to leave every hydrogen count where it was.
    A system the core cannot carry that way is skipped, never forced -- and skipping one does not cost
    the others.
    """
    candidates = []
    for atoms, bonds in _ring_systems(molecule):
        rewrite = tuple(b for b in bonds if molecule.order_of(*b) != 4)
        if not rewrite or not _conjugated(molecule, atoms):
            continue
        if sum(molecule.element_of(n) in _MOBILE and len(tuple(molecule.neighbors_of(n))) == 2
               for n in atoms) < 2:
            continue                       # one candidate nitrogen is one placement; there is no choice
        candidates.append((atoms, rewrite))
    if not candidates:
        return None, ()

    hydrogens = {n: molecule.implicit_h_of(n) for n in molecule.atom_numbers}
    work = molecule.copy()
    taken = []
    for atoms, rewrite in candidates:
        probe = work.copy()
        with probe.edit():
            for u, v in rewrite:
                probe.set_order(u, v, 4)
            for n in atoms:
                probe.set_hydrogens(n, hydrogens[n])
        # the count is restored explicitly above and re-read here: an aromatic pnictogen whose class the
        # ring decides answers `unknown`, and an unknown is not a placement anybody may choose.
        check = probe.copy()
        if check.kekule().unresolved:
            continue
        if any(check.implicit_h_of(n) != hydrogens[n] for n in atoms):
            continue
        work = probe
        taken.append((atoms, rewrite))
    return (work, tuple(taken)) if taken else (None, ())


# ---------------------------------------------------------------------------------------------------
# The amidines and guanidines.  No ring, so no kekuliser: which nitrogen may hold the double bond is a
# closed-form question about its valence.  Units sharing a nitrogen are decided together.
# ---------------------------------------------------------------------------------------------------

def _amidine_carbon(molecule: MoleculeContainer, n: int) -> bool:
    """An acyclic, neutral, non-radical sp2 carbon holding exactly one double bond.

    Carbon-local, so it can be answered before any nitrogen is classified: `_amidine_site` needs the
    carbon set and `_amidines` closes the loop between the two.
    """
    if molecule.element_of(n) != 6 or molecule.in_ring_of(n):
        return False
    if molecule.charge_of(n) or molecule.radical_of(n):
        return False
    h = molecule.implicit_h_of(n)
    if h is None:
        return False
    neighbors = tuple(molecule.neighbors_of(n))
    orders = [molecule.order_of(n, m) for m in neighbors]
    if any(o not in (1, 2) for o in orders) or orders.count(2) != 1:
        return False
    return len(neighbors) + h == 3


def _amidine_site(molecule: MoleculeContainer, n: int, carbon: int, carbons: set[int]) -> bool:
    """A nitrogen whose hydrogen count is decided by whether it takes this carbon's double bond.

    Non-radical, carrying at most one unit of charge, at most three heavy neighbours, and no multiple bond
    of its own anywhere else, which is what keeps a nitro or an azo group out.  A double bond to another
    amidine carbon is the one exception: that bond is mobile too, so the site set of a biguanide's bridge is
    the same whichever of the two carbons the drawing happened to put it on, which is what makes both
    drawings rank in one frame.

    A charged nitrogen is a site like any other, and its charge is placed with the double bond: which
    nitrogen of a guanidinium holds the `+` follows from which one holds the hydrogen.  Room is a question
    for the valence identity in `_choose_amidines` and not for a neighbour count, and it has to be: a
    tertiary nitrogen has no room for a hydrogen at all neutral, one as a cation, so gating on the count
    would make metformin's site set depend on where its `+` was drawn.  Four heavy neighbours is the count
    that decides, holding neither hydrogen nor double bond under any charge in range.
    """
    if molecule.element_of(n) != 7 or molecule.in_ring_of(n):
        return False
    if molecule.charge_of(n) not in (-1, 0, 1) or molecule.radical_of(n):
        return False
    if molecule.implicit_h_of(n) is None:
        return False
    neighbors = tuple(molecule.neighbors_of(n))
    if len(neighbors) > 3 or molecule.order_of(n, carbon) not in (1, 2):
        return False
    return all(molecule.order_of(n, m) == 1 or (molecule.order_of(n, m) == 2 and m in carbons)
               for m in neighbors if m != carbon)


def _amidines(molecule: MoleculeContainer) -> list[tuple[int, list[int]]]:
    """`(carbon, its mobile nitrogens)` per acyclic amidine or guanidine unit.

    Acyclic on both counts, the carbon and every nitrogen, so an amidine unit can never overlap a ring
    system `_aromatized` handed upward and the two answers cannot contradict each other.  A cyclic
    amidine is left for the ring path to reach through its ring.

    A unit is kept only when its one double bond lands on a site -- which is what an amide's C=O fails --
    and dropping a carbon can take a site away from its neighbour, so the two halves are iterated to a
    fixed point rather than computed once.  One site is a unit: the acceptor is forced, but a forced
    acceptor still denies a shared nitrogen to the unit next to it.
    """
    carbons = {n for n in molecule.atom_numbers if _amidine_carbon(molecule, n)}
    while True:
        sites = {c: sorted(m for m in molecule.neighbors_of(c) if _amidine_site(molecule, m, c, carbons))
                 for c in carbons}
        dropped = {c for c in carbons if not any(molecule.order_of(c, m) == 2 for m in sites[c])}
        if not dropped:
            return [(c, sites[c]) for c in sorted(carbons)]
        carbons -= dropped


def _amidine_systems(groups: list[tuple[int, list[int]]]) -> list[list[tuple[int, list[int]]]]:
    """The units merged where they share a nitrogen -- one system per connected component.

    A bridging nitrogen accepts one C=N and not two, so two units sharing one are not independent and a
    biguanide decided as two amidines puts a double bond on its bridge twice.
    """
    owners: dict[int, list[int]] = {}
    for carbon, sites in groups:
        for n in sites:
            owners.setdefault(n, []).append(carbon)
    lookup = dict(groups)
    seen: set[int] = set()
    out = []
    for start, _ in groups:
        if start in seen:
            continue
        stack, members = [start], []
        seen.add(start)
        while stack:
            c = stack.pop()
            members.append(c)
            for n in lookup[c]:
                for other in owners[n]:
                    if other not in seen:
                        seen.add(other)
                        stack.append(other)
        members.sort()
        out.append([(c, lookup[c]) for c in members])
    return out


def _choose_amidines(molecule: MoleculeContainer, unit: list[tuple[int, list[int]]],
                     ranks: dict[int, int]):
    """`({carbon: acceptor}, {site: charge})` for one system, `None` when it already has it, `'budget'`
    when too big.

    A nitrogen with charge `q` and `d` heavy neighbours carries `3 + q - d` hydrogens single-bonded and
    `2 + q - d` double-bonded, so with the charges dealt the total over the system is fixed whichever
    nitrogens accept and there is nothing to prove admissible -- only a count no nitrogen can carry to
    reject.  What has to be enforced is that each nitrogen accepts at most once, which makes an assignment
    a perfect matching of the carbons onto their sites.

    The candidate graph is a forest -- a cycle through it is a ring, and every atom here is acyclic -- so
    no two assignments share an acceptor set.  That is what leaves the key a strict total order: the
    sorted ranks of the nitrogens that KEEP their hydrogen, then the sorted ranks holding each charge,
    the same shape of key the ring path uses.
    """
    carbons = [c for c, _ in unit]
    all_sites = sorted({n for _, sites in unit for n in sites})
    degrees = {n: len(tuple(molecule.neighbors_of(n))) for n in all_sites}
    current = {c: next((n for n in sites if molecule.order_of(c, n) == 2), None) for c, sites in unit}
    # a nitrogen already holding two C=N is over-valent, and repairing one is `standardize()`'s job
    if len(set(current.values())) != len(current):
        return None
    charged = {n: molecule.charge_of(n) for n in all_sites}
    counts: dict[int, int] = {}
    for n in all_sites:
        counts[charged[n]] = counts.get(charged[n], 0) + 1
    order = tuple((charge, counts[charge]) for charge in sorted(counts))

    trials, rest = 1, len(all_sites)
    for _, count in order:
        trials *= comb(rest, count)
        rest -= count
    for _, sites in unit:
        trials *= len(sites)
    if trials > _TRIALS_MAX:
        return 'budget'

    best = None
    for charges in _deal(tuple(all_sites), order):
        free = {n: 3 + charges[n] - degrees[n] for n in all_sites}
        for choice in product(*(sites for _, sites in unit)):
            taken = frozenset(choice)
            if len(taken) != len(choice):
                continue                  # two carbons cannot share one acceptor
            if any(free[n] < (1 if n in taken else 0) for n in all_sites):
                continue                  # a nitrogen asked for a hydrogen it cannot carry
            key = [sorted(ranks[n] for n in all_sites if n not in taken)]
            key.extend(sorted(ranks[n] for n in all_sites if charges[n] == charge)
                       for charge, _ in order)
            if best is None or key < best[0]:
                best = (key, dict(zip(carbons, choice)), charges)
    if best is None or (best[1], best[2]) == (current, charged):
        return None
    return best[1], best[2]


# ---------------------------------------------------------------------------------------------------
# The pass.
# ---------------------------------------------------------------------------------------------------

def _place_rings(molecule: MoleculeContainer, work, systems, sites: list[int],
                 ranks: dict[int, int], lines: list[tuple[str, tuple[int, ...], str]]) -> bool:
    """Decide every ring system, aromatic as drawn or spelled aromatic for the purpose.  Did it move?"""
    target = work if work is not None else molecule
    if len(sites) < 2:
        return False                      # one site is one placement; zero is none

    # every read happens before the session opens: the container refuses a read while a journal is
    # pending, so the whole plan is computed first and the session is pure writes.
    rule = _RULE if work is None else _RULE_KEKULE
    plan: list[tuple[list[int], dict[int, tuple[int, int]]]] = []
    for group in _groups(target, sites):
        if len(group) < 2:
            continue
        answer = _choose(target, group, ranks)
        if answer is None:
            continue
        if answer == 'budget':
            lines.append((_RULE_BUDGET, tuple(group),
                          f'ring system {tuple(group)!r} has {len(group)} mobile sites, more '
                          f'placements than the {_TRIALS_MAX}-trial budget allows; it was left as '
                          f'drawn rather than half-searched'))
            continue
        plan.append((group, answer))
        lines.append((rule, tuple(group),
                      f'ring system {tuple(group)!r}: mobile hydrogen(s) placed on '
                      f'{tuple(n for n in group if answer[n][0])!r} and charge(s) on '
                      f'{tuple(n for n in group if answer[n][1])!r}, the canonical placement for this '
                      f'skeleton'))

    if not plan:
        return False

    if work is None:
        with molecule.edit():
            for _, placement in plan:
                for n, (hydrogens, charge) in placement.items():
                    molecule.set_charge(n, charge)
                    molecule.set_hydrogens(n, hydrogens)
                # no parity restore, deliberately: `set_hydrogens` and `set_charge` do not clear one,
                # only `delete_atom` does.  Pinned by
                # `test_a_stereocentre_is_not_touched_and_needs_no_parity_restore`.
        return True

    # The placement was decided on an aromatic spelling the caller never asked for, so the working copy
    # is kekulised and only the bonds this module itself spelled aromatic are read back.  A system whose
    # group did not move keeps the orders it was drawn with: the pass answers where a hydrogen goes, and
    # rewriting a Kekule form nobody asked about is not that answer.
    with work.edit():
        for _, placement in plan:
            for n, (hydrogens, charge) in placement.items():
                work.set_charge(n, charge)
                work.set_hydrogens(n, hydrogens)
    if work.kekule().unresolved:
        lines.append((_RULE_REFUSED, tuple(sorted(n for group, _ in plan for n in group)),
                      'the canonical placement has no Kekule form for the molecule as a whole, though '
                      'it had one for each system alone; nothing was written'))
        return False

    moved = frozenset(n for group, _ in plan for n in group)
    placed = {n: state for _, placement in plan for n, state in placement.items()}
    orders: list[tuple[int, int, int]] = []
    for atoms, rewrite in systems:
        if atoms.isdisjoint(moved):
            continue
        if any(work.implicit_h_of(n) != molecule.implicit_h_of(n)
               for n in atoms if n not in moved) \
                or any((work.implicit_h_of(n), work.charge_of(n)) != placed[n]
                       for n in atoms if n in moved):
            lines.append((_RULE_REFUSED, tuple(sorted(atoms)),
                          f'the aromatic round trip of ring system {tuple(sorted(atoms))!r} did not '
                          f'give back the hydrogen counts the placement asked for; nothing was written'))
            return False                  # all or nothing: half a plan is a corrupted molecule
        orders.extend((u, v, work.order_of(u, v)) for u, v in rewrite)

    # the plan's counts and not the working copy's: `kekule()` is allowed to repair, and the guard above
    # only proves it did not need to here.  Writing what was decided keeps the two readable side by side.
    with molecule.edit():
        for n in moved:
            hydrogens, charge = placed[n]
            molecule.set_charge(n, charge)
            molecule.set_hydrogens(n, hydrogens)
        for u, v, order in orders:
            molecule.set_order(u, v, order)
    return True


def _place_amidines(molecule: MoleculeContainer, groups: list[tuple[int, list[int]]],
                    ranks: dict[int, int],
                    lines: list[tuple[str, tuple[int, ...], str]]) -> bool:
    """Decide every acyclic amidine and guanidine.  Did anything move?"""
    if not groups:
        return False
    plan: list[tuple[list[tuple[int, list[int]]], dict[int, int], dict[int, int], dict[int, int]]] = []
    for unit in _amidine_systems(groups):
        answer = _choose_amidines(molecule, unit, ranks)
        if answer is None:
            continue
        sites = sorted({n for _, s in unit for n in s})
        if answer == 'budget':
            lines.append((_RULE_BUDGET, tuple(sites),
                          f'amidine system {tuple(sites)!r} has more placements than the '
                          f'{_TRIALS_MAX}-trial budget allows; it was left as drawn rather than '
                          f'half-searched'))
            continue
        accepted, charges = answer
        acceptors = frozenset(accepted.values())
        hydrogens = {n: (2 if n in acceptors else 3) + charges[n] - len(tuple(molecule.neighbors_of(n)))
                     for n in sites}
        if min(hydrogens.values()) < 0:
            continue                      # no nitrogen may be asked for a hydrogen it does not have
        plan.append((unit, accepted, charges, hydrogens))
        lines.append((_RULE_AMIDINE, tuple(sites),
                      f'amidine at atom(s) {tuple(c for c, _ in unit)!r}: the double bond(s) to '
                      f'{tuple(sites)!r} placed on {tuple(accepted[c] for c, _ in unit)!r} and '
                      f'charge(s) on {tuple(n for n in sites if charges[n])!r}, the canonical assignment '
                      f'for this skeleton, and the hydrogens follow it'))

    if not plan:
        return False
    with molecule.edit():
        for unit, accepted, charges, hydrogens in plan:
            for carbon, sites in unit:
                for n in sites:
                    molecule.set_order(carbon, n, 2 if n == accepted[carbon] else 1)
            for n, h in hydrogens.items():
                molecule.set_charge(n, charges[n])
                molecule.set_hydrogens(n, h)
    return True


def standardize_isomers(molecule: MoleculeContainer) -> bool:
    """Put every mobile hydrogen and charge where the canonical order says it goes.  Did it move?

    The stage that makes two tautomers of one compound store the same molecule.  Three shapes carry a
    mobile hydrogen and all three are decided here -- an aromatic ring system, a conjugated ring system
    stored in its Kekule form as a lactam is, and an acyclic amidine or guanidine.  The aromatic form is
    not a precondition: `CC1=CC=NC(=O)N1` and `CC1=NC(=O)NC=C1` are one compound and land on one molecule,
    so `canonicalize()` runs `thiele()` first for `thiele()`'s own sake and not for this pass.

    Moving a hydrogen from oxygen to nitrogen is *not* this pass: an enol and its ketone are a drawing to
    repair, and the repair is a `tautomer` row in `standardize_groups.tsv`.

    `molecule.log` gets one `INFO` record per group whose placement changed -- information and not a
    repair, since every form involved was a valid molecule -- and one `LOST` record per group left as it
    arrived, whether for exceeding the trial budget or because the placement could not be written whole.
    Never raises.
    """
    lines: list[tuple[str, tuple[int, ...], str]] = []
    groups = _amidines(molecule)
    work, systems = _aromatized(molecule)
    target = work if work is not None else molecule
    sites = _sites(target)

    # ONE frame for both paths, and that is what makes them independent -- their gates do not overlap,
    # an amidine's excluding every ring atom, but their answers would: a ring hydrogen and an amidine's
    # C=N are both mobile, so a frame holding either ranks the other's placements by where the drawing
    # happened to put it.  A guanidinyl pyrazole decided in two frames was not idempotent.  The double
    # bonds are stripped along with the hydrogens, being half of what two amidine spellings differ by,
    # and the ring systems are read on the aromatic copy where no Kekule order says anything either.
    ranks = _ranks(target, sites + [n for _, s in groups for n in s],
                   [(carbon, n) for carbon, s in groups for n in s])

    changed = _place_rings(molecule, work, systems, sites, ranks, lines)
    changed = _place_amidines(molecule, groups, ranks, lines) or changed

    with recording(molecule, stage='isomers') as log:
        for rule, atoms, message in lines:
            log.append(LogRecord(rule, atoms, message,
                                 LOST if rule in (_RULE_BUDGET, _RULE_REFUSED) else INFO))
    return changed
