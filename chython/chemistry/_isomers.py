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
  closed-form question about its valence (`_amidines`).

What is deliberately *not* here: moving a hydrogen from oxygen or sulfur onto nitrogen.  An enol and
its ketone are not two valid spellings of one drawing, they are a drawing to repair, and the repair is
a `tautomer` row in `standardize_groups.tsv` with a `why` -- which is why this module logs `INFO` and
never `REPAIRED`.
"""
from collections.abc import Iterable, Sequence
from itertools import combinations
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

    A neutral or anionic, non-radical N, P or As with exactly two heavy neighbours, both bonds
    aromatic, and at most one hydrogen -- nothing in the graph says whether it donates its lone pair.
    Three neighbours leaves no room for a double bond either way, so it is not a site; nor is an atom
    holding its hydrogen explicitly, which `implicify_hydrogens()` converts to this spelling first.

    Pure reads, and it must stay that way: the container refuses a read once an edit session is open.
    """
    out: list[int] = []
    for n in molecule.atom_numbers:
        if molecule.element_of(n) not in _MOBILE or molecule.radical_of(n):
            continue
        if molecule.charge_of(n) not in (0, -1):
            continue
        h = molecule.implicit_h_of(n)
        if h is None or h > 1:
            continue                      # an unknown count is not a placement anybody may choose
        neighbors = tuple(molecule.neighbors_of(n))
        if len(neighbors) != 2 or any(molecule.order_of(n, m) != 4 for m in neighbors):
            continue                      # three neighbours has no room either way; a non-aromatic
        out.append(n)                     # bond has already decided the question
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


def _admissible(molecule: MoleculeContainer, group: list[int],
                protons: frozenset[int], anions: frozenset[int]) -> bool:
    """Does this placement give a Kekule form for this group?  The kekuliser is the oracle.

    This group's systems and not the whole molecule's: one ring nobody can kekulise (`c1cccc1`) would
    otherwise make every placement of every other ring inadmissible.

    Asking the oracle is not enough on its own, because `kekule()` repairs by design: handed a ring that
    cannot carry the hydrogens it was given, it drops one and logs it, and handed a neutral nitrogen that
    has to be a cation it charges it.  A placement that kekulises only because the kekuliser rewrote it
    is a placement of a different molecule, so the counts are read back and compared with what was asked
    for.  1,2,4-triazol-3-one taught this: protonating both nitrogens flanking its lone ring carbon
    leaves that carbon no partner for a double bond, and the relaxation hid it by taking the hydrogens
    away.
    """
    work = molecule.copy()
    with work.edit():
        for n in group:
            work.set_charge(n, -1 if n in anions else 0)
            work.set_hydrogens(n, 1 if n in protons else 0)
    members = frozenset(group)
    if any(members.intersection(system) for system in work.kekule().unresolved):
        return False
    return all(work.implicit_h_of(n) == (1 if n in protons else 0)
               and work.charge_of(n) == (-1 if n in anions else 0) for n in group)


def _choose(molecule: MoleculeContainer, group: list[int], ranks: dict[int, int]):
    """The canonical placement for one group, `None` when nothing must move, `'budget'` when too big.

    Protons and charges are two candidate sets over the same sites, both counts read off the group
    rather than assumed.  The key is `(sorted proton ranks, sorted charge ranks)`, a strict total order
    because `atoms_order` is a permutation -- so no two placements share a key and no tie is left for an
    arbitrary rule to break.
    """
    held = frozenset(n for n in group if molecule.implicit_h_of(n))
    charged = frozenset(n for n in group if molecule.charge_of(n) == -1)
    k, q = len(held), len(charged)
    if k + q == 0 or k + q == len(group):
        return None                       # every site alike: there is no distribution to choose
    trials = 0
    for protons in combinations(group, k):
        trials += len(tuple(combinations([n for n in group if n not in protons], q)))
        if trials > _TRIALS_MAX:
            return 'budget'
    best = None
    for protons in combinations(group, k):
        rest = [n for n in group if n not in protons]
        for anions in combinations(rest, q):
            key = (sorted(ranks[n] for n in protons), sorted(ranks[n] for n in anions))
            if best is not None and key >= best[0]:
                continue                  # cheaper than the oracle, so it goes first
            if _admissible(molecule, group, frozenset(protons), frozenset(anions)):
                best = (key, frozenset(protons), frozenset(anions))
    if best is None or (best[1], best[2]) == (held, charged):
        return None
    return best[1], best[2]


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
# closed-form question about its valence.
# ---------------------------------------------------------------------------------------------------

def _amidine_site(molecule: MoleculeContainer, n: int, carbon: int) -> bool:
    """A nitrogen whose hydrogen count is decided by whether it takes this carbon's double bond.

    Neutral, non-radical, at most two heavy neighbours -- three leaves no room for the double bond -- and
    no multiple bond of its own anywhere else, which is what keeps a nitro or an azo group out.
    """
    if molecule.element_of(n) != 7 or molecule.in_ring_of(n):
        return False
    if molecule.charge_of(n) or molecule.radical_of(n) or molecule.implicit_h_of(n) is None:
        return False
    neighbors = tuple(molecule.neighbors_of(n))
    if len(neighbors) > 2 or molecule.order_of(n, carbon) not in (1, 2):
        return False
    return all(molecule.order_of(n, m) == 1 for m in neighbors if m != carbon)


def _amidines(molecule: MoleculeContainer) -> list[tuple[int, list[int]]]:
    """`(carbon, its mobile nitrogens)` per acyclic amidine or guanidine.

    Acyclic on both counts, the carbon and every nitrogen, so an amidine group can never overlap a ring
    system `_aromatized` handed upward and the two answers cannot contradict each other.  A cyclic
    amidine is left for the ring path to reach through its ring.
    """
    out = []
    for n in molecule.atom_numbers:
        if molecule.element_of(n) != 6 or molecule.in_ring_of(n):
            continue
        if molecule.charge_of(n) or molecule.radical_of(n):
            continue
        h = molecule.implicit_h_of(n)
        if h is None:
            continue
        neighbors = tuple(molecule.neighbors_of(n))
        orders = [molecule.order_of(n, m) for m in neighbors]
        if any(o not in (1, 2) for o in orders) or orders.count(2) != 1:
            continue
        if len(neighbors) + h != 3:
            continue                       # sp2, and that one double bond is the only one it has
        sites = sorted(m for m in neighbors if _amidine_site(molecule, m, n))
        if len(sites) < 2 or not any(molecule.order_of(n, m) == 2 for m in sites):
            continue                       # the C=N has to be one this group is allowed to move
        out.append((n, sites))
    return out


def _choose_amidine(molecule: MoleculeContainer, carbon: int, sites: list[int],
                    ranks: dict[int, int]) -> int | None:
    """Which nitrogen takes the C=N, or `None` when it already has it.

    A neutral nitrogen with `d` heavy neighbours carries `3 - d` hydrogens single-bonded and `2 - d`
    double-bonded, so the total over the group is fixed whichever one accepts and there is nothing to
    prove admissible.  The key mirrors the ring one -- the sorted ranks of the nitrogens that KEEP their
    hydrogen, minimised -- so one canonical order answers both paths.
    """
    current = next(n for n in sites if molecule.order_of(carbon, n) == 2)
    acceptor = min(sites, key=lambda a: sorted(ranks[n] for n in sites if n != a))
    return None if acceptor == current else acceptor


# ---------------------------------------------------------------------------------------------------
# The pass.
# ---------------------------------------------------------------------------------------------------

def _place_rings(molecule: MoleculeContainer, lines: list[tuple[str, tuple[int, ...], str]]) -> bool:
    """Decide every ring system, aromatic as drawn or spelled aromatic for the purpose.  Did it move?"""
    work, systems = _aromatized(molecule)
    target = work if work is not None else molecule

    sites = _sites(target)
    if len(sites) < 2:
        return False                      # one site is one placement; zero is none

    # every read happens before the session opens: the container refuses a read while a journal is
    # pending, so the whole plan is computed first and the session is pure writes.
    ranks = _ranks(target, sites)
    rule = _RULE if work is None else _RULE_KEKULE
    plan: list[tuple[list[int], frozenset[int], frozenset[int]]] = []
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
        protons, anions = answer
        plan.append((group, protons, anions))
        lines.append((rule, tuple(group),
                      f'ring system {tuple(group)!r}: mobile hydrogen(s) placed on '
                      f'{tuple(sorted(protons))!r} and charge(s) on {tuple(sorted(anions))!r}, the '
                      f'canonical placement for this skeleton'))

    if not plan:
        return False

    if work is None:
        with molecule.edit():
            for group, protons, anions in plan:
                for n in group:
                    molecule.set_charge(n, -1 if n in anions else 0)
                    molecule.set_hydrogens(n, 1 if n in protons else 0)
                # no parity restore, deliberately: `set_hydrogens` and `set_charge` do not clear one,
                # only `delete_atom` does.  Pinned by
                # `test_a_stereocentre_is_not_touched_and_needs_no_parity_restore`.
        return True

    # The placement was decided on an aromatic spelling the caller never asked for, so the working copy
    # is kekulised and only the bonds this module itself spelled aromatic are read back.  A system whose
    # group did not move keeps the orders it was drawn with: the pass answers where a hydrogen goes, and
    # rewriting a Kekule form nobody asked about is not that answer.
    with work.edit():
        for group, protons, anions in plan:
            for n in group:
                work.set_charge(n, -1 if n in anions else 0)
                work.set_hydrogens(n, 1 if n in protons else 0)
    if work.kekule().unresolved:
        lines.append((_RULE_REFUSED, tuple(sorted(n for group, _, _ in plan for n in group)),
                      'the canonical placement has no Kekule form for the molecule as a whole, though '
                      'it had one for each system alone; nothing was written'))
        return False

    moved = frozenset(n for group, _, _ in plan for n in group)
    placed = {n: (1 if n in protons else 0, -1 if n in anions else 0)
              for group, protons, anions in plan for n in group}
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


def _place_amidines(molecule: MoleculeContainer,
                    lines: list[tuple[str, tuple[int, ...], str]]) -> bool:
    """Decide every acyclic amidine and guanidine.  Did anything move?"""
    groups = _amidines(molecule)
    if not groups:
        return False

    # the double bond is stripped along with the hydrogens: it is half of what the two spellings differ
    # by, so a frame that kept it would rank the two drawings differently.
    ranks = _ranks(molecule, [n for _, sites in groups for n in sites],
                   [(carbon, n) for carbon, sites in groups for n in sites])
    plan: list[tuple[int, list[int], int, dict[int, int]]] = []
    for carbon, sites in groups:
        acceptor = _choose_amidine(molecule, carbon, sites, ranks)
        if acceptor is None:
            continue
        hydrogens = {n: (2 if n == acceptor else 3) - len(tuple(molecule.neighbors_of(n)))
                     for n in sites}
        if min(hydrogens.values()) < 0:
            continue                      # no nitrogen may be asked for a hydrogen it does not have
        plan.append((carbon, sites, acceptor, hydrogens))
        lines.append((_RULE_AMIDINE, tuple(sites),
                      f'amidine at atom {carbon}: the double bond to {tuple(sites)!r} placed on '
                      f'{acceptor}, the canonical acceptor for this skeleton, and the hydrogens '
                      f'follow it'))

    if not plan:
        return False
    with molecule.edit():
        for carbon, sites, acceptor, hydrogens in plan:
            for n in sites:
                molecule.set_order(carbon, n, 2 if n == acceptor else 1)
            for n, h in hydrogens.items():
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
    # rings first: an amidine's gates exclude every ring atom, so the two passes are independent and the
    # order is a convenience rather than a dependency.
    changed = _place_rings(molecule, lines)
    changed = _place_amidines(molecule, lines) or changed

    with recording(molecule, stage='isomers') as log:
        for rule, atoms, message in lines:
            log.append(LogRecord(rule, atoms, message,
                                 LOST if rule in (_RULE_BUDGET, _RULE_REFUSED) else INFO))
    return changed
