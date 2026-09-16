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
"""Which Kekule form a compound is stored in: the one whose double bonds sit inside the small rings.

A compound with more than one Kekule form has an aromatic form that depends on which one it is stored
in.  `thiele()` reads the drawing it is handed and refuses a candidate ring whose atom holds its double
bond outside the ring -- the rule that tells p-benzoquinone from benzene.  A candidate ring fused to a
ring `thiele()` does not consider can have its fusion carbons' double bonds absorbed by that
neighbour, and is then refused for looking exocyclically doubled: dibenzo[a,e]cyclooctatetraene's
eight-ring and biphenylene's four-ring each do it in one of their five Kekule forms.

So the form is picked before `thiele()` runs, and `thiele()`'s rule stays what it says.  A ring is
spelled alternating on a copy and the rest of its conjugated component handed back to `kekule()`; the
form kept is the one with the most double bonds inside a ring of 5 to 7 atoms, the window `thiele()`
considers (`core/_thiele.pxi`).  The kekuliser writes every order outside the ring, and a trial is
accepted only when it comes back as the same molecule -- `kekule()` repairs by design, so a hydrogen it
dropped or a charge it added means the trial failed and not that the form is better.

A local optimum and not a proved global one: rings are filled one at a time, each round taking the best
trial by score and then by canonical bytes.  That is order-independent per round, so a ring that CAN be
filled is filled whichever form arrived; two forms of equal score are still two forms.

A RING OUTSIDE THAT WINDOW HAS NO AROMATIC FORM TO COLLAPSE ITS ALTERNATIONS INTO ONE, so the phase it
is drawn in survives into storage: 1,2-dimethylcyclooctatetraene's two bond-shift drawings are one
compound.  An already-alternating ring of eight atoms or more -- or of four -- is therefore spelled both
ways as well, and the phase kept is the one with the smaller canonical form.  The score has to come out
equal, so no small ring's double bond is ever traded for a spelling, and the bytes have to be strictly
smaller, so the choice terminates.

WHAT IS OFFERED A PHASE IS A RING OF `mol.rings`, THE SSSR, and not every cycle: naphthalene's two
six-rings are inside the window and its ten-atom perimeter is not in the set at all, so both of its
Kekule spellings are left exactly as drawn for `thiele()` to decide.  The stereo rule reads wider -- an
alternating cycle is walked wherever it runs, a fused perimeter included -- so a bond can be
non-stereogenic in a form this stage has nothing to offer, which is dibenzo[a,e]cyclooctatetraene once
its benzo rings hold their own double bonds.

A configuration pins a double bond only where the constitution makes it stereogenic.  A cis/trans parity
written on a bond an alternating cycle can move describes the phase and not the compound's geometry
(`core/_stereo.pxi`, `SU_SHIFTABLE`), and `canonicalize()` drops those with `validate_stereo()` before
this stage runs.  Every other stated configuration is carried, and a form that cannot carry one is a
trial that failed.
"""
from ..core import INFO, LOST, LogRecord, MoleculeContainer, recording
from ._implicit import check_valence


__all__ = ['standardize_kekule']


#: Table-qualified, as every rule id must be.  The budget is a `LOST` record, the fill and the phase
#: `INFO` ones: both forms of a ring either touches were valid molecules, so nothing here is a repair.
_RULE = 'kekule-form:ring-filled'
_RULE_PHASE = 'kekule-form:ring-phase'
_RULE_BUDGET = 'kekule-form:budget'

#: The ring sizes `thiele()` considers, `core/_thiele.pxi`.  A double bond inside one of these is worth
#: something to the aromatiser and one anywhere else is not, so nothing else scores -- and a ring outside
#: the window is the one whose alternation nothing downstream collapses, so it is the one the phase stage
#: owns.
_SIZE_MIN = 5
_SIZE_MAX = 7

#: Kekulisations per molecule.  Two per deficient ring per round, and the gate lets almost nothing
#: through -- 16 of 212 molecules in the tree's own corpus, none of them filled -- so the cap is what a
#: pathological fused system pays rather than a budget ordinary input reaches.
_TRIALS_MAX = 512


def _edges(ring: tuple) -> list[tuple[int, int]]:
    return list(zip(ring, ring[1:] + ring[:1]))


def _scored(molecule: MoleculeContainer) -> set[frozenset]:
    """The bonds a double bond scores on: inside a ring of `_SIZE_MIN` to `_SIZE_MAX` atoms."""
    out = set()
    for ring in molecule.rings:
        if _SIZE_MIN <= len(ring) <= _SIZE_MAX:
            out.update(frozenset(e) for e in _edges(ring))
    return out


def _score(molecule: MoleculeContainer, scored: set[frozenset]) -> int:
    return sum(1 for b in molecule.bonds() if b.order == 2 and frozenset((b.n, b.m)) in scored)


def _doubles(molecule: MoleculeContainer) -> dict[int, int]:
    """Double bonds per atom, which is what says whether a bond is one a matching may move."""
    out = dict.fromkeys(molecule.atom_numbers, 0)
    for b in molecule.bonds():
        if b.order == 2:
            out[b.n] += 1
            out[b.m] += 1
    return out


def _state(molecule: MoleculeContainer) -> dict[int, tuple]:
    """Everything a trial has to give back unchanged, since another compound is not another form."""
    return {n: (molecule.implicit_h_of(n), molecule.charge_of(n), molecule.radical_of(n))
            for n in molecule.atom_numbers}


def _stated(molecule: MoleculeContainer) -> frozenset[tuple]:
    """Every configuration this molecule states.  Parity 0 is nobody's assertion, so it is not one.

    A parity on a unit the constitution does not make stereogenic is not one either: a cis/trans sign
    written on a bond an alternating cycle can move says which phase the ring was drawn in, which is the
    very thing this pass chooses (`core/_stereo.pxi`, `SU_SHIFTABLE`).
    """
    return frozenset((u['kind'], u['anchor'], u['refs'], u['parity'])
                     for u in molecule.stereo_units() if u['parity'] and u['stereogenic'])


def _deficient(molecule: MoleculeContainer, doubles: dict[int, int]) -> list[tuple]:
    """The rings worth a trial: an even candidate ring holding fewer double bonds than it has room for.

    Every ring atom must hold exactly one double bond -- an atom holding none has nothing to move into
    the ring, and two is a cumulene rather than this shape -- and every ring bond must be Kekule, an
    aromatic one meaning `thiele()` has already taken the ring.  An odd ring is not tried at all: which
    of its atoms donates a lone pair instead of a double bond is a second choice, and this pass makes
    the one it can prove.
    """
    out = []
    for ring in molecule.rings:
        if len(ring) % 2 or not _SIZE_MIN <= len(ring) <= _SIZE_MAX:
            continue
        if any(doubles[n] != 1 for n in ring):
            continue
        orders = [molecule.order_of(u, v) for u, v in _edges(ring)]
        if any(o not in (1, 2) for o in orders):
            continue
        if orders.count(2) < len(ring) // 2:
            out.append(ring)
    return out


def _alternating(molecule: MoleculeContainer, doubles: dict[int, int]) -> list[tuple]:
    """The rings with only a phase left to choose: even, no candidate for `thiele()`, already alternating.

    The same shape `_deficient` wants -- every ring atom holding exactly one double bond, every ring bond
    Kekule -- with the room already full, which for a ring of `n` atoms is `n // 2` double bonds: one per
    atom and no atom in two of them is a perfect matching of the cycle, so full means alternating.  Sizes
    `thiele()` considers are left out, their alternations being the aromatiser's to collapse.

    A RING OF THE SSSR, so a perimeter alternation is not one of these: naphthalene drawn with its fusion
    bond single alternates around all ten atoms and is still left alone.
    """
    out = []
    for ring in molecule.rings:
        if len(ring) % 2 or _SIZE_MIN <= len(ring) <= _SIZE_MAX:
            continue
        if any(doubles[n] != 1 for n in ring):
            continue
        orders = [molecule.order_of(u, v) for u, v in _edges(ring)]
        if any(o not in (1, 2) for o in orders):
            continue
        if orders.count(2) == len(ring) // 2:
            out.append(ring)
    return out


def _component(molecule: MoleculeContainer, ring: tuple, doubles: dict[int, int]) -> list[tuple[int, int]]:
    """The Kekule bonds of the conjugated component holding `ring`, so a trial cannot reach past it.

    Walked over the atoms carrying a double bond, which are exactly the ones a Kekule form is a perfect
    matching of.  An amide's nitrogen stops the walk and its C=O is left alone.
    """
    seen = set(ring)
    stack = list(ring)
    while stack:
        v = stack.pop()
        for w in molecule.neighbors_of(v):
            if w in seen or not doubles[w] or molecule.order_of(v, w) not in (1, 2):
                continue
            seen.add(w)
            stack.append(w)
    return [(b.n, b.m) for b in molecule.bonds()
            if b.order in (1, 2) and b.n in seen and b.m in seen]


def _filled(molecule: MoleculeContainer, ring: tuple, phase: int, doubles: dict[int, int]):
    """`ring` spelled alternating and the rest of its component re-kekulised, or `None`.

    Both phases of an even ring are alternations, and they are different forms of the whole molecule:
    naphthalene's fusion bond carries a double in one of its three forms and not in the others.
    """
    edges = _edges(ring)
    fixed = {frozenset(e) for e in edges}
    aromatic = [(u, v) for u, v in _component(molecule, ring, doubles) if frozenset((u, v)) not in fixed]
    state = _state(molecule)
    stated = _stated(molecule)

    work = molecule.copy()
    with work.edit():
        for i, (u, v) in enumerate(edges):
            work.set_order(u, v, 2 if i % 2 == phase else 1)
    result = work.kekule(aromatic_bonds=aromatic)
    if result.unresolved or any(record.severity != INFO for record in result.log):
        return None                        # a relaxation fired: this is a repair and not a form
    if _state(work) != state or check_valence(work):
        return None
    if not stated <= _stated(work):
        return None                        # a configuration the form cannot carry
    return work


def standardize_kekule(molecule: MoleculeContainer) -> bool:
    """Store the Kekule form whose double bonds sit inside the rings.  Did any order move?

    THE STAGE THAT MAKES `thiele()` INDEPENDENT OF THE DRAWING.  A compound with several Kekule forms
    has one aromatic form per form otherwise: `C1=CC2=CC=C3C=CC=CC3=CC=C2C=C1` and
    `C1=CC=C2C=CC3=CC=CC=C3C=CC2=C1` are dibenzo[a,e]cyclooctatetraene twice, and only the second
    aromatises, its benzo rings holding the double bonds the first one lends to the eight-ring.

    A RING NO CANDIDATE SIZE COVERS IS GIVEN A CANONICAL PHASE INSTEAD, since nothing downstream will
    collapse its two alternations: `C/C1=C/C=C\\C=C/C=C\\1/C` and `C/C1=C(\\C)/C=C\\C=C/C=C\\1` are
    1,2-dimethylcyclooctatetraene's two bond-shift drawings and share a form once this has run.

    Run on a kekulised molecule and before `thiele()`, which is where `canonicalize()` runs it.  On an
    aromatic molecule there is nothing to choose and the answer is `False`: a ring `thiele()` has taken
    is a ring this pass leaves alone.

    Neither a repair nor a refusal: every form involved kekulises and holds the same atoms, hydrogens
    and charges, so `molecule.log` gets one `INFO` record per ring filled or phased -- and one `LOST`
    record when a fused system has more trials than the budget allows, which leaves the molecule as it
    arrived.  Never raises.
    """
    if not molecule.rings:
        return False                       # a form is a ring's to choose

    scored = _scored(molecule)
    lines: list[tuple[str, tuple[int, ...], str]] = []
    current = molecule
    score = _score(current, scored)
    filled: list[tuple] = []
    trials = 0
    budget = False
    # each round strictly raises the score, which is bounded by `len(scored)`, so this terminates
    while True:
        doubles = _doubles(current)
        pool = []
        for ring in _deficient(current, doubles):
            if trials + 2 > _TRIALS_MAX:
                budget = True
                lines.append((_RULE_BUDGET, tuple(sorted(ring)),
                              f'ring {tuple(ring)!r} and the ones after it were not tried: this '
                              f'molecule reached the {_TRIALS_MAX}-kekulisation budget, so it keeps '
                              f'whichever form it had reached, which two drawings of it may not share'))
                break
            for phase in (0, 1):
                trials += 1
                work = _filled(current, ring, phase, doubles)
                if work is None:
                    continue
                filled_score = _score(work, scored)
                if filled_score > score:
                    pool.append((filled_score, work.canonical_bytes, ring, work))
        if budget or not pool:
            break
        # the best score, then the smallest canonical form: a round is a choice over a set and must not
        # depend on the order the rings came in
        score, _, ring, current = min(pool, key=lambda t: (-t[0], t[1]))
        filled.append(ring)

    # then the phase, which the score cannot decide: the rings this reaches hold no bond it counts.  Each
    # round strictly lowers the canonical form, which is what terminates it -- the same ring can be
    # flipped twice only by coming back to a form already left behind.
    phased: list[tuple] = []
    while not budget:
        doubles = _doubles(current)
        rings = _alternating(current, doubles)
        if not rings:
            break                          # ahead of the canonical form, which is not a free read
        best = current.canonical_bytes
        pool = []
        for ring in rings:
            if trials + 2 > _TRIALS_MAX:
                budget = True
                lines.append((_RULE_BUDGET, tuple(sorted(ring)),
                              f'ring {tuple(ring)!r} and the ones after it were not phased: this '
                              f'molecule reached the {_TRIALS_MAX}-kekulisation budget, so it keeps '
                              f'whichever form it had reached, which two drawings of it may not share'))
                break
            for phase in (0, 1):
                trials += 1
                work = _filled(current, ring, phase, doubles)
                # equal score and nothing else: a phase is a tie-break and never a reason to give up a
                # double bond `thiele()` would have used
                if work is None or _score(work, scored) != score:
                    continue
                if work.canonical_bytes < best:
                    pool.append((work.canonical_bytes, ring, work))
        if budget or not pool:
            break
        _, ring, current = min(pool, key=lambda t: t[0])
        phased.append(ring)

    if current is molecule:
        if lines:
            with recording(molecule, stage='kekule-form') as log:
                for rule, atoms, message in lines:
                    log.append(LogRecord(rule, atoms, message, LOST))
        return False

    # every read first: the container refuses a read while a session is open.  Only the orders are
    # written -- each trial gave back the hydrogens, charges and radicals it was handed, so there is
    # nothing else the chosen form differs by.
    orders = [(b.n, b.m, b.order) for b in current.bonds() if molecule.order_of(b.n, b.m) != b.order]
    with molecule.edit():
        for u, v, order in orders:
            molecule.set_order(u, v, order)

    with recording(molecule, stage='kekule-form') as log:
        for ring in filled:
            log.append(LogRecord(_RULE, tuple(sorted(ring)),
                                 f'ring {tuple(ring)!r} was given the double bonds its atoms held '
                                 f'outside it, the Kekule form of this molecule that most of its '
                                 f'small rings can be aromatic in', INFO))
        for ring in phased:
            log.append(LogRecord(_RULE_PHASE, tuple(sorted(ring)),
                                 f'ring {tuple(ring)!r} had its double bonds shifted by one bond: no '
                                 f'aromatic form collapses a ring this size, so the alternation with '
                                 f'the smaller canonical form is the one stored', INFO))
        for rule, atoms, message in lines:
            log.append(LogRecord(rule, atoms, message, LOST))
    return True
