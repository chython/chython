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
"""Bond orders for a ligand whose connectivity is known and whose orders are not.

A stated hydrogen count pins an atom's order sum and is the only thing that can force a multiple bond; an
`H_UNKNOWN` count is slack.  Assignment is all-or-nothing per connected fragment of the open subgraph,
deterministic, and writes bond orders only, in Kekule form -- never a charge, radical, count or bond.
"""
from collections import defaultdict
from collections.abc import Sequence
from typing import NamedTuple
from ._implicit import environment_of
from ..core import INFO, LOST, LogRecord, MoleculeContainer, REFUSED, recording
from ..core._core import valence_check, valence_implicit_h, valence_rules


__all__ = ['saturate']


#: One id per outcome a caller filters on; see the table in `saturate`'s docstring.
_RULE = 'saturate:orders'
_RULE_AMBIGUOUS = 'saturate:ambiguous'
_RULE_UNFORCED = 'saturate:unforced'
_RULE_NO_STATE = 'saturate:no-valence-state'
_RULE_GAP = 'saturate:collection-gap'
_RULE_AROMATIC = 'saturate:aromatic-bond'
_RULE_OVERSIZED = 'saturate:oversized'
_RULE_BUDGET = 'saturate:budget'
_RULE_UNSATISFIED = 'saturate:unsatisfied'

#: The highest order this pass will write.  Order 4 is aromatic and order 8 is dative; neither has a
#: valence row, so neither is derivable here.
_ORDER_MAX = 3

#: Open bonds per fragment.  Bounds the search, not the input: a fragment nothing pins is answered
#: before this cap is consulted, at any size.  The recursion below is one frame per bond.
_BONDS_MAX = 256

#: Search nodes per fragment.  Reaching it is a refusal with a record, never a truncated answer.
_NODES_MAX = 100_000

#: `(atomic number, charge, radical) -> {stated hydrogens or None: (order sums, ...)}`, built once
#: from the compiled valence collection.  `None` collects every order sum at the key regardless of
#: hydrogen count, which is what an `H_UNKNOWN` atom is allowed to take.
_STATES_CACHE: dict[tuple[int, int, bool], dict[int | None, tuple[int, ...]]] = {}


def _states() -> dict[tuple[int, int, bool], dict[int | None, tuple[int, ...]]]:
    """The valence collection indexed the way this pass asks it questions.

    The environment column is dropped here, because whether a row's environment is satisfied depends
    on the assignment being searched for.  So this is the optimistic filter and `_accepts` is the
    exact one, asked of a complete assignment inside the search and before anything is written.
    """
    if not _STATES_CACHE:
        collected: dict[tuple[int, int, bool], dict[int | None, set[int]]] = defaultdict(
            lambda: defaultdict(set))
        for z, charge, radical, bonds, hydrogens, _ in valence_rules():
            at = collected[(z, charge, radical)]
            at[hydrogens].add(bonds)
            at[None].add(bonds)
        for key, by_h in collected.items():
            _STATES_CACHE[key] = {h: tuple(sorted(bonds)) for h, bonds in by_h.items()}
    return _STATES_CACHE


class _State(NamedTuple):
    """Everything about one atom that the valence collection has to be asked, read once in phase 1.

    Read for every atom, frozen ones included: an atom this pass will not touch is still a neighbour.
    No field is ever written, which is why phase 6 may take them from here.
    """
    element: int
    charge: int
    radical: bool
    hydrogens: int | None                      #: None is `H_UNKNOWN`: slack, not zero


def _accepts(state: _State, order_sum: int, environment: Sequence[tuple[int, int]]) -> bool:
    """Does the collection accept this atom at this order sum in this neighbourhood?  Exactly.

    The one exact test, asked by both the search and the verdict -- accepting an assignment on the
    optimistic filter alone would let the verdict call an already-written atom a violation.  A stated
    hydrogen count is a complete question for `valence_check`; an `H_UNKNOWN` one becomes "does any
    row accept this order sum here", which is `valence_implicit_h`.
    """
    if state.hydrogens is None:
        return valence_implicit_h(state.element, state.charge, state.radical, order_sum,
                                  environment) is not None
    return valence_check(state.element, state.charge, state.radical, order_sum, state.hydrogens,
                         environment) == 'valid'


class _Site:
    """One connected fragment of the open subgraph, and its answer.

    A plain object rather than a tuple because it is built in three phases -- collected, solved,
    applied -- and a name on each field is what makes the phases readable.
    """
    __slots__ = ('atoms', 'bonds', 'solution', 'note')

    def __init__(self, atoms: tuple[int, ...], bonds: tuple[tuple[int, int, int], ...]):
        self.atoms = atoms
        #: `(low stable id, high stable id, headroom)`, sorted.  The sort is load-bearing, not
        #: tidiness: it is the search's tie-break, so walking the open bonds unsorted changes the
        #: answer on aromatic fragments.
        self.bonds = bonds
        self.solution: tuple[int, ...] | None = None
        #: `(rule, message, severity)` for whatever the search could not do, or could not promise
        self.note: tuple[str, str, str] | None = None


def _solve(site: _Site, needs: dict[int, tuple[int, ...]], state: dict[int, _State],
           neighbours: dict[int, list[tuple[int, int]]]) -> int:
    """Assign every open bond in one fragment, or refuse the fragment whole.  Returns nodes used.

    Depth-first over the fragment's bonds in sorted order, smallest increment first, with forward
    checking on both endpoints.  That makes the first solution the one that invents least
    unsaturation, and makes the enumeration a pure function of the stable ids.  The search stops at
    the second solution: uniqueness is a yes/no question.  A complete assignment is accepted only by
    the exact test `_accepts`, so the fragment can be refused whole instead of written and then
    complained about.
    """
    positions = {n: i for i, n in enumerate(site.atoms)}
    allowed = [frozenset(needs[n]) for n in site.atoms]
    #: per atom, the increment its still-unassigned open bonds could yet supply
    rest = [0] * len(site.atoms)
    for u, v, headroom in site.bonds:
        rest[positions[u]] += headroom
        rest[positions[v]] += headroom
    got = [0] * len(site.atoms)
    assigned = [0] * len(site.bonds)
    solutions: list[tuple[int, ...]] = []
    nodes = 0

    # What `_accepts` needs, per atom of the fragment in `site.atoms` order: the closed bonds'
    # contribution, which the search cannot change, and the open bonds as
    # `(position in the search, current order, neighbour's element)`.  Split once, not per candidate.
    at_depth = {(u, v): depth for depth, (u, v, _) in enumerate(site.bonds)}
    closed_sum: list[int] = []
    closed_env: list[tuple[tuple[int, int], ...]] = []
    open_env: list[tuple[tuple[int, int, int], ...]] = []
    for n in site.atoms:
        fixed_sum = 0
        fixed: list[tuple[int, int]] = []
        movable: list[tuple[int, int, int]] = []
        for m, order in neighbours[n]:
            depth = at_depth.get((min(n, m), max(n, m)))
            if depth is None:
                fixed_sum += order
                fixed.append((order, state[m].element))
            else:
                movable.append((depth, order, state[m].element))
        closed_sum.append(fixed_sum)
        closed_env.append(tuple(fixed))
        open_env.append(tuple(movable))

    def feasible(i: int) -> bool:
        low = got[i]
        return any(low <= want <= low + rest[i] for want in allowed[i])

    def exact() -> bool:
        """The collection asked about the assignment now on the table, environment column included."""
        for i, n in enumerate(site.atoms):
            order_sum = closed_sum[i]
            environment = list(closed_env[i])
            for depth, order, element in open_env[i]:
                order_sum += order + assigned[depth]
                environment.append((order + assigned[depth], element))
            if not _accepts(state[n], order_sum, environment):
                return False
        return True

    def walk(depth: int) -> None:
        nonlocal nodes
        nodes += 1
        if depth == len(site.bonds):
            # the environment half cannot be checked any earlier: a partial assignment has no
            # complete neighbourhood to ask about.
            if all(got[i] in allowed[i] for i in range(len(site.atoms))) and exact():
                solutions.append(tuple(assigned))
            return
        u, v, headroom = site.bonds[depth]
        i, j = positions[u], positions[v]
        rest[i] -= headroom
        rest[j] -= headroom
        for extra in range(headroom + 1):        # smallest first: never invent unsaturation
            got[i] += extra
            got[j] += extra
            if feasible(i) and feasible(j):
                assigned[depth] = extra
                walk(depth + 1)
            got[i] -= extra
            got[j] -= extra
            if len(solutions) > 1 or nodes > _NODES_MAX:
                break
        assigned[depth] = 0
        rest[i] += headroom
        rest[j] += headroom

    walk(0)
    if solutions:
        site.solution = solutions[0]
        if len(solutions) > 1:
            differing = sorted({n for k, (u, v, _) in enumerate(site.bonds)
                                if solutions[0][k] != solutions[1][k] for n in (u, v)})
            site.note = (_RULE_AMBIGUOUS,
                         f'atoms {tuple(differing)!r}: more than one assignment satisfies every '
                         'stated hydrogen count, and the one reported is the one that raises the '
                         'earliest bond least; an aromatic ring is ambiguous this way by '
                         'construction, so run thiele() if the Kekule choice is what differs, and '
                         'look at the fragment by hand if it is not', LOST)
        elif nodes > _NODES_MAX:
            # LOST, not REFUSED: this branch writes the assignment it found and only gives up on
            # proving it unique.  Its sibling below writes nothing, so the two must stay distinct.
            site.note = (_RULE_BUDGET,
                         f'atoms {site.atoms!r}: the search ran out of budget after {_NODES_MAX} '
                         'steps while checking whether the assignment is unique, so it is reported '
                         'as assigned but not as the only answer', LOST)
    elif nodes > _NODES_MAX:
        site.note = (_RULE_BUDGET,
                     f'atoms {site.atoms!r}: no assignment found within {_NODES_MAX} search steps; '
                     'every bond in the fragment is left as it was.  This is beyond the ligand this '
                     'pass is for', REFUSED)
    else:
        site.note = (_RULE_NO_STATE,
                     f'atoms {site.atoms!r}: no assignment of double and triple bonds puts every '
                     'atom of this fragment in a state the collection accepts, given the stated '
                     'hydrogen counts and charges and each atom\'s neighbourhood; every bond in it '
                     'is left as it was', REFUSED)
    return nodes


def saturate(molecule: MoleculeContainer) -> bool:
    """Raise bonds to double and triple until every atom's valence is satisfied.  Deterministically.

    Returns `True` when every atom ended in a state the valence collection calls valid -- success,
    not "changed", unlike `standardize()`.  An already fully ordered molecule returns `True` and
    writes nothing.  Everything not done is a `LogRecord` naming the atoms:

    ==============================  ============  =================================================
    rule                            severity      what it says
    ==============================  ============  =================================================
    ``saturate:orders``             info          these bonds were raised
    ``saturate:unforced``           info          nothing in this fragment demanded a multiple bond
    ``saturate:ambiguous``          lost          the answer is not unique; one of them is reported
    ``saturate:no-valence-state``   refused       this atom, or this fragment, has no assignment
    ``saturate:collection-gap``     lost          the collection describes nothing for this state
    ``saturate:aromatic-bond``      lost          an order-4 bond is here already; kekulise first
    ``saturate:oversized``          refused       the fragment's search is bigger than a ligand's
    ``saturate:budget``             lost/refused  the search was cut off: lost with an answer
                                                  written, refused with none
    ``saturate:unsatisfied``        lost          the finished atom is in a state no row accepts
    ==============================  ============  =================================================

    Intended for a ligand -- tens of atoms; amino acids and nucleotides are a template lookup.  Never
    alters a charge, a radical, a hydrogen count or the set of bonds, and never lowers an order.
    """
    states = _states()
    lines: list[tuple[str, tuple[int, ...], str, str]] = []

    # ------------------------------------------------------------------ phase 1: read, pure
    # The container refuses these reads once an edit session is open, so collect everything first.
    # Sorted, because the answer must not depend on arena order.
    atoms = sorted(molecule.atom_numbers)
    neighbours: dict[int, list[tuple[int, int]]] = {}
    state: dict[int, _State] = {}
    needs: dict[int, tuple[int, ...]] = {}
    frozen: set[int] = set()

    for n in atoms:
        order_sum, _, aromatic = environment_of(molecule, n)
        # for every atom, before any branch: a frozen atom is not searched but is still somebody's
        # neighbour, and the exact test reads the neighbours' elements
        state[n] = _State(molecule.element_of(n), molecule.charge_of(n), molecule.radical_of(n),
                          molecule.implicit_h_of(n))
        # order 8 contributes nothing to a valence and order 4 has no row; `environment_of` is this
        # package's one statement of both policies
        neighbours[n] = sorted((m, molecule.order_of(n, m)) for m in molecule.neighbors_of(n)
                               if 1 <= molecule.order_of(n, m) <= _ORDER_MAX)
        if aromatic:
            frozen.add(n)
            lines.append((_RULE_AROMATIC, (n,),
                          f'atom {n} carries {aromatic} aromatic bond(s), which no valence row '
                          'admits, so its bonds are left alone; kekulise before saturating', LOST))
            continue
        key = (state[n].element, state[n].charge, state[n].radical)
        if key not in states:
            frozen.add(n)
            lines.append((_RULE_GAP, (n,),
                          f'atom {n}: the valence collection describes no state for '
                          f'{molecule.atom(n).atomic_symbol} in charge {key[1]}'
                          f'{" as a radical" if key[2] else ""}, so its bonds are left alone -- a '
                          'gap in the collection, not a claim about the molecule', LOST))
            continue
        hydrogens = state[n].hydrogens                 # None is H_UNKNOWN: slack, not zero
        reachable = tuple(v - order_sum for v in states[key].get(hydrogens, ())
                          if v >= order_sum)
        if not reachable:
            frozen.add(n)
            lines.append((_RULE_NO_STATE, (n,),
                          f'atom {n}: no valence row accepts {molecule.atom(n).atomic_symbol} in '
                          f'charge {key[1]} with bond order sum {order_sum} or more and '
                          + (f'{hydrogens} hydrogen(s)' if hydrogens is not None
                             else 'any hydrogen count')
                          + ', so its bonds are left alone', REFUSED))
            continue
        needs[n] = reachable

    # ------------------------------------------------------------------ phase 2: propagate
    # A bond is OPEN while both its ends can still accept order.  Closing one lowers what its
    # neighbours can be handed, which can settle or starve them in turn, so this is a fixpoint rather
    # than a single sweep.  Atoms in sorted order, so the fixpoint is reached identically every run.
    open_bonds: set[tuple[int, int]] = set()
    for n in atoms:
        if n in frozen:
            continue
        for m, order in neighbours[n]:
            if m in frozen or n > m:
                continue
            if order < _ORDER_MAX:
                open_bonds.add((n, m))

    def headroom_of(u: int, v: int) -> int:
        return _ORDER_MAX - molecule.order_of(u, v)

    changed = True
    while changed:
        changed = False
        for n in atoms:
            if n in frozen:
                continue
            incident = [(u, v) for u, v in ((min(n, m), max(n, m)) for m, _ in neighbours[n])
                        if (u, v) in open_bonds]
            capacity = sum(headroom_of(u, v) for u, v in incident)
            reachable = tuple(v for v in needs[n] if v <= capacity)
            if not reachable:
                # starved: its neighbourhood cannot supply any state the collection accepts
                frozen.add(n)
                del needs[n]
                open_bonds.difference_update(incident)
                lines.append((_RULE_NO_STATE, (n,),
                              f'atom {n}: every valence state left to it needs more bond order than '
                              'its neighbours can accept, so its bonds are left as they were',
                              REFUSED))
                changed = True
                continue
            if reachable != needs[n]:
                needs[n] = reachable
                changed = True
            if not max(reachable) and incident:
                # settled at zero: nothing may be raised here, and closing its bonds is what lets
                # the rest fall apart into independent fragments
                open_bonds.difference_update(incident)
                changed = True

    # ------------------------------------------------------------------ phase 3: fragments
    # Sorted once and read from here on: `open_bonds` is the only set whose layout could otherwise
    # reach the answer, and no loop below may read an unordered container.
    opened = sorted(open_bonds)
    adjacency: dict[int, list[int]] = defaultdict(list)
    for u, v in opened:
        adjacency[u].append(v)
        adjacency[v].append(u)

    sites: list[_Site] = []
    seen: set[int] = set()
    for start in atoms:                                # sorted seeds, so fragment order is stable
        if start in seen or start not in adjacency:
            continue
        stack = [start]
        seen.add(start)
        members = {start}
        while stack:
            n = stack.pop()
            for m in sorted(adjacency[n]):
                if m not in members:
                    members.add(m)
                    seen.add(m)
                    stack.append(m)
        member_atoms = tuple(sorted(members))
        bonds = tuple((u, v, headroom_of(u, v)) for u, v in opened if u in members)
        sites.append(_Site(member_atoms, bonds))

    # ------------------------------------------------------------------ phase 4: solve each fragment
    for site in sites:
        if not any(min(needs[n]) for n in site.atoms):
            # nothing here demands a raise, so there is nothing for a search to satisfy and picking an
            # assignment would be choosing a compound rather than deriving it.  Answered before the
            # size cap deliberately: it needs no search and is right at any size.
            lines.append((_RULE_UNFORCED, site.atoms,
                          f'atoms {site.atoms!r}: no stated hydrogen count in this fragment demands '
                          'a multiple bond, so every bond in it stays single -- with no hydrogen '
                          'counts and no coordinates there is nothing to derive an order from',
                          INFO))
            continue
        if len(site.bonds) > _BONDS_MAX:
            lines.append((_RULE_OVERSIZED, site.atoms,
                          f'{len(site.bonds)} open bonds in one fragment is past the {_BONDS_MAX} '
                          'this pass accepts: it is for a ligand of tens of atoms, and a chain of '
                          'residues is a template lookup rather than a search.  Every bond in the '
                          'fragment is left as it was', REFUSED))
            continue
        _solve(site, needs, state, neighbours)
        if site.note is not None:
            rule, message, severity = site.note
            lines.append((rule, site.atoms, message, severity))

    # ------------------------------------------------------------------ phase 5: write, once
    writes: list[tuple[int, int, int]] = []
    for site in sites:
        if site.solution is None:
            continue
        raised = [(u, v, molecule.order_of(u, v) + extra)
                  for (u, v, _), extra in zip(site.bonds, site.solution) if extra]
        if raised:
            writes.extend(raised)
            lines.append((_RULE, tuple(sorted({n for u, v, _ in raised for n in (u, v)})),
                          'raised ' + ', '.join(f'{u}-{v} to order {order}'
                                                for u, v, order in raised), INFO))
    if writes:
        with molecule.edit():
            for u, v, order in writes:
                molecule.set_order(u, v, order)

    # ------------------------------------------------------------------ phase 6: verdict
    # Asked of the finished molecule and of every atom, not only the touched ones: an untouched atom
    # can still be the one the file got wrong.  Since `_accepts` is the same test the search used, a
    # fragment this pass wrote cannot be reported here -- only what it declined to assign.
    satisfied = True
    for n in atoms:
        order_sum, environment, aromatic = environment_of(molecule, n)
        if aromatic or n in frozen:
            satisfied = False                          # already reported above, with its reason
            continue
        if _accepts(state[n], order_sum, environment):
            continue
        satisfied = False
        hydrogens = state[n].hydrogens
        # the word comes from `valence_check`, which distinguishes a hole in the collection from a
        # claim about the molecule.  A slack atom has no count to check, and phase 1 would have frozen
        # it had its element been undescribed, so a violation is what is left.
        verdict = ('violation' if hydrogens is None
                   else valence_check(state[n].element, state[n].charge, state[n].radical,
                                      order_sum, hydrogens, environment))
        lines.append((_RULE_UNSATISFIED, (n,),
                      f'atom {n} ({molecule.atom(n).atomic_symbol}) ends with bond order sum '
                      f'{order_sum} and '
                      + ('an unstated hydrogen count' if hydrogens is None
                         else f'{hydrogens} hydrogen(s)')
                      + f', which the collection calls a {verdict}', LOST))

    with recording(molecule, stage='saturate') as log:
        for rule, touched, message, severity in lines:
            log.append(LogRecord(rule, touched, message, severity))
    return satisfied
