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
"""Pair a biradical, or move a dipole's charge, along an alternating-bond path.

Endpoints come from `tables/resonance.tsv`.  Only orders, charges, radicals and derived hydrogen
counts are written; total formal charge is conserved, an aromatic bond is never crossed, and every
atom a patch touches is valence-checked.  Deterministic: everything is visited by atom number.
"""
from collections.abc import Iterator, Sequence
from ._implicit import calc_implicit
from ._standardize import LogRecord
from ._tables import Endpoint, resonance_rules_by_role
from ..core import MoleculeContainer, recording
from ..core._core import valence_has_rules, valence_implicit_h


__all__ = ['fix_resonance']


# A flip is one bond of a path: `(u, v, new order)`, in path order from the start.
Flip = tuple[int, int, int]

# Edge expansions one path search may make before giving up with a log line.
_WALK_BUDGET = 200_000

# The core's charge span, as in `_standardize.py`.  A patch that would leave an atom outside it is
# refused whole rather than clipped.
_CHARGE_MIN, _CHARGE_MAX = -4, 8


# --------------------------------------------------------------------------------------------------
# the walk.  Knows about bond orders and nothing else
# --------------------------------------------------------------------------------------------------

def alternating_paths(molecule: MoleculeContainer, start: int, targets: frozenset[int],
                      allowed: frozenset[int], *, odd_length_only: bool = False,
                      minimum_length: int = 1,
                      budget: int = _WALK_BUDGET) -> Iterator[list[Flip]]:
    """Yield every simple path from `start` into `targets` whose bond orders can alternate.

    A path is yielded as `[(u, v, new order), ...]`: the first bond goes up by one, the second down,
    and so on.  Only atoms in `allowed` are entered, only orders 1, 2 and 3 are crossed, and only new
    orders in 1..3 are produced -- so aromatic (4) and dative (8) bonds are impassable in the walk
    itself, not merely in whatever chose `allowed`.  Nothing here reads a charge, a hydrogen count or
    a valence rule.  `odd_length_only` keeps the paths whose last flip is an increase;
    `minimum_length` is the shortest path accepted.  Depth-first, so not shortest-first, but the
    order is a function of the atom numbering alone.
    """
    if start in targets:                     # a zero-length path moves nothing
        return
    path: list[Flip] = []
    seen = {start}
    # `stack` holds `(from, to, depth, new order)`.  Pushed in descending atom order so that `pop`
    # takes the lowest first -- this function's determinism is these two `sorted` calls.
    stack: list[tuple[int, int, int, int]] = []
    for m in sorted(molecule.neighbors_of(start), reverse=True):
        if m not in allowed:
            continue
        order = molecule.order_of(start, m)
        if 1 <= order <= 2:                  # 3 + 1 is not a bond order
            stack.append((start, m, 0, order + 1))
    steps = 0
    while stack:
        steps += 1
        if steps > budget:                   # pragma: no cover - needs a pathological input
            return
        last, current, depth, order = stack.pop()
        if len(path) > depth:                # backtracking: drop the tail this frame replaces
            seen.difference_update(m for _, m, _ in path[depth:])
            del path[depth:]
        path.append((last, current, order))

        if current in targets and len(path) >= minimum_length \
                and (not odd_length_only or len(path) % 2):
            yield list(path)
            continue                         # a target is an endpoint, never an interior atom

        seen.add(current)
        delta = -1 if (depth + 1) % 2 else 1
        outgoing = []
        for m in molecule.neighbors_of(current):
            if m in seen or m not in allowed:
                continue
            order = molecule.order_of(current, m)
            if not 1 <= order <= 3:          # aromatic (4) and dative (8) are impassable
                continue
            if 1 <= order + delta <= 3:
                outgoing.append((current, m, depth + 1, order + delta))
        stack.extend(sorted(outgoing, reverse=True))


# --------------------------------------------------------------------------------------------------
# endpoint classification.  The charge-aware half, and all of it table-driven
# --------------------------------------------------------------------------------------------------

def _matches(molecule: MoleculeContainer, endpoint: Endpoint) -> set[int]:
    """The anchor atom of every embedding of one row, or an empty set if the screen says no."""
    if not endpoint.query.may_match(molecule):
        return set()
    return {mapping[endpoint.anchor] for mapping in endpoint.query.get_mapping(molecule)}


def _role(molecule: MoleculeContainer, rows: Sequence[Endpoint]) -> dict[int, str]:
    """`{atom: the id of the first row that claimed it}` over one role's rows, in file order."""
    out: dict[int, str] = {}
    for endpoint in rows:
        for n in _matches(molecule, endpoint):
            out.setdefault(n, endpoint.id)
    return out


def _carries_aromatic(molecule: MoleculeContainer, n: int) -> bool:
    return any(molecule.order_of(n, m) == 4 for m in molecule.neighbors_of(n))


class _Endpoints:
    """One molecule's classification: which atoms a path may use, and who may end one."""
    __slots__ = ('allowed', 'radicals', 'donors', 'acceptors', 'why')

    def __init__(self, allowed: frozenset[int], radicals: list[int], donors: list[int],
                 acceptors: list[int], why: dict[int, str]):
        self.allowed = allowed
        self.radicals = radicals
        self.donors = donors
        self.acceptors = acceptors
        self.why = why                       # atom -> id of the row that accepted it, for the log


def _classify(molecule: MoleculeContainer, refuse) -> _Endpoints:
    """Run the table over `molecule` and hand back the endpoint sets, logging every veto.

    A veto is logged only when the atom ALSO matched an accept row for the same role.  An atom no
    accept row wanted is not being refused, it is simply not an endpoint, and a record for it would
    bury the ones that mean something.
    """
    rows = resonance_rules_by_role()

    allowed = set()
    for endpoint in rows['path']:
        allowed |= _matches(molecule, endpoint)
    # an atom carrying an aromatic bond leaves the walk entirely: the valence collection has no
    # aromatic row, so there is nothing to check it against.  Kekulise first.
    allowed = frozenset(n for n in allowed if not _carries_aromatic(molecule, n))

    why: dict[int, str] = {}
    accepted: dict[str, list[int]] = {}
    for role, veto_role in (('radical', None), ('donor', 'veto_donor'),
                            ('acceptor', 'veto_acceptor')):
        claimed = _role(molecule, rows[role])
        vetoed = _role(molecule, rows[veto_role]) if veto_role else {}
        keep = []
        for n in sorted(claimed):
            if n in vetoed:
                refuse(vetoed[n], (n,),
                       f'atom {n} matched {claimed[n]} but is vetoed by {vetoed[n]}: '
                       f'{rows_comment(rows[veto_role], vetoed[n])}')
                continue
            if n not in allowed:
                refuse(claimed[n], (n,),
                       f'atom {n} matched {claimed[n]} but carries an aromatic bond or is outside '
                       'the organic set, so no valence row describes it; kekulise first')
                continue
            if molecule.implicit_h_of(n) is None:
                refuse(claimed[n], (n,),
                       f'atom {n} matched {claimed[n]} but its record does not state a hydrogen '
                       'count, so the state after a patch cannot be derived')
                continue
            keep.append(n)
            why[n] = claimed[n]
        accepted[role] = keep

    return _Endpoints(allowed, accepted['radical'], accepted['donor'], accepted['acceptor'], why)


def rows_comment(rows: Sequence[Endpoint], row_id: str) -> str:
    """The table's own words for one row, so a log line quotes the chemistry and not just an id."""
    for endpoint in rows:
        if endpoint.id == row_id:
            return endpoint.comment
    return ''                                                                # pragma: no cover


# --------------------------------------------------------------------------------------------------
# validating and writing one patch
# --------------------------------------------------------------------------------------------------

def _validate(molecule: MoleculeContainer, flips: Sequence[Flip], charges: dict[int, int],
              radicals: dict[int, bool]) -> str | None:
    """None when every atom the patch touches keeps a describable valence, else why not.

    Every atom, both ends and every interior one.  An interior atom's charge, radical state and
    bond-order sum are invariant along a path, so only its environment changes -- which is why
    `valence_implicit_h`, and not the coarse `valence_has_rules`, is the last word.
    """
    new_order: dict[tuple[int, int], int] = {}
    for u, v, order in flips:
        new_order[(u, v)] = new_order[(v, u)] = order
    touched = {u for u, _, _ in flips} | {v for _, v, _ in flips}

    for n in sorted(touched):
        order_sum = 0
        environment: list[tuple[int, int]] = []
        for m in molecule.neighbors_of(n):
            order = new_order.get((n, m), molecule.order_of(n, m))
            if order == 8:                   # dative: outside valence bookkeeping, as elsewhere
                continue
            if order == 4:                   # pragma: no cover - `allowed` excludes these atoms
                return f'atom {n} carries an aromatic bond, which no valence row describes'
            order_sum += order
            environment.append((order, molecule.element_of(m)))

        charge = molecule.charge_of(n) + charges.get(n, 0)
        if charge < _CHARGE_MIN or charge > _CHARGE_MAX:
            return (f'atom {n} would take charge {charge}, outside '
                    f'{_CHARGE_MIN}..{_CHARGE_MAX}')
        radical = radicals.get(n, molecule.radical_of(n))
        element = molecule.element_of(n)
        if not valence_has_rules(element, charge, radical, order_sum):
            return (f'atom {n} would become atomic number {element} with charge {charge}, '
                    f'{"a radical" if radical else "no radical"} and a bond order sum of '
                    f'{order_sum}, which the valence collection has no rule for at all')
        if valence_implicit_h(element, charge, radical, order_sum, environment) is None:
            return (f'atom {n} would become atomic number {element} with charge {charge} and a bond '
                    f'order sum of {order_sum} in the environment {tuple(environment)}, which no '
                    'valence row admits')
    return None


def _write(molecule: MoleculeContainer, flips: Sequence[Flip], charges: dict[int, int],
           radicals: dict[int, bool]) -> set[int]:
    """Write one validated patch in a single edit scope and re-derive the hydrogens.  All or none.

    Every read happens before the scope opens: a container with a pending journal refuses to be read
    from.  One scope, so the arena rebuilds its derived words and re-bases the stereo parities once
    for the whole patch rather than once per bond.
    """
    absolute = {n: molecule.charge_of(n) + delta for n, delta in charges.items()}
    with molecule.edit():
        for n, charge in absolute.items():
            molecule.set_charge(n, charge)
        for n, radical in radicals.items():
            molecule.set_radical(n, radical)
        for u, v, order in flips:
            molecule.set_order(u, v, order)
    touched = {u for u, _, _ in flips} | {v for _, v, _ in flips}
    touched |= set(absolute) | set(radicals)
    for n in sorted(touched):
        calc_implicit(molecule, n)
    return touched


def _charge_potential(molecule: MoleculeContainer, n: int, delta: int) -> tuple[int, int]:
    """`(is charged, is charged and not nitrogen)` for atom `n` after `delta`.

    Summed over a patch's two ends and compared lexicographically, this potential must strictly
    decrease for a dipole patch to be accepted; that is what makes the pass terminate and be
    idempotent.
    """
    charge = molecule.charge_of(n) + delta
    if not charge:
        return 0, 0
    return 1, 0 if molecule.element_of(n) == 7 else 1


# --------------------------------------------------------------------------------------------------
# the pass
# --------------------------------------------------------------------------------------------------

def _sweep(molecule: MoleculeContainer, refuse) -> set[int]:
    """One classification and one pass over its endpoints.  Returns the atoms written."""
    endpoints = _classify(molecule, refuse)
    written: set[int] = set()

    # radicals first: a pair joined by an odd-length alternating path each gain one bond order and go
    # closed-shell.  A path of length 1 is the point here (`[CH2][CH2]` is ethylene), which is why
    # `minimum_length` differs between the two loops.
    remaining = list(endpoints.radicals)
    while len(remaining) > 1:
        n = remaining.pop(0)
        for flips in alternating_paths(molecule, n, frozenset(remaining), endpoints.allowed,
                                       odd_length_only=True):
            end = flips[-1][1]
            radicals = {n: False, end: False}
            why = _validate(molecule, flips, {}, radicals)
            if why is not None:
                refuse(endpoints.why[n], _atoms_of(flips), f'refused: {why}')
                continue
            written |= _write(molecule, flips, {}, radicals)
            remaining.remove(end)
            _applied(refuse, endpoints.why[n], _atoms_of(flips),
                     f'paired the radicals on atoms {n} and {end} along {_render(flips)}')
            break
        else:
            refuse(endpoints.why[n], (n,),
                   f'refused: no alternating path of odd length from radical atom {n} to another '
                   'radical crosses only single, double and triple bonds between non-aromatic atoms')

    # then dipoles.  `minimum_length=2`: a donor bonded straight to an acceptor is a charge
    # annihilation, not the double-bond transfer this walk is for.
    acceptors = list(endpoints.acceptors)
    for n in endpoints.donors:
        if not acceptors:
            break
        found = False
        for flips in alternating_paths(molecule, n, frozenset(acceptors), endpoints.allowed,
                                       minimum_length=2):
            end = flips[-1][1]
            before = _charge_potential(molecule, n, 0)
            after = _charge_potential(molecule, n, 1)
            before = (before[0] + _charge_potential(molecule, end, 0)[0],
                      before[1] + _charge_potential(molecule, end, 0)[1])
            after = (after[0] + _charge_potential(molecule, end, -1)[0],
                     after[1] + _charge_potential(molecule, end, -1)[1])
            if after >= before:
                refuse(endpoints.why[n], _atoms_of(flips),
                       f'refused: moving charge from atom {n} to atom {end} would not reduce '
                       f'(charged atoms, charges off nitrogen) from {before}, so the pass would not '
                       'terminate')
                continue
            charges = {n: 1, end: -1}
            why = _validate(molecule, flips, charges, {})
            if why is not None:
                refuse(endpoints.why[n], _atoms_of(flips), f'refused: {why}')
                continue
            written |= _write(molecule, flips, charges, {})
            acceptors.remove(end)
            found = True
            _applied(refuse, endpoints.why[n], _atoms_of(flips),
                     f'moved one unit of charge from atom {n} to atom {end} along {_render(flips)}')
            break
        if not found:
            continue
    return written


def _atoms_of(flips: Sequence[Flip]) -> tuple[int, ...]:
    return tuple(sorted({u for u, _, _ in flips} | {v for _, v, _ in flips}))


def _render(flips: Sequence[Flip]) -> str:
    return ', '.join(f'{u}-{v} -> order {order}' for u, v, order in flips)


def _applied(refuse, rule: str, atoms: tuple[int, ...], message: str) -> None:
    """Same sink as a refusal.  One channel, so a report reads in the order things happened."""
    refuse(rule, atoms, message)


def fix_resonance(molecule: MoleculeContainer) -> bool:
    """Pair biradicals and move dipole charges into a neutral form, in place.  Did it change?

    A repair the caller asks for: nothing in the library calls this.  `molecule.log` gets a record per
    patch applied and per patch refused, each naming the `resonance.tsv` row that fired or vetoed the
    endpoint; a refusal is never an exception.  Runs to a fixed point -- each patch strictly decreases
    `_charge_potential`, so a second call is a no-op.
    """
    seen: set[LogRecord] = set()
    records: list[LogRecord] = []

    def refuse(rule: str, atoms: tuple[int, ...], message: str) -> None:
        record = LogRecord(rule, atoms, message)
        if record not in seen:               # a fixed point revisits its refusals; say each once
            seen.add(record)
            records.append(record)

    written: set[int] = set()
    # at most one round per atom.  Belt and braces -- a round that writes nothing ends the loop -- so
    # that a future rule which forgets to decrease the potential hangs a test and not the process.
    for _ in range(molecule.atom_count + 1):
        touched = _sweep(molecule, refuse)
        if not touched:
            break
        written |= touched

    with recording(molecule, stage='resonance') as log:
        log.extend(records)
    return bool(written)
