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
"""The rewrite-rule engine: match a repair rule, apply its patch, recompute what the patch invalidated.

A repair pass the caller asks for -- no reader or writer runs it.  The patch language is a signed
charge delta and an optional absolute radical flag per matched atom, plus a new order in {1,2,3,8}
per matched bond; no rule touches the atom or bond set, an isotope, an aromatic order or a hydrogen
count.  A patch is validated in full before any of it is written.
"""
from collections.abc import MutableSequence
from ._implicit import calc_implicit
from ._tables import Rule, groups_rules, metals_rules
from ..core import LogRecord, MoleculeContainer, recording


# `LogRecord` lives in `chython.core._log` because the SMIRKS patcher needs it and cannot import this
# package.  Re-exported, never re-defined: two NamedTuples with the same fields are two types.
__all__ = ['LogRecord', 'standardize']


# The core's charge span.  A patch that would leave an atom outside it is refused whole rather than
# clipped: clipping would silently change the formal charge the rule was written to produce.
_CHARGE_MIN, _CHARGE_MAX = -4, 8


def _apply(molecule: MoleculeContainer, rule: Rule, mapping: dict[int, int],
           log: MutableSequence) -> set[int] | None:
    """Apply one rule's patch at one match, or refuse it whole.  Returns the atoms written.

    All-or-nothing: every charge is computed and range-checked before the first write.  `None` means
    nothing was written and a log line says why.
    """
    written: set[int] = set()
    charges: list[tuple[int, int]] = []
    radicals: list[tuple[int, bool]] = []

    for number, delta, radical in rule.atom_fix:
        n = mapping[rule.numbers[number]]
        if delta:
            charge = molecule.charge_of(n) + delta
            if charge < _CHARGE_MIN or charge > _CHARGE_MAX:
                log.append(LogRecord(rule.id, tuple(sorted(mapping.values())),
                                     f'refused: atom {n} would take charge {charge}, outside '
                                     f'{_CHARGE_MIN}..{_CHARGE_MAX}; nothing was written'))
                return None
            charges.append((n, charge))
            written.add(n)
        if radical is not None and molecule.radical_of(n) != radical:
            radicals.append((n, radical))
            written.add(n)

    orders: list[tuple[int, int, int]] = []
    for a, b, order in rule.bonds_fix:
        u, v = mapping[rule.numbers[a]], mapping[rule.numbers[b]]
        if molecule.order_of(u, v) != order:
            orders.append((u, v, order))
            written.update((u, v))

    if not written:
        return None

    # one edit scope, so the derived words are rebuilt once rather than per write.
    with molecule.edit():
        for n, charge in charges:
            molecule.set_charge(n, charge)
        for n, radical in radicals:
            molecule.set_radical(n, radical)
        for u, v, order in orders:
            molecule.set_order(u, v, order)

    log.append(LogRecord(rule.id, tuple(sorted(mapping.values())), rule.why))
    return written


def _pass(molecule: MoleculeContainer, rules: tuple[Rule, ...],
          log: MutableSequence, fix_tautomers: bool = True) -> set[int]:
    """Run one rule table over `molecule`.  Returns every atom any patch wrote.

    The overlap policy is asymmetric on purpose: a match is rejected when any of its atoms has been
    seen, but what a match *contributes* to the seen set is only its site -- the matched atoms minus
    the rule's shared anchors.  Test the whole match, record the site.  Excluding anchors from the
    update is what lets several ligands share one metal (`Fe(CO)3` is three matches on one iron);
    excluding them from the test as well lets a match overlapping an already-repaired site by exactly
    its own wildcard slip through and undo it.  The set is per rule, not per table, so a rule listed
    twice can catch what the first copy's own dedupe rejected.  Cross-rule ordering is handled by the
    mutation itself: rules run in file order against the molecule as it now stands.
    """
    written: set[int] = set()
    for rule in rules:
        if not fix_tautomers and rule.tautomer:
            continue
        # the cheap screen: answers 'definitely not' from the feature words, with no embedding search
        if not rule.query.may_match(molecule):
            continue
        anchors = {rule.numbers[number] for number in rule.anchors}
        seen: set[int] = set()
        # materialised before the first patch: iterating lazily while mutating would search a graph
        # that is changing underneath the search
        for mapping in tuple(rule.query.get_mapping(molecule)):
            # the test is the whole match; the update below is the site.  See the docstring.
            if not set(mapping.values()).isdisjoint(seen):
                continue
            seen |= {molecule_atom for query_atom, molecule_atom in mapping.items()
                     if query_atom not in anchors}
            touched = _apply(molecule, rule, mapping, log)
            if touched:
                written |= touched
    return written


def standardize(molecule: MoleculeContainer, *, fix_hydrogens: bool = True,
                fix_tautomers: bool = True) -> bool:
    """Repair mis-drawn functional groups and metal-organic bonding in place.  Did anything change?

    Runs the functional-group rules and then the metal-organic ones, then recomputes the implicit
    hydrogen count of every atom a patch wrote -- charge, radical state and bond order all change what
    the valence collection gives an atom.  `fix_hydrogens=False` skips that recompute, for a caller about
    to kekulise anyway.  `molecule.log` takes a record per patch applied and per patch refused.

    `fix_tautomers=False` withholds the group rules whose repair displaces a hydrogen between heavy
    atoms; no rule writes a hydrogen count, so that is a `tautomer` column on the row rather than
    something the engine can infer.  Aromatic bonds are left alone: no valence row admits an order-4
    environment, so kekulise first or an aromatic atom whose charge a patch changed gets `H_UNKNOWN`.
    """
    with recording(molecule, stage='standardize') as lg:
        written = _pass(molecule, groups_rules(), lg, fix_tautomers)
        written |= _pass(molecule, metals_rules(), lg, fix_tautomers)
    if not written:
        return False
    if fix_hydrogens:
        for n in sorted(written):
            calc_implicit(molecule, n)
    return True
