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
"""Moving hydrogens between the graph and the count.

`implicify_hydrogens` folds an ordinary hydrogen atom into its neighbour's implicit count and is on the
deduplication path: `[CH4]` and `[H]C([H])([H])[H]` do not hash equal until it has run.
`explicify_hydrogens` is the other direction and is not on that path -- it is for consumers that need
every hydrogen to be an addressable vertex.  Both return a count of atoms moved, never a bool.
"""
from ..core import (H_IMPLICIT_MAX, INFO, LOST, REFUSED, REPAIRED, LogRecord, MoleculeContainer,
                    recording)


__all__ = ['explicify_hydrogens', 'implicify_hydrogens']


#: Table-qualified rule ids, as every log record's `rule` must be -- never a bare index.
_RULE = 'hydrogens:implicify'
_RULE_BRIDGE = 'hydrogens:implicify-bridging'
_RULE_UNKNOWN = 'hydrogens:implicify-unknown-count'
_RULE_FULL = 'hydrogens:implicify-count-full'
_RULE_EXPLICIT = 'hydrogens:explicify'
_RULE_EXPLICIT_UNKNOWN = 'hydrogens:explicify-unknown-count'


def _plan(molecule: MoleculeContainer, log):
    """Which hydrogens to fold, the counts that result, and the parities to restore.

    Pure reads, and it must stay that way: the container refuses a read once an edit session is open.
    Refused kinds are an isotope, a charge, a radical, a bridging hydride, H2, a non-single bond, a
    neighbour whose count is unknown and a neighbour already holding the largest count the field can
    take -- each carries something an implicit count cannot record.
    """
    doomed: list[int] = []
    counts: dict[int, int] = {}

    for n in molecule.atom_numbers:
        if molecule.element_of(n) != 1:
            continue
        if molecule.isotope_of(n) or molecule.charge_of(n) or molecule.radical_of(n):
            continue  # a label, a hydride, a proton, a radical -- not a count
        neighbors = tuple(molecule.neighbors_of(n))
        if len(neighbors) != 1:
            if len(neighbors) > 1:
                # a bridging hydride: a record, not an exception, and the molecule is left untouched.
                log.append(LogRecord(_RULE_BRIDGE, (n,),
                                     f'hydrogen {n} bridges {len(neighbors)} atoms; it is not any '
                                     f'one atom\'s hydrogen count and was left as an atom',
                                     REFUSED))
            continue
        other = neighbors[0]
        if molecule.element_of(other) == 1:
            continue  # H2: neither atom has a heavy neighbour to fold into
        if molecule.order_of(n, other) != 1:
            log.append(LogRecord(_RULE, (n, other),
                                 f'hydrogen {n} is bonded to {other} by order '
                                 f'{molecule.order_of(n, other)}; a hydrogen count cannot record '
                                 f'that, so it was left as an atom', REFUSED))
            continue
        if molecule.implicit_h_of(other) is None:
            log.append(LogRecord(_RULE_UNKNOWN, (n, other),
                                 f'atom {other} has an unknown implicit hydrogen count, so '
                                 f'folding hydrogen {n} into it would report a total that is not '
                                 f'known; left as an atom', LOST))
            continue
        current = counts.get(other, molecule.implicit_h_of(other))
        if current >= H_IMPLICIT_MAX:
            # A record that drew more explicit hydrogens on one atom than the count field holds.  The
            # decision is per hydrogen and taken HERE, before the session opens: `set_hydrogens` would
            # raise from inside `_implicify_apply`, which is a half-applied edit and not an answer.
            log.append(LogRecord(_RULE_FULL, (n, other),
                                 f'atom {other} would reach {current + 1} implicit hydrogens and the '
                                 f'count records at most {H_IMPLICIT_MAX}, so hydrogen {n} was left as '
                                 f'an atom', REFUSED))
            continue
        doomed.append(n)
        counts[other] = current + 1

    # Read every affected anchor's parity now.  `delete_atom` clears it, and it cannot be read back
    # from inside the session.
    parities: dict[int, int] = {}
    for n in counts:
        p = molecule.parity_of(n)
        if p:
            parities[n] = p
    return doomed, counts, parities


def implicify_hydrogens(molecule: MoleculeContainer) -> int:
    """Fold ordinary hydrogen atoms into their neighbours' implicit counts.  How many atoms went.

    Returns the number of hydrogen atoms removed, which is not the number of anchors touched: for
    methane the answer is 4 where the log holds one record.

    `molecule.log` gets one record per anchor folded (naming the atom and its new total) and one per
    hydrogen refused, whether or not anyone asked.
    """
    with recording(molecule, stage='implicify') as lg:
        doomed, counts, parities = _plan(molecule, lg)
        if not doomed:
            return 0
        _implicify_apply(molecule, doomed, counts, parities)
        for n, total in counts.items():
            lg.append(LogRecord(_RULE, (n,),
                                f'atom {n}: explicit hydrogen atom(s) folded into its count, now '
                                f'{total}', REPAIRED))
    return len(doomed)


def _implicify_apply(molecule, doomed, counts, parities):
    with molecule.edit():
        for n in doomed:
            molecule.delete_atom(n)
        for n, total in counts.items():
            molecule.set_hydrogens(n, total)
        for n, p in parities.items():
            # restore what `delete_atom` cleared, verbatim: without this, deleting the explicit
            # hydrogen of `F[C@]([H])(Cl)Br` racemises the centre.
            molecule.set_parity(n, p)


def explicify_hydrogens(molecule: MoleculeContainer) -> int:
    """Turn every implicit hydrogen count into hydrogen atoms.  How many atoms arrived.

    Returns the number of hydrogen atoms added.  The new atoms are unmapped -- an invented hydrogen has
    no counterpart to correspond to, so a caller needing mapped hydrogens numbers them itself -- and
    are uncharged, non-radical, non-isotopic with a count of zero; the anchor's count goes to zero.

    An atom whose implicit count is unknown gets nothing and a `LOST` record: there is no number to
    expand and inventing zero would answer a question the record never answered.  Severity is `INFO`,
    not `REPAIRED`: both spellings are true statements about the same compound.
    """
    with recording(molecule, stage='explicify') as lg:
        plan: list[tuple] = []
        for n in molecule.atom_numbers:
            h = molecule.implicit_h_of(n)
            if h is None:
                lg.append(LogRecord(_RULE_EXPLICIT_UNKNOWN, (n,),
                                    f'atom {n} has an unknown implicit hydrogen count, so no hydrogen '
                                    f'atoms could be made explicit for it', LOST))
            elif h:
                plan.append((n, h))
        if not plan:
            return 0

        with molecule.edit():
            for n, h in plan:
                for _ in range(h):
                    # `implicit_h=0` and not the default: the default is H_UNKNOWN, which would leave
                    # the molecule's count unanswerable right after a pass that stated it.
                    molecule.add_bond(n, molecule.add_atom(1, implicit_h=0), 1)
                molecule.set_hydrogens(n, 0)
            # no parity restore, deliberately: adding an atom does not clear one.

        for n, h in plan:
            lg.append(LogRecord(_RULE_EXPLICIT, (n,),
                                f'atom {n}: {h} implicit hydrogen(s) written out as atoms', INFO))
    return sum(h for _, h in plan)
