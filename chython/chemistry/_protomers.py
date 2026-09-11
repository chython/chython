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
"""Protomers: moving a proton between charged sites.  `neutralize()` and nothing else yet.

Reads `tables/acids.tsv` -- an `acid` row is a cation holding an implicit hydrogen, a `base` row an
anion that can take one -- and moves the proton from one to the other, leaving both ends neutral.  No
bond and no atom changes: a charge and an implicit count do.

NOTHING OVERSHOOTS ZERO.  A move is admissible when it takes every component it touches closer to
charge zero, or when it is a pair inside ONE component, which leaves that component's charge alone.
That single rule is what stops nitrate's second oxygen being protonated into `H2NO3+` while sulfate's
second still gives sulfuric acid, and it is why `keep_charge=False` is one branch rather than a
different algorithm.
"""
from ._implicit import environment_of
from ._tables import AcidRow, acids_rules_by_role
from ..core import INFO, REFUSED, LogRecord, MoleculeContainer, recording
from ..core._core import valence_check


__all__ = ['neutralize']


def _candidates(molecule: MoleculeContainer, role: str) -> dict[int, AcidRow]:
    """Site atom -> the first row that claimed it, in file order and ascending atom order.

    Every row of a role demands the same charge, so the first claim is as good as any and the id in
    the log is the most specific pattern the table had.
    """
    out: dict[int, AcidRow] = {}
    for row in acids_rules_by_role()[role]:
        for mapping in row.query.get_mapping(molecule):
            out.setdefault(mapping[row.anchor], row)
    return dict(sorted(out.items()))


def _admitted(molecule: MoleculeContainer, candidates: dict[int, AcidRow], delta: int,
              counts: dict[int, int], lines: list[LogRecord]) -> dict[int, AcidRow]:
    """The candidates whose hydrogen count is derivable and whose neutral form has a valence row.

    Fills `counts` with each admitted site's implicit hydrogen count.  The edit session below writes the
    new count and cannot read the old one -- a container read inside an open scope answers from the
    pre-scope arena and raises -- so every site is measured here, before anything is planned.
    """
    out: dict[int, AcidRow] = {}
    for n, row in candidates.items():
        hydrogens = molecule.implicit_h_of(n)
        if hydrogens is None:
            lines.append(LogRecord(row.id, (n,),
                                   f'atom {n} has no derivable implicit hydrogen count, so a proton '
                                   f'cannot be counted off or onto it', REFUSED))
            continue
        order_sum, environment, aromatic = environment_of(molecule, n)
        # an aromatic bond has no valence row, so the question cannot be put -- `check_valence` calls
        # that `unknown` rather than a violation, and a pyridinium must stay deprotonatable.
        if not aromatic and valence_check(molecule.element_of(n), 0, molecule.radical_of(n),
                                         order_sum, hydrogens + delta, environment) == 'violation':
            lines.append(LogRecord(row.id, (n,),
                                   f'atom {n} was left charged: neutral with {hydrogens + delta} '
                                   f'implicit hydrogen(s) is a valence violation', REFUSED))
            continue
        out[n] = row
        counts[n] = hydrogens
    return out


def _charges(molecule: MoleculeContainer, labels: dict[int, int]) -> dict[int, int]:
    """Component label -> its total formal charge."""
    out: dict[int, int] = {}
    for n, label in labels.items():
        out[label] = out.get(label, 0) + molecule.charge_of(n)
    return out


def neutralize(molecule: MoleculeContainer, *, keep_charge: bool = True) -> bool:
    """Move every proton the acid/base table can move from a cation onto an anion.  Anything moved?

    `[NH3+]CC(=O)[O-]` becomes `NCC(=O)O` and `C[NH3+].[Cl-]` becomes `CN.Cl`.  Charges and implicit
    hydrogen counts are the only things written; the atoms, the bonds and the components are untouched,
    which is what separates this from `split_salts`, which cuts a bond.

    `keep_charge=True` moves protons in PAIRS, so the total charge is exactly preserved and a record
    that cannot be balanced comes back partly neutral -- `[NH3+]CC[NH3+].[O-][N+](=O)[O-]` gives
    ethylenediamine's monocation beside nitric acid, since the second nitrate oxygen would take the
    nitrate past zero.  `keep_charge=False` lets a site act alone, as far as its own component's charge
    allows: `C[NH3+]` alone becomes `CN`.

    Sites are found by `tables/acids.tsv`, whose `h` primitive reads IMPLICIT hydrogens, so a molecule
    carrying explicit hydrogen atoms wants `implicify_hydrogens()` first.  A site whose hydrogen count
    is not derivable, and one whose neutral form no valence row accepts, is refused and logged.
    """
    lines: list[LogRecord] = []
    acids = _candidates(molecule, 'acid')
    bases = _candidates(molecule, 'base')
    if not acids and not bases:
        return False
    if keep_charge and not (acids and bases):
        return False    # nothing to pair with, so no decision was made and there is nothing to log

    counts: dict[int, int] = {}
    donors = _admitted(molecule, acids, -1, counts, lines)
    acceptors = _admitted(molecule, bases, 1, counts, lines)
    labels = molecule.component_labels()
    charges = _charges(molecule, labels)
    site_charges = {n: molecule.charge_of(n) for n in (*donors, *acceptors)}

    # a pair inside one component leaves its charge alone, so those go first and unconditionally.
    grouped: dict[int, tuple[list[int], list[int]]] = {}
    for n in donors:
        grouped.setdefault(labels[n], ([], []))[0].append(n)
    for n in acceptors:
        grouped.setdefault(labels[n], ([], []))[1].append(n)

    moves: list[tuple[int | None, int | None]] = []
    spare_donors: list[int] = []
    spare_acceptors: list[int] = []
    for label in sorted(grouped):
        inside, outside = grouped[label]
        paired = min(len(inside), len(outside))
        moves.extend(zip(inside[:paired], outside[:paired]))
        spare_donors.extend(inside[paired:])
        spare_acceptors.extend(outside[paired:])

    # what is left crosses a component boundary, and both ends must move toward zero.  One side of
    # every component is exhausted by now, so a spare donor and a spare acceptor are never the same
    # component and the two tests are independent.
    if keep_charge:
        available = list(spare_acceptors)
        for n in spare_donors:
            if charges[labels[n]] <= 0:
                lines.append(_stranded(n, donors[n], 'deprotonating', charges[labels[n]]))
                continue
            for m in available:
                if charges[labels[m]] < 0:
                    moves.append((n, m))
                    charges[labels[n]] -= 1
                    charges[labels[m]] += 1
                    available.remove(m)
                    break
            else:
                lines.append(LogRecord(donors[n].id, (n,),
                                       f'atom {n} stays at charge {site_charges[n]}: no anion '
                                       f'is left that could take its proton without being taken past '
                                       f'charge zero', REFUSED))
        for m in available:
            lines.append(LogRecord(acceptors[m].id, (m,),
                                   f'atom {m} stays at charge {site_charges[m]}: no cation is '
                                   f'left with a proton to give it', REFUSED))
    else:
        for n in spare_donors:
            if charges[labels[n]] > 0:
                moves.append((n, None))
                charges[labels[n]] -= 1
            else:
                lines.append(_stranded(n, donors[n], 'deprotonating', charges[labels[n]]))
        for m in spare_acceptors:
            if charges[labels[m]] < 0:
                moves.append((None, m))
                charges[labels[m]] += 1
            else:
                lines.append(_stranded(m, acceptors[m], 'protonating', charges[labels[m]]))

    if not moves:
        with recording(molecule, stage='neutralize') as log:
            log.extend(lines)
        return False

    with molecule.edit():
        for n, m in moves:
            if n is not None:
                molecule.set_charge(n, site_charges[n] - 1)
                molecule.set_hydrogens(n, counts[n] - 1)
            if m is not None:
                molecule.set_charge(m, site_charges[m] + 1)
                molecule.set_hydrogens(m, counts[m] + 1)

    with recording(molecule, stage='neutralize') as log:
        log.extend(lines)
        for n, m in moves:
            if n is None:
                log.append(LogRecord(acceptors[m].id, (m,),
                                     f'atom {m} took a proton and is neutral; keep_charge=False, so '
                                     f'no cation paid for it', INFO))
            elif m is None:
                log.append(LogRecord(donors[n].id, (n,),
                                     f'atom {n} gave its proton up and is neutral; keep_charge=False, '
                                     f'so no anion took it', INFO))
            else:
                log.append(LogRecord(donors[n].id, (n, m),
                                     f'the proton on {n} moved to {m} ({acceptors[m].id}); both ends '
                                     f'are neutral and the total charge is unchanged', INFO))
    return True


def _stranded(n: int, row: AcidRow, doing: str, charge: int) -> LogRecord:
    """The refusal for a site whose own component would be taken past charge zero."""
    return LogRecord(row.id, (n,),
                     f'atom {n} stays charged: its component is at charge {charge}, so {doing} it '
                     f'would take that component away from zero rather than toward it', REFUSED)
