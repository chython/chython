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
"""Salts: cutting the ionic bond (`split_salts`) and reading the record as compound plus what was drawn
beside it (`decompose_salts`).

Both read `tables/salts.tsv` and act on all 93 metals the core's `[M]` accepts.  `split_salts` is
all-or-nothing per cation atom and logs `REPAIRED`; `decompose_salts` changes nothing, logs nothing and
returns its answer.  A dative bond (order 8) is the exemption signal, never the trigger: `standardize()`
installs it to record coordination that must be preserved.
"""
from collections.abc import Iterable
from typing import NamedTuple
from ._hydrogens import implicify_hydrogens
from ._protomers import neutralize as _neutralize
from ._tables import SALT_ROLES, SaltRow, salts_rows_by_role, salts_species_keys
from ..core import REFUSED, REPAIRED, LogRecord, MoleculeContainer, recording


__all__ = ['SaltComposition', 'decompose_salts', 'split_salts']


#: The arena's charge domain, declared in `core/_atom_arena.pxi` and mirrored here as a refusal bound
#: rather than left to raise: an acceptor that would go to -5 refuses its cation and says so.
_CHARGE_MIN = -4

#: A `keep` item is a row id if it contains this, an element symbol otherwise.  Row ids are
#: table-qualified everywhere in this package, so the test is not a heuristic.
_ID_MARK = ':'

_SYMBOLS: dict[str, int] = {}


def _atomic_number(symbol: str, label: str) -> int:
    """`'Na'` -> 11, cached.

    Through a throwaway container because `chython.core` exposes no symbol table.  If core ever grows
    a public `atomic_number()`, this goes.
    """
    if symbol not in _SYMBOLS:
        probe = MoleculeContainer()
        try:
            _SYMBOLS[symbol] = probe.element_of(probe.add_atom(symbol))
        except ValueError:
            raise ValueError(f'{label} item {symbol!r} is none of: a row id (those carry a '
                             f'{_ID_MARK!r}, as in \'salts:water\'), a role '
                             f'({", ".join(SALT_ROLES)}), or an element symbol') from None
    return _SYMBOLS[symbol]


class _Keep:
    """`keep=` resolved into three sets, once per call.

    Row ids, roles and element symbols.  A `MoleculeContainer` is refused: this argument protects a
    CATION ATOM from being split, and a whole component is not a thing either pass here compares
    against.
    """
    __slots__ = ('ids', 'roles', 'elements')

    def __init__(self, items: Iterable):
        self.ids: set[str] = set()
        self.roles: set[str] = set()
        self.elements: set[int] = set()
        for item in items:
            if isinstance(item, MoleculeContainer):
                raise ValueError(
                    'a MoleculeContainer in `keep` names a whole COMPONENT, and this pass cuts bonds '
                    'inside one component rather than deleting components -- there is nothing for it to '
                    'compare against.  Name the element (keep=[\'Na\']), the row (keep=[\'salts:metal\']) '
                    'or the role (keep=[\'cation\']) instead')
            elif isinstance(item, str):
                if _ID_MARK in item:
                    self.ids.add(item)
                elif item in SALT_ROLES:
                    self.roles.add(item)
                else:
                    self.elements.add(_atomic_number(item, 'keep'))
            else:
                raise TypeError(f'keep items are row ids, roles or element symbols, not '
                                f'{type(item).__name__}')
        # ids are checked against the table, since a typo would otherwise be a silent no-op.
        unknown = self.ids - {row.id for rows in salts_rows_by_role().values() for row in rows}
        if unknown:
            raise ValueError(f'keep names no such row: {", ".join(sorted(unknown))}')

    def rows(self, rows: tuple[SaltRow, ...]):
        """The `cation`/`acceptor` rows still in play: the ones this keep set does not name."""
        return tuple(r for r in rows if r.id not in self.ids and r.role not in self.roles)


def _acceptor_atoms(molecule: MoleculeContainer, keep: _Keep) -> set[int]:
    """Every atom that some surviving `acceptor` row names as its `:1`.

    Built once per molecule, not per cation: the answer does not depend on which metal is asking.
    """
    out: set[int] = set()
    for row in keep.rows(salts_rows_by_role()['acceptor']):
        for mapping in row.query.get_mapping(molecule):
            out.add(mapping[row.anchor])
    return out


def _cation_atoms(molecule: MoleculeContainer, keep: _Keep) -> dict[int, SaltRow]:
    """Cation atom -> the first row that claimed it, in file order and ascending atom order."""
    out: dict[int, SaltRow] = {}
    for row in keep.rows(salts_rows_by_role()['cation']):
        for mapping in row.query.get_mapping(molecule):
            n = mapping[row.anchor]
            if molecule.element_of(n) not in keep.elements:
                out.setdefault(n, row)
    return dict(sorted(out.items()))


def split_salts(molecule: MoleculeContainer, *, keep: Iterable = ()) -> bool:
    """Cut every ionic cation-acceptor bond and move the charge onto the two ends.  Did anything cut?

    `CC(=O)O[Na]` becomes `CC(=O)[O-].[Na+]`.  The atom count does not change and no component is
    deleted; a component is never deleted by anything in this module -- `decompose_salts` reports
    instead.

    All-or-nothing per cation atom: a dative bond, a neighbour that is not an acceptor, an untabulated
    resulting charge or an implicit hydrogen on the cation refuses the whole atom and logs the reasons.
    `N[Pt](N)(Cl)Cl` therefore comes back intact rather than half-split, and so do ferrocene and the
    metal carbonyls.

    `keep=` takes row ids, roles (`'cation'`, `'acceptor'`) and element symbols; a `MoleculeContainer`
    is refused, this pass having no component to compare one against.
    """
    resolved = _Keep(keep)
    cations = _cation_atoms(molecule, resolved)
    if not cations:
        return False
    acceptors = _acceptor_atoms(molecule, resolved)

    # plan first, entirely in reads: the container refuses a read while a journal is pending, and the
    # all-or-nothing rule needs every condition on an atom known before any bond of it is touched.
    cuts: list[tuple[int, int]] = []
    charges: dict[int, int] = {}
    lines: list[LogRecord] = []
    for n, row in cations.items():
        neighbors = tuple(molecule.neighbors_of(n))
        if not neighbors:
            continue        # a lone ion: nothing to cut.  `decompose_salts` is the pass that counts it.
        reasons: list[str] = []
        hydrogens = molecule.implicit_h_of(n)
        if hydrogens is None:
            reasons.append('its implicit hydrogen count is unknown')
        elif hydrogens:
            reasons.append(f'it carries {hydrogens} implicit hydrogen(s)')
        for m in sorted(neighbors):
            order = molecule.order_of(n, m)
            if order == 8:
                reasons.append(f'the bond to {m} is dative (order 8), which standardize() installs to '
                               f'record coordination that must be preserved')
            elif m not in acceptors:
                reasons.append(f'atom {m} matches no acceptor row, so it cannot take the charge')
            elif molecule.charge_of(m) - 1 < _CHARGE_MIN:
                reasons.append(f'atom {m} is already at charge {molecule.charge_of(m)} and cannot go '
                               f'below {_CHARGE_MIN}')
        new_charge = molecule.charge_of(n) + len(neighbors)
        if new_charge not in row.charges:
            reasons.append(f'cutting {len(neighbors)} bond(s) would leave it at charge {new_charge}, '
                           f'which {row.id} does not tabulate '
                           f'({", ".join(str(c) for c in sorted(row.charges))})')
        if reasons:
            lines.append(LogRecord(row.id, (n,),
                                   f'atom {n} was not split: ' + '; '.join(reasons), REFUSED))
            continue
        for m in sorted(neighbors):
            cuts.append((n, m))
            charges[m] = molecule.charge_of(m) - 1
        charges[n] = new_charge

    if not cuts:
        with recording(molecule, stage='split-salts') as log:
            log.extend(lines)
        return False

    with molecule.edit():
        for n, m in cuts:
            molecule.delete_bond(n, m)
        for n, charge in charges.items():
            molecule.set_charge(n, charge)

    with recording(molecule, stage='split-salts') as log:
        log.extend(lines)
        for n, m in cuts:
            log.append(LogRecord(cations[n].id, (n, m),
                                 f'the bond between {n} and {m} was ionic, not covalent; cut, leaving '
                                 f'{n} at charge {charges[n]} and {m} at {charges[m]}', REPAIRED))
    return True


class SaltComposition(NamedTuple):
    """The inventory `decompose_salts()` returns.

    | field         | keyed by                | is                                               |
    | ---           | ---                     | ---                                              |
    | `parents`     | --                      | the compound itself, in its neutral drawing      |
    | `counterions` | row id, `'salts:tfa'`   | equivalents of each acid or base beside it       |
    | `solvates`    | row id, `'salts:water'` | equivalents of each solvent of crystallisation   |
    | `cations`     | element symbol, `'Na'`  | lone cation atoms, counted per element           |

    A `base` row counts into `counterions`: the two roles differ in which side of the salt a species
    came from, not in being what the compound was drawn beside.

    `cations` is separate from `counterions` because all 93 metals share one row id: keyed by `row.id`
    a sodium and a potassium salt would be indistinguishable.  The count is of ions, not of charge
    equivalents -- one `Ca` is `{'Ca': 1}`.
    """
    parents: tuple[MoleculeContainer, ...]
    counterions: dict[str, int]
    solvates: dict[str, int]
    cations: dict[str, int]


def decompose_salts(molecule: MoleculeContainer) -> SaltComposition:
    """Read `molecule` as a compound plus what was drawn beside it.  Changes nothing, logs nothing.

        smiles('NCC(=O)O.OC(=O)C(F)(F)F.O').decompose_salts()
        # SaltComposition(parents=(smiles('C(CN)(=O)O'),), counterions={'salts:tfa': 1},
        #                 solvates={'salts:water': 1}, cations={})

    ONE DRAWING PER COMPOUND, which is what makes the counts comparable across a corpus: the work runs
    on a copy that is hydrogen-implicified, salt-split, neutralized and aromatized, so `CC(=O)O[Na]`,
    `CC(=O)[O-].[Na+]` and `CC(=O)O.[Na+]` all report one `Na`.  `parents` holds that form and not the
    caller's drawing, and neutralizing first is why no conjugate base needs a row of its own.

    A TABULATED SPECIES IS ONLY A COUNTERION WHEN SOMETHING ELSE IS THERE TO BE THE COMPOUND.  Every
    component is a solvate row, a species row (`counterion` or `base`) or unmatched, and the three
    decide together: unmatched components are the parents when there are any; failing that the species
    rows are, with the solvates counted; failing that the solvates are.  So acetic acid answers itself
    rather than an empty `parents` and one equivalent of `salts:acetic`, and sodium chloride answers
    hydrochloric acid and one `Na`.

    A LONE CATION IS NEVER A PARENT, counted by element symbol ahead of the three -- keyed by row id all
    93 metals would be one bucket.  The count is of ions and not of charge equivalents: one `Ca` is
    `{'Ca': 1}`.

    Parents dedup by canonical bytes, so two drawn equivalents of one compound are one parent and two
    enantiomers are two.
    """
    probe = molecule.copy()
    # an explicit hydrogen is a key difference, so a solvate drawn with one would match no row; a
    # covalently drawn metal is one component until the ionic bond is cut; `keep_charge=False`
    # deliberately overrides what `canonicalize()` preserves, the net charge being part of the compound
    # but not part of an inventory of what it is beside; and `thiele()` because the writer writes what
    # is stored, so a Kekule toluene is otherwise not the tabulated one.
    implicify_hydrogens(probe)
    split_salts(probe)
    _neutralize(probe, keep_charge=False)
    probe.thiele()

    keys = salts_species_keys()
    cations: dict[str, int] = {}
    solvates: list[tuple[MoleculeContainer, SaltRow]] = []
    species: list[tuple[MoleculeContainer, SaltRow]] = []
    unmatched: list[MoleculeContainer] = []
    lone = _cation_atoms(probe, _Keep(()))       # `keep` is empty: nothing here deletes anything
    for component in probe.split():
        atoms = tuple(component.atom_numbers)
        if len(atoms) == 1 and atoms[0] in lone:
            symbol = component.atom(atoms[0]).atomic_symbol
            cations[symbol] = cations.get(symbol, 0) + 1
            continue
        row = keys.get(format(component, '!s'))
        if row is None:
            unmatched.append(component)
        elif row.role == 'solvate':
            solvates.append((component, row))
        else:
            species.append((component, row))

    if unmatched:
        parents, counted = unmatched, species + solvates
    elif species:
        parents, counted = [c for c, _ in species], solvates
    else:
        parents, counted = [c for c, _ in solvates], []

    counterion_counts: dict[str, int] = {}
    solvate_counts: dict[str, int] = {}
    for _, row in counted:
        bucket = solvate_counts if row.role == 'solvate' else counterion_counts
        bucket[row.id] = bucket.get(row.id, 0) + 1

    seen: set[bytes] = set()
    unique: list[MoleculeContainer] = []
    for component in parents:
        if component.canonical_bytes not in seen:
            seen.add(component.canonical_bytes)
            unique.append(component)
    return SaltComposition(tuple(unique), counterion_counts, solvate_counts, cations)
