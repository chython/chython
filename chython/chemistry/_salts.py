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
from collections.abc import Iterable, MutableSequence
from typing import NamedTuple
from ._hydrogens import implicify_hydrogens
from ._protomers import neutralize as _neutralize
from ._tables import SALT_CLASSES, SaltRow, salts_rows_by_klass, salts_species_keys
from ..core import REFUSED, REPAIRED, LogRecord, MoleculeContainer, recording


__all__ = ['ComponentRow', 'DEFAULT_STABILIZER_CLASSES', 'SaltComposition', 'decompose_salts', 'split_salts']


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
                             f'{_ID_MARK!r}, as in \'salts:water\'), a class '
                             f'({", ".join(SALT_CLASSES)}), or an element symbol') from None
    return _SYMBOLS[symbol]


class _Keep:
    """`keep=` resolved into three sets, once per call.

    Row ids, classes and element symbols.  A `MoleculeContainer` is refused: this argument protects a
    CATION ATOM from being split, and a whole component is not a thing either pass here compares
    against.
    """
    __slots__ = ('ids', 'classes', 'elements')

    def __init__(self, items: Iterable):
        self.ids: set[str] = set()
        self.classes: set[str] = set()
        self.elements: set[int] = set()
        for item in items:
            if isinstance(item, MoleculeContainer):
                raise ValueError(
                    'a MoleculeContainer in `keep` names a whole COMPONENT, and this pass cuts bonds '
                    'inside one component rather than deleting components -- there is nothing for it to '
                    'compare against.  Name the element (keep=[\'Na\']), the row (keep=[\'salts:metal\']) '
                    'or the class (keep=[\'metal_cation\']) instead')
            elif isinstance(item, str):
                if _ID_MARK in item:
                    self.ids.add(item)
                elif item in SALT_CLASSES:
                    self.classes.add(item)
                else:
                    self.elements.add(_atomic_number(item, 'keep'))
            else:
                raise TypeError(f'keep items are row ids, classes or element symbols, not '
                                f'{type(item).__name__}')
        # ids are checked against the table, since a typo would otherwise be a silent no-op.
        unknown = self.ids - {row.id for rows in salts_rows_by_klass().values() for row in rows}
        if unknown:
            raise ValueError(f'keep names no such row: {", ".join(sorted(unknown))}')

    def rows(self, rows: tuple[SaltRow, ...]):
        """The `embed` rows still in play: the ones this keep set does not name."""
        return tuple(r for r in rows if r.id not in self.ids and r.klass not in self.classes)


def _acceptor_atoms(molecule: MoleculeContainer, keep: _Keep) -> set[int]:
    """Every atom that some surviving `charge_acceptor` row names as its `:1`.

    Built once per molecule, not per cation: the answer does not depend on which metal is asking.
    """
    out: set[int] = set()
    for row in keep.rows(salts_rows_by_klass()['charge_acceptor']):
        for mapping in row.query.get_mapping(molecule):
            out.add(mapping[row.anchor])
    return out


def _cation_atoms(molecule: MoleculeContainer, keep: _Keep) -> dict[int, SaltRow]:
    """Cation atom -> the first row that claimed it, in file order and ascending atom order."""
    out: dict[int, SaltRow] = {}
    for row in keep.rows(salts_rows_by_klass()['metal_cation']):
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

    `keep=` takes row ids, classes (`'metal_cation'`, `'charge_acceptor'`) and element symbols; a
    `MoleculeContainer` is refused, this pass having no component to compare one against.
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


#: Narrow by design.  Widening is a deliberate act by a caller who has looked; a wide default silently
#: renames records.  Water, the mineral acids, the sulfonic acids and the C1-C2 carboxylic acids are the
#: species a hydrate or an `-HCl` is made of, and nothing else.  The organic solvents, the long and
#: polyfunctional carboxylic acids, the amine bases and the amino acids are all tabulated, all classified,
#: and all parents.
DEFAULT_STABILIZER_CLASSES = ('water', 'mineral_acid', 'sulfonic_acid', 'short_carboxylic_acid')

#: Which classes are a FORMER -- a species added to stabilize a compound.  `quaternary_ammonium` and
#: `metal_cation` are not: they are the cation half of a pair, not something added.
_ACID_CLASSES = frozenset({'mineral_acid', 'sulfonic_acid', 'short_carboxylic_acid', 'carboxylic_acid',
                           'aromatic_acid', 'fatty_acid'})
_BASE_CLASSES = frozenset({'amine_base', 'amino_acid'})
_SOLVENT_CLASSES = frozenset({'water', 'alcohol', 'hydrocarbon', 'halo_solvent', 'aprotic_solvent'})
_FORMER_CLASSES = _ACID_CLASSES | _BASE_CLASSES | _SOLVENT_CLASSES


class ComponentRow(NamedTuple):
    """One component of the normalized probe, featurized.

    `molecule` is the NORMALIZED form and `atoms` points back at the caller's drawing.  That pairing is
    the answer to "which one is it": the normalized copy costs chemistry -- charges, bond orders, explicit
    hydrogens -- and not identity, so a caller needing its own drawing has the atom numbers to rebuild it
    and a caller needing comparability has the form that compares.

    ONE DRAWN COMPONENT CAN YIELD TWO ROWS.  `CC(=O)O[Na]` is one component as drawn and two after the
    probe's `split_salts()`; the two rows' `atoms` partition the one drawn component.

    `is_lone_metal` is CHARGE-BLIND, and that is load-bearing: a neutral `[Na]` has residual charge 0, so
    keyed on charge it would read as droppable and sodium metal would be stripped off an acid.

    `equivalents` counts the rows that share this component's CANONICAL KEY **and** its ROLE.  A chloride
    on counter-ion duty and a free HCl copy share the canonical key but play different roles, so each
    counts only its own copies; `C[N+](C)(C)C.[Cl-].Cl` gives the on-duty chloride equivalents=1 and the
    stabilizer HCl equivalents=1.
    """
    atoms: tuple[int, ...]
    molecule: MoleculeContainer
    species: str | None
    klass: str | None
    equivalents: int
    heavy_atoms: int
    carbon_count: int
    ring_count: int
    charge: int
    residual_charge: int
    is_lone_metal: bool
    is_organometallic: bool
    role: str


class SaltComposition(NamedTuple):
    """The inventory `decompose_salts()` returns.

    | member            | is                                                                      |
    | ---               | ---                                                                     |
    | `components`      | every row, in canonical order                                           |
    | `parents`         | the rows whose `role` is `'parent'`                                     |
    | `stabilizers`     | the rows whose `role` is `'stabilizer'`                                 |
    | `tags`            | what the record IS -- never what this pass did with it                  |
    | `charge`          | the record's drawn charge                                               |
    | `residual_charge` | the record's intrinsic charge, after `neutralize(keep_charge=False)`     |
    """
    components: tuple[ComponentRow, ...]
    parents: tuple[ComponentRow, ...]
    stabilizers: tuple[ComponentRow, ...]
    tags: frozenset[str]
    charge: int
    residual_charge: int

    def equivalents_by_species(self) -> dict[str, int]:
        """`{species: equivalents}` over the stabilizer rows, untabulated rows omitted.

        `equivalents` is already keyed by canonical key AND role, so this is a regroup of the stabilizer
        rows: a species that appears twice as a stabilizer reports 2, and a copy on counter-ion duty is
        not a stabilizer and is not counted here.
        """
        return {row.species: row.equivalents for row in self.stabilizers if row.species is not None}

    def compose(self, rows: Iterable[ComponentRow]) -> MoleculeContainer:
        """`rows` as one molecule -- the sugar for putting a chosen selection back together.

        Over `union()` and never `substructure()`: cutting bonds destroys the frame a parity is stated in,
        so `substructure()` drops parities, wedges, stereo groups and CIP descriptors, while `split()` and
        `union()` preserve them.  The rows came from `split()`, so `remap=False` reassembles them at their
        own atom numbers and the result is the selection drawn directly.
        """
        rows = tuple(rows)
        if not rows:
            raise ValueError('compose() of no rows has no molecule; a selection of one component is the '
                             'smallest thing there is to compose')
        mine = {id(row.molecule) for row in self.components}
        if any(id(row.molecule) not in mine for row in rows):
            raise ValueError('compose() takes rows of this composition; a row from another record carries '
                             'atom numbers this one never issued, and `remap=False` would collide')
        out = rows[0].molecule.copy()
        for row in rows[1:]:
            out = out.union(row.molecule, remap=False)
        return out


def _organometallic(component: MoleculeContainer) -> bool:
    """A carbon bonded to a metal, or to boron.

    Boron is named with its reason: it carries the skeleton being transformed, like a Grignard, which is
    why a Molander trifluoroborate is a reactant and not a counterion.  A chemistry judgement, so it
    lives here and not on the atom view.
    """
    for atom in component.atoms():
        if atom.element == 6:
            for m in component.neighbors_of(atom.n):
                other = component.atom(m)
                if other.element == 5 or other.is_metal:
                    return True
    return False


def _counter_ion_duty(rows: list[dict]) -> None:
    """Mark the rows that are counter-ion to an intrinsic charge.  Sets `on_duty` in place.

    A component whose residual charge is 0 but whose DRAWN charge is nonzero and opposite in sign to the
    record's intrinsic charge stays, up to `|intrinsic|` equivalents of it.  Residual charge alone would
    let a chloride leave a quaternary ammonium; the record's charge balance alone would take the water
    off a sodium sulfonate monohydrate together with the sulfonate.  This clause does both.
    """
    intrinsic = sum(row['residual_charge'] for row in rows)
    if not intrinsic:
        return
    want = -1 if intrinsic > 0 else 1
    budget = abs(intrinsic)
    for row in rows:                                      # canonical order, so the choice is the same
        if budget <= 0:
            break
        charge = row['charge']
        if not row['residual_charge'] and charge and (charge > 0) == (want > 0):
            row['on_duty'] = True
            budget -= abs(charge)


def _undrawn_metal_duty(rows: list[dict]) -> None:
    """Mark every anion and every neutral acid on duty when a lone metal lacks its anion.  In place.

    TWO SPELLINGS OF ONE DRAWING ERROR, and the anion is absent from the drawing rather than from the
    compound either way: a metal drawn neutral states no charge at all, and a metal in a record whose drawn
    charges do not balance -- in either direction -- carries a charge the drawing does not account for.
    `CC(=O)O.[Na].CCBr`, `CC(=O)O.[Na+].CCBr` and `CC(=O)[O-].CC(=O)[O-].[Mg+].CCBr` all keep their acetate.

    ALL-OR-NOTHING, WITH NO BUDGET AND NO ORDERING: the missing charge is exactly the quantity the record
    does not state, so there is nothing to count equivalents against.  And `decompose_salts` reports rather
    than repairs, so it may not borrow `fix_salt_charges()`'s determination of that charge -- the answer for
    a magnesium drawn `[Mg+]` must not depend on whether the caller standardized first.
    """
    unpaid = sum(row['charge'] for row in rows) != 0
    if not any(row['is_lone_metal'] and (unpaid or not row['charge']) for row in rows):
        return
    for row in rows:
        if row['is_anion'] or (row['acid_sites'] and not row['charge']):
            row['on_duty'] = True


def _acid_sites(component: MoleculeContainer) -> tuple[SaltRow, ...]:
    """The `protic_acid` rows matching `component`, most acidic first.

    Rows in file order, which the loader pins to non-decreasing `order`, so the first hit is the most
    acidic row that fires.  Ties -- rows of equal `order` -- are left tied for the caller to break.
    """
    return tuple(row for row in salts_rows_by_klass()['protic_acid']
                 if next(row.query.get_mapping(component), None) is not None)


def _tags(build: list[dict], roles: list[str], guarded: bool) -> frozenset[str]:
    """What the record IS.  Never what the split did with it -- see `SaltComposition`.

    An ANION here is a component drawn negative whose residual charge is 0: a conjugate base, so a
    trifluoroborate is not one and its potassium is not a metal salt.  UN-IONIZED is drawn charge 0, the
    normalized probe having neutralized every acid.

    `acid_sites` and `is_anion` are read off the rows, the caller having computed them before the duty
    clauses that need them.
    """
    tags = set()
    if len(build) == 1:
        tags.add('single')
    if guarded:
        tags.add('stabilizer_only')
    if sum(row['charge'] for row in build) != 0:
        tags.add('charge_unbalanced')

    formers = {row['species'] for row in build if row['klass'] in _FORMER_CLASSES}
    if len(formers) > 1:
        tags.add('competing_formers')
    if len(build) > 1:
        if any(row['klass'] == 'water' for row in build):
            tags.add('hydrate')
        if any(row['klass'] in _SOLVENT_CLASSES - {'water'} for row in build):
            tags.add('solvate')

    if any(row['residual_charge'] for row, role in zip(build, roles) if role == 'parent'):
        tags.add('ion_pair')

    for row, role in zip(build, roles):
        others = [other for other in build if other is not row]
        partnered = any(other['is_anion'] or (other['acid_sites'] and not other['charge'])
                        for other in others)
        if row['klass'] in _ACID_CLASSES or row['acid_sites']:
            if any(not other['is_lone_metal'] for other in others):
                tags.add('acid_salt')
        if row['klass'] in _BASE_CLASSES and any(other['acid_sites'] for other in others):
            tags.add('base_salt')
        if row['is_lone_metal']:
            if role == 'parent':
                if partnered:
                    tags.add('metal_salt')
                elif not row['charge']:
                    tags.add('elemental_metal')
            if not row['charge'] and partnered:
                tags.add('charges_undrawn')

    parents = [row for row, role in zip(build, roles) if role == 'parent']
    if len(parents) > 1 and not any(row['klass'] in _FORMER_CLASSES or row['is_lone_metal']
                                    or row['on_duty'] for row in parents):
        tags.add('mixture')
    return frozenset(tags)


def decompose_salts(molecule: MoleculeContainer, *, classes=DEFAULT_STABILIZER_CLASSES,
                    max_atoms: int | None = None, discardable=()) -> SaltComposition:
    """Read `molecule` as its components, each featurized.  Changes nothing, logs nothing.

        smiles('NCC(=O)O.OC(=O)C(F)(F)F.O').decompose_salts().equivalents_by_species()
        # {'salts:tfa': 1, 'salts:water': 1}

    ONE DRAWING PER COMPONENT, which is what makes the answer comparable across a corpus: the work runs
    on a copy that is hydrogen-implicified, salt-split, neutralized and aromatized, so the charge-balanced
    spellings of sodium acetate (`CC(=O)[O-].[Na+]`, `CC(=O)O[Na]`) report the same rows.

    A STABILIZER IS ELIGIBLE, NEVER LEFTOVER.  A row leaves only when its residual charge is 0, it is
    neither a lone metal nor organometallic, nothing puts it on counter-ion duty (`_counter_ion_duty`,
    `_undrawn_metal_duty`), its class is named by `classes` or the row by `discardable`, and it fits
    `max_atoms`.  Subject to one guard: the parents must contain a row that is neither a lone metal nor a
    recognized solvent (`_SOLVENT_CLASSES`), and when they do not, every component is a parent -- so water
    alone answers water rather than nothing.

    | argument      | is                                                                          |
    | ---           | ---                                                                         |
    | `classes`     | which `klass` values may be a stabilizer                                    |
    | `max_atoms`   | heavy-atom ceiling on a stabilizer; `None` is no ceiling                     |
    | `discardable` | species keys, row ids or class names that are also eligible                  |
    """
    probe = molecule.copy()
    # an explicit hydrogen is a key difference, so a solvate drawn with one would match no row; a
    # covalently drawn metal is one component until the ionic bond is cut; and `thiele()` because the
    # writer writes what is stored, so a Kekule toluene is otherwise not the tabulated one.
    implicify_hydrogens(probe)
    split_salts(probe)
    probe.thiele()
    # the drawn charge is read HERE, before neutralizing, so neither number costs a second pass
    drawn = {n: probe.atom(n).charge for n in probe}
    # `keep_charge=False` deliberately overrides what `canonicalize()` preserves: the net charge is part
    # of the compound but not part of an inventory of what it is beside, and neutralizing first is why no
    # conjugate base needs a row of its own.
    _neutralize(probe, keep_charge=False)
    probe.thiele()

    keys = salts_species_keys()
    lone = _cation_atoms(probe, _Keep(()))       # `keep` is empty: nothing here deletes anything
    wanted = frozenset(classes)
    named = frozenset(discardable)
    build: list[dict] = []
    for component in probe.split():
        atoms = tuple(sorted(component.atom_numbers))
        key = format(component, '!s')
        row = keys.get(key)
        is_lone_metal = len(atoms) == 1 and component.atom(atoms[0]).is_metal
        if is_lone_metal and atoms[0] in lone:
            species, klass = lone[atoms[0]].id, lone[atoms[0]].klass
        elif row is not None:
            species, klass = row.id, row.klass
        else:
            species, klass = None, None
        build.append({'atoms': atoms, 'molecule': component, 'key': key, 'species': species,
                      'klass': klass, 'charge': sum(drawn[n] for n in atoms),
                      'residual_charge': sum(atom.charge for atom in component.atoms()),
                      'is_lone_metal': is_lone_metal, 'is_organometallic': _organometallic(component),
                      'on_duty': False})
    build.sort(key=lambda row: (row['molecule'].canonical_bytes, row['atoms']))
    for row in build:
        row['acid_sites'] = _acid_sites(row['molecule'])
        row['is_anion'] = row['charge'] < 0 and not row['residual_charge']
    _counter_ion_duty(build)
    _undrawn_metal_duty(build)

    def eligible(row: dict) -> bool:
        if row['residual_charge'] or row['is_lone_metal'] or row['is_organometallic'] or row['on_duty']:
            return False
        if row['klass'] not in wanted and not ({row['species'], row['klass'], row['key']} & named):
            return False
        return max_atoms is None or row['molecule'].atoms_count <= max_atoms

    roles = ['stabilizer' if eligible(row) else 'parent' for row in build]
    # the guard: without a parent that is neither a lone metal nor a recognized solvent, nothing here is
    # the compound, so every component is.  On-duty status is NOT excluded: a propane sulfonate on duty
    # for its Na+ is still the compound; a THF parent (aprotic_solvent) is only a solvent.
    guarded = False
    if not any(role == 'parent' and not row['is_lone_metal'] and row['klass'] not in _SOLVENT_CLASSES
               for row, role in zip(build, roles)):
        roles = ['parent'] * len(build)
        guarded = True

    # counts keyed by (canonical key, role): a component on duty and a free copy of the same structure
    # play different parts, so the on-duty row reports 1 and the stabilizer row reports its own count.
    counts: dict[tuple[str, str], int] = {}
    for row, role in zip(build, roles):
        k = (row['key'], role)
        counts[k] = counts.get(k, 0) + 1

    components = tuple(ComponentRow(row['atoms'], row['molecule'], row['species'], row['klass'],
                                    counts[(row['key'], role)], row['molecule'].atoms_count,
                                    row['molecule'].carbon_count, row['molecule'].rings_count,
                                    row['charge'], row['residual_charge'], row['is_lone_metal'],
                                    row['is_organometallic'], role)
                       for row, role in zip(build, roles))
    return SaltComposition(components,
                           tuple(row for row in components if row.role == 'parent'),
                           tuple(row for row in components if row.role == 'stabilizer'),
                           _tags(build, roles, guarded),
                           sum(row.charge for row in components),
                           sum(row.residual_charge for row in components))


#: Table-qualified, as every rule id must be.  Two: charging a metal the drawing left neutral and moving a
#: proton onto it are different claims about the drawing, and a consumer filtering the log can decline one.
_RULE_METAL = 'salts:metal-charge'
_RULE_PROTON = 'salts:charge-transfer'

#: The `<= 3` cut of §12.2's step 1.  An s- or p-block metal takes its group number as its ionic charge;
#: above that the electron count implies no charge at all -- Ti 4, Fe 8, Zn 12, Ag 11 are counts and not
#: charges -- and 0 is the f block's unknown, which refuses for the same reason.
_DETERMINATE_VALENCE = 3


def fix_salt_charges(molecule: MoleculeContainer, log: MutableSequence) -> set[int]:
    """Move the charges a salt was drawn without: `CC(=O)O.[Na]` is `CC(=O)[O-].[Na+]`.  Written ids.

    A `standardize()` stage and not a public pass, for the reasons `unite_organometallics()` gives plus
    one: run before that stage, this one would count a free halide as an anion equivalent that is about to
    stop being free.  `CC[Zn].[Cl-].CC(=O)O.[Na]` is the record that forces the order.

    THE STAGE WRITES CHARGES AND NOTHING ELSE -- no hydrogen count, no bond, no atom.  `standardize()`
    recomputes the implicit hydrogen count of every id returned here, through `calc_implicit`, which is the
    one derivation.  Every such id is a free metal or a site this stage has just charged, so none of them is
    the one atom class the ring decides: `calc_implicit` on an uncharged aromatic pnictogen answers
    `H_UNKNOWN`, while the `[n-]` a deprotonated tetrazole leaves has its class fixed by the charge and
    answers 0.

    ALL-OR-NOTHING PER RECORD, not per site, and both steps are planned entirely in reads: step 1 charges a
    metal on the strength of a balance step 2 may refuse.  `[Na]` alone is that case, and it comes back as
    `[Na]`.
    """
    metals = [atom.n for atom in molecule.atoms() if atom.is_metal and not atom.degree]
    if not metals:
        return set()

    charges = {n: molecule.charge_of(n) for n in metals}
    planned: dict[int, int] = {}
    total = sum(charges.values())
    if not any(charges.values()):
        # step 1: no metal carries a charge, so the drawing states nothing to override
        for n in metals:
            electrons = molecule.atom(n).valence_electrons
            if not 0 < electrons <= _DETERMINATE_VALENCE:
                if electrons == 0:
                    msg = (f'atom {n} is a free metal whose valence electron count is not known for this '
                           f'element; no ionic charge follows, so the whole record is left as drawn')
                else:
                    msg = (f'atom {n} is a free metal with {electrons} valence electrons, a count and not '
                           f'a charge; no ionic charge follows, so the whole record is left as drawn')
                log.append(LogRecord(_RULE_METAL, (n,), msg, REFUSED))
                return set()
            planned[n] = electrons
        total = sum(planned.values())

    # A: the anion equivalents already drawn, net of any cation that is not one of these free metals
    negative = positive = 0
    for atom in molecule.atoms():
        if atom.n in charges:
            continue
        if atom.charge < 0:
            negative -= atom.charge
        elif atom.charge > 0:
            positive += atom.charge
    anions = max(0, negative - positive)

    sites: list[tuple[int, int, int, str]] = []
    if total > anions:
        ranks = molecule.atoms_order
        best: dict[int, tuple[int, str]] = {}
        for row in salts_rows_by_klass()['protic_acid']:
            for mapping in row.query.get_mapping(molecule):
                n = mapping[row.anchor]
                if molecule.charge_of(n) or (n in best and best[n][0] <= row.order):
                    continue
                best[n] = (row.order, row.id)
        sites = sorted((order, ranks[n], n, row_id) for n, (order, row_id) in best.items())

    deficit = total - anions
    if deficit > 0 and len(sites) < deficit:
        log.append(LogRecord(_RULE_PROTON, tuple(metals),
                             f'the free metals carry {total:+d} against {anions} drawn anion '
                             f'equivalent(s) and {len(sites)} acidic site(s); the charge has nothing to '
                             f'sit on, so the whole record is left as drawn', REFUSED))
        return set()

    raises: dict[int, int] = {}
    if deficit < 0:
        # case 2: drawn anions outnumber the metals' charge, so every free metal below its group number
        # takes it -- one ionic charge per metal, the assumption step 1 makes -- and the resulting total
        # must equal the drawn anions exactly.  No drawn anions means nothing in the record states which
        # way a charge should move.
        if not anions:
            log.append(LogRecord(_RULE_METAL, tuple(metals),
                                 f'the free metals carry {total:+d} with no drawn anion equivalents; '
                                 f'the record states no direction for a charge move, so it is left '
                                 f'as drawn', REFUSED))
            return set()
        group_numbers: dict[int, int] = {}
        for n in metals:
            electrons = molecule.atom(n).valence_electrons
            if not 0 < electrons <= _DETERMINATE_VALENCE:
                if electrons == 0:
                    msg = (f'atom {n} is a free metal whose valence electron count is not known for this '
                           f'element; no ionic charge follows, so the whole record is left as drawn')
                else:
                    msg = (f'atom {n} is a free metal with {electrons} valence electrons, a count and not '
                           f'a charge; no single ionic charge follows, so the whole record is left as '
                           f'drawn')
                log.append(LogRecord(_RULE_METAL, (n,), msg, REFUSED))
                return set()
            group_numbers[n] = electrons
        held = {n: planned.get(n, charges[n]) for n in metals}
        raises = {n: group_numbers[n] for n in metals if held[n] < group_numbers[n]}
        reached = sum(raises.get(n, held[n]) for n in metals)
        if reached != anions:
            log.append(LogRecord(_RULE_METAL, tuple(metals),
                                 f'the free metals at their group numbers carry {reached:+d} against '
                                 f'{anions} drawn anion equivalent(s), so no assignment of ionic charges '
                                 f'balances the record; it is left as drawn', REFUSED))
            return set()

    taken = sites[:deficit] if deficit > 0 else []
    if not planned and not raises and not taken:
        return set()

    with molecule.edit():
        for n, charge in planned.items():
            molecule.set_charge(n, charge)
        for n, charge in raises.items():
            molecule.set_charge(n, charge)
        for _, _, n, _ in taken:
            molecule.set_charge(n, -1)   # sites are uncharged; the `best` build skips charged ones

    for n, charge in planned.items():
        log.append(LogRecord(_RULE_METAL, (n,), f'atom {n} is a free metal drawn neutral; charged '
                                                f'{charge:+d}, its group number', REPAIRED))
    for n, charge in raises.items():
        log.append(LogRecord(_RULE_METAL, (n,), f'atom {n} carries less charge than the drawn anions need; '
                                                f'raised to {charge:+d}, its group number', REPAIRED))
    for order, _, n, row_id in taken:
        log.append(LogRecord(_RULE_PROTON, (n,), f'atom {n} matches {row_id} at acidity rung {order} and '
                                                 f'is the site the metals take a proton from; charged -1',
                             REPAIRED))
    return set(planned) | set(raises) | {n for _, _, n, _ in taken}
