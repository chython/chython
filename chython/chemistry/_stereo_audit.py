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
"""`audit_stereo`: what a drawing's stereo labels state, as findings an agent can act on.

A reader: the molecule is not changed and nothing is logged.  Codes and messages are `tables/stereo_audit.tsv`.
"""
from collections.abc import Sequence
from typing import NamedTuple
from ._tables import stereo_audit_codes
from ..core import STEREO_ABS, STEREO_AND, STEREO_OR, STEREO_UNSPECIFIED, identity_memo


__all__ = ['Fix', 'Finding', 'StereoAudit', 'audit_stereo']

_SEVERITY_ORDER = {'error': 0, 'warning': 1, 'info': 2}
# Read-time log rules passed through as findings: the `read` rows of `docs/stereo.rst`, Index of stereo log rules.
_PASS_THROUGH = frozenset((
    'smiles:stereo-unknown', 'smiles:chiral-class-unsupported', 'smiles:stereo-no-unit',
    'smiles:stereo-unit-incomplete', 'smiles:stereo-too-many-dirs', 'smiles:stereo-ambiguous-h',
    'smiles:stereo-unsupported-unit', 'smiles:stereo-allene-bad-bonds', 'smiles:stereo-allene-overfilled',
    'smiles:stereo-one-sided', 'smiles:stereo-same-side', 'smiles:stereo-unused-dirs',
    'smiles:doubled-bond-direction', 'smiles:ring-bond-dir-conflict', 'smiles:cx-wavy-bond',
    'smiles:cx-wavy-bad-index', 'smiles:cx-cis-trans-unknown', 'smiles:cx-ctu-bad-index',
    'smiles:cx-wedge-not-applied', 'smiles:cx-wedge-bad-index', 'smiles:cx-wedge-contradictory',
    'smiles:cx-relative-flag-redundant', 'smiles:cx-relative-flag-unbacked', 'smiles:stereo-group-malformed',
    'smiles:stereo-group-empty', 'smiles:stereo-group-bad-index', 'smiles:stereo-group-duplicate',
    'smiles:stereo-group-renumbered', 'smiles:stereo-group-refused', 'smiles:cip-field-malformed',
    'smiles:cip-field-bad-descriptor', 'smiles:cip-undetermined', 'smiles:cip-bad-index', 'smiles:cip-duplicate',
    'smiles:cip-refused',
    'smarts:unresolved-stereo-dropped',
    'ctab:wedge-dropped', 'ctab:configuration-out-of-range', 'ctab:stereo-out-of-range',
    'ctab:stereo-group-dropped', 'ctab:chiral-flag-ignored', 'ctab:dimensionality-mismatch',
    'v2000:bad-bond-stereo', 'v2000:cis-or-trans-on-non-double', 'v2000:stereo-care',
    'v2000:chiral-flag-contradicts-collection',
    'v3000:chiral-flag-contradicts-collection', 'v3000:bond-bad-cfg', 'v3000:collection-unknown',
    'v3000:collection-no-group', 'v3000:collection-wrong-object', 'v3000:collection-atom-ref',
    'v3000:collection-bond-ref', 'v3000:collection-group-renumbered', 'v3000:collection-group-dropped',
    'v3000:collection-namespaces-merged', 'v3000:collection-duplicate-member',
    'v3000:collection-member-unconfigured',
    'sgroup:stereo-label-promoted', 'sgroup:stereo-label-not-promoted',
    'wedge:either-bond', 'wedge:drawing-field-disagree', 'wedge:reanchored', 'wedge:3d-wedge-overridden',
    'wedge:no-coords-read', 'wedge:hash-between-centres', 'wedge:non-stereogenic', 'wedge:on-double-bond',
    'wedge:undrawn-directions', 'wedge:ambiguous-drawing', 'wedge:ambiguous-parity-frame',
    'wedge:degenerate-geometry', 'wedge:collinear', 'wedge:same-side-substituents', 'wedge:no-partner-terminal',
    'wedge:bad-letter', 'wedge:ambiguous-bare-letter', 'wedge:wrong-ref-count', 'wedge:wrong-bond',
    'wedge:wrong-substituent', 'wedge:allene-chain-ends', 'wedge:allene-terminals', 'wedge:allene-one-direction',
    'wedge:allene-both-terminals', 'wedge:allene-same-direction', 'wedge:atropisomer-pivots',
    'wedge:atropisomer-both-pivots', 'wedge:atropisomer-same-direction', 'wedge:atropisomer-unsupported-mark',
    'wedge:helix-contradictory', 'wedge:hapto-not-bonded', 'wedge:unsupported-kind', 'wedge:unsettled',
    'cml:atom-parity-bad-value', 'cml:atom-parity-wrong-ref-count', 'cml:atom-parity-unknown-atom',
    'cml:atom-parity-phantom-twice', 'cml:bond-stereo-empty', 'cml:bond-stereo-unknown-letter',
    'cml:bond-stereo-unknown-convention', 'cml:bond-stereo-convention-value-bad',
    'cml:bond-stereo-cis-or-trans-unknown', 'cml:bond-stereo-wrong-ref-count', 'cml:bond-stereo-unknown-atom',
    'mrv:bond-stereo-empty', 'mrv:bond-stereo-bad-letter', 'mrv:bond-stereo-bad-value', 'mrv:bond-stereo-bad-code',
    'mrv:bond-stereo-foreign-dict', 'mrv:bond-stereo-cis-trans-unknown',
    'mmcif:chem-comp-stereo',
))
# Map core log severity strings to the audit severity vocabulary.
_LOG_SEVERITY = {'lost': 'warning', 'refused': 'error', 'repaired': 'info', 'info': 'info'}
# String name to core kind int, for the `args` of `set-group` and `merge-groups`.
_KIND = {'ABS': STEREO_ABS, 'AND': STEREO_AND, 'OR': STEREO_OR}
_KIND_NAME = {v: k for k, v in _KIND.items()}
_DEAD = ('stereo-audit:dead-unconfigured', 'stereo-audit:dead-not-stereogenic', 'stereo-audit:dead-restating-group')
_LABEL = {STEREO_ABS: 'abs', STEREO_AND: '&', STEREO_OR: 'o'}
_UNREALISABLE = 'stereo:bridgehead-unrealisable'
_STRENGTH = {STEREO_AND: 0, STEREO_OR: 1, STEREO_ABS: 2}
_SET_LABEL = ('stereo-audit:bridgehead-redundant-label', 'stereo-audit:bridgehead-half-labelled',
              'stereo-audit:bridgehead-merged-groups')


class Fix(NamedTuple):
    """One actionable suggestion: what to change, on which atoms, and what it costs.

    `action` is one of `'drop-label'`, `'set-sign'`, `'merge-groups'`, `'set-group'`, `'clear-sign'`.
    `atoms` keeps the owner pair for an axis; `_apply_fix` resolves it to the unit's anchor before
    calling `set_parity`.  `args` carries the action's parameters: `('AND', 1)` for `set-group`;
    `(1,)` or `(2,)` for `set-sign`; `()` otherwise.  `preserves_meaning` is True only when the
    fixed copy compares `==` and `str()` equal to the original.  `result` is `str()` of the fixed copy.
    """
    action: str
    atoms: tuple
    args: tuple
    preserves_meaning: bool
    result: str


class Finding(NamedTuple):
    """One labelling observation, with its severity, context and actionable fixes.

    `code` is the qualified rule id (`'stereo-audit:<code>'` or a pass-through rule id).
    `severity` is `'error'`, `'warning'` or `'info'`.  `atoms` are the atoms the observation is
    about.  `group` is `(kind, stored_id)` from `stereo_groups()`, or None.  `because` names the
    atoms or groups that motivate the observation.  `message` is the formatted template sentence.
    `fixes` is a sequence of `Fix` suggestions (`Fixes`), equal to the tuple of them.
    """
    code: str
    severity: str
    atoms: tuple
    group: tuple | None
    because: tuple
    message: str
    fixes: tuple


class Fixes(Sequence):
    """A finding's fixes, built on first read: each costs a fixed copy's `str()` and identity, and a finding
    that only reports its fixes need not pay for them until someone looks.  Fixes with an equal `result`
    keep the first.  Equal to, hashed and pickled as the plain tuple of `Fix`."""
    __slots__ = ('_pending', '_fixes')

    def __init__(self, fixes):
        self._pending, self._fixes = tuple(fixes), None

    def _built(self):
        if self._fixes is None:
            unique = {}
            with identity_memo():                   # a finding's fixes share most of their identities
                for fix in self._pending:
                    fix = fix() if callable(fix) else fix
                    unique.setdefault(fix.result, fix)
            self._fixes, self._pending = tuple(unique.values()), ()
        return self._fixes

    def __getitem__(self, i):
        return self._built()[i]

    def __len__(self):
        return len(self._built())

    def __eq__(self, other):
        return self._built() == (other._built() if isinstance(other, Fixes) else other)

    def __hash__(self):
        return hash(self._built())

    def __repr__(self):
        return repr(self._built())

    def __reduce__(self):
        return tuple, (self._built(),)


class StereoAudit(NamedTuple):
    """The full result of one `audit_stereo` call.

    `findings` is ordered by severity (error, warning, info) then by first atom then by code.
    `clean` is True when no finding is above `'info'` severity.
    """
    findings: tuple

    @property
    def clean(self):
        """No error or warning findings."""
        return all(f.severity == 'info' for f in self.findings)


def _first(atoms):
    """Sort key: the first atom number, or 0 for a molecule-level finding."""
    a = atoms[0] if atoms else 0
    return a[0] if isinstance(a, tuple) else a


def _units_text(xs):
    """Units as plain text: `2, 5`; an axis is its owner pair, `4-6`."""
    return ', '.join(f'{x[0]}-{x[1]}' if isinstance(x, tuple) else str(x) for x in xs)


def _group_text(group):
    """A `(kind, id)` key, or a tuple of them, as CXSMILES labels: `&1`, `o2`, `abs`."""
    if not group:
        return ''
    keys = group if isinstance(group[0], tuple) else (group,)
    return ', '.join(_LABEL[k] if k == STEREO_ABS else f'{_LABEL[k]}{g}' for k, g in keys)


def _finding(code, atoms, group=None, because=(), fixes=()):
    """Build one Finding from the table row for `code`.  `fixes` holds `Fix`es or `_Context.later` thunks."""
    row = stereo_audit_codes()[code]
    message = row.template.format(atoms=_units_text(atoms), because=_units_text(because), group=_group_text(group))
    return Finding(f'stereo-audit:{code}', row.severity, tuple(atoms), group,
                   tuple(because), message, Fixes(fixes))


def _anchors(units):
    """Axis owner pair to its anchor, for `set_parity`, which takes a single atom."""
    return {u['owners']: u['anchor'] for u in units if isinstance(u['owners'], tuple)}


def _apply_fix(molecule, action, atoms, args, anchors=None):
    """A copy of `molecule` with one fix applied through the edit session.

    An owner-pair element in `atoms` is resolved to its unit's anchor through `anchors` (`_anchors`,
    computed from `molecule` when None) before calling `set_parity`.  `set_stereo_group` accepts
    pairs natively, so no resolution is needed for label actions.
    """
    owners_to_anchor = _anchors(molecule.stereo_units()) if anchors is None else anchors
    m = molecule.copy()
    with m.edit() as e:
        if action == 'set-sign':
            owner = atoms[0]
            target = owners_to_anchor[owner] if isinstance(owner, tuple) else owner
            e.set_parity(target, args[0])
        elif action == 'clear-sign':
            for a in atoms:
                target = owners_to_anchor[a] if isinstance(a, tuple) else a
                e.set_parity(target, 0)
        elif action == 'drop-label':
            for a in atoms:
                # STEREO_UNSPECIFIED (0) clears the group membership entirely -- the spelling
                # `clean_stereo_groups` uses by zeroing the group byte at each slot.
                e.set_stereo_group(a, STEREO_UNSPECIFIED)
        elif action == 'set-group':
            kind = _KIND[args[0]]
            gid = args[1] if len(args) > 1 else 0
            for a in atoms:
                e.set_stereo_group(a, kind, gid)
        elif action == 'merge-groups':
            kind = _KIND[args[0]]
            gid = args[1]
            for a in atoms:
                e.set_stereo_group(a, kind, gid)
    return m


class _Context:
    """The input, its baseline copy, and the unit table, computed once per `audit_stereo` call.

    The audit compares MEANINGS: `==` and `str()` follow the stored labels, so two molecules are compared, and
    a fix's `result` written, after `canonicalize_stereo_groups()` on a copy (`canonical`).
    """
    __slots__ = ('molecule', 'baseline', 'canonical_baseline', 'baseline_str', 'units', 'anchors', '_canon')

    def __init__(self, molecule):
        # a snapshot, since fixes deferred by `later` run after the call returns and the caller may edit
        self.molecule = molecule.copy()
        # str() and == can write to a log (`stereo:bridgehead-unrealisable`); a copy writes to its own
        self.baseline = molecule.copy()
        self._canon = {}
        self.canonical_baseline = self.canonical(self.baseline)
        self.baseline_str = str(self.canonical_baseline)
        self.units = molecule.stereo_units()
        self.anchors = _anchors(self.units)

    def canonical(self, molecule):
        """`molecule` canonicalized, on a copy, memoised for the molecules this audit made."""
        key = id(molecule)
        hit = self._canon.get(key)
        if hit is None or hit[0] is not molecule:
            out = molecule.copy()
            out.canonicalize_stereo_groups()
            hit = self._canon[key] = (molecule, out)
        return hit[1]

    def same(self, a, b):
        """Do `a` and `b` state one meaning?"""
        return self.canonical(a) == self.canonical(b)

    def apply(self, action, atoms, args=(), molecule=None):
        """`_apply_fix` on `molecule`, the input by default, resolving axes through the unit table."""
        return _apply_fix(self.molecule if molecule is None else molecule, action, atoms, args, self.anchors)

    def fix(self, action, atoms, args=()):
        """One Fix: the canonical `str()` of the fixed copy, and whether its meaning is the baseline's --
        `==`, `hash` and `str()` of the two canonicalized."""
        fixed = self.canonical(self.apply(action, atoms, args))
        fixed_str = str(fixed)
        base = self.canonical_baseline
        same = fixed == base and hash(fixed) == hash(base) and fixed_str == self.baseline_str
        return Fix(action, tuple(atoms), tuple(args), same, fixed_str)

    def later(self, action, atoms, args=()):
        """`fix` deferred until the finding's fixes are read: for a finding that only reports them."""
        return lambda: self.fix(action, atoms, args)


def _dead(ctx, code, atoms, group):
    """A dead-* finding on `atoms` of `group`, or None when dropping the label moves the meaning.

    A label whose drop moves `==`, `hash` or `str()` is respelled onto another unit by the canonical
    form, not dead, and gets no dead-* finding.
    """
    fix = ctx.fix('drop-label', atoms)
    if fix.preserves_meaning:
        return _finding(code, atoms, group=group, fixes=(fix,))


def _units(ctx):
    """Per-unit findings: unknown signs and dead labels."""
    out = []
    molecule = ctx.molecule
    units_by_anchor = {u['anchor']: u for u in ctx.units}
    units_by_owners = {u['owners']: u for u in ctx.units if isinstance(u['owners'], tuple)}

    # unsigned stereogenic units; a unit the labels already state, signing it keeps the meaning.
    # An unsigned coupled bridgehead set is one choice: one finding on its lowest member, the rest `because`.
    stereogenic = [u for u in ctx.units if u['stereogenic']]
    unsigned = {u['anchor'] for u in stereogenic if u['parity'] == 0 and not u['coupled_implied']}
    for u in stereogenic:
        if u['anchor'] not in unsigned:
            continue
        partners = tuple(sorted(c for c in u['coupled'] if c in unsigned))
        if partners and partners[0] < u['anchor']:
            continue                    # named by the lowest member of its coupled set
        atoms = (u['owners'],)
        fixes = tuple(ctx.fix('set-sign', atoms, (p,)) for p in (1, 2))
        if not all(x.preserves_meaning for x in fixes):
            out.append(_finding('unknown-unit', atoms, because=partners, fixes=fixes))

    # dead group members
    for (kind, gid), members in molecule.dead_stereo_groups().items():
        for x in members:
            if kind != STEREO_ABS and molecule.stereo_is_live(x):
                continue                # a live member of a dead AND/OR group: `_groups` reports the group
            if isinstance(x, tuple):
                u = units_by_owners.get(x)
            else:
                u = units_by_anchor.get(x)
            parity = u['parity'] if u is not None else 0
            code = 'dead-unconfigured' if parity == 0 else 'dead-not-stereogenic'
            f = _dead(ctx, code, (x,), (kind, gid))
            if f is not None:
                out.append(f)
    return out


def _member_key(x):
    """Sort key over members: an atom number, an axis after the centre on its first owner."""
    return (x, 0) if isinstance(x, int) else (x[0], 1, x[1])


def _restates(ctx, members):
    """Is `members` a configuration that is its own mirror image?  With every label read as ABS, inverting
    them together restates the molecule and inverting any one of them does not: a cis 1,3-cyclobutane,
    whatever else the molecule carries.  One member is never one -- its sign then states nothing alone."""
    if len(members) < 2:
        return False
    labelled = tuple(x for ms in ctx.molecule.stereo_groups().values() for x in ms)
    fixed = ctx.apply('set-group', labelled, ('ABS',))
    anchors = ctx.anchors

    def inverted(xs):
        out = fixed.copy()
        with out.edit() as e:
            for x in xs:
                target = anchors[x] if isinstance(x, tuple) else x
                e.set_parity(target, 3 - fixed.parity_of(target))
        return out

    return ctx.same(inverted(members), fixed) and all(not ctx.same(inverted((x,)), fixed) for x in members)


def _groups(ctx):
    """Group findings: restating, racemate, split, absorbed groups and unlabelled signs.

    `live_stereo_groups()` keys are canonical ids; every member is mapped back to its stored key.
    """
    out = []
    molecule = ctx.molecule
    stored = molecule.stereo_groups()
    stored_of = {x: key for key, members in stored.items() for x in members}
    component_of = {a: c for c, atoms in enumerate(molecule.connected_components) for a in atoms}

    def component(x):
        return component_of[x if isinstance(x, int) else x[0]]

    # a dead group whose members are each live: inverting it restates the molecule
    dead_groups = molecule.dead_stereo_groups()
    for key, members in dead_groups.items():
        if key[0] == STEREO_ABS:
            continue                    # ABS has no inversion; `_units` reports its dead members
        restating = tuple(x for x in members if molecule.stereo_is_live(x))
        if restating:
            f = _dead(ctx, 'dead-restating-group', restating, key)
            if f is not None:
                out.append(f)

    live = molecule.live_stereo_groups()
    # live AND groups per component, each as its members in that component
    ands = {}
    for (kind, _), members in live.items():
        if kind == STEREO_AND:
            split = {}
            for x in members:
                split.setdefault(component(x), []).append(x)
            for c, part in split.items():
                ands.setdefault(c, []).append((part, len(members)))
    live_units = {}
    for u in ctx.units:
        if u['parity'] and molecule.stereo_is_live(u['owners']):
            live_units.setdefault(component(u['owners']), []).append(u['owners'])
    for c, parts in ands.items():
        if len(parts) == 1:
            [(part, size)] = parts
            if size == 1 and live_units.get(c) == part:
                out.append(_finding('and-single-member-racemate', tuple(part), group=stored_of.get(part[0])))
        else:
            atoms = tuple(sorted((x for part, _ in parts for x in part), key=_member_key))
            keys = tuple(sorted({stored_of[x] for x in atoms if x in stored_of}))
            gid = min((g for _, g in keys), default=1)
            fix = ctx.later('merge-groups', atoms, ('AND', gid))
            out.append(_finding('and-split', atoms, group=keys, fixes=(fix,)))

    # a stored AND or OR member read as ABS: another group fixes its phase, or with no other group in its
    # component the configuration is its own mirror image, which a group states as ABS does
    absolute = set(live.get((STEREO_ABS, 0), ()))
    for key, members in stored.items():
        if key[0] not in (STEREO_AND, STEREO_OR):
            continue
        absorbed = tuple(x for x in members if x in absolute)
        if not absorbed:
            continue
        fix = ctx.fix('set-group', absorbed, ('ABS',))
        if not fix.preserves_meaning:
            continue                    # a respelling of the same groups, not an absorbed phase
        here = {component(x) for x in absorbed}
        because = tuple(sorted((x for k, ms in stored.items() if k != key for x in ms if component(x) in here),
                               key=_member_key))
        if not because or _restates(ctx, absorbed):
            out.append(_finding('meso-absolute', absorbed, group=key, fixes=(fix,)))
        elif key[0] == STEREO_AND:
            out.append(_finding('and-absorbed', absorbed, group=key, because=because, fixes=(fix,)))

    # a sign without a group in a component where another unit carries a live AND or OR label
    dead = {x for ms in dead_groups.values() for x in ms}
    relative = {component(x) for (kind, _), ms in stored.items() if kind in (STEREO_AND, STEREO_OR)
                for x in ms if x not in dead}
    if relative:
        free = {kind: 1 + max((g for k, g in stored if k == kind), default=0) for kind in (STEREO_AND, STEREO_OR)}
        for u in ctx.units:
            x = u['owners']
            if not u['parity'] or x in stored_of or component(x) not in relative:
                continue
            if ctx.same(ctx.apply('clear-sign', (x,)), ctx.baseline):
                continue                # the sign states nothing
            fixes = (ctx.later('set-group', (x,), ('ABS',)),
                     ctx.later('set-group', (x,), ('AND', free[STEREO_AND])),
                     ctx.later('set-group', (x,), ('OR', free[STEREO_OR])))
            out.append(_finding('unlabelled-sign', (x,), fixes=fixes))
    return out


def _ring_systems(molecule):
    """Atom number to ring system id: rings sharing an atom, fused or spiro, are one system."""
    parent = {}

    def root(a):
        while parent[a] != a:
            a = parent[a]
        return a

    for ring in molecule.sssr:
        for a in ring:
            parent.setdefault(a, a)
        r = root(ring[0])
        for a in ring[1:]:
            q = root(a)
            if q != r:
                parent[q] = r
    return {a: root(a) for a in parent}


def _unit_atoms(u):
    """The anchor and owners of a unit."""
    x = u['owners']
    return {u['anchor'], *x} if isinstance(x, tuple) else {u['anchor'], x}


def _waits_for(ctx):
    """`dead-waits-for`: a sign that states nothing until an unsigned unit of its ring system is signed.

    Unit `u` is a candidate when clearing its sign leaves the copy `==` the baseline.  It waits for an
    unsigned stereogenic unit `c` of its ring system when signing `c` moves the meaning and clearing `u`
    in that signed copy moves it again.  A `c` the drawing already states, as the spiro centre of
    `|&1:1,&2:6|` on 2,6-dimethylspiro[3.3]heptane, is not waited for.
    """
    out = []
    system = _ring_systems(ctx.molecule)
    if not system:
        return out
    stored_of = {x: key for key, members in ctx.molecule.stereo_groups().items() for x in members}

    def systems(u):
        return {system[a] for a in _unit_atoms(u) if a in system}

    unsigned = [(c['owners'], systems(c)) for c in ctx.units
                if c['stereogenic'] and not c['parity'] and not c['coupled_implied']]
    unsigned = [(y, s) for y, s in unsigned if s]
    signed_copy = {}
    sign_fixes = {}
    for u in ctx.units:
        if not u['parity'] or not unsigned:
            continue
        here = systems(u)
        targets = [y for y, s in unsigned if s & here]
        x = u['owners']
        if not targets or not ctx.same(ctx.apply('clear-sign', (x,)), ctx.baseline):
            continue
        because, fixes = [], []
        for y in targets:
            if y not in signed_copy:
                signed = ctx.apply('set-sign', (y,), (1,))
                signed_copy[y] = None if ctx.same(signed, ctx.baseline) else signed   # None: the drawing states `y`
            signed = signed_copy[y]
            if signed is None or ctx.same(ctx.apply('clear-sign', (x,), molecule=signed), signed):
                continue
            if y not in sign_fixes:
                sign_fixes[y] = tuple(ctx.later('set-sign', (y,), (p,)) for p in (1, 2))
            because.append(y)
            fixes.extend(sign_fixes[y])
        if because:
            out.append(_finding('dead-waits-for', (x,), group=stored_of.get(x), because=because, fixes=fixes))
    return out


def _logs_unrealisable(molecule):
    """Does an uncached write of `molecule` log `stereo:bridgehead-unrealisable`?

    Asked of a copy with its log cleared: `copy()` carries the log, and the writer logs the record once per
    log, so a record inherited from the source would answer for the source.  `copy()` also carries the
    cached `str()` and identity, which therefore log nothing; the slot-order spec `'i'` is never cached and
    always reads the coupled sets.
    """
    probe = molecule.copy()
    del probe.log[:]
    format(probe, 'i')
    return any(r.rule == _UNREALISABLE for r in probe.log)


def _label_strength(key):
    """How strongly a stored `(kind, id)` states a coupled set: AND, OR, ABS, none -- the core's order."""
    return 3 if key is None else _STRENGTH[key[0]]


def _bridgeheads(ctx):
    """Coupled bridgehead findings: unrealisable signs, implied members, and labels that disagree or are missing.

    A coupled set is one element and takes its sign and group from its free member: the signed member with
    the strongest stored label -- AND, then OR, then ABS, then none -- the first in unit order among equals.
    The set is labelled on every member: a member with no label beside a labelled set is
    `bridgehead-half-labelled`, another group of the same kind is `bridgehead-merged-groups`, which the
    set joins into one, and a weaker label is `bridgehead-redundant-label`.  A set is unrealisable when the core logs
    `stereo:bridgehead-unrealisable` on the baseline copy and on a copy with every other set's signs
    cleared.  A one-member `set-sign` fix is kept when that copy, so flipped, logs nothing; a set no single
    flip realises gets `clear-sign` over the set.  Every fix's `result` is the whole molecule fixed.
    """
    out = []
    units = {u['anchor']: u for u in ctx.units}
    order = {u['anchor']: i for i, u in enumerate(ctx.units)}
    stored_of = {x: key for key, members in ctx.molecule.stereo_groups().items() for x in members}
    sets, seen = [], set()
    for u in ctx.units:
        a = u['anchor']
        if not u['coupled'] or a in seen:
            continue
        members = tuple(sorted((a, *u['coupled'])))
        seen.update(members)
        signed = sorted((x for x in members if units[x]['parity']), key=order.__getitem__)
        if signed:
            sets.append((members, signed))

    conflict = _logs_unrealisable(ctx.baseline)
    for members, signed in sets:
        if conflict and len(signed) > 1:
            others = [x for _, s in sets if s is not signed for x in s]
            alone = ctx.apply('clear-sign', others) if others else None     # this set judged alone
            if alone is None or _logs_unrealisable(alone):
                atoms = tuple(sorted(signed))
                fixes = []
                for x in atoms:
                    args = (3 - units[x]['parity'],)
                    if not _logs_unrealisable(ctx.apply('set-sign', (x,), args, molecule=alone)):
                        fixes.append(ctx.fix('set-sign', (x,), args))
                if not fixes:
                    fixes.append(ctx.fix('clear-sign', atoms))
                out.append(_finding('bridgehead-unrealisable', atoms, fixes=fixes))
                continue
        source = min(signed, key=lambda x: (_label_strength(stored_of.get(x)), order[x]))
        for x in members:
            if units[x]['coupled_implied']:
                out.append(_finding('bridgehead-implied', (x,), because=(source,)))
        key = stored_of.get(source)
        for x in members:
            if x == source or (x not in stored_of and key is None):
                continue
            other = stored_of.get(x)
            args = None
            if key is not None:
                args = (_KIND_NAME[key[0]],) if key[0] == STEREO_ABS else (_KIND_NAME[key[0]], key[1])
            if other == key:
                continue
            if other is None:
                code, fixes = 'bridgehead-half-labelled', (ctx.fix('set-group', (x,), args),)
            elif key is not None and key[0] != STEREO_ABS and other[0] == key[0]:
                rest = tuple(sorted(y for y, k in stored_of.items() if k == other))
                code, fixes = 'bridgehead-merged-groups', (ctx.fix('merge-groups', rest, (_KIND_NAME[key[0]], key[1])),)
            else:
                code = 'bridgehead-redundant-label'
                fixes = (ctx.fix('set-group', (x,), args),) if key is not None else ()
                fixes += (ctx.fix('drop-label', (x,)),)
            fixes = tuple(f for f in fixes if f.preserves_meaning)
            if fixes:
                out.append(_finding(code, (x,), group=other, because=(source,), fixes=fixes))
    return out


_STRAINED_FUSIONS = {(3, 6), (4, 5)}     # trans is made but strained; smaller ones cannot be built (core)


def _frame_parity(refs, frame):
    """Parity of the permutation taking `(exo, *frame)` onto `refs`; None when a frame atom is not a ref."""
    perm = []
    for f in frame:
        if f not in refs:
            return None
        perm.append(refs.index(f))
    perm.insert(0, ({0, 1, 2, 3} - set(perm)).pop())
    return sum(perm[i] > perm[j] for i in range(4) for j in range(i + 1, 4)) & 1


def _fused_strained(ctx):
    """`fused-strained` on each trans-fused 3/6 or 4/5 pair: two signed fusion atoms sharing a bond.

    The frame is the core's for a coupled pair: (exo, then the three bridges ordered by their atom at `a`), the
    shared bond being the zero bridge, and out,out -- cis -- where the parity xor matches OUT_OUT_XOR (1).
    """
    molecule, out = ctx.molecule, []
    tetra = {u['anchor']: u for u in ctx.units if u['kind'] == 0 and u['parity'] and isinstance(u['owners'], int)}
    rings = [tuple(r) for r in molecule.sssr]
    for a in tetra:
        for b in molecule.neighbors_of(a):
            if b <= a or b not in tetra:
                continue
            through = sorted((r for r in rings if a in r and b in r
                              and abs(r.index(a) - r.index(b)) in (1, len(r) - 1)), key=len)
            if len(through) < 2 or (len(through[0]), len(through[1])) not in _STRAINED_FUSIONS:
                continue
            ua, ub = tetra[a], tetra[b]
            if sum(x is not None for x in ua['refs']) < 3 or sum(x is not None for x in ub['refs']) < 3:
                continue
            bridges = [(b, a)]
            for r in through[:2]:
                i, j = r.index(a), r.index(b)
                step = 1 if (i - j) % len(r) == 1 else -1           # walk from `a` away from `b`
                bridges.append((r[(i + step) % len(r)], r[(j - step) % len(r)]))
            bridges.sort()
            pa = _frame_parity(ua['refs'], [f for f, _ in bridges])
            pb = _frame_parity(ub['refs'], [g for _, g in bridges])
            if pa is None or pb is None:
                continue
            if ((ua['parity'] - 1) ^ (ub['parity'] - 1)) != 1 ^ pa ^ pb:          # not out,out: trans
                out.append(_finding('fused-strained', (a, b),
                                    fixes=(ctx.later('set-sign', (b,), (3 - ub['parity'],)),
                                           ctx.later('set-sign', (a,), (3 - ua['parity'],)))))
    return out


def _limits(molecule):
    """Molecule-level limit findings: identity-approximate and truncated."""
    out = []
    if molecule.identity_approximate:
        out.append(_finding('identity-approximate', ()))
    if molecule.stereo_truncated:
        out.append(_finding('truncated', ()))
    return out


def _pass_through(molecule):
    """Read-time stereo log records, those `_PASS_THROUGH` names, forwarded as findings without fixes."""
    out = []
    for record in molecule.log:
        if record.rule in _PASS_THROUGH:
            out.append(Finding(record.rule, _LOG_SEVERITY[record.severity], record.atoms, None, (),
                               record.message, ()))
    return out


def audit_stereo(molecule):
    """One report of what `molecule`'s stereo labels state; see `docs/stereo.rst`, Auditing labels.

    A reader: the molecule is not changed and `molecule.log` is not written.  `str()` and `==` are
    computed on a copy taken here, so even a bridgehead molecule that logs on `str()` stays clean.
    Every finding carries the code it describes and a tuple of `Fix` suggestions.  `StereoAudit.clean`
    is True when no finding is above `'info'` severity.

    Every fix and comparison is a copy of one constitution, so the whole audit runs in one
    `identity_memo()` scope: each phase combination's identity is computed once across all of them.
    """
    with identity_memo():
        return _audit(molecule)


def _audit(molecule):
    ctx = _Context(molecule)
    bridgeheads = _bridgeheads(ctx)
    # `bridgehead-unrealisable` alone names the members of an impossible pair
    impossible = {x for f in bridgeheads if f.code == 'stereo-audit:bridgehead-unrealisable' for x in f.atoms}
    # `bridgehead-redundant-label` names the label of its whole set
    relabelled = {x for f in bridgeheads if f.code in _SET_LABEL for x in (*f.atoms, *f.because)}
    waits = [f for f in _waits_for(ctx) if not impossible.intersection(f.atoms)]
    waiting = {f.atoms[0] for f in waits}
    awaited = {(c,) for f in waits for c in f.because}
    out = waits + bridgeheads + _fused_strained(ctx)
    for f in _units(ctx) + _groups(ctx):
        if impossible.intersection(f.atoms):
            continue
        if f.code in _DEAD and waiting.issuperset(f.atoms):
            continue                    # replaced by `dead-waits-for`
        if f.code == 'stereo-audit:unknown-unit' and f.atoms in awaited:
            continue                    # `dead-waits-for` names it and carries its set-sign fixes
        if f.code == 'stereo-audit:unlabelled-sign' and relabelled.intersection(f.atoms):
            continue                    # the set's label is already reported
        out.append(f)
    out.extend(_limits(molecule))
    out.extend(_pass_through(molecule))
    out.sort(key=lambda f: (_SEVERITY_ORDER[f.severity], _first(f.atoms), f.code))
    return StereoAudit(tuple(out))
