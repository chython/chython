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
"""Tripos MOL2 reader: ``@<TRIPOS>``-tagged sections, one record per ``MOLECULE`` tag.

:func:`read_mol2` yields ``(molecule, log)`` per record and files a
:class:`~chython.formats.ctfile.FailedRecord` for one that will not parse; :func:`mol2_mol` reads a
single record from a string and raises :class:`Mol2ParseError` instead.  Bond type ``ar`` is stored as
order 4 and nothing kekulizes it.  Log prefixes: ``atom``, ``bond``, ``record``, ``unsupported``.

Every record's lines land on the molecule that record produced, under ``mol.log`` at stage ``read``, and
the yielded or caller-supplied list is a second copy of the same lines rather than the storage.
"""
from collections.abc import Generator
from io import StringIO
from pathlib import Path

from ._text import require_text
from .ctfile import FailedRecord
from ..core import LogRecord, LOST, MoleculeContainer, REPAIRED
from ..core._core import element_symbols


__all__ = ['Mol2ParseError', 'mol2', 'mol2_mol', 'read_mol2']


class Mol2ParseError(Exception):
    """The MOL2 record cannot be parsed: the structure is unreadable, not merely wrong.

    Raised only where a best-effort read would invent data (an ATOM line with fewer than the six
    required fields).  :func:`read_mol2` catches it and files a
    :class:`~chython.formats.ctfile.FailedRecord`; :func:`mol2_mol` lets it out.
    """


_ELEMENTS: frozenset = frozenset(element_symbols()[1:])   # index 0 is R, the fragment marker

_TAG = '@<TRIPOS>'

# Sections this reader reads.  Every other one is named in an `unsupported:` line.
_CLAIMED_SECTIONS = frozenset({'MOLECULE', 'ATOM', 'BOND'})

# A section whose loss is worth naming specifically: UNITY/SYBYL writers put formal charges in
# UNITY_ATOM_ATTR, so the user whose charges vanished needs to be told which section held them.
_SECTION_CONSEQUENCE: dict[str, str] = {
    'UNITY_ATOM_ATTR': 'formal charges written there are not read; every atom keeps the charge '
                       'column of the ATOM block',
}

# The Tripos charge types that store continuous partial charges rather than integer formal ones --
# the specification's list minus NO_CHARGES.  A record declaring one of these has a charge column
# that must not be rounded to a formal charge.  The declared type is upper-cased before the test,
# a lower-cased spelling being a writer's habit and not a statement that the charges are formal.

_PARTIAL_CHARGE_TYPES = frozenset({
    'DEL_RE', 'GASTEIGER', 'GAST_HUCK', 'HUCKEL', 'PULLMAN',
    'GAUSS80_CHARGES', 'AMPAC_CHARGES', 'MULLIKEN_CHARGES',
    'DICT_CHARGES', 'MMFF94_CHARGES', 'USER_CHARGES',
})

# MOL2 bond type strings -> internal order, or None for "skip this bond entirely".  A type not in
# this map is logged and stored as single rather than dropped.

_BOND_ORDER: dict[str, int | None] = {
    '1': 1,
    '2': 2,
    '3': 3,
    '4': 4,    # some writers use '4' for aromatic; not a Tripos type, so it is logged
    'ar': 4,   # standard aromatic
    'am': 1,   # amide: single bond in graph terms; logged as unsupported
    'un': 1,   # unknown order: single is the least-wrong default; logged as unsupported
    'du': None,  # dummy bond: skip
    'nc': None,  # not connected: skip
}

# What a legal Tripos bond type this library cannot represent costs the reader: all four of them.
_BOND_NOTE: dict[str, str] = {
    'am': 'bond type "am" (amide) is not modelled as a distinct order; stored as single',
    'un': 'bond type "un" (order unknown) has no representation; stored as single',
    'du': 'bond type "du" (dummy bond) is not modelled; the bond is not stored',
    'nc': 'bond type "nc" (not connected) is not modelled; the bond is not stored',
}

# A bond type the format does not define but writers emit anyway: the file is wrong and the reader
# copes, so the line takes no `unsupported:` prefix.
_NONSTANDARD_BOND_NOTE: dict[str, str] = {
    '4': 'bond type "4" is not a Tripos bond type; read as aromatic, which is how the writers '
         'that emit it spell "ar"',
}

# SYBYL types whose suffix states a coordination geometry rather than a hybridization.  The TSV's
# only vocabulary is chython's hybridization codes, so these rows carry 0 -- the same 0 hydrogen
# carries, which is why the reason is named here and not reconstructed from the value.

_GEOMETRY_ONLY_TYPES: dict[str, str] = {
    'Cr.th': 'tetrahedral', 'Cr.oh': 'octahedral', 'Co.oh': 'octahedral', 'Ru.oh': 'octahedral',
}


# --- intermediate atom / bond records ----------------------------------------------------- #

class _Atom:
    """One line from the ATOM block, in the file's own terms."""
    __slots__ = ('file_id', 'lineno', 'x', 'y', 'z', 'sybyl_type', 'element', 'hybridization',
                 'charge')

    def __init__(self, file_id, lineno, x, y, z, sybyl_type, element, hybridization, charge):
        self.file_id = file_id          # the file's atom id, or None when it was unreadable
        self.lineno = lineno            # 1-based line position inside the ATOM block
        self.x = x
        self.y = y
        self.z = z
        self.sybyl_type = sybyl_type
        self.element = element          # resolved element symbol, or None for pseudo-atoms
        self.hybridization = hybridization   # int 0-5 from the type; checked against the bonds
        self.charge = charge            # int formal charge (0 when partial or absent)


class _Bond:
    """One line from the BOND block, in the file's own terms."""
    __slots__ = ('file_id', 'lineno', 'a', 'b', 'order', 'sybyl_type')

    def __init__(self, file_id, lineno, a, b, order, sybyl_type):
        self.file_id = file_id
        self.lineno = lineno
        self.a = a        # 0-based index into the atom list
        self.b = b
        self.order = order       # None means "skip"
        self.sybyl_type = sybyl_type


# --- one name per thing ------------------------------------------------------------------- #

def _atom_ref(lineno: int, file_id: int | None) -> str:
    """How an atom is named in every log line: its ATOM-block line and the id the file gave it."""
    return f'atom line {lineno} (id {file_id})' if file_id is not None \
        else f'atom line {lineno} (id unreadable)'


def _bond_ref(lineno: int, file_id: int | None) -> str:
    """How a bond is named in every log line: its BOND-block line and the id the file gave it."""
    return f'bond line {lineno} (id {file_id})' if file_id is not None \
        else f'bond line {lineno} (id unreadable)'


def _section_tag(line: str) -> str | None:
    """The section name of a ``@<TRIPOS>`` line, upper-cased, or ``None`` for a content line.

    Leading whitespace is stripped before the test, so an indented tag starts a section.  The one
    place this question is asked: a second spelling of it lets a tag start a record and not a section.
    """
    stripped = line.lstrip()
    if stripped.upper().startswith(_TAG):
        return stripped[len(_TAG):].strip().upper()
    return None


# --- SYBYL type resolution ---------------------------------------------------------------- #

def _resolve_type(sybyl_type: str) -> tuple[str | None, int, str | None]:
    """Resolve a SYBYL type string (``C.ar``, ``N.3``, ``O.co2``) to ``(element, hybridization, note)``.

    The type is the file's only statement about hybridization; the mapping is
    ``chemistry/tables/sybyl_types.tsv``.  ``element`` is ``None`` for a pseudo-atom (LP, Du) with no
    nucleus and ``''`` for a type not recognised at all, whose atom cannot be stored.
    ``hybridization`` is a chython code 0-5, 0 when the type states none.  ``note`` is a reason and not a
    log line -- the caller prefixes it and names the atom -- or ``None``.
    """
    from ..chemistry._tables import sybyl_types
    table = sybyl_types()

    row = table.get(sybyl_type)
    if row is not None:
        geometry = _GEOMETRY_ONLY_TYPES.get(sybyl_type)
        note = None if geometry is None else (
            f'SYBYL type {sybyl_type!r} states {geometry} coordination, which chython has no '
            f'hybridization code for')
        # Recognised type.  An empty element means a pseudo-atom.
        return (None if row.element == '' else row.element), row.hybridization, note

    # Not in the table.  Try to split on the first dot and use the prefix as the element.
    if '.' in sybyl_type:
        prefix, _, suffix = sybyl_type.partition('.')
        # Normalise to title case for element lookup (Mol2 files sometimes uppercase all).
        candidate = prefix.title() if len(prefix) > 1 else prefix.upper()
        if candidate in _ELEMENTS:
            return candidate, 0, (f'SYBYL type {sybyl_type!r} is not in the type table; element '
                                  f'{candidate} was read from its prefix and the {suffix!r} tag '
                                  f'was not interpreted')
        # The part before the dot is not an element either.  Fall through.

    # Try the whole string as an element symbol (bare element, no dot).
    candidate = sybyl_type.title() if len(sybyl_type) > 1 else sybyl_type.upper()
    if candidate in _ELEMENTS:
        return candidate, 0, None   # a bare element symbol states no hybridization to lose

    # Completely unknown.
    return '', 0, (f'SYBYL type {sybyl_type!r} is not in the type table and its prefix is not an '
                   f'element symbol; the atom is not stored')


# --- section splitter --------------------------------------------------------------------- #

def _split_sections(lines: list[str]) -> list[tuple[str, list[str]]]:
    """Split a record's lines into ``[(section_name, [content_lines]), ...]``, in file order.

    ``section_name`` is the tag after ``@<TRIPOS>`` in upper case; lines before the first tag come
    back under the empty name, which is where bare MOLECULE content lands.  A list and not a dict: a
    dict cannot hold two sections of one name, and a record with two ATOM blocks has them.
    """
    sections: list[tuple[str, list[str]]] = [('', [])]
    for line in lines:
        stripped = line.rstrip('\r\n')
        name = _section_tag(stripped)
        if name is None:
            sections[-1][1].append(stripped)
        else:
            sections.append((name, []))
    return sections


def _split_records(lines: list[str]) -> list[list[str]]:
    """Split lines into records at every ``@<TRIPOS>MOLECULE`` tag, tag line excluded.

    Lines before the first tag belong to no record -- a MOL2 file may open with ``#`` comments.  Input
    with no tag at all is bare MOLECULE content, from a caller who stripped the tag, and is one record.
    """
    records: list[list[str]] = []
    current: list[str] | None = None
    for line in lines:
        if _section_tag(line) == 'MOLECULE':
            current = []
            records.append(current)
        elif current is not None:
            current.append(line)
    if not records and any(x.strip() for x in lines):
        return [list(lines)]
    return records


# --- MOLECULE section parser -------------------------------------------------------------- #

def _parse_molecule(lines: list[str],
                    log: list[str]) -> tuple[str, str, int | None, int | None]:
    """Parse MOLECULE content → ``(title, charge_type, num_atoms, num_bonds)``.

    The section is positional and is read positionally: line 1 the name, 2 the counts, 3 the molecule
    type, 4 the charge type.  No blank line may be dropped before indexing -- an empty name line is the
    commonest placeholder there is, and skipping it shifts the charge type off the end.
    ``num_atoms``/``num_bonds`` are ``None`` when the file stated no count, which is not a count of zero.
    Both are advisory: the ATOM and BOND blocks win and a discrepancy is logged.
    """
    title = lines[0].strip() if lines else ''

    num_atoms: int | None = None
    num_bonds: int | None = None
    if len(lines) > 1:
        parts = lines[1].split()
        try:
            num_atoms = int(parts[0])
        except (IndexError, ValueError):
            log.append(LogRecord('mol2:unparseable-counts', (),
                                 'record: counts line is not parseable; atom/bond count checks skipped'))
        if num_atoms is not None:
            if len(parts) > 1:
                try:
                    num_bonds = int(parts[1])
                except ValueError:
                    log.append(LogRecord('mol2:unparseable-counts', (),
                                         f'record: counts line bond field is {parts[1]!r}, not a number; '
                                         f'bond count check skipped'))
            # num_subst / num_feat / num_sets: a non-zero one is a construct the container drops, a
            # zero states nothing.
            extras = []
            for name, index in (('substructures', 2), ('features', 3), ('sets', 4)):
                if len(parts) > index:
                    try:
                        value = int(parts[index])
                    except ValueError:
                        log.append(LogRecord('mol2:unparseable-counts', (),
                                             f'record: counts line field for {name} is {parts[index]!r}, '
                                             f'not a number'))
                        continue
                    if value:
                        extras.append(f'{value} {name}')
            if extras:
                log.append(LogRecord('mol2:unsupported-subsystem-counts', (),
                                     f'unsupported: MOLECULE counts line states {", ".join(extras)}; the '
                                     f'substructure, feature and set model has no storage in chython', LOST))

    charge_type = lines[3].strip() if len(lines) > 3 else ''
    return title, charge_type, num_atoms, num_bonds


# --- ATOM section parser ------------------------------------------------------------------ #

def _parse_atoms(
    lines: list[str], charge_type: str, log: list[str]
) -> tuple[list[_Atom], int, dict[int, int], dict[int, str]]:
    """Parse ATOM section content lines → ``(atoms, total_atom_lines, id_to_index, dropped)``.

    ``total_atom_lines`` counts the block's non-blank lines, which is what the MOLECULE header claim is
    checked against -- never the atoms that were kept.  ``id_to_index`` maps the file's atom ids of the
    kept atoms to 0-based positions; an id is claimed by the first line stating it, a second one logged,
    since a bond names its endpoints by id.  ``dropped`` is ``{file_id: reason}`` for a claimed id whose
    atom was not stored, which is how the bond parser tells our limitation from the file's error.  Raises
    :class:`Mol2ParseError` only for a line with fewer than the six required fields (atom_id, atom_name,
    x, y, z, atom_type), which cannot even identify the atom.
    """
    use_partial = charge_type.upper() in _PARTIAL_CHARGE_TYPES
    partial_logged = False
    unstated_logged = False
    substructures: set = set()
    status_bits = False
    atoms: list[_Atom] = []
    total_atom_lines = 0
    claimed: dict[int, str] = {}
    id_to_index: dict[int, int] = {}
    dropped: dict[int, str] = {}

    for lineno, line in enumerate(lines, 1):
        if not line.strip():
            continue
        total_atom_lines += 1
        parts = line.split()
        if len(parts) < 6:
            raise Mol2ParseError(
                f'ATOM block line {lineno} has {len(parts)} fields (need at least 6): {line!r}')

        file_id: int | None
        try:
            file_id = int(parts[0])
        except ValueError:
            file_id = None
        ref = _atom_ref(lineno, file_id)
        if file_id is None:
            # No invented id: a line number in the id namespace collides with a real id.
            log.append(LogRecord('mol2:unreadable-atom-id', (),
                                 f'{ref}: atom id {parts[0]!r} is not an integer; no bond can reference '
                                 f'this atom'))
        elif file_id in claimed:
            log.append(LogRecord('mol2:duplicate-atom-id', (),
                                 f'{ref}: atom id {file_id} was already stated by {claimed[file_id]}; bonds '
                                 f'naming it bind to the first, and this atom is unreachable'))
            file_id = None   # stored, but unnameable: the first claim keeps the id
        else:
            claimed[file_id] = ref

        # Substructure annotation, reported once for the block and only when it names more than one
        # substructure: the subst columns are present in practically every file, and `1 LIG` on every
        # line of a one-residue ligand states nothing the container flattens.  The mandatory
        # ``atom_name`` is a property of the format pairing and gets no per-record line at all.
        if len(parts) >= 8:
            substructures.add((parts[6], parts[7]))
        if len(parts) >= 10:
            status_bits = True

        try:
            x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
        except ValueError:
            log.append(LogRecord('mol2:bad-coordinates', (),
                                 f'{ref}: coordinate fields are not numbers in {line!r}; stored as 0 0 0',
                                 REPAIRED))
            x = y = z = 0.0

        sybyl_type = parts[5]
        element, hybridization, note = _resolve_type(sybyl_type)
        if note is not None:
            # Every reason `_resolve_type` gives is a limit of our table, not a defect in the file.
            log.append(LogRecord('mol2:unknown-sybyl-type', (),
                                 f'unsupported: {ref}: {note}', LOST))

        # Pseudo-atom: no nucleus, so no graph atom.  The file is fine, hence `unsupported:`.
        if element is None:
            log.append(LogRecord('mol2:pseudo-atom', (),
                                 f'unsupported: {ref}: type {sybyl_type!r} is a pseudo-atom (no nucleus); '
                                 f'not stored', LOST))
            if file_id is not None:
                dropped[file_id] = f'a pseudo-atom of type {sybyl_type!r}'
            continue

        # Unrecognised type: `_resolve_type` has already said so, in its own words.
        if element == '':
            if file_id is not None:
                dropped[file_id] = f'an unrecognised SYBYL type {sybyl_type!r}'
            continue

        # Formal charge from the optional charge column.
        charge = 0
        if len(parts) >= 9:
            if use_partial:
                if not partial_logged:
                    log.append(LogRecord('mol2:partial-charges-discarded', (),
                                         f'unsupported: partial charges (charge_type={charge_type!r}) '
                                         f'are not stored; formal charges set to 0', LOST))
                    partial_logged = True
            else:
                if not charge_type and not unstated_logged:
                    # No charge type stated, so reading the column as formal is our assumption.
                    log.append(LogRecord('mol2:no-charge-type', (),
                                         'record: MOLECULE states no charge type; the ATOM block charge '
                                         'column is read as formal charges'))
                    unstated_logged = True
                try:
                    raw = float(parts[8])
                    charge = int(round(raw))
                except ValueError:
                    log.append(LogRecord('mol2:bad-charge', (),
                                         f'{ref}: charge field {parts[8]!r} is not a number; stored as 0',
                                         REPAIRED))

        if file_id is not None:
            id_to_index[file_id] = len(atoms)
        atoms.append(_Atom(file_id, lineno, x, y, z, sybyl_type, element, hybridization, charge))

    if len(substructures) > 1:
        log.append(LogRecord('mol2:substructure-annotation', (),
                             f'unsupported: the ATOM block assigns its atoms to {len(substructures)} '
                             f'substructures (the subst_id and subst_name columns); chython stores no '
                             f'residue annotation', LOST))
    if status_bits:
        log.append(LogRecord('mol2:status-bits', (),
                             'unsupported: the ATOM block states status_bit values; chython stores none',
                             LOST))

    return atoms, total_atom_lines, id_to_index, dropped


# --- BOND section parser ------------------------------------------------------------------ #

def _parse_bonds(lines: list[str], id_to_index: dict[int, int], dropped: dict[int, str],
                 log: list[str]) -> tuple[list[_Bond], int]:
    """Parse BOND section content lines → ``(bonds, total_bond_lines)``.

    A bond to an id in ``dropped`` -- seen in the ATOM block and not stored -- is our limitation and
    takes the ``unsupported:`` prefix; a bond to an id in neither map is the file's error and takes the
    plain ``bond`` prefix.  ``total_bond_lines`` counts the block's non-blank lines, which is what the
    MOLECULE header claim is checked against, a legal ``du``/``nc`` bond reducing the stored count
    without making the header wrong.  The bond type is resolved before the duplicate check, a ``du``
    bond over an already-bonded pair being a dummy bond rather than a duplicate.
    """
    bonds: list[_Bond] = []
    total_bond_lines = 0
    seen: set = set()

    for lineno, line in enumerate(lines, 1):
        if not line.strip():
            continue
        total_bond_lines += 1
        parts = line.split()
        if len(parts) < 4:
            log.append(LogRecord('mol2:malformed-bond-line', (),
                                 f'bond line {lineno}: line has {len(parts)} fields (need at least 4), '
                                 f'skipped'))
            continue

        file_id: int | None
        try:
            file_id = int(parts[0])
        except ValueError:
            file_id = None
        ref = _bond_ref(lineno, file_id)
        if file_id is None:
            log.append(LogRecord('mol2:unreadable-bond-id', (),
                                 f'{ref}: bond id {parts[0]!r} is not an integer'))

        sybyl_bond = parts[3]
        if sybyl_bond in _BOND_ORDER:
            order = _BOND_ORDER[sybyl_bond]
            note = _BOND_NOTE.get(sybyl_bond)
            if note is not None:
                log.append(LogRecord('mol2:unsupported-bond-type', (),
                                     f'unsupported: {ref}: {note}', LOST))
            nonstandard = _NONSTANDARD_BOND_NOTE.get(sybyl_bond)
            if nonstandard is not None:
                log.append(LogRecord('mol2:nonstandard-bond-type', (),
                                     f'{ref}: {nonstandard}', REPAIRED))
        else:
            log.append(LogRecord('mol2:unknown-bond-type', (),
                                 f'{ref}: bond type {sybyl_bond!r} is not a known type, stored as single',
                                 REPAIRED))
            order = 1

        if order is None:
            continue                     # du (dummy) or nc: the note above is the whole report

        try:
            aid_a = int(parts[1])
            aid_b = int(parts[2])
        except ValueError:
            log.append(LogRecord('mol2:unreadable-bond-atoms', (),
                                 f'{ref}: atom ids {parts[1]!r} / {parts[2]!r} are not integers, skipped'))
            continue

        idx_a = id_to_index.get(aid_a)
        idx_b = id_to_index.get(aid_b)
        if idx_a is None or idx_b is None:
            missing_id = aid_a if idx_a is None else aid_b
            reason = dropped.get(missing_id)
            if reason is not None:
                log.append(LogRecord('mol2:bond-to-unstored-atom', (),
                                     f'unsupported: {ref}: endpoint atom {missing_id} is {reason} and was '
                                     f'not stored; the bond is not stored', LOST))
            else:
                log.append(LogRecord('mol2:bond-missing-atom', (),
                                     f'{ref}: references atom id {missing_id} which is not in the ATOM '
                                     f'block, skipped'))
            continue
        if idx_a == idx_b:
            log.append(LogRecord('mol2:self-loop-bond', (),
                                 f'{ref}: self-loop on atom id {aid_a}, skipped'))
            continue

        key = (min(idx_a, idx_b), max(idx_a, idx_b))
        if key in seen:
            log.append(LogRecord('mol2:duplicate-bond', (),
                                 f'{ref}: duplicate of an earlier bond between the same two atoms, skipped'))
            continue
        seen.add(key)

        bonds.append(_Bond(file_id, lineno, idx_a, idx_b, order, sybyl_bond))

    return bonds, total_bond_lines


# --- molecule builder --------------------------------------------------------------------- #

def _build(title: str, atoms: list[_Atom], bonds: list[_Bond],
           log: list[str]) -> MoleculeContainer:
    """Build a :class:`~chython.core.MoleculeContainer` from the parsed intermediate lists.

    The ctfile pattern: atoms and bonds in one edit scope, coordinates in a second, then implicit
    hydrogen counts outside any scope, which is where the arena is sealed.
    """
    mol = MoleculeContainer()
    sids: list[int] = []

    with mol.edit():
        for a in atoms:
            ref = _atom_ref(a.lineno, a.file_id)
            full = {'charge': a.charge}
            for drop in ((), ('charge',)):
                kwargs = {k: v for k, v in full.items() if k not in drop}
                try:
                    sid = mol.add_atom(a.element, **kwargs)
                except ValueError as e:
                    reason = e
                    continue
                if drop:
                    log.append(LogRecord('mol2:atom-charge-dropped', (sid,),
                                         f'{ref} {a.element}: {reason}; dropped {", ".join(drop)}',
                                         REPAIRED))
                break
            else:
                log.append(LogRecord('mol2:atom-unstorable', (),
                                     f'{ref} {a.element}: cannot be stored even without charge; skipped',
                                     LOST))
                sids.append(-1)   # sentinel so bond indexing stays aligned
                continue
            sids.append(sid)

        n = len(sids)
        for bond in bonds:
            sa = sids[bond.a] if bond.a < n else -1
            sb = sids[bond.b] if bond.b < n else -1
            if sa == -1 or sb == -1:
                log.append(LogRecord('mol2:bond-skipped', (),
                                     f'{_bond_ref(bond.lineno, bond.file_id)}: one endpoint was not stored, '
                                     f'skipped'))
                continue
            mol.add_bond(sa, sb, bond.order)

    pairs = [(sid, a) for sid, a in zip(sids, atoms) if sid != -1]

    # Coordinates, only when at least one atom has a non-zero one.  MOL2 is a 3D format, so both
    # segments are filled: `SEG_XY` is what a depiction reads, `SEG_CONFORMERS` the stated geometry.
    if pairs and any(a.x or a.y or a.z for _, a in pairs):
        solid = any(a.z for _, a in pairs)
        with mol.edit():
            for sid, a in pairs:
                try:
                    mol.set_xy(sid, a.x, a.y)
                    if solid:
                        mol.set_xyz(sid, a.x, a.y, a.z)
                except ValueError as e:
                    log.append(LogRecord('mol2:coordinates-dropped', (sid,),
                                         f'coordinates for {_atom_ref(a.lineno, a.file_id)} dropped: {e}',
                                         LOST))

    # Implicit hydrogen counts via the chemistry layer (lazy import, one call per atom).
    from ..chemistry._implicit import calc_implicit
    for sid, _ in pairs:
        calc_implicit(mol, sid)

    # Hybridization is derived from the bonds, so the SYBYL type's claim is a check: where it
    # contradicts the bonds the record also wrote, the bonds win and the difference is logged.
    for sid, a in pairs:
        if a.hybridization and mol.hybridization_of(sid) != a.hybridization:
            log.append(LogRecord('mol2:hybridization-mismatch', (sid,),
                                 f'{_atom_ref(a.lineno, a.file_id)}: SYBYL type {a.sybyl_type!r} states '
                                 f'hybridization {a.hybridization} but its bonds give '
                                 f'{mol.hybridization_of(sid)}; the bonds are used'))

    if title:
        mol.set_title(title)

    return mol


# --- record parser ------------------------------------------------------------------------ #

def _parse_record(lines: list[str], log: list[str]) -> MoleculeContainer:
    """Parse one MOL2 record (the lines after its ``@<TRIPOS>MOLECULE`` tag).

    Raises :class:`Mol2ParseError` when the ATOM block is structurally unreadable.  Every other
    malformation is logged and the best-effort molecule is returned.

    Every line goes to this record's own list first, which is then absorbed onto the molecule's ``log``
    and extended onto the caller's.  ``mol.log`` is the storage and the absorb is unconditional; the
    private list is what keeps record 3's lines off molecule 4 whatever the caller's list already holds.
    """
    own: list = []
    try:
        mol = _read_sections(lines, own)
    finally:
        log.extend(own)   # even when the record raised: what was found before it is still the answer
    mol.log.absorb('read', own)
    return mol


def _read_sections(lines: list[str], log: list[str]) -> MoleculeContainer:
    """The record's sections, in the file's order, as a molecule.  Logs to *log* and nowhere else."""
    mol_lines: list[str] = []
    atom_lines: list[str] = []
    bond_lines: list[str] = []
    unclaimed: dict[str, int] = {}

    for name, content in _split_sections(lines):
        if name in ('', 'MOLECULE'):
            target = mol_lines
        elif name == 'ATOM':
            target = atom_lines
        elif name == 'BOND':
            target = bond_lines
        else:
            unclaimed[name] = unclaimed.get(name, 0) + 1
            continue
        if target and any(x.strip() for x in content):
            log.append(LogRecord('mol2:duplicate-section', (),
                                 f'record: a second {name or "MOLECULE"} section in one record; its lines '
                                 f'are read as a continuation of the first'))
        target.extend(content)

    # Every section this reader does not read, by name.
    for name, times in unclaimed.items():
        how_many = 'section is' if times == 1 else f'{times} sections are'
        consequence = _SECTION_CONSEQUENCE.get(name)
        if consequence is None:
            log.append(LogRecord('mol2:unsupported-section', (),
                                 f'unsupported: the {name!r} {how_many} not read', LOST))
        else:
            log.append(LogRecord('mol2:unsupported-section', (),
                                 f'unsupported: the {name!r} {how_many} not read: {consequence}', LOST))

    title, charge_type, num_atoms, num_bonds = _parse_molecule(mol_lines, log)

    atoms, total_atom_lines, id_to_index, dropped = _parse_atoms(atom_lines, charge_type, log)

    # Both count checks compare the header's claim against the block's line count, never against what
    # was stored: a pseudo-atom or a `du` bond reduces the stored count without making the header wrong.
    if num_atoms is not None and total_atom_lines != num_atoms:
        log.append(LogRecord('mol2:count-mismatch', (),
                             f'record: MOLECULE header claims {num_atoms} atoms but the ATOM block has '
                             f'{total_atom_lines} lines; the block is used'))

    bonds, total_bond_lines = _parse_bonds(bond_lines, id_to_index, dropped, log)

    if num_bonds is not None and total_bond_lines != num_bonds:
        log.append(LogRecord('mol2:count-mismatch', (),
                             f'record: MOLECULE header claims {num_bonds} bonds but the BOND block has '
                             f'{total_bond_lines} lines; the block is used'))

    return _build(title, atoms, bonds, log)


# --- public entry points ------------------------------------------------------------------ #

def mol2_mol(data, *, log: list[str] | None = None) -> MoleculeContainer:
    """Parse one MOL2 record from a string or a list of lines.

    A ``@<TRIPOS>MOLECULE`` header line is stripped, so bare content and a complete record are both
    accepted.  A string holding several records gives back the **first**, with a ``record:`` line naming
    how many there were -- the policy :func:`~chython.formats.ctfile.mol` follows for a multi-record
    SDF string.  ``log`` receives damage messages and so does the returned molecule's own ``log``.  Raises
    :class:`Mol2ParseError` when the ATOM block is structurally unreadable: this entry point was asked for
    one molecule and has no next record.
    """
    if log is None:
        log = []
    lines = data.split('\n') if isinstance(data, str) else list(data)
    records = _split_records([x.rstrip('\r\n') for x in lines]) or [[]]
    mol = _parse_record(records[0], log)
    if len(records) > 1:
        # Logged after the parse rather than before it: the sentence is about the molecule that came
        # back, so it goes on that molecule's own log too.
        extra = [LogRecord('mol2:multiple-records', (),
                           f'record: the input holds {len(records)} MOL2 records; the first is returned '
                           f'and read_mol2() is the call that yields them all')]
        mol.log.absorb('read', extra)
        log.extend(extra)
    return mol


def mol2(data, *, log=None):
    """Every molecule in a MOL2 document.

    :param data: MOL2 text.  ALWAYS text -- unlike :func:`read_mol2`, which opens a path, this reads
        what it is given, so a document is never mistaken for a filename.
    :param log: a list to append damage reports to.  Every record's own lines also land on the molecule
        it produced, under ``mol.log``, whether or not this is passed.

    Answers a list.  A record chython cannot build becomes a `FailedRecord` in its place rather than
    raising: one damaged record in a file does not hide the rest.
    """
    data = require_text(data, 'mol2')
    log = [] if log is None else log
    records = []
    # `_iter_records` yields `(record, its own log)`; the per-record logs are flattened into the
    # caller's, so the reader's own `mol2:parse-failure` sentence is the one a failure reports.
    for record, record_log in _iter_records(StringIO(data), list, True):
        log.extend(record_log)
        records.append(record)
    return records


def _read_one(lines: list[str], log: list[str], position: int):
    """One record as ``(molecule, log)``, or as ``(FailedRecord, log)`` when it will not parse.

    One handler for both a bug of ours and damage in the file: either way the caller gets a record to
    look at.
    """
    try:
        return _parse_record(lines, log), log
    except Exception as e:
        log.append(LogRecord('mol2:parse-failure', (),
                             f'record: this record could not be parsed and holds no molecule: {e}', LOST))
        return FailedRecord(position, lines, e), log


def _iter_records(stream, log_factory, owned: bool) -> Generator[tuple[object, list[str]],
                                                                 None, None]:
    """Yield one ``(molecule-or-FailedRecord, log)`` pair per ``@<TRIPOS>MOLECULE`` in *stream*."""
    position = -1
    record_lines: list[str] = []
    in_record = False
    try:
        for raw_line in stream:
            line = raw_line.rstrip('\r\n')
            if _section_tag(line) == 'MOLECULE':
                if in_record and record_lines:
                    position += 1
                    yield _read_one(record_lines, log_factory(), position)
                record_lines = []
                in_record = True
            elif in_record:
                record_lines.append(line)

        # Last record (no trailing MOLECULE tag to flush it).
        if in_record and record_lines:
            position += 1
            yield _read_one(record_lines, log_factory(), position)
    finally:
        if owned:
            stream.close()


def read_mol2(source, *, log_factory=None) -> Generator[tuple[object, list[str]], None, None]:
    """Yield ``(molecule, log)`` for every MOL2 record in *source*.

    *source* is a file path (``str`` or :class:`pathlib.Path`), a file-like object opened in text mode, or
    a string of MOL2 text; ``@<TRIPOS>MOLECULE`` lines separate the records.  A ``str`` with no
    ``@<TRIPOS>`` in it is a path and is opened before this call returns, so a misspelled filename raises
    ``FileNotFoundError`` at the call site for both spellings of a path.  A record that cannot be parsed
    does not end the iteration: its place is taken by a :class:`~chython.formats.ctfile.FailedRecord`
    carrying the lines and the error.  *log_factory* returns a fresh log list per record.
    """
    if log_factory is None:
        log_factory = list

    if isinstance(source, Path):
        return _iter_records(source.open(encoding='utf-8', errors='replace'), log_factory, True)
    if isinstance(source, str):
        if _TAG in source.upper():
            return _iter_records(StringIO(source), log_factory, True)
        return _iter_records(open(source, encoding='utf-8', errors='replace'), log_factory, True)
    # A file-like object the caller opened: iterate it and leave it open.
    return _iter_records(source, log_factory, False)
