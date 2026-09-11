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
"""Read the checked-in tables in `tables/` and hand back compiled rules.

Compilation is lazy and cached: `import chython.chemistry` must stay cheap.  This module and `_smarts.py`
may not import the passes beside them, which `test/test_dependency_direction.py` gates.  Every table is
lexed by the core's `read_smarts`, the standardization tables through `compile_smarts` for its numbering.
"""
from collections.abc import Sequence
from importlib.resources import files
from typing import Any, NamedTuple
from ._smarts import compile_smarts
from ..core import QueryContainer, read_smarts, read_smiles


__all__ = ['ACID_ROLES', 'AbbreviationRow', 'AcidRow', 'CRIPPEN_CATCH_ALLS', 'CRIPPEN_ROLES',
           'CRIPPEN_TYPES', 'CrippenRow', 'Endpoint', 'HBOND_ROLES', 'HBondRow',
           'MACCS_EXPECTATIONS', 'MACCS_KINDS', 'MACCS_PREDICATES', 'MACCS_UNSET_KEYS',
           'MaccsCorpusRow', 'MaccsRow', 'PHARMACOPHORE_ROLES', 'PharmacophoreRow', 'QedAlertRow',
           'RESONANCE_ROLES', 'ROTATABLE_ROLES', 'Rule',
           'RotatableRow', 'SALT_ROLES', 'SaltRow', 'SybylType', 'TPSA_CLASSES', 'TpsaRow',
           'abbreviation_row', 'abbreviations_by_label', 'abbreviations_rows',
           'acids_rules', 'acids_rules_by_role', 'acids_table_text',
           'covalent_radii', 'crippen_rules', 'crippen_rules_by_role', 'first_match',
           'groups_rules', 'hbond_rules', 'hbond_rules_by_role', 'maccs_corpus',
           'maccs_corpus_by_key', 'maccs_rules', 'maccs_rules_by_key', 'metals_rules',
           'pharmacophore_rules', 'pharmacophore_rules_by_role', 'qed_alerts', 'read_table',
           'resonance_rules', 'resonance_rules_by_role', 'resonance_table_text', 'rotatable_rules',
           'rotatable_rules_by_role', 'salts_rows', 'salts_rows_by_role',
           'salts_species_keys', 'salts_table_text', 'standardize_rules', 'sybyl_types',
           'tpsa_rules']


class Rule(NamedTuple):
    """One repair, compiled.

    `id` is the row's table-qualified identity (`groups:13`) and travels in every log record.
    `numbers` maps the patch's atom number to the core query's stable id.  `atom_fix` is
    `[(number, charge_delta, radical_or_None)]` in file order, which decides the order of a
    rejection message.  `anchors` names the atoms two matches of this rule may share (see
    `_anchors`).  `why` is the message a log record carries; `after` names earlier rows of the same
    table this one must follow; `examples` is one `IN>>OUT` per `,` alternative of the pattern.
    """
    id: str
    query: QueryContainer
    numbers: dict[int, int]
    atom_fix: tuple[tuple[int, int, bool | None], ...]
    bonds_fix: tuple[tuple[int, int, int], ...]
    tautomer: bool
    after: tuple[str, ...]
    examples: tuple[str, ...]
    why: str
    smarts: str
    anchors: frozenset[int]


_RADICAL = {'-': None, '0': False, '1': True}
_RULES_CACHE: dict[str, tuple[Rule, ...]] = {}


def read_table(name: str) -> list[dict[str, str]]:
    """Parse one TSV out of `tables/` into a list of row dicts, without compiling anything.

    `name` is the bare filename; the `tables/` prefix is added here.  Keep that prefix inside the
    one `joinpath` literal: `chython/test/test_packaging.py` scans for the whole relative path to
    check the table is declared in `[tool.setuptools.package-data]`, and an undeclared table is
    missing only from the wheel, never from a checkout.
    """
    try:
        # `encoding='utf-8'`: the tables are shipped bytes and the codec that reads them may not be the
        # user's locale -- `read_text` without it decodes `tables/protective.tsv` differently on a cp1252
        # host, and raises on one whose codec has no mapping for the byte.
        text = files(__package__).joinpath(f'tables/{name}').read_text(encoding='utf-8')
    except (FileNotFoundError, ModuleNotFoundError) as e:   # pragma: no cover - a packaging failure
        raise FileNotFoundError(
            f'tables/{name} is missing from {__package__}.  A table read at runtime must be named in '
            '[tool.setuptools.package-data]; an undeclared one yields a wheel that imports fine and '
            'fails here, on the first molecule.') from e
    header: list[str] | None = None
    rows: list[dict[str, str]] = []
    for number, line in enumerate(text.split('\n'), 1):
        if not line or line.startswith('#'):
            continue
        cells = line.split('\t')
        if header is None:
            header = cells
            continue
        if len(cells) != len(header):
            raise ValueError(f'{name}:{number} has {len(cells)} cells, not {len(header)}')
        rows.append(dict(zip(header, cells)))
    if header is None:
        raise ValueError(f'{name} has no header row')
    return rows


def _decode_atom_fix(cell: str, row_id: str) -> tuple[tuple[int, int, bool | None], ...]:
    if cell == '-':
        return ()
    out = []
    for entry in cell.split(';'):
        parts = entry.split(':')
        if len(parts) != 3:
            raise ValueError(f'{row_id}: atom_fix entry {entry!r} is not slot:delta:radical')
        number, delta, radical = parts
        if radical not in _RADICAL:
            raise ValueError(f'{row_id}: radical {radical!r} is not one of - 0 1')
        out.append((int(number), int(delta), _RADICAL[radical]))
    return tuple(out)


def _decode_bonds_fix(cell: str, row_id: str) -> tuple[tuple[int, int, int], ...]:
    if cell == '-':
        return ()
    out = []
    for entry in cell.split(';'):
        parts = entry.split(':')
        if len(parts) != 3:
            raise ValueError(f'{row_id}: bonds_fix entry {entry!r} is not a:b:order')
        a, b, order = (int(p) for p in parts)
        if order not in (1, 2, 3, 8):
            # order 4 is deliberately absent: no repair writes an aromatic bond.
            raise ValueError(f'{row_id}: bond order {order} is not one of 1 2 3 8')
        out.append((a, b, order))
    return tuple(out)


def _decode_after(cell: str, row_id: str, index: dict[str, int]) -> tuple[str, ...]:
    """The `after` column, checked against the rows already loaded from this table.

    `index` holds only the rows read so far, so an unknown id is either a typo or a forward
    reference -- and rules run in file order, so a forward reference cannot be satisfied.
    """
    if cell == '-':
        return ()
    out = []
    for entry in cell.split(';'):
        if entry not in index:
            raise ValueError(f'{row_id}: after names {entry!r}, which is not an earlier row of this '
                             'table.  Rules run in file order, so an obligation to follow a later '
                             'row cannot be met')
        out.append(entry)
    return tuple(out)


def _anchors(wildcards: dict[int, str], written: frozenset[int]) -> frozenset[int]:
    """The atoms two matches of one rule may share.

    The asymmetry between the two wildcard kinds is deliberate.  An `[A]` the patch *writes* is part
    of the site being repaired, so a second match sharing it would apply the same charge delta
    twice.  An `[M]` is shared even when written: that is ferrocene, where each ring must add `+1`
    to the same iron to reach `[Fe+2]` with two `[Cp-]`.
    """
    return frozenset(number for number, kind in wildcards.items()
                     if kind == 'metal' or number not in written)


def _compile(name: str) -> tuple[Rule, ...]:
    rules = []
    index: dict[str, int] = {}
    for row in read_table(name):
        row_id = row['id']
        query, numbers, wildcards = compile_smarts(row['smarts'])
        atom_fix = _decode_atom_fix(row['atom_fix'], row_id)
        bonds_fix = _decode_bonds_fix(row['bonds_fix'], row_id)
        # a patch addressing an undeclared atom is a typo; refuse at load time rather than as a
        # KeyError on the first molecule that happens to match.
        for number, _, _ in atom_fix:
            if number not in numbers:
                raise ValueError(f'{row_id}: atom_fix names atom {number}, which '
                                 f'{row["smarts"]!r} does not declare')
        for a, b, _ in bonds_fix:
            for number in (a, b):
                if number not in numbers:
                    raise ValueError(f'{row_id}: bonds_fix names atom {number}, which '
                                     f'{row["smarts"]!r} does not declare')
        rules.append(Rule(row_id, query, numbers, atom_fix, bonds_fix,
                          row['tautomer'] == '1',
                          _decode_after(row['after'], row_id, index),
                          () if row['examples'] == '-' else tuple(row['examples'].split(';')),
                          row['why'], row['smarts'],
                          _anchors(wildcards, frozenset(n for n, _, _ in atom_fix))))
        index[row_id] = len(index)
    return tuple(rules)


def groups_rules() -> tuple[Rule, ...]:
    """The 82 functional-group repairs, in file order."""
    if 'groups' not in _RULES_CACHE:
        _RULES_CACHE['groups'] = _compile('standardize_groups.tsv')
    return _RULES_CACHE['groups']


def metals_rules() -> tuple[Rule, ...]:
    """The 19 metal-organic repairs, in file order."""
    if 'metals' not in _RULES_CACHE:
        _RULES_CACHE['metals'] = _compile('standardize_metals.tsv')
    return _RULES_CACHE['metals']


def standardize_rules() -> tuple[Rule, ...]:
    """Both tables, groups first -- the order `standardize()` applies them in."""
    return groups_rules() + metals_rules()


# --- resonance.tsv ------------------------------------------------------------------------------ #

#: The closed vocabulary.  A row naming anything else is a load-time error, because a role the pass
#: never asks for is a row that silently does nothing.
RESONANCE_ROLES = ('path', 'radical', 'donor', 'acceptor', 'veto_donor', 'veto_acceptor')


class Endpoint(NamedTuple):
    """One row of `tables/resonance.tsv`, compiled.

    `anchor` is the query stable id of the atom mapped `:1` -- the endpoint itself.  A row may
    describe a neighbourhood (`veto_azide` names three atoms) and only the anchor is the endpoint,
    so the pass reads one atom out of every embedding rather than the whole match.
    """
    id: str
    role: str
    query: QueryContainer
    anchor: int
    probe: str
    comment: str


_RESONANCE_CACHE: dict[str, object] = {}


def _compile_resonance() -> tuple[Endpoint, ...]:
    out = []
    seen = set()
    for row in read_table('resonance.tsv'):
        row_id = row['id']
        if row_id in seen:
            raise ValueError(f'resonance.tsv: {row_id} appears twice; an id is a log record\'s '
                             'only handle on a pattern and must name one row')
        seen.add(row_id)
        if row['role'] not in RESONANCE_ROLES:
            raise ValueError(f'{row_id}: role {row["role"]!r} is not one of '
                             f'{", ".join(RESONANCE_ROLES)}')
        query = read_smarts(row['smarts'])
        anchors = [n for n, number in query.map_numbers().items() if number == 1]
        if len(anchors) != 1:
            raise ValueError(f'{row_id}: {row["smarts"]!r} maps {len(anchors)} atoms to `:1`, not '
                             'one.  The endpoint is the `:1` atom, so a row without exactly one has '
                             'no endpoint for the pass to read out of a match')
        out.append(Endpoint(row_id, row['role'], query, anchors[0], row['probe'], row['comment']))
    return tuple(out)


def resonance_rules() -> tuple[Endpoint, ...]:
    """Every row of `tables/resonance.tsv`, in file order."""
    if 'rows' not in _RESONANCE_CACHE:
        _RESONANCE_CACHE['rows'] = _compile_resonance()
    return _RESONANCE_CACHE['rows']                                           # type: ignore[return-value]


def resonance_rules_by_role() -> dict[str, tuple[Endpoint, ...]]:
    """The same rows grouped by role, every role present even when it has no rows.

    Grouped here rather than in the pass so that a role with no rows is a `()` the pass iterates
    over rather than a `KeyError` on the first molecule.
    """
    if 'by_role' not in _RESONANCE_CACHE:
        grouped: dict[str, list] = {role: [] for role in RESONANCE_ROLES}
        for endpoint in resonance_rules():
            grouped[endpoint.role].append(endpoint)
        _RESONANCE_CACHE['by_role'] = {role: tuple(rows) for role, rows in grouped.items()}
    return _RESONANCE_CACHE['by_role']                                        # type: ignore[return-value]


def resonance_table_text() -> str:
    """`tables/resonance.tsv` verbatim, for the gate that reads the table without compiling it."""
    return files(__package__).joinpath('tables/resonance.tsv').read_text(encoding='utf-8')


# --- salts.tsv ---------------------------------------------------------------------------------- #
#
# This table has two matching mechanisms: a `cation`/`acceptor` row is a SMARTS matched by embedding,
# a `counterion`/`base`/`solvate` row names a compound matched by string equality of its key.  `role`
# says which.

#: The closed vocabulary.  A row naming anything else is a load-time error, because a role the pass
#: never asks for is a row that silently does nothing.
SALT_ROLES = ('cation', 'acceptor', 'counterion', 'base', 'solvate')

#: Which roles are SMARTS matched by embedding.  The rest are SMILES matched by canonical equality.
_SMARTS_ROLES = frozenset({'cation', 'acceptor'})


class SaltRow(NamedTuple):
    """One row of `tables/salts.tsv`, compiled.

    Exactly one of `query` (a `cation` or `acceptor` row) and `key` (the other three roles) is set.
    `anchor` is the stable id of the query atom mapped `:1` -- the subject, since a row may name a
    whole neighbourhood.  `key` is `format(species, '!s')`, the stereo-free canonical SMILES a species
    row is matched by, `None` for a SMARTS row.  `charges` is the set of charges a cation may end up
    with, `cation` rows only.
    """
    id: str
    role: str
    pattern: str
    query: QueryContainer | None
    key: str | None
    anchor: int
    charges: frozenset[int]
    comment: str


_SALTS_CACHE: dict[str, object] = {}


def _compile_salts() -> tuple[SaltRow, ...]:
    # deferred only to keep the module import cheap; still `chython.core`, not a layering exception.
    from ..core import read_smiles

    out = []
    seen = set()
    for row in read_table('salts.tsv'):
        row_id = row['id']
        if row_id in seen:
            raise ValueError(f'salts.tsv: {row_id} appears twice; an id is a log record\'s only '
                             'handle on a row and must name one')
        seen.add(row_id)
        role = row['role']
        if role not in SALT_ROLES:
            raise ValueError(f'{row_id}: role {role!r} is not one of {", ".join(SALT_ROLES)}')

        query = None
        anchor = 0
        key = None
        if role in _SMARTS_ROLES:
            query = read_smarts(row['pattern'])
            anchors = [n for n, number in query.map_numbers().items() if number == 1]
            if len(anchors) != 1:
                raise ValueError(f'{row_id}: {row["pattern"]!r} maps {len(anchors)} atoms to `:1`, '
                                 'not one.  The subject of a SMARTS row is its `:1` atom, so a row '
                                 'without exactly one has no subject for the pass to read out')
            anchor = anchors[0]
        else:
            species = read_smiles(row['pattern'])
            species.thiele()
            key = format(species, '!s')

        if row['charges'] == '-':
            charges = frozenset()
            if role == 'cation':
                raise ValueError(f'{row_id}: a cation row must list the charges it may end up with; '
                                 '`-` would make the overcharge guard vacuous and split a metal '
                                 'carbonyl')
        else:
            if role != 'cation':
                raise ValueError(f'{row_id}: charges are the cation overcharge guard and mean nothing '
                                 f'for a {role} row; write `-`')
            charges = frozenset(int(c) for c in row['charges'].split(';'))

        out.append(SaltRow(row_id, role, row['pattern'], query, key, anchor, charges, row['comment']))
    return tuple(out)


def salts_rows() -> tuple[SaltRow, ...]:
    """Every row of `tables/salts.tsv`, in file order."""
    if 'rows' not in _SALTS_CACHE:
        _SALTS_CACHE['rows'] = _compile_salts()
    return _SALTS_CACHE['rows']                                               # type: ignore[return-value]


def salts_rows_by_role() -> dict[str, tuple[SaltRow, ...]]:
    """The same rows grouped by role, every role present even when it has no rows.

    Grouped here rather than in the pass so a role with no rows is a `()` the pass iterates over
    rather than a `KeyError` on the first molecule.
    """
    if 'by_role' not in _SALTS_CACHE:
        grouped: dict[str, list] = {role: [] for role in SALT_ROLES}
        for row in salts_rows():
            grouped[row.role].append(row)
        _SALTS_CACHE['by_role'] = {role: tuple(rows) for role, rows in grouped.items()}
    return _SALTS_CACHE['by_role']                                            # type: ignore[return-value]


def salts_species_keys() -> dict[str, SaltRow]:
    """The `counterion`, `base` and `solvate` rows keyed by `format(species, '!s')`.

    One dict lookup per component is the whole match: the key is a stereo-free canonical SMILES, so a
    component either IS a tabulated species or is not, and there is no candidate list to walk.
    """
    if 'by_key' not in _SALTS_CACHE:
        index: dict[str, SaltRow] = {}
        for row in salts_rows():
            if row.key is not None:
                if row.key in index:
                    raise ValueError(f'salts.tsv: {row.id} and {index[row.key].id} are the same '
                                     f'species ({row.key}); two rows for one compound make which id a '
                                     'record reports depend on file order')
                index[row.key] = row
        _SALTS_CACHE['by_key'] = index
    return _SALTS_CACHE['by_key']                                             # type: ignore[return-value]


def salts_table_text() -> str:
    """`tables/salts.tsv` verbatim, for the gates that read it without compiling it."""
    return files(__package__).joinpath('tables/salts.tsv').read_text(encoding='utf-8')


# --- acids.tsv ---------------------------------------------------------------------------------- #

#: The closed vocabulary.  A row naming anything else is a load-time error, because a role the pass
#: never asks for is a row that silently does nothing.
ACID_ROLES = ('acid', 'base')


class AcidRow(NamedTuple):
    """One row of `tables/acids.tsv`, compiled.

    `anchor` is the query stable id of the atom mapped `:1` -- the site the proton comes off or goes
    onto.  A row may describe a neighbourhood (`acids:nitrate` names four atoms) and only the anchor
    is the site, so the pass reads one atom out of every embedding.
    """
    id: str
    role: str
    query: QueryContainer
    anchor: int
    smarts: str
    probe: str
    comment: str


_ACIDS_CACHE: dict[str, object] = {}


def _compile_acids() -> tuple[AcidRow, ...]:
    out = []
    seen = set()
    for row in read_table('acids.tsv'):
        row_id = row['id']
        if row_id in seen:
            raise ValueError(f'acids.tsv: {row_id} appears twice; an id is a log record\'s only '
                             'handle on a pattern and must name one row')
        seen.add(row_id)
        if row['role'] not in ACID_ROLES:
            raise ValueError(f'{row_id}: role {row["role"]!r} is not one of {", ".join(ACID_ROLES)}')
        query = read_smarts(row['smarts'])
        anchors = [n for n, number in query.map_numbers().items() if number == 1]
        if len(anchors) != 1:
            raise ValueError(f'{row_id}: {row["smarts"]!r} maps {len(anchors)} atoms to `:1`, not '
                             'one.  The site is the `:1` atom, so a row without exactly one has no '
                             'site for the pass to move a proton off or onto')
        out.append(AcidRow(row_id, row['role'], query, anchors[0], row['smarts'], row['probe'],
                           row['comment']))
    return tuple(out)


def acids_rules() -> tuple[AcidRow, ...]:
    """Every row of `tables/acids.tsv`, in file order."""
    if 'rows' not in _ACIDS_CACHE:
        _ACIDS_CACHE['rows'] = _compile_acids()
    return _ACIDS_CACHE['rows']                                               # type: ignore[return-value]


def acids_rules_by_role() -> dict[str, tuple[AcidRow, ...]]:
    """The same rows grouped by role, every role present even when it has no rows.

    Grouped here rather than in the pass so that a role with no rows is a `()` the pass iterates
    over rather than a `KeyError` on the first molecule.
    """
    if 'by_role' not in _ACIDS_CACHE:
        grouped: dict[str, list] = {role: [] for role in ACID_ROLES}
        for row in acids_rules():
            grouped[row.role].append(row)
        _ACIDS_CACHE['by_role'] = {role: tuple(rows) for role, rows in grouped.items()}
    return _ACIDS_CACHE['by_role']                                            # type: ignore[return-value]


def acids_table_text() -> str:
    """`tables/acids.tsv` verbatim, for the gate that reads the table without compiling it."""
    return files(__package__).joinpath('tables/acids.tsv').read_text(encoding='utf-8')


# --- sybyl_types.tsv ---------------------------------------------------------------------------- #

class SybylType(NamedTuple):
    """One row of `tables/sybyl_types.tsv`, compiled.

    `element` is the element symbol, or `''` for a pseudo-atom (lone pair, dummy, wildcard alias)
    that must not become an atom in the graph.  `hybridization` is the V3 code: 0=not determined,
    1=sp3, 2=sp2, 3=sp, 4=aromatic, 5=cumulated.
    """
    element: str
    hybridization: int


_SYBYL_CACHE: dict[str, object] = {}


def _compile_sybyl_types() -> dict[str, SybylType]:
    out: dict[str, SybylType] = {}
    for row in read_table('sybyl_types.tsv'):
        out[row['type']] = SybylType(row['element'], int(row['hybridization']))
    return out


def sybyl_types() -> dict[str, SybylType]:
    """SYBYL type string → :class:`SybylType`, loaded lazily on first use.

    Keys are the type strings as they appear in MOL2 ATOM blocks (``'C.3'``, ``'N.ar'``,
    ``'O.co2'``).  Types absent here -- bare element symbols like ``'Br'`` -- are resolved by the
    MOL2 reader's fallback.
    """
    if 'rows' not in _SYBYL_CACHE:
        _SYBYL_CACHE['rows'] = _compile_sybyl_types()
    return _SYBYL_CACHE['rows']                                                # type: ignore[return-value]


# --- rotatable.tsv ------------------------------------------------------------------------------ #
#
# The `rotatable` row is symmetric in `:1` and `:2`, so the matcher offers each bond in both
# directions and the pass deduplicates with a set.

#: The closed vocabulary.  A row naming anything else is a load-time error, because a role the pass
#: never asks for is a row that silently does nothing.
ROTATABLE_ROLES = ('rotatable', 'exclude')


class RotatableRow(NamedTuple):
    """One row of `tables/rotatable.tsv`.

    `numbers` maps the pattern's atom number to the query's stable id, so the pass can ask which two
    atoms `:1` and `:2` matched.
    """
    id: str
    role: str
    pattern: str
    query: QueryContainer
    numbers: dict[int, int]
    description: str


_ROTATABLE_CACHE: dict[str, object] = {}


def rotatable_rules() -> tuple[RotatableRow, ...]:
    """Rows of `tables/rotatable.tsv`, compiled, in file order.  Loaded on first call."""
    if 'rows' not in _ROTATABLE_CACHE:
        rows = []
        for n, row in enumerate(read_table('rotatable.tsv'), 1):
            role = row['role']
            if role not in ROTATABLE_ROLES:
                raise ValueError(f'rotatable.tsv line {n}: role {role!r} is not one of '
                                 f'{ROTATABLE_ROLES}')
            pattern = row['pattern']
            query, numbers, _ = compile_smarts(pattern)
            # check `query.map_numbers()`, not `numbers`: `compile_smarts` auto-numbers every
            # unnumbered atom, so `1 in numbers` would be an unfailable check.
            if not {1, 2} <= set(query.map_numbers().values()):
                raise ValueError(f'rotatable.tsv line {n}: the pattern must map :1 and :2 onto the '
                                 f'two atoms whose bond is the subject; {pattern!r} does not')
            rows.append(RotatableRow(f'rotatable:{row["id"]}', role, pattern, query, numbers,
                                     row['description']))
        _ROTATABLE_CACHE['rows'] = tuple(rows)
    return _ROTATABLE_CACHE['rows']                                            # type: ignore[return-value]


def rotatable_rules_by_role() -> dict[str, tuple[RotatableRow, ...]]:
    """`rotatable_rules()` grouped by role, file order preserved inside each group."""
    if 'by_role' not in _ROTATABLE_CACHE:
        out: dict[str, list] = {r: [] for r in ROTATABLE_ROLES}
        for row in rotatable_rules():
            out[row.role].append(row)
        _ROTATABLE_CACHE['by_role'] = {k: tuple(v) for k, v in out.items()}
    return _ROTATABLE_CACHE['by_role']                                         # type: ignore[return-value]


# --- hbond.tsv ---------------------------------------------------------------------------------- #
#
# The subject atom is always `:1`.  An amide nitrogen is z1, the same as an amine nitrogen, so the
# amine acceptor rows (9, 10, 11) exclude amides and sulfonamides by naming every heavy neighbour and
# demanding an sp3-or-aromatic carbon.  Do not collapse them into one `[N;D1,D2,D3;z1:1]`.

#: The closed vocabulary.  A row naming anything else is a load-time error, because a role the pass
#: never asks for is a row that silently does nothing.
HBOND_ROLES = ('donor', 'acceptor')


class HBondRow(NamedTuple):
    """One row of `tables/hbond.tsv`, compiled.

    `numbers` maps the pattern's atom number to the query's stable id, so the pass can ask which
    atom `:1` matched without walking the query.
    """
    id: str
    role: str
    pattern: str
    query: QueryContainer
    numbers: dict[int, int]
    description: str


_HBOND_CACHE: dict[str, object] = {}


def hbond_rules() -> tuple[HBondRow, ...]:
    """Rows of `tables/hbond.tsv`, compiled, in file order.  Loaded on first call."""
    if 'rows' not in _HBOND_CACHE:
        rows = []
        for n, row in enumerate(read_table('hbond.tsv'), 1):
            role = row['role']
            if role not in HBOND_ROLES:
                raise ValueError(f'hbond.tsv line {n}: role {role!r} is not one of {HBOND_ROLES}')
            pattern = row['pattern']
            query, numbers, _ = compile_smarts(pattern)
            # check `query.map_numbers()`, not `numbers`: `compile_smarts` auto-numbers every
            # unnumbered atom, so `1 in numbers` would be an unfailable check.
            if 1 not in set(query.map_numbers().values()):
                raise ValueError(f'hbond.tsv line {n}: the pattern must map :1 onto the subject '
                                 f'atom; {pattern!r} does not')
            rows.append(HBondRow(f'hbond:{row["id"]}', role, pattern, query, numbers,
                                 row['description']))
        _HBOND_CACHE['rows'] = tuple(rows)
    return _HBOND_CACHE['rows']                                                 # type: ignore[return-value]


def hbond_rules_by_role() -> dict[str, tuple[HBondRow, ...]]:
    """`hbond_rules()` grouped by role, file order preserved inside each group."""
    if 'by_role' not in _HBOND_CACHE:
        out: dict[str, list] = {r: [] for r in HBOND_ROLES}
        for row in hbond_rules():
            out[row.role].append(row)
        _HBOND_CACHE['by_role'] = {k: tuple(v) for k, v in out.items()}
    return _HBOND_CACHE['by_role']                                              # type: ignore[return-value]


# --- pharmacophore.tsv -------------------------------------------------------------------------- #
#
# The subject atom is always `:1`.  `donor` and `acceptor` are deliberately absent: they are
# `hbond.tsv`'s answer, and a second SMARTS spelling would be a definition no test compares.  The
# load-time disjointness assertion below keeps them out.

#: The closed vocabulary.  A row naming anything else is a load-time error, because a role the pass
#: never asks for is a row that silently does nothing.
PHARMACOPHORE_ROLES = ('positive', 'negative', 'aromatic', 'hydrophobe')


class PharmacophoreRow(NamedTuple):
    """One row of `tables/pharmacophore.tsv`, compiled.

    `numbers` maps the pattern's atom number to the query's stable id, so the pass can ask which
    atom `:1` matched without walking the query.
    """
    id: str
    role: str
    pattern: str
    query: QueryContainer
    numbers: dict[int, int]
    description: str


_PHARMACOPHORE_CACHE: dict[str, object] = {}


def pharmacophore_rules() -> tuple[PharmacophoreRow, ...]:
    """Rows of `tables/pharmacophore.tsv`, compiled, in file order.  Loaded on first call."""
    if 'rows' not in _PHARMACOPHORE_CACHE:
        assert not set(PHARMACOPHORE_ROLES) & set(HBOND_ROLES), \
            "donor and acceptor are tables/hbond.tsv's; pharmacophore.tsv must not re-spell them"
        rows = []
        for n, row in enumerate(read_table('pharmacophore.tsv'), 1):
            role = row['role']
            if role not in PHARMACOPHORE_ROLES:
                raise ValueError(f'pharmacophore.tsv line {n}: role {role!r} is not one of '
                                 f'{PHARMACOPHORE_ROLES}')
            pattern = row['pattern']
            query, numbers, _ = compile_smarts(pattern)
            # check `query.map_numbers()`, not `numbers`: `compile_smarts` auto-numbers every
            # unnumbered atom, so `1 in numbers` would be an unfailable check.
            if 1 not in set(query.map_numbers().values()):
                raise ValueError(f'pharmacophore.tsv line {n}: the pattern must map :1 onto the '
                                 f'subject atom; {pattern!r} does not')
            rows.append(PharmacophoreRow(f'pharmacophore:{row["id"]}', role, pattern, query,
                                         numbers, row['description']))
        _PHARMACOPHORE_CACHE['rows'] = tuple(rows)
    return _PHARMACOPHORE_CACHE['rows']                                          # type: ignore[return-value]


def pharmacophore_rules_by_role() -> dict[str, tuple[PharmacophoreRow, ...]]:
    """`pharmacophore_rules()` grouped by role, file order preserved inside each group."""
    if 'by_role' not in _PHARMACOPHORE_CACHE:
        out: dict[str, list] = {r: [] for r in PHARMACOPHORE_ROLES}
        for row in pharmacophore_rules():
            out[row.role].append(row)
        _PHARMACOPHORE_CACHE['by_role'] = {k: tuple(v) for k, v in out.items()}
    return _PHARMACOPHORE_CACHE['by_role']                                       # type: ignore[return-value]


# --- tpsa.tsv ----------------------------------------------------------------------------------- #
#
# Ertl, Rohde, Selzer, J. Med. Chem. 2000, 43, 3714, Table 1.  Two classes: NO (the published TPSA)
# and SP (the optional sulfur/phosphorus extension).  First match wins in file order -- an epoxide
# oxygen also matches the ether row, so the epoxide row must precede it.  The subject atom is `:1`.

#: The closed vocabulary.  A row naming anything else is a load-time error, because a class the pass
#: never sums is a row that silently contributes nothing.
TPSA_CLASSES = ('NO', 'SP')


class TpsaRow(NamedTuple):
    """One row of `tables/tpsa.tsv`.  Ertl, Rohde, Selzer, J. Med. Chem. 2000, 43, 3714."""
    id: str
    element_class: str
    contribution: float
    pattern: str
    query: QueryContainer
    numbers: dict[int, int]
    description: str


_TPSA_CACHE: dict[str, object] = {}


def tpsa_rules() -> tuple[TpsaRow, ...]:
    """Rows of `tables/tpsa.tsv`, compiled, IN FILE ORDER -- which is MATCH ORDER.

    The patterns overlap and the first match wins, so re-sorting this table changes the descriptor.
    """
    if 'rows' not in _TPSA_CACHE:
        rows = []
        for n, row in enumerate(read_table('tpsa.tsv'), 1):
            cls = row['element_class']
            if cls not in TPSA_CLASSES:
                raise ValueError(f'tpsa.tsv line {n}: element_class {cls!r} is not one of '
                                 f'{TPSA_CLASSES}')
            pattern = row['pattern']
            query, numbers, _ = compile_smarts(pattern)
            # a pattern with no `:1` types nothing; refusing it here spares `first_match` the check.
            if 1 not in set(query.map_numbers().values()):
                raise ValueError(f'tpsa.tsv line {n}: the pattern maps no atom to 1, so it types '
                                 f'nothing')
            rows.append(TpsaRow(f"tpsa:{row['id']}", cls, float(row['contribution']), pattern,
                                query, numbers, row['description']))
        _TPSA_CACHE['rows'] = tuple(rows)
    return _TPSA_CACHE['rows']                                                   # type: ignore[return-value]


# --- the shared first-match resolver ------------------------------------------------------------ #

def first_match(rows: Sequence[Any], molecule: object) -> dict[int, Any]:
    """Resolve an overlapping, order-dependent atom-typing table against a molecule.

    Returns `{n: row}` keeping, per atom, the earliest row in `rows` that matches it.  Shared
    by `tables/tpsa.tsv` and `tables/crippen.tsv`; every row must carry `numbers` and have mapped its
    subject to 1, which both loaders check.  `rows` is `Sequence[Any]` because the two callers pass
    different NamedTuples and only `.numbers` and `.query` are read.
    """
    out = {}
    for row in rows:
        subject = row.numbers[1]
        for mapping in row.query.get_mapping(molecule):
            out.setdefault(mapping[subject], row)
    return out


# --- crippen.tsv -------------------------------------------------------------------------------- #
#
# Wildman, Crippen, J. Chem. Inf. Comput. Sci. 1999, 39, 868, Table 1.  A logP and a molar
# refractivity contribution over one 72-type inventory, resolved first-match-wins by `first_match`.

CRIPPEN_ROLES = ('heavy', 'hydrogen')

#: The only rows allowed to carry no probe: nothing reaches a catch-all until every row above it has
#: failed, so "it never fires" is a well-covered block rather than a defect.
CRIPPEN_CATCH_ALLS = ('CS', 'HS', 'NS', 'OS')

# FILE ORDER, WHICH IS MATCH ORDER -- not the paper's numbering; seven rows move, each a special case
# of the row that would otherwise swallow it (see crippen.tsv's header for the measurement).  `S2`
# appears twice on purpose: its two published alternatives are a charge test and an S=X bond test,
# which no single chython pattern expresses, so both rows carry type S2 and the same id.  Hence 73
# entries for 72 types, and `crippen_rules()` checks this as a sequence rather than a set.
CRIPPEN_TYPES = ('C8', 'C2', 'C1', 'C3', 'C4', 'C5', 'C26', 'C6', 'C7', 'C9', 'C10', 'C11', 'C12',
                 'C13', 'C14', 'C15', 'C16', 'C17', 'C18', 'C19', 'C20', 'C21', 'C22', 'C23', 'C24',
                 'C25', 'C27', 'CS', 'H1', 'H4', 'H2', 'H3', 'HS', 'N1', 'N2', 'N3', 'N4', 'N5',
                 'N6', 'N7', 'N8', 'N9', 'N10', 'N11', 'N12', 'N13', 'N14', 'NS', 'O1', 'O2', 'O3',
                 'O4', 'O5', 'O6', 'O12', 'O7', 'O8', 'O10', 'O11', 'O9', 'OS', 'F', 'Cl', 'Br',
                 'I', 'Hal', 'P', 'S2', 'S2', 'S1', 'S3', 'Me1', 'Me2')

_CRIPPEN_CACHE: dict[str, object] = {}


class CrippenRow(NamedTuple):
    """One row of `tables/crippen.tsv`.  Wildman, Crippen, J. Chem. Inf. Comput. Sci. 1999, 39, 868."""
    id: str
    type: str
    role: str
    logp: float
    mr: float
    mr_published: bool
    pattern: str
    query: QueryContainer
    numbers: dict[int, int]
    probe: str
    description: str


def crippen_rules() -> tuple[CrippenRow, ...]:
    """Rows of `tables/crippen.tsv`, compiled, IN FILE ORDER -- which is MATCH ORDER.

    The paper's types overlap deliberately, so file order is the disambiguation between them and
    re-sorting changes both descriptors.  The inventory check turns a dropped or mis-transcribed row
    into an error at first use rather than a quietly missing contribution.
    """
    if 'rows' not in _CRIPPEN_CACHE:
        rows = []
        for n, row in enumerate(read_table('crippen.tsv'), 1):
            type_, role, mr = row['type'], row['role'], row['mr']
            if role not in CRIPPEN_ROLES:
                raise ValueError(f'crippen.tsv line {n}: role {role!r} is not one of {CRIPPEN_ROLES}')
            probe = row['probe']
            if (probe == '-') != (type_ in CRIPPEN_CATCH_ALLS):
                raise ValueError(f'crippen.tsv line {n}: only the catch-alls {CRIPPEN_CATCH_ALLS} may '
                                 f'omit a probe, and each of them must')
            # a blank mr cell is `-`, stored as 0.0 and flagged: O7 and O9 carry a real published mr
            # of exactly 0, so an absent value must stay distinguishable from them.
            published = mr != '-'
            pattern = row['pattern']
            query, numbers, _ = compile_smarts(pattern)
            if 1 not in set(query.map_numbers().values()):
                raise ValueError(f'crippen.tsv line {n}: the pattern maps no atom to 1, so it types '
                                 f'nothing')
            rows.append(CrippenRow(f'crippen:{type_}', type_, role, float(row['logp']),
                                   float(mr) if published else 0.0, published, pattern,
                                   query, numbers, probe, row['description']))
        if tuple(r.type for r in rows) != CRIPPEN_TYPES:
            raise ValueError('crippen.tsv: the type column is not the inventory in file order')
        _CRIPPEN_CACHE['rows'] = tuple(rows)
    return _CRIPPEN_CACHE['rows']                                                # type: ignore[return-value]


def crippen_rules_by_role() -> dict[str, tuple[CrippenRow, ...]]:
    """`crippen_rules()` split by role, file order preserved inside each.

    The two roles resolve separately: a hydrogen row's subject is the carrier atom, which a heavy row
    would claim first in an unscoped pass.
    """
    if 'by_role' not in _CRIPPEN_CACHE:
        out: dict[str, list] = {r: [] for r in CRIPPEN_ROLES}
        for row in crippen_rules():
            out[row.role].append(row)
        _CRIPPEN_CACHE['by_role'] = {k: tuple(v) for k, v in out.items()}
    return _CRIPPEN_CACHE['by_role']                                             # type: ignore[return-value]


# --- maccs.tsv ---------------------------------------------------------------------------------- #
#
# The 166 published keys are not 166 SMARTS: some are a count, two count RINGS, one asks about isotopes
# and one about the record's fragment count.  Hence four kinds and a named-predicate registry.

#: The four kinds a key row can be.  `unset` is a key with no stated definition, not a key that is hard.
MACCS_KINDS = ('smarts', 'count', 'predicate', 'unset')

#: The registry a `predicate` row names.  `_maccs.py`'s `MACCS_PREDICATE_FNS` asserts equality with this
#: tuple AT MODULE IMPORT, so a name added to one alone raises on `import chython.chemistry`.
MACCS_PREDICATES = ('isotope', 'atomic_number_gt_103', 'charge', 'fragments_gt_1', 'ring_present',
                    'aromatic_rings_gt_1', 'six_rings_gt_1')

#: The keys that ship permanently zero.  Key 44's published description is the literal placeholder
#: `OTHER`: there is nothing to transcribe, a guessed pattern would be invented chemistry, and a
#: documented zero is the honest answer.  Asserted in both directions, so a second entry is a decision
#: taken here and not at transcription time.
MACCS_UNSET_KEYS = (44,)

#: The two answers a corpus row can state.
MACCS_EXPECTATIONS = ('set', 'unset')

_MACCS_CACHE: dict[str, object] = {}
_MACCS_CORPUS_CACHE: dict[str, object] = {}


class MaccsRow(NamedTuple):
    """One published MACCS key.  Durant, Leland, Henry, Nourse, JCICS 2002, 42, 1273.

    `query` is None for a `predicate` and an `unset` row; `count` is the number of DISTINCT matched atom
    sets at which the bit is set, which is how the published `> n` wordings are expressed.
    """
    id: str
    key: int
    kind: str
    pattern: str
    query: QueryContainer | None
    count: int
    predicate: str
    description: str


class MaccsCorpusRow(NamedTuple):
    """One acceptance row of `tables/maccs_corpus.tsv`."""
    key: int
    expectation: str
    smiles: str
    name: str


def maccs_rules() -> tuple[MaccsRow, ...]:
    """The 166 published MACCS keys, compiled, ordered by key number.  Loaded on first call."""
    if 'rows' not in _MACCS_CACHE:
        rows = []
        for n, row in enumerate(read_table('maccs.tsv'), 1):
            key, kind = int(row['key']), row['kind']
            pattern, predicate = row['pattern'], row['predicate']
            if kind not in MACCS_KINDS:
                raise ValueError(f'maccs.tsv line {n}: kind {kind!r} is not one of {MACCS_KINDS}')
            if kind == 'predicate':
                if predicate not in MACCS_PREDICATES:
                    raise ValueError(f'maccs.tsv line {n}: predicate {predicate!r} is not one of '
                                     f'{MACCS_PREDICATES}')
                query = None
            elif kind == 'unset':
                if key not in MACCS_UNSET_KEYS:
                    raise ValueError(f'maccs.tsv line {n}: key {key} may not ship unset; only '
                                     f'{MACCS_UNSET_KEYS} have no published definition')
                if pattern != '-' or predicate != '-':
                    raise ValueError(f'maccs.tsv line {n}: an unset row states no pattern and no '
                                     f'predicate')
                query = None
            else:
                # a MACCS pattern types no subject atom -- the engine counts distinct matched atom
                # SETS -- so neither the numbering nor the wildcards is kept.
                query, _, _ = compile_smarts(pattern)
            rows.append(MaccsRow(f'maccs:{key}', key, kind, pattern, query, int(row['count']),
                                 predicate, row['description']))
        rows.sort(key=lambda r: r.key)
        _MACCS_CACHE['rows'] = tuple(rows)
    return _MACCS_CACHE['rows']                                                  # type: ignore[return-value]


def maccs_rules_by_key() -> dict[int, MaccsRow]:
    """`maccs_rules()` keyed by published key number."""
    if 'by_key' not in _MACCS_CACHE:
        _MACCS_CACHE['by_key'] = {r.key: r for r in maccs_rules()}
    return _MACCS_CACHE['by_key']                                                # type: ignore[return-value]


def maccs_corpus() -> tuple[MaccsCorpusRow, ...]:
    """Rows of `tables/maccs_corpus.tsv`, in file order.  Loaded on first call.

    The SMILES are NOT compiled here: `_tables.py` may not import a pass, and a data file should not
    become 330 containers before anything asks.  The test parses them.
    """
    if 'rows' not in _MACCS_CORPUS_CACHE:
        rows = []
        for n, row in enumerate(read_table('maccs_corpus.tsv'), 1):
            expectation = row['expectation']
            if expectation not in MACCS_EXPECTATIONS:
                raise ValueError(f'maccs_corpus.tsv line {n}: expectation {expectation!r} is not '
                                 f'one of {MACCS_EXPECTATIONS}')
            key = int(row['key'])
            if key in MACCS_UNSET_KEYS:
                raise ValueError(f'maccs_corpus.tsv line {n}: key {key} has no published definition '
                                 f'and ships permanently unset; it cannot be exemplified')
            rows.append(MaccsCorpusRow(key, expectation, row['smiles'], row['name']))
        _MACCS_CORPUS_CACHE['rows'] = tuple(rows)
    return _MACCS_CORPUS_CACHE['rows']                                           # type: ignore[return-value]


def maccs_corpus_by_key() -> dict[int, tuple[MaccsCorpusRow, ...]]:
    """`maccs_corpus()` grouped by key number."""
    if 'by_key' not in _MACCS_CORPUS_CACHE:
        out: dict[int, list] = {}
        for row in maccs_corpus():
            out.setdefault(row.key, []).append(row)
        _MACCS_CORPUS_CACHE['by_key'] = {k: tuple(v) for k, v in out.items()}
    return _MACCS_CORPUS_CACHE['by_key']                                         # type: ignore[return-value]


# --- qed_alerts.tsv ----------------------------------------------------------------------------- #
#
# QED's eighth input is a count of structural alerts.  Every row carries a `probe`, a public compound the
# alert must match: an alert that matches nothing lowers no score and would sit here unnoticed.

_QED_ALERTS_CACHE: dict[str, object] = {}


class QedAlertRow(NamedTuple):
    """One structural alert QED counts.  Brenk et al., ChemMedChem 2008, 3, 435."""
    id: str
    name: str
    pattern: str
    query: QueryContainer
    probe: str
    description: str


def qed_alerts() -> tuple[QedAlertRow, ...]:
    """The QED structural alerts, compiled, in file order.  Loaded on first call."""
    if 'rows' not in _QED_ALERTS_CACHE:
        rows = []
        for n, row in enumerate(read_table('qed_alerts.tsv'), 1):
            pattern = row['pattern']
            # `compile_smarts` returns `(query, numbers, wildcards)`; an alert has no mapped atom and no
            # patch column, so only the query is kept -- unpacked rather than indexed, so that a
            # signature change fails here rather than later.
            query, _, _ = compile_smarts(pattern)
            rows.append(QedAlertRow(f'qed_alerts:{row["id"]}', row['name'], pattern, query,
                                    row['probe'], row['description']))
        _QED_ALERTS_CACHE['rows'] = tuple(rows)
    return _QED_ALERTS_CACHE['rows']                                             # type: ignore[return-value]


# --- abbreviations.tsv -------------------------------------------------------------------------- #
#
# The one table whose fragment is a MOLECULE and not a query: an abbreviation names a structure to graft,
# not a pattern to find, so it is lexed by `read_smiles` rather than by `read_smarts`.


class AbbreviationRow(NamedTuple):
    """One row of `tables/abbreviations.tsv`, compiled.

    `fragment` is the group as a molecule, with `marker` the stable id of the `*` standing for the bond
    to the rest of the structure and `attachment` the id of the atom that bond reaches.  Every hydrogen
    count in `fragment` is already the count of the ATTACHED group, which is what the marker buys.
    """
    id: str
    label: str
    smiles: str
    synonyms: tuple[str, ...]
    fragment: Any
    marker: int
    attachment: int


_ABBREVIATIONS_CACHE: dict[str, object] = {}


def abbreviations_rows() -> tuple[AbbreviationRow, ...]:
    """Rows of `tables/abbreviations.tsv`, compiled, in file order.  Loaded on first call."""
    if 'rows' not in _ABBREVIATIONS_CACHE:
        rows = []
        for n, row in enumerate(read_table('abbreviations.tsv'), 1):
            label = row['label']
            smiles = row['smiles']
            fragment = read_smiles(smiles)
            markers = [a.n for a in fragment.atoms() if a.is_r]
            if len(markers) != 1:
                raise ValueError(f'abbreviations.tsv line {n}: {label} has {len(markers)} `*` markers, '
                                 f'and a group attached anywhere but at one atom is a different fact')
            marker = markers[0]
            neighbors = list(fragment.neighbors_of(marker))
            if len(neighbors) != 1:
                raise ValueError(f'abbreviations.tsv line {n}: {label}\'s marker has '
                                 f'{len(neighbors)} bonds, so it names no attachment')
            if fragment.order_of(marker, neighbors[0]) != 1:
                raise ValueError(f'abbreviations.tsv line {n}: {label} attaches by a multiple bond; '
                                 f'the marker states the single bond a contracted group hangs by')
            synonyms = () if row['synonyms'] == '-' else tuple(row['synonyms'].split(','))
            rows.append(AbbreviationRow(f'abbreviations:{label}', label, smiles, synonyms,
                                        fragment, marker, neighbors[0]))
        _ABBREVIATIONS_CACHE['rows'] = tuple(rows)
    return _ABBREVIATIONS_CACHE['rows']                                          # type: ignore[return-value]


def abbreviations_by_label() -> dict[str, AbbreviationRow]:
    """Every spelling in `tables/abbreviations.tsv`, label and synonym alike, to its row.

    Case-folded keys are in the same dict, and a folded key that two rows would claim is a load-time
    error: `abbreviation_row` answers one row per spelling, so the table may not contain the collision.
    """
    if 'by_label' not in _ABBREVIATIONS_CACHE:
        exact: dict[str, AbbreviationRow] = {}
        folded: dict[str, AbbreviationRow] = {}
        for row in abbreviations_rows():
            for spelling in (row.label, *row.synonyms):
                if spelling in exact:
                    raise ValueError(f'abbreviations.tsv: {spelling!r} is claimed by both '
                                     f'{exact[spelling].label} and {row.label}')
                exact[spelling] = row
                key = spelling.casefold()
                if key in folded and folded[key] is not row:
                    raise ValueError(f'abbreviations.tsv: {spelling!r} folds onto a spelling of '
                                     f'{folded[key].label}, so a case-blind lookup has two answers')
                folded[key] = row
        _ABBREVIATIONS_CACHE['by_label'] = exact
        _ABBREVIATIONS_CACHE['folded'] = folded
    return _ABBREVIATIONS_CACHE['by_label']                                      # type: ignore[return-value]


def abbreviation_row(label: str) -> AbbreviationRow | None:
    """The row `label` names, or None.  Exact spelling first, then case-folded.

    Two passes rather than one folded lookup, and the order is the point: a file's own spelling wins
    where the table has it, so `Ts` cannot be answered by a row for `ts` that the table grows later.
    """
    exact = abbreviations_by_label()
    if label in exact:
        return exact[label]
    return _ABBREVIATIONS_CACHE['folded'].get(label.casefold())                  # type: ignore[union-attr]


# --- covalent_radii.tsv ------------------------------------------------------------------------- #
#
# One radius per element and no pattern, so there is no NamedTuple: a row is `z -> radius` and a
# caller wants the number.  The table's own header states why this is not `core`'s `atomic_radius`.

_COVALENT_RADII_CACHE: dict[str, dict[int, float]] = {}


def covalent_radii() -> dict[int, float]:
    """Atomic number -> single-bond covalent radius in Angstroms, for the elements the table covers.

    An element with no row is ABSENT rather than zero, so a caller that needs a radius asks and gets
    a `KeyError` or a miss instead of a threshold built from nothing.  `perceive_bonds` reads this.
    """
    if 'radii' not in _COVALENT_RADII_CACHE:
        radii: dict[int, float] = {}
        for row in read_table('covalent_radii.tsv'):
            z = int(row['z'])
            if z in radii:
                raise ValueError(f'covalent_radii.tsv: element {z} appears twice')
            radius = float(row['radius'])
            if radius <= 0.:
                raise ValueError(f'covalent_radii.tsv: element {z} states radius {radius}, which is '
                                 'not a length; an element the survey does not cover has no row')
            radii[z] = radius
        _COVALENT_RADII_CACHE['radii'] = radii
    return _COVALENT_RADII_CACHE['radii']
