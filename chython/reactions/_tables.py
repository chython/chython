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
"""Read the reaction corpus out of `tables/` and compose it into `ReactionTemplate`s.

`chemistry/_tables.py`'s conventions: TSV in `tables/`, a NamedTuple per row, a module-level cache, an
accessor, nothing loaded at import, rule ids qualified by their table (`reactions:13`).  A reaction row
names `functional.tsv` groups rather than repeating them.  This is the knowledge half of the package
and may not import the enumerators beside it.
"""
from importlib.resources import files
from typing import NamedTuple
from ..core import QueryContainer, ReactionTemplate, read_smarts, read_smirks


__all__ = ['FunctionalGroup', 'ProtectiveGroup', 'ReactionRule', 'Role', 'SLOT_STRIDE',
           'compose_smirks', 'functional_rules', 'protective_rules', 'read_table', 'reaction_rules',
           'roles']


#: Slot *i*'s map numbers are offset by `i * SLOT_STRIDE`, so a group's own numbers never move: adding
#: an atom to a functional group cannot invalidate the product side of a row that references it.
SLOT_STRIDE = 100


class FunctionalGroup(NamedTuple):
    """One row of `functional.tsv`: a named pattern, sealed.

    `query` answers the presence question and `smarts` is what the composer splices into a reactant
    side; both come from the one cell.

    `example` and `decoys` are the row's own acceptance test, run by `test/test_functional.py`: a
    pattern that can never match is otherwise invisible, since an absent group and an unmatchable one
    are the same missing key in `mol.functional_groups()`.
    """
    id: str
    name: str
    smarts: str
    query: QueryContainer
    example: str
    decoys: tuple[str, ...]
    description: str


class ReactionRule(NamedTuple):
    """One row of `reactions.tsv`, composed into one or two templates.

    `groups` is the reactant slots in order, so slot *i*'s map numbers are offset by `i * SLOT_STRIDE`.
    `template` is the intermolecular composition; `intramolecular` is the one the `ring_sizes` column
    asks for, or `None`.  Both carry this row's `id` as their `rule_id`, so a log record names the row
    rather than the composed string nobody wrote.

    `len(groups)` is not an arity: `react()` hands a row every input at once, so a two-slot row applies
    to two molecules, to one mixture, and -- through `intramolecular` -- to one molecule holding both.

    `probe` is the row's own acceptance test, `<reactants>>><product>`, run by `test/test_probes.py`:
    a row that composes cleanly and never fires is otherwise invisible.  EACH SIDE IS ONE RECORD --
    the reactant side finds its components wherever they are, so `CCBr.Oc1ccccc1` reads as one record
    and satisfies the intermolecular template, while a connected reactant side satisfies the
    intramolecular one.
    """
    id: str
    name: str
    groups: tuple[str, ...]
    product: str
    template: ReactionTemplate
    description: str
    probe: str = ''
    ring_atom: int = 0
    ring_sizes: tuple[int, ...] = ()
    intramolecular: ReactionTemplate | None = None

    @property
    def templates(self) -> tuple[ReactionTemplate, ...]:
        """Every template this row composed to, intermolecular first.  Order is outcome order."""
        if self.intramolecular is None:
            return (self.template,)
        return (self.template, self.intramolecular)


class ProtectiveGroup(NamedTuple):
    """One row of `protective.tsv`: a protecting group and the patch that removes it.

    `template` is composed from the row's own two SMARTS columns rather than from named slots -- a
    protecting group is single-use and does not belong in a table of shared patterns.  Otherwise it is
    a reaction row: one `read_smirks` call, a table-qualified `rule_id`, deletion by absence.

    `size` is the specificity, which is why it is a field the TSV does not carry.  Rows are served
    largest-first, so `hydroxyl_boc` (11 atoms) claims a Boc before `hydroxyl_tbu` (5) can take its
    tert-butyl half and leave a carbonate behind.  Computing it means a row can be added anywhere.
    """
    id: str
    name: str
    protects: tuple[str, ...]
    smarts: str
    product: str
    template: ReactionTemplate
    protected: str
    cleaved: str
    decoys: tuple[str, ...]
    description: str

    @property
    def size(self) -> int:
        """How many atoms the reactant side matches.  See the class docstring."""
        return len(self.template.reactants)


class Role(NamedTuple):
    """One row of `roles.tsv`: a coupling handle and the patch that cuts it and caps the cut.

    The cap is part of `product`, spelled `[#0:ROLE_CAP]`, so the attachment point exists inside the
    patch -- which is what lets the cut centre keep its configuration, and the site says `@=` to keep
    it: a cut takes no configuration away with the fragment that leaves, and a patch drops the
    configuration at its own reaction centre unless a template says otherwise.  `report=True` says
    which atom the marker landed on.
    """
    id: str
    name: str
    group: str
    product: str
    template: ReactionTemplate
    example: str
    decoys: tuple[str, ...]
    description: str


#: The map number every `roles.tsv` row gives its cap, so the enumerator finds the marker by number
#: instead of scanning the product for an R -- an input that already carries one stays distinguishable.
ROLE_CAP = 20

_FUNCTIONAL_CACHE: dict[str, FunctionalGroup] = {}

#: The composed corpus grouped by reaction name, or `None` before anything asked for it.
_RULES_CACHE: dict[str, tuple[ReactionRule, ...]] | None = None

#: `protective.tsv` composed and sorted, keyed by name, or `None`.  Separate from `_RULES_CACHE` so a
#: process that only deprotects does not compile the reaction corpus.
_PROTECTIVE_CACHE: dict[str, ProtectiveGroup] | None = None

#: `roles.tsv` composed and grouped by role name, or `None` before anything asked for it.
_ROLES_CACHE: dict[str, tuple[Role, ...]] | None = None


def read_table(name: str) -> list[dict[str, str]]:
    """Parse one TSV out of `tables/` into a list of row dicts, without compiling anything.

    A copy of `chemistry/_tables.py`'s function and not an import: that one resolves `tables/` against
    its own `__package__`, and `chemistry` is a sibling this package may not import.  The `tables/`
    prefix belongs inside the `joinpath` call -- a package-data pattern is matched against the path
    relative to the package, so a bare filename ships nothing.
    """
    try:
        # `encoding='utf-8'`: the tables are shipped bytes and the codec that reads them may not be the
        # user's locale -- without it a table holding a non-ASCII byte decodes differently on a cp1252
        # host, and raises on one whose codec has no mapping for that byte.
        text = files(__package__).joinpath(f'tables/{name}').read_text(encoding='utf-8')
    except (FileNotFoundError, ModuleNotFoundError) as e:   # pragma: no cover - a packaging failure
        raise FileNotFoundError(
            f'tables/{name} is missing from {__package__}.  A table read at runtime must be named in '
            '[tool.setuptools.package-data]; an undeclared one yields a wheel that imports fine and '
            'fails here, on the first enumeration.') from e
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


def _offset_map_numbers(smarts: str, offset: int) -> str:
    """Add `offset` to every `:N` map number of one SMARTS string.

    A string rewrite, because the composed side goes through `read_smirks` whole.  The scan tracks
    brackets rather than using a regex: `:` is also the aromatic bond token, and the only `:` a map
    number can follow is one inside a bracket.
    """
    if not offset:
        return smarts
    out = []
    depth = 0
    i = 0
    n = len(smarts)
    while i < n:
        c = smarts[i]
        if c == '[':
            depth += 1
        elif c == ']':
            depth -= 1
        elif c == ':' and depth > 0:
            j = i + 1
            while j < n and smarts[j].isdigit():
                j += 1
            if j > i + 1:                     # a map number, not an aromatic bond inside a bracket
                out.append(':%d' % (int(smarts[i + 1:j]) + offset))
                i = j
                continue
        out.append(c)
        i += 1
    return ''.join(out)


def _slots(groups: tuple[str, ...]) -> list[str]:
    """Each named group's SMARTS, offset by its slot."""
    known = functional_rules()
    parts = []
    for slot, name in enumerate(groups):
        group = known.get(name)
        if group is None:
            raise ValueError(f'unknown functional group {name!r}; functional.tsv names '
                             f'{len(known)} groups')
        parts.append(_offset_map_numbers(group.smarts, slot * SLOT_STRIDE))
    return parts


def compose_smirks(groups: tuple[str, ...], product: str, *, intramolecular: bool = False) -> str:
    """The SMIRKS one reaction row means: its groups, offset by slot, then `>>` and the product side.

    Public because printing it is how a row author sees the numbering their product side must address.

    `intramolecular=False` composes `(A).(B)`, demanding the groups be in DIFFERENT molecule
    components; `True` composes `(A.B)`, demanding one.  Never a bare `.`, which says only "not bonded"
    and so also matches the groups bonded together -- a cyclization with no ring size stated.  A
    one-slot row gets no grouping, there being nothing to constrain.
    """
    parts = _slots(groups)
    if len(parts) < 2:
        reactants = ''.join(parts)
    elif intramolecular:
        reactants = '(%s)' % '.'.join(parts)
    else:
        reactants = '.'.join('(%s)' % part for part in parts)
    return '%s>>%s' % (reactants, product)


def _parse_ring_sizes(cell: str, row_id: str) -> tuple[int, tuple[int, ...]]:
    """`'1:5,6,7'` -> `(1, (5, 6, 7))`, and `''` -> `(0, ())`.

    Refused at load rather than tolerated: every way of getting the cell wrong composes a template that
    reads and never fires, which is invisible.
    """
    if not cell:
        return 0, ()
    atom, _, sizes = cell.partition(':')
    if not sizes:
        raise ValueError(f'{row_id}: ring_sizes is {cell!r}; the column is `<atom>:<sizes>`, as in '
                         '`1:5,6,7` -- "the ring through product atom :1 is 5-, 6- or 7-membered"')
    try:
        parsed = (int(atom), tuple(int(size) for size in sizes.split(',')))
    except ValueError as e:
        raise ValueError(f'{row_id}: ring_sizes is {cell!r}, which is not `<atom>:<sizes>`') from e
    for size in parsed[1]:
        if not 3 <= size <= 14:
            raise ValueError(f'{row_id}: ring_sizes asks for a {size}-membered ring; the `r` primitive '
                             'spans 3-14, so this row would compose a template that never fires')
    return parsed


def _attach_ring_sizes(product: str, atom: int, sizes: tuple[int, ...]) -> str:
    """Add `;r<a>,r<b>,...` to the product-side atom numbered `atom`.

    The disjunction goes last inside the bracket, immediately before the map number, which is safe only
    because a map number is not part of the AND/OR grammar: a high AND after a `,` list would bind to
    its last alternative alone, but `:1` is an atom-level field, so `[N;D1;z1;x0;r5,r6,r7:101]`
    constrains the whole atom and not just its `r7` branch.  Pinned by `test/test_tables.py`.
    """
    ring = ';' + ','.join('r%d' % size for size in sizes)
    out = []
    depth = 0
    i = 0
    found = False
    while i < len(product):
        c = product[i]
        if c == '[':
            depth += 1
        elif c == ']':
            depth -= 1
        elif c == ':' and depth > 0:
            j = i + 1
            while j < len(product) and product[j].isdigit():
                j += 1
            if j > i + 1:
                if int(product[i + 1:j]) == atom:
                    out.append(ring)
                    found = True
                out.append(product[i:j])
                i = j
                continue
        out.append(c)
        i += 1
    if not found:
        raise ValueError(f'ring_sizes names product atom :{atom}, which {product!r} does not mention. '
                         'The ring constraint has to land on an atom the product side states, and on '
                         'one that lies on the new ring -- otherwise the template never fires.')
    return ''.join(out)


def functional_rules() -> dict[str, FunctionalGroup]:
    """Every row of `functional.tsv` by name, compiled on first use.

    Keyed by name because both consumers look a group up that way: the composer resolving a row's
    slots, and `mol.functional_groups()` reporting under the same names.  THE TABLE AND THE PASS ARE
    TWO QUESTIONS: this reads the rows, `functional_groups(molecule)` asks a molecule which of them
    it carries.
    """
    if not _FUNCTIONAL_CACHE:
        for row in read_table('functional.tsv'):
            name = row['name']
            if name in _FUNCTIONAL_CACHE:
                raise ValueError(f'functional.tsv: {name!r} is defined twice')
            row_id = 'functional:%s' % row['id']
            try:
                query = read_smarts(row['smarts'])
            except Exception as e:
                raise ValueError(f'{row_id}: {row["smarts"]!r} does not read as SMARTS: {e}') from e
            _FUNCTIONAL_CACHE[name] = FunctionalGroup(row_id, name, row['smarts'], query,
                                                      row['example'],
                                                      tuple(filter(None, row['decoys'].split(','))),
                                                      row['description'])
    return _FUNCTIONAL_CACHE


def reaction_rules() -> dict[str, tuple[ReactionRule, ...]]:
    """The whole corpus -- every row of `reactions.tsv` -- composed and cached, GROUPED BY NAME.

    The one accessor; a caller wanting a subset asks for it by `reaction=` name, and the name is
    therefore the key.  Grouped and not one row deep because A REACTION NAME NAMES A ROW FAMILY: 294
    rows under 72 names, `amidation` being three of them, one per way the acid is activated.  Same shape
    as `roles()` for the same reason, while `functional_rules()` is one row deep because a group name
    there names exactly one row.

    Insertion order is table order, inside a family and across them, so iterating `.values()` reads the
    table in the order it was written.  Caching is what keeps composition free: one lex per row per
    process, rather than a SMARTS lex inside `react()`'s loop.
    """
    global _RULES_CACHE
    if _RULES_CACHE is not None:
        return _RULES_CACHE
    out = []
    for row in read_table('reactions.tsv'):
        row_id = 'reactions:%s' % row['id']
        groups = tuple(row['groups'].split(','))
        ring_atom, ring_sizes = _parse_ring_sizes(row.get('ring_sizes', ''), row_id)
        template = _read(compose_smirks(groups, row['product']), row_id, row['name'])

        intramolecular = None
        if ring_sizes:
            if len(groups) < 2:
                raise ValueError(f'{row_id} ({row["name"]}) has one slot and a ring_sizes column. A '
                                 'one-slot row is already one molecule, so there is no intermolecular '
                                 'reading to distinguish it from: write the `r` into `product`.')
            product = _attach_ring_sizes(row['product'], ring_atom, ring_sizes)
            intramolecular = _read(compose_smirks(groups, product, intramolecular=True),
                                   row_id, row['name'])

        out.append(ReactionRule(row_id, row['name'], groups, row['product'], template,
                                row['description'], row['probe'], ring_atom, ring_sizes,
                                intramolecular))
    grouped: dict[str, list[ReactionRule]] = {}
    for rule in out:
        grouped.setdefault(rule.name, []).append(rule)
    _RULES_CACHE = {name: tuple(family) for name, family in grouped.items()}
    return _RULES_CACHE


#: The elements a `protects` value is allowed to unmask, checked against the atoms each row keeps.
#: `carbonyl` and `carboxyl` keep the carbon and rebuild the oxygen, so both name C.
PROTECTS_ELEMENTS = {'hydroxyl': 'O', 'diol': 'O', 'amine': 'N', 'thiol': 'S',
                     'carbonyl': 'C', 'carboxyl': 'C'}


def protective_rules() -> dict[str, ProtectiveGroup]:
    """Every row of `protective.tsv`, composed and cached, BY NAME, MOST SPECIFIC FIRST.

    Keyed by name because `deprotect()` selects by name and a duplicate is already refused below for
    that reason.  One row deep, not a family: a protecting group's name names one row.

    `.values()` is sorted by reactant atom count descending, table order as the tiebreak, and the order is
    load-bearing: `hydroxyl_tbu` applied to a Boc-protected alcohol yields a carbonate instead, because
    a Boc contains a tert-butyl ether's worth of atoms.  Serving the bigger pattern first lets the
    specific rule consume the site before the general one is offered it.

    Size is a proxy for specificity, not a definition -- two equal-size patterns can overlap.  What
    makes that safe is that overlap is resolved by CLAIM (`_enumerate._claims`): the sort decides who is
    asked first, the claim decides who wins.
    """
    global _PROTECTIVE_CACHE
    if _PROTECTIVE_CACHE is not None:
        return _PROTECTIVE_CACHE
    out = []
    seen = set()
    for row in read_table('protective.tsv'):
        row_id = 'protective:%s' % row['id']
        name = row['name']
        if name in seen:
            raise ValueError(f'protective.tsv: {name!r} is defined twice; `deprotect()` selects rules '
                             'by name, so a duplicate would silently shadow a pattern')
        seen.add(name)
        protects = tuple(row['protects'].split(','))
        for what in protects:
            if what not in PROTECTS_ELEMENTS:
                raise ValueError(f'{row_id} ({name}) protects {what!r}; the column takes '
                                 f'{"|".join(sorted(PROTECTS_ELEMENTS))}, comma-joined')
        template = _read('%s>>%s' % (row['smarts'], row['product']), row_id, name)
        if not template.deleted_atoms:
            raise ValueError(f'{row_id} ({name}) deletes nothing, so it is not a deprotection. A '
                             'protecting group leaves by being absent from the product side.')
        out.append(ProtectiveGroup(row_id, name, protects, row['smarts'], row['product'], template,
                                   row['protected'], row['cleaved'],
                                   tuple(filter(None, row['decoys'].split(','))), row['description']))
    # `sorted` is stable, so table order IS the tiebreak and no row needs to carry its own index; dict
    # insertion order then carries the sort, which is why `.values()` needs no re-sorting.
    _PROTECTIVE_CACHE = {rule.name: rule for rule in sorted(out, key=lambda rule: -rule.size)}
    return _PROTECTIVE_CACHE


def _read(smirks: str, row_id: str, name: str) -> ReactionTemplate:
    """One `read_smirks` call, with the composed string in the failure message.

    Nobody wrote the composed string, so a row's mistake is only findable if the message shows both it
    and the row.
    """
    try:
        return read_smirks(smirks, rule_id=row_id)
    except Exception as e:
        raise ValueError(f'{row_id} ({name}) composes to {smirks!r}, which does not read: {e}') from e


def _cap_in_product(product: str, row_id: str) -> None:
    """Refuse a `product` that states no cap.

    At load, the way `_attach_ring_sizes` refuses a ring atom nobody stated: a row whose patch cuts a
    handle and leaves no attachment point yields a fragment nothing can be coupled to, and the
    enumerator would fail on the missing map number one match later.

    A substring test is the whole of it -- a map number is stated inside a bracket, so the closing
    bracket always follows it -- and the number is `ROLE_CAP` table-wide, so an unmapped marker and a
    differently numbered one are both refused: `where[ROLE_CAP]` names the marker without a scan.
    """
    if f'#0:{ROLE_CAP}]' not in product:
        raise ValueError(f'{row_id}: {product!r} states no cap. The patch has to build the attachment '
                         f'point it leaves, spelled `[#0:{ROLE_CAP}]`.')


def _keep_at_the_site(template: ReactionTemplate, row_id: str) -> None:
    """Refuse a row whose cut site does not state `@=`.

    A cut takes no configuration away with the fragment that leaves, and a patch drops the
    configuration at its own reaction centre -- so the site of every row in this table states `@=`,
    including the 57 whose group can never match a stereogenic atom.  Uniform because the statement is
    about the CUT and not about the row: a group whose SMARTS is later widened would otherwise start
    losing configurations quietly.

    Asked of the compiled template rather than of the string: the site is whatever the cap is bonded to,
    which the patch's own bond list answers exactly and a substring test only guesses at.
    """
    cap = next(sid for sid, number in template.product_map_numbers.items() if number == ROLE_CAP)
    site = next(u if v == cap else v for u, v in template.product_bonds if cap in (u, v))
    if site not in template.product_stereo_keep:
        raise ValueError(f'{row_id}: the atom the cap hangs off states no `@=`, so a cut at a '
                         f'configured centre would come back unconfigured. Every row of this table '
                         f'keeps what it cut.')


def roles() -> dict[str, tuple[Role, ...]]:
    """Every row of `roles.tsv`, grouped by role name, composed on first use.

    Grouped rather than flat because a role IS the unit a caller selects:
    `sticky_fragments('aryl_halide')` wants all three halides at once.
    """
    global _ROLES_CACHE
    if _ROLES_CACHE is not None:
        return _ROLES_CACHE
    known = functional_rules()
    out: dict[str, list[Role]] = {}
    for row in read_table('roles.tsv'):
        row_id = 'roles:%s' % row['id']
        group = known.get(row['group'])
        if group is None:
            raise ValueError(f'{row_id} ({row["name"]}) names group {row["group"]!r}; '
                             f'functional.tsv names {len(known)} groups')
        _cap_in_product(row['product'], row_id)
        template = _read('%s>>%s' % (group.smarts, row['product']), row_id, row['name'])
        _keep_at_the_site(template, row_id)
        out.setdefault(row['name'], []).append(
            Role(row_id, row['name'], row['group'], row['product'], template,
                 row['example'], tuple(filter(None, row['decoys'].split(','))), row['description']))
    _ROLES_CACHE = {name: tuple(rows) for name, rows in out.items()}
    return _ROLES_CACHE
