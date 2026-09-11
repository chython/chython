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
"""Legacy PDB reader: fixed columns, the compatibility half of the pair.

The limits are structural -- 99,999 atoms (five-column serial), 62 chains (one column) -- and new
depositions get no legacy file from 2027-07-21.  The element comes from columns 77-78 and nowhere else,
never from the atom name.  Bonds come only from ``CONECT``, ``SSBOND`` and ``LINK``, none of which
states an order; since all three sit after the last ``ENDMDL`` and apply to every model, records are
materialised before being yielded where the mmCIF reader streams.  Every unread record type is named in
one aggregated ``unsupported:`` line.
"""
from collections.abc import Iterable, Iterator
from pathlib import Path

from ._records import PDBAtom, PDBBond, PDBRecord, normalize_element, range_messages
from .._text import require_text
from ...core import LogRecord, LOST, REPAIRED, REFUSED


__all__ = ['pdb', 'read_pdb']


#: The deck's own bookkeeping: nothing in them to read or to lose.
_STRUCTURAL = frozenset(('END', 'MASTER'))

#: The B-factor field is six columns at two decimals, so a value over 999.99 does not fit it -- stored
#: and logged, since a file carrying one was written outside the format.
_B_FACTOR_LIMIT = 999.99

#: A blank symmetry operator, or the identity one, in the `SSBOND`/`LINK` operator fields.
_IDENTITY_SYMMETRY = frozenset(('', '1555', '1_555'))

#: The two component ids wwPDB issues for water.  Legacy PDB has no entity table, so the residue name
#: read against this registry is the only statement available.
_WATER = frozenset(('HOH', 'DOD'))


def _field(line: str, start: int, end: int) -> str:
    """Columns *start*..*end*, 1-based and inclusive, stripped.  ``''`` past the end of the line.

    The strip makes the fixed-column reader line-ending agnostic, so nothing here normalises a CRLF: a
    stray CR sits past column 80 and comes off with the padding.  (The STAR lexer splits on whitespace
    and so must strip it explicitly.)
    """
    return line[start - 1:end].strip()


class _Statement:
    """One stated bond, unresolvable until a record's atoms are known.

    ``CONECT`` names serials and ``SSBOND``/``LINK`` name residues; both sit after the last ``ENDMDL``,
    so resolution is deferred and then run once per record.
    """
    __slots__ = ('kind', 'first', 'second')

    def __init__(self, kind: str, first, second):
        self.kind = kind
        self.first = first
        self.second = second


class _Damage:
    """Per-atom damage counted rather than logged line by line, keeping the first line number.

    A file whose every line is truncated is one broken writer, not a hundred thousand findings.
    """
    __slots__ = ('counts', 'first', 'elements')

    def __init__(self):
        self.counts: dict[str, int] = {}
        self.first: dict[str, int] = {}
        self.elements: dict[str, int] = {}

    def hit(self, what: str, lineno: int) -> None:
        self.counts[what] = self.counts.get(what, 0) + 1
        self.first.setdefault(what, lineno)

    def element(self, token: str) -> None:
        self.elements[token] = self.elements.get(token, 0) + 1


def _integer(text: str, what: str, lineno: int, sink) -> int | None:
    if not text:
        return None
    try:
        return int(text)
    except ValueError:
        if isinstance(sink, _Damage):
            # One counter per field name: a file whose every sequence number is bad is one writer.
            art = 'an' if what[0].lower() in 'aeioux' else 'a'
            sink.hit(f'{art} {what} that is not an integer', lineno)
        else:
            sink.append(LogRecord('pdb:unreadable-integer', (),
                                  f'atom: {what} {text!r} on line {lineno} is not an integer; not '
                                  f'stored',
                                  LOST))
        return None


def _number(text: str, what: str, lineno: int, sink) -> float | None:
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        if isinstance(sink, _Damage):
            # One counter per field name, same as _integer above.
            art = 'an' if what[0].lower() in 'aeioux' else 'a'
            sink.hit(f'{art} {what} that is not a number', lineno)
        else:
            sink.append(LogRecord('pdb:unreadable-number', (),
                                  f'atom: {what} {text!r} on line {lineno} is not a number; not stored',
                                  LOST))
        return None


def _charge(text: str, lineno: int, sink) -> int:
    """Columns 79-80.  The format spells it ``1+``; writers also spell it ``+1``."""
    if not text:
        return 0
    if text[-1] in '+-':
        text = text[-1] + text[:-1]
    try:
        return int(text)
    except ValueError:
        if isinstance(sink, _Damage):
            sink.hit('a formal charge that is not a charge', lineno)
        else:
            sink.append(LogRecord('pdb:unreadable-charge', (),
                                  f'atom: formal charge {text!r} on line {lineno} is not a charge; '
                                  f'read as 0',
                                  REPAIRED))
        return 0


def _residue_name(line: str, lineno: int, damage: _Damage) -> str:
    """The residue name, columns 18-20, widened to 21 for a writer that used the blank column.

    Column 21 is blank in the format and column 22 is the chain id.  A four-character component name
    puts its last character in 21, so the field is 18-21 whenever 21 is not blank -- and stops there,
    because chasing the run further would eat the chain id.
    """
    if line[20:21].strip():
        name = _field(line, 18, 21)
        damage.hit(f'a four-character residue name in columns 18-21, one wider than the three the '
                   f'format gives it (read as {name!r}; the chain id is still column 22)', lineno)
        return name
    return _field(line, 18, 20)


def _atom(line: str, lineno: int, damage: _Damage, log: list) -> PDBAtom | None:
    """One ``ATOM``/``HETATM`` line.  ``None`` only when there is not even a serial field to store."""
    if len(line.rstrip()) < 11:
        damage.hit('an ATOM/HETATM line that ends before its serial field, so no atom is stored',
                   lineno)
        return None
    serial = _integer(_field(line, 7, 11), 'atom serial', lineno, damage)
    if len(line.rstrip()) < 54:
        damage.hit('an ATOM/HETATM line that ends inside or before its coordinate fields, so at least '
                   'one coordinate is not stored', lineno)

    element_field = _field(line, 77, 78)
    if element_field:
        element, isotope, message = normalize_element(element_field)
        if element is None:
            damage.element(element_field)
    else:
        element, isotope = None, 0
        damage.hit('an atom with nothing in the element columns 77-78, so its element is not stored '
                   '(and is not guessed from its atom name)', lineno)

    name = _residue_name(line, lineno, damage)
    sequence = _integer(_field(line, 23, 26), 'residue sequence number', lineno, damage)

    occupancy = _number(_field(line, 55, 60), 'occupancy', lineno, damage)
    b_factor = _number(_field(line, 61, 66), 'B factor', lineno, damage)
    for what, _ in range_messages(occupancy, b_factor,
                                  f'line {lineno}' if serial is None else f'atom {serial}'):
        damage.hit(what, lineno)
    if b_factor is not None and b_factor > _B_FACTOR_LIMIT:
        damage.hit(f'a B factor over {_B_FACTOR_LIMIT}, which does not fit the six columns the format '
                   f'gives it at the two decimals it asks for', lineno)

    return PDBAtom(
        element,
        _number(_field(line, 31, 38), 'x coordinate', lineno, damage),
        _number(_field(line, 39, 46), 'y coordinate', lineno, damage),
        _number(_field(line, 47, 54), 'z coordinate', lineno, damage),
        isotope=isotope,
        charge=_charge(_field(line, 79, 80), lineno, damage),
        serial=serial,
        atom_name=_field(line, 13, 16) or None,
        residue_name=name or None,
        chain=_field(line, 22, 22) or None,
        auth_chain=_field(line, 22, 22) or None,
        auth_seq=sequence,
        residue_seq=sequence,
        ins_code=_field(line, 27, 27) or None,
        alt_loc=_field(line, 17, 17) or None,
        occupancy=occupancy,
        b_factor=b_factor,
        hetatm=line[:6].strip().upper() == 'HETATM',
        entity_type='water' if name.upper() in _WATER else None)


def _conect(line: str, lineno: int, log: list) -> list[_Statement]:
    """``CONECT``: one serial and up to four partners.  No order, in any writer."""
    serial = _integer(_field(line, 7, 11), 'CONECT serial', lineno, log)
    if serial is None:
        log.append(LogRecord('pdb:conect-no-serial', (),
                             f'bond: CONECT on line {lineno} names no atom serial; the record is '
                             f'not read'))
        return []
    statements = []
    for start in (12, 17, 22, 27):
        partner = _integer(_field(line, start, start + 4), 'CONECT partner serial', lineno, log)
        if partner is not None:
            statements.append(_Statement('conect', serial, partner))
    if _field(line, 32, 80):
        log.append(LogRecord('pdb:conect-obsolete-fields', (),
                             f'unsupported: CONECT on line {lineno} carries the obsolete '
                             f'hydrogen-bond and salt-bridge fields past column 31; they are not '
                             f'modelled',
                             LOST))
    return statements


def _residue_partner(line: str, name_columns: tuple[int, int], chain_column: int,
                     sequence_columns: tuple[int, int], ins_column: int,
                     atom_name: str | None) -> tuple:
    """The residue-and-atom tuple an ``SSBOND`` or ``LINK`` partner names."""
    sequence = _field(line, *sequence_columns)
    try:
        number = int(sequence)
    except ValueError:
        number = None            # an unparseable sequence number cannot match any atom, and says so
    return (_field(line, chain_column, chain_column) or None,
            _field(line, *name_columns) or None,
            number,
            _field(line, ins_column, ins_column) or None,
            atom_name or None)


def _symmetry(line: str, kind: str, lineno: int, log: list) -> bool:
    """Whether both partners sit in the deposited coordinates rather than in a symmetry image.

    A bond to an image would join an atom to a copy that is not in the file, so it is named, not built.
    """
    operators = [operator for operator in (_field(line, 60, 65), _field(line, 67, 72))
                 if operator not in _IDENTITY_SYMMETRY]
    if operators:
        log.append(LogRecord('pdb:symmetry-bond', (),
                             f'unsupported: {kind} on line {lineno} joins a symmetry image under '
                             f'operator {", ".join(operators)}, which is not among the coordinates '
                             f'in the file; no bond built',
                             LOST))
        return False
    return True


def _ssbond(line: str, lineno: int, log: list) -> list[_Statement]:
    """``SSBOND``: a disulfide between two named cysteines, SG to SG."""
    if not _symmetry(line, 'SSBOND', lineno, log):
        return []
    return [_Statement('ssbond',
                       _residue_partner(line, (12, 14), 16, (18, 21), 22, 'SG'),
                       _residue_partner(line, (26, 28), 30, (32, 35), 36, 'SG'))]


def _link(line: str, lineno: int, log: list) -> list[_Statement]:
    """``LINK``: an inter-residue or metal-ligand bond, with both atoms named."""
    if not _symmetry(line, 'LINK', lineno, log):
        return []
    return [_Statement('link',
                       _residue_partner(line, (18, 20), 22, (23, 26), 27, _field(line, 13, 16)),
                       _residue_partner(line, (48, 50), 52, (53, 56), 57, _field(line, 43, 46)))]


def _resolve(record: PDBRecord, statements: list[_Statement], log: list) -> None:
    """Turn the stated bonds into bonds over *record*'s own atoms."""
    by_serial: dict[int, list[int]] = {}
    for index, atom in enumerate(record.atoms):
        if atom.serial is not None:
            by_serial.setdefault(atom.serial, []).append(index)

    # Reported whether or not the file states a bond to resolve: a repeated serial is a broken writer
    # either way, and a file with no CONECT records has nothing else to notice it by.
    duplicated = sorted(serial for serial, found in by_serial.items() if len(found) > 1)
    if duplicated:
        log.append(LogRecord('pdb:duplicate-serial', (),
                             f'atom: {len(duplicated)} atom serial(s) occur more than once '
                             f'({duplicated[0]} first); a CONECT naming one is resolved to the atom '
                             f'that came first',
                             REPAIRED))

    if not statements:
        return
    by_residue: dict[tuple, list[int]] = {}
    for index, atom in enumerate(record.atoms):
        by_residue.setdefault((atom.chain, atom.residue_name, atom.residue_seq, atom.ins_code,
                               atom.atom_name), []).append(index)

    missing: dict[str, int] = {}
    missing_first: dict[str, int] = {}
    ambiguous: dict[str, int] = {}
    # How many times each pair was stated, and by which record type.  A well-formed file states each
    # CONECT pair twice, once from each atom, so the allowance is two there and one elsewhere.
    stated: dict[tuple, int] = {}
    allowance: dict[tuple, int] = {}
    seen: dict[tuple, PDBBond] = {}
    for statement in statements:
        if statement.kind == 'conect':
            first = by_serial.get(statement.first)
            second = by_serial.get(statement.second)
            if first is None or second is None:
                # Counted rather than logged per record -- a file that lost a chain names every serial
                # of it -- with the first absent serial named in the aggregate line.
                absent = statement.first if first is None else statement.second
                missing['conect'] = missing.get('conect', 0) + 1
                missing_first.setdefault('conect', absent)
                continue
            pair = (first[0], second[0])
        else:
            first = by_residue.get(statement.first)
            second = by_residue.get(statement.second)
            if first is None or second is None:
                missing[statement.kind] = missing.get(statement.kind, 0) + 1
                continue
            if len(first) > 1 or len(second) > 1:
                # SSBOND and LINK name a residue and an atom name, and neither distinguishes one
                # alternate conformer from another, so the file does not say which is meant.
                ambiguous[statement.kind] = ambiguous.get(statement.kind, 0) + 1
            pair = (first[0], second[0])
        if pair[0] == pair[1]:
            log.append(LogRecord('pdb:self-bond', (),
                                 f'bond: {statement.kind.upper()} joins an atom to itself; no bond '
                                 f'built',
                                 REFUSED))
            continue
        bond = PDBBond(pair[0], pair[1], 1, stated_order=False, source=statement.kind)
        stated[bond.key] = stated.get(bond.key, 0) + 1
        allowance[bond.key] = max(allowance.get(bond.key, 1),
                                  2 if statement.kind == 'conect' else 1)
        if bond.key in seen:
            continue
        seen[bond.key] = bond
        record.bonds.append(bond)
    repeated = sum(count - allowance[key] for key, count in stated.items()
                   if count > allowance[key])

    for kind, count in sorted(missing.items()):
        if kind == 'conect':
            log.append(LogRecord('pdb:absent-serial', (),
                                 f'bond: {count} CONECT record(s) name an atom serial that is not '
                                 f'in this record ({missing_first[kind]} first); no bond built for '
                                 f'those'))
        else:
            log.append(LogRecord('pdb:absent-residue', (),
                                 f'bond: {count} {kind.upper()} record(s) name a residue or atom '
                                 f'that is not in this record; no bond built for those'))
    for kind, count in sorted(ambiguous.items()):
        log.append(LogRecord('pdb:ambiguous-conformer', (),
                             f'bond: {count} {kind.upper()} record(s) name an atom that this record '
                             f'holds in more than one alternate conformer, and the record type has no '
                             f'altLoc field; the bond is built to the conformer that came first'))
    if repeated:
        # Past the reciprocal statement the format asks for, a repeated pair means a double bond to some
        # writers and a duplicate to others, so it is read as one single bond and reported.
        log.append(LogRecord('pdb:repeated-pair', (),
                             f'bond: {repeated} stated bond(s) restate a pair beyond the reciprocal '
                             f'CONECT the format asks for; the repetition is not read as a bond order'))
    if record.bonds:
        log.append(LogRecord('pdb:bonds-no-order', (),
                             f'bond: {len(record.bonds)} bond(s) come from CONECT, SSBOND or LINK '
                             f'records, none of which states an order; all are stored as single '
                             f'bonds'))


def _report(damage: _Damage, log: list) -> None:
    """The per-atom damage counters, one line each."""
    for what, count in sorted(damage.counts.items(), key=lambda item: damage.first[item[0]]):
        log.append(LogRecord('pdb:atom-damage', (),
                             f'atom: {count} line(s) hold {what} (first on line {damage.first[what]})'))
    if damage.elements:
        named = ', '.join(f'{token!r} ({count})'
                          for token, count in sorted(damage.elements.items()))
        log.append(LogRecord('pdb:element-not-symbol', (),
                             f'atom: element columns 77-78 hold something that is not an element '
                             f'symbol: {named}; the element is not stored for those atoms',
                             LOST))


def _iter_lines(source: str | Path | Iterable[str]) -> Iterator[str]:
    """Lines from *source*: a ``str`` is the file's text, a :class:`~pathlib.Path` is a path."""
    if isinstance(source, Path):
        with source.open(encoding='utf8', errors='replace') as f:
            yield from f
    elif isinstance(source, str):
        yield from source.splitlines()
    else:
        yield from source


def read_pdb(source: str | Path | Iterable[str], *, log: list | None = None) \
        -> Iterator[PDBRecord]:
    """Yield a :class:`PDBRecord` per ``MODEL`` in *source*, or one record for a file with none.

    Never raises on a file that is merely wrong -- a truncated line, a negative occupancy, a ``CONECT``
    naming an absent serial -- which is stored as far as it can be, and logged.
    """
    log = [] if log is None else log
    # Damage found while scanning belongs to the file, not to one model: given to the caller once and
    # prepended to every record's own log, so a record reads on its own without repeating it per model.
    file_log: list = []
    damage = _Damage()
    statements: list[_Statement] = []
    unread: dict[str, int] = {}
    records: list[PDBRecord] = []
    current: PDBRecord | None = None
    entry: str | None = None
    title_parts: list[str] = []
    model: int | None = None

    def record_for(model_number: int | None) -> PDBRecord:
        nonlocal current
        if current is None:
            current = PDBRecord(model=model_number)
            records.append(current)
        return current

    for lineno, raw in enumerate(_iter_lines(source), 1):
        line = raw.rstrip('\n')             # a CRLF's CR needs nothing here; see `_field`
        tag = line[:6].strip().upper()

        if tag in ('ATOM', 'HETATM'):
            atom = _atom(line, lineno, damage, file_log)
            if atom is not None:
                atom.model = model
                record_for(model).atoms.append(atom)
        elif tag == 'MODEL':
            model = _integer(_field(line, 11, 14), 'model serial', lineno, file_log)
            current = None
            record_for(model)
        elif tag == 'ENDMDL':
            current = None
        elif tag == 'CONECT':
            statements.extend(_conect(line, lineno, file_log))
        elif tag == 'SSBOND':
            statements.extend(_ssbond(line, lineno, file_log))
        elif tag == 'LINK':
            statements.extend(_link(line, lineno, file_log))
        elif tag == 'HEADER':
            entry = _field(line, 63, 66) or None
        elif tag == 'TITLE':
            title_parts.append(_field(line, 11, 80))
        elif tag in _STRUCTURAL or not tag:
            continue
        else:
            unread[tag] = unread.get(tag, 0) + 1

    title = ' '.join(part for part in title_parts if part) or None
    _report(damage, file_log)
    if unread:
        named = ', '.join(f'{tag} ({count})' for tag, count in sorted(unread.items()))
        file_log.append(LogRecord('pdb:unread-records', (),
                                  f'unsupported: {sum(unread.values())} record(s) of {len(unread)} '
                                  f'type(s) are not modelled: {named}',
                                  LOST))
    if not records:
        file_log.append(LogRecord('pdb:no-atoms', (),
                                  'record: the file states no ATOM or HETATM records; no atoms read',
                                  LOST))
    log.extend(file_log)

    for record in records:
        record.entry_id = entry
        record.title = title
        own: list = []
        if len(records) > 1:
            own.append(LogRecord('pdb:model-selected', (),
                                 f'record: the file states {len(records)} models; this record holds '
                                 f'model {record.model}'))
        _resolve(record, statements, own)
        if not record.bonds and record.atoms:
            own.append(LogRecord('pdb:no-connectivity', (),
                                 f'bond: the file states no connectivity for this record; '
                                 f'{len(record.atoms)} atom(s) are unbonded'))
        else:
            unbonded = record.unbonded_count()
            if unbonded:
                own.append(LogRecord('pdb:unbonded-atoms', (),
                                     f'bond: {unbonded} atom(s) are touched by no stated bond'))
        record.log = file_log + own
        log.extend(own)
        yield record


def pdb(text: str, *, log: list | None = None) -> list[PDBRecord]:
    """Every record in the legacy PDB *text*, as a list.  The string form of :func:`read_pdb`.

    *text* is always the file's text and never a path: :func:`read_pdb` is the file reader.
    """
    return list(read_pdb(require_text(text, 'pdb'), log=log))
