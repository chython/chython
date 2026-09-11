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
"""PDBx/mmCIF reader: the archive's canonical format and the release-critical half of the pair.

Read, per category: ``_atom_site``, ``_entity``, ``_chem_comp_bond`` (orders), ``_struct_conn``,
``_entry.id`` and ``_struct.title``; every other category, and every unread ``_atom_site`` item, is
named in one aggregated ``unsupported:`` line.  Bonds are only ever the ones those tables state -- no
cutoff, no template, no saturation.  Alternate conformers never bond to each other.  One record per
``pdbx_PDB_model_num``, an NMR ensemble being alternative structures rather than duplicated atoms.
"""
from collections.abc import Iterable, Iterator
from pathlib import Path

from ._records import PDBAtom, PDBBond, PDBRecord, normalize_element, range_messages
from .._text import require_text
from ._star import StarLoop, is_null, parse_star
from ...core import LogRecord, LOST, REPAIRED, REFUSED


__all__ = ['mmcif', 'read_mmcif']


# --------------------------------------------------------------------------- vocabulary

#: `_chem_comp_bond.value_order` / `_struct_conn.pdbx_value_order` spellings chython can hold.  A stated
#: AROM is chython's order 4, stored as such and never kekulised.
_ORDERS = {'sing': 1, 'doub': 2, 'trip': 3, 'arom': 4}

#: Orders the format states and chython has no bond for.  Each is stored as a single bond so the
#: connectivity survives, and named in the log so the approximation is not silent.
_UNHELD_ORDERS = {'quad': 'quadruple', 'poly': 'polymeric', 'delo': 'delocalised', 'pi': 'pi'}

#: `_struct_conn.conn_type_id` values that are not covalent connectivity at all.  Skipped, counted.
_NON_COVALENT_CONN = {'hydrog': 'hydrogen bond', 'saltbr': 'salt bridge',
                      'mismat': 'mismatched base pair'}

#: `_struct_conn.conn_type_id` values that are a covalent single bond whose order the row may state in
#: `pdbx_value_order`.  `covale_base`/`covale_phosphate`/`covale_sugar` are the nucleic-acid link
#: spellings and `modres` links a modified residue to its parent component -- all ordinary dictionary
#: values.  The row's own spelling is kept on the bond, so a caller can filter for `modres`.
_COVALENT_CONN = frozenset(('covale', 'covale_base', 'covale_phosphate', 'covale_sugar', 'modres'))

#: The identity operator.  A `_struct_conn` partner under any other operator is an atom of a symmetry
#: image that is not in the deposited coordinates, so the bond has no second end to attach to here.
_IDENTITY_SYMMETRY = '1_555'

#: `_atom_site` items this reader reads.  Anything else in the category is reported, so adding a read
#: item means adding it here too, which is what stops the two drifting apart silently.
_ATOM_SITE_READ = frozenset((
    '_atom_site.group_pdb', '_atom_site.id', '_atom_site.type_symbol', '_atom_site.label_atom_id',
    '_atom_site.label_alt_id', '_atom_site.label_comp_id', '_atom_site.label_asym_id',
    '_atom_site.label_entity_id', '_atom_site.label_seq_id', '_atom_site.pdbx_pdb_ins_code',
    '_atom_site.cartn_x', '_atom_site.cartn_y', '_atom_site.cartn_z', '_atom_site.occupancy',
    '_atom_site.b_iso_or_equiv', '_atom_site.pdbx_formal_charge', '_atom_site.auth_seq_id',
    '_atom_site.auth_asym_id', '_atom_site.pdbx_pdb_model_num'))
# `auth_comp_id` and `auth_atom_id` are deliberately absent: nothing keeps them per atom, so they belong
# in the unread-item line.  A `_struct_conn` partner is still searched for by them.

#: Categories this reader consumes.  Everything else present in the block goes into the aggregated
#: `unsupported:` line.
_CATEGORIES_READ = frozenset(('_atom_site', '_chem_comp_bond', '_struct_conn', '_entity', '_entry',
                              '_struct'))


# --------------------------------------------------------------------------- damage, counted

class _Damage:
    """One line per *kind* of malformation, however many rows carry it.

    The first row's message is kept verbatim, since it names a row, and the count of the rest is
    appended to it.  The legacy reader's class of the same name reports a different message shape and
    is deliberately not shared.
    """
    __slots__ = ('_first', '_counts')

    def __init__(self):
        self._first: dict[str, str] = {}
        self._counts: dict[str, int] = {}

    def hit(self, what: str, message: str) -> None:
        """Record one occurrence of *what*, keeping *message* if it is the first."""
        self._counts[what] = self._counts.get(what, 0) + 1
        self._first.setdefault(what, message)

    def report(self, log: list) -> None:
        """Append one line per kind, in the order the kinds were first seen."""
        for what, message in self._first.items():
            more = self._counts[what] - 1
            log.append(LogRecord('mmcif:damage-report', (),
                                 message if not more else f'{message} (and {more} more row(s))'))


# --------------------------------------------------------------------------- field readers

def _text(value) -> str | None:
    """A CIF value as text, with both nulls collapsed to ``None`` -- they stay distinct in the block."""
    if value is None or is_null(value):
        return None
    return value


def _integer(value, what: str, where: str, damage: _Damage) -> int | None:
    """An integer field, or ``None`` with a counted line when the text is not one."""
    text = _text(value)
    if text is None:
        return None
    try:
        return int(text)
    except ValueError:
        damage.hit(f'{what} is not an integer',
                   f'atom: {what} {text!r} on {where} is not an integer; not stored')
        return None


def _number(value, what: str, where: str, damage: _Damage) -> float | None:
    """A float field, or ``None`` with a counted line.  Trailing esd in parentheses is stripped."""
    text = _text(value)
    if text is None:
        return None
    try:
        return float(text)
    except ValueError:
        head = text.split('(', 1)[0]
        try:
            number = float(head)
        except ValueError:
            damage.hit(f'{what} is not a number',
                       f'atom: {what} {text!r} on {where} is not a number; not stored')
            return None
        damage.hit(f'{what} carries an uncertainty',
                   f'atom: {what} {text!r} on {where} carries an uncertainty; read as {number}')
        return number


def _charge(value, where: str, damage: _Damage) -> int:
    """A formal charge.  ``2-`` and ``-2`` are both written by real files and both are read."""
    text = _text(value)
    if text is None:
        return 0
    text = text.strip()
    if text[-1:] in '+-' and len(text) > 1:      # trailing-sign spelling
        text = text[-1] + text[:-1]
    try:
        return int(text)
    except ValueError:
        damage.hit('formal charge is not a charge',
                   f'atom: formal charge {_text(value)!r} on {where} is not a charge; read as 0')
        return 0


def _order(value, damage: _Damage, where: str) -> tuple:
    """``(order, stated)`` for a ``value_order`` field.

    A field stating nothing is silent here and counted by the caller -- it is the commonest omission
    in a ``_struct_conn`` row, so a line each would bury the log.
    """
    text = _text(value)
    if text is None:
        return 1, False
    key = text.strip().lower()
    if key in _ORDERS:
        return _ORDERS[key], True
    if key in _UNHELD_ORDERS:
        damage.hit(f'{key} bond order',
                   f'unsupported: {_UNHELD_ORDERS[key]} bond order stated on {where} is not '
                   f'modelled; the bond is stored as a single bond')
        return 1, False
    damage.hit(f'bond order {key!r} is not a spelling',
               f'bond: bond order {text!r} on {where} is not a value_order spelling; the bond is '
               f'stored as a single bond')
    return 1, False


# --------------------------------------------------------------------------- atoms

def _entity_types(block) -> dict[str, str]:
    """``label_entity_id`` -> ``_entity.type``, empty when the file states no entities."""
    loop = block.loop('_entity.id')
    if loop is None:
        one = _text(block.get('_entity.id'))
        if one is None:
            return {}
        return {one: (_text(block.get('_entity.type')) or '').lower()}
    types = {}
    for row in loop.rows:
        key = _text(loop.value(row, '_entity.id'))
        if key is None:
            continue
        types[key] = (_text(loop.value(row, '_entity.type')) or '').lower()
    return types


def _atom_site_rows(block):
    """The ``_atom_site`` rows as ``(loop, rows)``, tolerating a single atom written as scalars.

    A one-atom block may state ``_atom_site.id 1`` as plain items rather than a ``loop_``; the format
    permits it and a reader that only looks for the loop finds no atoms at all.
    """
    loop = block.loop('_atom_site.id')
    if loop is not None:
        return loop, loop.rows
    scalars = {t: v for t, v in block.items.items() if t.startswith('_atom_site.')}
    if not scalars:
        return None, []
    tags = sorted(scalars)
    return StarLoop(tags, [[scalars[t] for t in tags]]), None


def _read_atoms(block, log: list) -> list[PDBAtom]:
    """Every ``_atom_site`` row as a :class:`PDBAtom`, in file order."""
    loop, rows = _atom_site_rows(block)
    if loop is None:
        log.append(LogRecord('mmcif:no-atom-site', (),
                             'atom: the block states no _atom_site rows; no atoms read',
                             LOST))
        return []
    if rows is None:
        rows = loop.rows
    entities = _entity_types(block)
    damage = _Damage()
    atoms = []
    for number, row in enumerate(rows, 1):
        serial = _integer(loop.value(row, '_atom_site.id'), 'atom id', f'_atom_site row {number}',
                          damage)
        where = f'_atom_site row {number}' if serial is None else f'atom {serial}'

        symbol = _text(loop.value(row, '_atom_site.type_symbol')) or ''
        element, isotope, message = normalize_element(symbol)
        if message is not None:
            damage.hit(f'element field {symbol!r}', f'{message} ({where})')
        elif element is None:
            damage.hit('no element stated',
                       f'atom: no element stated on {where}; element not stored')

        group = (_text(loop.value(row, '_atom_site.group_pdb')) or '').upper()
        entity = _text(loop.value(row, '_atom_site.label_entity_id'))
        occupancy = _number(loop.value(row, '_atom_site.occupancy'), 'occupancy', where, damage)
        b_factor = _number(loop.value(row, '_atom_site.b_iso_or_equiv'), 'B factor', where, damage)
        for what, message in range_messages(occupancy, b_factor, where):
            damage.hit(what, message)
        atoms.append(PDBAtom(
            element,
            _number(loop.value(row, '_atom_site.cartn_x'), 'x coordinate', where, damage),
            _number(loop.value(row, '_atom_site.cartn_y'), 'y coordinate', where, damage),
            _number(loop.value(row, '_atom_site.cartn_z'), 'z coordinate', where, damage),
            isotope=isotope,
            charge=_charge(loop.value(row, '_atom_site.pdbx_formal_charge'), where, damage),
            serial=serial,
            atom_name=_text(loop.value(row, '_atom_site.label_atom_id')),
            residue_name=_text(loop.value(row, '_atom_site.label_comp_id')),
            chain=_text(loop.value(row, '_atom_site.label_asym_id')),
            auth_chain=_text(loop.value(row, '_atom_site.auth_asym_id')),
            auth_seq=_integer(loop.value(row, '_atom_site.auth_seq_id'), 'auth_seq_id', where,
                              damage),
            residue_seq=_integer(loop.value(row, '_atom_site.label_seq_id'), 'label_seq_id', where,
                                 damage),
            ins_code=_text(loop.value(row, '_atom_site.pdbx_pdb_ins_code')),
            alt_loc=_text(loop.value(row, '_atom_site.label_alt_id')),
            occupancy=occupancy,
            b_factor=b_factor,
            hetatm=group == 'HETATM',
            entity_type=entities.get(entity) if entity is not None else None,
            model=_integer(loop.value(row, '_atom_site.pdbx_pdb_model_num'), 'model number', where,
                           damage)))
    damage.report(log)
    _report_atom_site_items(loop, log)
    return atoms


def _report_atom_site_items(loop, log: list) -> None:
    """Name every ``_atom_site`` item the reader did not read, in one line."""
    extra = sorted(t for t in loop.tags if t not in _ATOM_SITE_READ)
    if extra:
        log.append(LogRecord('mmcif:unread-atom-site-items', (),
                             f'unsupported: {len(extra)} _atom_site item(s) not modelled: '
                             + ', '.join(extra),
                             LOST))


# --------------------------------------------------------------------------- bonds

def _alt_compatible(one: str | None, other: str | None) -> bool:
    """Whether two atoms may be bonded, given their alternate-conformer ids.

    An atom with no alt id belongs to every conformer; two atoms with different alt ids describe the
    same place and are never bonded to each other.
    """
    return one is None or other is None or one == other


def _component_bonds(block, log: list) -> dict[str, list]:
    """``comp_id`` -> ``[(atom_name_1, atom_name_2, order, stated_order), ...]``."""
    loop = block.loop('_chem_comp_bond.comp_id')
    if loop is None:
        return {}
    table: dict[str, list] = {}
    damage = _Damage()
    stereo = unstated = 0
    for number, row in enumerate(loop.rows, 1):
        comp = _text(loop.value(row, '_chem_comp_bond.comp_id'))
        one = _text(loop.value(row, '_chem_comp_bond.atom_id_1'))
        other = _text(loop.value(row, '_chem_comp_bond.atom_id_2'))
        if comp is None or one is None or other is None:
            damage.hit('a row naming no component and two atoms',
                       f'bond: _chem_comp_bond row {number} does not name a component and two '
                       f'atoms; the row is not read')
            continue
        value = loop.value(row, '_chem_comp_bond.value_order')
        if _text(value) is None and (_text(loop.value(row, '_chem_comp_bond.pdbx_aromatic_flag'))
                                     or '').upper() == 'Y':
            # The aromatic flag is the only order statement on this row, so it applies; where a
            # value_order is stated too, its Kekule order is strictly more information.
            order, stated = 4, True
        else:
            if _text(value) is None:
                unstated += 1
            order, stated = _order(value, damage, f'_chem_comp_bond row {number}')
        if (_text(loop.value(row, '_chem_comp_bond.pdbx_stereo_config')) or 'N').upper() != 'N':
            stereo += 1
        table.setdefault(comp, []).append((one, other, order, stated))
    damage.report(log)
    if unstated:
        log.append(LogRecord('mmcif:chem-comp-no-order', (),
                             f'bond: {unstated} _chem_comp_bond row(s) state no bond order at all; '
                             f'each is stored as a single bond'))
    if stereo:
        log.append(LogRecord('mmcif:chem-comp-stereo', (),
                             f'unsupported: {stereo} _chem_comp_bond row(s) state a bond stereo '
                             f'configuration that is not modelled',
                             LOST))
    return table


def _residue_atoms(record: PDBRecord) -> dict[tuple, dict[str, list]]:
    """``residue_key`` -> ``atom_name`` -> ``[(index, alt_loc), ...]`` over one record."""
    residues: dict[tuple, dict[str, list]] = {}
    for index, atom in enumerate(record.atoms):
        if atom.atom_name is None:
            continue
        residues.setdefault(atom.residue_key, {}).setdefault(atom.atom_name, []) \
            .append((index, atom.alt_loc))
    return residues


def _add_component_bonds(record: PDBRecord, table: dict[str, list], seen: dict, log: list) -> None:
    """Every ``_chem_comp_bond`` row, applied to every instance of its component in this record."""
    if not table:
        return
    absent = 0
    for key, by_name in _residue_atoms(record).items():
        rows = table.get(key[2])           # residue_name
        if not rows:
            continue
        for one, other, order, stated in rows:
            first, second = by_name.get(one), by_name.get(other)
            if not first or not second:
                absent += 1
                continue
            for a, alt_a in first:
                for b, alt_b in second:
                    if a == b or not _alt_compatible(alt_a, alt_b):
                        continue
                    _append_bond(record, PDBBond(a, b, order, stated_order=stated,
                                                 source='chem_comp_bond'), seen, log)
    if absent:
        # Routine rather than damage: the component definition lists every atom including hydrogens,
        # and a crystal structure states coordinates for a subset of them.
        log.append(LogRecord('mmcif:chem-comp-absent-atom', (),
                             f'bond: {absent} _chem_comp_bond row(s) name an atom that has no '
                             f'coordinates in this model; no bond built for those'))


def _atom_indexes(record: PDBRecord) -> tuple:
    """Two lookups for a ``_struct_conn`` partner: by label identifiers and by auth identifiers.

    Neither is redundant -- ``label_seq_id`` is null for every non-polymer, so a ligand or metal
    partner is findable only by its auth numbering.
    """
    label: dict[tuple, list] = {}
    auth: dict[tuple, list] = {}
    for index, atom in enumerate(record.atoms):
        label.setdefault((atom.chain, atom.residue_name, atom.residue_seq, atom.ins_code,
                          atom.atom_name), []).append((index, atom.alt_loc))
        auth.setdefault((atom.auth_chain, atom.residue_name, atom.auth_seq, atom.ins_code,
                         atom.atom_name), []).append((index, atom.alt_loc))
    return label, auth


def _partner(loop, row, side: str, label: dict, auth: dict) -> tuple:
    """``(candidates, why_not)`` for one ``_struct_conn`` partner, as ``(index, alt_loc)`` pairs.

    A row naming an alternate-conformer id means the bond exists only in that conformer, so candidates
    are narrowed to it; a row naming none leaves every conformer a candidate for
    :func:`_alt_compatible`.  ``why_not`` distinguishes the three ways this comes back empty -- no
    sequence number stated, a stated one that misses, a conformer the model does not contain -- and is
    never ``None`` on an empty list, since the caller interpolates it into a message.
    """
    name = _text(loop.value(row, f'_struct_conn.ptnr{side}_label_atom_id')) \
        or _text(loop.value(row, f'_struct_conn.ptnr{side}_auth_atom_id'))
    comp = _text(loop.value(row, f'_struct_conn.ptnr{side}_label_comp_id')) \
        or _text(loop.value(row, f'_struct_conn.ptnr{side}_auth_comp_id'))
    ins = _text(loop.value(row, f'_struct_conn.pdbx_ptnr{side}_pdb_ins_code'))
    seq = _text(loop.value(row, f'_struct_conn.ptnr{side}_label_seq_id'))
    chain = _text(loop.value(row, f'_struct_conn.ptnr{side}_label_asym_id'))
    found = None
    searched = False
    if seq is not None and chain is not None:
        searched = True
        try:
            found = label.get((chain, comp, int(seq), ins, name))
        except ValueError:
            found = None
    if not found:
        auth_seq = _text(loop.value(row, f'_struct_conn.ptnr{side}_auth_seq_id'))
        auth_chain = _text(loop.value(row, f'_struct_conn.ptnr{side}_auth_asym_id'))
        if auth_seq is not None and auth_chain is not None:
            searched = True
            try:
                found = auth.get((auth_chain, comp, int(auth_seq), ins, name))
            except ValueError:
                found = None
    if not found:
        return [], ('names an atom that is not in this model' if searched else
                    f'states neither a label nor an auth sequence number for partner {side}, so it '
                    f'names no atom this model can be searched for')
    alt = _text(loop.value(row, f'_struct_conn.pdbx_ptnr{side}_label_alt_id'))
    if alt is None:
        return found, None
    narrowed = [pair for pair in found if pair[1] is None or pair[1] == alt]
    if not narrowed:
        return [], (f'names alternate conformer {alt!r} of partner {side}, which this model does not '
                    f'contain')
    return narrowed, None


def _add_struct_conn(record: PDBRecord, block, seen: dict, log: list) -> None:
    """Every inter-component bond ``_struct_conn`` states: disulfide, covalent link, metal contact."""
    loop = block.loop('_struct_conn.conn_type_id')
    if loop is None:
        return
    label, auth = _atom_indexes(record)
    damage = _Damage()
    non_covalent: dict[str, int] = {}
    symmetry: dict[str, int] = {}
    metal = unstated = 0
    for number, row in enumerate(loop.rows, 1):
        kind = (_text(loop.value(row, '_struct_conn.conn_type_id')) or '').lower()
        where = f'_struct_conn row {number}'
        if kind in _NON_COVALENT_CONN:
            non_covalent[kind] = non_covalent.get(kind, 0) + 1
            continue
        one = _text(loop.value(row, '_struct_conn.ptnr1_symmetry')) or _IDENTITY_SYMMETRY
        other = _text(loop.value(row, '_struct_conn.ptnr2_symmetry')) or _IDENTITY_SYMMETRY
        if one != _IDENTITY_SYMMETRY or other != _IDENTITY_SYMMETRY:
            operator = one if one != _IDENTITY_SYMMETRY else other
            symmetry[operator] = symmetry.get(operator, 0) + 1
            continue
        first, why_first = _partner(loop, row, '1', label, auth)
        second, why_second = _partner(loop, row, '2', label, auth)
        if not first or not second:
            log.append(LogRecord('mmcif:struct-conn-no-partner', (),
                                 f'bond: {where} {why_first or why_second}; no bond built'))
            continue
        # The pairs are worked out before the row's order is, so a row whose every pair is discarded is
        # not counted among the bonds the aggregate lines below report.
        pairs = []
        self_referential = False
        for a, alt_a in first:
            for b, alt_b in second:
                if a == b:
                    self_referential = True
                elif _alt_compatible(alt_a, alt_b):
                    pairs.append((a, b))
        if not pairs:
            if self_referential:
                log.append(LogRecord('mmcif:self-bond', (),
                                     f'bond: {where} joins an atom to itself; no bond built',
                                     REFUSED))
            continue
        if kind == 'metalc':
            order, stated = 8, True
            metal += 1
        elif kind == 'disulf':
            order, stated = 1, True      # the connection type states the order
        elif kind in _COVALENT_CONN:
            value = loop.value(row, '_struct_conn.pdbx_value_order')
            if _text(value) is None:
                unstated += 1
            order, stated = _order(value, damage, where)
        else:
            damage.hit(f'connection type {kind!r}',
                       f'bond: {where} states connection type {kind!r}, which is not a connection '
                       f'type this reader knows; the bond is stored as a single bond')
            order, stated = 1, False
        for a, b in pairs:
            _append_bond(record, PDBBond(a, b, order, stated_order=stated, source='struct_conn',
                                         conn_type=kind), seen, log)
    damage.report(log)
    for kind, count in sorted(non_covalent.items()):
        log.append(LogRecord('mmcif:non-covalent-conn', (),
                             f'unsupported: {count} _struct_conn row(s) state a '
                             f'{_NON_COVALENT_CONN[kind]}, which is not a bond order chython holds; '
                             f'no bond built for those',
                             LOST))
    for operator, count in sorted(symmetry.items()):
        log.append(LogRecord('mmcif:symmetry-conn', (),
                             f'unsupported: {count} _struct_conn row(s) join an atom under symmetry '
                             f'operator {operator}, whose image is not among the deposited '
                             f'coordinates; no bond built for those',
                             LOST))
    if unstated:
        log.append(LogRecord('mmcif:struct-conn-no-order', (),
                             f'bond: {unstated} _struct_conn covalent link(s) state no bond order; '
                             f'each is stored as a single bond'))
    if metal:
        log.append(LogRecord('mmcif:metal-conn-no-direction', (),
                             f'bond: {metal} metal coordination bond(s) stored as coordination bonds '
                             f'in the order the file names the partners; _struct_conn states no donor '
                             f'direction'))


def _append_bond(record: PDBRecord, bond: PDBBond, seen: dict, log: list) -> None:
    """Add *bond* unless the pair is already bonded; log a second statement of a different order."""
    key = bond.key
    known = seen.get(key)
    if known is not None:
        if known.order != bond.order:
            log.append(LogRecord('mmcif:duplicate-bond', (),
                                 f'bond: atoms {key[0]} and {key[1]} are bonded twice, as order '
                                 f'{known.order} by {known.source} and as order {bond.order} by '
                                 f'{bond.source}; the first statement is the one kept',
                                 REPAIRED))
        return
    seen[key] = bond
    record.bonds.append(bond)


# --------------------------------------------------------------------------- records

def _split_models(atoms: list[PDBAtom], block) -> list[PDBRecord]:
    """One record per ``pdbx_PDB_model_num``, in the order the models first appear."""
    entry = _text(block.get('_entry.id')) or block.name
    title = _text(block.get('_struct.title'))
    records: dict[int | None, PDBRecord] = {}
    for atom in atoms:
        record = records.get(atom.model)
        if record is None:
            record = records[atom.model] = PDBRecord(entry_id=entry, title=title, model=atom.model)
        record.atoms.append(atom)
    if not records:
        return [PDBRecord(entry_id=entry, title=title)]
    return list(records.values())


def _is_polymer_monomer(atom: PDBAtom) -> bool:
    """Whether *atom* belongs to a polymer, by what the file states and by nothing weaker.

    ``_entity`` is the direct statement and wins where the file carries it; it is not mandatory mmCIF.
    The fallback is a dictionary fact rather than an inference: ``_atom_site.label_seq_id`` is defined
    only for a polymer entity, so a non-polymer row carries the *inapplicable* null there.
    """
    if atom.entity_type is not None:
        return atom.entity_type == 'polymer'
    return atom.residue_seq is not None


def _report_polymer_linkage(record: PDBRecord, log: list) -> None:
    """Name the one kind of bond an archive entry does not state as an atom pair.

    The bond joining one polymer monomer to the next is stated as a *sequence* -- ``_entity_poly_seq``
    plus the component's ``_chem_comp.type`` -- and never as a pair of atom names, so building it needs
    the per-component attachment-point table and a separate explicit pass.  It is therefore not built,
    and the log says how many monomers are bonded only within themselves.
    """
    monomers = {atom.residue_key for atom in record.atoms if _is_polymer_monomer(atom)}
    if len(monomers) > 1:
        log.append(LogRecord('mmcif:polymer-linkage', (),
                             f'unsupported: the file states its polymer linkage as a sequence and '
                             f'not as pairs of atoms, so no bond is built between consecutive '
                             f'monomers; {len(monomers)} polymer monomer(s) are bonded only within '
                             f'themselves',
                             LOST))


def _report_categories(block, log: list) -> None:
    """Name every category present and not read, in one line."""
    extra = sorted(block.categories() - _CATEGORIES_READ)
    if extra:
        log.append(LogRecord('mmcif:unread-categories', (),
                             f'unsupported: {len(extra)} mmCIF category(ies) not modelled: '
                             + ', '.join(extra),
                             LOST))


def _read_block(block, log: list) -> list[PDBRecord]:
    """Every record one ``data_`` block yields."""
    block_log: list = []
    atoms = _read_atoms(block, block_log)
    table = _component_bonds(block, block_log)
    _report_categories(block, block_log)

    records = _split_models(atoms, block)
    # Block-level damage belongs to every record the block yielded, a record being self-contained, and
    # to the caller's flat list exactly once -- hence two separate extends.
    log.extend(block_log)
    for record in records:
        own: list = []
        if len(records) > 1:
            own.append(LogRecord('mmcif:model-selected', (),
                                 f'record: the block states {len(records)} models; this record holds '
                                 f'model {record.model}'))
        seen: dict = {}
        _add_component_bonds(record, table, seen, own)
        _add_struct_conn(record, block, seen, own)
        _report_polymer_linkage(record, own)
        if not record.bonds and record.atoms:
            own.append(LogRecord('mmcif:no-connectivity', (),
                                 f'bond: the file states no connectivity for this record; '
                                 f'{len(record.atoms)} atom(s) are unbonded'))
        else:
            unbonded = record.unbonded_count()
            if unbonded:
                own.append(LogRecord('mmcif:unbonded-atoms', (),
                                     f'bond: {unbonded} atom(s) are touched by no stated bond'))
        record.log = block_log + own
        log.extend(own)
    return records


# --------------------------------------------------------------------------- public readers

def read_mmcif(source: Path | Iterable[str], *, log: list | None = None) \
        -> Iterator[PDBRecord]:
    """Yield a :class:`PDBRecord` per model per ``data_`` block in *source*.

    *source* is a :class:`~pathlib.Path`, an open file object, or any iterable of lines -- read once,
    line by line, so a 25 MB entry is never materialised whole.  A ``str`` is always the file's text;
    :func:`mmcif` is the name to reach for there.  Never raises on a file that is merely wrong.
    """
    log = [] if log is None else log
    for block in parse_star(source, log=log):
        yield from _read_block(block, log)


def mmcif(text: str, *, log: list | None = None) -> list[PDBRecord]:
    """Every record in the mmCIF *text*, as a list.  The string form of :func:`read_mmcif`.

    *text* is always the file's text and never a path: :func:`read_mmcif` is the file reader.
    """
    return list(read_mmcif(require_text(text, 'mmcif'), log=log))
