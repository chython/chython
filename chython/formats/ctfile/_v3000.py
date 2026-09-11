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
"""MDL V3000: CTAB in, CTAB out.  Reading recovers; writing stays inside the specification.

The trap is the index alphabet: a V3000 atom index is an arbitrary positive integer and not a
position -- a file may number its atoms 1, 7, 100 -- so every reference the file makes (bond
endpoints, S-group atom lists, collection members) is resolved through a map and stored as a 0-based
position.  Treating an index as a position works on almost every file, then builds a different one.
"""

from ._ctab import (Ctab, CtabAtom, CtabBond, order_from_bond_type, LABEL_ELEMENT,
                    STEREO_FROM_COLLECTION, STEREO_TO_COLLECTION, WEDGE_FROM_V3000,
                    WEDGE_TO_V3000)
from ._errors import MalformedCtfile, UnsupportedCtfile
from ._hydrogens import (MRV_IMPLICIT_H, ZERO_VALENCE, apply_mrv_implicit_h, implicit_h_records,
                         valence_for_write)
from ._sgroup import (NO_INDEX, UNSUPPORTED, SGroup, checked_index, format_fielddisp,
                      merge_log, normalize_indices, parse_fielddisp, resolve_output)
from ._tokens import emit_v30, join_continuations, parse_list, quote_value, tokenize
from ...core import LogRecord, LOST, R_INDEX_MAX, REPAIRED
from ...core._core import element_symbols
from ...core.wedge import wedge_in_file_order, wedges_for_write


__all__ = ['parse_v3000', 'emit_v3000', 'parse_ctab_block', 'V3000_STAMP']


V3000_STAMP = 'V3000'
_SYMBOLS = element_symbols()
_SYMBOL_SET = frozenset(_SYMBOLS[1:])
# Case-insensitive lookup for the recovery below.  Real files write CL, BR, and occasionally cl.
_SYMBOL_FOLD = {s.upper(): s for s in _SYMBOLS[1:]}

# Atom-type tokens that are valid CTfile and are CONSTRAINTS: each stands for a set of elements, which
# a molecule cannot hold, so each is refused by name -- "this is a query file" is actionable, "unknown
# element 'Q'" is not.  `LP`, `Pol`, `Mod` and `Dummy` are NOT here: they name one atom each rather
# than a set, so they read as the marker with the token as its alias, like any other label.
_QUERY_TYPES = {
    'A': 'any heavy atom', 'AH': 'any atom', 'Q': 'any heteroatom', 'QH': 'any heteroatom or H',
    'M': 'any metal', 'MH': 'any metal or H', 'X': 'any halogen', 'XH': 'any halogen or H',
    'L': 'atom list',
}
# Hydrogen isotope aliases, not query features: a file that writes D means a deuterium.
_H_ISOTOPES = {'D': 2, 'T': 3}

# Bond keywords whose value this release models.  Anything else on a bond line is logged by name and
# dropped: a query keyword silently kept would make the molecule claim a constraint it does not hold.
_BOND_MODELLED = frozenset(('CFG', 'TOPO', 'RXCTR'))

# S-group keywords this module turns into structure.  Everything else rides in `SGroup.fields` as a
# list of raw values -- a list, because CSTATE, BRKXYZ and FIELDDATA all legitimately repeat on one
# record and a dict of strings would keep the last one only.
_SGROUP_ATOM_KEYS = frozenset(('ATOMS', 'PATOMS'))
_SGROUP_BOND_KEYS = frozenset(('XBONDS', 'CBONDS'))

# Every V3000 keyword whose value contains an index, and what it indexes.  `SGroup.fields` is
# re-emitted verbatim while the writer *regenerates* atom, bond and S-group numbering, so a keyword
# holding indices that rides through untranslated names different atoms on the way out -- silently,
# and only on files not already numbered `1..n` in order.
#
# So an index-valued keyword has two permitted fates and passing through is not one of them:
# translated into stored references (`ATOMS`, `PATOMS`, `CBONDS`, `XBONDS`, `CSTATE`, `PARENT` and the
# collection `ATOMS` all are), or dropped by name with a log line.  `ENDPTS` (bond block) and
# `ATTCHORD` (atom block) are listed for completeness; those blocks drop what they do not model.
_INDEX_VALUED = {
    'SAP': 'atom',        # superatom attachment point: (3 atom leaving-atom id)
    'XBHEAD': 'bond',     # a multiple group's crossing bonds
    'XBCORR': 'bond',     # ... and the pairs correlating them
    'ENDPTS': 'atom',     # bond block: a variable-attachment bond's candidate atoms
    'ATTCHORD': 'atom',   # atom block: R-group template attachment order
    'MEMBERS': 'atom',    # collection block
    'OBJ3DS': 'object',
}

# The three enhanced-stereo collection prefixes, and the tag this library maps each onto.
_COLLECTION_PREFIXES = (('STEABS', 'ABS'), ('STERAC', 'RAC'), ('STEREL', 'REL'))

#: The highest group id the arena stores beside a kind in one byte, mirrored from the core's
#: `STEREO_GROUP_MAX` (`core/_molecule_arena.pxi`), which a Python layer cannot read from a `DEF`.
#: A file's id above it is renumbered by `_parse_collections`, the id being a label and not data.
_MAX_STEREO_GROUP = 63


def _bonds_keyword(sgroup_type):
    """Which bond keyword a record of this type carries its bond list in.

    ``DAT`` attaches its datum to bonds with ``CBONDS``; a superatom or a polymer names its
    *crossing* bonds with ``XBONDS``.  One stored list serves both, since the type says which it is.
    """
    return 'CBONDS' if sgroup_type == 'DAT' else 'XBONDS'


def _split_kv(token):
    """``KEY=value`` -> ``(KEY, value)`` with the value unquoted; a positional token -> ``(None, t)``."""
    i = token.find('=')
    if i <= 0:
        return None, token
    key = token[:i]
    if not key.replace('_', '').isalnum():
        return None, token
    return key.upper(), _unquote(token[i + 1:])


def _unquote(value):
    if len(value) >= 2 and value[0] == '"' and value[-1] == '"':
        return value[1:-1].replace('""', '"')
    return value


def _int(value, default=0):
    try:
        return int(value)
    except ValueError:
        try:  # a writer that puts 1.0 where the spec says an integer
            return int(float(value))
        except ValueError:
            return default


def _float(value, default=0.0):
    try:
        return float(value)
    except ValueError:
        return default


def _resolve_element(token, index, log):
    """``(element, isotope, label, r_index)`` for an atom-type token.  Raises for a query or template atom.

    `label` is the token itself where it names no element, as in V2000; V3000 has no alias line, so
    this field is the only place such a text can arrive from.
    """
    token = _unquote(token)
    if token in _SYMBOL_SET:
        return token, 0, None, 0
    if token in _H_ISOTOPES:
        return 'H', _H_ISOTOPES[token], None, 0
    upper = token.upper()
    if upper in ('R', 'R#', '*') or (upper.startswith('R') and upper[1:].isdigit()):
        # `R#` takes its group from an `RGROUPS=` keyword on the same line; absent one it stays 0.
        # `*` is the same attachment point spelt differently.  Upper-cased for the same reason
        # `_SYMBOL_FOLD` exists: real files case-fold whole records.  `RB`, `RU` and `RN` reach
        # this line only as themselves, since `_SYMBOL_SET` and `_H_ISOTOPES` are tested first and no
        # element is `R` followed by digits.
        if token != upper:
            log.append(LogRecord('v3000:atom-type-folded', (),
                                 f'atom {index}: atom type {token!r} read as {upper!r}', REPAIRED))
        group = int(upper[1:]) if upper[1:].isdigit() else 0
        if group > R_INDEX_MAX:
            log.append(LogRecord('v3000:r-index-too-wide', (),
                                 f'atom {index}: R index {group} is past R_INDEX_MAX ({R_INDEX_MAX}), '
                                 f'so the marker is left unindexed', LOST))
            group = 0
        return 'R', 0, None, group
    if token in _QUERY_TYPES:
        raise UnsupportedCtfile(f'atom {index} is a query atom ({_QUERY_TYPES[token]}); this is a '
                                f'query CTAB and chython reads structure CTABs only')
    if token.startswith('[') or upper.startswith('NOT'):
        raise UnsupportedCtfile(f'atom {index} is an atom list; this is a query CTAB and chython '
                                f'reads structure CTABs only')
    folded = _SYMBOL_FOLD.get(upper)
    if folded is not None:
        log.append(LogRecord('v3000:element-folded', (), f'atom {index}: element {token!r} read as {folded}', REPAIRED))
        return folded, 0, None, 0
    log.append(LogRecord('v3000:atom-type-as-label', (),
                         f'atom {index}: atom type {token!r} names no element, so it is read as the '
                         f'display label it is: the atom is kept as the marker {LABEL_ELEMENT!r} '
                         f'carrying {token!r} as its alias. Whatever the label abbreviates is not in '
                         f'the structure', LOST))
    return LABEL_ELEMENT, 0, token, 0


def parse_v3000(lines, log=None):
    """Parse a whole V3000 molfile -- four header lines then the ``M  V30`` block -- into a `Ctab`.

    `lines` is the record's physical lines with line endings already removed or not; either is fine.
    """
    log = [] if log is None else log
    if len(lines) < 4:
        raise MalformedCtfile(f'molfile has {len(lines)} lines, fewer than the four header lines')

    body = []
    for line in lines[4:]:
        if line.startswith('M  END'):
            break
        body.append(line)
    else:
        # Not fatal: a record truncated at the end of a file usually still has a complete CTAB.
        log.append(LogRecord('v3000:no-m-end', (), 'no M  END; CTAB read to end of record'))

    ctab = parse_ctab_block(join_continuations(body, log), log)
    # `Ctab.build` starts its own log from this one, so without this line every V3000 recovery is
    # invisible to a caller reading the build log.
    ctab.log = log
    ctab.title = lines[0].rstrip()
    ctab.program = lines[1].rstrip()
    ctab.comment = lines[2].rstrip()
    # Columns 20-22 of the program line are the dimensionality stamp.  Informational only: a non-zero
    # z is what makes a structure 3D, since writers emitting 3D under a "2D" stamp are common.
    ctab.dimensionality = lines[1][20:22].strip() if len(lines[1]) > 20 else ''
    return ctab


def parse_ctab_block(logical, log=None):
    """Parse joined, prefix-stripped V3000 logical lines into a `Ctab`."""
    log = [] if log is None else log
    ctab = Ctab()
    index_of = {}      # the file's atom index -> position in ctab.atoms
    bond_of = {}       # the file's bond index -> (position a, position b)
    sgroup_lines = []
    collection_lines = []
    declared = None

    i = 0
    n = len(logical)
    while i < n:
        line = logical[i].strip()
        i += 1
        if not line:
            continue
        upper = line.upper()
        if upper.startswith('BEGIN CTAB'):
            continue
        if upper.startswith('END CTAB'):
            break
        if upper.startswith('COUNTS'):
            declared = tokenize(line)[1:]
            continue
        if upper.startswith('BEGIN '):
            block = upper[6:].split()[0]
            body = []
            while i < n:
                inner = logical[i].strip()
                i += 1
                if inner.upper().startswith('END ' + block):
                    break
                body.append(inner)
            else:
                log.append(LogRecord('v3000:block-unclosed', (), f'{block} block not closed; read to end of CTAB'))
            if block == 'ATOM':
                _parse_atoms(body, ctab, index_of, log)
            elif block == 'BOND':
                _parse_bonds(body, ctab, index_of, bond_of, log)
            elif block == 'SGROUP':
                sgroup_lines = body
            elif block == 'COLLECTION':
                collection_lines = body
            elif block in ('TEMPLATE', 'OBJ3D', 'RGROUP'):
                log.append(LogRecord(
                    'v3000:unsupported-block', (), f'unsupported: {block} block ignored ({len(body)} line(s))', LOST))
            else:
                log.append(LogRecord(
                    'v3000:unknown-block', (), f'unsupported: unknown block {block} ignored ({len(body)} line(s))',
                    LOST))
            continue
        if upper.startswith('LINKNODE'):
            log.append(LogRecord('v3000:linknode', (), 'unsupported: LINKNODE ignored', LOST))
            continue
        log.append(LogRecord('v3000:unrecognised-line', (), f'unrecognised V3000 line ignored: {line!r:.60}', LOST))
    else:
        # No END CTAB.  Everything before the truncation parsed, so the CTAB is still usable.
        log.append(LogRecord('v3000:no-end-ctab', (), 'END CTAB missing; CTAB read to end of input'))

    # S-groups and collections last: both reference atoms and bonds by index, so they need the maps
    # complete whatever order the file put its blocks in.  Files with SGROUP before ATOM exist.
    if sgroup_lines:
        _parse_sgroups(sgroup_lines, ctab, index_of, bond_of, log)
        # After the S-groups exist and not before: a stated hydrogen count arrives as one of them.
        apply_mrv_implicit_h(ctab, log)
    if collection_lines:
        _parse_collections(collection_lines, ctab, index_of, log)

    if declared:
        _check_counts(declared, ctab, log)
    if not ctab.atoms:
        raise MalformedCtfile('CTAB has no atom block')
    return ctab


def _check_counts(declared, ctab, log):
    """Compare the COUNTS line with what the blocks actually held.

    The blocks win: a wrong count is a writer bug and the blocks are the data, and a reader trusting
    the count truncates a file whose atom block is longer than it claims.
    """
    if len(declared) >= 1:
        na = _int(declared[0], -1)
        if na >= 0 and na != len(ctab.atoms):
            log.append(LogRecord(
                'v3000:counts-mismatch', (), f'COUNTS says {na} atoms, atom block has {len(ctab.atoms)}'))
    if len(declared) >= 2:
        nb = _int(declared[1], -1)
        if nb >= 0 and nb != len(ctab.bonds):
            log.append(LogRecord(
                'v3000:counts-mismatch', (), f'COUNTS says {nb} bonds, bond block has {len(ctab.bonds)}'))
    if len(declared) >= 5:
        ctab.chiral = _int(declared[4]) == 1


def _parse_atoms(body, ctab, index_of, log):
    for line in body:
        tokens = tokenize(line, log)
        if len(tokens) < 6:
            log.append(LogRecord(
                'v3000:atom-short', (), f'atom line has {len(tokens)} field(s), needs 6: {line!r:.60}'))
            if len(tokens) < 2:
                continue
        file_index = _int(tokens[0], -1)
        if file_index < 0:
            log.append(LogRecord(
                'v3000:atom-bad-index', (), f'atom line with unreadable index ignored: {line!r:.60}', LOST))
            continue
        if file_index in index_of:
            log.append(LogRecord(
                'v3000:atom-repeated', (), f'atom index {file_index} repeated; the second is ignored', LOST))
            continue
        symbol, isotope, label, r_index = _resolve_element(tokens[1], file_index, log)
        atom = CtabAtom(symbol)
        atom.isotope = isotope
        atom.label = label
        atom.r_index = r_index
        atom.file_index = file_index
        if len(tokens) > 4:
            atom.x = _float(tokens[2])
            atom.y = _float(tokens[3])
            atom.z = _float(tokens[4])
        if len(tokens) > 5:
            atom.map_number = _int(tokens[5])

        unknown = []
        for token in tokens[6:]:
            key, value = _split_kv(token)
            if key is None:
                log.append(LogRecord('v3000:atom-positional-field', (),
                                     f'atom {file_index}: positional field {value!r:.20} after the sixth '
                                     f'ignored', LOST))
                continue
            if key == 'CHG':
                atom.charge = _int(value)
            elif key == 'RAD':
                rad = _int(value)
                # The core stores one radical bit.  A doublet is that bit; a singlet or triplet
                # diradical is two unpaired electrons on one atom, so it is set and reported.
                atom.radical = rad != 0
                if rad in (1, 3):
                    kind = 'singlet' if rad == 1 else 'triplet'
                    log.append(LogRecord('v3000:rad-diradical', (),
                                         f'atom {file_index}: RAD={rad} ({kind} diradical) stored as a '
                                         f'single radical', REPAIRED))
                elif rad not in (0, 2):
                    log.append(LogRecord(
                        'v3000:rad-out-of-range', (),
                        f'atom {file_index}: RAD={rad} out of range, stored as a radical', REPAIRED))
            elif key == 'MASS':
                atom.isotope = _int(value)
            elif key == 'CFG':
                atom.parity = _int(value)
            elif key == 'VAL':
                val = _int(value)
                # VAL=-1 is the spec's spelling of "valence zero", which is a real statement about,
                # say, a bare metal ion, and is not the same as VAL absent.
                atom.valence = 0 if val == -1 else val
            elif key == 'HCOUNT':
                log.append(LogRecord('v3000:hcount', (),
                                     f'unsupported: atom {file_index}: HCOUNT={value} is a query field, ignored; '
                                     f'hydrogens computed from valence rules', LOST))
            elif key == 'RGROUPS':
                # `(N i1 ... iN)`, decoded by the same helper `ATOMS=` uses.  One atom carries one
                # group here, so a member of several keeps the first and says so.
                members = parse_list(value, log)
                if len(members) > 1:
                    log.append(LogRecord('v3000:rgroups-multiple', (),
                                         f'atom {file_index}: RGROUPS names {len(members)} groups; the '
                                         f'first is kept', LOST))
                group = _int(members[0], -1) if members else -1
                if 0 <= group <= R_INDEX_MAX:
                    if atom.r_index and atom.r_index != group:
                        # A record that states an index twice.  The property line is the format's override,
                        # so it wins, and a reader saying which one it took costs one line.
                        log.append(LogRecord('v3000:rgroups-overrides-type-token', (),
                                             f'atom {file_index}: the type token reads R{atom.r_index} and '
                                             f'RGROUPS names group {group}; the property line wins',
                                             REPAIRED))
                    atom.r_index = group
                else:
                    log.append(LogRecord('v3000:rgroups-out-of-range', (),
                                         f'atom {file_index}: RGROUPS value {(members[0] if members else "")!r} '
                                         f'is not an R index in 0-{R_INDEX_MAX}; the marker is left '
                                         f'unindexed', LOST))
            else:
                unknown.append(key)
        if unknown:
            log.append(LogRecord(
                'v3000:atom-unknown-keywords', (),
                f'unsupported: atom {file_index}: ignored keyword(s) {", ".join(sorted(set(unknown)))}', LOST))
        index_of[file_index] = len(ctab.atoms)
        ctab.atoms.append(atom)


def _parse_bonds(body, ctab, index_of, bond_of, log):
    for line in body:
        tokens = tokenize(line, log)
        if len(tokens) < 4:
            log.append(LogRecord(
                'v3000:bond-short', (), f'bond line has {len(tokens)} field(s), needs 4: {line!r:.60}'))
            continue
        file_index = _int(tokens[0], -1)
        type_ = _int(tokens[1], 1)
        a = _int(tokens[2], -1)
        b = _int(tokens[3], -1)
        if a not in index_of or b not in index_of:
            log.append(LogRecord('v3000:bond-unknown-atom', (),
                                 f'bond {file_index} references unknown atom index '
                                 f'{a if a not in index_of else b}, dropped', LOST))
            continue
        # The same translation V2000 uses; putting the file's number straight into `order` turns a
        # vendor's coordination bond (`M  V30 1 9 3 7`) into a single one.  The endpoint check stays
        # above it: a bond naming an atom that does not exist is dropped whatever its type.
        order = order_from_bond_type(type_, f'bond {file_index}', log)
        bond = CtabBond(index_of[a], index_of[b], order)
        for token in tokens[4:]:
            key, value = _split_kv(token)
            if key is None:
                continue
            if key == 'CFG':
                cfg = _int(value)
                if cfg in WEDGE_FROM_V3000:
                    bond.wedge = WEDGE_FROM_V3000[cfg]
                else:
                    log.append(LogRecord(
                        'v3000:bond-bad-cfg', (), f'bond {file_index}: CFG={cfg} not a V3000 wedge code, ignored',
                        LOST))
            elif key == 'TOPO':
                bond.topology = _int(value)
            elif key == 'RXCTR':
                bond.reacting_center = _int(value)
            elif key not in _BOND_MODELLED:
                log.append(LogRecord(
                    'v3000:bond-unknown-keyword', (), f'unsupported: bond {file_index}: ignored keyword {key}', LOST))
        bond_of[file_index] = (bond.a, bond.b)
        ctab.bonds.append(bond)


def _parse_sgroups(body, ctab, index_of, bond_of, log):
    for line in body:
        tokens = tokenize(line, log)
        if len(tokens) < 2:
            log.append(LogRecord('v3000:sgroup-short', (), f'sgroup line too short, ignored: {line!r:.60}', LOST))
            continue
        index = _int(tokens[0], NO_INDEX)
        stype = _unquote(tokens[1]).upper()
        if len(tokens) > 2 and '=' not in tokens[2]:
            # V3000 states both numbers as unbounded integers, so both can arrive outside what the
            # model holds.  The external one is checked here, where the token's presence is still
            # known; `index` needs the whole Ctab, so `normalize_indices` does it after the loop.
            ext = checked_index(_int(tokens[2], NO_INDEX), f'sgroup {index} external number', log)
        else:
            ext = NO_INDEX
        sg = SGroup(stype, index, ext)
        bonds_key = _bonds_keyword(stype)
        start = 3 if (len(tokens) > 2 and '=' not in tokens[2]) else 2

        for token in tokens[start:]:
            key, value = _split_kv(token)
            if key is None:
                sg.log.append(LogRecord(
                    'v3000:sgroup-positional-field', (), f'positional field {value!r:.20} ignored', LOST))
                continue
            if key in _SGROUP_ATOM_KEYS:
                target = sg.atoms if key == 'ATOMS' else sg.patoms
                for item in parse_list(value, sg.log):
                    idx = _int(item, -1)
                    if idx in index_of:
                        target.append(index_of[idx])
                    else:
                        sg.log.append(LogRecord(
                            'v3000:sgroup-atom-ref', (), f'{key} references unknown atom index {item}', LOST))
            elif key in _SGROUP_BOND_KEYS:
                if key != bonds_key:
                    # A DAT with XBONDS, or a superatom with CBONDS.  Merging it into the modelled bond
                    # list would emit it under the other keyword, and keeping it verbatim would leave
                    # bond indices the writer renumbers (see `_INDEX_VALUED`).
                    sg.log.append(LogRecord('v3000:sgroup-wrong-bond-key', (),
                                            f'{UNSUPPORTED}{key} on a {stype} group is not modelled; '
                                            f'it will not be re-emitted', LOST))
                    continue
                for item in parse_list(value, sg.log):
                    idx = _int(item, -1)
                    if idx in bond_of:
                        sg.bonds.append(bond_of[idx])
                    else:
                        sg.log.append(LogRecord(
                            'v3000:sgroup-bond-ref', (), f'{key} references unknown bond index {item}', LOST))
            elif key == 'CSTATE':
                # `CSTATE=(4 <bond> x y z)`: the leading value is a bond index, which the writer
                # renumbers, so it is split into (endpoint pair, vector tail).  An unresolvable index
                # keeps the whole value as text.
                items = parse_list(value, sg.log)
                idx = _int(items[0], -1) if items else -1
                if idx in bond_of:
                    sg.cstates.append((bond_of[idx], ' '.join(items[1:])))
                else:
                    sg.cstates.append((None, ' '.join(items)))
                    sg.log.append(LogRecord('v3000:sgroup-cstate-ref', (),
                                            f'CSTATE references unknown bond index '
                                            f'{items[0] if items else "(empty)"}, kept verbatim', LOST))
            elif key == 'FIELDNAME':
                sg.name = value
            elif key == 'FIELDDATA':
                # Bytes, not text: latin-1 round-trips every byte, so this decodes nothing and loses
                # nothing.  `SGroup.field_data` is where a caller asks for a string.
                sg.data.append(value.encode('latin-1', errors='backslashreplace'))
            elif key == 'FIELDDISP':
                sg.disp = parse_fielddisp(value, sg.log)
                if sg.disp is None:
                    sg.fields.setdefault(key, []).append(value)
            elif key == 'PARENT':
                sg.parent = checked_index(_int(value, NO_INDEX), f'sgroup {index} PARENT', log)
            elif key == 'SUBTYPE':
                sg.subtype = value
            elif key in _INDEX_VALUED:
                # Dropped rather than kept: the value names atoms or bonds by the file's numbering,
                # which the writer regenerates, so keeping it verbatim would point at other atoms.
                # Losing a superatom's attachment point is a declared gap; moving it is a wrong file.
                sg.log.append(LogRecord(
                    'v3000:sgroup-index-valued', (),
                    f'{UNSUPPORTED}{key} references {_INDEX_VALUED[key]} indices and is not modelled; '
                    f'it will not be re-emitted', LOST))
            else:
                sg.fields.setdefault(key, []).append(value)

        merge_log(log, sg.log, f'sgroup {index} {stype}: ')
        ctab.sgroups.append(sg)
    normalize_indices(ctab.sgroups, log)


def _parse_collections(body, ctab, index_of, log):
    """Read the collection block into ``ctab.groups`` as ``position -> (kind, group)``.

    **A group id is a label, so an out-of-range one is renumbered rather than dropped.** What a
    collection states is which atoms share a group and of which kind; the number naming it carries
    nothing further, which is why the stored id is opaque (`core/_stereo.pxi`, ruling F79). The arena
    holds 1..63 beside the kind in one byte and files exceed it -- ``MDLV30/STERAC1384`` occurs in the
    wild -- so the id is mapped to a free one of its own kind, in file order, consistently for
    every line naming it. Only a record already holding 63 groups of that kind has nothing free left.
    """
    parsed = []
    taken = {}                                  # kind -> the in-range ids the file itself states
    for line in body:
        tokens = tokenize(line, log)
        if not tokens:
            continue
        name = _unquote(tokens[0]).upper()
        if not name.startswith('MDLV30/'):
            log.append(LogRecord('v3000:collection-unknown', (), f'unsupported: collection {name!r:.30} ignored', LOST))
            continue
        tag = name[7:]
        matched = next((p for p in _COLLECTION_PREFIXES if tag.startswith(p[0])), None)
        if matched is None:
            log.append(LogRecord('v3000:collection-unknown', (), f'unsupported: collection MDLV30/{tag} ignored', LOST))
            continue
        prefix, suffix = matched
        kind = STEREO_FROM_COLLECTION[suffix]
        if suffix == 'ABS':
            group = 0  # ABS is one bucket, not a numbered group
        else:
            group = _int(tag[len(prefix):], 0)
            if group < 1:
                log.append(LogRecord(
                    'v3000:collection-no-group', (), f'collection MDLV30/{tag} has no group number, read as group 1',
                    REPAIRED))
                group = 1
            if group <= _MAX_STEREO_GROUP:
                taken.setdefault(kind, set()).add(group)
        positions = []
        for token in tokens[1:]:
            key, value = _split_kv(token)
            if key != 'ATOMS':
                continue
            for item in parse_list(value, log):
                idx = _int(item, -1)
                if idx in index_of:
                    positions.append(index_of[idx])
                else:
                    log.append(LogRecord(
                        'v3000:collection-atom-ref', (),
                        f'collection MDLV30/{tag} references unknown atom index {item}', LOST))
        parsed.append((kind, group, tag, positions))

    renumbered = {}
    for kind, group, tag, positions in parsed:
        if group > _MAX_STEREO_GROUP:
            key = (kind, group)
            if key not in renumbered:
                free = taken.setdefault(kind, set())
                new = next((i for i in range(1, _MAX_STEREO_GROUP + 1) if i not in free), None)
                if new is None:
                    log.append(LogRecord(
                        'v3000:collection-group-dropped', (),
                        f'collection MDLV30/{tag} group {group} is outside 1..{_MAX_STEREO_GROUP} and all '
                        f'{_MAX_STEREO_GROUP} ids of its kind are taken, dropped', LOST))
                    renumbered[key] = None
                else:
                    log.append(LogRecord(
                        'v3000:collection-group-renumbered', (),
                        f'collection MDLV30/{tag} group {group} is outside 1..{_MAX_STEREO_GROUP}, '
                        f'renumbered to {new}', REPAIRED))
                    free.add(new)
                    renumbered[key] = new
            group = renumbered[key]
            if group is None:
                continue
        for position in positions:
            ctab.groups[position] = (kind, group)


def _coord(value):
    """A coordinate in V3000's free format.

    Trailing zeros are trimmed, as reference writers do, which also keeps the line short enough to
    avoid a continuation; the value is unchanged.
    """
    text = f'{value:.4f}'.rstrip('0').rstrip('.')
    return text if text and text != '-0' else '0'


def emit_v3000(mol, sgroups=None, *, title=None, program='', comment='',
               log=None):
    """Render `mol` as V3000 molfile lines (no line endings).  Returns ``(lines, log)``.

    One extension to the spec is emitted and nothing else: the ``MRV_IMPLICIT_H`` data S-group, for an
    atom whose hydrogen count the valence rules would not reproduce.  In particular the atom ``CFG``
    keyword is left out although several tools write it -- the spec makes bond ``CFG``, the wedge,
    what defines stereo and atom ``CFG`` informational, and an informational field computed from a
    frame-relative parity can disagree with the wedges beside it.

    Bond orders are written as stored, aromatic order 4 as bond type 4, exactly as in V2000: type 4 is
    how a CTfile spells an aromatic bond and this reader reads it back.  Same census as V2000: of the
    five toolkits measured in ``chython/formats/test/``, three write type 4 for an aromatic ring by
    default and all five read it back as one.  Call ``kekule()`` first for an alternating file, which
    is what the other two write.

    An aromatic bond does change ``VAL=``, which is not written for an atom holding one: ``VAL`` is a
    total valence, and reaching it from a face-value aromatic order needs the valence tables and the
    aromatic classifier, which a format module has not.  The hydrogen count is then carried by the
    ``MRV_IMPLICIT_H`` group alone.  See ``valence_for_write``.
    """
    out = [] if log is None else log
    # `None` for either means the caller did not say, so the molecule answers -- see `resolve_output`.
    title, sgroups = resolve_output(mol, title, sgroups, log=out)
    sids = list(mol.atom_numbers)
    position = {sid: i + 1 for i, sid in enumerate(sids)}

    wedges, _ = wedges_for_write(mol, out)
    wedge_of = {}
    for narrow, wide, code in wedges:
        wedge_of[(narrow, wide)] = code

    bonds = list(mol.bonds())
    groups = mol.canonical_stereo_groups() if mol.has_stereo_groups else {}
    # A hydrogen count the valence rules would not reproduce is stated as an MRV_IMPLICIT_H data
    # S-group, and in `VAL=` too when every bond on the atom has an integral order.  Re-derived rather
    # than passed through, so repeated round trips accumulate nothing.  Counted in COUNTS, hence here.
    implicit_records, implicit_log = implicit_h_records(mol, sgroups)
    out.extend(implicit_log)
    records = [r for r in (sgroups if sgroups is not None else ())
               if not (r.is_data() and r.name == MRV_IMPLICIT_H)]
    records += implicit_records

    # The chiral flag says the whole structure is one known enantiomer, which any AND or OR collection
    # contradicts, so the two are never both asserted.
    relative = any(kind != 1 for kind, _ in groups)
    configured = any(mol.parity_of(sid) for sid in sids)
    chiral = 1 if configured and not relative else 0

    lines = [title[:80], program[:80], comment[:80],
             f'  0  0  0  0  0  0            999 {V3000_STAMP}']
    body = [f'COUNTS {len(sids)} {len(bonds)} {len(records)} 0 {chiral}']

    has_xyz = mol.has_3d
    has_xy = mol.has_coordinates
    body.append('BEGIN ATOM')
    for sid in sids:
        atom = mol.atom(sid)
        # As in `emit_v2000`: the conformer wins over the depiction, one atom line holding one triple.
        if has_xyz:
            x, y, zc = mol.xyz_of(sid)
        else:
            x, y = mol.xy_of(sid) if has_xy else (0.0, 0.0)
            zc = 0.0
        fields = [str(position[sid]), atom.atomic_symbol, _coord(x), _coord(y), _coord(zc),
                  str(atom.map_number)]
        if atom.charge:
            fields.append(f'CHG={atom.charge}')
        if atom.is_radical:
            fields.append('RAD=2')
        if atom.isotope:
            fields.append(f'MASS={atom.isotope}')
        # No `VAL=` for an atom holding `H_UNKNOWN` -- the field implies a definite hydrogen count --
        # nor for one holding an aromatic bond, whose total valence is not this layer's to compute.
        valence = valence_for_write(mol, sid)
        if valence is not None:
            # V3000 spells a stated zero valence -1, not 15.  Emitted only when the valence tables
            # would give this atom a different hydrogen count than it is carrying.
            fields.append(f'VAL={-1 if valence == ZERO_VALENCE else valence}')
        # `R#` with the group in `RGROUPS=` is the format's spelling for an indexed marker.
        if atom.is_r and atom.r_index:
            fields[1] = 'R#'
            fields.append(f'RGROUPS=(1 {atom.r_index})')
        body.append(' '.join(fields))
    body.append('END ATOM')

    body.append('BEGIN BOND')
    bond_position = {}
    for i, bond in enumerate(bonds, 1):
        # `a`, `b` are the endpoints in file order, which the wedge may reverse; the bond's own
        # endpoints stay `bond.n`, `bond.m`, which is what `bond_position` is keyed by.  `bond.order`
        # goes straight into the field, as in V2000: the stored orders and the format's bond types are
        # the same numbers, so a mixed-representation molecule writes a block of mixed types.
        a, b, code = wedge_in_file_order(wedge_of, bond.n, bond.m)
        fields = [str(i), str(bond.order), str(position[a]), str(position[b])]
        if code:
            fields.append(f'CFG={WEDGE_TO_V3000[code]}')
        body.append(' '.join(fields))
        bond_position[(bond.n, bond.m)] = i
        bond_position[(bond.m, bond.n)] = i
    body.append('END BOND')

    if records:
        body.append('BEGIN SGROUP')
        # An S-group index is regenerated as the record's position here, so `PARENT` is translated like
        # any other index: a file numbering its groups 1, 2, 5 would else name group 5 of three.
        renumber = {record.index: i for i, record in enumerate(records, 1)
                    if record.index != NO_INDEX}
        for i, record in enumerate(records, 1):
            body.append(_emit_sgroup(record, i, position, bond_position, out, renumber))
        body.append('END SGROUP')

    if groups:
        body.append('BEGIN COLLECTION')
        for (kind, group), members in sorted(groups.items()):
            tag = STEREO_TO_COLLECTION[kind]
            if kind != 1:
                tag = f'{tag}{group}'
            listed = ' '.join(str(position[m]) for m in sorted(members) if m in position)
            body.append(f'MDLV30/{tag} ATOMS=({len(members)} {listed})')
        body.append('END COLLECTION')

    body.insert(0, 'BEGIN CTAB')
    body.append('END CTAB')
    for content in body:
        lines.extend(emit_v30(content))
    lines.append('M  END')
    return lines, out


def _emit_sgroup(record, index, position, bond_position, log, renumber=None):
    """One S-group as a logical V3000 line, in a fixed keyword order.

    The order is canonical rather than as-read, so two files stating the same S-groups produce
    byte-identical output; values themselves are preserved exactly.

    `renumber` maps each record's own S-group index to the number it is being written under, which is
    what `PARENT` is translated through.  A parent naming a group no longer being written is dropped
    and reported rather than invented.
    """
    ext = record.ext_index if record.ext_index != NO_INDEX else index
    fields = [str(index), record.type, str(ext)]
    if record.parent != NO_INDEX:
        parent = (renumber or {}).get(record.parent, record.parent if renumber is None else None)
        if parent is None:
            log.append(LogRecord('v3000:sgroup-parent-lost', (),
                                 f'sgroup {index} {record.type}: PARENT={record.parent} names a group that is '
                                 f'not being written, reference dropped', LOST))
        else:
            fields.append(f'PARENT={parent}')
    if record.subtype:
        fields.append(f'SUBTYPE={quote_value(record.subtype)}')

    atoms = [position[a] for a in record.atoms if a in position]
    if len(atoms) != len(record.atoms):
        log.append(LogRecord('v3000:sgroup-atom-lost', (),
                             f'sgroup {index} {record.type}: '
                             f'{len(record.atoms) - len(atoms)} atom reference(s) no longer in the molecule', LOST))
    if atoms:
        fields.append(f'ATOMS=({len(atoms)} {" ".join(str(a) for a in atoms)})')
    patoms = [position[a] for a in record.patoms if a in position]
    if patoms:
        fields.append(f'PATOMS=({len(patoms)} {" ".join(str(a) for a in patoms)})')

    if record.bonds:
        # A stored bond reference is an endpoint pair, and a pair can name two live atoms with no bond
        # between them: the arena keeps atom references correct across an edit but cannot see
        # `delete_bond`.
        numbers = []
        for pair in record.bonds:
            number = bond_position.get(pair)
            if number is None:
                log.append(LogRecord('v3000:sgroup-bond-lost', (),
                                     f'sgroup {index} {record.type}: bond {pair[0]}-{pair[1]} no longer '
                                     f'exists, reference dropped', LOST))
            else:
                numbers.append(number)
        if numbers:
            key = _bonds_keyword(record.type)
            fields.append(f'{key}=({len(numbers)} {" ".join(str(x) for x in sorted(numbers))})')

    for pair, tail in record.cstates:
        tail_items = tail.split()
        if pair is None:
            # Never resolved to a bond, so it stays text, with its own leading count added back.
            values = tail_items
        else:
            number = bond_position.get(pair)
            if number is None:
                log.append(LogRecord('v3000:sgroup-cstate-lost', (),
                                     f'sgroup {index} {record.type}: CSTATE bond {pair[0]}-{pair[1]} no '
                                     f'longer exists, state dropped', LOST))
                continue
            values = [str(number)] + tail_items
        fields.append(f'CSTATE=({len(values)} {" ".join(values)})')

    if record.name:
        fields.append(f'FIELDNAME={quote_value(record.name)}')
    if record.disp is not None:
        fields.append(f'FIELDDISP={quote_value(format_fielddisp(record.disp, log))}')
    for datum in record.data:
        fields.append(f'FIELDDATA={quote_value(datum.decode("latin-1"))}')
    for key in sorted(record.fields):
        for value in record.fields[key]:
            fields.append(f'{key}={quote_value(value)}')
    return ' '.join(fields)
