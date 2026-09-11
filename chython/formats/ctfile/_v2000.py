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
"""MDL V2000: the fixed-column CTAB and the ``M  `` properties block.

Fields are identified by column, so a short line returns a wrong-but-plausible substring: every
field is read through :func:`_field`, and the positions live in named constants.  The blocks outrank
the counts line, and the properties block outranks the atom block -- ``M  CHG`` supersedes ``ccc``
and ``M  ISO`` supersedes ``dd``, per the spec -- so it is applied after the atom block is read.
"""

from ._ctab import (Ctab, CtabAtom, CtabBond, order_from_bond_type, LABEL_ELEMENT,
                    WEDGE_FROM_V2000, WEDGE_TO_V2000)
from ._errors import MalformedCtfile, UnsupportedCtfile
from ._hydrogens import (MRV_IMPLICIT_H, ZERO_VALENCE, apply_mrv_implicit_h, implicit_h_records,
                         valence_for_write)
from ._sgroup import (NO_INDEX, SGroup, SGroupStore, checked_index, format_fielddisp,
                      merge_log, normalize_indices, parse_fielddisp, resolve_output)
from ...core.wedge import wedge_in_file_order, wedges_for_write
from ...core import INFO, LogRecord, LOST, R_INDEX_MAX, REPAIRED, WEDGE_EITHER, WEDGE_NONE
from ...core._core import element_symbols


__all__ = ['parse_v2000', 'emit_v2000', 'V2000_STAMP']


V2000_STAMP = 'V2000'
_SYMBOLS = element_symbols()

# The counts line.  `lll` (atom lists), `sss` (stext) and `mmm` (the always-999 properties count)
# describe blocks V2000 has not used since the 1990s, or restate what the block itself says.
_COUNTS_ATOMS = (0, 3)
_COUNTS_BONDS = (3, 6)
_COUNTS_CHIRAL = (12, 15)
_COUNTS_VERSION = (33, 39)

# The atom line.  Named because `line[36:39]` at the point of use is unreviewable.
_ATOM_X = (0, 10)
_ATOM_Y = (10, 20)
_ATOM_Z = (20, 30)
_ATOM_SYMBOL = (31, 34)
_ATOM_MASS_DIFF = (34, 36)
_ATOM_CHARGE = (36, 39)
_ATOM_PARITY = (39, 42)
_ATOM_HCOUNT = (42, 45)      # a query field: "n-1 or more hydrogens".  Read, logged, not applied.
_ATOM_STEREO_CARE = (45, 48)
_ATOM_VALENCE = (48, 51)
_ATOM_MAP = (60, 63)

# The bond line.
_BOND_A = (0, 3)
_BOND_B = (3, 6)
_BOND_TYPE = (6, 9)
_BOND_STEREO = (9, 12)
_BOND_TOPOLOGY = (15, 18)
_BOND_CENTER = (18, 21)

# `ccc` is a code, not a number: 1..3 are +3..+1, 5..7 are -1..-3, and 4 is not a charge at all but
# a doublet radical, so it sets the radical flag.
_CHARGE_CODES = {0: (0, False), 1: (3, False), 2: (2, False), 3: (1, False), 4: (0, True),
                 5: (-1, False), 6: (-2, False), 7: (-3, False)}
_CHARGE_TO_CODE = {0: 0, 3: 1, 2: 2, 1: 3, -1: 5, -2: 6, -3: 7}

# Bond types are not declared here: `order_from_bond_type` in `_ctab` is the one translation both
# versions call, and the type table lives in its docstring.

# Query atom symbols.  `A` and `Q` are any-atom shorthands, `L` heads an atom list.  Each is a
# valid CTfile statement that a molecule cannot hold, so each is refused by name rather than guessed
# at.  `R`, `R#`, `R<n>` and `*` are NOT query symbols here -- they are the marker element 0.
_QUERY_SYMBOLS = {'A': 'any atom', 'Q': 'any heteroatom', 'L': 'an atom list'}

_SGROUP_TEXT_LABEL = {'MUL': 'MULT'}  # everything else spells its SMT text LABEL= in V3000

# V2000 writes an S-group number in three columns, so this format's ceiling is 999, below the 65534
# the record model holds (V3000's).  Reachable: read a V3000 file numbering its groups 1..1500 and
# write it as V2000.  A wider number would push every following column right.
_SGROUP_NUMBER_MAX = 999


def _field(line, span):
    """The text of `span` in `line`, or ``''`` when the line stops short.

    A truncated trailing field is the common damage in a column format, and ``''`` reads as
    "not stated" where a plain slice would give a wrong-but-plausible substring.
    """
    start, end = span
    if start >= len(line):
        return ''
    return line[start:end]


def _int(text, default=0):
    text = text.strip()
    if not text:
        return default
    try:
        return int(text)
    except ValueError:
        try:  # `1.0` in an integer column; a float formatter reached a field it should not have
            return int(float(text))
        except ValueError:
            return default


def _float(text, default=None):
    text = text.strip()
    if not text:
        return default
    try:
        return float(text)
    except ValueError:
        return default


def _resolve_element(symbol, log, line_number):
    """``(element, isotope, label, r_index)`` for a V2000 atom symbol, or a refusal.

    Recovers the case-folding that fixed-column writers produce (`CL`, `br`) with a log line, since
    the symbol column is left-justified text and some writers upper-case the whole record.

    `label` is the text itself for the last case, where it names no element: the reference lists the
    tokens this field takes and says nothing about free text, and implementations differ over what a
    label there means.  A record is not lost over it -- see :data:`LABEL_ELEMENT`.
    """
    symbol = symbol.strip()
    if not symbol:
        raise MalformedCtfile(f'atom line {line_number}: no element symbol')
    upper = symbol.upper()
    if upper in ('R', 'R#', '*') or (upper.startswith('R') and upper[1:].isdigit()):
        # `R#` takes its index from an `M  RGP` line later in the record; absent one it stays 0, an
        # unindexed marker rather than a refusal.  `*` is the same attachment point spelt differently.
        # Upper-cased because this column is fixed-width text and writers case-fold whole records; the
        # fold is logged the same way `cl` and `br` are, and `RB`, `RU` and `RN` fall through to the
        # element path since no element is `R` followed by digits.
        if symbol != upper:
            log.append(LogRecord('v2000:symbol-case-folded', (),
                                 f'atom line {line_number}: symbol {symbol!r} read as {upper!r}',
                                 REPAIRED))
        index = int(upper[1:]) if upper[1:].isdigit() else 0
        if index > R_INDEX_MAX:
            log.append(LogRecord('v2000:r-index-too-wide', (),
                                 f'atom line {line_number}: R index {index} is past R_INDEX_MAX '
                                 f'({R_INDEX_MAX}), so the marker is left unindexed', LOST))
            index = 0
        return 'R', 0, None, index
    if symbol in _QUERY_SYMBOLS:
        raise UnsupportedCtfile(f'atom line {line_number}: {symbol!r} is {_QUERY_SYMBOLS[symbol]}, '
                                f'which a molecule cannot represent. Read this file with a query '
                                f'reader')
    if symbol == 'D':
        log.append(LogRecord('v2000:d-as-hydrogen', (),
                             f'atom line {line_number}: D read as hydrogen isotope 2', REPAIRED))
        return 'H', 2, None, 0
    if symbol == 'T':
        log.append(LogRecord('v2000:t-as-hydrogen', (),
                             f'atom line {line_number}: T read as hydrogen isotope 3', REPAIRED))
        return 'H', 3, None, 0
    if symbol in _SYMBOLS:
        return symbol, 0, None, 0
    folded = symbol.capitalize()
    if folded in _SYMBOLS:
        log.append(LogRecord('v2000:symbol-case-folded', (),
                             f'atom line {line_number}: symbol {symbol!r} read as {folded!r}',
                             REPAIRED))
        return folded, 0, None, 0
    log.append(LogRecord('v2000:symbol-as-label', (),
                         f'atom line {line_number}: symbol {symbol!r} names no element, so it is read '
                         f'as the display label it is: the atom is kept as the marker '
                         f'{LABEL_ELEMENT!r} carrying {symbol!r} as its alias. Whatever the label '
                         f'abbreviates is not in the structure', LOST))
    return LABEL_ELEMENT, 0, symbol, 0


def parse_v2000(lines, log=None):
    """Parse a V2000 record into a :class:`~chython.formats.ctfile._ctab.Ctab`.

    `lines` is the whole record: title, program, comment, counts, then the blocks.  A missing
    trailing ``M  END`` is reported rather than fatal.
    """
    out = [] if log is None else log
    if len(lines) < 4:
        raise MalformedCtfile(f'V2000 record has {len(lines)} lines; the header alone needs 4')

    ctab = Ctab()
    ctab.log = out
    ctab.title = lines[0].rstrip()
    ctab.program = lines[1].rstrip()
    ctab.comment = lines[2].rstrip()
    # Columns 20-22 of the program line, which is the only place a V2000 record says whether its
    # coordinates are a 2D drawing or a 3D structure.  Absent in most real files.
    ctab.dimensionality = _field(lines[1], (20, 22)).strip()

    counts = lines[3]
    declared_atoms = _int(_field(counts, _COUNTS_ATOMS), -1)
    declared_bonds = _int(_field(counts, _COUNTS_BONDS), -1)
    ctab.chiral = _int(_field(counts, _COUNTS_CHIRAL)) == 1
    if declared_atoms < 0:
        raise MalformedCtfile(f'V2000 counts line: no atom count in {counts!r:.40}')
    if declared_bonds < 0:
        out.append(LogRecord('v2000:missing-bond-count', (),
                             f'V2000 counts line: no bond count in {counts!r:.40}, read as 0',
                             REPAIRED))
        declared_bonds = 0

    # The counts line is a hint, not a contract: it is wrong in real files, and a reader trusting it
    # reads the first bond line as an atom.  So each block is bounded by the declared count *and*
    # ended by the shape of the line -- an atom line runs to at least column 31 before its symbol, a
    # bond line is 21 columns.  A full-width atom line with a blank symbol still raises.
    i = 4
    end = len(lines)
    while i < end and len(ctab.atoms) < declared_atoms:
        line = lines[i]
        if _ends_a_block(line):
            break
        ctab.atoms.append(_parse_atom(line, len(ctab.atoms) + 1, out))
        i += 1
    if len(ctab.atoms) < declared_atoms:
        out.append(LogRecord('v2000:atom-count-mismatch', (),
                             f'counts line declares {declared_atoms} atoms but the atom block ends after '
                             f'{len(ctab.atoms)}; read what is there. Every number this record states about an '
                             f'atom past {len(ctab.atoms)} -- bonds, S-groups, collections -- refers to an atom '
                             f'that is not here and will be dropped',
                             REPAIRED))

    bond_lines = 0
    while i < end and bond_lines < declared_bonds:
        line = lines[i]
        if _ends_a_block(line, atom_line=False):
            break
        bond_lines += 1
        bond = _parse_bond(line, bond_lines, len(ctab.atoms), out)
        if bond is not None:
            ctab.bonds.append(bond)
        i += 1
    if bond_lines < declared_bonds:
        out.append(LogRecord('v2000:bond-count-mismatch', (),
                             f'counts line declares {declared_bonds} bonds but the bond block ends after '
                             f'{bond_lines}; read what is there',
                             REPAIRED))

    _parse_properties(lines[i:], ctab, out)
    return ctab


#: Column an atom line has reached by the time its symbol starts.  A bond line is 21 columns wide, so
#: both are recognisable by width alone.
_ATOM_LINE_WIDTH = 31


def _ends_a_block(line, atom_line=True):
    """Whether `line` is the start of the next block rather than another line of this one."""
    if not line.strip():
        return True  # blank line inside a block: nothing follows it that this block can use
    if line[:3] in ('M  ', 'A  ', 'V  ', 'G  ', 'S  ') or line.startswith('$$$$'):
        return True
    return atom_line and len(line.rstrip()) < _ATOM_LINE_WIDTH


def _parse_atom(line, number, log):
    atom = CtabAtom()
    atom.file_index = number
    symbol, isotope, label, r_index = _resolve_element(_field(line, _ATOM_SYMBOL), log, number)
    atom.element = symbol
    atom.isotope = isotope
    atom.label = label
    atom.r_index = r_index

    x = _float(_field(line, _ATOM_X))
    y = _float(_field(line, _ATOM_Y))
    z = _float(_field(line, _ATOM_Z))
    if x is None or y is None:
        log.append(LogRecord('v2000:bad-coordinates', (),
                             f'atom line {number}: unreadable coordinates {_field(line, (0, 30))!r}, '
                             f'placed at the origin',
                             REPAIRED))
        x = y = 0.0
    atom.x = x
    atom.y = y
    atom.z = z or 0.0

    code = _int(_field(line, _ATOM_CHARGE))
    if code in _CHARGE_CODES:
        atom.charge, atom.radical = _CHARGE_CODES[code]
    else:
        log.append(LogRecord('v2000:bad-charge-code', (),
                             f'atom line {number}: charge code {code} is not one of 0-7, read as neutral',
                             REPAIRED))

    atom.mass_diff = _int(_field(line, _ATOM_MASS_DIFF))
    atom.parity = _int(_field(line, _ATOM_PARITY))
    atom.map_number = _int(_field(line, _ATOM_MAP))

    valence = _int(_field(line, _ATOM_VALENCE))
    if valence == ZERO_VALENCE:
        atom.valence = 0
    elif valence:
        atom.valence = valence

    hcount = _int(_field(line, _ATOM_HCOUNT))
    if hcount:
        # `hhh` is a query field: 1 means "0 or more hydrogens", 2 means "1 or more".  Not a count,
        # so reading it as one would disagree with every other tool.  Same as V3000 `HCOUNT=`.
        log.append(LogRecord('v2000:query-hcount', (),
                             f'unsupported: atom line {number}: query hydrogen field hhh={hcount} ignored; it states a '
                             f'minimum, not a count',
                             LOST))
    if _int(_field(line, _ATOM_STEREO_CARE)):
        log.append(LogRecord('v2000:stereo-care', (),
                             f'unsupported: atom line {number}: stereo care box is a query field, ignored',
                             LOST))
    return atom


def _parse_bond(line, number, atom_count, log):
    a = _int(_field(line, _BOND_A), 0)
    b = _int(_field(line, _BOND_B), 0)
    if not a or not b:
        log.append(LogRecord('v2000:bad-bond-atom-numbers', (),
                             f'bond line {number}: unreadable atom numbers {_field(line, (0, 6))!r}, dropped',
                             LOST))
        return None
    if not (1 <= a <= atom_count and 1 <= b <= atom_count):
        log.append(LogRecord('v2000:bond-atom-out-of-range', (),
                             f'bond line {number}: atom number out of 1..{atom_count}, dropped',
                             LOST))
        return None

    order = order_from_bond_type(_int(_field(line, _BOND_TYPE), 1), f'bond line {number}', log)
    bond = CtabBond(a - 1, b - 1, order)
    stereo = _int(_field(line, _BOND_STEREO))
    if stereo in WEDGE_FROM_V2000:
        bond.wedge = WEDGE_FROM_V2000[stereo]
    elif stereo == 3:
        # `3` on a double bond is "cis or trans, either", the counterpart of the `4` that means
        # "up or down, either" on a single bond; the core has one code for both.
        bond.wedge = WEDGE_EITHER
        if order != 2:
            log.append(LogRecord('v2000:cis-or-trans-on-non-double', (),
                                 f'bond line {number}: stereo 3 (cis or trans) on a bond of type {order}'))
    else:
        log.append(LogRecord('v2000:bad-bond-stereo', (),
                             f'bond line {number}: bond stereo {stereo} is not 0, 1, 3, 4 or 6, ignored',
                             LOST))

    bond.topology = _int(_field(line, _BOND_TOPOLOGY))
    bond.reacting_center = _int(_field(line, _BOND_CENTER))
    return bond


def _parse_properties(lines, ctab, log):
    """The ``M  `` / ``A  `` / ``V  `` block, and the S-group assembly that goes with it.

    A V2000 S-group is spread over as many as a dozen lines sharing only an S-group number --
    ``STY`` declares it, ``SAL`` gives it atoms, ``SDT`` its field definition, ``SED``/``SCD`` its
    data, ``SDD`` its display position -- and they may arrive in any order, including data before
    declaration, so a record is created on first mention by whichever line mentions it.
    """
    records = {}          # sgroup number -> SGroup
    data_parts = {}       # sgroup number -> [str], SCD lines awaiting their SED
    atom_count = len(ctab.atoms)
    bond_count = len(ctab.bonds)
    ended = False
    alias_for = None
    abbreviation_for = None   # the atom pair of a `G  aaappp` line whose text line has not arrived
    abbreviations = []        # (atom pair, text) per G line, judged once the S-group records are in

    def record(number):
        try:
            return records[number]
        except KeyError:
            sg = records[number] = SGroup('GEN', number)
            log.append(LogRecord('v2000:sgroup-predeclared', (),
                                 f'sgroup {number} used before M  STY declared it; read as GEN',
                                 REPAIRED))
            return sg

    for line in lines:
        line = line.rstrip('\r\n')
        if alias_for is not None:
            # The line after `A  aaa` is free alias text and may look like a property line.
            ctab.aliases[alias_for] = line.rstrip()
            alias_for = None
            continue
        if abbreviation_for is not None:
            # The line after `G  aaappp` is the abbreviation text and may look like a property line.
            abbreviations.append((abbreviation_for, line.strip()))
            abbreviation_for = None
            continue
        if line.startswith('M  END'):
            ended = True
            break
        if not line.strip():
            continue

        if line.startswith('A  '):
            index = _int(_field(line, (3, 6)), 0)
            if 1 <= index <= atom_count:
                alias_for = index - 1
            else:
                log.append(LogRecord('v2000:alias-atom-out-of-range', (),
                                     f'atom alias for atom {index}, which is out of 1..{atom_count}',
                                     LOST))
                alias_for = None
                continue
            continue
        if line.startswith('G  '):
            # A superatom's display label in the obsolete spelling: `aaa` is an atom of the contracted
            # group and `ppp` the atom outside that the crossing bond reaches, and BOTH are in the atom
            # block already.  So a G line contracts nothing and expands to nothing -- it is the fact a
            # `SUP` record spells with `SAL` + `SMT`, which is where V3 keeps it.  Reading it as an
            # ALIAS would be wrong: the label would land on an atom whose group is already drawn, and
            # `expand_abbreviations` would graft a second copy of it.
            abbreviation_for = (_int(_field(line, (3, 6)), 0), _int(_field(line, (6, 9)), 0))
            continue
        if line.startswith('V  '):
            index = _int(_field(line, (3, 6)), 0)
            if 1 <= index <= atom_count:
                # An atom value: free text attached to an atom, so stored as an alias.
                ctab.aliases[index - 1] = line[7:].rstrip()
            continue
        if not line.startswith('M  '):
            log.append(LogRecord('v2000:unrecognised-props-line', (),
                                 f'unrecognised properties line: {line!r:.60}',
                                 LOST))
            continue

        tag = line[3:6]
        if tag in ('CHG', 'ISO', 'RAD', 'RGP'):
            _parse_atom_property(line, tag, ctab, log)
        elif tag == 'STY':
            for number, value in _pairs(line, log):
                if number in records:
                    records[number].type = value
                else:
                    records[number] = SGroup(value, number)
        elif tag == 'SST':
            for number, value in _pairs(line, log):
                record(number).subtype = value
        elif tag == 'SLB':
            for number, value in _pairs(line, log):
                record(number).ext_index = checked_index(_int(value, NO_INDEX),
                                                         f'sgroup {number} SLB', log)
        elif tag == 'SPL':
            for number, value in _pairs(line, log):
                record(number).parent = checked_index(_int(value, NO_INDEX),
                                                      f'sgroup {number} SPL', log)
        elif tag in ('SAL', 'SPA'):
            sg = record(_int(_field(line, (7, 10)), 0))
            target = sg.atoms if tag == 'SAL' else sg.patoms
            for value in _list(line, log):
                if 1 <= value <= atom_count:
                    target.append(value - 1)
                else:
                    sg.log.append(LogRecord('v2000:sgroup-atom-ref-out-of-range', (),
                                            f'M  {tag} references atom {value}, out of 1..{atom_count}',
                                            LOST))
        elif tag == 'SBL':
            sg = record(_int(_field(line, (7, 10)), 0))
            for value in _list(line, log):
                if 1 <= value <= bond_count:
                    bond = ctab.bonds[value - 1]
                    sg.bonds.append((bond.a, bond.b))
                else:
                    sg.log.append(LogRecord('v2000:sgroup-bond-ref-out-of-range', (),
                                            f'M  SBL references bond {value}, out of 1..{bond_count}',
                                            LOST))
        elif tag == 'SBV':
            # The V2000 spelling of V3000's CSTATE: an S-group bond and a vector along it.  Its first
            # value is a bond number, which is a position, so it becomes an endpoint pair.
            sg = record(_int(_field(line, (7, 10)), 0))
            index = _int(_field(line, (10, 14)), 0)
            tail = line[14:].rstrip()
            if 1 <= index <= bond_count:
                bond = ctab.bonds[index - 1]
                sg.cstates.append(((bond.a, bond.b), ' '.join(tail.split())))
            else:
                sg.cstates.append((None, ' '.join((str(index) + ' ' + tail).split())))
                sg.log.append(LogRecord('v2000:sgroup-bond-ref-out-of-range', (),
                                        f'M  SBV references bond {index}, out of 1..{bond_count}',
                                        LOST))
        elif tag == 'SMT':
            sg = record(_int(_field(line, (7, 10)), 0))
            text = line[11:].rstrip()
            key = _SGROUP_TEXT_LABEL.get(sg.type, 'LABEL')
            sg.fields.setdefault(key, []).append(text)
        elif tag == 'SDT':
            _parse_sdt(line, record(_int(_field(line, (7, 10)), 0)))
        elif tag == 'SDD':
            sg = record(_int(_field(line, (7, 10)), 0))
            # V2000 puts the display coordinates in two F10.4 columns after the S-group number, then
            # styling -- byte-for-byte a V3000 FIELDDISP value, so one model serves both.
            sg.disp = parse_fielddisp(line[11:], sg.log)
            if sg.disp is None:
                sg.fields.setdefault('FIELDDISP', []).append(line[11:].rstrip())
        elif tag in ('SED', 'SCD'):
            number = _int(_field(line, (7, 10)), 0)
            sg = record(number)
            # Trailing spaces are dropped from each part: a writer pads a data line to its full
            # 69-column width, so keeping the padding injects spaces into a continued datum.
            data_parts.setdefault(number, []).append(line[11:].rstrip())
            if tag == 'SED':
                sg.data.append(''.join(data_parts.pop(number)).encode('latin-1',
                                                                      errors='backslashreplace'))
        elif tag == 'ALS':
            raise UnsupportedCtfile('M  ALS declares an atom list, which a molecule cannot '
                                    'represent. Read this file with a query reader')
        elif tag in ('SDS', 'SCN', 'SAP', 'SCL', 'SNC', 'SPS', 'CRS', 'MRV', 'LOG', 'APO',
                     'SBT', 'SGD', 'SMS', 'SPM', 'REG', 'SEQ', 'SUP'):
            # Known S-group and query keywords with no model here.  Named individually rather than
            # wildcarded so the log says which one, and so an unheard-of keyword still reaches the
            # catch-all below.
            log.append(LogRecord('v2000:unmodelled-property', (),
                                 f'unsupported: M  {tag} is not modelled; it will not be re-emitted',
                                 LOST))
        else:
            # A declared fidelity gap: an unrecognised field is normally preserved, but the molecule
            # has no segment for an opaque V2000 property line, so there is nowhere to keep this one.
            # Closing it needs a core segment.  V3000 has no such gap -- an unknown SGROUP keyword
            # rides in `fields`.
            log.append(LogRecord('v2000:unrecognised-property', (),
                                 f'unsupported: unrecognised property M  {tag} dropped, not re-emitted: {line!r:.60}',
                                 LOST))

    if data_parts:
        for number, parts in data_parts.items():
            # SCD lines with no closing SED.  The datum is real and is kept.
            records[number].data.append(''.join(parts).encode('latin-1',
                                                              errors='backslashreplace'))
            log.append(LogRecord('v2000:scd-not-closed', (),
                                 f'sgroup {number}: M  SCD data was not closed by an M  SED line'))
    if not ended:
        log.append(LogRecord('v2000:no-m-end', (),
                             'no M  END; properties block read to the end of the record'))
    for (index, attachment), text in abbreviations:
        covered = any(sg.type == 'SUP' and index - 1 in sg.atoms for sg in records.values())
        log.append(LogRecord('v2000:group-abbreviation', (),
                             f'G  line: {text!r} labels the group at atom {index}, bonded to atom '
                             f'{attachment}' + (', which a SUP S-group already states' if covered else
                                                '; no SUP S-group gives its extent, so the label has '
                                                'no atoms to hold it and is dropped'),
                             INFO if covered else LOST))

    for number in sorted(records):
        sg = records[number]
        merge_log(log, sg.log, f'sgroup {number} {sg.type}: ')
        ctab.sgroups.append(sg)
    normalize_indices(ctab.sgroups, log)
    apply_mrv_implicit_h(ctab, log)


def _parse_atom_property(line, tag, ctab, log):
    """``M  CHG`` / ``M  ISO`` / ``M  RAD`` / ``M  RGP``: a count, then that many (atom, value) pairs.

    These override the atom block, per the spec: ``ccc`` cannot hold a charge past 3 and ``dd``
    cannot hold an absolute mass number at all, so a writer puts a placeholder there and the truth
    here.
    """
    count = _int(_field(line, (6, 9)), -1)
    if count < 0:
        log.append(LogRecord('v2000:missing-entry-count', (),
                             f'M  {tag}: no entry count in {line!r:.40}',
                             LOST))
        return
    for i in range(count):
        base = 10 + i * 8
        index = _int(_field(line, (base, base + 3)), 0)
        value = _int(_field(line, (base + 4, base + 7)), 0)
        if not 1 <= index <= len(ctab.atoms):
            log.append(LogRecord('v2000:atom-ref-out-of-range', (),
                                 f'M  {tag} references atom {index}, out of 1..{len(ctab.atoms)}',
                                 LOST))
            continue
        atom = ctab.atoms[index - 1]
        if tag == 'CHG':
            atom.charge = value
        elif tag == 'ISO':
            atom.isotope = value
        elif tag == 'RGP':
            # The R group number.  Out of the field's range it is dropped rather than truncated: a
            # wrong index names a different fragment, which is worse than an unindexed marker.
            if 0 <= value <= R_INDEX_MAX:
                if atom.r_index and atom.r_index != value:
                    # A record that states an index twice.  The property line is the format's override,
                    # so it wins, and a reader saying which one it took costs one line.
                    log.append(LogRecord('v2000:rgp-overrides-symbol', (),
                                         f'atom {index}: the symbol column reads R{atom.r_index} and '
                                         f'M  RGP names group {value}; the property line wins',
                                         REPAIRED))
                atom.r_index = value
            else:
                log.append(LogRecord('v2000:r-index-too-wide', (),
                                     f'atom {index}: M  RGP {value} is past R_INDEX_MAX '
                                     f'({R_INDEX_MAX}), so the marker is left unindexed', LOST))
        else:
            # `M  RAD` states a multiplicity: 1 singlet (a carbene), 2 doublet, 3 triplet.  The core
            # has one radical bit, so 1 and 3 are recorded as a radical too rather than as closed
            # shell.
            atom.radical = value != 0
            if value in (1, 3):
                atom.radical = True
                log.append(LogRecord('v2000:rad-two-electron', (),
                                     f'atom {index}: M  RAD {value} is a two-electron state; stored as a '
                                     f'single radical, which is the closest the core can hold',
                                     REPAIRED))


def _parse_sdt(line, sg):
    """``M  SDT``: the DAT field definition -- name, type, units, and a query the reader ignores.

    Fixed columns: 30 for the name, 2 for the type, 20 for units or format, then the query operator
    and its data.  Only the name has a model; the rest is kept for the round trip.
    """
    sg.name = line[11:41].strip()
    for key, span in (('FIELDTYPE', (41, 43)), ('FIELDINFO', (43, 63)),
                      ('QUERYTYPE', (63, 65)), ('QUERYOP', (65, 85))):
        value = _field(line, span).strip()
        if value:
            sg.fields.setdefault(key, []).append(value)


def _pairs(line, log):
    """``M  XXX  nn8 aaa vvv aaa vvv ...`` -- a count then that many 8-column pairs."""
    count = _int(_field(line, (6, 9)), -1)
    if count < 0:
        log.append(LogRecord('v2000:missing-entry-count', (),
                             f'{line[:6]}: no entry count in {line!r:.40}',
                             LOST))
        return
    for i in range(count):
        base = 10 + i * 8
        number = _int(_field(line, (base, base + 3)), 0)
        yield number, _field(line, (base + 4, base + 7)).strip()


def _list(line, log):
    """``M  XXX sss nn8 aaa aaa ...`` -- an S-group number, a count, then 4-column entries."""
    count = _int(_field(line, (10, 13)), -1)
    if count < 0:
        log.append(LogRecord('v2000:missing-entry-count', (),
                             f'{line[:6]}: no entry count in {line!r:.40}',
                             LOST))
        return
    for i in range(count):
        base = 14 + i * 4
        yield _int(_field(line, (base, base + 3)), 0)


def _coord(value):
    """A coordinate in the V2000 10-character F10.4 column.

    ``%10.4f`` needs 11 characters at 100000 and at -10000, and one character of overflow shifts
    every following field on the line, which a fixed-column reader misreads as chemistry.  Hence the
    refusal, naming V3000's free-format coordinates as the fix.
    """
    text = f'{value:10.4f}'
    if len(text) > 10:
        raise MalformedCtfile(f'coordinate {value} does not fit the V2000 10-character column; '
                              f'write this structure as V3000')
    return text


def emit_v2000(mol, sgroups=None, *, title=None, program='', comment='', log=None):
    """Render `mol` as a V2000 record: a list of lines with no trailing newlines.

    Spec-conformant output, with one vendor extension.  Three points a writer can get backwards:

    * a charge is written **both** in the atom line's ``ccc`` and in ``M  CHG``, because the spec says
      the properties block wins and many readers look at only one of the two.  A charge past the reach
      of ``ccc`` writes 0 in the column and the truth in ``M  CHG``;
    * bond orders are written as stored: 1, 2 and 3 for a Kekule bond and 4 for an aromatic one, which
      is how a CTfile spells one.  Caveat to weigh when choosing a representation: the spec lists type
      4 among the *query* bond types, so some consumers read it as an aromatic query rather than as a
      delocalised bond.  Call ``kekule()`` first for an alternating file;
    * a hydrogen count the valence rules would not reproduce is written as the ``MRV_IMPLICIT_H`` data
      S-group, and **also** in ``vvv`` where every bond on the atom has an integral order.  An atom
      holding an aromatic bond gets the S-group alone: ``vvv`` is a *total* valence, and reaching it
      from a face-value aromatic order needs the valence tables and the aromatic classifier, which a
      format module has not.  See ``valence_for_write``;
    * an R-atom's index travels in ``M  RGP`` because that is the spelling the format specifies.
      ``atomic_symbol`` answers ``R7``, which would fit the three-character symbol column, but
      ``R#`` plus ``M  RGP`` is the form other readers expect.  An unindexed R writes the bare ``R``.
    """
    out = [] if log is None else log
    # `None` for either means the caller did not say, so the molecule answers -- see `resolve_output`.
    title, sgroups = resolve_output(mol, title, sgroups, log=out)
    sids = list(mol.atom_numbers)
    if len(sids) > 999:
        raise MalformedCtfile(f'{len(sids)} atoms will not fit the V2000 3-character count field; '
                              f'write this structure as V3000')
    position = {sid: i + 1 for i, sid in enumerate(sids)}

    wedges, _ = wedges_for_write(mol, out)
    wedge_of = {(narrow, wide): code for narrow, wide, code in wedges}

    bonds = list(mol.bonds())
    if len(bonds) > 999:
        raise MalformedCtfile(f'{len(bonds)} bonds will not fit the V2000 3-character count field; '
                              f'write this structure as V3000')

    groups = mol.canonical_stereo_groups() if mol.has_stereo_groups else {}
    relative = any(kind != 1 for kind, _ in groups)
    configured = any(mol.parity_of(sid) for sid in sids)
    chiral = 1 if configured and not relative else 0
    if groups and relative:
        # V2000's chiral flag is one bit for the whole record, so an AND or OR collection cannot be
        # written -- and writing the atoms without it would state a single known enantiomer.
        out.append(LogRecord('v2000:enhanced-stereo-not-written', (),
                             'V2000 has no enhanced stereo groups; the AND/OR collections in this structure '
                             'are not written. Write V3000 to keep them',
                             LOST))

    lines = [title[:80], program[:80], comment[:80],
             f'{len(sids):3d}{len(bonds):3d}  0  0{chiral:3d}  0            999 {V2000_STAMP}']

    # The geometry wins over the depiction: a CTAB atom line has one coordinate triple.
    # `SEG_CONFORMERS` is what the file stated, while `SEG_XY` may be a layout engine's projection of
    # it, and writing the projection would flatten a 3D record on every round trip.  A molecule with
    # no conformer writes its `xy` and a z of zero.
    has_xyz = mol.has_3d
    has_xy = mol.has_coordinates
    charges = []
    isotopes = []
    radicals = []
    rgroups = []
    for sid in sids:
        atom = mol.atom(sid)
        if has_xyz:
            x, y, zc = mol.xyz_of(sid)
        else:
            x, y = mol.xy_of(sid) if has_xy else (0.0, 0.0)
            zc = 0.0
        code = _CHARGE_TO_CODE.get(atom.charge, 0)
        # `vvv` 0 means "not stated", which is what an atom holding `H_UNKNOWN` gets, and also one
        # holding an aromatic bond -- whose total valence exists but is not this layer's to compute.
        valence = valence_for_write(mol, sid)
        # `R#` plus an `M  RGP` entry is the format's spelling for an indexed marker.  `atomic_symbol`
        # answers `R7`, which this column would hold -- but the spelling other readers expect is the
        # one with the group number in the properties block.
        if atom.is_r and atom.r_index:
            symbol = 'R#'
            rgroups.append((position[sid], atom.r_index))
        else:
            symbol = atom.atomic_symbol
        lines.append(f'{_coord(x)}{_coord(y)}{_coord(zc)} '
                     f'{symbol:<3s} 0{code:3d}  0  0  0{valence or 0:3d}'
                     f'  0  0  0{atom.map_number:3d}  0  0')
        if atom.charge:
            charges.append((position[sid], atom.charge))
        if atom.isotope:
            isotopes.append((position[sid], atom.isotope))
        if atom.is_radical:
            radicals.append((position[sid], 2))

    bond_position = {}
    for i, bond in enumerate(bonds, 1):
        # `bond.order` goes straight into the column: the stored orders and the format's bond types
        # are the same four numbers for 1, 2, 3 and 4.  A molecule of mixed representation therefore
        # writes a block of mixed types rather than one type chosen for the whole record.
        a, b, code = wedge_in_file_order(wedge_of, bond.n, bond.m)
        stereo = WEDGE_TO_V2000.get(code or WEDGE_NONE, 0)
        if code == WEDGE_EITHER and bond.order == 2:
            stereo = 3  # "cis or trans" rather than "up or down"
        lines.append(f'{position[a]:3d}{position[b]:3d}{bond.order:3d}{stereo:3d}  0  0  0')
        bond_position[(bond.n, bond.m)] = i
        bond_position[(bond.m, bond.n)] = i

    for tag, entries in (('CHG', charges), ('ISO', isotopes), ('RAD', radicals), ('RGP', rgroups)):
        # Eight entries per line is the format's own limit, not a wrapping preference.
        for chunk in (entries[i:i + 8] for i in range(0, len(entries), 8)):
            body = ''.join(f' {n:3d} {v:3d}' for n, v in chunk)
            lines.append(f'M  {tag}{len(chunk):3d}{body}')

    for position_index, text in sorted((position[sid], text)
                                       for sid, text in (sgroups.aliases.items()
                                                         if sgroups else ())
                                       if sid in position):
        lines.append(f'A  {position_index:3d}')
        lines.append(text)

    # The stated-hydrogen S-groups, re-derived and merged with the record's own.  An MRV_IMPLICIT_H
    # that came in with the file is dropped in favour of the fresh one: the stored count is the truth,
    # so re-deriving keeps repeated round trips idempotent.
    implicit_records, implicit_log = implicit_h_records(mol, sgroups)
    out.extend(implicit_log)
    kept = [r for r in (sgroups.records if sgroups else ())
            if not (r.is_data() and r.name == MRV_IMPLICIT_H)]
    if kept or implicit_records:
        merged = SGroupStore(kept + implicit_records)
        lines.extend(_emit_sgroups(merged, position, bond_position, out))
    lines.append('M  END')
    return lines, out


def _emit_sgroups(store, position, bond_position, log):
    """The S-group half of the properties block, one keyword at a time.

    Grouped by keyword, not by S-group -- ``M  STY`` for every record, then every record's ``M  SAL``,
    and so on.  That is what the spec's examples do, and the only order in which the
    8-entries-per-line packing of ``STY``, ``SLB`` and ``SST`` can be filled.
    """
    lines = []
    records = [r for r in store.records if r.type]
    if not records:
        return lines
    # A record keeps the number it came in with when this format can write it, and is given the lowest
    # free one when it cannot -- either because it has none, or because it is too wide for three
    # columns.  Not the record's *position*, which can collide with a number another record states.
    numbers = {}
    taken = {r.index for r in records if 0 <= r.index <= _SGROUP_NUMBER_MAX}
    free = (x for x in range(1, _SGROUP_NUMBER_MAX + 1) if x not in taken)
    writable = []
    for record in records:
        if 0 <= record.index <= _SGROUP_NUMBER_MAX:
            numbers[id(record)] = record.index
        else:
            number = next(free, None)
            if number is None:
                # More than 999 S-groups: no number left to give, and a record with no number cannot
                # be referred to by any of its own lines.
                log.append(LogRecord('v2000:sgroup-count-exceeded', (),
                                     f'sgroup {record.type}: more than {_SGROUP_NUMBER_MAX} S-groups, '
                                     f'record dropped',
                                     LOST))
                continue
            if record.index != NO_INDEX:
                log.append(LogRecord('v2000:sgroup-number-out-of-range', (),
                                     f'sgroup number {record.index} does not fit V2000\'s three columns, '
                                     f'written as {number}',
                                     REPAIRED))
            numbers[id(record)] = number
        if len(record.type) > 3:
            # Truncating changes the type and dropping loses the record, so the wide line is written
            # and the log says so.  V3000 has no such limit.
            log.append(LogRecord('v2000:sgroup-type-too-wide', (),
                                 f'sgroup {numbers[id(record)]}: type {record.type!r} is wider than V2000\'s '
                                 f'three columns, the M  STY line will not be re-readable'))
        writable.append(record)
    records = writable
    if not records:
        return lines
    # What a written number resolves to, for the keywords that name one instead of being one.  Keyed
    # on the stated index, which is the alphabet `parent` speaks.
    renumber = {r.index: numbers[id(r)] for r in records if r.index != NO_INDEX}

    def packed(tag, entries):
        for chunk in (entries[i:i + 8] for i in range(0, len(entries), 8)):
            body = ''.join(f' {n:3d} {v:>3s}' for n, v in chunk)
            lines.append(f'M  {tag}{len(chunk):3d}{body}')

    labels = []
    for r in records:
        if r.ext_index == NO_INDEX:
            continue
        elif 0 <= r.ext_index <= _SGROUP_NUMBER_MAX:
            labels.append((numbers[id(r)], str(r.ext_index)))
        else:  # an external label is a label, so there is no free one to substitute
            log.append(LogRecord('v2000:sgroup-ext-index-out-of-range', (),
                                 f'sgroup {numbers[id(r)]} {r.type}: external number {r.ext_index} does not '
                                 f'fit V2000\'s three columns, dropped',
                                 LOST))

    parents = []
    for r in records:
        if r.parent == NO_INDEX:
            continue
        parent = renumber.get(r.parent)
        if parent is None:
            log.append(LogRecord('v2000:sgroup-parent-not-written', (),
                                 f'sgroup {numbers[id(r)]} {r.type}: parent {r.parent} names a group that is '
                                 f'not being written, reference dropped',
                                 LOST))
        else:
            parents.append((numbers[id(r)], str(parent)))

    packed('STY', [(numbers[id(r)], r.type) for r in records])
    packed('SST', [(numbers[id(r)], r.subtype) for r in records if r.subtype])
    packed('SLB', labels)
    packed('SPL', parents)

    for record in records:
        number = numbers[id(record)]
        for tag, refs in (('SAL', record.atoms), ('SPA', record.patoms)):
            live = [position[a] for a in refs if a in position]
            if len(live) != len(refs):
                log.append(LogRecord('v2000:sgroup-dead-atom-refs', (),
                                     f'sgroup {number} {record.type}: '
                                     f'{len(refs) - len(live)} atom reference(s) no longer in the molecule',
                                     LOST))
            # 15 entries per line, which is what the 4-column fields and an 80-character line allow.
            for chunk in (live[i:i + 15] for i in range(0, len(live), 15)):
                lines.append(f'M  {tag}{number:4d}{len(chunk):3d}'
                             + ''.join(f'{n:4d}' for n in chunk))

        live_bonds = []
        for pair in record.bonds:
            index = bond_position.get(pair)
            if index is None:
                log.append(LogRecord('v2000:sgroup-dead-bond-ref', (),
                                     f'sgroup {number} {record.type}: bond {pair[0]}-{pair[1]} no longer '
                                     f'exists, reference dropped',
                                     LOST))
            else:
                live_bonds.append(index)
        for chunk in (sorted(live_bonds)[i:i + 15] for i in range(0, len(live_bonds), 15)):
            lines.append(f'M  SBL{number:4d}{len(chunk):3d}' + ''.join(f'{n:4d}' for n in chunk))

        for pair, tail in record.cstates:
            if pair is None:
                continue  # never resolved to a bond on the way in; there is no number to write
            index = bond_position.get(pair)
            if index is None:
                log.append(LogRecord('v2000:sgroup-dead-sbv-bond', (),
                                     f'sgroup {number} {record.type}: M  SBV bond {pair[0]}-{pair[1]} no '
                                     f'longer exists, vector dropped',
                                     LOST))
                continue
            lines.append(f'M  SBV{number:4d}{index:4d} {tail}'.rstrip())

        for key in ('LABEL', 'MULT'):
            for text in record.fields.get(key, ()):
                lines.append(f'M  SMT{number:4d} {text}')
        if record.name or any(k in record.fields for k in ('FIELDTYPE', 'FIELDINFO')):
            field_type = _first(record.fields, 'FIELDTYPE')
            info = _first(record.fields, 'FIELDINFO')
            query_type = _first(record.fields, 'QUERYTYPE')
            query_op = _first(record.fields, 'QUERYOP')
            lines.append(f'M  SDT{number:4d} {record.name:<30s}{field_type:<2s}{info:<20s}'
                         f'{query_type:<2s}{query_op:<20s}'.rstrip())
        if record.disp is not None:
            lines.append(f'M  SDD{number:4d} {format_fielddisp(record.disp, log)}')
        for datum in record.data:
            text = datum.decode('latin-1')
            # 69 characters per line is the format's limit; the last part is SED, the rest SCD.
            chunks = [text[i:i + 69] for i in range(0, len(text), 69)] or ['']
            for chunk in chunks[:-1]:
                lines.append(f'M  SCD{number:4d} {chunk}')
            lines.append(f'M  SED{number:4d} {chunks[-1]}')
    return lines


def _first(fields, key):
    values = fields.get(key)
    return values[0] if values else ''
