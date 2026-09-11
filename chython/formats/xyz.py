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
"""XYZ format reader: atoms and coordinates out of a string, as :class:`XYZFrame` objects.

The format states no bonds, so ``xyz(text)`` returns frames and not a
:class:`~chython.core.MoleculeContainer`; perception is a separate, explicitly invoked pass.  Every
frame in the string comes back -- a trajectory is not read first-frame-only.  Coordinates are
Angstroms as written, all three of them.  Log prefixes: ``atom``, ``record``, ``unsupported``.
"""

from __future__ import annotations

from re import findall, IGNORECASE

from ._text import require_text
from ..core._core import element_symbols
from ..core._log import LOST, LogRecord, REPAIRED


__all__ = ['XYZAtom', 'XYZFrame', 'build_molecule', 'xyz', 'xyz_conformers']


_SYMBOLS = element_symbols()          # ('R', 'H', 'He', ..., 'Og')
_VALID = frozenset(_SYMBOLS[1:])      # all 118 element symbols, proper case
_SYM_UPPER = {s.upper(): s for s in _VALID}   # 'CL' -> 'Cl', 'FE' -> 'Fe', etc.


class XYZAtom:
    """One atom line from an XYZ file, with x, y, z in Angstroms as written.

    ``element`` is the normalized symbol; an unrecognized token is kept raw for the perception pass to
    deal with.  ``isotope`` is 0 unless the symbol implied one -- ``D`` and ``T`` normalize to ``H``
    with isotope 2 and 3.
    """
    __slots__ = ('element', 'isotope', 'x', 'y', 'z')

    def __init__(self, element: str, isotope: int, x: float, y: float, z: float):
        self.element = element
        self.isotope = isotope
        self.x = x
        self.y = y
        self.z = z

    def __repr__(self):
        iso = f'/{self.isotope}' if self.isotope else ''
        return f'XYZAtom({self.element}{iso}, {self.x:.4f}, {self.y:.4f}, {self.z:.4f})'

    def __eq__(self, other):
        if not isinstance(other, XYZAtom):
            return NotImplemented
        return (self.element == other.element and self.isotope == other.isotope
                and self.x == other.x and self.y == other.y and self.z == other.z)

    def __hash__(self):
        # By value, matching __eq__; safe because the parser fills these fields once.
        return hash((self.element, self.isotope, self.x, self.y, self.z))


class XYZFrame:
    """One frame from an XYZ file: atoms, 3D coordinates, and parse log.

    ``stated_count`` is what the count line said and is never revised to match what was found;
    ``len(atoms)`` is how many atom lines parsed.  ``title`` is the comment line verbatim, which may
    carry ``charge=N``/``radical=N`` pairs from older chython tools or extended-XYZ
    ``Properties=``/``Lattice=`` pairs; both are stored unchanged and the extended fields are also
    reported as ``unsupported:``.  ``log`` holds this frame's messages, and holds them because a frame is
    not a container: XYZ states no bond, so there is no molecule to write to until a later pass builds
    one.  A reader that returns a :class:`~chython.core.MoleculeContainer` writes to ``mol.log`` instead.
    """
    __slots__ = ('title', 'atoms', 'log', 'stated_count')

    def __init__(self):
        self.title: str = ''
        self.atoms: list = []
        self.log: list = []
        self.stated_count: int = 0

    def __len__(self):
        return len(self.atoms)

    def __repr__(self):
        return (f'XYZFrame(stated={self.stated_count}, atoms={len(self.atoms)}, '
                f'title={self.title!r:.30})')


def _normalize_element(raw: str) -> tuple:
    """``(symbol, isotope, LogRecord_or_None)``, the record ready to append.

    The chain, in order: direct match; case fold ('CL' -> 'Cl'); bare integer as an atomic number;
    'D'/'T' as hydrogen with isotope 2/3; trailing non-alphabetic characters stripped and the first two
    steps retried; anything else stored raw with an ``atom:`` message.
    """
    # 1. Direct match
    if raw in _VALID:
        return raw, 0, None

    # 2. Case normalization
    upper = raw.upper()
    sym = _SYM_UPPER.get(upper)
    if sym is not None:
        return sym, 0, LogRecord('xyz:symbol-case-corrected', (),
                                 f'atom: element symbol {raw!r} corrected to {sym!r}',
                                 REPAIRED)

    # 3. Bare atomic number ('6' → 'C')
    if raw.isdigit():
        z = int(raw)
        if 1 <= z <= 118:
            sym = _SYMBOLS[z]
            return sym, 0, LogRecord('xyz:atomic-number-to-symbol', (),
                                     f'atom: atomic number {z} stored as element {sym!r}',
                                     REPAIRED)
        return raw, 0, LogRecord('xyz:atomic-number-out-of-range', (),
                                 f'atom: atomic number {raw!r} out of range 1-118',
                                 LOST)

    # 4. Hydrogen isotope aliases
    if raw == 'D':
        return 'H', 2, LogRecord('xyz:deuterium-to-hydrogen', (),
                                 "atom: 'D' treated as hydrogen (deuterium, isotope 2)",
                                 REPAIRED)
    if raw == 'T':
        return 'H', 3, LogRecord('xyz:tritium-to-hydrogen', (),
                                 "atom: 'T' treated as hydrogen (tritium, isotope 3)",
                                 REPAIRED)

    # 5. Trailing non-alphabetic characters stripped ('C1' → 'C', 'O-' → 'O')
    prefix = ''
    for ch in raw:
        if ch.isalpha():
            prefix += ch
        else:
            break
    if prefix and prefix != raw:
        if prefix in _VALID:
            return prefix, 0, LogRecord('xyz:symbol-trailing-stripped', (),
                                        f'atom: trailing characters stripped from {raw!r},'
                                        f' stored as {prefix!r}',
                                        REPAIRED)
        upper_prefix = prefix.upper()
        sym2 = _SYM_UPPER.get(upper_prefix)
        if sym2 is not None:
            return sym2, 0, LogRecord('xyz:symbol-trailing-stripped', (),
                                      f'atom: trailing characters stripped and case corrected '
                                      f'from {raw!r}, stored as {sym2!r}',
                                      REPAIRED)

    # 6. Unknown
    return raw, 0, LogRecord('xyz:symbol-unrecognized', (),
                             f'atom: unrecognized element symbol {raw!r}',
                             LOST)


def _looks_like_atom_line(line: str) -> bool:
    """True when *line* looks like an XYZ atom line: 4+ tokens whose 2nd, 3rd and 4th are floats.

    A count line is one bare integer and an atom line is at least four tokens with three parseable
    coordinates, so the two shapes cannot be confused and a one-line lookahead is exact.  Extra
    trailing columns do not matter here, as they do not to the atom parser.
    """
    parts = line.split()
    if len(parts) < 4:
        return False
    try:
        float(parts[1])
        float(parts[2])
        float(parts[3])
        return True
    except ValueError:
        return False


def _looks_like_truncated_atom_line(line: str) -> bool:
    """True when *line* is an atom line that stops mid-coordinate: 2-3 tokens, the rest floats.

    The shape a killed job leaves on its last line, recognised past the stated count as well as inside
    it so the last atom of a half-written trajectory is reported.  Still unambiguous against a count
    line, which is a single token.
    """
    parts = line.split()
    if not 2 <= len(parts) <= 3:
        return False
    try:
        for part in parts[1:]:
            float(part)
    except ValueError:
        return False
    return True


def _parse_count(line: str) -> int | None:
    """The atom count *line* states, or ``None`` when it does not state one.

    A count line is a bare non-negative integer: ``'3 atoms'``, ``'3.0'``, ``'-3'`` and a count with a
    trailing comment are not, and no part of them is guessed at, a wrong guess inventing a frame
    boundary.
    """
    try:
        count = int(line)
    except ValueError:
        return None
    return count if count >= 0 else None


def _parse_extended_xyz(title: str, frame_log: list, caller_log: list) -> None:
    """Detect and log extended-XYZ (ASE/libAtoms) comment fields.

    The dialect carries key=value pairs in the comment line and is recognised by ``Properties=`` or
    ``Lattice=``; every field gets an ``unsupported:`` line.  Atom lines are read as symbol x y z
    regardless, so a non-standard ``Properties`` column order -- position before species -- is read
    wrongly and says so.
    """
    # key=value, the value quoted or unquoted non-whitespace
    pairs = findall(r'([A-Za-z_][A-Za-z0-9_]*)=("(?:[^"\\]|\\.)*"|[^\s]+)', title)
    if not pairs:
        return  # plain comment, nothing to log

    keys = {k for k, _ in pairs}
    if 'Properties' not in keys and 'Lattice' not in keys:
        return  # key=value in comment but not the extended-XYZ dialect

    for key, value in pairs:
        if key == 'Properties':
            # Standard layout: species:S:1:pos:R:3, possibly with more columns.
            v = value.strip('"')
            cols = v.split(':')
            if len(cols) < 6 or cols[0].lower() != 'species' or cols[3].lower() != 'pos':
                msg = (f'unsupported: extended XYZ Properties column order {value!r} not stored; '
                       f'atom lines parsed as symbol x y z anyway')
                frame_log.append(LogRecord('xyz:extended-xyz-field', (), msg, LOST))
                caller_log.append(LogRecord('xyz:extended-xyz-field', (), msg, LOST))
            else:
                # Standard order, so the atom lines are read correctly; the metadata is still lost.
                msg = f'unsupported: extended XYZ field {key!r} not stored'
                frame_log.append(LogRecord('xyz:extended-xyz-field', (), msg, LOST))
                caller_log.append(LogRecord('xyz:extended-xyz-field', (), msg, LOST))
        else:
            msg = f'unsupported: extended XYZ field {key!r} not stored'
            frame_log.append(LogRecord('xyz:extended-xyz-field', (), msg, LOST))
            caller_log.append(LogRecord('xyz:extended-xyz-field', (), msg, LOST))


def xyz(text: str, *, log: list | None = None) -> list:
    """Parse all XYZ frames from *text* and return a list of :class:`XYZFrame` objects.

    Each frame holds the symbols and 3D coordinates the file stated; no bond is perceived here.  *log*
    is an optional caller-supplied list receiving every frame's messages flat.

    Nothing raises on chemically wrong or malformed input.  *text* must be a ``str`` -- ``bytes`` raises
    ``TypeError`` at the call site rather than guessing an encoding.
    """
    text = require_text(text, 'xyz')
    log = [] if log is None else log
    frames: list = []
    lines = text.splitlines()
    n = len(lines)
    i = 0

    while i < n:
        line = lines[i].strip()

        # A blank separator between frames is common in trajectory writers.
        if not line:
            i += 1
            continue

        # No count line, no frame: everything up to the next one goes unread, with a line saying how
        # many atom lines a bad header ('3 atoms', '3.0', a BOM glued to the digits) cost.  The count is
        # not reconstructed -- that would be inventing a frame.
        count = _parse_count(line)
        if count is None:
            lost = 1 if _looks_like_atom_line(line) or _looks_like_truncated_atom_line(line) else 0
            i += 1
            while i < n:
                candidate = lines[i].strip()
                if candidate:
                    if _parse_count(candidate) is not None:
                        break  # the next frame starts here; leave it for the outer loop
                    if _looks_like_atom_line(candidate) or _looks_like_truncated_atom_line(candidate):
                        lost += 1
                i += 1
            log.append(LogRecord('xyz:no-count-line', (),
                                 f'record: {line!r} does not state an atom count; the block it starts is not '
                                 f'read as a frame ({lost} atom line(s) discarded)',
                                 LOST))
            continue

        # Comment line: always the line immediately after the count, even when blank.
        i += 1
        title = lines[i] if i < n else ''
        i += 1

        frame = XYZFrame()
        frame.stated_count = count
        frame.title = title

        # Exactly `count` non-blank lines, stopping early at the next frame's count line, which is the
        # only frame boundary this format has.
        atoms_found = 0
        while atoms_found < count and i < n:
            atom_line = lines[i].strip()

            # A blank line inside the block consumes no slot.
            if not atom_line:
                i += 1
                continue

            # A count line here is the next frame's, not a malformed atom.  One definition of "count
            # line" for the whole reader, so a shape rejected as a header cannot end a frame either.
            if _parse_count(atom_line) is not None:
                break  # do not consume the next frame's count

            i += 1
            atoms_found += 1

            parts = atom_line.split()
            if len(parts) < 4:
                msg = f'atom: fewer than 4 fields on atom line {atoms_found}: {atom_line!r}'
                frame.log.append(LogRecord('xyz:atom-line-too-short', (), msg, LOST))
                log.append(LogRecord('xyz:atom-line-too-short', (), msg, LOST))
                continue

            raw_sym = parts[0]
            try:
                x, y, z_coord = float(parts[1]), float(parts[2]), float(parts[3])
            except ValueError:
                msg = f'atom: non-numeric coordinate on atom line {atoms_found}: {atom_line!r}'
                frame.log.append(LogRecord('xyz:non-numeric-coordinate', (), msg, LOST))
                log.append(LogRecord('xyz:non-numeric-coordinate', (), msg, LOST))
                continue

            sym, isotope, norm_msg = _normalize_element(raw_sym)
            if norm_msg is not None:
                frame.log.append(LogRecord(*norm_msg))
                log.append(LogRecord(*norm_msg))

            frame.atoms.append(XYZAtom(sym, isotope, x, y, z_coord))

        # Over-stated count: the truncation a job killed mid-write leaves behind.
        if atoms_found != count:
            msg = (f'record: count stated {count} atoms; '
                   f'{atoms_found} found before end of frame')
            frame.log.append(LogRecord('xyz:count-overstated', (), msg, LOST))
            log.append(LogRecord('xyz:count-overstated', (), msg, LOST))

        # Under-stated count: only the stated N atoms are kept, since the stated count is what
        # resynchronises every later frame, and the surplus lines are counted and reported.  Both
        # atom-line shapes are counted; neither can be confused with the next frame's header.
        elif i < n:
            surplus = 0
            truncated = []
            while i < n:
                peek = lines[i].strip()
                if not peek:
                    i += 1
                    continue
                if _looks_like_atom_line(peek):
                    surplus += 1
                elif _looks_like_truncated_atom_line(peek):
                    surplus += 1
                    truncated.append(peek)
                else:
                    break
                i += 1
            # Reported lines are consumed, or the outer loop would meet them again and report the same
            # lines a second time under the other reading.
            if surplus:
                msg = (f'record: count stated {count} atoms; '
                       f'{surplus} surplus atom line(s) not stored')
                frame.log.append(LogRecord('xyz:count-understated', (), msg, LOST))
                log.append(LogRecord('xyz:count-understated', (), msg, LOST))
                for peek in truncated:
                    # Same damage as inside the count, so the same prefix and wording.
                    msg = f'atom: fewer than 4 fields on surplus atom line: {peek!r}'
                    frame.log.append(LogRecord('xyz:atom-line-too-short', (), msg, LOST))
                    log.append(LogRecord('xyz:atom-line-too-short', (), msg, LOST))

        # After atom parsing, so the log reads atom issues first and metadata second.
        _parse_extended_xyz(title, frame.log, log)

        frames.append(frame)

    return frames


def xyz_conformers(molecule, frames, *, log: list | None = None) -> int:
    """Store each :class:`XYZFrame` as one conformer of `molecule`; return how many landed.

    A frame is a state of a molecule the caller already has, so the topology comes from a MOL, an SDF,
    a SMILES, :func:`build_molecule` or :func:`chython.formats.pdb.build_molecule`.  Atoms are matched
    POSITIONALLY against ``molecule.atom_numbers``, and a frame whose length or element sequence
    disagrees is logged and skipped while the rest still land.  The return value is what a caller
    compares against ``len(frames)`` without reading the log.

    Frames APPEND: a molecule carrying no conformer takes the first frame as model 0, and one already
    carrying models keeps them and gains the frames after them.  Each conformer's ``ext_index`` is the
    frame's ordinal in `frames`, stored verbatim.
    """
    log = [] if log is None else log
    own: list = []
    numbers = molecule.atom_numbers
    stored = _store_frames(molecule, frames, numbers,
                           [_SYMBOLS[molecule.element_of(n)] for n in numbers], own)
    molecule.log.absorb('read', own)
    log.extend(own)
    return stored


def _store_frames(molecule, frames, numbers: list, symbols: list, own: list) -> int:
    """Each frame as one appended conformer, positionally; how many landed.  Records into `own`.

    `symbols` is what each position is expected to state, which is not always the molecule's own
    element: :func:`build_molecule` stores an unreadable symbol as the R marker and still has to match
    the frames against the symbol the file wrote.
    """
    stored = 0
    for ordinal, frame in enumerate(frames):
        if len(frame.atoms) != len(numbers):
            own.append(LogRecord('xyz:frame-atom-count', (),
                                 f'record: frame {ordinal} holds {len(frame.atoms)} atom(s) where the '
                                 f'molecule holds {len(numbers)}, so it is not stored', LOST))
            continue
        mismatch = next((i for i, (a, s) in enumerate(zip(frame.atoms, symbols)) if a.element != s), -1)
        if mismatch >= 0:
            own.append(LogRecord('xyz:frame-element-mismatch', (),
                                 f'record: frame {ordinal} states {frame.atoms[mismatch].element!r} at '
                                 f'position {mismatch} where the molecule holds {symbols[mismatch]!r}, '
                                 f'so it is not stored', LOST))
            continue
        try:
            # ALL-OR-NOTHING PER FRAME, and the edit scope is what enforces it: leaving the `with` by
            # exception discards the journal, so a coordinate the container refuses drops the whole
            # frame rather than half of one.
            with molecule.edit():
                model = molecule.add_conformer(ext_index=ordinal)
                for n, atom in zip(numbers, frame.atoms):
                    molecule.set_xyz(n, atom.x, atom.y, atom.z, model=model)
        except ValueError as e:
            own.append(LogRecord('xyz:frame-refused', (),
                                 f'record: frame {ordinal} holds a coordinate the container refuses '
                                 f'({e}), so it is not stored', LOST))
            continue
        stored += 1
    return stored


def build_molecule(frames, *, log: list | None = None):
    """One or more :class:`XYZFrame` objects as a :class:`~chython.core.MoleculeContainer`.

    `frames` is a frame or a sequence of them.  The FIRST states the atoms; every frame states one
    model, in the order read, each carrying its ordinal as its ``ext_index``.  A frame whose atom count
    or element sequence disagrees with the first is logged and skipped, the way :func:`xyz_conformers`
    treats one that disagrees with its molecule.

    WHAT THIS DOES NOT DO IS PERCEIVE.  The molecule arrives with atoms, coordinates and NO bond, so
    two explicit calls follow it: ``chython.chemistry.perceive_bonds()`` for the connectivity the
    distances imply, then ``chython.chemistry.saturate()`` for the orders the hydrogen counts force.

    Every atom states ZERO implicit hydrogens, which is the format's own statement: an XYZ record lists
    every atom, hydrogens included, so a hydrogen not written is a hydrogen not there.  A symbol the
    reader could not resolve becomes the R marker rather than a dropped row -- the coordinate is a fact
    the file stated and the element is the part that is missing.  ``charge=``/``radical=`` in the
    comment line is reported and not applied: it names no atom, and nothing here guesses which one it
    meant.

    *log* is an optional list receiving this pass's own findings; the molecule's ``log`` receives those
    and the reader's findings for every frame stored on it.
    """
    from ..core._core import MoleculeContainer

    log = [] if log is None else log
    if isinstance(frames, XYZFrame):
        frames = [frames]
    else:
        frames = list(frames)
    own: list = []
    molecule = MoleculeContainer()
    if not frames:
        molecule.log.absorb('read', own)
        return molecule

    symbols = [atom.element for atom in frames[0].atoms]
    numbers = []
    for position, atom in enumerate(frames[0].atoms):
        if atom.element in _VALID:
            numbers.append(molecule.add_atom(atom.element, isotope=atom.isotope, implicit_h=0))
            continue
        # The reader already said the symbol is not an element; this says what became of the atom.
        own.append(LogRecord('xyz:symbol-not-an-element', (),
                             f'atom: {atom.element!r} at position {position} is not an element, so the '
                             'atom is stored as the R marker and keeps its coordinate', LOST))
        numbers.append(molecule.add_atom('R', implicit_h=0))

    for frame in frames:
        if 'charge=' in frame.title or 'radical=' in frame.title:
            own.append(LogRecord('xyz:title-charge-not-applied', (),
                                 f'record: comment line {frame.title!r} states a charge or a radical '
                                 'count for the record as a whole; it names no atom and is not '
                                 'applied', LOST))

    _store_frames(molecule, frames, numbers, symbols, own)
    molecule.log.absorb('read', [record for frame in frames for record in frame.log] + own)
    log.extend(own)
    return molecule
