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
"""The S-group record: a CTfile Sgroup as this library models it.

``DAT`` is fully modelled -- ``FIELDNAME``, ``FIELDDATA``, ``FIELDDISP`` are parsed and re-emitted.
Every other type (``SUP`` ``MUL`` ``SRU`` ``MON`` ``COP`` ``GEN`` ...) is preserved rather than
interpreted: references are translated, every keyword rides in ``fields``.  ``fields`` is keyed by
keyword, so keyword order across the line is *not* preserved -- the writer emits a canonical order.
"""

from itertools import chain

from ...core import LogRecord, LOST, REPAIRED


__all__ = ['SGroup', 'SGroupStore', 'SGROUP_TYPES', 'NO_INDEX', 'DISP_MAX', 'DISP_MIN',
           'resolve_output', 'UNSUPPORTED', 'merge_log',
           'add_data_sgroup', 'data_sgroups', 'FIELDDISP_TAIL']

#: The stable prefix for log lines that report an unmodelled construct.  A caller filters on
#: ``startswith(UNSUPPORTED)``, so writers and the ratchet import this one definition.
UNSUPPORTED = 'unsupported: '


def merge_log(log, sub, prefix):
    """Merge a per-sgroup log into the record's log without burying the `unsupported: ` marker.

    A caller filters on `startswith`, so an entry carrying the marker keeps it in front of the
    location this adds rather than in the middle.
    """
    for entry in sub:
        if str(entry).startswith(UNSUPPORTED):
            log.append(LogRecord('sgroup:merged', (), f'{UNSUPPORTED}{prefix}{str(entry)[len(UNSUPPORTED):]}', LOST))
        else:
            log.append(LogRecord('sgroup:merged', (), f'{prefix}{str(entry)}', getattr(entry, 'severity', 'info')))


# CTfile Sgroup types.  Everything outside this set is still read and preserved -- the set exists
# to decide *display* order on write and to spell the DAT special case, not to gate acceptance.
SGROUP_TYPES = frozenset((
    'SUP', 'MUL', 'SRU', 'MON', 'COP', 'CRO', 'MOD', 'GRA', 'COM', 'MER', 'FOR', 'MIX', 'ANY',
    'GEN', 'DAT',
))

# The sentinel for "no Sgroup number".  0xFFFF rather than 0, because V2000 `M  STY` writes the
# number in a 3-char field where `  0` is representable and nothing in CTfile forbids it.  Also the
# arena's `parent`/`ext_index` sentinel, declared here because this is where the field lives.
NO_INDEX = 0xFFFF

# The largest Sgroup number distinguishable from "no number".  The sentinel is a value inside the
# format's own domain -- a V3000 `Sgroup 65535` is legal, the keyword takes an unbounded integer --
# and the number lives in a uint16 arena slot, so numbers outside 0..65534 are renumbered on read
# and reported; see `normalize_indices`.
INDEX_MAX = NO_INDEX - 1

# The positional order of the string slots, fixed so a fixed-width record can hold them.  Slots 0-4
# are single-valued; slot 5 onward is FIELDDATA, which is the only multi-valued one.  Non-DAT types
# reuse slot 0 for SUBSCRIPT/LABEL and slot 1 for CLASS.
STRING_SLOTS = ('FIELDNAME', 'FIELDINFO', 'FIELDTYPE', 'QUERYTYPE', 'QUERYOP')

# What an F10.4 field can hold.  Asymmetric because the sign costs a column: `99999.9999` and
# `-9999.9999` are ten characters, `-10000.0000` is eleven.  A wider value on read means the file was
# not written to the fixed layout, so the anchor is kept as opaque text; on write one extra character
# makes the *next* field start early, so a `FIELDDISP` of `1e9, 1e9` reparses with a changed y.
DISP_MAX = 99999.9999
DISP_MIN = -9999.9999


class SGroup:
    """One CTfile Sgroup.

    Atom and bond references are held in whatever alphabet the owner is using:

    * inside a :class:`~chython.formats.ctfile._ctab.Ctab` they are **0-based positions** in the
      Ctab's atom and bond lists;
    * inside an :class:`SGroupStore` they are **stable ids** of a built molecule.

    Bonds are held as ``(a, b)`` endpoint pairs in both alphabets, never as bond indices: a bond index
    is positional and does not survive an edit or a reordered write.  The file's 1-based indices are
    translated in and regenerated out, so a round trip preserves the set of referenced bonds but not
    their numbers -- CTfile bond numbers are positions, not identities.
    """
    __slots__ = ('type', 'index', 'ext_index', 'parent', 'subtype', 'atoms', 'patoms', 'bonds',
                 'cstates', 'fields', 'data', 'name', 'disp', 'log')

    def __init__(self, type='DAT', index=NO_INDEX, ext_index=NO_INDEX, parent=NO_INDEX):
        self.type = type
        self.index = index
        self.ext_index = ext_index
        self.parent = parent
        self.subtype = ''
        self.atoms = []          # atom references
        self.patoms = []         # SPA / PATOMS= -- the parent-atom subset of a MUL group
        self.bonds = []          # list of (a, b) endpoint pairs
        # CSTATE, as [((a, b), tail), ...].  Structured rather than verbatim because its first value
        # is a BOND INDEX -- a position in the bond block -- so passing the string through after a
        # renumber would point it at a different bond.  The vector tail rides along as text.
        self.cstates = []
        self.fields = {}         # every keyword this release does not model: key -> list of values
        self.data = []           # FIELDDATA, list[bytes], in file order
        self.name = ''           # FIELDNAME
        self.disp = None         # FIELDDISP parsed as (x, y, rest) or None
        self.log = []

    @property
    def field_data(self):
        """``FIELDDATA`` joined and decoded for display, with undecodable bytes replaced.

        The bytes in :attr:`data` are the truth and are what the writer emits; this accessor is
        explicitly lossy (``errors='replace'``) so a caller can print a datum without deciding an
        encoding.
        """
        return b'\n'.join(self.data).decode('utf8', errors='replace')

    def is_data(self):
        return self.type == 'DAT'

    def translate(self, atom_map, log=None):
        """Return a copy with every atom reference mapped through `atom_map`.

        Used at both boundaries: positions to stable ids after a build, and stable ids to positions
        before a write.  A reference with no image in `atom_map` is dropped and reported, but a record
        whose atom list empties this way stays present.  Bond pairs die when either endpoint dies,
        matching the arena's carry rule.
        """
        out = SGroup(self.type, self.index, self.ext_index, self.parent)
        out.subtype = self.subtype
        out.fields = dict(self.fields)
        out.data = list(self.data)
        out.name = self.name
        out.disp = self.disp
        out.log = list(self.log)

        lost = 0
        for a in self.atoms:
            if a in atom_map:
                out.atoms.append(atom_map[a])
            else:
                lost += 1
        for a in self.patoms:
            if a in atom_map:
                out.patoms.append(atom_map[a])
        for a, b in self.bonds:
            if a in atom_map and b in atom_map:
                out.bonds.append((atom_map[a], atom_map[b]))
            else:
                lost += 1
        for pair, tail in self.cstates:
            if pair is None:
                out.cstates.append((None, tail))  # unresolved on read, kept as text
            elif pair[0] in atom_map and pair[1] in atom_map:
                out.cstates.append(((atom_map[pair[0]], atom_map[pair[1]]), tail))
            else:
                # Demoted to unresolved rather than removed, as the arena requires: the vector tails
                # live in a run whose length derives from the number of cstates, so dropping an entry
                # would shift every later pair onto the tail before it.
                out.cstates.append((None, tail))
                lost += 1
        if lost and log is not None:
            log.append(LogRecord('sgroup:refs-dropped', (),
                                 f'sgroup {self.index} {self.type}: {lost} reference(s) dropped, '
                                 f'atom no longer in the molecule', LOST))
        return out

    def __repr__(self):
        bits = [self.type, f'index={self.index}']
        if self.name:
            bits.append(f'name={self.name!r}')
        if self.data:
            bits.append(f'data={self.field_data!r:.30}')
        bits.append(f'atoms={self.atoms}')
        return f'SGroup({", ".join(bits)})'


class SGroupStore:
    """The S-groups of one molecule, keyed in the molecule's own stable ids.

    The molecule is the storage -- the core holds the records in persistent arena segments -- and this
    class is the view the emitters take and a caller edits, reached through :meth:`to_molecule` and
    :meth:`from_molecule`.

    This class is the encoding boundary: records here hold ``str`` for the keyword-ish fields
    (``type``, ``subtype``, ``name``, ``fields``) and ``bytes`` for ``data``, which a FIELDDATA value
    is not required to be text at all; the arena holds all of them as ``bytes``, never decoding.
    """
    __slots__ = ('records', 'aliases', 'log')

    def __init__(self, records=(), aliases=None, log=None):
        self.records = list(records)
        # Atom aliases (V2000 `A  <n>` lines, MRV mrvAlias): display labels with no place in the
        # chemistry model, kept verbatim rather than modelled.
        self.aliases = dict(aliases) if aliases else {}
        self.log = list(log) if log else []

    def __len__(self):
        return len(self.records)

    def __iter__(self):
        return iter(self.records)

    def __bool__(self):
        return bool(self.records) or bool(self.aliases)

    def data_records(self):
        """Just the DAT records, in file order."""
        return [r for r in self.records if r.is_data()]

    def by_name(self, name):
        """DAT records whose ``FIELDNAME`` is `name`."""
        return [r for r in self.records if r.is_data() and r.name == name]

    def translate(self, atom_map):
        """A new store with every reference mapped, per :meth:`SGroup.translate`."""
        log = []
        records = [r.translate(atom_map, log) for r in self.records]
        aliases = {atom_map[k]: v for k, v in self.aliases.items() if k in atom_map}
        return SGroupStore(records, aliases, chain(self.log, log))

    def to_molecule(self, mol, log=None):
        """Write every record and alias onto `mol`, replacing whatever it held.  Returns `mol`.

        References must already be in `mol`'s stable ids -- :meth:`translate` is what puts them there.

        Two setters, because the core's are narrow: ``set_sgroups`` replaces the records and leaves
        aliases alone, ``set_aliases`` the reverse.  Each is a full blob rebuild, so the empty case is
        skipped rather than paid.

        A dangling ``PARENT`` is dropped and logged, never raised: the core refuses a parent naming an
        index no record carries, and a file stating ``PARENT=3`` with no Sgroup 3 is ordinary damage
        that must not cost the record.  The emitter drops the same reference.
        """
        if self.records or mol.sgroups:
            numbered = {r.index for r in self.records if r.index != NO_INDEX}
            records = []
            for r in self.records:
                d = self._to_dict(r)
                if d['parent'] != NO_INDEX and d['parent'] not in numbered:
                    if log is not None:
                        log.append(LogRecord('sgroup:parent-dropped', (),
                                             f'sgroup {r.index} {r.type}: PARENT={d["parent"]} dropped, no '
                                             f'sgroup in this record carries that number', LOST))
                    d['parent'] = NO_INDEX
                records.append(d)
            mol.set_sgroups(records)
        if self.aliases or mol.aliases:
            mol.set_aliases(self.aliases)
        return mol

    @classmethod
    def from_molecule(cls, mol):
        """Read `mol`'s records and aliases back out as a store, keyed in its stable ids.

        Carries the molecule's own ``sgroup_log`` in: that log is where the core records a reference it
        had to drop across an edit, which cannot be written into the records themselves.
        """
        # Alias text is decoded here, because the core hands aliases back as `bytes` and this store's
        # alphabet is `str`.  Skipping the decode puts a `bytes` into the `A  <n>` line and the join
        # of the whole record fails.
        return cls((cls._from_dict(d) for d in mol.sgroups),
                   {n: t.decode('utf8', errors='replace') if isinstance(t, bytes) else t
                    for n, t in mol.aliases.items()},
                   [x.decode('utf8', errors='replace') if isinstance(x, bytes) else x
                    for x in mol.sgroup_log])

    @staticmethod
    def _to_dict(r):
        """One :class:`SGroup` as the core's record dict.

        ``fields`` flattens: this model keys it by keyword to a list of values, the arena stores a flat
        run of ``(key, value)`` pairs.  Flattening in keyword order and regrouping on the way back
        preserves each keyword's values in order, which is what :class:`SGroup` claims survives.

        EVERY STRING CROSSING INTO THE ARENA IS ENCODED HERE, the CSTATE vector tail included: the
        blob run is bytes, and this side of the boundary is where the model's `str` becomes them.
        """
        disp, tail = (None, b'')
        if r.disp is not None:
            x, y, tail = r.disp
            disp = (x, y)
        return {'type': r.type.encode('utf8'), 'subtype': r.subtype.encode('utf8'),
                'name': r.name.encode('utf8'),
                'disp': disp,
                'disp_tail': tail.encode('utf8') if isinstance(tail, str) else tail,
                'index': r.index, 'ext_index': r.ext_index, 'parent': r.parent,
                'atoms': tuple(r.atoms), 'patoms': tuple(r.patoms), 'bonds': tuple(r.bonds),
                'cstates': tuple((pair, t.encode('utf8') if isinstance(t, str) else t)
                                 for pair, t in r.cstates),
                'data': tuple(r.data),
                'fields': tuple((k.encode('utf8'), v.encode('utf8') if isinstance(v, str) else v)
                                for k, vs in r.fields.items() for v in vs),
                'log': tuple(str(x).encode('utf8') for x in r.log)}

    @staticmethod
    def _from_dict(d):
        """The core's record dict as one :class:`SGroup`.  The inverse of :meth:`_to_dict`."""
        out = SGroup(d['type'].decode('utf8', errors='replace'), d['index'], d['ext_index'],
                     d['parent'])
        out.subtype = d['subtype'].decode('utf8', errors='replace')
        out.name = d['name'].decode('utf8', errors='replace')
        if d['disp'] is None:
            out.disp = None
        else:
            out.disp = (d['disp'][0], d['disp'][1],
                        d['disp_tail'].decode('utf8', errors='replace'))
        out.atoms = list(d['atoms'])
        out.patoms = list(d['patoms'])
        out.bonds = [tuple(b) for b in d['bonds']]
        out.cstates = [(pair, tail.decode('utf8', errors='replace') if isinstance(tail, bytes)
                        else tail) for pair, tail in d['cstates']]
        out.data = list(d['data'])
        out.log = [x.decode('utf8', errors='replace') for x in d['log']]
        fields = {}
        for key, value in d['fields']:
            fields.setdefault(key.decode('utf8', errors='replace'), []).append(
                value.decode('utf8', errors='replace'))
        out.fields = fields
        return out

    def next_index(self):
        """The lowest Sgroup number not already used, for a record being added.

        Numbers are the file's, not ours, so this only ever runs for a record built in memory.
        """
        used = {r.index for r in self.records if r.index != NO_INDEX}
        n = 1
        while n in used:
            n += 1
        return n

    def __repr__(self):
        return f'SGroupStore({len(self.records)} records, {len(self.aliases)} aliases)'


def checked_index(value, what, log):
    """`value` as a *reference* to an Sgroup number, or :data:`NO_INDEX` with a report.

    For ``PARENT`` and the external number, which point at a number rather than being one.  Call it at
    the parse site and nowhere later: the field spells absence as 65535, so a keyword stating 65535 and
    an omitted keyword are only distinguishable here.
    """
    if 0 <= value <= INDEX_MAX:
        return value
    log.append(LogRecord(
        'sgroup:index-out-of-range', (), f'{what}={value} outside 0..{INDEX_MAX}, reference dropped', LOST))
    return NO_INDEX


def normalize_indices(sgroups, log):
    """Bring every Sgroup number in `sgroups` inside ``0..INDEX_MAX``, in place.

    A record whose own number is out of domain is renumbered, not demoted to :data:`NO_INDEX`: two
    records stating 70000 would both become unnumbered and the V2000 writer would then emit them under
    the same regenerated number.  The replacement is the lowest number no other record claims.

    A ``PARENT`` naming an out-of-domain group is not rescued -- that reference is dropped where it is
    read, by :func:`checked_index`, and reported.
    """
    used = {sg.index for sg in sgroups if 0 <= sg.index <= INDEX_MAX}
    n = 1
    for sg in sgroups:
        if 0 <= sg.index <= INDEX_MAX:
            continue
        while n in used:
            n += 1
        used.add(n)
        log.append(LogRecord(
            'sgroup:renumbered', (), f'sgroup number {sg.index} outside 0..{INDEX_MAX}, renumbered to {n}', REPAIRED))
        sg.index = n


def parse_fielddisp(value, log=None):
    """Parse a ``FIELDDISP`` value into ``(x, y, rest)``.

    ``FIELDDISP`` is a fixed-layout string, not free text::

        "   49.5979   -3.8125    DA    ALL  1       5"
          x (F10.4)  y (F10.4)  <-------- display styling -------->

    The x/y are display coordinates in the molecule's own frame, which is why they are parsed out
    rather than passed through: re-emitting the string unchanged after the molecule moved anchors the
    label at empty space.  The styling tail is round-tripped verbatim.

    Returns ``None`` when the value is too short, the coordinates do not parse, or they fall outside
    what the layout can express (see :data:`DISP_MAX`); the caller then keeps the raw string in
    ``fields``, so the anchor survives as opaque text rather than being clamped.
    """
    if value is None or len(value) < 20:
        if log is not None and value:
            log.append(LogRecord('sgroup:fielddisp-short', (), f'FIELDDISP too short to parse: {value!r:.40}', LOST))
        return None
    try:
        x = float(value[0:10])
        y = float(value[10:20])
    except ValueError:
        if log is not None:
            log.append(LogRecord(
                'sgroup:fielddisp-bad-coords', (), f'FIELDDISP coordinates unparseable: {value[:20]!r}', LOST))
        return None
    if not (DISP_MIN <= x <= DISP_MAX and DISP_MIN <= y <= DISP_MAX):
        if log is not None:
            log.append(LogRecord(
                'sgroup:fielddisp-out-of-range', (),
                f'FIELDDISP coordinates outside F10.4 range, kept as text: {value[:20]!r}', LOST))
        return None
    return (x, y, value[20:])


def format_fielddisp(disp, log=None):
    """Render ``(x, y, rest)`` back into a ``FIELDDISP`` value.

    A coordinate wider than the field is clamped to the field's limit and reported: emitting the wide
    number would shift every column after it, and omitting ``FIELDDISP`` would delete the anchor.
    """
    x, y, rest = disp
    if not (DISP_MIN <= x <= DISP_MAX and DISP_MIN <= y <= DISP_MAX):
        if log is not None:
            log.append(LogRecord(
                'sgroup:fielddisp-clamped', (), f'FIELDDISP anchor ({x:g}, {y:g}) does not fit F10.4, clamped',
                REPAIRED))
        x = min(max(x, DISP_MIN), DISP_MAX)
        y = min(max(y, DISP_MIN), DISP_MAX)
    return f'{x:10.4f}{y:10.4f}{rest}'


#: The FIELDDISP tail written when a caller states an anchor and nothing else: absolute placement,
#: attached, all occurrences.  Copied verbatim rather than composed -- the layout is fixed-width, and a
#: reader that understands only the fixed spelling still finds the anchor.
FIELDDISP_TAIL = '    DA    ALL  1       5'


def add_data_sgroup(molecule, name, data, *, atoms=(), bonds=(), disp=None, log=None):
    """Attach one ``DAT`` S-group to `molecule` and return it.

    The whole of a data label in one call: ``FIELDNAME``, ``FIELDDATA``, the atoms and bonds it labels,
    and a ``FIELDDISP`` anchor.  The record is APPENDED -- the core's ``set_sgroups`` replaces the whole
    set, which is the rebuild this saves the caller.  Both emitters write it, so `mol(m)` and
    `mol(m, version=3000)` both carry the label.

    :param name: ``FIELDNAME``.
    :param data: ``FIELDDATA`` -- ``str``, ``bytes``, or a sequence of either for a multi-value datum,
        in file order.
    :param atoms: the atom numbers labelled; empty for a datum about the whole molecule.
    :param bonds: ``(n, m)`` endpoint pairs.  A bond label lists its endpoints in `atoms` as well,
        which is what puts the auto anchor on the bond instead of beside it.
    :param disp: ``None`` (default) computes the anchor from the referenced atoms' coordinates;
        ``(x, y)`` states it with :data:`FIELDDISP_TAIL`; ``(x, y, tail)`` states both; ``False``
        writes no ``FIELDDISP``.
    :param log: a list to append reports to.

    The returned record is a snapshot: the molecule is the storage, so editing it afterwards changes
    nothing -- call again, or go through ``set_sgroups``.
    """
    log = [] if log is None else log
    atoms = [int(n) for n in atoms]
    bonds = [(int(n), int(m)) for n, m in bonds]
    numbers = set(molecule.atom_numbers)
    for n in atoms:
        if n not in numbers:
            raise ValueError(f'atom {n} is not in this molecule')
    for n, m in bonds:
        if molecule.order_of(n, m) is None:
            raise ValueError(f'there is no bond {n}-{m} in this molecule')

    store = SGroupStore.from_molecule(molecule)
    record = SGroup('DAT', index=store.next_index())
    record.name = name
    record.data = [x if isinstance(x, bytes) else str(x).encode('utf8')
                   for x in ((data,) if isinstance(data, (str, bytes)) else data)]
    record.atoms = atoms
    record.bonds = bonds
    record.disp = _resolve_disp(molecule, atoms, disp, log)
    store.records.append(record)
    store.to_molecule(molecule, log)
    return record


def _resolve_disp(molecule, atoms, disp, log):
    """The record's ``(x, y, tail)``, or ``None``."""
    if disp is False:
        return None
    if disp is not None:
        if len(disp) == 2:
            return (float(disp[0]), float(disp[1]), FIELDDISP_TAIL)
        x, y, tail = disp
        return (float(x), float(y), tail)
    if not atoms:
        log.append(LogRecord('sgroup:no-anchor', (),
                             'FIELDDISP not written: the record references no atom to anchor to',
                             LOST))
        return None
    if not molecule.has_coordinates:
        log.append(LogRecord('sgroup:no-anchor', tuple(atoms),
                             'FIELDDISP not written: the molecule states no coordinates; clean2d() '
                             'or an explicit disp= gives the label an anchor', LOST))
        return None
    points = [molecule.xy_of(n) for n in atoms]
    return (sum(p[0] for p in points) / len(points),
            sum(p[1] for p in points) / len(points), FIELDDISP_TAIL)


def data_sgroups(molecule, name=None):
    """`molecule`'s ``DAT`` records, or just those whose ``FIELDNAME`` is `name`.

    Snapshots, for the reason :func:`add_data_sgroup` gives: the arena is the storage.
    """
    store = SGroupStore.from_molecule(molecule)
    return store.data_records() if name is None else store.by_name(name)


def resolve_output(mol, title, sgroups, log=None):
    """What a writer should emit for `title` and `sgroups`, given what the caller did or did not say.

    ``None`` means the caller did not say, and the answer is then what the molecule holds; anything
    else is an instruction and wins.  The sentinel is ``None`` and not a falsy test, because
    ``sgroups=SGroupStore()`` must be able to say "write none".

    Title and S-groups resolve together because in the arena the title is handle 0 of the very blob the
    S-group records live in.  The title comes back as ``str``, which is what
    :attr:`MoleculeContainer.title` now is; a caller-supplied buffer is decoded with
    ``surrogateescape``, so nothing is lost and nothing is logged.
    """
    if title is None:
        title = mol.title
    elif isinstance(title, (bytes, bytearray, memoryview)):
        # A caller may hand a raw name line.  `surrogateescape` and not `replace`: the byte survives to
        # the writer, which re-encodes it.  Only XML cannot carry it -- see `xml_text`.
        title = bytes(title).decode('utf8', 'surrogateescape')
    if sgroups is None:
        sgroups = SGroupStore.from_molecule(mol)
    return title, sgroups
