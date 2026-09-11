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
"""``Ctab`` -- the one intermediate every MDL and MRV path goes through.

Index alphabet inside a ``Ctab``: **0-based positions** in ``atoms`` and ``bonds``.  Never a stable
id and never the file's 1-based number.  V3000 atom indices in particular are arbitrary positive
integers rather than positions, so the parser maps them; the writer regenerates them.
"""

from ...core import (LogRecord, LOST, MoleculeContainer, REPAIRED, STEREO_ABS, STEREO_AND,
                     STEREO_OR, WEDGE_DOWN, WEDGE_EITHER, WEDGE_NONE, WEDGE_UP)
from ._errors import MalformedCtfile, UnsupportedCtfile
from ._hydrogens import MOLFILE_CHANNELS, calc_implicit
from ._sgroup import SGroupStore
from ...core.wedge import assign_parities


__all__ = ['Ctab', 'CtabAtom', 'CtabBond', 'WEDGE_FROM_V2000', 'WEDGE_FROM_V3000',
           'WEDGE_TO_V2000', 'WEDGE_TO_V3000', 'STEREO_FROM_COLLECTION', 'STEREO_TO_COLLECTION',
           'STORABLE_ORDERS', 'QUERY_BOND_TYPES', 'LABEL_ELEMENT', 'order_from_bond_type']


# What an atom becomes when the file's symbol field holds free text -- a drawn label like `Me`, a
# polymer bead `Pol`, a registry identifier -- and so names no element: the marker, element 0, with
# the text as its alias.  One element for every unnamed atom type, so nothing downstream has to know
# which spelling a file used, and no borrowed element: a marker's valence, formula contribution and
# hydrogen count are all zero by definition, where a borrowed `C` reads as chemistry to anything that
# does not also consult the alias.  A label that abbreviates a real fragment becomes one only when a
# caller asks for it, by an explicit pass over the aliases after reading.
LABEL_ELEMENT = 'R'


# The wedge encodings of the two MDL versions disagree, and mixing them up is the classic MDL bug:
# V2000 spells "down" 6 while V3000 spells it 3, and V3000's 2 is V2000's 4.  Both are declared here
# once, as tables, so no parser or writer restates the numbers (RULES.md section 6).
WEDGE_FROM_V2000 = {0: WEDGE_NONE, 1: WEDGE_UP, 4: WEDGE_EITHER, 6: WEDGE_DOWN}
WEDGE_FROM_V3000 = {0: WEDGE_NONE, 1: WEDGE_UP, 2: WEDGE_EITHER, 3: WEDGE_DOWN}
WEDGE_TO_V2000 = {WEDGE_NONE: 0, WEDGE_UP: 1, WEDGE_EITHER: 4, WEDGE_DOWN: 6}
WEDGE_TO_V3000 = {WEDGE_NONE: 0, WEDGE_UP: 1, WEDGE_EITHER: 2, WEDGE_DOWN: 3}

# Enhanced stereo collections.  Note the direction, which is easy to invert: CTfile STERAC is a
# *racemate* -- both enantiomers present -- which is the AND semantics, and STEREL is one enantiomer
# of unknown absolute configuration, which is OR.  Same partition as CXSMILES `|&N:|` and `|oN:|`.
STEREO_FROM_COLLECTION = {'ABS': STEREO_ABS, 'RAC': STEREO_AND, 'REL': STEREO_OR}
STEREO_TO_COLLECTION = {STEREO_ABS: 'STEABS', STEREO_AND: 'STERAC', STEREO_OR: 'STEREL'}

# Bond orders a CTfile can state and this reader stores as stated.  Order 4 is in the set: a file
# that says aromatic is stored aromatic and one that says Kekule is stored Kekule.  Nothing here
# kekulises; that is an explicit operation a caller asks for.
STORABLE_ORDERS = frozenset((1, 2, 3, 4, 8))

# The bond types a molecule cannot hold, being CONSTRAINTS rather than bonds: "either of two orders"
# has nowhere to go in a structure container.  A record containing one is refused by name -- a refusal
# at the answer boundary -- so the caller can reach for a query reader instead of guessing why.
QUERY_BOND_TYPES = {5: 'single or double', 6: 'single or aromatic', 7: 'double or aromatic'}


def order_from_bond_type(type_, where, log):
    """A CTfile bond type as a storable chython order.  Raises for a query type.

    `where` names the bond the way the calling format numbers them (``'bond line 7'`` for V2000,
    ``'bond 7'`` for V3000), so the message points at something the caller can find in the file.

    The types, and why each is read as it is:

    - **1, 2, 3, 4** -- stated as stated.  Aromatic included; see :data:`STORABLE_ORDERS`.
    - **5, 6, 7** -- :data:`QUERY_BOND_TYPES`, refused with ``UnsupportedCtfile``.
    - **8** -- the spec's query "any bond", read as chython's order 8 with the caveat logged: chython's
      writer puts its internal order straight into this column, so every dative or ionic contact it has
      written is on disk as type 8.  A genuine query file is caught by its query atoms.
    - **9** -- the coordination bond, which *is* chython's order 8.  Exact.
    - **10** -- V3000's hydrogen bond.  chython has no order for one; read as single it would join the
      valence arithmetic, while order 8 is excluded from valence, degree and heteroatom counts and is
      removable by name, so it is stored as a non-covalent contact and prefixed ``unsupported: ``
      because the *kind* of contact is lost.
    - **anything else** -- not a CTfile bond type at all; read as single, logged.

    The writer states type 8 silently: the caveat is about provenance, and a file this library wrote did
    not come from a query editor.  The accepted cost is that a conforming third-party reader is entitled
    to read type 8 as the query it is defined to be, so such a file is not portable.
    """
    if type_ in STORABLE_ORDERS:
        if type_ == 8:
            log.append(LogRecord('ctab:type-8-as-dative', (),
                                 f'{where}: type 8 read as chython order 8, a dative or ionic contact. The '
                                 f'CTfile specification calls type 8 the query "any bond"; it is read this way '
                                 f'because chython writes order 8 into this column. If this file came from a '
                                 f'query editor the bond is a wildcard and this molecule is not what it means'))
        return type_
    if type_ in QUERY_BOND_TYPES:
        raise UnsupportedCtfile(f'{where}: type {type_} is a query bond '
                                f'({QUERY_BOND_TYPES[type_]}), which a molecule cannot represent. '
                                f'Read this file with a query reader')
    if type_ == 9:
        log.append(LogRecord('ctab:coordination-bond', (),
                             f'{where}: coordination bond type 9 read as chython order 8', REPAIRED))
        return 8
    if type_ == 10:
        log.append(LogRecord('ctab:hydrogen-bond', (),
                             f'unsupported: {where}: hydrogen bond type 10 stored as chython order 8, a '
                             f'non-covalent contact -- there is no hydrogen-bond order, and reading it as '
                             f'single would add it to the valence of both atoms. That it was a hydrogen bond '
                             f'is not preserved', LOST))
        return 8
    log.append(LogRecord('ctab:unknown-bond-type', (),
                         f'{where}: type {type_} is not a CTfile bond type, read as single', REPAIRED))
    return 1


class CtabAtom:
    """One atom line's worth of information, in the file's own terms.

    ``valence`` is the file's explicit total-valence statement (V2000 ``vvv``, V3000 ``VAL=``) with
    ``None`` for "not stated" and ``0`` for the spec's zero-valence marker, kept separate from anything
    derived.

    ``parity`` is the atom block's own stereo statement -- V2000 ``sss``, V3000 ``CFG=`` -- where 1 is
    odd, 2 is even, 3 is "either" and 0 is unstated.  A real stereo source rather than a note in the
    margin: :meth:`Ctab.build` passes it to :func:`chython.core.wedge.assign_parities` as ``stated``,
    so a record whose layout says nothing still comes back configured, the drawing winning where both
    speak and the disagreement logged.  Read direction only -- the writers state the configuration in
    the wedge and decline this column, one statement being better than two that can contradict.

    ``label`` is the symbol field's text when it names no element, with ``element`` then holding
    :data:`LABEL_ELEMENT`; ``None`` for every atom whose element the file did state.
    :meth:`Ctab.build` files the text as the atom's alias.

    ``r_index`` is the R group number, from the ``R<n>`` symbol column or an ``M  RGP`` line; 0 for
    an unindexed marker.  Meaningful only when ``element`` is ``'R'``.
    """
    __slots__ = ('element', 'charge', 'isotope', 'radical', 'map_number', 'x', 'y', 'z', 'parity',
                 'valence', 'stated_h', 'mass_diff', 'file_index', 'label', 'r_index')

    def __init__(self, element='C', x=0.0, y=0.0, z=0.0):
        self.element = element
        self.charge = 0
        self.isotope = 0
        self.radical = False
        self.map_number = 0
        self.x = x
        self.y = y
        self.z = z
        self.parity = 0        # the atom block's stereo field; a read channel, see the docstring
        self.valence = None    # explicit total valence, None = unstated
        self.stated_h = None   # an explicit H count (HCOUNT=, MRV_IMPLICIT_H), None = unstated
        self.mass_diff = 0     # V2000 `dd`, superseded by M ISO
        self.file_index = 0    # the number this atom carried in the file, for diagnostics
        self.label = None      # free text where an element symbol belongs; see the class docstring
        self.r_index = 0       # the R group number, from R<n> or M  RGP; 0 for an unindexed marker

    def __repr__(self):
        return f'CtabAtom({self.element}, charge={self.charge}, xy=({self.x:.4f}, {self.y:.4f}))'


class CtabBond:
    """One bond line.  ``a`` and ``b`` are 0-based atom positions; ``a`` is the wedge's narrow end.

    ``order`` is the file's bond type verbatim, 4 for aromatic included, and ``Ctab.build`` stores it as
    it stands.  The query types 5, 6 and 7 never reach here -- the parser refuses the record by name --
    while type 8 does, as chython's order 8; see :func:`order_from_bond_type`.

    ``configuration`` is a **non-geometric** double-bond descriptor: ``(letter, refs)``, where `letter`
    is ``'C'`` or ``'T'`` and `refs` is the four 0-based atom positions the document measured it over,
    empty when it names none.  No CTAB carries it in either direction -- MDL states a double bond's
    configuration in the coordinates and nowhere else, while CML and MRV spell the letter out.
    ``Ctab.build`` hands it to :func:`chython.core.wedge.assign_parities` as ``configurations``, ranked
    exactly as the atom parity field is: consulted where the layout is silent, overridden where it is
    not, logged where the two disagree.  A writer DERIVES this from ``mol.parity_of`` via
    :func:`chython.core.wedge.cis_trans_letter` and must never replay what a reader stored here, so a
    self-contradicting input file (letter ``C``, coordinates ``T``) does not round-trip byte for byte.
    """
    __slots__ = ('a', 'b', 'order', 'wedge', 'topology', 'reacting_center', 'configuration')

    def __init__(self, a, b, order=1, wedge=WEDGE_NONE):
        self.a = a
        self.b = b
        self.order = order
        self.wedge = wedge
        self.topology = 0
        self.reacting_center = 0
        self.configuration = None   # ('C'|'T', refs) or None; see the class docstring

    def __repr__(self):
        return f'CtabBond({self.a}-{self.b}, order={self.order}, wedge={self.wedge})'


class Ctab:
    """A parsed connection table, version- and surface-independent."""
    __slots__ = ('title', 'program', 'comment', 'dimensionality', 'chiral', 'atoms', 'bonds',
                 'sgroups', 'groups', 'aliases', 'log', 'meta', 'unknown_hydrogens', 'channels')

    def __init__(self):
        self.title = ''
        self.program = ''
        self.comment = ''
        self.dimensionality = ''
        self.chiral = False
        self.atoms = []
        self.bonds = []
        self.sgroups = []
        # atom position -> (stereo kind, group id).  Populated from V3000 COLLECTION, V2000's
        # enhanced-stereo Sgroup encoding, or MRV's @mrvStereoGroup -- one model, three spellings.
        self.groups = {}
        self.aliases = {}
        self.log = []
        #: Record metadata the dialect that filled this table collected -- an SDF's data fields, a CML
        #: `<propertyList>`.  `build` copies it onto `mol.meta`, which is the one place it lives.
        self.meta = {}
        # Set by `build`: stable ids whose implicit hydrogen count nothing determined.
        self.unknown_hydrogens = ()
        # How the dialect that filled this table spells a stated hydrogen count and a stated total
        # valence.  Read only by `build`, to compose the advice half of an unknown-count log line, so
        # that it never names a channel the file in front of the caller does not have.
        self.channels = MOLFILE_CHANNELS

    def __len__(self):
        return len(self.atoms)

    def bond_index(self):
        """``{(a, b): position}`` for every bond, both orientations.

        The file's S-group bond references are 1-based positions into the bond block; this is the
        table that turns them into the endpoint pairs :class:`SGroup` holds.
        """
        out = {}
        for i, b in enumerate(self.bonds):
            out[(b.a, b.b)] = i
            out[(b.b, b.a)] = i
        return out

    def build(self, *, ignore_stereo=False):
        """Build a core ``MoleculeContainer``.

        Returns ``(molecule, store, log)`` where `store` is an :class:`SGroupStore` keyed in the
        molecule's stable ids and `log` is every recovery this build made, appended to whatever the
        parser already recorded.  Also sets :attr:`unknown_hydrogens`.

        A bad structure is not a reason to fail: an illegal valence, a nonsense charge and an
        underivable hydrogen count are all stored and logged.  What still raises is a record that cannot
        be *parsed* or whose atoms cannot be stored at all.

        The order of operations is not arbitrary and each step depends on the one before:

        1. atoms and bonds -- **at the orders the file states**, aromatic bond type 4 included, so a
           Kekule drawing is stored Kekule and an aromatic one aromatic, in either direction;
        2. coordinates and wedges -- the as-drawn direction, stored before anything interprets it;
        3. implicit hydrogens;
        4. stereo groups, then parities -- parities last, perception needing the finished constitution,
           and each computed in the frame the core names.
        """
        log = list(self.log)
        mol = MoleculeContainer()
        sids = []
        # One edit scope for the whole constitution.  Outside a scope every `add_atom` applies its
        # journal at once, so an N-atom record costs N buffer copies and reading an SDF is quadratic in
        # its largest molecule.  It is also why the steps that read the finished constitution run
        # outside this block: `neighbors_of`, `order_of` and the geometry accessors all require a clean
        # arena.  Core validation happens before anything is journalled, so the degradation loop below
        # still sees its `ValueError` at the call.
        z = {}
        aromatic = []
        with mol.edit():
            for i, a in enumerate(self.atoms):
                # A property outside the core's domain -- a charge past +8, a mass number past 65535 --
                # must not cost the atom: dropping it renumbers every atom after it and invalidates
                # every bond, collection and S-group reference in the record.  So the atom goes in
                # either way and the properties are surrendered one at a time, least useful first.
                full = {'charge': a.charge, 'isotope': a.isotope, 'radical': a.radical,
                        'map_number': a.map_number}
                for drop in ((), ('charge',), ('charge', 'isotope'),
                             ('charge', 'isotope', 'map_number')):
                    kwargs = {k: v for k, v in full.items() if k not in drop}
                    try:
                        sid = mol.add_atom(a.element, **kwargs)
                    except ValueError as e:
                        reason = e
                        continue
                    if drop:
                        log.append(LogRecord('ctab:atom-property-dropped', (),
                                             f'atom {i + 1} {a.element}: {reason}; '
                                             f'dropped {", ".join(drop)}', LOST))
                    break
                else:
                    raise MalformedCtfile(f'atom {i + 1} {a.element} cannot be stored: {reason}')
                sids.append(sid)
                if a.element == 'R' and a.r_index:
                    mol.set_r_index(sid, a.r_index)
                elif a.r_index:
                    log.append(LogRecord('ctab:rgp-on-a-non-marker', (),
                                         f'atom {i + 1} {a.element}: M  RGP assigns an R group to an '
                                         f'atom that is not a marker; the assignment is dropped', LOST))

            seen = set()
            n = len(sids)
            for i, b in enumerate(self.bonds):
                if not (0 <= b.a < n and 0 <= b.b < n):
                    log.append(LogRecord('ctab:bond-out-of-range', (),
                                         f'bond {i + 1} references an atom outside the block, dropped', LOST))
                    continue
                if b.a == b.b:
                    log.append(LogRecord('ctab:self-loop', (),
                                         f'self-loop bond {i + 1} dropped', LOST))
                    continue
                key = (b.a, b.b) if b.a < b.b else (b.b, b.a)
                if key in seen:
                    log.append(LogRecord('ctab:duplicate-bond', (),
                                         f'duplicate bond {i + 1} dropped', LOST))
                    continue
                seen.add(key)
                order = b.order
                if order == 4:
                    aromatic.append(i + 1)
                elif order not in STORABLE_ORDERS:
                    log.append(LogRecord('ctab:bond-order-read-as-single', (),
                                         f'bond {i + 1} order {order} read as single', REPAIRED))
                    order = 1
                mol.add_bond(sids[b.a], sids[b.b], order)
                if b.reacting_center:
                    # Parsed and not stored: the arena has no per-bond field for it yet.  NOT prefixed
                    # `unsupported:` -- MDL models this and so will we, so it is a gap in our storage
                    # rather than a construct to send the caller to another tool for.
                    log.append(LogRecord('ctab:reacting-centre-not-stored', (),
                                         f'bond {i + 1}: reacting-centre code {b.reacting_center} not stored',
                                         LOST))
                if b.topology:
                    # Permanently ours to refuse: TOPO states "ring bond only" or "chain bond only",
                    # a QUERY constraint a structure record has nowhere to put.
                    log.append(LogRecord('ctab:topology-not-stored', (),
                                         f'unsupported: bond {i + 1}: TOPO/topology {b.topology} is a '
                                         f'query constraint, not stored', LOST))

            # Coordinates.  A file with every coordinate at the origin has no layout -- writing the
            # segment anyway would claim one and make `has_coordinates` a lie.
            if any(a.x or a.y or a.z for a in self.atoms):
                # Both segments on a 3D record: `SEG_XY` is the depiction and `SEG_CONFORMERS` the
                # geometry, and a V3000 3D block states one thing serving as both.  Filling only the
                # conformer leaves the record undepictable; filling only `xy` throws the geometry away.
                solid = any(a.z for a in self.atoms)
                for pos, (sid, a) in enumerate(zip(sids, self.atoms)):
                    try:
                        mol.set_xy(sid, a.x, a.y)
                        if solid:
                            mol.set_xyz(sid, a.x, a.y, a.z)
                    except ValueError as e:
                        log.append(LogRecord('ctab:coordinates-dropped', (),
                                             f'atom {pos + 1} coordinates dropped: {e}', LOST))
                if solid:
                    # Parity assignment runs inside this edit session, where `xyz_of` refuses to answer
                    # (every geometry accessor requires a clean arena), so the reader keeps its own
                    # copy of z for the pre-seal question.
                    for sid, a in zip(sids, self.atoms):
                        z[sid] = a.z

            # Wedges, verbatim, before any interpretation.  The narrow end is the bond line's first
            # atom -- CTfile puts the wedge's point at atom 1 -- and that is the whole reason a
            # bond's atom order is information rather than an implementation detail.
            for i, b in enumerate(self.bonds):
                if b.wedge != WEDGE_NONE and b.a != b.b and 0 <= b.a < n and 0 <= b.b < n:
                    try:
                        mol.set_wedge(sids[b.a], sids[b.b], b.wedge)
                    except (KeyError, ValueError) as e:
                        log.append(LogRecord('ctab:wedge-dropped', (),
                                             f'wedge on bond {i + 1} dropped: {e}', LOST))

        # Nothing is logged for an aromatic bond: it is stored as the file drew it, and
        # `mol.aromatic_bond_count` answers "what representation is this molecule in" from the bonds
        # themselves, so no log line here could be more trustworthy than asking.

        # Implicit hydrogens.  Building an atom derives no count and an unset one stores as
        # `H_UNKNOWN`, so a record read without this step comes out with every count unknown.
        explicit_h = {}
        for i, a in enumerate(self.atoms):
            if a.stated_h is not None:
                explicit_h[sids[i]] = a.stated_h
        hydrogens = calc_implicit(mol, stated=explicit_h,
                                  valences={sids[i]: a.valence for i, a in enumerate(self.atoms)
                                            if a.valence is not None},
                                  channels=self.channels)
        log.extend(hydrogens.log)
        # Published as one attribute: the molecule answers "how many" in constant time
        # (`mol.unknown_h_count`), and this slot answers "which ones", which is what a repair needs.
        self.unknown_hydrogens = hydrogens.unknown

        # Enhanced stereo groups.  Written before parities because a group is a statement about a
        # centre that exists whether or not its configuration is known.
        if self.groups:
            with mol.edit():
                for pos, (kind, group) in self.groups.items():
                    if pos >= len(sids):
                        log.append(LogRecord('ctab:stereo-out-of-range', (),
                                             f'stereo collection references atom {pos + 1}, out of range'))
                        continue
                    try:
                        mol.set_stereo_group(sids[pos], kind, group)
                    except ValueError as e:
                        log.append(LogRecord('ctab:stereo-group-dropped', (),
                                             f'stereo group on atom {pos + 1} dropped: {e}', LOST))

        if not ignore_stereo:
            # The atom parity field is measured in the file's own atom-block order, so the positions
            # travel with it.  The stated double-bond configurations travel keyed by the bond's two
            # stable ids, low first, that being the only name for a bond the core recognises; a
            # reference frame is translated with them, and a bond whose frame names an out-of-range
            # position is dropped whole, half a frame being a different statement rather than a
            # narrower one.
            configurations = {}
            count = len(sids)
            for i, b in enumerate(self.bonds):
                if b.configuration is None or not (0 <= b.a < count and 0 <= b.b < count) \
                        or b.a == b.b:
                    continue
                letter, refs = b.configuration
                if any(not 0 <= r < count for r in refs):
                    log.append(LogRecord('ctab:configuration-out-of-range', (),
                                         f'bond {i + 1}: stated configuration {letter} references an atom '
                                         f'outside the atom block, dropped', LOST))
                    continue
                x, y = sids[b.a], sids[b.b]
                configurations[(x, y) if x < y else (y, x)] = (letter, tuple(sids[r] for r in refs))

            assign_parities(mol, z, log,
                            stated={sids[i]: a.parity for i, a in enumerate(self.atoms)
                                    if a.parity in (1, 2)},
                            positions={sid: i for i, sid in enumerate(sids)},
                            configurations=configurations)

        # S-groups, translated from Ctab positions into the molecule's stable ids.  Bond references
        # arrive as endpoint pairs already (the parsers translate the file's 1-based indices), so
        # this is one uniform atom-position mapping.
        atom_map = {i: sid for i, sid in enumerate(sids)}
        # A symbol field holding free text is filed here: the label is the file's only statement about
        # that atom, and a display label is what an alias is.  `setdefault`, so an `A  <n>` line --
        # which says "label" outright -- outranks one recovered from the element column.
        aliases = dict(self.aliases)
        for i, a in enumerate(self.atoms):
            if a.label is not None:
                aliases.setdefault(i, a.label)
        store = SGroupStore(self.sgroups, aliases).translate(atom_map)
        log.extend(store.log)
        store.log = []

        # Onto the molecule as well as into the returned store: a title or an S-group handed back only
        # BESIDE the molecule is dropped by the first consumer that keeps just the molecule, which is
        # most of them.  Both calls are skipped when there is nothing to say, each being a full blob
        # rebuild that most records do not need.
        if self.title:
            mol.set_title(self.title)
        if self.meta:
            mol.meta.update(self.meta)
        store.to_molecule(mol, log)

        # The log goes onto the molecule for the same reason, and is still RETURNED, because a caller
        # collecting into its own `log=` list may still do so.  `absorb` and not `extend`: it stamps
        # the `'read'` stage on every record that did not name one, so an `edit:sgroup` line the core
        # wrote during `to_molecule` keeps the stage `mol.sgroup_log` filters on.  The `if` asks
        # whether there is anything to say, never whether to say it -- `mol.log` builds its storage on
        # first touch, and a clean record should not pay for an empty one.
        if log:
            mol.log.absorb('read', log)

        return mol, store, log

    def __repr__(self):
        return (f'Ctab({len(self.atoms)} atoms, {len(self.bonds)} bonds, '
                f'{len(self.sgroups)} sgroups, title={self.title!r:.20})')
