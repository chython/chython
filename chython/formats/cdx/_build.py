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
"""The object tree to a :class:`~chython.formats.ctfile._ctab.Ctab`: a CDX document read as a molfile.

A document without a reaction step is **one record**: every structure fragment on every page, as
components, the way a molfile holds disconnected pieces.  A document with one is a reaction, each
participant its own record (:func:`read_reaction`).  The query vocabulary is refused where V3000 refuses
it:

| CDX | Read as |
|---|---|
| node, ``Element`` absent | carbon |
| ``GenericNickname`` R, R<n>, ``*`` | the marker, R index from R<n> |
| ``GenericNickname`` A, AH, Q, QH, M, MH, X, XH; ``ElementList`` | refused, :class:`UnsupportedCdx` |
| ``Nickname`` / ``Fragment`` with an inner fragment | expanded, plus a SUP S-group carrying the label |
| ``Nickname`` without one, ``Formula``, any other label | the marker with the label as alias, logged |
| ``ExternalConnectionPoint`` outside an abbreviation | the marker; ``PolymerBead`` aliased ``Pol`` |
| ``Order`` 1, 2, 3, 1.5 | 1, 2, 3, 4 |
| ``Order`` dative, ionic | 8 |
| ``Order`` hydrogen, threecenter | 8, logged: the kind of contact is not kept |
| ``Order`` 4, 5, 6, half orders but 1.5 | 1, logged: no chython order |
| ``Order`` with two or more values, ``any`` | refused, :class:`UnsupportedCdx` |

Coordinates: ``x / 30``, ``-y / 30`` -- CDX's y axis points down, and 30 points per unit is the scale of
the CDX conversions of V3000 files, so such a pair reads to the same coordinates.
"""

from ._errors import MalformedCdx, UnsupportedCdx
from ..ctfile._ctab import Ctab, CtabAtom, CtabBond, LABEL_ELEMENT
from ..ctfile._hydrogens import StatedChannels
from ..ctfile._sgroup import SGroup, normalize_indices
from ...core import (LogRecord, LOST, REPAIRED, STEREO_ABS, STEREO_AND, STEREO_OR, WEDGE_DOWN,
                     WEDGE_EITHER, WEDGE_UP)
from ...core._core import element_symbols, R_INDEX_MAX


__all__ = ['read_document', 'read_reaction', 'has_reaction', 'CDX_CHANNELS', 'POINTS_PER_UNIT']


#: Points per chython coordinate unit.
POINTS_PER_UNIT = 30.0

CDX_CHANNELS = StatedChannels(count='a `NumHydrogens` attribute', valence=None)

_SYMBOLS = element_symbols()

#: `GenericNickname` values that are constraints over a set of elements, refused by name as V3000 does.
_QUERY_NICKNAMES = {
    'A': 'any heavy atom', 'AH': 'any atom', 'Q': 'any heteroatom', 'QH': 'any heteroatom or H',
    'M': 'any metal', 'MH': 'any metal or H', 'X': 'any halogen', 'XH': 'any halogen or H',
}

_ORDERS = {'1': 1, '2': 2, '3': 3, '1.5': 4, 'dative': 8, 'ionic': 8}
#: Orders a molecule cannot hold as such, and what each is stored as.
_LOSSY_ORDERS = {'hydrogen': 8, 'threecenter': 8, '4': 1, '5': 1, '6': 1, '0.5': 1, '2.5': 1, '3.5': 1,
                 '4.5': 1, '5.5': 1}

# `Display`: (wedge, narrow end at B?).  Hash, Dash, Bold and the hollow wedges are drawing styles.
_WEDGES = {'WedgeBegin': (WEDGE_UP, True), 'WedgeEnd': (WEDGE_UP, False),
           'WedgedHashBegin': (WEDGE_DOWN, True), 'WedgedHashEnd': (WEDGE_DOWN, False),
           'Wavy': (WEDGE_EITHER, True), 'WavyWedgeBegin': (WEDGE_EITHER, True),
           'WavyWedgeEnd': (WEDGE_EITHER, False)}

_ENHANCED = {'Absolute': STEREO_ABS, 'And': STEREO_AND, 'Or': STEREO_OR}

#: Node attributes read here.
_NODE_READ = frozenset(('p', 'Element', 'NodeType', 'GenericNickname', 'Charge', 'Isotope', 'Radical',
                        'NumHydrogens', 'EnhancedStereoType', 'EnhancedStereoGroupNum',
                        'ExternalConnectionType', 'ExternalConnectionNum', 'BondOrdering', 'AltGroupID',
                        'Attachments'))
#: Node attributes carrying no chemistry, or chemistry derived from what is read: a reader honouring
#: them builds the same molecule.  `AS` is ChemDraw's computed CIP label, `Geometry` its perceived
#: coordination shape.
_NODE_INERT = frozenset(('Z', 'AS', 'Geometry', 'AbnormalValence', 'LabelDisplay', 'LabelJustification',
                         'LabelAlignment', 'color', 'bgcolor', 'AtomNumber', 'AtomID', 'Warning',
                         'IgnoreWarnings', 'ChemicalWarning', 'Visible', 'LabelFont', 'LabelSize',
                         'LabelFace', 'CaptionFont', 'CaptionSize', 'CaptionFace', 'LineWidth',
                         'MarginWidth', 'HashSpacing', 'BoldWidth', 'BondLength', 'BondSpacing',
                         'ChainAngle', 'LabelLineHeight', 'CaptionLineHeight', 'LabelJustification',
                         'BoundingBox', 'SupersededBy', 'fgalpha', 'bgalpha', 'alpha', 'NeedsClean',
                         'ShowAtomQuery', 'ShowAtomStereo', 'ShowAtomEnhancedStereo', 'ShowAtomNumber',
                         'ShowAtomID', 'ShowResidueID', 'ShowTerminalCarbonLabels',
                         'ShowNonTerminalCarbonLabels', 'HideImplicitHydrogens', 'ShowEnhancedStereo'))
_BOND_READ = frozenset(('B', 'E', 'Order', 'Display', 'Topology', 'RxnParticipation'))
_BOND_INERT = frozenset(('Z', 'BS', 'Display2', 'DoublePosition', 'BondCircularOrdering', 'BeginAttach',
                         'EndAttach', 'color', 'bgcolor', 'LineWidth', 'BoldWidth', 'HashSpacing',
                         'MarginWidth', 'BondSpacing', 'BondSpacingAbs', 'BondLength', 'LabelFont',
                         'LabelSize', 'LabelFace', 'Warning', 'IgnoreWarnings', 'ChemicalWarning',
                         'Visible', 'SupersededBy', 'fgalpha', 'bgalpha', 'alpha', 'ShowBondQuery',
                         'ShowBondStereo', 'ShowBondRxn', 'ShowBondDisplay', 'BoundingBox', 'CrossingBonds',
                         'NeedsClean', 'Connectivity'))

# CTfile reacting-centre codes for `RxnParticipation`, so `Ctab.build` reports them as it does RXCTR.
_RXN = {'ReactionCenter': 1, 'MakeOrBreak': 4, 'ChangeType': 8, 'MakeAndChange': 12,
        'NotReactionCenter': -1, 'NoChange': 2}
_TOPOLOGY = {'Ring': 1, 'Chain': 2}

# Containers walked for structure.  Any other object -- text, graphics, arrows, annotations -- holds no
# fragment a molecule is made of, and is counted into one log line rather than one line each.
_CONTAINERS = frozenset(('CDXML', 'page', 'group'))
# Document-wide tables: infrastructure, not a drawing object.
_DOCUMENT = frozenset(('colortable', 'fonttable'))


def read_document(root, log):
    """The `Ctab` of a CDX document tree: every structure in it, as one record."""
    state = _State(log)
    _walk(root, state)
    if state.steps:
        log.append(LogRecord('cdx:reaction-read-as-molecule', (),
                             f'{state.steps} reaction step(s) in the document; every structure is read as '
                             f'one molecule and the reaction roles are not kept', LOST))
    return _finish(state)


def has_reaction(root):
    """Does the document hold a reaction step?"""
    return any(o.kind == 'step' for o in root.walk())


#: Step attributes listing the objects on each side, in `ReactionContainer` order.
_SIDES = ('ReactionStepReactants', 'ReactionStepProducts', 'ReactionStepObjectsAboveArrow',
          'ReactionStepObjectsBelowArrow')
_MAPS = ('ReactionStepAtomMap', 'ReactionStepAtomMapManual', 'ReactionStepAtomMapAuto')


def read_reaction(root, log, step=0):
    """``(reactants, products, agents, meta)`` of one reaction step, each side a list of Ctabs.

    `step` indexes the document's steps in document order, across schemes; negative counts from the
    end.  The steps not read are logged.  A side lists object ids: a fragment
    or a group of fragments is one molecule, a text object above or below the arrow is the conditions,
    filed in `meta`.  ``ReactionStepAtomMap`` pairs a reactant node with a product node; pair *k* becomes
    map number *k* on both.  A structure no side names is not part of the reaction, and is logged.
    """
    steps = [o for o in root.walk() if o.kind == 'step']
    if not -len(steps) <= step < len(steps):
        raise IndexError(f'step {step} requested; the document holds {len(steps)} reaction step(s)')
    index_of_step = step % len(steps)
    step = steps[index_of_step]
    if len(steps) > 1:
        log.append(LogRecord('cdx:steps-dropped', (),
                             f'{len(steps)} reaction steps in the document; step {index_of_step} is read and '
                             f'the other {len(steps) - 1} are not', LOST))
    index = {}
    for o in root.walk():
        if o.id:
            index.setdefault(o.id, o)
    maps = {}
    for name in _MAPS:
        ids = _ids(step.attributes.get(name))
        if len(ids) % 2:
            log.append(LogRecord('cdx:atom-map-odd', (), f'{name} lists {len(ids)} ids, not pairs; the last is '
                                                         f'dropped', REPAIRED))
        for a, b in zip(ids[::2], ids[1::2]):
            if a not in maps and b not in maps:
                number = len(maps) // 2 + 1
                maps[a] = maps[b] = number
    sides = []
    texts = {}
    used = set()
    for name in _SIDES:
        side = []
        for oid in _ids(step.attributes.get(name)):
            obj = index.get(oid)
            if obj is None:
                log.append(LogRecord('cdx:step-reference-dangling', (),
                                     f'{name} names object {oid}, which the document does not hold', LOST))
            elif obj.kind in ('fragment', 'group'):
                used.update(o.id for o in obj.walk() if o.kind == 'fragment')
                state = _State([], maps)        # each molecule's own log, so it lands on that molecule
                if obj.kind == 'fragment':
                    _fragment(obj, state, 0)
                else:
                    _walk(obj, state)
                if state.ctab.atoms:
                    side.append(_finish(state))
                else:
                    log.append(LogRecord('cdx:empty-component', (),
                                         f'{name} object {oid} holds no atom; not a molecule', LOST))
            elif obj.kind == 't':
                if obj.text:
                    texts.setdefault(name, []).append(obj.text)
            else:
                log.append(LogRecord('cdx:step-object-not-read', (),
                                     f'{name} object {oid} is a <{obj.kind}>, not a structure; not read', LOST))
        sides.append(side)
    stray = [o.id for o in _structures(root) if o.id not in used]
    if stray:
        log.append(LogRecord('cdx:structure-outside-reaction', (),
                             f'{len(stray)} fragment(s) on the page are named by no side of the reaction step '
                             f'and are not read', LOST))
    meta = {}
    for name, label in (('ReactionStepObjectsAboveArrow', 'above arrow'),
                        ('ReactionStepObjectsBelowArrow', 'below arrow')):
        if name in texts:
            meta[label] = '\n'.join(texts[name])
    reactants, products, above, below = sides
    return reactants, products, above + below, meta


def _structures(obj):
    """The top-level fragments of a document: reached through pages and groups, not inside a node."""
    for child in obj.children:
        if child.kind == 'fragment':
            yield child
        elif child.kind in _CONTAINERS:
            yield from _structures(child)


def _finish(state):
    """Log what the walk collected and hand back the Ctab."""
    log = state.log
    ctab = state.ctab
    if state.altgroups:
        log.append(LogRecord('cdx:rgroup-definitions', (),
                             f'unsupported: {state.altgroups} R-group definition(s) not read; the R atoms '
                             f'are kept as markers', LOST))
    if state.drawing:
        log.append(LogRecord('cdx:drawing-objects', (),
                             f'{sum(state.drawing.values())} drawing object(s) carry no structure and are '
                             f'not read: ' + ', '.join(f'{k} x{v}' for k, v in sorted(state.drawing.items()))))
    for (kind, name), count in sorted(state.unread.items()):
        log.append(LogRecord('cdx:attribute-not-read', (),
                             f'unsupported: {kind} attribute {name} on {count} object(s) not read', LOST))
    if not ctab.atoms:
        raise MalformedCdx('the document holds no structure: no atom in any fragment')
    if ctab.sgroups:
        normalize_indices(ctab.sgroups, log)
    return ctab


class _State:
    __slots__ = ('ctab', 'log', 'maps', 'steps', 'altgroups', 'drawing', 'unread', 'groups', 'free')

    def __init__(self, log, maps=None):
        self.ctab = ctab = Ctab()
        ctab.channels = CDX_CHANNELS
        ctab.log = log
        self.log = log
        self.maps = {} if maps is None else maps  # node id -> map number, from a reaction step
        self.steps = 0
        self.altgroups = 0
        self.drawing = {}
        self.unread = {}
        self.groups = {}     # (kind, file group number) -> chython group id
        self.free = {}       # kind -> next free group id

    def unread_attributes(self, obj, read, inert):
        for name in obj.attributes:
            if name not in read and name not in inert:
                key = (obj.kind, name)
                self.unread[key] = self.unread.get(key, 0) + 1

    def group(self, kind, number):
        """A chython group id for the file's group `number` of `kind`: kept in 1..63, else the next free."""
        key = (kind, number)
        if key not in self.groups:
            taken = {g for (k, _), g in self.groups.items() if k == kind}
            if 1 <= number <= 63 and number not in taken:
                self.groups[key] = number
            else:
                free = next((g for g in range(1, 64) if g not in taken), None)
                if free is None:
                    self.log.append(LogRecord('cdx:stereo-groups-exhausted', (),
                                              f'stereo group {number}: all 63 ids of its kind are taken; '
                                              f'its centres join group 63', LOST))
                    free = 63
                self.groups[key] = free
        return self.groups[key]


def _walk(obj, state):
    for child in obj.children:
        kind = child.kind
        if kind == 'fragment':
            _fragment(child, state, 0)
        elif kind in _CONTAINERS:
            _walk(child, state)
        elif kind == 'scheme':
            state.steps += sum(1 for c in child.children if c.kind == 'step')
        elif kind == 'altgroup':
            state.altgroups += 1
        elif kind in _DOCUMENT:
            continue
        elif kind:
            state.drawing[kind] = state.drawing.get(kind, 0) + 1


def _fragment(fragment, state, depth):
    """Add a fragment's atoms and bonds to the Ctab.  Returns ``{node id: position}`` and the connection
    points ``[(point id, attaching position)]`` an enclosing abbreviation pairs its bonds with."""
    ctab = state.ctab
    log = state.log
    where = {}          # node id -> position of the atom standing for it
    expansions = {}     # nickname node id -> (inner fragment's points, its ConnectionOrder)
    points = []         # this fragment's own connection point node ids, for an enclosing nickname
    nodes = [c for c in fragment.children if c.kind == 'n']
    bonds = [c for c in fragment.children if c.kind == 'b']
    for child in fragment.children:
        if child.kind == 'fragment':
            inner, _ = _fragment(child, state, depth)
            where.update(inner)
        elif child.kind not in ('n', 'b'):
            state.drawing[child.kind] = state.drawing.get(child.kind, 0) + 1

    for node in nodes:
        state.unread_attributes(node, _NODE_READ, _NODE_INERT)
        node_type = node.attributes.get('NodeType', 'Element')
        inner = next((c for c in node.children if c.kind == 'fragment'), None)
        label = next((c.text for c in node.children if c.kind == 't' and c.text), None)
        if node_type == 'ExternalConnectionPoint' and depth:
            points.append(node.id)
            continue
        if node_type in ('Nickname', 'Fragment') and inner is not None:
            start = len(ctab.atoms)
            inner_where, inner_points = _fragment(inner, state, depth + 1)
            where.update(inner_where)
            order = _ids(inner.attributes.get('ConnectionOrder'))
            expansions[node.id] = (inner_points, order, _ids(node.attributes.get('BondOrdering')))
            sg = SGroup('SUP', len(ctab.sgroups) + 1)
            sg.atoms = list(range(start, len(ctab.atoms)))
            if label:
                sg.fields['LABEL'] = [label]
            ctab.sgroups.append(sg)
            continue
        atom = _atom(node, node_type, label, state)
        atom.file_index = node.id
        atom.map_number = state.maps.get(node.id, 0)
        if node.id in where:
            log.append(LogRecord('cdx:duplicate-id', (),
                                 f'node id {node.id} appears twice; the later one is kept for bond references',
                                 REPAIRED))
        where[node.id] = len(ctab.atoms)
        ctab.atoms.append(atom)

    # Connection points: the atom each is bonded to inside this fragment is where an outer bond lands.
    attach = {}
    for bond in bonds:
        b, e = _id(bond.attributes.get('B')), _id(bond.attributes.get('E'))
        if b in points and e in where:
            attach[b] = where[e]
        elif e in points and b in where:
            attach[e] = where[b]
    own_points = [(p, attach.get(p)) for p in points]

    # Outer bonds of each expanded abbreviation, paired with its points: `BondOrdering` on the node
    # against `ConnectionOrder` on the inner fragment, each falling back to document order.
    redirect = {}       # (nickname id, bond id) -> attaching position
    for nid, (inner_points, order, ordering) in expansions.items():
        incident = [bond.id for bond in bonds if nid in (_id(bond.attributes.get('B')), _id(bond.attributes.get('E')))]
        ranked_bonds = [x for x in ordering if x in incident] + [x for x in incident if x not in ordering]
        by_id = dict(inner_points)
        ranked_points = [p for p in order if p in by_id] + [p for p, _ in inner_points if p not in order]
        if len(ranked_points) < len(ranked_bonds):
            log.append(LogRecord('cdx:abbreviation-unattached', (),
                                 f'abbreviation node {nid}: {len(ranked_bonds)} bond(s) and '
                                 f'{len(ranked_points)} attachment point(s); the excess bonds attach to its '
                                 f'first atom', REPAIRED))
        first = min((v for k, v in where.items() if k != nid), default=None)
        for i, bid in enumerate(ranked_bonds):
            target = by_id.get(ranked_points[i]) if i < len(ranked_points) else None
            redirect[(nid, bid)] = target if target is not None else first

    for bond in bonds:
        state.unread_attributes(bond, _BOND_READ, _BOND_INERT)
        b, e = _id(bond.attributes.get('B')), _id(bond.attributes.get('E'))
        if b in points or e in points:
            continue        # the attachment stub of this abbreviation; an outer bond replaces it
        pa = redirect.get((b, bond.id), where.get(b))
        pb = redirect.get((e, bond.id), where.get(e))
        if pa is None or pb is None:
            log.append(LogRecord('cdx:bond-dangling', (),
                                 f'bond {bond.id} references node {b if pa is None else e}, which is no atom of '
                                 f'its fragment; dropped', LOST))
            continue
        _bond(bond, pa, pb, state)
    return where, own_points


def _atom(node, node_type, label, state):
    a = node.attributes
    log = state.log
    atom = CtabAtom()
    x, y = _point(a.get('p'), node.id, log)
    atom.x, atom.y = x / POINTS_PER_UNIT, -y / POINTS_PER_UNIT
    if node_type in ('Element', 'Unspecified', ''):
        number = _int(a.get('Element'), 6, 'Element', node.id, log)
        if 0 < number < len(_SYMBOLS):
            atom.element = _SYMBOLS[number]
        else:
            log.append(LogRecord('cdx:element-out-of-range', (),
                                 f'node {node.id}: Element {number} names no element; read as the marker',
                                 LOST))
            atom.element = LABEL_ELEMENT
    elif node_type == 'GenericNickname':
        name = a.get('GenericNickname', '')
        upper = name.upper()
        if upper in _QUERY_NICKNAMES:
            raise UnsupportedCdx(f'node {node.id} is a query atom ({_QUERY_NICKNAMES[upper]}); this is a '
                                 f'query document and chython reads structures only')
        atom.element = LABEL_ELEMENT
        if upper in ('R', '*') or (upper.startswith('R') and upper[1:].isdigit()):
            group = int(upper[1:]) if upper[1:].isdigit() else 0
            if group > R_INDEX_MAX:
                log.append(LogRecord('cdx:r-index-too-wide', (),
                                     f'node {node.id}: R index {group} is past R_INDEX_MAX ({R_INDEX_MAX}), so '
                                     f'the marker is left unindexed', LOST))
                group = 0
            atom.r_index = group
        else:
            atom.label = name or label or '?'
            _label_lost(node, atom.label, log)
    elif node_type in ('ElementList', 'ElementListNickname'):
        raise UnsupportedCdx(f'node {node.id} is an atom list; this is a query document and chython reads '
                             f'structures only')
    elif node_type == 'ExternalConnectionPoint':
        atom.element = LABEL_ELEMENT
        if a.get('ExternalConnectionType') == 'PolymerBead':
            atom.label = 'Pol'
    elif node_type in ('MultiAttachment', 'VariableAttachment'):
        # The bond a position-variation or delocalised attachment draws lands on this point, as V3000's
        # ENDPTS bond lands on its `*` atom; the list of atoms it stands for is not kept.
        atom.element = LABEL_ELEMENT
        log.append(LogRecord('cdx:attachment-list', (),
                             f'unsupported: node {node.id}: {node_type} attachment list not read; the point is '
                             f'kept as the marker', LOST))
    elif node_type == 'LinkNode':
        number = _int(a.get('Element'), 6, 'Element', node.id, log)
        atom.element = _SYMBOLS[number] if 0 < number < len(_SYMBOLS) else 'C'
        log.append(LogRecord('cdx:link-node', (), f'unsupported: node {node.id}: link node repeat count '
                                                  f'not read', LOST))
    elif node_type in ('AnonymousAlternativeGroup', 'NamedAlternativeGroup'):
        # An R-group atom.  Its group number is its label, `R4`, as V3000's `RGROUPS=` states it; the
        # `altgroup` it points at holds the members, not the name.
        atom.element = LABEL_ELEMENT
        name = (label or '').strip().upper()
        if name.startswith('R') and name[1:].isdigit() and int(name[1:]) <= R_INDEX_MAX:
            atom.r_index = int(name[1:])
        elif label:
            atom.label = label
    else:
        # Nickname or Fragment without an inner fragment, Formula, Monomer: a label.
        atom.element = LABEL_ELEMENT
        atom.label = label or a.get('GenericNickname') or node_type
        _label_lost(node, atom.label, log)
    atom.charge = _int(a.get('Charge'), 0, 'Charge', node.id, log)
    atom.isotope = _int(a.get('Isotope'), 0, 'Isotope', node.id, log)
    radical = a.get('Radical', 'None')
    atom.radical = radical not in ('None', '0', '')
    if radical in ('Singlet', 'Triplet'):
        log.append(LogRecord('cdx:diradical', (), f'node {node.id}: Radical={radical} (a diradical) stored as '
                                                  f'a single radical', REPAIRED))
    if 'NumHydrogens' in a:
        atom.stated_h = _int(a['NumHydrogens'], None, 'NumHydrogens', node.id, log)
    kind = _ENHANCED.get(a.get('EnhancedStereoType'))
    if kind is not None:
        if kind == STEREO_ABS:
            group = 0
        else:
            group = state.group(kind, _int(a.get('EnhancedStereoGroupNum'), 1, 'EnhancedStereoGroupNum',
                                           node.id, log))
        state.ctab.groups[len(state.ctab.atoms)] = (kind, group)
    return atom


def _label_lost(node, label, log):
    log.append(LogRecord('cdx:label-as-marker', (),
                         f'node {node.id}: label {label!r} names no element and carries no fragment, so it is '
                         f'read as the marker {LABEL_ELEMENT!r} with {label!r} as its alias. Whatever the label '
                         f'abbreviates is not in the structure', LOST))


def _bond(bond, pa, pb, state):
    a = bond.attributes
    log = state.log
    tokens = a.get('Order', '1').split()
    if len(tokens) != 1 or tokens[0] == 'any':
        raise UnsupportedCdx(f'bond {bond.id}: Order {a["Order"]!r} is a query bond, which a molecule cannot '
                             f'represent')
    token = tokens[0]
    order = _ORDERS.get(token)
    if order is None:
        order = _LOSSY_ORDERS.get(token)
        if order is None:
            log.append(LogRecord('cdx:unknown-bond-order', (), f'bond {bond.id}: Order {token!r} is not a CDX '
                                                               f'order, read as single', REPAIRED))
            order = 1
        else:
            stored = 'a non-covalent contact (order 8)' if order == 8 else 'single'
            log.append(LogRecord('cdx:bond-order-lost', (), f'unsupported: bond {bond.id}: Order {token} has no '
                                                            f'chython order, stored as {stored}', LOST))
    wedge, narrow_at_b = _WEDGES.get(a.get('Display', ''), (None, True))
    if wedge is not None and not narrow_at_b:
        pa, pb = pb, pa
    out = CtabBond(pa, pb, order)
    if wedge is not None:
        out.wedge = wedge
    out.reacting_center = _RXN.get(a.get('RxnParticipation'), 0)
    out.topology = _TOPOLOGY.get(a.get('Topology'), 0)
    state.ctab.bonds.append(out)


def _point(value, oid, log):
    if not value:
        return 0.0, 0.0
    try:
        x, y = (float(t) for t in value.split()[:2])
    except ValueError:
        log.append(LogRecord('cdx:bad-position', (), f'node {oid}: position {value!r:.30} unreadable; placed at '
                                                     f'the origin', REPAIRED))
        return 0.0, 0.0
    return x, y


def _int(value, default, name, oid, log):
    if value is None:
        return default
    try:
        return int(value)
    except ValueError:
        try:
            return int(float(value))
        except ValueError:
            log.append(LogRecord('cdx:bad-integer', (), f'node {oid}: {name}={value!r:.20} is not a number; '
                                                        f'the default is used', REPAIRED))
            return default


def _id(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _ids(value):
    if not value:
        return []
    out = []
    for t in value.split():
        try:
            out.append(int(t))
        except ValueError:
            pass
    return out
