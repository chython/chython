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
"""A molecule to the CDX object tree, and the tree to either spelling.

One fragment per connected component, one id counter for the document.  What a CDX document has no
spelling for is logged, never written approximately:

| Molecule | CDX |
|---|---|
| order 1, 2, 3, 4, 8 | ``Order`` omitted, 2, 3, 1.5, dative -- from the ligand when one end is a metal |
| wedge up / down / either, narrow end first | ``Display`` WedgeBegin / WedgedHashBegin / Wavy, ``B`` narrow |
| marker R<n>, ``Pol`` alias, other alias, bare | ``GenericNickname`` R<n>; polymer bead; label ``Nickname``; R |
| stereo group on an atom, an atropisomer's anchor, an axis's middle atom | ``EnhancedStereoType``, group number |
| stereo group on a double bond | not written: CDX has no bond-level enhanced stereo |
| hydrogen count | ``NumHydrogens`` on every atom but a carbon the valence rules reproduce |
| map numbers | a reaction step's ``ReactionStepAtomMap``; not written for a molecule |
| S-groups, data fields | not written yet |
"""

from ._binary import serialize_cdx
from ._build import POINTS_PER_UNIT
from ._cdxml import serialize_cdxml
from ._tree import CdxObject
from ..ctfile._ctab import axis_spelling
from ..ctfile._hydrogens import _needs_statement
from ..ctfile._sgroup import groups_for_write
from ...core import (LogRecord, LOST, STEREO_ABS, STEREO_AND, STEREO_OR, SU_ATROPISOMER, SU_HELICAL,
                     WEDGE_DOWN, WEDGE_EITHER, WEDGE_UP)
from ...core.reaction import ReactionContainer
from ...core.wedge import wedge_in_file_order, wedges_for_write


__all__ = ['write_cdx', 'write_cdxml', 'molecule_tree', 'reaction_tree']


_ORDERS = {1: None, 2: '2', 3: '3', 4: '1.5', 8: 'dative'}
_DISPLAY = {WEDGE_UP: 'WedgeBegin', WEDGE_DOWN: 'WedgedHashBegin', WEDGE_EITHER: 'Wavy'}
_ENHANCED = {STEREO_ABS: 'Absolute', STEREO_AND: 'And', STEREO_OR: 'Or'}


def write_cdxml(data, log):
    return serialize_cdxml(_tree(data, log), log)


def write_cdx(data, log):
    return serialize_cdx(_tree(data, log), log)


def _tree(data, log):
    return reaction_tree(data, log) if isinstance(data, ReactionContainer) else molecule_tree(data, log)


def _document():
    counter = iter(range(1, 1 << 31))
    root = CdxObject('CDXML', next(counter), {'CreationProgram': 'chython'})
    colors = CdxObject('colortable')
    colors.children = [CdxObject('color', 0, {'r': '1', 'g': '1', 'b': '1'}),
                       CdxObject('color', 0, {'r': '0', 'g': '0', 'b': '0'})]
    page = CdxObject('page', next(counter))
    root.children = [colors, page]
    return root, page, counter


def molecule_tree(mol, log):
    """The document tree of one molecule: one fragment per component, straight on the page."""
    root, page, counter = _document()
    page.children.extend(_molecule(mol, counter, log, (0.0, 0.0), group=False)[0])
    length = _bond_length([mol])
    if length:
        root.attributes['BondLength'] = length
    return root


def reaction_tree(reaction, log):
    """The document tree of one reaction: reactants, an arrow with the agents above it, products.

    Each molecule is shifted onto a left-to-right row, so the page shows a reaction rather than every
    molecule drawn over the others; a molecule's own coordinates are not changed.  A molecule of several
    components is a group, which a reader takes as one molecule.  Map numbers become the step's
    ``ReactionStepAtomMap``: a reactant node and a product node sharing a number are one pair.
    """
    start = len(log)
    root, page, counter = _document()
    gap, arrow_min = 1.0, 3.0
    cursor = 0.0
    nodes = {}
    placed = []

    def put(mol, x0, y_mode):
        box = _box(mol)
        ox = x0 - box[0]
        oy = -(box[2] + box[3]) / 2 if y_mode == 'centre' else 0.5 - box[2]
        objects, ids, cid = _molecule(mol, counter, log, (ox, oy), group=len(mol.connected_components) > 1,
                                      in_reaction=True)
        page.children.extend(objects)
        placed.append((mol, ids))
        return cid, box[1] - box[0]

    reactant_ids = []
    for mol in reaction.reactants:
        cid, width = put(mol, cursor, 'centre')
        reactant_ids.append(cid)
        nodes.setdefault('r', []).append(placed[-1][1])
        cursor += width + gap
    agents_width = sum(_box(m)[1] - _box(m)[0] for m in reaction.agents)
    agents_width += gap * max(len(reaction.agents) - 1, 0)
    tail = cursor
    head = tail + max(arrow_min, agents_width + 2 * gap)
    agent_ids = []
    x = tail + (head - tail - agents_width) / 2
    for mol in reaction.agents:
        cid, width = put(mol, x, 'above')
        agent_ids.append(cid)
        x += width + gap
    cursor = head + gap
    product_ids = []
    for mol in reaction.products:
        cid, width = put(mol, cursor, 'centre')
        product_ids.append(cid)
        nodes.setdefault('p', []).append(placed[-1][1])
        cursor += width + gap

    above, below = [], []
    meta = dict(reaction.meta)
    for key, target, y in (('above arrow', above, 0.3), ('below arrow', below, -0.6)):
        text = meta.pop(key, None)
        if text:
            oid = next(counter)
            page.children.append(CdxObject('t', oid, {'p': _points(tail + gap / 2, y)}, text=str(text)))
            target.append(oid)
    if meta:
        log.append(LogRecord('cdx:reaction-meta-unwritten', (),
                             f'not written: {len(meta)} reaction data field(s) -- no CDX spelling is implemented '
                             f'for them yet', LOST))

    arrow_id = next(counter)
    tail_p, head_p = _points(tail, 0.0), _points(head, 0.0)
    page.children.append(CdxObject('arrow', arrow_id, {
        'ArrowheadHead': 'Full', 'ArrowheadType': 'Solid',
        'Head3D': head_p + ' 0', 'Tail3D': tail_p + ' 0'}))

    pairs = []
    reactant_nodes = {}
    for ids, mol in zip(nodes.get('r', []), reaction.reactants):
        for sid, oid in ids.items():
            number = mol.atom(sid).map_number
            if number:
                reactant_nodes.setdefault(number, oid)
    for ids, mol in zip(nodes.get('p', []), reaction.products):
        for sid, oid in ids.items():
            number = mol.atom(sid).map_number
            if number and number in reactant_nodes:
                pairs.extend((reactant_nodes.pop(number), oid))

    step = {}
    for key, ids in (('ReactionStepReactants', reactant_ids), ('ReactionStepProducts', product_ids),
                     ('ReactionStepArrows', [arrow_id]), ('ReactionStepObjectsAboveArrow', agent_ids + above),
                     ('ReactionStepObjectsBelowArrow', below), ('ReactionStepAtomMap', pairs)):
        if ids:
            step[key] = ' '.join(map(str, ids))
    scheme = CdxObject('scheme', next(counter))
    scheme.children.append(CdxObject('step', next(counter), step))
    page.children.append(scheme)
    length = _bond_length(list(reaction.molecules()))
    if length:
        root.attributes['BondLength'] = length
    # One line for the reaction, not one per molecule.
    flat = [r for r in log[start:] if r.rule == 'cdx:no-coordinates']
    if len(flat) > 1:
        log[start:] = [r for r in log[start:] if r.rule != 'cdx:no-coordinates']
        log.append(LogRecord('cdx:no-coordinates', (),
                             f'{len(flat)} molecules have no 2D coordinates; their atoms are written at their row '
                             f'position, which a drawing program shows as one point each. Lay them out first',
                             LOST))
    return root


def _molecule(mol, counter, log, offset, group, in_reaction=False):
    """``(objects, {stable id: node id}, the id a reaction step names)`` for one molecule.

    `offset` is added to every coordinate, in chython units.  `group` wraps the fragments in a ``group``,
    which is what makes a several-component molecule one reaction participant.
    """
    sids = list(mol.atom_numbers)
    if sids and not mol.has_coordinates:
        log.append(LogRecord('cdx:no-coordinates', (),
                             'the molecule has no 2D coordinates; every atom is written at the origin, which a '
                             'drawing program shows as one point. Lay the molecule out first', LOST))
    ox, oy = offset
    xy = {}
    for sid in sids:
        x, y = mol.xy_of(sid) if mol.has_coordinates else (0.0, 0.0)
        xy[sid] = (x + ox, y + oy)

    wedges, _ = wedges_for_write(mol, log)
    wedge_of = {(narrow, wide): code for narrow, wide, code in wedges}
    enhanced = _enhanced(mol, log)
    aliases = mol.aliases
    lost = [name for name, present in (
        ('map numbers', not in_reaction and any(mol.atom(s).map_number for s in sids)),
        ('S-groups', bool(mol.sgroups)),
        ('data fields', bool(mol.meta))) if present]
    if lost:
        log.append(LogRecord('cdx:not-written', (), f'not written: {", ".join(lost)} -- no CDX spelling is '
                                                    f'implemented for them yet', LOST))

    holder = CdxObject('group', next(counter)) if group else None
    objects = []
    node_id = {}
    for component in mol.connected_components:
        fragment = CdxObject('fragment', next(counter))
        objects.append(fragment)
        for sid in sids:
            if sid in component:
                node_id[sid] = next(counter)
                fragment.children.append(_node(mol, sid, node_id[sid], xy[sid], aliases, enhanced))
        for bond in mol.bonds():
            if bond.n not in component:
                continue
            a, b, code = wedge_in_file_order(wedge_of, bond.n, bond.m)
            if bond.order == 8 and code is None and mol.atom(a).is_metal and not mol.atom(b).is_metal:
                a, b = b, a     # a dative bond starts at the donor: the ligand, not the metal
            attributes = {'B': str(node_id[a]), 'E': str(node_id[b])}
            order = _ORDERS.get(bond.order)
            if bond.order not in _ORDERS:
                log.append(LogRecord('cdx:bond-order-unwritten', (),
                                     f'bond {a}-{b}: order {bond.order} has no CDX spelling, written single', LOST))
            if order:
                attributes['Order'] = order
            if code in _DISPLAY:
                attributes['Display'] = _DISPLAY[code]
            fragment.children.append(CdxObject('b', next(counter), attributes))
    if holder is not None:
        holder.children = objects
        return [holder], node_id, holder.id
    return objects, node_id, objects[0].id if objects else 0


def _box(mol):
    """``(min x, max x, min y, max y)`` in chython units; zeros for a molecule without a layout."""
    if not mol.has_coordinates or not mol.atom_count:
        return 0.0, 0.0, 0.0, 0.0
    xs, ys = zip(*(mol.xy_of(s) for s in mol.atom_numbers))
    return min(xs), max(xs), min(ys), max(ys)


def _bond_length(molecules):
    lengths = []
    for mol in molecules:
        if mol.has_coordinates:
            for b in mol.bonds():
                (x1, y1), (x2, y2) = mol.xy_of(b.n), mol.xy_of(b.m)
                lengths.append(((x1 - x2) ** 2 + (y1 - y2) ** 2) ** .5)
    mean = sum(lengths) / len(lengths) if lengths else 0
    return _number(mean * POINTS_PER_UNIT) if mean else None


def _points(x, y):
    return f'{_number(x * POINTS_PER_UNIT)} {_number(-y * POINTS_PER_UNIT)}'


def _node(mol, sid, oid, xy, aliases, enhanced):
    atom = mol.atom(sid)
    x, y = xy
    attributes = {'p': f'{_number(x * POINTS_PER_UNIT)} {_number(-y * POINTS_PER_UNIT)}'}
    node = CdxObject('n', oid, attributes)
    if atom.is_r:
        alias = aliases.get(sid)
        if isinstance(alias, bytes):
            alias = alias.decode('utf-8', errors='replace')
        if atom.r_index:
            attributes['NodeType'] = 'GenericNickname'
            attributes['GenericNickname'] = f'R{atom.r_index}'
        elif alias == 'Pol':
            attributes['NodeType'] = 'ExternalConnectionPoint'
            attributes['ExternalConnectionType'] = 'PolymerBead'
        elif alias:
            attributes['NodeType'] = 'Nickname'
            node.children.append(CdxObject('t', 0, text=alias))
        else:
            attributes['NodeType'] = 'GenericNickname'
            attributes['GenericNickname'] = 'R'
    else:
        if atom.atomic_symbol != 'C':
            attributes['Element'] = str(atom.element)
        hydrogens = mol.implicit_h_of(sid)
        if hydrogens is not None and (atom.atomic_symbol != 'C' or _needs_statement(mol, sid)):
            attributes['NumHydrogens'] = str(hydrogens)
    if atom.charge:
        attributes['Charge'] = str(atom.charge)
    if atom.isotope:
        attributes['Isotope'] = str(atom.isotope)
    if atom.is_radical:
        attributes['Radical'] = 'Doublet'
    if sid in enhanced:
        kind, group = enhanced[sid]
        attributes['EnhancedStereoType'] = _ENHANCED[kind]
        if kind != STEREO_ABS:
            attributes['EnhancedStereoGroupNum'] = str(group)
    return node


def _enhanced(mol, log):
    """``{stable id: (kind, group)}``: each collection member as the node CDX can state it on."""
    out = {}
    for (kind, group), members in groups_for_write(mol).items():
        for member in members:
            if not isinstance(member, tuple):
                out[member] = (kind, group)
                continue
            anchor = mol.stereo_group_anchor_of(member)
            unit = None if anchor is None else mol.unit_of(anchor)
            if unit is not None and unit['kind'] in (SU_ATROPISOMER, SU_HELICAL):
                out[anchor] = (kind, group)        # its anchor pivot, as V3000 writes it
                continue
            spelling = axis_spelling(mol, member[0], member[1])
            if spelling is not None and spelling[0] == 'atom':
                out[spelling[1]] = (kind, group)   # an allene's middle carbon, where ChemDraw puts it
                continue
            log.append(LogRecord('cdx:bond-stereo-group-unwritten', member,
                                 f'stereo group {kind}:{group} on the axis {member[0]}-{member[1]} not written: '
                                 f'CDX states enhanced stereo on atoms only', LOST))
    return out


def _number(v):
    text = f'{v:.2f}'.rstrip('0').rstrip('.')
    return '0' if text in ('', '-0') else text
