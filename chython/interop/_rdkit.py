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
"""
RDKit conversion, both directions.  Reached through `interop.rdkit`.

Nothing kekulizes: an order-4 bond goes out as `BondType.AROMATIC` and comes back off
`GetBondType()`, never off `GetIsAromatic()` (RDKit's perception of it).  An import's losses go to the
returned container's `.log`, an export's to the caller's optional `log` list -- see `_records.py`;
`ToolkitError` means RDKit refused and there is no molecule left.

A reaction converts as a `rdChemReactions.ChemicalReaction`, the one RDKit object that holds three
sides and their atom-atom mapping at once, so `map_number` survives a round trip through it.
"""
from ..core import (H_IMPLICIT_MAX, LogRecord, MoleculeContainer as V3Molecule,
                    ReactionContainer as V3Reaction, STEREO_ABS, STEREO_AND, STEREO_OR,
                    SU_CIS_TRANS, SU_TETRA)
from ..exceptions import ToolkitError, UnconvertibleType
from ._records import deliver, mirror
from ._stereo import cis_trans_frame, kind_name, set_parity_by_probe

# A dative bond points donor -> acceptor in RDKit and chython's graph is undirected, so the direction
# is reconstructed on the way out: main-group ligand first, metal second.
_MAIN_GROUP = frozenset({1, 2, 5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18,
                         32, 33, 34, 35, 36, 51, 52, 53, 54})

# Order 8 is chython's "unspecified or dative"; DATIVE is the closest RDKit type surviving a round trip.
_ORDERS = None  # built on first use, since naming a BondType means importing RDKit


def _order_tables():
    """The two bond-order tables, built once, on the first conversion.

    Lazy because a module-level dict keyed by `BondType` would import RDKit at import time, which
    the dispatch test forbids.
    """
    global _ORDERS
    if _ORDERS is None:
        from rdkit.Chem import BondType

        out = {1: BondType.SINGLE, 2: BondType.DOUBLE, 3: BondType.TRIPLE, 4: BondType.AROMATIC,
               8: BondType.DATIVE}
        # Every RDKit type chython can read, and only those.  Anything else is reported and stored as
        # order 8, not as a single bond, which would be a statement the source never made.
        back = {BondType.SINGLE: 1, BondType.DOUBLE: 2, BondType.TRIPLE: 3, BondType.AROMATIC: 4,
                BondType.ZERO: 8, BondType.UNSPECIFIED: 8, BondType.DATIVE: 8,
                BondType.DATIVEONE: 8, BondType.DATIVEL: 8, BondType.DATIVER: 8}
        _ORDERS = (out, back)
    return _ORDERS


def _note(log, line):
    """Append one loss to a record list -- the caller's on export, the import's own on the way in."""
    if log is not None:
        log.append(LogRecord('rdkit:note', (), line))


def _bump(counts, name):
    """Count one loss of a named kind, so the log gets one line per kind and not one per atom."""
    counts[name] = counts.get(name, 0) + 1


def to_rdkit(mol, /, *, log=None, keep_mapping=True, keep_numbers=False, keep_hydrogens=True,
             keep_coordinates=None, absolute=False):
    """
    Export a chython container as an RDKit `Mol`, or a reaction as a `ChemicalReaction`.

    Atoms are emitted in the molecule's own iteration order, which is part of the interface:
    `depict/layout/molecule.py` maps coordinates back by position.  Aromatic bonds go out as
    `BondType.AROMATIC`; nothing here kekulizes.

    :param log: optional list receiving one human-readable line per reportable loss.
    :param keep_mapping: carry the atom-atom mapping across -- RDKit's atom map number is set from
        `map_number`, and an atom that carries none gets none.  An unmapped molecule therefore
        exports clean, which is what makes `MolToSmiles` of the result a plain SMILES.
    :param keep_numbers: put chython's stable atom *ids* in RDKit's map field instead.  A different
        question from the one above, and the reason the two are separate flags: this is a label to
        match results back on, not a mapping.  RDKit has ONE integer per atom, so this wins the field
        and a `map_number` it displaces is logged rather than guessed at.
    :param keep_hydrogens: state the implicit-H count on the RDKit atom (`SetNumExplicitHs` plus
        `SetNoImplicit`, so the stated count is final) rather than letting RDKit perceive one.
    :param keep_coordinates: export the 2D layout as a conformer.  `None` exports it only when a
        layout exists.
    :param absolute: add an ABS stereo group for stereocentres not in an AND/OR group.
    """
    if isinstance(mol, V3Reaction):
        return _reaction_to_rdkit(mol, log=log, keep_mapping=keep_mapping, keep_numbers=keep_numbers,
                                  keep_hydrogens=keep_hydrogens, keep_coordinates=keep_coordinates,
                                  absolute=absolute)

    from rdkit.Chem import Atom, Conformer, RWMol, SanitizeFlags, SanitizeMol
    from rdkit.Chem.rdmolops import AssignStereochemistry, FastFindRings, \
        SetDoubleBondNeighborDirections

    if not isinstance(mol, V3Molecule):
        raise UnconvertibleType(f'{type(mol).__name__} has no RDKit form: only a molecule and a '
                                f'reaction container convert, and a query is neither')

    forward, _ = _order_tables()
    # (number, atomic number, charge, isotope, radical, implicit H or None, map number, xy)
    atoms = [(a.n, a.element, a.charge, a.isotope, a.is_radical, a.implicit_h,
              a.map_number, a.xy) for a in mol.atoms()]
    bonds = [(b.n, b.m, b.order) for b in mol.bonds()]

    # Aromatic from the stored bond order, not from `hybridization`, which is a derived cache.
    aromatic = {n for x, y, o in bonds if o == 4 for n in (x, y)}

    rw = RWMol()
    index = {}  # chython number -> RDKit index.  Insertion order is the molecule's iteration order.
    unknown_h = 0
    shadowed_maps = 0
    for n, number, charge, isotope, radical, hydrogens, map_number, _ in atoms:
        ra = Atom(number)
        if charge:
            ra.SetFormalCharge(charge)
        if isotope:
            ra.SetIsotope(isotope)
        if radical:
            ra.SetNumRadicalElectrons(1)
        if keep_numbers:
            ra.SetAtomMapNum(n)
            if keep_mapping and map_number and map_number != n:
                shadowed_maps += 1
        elif keep_mapping and map_number:
            ra.SetAtomMapNum(map_number)
        if keep_hydrogens:
            if hydrogens is None:
                # H_UNKNOWN: RDKit cannot spell "nobody stated a count" and a stated zero is a
                # different molecule, so the count is left for RDKit to perceive, and logged.
                unknown_h += 1
            else:
                ra.SetNumExplicitHs(hydrogens)
                ra.SetNoImplicit(True)
        if n in aromatic:
            ra.SetIsAromatic(True)
        index[n] = rw.AddAtom(ra)

    if unknown_h:
        _note(log, f'{unknown_h} atom(s) have no implicit hydrogen count (H_UNKNOWN); RDKit has no '
                   f'spelling for that and will perceive one instead')
    if shadowed_maps:
        _note(log, f'{shadowed_maps} atom(s) carry a map number different from their atom number; '
                   f'keep_numbers wrote the atom number, so the mapping is not in the RDKit '
                   f'molecule')

    elements = {a[0]: a[1] for a in atoms}
    for n, m, order in bonds:
        if order == 8 and elements[n] not in _MAIN_GROUP:
            n, m = m, n  # a dative bond points from the donor to the acceptor
        rb = rw.AddBond(index[n], index[m], forward[order])
        if order == 4:
            rw.GetBondWithIdx(rb - 1).SetIsAromatic(True)

    reverse = {v: k for k, v in index.items()}
    _v3_export_stereo(mol, rw, index, reverse, log)
    _export_stereo_groups(mol, rw, index, absolute, log)

    if keep_coordinates is None:
        keep_coordinates = any(xy and (xy[0] or xy[1]) for *_, xy in atoms)
    if keep_coordinates:
        conf = Conformer(len(atoms))
        for n, *_, xy in atoms:
            x, y = xy or (0., 0.)
            conf.SetAtomPosition(index[n], (x, y, 0.))
        conf.Set3D(False)
        rw.AddConformer(conf, assignId=True)

    # THE GEOMETRY CONFORMERS FOLLOW THE LAYOUT, so the layout keeps id 0, which is what
    # `GetConformer()` returns.  `keep_coordinates` does not gate them: that flag is about a depiction,
    # not a geometry.  One RDKit conformer per model, in model order.
    for model in mol.conformers:
        conf = Conformer(len(atoms))
        for n, *_ in atoms:
            conf.SetAtomPosition(index[n], model.xyz_of(n))
        conf.Set3D(True)
        rw.AddConformer(conf, assignId=True)

    # Not plain `SanitizeMol(rw)`: its default ops rewrite what was just written -- KEKULIZE and
    # SETAROMATICITY rewrite the representation, CLEANUP rewrites nitro and azide, CLEANUPCHIRALITY
    # drops the chiral tags, FINDRADICALS recomputes the radical counts, ADJUSTHS moves hydrogens
    # between the explicit and implicit sides.  Only the four cache-filling ops are safe.
    ops = (SanitizeFlags.SANITIZE_PROPERTIES | SanitizeFlags.SANITIZE_SYMMRINGS
           | SanitizeFlags.SANITIZE_SETCONJUGATION | SanitizeFlags.SANITIZE_SETHYBRIDIZATION)
    try:
        SanitizeMol(rw, sanitizeOps=ops)
    except Exception as e:
        # A valence RDKit will not accept is not a reason to hand back nothing; only the caches are
        # lost, and that is logged.
        try:
            rw.UpdatePropertyCache(strict=False)
            FastFindRings(rw)
        except Exception as e2:  # pragma: no cover -- RDKit refusing even the permissive path
            raise ToolkitError(f'RDKit refused the molecule: {e2}') from e2
        _note(log, f'RDKit sanitization failed ({e}); valence and ring caches are approximate')
    else:
        AssignStereochemistry(rw, cleanIt=False, force=True, flagPossibleStereoCenters=True)
        # RDKit derives the SMILES bond directions from the stereo atoms only for single-fragment
        # molecules; without this any salt or solvate loses its cis/trans marks on SMILES export.
        SetDoubleBondNeighborDirections(rw)
    return rw


def _v3_export_stereo(mol, rw, index, reverse, log):
    """Tetrahedral and cis/trans configurations of a V3 molecule onto the RDKit molecule.

    Two measured calibrations: RDKit lists an atom's directions in bond order with the implicit
    hydrogen appended *last*, and `CHI_TETRAHEDRAL_CCW` is SMILES `@` over that list; and
    `translate_stereo` answers 2 for an odd permutation parity, which is `@`, hence CCW.
    """
    from rdkit.Chem import BondStereo, ChiralType

    dropped = {}
    wedged_unsigned = 0
    wedged = {n for n, m, _ in mol.wedges()}
    for u in mol.stereo_units():
        kind = u['kind']
        anchor = u['anchor']
        if not u['parity']:
            # A drawn centre with no derived parity: RDKit takes tags, not wedges, so the
            # configuration in the drawing does not reach it.
            if u['stereogenic'] and anchor in wedged:
                wedged_unsigned += 1
            continue
        if kind == SU_TETRA:
            ra = rw.GetAtomWithIdx(index[anchor])
            order = [reverse[x.GetIdx()] for x in ra.GetNeighbors()]
            order += [None] * (4 - len(order))
            try:
                parity = mol.translate_stereo(anchor, tuple(order))
            except (KeyError, ValueError) as e:
                _bump(dropped, f'tetrahedral centre RDKit will not frame ({e})')
                continue
            ra.SetChiralTag(ChiralType.CHI_TETRAHEDRAL_CCW if parity == 2 else
                            ChiralType.CHI_TETRAHEDRAL_CW)
        elif kind == SU_CIS_TRANS:
            frame = cis_trans_frame(mol, anchor, u['refs'])
            if frame is None:  # pragma: no cover -- perception refuses such a unit
                _bump(dropped, 'cis/trans bond that cannot be framed on the graph')
                continue
            near, far, partner = frame.near, frame.far, frame.partner
            try:
                parity = mol.translate_stereo(anchor, frame.order)
            except (KeyError, ValueError) as e:
                _bump(dropped, f'cis/trans bond RDKit will not frame ({e})')
                continue
            rb = rw.GetBondBetweenAtoms(index[anchor], index[partner])
            if rb.GetBeginAtomIdx() == index[anchor]:
                rb.SetStereoAtoms(index[near], index[far])
            else:
                rb.SetStereoAtoms(index[far], index[near])
            # parity 2 is "the two chosen directions are on the same side" = RDKit STEREOZ over the
            # same two stereo atoms.
            rb.SetStereo(BondStereo.STEREOZ if parity == 2 else BondStereo.STEREOE)
        else:
            _bump(dropped, kind_name(kind))

    for name, count in dropped.items():
        _note(log, f'{count} {name} configuration(s) dropped: RDKit has no form for them')
    if wedged_unsigned:
        _note(log, f'{wedged_unsigned} stereocentre(s) are drawn with a wedge but carry no derived '
                   f'parity; RDKit takes configurations and not wedges, so those are not exported')


def _export_stereo_groups(mol, rw, index, absolute, log):
    """AND / OR / ABS stereo groups onto `Chem.StereoGroup`s.

    `absolute=True` additionally names the centres nothing else claims.
    """
    from rdkit.Chem import CreateStereoGroup, StereoGroupType

    signed = _configured_tetrahedral(mol)
    rac, rel, ast = {}, {}, []
    for (kind, gid), members in mol.stereo_groups().items():
        members = [index[n] for n in members if n in signed]
        if not members:
            continue
        if kind == STEREO_AND:
            rac.setdefault(gid, []).extend(members)
        elif kind == STEREO_OR:
            rel.setdefault(gid, []).extend(members)
        elif kind == STEREO_ABS:
            ast.extend(members)
    claimed = {n for ms in (*rac.values(), *rel.values()) for n in ms} | set(ast)
    if absolute:
        ast.extend(index[n] for n in signed if index[n] not in claimed)

    groups = []
    if ast:
        groups.append((StereoGroupType.STEREO_ABSOLUTE, 0, sorted(set(ast))))
    groups.extend((StereoGroupType.STEREO_AND, gid, rac[gid]) for gid in sorted(rac))
    groups.extend((StereoGroupType.STEREO_OR, gid, rel[gid]) for gid in sorted(rel))
    if groups:
        sgs = []
        for gt, gid, members in groups:
            sg = CreateStereoGroup(gt, rw, members, [], gid)
            sg.SetWriteId(gid)  # otherwise RDKit renumbers the groups from one on write
            sgs.append(sg)
        rw.SetStereoGroups(sgs)


def _configured_tetrahedral(mol):
    """Atom numbers whose tetrahedral configuration was actually written to the RDKit molecule.

    RDKit drops a stereo group over a centre with no chiral tag, so group membership is filtered by
    what the export produced.
    """
    return {u['anchor'] for u in mol.stereo_units() if u['parity'] and u['kind'] == SU_TETRA}


def _reaction_to_rdkit(rxn, /, *, log=None, **kwargs):
    """A `ReactionContainer` as a `rdChemReactions.ChemicalReaction`.

    THE THREE SIDES GO TO THE THREE TEMPLATE LISTS, in `molecules()` order, so a round trip through
    here does not move a catalyst onto the left.  RDKit's own word for the middle side is "agent",
    which is chython's, so nothing is renamed on the way.

    A CHEMICAL REACTION AND NOT THREE LISTS OF `Mol`, because the atom-atom mapping is the point: this
    is the one RDKit object whose SMILES writer emits `:n` labels across an arrow, so `keep_mapping`
    means something here that it cannot mean for a lone molecule.  It is RDKit's transform type as
    well as its record type, and a record put in it is not thereby a transform -- `Initialize()` is
    left to a caller who wants to RUN it.

    `log` is shared by every side, so one list reports the whole record; each molecule's losses are
    the same ones `to_rdkit` reports of it alone.
    """
    from rdkit.Chem.rdChemReactions import ChemicalReaction

    rr = ChemicalReaction()
    for add, side in ((rr.AddReactantTemplate, rxn.reactants), (rr.AddAgentTemplate, rxn.agents),
                      (rr.AddProductTemplate, rxn.products)):
        for m in side:
            add(to_rdkit(m, log=log, **kwargs))
    return rr


def from_rdkit(data, /, *, log=None):
    """
    Import an RDKit `Mol` or `RWMol` as a chython molecule, or a `ChemicalReaction` as a reaction.

    Bond order comes from `GetBondType()`, so a kekulized molecule stays kekulized and an aromatic
    one stays aromatic; `GetIsAromatic()` is deliberately not consulted.  RDKit atom map numbers
    become V3 `map_number`s, not stable ids.

    EVERY RECORD LANDS ON THE RETURNED CONTAINER'S `.log`, in stage `'interop'`, with nothing passed
    in: this direction produces the container that owns them.  `log`, when given, receives a copy of
    the same records, which is how a caller reading a whole file keeps one sequence for the batch.
    """
    # a `ChemicalReaction` first, because it is recognised by the templates it holds and has no
    # `GetAtoms` of its own -- so the molecule test below would refuse it.
    if hasattr(data, 'GetReactants') and hasattr(data, 'GetProducts'):
        return _reaction_from_rdkit(data, log=log)
    if not hasattr(data, 'GetAtoms') or not hasattr(data, 'GetBonds'):
        raise UnconvertibleType(f'{type(data).__name__} is neither a chython container nor an RDKit '
                                f'Mol/RWMol/ChemicalReaction')
    from rdkit.Chem import StereoGroupType

    _, back = _order_tables()
    # Collected here and handed to `deliver` at the end, because the first record is emitted before
    # there is a molecule to put it on.
    records = []
    if not _valence_cache(data):
        # An `RWMol` built atom by atom has no valence cache, and `GetNumImplicitHs` on it raises a
        # precondition violation rather than answering.  Fill it permissively and report.
        try:
            data.UpdatePropertyCache(strict=False)
        except Exception as e:
            raise ToolkitError(f'RDKit molecule has no valence information and refused to compute '
                               f'it: {e}') from e
        _note(records, 'the RDKit molecule had no valence cache; its implicit hydrogen counts are '
                       'RDKit\'s perception rather than a stated count')

    mol = V3Molecule()
    number = {}  # RDKit index -> chython stable id
    clamped_charges = clamped_h = radicals = dropped_maps = 0
    unknown_orders = {}
    try:
        with mol.edit():
            for ra in data.GetAtoms():
                charge = ra.GetFormalCharge()
                if charge < -4 or charge > 8:
                    clamped_charges += 1
                    charge = -4 if charge < -4 else 8
                electrons = ra.GetNumRadicalElectrons()
                if electrons > 1:
                    radicals += 1
                hydrogens = ra.GetNumExplicitHs() + ra.GetNumImplicitHs()
                if hydrogens > H_IMPLICIT_MAX:
                    clamped_h += 1
                    hydrogens = H_IMPLICIT_MAX
                map_number = ra.GetAtomMapNum()
                if map_number > 9999:
                    dropped_maps += 1
                    map_number = 0
                number[ra.GetIdx()] = mol.add_atom(
                    ra.GetAtomicNum(), charge=charge, isotope=ra.GetIsotope(),
                    radical=bool(electrons), map_number=map_number, implicit_h=hydrogens)
            for rb in data.GetBonds():
                bt = rb.GetBondType()
                order = back.get(bt)
                if order is None:
                    unknown_orders[str(bt)] = unknown_orders.get(str(bt), 0) + 1
                    order = 8
                mol.add_bond(number[rb.GetBeginAtomIdx()], number[rb.GetEndAtomIdx()], order)
    except Exception as e:
        raise ToolkitError(f'RDKit molecule could not be read: {e}') from e

    if clamped_charges:
        _note(records, f'{clamped_charges} formal charge(s) outside -4..8 clamped to the range '
                       f'chython stores')
    if radicals:
        _note(records, f'{radicals} atom(s) carry more than one radical electron; chython stores a '
                       f'single radical flag, so the count is lost')
    if clamped_h:
        _note(records, f'{clamped_h} hydrogen count(s) above {H_IMPLICIT_MAX} clamped: 15 is the '
                       f'"not stated" sentinel and cannot also be a count')
    if dropped_maps:
        _note(records, f'{dropped_maps} atom map number(s) above 9999 dropped')
    for name, count in unknown_orders.items():
        _note(records, f'{count} bond(s) of type {name} stored as order 8 (unspecified): chython has '
                       f'no order for that type')

    _import_stereo(data, mol, number, records)
    _import_stereo_groups(data, mol, number, StereoGroupType, records)
    _import_conformers(data, mol, number, records)
    deliver(mol, records, log)
    return mol


def _valence_cache(data):
    """Has this molecule's implicit valence been computed?  Asked of one atom, because RDKit computes
    the cache for the whole molecule at once and an empty molecule needs nothing."""
    if not data.GetNumAtoms():
        return True
    try:
        data.GetAtomWithIdx(0).GetNumImplicitHs()
    except RuntimeError:
        return False
    return True


def _import_stereo(data, mol, number, log):
    """RDKit chiral tags and double-bond stereo into V3 parities.

    The write is done by probing rather than by arithmetic: write parity 1, ask `translate_stereo`
    what that looks like in RDKit's order, keep it or flip it.  Restating the permutation-parity rule
    here would risk a sign error, i.e. a silently mirrored molecule.
    """
    from rdkit.Chem import BondStereo, ChiralType

    tags = {}
    for ra in data.GetAtoms():
        tag = ra.GetChiralTag()
        if tag == ChiralType.CHI_TETRAHEDRAL_CW or tag == ChiralType.CHI_TETRAHEDRAL_CCW:
            tags[number[ra.GetIdx()]] = ([number[x.GetIdx()] for x in ra.GetNeighbors()],
                                         2 if tag == ChiralType.CHI_TETRAHEDRAL_CCW else 1)
        elif tag != ChiralType.CHI_UNSPECIFIED:
            _note(log, f'atom {ra.GetIdx()} chiral tag {tag} is not tetrahedral and was dropped')

    cis = {BondStereo.STEREOZ, BondStereo.STEREOCIS}
    trans = {BondStereo.STEREOE, BondStereo.STEREOTRANS}
    signs = []
    for rb in data.GetBonds():
        s = rb.GetStereo()
        if s in cis or s in trans:
            a, b = rb.GetStereoAtoms()
            signs.append((number[rb.GetBeginAtomIdx()], number[rb.GetEndAtomIdx()],
                          number[a], number[b], 2 if s in cis else 1))
        elif s == BondStereo.STEREOANY:
            _note(log, f'bond {rb.GetBeginAtomIdx()}-{rb.GetEndAtomIdx()} is marked as unknown '
                       f'geometry; chython has no "either" double bond and it was dropped')

    if not tags and not signs:
        return

    units = {u['anchor']: u for u in mol.stereo_units()}
    dropped = 0
    for n, (env, want) in tags.items():
        u = units.get(n)
        if u is None or u['kind'] != SU_TETRA:
            dropped += 1
            continue
        order = tuple(env + [None] * (4 - len(env)))
        if not set_parity_by_probe(mol, n, order, want):
            dropped += 1
    for n, m, a, b, want in signs:
        for anchor, near, far in ((n, a, b), (m, b, a)):
            u = units.get(anchor)
            if u is None or u['kind'] != SU_CIS_TRANS:
                continue
            order = _bond_order_for(u['refs'], near, far)
            if order is not None and set_parity_by_probe(mol, anchor, order, want):
                break
        else:
            dropped += 1
    if dropped:
        _note(log, f'{dropped} configuration(s) from RDKit could not be placed on a chython stereo '
                   f'unit and were dropped')


def _bond_order_for(refs, near, far):
    """`refs` permuted so that `near` and `far` head the two halves, or None if they are not in it.

    `translate_stereo` requires each half of the order to map onto one stored half, and allows the
    two halves to be exchanged wholesale.
    """
    a, b = refs[:2], refs[2:]
    if near in a and far in b:
        head, tail = a, b
    elif near in b and far in a:
        head, tail = b, a  # the legal wholesale exchange of the two halves
    else:
        return None
    return (near, head[1] if head[0] == near else head[0],
            far, tail[1] if tail[0] == far else tail[0])


def _import_stereo_groups(data, mol, number, StereoGroupType, log):
    """`Chem.StereoGroup`s into V3's stored AND / OR / ABS marks.

    RDKit group ids are free integers and V3's are 1..63 per kind, so an id is kept whenever it fits
    and only a clash or an out-of-range value moves one (which is logged).
    """
    groups = data.GetStereoGroups()
    if not groups:
        return
    kinds = {StereoGroupType.STEREO_AND: STEREO_AND, StereoGroupType.STEREO_OR: STEREO_OR,
             StereoGroupType.STEREO_ABSOLUTE: STEREO_ABS}
    marked = set()
    taken = {STEREO_AND: set(), STEREO_OR: set()}
    renumbered = overflow = 0
    with mol.edit():
        for sg in groups:
            kind = kinds.get(sg.GetGroupType())
            if kind is None:  # pragma: no cover -- RDKit has exactly these three today
                _note(log, f'stereo group of type {sg.GetGroupType()} dropped: chython stores '
                           f'absolute, AND and OR only')
                continue
            if kind is STEREO_ABS:
                gid = 0  # V3 spells the absolute kind with no id; there is only ever one of it
            else:
                used = taken[kind]
                gid = sg.GetReadId()
                if gid < 1 or gid > 63 or gid in used:
                    gid = next((i for i in range(1, 64) if i not in used), 0)
                    if not gid:
                        overflow += 1
                        continue
                    renumbered += 1
                used.add(gid)
            for ra in sg.GetAtoms():
                n = number[ra.GetIdx()]
                if n in marked:  # RDKit allows an atom in two groups at once; the first one wins
                    continue
                marked.add(n)
                mol.set_stereo_group(n, kind, gid)
    if renumbered:
        _note(log, f'{renumbered} stereo group id(s) renumbered: chython stores 1..63 per kind')
    if overflow:
        _note(log, f'{overflow} stereo group(s) dropped: chython stores at most 63 per kind')


def _import_conformers(data, mol, number, log):
    """The first 2D conformer becomes the layout and EVERY 3D one becomes a model.

    One layout because a molecule has one drawing; N models because the arena holds N.  Neither
    substitutes for the other: the xy of a 3D conformer is a projection rather than a layout, and the z
    of a 2D one is zero.  A molecule with both gets both, `set_xy` and `set_xyz` being independent by
    design.  Each model carries the RDKit conformer id as its `ext_index`.
    """
    plane = None
    planes = 0
    solids = []
    for c in data.GetConformers():
        if c.Is3D():
            solids.append(c)
        else:
            planes += 1
            if plane is None:
                plane = c
    if planes > 1:
        _note(log, f'{planes - 1} of {planes} 2D conformer(s) dropped: a molecule has one layout')
    if plane is None and not solids:
        return
    flat = plane.GetPositions() if plane is not None else None
    # ONE SESSION FOR ALL OF IT: adds and sets may share a scope, so N models cost one respan.
    with mol.edit():
        if flat is not None:
            for ra in data.GetAtoms():
                mol.set_xy(number[ra.GetIdx()], flat[ra.GetIdx()][0], flat[ra.GetIdx()][1])
        for c in solids:
            full = c.GetPositions()
            model = mol.add_conformer(ext_index=c.GetId())
            for ra in data.GetAtoms():
                x, y, z = full[ra.GetIdx()]
                mol.set_xyz(number[ra.GetIdx()], x, y, z, model=model)


def _reaction_from_rdkit(data, /, *, log=None):
    """A `ChemicalReaction` as a `ReactionContainer`, one side at a time.

    RDKit atom map numbers become `map_number`s, as they do for a lone molecule -- which is what makes
    the mapping the thing that survives the round trip, and the reason a reaction goes through this
    type rather than through three separate molecule conversions the caller reassembles.

    A TEMPLATE IS READ AS A RECORD.  A reaction parsed from SMARTS holds query features chython has no
    atom for, and each of those raises out of the molecule conversion; nothing here inspects the
    templates first to refuse more politely, because the honest answer is the one the molecule
    conversion already gives.

    EACH MOLECULE KEEPS ITS OWN RECORDS AND THE REACTION GETS A COPY stamped with `'agents[0]'` and
    the like, since a record's `atoms` are stable ids in one molecule and mean nothing pooled.
    """
    reactants = [from_rdkit(m, log=log) for m in data.GetReactants()]
    products = [from_rdkit(m, log=log) for m in data.GetProducts()]
    agents = [from_rdkit(m, log=log) for m in data.GetAgents()]
    rxn = V3Reaction(reactants, products, agents)
    for side, molecules in (('reactants', reactants), ('agents', agents), ('products', products)):
        for i, m in enumerate(molecules):
            mirror(rxn, f'{side}[{i}]', m.log)
    return rxn


__all__ = ['from_rdkit', 'to_rdkit']
