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
OpenBabel conversion, both directions.  Reached through `interop.openbabel`.

Never call EndModify(): after CloneData() attaches an OBTetrahedralStereo or OBCisTransStereo, it
re-perceives the whole structure and discards every stereo object just set.  SetAromaticPerceived(True)
plus SetChiralityPerceived(True) stop OBMol perceiving on its own.  CloneData is the Python entry point
for attaching stereo data; SetData() is C++-only.
"""
from ..core import LogRecord, LOST, REPAIRED, SU_ALLENE, SU_CIS_TRANS, SU_TETRA
from ..exceptions import UnconvertibleType
from ._records import deliver
from ._stereo import cis_trans_partner, kind_name, set_parity_by_probe


def to_openbabel(mol, /, *, log=None):
    """
    Export a `chython.core.MoleculeContainer` as an OpenBabel OBMol.

    Atoms are emitted in the molecule's own iteration order, so OBMol's 1-based index i is the i-th
    atom of `mol.atoms()`; that is part of the interface, since ``depict/layout/molecule.py`` maps 2D
    coordinates back by position.  Nothing kekulizes: order 4 goes out as bond order 1 with the
    aromatic flag set.  Losses: dative bonds (order 8) become single bonds, H_UNKNOWN becomes a stated
    zero, allene/cumulene stereo is not exported (OBMol 3.1.0 exposes no OBExtendedTetrahedralStereo
    in its Python bindings) and neither are atropisomers.

    :param log: optional list; one human-readable str is appended per reportable loss.
    :raises UnconvertibleType: if mol is not a chython MoleculeContainer.
    """
    from openbabel import openbabel as ob
    from ..core import MoleculeContainer as V3

    if not isinstance(mol, V3):
        raise UnconvertibleType(
            f'to_openbabel expects a chython MoleculeContainer, got {type(mol).__qualname__}'
        )

    # OBStereo.ImplicitRef is -2 as a signed C long; mask to unsigned for Python comparisons.
    implicit_ref = ob.OBStereo.ImplicitRef & 0xFFFFFFFFFFFFFFFF

    ob_mol = ob.OBMol()
    # Do not call BeginModify()/EndModify(): EndModify() re-perceives the whole structure and destroys
    # the stereo data attached via CloneData().  SetAromaticPerceived(True) is the right guard.
    idx = {}  # chython atom id → OBMol 1-based GetIdx() (for AddBond)
    ids = {}  # chython atom id → OBMol atom GetId() (for stereo refs)

    for atom in mol.atoms():
        n = atom.n
        oa = ob_mol.NewAtom()
        oa.SetAtomicNum(atom.element)
        if atom.charge:
            oa.SetFormalCharge(atom.charge)
        if atom.isotope:
            oa.SetIsotope(atom.isotope)
        if atom.is_radical:
            oa.SetSpinMultiplicity(2)
        h = atom.implicit_h
        if h is None:
            # OBMol cannot distinguish "zero" from "unstated": both read as 0 ImplicitHCount.
            if log is not None:
                log.append(LogRecord('openbabel:h-unknown-treated-as-zero', (n,),
                                     f'atom {n} (element {atom.element}): implicit hydrogen count is '
                                     f'H_UNKNOWN; OpenBabel cannot distinguish "zero" from "unstated", '
                                     f'treated as zero in OBMol', REPAIRED))
            # leave OBMol at its default (0); do not call SetImplicitHCount
        else:
            oa.SetImplicitHCount(h)
        idx[n] = oa.GetIdx()
        ids[n] = oa.GetId()

    for bond in mol.bonds():
        n, m, order = bond.n, bond.m, bond.order
        if order == 8:
            if log is not None:
                log.append(LogRecord('openbabel:dative-bond-stored-as-single', (n, m),
                                     f'bond {n}-{m}: dative/order-8 bond has no OpenBabel equivalent; '
                                     f'stored as single bond in OBMol', REPAIRED))
            ob_order = 1
        elif order == 4:
            ob_order = 1  # aromatic: single + aromatic flag set below
        else:
            ob_order = order
        ob_mol.AddBond(idx[n], idx[m], ob_order)
        if order == 4:
            ob_mol.GetBond(idx[n], idx[m]).SetAromatic(True)
            ob_mol.GetAtom(idx[n]).SetAromatic(True)
            ob_mol.GetAtom(idx[m]).SetAromatic(True)

    ob_mol.SetAromaticPerceived(True)

    _export_v3_stereo(mol, ob_mol, ids, implicit_ref, log)

    ob_mol.SetChiralityPerceived(True)
    return ob_mol


def _ob_ref(ids, chy_id, implicit_ref):
    """Map a chython atom number (or None for an unnamed direction) to an OBMol atom id."""
    return implicit_ref if chy_id is None else ids[chy_id]


def _export_v3_stereo(mol, ob_mol, ids, implicit_ref, log):
    from openbabel import openbabel as ob

    for unit in mol.stereogenic_units():
        kind = unit['kind']
        anchor = unit['anchor']
        parity = unit['parity']
        if parity == 0:
            continue  # no configuration stored for this unit

        refs = unit['refs']

        if kind == SU_TETRA:
            # parity 2 (odd, `@`) = AntiClockwise, parity 1 (even, `@@`) = Clockwise; OB takes
            # from_or_towards = refs[0] and MakeRefs(refs[1], refs[2], refs[3]).
            config = ob.OBTetrahedralConfig()
            config.center = ids[anchor]
            config.from_or_towards = _ob_ref(ids, refs[0], implicit_ref)
            config.view = ob.OBStereo.ViewFrom
            config.refs = ob.OBStereo.MakeRefs(
                _ob_ref(ids, refs[1], implicit_ref),
                _ob_ref(ids, refs[2], implicit_ref),
                _ob_ref(ids, refs[3], implicit_ref),
            )
            config.winding = ob.OBStereo.AntiClockwise if parity == 2 else ob.OBStereo.Clockwise
            config.specified = True
            ts = ob.OBTetrahedralStereo(ob_mol)
            ts.SetConfig(config)
            ob_mol.CloneData(ts)  # SetData is C++-only; CloneData is the Python entry point

        elif kind == SU_CIS_TRANS:
            # refs = (anchor_dir0, anchor_dir1, partner_dir0, partner_dir1); parity 2 = cis = refs[0]
            # and refs[2] on the same side.  In OBMol ShapeU, positions 0 and 3 are the same side and
            # 0 and 2 are opposite, so cis is MakeRefs(a0, a1, b1, b0) (partner pair swapped) and trans
            # is MakeRefs(a0, a1, b0, b1).  Only the partner lookup is shared with `_stereo`.
            partner = cis_trans_partner(mol, anchor, refs)
            if partner is None:
                if log is not None:
                    log.append(LogRecord('openbabel:cis-trans-frame-not-found', (anchor,),
                                         f'atom {anchor}: cis/trans stereo cannot be framed on the graph; '
                                         f'not exported', LOST))
                continue

            a0 = _ob_ref(ids, refs[0], implicit_ref)
            a1 = _ob_ref(ids, refs[1], implicit_ref)
            b0 = _ob_ref(ids, refs[2], implicit_ref)
            b1 = _ob_ref(ids, refs[3], implicit_ref)

            config = ob.OBCisTransConfig()
            config.begin = ids[anchor]
            config.end = ids[partner]
            if parity == 2:  # cis: anchor_dir0 same side as partner_dir0 → positions 0 and 3
                config.refs = ob.OBStereo.MakeRefs(a0, a1, b1, b0)
            else:            # trans: anchor_dir0 opposite side from partner_dir0 → positions 0 and 2
                config.refs = ob.OBStereo.MakeRefs(a0, a1, b0, b1)
            config.shape = ob.OBStereo.ShapeU
            config.specified = True
            cts = ob.OBCisTransStereo(ob_mol)
            cts.SetConfig(config)
            ob_mol.CloneData(cts)

        elif kind == SU_ALLENE:
            if log is not None:
                log.append(LogRecord('openbabel:allene-not-representable', (anchor,),
                                     f'atom {anchor}: allene/cumulene stereo not exported -- OBMol 3.1.0 has '
                                     f'no OBExtendedTetrahedralStereo in its Python bindings', LOST))

        elif log is not None:
            # atropisomers, and whatever kind the core grows next: named from the core's own table so
            # a new one cannot arrive unreported
            log.append(LogRecord('openbabel:stereo-kind-not-representable', (anchor,),
                                 f'atom {anchor}: {kind_name(kind)} stereo has no OBMol representation; '
                                 f'not exported', LOST))


def from_openbabel(data, /, *, log=None):
    """
    Import an OpenBabel OBMol as a V3 chython MoleculeContainer.

    Stores what OpenBabel reported; aromaticity is not re-perceived (order 4 when OBMol marks the bond
    IsAromatic(), the numeric order otherwise).  Implicit hydrogen counts come literally from
    GetImplicitHCount(), which is always numeric (0 when unset), so H_UNKNOWN is not recoverable and
    nothing is logged -- "unset" and "zero" are indistinguishable.  Tetrahedral stereo comes through
    OBStereoFacade; cis/trans through OBCisTransStereo.IsCis() on each side's primary substituent,
    which is symmetric in ``begin``/``end``.  Allene stereo is not imported (no binding for it).

    EVERY RECORD LANDS ON THE RETURNED MOLECULE'S `.log`, in stage `'interop'`, with nothing passed
    in; `log`, when given, receives a copy of the same records.

    :param log: optional list receiving a copy of the records put on the molecule.
    :raises UnconvertibleType: if data is not an openbabel.OBMol.
    """
    from openbabel import openbabel as ob
    if not isinstance(data, ob.OBMol):
        raise UnconvertibleType(
            f'from_openbabel expects an openbabel.OBMol, got {type(data).__qualname__}'
        )

    from ..core import MoleculeContainer

    implicit_ref = ob.OBStereo.ImplicitRef & 0xFFFFFFFFFFFFFFFF

    mol = MoleculeContainer()
    ob_id_to_chy = {}  # OBMol atom GetId() → chython atom number
    records = []

    with mol.edit():
        for i in range(1, data.NumAtoms() + 1):
            atom = data.GetAtom(i)
            sid = mol.add_atom(
                atom.GetAtomicNum(),
                charge=atom.GetFormalCharge(),
                isotope=atom.GetIsotope(),
                radical=atom.GetSpinMultiplicity() != 0,
                implicit_h=atom.GetImplicitHCount(),
            )
            ob_id_to_chy[atom.GetId()] = sid

        for i in range(data.NumBonds()):
            bond = data.GetBond(i)
            n_chy = ob_id_to_chy[data.GetAtom(bond.GetBeginAtomIdx()).GetId()]
            m_chy = ob_id_to_chy[data.GetAtom(bond.GetEndAtomIdx()).GetId()]
            order = 4 if bond.IsAromatic() else bond.GetBondOrder()
            mol.add_bond(n_chy, m_chy, order)

    # Stereo must be configured outside the edit scope: set_parity() calls _require_clean().
    facade = ob.OBStereoFacade(data)

    for i in range(1, data.NumAtoms() + 1):
        atom = data.GetAtom(i)
        ob_id = atom.GetId()
        if not facade.HasTetrahedralStereo(ob_id):
            continue
        td = facade.GetTetrahedralStereo(ob_id)
        cfg = td.GetConfig()
        if not cfg.specified:
            continue

        center_chy = ob_id_to_chy.get(ob_id)
        if center_chy is None:
            continue

        def _to_chy(oid, _ir=implicit_ref, _m=ob_id_to_chy):
            return None if oid == _ir else _m.get(oid)

        # OBMol order is (from_or_towards, refs[0], refs[1], refs[2]); AntiClockwise = `@` = parity 2,
        # Clockwise = `@@` = parity 1.
        order_tuple = (
            _to_chy(cfg.from_or_towards),
            _to_chy(cfg.refs[0]),
            _to_chy(cfg.refs[1]),
            _to_chy(cfg.refs[2]),
        )
        target_parity = 2 if cfg.winding == ob.OBStereo.AntiClockwise else 1

        # a rejected frame leaves the centre cleared and answers False
        if not set_parity_by_probe(mol, center_chy, order_tuple, target_parity):
            records.append(LogRecord('openbabel:tetrahedral-not-placeable', (center_chy,),
                                     f'atom {center_chy}: could not set tetrahedral parity from OBMol; center '
                                     f'cleared', LOST))

    # OBMol reports cis/trans on both terminals, so the bonds are deduplicated here.
    seen_ct_bonds = set()
    for i in range(1, data.NumAtoms() + 1):
        atom = data.GetAtom(i)
        ob_id = atom.GetId()
        if not facade.HasCisTransStereo(ob_id):
            continue
        ct = facade.GetCisTransStereo(ob_id)
        cfg = ct.GetConfig()
        if not cfg.specified:
            continue

        begin_chy = ob_id_to_chy.get(cfg.begin)
        end_chy = ob_id_to_chy.get(cfg.end)
        if begin_chy is None or end_chy is None:
            continue

        bond_key = (min(begin_chy, end_chy), max(begin_chy, end_chy))
        if bond_key in seen_ct_bonds:
            continue
        seen_ct_bonds.add(bond_key)

        # the primary substituent is the first non-implicit ref on each side
        def _primary(refs_slice, _ir=implicit_ref):
            for r in refs_slice:
                if r != _ir:
                    return r
            return None

        begin_primary_ob = _primary([cfg.refs[0], cfg.refs[1]])
        end_primary_ob = _primary([cfg.refs[2], cfg.refs[3]])
        if begin_primary_ob is None or end_primary_ob is None:
            continue  # cannot determine geometry without at least one ref per side

        # IsCis is symmetric: True iff the two primaries are on the same side, which is parity 2.
        is_cis = ct.IsCis(begin_primary_ob, end_primary_ob)
        target_parity = 2 if is_cis else 1

        cb = mol.chiral_bonds()
        unit = cb.get(bond_key)
        if unit is None:
            continue  # V3 does not recognise this as a stereogenic CT bond

        anchor_chy = unit['anchor']
        try:
            mol.set_parity(anchor_chy, target_parity)
        except (KeyError, ValueError):
            records.append(LogRecord('openbabel:cis-trans-not-placeable', (begin_chy, end_chy),
                                     f'bond {begin_chy}-{end_chy}: could not set cis/trans parity; cleared',
                                     LOST))

    records.append(LogRecord('openbabel:coordinates-not-imported', (),
                             'coordinates are not imported; the source molecule\'s layout is dropped', LOST))
    deliver(mol, records, log)
    return mol
