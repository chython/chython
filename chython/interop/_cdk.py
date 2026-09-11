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
CDK conversion, both directions.  Reached through `interop.cdk`.

The JVM comes from `_java.py:get_cdk`, which starts it once and caches the package -- do not start a
second one here.  `CDK_PATH` selects the jar.  Parity conventions and the aromatic/dative/H_UNKNOWN
bond mappings are documented at the functions that apply them.
"""
from functools import cache
from ..core import LogRecord, LOST, REPAIRED, SU_ALLENE, SU_CIS_TRANS, SU_TETRA
from ..exceptions import DirectionNotImplemented, UnconvertibleType
from ._records import deliver
from ._stereo import cis_trans_frame, kind_name, set_parity_by_probe


@cache
def _cdk_handles():
    """Resolve and cache CDK/JPype handles used by both directions."""
    from jpype import JArray, JClass
    from ._java import get_cdk

    cdk = get_cdk()
    order = JClass('org.openscience.cdk.interfaces.IBond$Order')
    return {
        'builder': cdk.silent.SilentChemObjectBuilder.getInstance(),
        'Atom': JClass('org.openscience.cdk.Atom'),
        'IAtomArr': JArray(JClass('org.openscience.cdk.interfaces.IAtom')),
        'IBondArr': JArray(JClass('org.openscience.cdk.interfaces.IBond')),
        'Integer': JClass('java.lang.Integer'),
        'TetrahedralChirality': JClass('org.openscience.cdk.stereo.TetrahedralChirality'),
        'DoubleBondStereochemistry': JClass('org.openscience.cdk.stereo.DoubleBondStereochemistry'),
        'ExtendedTetrahedral': JClass('org.openscience.cdk.stereo.ExtendedTetrahedral'),
        'THStereo': JClass('org.openscience.cdk.interfaces.ITetrahedralChirality$Stereo'),
        'DBConf': JClass('org.openscience.cdk.interfaces.IDoubleBondStereochemistry$Conformation'),
        # aromatic (order 4) has no CDK bond order: UNSET + the aromatic flag is CDK's kekule-free
        # model.  Dative (order 8) also maps to UNSET, but *without* the aromatic flag.
        'order_map': {1: order.SINGLE, 2: order.DOUBLE, 3: order.TRIPLE, 4: order.UNSET, 8: order.UNSET},
    }


def to_cdk(mol, /, *, log=None):
    """
    Export a chython container as a CDK `IAtomContainer`.

    Atoms are emitted in the molecule's own iteration order, so CDK's 0-based index i is the i-th
    atom; that is part of the interface, since `depict/layout/molecule.py` maps coordinates back by
    position.  CDK has no aromatic bond order, so order 4 is written as `Order.UNSET` with the
    aromatic flag on the bond and its atoms, never kekulized; order 8 (dative) is `Order.UNSET`
    without the flag and is reported as a loss.  H_UNKNOWN maps to Java null without loss.
    Tetrahedral, cis/trans and allene stereo are exported; atropisomers are not (no CDK form).

    :param log: optional list receiving one human-readable line per reportable loss.
    """
    from ..core import MoleculeContainer as CoreMoleculeContainer

    if isinstance(mol, CoreMoleculeContainer):
        return _to_cdk_v3(mol, log=log)

    raise UnconvertibleType(
        f'CDK export is not supported for {type(mol).__name__}; only MoleculeContainer'
    )


def from_cdk(data, /, *, log=None):
    """
    Import a CDK `IAtomContainer` as a chython V3 `MoleculeContainer`.

    Aromaticity is taken from what CDK says (UNSET + aromatic flag), not re-perceived.  An UNSET bond
    without the aromatic flag becomes a dative bond (order 8) and is reported as a loss.

    EVERY RECORD LANDS ON THE RETURNED MOLECULE'S `.log`, in stage `'interop'`, with nothing passed
    in; `log`, when given, receives a copy of the same records.

    :param log: optional list receiving a copy of the records put on the molecule.
    """
    try:
        from jpype import JClass
        IAtomContainer = JClass('org.openscience.cdk.interfaces.IAtomContainer')
    except (ImportError, TypeError, RuntimeError):
        # THREE ERRORS FOR ONE MISSING TOOLKIT, and only the first is an `ImportError`: jpype absent
        # raises that, a name the classpath cannot resolve raises `TypeError`, and no JVM at all raises
        # `JVMNotRunning`, a `RuntimeError`.  An environment with jpype installed and no CDK jar is the
        # common one -- it is what CI has -- and there the refusal has to be this exception and not
        # whichever one jpype picked.
        raise DirectionNotImplemented(
            'interop.cdk import requires jpype and CDK; install them and set CDK_PATH'
        )
    if not isinstance(data, IAtomContainer):
        raise DirectionNotImplemented(
            f'interop.cdk import direction does not accept {type(data).__name__}; '
            f'pass a CDK IAtomContainer'
        )
    return _from_cdk(data, log=log)


def _to_cdk_v3(mol, *, log=None):
    """Export a V3 MoleculeContainer to a CDK IAtomContainer."""
    c = _cdk_handles()
    Atom = c['Atom']
    Integer = c['Integer']
    order_map = c['order_map']
    IAtomArr = c['IAtomArr']
    IBondArr = c['IBondArr']

    cdk_mol = c['builder'].newAtomContainer()

    idx = {}       # atom number -> 0-based CDK index
    cdk_atoms = []
    radicals = []

    for i, atom in enumerate(mol.atoms()):
        n = atom.n
        oa = Atom(Integer(atom.element))
        if atom.charge:
            oa.setFormalCharge(Integer(atom.charge))
        if atom.isotope:
            oa.setMassNumber(Integer(atom.isotope))
        if atom.is_radical:
            radicals.append(i)
        ih = atom.implicit_h
        if ih is not None:
            oa.setImplicitHydrogenCount(Integer(ih))
        else:
            oa.setImplicitHydrogenCount(None)  # H_UNKNOWN -> CDK null (preserves the unknown state)
        idx[n] = i
        cdk_atoms.append(oa)

    cdk_mol.setAtoms(IAtomArr(cdk_atoms))
    for i in radicals:
        cdk_mol.addSingleElectron(i)

    bond_of = {}  # (n, m) and (m, n) -> CDK IBond, for stereo lookup
    for bond in mol.bonds():
        n, m, order = bond.n, bond.m, bond.order
        if order == 8 and log is not None:
            log.append(LogRecord('cdk:dative-bond-written-as-unset', (n, m),
                                 f'atom {n}-atom {m}: dative bond (order 8) written as UNSET; geometry is lost',
                                 LOST))
        cdk_mol.addBond(idx[n], idx[m], order_map[order])
        bd = cdk_mol.getBond(cdk_mol.getBondCount() - 1)
        bond_of[n, m] = bond_of[m, n] = bd
        if order == 4:
            bd.setIsAromatic(True)
            cdk_atoms[idx[n]].setIsAromatic(True)
            cdk_atoms[idx[m]].setIsAromatic(True)

    TH = c['TetrahedralChirality']
    THStereo = c['THStereo']
    DB = c['DoubleBondStereochemistry']
    DBConf = c['DBConf']
    ET = c['ExtendedTetrahedral']

    # adjacency for allene terminal discovery and cis/trans other-terminal lookup
    adj = {}
    for bond in mol.bonds():
        n, m = bond.n, bond.m
        adj.setdefault(n, []).append(m)
        adj.setdefault(m, []).append(n)

    for unit in mol.stereo_units():
        parity = unit['parity']
        if parity == 0:
            continue
        kind = unit['kind']
        anchor = unit['anchor']
        refs = unit['refs']

        if kind == SU_TETRA:
            _export_tetrahedral(
                anchor, refs, parity, idx, cdk_atoms, cdk_mol, TH, THStereo, IAtomArr
            )
        elif kind == SU_CIS_TRANS:
            _export_cis_trans(
                anchor, refs, mol, idx, cdk_mol, bond_of, DB, DBConf, IBondArr, log
            )
        elif kind == SU_ALLENE:
            _export_allene(
                anchor, refs, parity, idx, cdk_atoms, adj, cdk_mol, ET, THStereo, IAtomArr
            )
        elif log is not None:
            # atropisomers, and whatever kind the core grows next: reported, never silently dropped
            log.append(LogRecord('cdk:stereo-kind-not-representable', (anchor,),
                                 f'{kind_name(kind)} stereo at atom {anchor} has no CDK representation; dropped',
                                 LOST))

    return cdk_mol


def _export_tetrahedral(anchor, refs, parity, idx, cdk_atoms, cdk_mol, TH, THStereo, IAtomArr):
    """Add a CDK TetrahedralChirality element for one V3 tetrahedral unit."""
    focus_cdk = cdk_atoms[idx[anchor]]
    ligands = []
    for r in refs:
        if r is None:
            ligands.append(focus_cdk)   # implicit H / lone pair: CDK convention = focus itself
        else:
            ligands.append(cdk_atoms[idx[r]])
    # parity 2 (odd) = anticlockwise = @; parity 1 (even) = clockwise = @@
    winding = THStereo.ANTI_CLOCKWISE if parity == 2 else THStereo.CLOCKWISE
    cdk_mol.addStereoElement(TH(focus_cdk, IAtomArr(ligands), winding))


def _export_cis_trans(anchor, refs, mol, idx, cdk_mol, bond_of, DB, DBConf, IBondArr, log=None):
    """Add a CDK DoubleBondStereochemistry element for one V3 cis/trans unit.

    Framing comes from `_stereo.cis_trans_frame`, shared with `_rdkit` and `_cdpkit`.  The parity is
    read in the frame actually written, not in the unit's own refs order: the two coincide only when
    the near reference happens to be `refs[0]`, which the core does not promise.
    """
    frame = cis_trans_frame(mol, anchor, refs)
    if frame is None:
        if log is not None:
            log.append(LogRecord('cdk:cis-trans-frame-not-found', (anchor,),
                                 f'cis/trans stereo at atom {anchor} cannot be framed on the graph; dropped',
                                 LOST))
        return

    try:
        nb_bond_anchor = bond_of[anchor, frame.near]
        nb_bond_other = bond_of[frame.partner, frame.far]
        double_bond = bond_of[anchor, frame.partner]
    except KeyError:
        if log is not None:
            log.append(LogRecord('cdk:cis-trans-framing-bond-absent', (anchor,),
                                 f'cis/trans stereo at atom {anchor}: a framing bond is absent from the CDK '
                                 f'molecule; dropped', LOST))
        return

    try:
        p = mol.translate_stereo(anchor, frame.order)
    except (KeyError, ValueError) as e:
        if log is not None:
            log.append(LogRecord('cdk:cis-trans-config-not-readable', (anchor,),
                                 f'cis/trans stereo at atom {anchor}: chython will not read the configuration '
                                 f"in CDK's frame ({e}); dropped", LOST))
        return
    if p == 0:
        return

    # parity 1 (even) = OPPOSITE (trans); parity 2 (odd) = TOGETHER (cis)
    conf = DBConf.OPPOSITE if p == 1 else DBConf.TOGETHER
    cdk_mol.addStereoElement(DB(double_bond, IBondArr([nb_bond_anchor, nb_bond_other]), conf))


def _export_allene(anchor, refs, parity, idx, cdk_atoms, adj, cdk_mol, ET, THStereo, IAtomArr):
    """Add a CDK ExtendedTetrahedral element for one V3 allene unit.

    V3 refs layout: (near_sub0, near_sub1, far_sub0, far_sub1).
    anchor = central sp-carbon.
    CDK peripherals: [near_sub0 or t1, t1, t2, far_sub0 or t2].
    Winding: parity 2 → ANTI_CLOCKWISE, parity 1 → CLOCKWISE.
    """
    refs_0, refs_1, refs_2, refs_3 = refs
    focus_cdk = cdk_atoms[idx[anchor]]

    terminals = adj.get(anchor, [])
    if len(terminals) < 2:
        return

    # Identify near terminal (t1: has refs_0 or refs_1 as neighbor) vs far (t2: has refs_2/refs_3)
    t1, t2 = None, None
    for term in terminals:
        term_nbs = set(adj.get(term, []))
        term_nbs.discard(anchor)
        near_match = (refs_0 is not None and refs_0 in term_nbs) or \
                     (refs_1 is not None and refs_1 in term_nbs)
        far_match  = (refs_2 is not None and refs_2 in term_nbs) or \
                     (refs_3 is not None and refs_3 in term_nbs)
        if near_match and t1 is None:
            t1 = term
        elif far_match and t2 is None:
            t2 = term

    if t1 is None or t2 is None:
        t1, t2 = terminals[0], terminals[1]

    t1_cdk = cdk_atoms[idx[t1]]
    t2_cdk = cdk_atoms[idx[t2]]

    p0 = cdk_atoms[idx[refs_0]] if refs_0 is not None else t1_cdk
    p3 = cdk_atoms[idx[refs_2]] if refs_2 is not None else t2_cdk

    winding = THStereo.ANTI_CLOCKWISE if parity == 2 else THStereo.CLOCKWISE
    cdk_mol.addStereoElement(ET(focus_cdk, IAtomArr([p0, t1_cdk, t2_cdk, p3]), winding))


def _from_cdk(data, *, log=None):
    """Import a CDK IAtomContainer as a V3 MoleculeContainer."""
    from ..core import MoleculeContainer
    from jpype import JClass

    Integer = JClass('java.lang.Integer')
    ITetrahedralChirality = JClass('org.openscience.cdk.interfaces.ITetrahedralChirality')
    IDoubleBondStereochemistry = JClass(
        'org.openscience.cdk.interfaces.IDoubleBondStereochemistry'
    )
    try:
        ExtendedTetrahedral = JClass('org.openscience.cdk.stereo.ExtendedTetrahedral')
        has_et = True
    except Exception:
        has_et = False

    mol = MoleculeContainer()
    records = []
    n_atoms = data.getAtomCount()

    # atom index (0-based CDK) -> V3 atom number
    atom_map = {}

    for i in range(n_atoms):
        a = data.getAtom(i)
        an = int(a.getAtomicNumber() or 0)
        if an == 0:
            an = 6  # fallback: treat unknown atomic number as carbon

        charge = int(a.getFormalCharge() or 0)

        iso_obj = a.getMassNumber()
        isotope = int(iso_obj) if iso_obj is not None else 0

        radical = data.getConnectedSingleElectronsCount(a) > 0

        # Java null implicit H count -> H_UNKNOWN; only CDK preserves the unknown state
        ih_obj = a.getImplicitHydrogenCount()
        implicit_h = int(ih_obj) if ih_obj is not None else None

        sid = mol.add_atom(an, charge=charge, isotope=isotope,
                           radical=radical, implicit_h=implicit_h)
        atom_map[i] = sid

    n_bonds = data.getBondCount()
    for i in range(n_bonds):
        b = data.getBond(i)
        atoms_in_bond = [b.getAtom(j) for j in range(b.getAtomCount())]
        if len(atoms_in_bond) != 2:
            continue
        ia0 = int(data.indexOf(atoms_in_bond[0]))
        ia1 = int(data.indexOf(atoms_in_bond[1]))
        if ia0 < 0 or ia1 < 0:
            continue
        sid0 = atom_map[ia0]
        sid1 = atom_map[ia1]

        order_obj = b.getOrder()
        order_name = str(order_obj.name()) if order_obj is not None else 'UNSET'
        is_aromatic = bool(b.isAromatic())

        if order_name == 'SINGLE':
            order = 1
        elif order_name == 'DOUBLE':
            order = 2
        elif order_name == 'TRIPLE':
            order = 3
        elif order_name == 'UNSET':
            if is_aromatic:
                order = 4
            else:
                order = 8
                records.append(LogRecord('cdk:unset-bond-imported-as-dative', (sid0, sid1),
                                         f'bond {sid0}-{sid1}: UNSET non-aromatic bond imported as dative (order 8)',
                                         REPAIRED))
        else:
            order = 1  # unknown bond type: treat as single
        mol.add_bond(sid0, sid1, order)

    for se in data.stereoElements():
        # a malformed element is skipped and reported, never raised: it must not cost the caller the
        # rest of the molecule
        try:
            if isinstance(se, ITetrahedralChirality):
                _import_tetrahedral(se, data, atom_map, mol, records)
            elif isinstance(se, IDoubleBondStereochemistry):
                _import_cis_trans(se, data, atom_map, mol, records)
            elif has_et and isinstance(se, ExtendedTetrahedral):
                _import_allene(se, data, atom_map, mol, records)
            else:
                records.append(LogRecord('cdk:stereo-element-no-equivalent', (),
                                         f'CDK stereo element {type(se).__name__} has no chython equivalent; '
                                         f'dropped', LOST))
        except Exception as e:
            records.append(LogRecord('cdk:stereo-element-read-error', (),
                                     f'a CDK stereo element could not be read ({type(e).__name__}: {e}); dropped',
                                     LOST))

    records.append(LogRecord('cdk:coordinates-not-imported', (),
                             'coordinates are not imported; the source molecule\'s layout is dropped', LOST))
    deliver(mol, records, log)
    return mol


def _import_tetrahedral(se, data, atom_map, mol, records):
    """Set V3 tetrahedral parity from a CDK TetrahedralChirality.  True when one was placed."""
    focus = se.getChiralAtom()
    focus_idx = int(data.indexOf(focus))
    anchor = atom_map[focus_idx]

    ligands = se.getLigands()
    order = []
    for lig in ligands:
        li = int(data.indexOf(lig))
        if li == focus_idx:
            order.append(None)   # focus itself = implicit H placeholder
        else:
            order.append(atom_map[li])
    order = tuple(order)

    stereo = se.getStereo()
    # ANTI_CLOCKWISE → parity 2; CLOCKWISE → parity 1
    cdk_parity = 2 if str(stereo.name()) == 'ANTI_CLOCKWISE' else 1

    if _set_parity_via_probe(mol, anchor, order, cdk_parity):
        return True
    records.append(LogRecord('cdk:tetrahedral-not-placeable', (anchor,),
                             f'atom {anchor}: CDK reports a tetrahedral configuration that chython cannot place '
                             f'on a stereo unit here; dropped', LOST))
    return False


def _import_cis_trans(se, data, atom_map, mol, records):
    """Set V3 cis/trans parity from a CDK DoubleBondStereochemistry.  True when one was placed."""
    db = se.getStereoBond()
    db_atoms = [db.getAtom(j) for j in range(db.getAtomCount())]
    if len(db_atoms) != 2:
        return
    ia0 = int(data.indexOf(db_atoms[0]))
    ia1 = int(data.indexOf(db_atoms[1]))
    anchor_cdk = ia0
    other_cdk = ia1

    anchor = atom_map[anchor_cdk]
    other_terminal = atom_map[other_cdk]

    nb_bonds = se.getBonds()
    if len(nb_bonds) < 2:
        return

    def _other_atom_idx(bond, given_cdk_idx):
        for j in range(bond.getAtomCount()):
            ai = int(data.indexOf(bond.getAtom(j)))
            if ai != given_cdk_idx:
                return ai
        return -1

    sub_anchor_cdk = _other_atom_idx(nb_bonds[0], anchor_cdk)
    sub_other_cdk  = _other_atom_idx(nb_bonds[1], other_cdk)
    if sub_anchor_cdk < 0 or sub_other_cdk < 0:
        return

    sub_anchor = atom_map[sub_anchor_cdk]
    sub_other  = atom_map[sub_other_cdk]

    conf = se.getStereo()
    # TOGETHER = cis → parity 2; OPPOSITE = trans → parity 1
    cdk_parity = 2 if str(conf.name()) == 'TOGETHER' else 1

    # either terminal may anchor the unit, so both framings are tried
    order_a = (sub_anchor, None, sub_other, None)
    order_b = (sub_other, None, sub_anchor, None)
    for unit_anchor, order in [(anchor, order_a), (other_terminal, order_b)]:
        if _set_parity_via_probe(mol, unit_anchor, order, cdk_parity):
            return True
    records.append(LogRecord('cdk:cis-trans-not-placeable', (anchor, other_terminal),
                             f'bond {anchor}-{other_terminal}: CDK reports a double-bond configuration that '
                             f'chython cannot place on a stereo unit at either terminal; dropped', LOST))
    return False


def _import_allene(se, data, atom_map, mol, records):
    """Set V3 allene parity from a CDK ExtendedTetrahedral.

    CDK uses `.peripherals()` (not getPeripherals) and `.winding()` (not getStereo).
    Peripherals layout: [p0, t1, t2, p3] where t1/t2 are terminal atoms and p0/p3 are
    their substituents (terminal atom serves as implicit-H placeholder if p == terminal).
    """
    focus = se.getFocus()
    focus_idx = int(data.indexOf(focus))
    anchor = atom_map[focus_idx]

    peripherals = se.peripherals()
    if len(peripherals) < 4:
        return
    p_idxs = [int(data.indexOf(p)) for p in peripherals]
    t1_cdk = p_idxs[1]
    t2_cdk = p_idxs[2]
    p0_cdk = p_idxs[0]
    p3_cdk = p_idxs[3]

    # a peripheral equal to its terminal atom is the implicit-H placeholder -> None in V3
    p0 = None if p0_cdk == t1_cdk else atom_map[p0_cdk]
    p3 = None if p3_cdk == t2_cdk else atom_map[p3_cdk]

    winding = se.winding()
    cdk_parity = 2 if str(winding.name()) == 'ANTI_CLOCKWISE' else 1

    order = (p0, None, p3, None)
    if _set_parity_via_probe(mol, anchor, order, cdk_parity):
        return True
    records.append(LogRecord('cdk:allene-not-placeable', (anchor,),
                             f'atom {anchor}: CDK reports an allene configuration that chython cannot place on a '
                             f'stereo unit here; dropped', LOST))
    return False


def _set_parity_via_probe(mol, anchor, order, cdk_parity):
    """Write the stored parity that reads back as CDK's in `order`.  See `_stereo.set_parity_by_probe`.

    A frame the core rejects clears the centre and returns `False` rather than leaving half a
    configuration behind, which is indistinguishable from a real one.
    """
    return set_parity_by_probe(mol, anchor, order, cdk_parity)
