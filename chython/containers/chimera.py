# -*- coding: utf-8 -*-
#
#  Copyright 2025, 2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from functools import cache
from typing import TYPE_CHECKING
from .._java import get_cdk


if TYPE_CHECKING:
    from chython import MoleculeContainer


@cache
def _cdk_classes():
    """Resolve and cache the CDK/JPype handles used by the direct builder."""
    from jpype import JArray, JClass

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
        'THStereo': JClass('org.openscience.cdk.interfaces.ITetrahedralChirality$Stereo'),
        'DBConf': JClass('org.openscience.cdk.interfaces.IDoubleBondStereochemistry$Conformation'),
        # aromatic (order 4) has no CDK bond order; UNSET + aromatic flag mirrors CDK's own kekule-free model
        'order_map': {1: order.SINGLE, 2: order.DOUBLE, 3: order.TRIPLE, 4: order.UNSET, 8: order.UNSET},
    }


class Chimera:
    __slots__ = ()

    def to_cdk(self: 'MoleculeContainer'):
        """
        Convert molecule to CDK IAtomContainer object.

        Atoms are built directly from the graph, so atom order is preserved:
        `mol.getAtom(i)` (0-based) corresponds to the i-th atom yielded by
        `self.atoms()`. Aromatic, tetrahedral and cis-trans stereo are set
        directly. Allene/cumulene stereo is not exported: CDK cannot round-trip
        it here, mirroring `to_rdkit`.

        The build is ~2x slower than the SMILES round-trip (every atom and bond
        crosses the JPype/JVM boundary), but preserves atom order without a
        separate `smiles_atoms_order` lookup.
        """
        c = _cdk_classes()
        Atom = c['Atom']
        Integer = c['Integer']
        order_map = c['order_map']

        mol = c['builder'].newAtomContainer()
        idx = {}  # chython atom number -> position in the atom array (0-based)
        atoms = []
        radicals = []
        for i, (n, a) in enumerate(self.atoms()):
            oa = Atom(a.atomic_number)
            if a.charge:
                oa.setFormalCharge(Integer(a.charge))
            if a.isotope:
                oa.setMassNumber(Integer(a.isotope))
            if a.is_radical:
                radicals.append(i)
            if a.implicit_hydrogens is not None:
                oa.setImplicitHydrogenCount(Integer(a.implicit_hydrogens))
            if a.hybridization == 4:
                oa.setIsAromatic(True)
            idx[n] = i
            atoms.append(oa)
        mol.setAtoms(c['IAtomArr'](atoms))
        for i in radicals:  # after setAtoms: a radical is one unpaired electron on that atom index
            mol.addSingleElectron(i)

        bond_of = {}  # (n, m) and (m, n) -> IBond, for stereo neighbour lookup
        for n, m, b in self.bonds():
            mol.addBond(idx[n], idx[m], order_map[b.order])
            bd = mol.getBond(mol.getBondCount() - 1)
            bond_of[n, m] = bond_of[m, n] = bd
            if b.order == 4:
                bd.setIsAromatic(True)

        TH = c['TetrahedralChirality']
        THStereo = c['THStereo']
        IAtomArr = c['IAtomArr']
        for n, a in self.atoms():
            if a.stereo is None or n not in self.stereogenic_tetrahedrons:
                continue  # allenes are not supported
            env = list(self._bonds[n])
            s = self._translate_tetrahedron_sign(n, env)
            ligands = [atoms[idx[x]] for x in env]
            if len(ligands) < 4:  # implicit hydrogen is encoded as the central atom itself
                ligands.append(atoms[idx[n]])
            winding = THStereo.ANTI_CLOCKWISE if s else THStereo.CLOCKWISE
            mol.addStereoElement(TH(atoms[idx[n]], IAtomArr(ligands), winding))

        DB = c['DoubleBondStereochemistry']
        DBConf = c['DBConf']
        IBondArr = c['IBondArr']
        for n, m, b in self.bonds():
            if b.stereo is None:
                continue
            nm = self._stereo_cis_trans_centers.get(n)
            if nm is None or n not in nm or m not in nm:
                continue
            n1, m1, *_ = self.stereogenic_cis_trans[nm]
            neighbours = IBondArr([bond_of[n, n1], bond_of[m, m1]])
            conf = DBConf.TOGETHER if b.stereo else DBConf.OPPOSITE  # True -> cis/Z
            mol.addStereoElement(DB(bond_of[n, m], neighbours, conf))

        return mol

    def to_openbabel(self: 'MoleculeContainer'):
        """
        Convert molecule to OpenBabel OBMol object.

        Atoms are built directly from the graph, so atom order is preserved:
        OBMol atom index `n` (1-based) corresponds to the n-th atom yielded by
        `self.atoms()`. Aromatic, tetrahedral and cis-trans stereo are set
        directly. Allene/cumulene stereo is not exported: OpenBabel cannot carry
        it (no extended-tetrahedral support in the bundled build), mirroring
        `to_rdkit`.
        """
        from openbabel import openbabel

        # OBMol atom ids default to idx-1; we reuse that as the stable stereo ref.
        implicit_ref = openbabel.OBStereo.ImplicitRef & 0xFFFFFFFFFFFFFFFF

        mol = openbabel.OBMol()
        # NB: do NOT wrap the build in BeginModify()/EndModify() - EndModify
        # re-perceives the structure and discards the stereo we set below.
        idx = {}  # chython atom number -> OBMol atom index (1-based, for bonds)
        ids = {}  # chython atom number -> OBAtom id (stereo refs are keyed by id)
        for n, a in self.atoms():
            oa = mol.NewAtom()
            oa.SetAtomicNum(a.atomic_number)
            if a.charge:
                oa.SetFormalCharge(a.charge)
            if a.isotope:
                oa.SetIsotope(a.isotope)
            if a.is_radical:
                oa.SetSpinMultiplicity(2)
            if a.implicit_hydrogens is not None:
                oa.SetImplicitHCount(a.implicit_hydrogens)
            if a.hybridization == 4:
                oa.SetAromatic(True)
            idx[n] = oa.GetIdx()
            ids[n] = oa.GetId()

        for n, m, b in self.bonds():
            # aromatic (order 4) has no numeric bond order: add as single + aromatic flag
            mol.AddBond(idx[n], idx[m], 1 if b.order == 4 else b.order)
            if b.order == 4:
                mol.GetBond(idx[n], idx[m]).SetAromatic(True)
        mol.SetAromaticPerceived(True)

        for n, a in self.atoms():
            if a.stereo is None:
                continue
            if n not in self.stereogenic_tetrahedrons:
                continue  # allenes are not supported
            env = list(self._bonds[n])
            s = self._translate_tetrahedron_sign(n, env)
            refs = [ids[x] for x in env]
            if len(refs) < 4:  # implicit hydrogen
                refs.append(implicit_ref)
            config = openbabel.OBTetrahedralConfig()
            config.center = ids[n]
            config.from_or_towards = refs[0]
            config.view = openbabel.OBStereo.ViewFrom
            config.refs = openbabel.OBStereo.MakeRefs(refs[1], refs[2], refs[3])
            config.winding = openbabel.OBStereo.AntiClockwise if s else openbabel.OBStereo.Clockwise
            config.specified = True
            ts = openbabel.OBTetrahedralStereo(mol)
            ts.SetConfig(config)
            mol.CloneData(ts)  # SetData is C++ only; CloneData is the python entry point

        for n, m, b in self.bonds():
            if b.stereo is None:
                continue
            nm = self._stereo_cis_trans_centers.get(n)
            if nm is None or n not in nm or m not in nm:
                continue
            n1, m1, *_ = self.stereogenic_cis_trans[nm]
            config = openbabel.OBCisTransConfig()
            config.begin = ids[n]
            config.end = ids[m]
            # ShapeU refs: [begin-neighbor, ?, ?, ?]; cis puts the end neighbor
            # at index 3 (same side), trans at index 2 (opposite side).
            if b.stereo:  # cis
                config.refs = openbabel.OBStereo.MakeRefs(ids[n1], implicit_ref, implicit_ref, ids[m1])
            else:  # trans
                config.refs = openbabel.OBStereo.MakeRefs(ids[n1], implicit_ref, ids[m1], implicit_ref)
            config.shape = openbabel.OBStereo.ShapeU
            config.specified = True
            cts = openbabel.OBCisTransStereo(mol)
            cts.SetConfig(config)
            mol.CloneData(cts)

        mol.SetChiralityPerceived(True)
        return mol

    def to_indigo(self: 'MoleculeContainer'):
        """
        Convert molecule to Indigo molecule object.

        Atoms are built directly from the graph, so atom order is preserved:
        `mol.iterateAtoms()` yields atoms in the same order as `self.atoms()`.
        Aromatic and tetrahedral stereo are set directly.

        Cis-trans and allene stereo are NOT exported: Indigo derives them from 2D
        coordinates (via `markStereobonds`), which a freshly parsed molecule does
        not have. This mirrors `to_rdkit`, which also drops allenes.
        """
        from indigo import Indigo

        mol = Indigo().createMolecule()
        mapping = {}
        for n, a in self.atoms():
            ia = mol.addAtom(a.atomic_symbol)
            if a.charge:
                ia.setCharge(a.charge)
            if a.isotope:
                ia.setIsotope(a.isotope)
            if a.is_radical:
                ia.setRadical(2)  # 2 == doublet in Indigo
            if a.implicit_hydrogens is not None:
                ia.setImplicitHCount(a.implicit_hydrogens)
            mapping[n] = ia.index()

        for n, m, b in self.bonds():
            # Indigo accepts aromatic (order 4) bonds directly
            mol.getAtom(mapping[n]).addBond(mol.getAtom(mapping[m]), b.order)

        for n, a in self.atoms():
            if a.stereo is None:
                continue
            if n not in self.stereogenic_tetrahedrons:
                continue  # allenes are not supported
            env = list(self._bonds[n])
            s = self._translate_tetrahedron_sign(n, env)
            pyramid = [mapping[x] for x in env]
            while len(pyramid) < 4:  # implicit hydrogen encoded as -1
                pyramid.append(-1)
            if s:  # swap two refs to flip handedness for chython's sign convention
                pyramid[1], pyramid[2] = pyramid[2], pyramid[1]
            mol.getAtom(mapping[n]).addStereocenter(Indigo.ABS, pyramid[0], pyramid[1], pyramid[2], pyramid[3])

        mol.aromatize()
        return mol


__all__ = ['Chimera']
