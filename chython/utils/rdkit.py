# -*- coding: utf-8 -*-
#
#  Copyright 2019-2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from rdkit.Chem import AssignStereochemistry, Atom, BondStereo, BondType, ChiralType, Conformer, RWMol, SanitizeMol
from ..containers import MoleculeContainer
from ..exceptions import IsChiral, NotChiral, ValenceError
from ..periodictable import Element


def from_rdkit_molecule(data):
    """
    RDKit molecule object to MoleculeContainer converter
    """
    mol = MoleculeContainer()
    parsed_mapping = mol._parsed_mapping
    mol_conformers = mol._conformers
    bonds = mol._bonds

    atoms, mapping = [], []
    tetrahedron_stereo = []
    for a in data.GetAtoms():
        e = Element.from_symbol(a.GetSymbol())
        isotope = a.GetIsotope()
        if isotope:
            e = e(isotope)
        else:
            e = e()
        atom = {'atom': e, 'charge': a.GetFormalCharge()}

        radical = a.GetNumRadicalElectrons()
        if radical:
            atom['is_radical'] = True

        atoms.append(atom)
        mapping.append(a.GetAtomMapNum())
        tetrahedron_stereo.append(a.GetChiralTag())

    conformers = []
    c = data.GetConformers()
    if c:
        for atom, (x, y, _) in zip(atoms, c[0].GetPositions()):
            atom['xy'] = (x, y)
        for c in c:
            if c.Is3D():
                conformers.append(c.GetPositions())

    new_map = []
    for a, n in zip(atoms, mapping):
        a = mol.add_atom(**a)
        new_map.append(a)
        parsed_mapping[a] = n

    stereo = []
    for b in data.GetBonds():
        n, m = new_map[b.GetBeginAtomIdx()], new_map[b.GetEndAtomIdx()]
        mol.add_bond(n, m, _rdkit_bond_map[b.GetBondType()])
        s = b.GetStereo()
        if s == _cis:
            nn, nm = b.GetStereoAtoms()
            stereo.append((mol.add_cis_trans_stereo, n, m, new_map[nn], new_map[nm], True))
        elif s == _trans:
            nn, nm = b.GetStereoAtoms()
            stereo.append((mol.add_cis_trans_stereo, n, m, new_map[nn], new_map[nm], False))

    for n, s in zip(new_map, tetrahedron_stereo):
        if s == _chiral_cw:
            env = bonds[n]
            env = [x for x in new_map if x in env]
            stereo.append((mol.add_atom_stereo, n, env, False))
        elif s == _chiral_ccw:
            env = bonds[n]
            env = [x for x in new_map if x in env]
            stereo.append((mol.add_atom_stereo, n, env, True))

    while stereo:
        fail_stereo = []
        old_stereo = len(stereo)
        for f, *args in stereo:
            try:
                f(*args, clean_cache=False)
            except NotChiral:
                fail_stereo.append((f, *args))
            except IsChiral:
                pass
            except ValenceError:
                mol.flush_cache()
                break
        else:
            stereo = fail_stereo
            if len(stereo) == old_stereo:
                break
            mol.flush_stereo_cache()
            continue
        break

    for c in conformers:
        mol_conformers.append({k: tuple(v) for k, v in zip(new_map, c)})
    return mol


def to_rdkit_molecule(data: MoleculeContainer):
    """
    MoleculeContainer to RDKit molecule object converter.

    Note: implicit hydrogens data omitted.
    """
    mol = RWMol()
    mapping = {}
    bonds = data._bonds

    for n, a in data.atoms():
        ra = Atom(a.atomic_number)
        ra.SetAtomMapNum(n)
        if a.charge:
            ra.SetFormalCharge(a.charge)
        if a.isotope:
            ra.SetIsotope(a.isotope)
        if a.is_radical:
            ra.SetNumRadicalElectrons(1)
        mapping[n] = mol.AddAtom(ra)

    for n, m, b in data.bonds():
        mol.AddBond(mapping[n], mapping[m], _bond_map[b.order])

    for n in data._atoms_stereo:
        ra = mol.GetAtomWithIdx(mapping[n])
        env = bonds[n]
        s = data._translate_tetrahedron_sign(n, [x for x in mapping if x in env])
        ra.SetChiralTag(_chiral_ccw if s else _chiral_cw)

    for nm, s in data._cis_trans_stereo.items():
        n, m = nm
        if m in bonds[n]:  # cumulenes unsupported
            nn, nm, *_ = data._stereo_cis_trans[nm]
            b = mol.GetBondBetweenAtoms(mapping[n], mapping[m])
            b.SetStereoAtoms(mapping[nn], mapping[nm])
            b.SetStereo(_cis if s else _trans)

    conf = Conformer()
    for n, a in data.atoms():
        conf.SetAtomPosition(mapping[n], (a.x, a.y, 0))
    conf.Set3D(False)
    mol.AddConformer(conf, assignId=True)

    for c in data._conformers:
        conf = Conformer()
        for n, xyz in c.items():
            conf.SetAtomPosition(mapping[n], xyz)
        mol.AddConformer(conf, assignId=True)

    SanitizeMol(mol)
    AssignStereochemistry(mol, flagPossibleStereoCenters=True, force=True)
    return mol


_rdkit_bond_map = {BondType.SINGLE: 1, BondType.DOUBLE: 2, BondType.TRIPLE: 3, BondType.AROMATIC: 4}
_bond_map = {1: BondType.SINGLE, 2: BondType.DOUBLE, 3: BondType.TRIPLE, 4: BondType.AROMATIC}

_chiral_cw = ChiralType.CHI_TETRAHEDRAL_CW
_chiral_ccw = ChiralType.CHI_TETRAHEDRAL_CCW
_trans = BondStereo.STEREOE
_cis = BondStereo.STEREOZ

__all__ = ['from_rdkit_molecule', 'to_rdkit_molecule']
