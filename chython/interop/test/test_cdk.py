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
CDK converter tests -- both directions.

Skipped when jpype or the CDK jar is absent; the `cdk` fixture (conftest.py) handles that.
"""
import pytest
from chython.interop import cdk as interop_cdk
from chython.interop._cdk import to_cdk, from_cdk
from chython.core import MoleculeContainer as V3Molecule
from chython.exceptions import DirectionNotImplemented, UnconvertibleType

from .conftest import requires_jpype


@requires_jpype
def test_atom_order_is_preserved_in_to_cdk(cdk):
    """CDK atom index i must correspond to the i-th atom in chython's iteration order."""
    mol = V3Molecule()
    c = mol.add_atom(6)
    n = mol.add_atom(7)
    o = mol.add_atom(8)
    mol.add_bond(c, n, 1)
    mol.add_bond(n, o, 2)
    cmol = to_cdk(mol)
    atoms = [mol.atoms().__next__().n]
    numbers = [a.n for a in mol.atoms()]
    for i, sid in enumerate(numbers):
        expected_an = next(a.element for a in mol.atoms() if a.n == sid)
        assert int(cmol.getAtom(i).getAtomicNumber()) == expected_an


@requires_jpype
def test_constitution_elements(cdk):
    """Elements survive the round-trip (CHN molecule: alanine-like)."""
    mol = V3Molecule()
    c1 = mol.add_atom(6, implicit_h=3)   # CH3
    c2 = mol.add_atom(6, implicit_h=1)   # CH
    n  = mol.add_atom(7, implicit_h=2)   # NH2
    o1 = mol.add_atom(8, implicit_h=0)   # C=O
    o2 = mol.add_atom(8, implicit_h=1)   # OH
    mol.add_bond(c1, c2, 1)
    mol.add_bond(c2, n,  1)
    mol.add_bond(c2, o1, 2)
    mol.add_bond(o1, o2, 1)

    cmol = to_cdk(mol)
    back = from_cdk(cmol)

    orig_an = sorted(a.element for a in mol.atoms())
    back_an = sorted(a.element for a in back.atoms())
    assert orig_an == back_an


@requires_jpype
def test_constitution_charges(cdk):
    """Formal charges survive the round-trip."""
    mol = V3Molecule()
    n_plus = mol.add_atom(7, charge=1, implicit_h=4)
    cl_minus = mol.add_atom(17, charge=-1, implicit_h=0)
    cmol = to_cdk(mol)
    back = from_cdk(cmol)
    charges = {a.element: a.charge for a in back.atoms()}
    assert charges[7] == 1
    assert charges[17] == -1


@requires_jpype
def test_constitution_isotopes(cdk):
    """Isotope labels survive the round-trip."""
    mol = V3Molecule()
    c13 = mol.add_atom(6, isotope=13, implicit_h=4)
    h2  = mol.add_atom(1, isotope=2)
    mol.add_bond(c13, h2, 1)
    cmol = to_cdk(mol)
    back = from_cdk(cmol)
    isos = {a.element: a.isotope for a in back.atoms()}
    assert isos[6] == 13
    assert isos[1] == 2


@requires_jpype
def test_constitution_radical(cdk):
    """Radical flag survives the round-trip."""
    mol = V3Molecule()
    c_rad = mol.add_atom(6, radical=True, implicit_h=3)
    cmol = to_cdk(mol)
    back = from_cdk(cmol)
    atoms = list(back.atoms())
    assert atoms[0].is_radical


@requires_jpype
def test_constitution_bond_orders(cdk):
    """Bond orders 1, 2, 3 survive the round-trip."""
    mol = V3Molecule()
    c1 = mol.add_atom(6, implicit_h=3)
    c2 = mol.add_atom(6, implicit_h=1)
    c3 = mol.add_atom(6, implicit_h=0)
    c4 = mol.add_atom(6, implicit_h=1)
    mol.add_bond(c1, c2, 1)
    mol.add_bond(c2, c3, 2)
    mol.add_bond(c3, c4, 3)
    cmol = to_cdk(mol)
    back = from_cdk(cmol)
    orders = sorted(b.order for b in back.bonds())
    assert orders == [1, 2, 3]


@requires_jpype
def test_constitution_aromatic_bond(cdk):
    """Aromatic bonds (order 4) survive the round-trip via UNSET+aromatic flag."""
    mol = V3Molecule()
    c1 = mol.add_atom(6, implicit_h=1)
    c2 = mol.add_atom(6, implicit_h=1)
    c3 = mol.add_atom(6, implicit_h=1)
    c4 = mol.add_atom(6, implicit_h=1)
    c5 = mol.add_atom(6, implicit_h=1)
    c6 = mol.add_atom(6, implicit_h=1)
    mol.add_bond(c1, c2, 4)
    mol.add_bond(c2, c3, 4)
    mol.add_bond(c3, c4, 4)
    mol.add_bond(c4, c5, 4)
    mol.add_bond(c5, c6, 4)
    mol.add_bond(c6, c1, 4)
    cmol = to_cdk(mol)
    back = from_cdk(cmol)
    orders = sorted(b.order for b in back.bonds())
    assert orders == [4, 4, 4, 4, 4, 4]


@requires_jpype
def test_h_unknown_is_preserved_export(cdk):
    """H_UNKNOWN (implicit_h=None) maps to CDK null without loss."""
    mol = V3Molecule()
    c = mol.add_atom(6)   # implicit_h=None -> H_UNKNOWN
    cmol = to_cdk(mol)
    assert cmol.getAtom(0).getImplicitHydrogenCount() is None


@requires_jpype
def test_h_unknown_is_preserved_import(cdk):
    """CDK null implicit H maps back to H_UNKNOWN (implicit_h_of returns None)."""
    mol = V3Molecule()
    c = mol.add_atom(6)   # H_UNKNOWN
    cmol = to_cdk(mol)
    back = from_cdk(cmol)
    atom = next(iter(back.atoms()))
    assert atom.implicit_h is None


@requires_jpype
def test_h_zero_is_preserved(cdk):
    """Explicit H count of 0 is not silently changed to H_UNKNOWN."""
    mol = V3Molecule()
    c = mol.add_atom(6, implicit_h=0)
    cmol = to_cdk(mol)
    ih = cmol.getAtom(0).getImplicitHydrogenCount()
    assert ih is not None
    assert int(ih) == 0
    back = from_cdk(cmol)
    atom = next(iter(back.atoms()))
    assert atom.implicit_h == 0


@requires_jpype
def test_aromatic_bond_carries_aromatic_flag(cdk):
    """An order-4 bond is written as CDK UNSET + aromatic flag; the flag is present."""
    mol = V3Molecule()
    c1 = mol.add_atom(6, implicit_h=1)
    c2 = mol.add_atom(6, implicit_h=1)
    mol.add_bond(c1, c2, 4)
    cmol = to_cdk(mol)
    bond = cmol.getBond(0)
    assert bond.isAromatic()
    assert str(bond.getOrder().name()) == 'UNSET'


@requires_jpype
def test_single_bond_has_no_aromatic_flag(cdk):
    """A single bond (order 1) does NOT carry the aromatic flag."""
    mol = V3Molecule()
    c1 = mol.add_atom(6, implicit_h=3)
    c2 = mol.add_atom(6, implicit_h=3)
    mol.add_bond(c1, c2, 1)
    cmol = to_cdk(mol)
    bond = cmol.getBond(0)
    assert not bond.isAromatic()
    assert str(bond.getOrder().name()) == 'SINGLE'


@requires_jpype
def test_dative_bond_reports_loss(cdk):
    """Dative bonds (order 8) are written as UNSET and logged as a loss."""
    mol = V3Molecule()
    n = mol.add_atom(7, implicit_h=3)
    b = mol.add_atom(5, implicit_h=3)
    mol.add_bond(n, b, 8)
    losses = []
    cmol = to_cdk(mol, log=losses)
    assert len(losses) == 1
    assert '8' in losses[0] or 'dative' in losses[0]
    bond = cmol.getBond(0)
    assert not bond.isAromatic()
    assert str(bond.getOrder().name()) == 'UNSET'


@requires_jpype
def test_tetrahedral_stereo_survives_roundtrip(cdk):
    """An alanine tetrahedral centre survives V3 -> CDK -> V3.

    Only that a parity is set, not its absolute sense: that depends on slot ordering.
    """
    v3 = V3Molecule()
    n_atom = v3.add_atom(7, implicit_h=2)
    c_atom = v3.add_atom(6, implicit_h=1)
    c_me   = v3.add_atom(6, implicit_h=3)
    c_coo  = v3.add_atom(6, implicit_h=0)
    o_dbl  = v3.add_atom(8, implicit_h=0)
    o_oh   = v3.add_atom(8, implicit_h=1)
    v3.add_bond(n_atom, c_atom, 1)
    v3.add_bond(c_atom, c_me,   1)
    v3.add_bond(c_atom, c_coo,  1)
    v3.add_bond(c_coo,  o_dbl,  2)
    v3.add_bond(c_coo,  o_oh,   1)

    v3.set_parity(c_atom, 2)

    cmol = to_cdk(v3)
    back = from_cdk(cmol)

    units = [u for u in back.stereo_units() if u['kind'] == 0]
    assert len(units) >= 1
    assert any(u['parity'] != 0 for u in units)


@requires_jpype
def test_tetrahedral_implicit_h_roundtrip(cdk):
    """Tetrahedral centre with implicit H (3 heavy neighbours + 1 H) survives."""
    v3 = V3Molecule()
    c  = v3.add_atom(6, implicit_h=1)  # the chiral centre
    f  = v3.add_atom(9,  implicit_h=0)
    cl = v3.add_atom(17, implicit_h=0)
    br = v3.add_atom(35, implicit_h=0)
    v3.add_bond(c, f,  1)
    v3.add_bond(c, cl, 1)
    v3.add_bond(c, br, 1)
    v3.set_parity(c, 1)

    cmol = to_cdk(v3)
    back = from_cdk(cmol)

    units = [u for u in back.stereo_units() if u['kind'] == 0]
    assert len(units) == 1
    assert units[0]['parity'] != 0


@requires_jpype
def test_cis_trans_stereo_survives_roundtrip(cdk):
    """E-but-2-ene cis/trans stereo survives V3 -> CDK -> V3."""
    v3 = V3Molecule()
    c1 = v3.add_atom(6, implicit_h=3)
    c2 = v3.add_atom(6, implicit_h=1)
    c3 = v3.add_atom(6, implicit_h=1)
    c4 = v3.add_atom(6, implicit_h=3)
    v3.add_bond(c1, c2, 1)
    v3.add_bond(c2, c3, 2)
    v3.add_bond(c3, c4, 1)
    # parity 1 = even = OPPOSITE (trans) for the lower-indexed terminal.
    v3.set_parity(c2, 1)

    cmol = to_cdk(v3)
    back = from_cdk(cmol)

    units = [u for u in back.stereo_units() if u['kind'] == 1]
    assert len(units) == 1
    assert units[0]['parity'] != 0


@requires_jpype
def test_allene_stereo_survives_roundtrip(cdk):
    """Chiral allene stereo survives V3 -> CDK -> V3 via ExtendedTetrahedral."""
    # (R)-1,3-dimethylallene: MeHC=C=CHMe
    v3 = V3Molecule()
    c_near = v3.add_atom(6, implicit_h=1)   # terminal near (CHMe)
    c_cen  = v3.add_atom(6, implicit_h=0)   # allene centre
    c_far  = v3.add_atom(6, implicit_h=1)   # terminal far (CHMe)
    me_near = v3.add_atom(6, implicit_h=3)
    me_far  = v3.add_atom(6, implicit_h=3)
    v3.add_bond(me_near, c_near, 1)
    v3.add_bond(c_near,  c_cen,  2)
    v3.add_bond(c_cen,   c_far,  2)
    v3.add_bond(c_far,   me_far, 1)
    v3.set_parity(c_cen, 2)

    cmol = to_cdk(v3)

    from jpype import JClass
    ET = JClass('org.openscience.cdk.stereo.ExtendedTetrahedral')
    et_elements = [se for se in cmol.stereoElements() if isinstance(se, ET)]
    assert len(et_elements) == 1

    back = from_cdk(cmol)
    units = [u for u in back.stereo_units() if u['kind'] == 2]
    assert len(units) == 1
    assert units[0]['parity'] != 0


@requires_jpype
def test_from_cdk_rejects_non_container():
    """from_cdk raises DirectionNotImplemented for a non-IAtomContainer."""
    with pytest.raises(DirectionNotImplemented):
        from_cdk(object())


@requires_jpype
def test_from_cdk_rejects_string():
    """from_cdk raises DirectionNotImplemented for a plain string."""
    with pytest.raises(DirectionNotImplemented):
        from_cdk('CC')


@requires_jpype
def test_dispatch_export_calls_to_cdk(cdk):
    """interop.cdk(V3Molecule()) reaches the export direction and returns an IAtomContainer."""
    from jpype import JClass
    IAtomContainer = JClass('org.openscience.cdk.interfaces.IAtomContainer')
    mol = V3Molecule()
    mol.add_atom(6)
    result = interop_cdk(mol)
    assert isinstance(result, IAtomContainer)


@requires_jpype
def test_dispatch_import_calls_from_cdk(cdk):
    """interop.cdk(IAtomContainer) reaches the import direction and returns a V3 molecule."""
    mol = V3Molecule()
    mol.add_atom(6)
    cmol = to_cdk(mol)
    back = interop_cdk(cmol)
    assert isinstance(back, V3Molecule)
