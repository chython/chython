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
Tests for the OpenBabel converter, both directions.
"""
import pytest

from chython.core import MoleculeContainer, H_UNKNOWN
from chython.core._core import read_smiles, write_smiles
from chython.exceptions import UnconvertibleType

from .conftest import requires_openbabel


def _v3_mol(*args_unused):
    """Return a V3 MoleculeContainer built manually (methanol: C-O)."""
    mol = MoleculeContainer()
    with mol.edit():
        c = mol.add_atom(6, implicit_h=3)
        o = mol.add_atom(8, implicit_h=1)
        mol.add_bond(c, o, 1)
    return mol


def _read(smiles_str):
    return read_smiles(smiles_str)


@requires_openbabel
def test_atom_order_preserved():
    """OBMol index i must correspond to the i-th atom in mol.atoms() iteration order."""
    from chython.interop._openbabel import to_openbabel

    mol = _read('CCO')  # C-C-O: three atoms, clear element sequence
    original_elements = [a.element for a in mol.atoms()]

    ob_mol = to_openbabel(mol)
    ob_elements = [ob_mol.GetAtom(i).GetAtomicNum() for i in range(1, ob_mol.NumAtoms() + 1)]

    assert ob_elements == original_elements, (
        f'atom order not preserved: original {original_elements}, OBMol {ob_elements}'
    )


@requires_openbabel
def test_constitution_round_trip_ethanol():
    """Ethanol constitution survives a V3 → OBMol → V3 round-trip."""
    from chython.interop._openbabel import to_openbabel, from_openbabel

    mol = _read('CCO')
    ob_mol = to_openbabel(mol)
    mol2 = from_openbabel(ob_mol)

    original = [(a.element, a.implicit_h) for a in mol.atoms()]
    result = [(a.element, a.implicit_h) for a in mol2.atoms()]
    assert original == result

    bonds_in = sorted((b.order, frozenset({b.n, b.m})) for b in mol.bonds())
    bonds_out = sorted((b.order, frozenset({b.n, b.m})) for b in mol2.bonds())
    assert len(bonds_in) == len(bonds_out)
    assert all(o1 == o2 for (o1, _), (o2, _) in zip(bonds_in, bonds_out))


@requires_openbabel
def test_isotope_and_charge_round_trip():
    """Isotope and formal charge survive V3 → OBMol → V3."""
    from chython.interop._openbabel import to_openbabel, from_openbabel

    mol = MoleculeContainer()
    with mol.edit():
        sid = mol.add_atom(6, isotope=13, charge=1, implicit_h=3)

    ob_mol = to_openbabel(mol)
    oa = ob_mol.GetAtom(1)
    assert oa.GetIsotope() == 13
    assert oa.GetFormalCharge() == 1

    mol2 = from_openbabel(ob_mol)
    atom2 = next(iter(mol2.atoms()))
    assert atom2.isotope == 13
    assert atom2.charge == 1


@requires_openbabel
def test_radical_round_trip():
    """SpinMultiplicity round-trips; non-zero OBMol spin becomes radical=True in V3."""
    from chython.interop._openbabel import to_openbabel, from_openbabel

    mol = MoleculeContainer()
    with mol.edit():
        mol.add_atom(6, radical=True, implicit_h=3)

    ob_mol = to_openbabel(mol)
    assert ob_mol.GetAtom(1).GetSpinMultiplicity() != 0

    mol2 = from_openbabel(ob_mol)
    assert next(iter(mol2.atoms())).is_radical


@requires_openbabel
def test_aromatic_bond_order_4_sets_aromatic_flag():
    """A stored order-4 bond comes out as OBMol's aromatic flag plus bond order 1."""
    from chython.interop._openbabel import to_openbabel

    mol = _read('c1ccccc1')  # benzene, all bonds stored as order 4
    ob_mol = to_openbabel(mol)

    for i in range(ob_mol.NumBonds()):
        b = ob_mol.GetBond(i)
        assert b.IsAromatic(), f'benzene bond {i} not marked aromatic in OBMol'
        assert b.GetBondOrder() == 1, (
            f'benzene bond {i}: expected bond order 1 (aromatic as flag), got {b.GetBondOrder()}'
        )


@requires_openbabel
def test_aromatic_bond_imports_as_order_4():
    """OBMol aromatic bonds must be stored as order 4 in V3, not re-perceived."""
    from chython.interop._openbabel import to_openbabel, from_openbabel
    from openbabel import openbabel as ob

    mol = _read('c1ccccc1')
    ob_mol = to_openbabel(mol)
    mol2 = from_openbabel(ob_mol)

    bond_orders = sorted(b.order for b in mol2.bonds())
    assert all(o == 4 for o in bond_orders), (
        f'benzene bonds not stored as order 4 after import: {bond_orders}'
    )


@requires_openbabel
def test_h_unknown_is_logged():
    """An atom with H_UNKNOWN must produce a log entry when exported."""
    from chython.interop._openbabel import to_openbabel

    mol = MoleculeContainer()
    with mol.edit():
        c_id = mol.add_atom(6)           # implicit_h not given → H_UNKNOWN
        o_id = mol.add_atom(8, implicit_h=1)
        mol.add_bond(c_id, o_id, 1)

    assert mol.implicit_h_of(c_id) is None, 'add_atom without implicit_h should give H_UNKNOWN'

    log = []
    to_openbabel(mol, log=log)
    assert any('H_UNKNOWN' in entry for entry in log), (
        f'expected H_UNKNOWN log entry, got: {log}'
    )


@requires_openbabel
def test_h_unknown_not_logged_when_stated():
    """No H_UNKNOWN log entry when implicit_h is explicitly given."""
    from chython.interop._openbabel import to_openbabel

    mol = MoleculeContainer()
    with mol.edit():
        c_id = mol.add_atom(6, implicit_h=4)
        o_id = mol.add_atom(8, implicit_h=0)
        mol.add_bond(c_id, o_id, 2)

    log = []
    to_openbabel(mol, log=log)
    assert not any('H_UNKNOWN' in entry for entry in log), (
        f'unexpected H_UNKNOWN log entry: {log}'
    )


@requires_openbabel
def test_dative_bond_is_logged_and_stored_as_single():
    """A dative bond (order 8) must be logged and stored as a single bond in OBMol."""
    from chython.interop._openbabel import to_openbabel

    mol = MoleculeContainer()
    with mol.edit():
        n_id = mol.add_atom(7, implicit_h=0)
        fe_id = mol.add_atom(26, implicit_h=0)
        mol.add_bond(n_id, fe_id, 8)

    log = []
    ob_mol = to_openbabel(mol, log=log)

    assert any('dative' in entry or 'order-8' in entry for entry in log), (
        f'expected dative bond log entry, got: {log}'
    )
    assert ob_mol.GetBond(0).GetBondOrder() == 1


@requires_openbabel
def test_allene_stereo_is_logged():
    """Allene/cumulene stereo with a configured parity must produce a log entry."""
    from chython.interop._openbabel import to_openbabel

    # F-CH=C=CH-F
    mol = MoleculeContainer()
    with mol.edit():
        f1 = mol.add_atom(9, implicit_h=0)
        c1 = mol.add_atom(6, implicit_h=1)
        ca = mol.add_atom(6, implicit_h=0)
        c2 = mol.add_atom(6, implicit_h=1)
        f2 = mol.add_atom(9, implicit_h=0)
        mol.add_bond(f1, c1, 1)
        mol.add_bond(c1, ca, 2)
        mol.add_bond(ca, c2, 2)
        mol.add_bond(c2, f2, 1)

    units = mol.stereogenic_units()
    allene_units = [u for u in units if u['kind'] == 2]
    if not allene_units:
        pytest.skip('molecule not recognised as allene stereocentre')

    anchor = allene_units[0]['anchor']
    mol.set_parity(anchor, 2)

    log = []
    to_openbabel(mol, log=log)
    assert any('allene' in entry or 'cumulene' in entry for entry in log), (
        f'expected allene log entry, got: {log}'
    )


@requires_openbabel
def test_allene_no_log_when_unconfigured():
    """No log entry for an allene with parity 0 (unconfigured)."""
    from chython.interop._openbabel import to_openbabel

    mol = MoleculeContainer()
    with mol.edit():
        f1 = mol.add_atom(9, implicit_h=0)
        c1 = mol.add_atom(6, implicit_h=1)
        ca = mol.add_atom(6, implicit_h=0)
        c2 = mol.add_atom(6, implicit_h=1)
        f2 = mol.add_atom(9, implicit_h=0)
        mol.add_bond(f1, c1, 1)
        mol.add_bond(c1, ca, 2)
        mol.add_bond(ca, c2, 2)
        mol.add_bond(c2, f2, 1)

    # parity stays 0 -- nothing to export, nothing to log
    log = []
    to_openbabel(mol, log=log)
    assert not any('allene' in entry or 'cumulene' in entry for entry in log), (
        f'unexpected allene log entry for unconfigured stereo: {log}'
    )


@requires_openbabel
def test_tetrahedral_stereo_round_trip_smiles():
    """[C@@H](F)(Cl)Br survives V3 -> OBMol -> V3 with parity intact."""
    from chython.interop._openbabel import to_openbabel, from_openbabel

    smiles_in = '[C@@H](F)(Cl)Br'  # parity 2 = @ = AntiClockwise
    mol = _read(smiles_in)

    units = mol.stereogenic_units()
    assert units, 'expected a TH stereocentre'
    assert units[0]['parity'] == 2, f'expected parity 2, got {units[0]["parity"]}'

    ob_mol = to_openbabel(mol)
    mol2 = from_openbabel(ob_mol)

    units2 = mol2.stereogenic_units()
    assert units2, 'stereocentre lost on round-trip'
    assert units2[0]['parity'] == 2, (
        f'parity changed: expected 2, got {units2[0]["parity"]}'
    )

    assert write_smiles(mol) == write_smiles(mol2), (
        f'SMILES mismatch: {write_smiles(mol)!r} vs {write_smiles(mol2)!r}'
    )


@requires_openbabel
def test_tetrahedral_stereo_opposite_round_trip_smiles():
    """[C@H](F)(Cl)Br (parity 1) also round-trips correctly."""
    from chython.interop._openbabel import to_openbabel, from_openbabel

    mol = _read('[C@H](F)(Cl)Br')
    units = mol.stereogenic_units()
    assert units and units[0]['parity'] == 1

    mol2 = from_openbabel(to_openbabel(mol))
    units2 = mol2.stereogenic_units()
    assert units2 and units2[0]['parity'] == 1, (
        f'parity changed: expected 1, got {units2[0]["parity"] if units2 else "none"}'
    )
    assert write_smiles(mol) == write_smiles(mol2)


@requires_openbabel
def test_two_distinct_stereocentres_preserved():
    """Both TH centres in (R,R)-2,3-dichlorobutane survive the round-trip."""
    from chython.interop._openbabel import to_openbabel, from_openbabel

    mol = _read('C[C@@H](Cl)[C@@H](Cl)C')  # (R,R)
    units = {u['anchor']: u['parity'] for u in mol.stereogenic_units()}
    assert len(units) == 2

    mol2 = from_openbabel(to_openbabel(mol))
    units2 = {u['anchor']: u['parity'] for u in mol2.stereogenic_units()}
    assert len(units2) == 2

    assert write_smiles(mol) == write_smiles(mol2), (
        f'SMILES mismatch: {write_smiles(mol)!r} vs {write_smiles(mol2)!r}'
    )


@requires_openbabel
def test_trans_alkene_round_trip():
    """(E)-1,2-difluoroethylene survives V3 -> OBMol -> V3."""
    from chython.interop._openbabel import to_openbabel, from_openbabel

    mol = _read('F/C=C/F')  # trans; parity 1
    units = mol.stereogenic_units()
    assert units and units[0]['parity'] == 1

    mol2 = from_openbabel(to_openbabel(mol))
    units2 = mol2.stereogenic_units()
    assert units2 and units2[0]['parity'] == 1, (
        f'trans parity changed: expected 1, got {units2[0]["parity"] if units2 else "none"}'
    )
    assert write_smiles(mol) == write_smiles(mol2)


@requires_openbabel
def test_cis_alkene_round_trip():
    """(Z)-1,2-difluoroethylene survives V3 → OBMol → V3."""
    from chython.interop._openbabel import to_openbabel, from_openbabel

    mol = _read(r'F/C=C\F')  # cis; parity 2
    units = mol.stereogenic_units()
    assert units and units[0]['parity'] == 2

    mol2 = from_openbabel(to_openbabel(mol))
    units2 = mol2.stereogenic_units()
    assert units2 and units2[0]['parity'] == 2, (
        f'cis parity changed: expected 2, got {units2[0]["parity"] if units2 else "none"}'
    )
    assert write_smiles(mol) == write_smiles(mol2)


@requires_openbabel
def test_cis_and_trans_distinguished():
    """The two isomers must produce distinct canonical SMILES after round-trip."""
    from chython.interop._openbabel import to_openbabel, from_openbabel

    mol_cis = _read(r'F/C=C\F')
    mol_trans = _read('F/C=C/F')

    smi_cis = write_smiles(from_openbabel(to_openbabel(mol_cis)))
    smi_trans = write_smiles(from_openbabel(to_openbabel(mol_trans)))

    assert smi_cis != smi_trans, (
        f'cis and trans round-tripped to same SMILES: {smi_cis!r}'
    )


@requires_openbabel
def test_trans_but2ene_round_trip():
    """(E)-but-2-ene (two methyl groups) round-trips correctly."""
    from chython.interop._openbabel import to_openbabel, from_openbabel

    mol = _read('C/C=C/C')  # trans-2-butene, parity 1
    units = mol.stereogenic_units()
    if not units:
        pytest.skip('but-2-ene not parsed as stereocentre')

    mol2 = from_openbabel(to_openbabel(mol))
    assert write_smiles(mol) == write_smiles(mol2)


@requires_openbabel
def test_from_openbabel_wrong_type_raises():
    """from_openbabel raises UnconvertibleType for anything that is not an OBMol."""
    from chython.interop._openbabel import from_openbabel

    with pytest.raises(UnconvertibleType):
        from_openbabel('not an OBMol')

    with pytest.raises(UnconvertibleType):
        from_openbabel(42)

    with pytest.raises(UnconvertibleType):
        from_openbabel(None)


@requires_openbabel
def test_to_openbabel_wrong_type_raises():
    """to_openbabel raises UnconvertibleType for non-chython containers."""
    from chython.interop._openbabel import to_openbabel

    with pytest.raises(UnconvertibleType):
        to_openbabel('CCO')

    with pytest.raises(UnconvertibleType):
        to_openbabel(42)

    with pytest.raises(UnconvertibleType):
        to_openbabel(None)


@requires_openbabel
def test_empty_molecule_round_trip():
    """An empty MoleculeContainer must survive the round-trip without error."""
    from chython.interop._openbabel import to_openbabel, from_openbabel

    mol = MoleculeContainer()
    ob_mol = to_openbabel(mol)
    assert ob_mol.NumAtoms() == 0

    mol2 = from_openbabel(ob_mol)
    assert sum(1 for _ in mol2.atoms()) == 0
