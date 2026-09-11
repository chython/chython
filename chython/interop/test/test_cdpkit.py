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
CDPKit exporter tests, export direction only.  Skipped when ``CDPL`` is absent.

Stereo is checked against an InChI oracle on non-aromatic molecules only: CDPKit's InChI writer
derives hydrogen counts from raw bond orders, and returns ``None`` for some simple topologies.
"""
import pytest

from .conftest import requires_cdpkit

# Optional InChI oracle; libinchi must be bundled and loadable.
try:
    from chython.core import inchi_library_loaded as _inchi_loaded, molecule_to_inchi
    _INCHI = _inchi_loaded()
except ImportError:
    _INCHI = False
    molecule_to_inchi = None

needs_inchi = pytest.mark.skipif(not _INCHI, reason='libinchi not loaded; InChI oracle absent')


def _v3_methanol():
    """CH3OH."""
    from chython.core import MoleculeContainer
    mol = MoleculeContainer()
    c = mol.add_atom(6, implicit_h=3)
    o = mol.add_atom(8, implicit_h=1)
    mol.add_bond(c, o, 1)
    return mol, c, o


def _v3_benzene():
    """Benzene ring with aromatic bonds (order 4)."""
    from chython.core import MoleculeContainer
    mol = MoleculeContainer()
    atoms = [mol.add_atom(6, implicit_h=1) for _ in range(6)]
    for i in range(6):
        mol.add_bond(atoms[i], atoms[(i + 1) % 6], 4)
    return mol, atoms


def _v3_chfclbr(parity):
    """
    CHFClBr with a given parity: 1 -> InChI /m1 (S), 2 -> /m0 (R).

    Refs are (F, Cl, Br, None); the implicit H maps to the centre atom in CDPKit.
    """
    from chython.core import MoleculeContainer
    mol = MoleculeContainer()
    c = mol.add_atom(6, implicit_h=1)
    f = mol.add_atom(9, implicit_h=0)
    cl = mol.add_atom(17, implicit_h=0)
    br = mol.add_atom(35, implicit_h=0)
    mol.add_bond(c, f, 1)
    mol.add_bond(c, cl, 1)
    mol.add_bond(c, br, 1)
    mol.set_parity(c, parity)
    return mol


def _v3_but2ene(parity):
    """But-2-ene CH3-CH=CH-CH3; parity 1 -> E (TRANS), 2 -> Z (CIS)."""
    from chython.core import MoleculeContainer
    mol = MoleculeContainer()
    c0 = mol.add_atom(6, implicit_h=3)
    c1 = mol.add_atom(6, implicit_h=1)
    c2 = mol.add_atom(6, implicit_h=1)
    c3 = mol.add_atom(6, implicit_h=3)
    mol.add_bond(c0, c1, 1)
    mol.add_bond(c1, c2, 2)
    mol.add_bond(c2, c3, 1)
    mol.set_parity(c1, parity)
    return mol


def _v3_allene():
    """
    1,3-Dichloroallene ClHC=C=CHCl and its centre atom; kind 2, unrepresentable in CDPKit.

    Dichloro rather than the parent allene: identical termini are not stereogenic.
    """
    from chython.core import MoleculeContainer
    mol = MoleculeContainer()
    c0 = mol.add_atom(6, implicit_h=1)      # ClCH= (left terminal)
    cl0 = mol.add_atom(17, implicit_h=0)
    c1 = mol.add_atom(6, implicit_h=0)      # =C= (allene centre)
    c2 = mol.add_atom(6, implicit_h=1)      # =CHCl (right terminal)
    cl2 = mol.add_atom(17, implicit_h=0)
    mol.add_bond(c0, cl0, 1)
    mol.add_bond(c0, c1, 2)
    mol.add_bond(c1, c2, 2)
    mol.add_bond(c2, cl2, 1)
    return mol, c1


@requires_cdpkit
def test_from_cdpkit_raises_direction_not_implemented():
    """
    The import direction raises ``DirectionNotImplemented`` and names itself as export-only.

    Export-only is permanent, not a stub, which is what ``match='export only'`` tests.
    """
    from chython.interop import cdpkit
    from chython.exceptions import DirectionNotImplemented

    with pytest.raises(DirectionNotImplemented, match='export only'):
        cdpkit(object())


@requires_cdpkit
def test_to_cdpkit_returns_basic_molecule():
    """``to_cdpkit`` returns a CDPKit ``BasicMolecule``."""
    import CDPL.Chem as Chem
    from chython.interop import cdpkit

    mol, *_ = _v3_methanol()
    result = cdpkit(mol)
    assert isinstance(result, Chem.BasicMolecule)


@requires_cdpkit
def test_atom_order_matches_iteration_order_v3():
    """CDPKit index *i* corresponds to the *i*-th atom in ``mol.atom_numbers``."""
    import CDPL.Chem as Chem
    from chython.interop import cdpkit

    mol, c_num, o_num = _v3_methanol()
    cdp = cdpkit(mol)

    nums = list(mol.atom_numbers)

    a0 = cdp.getAtom(0)
    a1 = cdp.getAtom(1)
    if nums[0] == c_num:
        assert Chem.getType(a0) == 6   # carbon
        assert Chem.getType(a1) == 8   # oxygen
    else:
        assert Chem.getType(a0) == 8
        assert Chem.getType(a1) == 6


@requires_cdpkit
def test_element_types_are_transferred():
    """Atomic numbers are stored on CDPKit atoms as the atom type."""
    import CDPL.Chem as Chem
    from chython.interop import cdpkit

    mol, c_num, o_num = _v3_methanol()
    cdp = cdpkit(mol)

    types = {Chem.getType(cdp.getAtom(i)) for i in range(cdp.numAtoms)}
    assert 6 in types   # carbon
    assert 8 in types   # oxygen


@requires_cdpkit
def test_formal_charge_is_transferred():
    """Positive and negative formal charges must survive the round-trip."""
    import CDPL.Chem as Chem
    from chython.interop import cdpkit
    from chython.core import MoleculeContainer

    mol = MoleculeContainer()
    mol.add_atom(7, charge=1, implicit_h=4)   # NH4+
    cdp = cdpkit(mol)
    assert Chem.getFormalCharge(cdp.getAtom(0)) == 1

    mol2 = MoleculeContainer()
    mol2.add_atom(8, charge=-1, implicit_h=0)  # O-
    cdp2 = cdpkit(mol2)
    assert Chem.getFormalCharge(cdp2.getAtom(0)) == -1


@requires_cdpkit
def test_isotope_is_transferred():
    """Non-zero isotope labels must be present in the CDPKit atom."""
    import CDPL.Chem as Chem
    from chython.interop import cdpkit
    from chython.core import MoleculeContainer

    mol = MoleculeContainer()
    mol.add_atom(6, isotope=13, implicit_h=4)  # 13CH4
    cdp = cdpkit(mol)
    assert Chem.getIsotope(cdp.getAtom(0)) == 13


@requires_cdpkit
def test_radical_is_transferred():
    """Radical atoms must carry a DOUBLET radical type in CDPKit."""
    import CDPL.Chem as Chem
    from chython.interop import cdpkit
    from chython.core import MoleculeContainer

    mol = MoleculeContainer()
    mol.add_atom(6, radical=True, implicit_h=3)  # methyl radical
    cdp = cdpkit(mol)
    assert Chem.getRadicalType(cdp.getAtom(0)) == Chem.RadicalType.DOUBLET


@requires_cdpkit
def test_implicit_hydrogen_count_is_transferred():
    """Explicit implicit-H values must be set on CDPKit atoms."""
    import CDPL.Chem as Chem
    from chython.interop import cdpkit
    from chython.core import MoleculeContainer

    mol = MoleculeContainer()
    mol.add_atom(6, implicit_h=4)  # methane
    cdp = cdpkit(mol)
    assert Chem.getImplicitHydrogenCount(cdp.getAtom(0)) == 4


@requires_cdpkit
def test_h_unknown_appends_to_log():
    """Atoms with an unknown implicit-H count trigger exactly one log message per molecule."""
    from chython.interop import cdpkit
    from chython.core import MoleculeContainer

    mol = MoleculeContainer()
    mol.add_atom(6)   # implicit_h defaults to H_UNKNOWN
    mol.add_atom(6)   # a second unknown atom must not add a second log entry

    log = []
    cdpkit(mol, log=log)

    assert len(log) == 1
    assert 'hydrogen' in str(log[0]).lower()


@requires_cdpkit
def test_h_unknown_no_log_when_log_is_none():
    """No-op when ``log=None`` (the default): the call must not crash."""
    from chython.interop import cdpkit
    from chython.core import MoleculeContainer

    mol = MoleculeContainer()
    mol.add_atom(6)   # H_UNKNOWN — log=None by default
    cdpkit(mol)       # must not raise


@requires_cdpkit
def test_aromatic_bonds_set_aromaticity_flag():
    """
    Order-4 bonds in chython mark both atoms and bonds aromatic in CDPKit.

    Checked at the flag level: CDPKit's InChI writer derives hydrogen counts from raw bond orders.
    """
    import CDPL.Chem as Chem
    from chython.interop import cdpkit

    mol, _ = _v3_benzene()
    cdp = cdpkit(mol)

    for i in range(cdp.numAtoms):
        assert Chem.getAromaticityFlag(cdp.getAtom(i)), f'atom {i} is not marked aromatic'

    for i in range(cdp.numBonds):
        b = cdp.getBond(i)
        assert Chem.getAromaticityFlag(b), f'bond {i} is not marked aromatic'


@requires_cdpkit
def test_aromatic_smiles_contains_lowercase():
    """CDPKit should generate a SMILES with lowercase aromatic atoms for benzene."""
    import CDPL.Chem as Chem
    from chython.interop import cdpkit

    mol, _ = _v3_benzene()
    cdp = cdpkit(mol)
    smi = Chem.generateSMILES(cdp)

    assert smi is not None
    assert 'c' in smi, f'expected aromatic SMILES, got {smi!r}'


@needs_inchi
@requires_cdpkit
def test_tetrahedral_stereo_parity1_matches_inchi():
    """
    Parity 1 on a V3 stereocentre maps to the same stereo descriptor as chython's own InChI.

    CHFClBr: all substituents are single atoms, so CDPKit's InChI writer works on it.
    """
    import CDPL.Chem as Chem
    from chython.interop import cdpkit

    mol = _v3_chfclbr(1)
    cdp = cdpkit(mol)
    inchi_cdp = Chem.generateINCHI(cdp)
    inchi_chy = molecule_to_inchi(mol)

    assert inchi_cdp is not None, 'CDPKit returned None — topology not supported'
    assert inchi_cdp == inchi_chy


@needs_inchi
@requires_cdpkit
def test_tetrahedral_stereo_parity2_matches_inchi():
    """Parity 2 (the other enantiomer) also matches."""
    import CDPL.Chem as Chem
    from chython.interop import cdpkit

    mol = _v3_chfclbr(2)
    cdp = cdpkit(mol)
    inchi_cdp = Chem.generateINCHI(cdp)
    inchi_chy = molecule_to_inchi(mol)

    assert inchi_cdp is not None
    assert inchi_cdp == inchi_chy


@needs_inchi
@requires_cdpkit
def test_tetrahedral_enantiomers_differ():
    """The two parities must produce distinct InChIs (different /m layer)."""
    import CDPL.Chem as Chem
    from chython.interop import cdpkit

    cdp1 = cdpkit(_v3_chfclbr(1))
    cdp2 = cdpkit(_v3_chfclbr(2))

    inchi1 = Chem.generateINCHI(cdp1)
    inchi2 = Chem.generateINCHI(cdp2)

    assert inchi1 is not None and inchi2 is not None
    assert inchi1 != inchi2


@needs_inchi
@requires_cdpkit
def test_cistrans_parity1_is_E():
    """Parity 1 on a cis/trans unit in V3 is TRANS (E); InChI layer ends with ``/b..+``."""
    import CDPL.Chem as Chem
    from chython.interop import cdpkit

    mol = _v3_but2ene(1)
    cdp = cdpkit(mol)
    inchi_cdp = Chem.generateINCHI(cdp)
    inchi_chy = molecule_to_inchi(mol)

    assert inchi_cdp is not None, 'CDPKit returned None'
    assert inchi_cdp == inchi_chy
    assert inchi_cdp.endswith('+'), f'expected /b..+ for E isomer, got {inchi_cdp!r}'


@needs_inchi
@requires_cdpkit
def test_cistrans_parity2_is_Z():
    """Parity 2 is CIS (Z); InChI layer ends with ``/b..-``."""
    import CDPL.Chem as Chem
    from chython.interop import cdpkit

    mol = _v3_but2ene(2)
    cdp = cdpkit(mol)
    inchi_cdp = Chem.generateINCHI(cdp)
    inchi_chy = molecule_to_inchi(mol)

    assert inchi_cdp is not None
    assert inchi_cdp == inchi_chy
    assert inchi_cdp.endswith('-'), f'expected /b..- for Z isomer, got {inchi_cdp!r}'


@needs_inchi
@requires_cdpkit
def test_cistrans_isomers_differ():
    """E and Z must give different InChIs."""
    import CDPL.Chem as Chem
    from chython.interop import cdpkit

    i1 = Chem.generateINCHI(cdpkit(_v3_but2ene(1)))
    i2 = Chem.generateINCHI(cdpkit(_v3_but2ene(2)))

    assert i1 is not None and i2 is not None
    assert i1 != i2


@needs_inchi
@requires_cdpkit
def test_cistrans_with_a_hydrogen_on_one_side_is_exported():
    """A double bond with an implicit hydrogen on one side keeps its configuration."""
    import CDPL.Chem as Chem
    from chython.core import read_smiles
    from chython.interop import cdpkit

    inchis = []
    for smi in (r'C/C=C/Cl', r'C/C=C\Cl'):
        log = []
        cdp = cdpkit(read_smiles(smi), log=log)
        assert not [line for line in log if 'cis/trans' in line], log
        inchi = Chem.generateINCHI(cdp)
        assert inchi is not None
        assert '/b' in inchi, f'configuration was dropped for {smi}: {inchi!r}'
        inchis.append(inchi)
    assert inchis[0] != inchis[1], 'E and Z came out identical'


@requires_cdpkit
def test_an_unrepresentable_kind_is_always_reported():
    """Whatever stereo kind the core grows next is reported, not dropped in silence."""
    import CDPL.Chem as Chem
    from chython.core import read_smiles
    from chython.interop._cdpkit import _to_cdpkit_v3

    class _UnknownKind:
        """A molecule reporting one unit of a kind this converter has no branch for."""

        def __init__(self, real):
            self._real = real

        def __getattr__(self, name):
            return getattr(self._real, name)

        def stereo_units(self):
            return [{'kind': 99, 'parity': 1, 'anchor': self._real.atom_numbers[0],
                     'refs': (None, None, None, None), 'stereogenic': True}]

    log = []
    _to_cdpkit_v3(_UnknownKind(read_smiles('CCO')), log, Chem)
    assert any('stereo kind 99' in line and 'skipped' in line for line in log), log


@requires_cdpkit
def test_allene_stereo_logs_loss():
    """
    Allene stereo (kind 2) cannot be encoded in CDPKit; the loss is logged, not raised.

    ``to_cdpkit`` must still return a molecule.
    """
    import CDPL.Chem as Chem
    from chython.interop import cdpkit

    mol, centre = _v3_allene()
    mol.set_parity(centre, 1)

    log = []
    cdp = cdpkit(mol, log=log)

    assert isinstance(cdp, Chem.BasicMolecule), 'expected a molecule despite unsupported stereo'
    assert log, 'expected at least one log entry for unsupported allene stereo'
    assert any('allene' in str(entry).lower() for entry in log), f'no allene mention in log: {log}'
