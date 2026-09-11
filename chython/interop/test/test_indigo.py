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
Tests for the Indigo converter, both directions, with an InChI oracle.

Import uses `stereo_units()` and not `stereogenic_units()`, so an atom perceived as topologically
symmetric before any parity is set is still found.
"""
from pathlib import Path
from importlib.util import find_spec

import pytest

from .conftest import requires_indigo


# Helpers loaded lazily so the module imports cleanly without Indigo installed.

def _try_load_libinchi():
    """Load the bundled libinchi from any candidate path; idempotent, returns availability."""
    from chython.core._core import ich_load_library, inchi_library_loaded
    if inchi_library_loaded():
        return True
    # Normally already loaded by `core/__init__.py` from beside itself.  The search below is for a
    # worktree, which carries the source but not the git-ignored binary, so it also looks in the main
    # checkout the worktree branched from.
    names = ('libinchi.dylib', 'libinchi.so', 'libinchi.dll')
    _core = Path(__file__).parent.parent.parent / 'core'
    _repo = Path(__file__).parent.parent.parent.parent           # repo (or worktree) root
    # typically .claude/worktrees/<hash>/ -> the real checkout sits three levels above it
    _common = _repo.parent.parent.parent / 'chython' / 'chython' / 'core'
    candidates = [d / n for d in (_core, _repo / 'chython' / 'core', _common) for n in names]
    for c in candidates:
        if c.exists():
            try:
                ich_load_library(str(c))
                if inchi_library_loaded():
                    return True
            except Exception:
                pass
    return inchi_library_loaded()


def _inchi(mol):
    """Compute an InChI, skipping the test if libinchi is not available."""
    from chython.core._core import molecule_to_inchi, inchi_library_loaded
    if not inchi_library_loaded() and not _try_load_libinchi():
        pytest.skip('libinchi not loaded; install or set INCHI_PATH')
    return molecule_to_inchi(mol)


def _require_inchi():
    """Skip the calling test if InChI is unavailable."""
    from chython.core._core import inchi_library_loaded
    if not inchi_library_loaded() and not _try_load_libinchi():
        pytest.skip('libinchi not loaded; install or set INCHI_PATH')


@requires_indigo
def test_basic_export_and_import():
    """Ethanol survives a round-trip as the same InChI."""
    from chython.core._core import read_smiles
    from chython.interop._indigo import to_indigo, from_indigo

    mol = read_smiles('CCO')
    ig = to_indigo(mol)
    mol_rt = from_indigo(ig)
    assert _inchi(mol) == _inchi(mol_rt)


@requires_indigo
@pytest.mark.parametrize('smi', [
    '[NH4+]',             # positive charge
    '[O-]',               # negative charge
    '[13CH4]',            # isotope
    '[CH3][13CH3]',       # isotope in chain
])
def test_constitution_properties(smi):
    """Charge and isotope survive the round-trip."""
    from chython.core._core import read_smiles
    from chython.interop._indigo import to_indigo, from_indigo

    mol = read_smiles(smi)
    ig = to_indigo(mol)
    mol_rt = from_indigo(ig)
    assert _inchi(mol) == _inchi(mol_rt)


@requires_indigo
def test_radical_round_trip():
    """
    A mono-radical is exported as Indigo doublet (102) and imported back as a radical.

    102 and not 2: Indigo rejects `setRadical(2)` with 'Unknown radical type'.
    """
    from chython.core._core import MoleculeContainer
    from chython.interop._indigo import to_indigo, from_indigo

    mol = MoleculeContainer()
    n = mol.add_atom(6, radical=True, implicit_h=3)  # methyl radical
    ig = to_indigo(mol)

    ig_radical = ig.getAtom(0).radical()
    assert ig_radical == 102, f'Expected Indigo doublet (102), got {ig_radical}'

    mol_rt = from_indigo(ig)
    rt_atoms = list(mol_rt.atoms())
    assert rt_atoms[0].is_radical, 'Radical flag lost on import'


@requires_indigo
def test_atom_order_is_interface():
    """
    `iterateAtoms()` on the exported Indigo molecule follows the V3 iteration order.

    `depict/layout/molecule.py` maps 2D coordinates back by position.
    """
    from chython.core._core import read_smiles
    from chython.interop._indigo import to_indigo

    mol = read_smiles('CCN(CC)CC')  # triethylamine: C, C, N, C, C, C, C
    ig = to_indigo(mol)

    v3_elements = [a.element for a in mol.atoms()]
    ig_symbols = [a.symbol() for a in ig.iterateAtoms()]
    from chython.core._core import element_symbols
    syms = element_symbols()
    v3_symbols = [syms[e] for e in v3_elements]

    assert v3_symbols == ig_symbols, (
        f'Atom order mismatch: V3 {v3_symbols} vs Indigo {ig_symbols}'
    )


@requires_indigo
def test_aromatic_bonds_stay_aromatic():
    """Stored aromatic bonds (order 4) go out as aromatic, not kekulized."""
    from chython.core._core import read_smiles
    from chython.interop._indigo import to_indigo

    mol = read_smiles('c1ccccc1')  # benzene: all bonds stored as order 4
    ig = to_indigo(mol)

    orders = [b.bondOrder() for b in ig.iterateBonds()]
    assert all(o == 4 for o in orders), (
        f'Expected all aromatic (order 4) bonds; got {orders}'
    )


@requires_indigo
def test_kekule_bonds_stay_kekule():
    """A molecule stored in Kekule form (alternating 1/2 bonds) goes out in Kekule form."""
    from chython.core._core import MoleculeContainer
    from chython.interop._indigo import to_indigo

    mol = MoleculeContainer()
    sids = [mol.add_atom(6) for _ in range(6)]
    for i in range(6):
        order = 2 if i % 2 == 0 else 1
        mol.add_bond(sids[i], sids[(i + 1) % 6], order)

    ig = to_indigo(mol)
    orders = sorted(b.bondOrder() for b in ig.iterateBonds())
    assert orders.count(1) == 3 and orders.count(2) == 3, (
        f'Expected 3 single + 3 double bonds; got {orders}'
    )


@requires_indigo
def test_h_unknown_is_logged_not_silent():
    """
    H_UNKNOWN is logged, never silently 0: unset and stated zero are different molecules.

    The count in the Indigo object is left to Indigo's valence perception rather than forced.
    """
    from chython.core._core import MoleculeContainer, H_UNKNOWN
    from chython.interop._indigo import to_indigo

    mol = MoleculeContainer()
    n = mol.add_atom(6, implicit_h=H_UNKNOWN)
    log = []
    to_indigo(mol, log=log)

    assert len(log) == 1, f'Expected exactly one log line; got {log}'
    assert 'H_UNKNOWN' in log[0], f'Log line does not mention H_UNKNOWN: {log[0]}'


@requires_indigo
def test_stated_zero_h_not_logged():
    """An explicitly stated zero implicit H count is representable in Indigo, so it is not logged."""
    from chython.core._core import MoleculeContainer
    from chython.interop._indigo import to_indigo

    mol = MoleculeContainer()
    mol.add_atom(6, implicit_h=0)
    log = []
    to_indigo(mol, log=log)

    assert log == [], f'Unexpected log line for stated-zero H: {log}'


@requires_indigo
def test_tetrahedral_stereo_export_enantiomers_differ():
    """Enantiomers must map to different Indigo canonical SMILES."""
    from chython.core._core import read_smiles
    from chython.interop._indigo import to_indigo

    mol_cw = read_smiles('[C@@H](F)(Cl)Br')
    mol_ccw = read_smiles('[C@H](F)(Cl)Br')
    ig_cw = to_indigo(mol_cw)
    ig_ccw = to_indigo(mol_ccw)

    assert ig_cw.canonicalSmiles() != ig_ccw.canonicalSmiles(), (
        '@@H and @H produced the same Indigo SMILES -- the parity mapping did not distinguish them'
    )


@requires_indigo
def test_tetrahedral_stereo_round_trip_inchi_cw():
    """InChI is preserved after V3 → Indigo → V3 for (R) enantiomer."""
    from chython.core._core import read_smiles
    from chython.interop._indigo import to_indigo, from_indigo

    mol = read_smiles('[C@@H](F)(Cl)Br')
    mol_rt = from_indigo(to_indigo(mol))
    assert _inchi(mol) == _inchi(mol_rt)


@requires_indigo
def test_tetrahedral_stereo_round_trip_inchi_ccw():
    """InChI is preserved after V3 → Indigo → V3 for (S) enantiomer."""
    from chython.core._core import read_smiles
    from chython.interop._indigo import to_indigo, from_indigo

    mol = read_smiles('[C@H](F)(Cl)Br')
    mol_rt = from_indigo(to_indigo(mol))
    assert _inchi(mol) == _inchi(mol_rt)


@requires_indigo
def test_tetrahedral_stereo_round_trip_inchi_multiple_centers():
    """
    A molecule with multiple tetrahedral stereocenters survives the round-trip.

    Exercises the `stereo_units()` path: one of the five reads as non-stereogenic before parity.
    """
    from chython.core._core import read_smiles
    from chython.interop._indigo import to_indigo, from_indigo

    mol = read_smiles('[C@H](Cl)([C@@H](Cl)[C@@H]([C@@H](C(O)=O)Cl)Cl)[C@@H](Cl)C(=O)O')
    mol_rt = from_indigo(to_indigo(mol))
    assert _inchi(mol) == _inchi(mol_rt)


@requires_indigo
def test_tetrahedral_stereo_enantiomers_still_differ_after_round_trip():
    """Enantiomers must still be distinguishable after a full round-trip."""
    from chython.core._core import read_smiles
    from chython.interop._indigo import to_indigo, from_indigo

    mol_cw = read_smiles('[C@@H](F)(Cl)Br')
    mol_ccw = read_smiles('[C@H](F)(Cl)Br')
    mol_cw_rt = from_indigo(to_indigo(mol_cw))
    mol_ccw_rt = from_indigo(to_indigo(mol_ccw))

    assert _inchi(mol_cw_rt) != _inchi(mol_ccw_rt), (
        'Round-tripped enantiomers have the same InChI -- the import parity did not distinguish them'
    )


@requires_indigo
def test_cis_trans_without_coordinates_is_logged():
    """Cis/trans stereo is logged as a loss when the molecule has no 2D layout."""
    from chython.core._core import read_smiles
    from chython.interop._indigo import to_indigo

    mol = read_smiles('C/C=C/C')  # trans-but-2-ene, no coordinates
    log = []
    to_indigo(mol, log=log)

    assert any('cis/trans' in line for line in log), (
        f'Expected a cis/trans loss log line; got {log}'
    )


@requires_indigo
def test_cis_trans_no_log_for_achiral_molecule():
    """No cis/trans log line for a molecule with no geometric stereo."""
    from chython.core._core import read_smiles
    from chython.interop._indigo import to_indigo

    mol = read_smiles('CC=CC')  # but-2-ene, no stereo specified
    log = []
    to_indigo(mol, log=log)

    ct_lines = [rec for rec in log if 'cis/trans' in rec]
    assert ct_lines == [], f'Unexpected cis/trans log for achiral molecule: {ct_lines}'


@requires_indigo
def test_allene_stereo_without_coordinates_is_logged():
    """An allene configuration Indigo cannot hold is reported, not dropped in silence."""
    from chython.core._core import read_smiles
    from chython.interop._indigo import to_indigo

    mol = read_smiles('NC(Br)=[C@]=C(O)C')  # 1-amino-1-bromo-3-hydroxy-buta-1,2-diene, chiral allene
    log = []
    to_indigo(mol, log=log)

    assert any('allene' in line for line in log), (
        f'An allene configuration was dropped with no log entry; got {log}'
    )


@requires_indigo
def test_allene_stereo_is_logged_even_with_coordinates():
    """
    Coordinates do not rescue an allene; the loss must still be reported.

    Indigo's `markStereobonds` reads 2D geometry for double-bond stereo and has no allene form to read
    it into, so gating the report on `has_coordinates` would be wrong.
    """
    from chython.core._core import read_smiles
    from chython.interop._indigo import to_indigo

    mol = read_smiles('NC(Br)=[C@]=C(O)C')
    for i, n in enumerate(mol.atom_numbers):
        mol.set_xy(n, float(i), 0.0)
    assert mol.has_coordinates

    log = []
    to_indigo(mol, log=log)

    assert any('allene' in line for line in log), (
        f'An allene configuration was dropped with no log entry once coordinates existed; got {log}'
    )


@requires_indigo
def test_no_allene_log_without_an_allene():
    """Negative control: a plain cis/trans molecule must not be reported as an allene loss."""
    from chython.core._core import read_smiles
    from chython.interop._indigo import to_indigo

    mol = read_smiles('C/C=C/C')
    log = []
    to_indigo(mol, log=log)

    assert not any('allene' in line for line in log), (
        f'Spurious allene loss reported for a cis/trans molecule: {log}'
    )


@requires_indigo
def test_dative_bond_is_logged_and_held():
    """
    Dative bonds (order 8) are passed to the Indigo object and a loss is logged.

    Indigo holds order-8 bonds internally but loses them on SMILES serialisation.
    """
    from chython.core._core import MoleculeContainer
    from chython.interop._indigo import to_indigo

    mol = MoleculeContainer()
    n1 = mol.add_atom(7)   # N
    n2 = mol.add_atom(26)  # Fe
    mol.add_bond(n1, n2, 8)

    log = []
    ig = to_indigo(mol, log=log)

    orders = [b.bondOrder() for b in ig.iterateBonds()]
    assert 8 in orders, f'Dative bond not passed to Indigo; got orders {orders}'

    assert any('dative' in line for line in log), (
        f'Expected a dative bond log line; got {log}'
    )


@requires_indigo
def test_from_indigo_raises_for_non_indigo_object():
    """from_indigo raises UnconvertibleType for any non-Indigo argument."""
    from chython.interop._indigo import from_indigo
    from chython.exceptions import UnconvertibleType

    with pytest.raises(UnconvertibleType):
        from_indigo(object())

    with pytest.raises(UnconvertibleType):
        from_indigo('a SMILES string')

    with pytest.raises(UnconvertibleType):
        from_indigo(42)


@requires_indigo
def test_from_indigo_accepts_indigo_molecule():
    """from_indigo does not raise UnconvertibleType for an actual Indigo object."""
    from indigo import Indigo
    from chython.interop._indigo import from_indigo

    ig = Indigo()
    mol = ig.loadMolecule('CCO')
    result = from_indigo(mol)
    assert result is not None


@requires_indigo
def test_aromatic_round_trip_benzene():
    """Benzene: aromatic bonds survive both directions."""
    from chython.core._core import read_smiles
    from chython.interop._indigo import to_indigo, from_indigo

    mol = read_smiles('c1ccccc1')
    ig = to_indigo(mol)
    mol_rt = from_indigo(ig)

    bonds_rt = [mol_rt.order_of(b.n, b.m) for b in mol_rt.bonds()]
    assert all(o == 4 for o in bonds_rt), (
        f'Expected all aromatic bonds in round-trip; got {bonds_rt}'
    )


@requires_indigo
def test_aromatic_round_trip_inchi():
    """Naphthalene InChI is preserved through the round-trip."""
    from chython.core._core import read_smiles
    from chython.interop._indigo import to_indigo, from_indigo

    mol = read_smiles('c1ccc2ccccc2c1')
    mol_rt = from_indigo(to_indigo(mol))
    assert _inchi(mol) == _inchi(mol_rt)


@requires_indigo
def test_stereo_corpus_tetrahedral_only():
    """
    Every stereo-corpus record whose configured units are all `SU_TETRA` round-trips by InChI.

    "Tetrahedral only" means every *configured* unit, not merely "no cis/trans": an allene record
    loses a `/t` entry through Indigo, which `test_allene_is_a_logged_loss` pins separately.
    """
    sdf_path = Path(__file__).parent.parent.parent.parent / 'test' / 'stereo.sdf'
    if not sdf_path.exists():
        pytest.skip('test/stereo.sdf not found')

    from chython.core import SU_TETRA
    from chython.formats.ctfile import SDFRead
    from chython.core._core import inchi_library_loaded
    from chython.interop._indigo import to_indigo, from_indigo

    if not inchi_library_loaded():
        pytest.skip('libinchi not loaded')

    matched = 0
    total_th_only = 0
    failures = []

    with SDFRead(str(sdf_path)) as reader:
        for mol in reader:
            units = mol.stereo_units()
            has_other_kind = any(u['kind'] != SU_TETRA and u['parity'] != 0 for u in units)
            if has_other_kind:
                continue  # cis/trans, allene, atropisomer: not this test's subject
            has_h_unknown = any(a.implicit_h is None for a in mol.atoms())
            if has_h_unknown:
                continue  # H_UNKNOWN loss will change InChI; skip to keep test clean

            total_th_only += 1
            inchi_orig = _inchi(mol)
            mol_rt = from_indigo(to_indigo(mol))
            inchi_rt = _inchi(mol_rt)
            if inchi_orig == inchi_rt:
                matched += 1
            else:
                failures.append((inchi_orig, inchi_rt))

    assert matched == total_th_only, (
        f'Only {matched}/{total_th_only} tetrahedral-only records round-tripped cleanly.\n'
        f'First failure:\n  orig: {failures[0][0]}\n  rt:   {failures[0][1]}' if failures else ''
    )
    # The corpus must contain a meaningful number of tetrahedral-only records.
    assert total_th_only >= 150, (
        f'Expected at least 150 tetrahedral-only records; found only {total_th_only}'
    )


@requires_indigo
def test_allene_is_a_logged_loss():
    """
    A configured allene loses its configuration through Indigo, and `to_indigo` says so.

    Indigo has no allene representation at all; InChI's `/t` layer is what makes the drop visible.
    Penta-2,3-diene, with the parity set directly: the SMILES reader does not configure allene units.
    """
    from chython.core import SU_ALLENE
    from chython.core._core import read_smiles
    from chython.interop._indigo import to_indigo, from_indigo

    mol = read_smiles('C/C=C=C/C')
    allene = next(u for u in mol.stereo_units() if u['kind'] == SU_ALLENE)
    assert allene['stereogenic'] and allene['parity'] == 0
    mol.set_parity(allene['anchor'], 1)

    log = []
    ig_mol = to_indigo(mol, log=log)
    assert any('allene' in entry for entry in log), log

    inchi_orig = _inchi(mol)
    assert '/t' in inchi_orig, inchi_orig  # InChI does hold the configuration
    assert '/t' not in _inchi(from_indigo(ig_mol))  # ... and Indigo cannot carry it


@requires_indigo
def test_the_index_reaches_the_indigo_object():
    """Indigo carries a pseudo-atom label verbatim, so the index is not a fact it cannot hold."""
    from chython.core._core import read_smiles
    from chython.interop._indigo import to_indigo

    mol = read_smiles('[R7]C')
    assert 'R7' in to_indigo(mol).smiles()


@requires_indigo
def test_the_export_says_the_marker_does_not_come_back():
    """``to_indigo`` logs, once per molecule, that the marker does not survive the round trip.

    The write direction is the whole assertion: one ``LOST`` record naming the marker, whatever the
    molecule's R count.
    """
    from chython.core._core import read_smiles
    from chython.interop._indigo import to_indigo

    mol = read_smiles('[R7]C')
    log = []
    to_indigo(mol, log=log)
    assert any('marker' in str(x) for x in log), log
