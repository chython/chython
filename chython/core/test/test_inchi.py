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
"""Tests for the native InChI binding (_inchi.pxi).

All tests are skipped when libinchi is not loaded (e.g. CI without the binary),
or when the InChI functions are not compiled into this build (older .so).
"""
import pytest
from chython.core import MoleculeContainer, read_smiles

try:
    from chython.core import (
        inchi,
        inchi_library_loaded,
        inchikey,
        molecule_to_inchi,
        molecule_to_inchikey,
        inchi_to_molecule,
    )
    _INCHI_AVAILABLE = inchi_library_loaded()
except ImportError:
    # InChI functions are not compiled into this build (.so predates _inchi.pxi).
    _INCHI_AVAILABLE = False

    def inchi_library_loaded():
        return False

    molecule_to_inchi = molecule_to_inchikey = inchi_to_molecule = None
    inchi = inchikey = None

pytestmark = pytest.mark.skipif(
    not _INCHI_AVAILABLE,
    reason='libinchi not loaded or InChI functions not compiled — skipping InChI tests'
)


# ---- helpers ---------------------------------------------------------------- #

def _build_ethanol():
    """Return ethanol (CCO) as a MoleculeContainer.  Atoms: C(1), C(2), O(3)."""
    mol = MoleculeContainer()
    with mol:
        c1 = mol.add_atom(6, implicit_h=3)   # CH3
        c2 = mol.add_atom(6, implicit_h=2)   # CH2
        o  = mol.add_atom(8, implicit_h=1)   # OH
        mol.add_bond(c1, c2, 1)
        mol.add_bond(c2, o, 1)
    return mol, c1, c2, o


def _build_benzene():
    """Return benzene (6 carbons in a Kekulé ring) as a MoleculeContainer."""
    mol = MoleculeContainer()
    with mol:
        sids = [mol.add_atom(6, implicit_h=1) for _ in range(6)]
        for i in range(6):
            mol.add_bond(sids[i], sids[(i + 1) % 6], 2 if i % 2 == 0 else 1)
    return mol


def _build_alanine():
    """Build alanine NH2-CH(CH3)-COOH (6 heavy atoms).

    Returns (mol, ca_sid) where ca_sid is the alpha carbon stable id.
    """
    mol = MoleculeContainer()
    with mol:
        n   = mol.add_atom(7, implicit_h=2)
        ca  = mol.add_atom(6, implicit_h=1)
        cme = mol.add_atom(6, implicit_h=3)
        cc  = mol.add_atom(6)
        ok  = mol.add_atom(8)
        oh  = mol.add_atom(8, implicit_h=1)
        mol.add_bond(n,  ca,  1)
        mol.add_bond(ca, cme, 1)
        mol.add_bond(ca, cc,  1)
        mol.add_bond(cc, ok,  2)
        mol.add_bond(cc, oh,  1)
    return mol, ca


# ---- library load test ------------------------------------------------------ #

def test_lib_loaded():
    """inchi_library_loaded() returns True (the skipif guard ensures this)."""
    assert inchi_library_loaded()


# ---- forward direction ------------------------------------------------------- #

def test_ethanol_forward():
    mol, *_ = _build_ethanol()
    inchi = molecule_to_inchi(mol)
    assert inchi.startswith('InChI=1S/')
    assert 'C2H6O' in inchi


def test_benzene_forward():
    mol = _build_benzene()
    inchi = molecule_to_inchi(mol)
    assert inchi.startswith('InChI=1S/')
    assert 'C6H6' in inchi


def test_ethanol_forward_known_inchi():
    """Ethanol InChI must equal the IUPAC standard value."""
    mol, *_ = _build_ethanol()
    assert molecule_to_inchi(mol) == 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3'


def test_nonstd_inchi_option():
    """GetINCHI (non-standard call) with options should still return an InChI string.

    Note: libinchi may return 'InChI=1S/' even when GetINCHI is called if the
    molecule happens to satisfy standard-InChI constraints.  We only verify that
    the function completes without error and returns something starting with 'InChI='.
    """
    mol, *_ = _build_ethanol()
    inchi_opt = molecule_to_inchi(mol, standard=False, options='-SNon')
    assert inchi_opt.startswith('InChI=')


def test_empty_molecule_raises():
    """molecule_to_inchi on an empty molecule should raise (no atoms)."""
    mol = MoleculeContainer()
    with pytest.raises(ValueError):
        molecule_to_inchi(mol)


# ---- InChIKey --------------------------------------------------------------- #

def test_inchikey_ethanol():
    mol, *_ = _build_ethanol()
    key = molecule_to_inchikey(mol)
    # InChIKey is always 27 chars: XXXXXXXXXXXXXX-XXXXXXXXXX-N
    assert len(key) == 27
    assert key.count('-') == 2
    assert key == 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N'


def test_inchikey_benzene():
    mol = _build_benzene()
    assert molecule_to_inchikey(mol) == 'UHOVQNZJYSORNB-UHFFFAOYSA-N'


# ---- reverse direction ------------------------------------------------------- #

def test_ethanol_reverse():
    inchi = 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3'
    mol = inchi_to_molecule(inchi)
    # Should have 3 heavy atoms
    sids = list(mol.atoms())
    assert len(sids) == 3


def test_ethanol_round_trip():
    """molecule → InChI → molecule → InChI must produce identical strings."""
    mol, *_ = _build_ethanol()
    inchi1 = molecule_to_inchi(mol)
    mol2 = inchi_to_molecule(inchi1)
    inchi2 = molecule_to_inchi(mol2)
    assert inchi1 == inchi2


def test_reverse_bad_inchi_raises():
    with pytest.raises((ValueError, Exception)):
        inchi_to_molecule('InChI=1S/utter/garbage')


# ---- stereo round-trip: alanine enantiomers ---------------------------------- #
#
# Correct alanine InChI strings (generated by molecule_to_inchi, verified correct):
# parity=1 → /t2-/m0/s1  (R-alanine)
# parity=2 → /t2-/m1/s1  (S-alanine)
#
# R-alanine InChIKey: QNAYBMKLOCPYGJ-REOHCLBHSA-N
#
_ALANINE_NOSTEREO = 'InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)'
_ALANINE_R_INCHI  = 'InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1'
_ALANINE_S_INCHI  = 'InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m1/s1'
_ALANINE_R_KEY    = 'QNAYBMKLOCPYGJ-REOHCLBHSA-N'


def test_forward_stereo_alanine():
    """Setting parity 1 vs 2 must yield two different InChI strings."""
    mol, ca = _build_alanine()
    mol.set_parity(ca, 1)
    inchi1 = molecule_to_inchi(mol)
    mol.set_parity(ca, 2)
    inchi2 = molecule_to_inchi(mol)
    assert inchi1 != inchi2
    # One must be R, other S
    assert {inchi1, inchi2} == {_ALANINE_R_INCHI, _ALANINE_S_INCHI}


def test_forward_stereo_r_alanine_key():
    """R-alanine (parity=1) must produce the canonical CAS InChIKey."""
    mol, ca = _build_alanine()
    mol.set_parity(ca, 1)
    key = molecule_to_inchikey(mol)
    assert key == _ALANINE_R_KEY


def test_reverse_stereo_r_alanine():
    """Parse R-alanine InChI and verify round-trip InChI matches."""
    mol = inchi_to_molecule(_ALANINE_R_INCHI)
    rt = molecule_to_inchi(mol)
    assert rt == _ALANINE_R_INCHI


def test_reverse_stereo_s_alanine():
    """Parse S-alanine InChI and verify round-trip InChI matches."""
    mol = inchi_to_molecule(_ALANINE_S_INCHI)
    rt = molecule_to_inchi(mol)
    assert rt == _ALANINE_S_INCHI


def test_reverse_stereo_parity_set():
    """After parsing a stereo InChI, at least one stereo unit should have parity != 0."""
    mol = inchi_to_molecule(_ALANINE_R_INCHI)
    units = mol.stereo_units()
    configured = [u for u in units if u['parity'] != 0]
    assert configured, 'no stereo unit with parity set after parsing stereo InChI'


def test_reverse_enantiomers_differ():
    """R and S alanine parsed from InChI must give different InChIKeys."""
    mol_r = inchi_to_molecule(_ALANINE_R_INCHI)
    mol_s = inchi_to_molecule(_ALANINE_S_INCHI)
    assert molecule_to_inchikey(mol_r) != molecule_to_inchikey(mol_s)


def test_forward_stereo_nostereo_no_t_layer():
    """Molecule without stereo parity set should produce InChI with no /t layer."""
    mol, ca = _build_alanine()
    # do not set parity: mol.set_parity(ca, 0) is default
    inchi = molecule_to_inchi(mol)
    assert '/t' not in inchi


# ---- aromatic-input invariants ---------------------------------------------- #
#
# Input fidelity: molecule_to_inchi must NEVER mutate the caller's molecule.
# This is true even when an internal copy is made for kekulisation.
#
# Aromatic-bond molecules (HE_AROMATIC stored / order 4) require a kekuliser that
# is not present in this build yet (the arom module registers it via
# _ich_set_kekule_fn).  Until then, molecule_to_inchi raises ValueError if given
# a molecule with aromatic bonds.
#
# The tests below cover:
#   1. Kekule benzene: correct InChIKey, bond orders unchanged after the call.
#   2. No aromatic kekuliser registered yet: molecule_to_inchi raises ValueError
#      for a molecule with an order-4 bond in the journal (before it is committed).
# Once the arom module lands and registers its kekuliser, the aromatic-storage
# tests in test_inchi_aromatic.py will be enabled.

def test_kekule_benzene_inchikey_and_bond_orders_preserved():
    """Kekule benzene: InChIKey is correct and bond orders are unchanged after the call.

    Bond order preservation is the input-fidelity invariant: molecule_to_inchi must
    never convert aromatic bonds to Kekule (or vice versa) in the caller's molecule.
    """
    mol = _build_benzene()

    # Collect bond orders before the call (Bond.n, Bond.m, Bond.order).
    bonds_before = {(b.n, b.m): b.order for b in mol.bonds()}

    key = molecule_to_inchikey(mol)
    assert key == 'UHOVQNZJYSORNB-UHFFFAOYSA-N', f'unexpected InChIKey: {key}'

    # Bond orders must be exactly the same after the call.
    bonds_after = {(b.n, b.m): b.order for b in mol.bonds()}
    assert bonds_before == bonds_after, 'molecule_to_inchi mutated the input bond orders'


def test_forward_does_not_mutate_charge_or_isotope():
    """Forward call must not change atom properties (charge, isotope) on the input."""
    mol, c1, c2, o = _build_ethanol()

    def _snapshot(m):
        return {a.n: (a.charge, a.isotope) for a in m.atoms()}

    before = _snapshot(mol)
    molecule_to_inchi(mol)
    after = _snapshot(mol)
    assert before == after, 'molecule_to_inchi mutated atom properties'


def test_pyridine_inchikey():
    """Pyridine (Kekule form) must produce the correct InChIKey."""
    mol = MoleculeContainer()
    with mol:
        n  = mol.add_atom(7, implicit_h=0)
        c2 = mol.add_atom(6, implicit_h=1)
        c3 = mol.add_atom(6, implicit_h=1)
        c4 = mol.add_atom(6, implicit_h=1)
        c5 = mol.add_atom(6, implicit_h=1)
        c6 = mol.add_atom(6, implicit_h=1)
        # Kekule: N=C-C=C-C=C (ring)
        mol.add_bond(n,  c2, 2)
        mol.add_bond(c2, c3, 1)
        mol.add_bond(c3, c4, 2)
        mol.add_bond(c4, c5, 1)
        mol.add_bond(c5, c6, 2)
        mol.add_bond(c6, n,  1)
    key = molecule_to_inchikey(mol)
    assert key == 'JUJWROOIHBZHMG-UHFFFAOYSA-N', f'unexpected pyridine InChIKey: {key}'


def test_naphthalene_inchikey():
    """Naphthalene (Kekule, fused ring) must produce the correct InChIKey."""
    # Standard InChIKey for naphthalene: UFWIBTONFRDIAS-UHFFFAOYSA-N
    inchi_naphthalene = 'InChI=1S/C10H8/c1-2-6-10-8-4-3-7-9(10)5-1/h1-8H'
    mol = inchi_to_molecule(inchi_naphthalene)
    key = molecule_to_inchikey(mol)
    assert key == 'UFWIBTONFRDIAS-UHFFFAOYSA-N', f'unexpected naphthalene InChIKey: {key}'


# ---- bond-kind stereo: cis/trans and allene --------------------------------- #
#
# Both kinds were DEAD before this suite existed: all four `translate_stereo` call sites in
# `_inchi.pxi` passed InChI's own `(X, A, B, Y)` neighbour tuple, whose middle two entries are the
# CHAIN atoms and not the unit's refs at all, so every call raised
# `ValueError: order is not a permutation of the unit refs`.  The tests below cover both kinds in
# both directions and, crucially, pin the two SIGNS -- the round trip alone cannot, because it is
# insensitive to a global flip by construction (`_ich_inchi_to_chython` is the exact inverse of
# `_ich_chython_to_inchi` for any value of the constant).
#
# THE ABSOLUTE ASSERTIONS are the ones that pin the signs.  They were measured on 2026-09-02
# against libinchi built from the bundled INCHI submodule (11a8798), which perceives both kinds
# from coordinates:
#
#   * cis/trans -- a 2D but-2-ene with the two methyls on OPPOSITE sides perceives as `/b4-3+`,
#     and a zero-coordinate 0D record of EVEN reproduces exactly that string.  `/b4-3+` is also
#     the published standard InChI of (E)-but-2-ene, so this assertion has an authority outside
#     this repo.
#   * allene -- the two mirror conformers of 1,3-dibromo-1,3-difluoroallene perceive as
#     `/t1-/m0/s1` and `/t1-/m1/s1`, and the `/m0` one is the conformer whose signed volume over
#     the unit's OWN ref order `(a1, a2, b1, b2)` is POSITIVE, which is the core's even (parity 1)
#     under the same signed-volume rule the tetrahedral convention uses.
#
# Flipping either constant in `_inchi.pxi` flips the corresponding absolute assertions and leaves
# every round-trip test green -- which is the whole reason the absolute ones are here.

_E_BUT_2_ENE = 'InChI=1S/C4H8/c1-3-4-2/h3-4H,1-2H3/b4-3+'   # methyls TRANS
_Z_BUT_2_ENE = 'InChI=1S/C4H8/c1-3-4-2/h3-4H,1-2H3/b4-3-'   # methyls CIS
_ALLENE_M0 = 'InChI=1S/C3Br2F2/c4-2(6)1-3(5)7/t1-/m0/s1'
_ALLENE_M1 = 'InChI=1S/C3Br2F2/c4-2(6)1-3(5)7/t1-/m1/s1'


def _build_but_2_ene(parity=0):
    """but-2-ene CH3-CH=CH-CH3.  Returns (mol, anchor_sid), anchor on the first sp2 carbon.

    refs of the SU_CIS_TRANS unit are (C1-methyl, None | C4-methyl, None): one named heavy
    substituent and one implicit hydrogen per terminal, so refs[0]/refs[2] are the methyls.
    """
    mol = MoleculeContainer()
    with mol:
        c1 = mol.add_atom(6, implicit_h=3)
        c2 = mol.add_atom(6, implicit_h=1)
        c3 = mol.add_atom(6, implicit_h=1)
        c4 = mol.add_atom(6, implicit_h=3)
        mol.add_bond(c1, c2, 1)
        mol.add_bond(c2, c3, 2)
        mol.add_bond(c3, c4, 1)
    if parity:
        mol.set_parity(c2, parity)
    return mol, c2


def _build_allene(parity=0):
    """1,3-dibromo-1,3-difluoroallene FBrC=C=CBrF.  Returns (mol, centre_sid).

    Tetrasubstituted, so both terminals carry two NAMED substituents and the unit's refs hold no
    unnamed direction -- the cleanest possible axial frame.
    """
    mol = MoleculeContainer()
    with mol:
        ca = mol.add_atom(6)
        cb = mol.add_atom(6)
        cc = mol.add_atom(6)
        f1 = mol.add_atom(9)
        br1 = mol.add_atom(35)
        f2 = mol.add_atom(9)
        br2 = mol.add_atom(35)
        mol.add_bond(ca, cb, 2)
        mol.add_bond(cb, cc, 2)
        mol.add_bond(ca, f1, 1)
        mol.add_bond(ca, br1, 1)
        mol.add_bond(cc, f2, 1)
        mol.add_bond(cc, br2, 1)
    if parity:
        mol.set_parity(cb, parity)
    return mol, cb


def _unit_of(mol, anchor_sid):
    for u in mol.stereo_units():
        if u['anchor'] == anchor_sid:
            return u
    raise AssertionError(f'no stereo unit anchored at {anchor_sid}')


# ---- the frame itself ------------------------------------------------------- #

def test_cis_trans_unit_is_perceived_with_one_named_atom_per_pair():
    """The fixture really is a bond-kind unit whose refs[0]/refs[2] are the two methyls."""
    mol, anchor = _build_but_2_ene()
    u = _unit_of(mol, anchor)
    assert u['n_refs'] == 4
    assert u['refs'][0] is not None and u['refs'][2] is not None, 'F26: pair slot 0 is named'
    # the two unnamed directions are the implicit hydrogens, one per terminal
    assert u['refs'][1] is None and u['refs'][3] is None


def test_allene_unit_is_perceived_with_four_named_refs():
    mol, centre = _build_allene()
    u = _unit_of(mol, centre)
    assert u['n_refs'] == 4
    assert all(r is not None for r in u['refs']), 'tetrasubstituted allene names all four'
    assert u['unnamed_mask'] == 0


def test_export_passes_zero_coordinates_so_the_0d_records_are_what_count():
    """inchi_api.h: 0D parities are honoured only when all atom coordinates are zero.

    `_ich_fill_atoms` memsets every ICH_Atom and never writes x/y/z, so the condition holds
    structurally -- but the observable proof is that the two parities produce two DIFFERENT
    strings.  Were libinchi perceiving from geometry instead, the 0D records would be ignored
    and both parities would collapse to the same output.
    """
    for build in (_build_but_2_ene, _build_allene):
        mol, _ = build()
        assert not mol.has_coordinates, 'fixture must carry no coordinates'
        even, _ = build(1)
        odd, _ = build(2)
        assert molecule_to_inchi(even) != molecule_to_inchi(odd), (
            f'{build.__name__}: the 0D stereo record had no effect on the output')


# ---- cis/trans, export ------------------------------------------------------ #

def test_cis_trans_export_does_not_raise():
    """The original defect: this call raised ValueError out of translate_stereo."""
    mol, _ = _build_but_2_ene(1)
    assert '/b' in molecule_to_inchi(mol)


def test_cis_trans_export_gives_the_two_spellings():
    even, _ = _build_but_2_ene(1)
    odd, _ = _build_but_2_ene(2)
    assert molecule_to_inchi(even) != molecule_to_inchi(odd)


def test_cis_trans_even_is_trans():
    """ABSOLUTE, and the assertion ICH_CIS_TRANS_FLIP is answerable to.

    Core parity 1 (even) means refs[0] and refs[2] -- here the two methyls -- are TRANS, so the
    export must be the published (E)-but-2-ene string.
    """
    even, _ = _build_but_2_ene(1)
    assert molecule_to_inchi(even) == _E_BUT_2_ENE
    odd, _ = _build_but_2_ene(2)
    assert molecule_to_inchi(odd) == _Z_BUT_2_ENE


def test_a_cis_trans_units_refs_0_is_bonded_to_its_anchor_and_refs_2_is_not():
    """WHY THE CIS/TRANS EXPORT ORIENTS NOTHING while the allene export must.

    Both records need `neighbor[0]` (X) bonded to `neighbor[1]` (A).  The allene anchors on the
    CENTRE, so `refs[0:2]` is whichever terminal's pair perception happened to fill first and the
    export has to orient A/B by adjacency to X -- that is
    `test_allene_terminals_are_oriented_by_the_refs_not_by_edge_order`.  The cis/trans export just
    writes `neighbor[0] = refs[0]`, `neighbor[1] = anchor` with no adjacency check, and this is the
    assumption that makes that safe: a cis/trans unit ANCHORS ON A TERMINAL, and it is the same
    terminal `refs[0:2]` was filled from (`_stereo.pxi` pass 2 walks from the lower-indexed
    terminal, anchors there, and fills `refs[0]` from `_terminal_pair` on it).

    So the asymmetry between the two branches is not an oversight in one of them, and the ONLY
    thing holding it up is a same-source guarantee in another file -- which is what this asserts
    directly, in both bond orders and with both one and two named substituents per terminal, rather
    than leaving it to be inferred from a green round trip.
    """
    def dichlorobutene(reverse):
        mol = MoleculeContainer()
        with mol:
            c1 = mol.add_atom(6, implicit_h=3)
            c2 = mol.add_atom(6)
            c3 = mol.add_atom(6)
            c4 = mol.add_atom(6, implicit_h=3)
            l2 = mol.add_atom(17)
            l3 = mol.add_atom(17)
            if reverse:
                mol.add_bond(c3, c4, 1)
                mol.add_bond(c2, c3, 2)
                mol.add_bond(c1, c2, 1)
            else:
                mol.add_bond(c1, c2, 1)
                mol.add_bond(c2, c3, 2)
                mol.add_bond(c3, c4, 1)
            mol.add_bond(c2, l2, 1)
            mol.add_bond(c3, l3, 1)
        return mol, c2

    cases = [('but-2-ene', _build_but_2_ene()),
             ('2,3-dichlorobut-2-ene', dichlorobutene(False)),
             ('2,3-dichlorobut-2-ene, chain reversed', dichlorobutene(True))]
    for label, (mol, anchor) in cases:
        u = _unit_of(mol, anchor)
        assert u['kind'] == 1, f'{label}: fixture is not a cis/trans unit'
        near = set(mol.neighbors_of(anchor))
        assert u['refs'][0] in near, f'{label}: refs[0] is not on the anchor terminal'
        assert u['refs'][2] not in near, f'{label}: refs[2] is on the anchor terminal too'
        if u['refs'][1] is not None:
            assert u['refs'][1] in near, f'{label}: refs[1] left the anchor pair'
        if u['refs'][3] is not None:
            assert u['refs'][3] not in near, f'{label}: refs[3] joined the anchor pair'


# ---- cis/trans, import ------------------------------------------------------ #

def test_cis_trans_import_does_not_raise():
    """The mirror defect on the reverse path; this raised ValueError too."""
    assert inchi_to_molecule(_E_BUT_2_ENE) is not None


def test_cis_trans_import_sets_a_parity():
    for s in (_E_BUT_2_ENE, _Z_BUT_2_ENE):
        mol = inchi_to_molecule(s)
        assert any(mol.parity_of(sid) != 0 for sid in mol.atom_numbers), f'no parity configured importing {s}'


def test_cis_trans_import_distinguishes_e_from_z():
    e = inchi_to_molecule(_E_BUT_2_ENE)
    z = inchi_to_molecule(_Z_BUT_2_ENE)
    assert molecule_to_inchi(e) != molecule_to_inchi(z)


def test_cis_trans_round_trip_preserves_configuration():
    for s in (_E_BUT_2_ENE, _Z_BUT_2_ENE):
        assert molecule_to_inchi(inchi_to_molecule(s)) == s


def test_cis_trans_round_trip_from_the_core_side():
    """Export, import, re-export: the configuration must survive both crossings."""
    for parity in (1, 2):
        mol, _ = _build_but_2_ene(parity)
        first = molecule_to_inchi(mol)
        assert molecule_to_inchi(inchi_to_molecule(first)) == first


def test_cis_trans_round_trip_with_two_named_substituents_per_terminal():
    """2,3-dichlorobut-2-ene: every terminal has TWO named substituents.

    InChI picks one of them as its X (or Y) and it need not be the one at the pair's slot 0.  When
    it picks the other, that is a within-pair transposition and the parity flips --
    `_ich_bond_order_from_refs` locates X and Y inside `refs` rather than assuming slot 0, so
    `translate_stereo` charges the flip correctly.
    """
    def build(parity):
        mol = MoleculeContainer()
        with mol:
            c1 = mol.add_atom(6, implicit_h=3)
            c2 = mol.add_atom(6)
            c3 = mol.add_atom(6)
            c4 = mol.add_atom(6, implicit_h=3)
            l2 = mol.add_atom(17)
            l3 = mol.add_atom(17)
            mol.add_bond(c1, c2, 1)
            mol.add_bond(c2, c3, 2)
            mol.add_bond(c3, c4, 1)
            mol.add_bond(c2, l2, 1)
            mol.add_bond(c3, l3, 1)
        mol.set_parity(c2, parity)
        return mol

    seen = set()
    for parity in (1, 2):
        s = molecule_to_inchi(build(parity))
        assert molecule_to_inchi(inchi_to_molecule(s)) == s, f'parity {parity} lost on round trip'
        seen.add(s)
    assert len(seen) == 2, 'the two parities must be two different strings'


# ---- allene, export -------------------------------------------------------- #

def test_allene_export_does_not_raise():
    """The original defect on the axial path."""
    mol, _ = _build_allene(1)
    assert '/t' in molecule_to_inchi(mol)


def test_allene_export_gives_the_two_spellings():
    even, _ = _build_allene(1)
    odd, _ = _build_allene(2)
    assert molecule_to_inchi(even) != molecule_to_inchi(odd)


def test_allene_even_is_the_positive_volume_conformer():
    """ABSOLUTE, and the assertion ICH_ALLENE_FLIP is answerable to.

    libinchi perceives `/t1-/m0/s1` for the conformer whose signed volume over the unit's own ref
    order (a1, a2, b1, b2) is POSITIVE, and the core calls that handedness even (parity 1).  So
    even must export as `/m0`.
    """
    even, _ = _build_allene(1)
    assert molecule_to_inchi(even) == _ALLENE_M0
    odd, _ = _build_allene(2)
    assert molecule_to_inchi(odd) == _ALLENE_M1


def test_allene_terminals_are_oriented_by_the_refs_not_by_edge_order():
    """`_ich_find_allene_terminals` returns terminals in CSR edge order, which is unrelated to
    which pair perception put first -- the anchor is the CENTRE, so refs[0:2] is simply one
    terminal's pair, not "the anchor's".  InChI requires neighbor[0] (X) to be bonded to
    neighbor[1] (A), so the export orients A/B by adjacency to X.  Adding the centre's two chain
    bonds in the opposite order must therefore give the SAME InChI for the same parity.
    """
    def build(parity, reverse):
        mol = MoleculeContainer()
        with mol:
            ca = mol.add_atom(6)
            cb = mol.add_atom(6)
            cc = mol.add_atom(6)
            f1 = mol.add_atom(9)
            br1 = mol.add_atom(35)
            f2 = mol.add_atom(9)
            br2 = mol.add_atom(35)
            if reverse:
                mol.add_bond(cb, cc, 2)
                mol.add_bond(ca, cb, 2)
            else:
                mol.add_bond(ca, cb, 2)
                mol.add_bond(cb, cc, 2)
            mol.add_bond(ca, f1, 1)
            mol.add_bond(ca, br1, 1)
            mol.add_bond(cc, f2, 1)
            mol.add_bond(cc, br2, 1)
        mol.set_parity(cb, parity)
        return mol

    for parity in (1, 2):
        a = molecule_to_inchi(build(parity, False))
        b = molecule_to_inchi(build(parity, True))
        assert a == b, f'parity {parity}: chain-bond order changed the axial sign ({a} vs {b})'


# ---- allene, import -------------------------------------------------------- #

def test_allene_import_does_not_raise():
    assert inchi_to_molecule(_ALLENE_M0) is not None


def test_allene_import_sets_a_parity_on_the_centre():
    mol = inchi_to_molecule(_ALLENE_M0)
    assert any(mol.parity_of(sid) != 0 for sid in mol.atom_numbers), 'no parity configured on allene import'


def test_allene_import_distinguishes_the_enantiomers():
    m0 = inchi_to_molecule(_ALLENE_M0)
    m1 = inchi_to_molecule(_ALLENE_M1)
    assert molecule_to_inchi(m0) != molecule_to_inchi(m1)


def test_allene_round_trip_preserves_configuration():
    for s in (_ALLENE_M0, _ALLENE_M1):
        assert molecule_to_inchi(inchi_to_molecule(s)) == s


def test_allene_round_trip_from_the_core_side():
    for parity in (1, 2):
        mol, _ = _build_allene(parity)
        first = molecule_to_inchi(mol)
        assert molecule_to_inchi(inchi_to_molecule(first)) == first


# ---- a stereogenic bond whose parity is explicitly undefined ---------------- #

# libinchi 1.07.5 renders an undefined bond configuration as "?" in the /b layer.  That state is
# NOT the same as "no stereo here", and the container does distinguish the two: the unit is
# perceived with `stereogenic=True` while `parity` stays 0.  `translate_stereo` returns 0 for both
# an undefined and an unset parity, so the distinction lives in the STEREOGENICITY MARK, not in
# the parity value.
#
# One asymmetry is worth recording, because it looks like a fidelity bug and is not ours.  Feeding
# "b3-1?,5-4+" back out drops the "?".  The mechanism, measured rather than assumed:
#
#   * a double-bond terminal left as a bare FREE VALENCE (carbon, one substituent, one H, radical
#     flag clear) is accepted by libinchi 1.07.5's half_stereo_bond_parity() and printed as "?";
#   * the SAME terminal carrying a doublet RADICAL is rejected outright by
#     bCanAtomHaveAStereoBond(), so no "?" is printed;
#   * `inchi_to_molecule` reads InChI's under-valent carbon as a radical -- which is the correct
#     chemical reading of that species -- so re-export takes the second path.
#
# The two forms are indistinguishable in the InChI string, so giving them different /b layers is a
# libinchi defect (upstream github #263, fixed after 1.07.5 by rejecting the free valence too).
# Our output is already what the corrected library produces, so there is nothing to repair here and
# nothing to represent: the "?" is absent because the bond is genuinely not stereogenic on a
# radical, and the container still holds the stereogenic-but-undefined state either way.
#
# The assertions below therefore pin OUR behaviour only.  They deliberately do NOT assert that
# libinchi prints "?" for the free-valence form, since that is exactly what upstream is changing.

_UNDEF_BOND = 'InChI=1S/C5H7/c1-3-5-4-2/h1,3-5H,2H3/b3-1?,5-4+'


def test_undefined_bond_parity_imports_as_stereogenic_but_unset():
    """"?" is representable: the unit exists, is stereogenic, and carries no parity."""
    mol = inchi_to_molecule(_UNDEF_BOND)
    bond_units = [u for u in mol.stereo_units() if u['kind'] == 1]
    assert len(bond_units) == 2, 'both double bonds are perceived as bond-kind units'
    undefined = [u for u in bond_units if u['parity'] == 0]
    assert len(undefined) == 1
    assert undefined[0]['stereogenic'], 'stereogenic flag must survive an undefined parity'


def test_undefined_bond_parity_does_not_swallow_the_defined_bond():
    """The "?" on one bond must not disturb the sign of the other."""
    assert molecule_to_inchi(inchi_to_molecule(_UNDEF_BOND)).endswith('/b5-4+')
    assert molecule_to_inchi(inchi_to_molecule(_UNDEF_BOND.replace('5-4+', '5-4-'))).endswith('/b5-4-')


def test_undefined_bond_parity_round_trip_reaches_a_fixed_point():
    """Re-export is stable: one pass normalises, further passes change nothing."""
    first = molecule_to_inchi(inchi_to_molecule(_UNDEF_BOND))
    assert molecule_to_inchi(inchi_to_molecule(first)) == first


class TestInchiFacade:
    """`inchi()` in both directions; `inchikey()` one way."""

    def test_a_molecule_becomes_an_inchi_string(self):
        assert inchi(read_smiles('CCO')).startswith('InChI=')

    def test_an_inchi_string_becomes_a_molecule(self):
        assert str(inchi('InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3')) == 'C(C)O'

    def test_a_round_trip_is_stable(self):
        text = inchi(read_smiles('c1ccccc1'))
        assert inchi(inchi(text)) == text

    def test_the_options_reach_the_writer(self):
        """`-SNon` drops the stereo layer, which is the observable proof the flags are passed on.

        `standard=` reaches the library too, but for most structures `GetINCHI` with no options
        answers what `GetStdINCHI` does, so it is not what an assertion can be built on.
        """
        m = read_smiles('C[C@H](N)O')
        assert '/t' in inchi(m)
        assert '/t' not in inchi(m, options='-SNon')
        assert '/t' not in inchi(m, standard=False, options='-SNon')

    def test_the_option_prefix_is_the_one_this_platform_reads(self):
        """What makes the assertion above hold off Linux, checked where the rewrite happens.

        `inchi_api.h` prefixes szOptions with `/` on Windows and `-` elsewhere, and libinchi IGNORES an
        option it does not recognise -- so `-SNon` there returns a stereo layer and nothing says why.
        Both spellings go in and the platform's comes out, which is a statement every host can check.
        """
        from sys import platform

        from chython.core._core import _ich_platform_options

        prefix = '/' if platform == 'win32' else '-'
        assert _ich_platform_options('-SNon') == prefix + 'SNon'
        assert _ich_platform_options('/SNon') == prefix + 'SNon'
        assert _ich_platform_options('-SNon -DoNotAddH') == f'{prefix}SNon {prefix}DoNotAddH'

    def test_a_string_that_is_not_an_inchi_is_refused_by_name(self):
        with pytest.raises(ValueError, match='InChI='):
            inchi('CCO')

    def test_an_inchikey_is_its_own_name_because_it_cannot_be_read_back(self):
        key = inchikey(read_smiles('CCO'))
        assert len(key) == 27 and key.count('-') == 2
        with pytest.raises(ValueError, match='InChI='):
            inchi(key)

    def test_anything_that_is_neither_is_refused_by_type(self):
        with pytest.raises(TypeError, match='inchi'):
            inchi(42)
