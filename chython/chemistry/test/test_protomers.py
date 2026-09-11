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
"""`neutralize()`: moving a proton from a tabulated cation onto a tabulated anion.

The invariant under test is that NOTHING OVERSHOOTS ZERO -- a component is only ever taken closer to
neutral -- which is why nitrate takes one proton and sulfate two, and why a lone cation only moves
under `keep_charge=False`.  The refusals matter as much as the moves: what stays charged is a decision
the pass has to state.
"""
# `__all__` and not the package object: `from ... import chemistry` would execute the facade, which
# `test_dependency_direction.py` ratchets against.
from .. import __all__ as CHEMISTRY_ALL, canonicalize, implicify_hydrogens, neutralize
from ...core import INFO, REFUSED, read_smiles as smiles


def _charge(molecule):
    return sum(molecule.charge_of(n) for n in molecule)


# what moves


def test_a_zwitterion_becomes_its_neutral_form():
    """One proton crosses the molecule and both ends go neutral -- unless the cation has none to give."""
    for s, expect in [
            ('[NH3+]CC(=O)[O-]', 'C(CN)(=O)O'),                       # glycine
            ('[NH3+]CCS(=O)(=O)[O-]', 'S(CCN)(O)(=O)=O'),             # taurine
            ('C[N+](C)(C)CCC([O-])=O', 'O=C([O-])CC[N+](C)(C)C')]:    # a betaine: untouched
        m = smiles(s)
        neutralize(m)
        assert format(m) == expect, s


def test_a_salt_becomes_the_free_acid_and_the_free_base():
    for s, expect in [
            ('C[NH3+].[Cl-]', 'CN.Cl'),
            ('CC(=O)[O-].[NH4+]', 'C(C)(=O)O.N'),
            ('[NH4+].[OH-]', 'N.O'),
            ('c1cc[nH+]cc1.CC(=O)[O-]', 'c1ccccn1.C(C)(=O)O')]:
        m = smiles(s)
        assert neutralize(m) is True, s
        assert format(m) == expect, s


def test_an_aromatic_cation_is_deprotonatable():
    """An aromatic bond has no valence row, so the post-move question cannot be put and must not be
    read as a violation: a pyridinium that could not give its proton up is the failure this pins."""
    m = smiles('c1cc[nH+]cc1.[Cl-]')
    assert neutralize(m) is True
    assert format(m) == 'c1ccccn1.Cl'


def test_a_stereocentre_survives():
    """Only a charge and an implicit count are written, so alanine keeps its configuration."""
    zwitterion, neutral = smiles('C[C@H]([NH3+])C(=O)[O-]'), smiles('C[C@H](N)C(=O)O')
    neutralize(zwitterion)
    canonicalize(zwitterion)
    canonicalize(neutral)
    assert zwitterion == neutral
    assert zwitterion != smiles('C[C@@H](N)C(=O)O')


# the overshoot guard


def test_nitrate_takes_one_proton_and_sulfate_two():
    """Both nitrate oxygens match `acids:nitrate` and only one may be protonated: the second would take
    the nitrate to +1.  Sulfate has two charges to spend and gets both protons."""
    m = smiles('[NH3+]CC[NH3+].[O-][N+](=O)[O-]')
    assert neutralize(m) is True
    assert format(m) == '[N+](=O)(O)[O-].C([NH3+])CN'

    m = smiles('[NH3+]CC[NH3+].[O-]S(=O)(=O)[O-]')
    assert neutralize(m) is True
    assert format(m) == 'O=S(O)(=O)O.C(N)CN'


def test_the_total_charge_is_preserved():
    """The contract of `keep_charge=True`: protons move in pairs, so the net charge cannot drift."""
    for s in ['[NH3+]CC(=O)[O-]', 'C[NH3+].[Cl-]', '[NH3+]CC[NH3+].[O-][N+](=O)[O-]',
              'CC(=O)[O-].[O-]C(C)=O.[NH4+]', '[NH4+].[OH-]', 'c1cc[nH+]cc1.CC(=O)[O-]']:
        m = smiles(s)
        before = _charge(m)
        neutralize(m)
        assert _charge(m) == before, s


def test_the_graph_never_changes():
    """No atom and no bond is touched -- what separates this from `split_salts`, which cuts a bond."""
    for s in ['[NH3+]CC(=O)[O-]', 'C[NH3+].[Cl-]', '[NH4+].[OH-]', '[NH3+]CC[NH3+].[O-]S(=O)(=O)[O-]']:
        m = smiles(s)
        atoms, bonds = len(m), m.bond_count
        neutralize(m)
        assert (len(m), m.bond_count) == (atoms, bonds), s


def test_a_second_call_moves_nothing():
    m = smiles('[NH3+]CC(=O)[O-]')
    assert neutralize(m) is True
    assert neutralize(m) is False


# what stays charged


def test_a_cation_with_no_proton_keeps_its_counterion():
    """A quaternary ammonium and a metal are not `acid` rows: nothing there can pay for the anion."""
    for s in ['C[N+](C)(C)CC(=O)[O-]', '[Na+].CC(=O)[O-]', 'C[N+](C)(C)C.[Cl-]']:
        m = smiles(s)
        before = m.canonical_bytes
        assert neutralize(m) is False, s
        assert m.canonical_bytes == before, f'{s} was modified'


def test_a_lone_ion_needs_keep_charge_off():
    """One side alone cannot be paired, so `keep_charge=True` declines without deciding anything."""
    for s in ['C[NH3+]', 'CC(=O)[O-]', 'c1cc[nH+]cc1']:
        m = smiles(s)
        assert neutralize(m) is False, s
        assert neutralize(m, keep_charge=False) is True, s
        assert _charge(m) == 0, s


def test_an_unbalanced_record_comes_back_partly_neutral():
    """All-or-nothing is per site, not per record: the pair that can be made is made and the leftover
    is reported.  A dication with one chloride keeps one of its two charges."""
    m = smiles('[NH3+]CC[NH3+].[Cl-]')
    log = m.log
    assert neutralize(m) is True
    assert format(m) == 'C([NH3+])CN.Cl'
    assert _charge(m) == 1
    assert len(log.refused()) == 1
    assert 'no anion is left' in log.refused()[0]


def test_a_neutral_form_no_valence_row_accepts_is_refused():
    """The valence question goes to the shared collection, not to an `h` primitive in the table: a bare
    `[O-]` would become a one-hydrogen neutral oxygen, which no row allows."""
    m = smiles('[O-]')
    log = m.log
    assert neutralize(m, keep_charge=False) is False
    assert format(m) == '[O-]'
    assert len(log.refused()) == 1
    assert 'valence violation' in log.refused()[0]


def test_a_site_taken_past_zero_alone_is_refused():
    """`keep_charge=False` is the same rule with one end: nitric acid's remaining oxygen is not
    protonated, because its component is already at zero."""
    m = smiles('[O-][N+](=O)[O-]')
    log = m.log
    assert neutralize(m, keep_charge=False) is True
    assert _charge(m) == 0
    assert len(log.refused()) == 1
    assert 'away from zero' in log.refused()[0]


# hydrogens, log and registration


def test_explicit_hydrogens_hide_the_site():
    """`acids.tsv` reads IMPLICIT hydrogens, so a cation drawn with hydrogen atoms is invisible until
    they are folded in.  The docstring says `implicify_hydrogens()` first, and this is why."""
    m = smiles('[H][N+]([H])([H])C.[Cl-]')
    assert neutralize(m) is False
    implicify_hydrogens(m)
    assert neutralize(m) is True
    assert format(m) == 'CN.Cl'


def test_the_log_names_the_stage_and_the_row():
    m = smiles('CC(=O)[O-].[NH4+]')
    log = m.log
    neutralize(m)
    records = log.by_stage('neutralize')
    assert len(records) == 1
    assert records[0].rule == 'acids:ammonium'
    assert records[0].severity == INFO
    assert all(r.rule.startswith('acids:') for r in log)


def test_a_refusal_and_a_move_are_one_record_type():
    m = smiles('[NH3+]CC[NH3+].[Cl-]')
    log = m.log
    neutralize(m)
    assert len({type(r) for r in log}) == 1
    assert {r.severity for r in log} == {INFO, REFUSED}


def test_the_container_method_is_the_pass():
    """Registration is by injection: `MoleculeContainer` is a `cdef class` and cannot be extended."""
    a, b = smiles('[NH3+]CC(=O)[O-]'), smiles('[NH3+]CC(=O)[O-]')
    assert a.neutralize() == neutralize(b)
    assert a == b
    assert 'neutralize' in CHEMISTRY_ALL


def test_canonicalize_runs_it():
    """A zwitterion and its neutral drawing are one compound and must share a key."""
    for zwitterion, neutral in [('[NH3+]CC(=O)[O-]', 'NCC(=O)O'), ('CC(=O)[O-].[NH4+]', 'CC(=O)O.N')]:
        a, b = smiles(zwitterion), smiles(neutral)
        canonicalize(a)
        canonicalize(b)
        assert a == b, zwitterion


def test_canonicalize_leaves_a_charge_it_cannot_pair():
    """The net charge is part of the compound, so betaine and sodium acetate keep theirs."""
    for s in ['C[N+](C)(C)CC(=O)[O-]', '[Na+].CC(=O)[O-]']:
        m = smiles(s)
        canonicalize(m)
        assert _charge(m) == 0, s
        assert any(m.charge_of(n) for n in m), f'{s} lost its charges'
