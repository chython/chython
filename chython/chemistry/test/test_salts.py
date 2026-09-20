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
"""`split_salts`, `decompose_salts` and the table behind them.

`split_salts` acts on all 93 `[M]` metals, kept safe by an all-or-nothing test per cation atom, so the
refusals matter as much as the splits.  `decompose_salts` returns one `ComponentRow` per component: a
stabilizer is eligible, never leftover -- the guard fires when no parent is a non-lone-metal outside the
recognized solvents, so water alone answers water and a THF monohydrate reports both as parents.
"""
from pytest import mark, raises
# `__all__` and not the package object: `from ... import chemistry` would execute the facade, which
# `test_dependency_direction.py` ratchets against.
from .. import (ComponentRow, DEFAULT_STABILIZER_CLASSES, __all__ as CHEMISTRY_ALL, decompose_salts,
                implicify_hydrogens, neutralize, split_salts)
from .. import _salts
from .._tables import (SALT_CLASSES, SALT_MATCHES, salts_rows, salts_rows_by_klass,
                       salts_species_keys, salts_table_text)
from ...core import REFUSED, REPAIRED, MoleculeContainer, read_smiles as smiles


# split_salts


def test_the_metals_chython_two_left_alone_now_split():
    """Zn, Al, Ag, Bi and Fe carboxylates split, and to the right metal charge.

    The three acetates of the aluminium salt print as ONE string three times: the canonical order is
    per component, so two isomorphic components cannot come out spelled differently.
    """
    for s, expect in [
            ('CC(=O)O[Na]', 'C(C)([O-])=O.[Na+]'),
            ('CC(=O)O[Zn]OC(C)=O', 'C(C)([O-])=O.C(C)([O-])=O.[Zn+2]'),
            ('[Ag]OC(C)=O', 'C(C)([O-])=O.[Ag+]'),
            ('[Al](OC(C)=O)(OC(C)=O)OC(C)=O', 'C(C)([O-])=O.C(C)([O-])=O.C(C)([O-])=O.[Al+3]'),
            ('[Bi](OC(C)=O)(OC(C)=O)OC(C)=O', 'C(C)([O-])=O.C(C)([O-])=O.C(C)([O-])=O.[Bi+3]'),
            ('CC(=O)O[Zn]Cl', 'C(C)([O-])=O.[Cl-].[Zn+2]'),            # two different acceptor rows
            ('[Na][Cl]', '[Na+].[Cl-]'),
            ('OP(=O)(O[Na])O[Na]', '[O-]P(=O)(O)[O-].[Na+].[Na+]')]:
        m = smiles(s)
        assert split_salts(m) is True, s
        assert format(m) == expect, s


def test_an_alkoxide_drawn_with_the_metal_bonded_splits():
    """The drawn O-metal bond states the alkoxide, so there is nothing left to disambiguate and
    `salts:alkoxide` reads it.  Silicon and boron are outside `[M]`, so a silyl ether is untouched."""
    for s, expect in [('CCO[Na]', 'C(C)[O-].[Na+]'),                          # sodium ethoxide
                      ('CC(C)(C)O[K]', 'C(C)([O-])(C)C.[K+]'),                # potassium tert-butoxide
                      ('CCS[K]', 'C(C)[S-].[K+]'),                            # potassium ethanethiolate
                      ('C[Si](C)(C)O[Na]', 'C[Si]([O-])(C)C.[Na+]'),          # sodium trimethylsilanolate
                      ('CCO[Mg]OCC', 'C(C)[O-].C(C)[O-].[Mg+2]')]:
        m = smiles(s)
        assert split_salts(m) is True, s
        assert format(m) == expect, s
    for s in ['CCOCC', 'CSC', 'CCO[Si](C)(C)C', 'CCOB(OCC)OCC']:              # no metal, no acceptor
        m = smiles(s)
        before = m.canonical_bytes
        assert split_salts(m) is False, s
        assert m.canonical_bytes == before, s


def test_the_atom_count_never_changes():
    """The atom count is a contract of this pass; nothing in this module deletes a component."""
    for s in ['CC(=O)O[Na]', 'CC(=O)O[Zn]OC(C)=O', '[Al](OC(C)=O)(OC(C)=O)OC(C)=O', '[Na][Cl]']:
        m = smiles(s)
        n = len(m)
        split_salts(m)
        assert len(m) == n, s


def test_a_ligand_that_is_not_an_acceptor_refuses_the_whole_atom():
    """The all-or-nothing test is per atom: cisplatin's two chlorides are acceptor rows and its two
    ammines are not, so a per-bond rule would half-split it.  Nothing is cut."""
    for s in ['N[Pt](N)(Cl)Cl',                                       # cisplatin
              '[Fe]1(C2C=CC=C2)C2C=CC=C12',                           # ferrocene
              '[Fe](C#[O+])(C#[O+])(C#[O+])(C#[O+])C#[O+]']:          # iron pentacarbonyl
        m = smiles(s)
        before = m.canonical_bytes
        log = m.log
        assert split_salts(m) is False, s
        assert m.canonical_bytes == before, f'{s} was modified'
        assert len(log.refused()) == 1
        assert 'matches no acceptor row' in log.refused()[0]


def test_a_dative_bond_exempts_the_atom_and_does_not_trigger_it():
    """Order 8 exempts, it does not trigger: `standardize()` installs it to record coordination that
    must be preserved, so "order 8 means salt" would undo the pass that ran before."""
    m = smiles('[Fe]~N(C)(C)C')
    log = m.log
    assert split_salts(m) is False
    assert len(log.refused()) == 1
    assert 'dative (order 8)' in log.refused()[0]
    assert 'must be preserved' in log.refused()[0]


def test_an_untabulated_resulting_charge_refuses_the_atom_before_anything_is_cut():
    """A resulting charge outside the tabulated set refuses the whole atom, and the allowed set is
    named in the message.  Nothing is left half-split."""
    m = smiles('[W](OC(C)=O)(OC(C)=O)(OC(C)=O)(OC(C)=O)(OC(C)=O)OC(C)=O')
    before = m.canonical_bytes
    log = m.log
    assert split_salts(m) is False
    assert m.canonical_bytes == before, 'the molecule was left half-split'
    assert 'would leave it at charge 6' in log.refused()[0]
    assert '1, 2, 3, 4' in log.refused()[0], 'the message does not name the tabulated set'


def test_an_implicit_hydrogen_on_the_cation_refuses_it():
    """Cutting bonds and adding charge with an uncounted hydrogen on the atom would report a valence
    nobody stated."""
    m = smiles('[AlH2]OC(C)=O')
    log = m.log
    assert split_salts(m) is False
    assert 'carries 2 implicit hydrogen(s)' in log.refused()[0]


def test_all_or_nothing_is_per_atom_and_not_per_molecule():
    """One molecule, one metal that splits and one that refuses, both answers in the same log."""
    m = smiles('CC(=O)O[Na].N[Pt](N)(Cl)Cl')
    log = m.log
    assert split_salts(m) is True
    assert format(m) == 'N[Pt](Cl)(N)Cl.C(C)([O-])=O.[Na+]'
    assert len(log.repaired()) == 1
    assert len(log.refused()) == 1


def test_an_already_ionic_salt_and_an_inner_salt_are_left_alone():
    """Nothing to cut, and an inner salt never enters: it has no `[M]`, so no cation row claims it."""
    for s in ['CC(=O)[O-].[Na+]', '[NH3+]CC(=O)[O-]', 'C[Si](C)(C)Cl', 'CC(=O)O']:
        m = smiles(s)
        before = m.canonical_bytes
        assert split_salts(m) is False, s
        assert m.canonical_bytes == before


def test_splitting_is_idempotent():
    """The second call has nothing left to cut and says so."""
    m = smiles('CC(=O)O[Zn]OC(C)=O')
    assert split_salts(m) is True
    assert split_salts(m) is False


def test_split_salts_logs_repaired_because_the_drawing_was_wrong():
    """A covalent bond from sodium to a carboxylate oxygen is a mis-drawing, so REPAIRED."""
    m = smiles('CC(=O)O[Na]')
    log = m.log
    split_salts(m)
    assert len(log) == 1
    assert log[0].severity == REPAIRED
    assert log[0].rule == 'salts:metal'
    assert log[0].atoms == (5, 4)


def test_the_stage_is_the_pass_name():
    m = smiles('CC(=O)O[Na]')
    split_salts(m)
    assert m.log[0].stage == 'split-salts'


# the keep surface: row id, class or element symbol


def test_keep_protects_a_cation_from_splitting_too():
    for keep in [['Na'], ['salts:metal']]:
        m = smiles('CC(=O)O[Na]')
        assert split_salts(m, keep=keep) is False, keep
        assert format(m) == '[Na]OC(=O)C'
    # ... and naming a different element does not
    m = smiles('CC(=O)O[Na]')
    assert split_salts(m, keep=['K']) is True


def test_keeping_an_acceptor_row_refuses_the_cation_that_needed_it():
    """`keep=['salts:oxo-acid']` means "carboxylate O is not an acceptor", so the sodium fails the
    all-or-nothing test rather than half-splitting."""
    m = smiles('CC(=O)O[Na]')
    log = m.log
    assert split_salts(m, keep=['salts:oxo-acid']) is False
    assert 'matches no acceptor row' in log.refused()[0]


def test_a_class_name_reaches_every_row_of_that_class():
    """A class name saves a caller from naming rows one by one and going stale each time the table
    grows.  Only the two embed classes name an atom; a class with no atom leaves the split alone."""
    for keep in [['metal_cation'], ['charge_acceptor']]:
        m = smiles('CC(=O)O[Na]')
        assert split_salts(m, keep=keep) is False, keep
    # a class that names no atom leaves the split alone
    m = smiles('CC(=O)O[Na]')
    assert split_salts(m, keep=['water']) is True


def test_there_is_no_strip_keyword():
    """Nothing here deletes a component, so there is nothing for a second selector to protect against."""
    with raises(TypeError):
        split_salts(smiles('CC(=O)O[Na]'), strip=['salts:water'])


def test_a_molecule_in_keep_is_refused_with_the_reason():
    """There is no component to compare one against here, so it is refused rather than silently
    ignored -- a caller who wrote it believes something is protected."""
    with raises(ValueError) as e:
        split_salts(smiles('CC(=O)O[Na]'), keep=[smiles('Cl')])
    assert 'names a whole COMPONENT' in str(e.value)
    assert "keep=['Na']" in str(e.value), 'the refusal does not say what to write instead'
    assert "keep=['metal_cation']" in str(e.value)


def test_a_keep_item_of_the_wrong_kind_is_refused_at_the_boundary():
    with raises(TypeError) as e:
        split_salts(smiles('CC(=O)O[Na]'), keep=[11])
    assert 'row ids, classes or element symbols' in str(e.value)

    with raises(ValueError) as e:
        split_salts(smiles('CC(=O)O[Na]'), keep=['Unobtainium'])
    assert 'a row id' in str(e.value) and 'a class' in str(e.value)


def test_a_typo_in_a_row_id_is_refused():
    """A typo must not be a silent no-op: a caller who wrote it believes something is protected."""
    with raises(ValueError) as e:
        split_salts(smiles('CC(=O)O[Na]'), keep=['salts:not-a-row'])
    assert 'no such row' in str(e.value)
    with raises(ValueError):
        split_salts(smiles('CC(=O)O[Na]'), keep=['salts:METAL'])     # case matters, and so does the check


# getting the pieces, which is composition and not a third pass


def test_split_alone_misses_a_salt_that_was_drawn_covalently():
    """`split()` reports the components a molecule already has, so `CC(=O)O[Na]` is one and a corpus
    filtered on `len(m.split()) == 1` keeps it.  `split()`'s own docstring must point at `split_salts`,
    since that pointer is how a caller finds the pass that fixes it."""
    m = smiles('CC(=O)O[Na]')
    assert len(m.split()) == 1, 'split() started cutting ionic bonds, which is not its job'
    assert 'split_salts' in MoleculeContainer.split.__doc__, \
        'split() no longer points at the pass that fixes this, so nobody will find it'

    working = m.copy()
    split_salts(working)
    parts = working.split()
    assert [format(p) for p in parts] == ['C(C)([O-])=O', '[Na+]']
    assert max(parts, key=len).canonical_bytes == smiles('CC(=O)[O-]').canonical_bytes
    assert m.canonical_bytes == smiles('CC(=O)O[Na]').canonical_bytes, 'the copy did not protect it'


def test_split_already_hands_back_new_molecules():
    """New objects, stereo intact, source numbering."""
    m = smiles('C[C@H](N)C(=O)[O-].[Na+]')
    parts = m.split()
    assert isinstance(parts, list) and len(parts) == 2
    assert all(isinstance(p, MoleculeContainer) for p in parts)
    assert all(p is not m for p in parts)
    parent = max(parts, key=len)
    assert parent.canonical_bytes == smiles('C[C@H](N)C(=O)[O-]').canonical_bytes, 'stereo was lost'
    assert parent.canonical_bytes != smiles('C[C@@H](N)C(=O)[O-]').canonical_bytes


def test_an_inner_salt_is_one_piece_either_way():
    """Glycine has no cation to cut, so the composition agrees with a bare `split()` on it."""
    m = smiles('[NH3+]CC(=O)[O-]')
    assert split_salts(m) is False
    assert len(m.split()) == 1


# decompose_salts: the default split


def test_decompose_salts_changes_nothing_and_logs_nothing():
    for s in ('NCC(=O)O.OC(=O)C(F)(F)F', 'CC(=O)O[Na].O', 'N[Pt](N)(Cl)Cl', 'Cl'):
        m = smiles(s)
        log, before = m.log, bytes(m.canonical_bytes)
        m.decompose_salts()
        assert bytes(m.canonical_bytes) == before, s
        assert not len(log), s


def test_the_salt_surface_is_two_methods():
    """One that edits and one that reports.  The third, which deleted components, asked the same
    question as the reporter and answered it destructively; a consumer reads `parents` instead."""
    mol = smiles('CC(=O)O[Na].O')
    assert callable(mol.split_salts) and callable(mol.decompose_salts)
    for gone in ('split_ionic', 'strip_salts', 'salt_composition'):
        assert not hasattr(mol, gone), f'{gone} still resolves'
        assert gone not in CHEMISTRY_ALL, f'{gone} is still on the package'
    assert set(_salts.__all__) == {'ComponentRow', 'DEFAULT_STABILIZER_CLASSES',
                                   'SaltComposition', 'decompose_salts', 'split_salts'}


def test_two_enantiomers_drawn_together_are_two_parents():
    """Parents are per component, not per canonical constitution: a racemate drawn as two components is
    two compounds, and collapsing them would report a stoichiometry that was not drawn."""
    r = smiles('N[C@@H](C)C(=O)O.N[C@H](C)C(=O)O').decompose_salts()
    assert len(r.parents) == 2


#: (record, parents, stabilizer row ids).  A parent is named by the spelling the NORMALIZED probe holds --
#: neutralized, so a drawn acetate is acetic acid here and an intrinsic `[Na+]` is still `[Na+]`.  Public
#: compounds throughout.
COMPOSITIONS = (
    ('CC(=O)Oc1ccccc1C(=O)O.O', ('CC(=O)Oc1ccccc1C(=O)O',), ('salts:water',)),
    ('CCN.Cl', ('CCN',), ('salts:hcl',)),
    ('c1ccncc1.OC(=O)C(F)(F)F', ('c1ccncc1',), ('salts:tfa',)),
    ('CC(=O)[O-].[Na+]', ('CC(=O)O', '[Na+]'), ()),
    ('CC(=O)O[Na]', ('CC(=O)O', '[Na+]'), ()),
    ('CC(=O)O.[Na+]', ('CC(=O)O', '[Na+]'), ()),
    ('CC(=O)O.[Na]', ('CC(=O)O', '[Na]'), ()),
    ('C[N+](C)(C)C.[Cl-]', ('C[N+](C)(C)C', 'Cl'), ()),
    ('CC[B-](F)(F)F.[K+]', ('CC[B-](F)(F)F', '[K+]'), ()),
    ('CCCS([O-])(=O)=O.[Na+].O', ('CCCS(O)(=O)=O', '[Na+]'), ('salts:water',)),
    ('[Na]', ('[Na]',), ()),
    ('CCBr.[Zn]', ('CCBr', '[Zn]'), ()),
    ('c1ccccc1.CC(=O)Oc1ccccc1C(=O)O', ('c1ccccc1', 'CC(=O)Oc1ccccc1C(=O)O'), ()),
    ('CCO.CCCO', ('CCO', 'CCCO'), ()),
    ('O.C1CCOC1', ('O', 'C1CCOC1'), ()),
    ('CC(=O)O', ('CC(=O)O',), ()),
    ('O', ('O',), ()),
)


def _canonical(spelling):
    mol = smiles(spelling)
    mol.thiele()
    return format(mol, '!s')


@mark.parametrize('spelling, parents, stabilizers', COMPOSITIONS)
def test_the_default_split_is_the_table(spelling, parents, stabilizers):
    answer = smiles(spelling).decompose_salts()
    assert sorted(format(row.molecule, '!s') for row in answer.parents) == \
        sorted(_canonical(p) for p in parents)
    assert sorted(row.species for row in answer.stabilizers) == sorted(stabilizers)


def test_every_component_gets_exactly_one_row_and_one_role():
    answer = smiles('CCCS([O-])(=O)=O.[Na+].O').decompose_salts()
    assert len(answer.components) == 3
    assert {row.role for row in answer.components} == {'parent', 'stabilizer'}
    assert len(answer.parents) + len(answer.stabilizers) == len(answer.components)
    # `atoms` partitions the caller's heavy atoms: 7 in the sulfonate, the sodium, the water
    assert sorted(n for row in answer.components for n in row.atoms) == list(range(1, 10))


def test_the_row_carries_the_shape_and_both_charges():
    answer = smiles('CC(=O)[O-].[Na+]').decompose_salts()
    acetate, sodium = sorted(answer.components, key=lambda row: row.heavy_atoms, reverse=True)
    assert (acetate.heavy_atoms, acetate.carbon_count, acetate.ring_count) == (4, 2, 0)
    assert (acetate.charge, acetate.residual_charge) == (-1, 0)
    assert (sodium.charge, sodium.residual_charge) == (1, 1)
    assert sodium.is_lone_metal and not acetate.is_lone_metal
    assert sodium.species == 'salts:metal' and sodium.klass == 'metal_cation'


def test_is_lone_metal_is_charge_blind_across_all_four_spellings():
    for spelling in ('CC(=O)[O-].[Na+]', 'CC(=O)O[Na]', 'CC(=O)O.[Na+]', 'CC(=O)O.[Na]'):
        answer = smiles(spelling).decompose_salts()
        assert [row.is_lone_metal for row in answer.components].count(True) == 1, spelling
        assert not answer.stabilizers, spelling


def test_a_counter_ion_of_an_intrinsic_charge_stays():
    # the chloride is drawn -1 and the record's intrinsic charge is +1, so it is on duty
    answer = smiles('C[N+](C)(C)C.[Cl-]').decompose_salts()
    assert len(answer.parents) == 2 and not answer.stabilizers
    # the same chloride beside an ammonium that HAS a neutral form is a stabilizer
    answer = smiles('CC[NH3+].[Cl-]').decompose_salts()
    assert [row.species for row in answer.stabilizers] == ['salts:hcl']


#: A lone metal states no charge when it is drawn neutral, and one the drawing does not account for when the
#: record's drawn charges do not balance either way; the anion it belongs to is then on counter-ion duty and
#: never leftover -- with a third component present, the guard is not what saves it.  Water is on no duty: it
#: is neither an anion nor a `protic_acid` site.  `[Mg+]` with two acetates is the record that pins the
#: negative direction: `fix_salt_charges` calls that magnesium under-charged, so the reader must not answer it
#: one way before `standardize()` and another way after.  The last four rows are the controls: a balanced
#: sodium chloride beside a free acid, a record with no metal at all, a hydrate of a charged metal, and an
#: amine hydrochloride, whose HCl leaves because the record balances.
UNDRAWN_METAL_DUTY = (
    ('Cl.[Na].CCBr', 'Cl', ()),
    ('CC(=O)[O-].[Na].CCBr', 'C(C)(=O)O', ()),
    ('[O-]S(=O)(=O)C.[Na].CCBr', 'CS(=O)(O)=O', ()),
    ('CC(=O)O.[Na].O.CCBr', 'C(C)(=O)O', ('salts:water',)),
    ('CC(=O)O.[Na+].CCBr', 'C(C)(=O)O', ()),
    ('Cl.[Na+].CCBr', 'Cl', ()),
    ('CC(=O)O.[Mg+2].CCBr', 'C(C)(=O)O', ()),
    ('CC(=O)[O-].CC(=O)[O-].[Mg+].CCBr', 'C(C)(=O)O', ()),
    ('[Na+].[Cl-].[Cl-].CCBr', 'Cl', ()),
    ('[Na+].[Cl-].CC(=O)O.CCBr', 'Cl', ('salts:acetic',)),
    ('CCN.Cl.CCBr', 'C(C)N', ('salts:hcl',)),
    ('O.[Na+].CCBr', 'C(C)Br', ('salts:water',)),
    ('CC[NH3+].[Cl-]', 'C(C)N', ('salts:hcl',)),
)


@mark.parametrize('spelling, parent, stabilizers', UNDRAWN_METAL_DUTY)
def test_an_anion_whose_metal_lacks_its_charge_is_never_leftover(spelling, parent, stabilizers):
    answer = smiles(spelling).decompose_salts()
    assert parent in [str(row.molecule) for row in answer.parents], spelling
    assert tuple(row.species for row in answer.stabilizers) == stabilizers, spelling


def test_an_under_charged_metal_answers_the_same_before_and_after_standardize():
    """The records `fix_salt_charges` repairs must not partition one way before it and another way after.

    The metal's own drawn charge is what the repair moves, so the comparison is over the partition: which
    components are parents, which leave, and every parent that is not the metal.
    """
    for spelling in ('CC(=O)[O-].CC(=O)[O-].[Mg+].CCBr', 'CC(=O)O.[Na].CCBr', 'CC(=O)O.[Na+].CCBr'):
        mol = smiles(spelling)
        drawn = smiles(spelling).decompose_salts()
        mol.standardize()
        repaired = mol.decompose_salts()
        assert ([str(row.molecule) for row in drawn.parents if not row.is_lone_metal]
                == [str(row.molecule) for row in repaired.parents if not row.is_lone_metal]), spelling
        assert len(drawn.parents) == len(repaired.parents), spelling
        assert ([row.species for row in drawn.stabilizers]
                == [row.species for row in repaired.stabilizers]), spelling


def test_a_sodium_sulfonate_monohydrate_keeps_its_pair_and_loses_its_water():
    answer = smiles('CCCS([O-])(=O)=O.[Na+].O').decompose_salts()
    assert [row.species for row in answer.stabilizers] == ['salts:water']
    assert len(answer.parents) == 2


def test_an_organometallic_is_never_a_stabilizer():
    answer = smiles('CC[B-](F)(F)F.[K+]').decompose_salts()
    borate = next(row for row in answer.components if row.is_organometallic)
    assert borate.role == 'parent' and borate.residual_charge == -1


def test_the_guard_makes_every_component_a_parent_when_all_are_formers():
    for spelling in ('O', 'O.O', 'O.Cl'):
        answer = smiles(spelling).decompose_salts()
        assert answer.parents and not answer.stabilizers, spelling


def test_classes_and_max_atoms_and_discardable_widen_the_split():
    record = 'CC(=O)Oc1ccccc1C(=O)O.c1ccccc1'
    assert not smiles(record).decompose_salts().stabilizers
    widened = smiles(record).decompose_salts(classes=DEFAULT_STABILIZER_CLASSES + ('hydrocarbon',))
    assert [row.species for row in widened.stabilizers] == ['salts:benzene']
    capped = smiles(record).decompose_salts(classes=DEFAULT_STABILIZER_CLASSES + ('hydrocarbon',),
                                            max_atoms=5)
    assert not capped.stabilizers
    named = smiles(record).decompose_salts(discardable=('salts:benzene',))
    assert [row.species for row in named.stabilizers] == ['salts:benzene']


def test_equivalents_counts_repeated_components():
    answer = smiles('CCN.Cl.Cl').decompose_salts()
    assert {row.equivalents for row in answer.stabilizers} == {2}
    assert answer.equivalents_by_species() == {'salts:hcl': 2}


def test_equivalents_counts_by_role_not_just_by_structure():
    """A component on counter-ion duty and a free copy share the canonical key but not the role.

    `C[N+](C)(C)C.[Cl-].Cl`: `[Cl-]` is on duty for TMA's intrinsic +1 charge, the second `Cl` is a
    free stabilizer.  Both normalize to `Cl`, but their roles differ, so each has equivalents=1.
    """
    answer = smiles('C[N+](C)(C)C.[Cl-].Cl').decompose_salts()
    assert answer.equivalents_by_species() == {'salts:hcl': 1}
    assert any(row.species == 'salts:hcl' and row.role == 'parent' for row in answer.components)


def test_the_record_charges_are_both_reported():
    answer = smiles('CC(=O)[O-].[Na+]').decompose_salts()
    assert (answer.charge, answer.residual_charge) == (0, 1)
    answer = smiles('CC(=O)O.[Na+]').decompose_salts()
    assert (answer.charge, answer.residual_charge) == (1, 1)


# the container methods, and the table


def test_the_two_are_container_methods_and_agree_with_the_functions():
    """The molecule's facade needs methods, which is why the injection hook exists."""
    for s in ['CC(=O)O[Na]', 'N[Pt](N)(Cl)Cl']:
        a, b = smiles(s), smiles(s)
        assert a.split_salts() == split_salts(b)
        assert a.canonical_bytes == b.canonical_bytes
    for s in ['CC(=O)O[Na]', 'NCC(=O)O.OC(=O)C(F)(F)F', 'Cl']:
        a, b = smiles(s), smiles(s)
        assert a.decompose_salts() == decompose_salts(b)


def test_the_container_method_forwards_every_keyword():
    record = 'CC(=O)Oc1ccccc1C(=O)O.c1ccccc1'
    assert not smiles(record).decompose_salts().stabilizers
    assert [row.species for row in
            smiles(record).decompose_salts(
                classes=DEFAULT_STABILIZER_CLASSES + ('hydrocarbon',)).stabilizers] == ['salts:benzene']
    assert not smiles(record).decompose_salts(classes=DEFAULT_STABILIZER_CLASSES + ('hydrocarbon',),
                                              max_atoms=5).stabilizers
    assert [row.species for row in
            smiles(record).decompose_salts(discardable=('salts:benzene',)).stabilizers] == ['salts:benzene']


def test_the_container_method_and_the_function_agree():
    mol = smiles('CCN.Cl.O')
    assert mol.decompose_salts().tags == decompose_salts(mol).tags


def test_a_species_row_compiles_to_a_stereo_free_canonical_key():
    keys = salts_species_keys()
    assert keys[format(smiles('OC(=O)C(F)(F)F'), '!s')].id == 'salts:tfa'


def test_the_key_names_a_constitution_and_not_a_stereoisomer():
    """One tartrate row, both enantiomers and the meso form."""
    row = salts_species_keys()[format(smiles('OC(C(O)C(O)=O)C(O)=O'), '!s')]
    assert row.id == 'salts:tartaric'
    for s in ['O[C@@H]([C@H](O)C(O)=O)C(O)=O', 'O[C@H]([C@H](O)C(O)=O)C(O)=O']:
        assert format(smiles(s), '!s') == row.key


def test_the_widened_corpus_recognizes_each_new_former():
    """A round trip through the key index covers every row at once and names the one that broke."""
    for pattern, expected in [('OC(=O)CC(O)(CC(O)=O)C(O)=O', 'salts:citric'),
                              ('OC(=O)CCC(O)=O', 'salts:succinic'),
                              ('OS(=O)(=O)c1ccccc1', 'salts:besylic'),
                              ('OCC(N)(CO)CO', 'salts:tromethamine'),
                              ('CNCC(O)C(O)C(O)C(O)CO', 'salts:meglumine'),
                              ('C[N+](C)(C)CCO', 'salts:choline'),
                              ('CN1CCCC1=O', 'salts:nmp')]:
        probe = smiles(pattern)
        probe.thiele()
        assert salts_species_keys()[format(probe, '!s')].id == expected, pattern


def test_every_row_of_the_table_round_trips_through_its_own_key():
    """The table is its own fixture: a row whose pattern does not key back to it is unreachable."""
    for row in salts_rows():
        if row.key is None:
            continue
        probe = smiles(row.pattern)
        probe.thiele()
        assert salts_species_keys()[format(probe, '!s')].id == row.id, row.id


def test_pyridine_is_tabulated_as_an_amine_base_not_a_solvent():
    """Pyridine under a solvate class would read pyridine hydrochloride as hydrochloric acid."""
    assert salts_species_keys()[format(smiles('c1ccncc1'), '!s')].klass == 'amine_base'


def test_a_kekule_drawing_still_keys_to_its_row():
    """`thiele()` runs on the copy before keying, so a Kekule toluene keys to the tabulated aromatic row."""
    for s in ['NCC(=O)O.Cc1ccccc1', 'NCC(=O)O.C1=CC=CC=C1C']:
        answer = smiles(s).decompose_salts()
        assert any(row.species == 'salts:toluene' for row in answer.components), s


def test_an_explicit_hydrogen_does_not_hide_a_stabilizer():
    """`implicify_hydrogens()` runs on the copy before keying, so an H-atom drawing still keys."""
    a = smiles('NCC(=O)O.[H]OC(=O)C(F)(F)F').decompose_salts()
    b = smiles('NCC(=O)O.OC(=O)C(F)(F)F').decompose_salts()
    assert (sorted(row.species for row in a.stabilizers) ==
            sorted(row.species for row in b.stabilizers) == ['salts:tfa'])


def test_an_isotope_labelled_solvate_does_not_match_the_unlabelled_row():
    """No `clean_isotopes()`: D2O is not water, so it is an unrecognized parent, not a hydrate."""
    answer = smiles('NCC(=O)O.[2H]O[2H]').decompose_salts()
    assert not any(row.species == 'salts:water' for row in answer.components)
    assert len(answer.parents) == 2


def test_a_coordination_complex_is_one_parent():
    """The `split_salts` all-or-nothing test reaches here: cisplatin comes back intact as one component."""
    answer = smiles('N[Pt](N)(Cl)Cl').decompose_salts()
    assert len(answer.parents) == 1 and not answer.stabilizers


def test_two_drawn_equivalents_report_one_row_per_drawn_copy():
    """Each drawn copy gets its own row; `equivalents` says how many copies play the same role."""
    answer = smiles('OC(=O)c1ccccc1O.OC(=O)c1ccccc1O').decompose_salts()
    assert len(answer.parents) == 2
    assert all(row.equivalents == 2 for row in answer.parents)


def test_no_two_rows_share_a_key():
    assert len(salts_species_keys()) == sum(1 for row in salts_rows() if row.key is not None)


def test_the_ionic_conjugates_are_gone_because_neutralize_reaches_them():
    ids = {row.id for row in salts_rows()}
    assert not ids & {'salts:acetate', 'salts:chloride', 'salts:tosylate', 'salts:ammonium'}
    assert 'salts:hf' in ids                       # fluoride's neutral twin, which was missing


def test_a_row_no_longer_carries_a_keep_flag():
    assert not hasattr(salts_rows()[0], 'keep')
    for line in salts_table_text().splitlines():
        if line.startswith('id\t'):
            assert line.split('\t') == ['id', 'match', 'klass', 'pattern', 'charges', 'order', 'comment']
            break
    else:
        raise AssertionError('salts.tsv has no column header')


def test_every_row_declares_a_known_match_and_class():
    for row in salts_rows():
        assert row.match in SALT_MATCHES, row.id
        assert row.klass in SALT_CLASSES, row.id


def test_embed_rows_carry_a_query_and_whole_rows_carry_a_key():
    for row in salts_rows():
        if row.match == 'embed':
            assert row.query is not None and row.key is None, row.id
            assert row.heavy_atoms == 0, row.id
        else:
            assert row.key is not None and row.query is None, row.id
            assert row.heavy_atoms > 0, row.id


def test_heavy_atoms_counts_the_species():
    rows = {row.id: row for row in salts_rows()}
    assert rows['salts:water'].heavy_atoms == 1
    assert rows['salts:acetic'].heavy_atoms == 4
    assert rows['salts:tartaric'].heavy_atoms == 10


def test_classes_are_grouped_and_none_is_missing():
    grouped = salts_rows_by_klass()
    assert set(grouped) == set(SALT_CLASSES)
    assert sum(len(rows) for rows in grouped.values()) == len(salts_rows())
    assert len(grouped['metal_cation']) == 1
    assert len(grouped['charge_acceptor']) == 6
    assert len(grouped['metal_protic']) == 2
    assert len(grouped['water']) == 1


def test_the_class_census_matches_the_table():
    expected = {'mineral_acid': 12, 'sulfonic_acid': 10, 'short_carboxylic_acid': 6,
                'carboxylic_acid': 18, 'aromatic_acid': 6, 'fatty_acid': 4, 'amino_acid': 4,
                'amine_base': 20, 'quaternary_ammonium': 1, 'water': 1, 'alcohol': 5,
                'hydrocarbon': 6, 'halo_solvent': 3, 'aprotic_solvent': 14}
    grouped = salts_rows_by_klass()
    assert {k: len(grouped[k]) for k in expected} == expected
    assert sum(expected.values()) == 110


def test_only_the_metal_cation_row_carries_charges():
    for row in salts_rows():
        assert bool(row.charges) == (row.klass == 'metal_cation'), row.id


def test_only_a_runged_row_carries_an_order():
    """One rung scale over the two classes `fix_salt_charges` ranks together, and no rung anywhere else."""
    for row in salts_rows():
        assert bool(row.order) == (row.klass in ('protic_acid', 'metal_protic')), row.id


def test_the_metal_row_covers_every_metal_the_core_calls_one():
    """93 metals: `[M]`'s membership is the core's business and this row inherits it, so a change
    there shows up here as a diff."""
    row, = salts_rows_by_klass()['metal_cation']
    assert row.pattern == '[M;*:1]'
    hits = 0
    for z in range(1, 119):
        probe = MoleculeContainer()
        probe.add_atom(z)
        if row.query.is_substructure(probe):
            hits += 1
    assert hits == 93
    # and `*` is what lets it see a charged one
    charged = smiles('[Na+]')
    assert row.query.is_substructure(charged)


def test_a_log_gives_one_filterable_record_type():
    m = smiles('CC(=O)O[Na].N[Pt](N)(Cl)Cl')
    log = m.log
    split_salts(m)
    assert len({type(r) for r in log}) == 1
    assert all(r.rule.startswith('salts:') for r in log)
    assert {r.severity for r in log} == {REPAIRED, REFUSED}


def test_no_log_is_the_same_answer_as_a_log():
    for s in ['CC(=O)O[Na]', 'N[Pt](N)(Cl)Cl', 'CN.Cl', '[Na+].[Cl-]', 'CC(=O)O[Na].O']:
        a, b = smiles(s), smiles(s)
        assert split_salts(a) == split_salts(b)
        assert a.decompose_salts() == b.decompose_salts()


def test_a_record_is_substring_matchable():
    """A caller filtering a log by text, rather than by rule id, is not broken."""
    m = smiles('CC(=O)O[Na]')
    log = m.log
    split_salts(m)
    assert len(log) == 1
    assert 'was ionic' in log[0]


#: Every `protic_acid` row, with one public compound each row must fire on.  The rung is the pass's
#: proton-choice order, so it is pinned here beside the pattern rather than left to file position.
ACID_LADDER = (
    ('salts:sulfonic-oh', 1, 'CS(=O)(=O)O'),
    ('salts:oxo-acid-oh', 2, 'OP(=O)(O)O'),
    ('salts:nitric-oh', 2, 'O[N+](=O)[O-]'),
    ('salts:acyl-sulfonamide', 2, 'O=C1NS(=O)(=O)c2ccccc21'),
    ('salts:sulfinic-oh', 2, 'CS(=O)O'),
    ('salts:thio-acid-sh', 2, 'CCOC(=S)S'),
    ('salts:hydrogen-halide', 3, 'Cl'),
    ('salts:carboxylic-oh', 4, 'CC(=O)O'),
    ('salts:tetrazole-1h', 4, 'c1nnn[nH]1'),
    ('salts:tetrazole-2h', 4, 'c1nn[nH]n1'),
    ('salts:thiophenol-sh', 4, 'Sc1ccccc1'),
    ('salts:phenol-oh', 5, 'Oc1ccccc1'),
    ('salts:imide', 5, 'O=C1CCC(=O)N1'),
    ('salts:sulfonamide', 5, 'Cc1ccc(cc1)S(N)(=O)=O'),
    ('salts:thiol-sh', 5, 'CCS'),
)

#: Every `metal_protic` row, on the same rung scale and below every acid.  A separate ladder because a
#: separate class: these fire for a free s-block metal and for nothing else in the library.
METAL_PROTIC_LADDER = (
    ('salts:water-oh', 6, 'O'),
    ('salts:alcohol-oh', 6, 'CCO'),
)

#: Nothing here is a salt-forming acid, and no `protic_acid` row may fire on any of it.  Every azole people
#: draw is in the list: only tetrazole sits in the salt-forming range, so every other one is a refusal by
#: design.  Water and the alcohols are here too -- they are `metal_protic`, which is a different question
#: and a different class, so a row of this table must not call either an acid.
NOT_ACIDS = (
    'c1cc[nH]c1', 'c1cnc[nH]1', 'c1cn[nH]c1', 'c1cn[nH]n1', 'c1c[nH]nn1', 'c1nc[nH]n1',
    'c1ccc2[nH]nnc2c1', 'c1ccc2[nH]ccc2c1', 'Cn1c(=O)c2[nH]cnc2n(C)c1=O',
    'CC(N)=O', 'c1ccccc1C(N)=O', 'CNC(C)=O', 'NC(N)=O', 'Nc1ccccc1', 'NO', 'CC=NO', 'CC(=O)NO',
    'CCO', 'O', 'CO', 'CC(C)(C)O', 'OC1CCCCC1', 'OCC1OC(O)C(O)C(O)C1O',
    'CC(C)=O', 'CS(C)=O', 'CC(=O)OC', 'COS(C)(=O)=O', 'C1CCOC1', 'c1ccncc1',
    'c1ccccc1[N+](=O)[O-]', 'C[N+](C)(C)[O-]', 'Clc1ccccc1', 'CCCl',
    'CSC', 'CCSSCC', 'CSc1ccccc1', 'CC(=O)SC',
)

#: What a `metal_protic` row must not fire on.  Every other O-H the table already places -- a phenol and a
#: carboxylic acid are acids and rank above these -- plus the ethers, which hold no hydrogen at all.
NOT_METAL_PROTIC = (
    'Oc1ccccc1', 'CC(=O)O', 'CS(=O)(=O)O', 'ON=O', 'CCOCC', 'C1CCOC1', 'COC', 'C[Si](C)(C)OC',
    'CCS', 'CC(=O)OC', 'CC(C)=O', 'NO', 'OO',
)


def test_the_acid_ladder_is_the_table():
    rows = salts_rows_by_klass()['protic_acid']
    assert [(row.id, row.order) for row in rows] == [(i, o) for i, o, _ in ACID_LADDER]
    assert [row.order for row in rows] == sorted(row.order for row in rows)


def test_every_acid_row_fires_on_its_compound():
    rows = {row.id: row for row in salts_rows_by_klass()['protic_acid']}
    for row_id, _, spelling in ACID_LADDER:
        mol = smiles(spelling)
        mol.standardize()
        mol.thiele()
        assert next(rows[row_id].query.get_mapping(mol), None) is not None, row_id


def test_no_acid_row_fires_on_a_non_acid():
    rows = salts_rows_by_klass()['protic_acid']
    for spelling in NOT_ACIDS:
        mol = smiles(spelling)
        mol.standardize()
        mol.thiele()
        fired = [row.id for row in rows if next(row.query.get_mapping(mol), None) is not None]
        assert not fired, (spelling, fired)


def test_the_metal_protic_ladder_is_the_table():
    rows = salts_rows_by_klass()['metal_protic']
    assert [(row.id, row.order) for row in rows] == [(i, o) for i, o, _ in METAL_PROTIC_LADDER]
    acids = salts_rows_by_klass()['protic_acid']
    assert min(row.order for row in rows) > max(row.order for row in acids), \
        'a metal_protic rung must sit below every acid: the metal spends itself on the acid first'


def test_every_metal_protic_row_fires_on_its_compound():
    rows = {row.id: row for row in salts_rows_by_klass()['metal_protic']}
    for row_id, _, spelling in METAL_PROTIC_LADDER:
        mol = smiles(spelling)
        mol.standardize()
        mol.thiele()
        assert next(rows[row_id].query.get_mapping(mol), None) is not None, row_id


def test_no_metal_protic_row_fires_on_an_acid_or_an_ether():
    rows = salts_rows_by_klass()['metal_protic']
    for spelling in NOT_METAL_PROTIC:
        mol = smiles(spelling)
        mol.standardize()
        mol.thiele()
        fired = [row.id for row in rows if next(row.query.get_mapping(mol), None) is not None]
        assert not fired, (spelling, fired)


def test_a_nitrophenol_is_a_phenol_and_not_a_nitric_acid():
    rows = {row.id: row for row in salts_rows_by_klass()['protic_acid']}
    mol = smiles('c1cc(O)ccc1[N+](=O)[O-]')
    mol.standardize()
    mol.thiele()
    assert next(rows['salts:phenol-oh'].query.get_mapping(mol), None) is not None
    assert next(rows['salts:nitric-oh'].query.get_mapping(mol), None) is None


# tags


#: (record, tags that must be present).  Containment and not equality: §14's table names the tags each
#: record is interesting for, and a derived tag set is free to carry more -- `CC(=O)O.[Na+]` is both
#: `charge_unbalanced` and an `ion_pair`.  Exclusivity is asserted separately, where it is a claim.
TAGGINGS = (
    ('CC(=O)Oc1ccccc1C(=O)O.O', {'hydrate'}),
    ('CCN.Cl', {'acid_salt'}),
    ('c1ccncc1.OC(=O)C(F)(F)F', {'acid_salt'}),
    ('CCN.CCS', {'acid_salt'}),                           # a thiol is an acid site, so this is a salt
    ('CCN.Sc1ccccc1', {'acid_salt'}),
    ('CCN.CCOC(=S)S', {'acid_salt'}),                     # xanthic acid
    ('CCN.CS(=O)O', {'acid_salt'}),                       # sulfinic: the sulfonic row is z5, this S is z2
    ('CCN.CCO', {'solvate'}),                             # and an alcohol is still no acid site
    ('CCS.[Na].CCBr', {'metal_salt', 'acid_salt'}),       # the thiol is what the lone sodium owes to
    ('CC(=O)[O-].[Na+]', {'metal_salt', 'ion_pair'}),
    ('CC(=O)O[Na]', {'metal_salt', 'ion_pair'}),
    ('CC(=O)O.[Na+]', {'metal_salt', 'charge_unbalanced'}),
    ('CC(=O)O.[Na]', {'metal_salt', 'charges_undrawn'}),
    ('C[N+](C)(C)C.[Cl-]', {'acid_salt', 'ion_pair'}),
    ('CC[B-](F)(F)F.[K+]', {'ion_pair'}),
    ('CCCS([O-])(=O)=O.[Na+].O', {'metal_salt', 'hydrate'}),
    ('[Na]', {'single', 'elemental_metal'}),
    ('CCBr.[Zn]', {'elemental_metal'}),
    ('c1ccccc1.CC(=O)Oc1ccccc1C(=O)O', {'solvate'}),
    ('CCO.CCCO', {'solvate'}),
    ('O.C1CCOC1', {'hydrate', 'solvate', 'competing_formers', 'stabilizer_only'}),
    ('CC(=O)O', {'single'}),
    ('O', {'single', 'stabilizer_only'}),
)


@mark.parametrize('spelling, tags', TAGGINGS)
def test_the_tags_are_the_table(spelling, tags):
    assert tags <= smiles(spelling).decompose_salts().tags, spelling


def test_a_borate_potassium_pair_is_not_a_metal_salt():
    tags = smiles('CC[B-](F)(F)F.[K+]').decompose_salts().tags
    assert 'metal_salt' not in tags and 'elemental_metal' not in tags


def test_mixture_is_what_is_left_when_the_table_recognizes_nothing():
    tags = smiles('CCBr.c1ccc2c(c1)cccc2C#N').decompose_salts().tags
    assert 'mixture' in tags
    for spelling in ('CC(=O)[O-].[Na+]', 'CCN.Cl', 'CC(=O)Oc1ccccc1C(=O)O.O',
                     'c1ccccc1.CC(=O)Oc1ccccc1C(=O)O'):
        assert 'mixture' not in smiles(spelling).decompose_salts().tags, spelling


def test_a_tag_does_not_flip_when_classes_widen():
    record = 'CC(=O)Oc1ccccc1C(=O)O.c1ccccc1'
    narrow = smiles(record).decompose_salts().tags
    wide = smiles(record).decompose_salts(classes=DEFAULT_STABILIZER_CLASSES + ('hydrocarbon',)).tags
    assert 'solvate' in narrow and narrow == wide


def test_a_solvent_tag_needs_the_component_drawn_neutral():
    """`klass` is read off the neutralized probe and the solvent tags are not.

    Sodium hydroxide keys as water, because a conjugate is not a row and that is what makes the inventory
    small; it is the counterion all the same, and `hydrate` is a claim about solvent of crystallization.
    The roles are untouched -- both components stay parents and nothing is stripped.
    """
    for spelling in ('[OH-].[Na+]', 'CC[O-].[Na+]', 'CC(C)(C)[O-].[K+]', 'C[O-].[Na+]'):
        answer = smiles(spelling).decompose_salts()
        assert 'hydrate' not in answer.tags and 'solvate' not in answer.tags, spelling
        assert {'metal_salt', 'ion_pair'} <= answer.tags, spelling
        assert len(answer.parents) == 2 and not answer.stabilizers, spelling
    # and the drawn-neutral water is still the hydrate it was
    for spelling in ('CC(=O)Oc1ccccc1C(=O)O.O', 'CCCS([O-])(=O)=O.[Na+].O', 'O.[Na]'):
        assert 'hydrate' in smiles(spelling).decompose_salts().tags, spelling


def test_an_attachment_point_is_a_parent_like_any_other():
    """A supported-synthesis drawing keys and desalts like anything else: the `[R]` is skeleton.

    `[R]` alone is `single` and no more -- it is neither a lone metal nor a recognized solvent, so the
    parent guard has nothing to promote back.
    """
    answer = smiles('[R]CCN.Cl').decompose_salts()
    assert [format(row.molecule, '!s') for row in answer.parents] == ['C(N)C[R]']
    assert [row.klass for row in answer.stabilizers] == ['mineral_acid']
    assert answer.tags == frozenset({'acid_salt'})
    assert 'hydrate' in smiles('[R]CCN.O').decompose_salts().tags
    assert smiles('[R]').decompose_salts().tags == frozenset({'single'})


def test_the_three_records_the_default_keeps_whole():
    api = 'CC(=O)Oc1ccccc1C(=O)O'                                  # aspirin: a carboxylic acid site
    for record, tag in ((f'{api}.Cc1ccccc1', 'solvate'),
                        (f'{api}.NCCCCC(N)C(=O)O', 'base_salt'),
                        (f'{api}.CCN(CC)CC', 'base_salt')):
        answer = smiles(record).decompose_salts()
        assert tag in answer.tags, record
        assert len(answer.parents) == 2 and not answer.stabilizers, record


# compose: putting a selection back together


def _normalized(spelling):
    """The form `decompose_salts()` compares against: five steps in the same order as the function body."""
    mol = smiles(spelling)
    implicify_hydrogens(mol)
    split_salts(mol)
    mol.thiele()
    neutralize(mol, keep_charge=False)
    mol.thiele()
    return mol


#: Each record is a stereo-carrying parent drawn beside a water.  `compose(parents)` must equal the parent
#: drawn dry -- so the four kinds of stereo chython stores each get a row.
COMPOSABLE = (
    ('C[C@H](N)C(=O)O.O', 'C[C@H](N)C(=O)O'),                       # a plain tetrahedral centre
    ('C[C@@H](O)c1ccccc1.O', 'C[C@@H](O)c1ccccc1'),                 # a second, on an aromatic
    ('[13CH3]C(=O)O.O', '[13CH3]C(=O)O'),                           # an isotope
    ('C/C=C/C(=O)O.O', 'C/C=C/C(=O)O'),                             # an E/Z bond
)


@mark.parametrize('record, dry', COMPOSABLE)
def test_compose_of_the_parents_is_the_parent_drawn_dry(record, dry):
    answer = smiles(record).decompose_salts()
    assert answer.compose(answer.parents).canonical_bytes == _normalized(dry).canonical_bytes


def test_compose_takes_any_selection_including_a_stabilizer_back():
    answer = smiles('CCN.Cl.O').decompose_salts()
    assert answer.compose(answer.components).canonical_bytes == _normalized('CCN.Cl.O').canonical_bytes
    water = next(row for row in answer.stabilizers if row.species == 'salts:water')
    rebuilt = answer.compose(answer.parents + (water,))
    assert rebuilt.canonical_bytes == _normalized('CCN.O').canonical_bytes


def test_compose_of_one_row_is_that_row():
    answer = smiles('CCN.Cl').decompose_salts()
    only = answer.parents[0]
    assert answer.compose((only,)).canonical_bytes == only.molecule.canonical_bytes


def test_compose_refuses_an_empty_selection_and_a_foreign_row():
    answer = smiles('CCN.Cl').decompose_salts()
    other = smiles('CCO').decompose_salts()
    with raises(ValueError):
        answer.compose(())
    with raises(ValueError):
        answer.compose(other.components)


def test_compose_unions_two_stereo_parents():
    # both parents carry a tetrahedral centre, so the union loop runs with parities in play
    answer = smiles('C[C@H](N)C(=O)O.C[C@@H](O)c1ccccc1.O').decompose_salts()
    assert len(answer.parents) == 2                          # or the copy path would carry the test
    assert answer.compose(answer.parents).canonical_bytes == \
        _normalized('C[C@H](N)C(=O)O.C[C@@H](O)c1ccccc1').canonical_bytes


# standardize integration


def test_standardize_runs_the_stage_and_recomputes_the_hydrogen_counts():
    mol = smiles('CC(=O)O.[Na]')
    assert mol.standardize()
    assert mol.canonical_bytes == smiles('CC(=O)[O-].[Na+]').canonical_bytes
    # the deprotonated oxygen's count comes from calc_implicit, not from the stage
    oxygen = next(atom for atom in mol.atoms() if atom.element == 8 and atom.charge == -1)
    assert oxygen.implicit_h == 0


def test_a_deprotonated_azole_keeps_the_ring_hydrogens():
    """The recomputation never lands on the one atom class the ring decides.

    `calc_implicit` on an uncharged aromatic pnictogen answers `H_UNKNOWN`, so a stage recomputing one
    would erase a stated N-H.  The ids this stage returns are free metals and sites it has just charged,
    and a charged `[n-]` has its class fixed by the charge.
    """
    mol = smiles('c1nnn[nH]1.[Na]')
    assert mol.standardize()
    assert mol.canonical_bytes == smiles('c1nnn[n-]1.[Na+]').canonical_bytes
    assert sum(atom.total_h for atom in mol.atoms()) == 1         # the ring CH, and nothing lost
    assert not mol.check_valence()


def test_the_stage_order_is_forced_by_the_zinc_record():
    mol = smiles('CC[Zn].[Cl-].CC(=O)O.[Na]')
    assert mol.standardize()
    assert mol.canonical_bytes == smiles('CC[Zn]Cl.CC(=O)[O-].[Na+]').canonical_bytes


def test_standardize_reproduces_the_four_spellings():
    expected = {'CC(=O)O.[Na]': 'CC(=O)[O-].[Na+]',
                'CC(=O)O.[Na+]': 'CC(=O)[O-].[Na+]',
                'CC(=O)[O-].[Na+]': 'CC(=O)[O-].[Na+]',
                'CC(=O)O[Na]': 'CC(=O)O[Na]'}          # covalent stays covalent: split_salts() owns that
    for drawn, after in expected.items():
        mol = smiles(drawn)
        mol.standardize()
        assert mol.canonical_bytes == smiles(after).canonical_bytes, drawn
        assert not mol.check_valence(), drawn


def test_canonicalize_does_not_undo_the_repair():
    mol = smiles('CC(=O)O.[Na]')
    mol.canonicalize()
    assert mol.canonical_bytes == smiles('CC(=O)[O-].[Na+]').canonical_bytes


def test_a_free_s_block_metal_beside_water_or_an_alcohol_is_the_alcoholate():
    """Through `standardize()` and through `canonicalize()`, which reaches this stage by running it: a
    hydroxide's own hydrogen comes back from `calc_implicit` like any other."""
    for drawn, after in (('CCO.[Na]', 'CC[O-].[Na+]'),
                         ('O.[Na]', '[OH-].[Na+]'),
                         ('CC(C)(C)O.[K]', 'CC(C)(C)[O-].[K+]'),
                         ('CCS.[Na]', 'CC[S-].[Na+]'),
                         ('C[Si](C)(C)O.[Na]', 'C[Si](C)(C)[O-].[Na+]'),
                         ('CCO.[Zn]', 'CCO.[Zn]')):        # not the s block: left as drawn
        for pipeline in ('standardize', 'canonicalize'):
            mol = smiles(drawn)
            getattr(mol, pipeline)()
            assert mol.canonical_bytes == smiles(after).canonical_bytes, (drawn, pipeline)
            assert not mol.check_valence(), (drawn, pipeline)
    mol = smiles('O.[Na]')
    mol.standardize()
    assert next(atom for atom in mol.atoms() if atom.element == 8).total_h == 1


def test_the_log_carries_one_stage_name_for_one_call():
    mol = smiles('CC(=O)O.[Na]')
    mol.standardize()
    assert {record.stage for record in mol.log} == {'standardize'}


# fix_salt_charges


#: §12.2's worked table, one row per case.  `None` as the answer means the record comes back as drawn and
#: the refusal is logged -- the invariant being that a refusal writes nothing at all, not that it writes
#: something harmless.
CHARGE_FIXES = (
    ('CC(=O)O.[Na]', 'CC(=O)[O-].[Na+]'),                   # step 1, then case 3
    ('CC(=O)O.[Na+]', 'CC(=O)[O-].[Na+]'),                  # case 3
    ('CC(=O)[O-].[Na+]', 'CC(=O)[O-].[Na+]'),               # case 1, nothing written
    ('CC(=O)O.CC(=O)O.[Mg]', 'CC(=O)[O-].CC(=O)[O-].[Mg+2]'),
    ('O=S(=O)([O-])[O-].[Mg+2]', 'O=S(=O)([O-])[O-].[Mg+2]'),   # case 1
    ('CS(=O)(=O)O.[Na]', 'CS(=O)(=O)[O-].[Na+]'),           # rung 1
    ('Cl.[Na]', '[Cl-].[Na+]'),                              # rung 3
    ('Oc1ccccc1.[Na]', '[O-]c1ccccc1.[Na+]'),                # rung 5
    ('CCS.[Na]', 'CC[S-].[Na+]'),                            # rung 5: a thiol, for any metal at all
    ('CCO.[Na]', 'CC[O-].[Na+]'),                            # rung 6: only an s-block metal opens this
    ('O.[Na]', '[OH-].[Na+]'),
    ('O.O.[Mg]', '[OH-].[OH-].[Mg+2]'),                      # two equivalents, two sites
    ('CC(=O)O.CCO.[Na]', 'CC(=O)[O-].CCO.[Na+]'),            # rung 4 beats rung 6: the acid spends it
    ('CCO.[Sc]', None),                                      # d block: an alcohol is no site for it
    ('CCO.[Al]', None),                                      # p block: the same
    ('CCO.[Be]', None),                                      # the s-block metal that reduces neither
    ('CC(=O)O.[Mg]', None),                                  # case 4: one equivalent, two wanted
    ('CC(=O)O.[Zn]', None),                                  # step 1 refused
    ('CCCCCC.[Na+]', None),                                  # case 4: nothing to sit on
    ('[Na]', None),                                          # case 4: a lone metal is the metal
    ('[Ce]', None),                                          # the f block: valence_electrons is 0
    ('[Al+].[O-]S(=O)(=O)[O-]', None),                       # case 2: +3 against two anion equivalents
    ('CCN.Cl', 'CCN.Cl'),                                    # no free metal: the stage never looks
)


def _fixed(spelling):
    """Run the stage alone, against a list, and answer (molecule, written ids, log).

    Mirrors what `standardize()` does: run the stage, then recompute implicit H counts for every
    written atom via `calc_implicit` -- the same step `standardize()` takes after each stage.
    """
    from .. import calc_implicit
    mol = smiles(spelling)
    log = []
    written = _salts.fix_salt_charges(mol, log)
    for n in sorted(written):
        calc_implicit(mol, n)
    return mol, written, log


@mark.parametrize('spelling, expected', CHARGE_FIXES)
def test_fix_salt_charges_is_the_worked_table(spelling, expected):
    mol, written, log = _fixed(spelling)
    if expected is None:
        assert mol.canonical_bytes == smiles(spelling).canonical_bytes, spelling
        assert not written, spelling
        assert log and all(record.severity == REFUSED for record in log), spelling
    else:
        assert mol.canonical_bytes == smiles(expected).canonical_bytes, spelling


def test_a_refusal_writes_nothing_at_all():
    for spelling in ('[Na]', 'CC(=O)O.[Zn]', 'CC(=O)O.[Mg]'):
        mol, written, _ = _fixed(spelling)
        assert format(mol) == format(smiles(spelling)), spelling
        assert written == set(), spelling


def test_the_written_ids_are_the_atoms_whose_charge_moved():
    mol, written, _ = _fixed('CC(=O)O.[Na]')
    assert len(written) == 2
    assert {mol.atom(n).charge for n in written} == {-1, 1}


def test_the_chosen_site_and_the_moved_charge_are_named_in_the_log():
    _, _, log = _fixed('CC(=O)O.[Na]')
    assert {record.rule for record in log} == {'salts:metal-charge', 'salts:charge-transfer'}
    assert all(record.severity == REPAIRED for record in log)
    transfer = next(record for record in log if record.rule == 'salts:charge-transfer')
    assert 'salts:carboxylic-oh' in transfer.message


def test_a_metal_protic_site_is_named_in_the_log_like_any_other():
    _, _, log = _fixed('CCO.[Na]')
    transfer = next(record for record in log if record.rule == 'salts:charge-transfer')
    assert 'salts:alcohol-oh' in transfer.message and 'rung 6' in transfer.message


def test_the_alcohol_is_a_site_for_the_s_block_and_for_no_other_metal():
    """The ruling: a free group 1 or 2 metal does not stand beside a hydroxyl.  Beryllium is the one
    s-block metal outside it, and the d and p blocks are outside it as blocks."""
    for symbol in ('Li', 'Na', 'K', 'Cs', 'Mg', 'Ca', 'Ba'):
        # one alcohol per equivalent: the stage is all-or-nothing, so a group 2 metal needs two
        charge = smiles('[%s]' % symbol).atom(1).valence_electrons
        mol, written, log = _fixed('.'.join(['CCO'] * charge + ['[%s]' % symbol]))
        assert written, symbol
        assert all(record.severity == REPAIRED for record in log), symbol
        assert sum(atom.charge for atom in mol.atoms() if atom.charge > 0) == charge, symbol
    for symbol in ('Be', 'Sc', 'Al', 'Zn', 'Fe', 'Ce'):
        mol, written, log = _fixed('.'.join(['CCO'] * 4 + ['[%s]' % symbol]))
        assert not written, symbol
        assert all(record.severity == REFUSED for record in log), symbol


def test_a_thiol_is_an_acid_for_every_metal_and_not_a_metal_protic_site():
    """The sulfur rows are `protic_acid`, so they need no s-block gate -- what gates them is the same
    determinate-valence step every acid goes through."""
    mol, written, log = _fixed('CCS.CCS.CCS.[Al]')
    assert mol.canonical_bytes == smiles('CC[S-].CC[S-].CC[S-].[Al+3]').canonical_bytes
    assert all(record.severity == REPAIRED for record in log)
    # and the alcohol is not, for the same metal and the same count
    mol, written, log = _fixed('CCO.CCO.CCO.[Al]')
    assert not written and all(record.severity == REFUSED for record in log)


def test_the_ladder_picks_the_more_acidic_of_two_sites():
    # a phenol and a carboxylic acid on one molecule: rung 4 beats rung 5
    mol, written, log = _fixed('OC(=O)c1ccccc1O.[Na]')
    assert mol.canonical_bytes == smiles('[O-]C(=O)c1ccccc1O.[Na+]').canonical_bytes
    transfer = next(record for record in log if record.rule == 'salts:charge-transfer')
    assert 'salts:carboxylic-oh' in transfer.message and 'rung 4' in transfer.message


def test_the_n_h_rows_fire_and_the_azoles_refuse():
    for spelling, expected in (('O=S1(=O)NC(=O)c2ccccc21.[Na]',
                                'O=S1(=O)[N-]C(=O)c2ccccc21.[Na+]'),
                               ('O=C1NC(=O)NC1(c1ccccc1)c1ccccc1.[Na]',
                                'O=C1[N-]C(=O)NC1(c1ccccc1)c1ccccc1.[Na+]'),
                               ('[nH]1nnnc1c1ccccc1.[Na]', '[n-]1nnnc1c1ccccc1.[Na+]')):
        mol, _, _ = _fixed(spelling)
        assert mol.canonical_bytes == smiles(expected).canonical_bytes, spelling
    for spelling in ('c1cc[nH]c1.[Na]', 'c1cnc[nH]1.[Na]', 'c1cn[nH]c1.[Na]'):
        mol, written, log = _fixed(spelling)
        assert mol.canonical_bytes == smiles(spelling).canonical_bytes, spelling
        assert not written and log, spelling


def test_case_2_raises_every_metal_to_its_group_number():
    # one ionic charge per metal: two magnesiums take +2 each, and four anion equivalents want exactly that
    mol, written, _ = _fixed('O=S(=O)([O-])[O-].O=S(=O)([O-])[O-].[Mg+].[Mg+]')
    assert sorted(mol.atom(n).charge for n in written if mol.atom(n).is_metal) == [2, 2]


def test_case_2_refuses_when_the_group_numbers_fall_short():
    # sodium at +1 against two chlorides: no second charge is available to it
    mol, written, log = _fixed('[Cl-].[Cl-].[Na]')
    assert format(mol) == format(smiles('[Cl-].[Cl-].[Na]'))
    assert written == set()
    assert len(log) == 1 and log[0].severity == REFUSED
    assert '+1' in log[0].message and '2 drawn anion' in log[0].message


def test_case_2_refuses_when_the_group_numbers_overshoot():
    # aluminium's group number is +3 and sulfate draws two anion equivalents; over and under are one refusal
    mol, written, log = _fixed('[Al+].[O-]S(=O)(=O)[O-]')
    assert format(mol) == format(smiles('[Al+].[O-]S(=O)(=O)[O-]'))
    assert written == set()
    assert len(log) == 1 and log[0].severity == REFUSED
    assert '+3' in log[0].message and '2 drawn anion' in log[0].message


def test_case_2_refuses_a_count_that_is_not_a_charge():
    """The raise target is the group number, so a metal whose valence electron count states no single
    ionic charge refuses the record: zinc's count is 12, and the f block states none at all."""
    for spelling, phrase in (('[Zn+].[O-]C(=O)C(=O)[O-]', 'a count and not a charge'),
                             ('[Ce+].[O-]C(=O)C(=O)[O-]', 'not known')):
        mol, written, log = _fixed(spelling)
        assert format(mol) == format(smiles(spelling)), spelling
        assert written == set(), spelling
        assert len(log) == 1 and log[0].severity == REFUSED, spelling
        assert phrase in log[0].message, spelling


def test_a_refused_case_2_leaves_no_hydrogen_behind():
    """A refusal writes no charge, and `standardize()` derives a hydrogen count only where one moved, so
    the brutto formula and the metal's implicit count come back exactly as drawn."""
    mol = smiles('[Al+].[O-]S(=O)(=O)[O-]')
    metal = next(atom.n for atom in mol.atoms() if atom.is_metal)
    brutto, hydrogens = mol.brutto, mol.implicit_h_of(metal)
    mol.standardize()
    assert mol.brutto == brutto
    assert mol.implicit_h_of(metal) == hydrogens
    assert mol.canonical_bytes == smiles('[Al+].[O-]S(=O)(=O)[O-]').canonical_bytes


#: A free metal drawn with a hydrogen, one record per path that would have written its charge: the four
#: acidity rungs an alcohol, water, a carboxylic acid and a benzylic alcohol reach, and the case-2 raise.
#: `[H-].[Na+].CCO` is the control -- drawn ionically, the hydride is a component and nothing is at risk.
HYDRIDE_RECORDS = ('[NaH]', '[NaH].CCO', '[NaH].O', '[NaH].CC(=O)O', '[NaH].c1ccccc1CO', '[KH].CCO',
                   '[AlH2].CCO', '[MgH+].[O-]S(=O)(=O)[O-]')


@mark.parametrize('spelling', HYDRIDE_RECORDS)
def test_an_implicit_hydrogen_on_the_free_metal_refuses_the_transfer(spelling):
    """Charging a free metal drops nothing it was drawn with.

    `[NaH].CCO` is sodium hydride in ethanol AS DRAWN.  Sodium ethoxide and H2 is a reaction, not a
    repair, and the hydride has nowhere to go, so the record is left alone -- the guard `split_salts()`
    already applies before cutting a bond, on the path that writes a charge without cutting one.
    """
    mol = smiles(spelling)
    metal = next(atom.n for atom in mol.atoms() if atom.is_metal)
    brutto, hydrogens = mol.brutto, mol.implicit_h_of(metal)
    log = mol.log
    mol.standardize()
    assert mol.brutto == brutto, spelling
    assert mol.implicit_h_of(metal) == hydrogens, spelling
    assert mol.canonical_bytes == smiles(spelling).canonical_bytes, spelling
    assert 'implicit hydrogen' in log.refused()[0], spelling


def test_the_ionically_drawn_hydride_is_not_the_guard_s_business():
    """`[H-]` is a component and not a count on the metal, so the record standardizes as any other."""
    mol = smiles('[H-].[Na+].CCO')
    brutto = mol.brutto
    mol.standardize()
    assert mol.brutto == brutto
    assert mol.canonical_bytes == smiles('[H-].[Na+].CCO').canonical_bytes
    assert not mol.log.refused()


def test_case_2_requires_drawn_anions():
    mol, written, log = _fixed('[Na-]')
    assert mol.canonical_bytes == smiles('[Na-]').canonical_bytes
    assert written == set()
    assert len(log) == 1 and log[0].severity == REFUSED


def test_step_1_refusal_message_names_the_cut():
    # `[Ce]`: electrons == 0, count not known for the f block
    _, _, ce_log = _fixed('[Ce]')
    assert len(ce_log) == 1 and ce_log[0].severity == REFUSED
    # `CC(=O)O.[Zn]`: electrons == 12, a valence electron count and not a charge
    _, _, zn_log = _fixed('CC(=O)O.[Zn]')
    assert len(zn_log) == 1 and zn_log[0].severity == REFUSED
    # each message names its own cut and not the other
    assert 'not known' in ce_log[0].message and 'not known' not in zn_log[0].message
    assert 'a count' in zn_log[0].message and 'a count' not in ce_log[0].message
