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
refusals matter as much as the splits.  `decompose_salts` reads the same table the other way round: a
tabulated species is a counterion only when something else is there to be the compound, so the record
that is nothing but a salt former answers with that former as its own parent.
"""
from pytest import raises
# `__all__` and not the package object: `from ... import chemistry` would execute the facade, which
# `test_dependency_direction.py` ratchets against.
from .. import SALT_ROLES, __all__ as CHEMISTRY_ALL, decompose_salts, split_salts
from .. import _salts
from .._tables import salts_rows, salts_rows_by_role, salts_species_keys, salts_table_text
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


# the keep surface: row id, role or element symbol


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


def test_a_role_name_reaches_every_row_of_that_role():
    """A role name saves a caller from naming rows one by one and going stale each time the table
    grows.  Only the two SMARTS roles mean anything here: nothing else names an atom."""
    for keep in [['cation'], ['acceptor']]:
        m = smiles('CC(=O)O[Na]')
        assert split_salts(m, keep=keep) is False, keep
    # a role that names no atom leaves the split alone
    m = smiles('CC(=O)O[Na]')
    assert split_salts(m, keep=['solvate']) is True


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
    assert "keep=['cation']" in str(e.value)


def test_a_keep_item_of_the_wrong_kind_is_refused_at_the_boundary():
    with raises(TypeError) as e:
        split_salts(smiles('CC(=O)O[Na]'), keep=[11])
    assert 'row ids, roles or element symbols' in str(e.value)

    with raises(ValueError) as e:
        split_salts(smiles('CC(=O)O[Na]'), keep=['Unobtainium'])
    assert 'a row id' in str(e.value) and 'a role' in str(e.value)


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


# decompose_salts: the relational order


def test_a_record_that_is_its_own_salt_former_is_its_own_parent():
    """Acetic acid is acetic acid.  Nothing is beside it, so nothing is counted beside it."""
    r = smiles('CC(=O)O').decompose_salts()
    assert [str(p) for p in r.parents] == ['C(C)(=O)O']
    assert (r.counterions, r.solvates, r.cations) == ({}, {}, {})

    r = smiles('Cl').decompose_salts()
    assert [str(p) for p in r.parents] == ['Cl']
    assert (r.counterions, r.solvates, r.cations) == ({}, {}, {})


def test_a_salt_reports_the_acid_as_the_parent_and_the_metal_as_a_cation():
    for s in ['CC(=O)O[Na]', 'CC(=O)[O-].[Na+]', 'CC(=O)O.[Na+]']:
        r = smiles(s).decompose_salts()
        assert [str(p) for p in r.parents] == ['C(C)(=O)O'], s
        assert r.cations == {'Na': 1} and r.counterions == {}, s
    # and sodium chloride is hydrochloric acid beside one sodium, rather than no compound at all
    r = smiles('[Na+].[Cl-]').decompose_salts()
    assert [str(p) for p in r.parents] == ['Cl'] and r.cations == {'Na': 1}


def test_a_compound_beside_a_former_reports_the_compound():
    r = smiles('NCC(=O)O.OC(=O)C(F)(F)F.O').decompose_salts()
    assert [str(p) for p in r.parents] == ['C(CN)(=O)O']
    assert r.counterions == {'salts:tfa': 1} and r.solvates == {'salts:water': 1}


def test_a_solvate_is_counted_even_when_every_component_is_tabulated():
    """The middle rung: acetic acid monohydrate is acetic acid, with the water counted."""
    r = smiles('CC(=O)O.O').decompose_salts()
    assert [str(p) for p in r.parents] == ['C(C)(=O)O'] and r.solvates == {'salts:water': 1}


def test_a_record_of_nothing_but_solvates_is_the_solvate():
    assert [str(p) for p in smiles('O.O').decompose_salts().parents] == ['O']
    assert [str(p) for p in smiles('CCO.O').decompose_salts().parents] == ['C(C)O', 'O']


def test_two_drawn_equivalents_of_one_compound_are_one_parent():
    r = smiles('OC(=O)c1ccccc1O.OC(=O)c1ccccc1O').decompose_salts()
    assert len(r.parents) == 1


def test_two_enantiomers_drawn_together_are_two_parents():
    """Parents dedup by canonical bytes, not by the constitution key: a racemate drawn as two
    components is two compounds, and collapsing them would report a stoichiometry that was not drawn."""
    r = smiles('N[C@@H](C)C(=O)O.N[C@H](C)C(=O)O').decompose_salts()
    assert len(r.parents) == 2


def test_decompose_salts_changes_nothing_and_logs_nothing():
    for s in ('NCC(=O)O.OC(=O)C(F)(F)F', 'CC(=O)O[Na].O', 'N[Pt](N)(Cl)Cl', 'Cl'):
        m = smiles(s)
        log, before = m.log, bytes(m.canonical_bytes)
        m.decompose_salts()
        assert bytes(m.canonical_bytes) == before, s
        assert not len(log), s


def test_the_salt_surface_is_two_methods():
    """One that edits and one that reports.  The third, which deleted components, asked the same
    question as the reporter and answered it destructively; a consumer takes `parents[0]` instead."""
    mol = smiles('CC(=O)O[Na].O')
    assert callable(mol.split_salts) and callable(mol.decompose_salts)
    for gone in ('split_ionic', 'strip_salts', 'salt_composition'):
        assert not hasattr(mol, gone), f'{gone} still resolves'
        assert gone not in CHEMISTRY_ALL, f'{gone} is still on the package'
    assert set(_salts.__all__) == {'SaltComposition', 'decompose_salts', 'split_salts'}


#  smiles                                          parent smiles, counterions, solvates, cations
COMPOSITIONS = [
    ('NCC(=O)O',                                    'NCC(=O)O', {}, {}, {}),
    ('NCC(=O)O.OC(=O)C(F)(F)F',                     'NCC(=O)O', {'salts:tfa': 1}, {}, {}),
    ('[NH3+]CC(=O)O.[O-]C(=O)C(F)(F)F',             'NCC(=O)O', {'salts:tfa': 1}, {}, {}),
    ('NCC(=O)O.OC(=O)C(F)(F)F.OC(=O)C(F)(F)F',      'NCC(=O)O', {'salts:tfa': 2}, {}, {}),
    ('NCCCCN.Cl.Cl.O',                              'NCCCCN', {'salts:hcl': 2},
     {'salts:water': 1}, {}),
    # a record whose charges do not balance is still read: the counterion is named and the parent comes
    # back neutral, since `keep_charge=False` takes each component as close to neutral as it goes.
    ('CN(C)CCOC(c1ccccc1)c1ccccc1.[O-]S(=O)(=O)c1ccc(C)cc1',
     'CN(C)CCOC(c1ccccc1)c1ccccc1', {'salts:tosylic': 1}, {}, {}),
    ('[NH2]CC(=O)[O-].[K+]',                        'NCC(=O)O', {}, {}, {'K': 1}),
    # a quaternary ammonium has no proton to give, so it keeps its charge and the chloride still counts
    ('C[N+](C)(C)C.[Cl-]',                          'C[N+](C)(C)C', {'salts:hcl': 1}, {}, {}),
    ('C[N+](C)(C)CC(=O)[O-]',                       'C[N+](C)(C)CC(=O)[O-]', {}, {}, {}),
]


def test_the_composition_of_a_record():
    for s, parent, counterions, solvates, cations in COMPOSITIONS:
        r = smiles(s).decompose_salts()
        assert len(r.parents) == 1, s
        assert r.parents[0].canonical_bytes == smiles(parent).canonical_bytes, s
        assert r.counterions == counterions, s
        assert r.solvates == solvates, s
        assert r.cations == cations, s


def test_an_explicit_hydrogen_does_not_hide_a_counterion():
    """The key is taken from a molecule with implicit hydrogens, where an explicit one is a difference;
    the implicification inside the pass is what keeps a hydrogen-atom drawing readable."""
    a = smiles('NCC(=O)O.[H]OC(=O)C(F)(F)F').decompose_salts()
    b = smiles('NCC(=O)O.OC(=O)C(F)(F)F').decompose_salts()
    assert a.counterions == b.counterions == {'salts:tfa': 1}


def test_a_kekule_drawing_of_a_solvate_still_matches():
    """`thiele()` runs on the copy, so a Kekule toluene keys to the tabulated aromatic one."""
    for s in ['NCC(=O)O.Cc1ccccc1', 'NCC(=O)O.C1=CC=CC=C1C']:
        assert smiles(s).decompose_salts().solvates == {'salts:toluene': 1}, s


def test_a_counterion_row_names_a_constitution_and_not_a_stereoisomer():
    """One tartrate row covers L, D, meso and undefined, the key being taken with stereo disabled."""
    for s in ['CN.O[C@H]([C@@H](O)C(O)=O)C(O)=O',
              'CN.O[C@@H]([C@H](O)C(O)=O)C(O)=O',
              'CN.OC(C(O)C(O)=O)C(O)=O']:
        r = smiles(s).decompose_salts()
        assert [str(p) for p in r.parents] == ['CN'], s
        assert r.counterions == {'salts:tartaric': 1}, s


def test_the_compound_of_interest_keeps_its_own_stereo():
    """Stereo is dropped from the KEY and not from the molecule, so a parent comes back configured."""
    r = smiles('C[C@H](N)C(=O)O.Cl').decompose_salts()
    assert r.counterions == {'salts:hcl': 1}
    assert r.parents[0].canonical_bytes == smiles('C[C@H](N)C(=O)O').canonical_bytes
    assert r.parents[0].canonical_bytes != smiles('C[C@@H](N)C(=O)O').canonical_bytes


def test_a_labelled_solvate_is_a_different_species():
    """No `clean_isotopes()`, matching how the table itself was loaded: D2O is not water, so it is
    reported as an unrecognized component rather than silently counted as a hydrate."""
    r = smiles('NCC(=O)O.[2H]O[2H]').decompose_salts()
    assert r.solvates == {}
    assert len(r.parents) == 2


def test_a_coordination_complex_is_a_parent_and_not_a_composition():
    """The all-or-nothing test inside `split_salts` reaches here: a dative bond exempts the whole atom,
    so cisplatin is one compound rather than a platinum and two chlorides."""
    r = smiles('N[Pt](N)(Cl)Cl').decompose_salts()
    assert len(r.parents) == 1
    assert r.counterions == {} and r.cations == {}


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


def test_a_named_base_is_a_base_and_not_a_counterion():
    by_role = salts_rows_by_role()
    assert {r.id for r in by_role['base']} >= {'salts:ammonia', 'salts:tromethamine', 'salts:choline'}
    assert 'salts:ammonia' not in {r.id for r in by_role['counterion']}


def test_a_base_counts_where_a_counterion_does():
    """The two roles differ in which side of the salt a species came from, not in being beside the
    compound, so `counterions` holds both and a consumer reads one field."""
    r = smiles('OC(=O)c1ccc(cc1)C(=O)Nc1ccccc1.OCC(N)(CO)CO').decompose_salts()
    assert r.counterions == {'salts:tromethamine': 1}


def test_a_species_that_is_both_solvent_and_base_is_tabulated_as_the_base():
    """Pyridine under `solvate` would read pyridine hydrochloride as hydrochloric acid, the species rung
    beating the solvate one.  Under `base` the record answers with both candidates, which is visible."""
    assert salts_species_keys()[format(smiles('c1ccncc1'), '!s')].role == 'base'
    r = smiles('c1ccncc1.Cl').decompose_salts()
    assert {str(p) for p in r.parents} == {'c1ccccn1', 'Cl'}
    assert r.counterions == {} and r.solvates == {}


def test_no_two_rows_share_a_key():
    assert len(salts_species_keys()) == sum(len(salts_rows_by_role()[r])
                                            for r in ('counterion', 'base', 'solvate'))


def test_the_ionic_conjugates_are_gone_because_neutralize_reaches_them():
    ids = {row.id for row in salts_rows()}
    assert not ids & {'salts:acetate', 'salts:chloride', 'salts:tosylate', 'salts:ammonium'}
    assert 'salts:hf' in ids                       # fluoride's neutral twin, which was missing


def test_a_row_no_longer_carries_a_keep_flag():
    assert not hasattr(salts_rows()[0], 'keep')
    for line in salts_table_text().splitlines():
        if line.startswith('id\t'):
            assert line.split('\t') == ['id', 'role', 'pattern', 'charges', 'comment']
            break
    else:
        raise AssertionError('salts.tsv has no column header')


def test_every_row_is_internally_consistent():
    """The loaded row shape the passes index into.  Exactly one of `query`/`key` per row is what makes
    `role` load-bearing; load-time checks themselves live in `_tables.py`."""
    ids = set()
    for row in salts_rows():
        assert row.id not in ids, f'{row.id} appears twice'
        ids.add(row.id)
        assert row.id.startswith('salts:'), f'{row.id} is not table-qualified'
        assert row.role in SALT_ROLES
        assert (row.query is None) != (row.key is None), \
            f'{row.id} has both a query and a key, or neither'
        if row.query is not None:
            assert row.anchor in row.query.map_numbers()
        assert bool(row.charges) == (row.role == 'cation'), \
            f'{row.id}: charges are the cation overcharge guard and mean nothing elsewhere'
        assert row.comment, f'{row.id} has no comment'


def test_the_metal_row_covers_every_metal_the_core_calls_one():
    """93 metals: `[M]`'s membership is the core's business and this row inherits it, so a change
    there shows up here as a diff."""
    row, = salts_rows_by_role()['cation']
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
