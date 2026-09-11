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
"""`mol.react()`, `mol @ mol` and `mol.functional_groups()` -- the whole enumeration surface.

Every substrate is a catalogue compound, by rule: no real screening scaffold may leak into a public
repository.  A product is compared as a MOLECULE and never as its SMILES string, since `==` is the
canonical form and a string would pin the writer's atom order.  The ORDER rows come out in is table
order and is pinned nowhere.
"""
import pytest
from .._enumerate import EnumeratedReaction, functional_group_hits, functional_groups
from .._tables import functional_rules
from ...core import H_UNKNOWN, read_smiles as smiles


def products(enumerated):
    """`{name: {product molecule}}` for one enumeration, so an assertion names a reaction and a product.

    Every row in this corpus makes one product molecule per outcome, byproducts being deleted rather
    than emitted; a row that ever emits two fails here rather than passing a half-assertion.
    """
    out = {}
    for row in enumerated:
        assert isinstance(row, EnumeratedReaction)
        assert len(row.reaction.products) == 1, f'{row.name} made {len(row.reaction.products)}'
        out.setdefault(row.name, set()).add(row.reaction.products[0])
    return out


def skeletons(molecules):
    """The same molecules with every implicit hydrogen count erased to `H_UNKNOWN`.

    A patch reaching an aromatic atom whose class only the ring decides leaves that count `H_UNKNOWN`;
    the repair is `kekule()` then `calc_implicit()`, and `calc_implicit` lives in the sibling package
    `chython.chemistry`, which this one may not import.  `chython/chemistry/test/
    test_reaction_hydrogen_repair.py` tests it end to end.  So a coupling is compared as a skeleton --
    erased on both sides, symmetric rather than lenient -- while every aliphatic outcome here still uses
    plain `==`, hydrogens included.
    """
    out = set()
    for molecule in molecules:
        molecule = molecule.copy()
        for n in molecule.atom_numbers:
            molecule.set_hydrogens(n, H_UNKNOWN)
        out.add(molecule)
    return out


# --- functional_groups ---------------------------------------------------------------------------

def test_a_group_absent_from_the_molecule_is_absent_from_the_dict():
    """Which is what makes `name in mol.functional_groups()` the presence test, with no zero to skip.

    Stated as two properties rather than as acetic acid's whole dict: the corpus grows, and a row added
    to `functional.tsv` must not be able to fail a test about the shape of the answer.  Every row's own
    reach is `test_functional.py`'s business.
    """
    found = smiles('CC(=O)O').functional_groups()
    assert 'carboxylic_acid' in found and 'ester' not in found
    assert all(count > 0 for count in found.values()), found


def test_the_count_is_distinct_matches_and_not_automorphic_ones():
    """Terephthalic acid has two acid groups.  Four would mean the symmetry was counted twice."""
    assert smiles('OC(=O)c1ccc(C(=O)O)cc1').functional_groups()['carboxylic_acid'] == 2


@pytest.mark.parametrize('name, smi, sites', [
    # A pattern that writes out equivalent atoms admits their permutations as separate mappings, so
    # counting mappings reports a factorial: one CF3 is 3! and one sulfonyl is 2!.
    ('trifluoromethyl', 'OCC(F)(F)F', 1),
    ('difluoromethyl', 'FC(F)c1ccccc1', 1),
    ('sulfonamide', 'NS(=O)(=O)c1ccccc1', 1),
    ('sulfone', 'CS(=O)(=O)C', 1),
    ('secondary_amine', 'CCNCC', 1),
    ('1_3_diketone', 'CC(=O)CC(=O)C', 1),
    # and the other side of the rule: these really are two sites, and stay 2
    ('secondary_amine', 'C1CNCCN1', 2),
    ('vicinal_diol', 'OCC(O)CO', 2),
    ('carboxylic_acid', 'OC(=O)c1ccc(C(=O)O)cc1', 2),
    ('arene_ch', 'c1ccccc1', 6),
])
def test_the_count_is_sites_and_not_mappings(name, smi, sites):
    """How many sets of atoms the pattern covers -- not how many ways it can be laid over them.

    The distinction is invisible to `name in mol.functional_groups()` and decides every count filter
    built on the answer, which is why it gets a table of both directions rather than one example.
    """
    assert smiles(smi).functional_groups()[name] == sites


@pytest.mark.parametrize('smi, name', [
    ('NC=O', 'primary_amide'),          # formamide
    ('CNC=O', 'secondary_amide'),       # N-methylformamide
    ('NC=S', 'thioamide'),              # thioformamide
    ('CC(=O)N', 'primary_amide'),       # and the D3 carbonyl the rows always reached
    ('CNC(=O)C', 'secondary_amide'),
    ('CC(N)=S', 'thioamide'),
])
def test_a_formamide_carbonyl_is_reached(smi, name):
    """Hydrogens do not count toward `D`, so a formamide's carbonyl carbon is `D2` and not `D3`.

    The three rows wrote `D3`, which reported formamide, N-methylformamide and thioformamide as carrying
    no group whatsoever -- the failure mode a `D` on a carbon that may bear a hydrogen always has.
    """
    assert name in smiles(smi).functional_groups()


def test_one_molecule_can_carry_several_groups():
    found = smiles('OCc1ccccc1Br').functional_groups()
    assert found['primary_alcohol'] == 1
    assert found['aryl_bromide'] == 1


def test_the_method_is_not_cached_across_an_edit():
    """The enumerators read this, so a stale count silently changes which templates are tried."""
    molecule = smiles('CCO')
    assert molecule.functional_groups() == {'primary_alcohol': 1}

    # propane: the alcohol is gone, and a cache would still be reporting it
    oxygen, = molecule.atoms_of_element(8)
    with molecule.edit() as e:
        e.set_element(oxygen, 'C')
    assert molecule.functional_groups() == {}

    # and back, because a cache invalidated once is not a cache invalidated
    with molecule.edit() as e:
        e.set_element(oxygen, 'O')
    assert molecule.functional_groups() == {'primary_alcohol': 1}


# --- the same answer, with the row ids -----------------------------------------------------------

def test_a_hit_carries_the_row_id_the_name_and_the_count():
    """An id is storage identity, so it comes back with the hit rather than through a second table
    lookup a caller has to know to make."""
    hits = functional_group_hits(smiles('OC(=O)c1ccc(C(=O)O)cc1'))
    by_name = {hit.name: hit for hit in hits}
    assert by_name['carboxylic_acid'].count == 2
    assert by_name['carboxylic_acid'].id == functional_rules()['carboxylic_acid'].id
    assert by_name['carboxylic_acid'].id.startswith('functional:')


def test_the_dict_form_is_the_same_answer_folded():
    mol = smiles('OCc1ccccc1Br')
    assert functional_groups(mol) == {hit.name: hit.count for hit in functional_group_hits(mol)}


def test_hits_come_in_table_order_so_a_stored_id_array_is_stable():
    """What lets a consumer store the ids as a sorted array and compare two molecules by set arithmetic
    rather than by dict merge."""
    order = [group.id for group in functional_rules().values()]
    hits = [hit.id for hit in functional_group_hits(smiles('OCc1ccccc1Br'))]
    assert hits == [i for i in order if i in set(hits)]


def test_the_container_method_is_the_function():
    mol = smiles('CC(=O)Nc1ccc(O)cc1')                      # paracetamol
    assert mol.functional_group_hits() == functional_group_hits(mol)


# --- react / @ -----------------------------------------------------------------------------------

def test_amidation():
    found = products(smiles('CC(=O)O') @ smiles('CCN'))
    assert found == {'amidation': {smiles('CC(=O)NCC')}}


def test_the_argument_order_does_not_decide_the_slot_order():
    """A row states chemical roles; the caller states what is in the flask.

    Nothing reorders the caller's molecules: both go to the matcher at once and the reactant side finds
    the acid where the acid is.
    """
    forward = products(smiles('CC(=O)O') @ smiles('CCN'))
    backward = products(smiles('CCN') @ smiles('CC(=O)O'))
    assert forward == backward == {'amidation': {smiles('CC(=O)NCC')}}


def test_a_mixture_in_one_container_is_the_same_reaction():
    """Two components of one input, rather than two inputs.

    The `(A).(B)` grouping asks about COMPONENTS, so a two-slot row applies to a single container
    holding a reagent mixture -- an ordinary SDF record.
    """
    mixture = smiles('CC(=O)O.CCN')
    assert products(mixture.react(reaction='amidation')) == {'amidation': {smiles('CC(=O)NCC')}}


def test_suzuki():
    """Both leaving groups go by absence: the bromine, and the boron with its two hydroxyls."""
    found = products(smiles('Brc1ccccc1') @ smiles('OB(O)c1ccc(C)cc1'))
    assert set(found) == {'suzuki'}
    assert skeletons(found['suzuki']) == skeletons([smiles('Cc1ccc(-c2ccccc2)cc1')])


def test_one_row_serves_two_sites_of_one_molecule():
    """4-bromoiodobenzene has two coupling sites, so ONE suzuki row matches twice and both couple.

    `aryl_halide` is `[Cl,Br,I;D1]`, so the two outcomes are two MATCHES of one row and not one match
    each of two -- which is what the `rule_id` assertion pins.
    """
    outcomes = list(smiles('Brc1ccc(I)cc1') @ smiles('OB(O)c1ccccc1'))
    assert {row.rule_id for row in outcomes} == {'reactions:8'}, 'one row, matched twice'
    found = products(outcomes)
    assert skeletons(found['suzuki']) == skeletons([smiles('Brc1ccc(-c2ccccc2)cc1'),
                                                   smiles('Ic1ccc(-c2ccccc2)cc1')])


def test_every_input_has_to_be_touched():
    """Three molecules and a two-component row: refused, because the answer would ignore an input.

    The test is on the OUTCOME: a template reports as `reactants` exactly the inputs its match touched,
    so an outcome naming two of three inputs is dropped -- silently, being a statement about the corpus
    and not about the caller.
    """
    acid, amine, spectator = smiles('CC(=O)O'), smiles('CCN'), smiles('Cc1ccccc1')
    assert products(acid.react(amine, spectator)) == {}
    assert products(acid.react(amine)) == {'amidation': {smiles('CC(=O)NCC')}}
    # and the spectator is refused for being untouched, not for being unreactive: an alcohol that
    # WOULD react on its own with nobody to react with is refused just the same
    assert products(acid.react(amine, smiles('CCO'))) == {}


def test_an_untouched_input_is_refused_but_an_untouched_component_is_not():
    """The distinction the every-input rule turns on, and it is easy to conflate.

    An untouched INPUT means the caller asked about something the answer ignores.  An untouched
    COMPONENT of a touched input comes out -- the salt rule, which keeps a counter-ion from vanishing.
    """
    salted = smiles('CC(=O)O.[Na+].[Cl-]')
    outcomes = list(salted.react(smiles('CCN'), reaction='amidation'))
    assert len(outcomes) == 1
    assert len(outcomes[0].reaction.reactants) == 2, 'two inputs, both touched'
    assert len(outcomes[0].reaction.products) == 3, 'the sodium and the chloride came through'


# --- intramolecular ------------------------------------------------------------------------------

def test_an_intramolecular_reaction_needs_no_second_molecule():
    """A lactam out of one amino acid, from the SAME ROW as the intermolecular amidation."""
    found = products(smiles('NCCCCCC(=O)O').react(reaction='amidation'))
    assert found == {'amidation': {smiles('O=C1CCCCCN1')}}, 'epsilon-caprolactam'


def test_the_ring_size_is_what_makes_it_a_reaction():
    """5, 6 and 7 close; 4 and 12 do not.  The product-side `r` is doing this and nothing else is.

    Which is why the intramolecular template is a separate composition rather than an unconstrained `.`
    join: `.` says only "not bonded", so it fires on a strained four-ring as happily as on a six.
    """
    closes = {4: 'NCCC(=O)O', 5: 'NCCCC(=O)O', 6: 'NCCCCC(=O)O', 7: 'NCCCCCC(=O)O',
              12: 'NCCCCCCCCCCC(=O)O'}
    fired = {size for size, s in closes.items() if products(smiles(s).react(reaction='amidation'))}
    assert fired == {5, 6, 7}, 'reactions:1 names 1:5,6,7'


def test_the_two_readings_never_both_fire():
    """The two groupings are complementary, so a row carries both templates with no deduplication."""
    together = list(smiles('NCCCCCC(=O)O').react(reaction='amidation'))
    apart = list(smiles('CC(=O)O').react(smiles('CCN'), reaction='amidation'))
    assert len(together) == 1 and len(apart) == 1


def test_a_row_with_no_ring_sizes_has_no_intramolecular_reading():
    """An empty `ring_sizes` cell says the corpus does not know which sizes close, not that none do."""
    tethered = smiles('Brc1ccccc1CCc1ccccc1B(O)O')
    assert 'aryl_bromide' in functional_groups(tethered)
    assert 'aryl_boronic_acid' in functional_groups(tethered)
    assert products(tethered.react(reaction='suzuki')) == {}


def test_selecting_one_reaction():
    both = products(smiles('Brc1ccccc1') @ smiles('NCC'))
    assert 'buchwald_hartwig' in both
    only = products(smiles('Brc1ccccc1').react(smiles('NCC'), reaction='buchwald_hartwig'))
    assert set(only) == {'buchwald_hartwig'}


def test_an_unknown_reaction_name_is_refused():
    """"No such reaction" and "does not apply here" are the same empty generator; only one is a typo."""
    with pytest.raises(ValueError) as exc:
        list(smiles('CC(=O)O').react(smiles('CCN'), reaction='suzukii'))
    assert 'suzukii' in str(exc.value)


def test_a_partner_with_no_matching_group_enumerates_nothing():
    assert products(smiles('CCCC') @ smiles('CCCC')) == {}


def test_an_outcome_names_the_row_it_came_from():
    """`name` is the chemistry, `rule_id` the spelling: only the id says which row earned the hit."""
    outcomes = list(smiles('CC(=O)O').react(smiles('CCN'), reaction='amidation'))
    assert {row.rule_id for row in outcomes} == {'reactions:1'}


# --- the single-molecule rows, through the one method --------------------------------------------
#
# `react()` with no partner: one table and one method, where the corpus once had three of each.

def test_oxidation_of_a_primary_alcohol():
    found = products(smiles('CCO').react(reaction='alcohol_to_aldehyde'))
    assert found == {'alcohol_to_aldehyde': {smiles('CC=O')}}


def test_reduction_of_a_ketone():
    found = products(smiles('CC(=O)C').react(reaction='ketone_to_alcohol'))
    assert found == {'ketone_to_alcohol': {smiles('CC(O)C')}}


def test_reduction_of_a_nitroarene_deletes_both_oxygens():
    """`[A:1]` alone: the nitrogen inherits its element, the oxygens are gone by absence.

    Also the narrowest test that a product side is EXPLICIT-ONLY about charge: the reactant states `+`
    on the nitrogen, the product does not, so the product nitrogen is neutral.
    """
    found = products(smiles('[O-][N+](=O)c1ccccc1').react(reaction='nitro_to_amine'))
    assert found == {'nitro_to_amine': {smiles('Nc1ccccc1')}}


def test_an_oxidation_keeps_the_substituents_it_matched():
    """The N-oxide, not `[NH3+][O-]`.

    `tertiary_amine` numbers its three substituents so a template can address them, which makes
    restating them the product side's job -- deletion being by absence.  The row that did not restate
    them returned a molecule of two atoms, and the only test naming this reaction asserted its name.
    """
    found = products(smiles('CN(C)C').react(reaction='nitrogen_oxidation'))
    assert found == {'nitrogen_oxidation': {smiles('C[N+](C)(C)[O-]')}}


@pytest.mark.parametrize('smi, what', [
    ('CN(C)C=O', 'a tertiary amide'),
    ('O=C(OC(C)(C)C)N1CCCCC1', 'an N-Boc amine'),
    ('CN(C)S(=O)(=O)C', 'a sulfonamide'),
])
def test_nothing_but_an_amine_is_offered_for_n_oxidation(smi, what):
    """`tertiary_amine` states its three substituents, so DMF is not a tertiary amine.

    The bare `[N;D3;z1;x0]` matched every tertiary amide and every carbamate, which reported an N-Boc
    amine as an amine and offered to oxidize DMF.
    """
    assert 'tertiary_amine' not in smiles(smi).functional_groups(), what
    assert not products(smiles(smi).react(reaction='nitrogen_oxidation'))


def test_a_transformation_creates_the_atom_it_needs():
    found = products(smiles('CCO').react())
    assert found['appel'] == {smiles('CCBr')}
    assert found['appel_chloride'] == {smiles('CCCl')}


def test_the_transition_state_of_an_outcome_holds_what_left_and_what_arrived():
    """The reactor numbers a pair and leaves a leaving or an arriving atom at 0, so the ML view of an
    outcome has to place an unmapped atom rather than count it.

    Buchwald-Hartwig on aziridine and bromobenzene.  The aryl carbon has three heavy neighbours before
    and three after -- two ring carbons and the bromine, then two ring carbons and the nitrogen -- and
    the C-Br bond has to be in the union for that to be true.  Without it, every aryl halide of one
    ring gives the same transition state.
    """
    out = next(iter(smiles('N1CC1') @ smiles('c1ccccc1Br')))
    view = out.reaction.modeling_view()
    assert view.unmapped == {'reactants': 1, 'products': 0}
    aryl = next(n for n, state in view.states.items() if state[:3] == (6, 0, 3))
    assert view.states[aryl][3:] == (0, 3), 'three heavy neighbours after as well'
    leaving = [n for n, state in view.states.items() if state[0] == 35]
    assert len(leaving) == 1, 'the bromine that left is one row of the union'
    assert view.union_bonds[tuple(sorted((aryl, leaving[0])))] == (1, 0), 'C-Br broken'


def test_a_created_atom_reaches_the_transition_state_too():
    """The mirror: pyridine N-oxidation creates the oxygen, so the reactor leaves it at 0.

    Dropped, this record's transition state has no bond change at all -- an oxidation that looks inert.
    """
    out = next(iter(smiles('c1ccccn1').react(reaction='nitrogen_oxidation')))
    view = out.reaction.modeling_view()
    assert view.unmapped == {'reactants': 0, 'products': 1}
    arriving = [n for n, state in view.states.items() if state[0] == 8]
    assert len(arriving) == 1, 'the oxide oxygen is one row of the union'
    assert [orders for orders in view.union_bonds.values() if orders[0] != orders[1]] == [(0, 1)]


def test_a_ring_bond_the_product_omits_is_deleted():
    """Epoxide hydrolysis: the product names three atoms in a chain, so the closing bond goes.

    Deletion by absence applies to a BOND and not only to an atom.
    """
    found = products(smiles('C1OC1c1ccccc1').react())
    assert found['epoxide_opening'] == {smiles('OCC(O)c1ccccc1')}


def test_no_partner_reaches_every_kind_of_single_molecule_row():
    """One call, and an oxidation, an interconversion and a reduction all come back out of it.

    Nothing in the answer says which was which; a caller who wants one asks for it by name.
    """
    found = products(smiles('CCO').react())
    assert {'alcohol_to_aldehyde', 'appel', 'appel_chloride'} <= set(found)

    # and a reduction from the same one method, on a substrate that has one
    assert 'ketone_to_alcohol' in products(smiles('CC(=O)C').react())


def test_a_named_reaction_is_the_only_scope_there_is():
    """`reaction=` scopes to one row family, which is finer than a whole heading of the corpus.

    What is genuinely gone is "every oxidation of this molecule" as a single question: the taxonomy was
    a filing decision nothing verified, so it is a `#` banner in the TSV and not a column.
    """
    alcohol = smiles('CCO')
    scoped = set(products(alcohol.react(reaction='alcohol_to_aldehyde')))
    assert scoped < set(products(alcohol.react())), 'a scope, not the whole corpus'
    assert scoped == {'alcohol_to_aldehyde'}


def test_a_call_with_partners_never_yields_a_one_slot_row():
    """A one-slot row touches one input, and every input must be touched -- so `acid.react(amine)` never
    reports the acid's own chlorination beside the amide."""
    paired = set(products(smiles('Brc1ccccc1').react(smiles('NCC'))))
    assert paired and not paired & set(products(smiles('Brc1ccccc1').react()))


def test_a_substrate_with_no_group_enumerates_nothing():
    assert products(smiles('CCCC').react()) == {}


# --- the prefilter is a prefilter ----------------------------------------------------------------

def test_present_groups_do_not_guarantee_a_product():
    """Both groups present, in one molecule, with an intramolecular template composed -- and no product.

    Glycine satisfies the multiset prefilter and `reactions:1` has an intramolecular template to try;
    the ring size refuses it (only a four-ring closes, the row names 5, 6, 7).  The enumeration is
    decided by the MATCH, so a non-empty `functional_groups()` is never a promise.
    """
    glycine = smiles('NCC(=O)O')
    present = functional_groups(glycine)
    assert 'carboxylic_acid' in present and 'primary_amine' in present
    assert products(glycine.react(reaction='amidation')) == {}


def test_the_prefilter_is_a_multiset():
    """A row naming one group twice needs two matches of it, from wherever in the inputs they come.

    Written against `_run` with a fixture rule, the corpus having no two-of-a-kind row yet.
    """
    from .._enumerate import _run
    from .._tables import ReactionRule, compose_smirks
    from ...core import read_smirks

    doubled = ('carboxylic_acid', 'carboxylic_acid')
    smirks = compose_smirks(doubled, '[A:1](=[A:2])-[A:101](=[A:102])')
    rule = ReactionRule('test:1', 'anhydride', doubled, '[A:1](=[A:2])-[A:101](=[A:102])',
                        read_smirks(smirks, rule_id='test:1'), 'two acids')
    corpus = {rule.name: (rule,)}
    assert not list(_run((smiles('CC(=O)O'),), corpus))
    assert list(_run((smiles('CC(=O)O'), smiles('CCC(=O)O')), corpus))


def test_a_carried_racemic_centre_stays_one_group():
    # Reported as `|&1:4|` in, `|&2:14|` out.  Needs the reactor's mapping: without it the writer has
    # nothing to merge on.
    mol = smiles('C[C@H](CCN)CC(O)=O |&1:1,r|')
    rxn = next(o.reaction for o in mol.react(reaction='amidation'))
    assert str(rxn) == 'O=C(O)C[C@@H](CCN)C>>O=C1NCC[C@@H](C)C1 |&1:4,14|'
