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
"""`protective.tsv`, `protective_groups()` and `deprotect()`.

Every row carries its own acceptance test: the `protected` and `cleaved` columns are parametrized over
below, so a row cannot enter the table without a worked example nor stay once the example stops holding.
What no test can check is whether the row is right about chemistry -- whether a reagent cleaving its group
exists and leaves the rest alone.  That is chemist review, and the `decoys` column is what makes it
reviewable.
"""
import pytest
from itertools import combinations
from .._enumerate import _claims, deprotect, protective_group_hits, protective_groups
from .._tables import PROTECTS_ELEMENTS, protective_rules, read_table
from ...core import H_UNKNOWN, read_smiles


BY_NAME = protective_rules()
RULES = tuple(BY_NAME.values())


def skeleton(molecule):
    """The molecule with every implicit hydrogen count erased.

    A patch reaching an aromatic atom whose class only the ring decides leaves that count `H_UNKNOWN`, and
    the repair needs `chython.chemistry.calc_implicit` -- a sibling package this one may not import.  So
    outcomes are compared with the count erased on both sides: symmetric, and blind to hydrogen
    bookkeeping and nothing else.  `chython/chemistry/test/` tests the repair end to end.
    """
    molecule = molecule.copy()
    for n in molecule.atom_numbers:
        molecule.set_hydrogens(n, H_UNKNOWN)
    return molecule


def stripped(molecule, *names, **kwargs):
    """The products of the first deprotection of `molecule`, as a set of skeletons."""
    outcome = next(deprotect(molecule, names, **kwargs), None)
    if outcome is None:
        return None
    return {skeleton(p) for p in outcome.reaction.products}


# --- every row, against its own worked example ----------------------------------------------------

@pytest.mark.parametrize('rule', RULES, ids=lambda rule: rule.name)
def test_every_row_cleaves_its_own_example(rule):
    """`protected` -> `cleaved`, through the row's own template.

    Applied directly rather than through `deprotect()`, so a failure here is the row's and not the claim
    walk's: the pattern matches, the product side keeps the atom it meant to, and it restated every bond
    between kept atoms -- a forgotten one hands back fragments instead of a molecule.
    """
    molecule = read_smiles(rule.protected)
    want = skeleton(read_smiles(rule.cleaved))
    outcomes = list(rule.template(molecule))
    assert outcomes, f'{rule.id} does not match its own `protected` example {rule.protected!r}'
    got = [{skeleton(p) for p in rxn.products} for rxn in outcomes]
    assert any(want in products for products in got), (
        f'{rule.id} turns {rule.protected!r} into '
        + ' | '.join('.'.join(str(p) for p in rxn.products) for rxn in outcomes)
        + f', not into {rule.cleaved!r}')


@pytest.mark.parametrize('rule', [r for r in RULES if r.decoys], ids=lambda rule: rule.name)
def test_no_row_matches_its_decoys(rule):
    """A decoy is a thing this pattern once caught and should not: a widening that went too far, recorded
    next to the pattern it went too far in, so the next widening has a fixed point to test against."""
    for decoy in rule.decoys:
        molecule = read_smiles(decoy)
        assert not next(rule.template.reactants.get_mapping(molecule), None), (
            f'{rule.id} matches its own decoy {decoy!r}')


@pytest.mark.parametrize('rule', RULES, ids=lambda rule: rule.name)
def test_a_row_reveals_the_element_it_says_it_does(rule):
    """The `protects` column, checked against the atoms the row keeps.

    `:1` is the atom being unmasked, so the FIRST category names its element -- which also pins the order
    of a two-part cell (`hydroxyl,amine`, not `amine,hydroxyl`) -- and every further category has to appear
    somewhere among the kept atoms.
    """
    molecule = read_smiles(rule.protected)
    mapping = next(rule.template.reactants.get_mapping(molecule))
    numbers = rule.template.reactant_map_numbers        # query index -> map number
    kept = {numbers[index]: molecule.atom(atom).atomic_symbol
            for index, atom in mapping.items() if index in numbers}

    assert kept[1] == PROTECTS_ELEMENTS[rule.protects[0]], (
        f'{rule.id} protects {rule.protects[0]!r} but its :1 atom is {kept[1]}')
    for what in rule.protects[1:]:
        assert PROTECTS_ELEMENTS[what] in kept.values(), (
            f'{rule.id} protects {what!r} and keeps no {PROTECTS_ELEMENTS[what]}')


# --- the table as a file --------------------------------------------------------------------------

def test_the_table_reads_and_ids_are_one_space():
    """103 rows numbered 1..103, which is what a positional id in a ported table has to be."""
    rows = read_table('protective.tsv')
    assert [int(row['id']) for row in rows] == list(range(1, len(rows) + 1))
    assert len(RULES) == len(rows) == 103


def test_a_protecting_group_is_named_once():
    """`deprotect()` selects by name, so a duplicate would silently shadow a pattern."""
    names = [rule.name for rule in RULES]
    assert len(set(names)) == len(names)


def test_a_rows_rule_id_is_table_qualified():
    for rule in RULES:
        table, _, number = rule.id.partition(':')
        assert table == 'protective', rule.id
        assert number.isdigit(), rule.id


def test_every_row_deletes_something():
    """The loader refuses a row that deletes nothing, and that is also the strip loop's termination
    argument: every pass removes an atom, so the molecule shrinks and the walk cannot cycle."""
    for rule in RULES:
        assert rule.template.deleted_atoms, rule.id


def test_a_name_says_what_it_protects():
    """The name prefix and the `protects` cell are the same claim, so they must agree -- which catches the
    copy-paste that gives a new thiol row a hydroxyl's `protects`."""
    for rule in RULES:
        assert rule.name.startswith(rule.protects[0] + '_'), rule.id


# --- specificity ----------------------------------------------------------------------------------

def test_the_rules_come_out_most_specific_first():
    sizes = [rule.size for rule in RULES]
    assert sizes == sorted(sizes, reverse=True)


def test_the_general_patterns_land_last():
    """The general patterns, named explicitly: each is a subset of a dozen others, so each is what a lost
    sort breaks first."""
    tail = {rule.name for rule in RULES[-6:]}
    for name in ('hydroxyl_methyl', 'hydroxyl_ethyl', 'hydroxyl_acyl', 'hydroxyl_allyl'):
        assert name in tail, f'{name} is a general pattern and should sort to the end'
    assert BY_NAME['hydroxyl_tbu'].size < BY_NAME['hydroxyl_boc'].size


def test_reversing_the_table_does_not_change_the_answer():
    """Sorting by pattern size makes specificity arithmetic, so destroying the file order changes
    nothing -- the general patterns need not be kept at the bottom by hand."""
    molecule = read_smiles('CC(C)OC(=O)OC(C)(C)C')      # Boc on isopropanol
    forward = [(claim.rule.name, claim.atoms) for claim in _claims(molecule, RULES)]
    shuffled = tuple(sorted(reversed(RULES), key=lambda rule: -rule.size))
    assert [(claim.rule.name, claim.atoms) for claim in _claims(molecule, shuffled)] == forward
    assert [name for name, _ in forward] == ['hydroxyl_boc']


def test_without_the_sort_a_boc_becomes_a_carbonate():
    """What the sort is FOR, measured rather than asserted.

    A Boc contains a tert-butyl ether's worth of atoms, so `hydroxyl_tbu` offered the site first strips
    the tert-butyl off the carbamate and leaves a carbonate -- not a deprotection of anything.
    """
    molecule = read_smiles('CC(C)OC(=O)OC(C)(C)C')
    wrong = {skeleton(p) for rxn in BY_NAME['hydroxyl_tbu'].template(molecule)
             for p in rxn.products}
    assert skeleton(read_smiles('C(O)(OC(C)C)=O')) in wrong
    # and the sorted walk never offers it that site
    assert protective_groups(molecule) == {'hydroxyl_boc': 1}


# --- the report -----------------------------------------------------------------------------------

def test_the_report_counts_sites_and_not_readings():
    """A tert-butyl's three methyls make six mappings of one site; the claim walk collapses them."""
    assert protective_groups(read_smiles('CC(C)OC(C)(C)C')) == {'hydroxyl_tbu': 1}


def test_two_sites_of_one_rule_are_two():
    """A bis-Boc diamine carries two Boc groups, which is a different fact from carrying one."""
    assert protective_groups(read_smiles('O=C(OC(C)(C)C)NCCCNC(=O)OC(C)(C)C')) == {'amine_boc': 2}


def test_two_groups_on_one_atom_are_two():
    """N,N-di-Boc: one nitrogen, two masks, and the second is not a shadow of the first.

    THIS IS WHY A CLAIM IS OVER THE DELETED ATOMS AND NOT THE WHOLE MATCH.  Both Boc groups match with the
    same nitrogen at `:1` and share no other atom, so a claim over the whole match drops the second as an
    overlap -- and then something smaller reaches the orphaned tert-butyl and reports an ether.
    """
    molecule = read_smiles('c1ccccc1CN(C(=O)OC(C)(C)C)C(=O)OC(C)(C)C')
    molecule.thiele()
    assert protective_groups(molecule) == {'amine_boc': 2}
    outcome = next(deprotect(molecule))
    assert outcome.names == ('amine_boc', 'amine_boc')
    assert {skeleton(p) for p in outcome.reaction.products} == {skeleton(read_smiles('c1ccccc1CN'))}


def test_two_rows_sharing_a_revealed_atom_both_claim_and_both_come_off():
    """`hydroxyl_amine_acetone`'s own example is an N-acetylated amino alcohol behind an acetonide.

    The two rows overlap on exactly the atom one of them reveals and both are claimed in ONE pass, no
    cascade needed -- again the deleted-atom claim.  So the full strip opens the ring AND removes the
    acetamide, which is why the row's `cleaved` column is checked against a single-rule strip instead.
    """
    molecule = read_smiles('N1(C(C)=O)C(C)(OC(C1)C)C')
    outcome = next(deprotect(molecule))
    assert outcome.names == ('hydroxyl_amine_acetone', 'amine_acyl')
    assert {skeleton(p) for p in outcome.reaction.products} == {skeleton(read_smiles('OC(C)CN'))}


def test_a_revealed_atom_must_keep_a_substituent():
    """Three shapes where two rows together would otherwise dissolve the substrate.

    Deprotection reveals a functional group ON something.  Each row guards its own site with a degree
    primitive, but two rows together can consume every neighbour of the atom they reveal, and then the
    "product" is a bare heteroatom -- water here, or ammonia.  The first is the sharpest: an
    ArCH2-O-C(=O)CH3 ester is a protected alcohol if the acetyl is the mask and a protected acid if the
    benzyl is, and exactly one can be true.
    """
    dmab = read_smiles('C(C)(C)CC(Nc1ccc(cc1)COC(C)=O)=C1C(=O)CC(CC1=O)(C)C')
    assert protective_groups(dmab) == {'hydroxyl_dmab_enamine': 1}
    assert stripped(dmab) == {skeleton(read_smiles('CC(=O)O'))}, 'stripped past acetic acid to water'

    ether = read_smiles('CCOC(C)(C)C')                                   # tert-butyl ethyl ether
    assert protective_groups(ether) == {'hydroxyl_tbu': 1}
    assert stripped(ether) == {skeleton(read_smiles('CCO'))}, 'stripped past ethanol to water'

    phth = read_smiles('O=C1N(Cc2ccccc2)C(=O)c2ccccc21')                 # N-benzylphthalimide
    phth.thiele()
    assert protective_groups(phth) == {'amine_phth': 1}
    assert stripped(phth) == {skeleton(read_smiles('c1ccccc1CN'))}, 'stripped past benzylamine to ammonia'


def test_an_unprotected_molecule_reports_nothing_and_deprotects_to_nothing():
    """Absent rather than zero, which is what makes `name in mol.protective_groups()` the presence test."""
    molecule = read_smiles('c1ccccc1CCN')
    assert protective_groups(molecule) == {}
    assert next(deprotect(molecule), None) is None


def test_a_protective_hit_carries_its_row_id():
    """The report with each row's id beside its count, for a consumer that persists ids."""
    molecule = read_smiles('O=C(OC(C)(C)C)NCCCNC(=O)OC(C)(C)C')
    hits = protective_group_hits(molecule)
    assert [(hit.name, hit.count) for hit in hits] == [('amine_boc', 2)]
    assert hits[0].id == BY_NAME['amine_boc'].id
    assert hits[0].id.startswith('protective:')
    assert protective_groups(molecule) == {hit.name: hit.count for hit in hits}


def test_the_report_and_the_strip_name_the_same_groups():
    """The two read one claim walk, so they cannot disagree by construction."""
    from collections import Counter
    molecule = read_smiles('O=C(OC(C)(C)C)NCCOC(C)=O')
    outcome = next(deprotect(molecule))
    assert Counter(outcome.names) == Counter(protective_groups(molecule))


# --- the strip ------------------------------------------------------------------------------------

def test_the_default_is_one_outcome_and_it_is_the_full_strip():
    molecule = read_smiles('O=C(OC(C)(C)C)NCCOC(C)=O')      # Boc-amine plus an acetate
    outcomes = list(deprotect(molecule))
    assert len(outcomes) == 1
    assert outcomes[0].names == ('amine_boc', 'hydroxyl_acyl')
    assert {skeleton(p) for p in outcomes[0].reaction.products} == {skeleton(read_smiles('NCCO'))}


def test_the_reaction_has_the_untouched_molecule_as_its_reactant():
    """A deprotection is a reaction and never an in-place mutation: the reactant is what went in, and the
    molecule handed to `deprotect()` still carries its protecting group afterwards."""
    molecule = read_smiles('CC(C)OC(=O)OC(C)(C)C')
    before = str(molecule)
    outcome = next(deprotect(molecule))
    assert len(outcome.reaction.reactants) == 1
    assert skeleton(outcome.reaction.reactants[0]) == skeleton(molecule)
    assert str(molecule) == before, 'deprotect() mutated its argument'


def test_a_counter_ion_survives_the_strip():
    """An untouched component of a touched input comes out as its own product.

    A template reports products SPLIT, so a multi-pass strip has to reunion them between passes; taking
    only the largest piece loses the salt on the second pass and not the first.
    """
    molecule = read_smiles('O=C(OC(C)(C)C)NCCCNC(=O)OC(C)(C)C.Cl')
    products = stripped(molecule)
    assert skeleton(read_smiles('Cl')) in products
    assert skeleton(read_smiles('NCCCN')) in products


def test_the_full_strip_takes_every_site_of_a_rule():
    """The FULL strip is exhaustive: a bis-Boc diamine loses both, in one outcome and not two.

    That is what "full" means and not a law about reagents -- `partial=True` offers the mono-Boc too.
    """
    outcome = next(deprotect(read_smiles('O=C(OC(C)(C)C)NCCCNC(=O)OC(C)(C)C')))
    assert outcome.names == ('amine_boc', 'amine_boc')
    assert outcome.rule_ids == (BY_NAME['amine_boc'].id,) * 2
    assert {skeleton(p) for p in outcome.reaction.products} == {skeleton(read_smiles('NCCCN'))}


def test_a_multi_atom_keep_comes_back_as_one_molecule():
    """The acetonide rows keep four atoms, so their product side restates three bonds; a forgotten one
    severs the molecule, so this is checked as connectivity and not only as equality."""
    products = stripped(read_smiles('CC1COC(C)(C)O1'))
    assert products == {skeleton(read_smiles('CC(O)CO'))}
    assert len(products) == 1


# --- partial --------------------------------------------------------------------------------------

def test_partial_yields_every_subset_largest_first():
    """2^k - 1 outcomes for k SITES, the full strip first, narrowing from there."""
    molecule = read_smiles('O=C(OC(C)(C)C)NCCOC(C)=O')      # two sites, two rules
    outcomes = list(deprotect(molecule, partial=True))
    assert len(outcomes) == 3 == 2 ** 2 - 1
    assert [len(o.names) for o in outcomes] == [2, 1, 1]
    assert outcomes[0].names == next(deprotect(molecule)).names
    assert {o.names for o in outcomes[1:]} == {('amine_boc',), ('hydroxyl_acyl',)}


def test_partial_over_one_rule_with_two_sites_offers_the_mono():
    """Two sites of ONE rule are two choices, chemistry not going to completion on request.

    The two sites here are symmetry-equivalent, so their subsets give the same products and are deduped by
    product multiset -- two outcomes, not three.
    """
    outcomes = list(deprotect(read_smiles('O=C(OC(C)(C)C)NCCCNC(=O)OC(C)(C)C'), partial=True))
    assert [o.names for o in outcomes] == [('amine_boc', 'amine_boc'), ('amine_boc',)]
    assert {skeleton(p) for p in outcomes[1].reaction.products} == \
           {skeleton(read_smiles('O=C(OC(C)(C)C)NCCCN'))}


def test_a_di_boc_amine_can_give_up_one_boc():
    """R-N(Boc)2 -> R-NHBoc, the case a rule-level enumeration cannot express: one nitrogen, two masks,
    the second harder than the first."""
    molecule = read_smiles('c1ccccc1CN(C(=O)OC(C)(C)C)C(=O)OC(C)(C)C')
    molecule.thiele()
    outcomes = list(deprotect(molecule, partial=True))
    assert [o.names for o in outcomes] == [('amine_boc', 'amine_boc'), ('amine_boc',)]
    assert {skeleton(p) for p in outcomes[1].reaction.products} == \
           {skeleton(read_smiles('O=C(OC(C)(C)C)NCc1ccccc1'))}


def test_partial_is_lazy():
    """The subsets are generated and not enumerated, which is what makes 2^k safe to offer with no cap."""
    molecule = read_smiles('O=C(OC(C)(C)C)NCCOC(C)=O')
    walker = deprotect(molecule, partial=True)
    assert next(walker).names == ('amine_boc', 'hydroxyl_acyl')     # and the rest is never built


def test_a_partial_subset_does_not_reopen_a_shadow():
    """THE CORRECTNESS ARGUMENT FOR PARTIAL DEPROTECTION.

    A molecule carrying a Boc AND a real tert-butyl ether gives `hydroxyl_tbu` one site it owns and one the
    Boc took from it.  Claims are computed over the WHOLE table on every pass and the subset only chooses
    which to act on; claiming with the subset instead re-offers the Boc's tert-butyl half as a carbonate.
    """
    molecule = read_smiles('CC(C)OC(=O)OC(C)(C)C.CCOC(C)(C)C')
    outcomes = {o.names: {skeleton(p) for p in o.reaction.products}
                for o in deprotect(molecule, partial=True)}
    assert set(outcomes) == {('hydroxyl_boc', 'hydroxyl_tbu'), ('hydroxyl_boc',), ('hydroxyl_tbu',)}
    only_tbu = outcomes[('hydroxyl_tbu',)]
    assert skeleton(read_smiles('CCO')) in only_tbu
    assert skeleton(read_smiles('CC(C)OC(=O)OC(C)(C)C')) in only_tbu, 'the Boc was not left whole'
    assert skeleton(read_smiles('C(O)(OC(C)C)=O')) not in only_tbu, 'the Boc became a carbonate'


def test_asking_for_a_group_that_is_only_a_shadow_yields_nothing():
    """There is no tert-butyl ether in a Boc: nothing, rather than a wrong answer or a refusal."""
    assert next(deprotect(read_smiles('CC(C)OC(=O)OC(C)(C)C'), ('hydroxyl_tbu',)), None) is None


# --- selection ------------------------------------------------------------------------------------

def test_names_selects_rows():
    molecule = read_smiles('O=C(OC(C)(C)C)NCCOC(C)=O')
    assert stripped(molecule, 'amine_boc') == {skeleton(read_smiles('NCCOC(C)=O'))}
    assert stripped(molecule, 'hydroxyl_acyl') == {skeleton(read_smiles('OCCNC(=O)OC(C)(C)C'))}


def test_protects_selects_by_what_comes_off():
    """An orthogonal question to `names`, because a reagent class cuts across the six categories."""
    molecule = read_smiles('O=C(OC(C)(C)C)NCCOC(C)=O')
    assert next(deprotect(molecule, protects='amine')).names == ('amine_boc',)
    assert next(deprotect(molecule, protects=('hydroxyl',))).names == ('hydroxyl_acyl',)
    assert next(deprotect(molecule, protects=('amine', 'hydroxyl'))).names == \
        ('amine_boc', 'hydroxyl_acyl')


def test_an_unknown_name_says_so_and_lists_the_alternatives():
    """Rather than yielding nothing, which is what a typo and an absent group look like together."""
    with pytest.raises(ValueError) as exc:
        list(deprotect(read_smiles('CCO'), ('amine_bok',)))
    assert 'amine_bok' in str(exc.value) and 'amine_boc' in str(exc.value)


def test_an_unknown_protects_says_so():
    with pytest.raises(ValueError) as exc:
        list(deprotect(read_smiles('CCO'), protects='alcohol'))
    assert 'alcohol' in str(exc.value) and 'hydroxyl' in str(exc.value)


def test_the_protects_column_takes_only_the_six():
    """The vocabulary is closed, so a new row cannot invent a seventh category by typing it."""
    for rule in RULES:
        for what in rule.protects:
            assert what in PROTECTS_ELEMENTS, rule.id


# --- the key structure the accessor hands out -----------------------------------------------------

def test_protective_rules_is_keyed_by_the_name_deprotect_selects_by():
    """The name is already the primary key -- a duplicate is refused at load for exactly this reason.

    The tuple carried only the sort, and dict insertion order carries it just as well, so `.values()`
    is the specificity order and nothing that iterates loses it.
    """
    rules = protective_rules()
    assert isinstance(rules, dict)
    assert all(name == rule.name for name, rule in rules.items())
    assert len(rules) == len(read_table('protective.tsv'))


def test_the_specificity_order_survives_the_key():
    """MOST SPECIFIC FIRST is load-bearing -- `hydroxyl_tbu` offered a Boc yields a carbonate -- and it
    is the values' order, not something a caller has to re-sort."""
    sizes = [rule.size for rule in protective_rules().values()]
    assert sizes == sorted(sizes, reverse=True)
    assert protective_rules()['hydroxyl_boc'].size > protective_rules()['hydroxyl_tbu'].size


# --- caching and laziness -------------------------------------------------------------------------

def test_the_templates_are_cached():
    assert protective_rules() is protective_rules()


def test_nothing_loads_at_import():
    """A process that only writes SMILES must not compile 103 templates."""
    from subprocess import run
    from sys import executable
    script = ('import chython.reactions._tables as t;'
              'print(t._PROTECTIVE_CACHE, len(t._FUNCTIONAL_CACHE), t._RULES_CACHE)')
    result = run([executable, '-c', script], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout.split() == ['None', '0', 'None'], result.stdout


# --- the container methods ------------------------------------------------------------------------

def test_the_methods_are_on_the_container():
    """Registered by injection, because `MoleculeContainer` is a `cdef class`."""
    molecule = read_smiles('CC(C)OC(=O)OC(C)(C)C')
    assert molecule.protective_groups() == {'hydroxyl_boc': 1}
    assert next(molecule.deprotect()).names == ('hydroxyl_boc',)
    assert next(molecule.deprotect('hydroxyl_boc')).names == ('hydroxyl_boc',)
    assert next(molecule.deprotect(protects='hydroxyl')).names == ('hydroxyl_boc',)
    assert len(list(molecule.deprotect(partial=True))) == 1


def test_deprotect_is_always_an_iterator():
    """So `partial` widens the answer set without changing its shape, and one answer is one `next`."""
    molecule = read_smiles('c1ccccc1CCN')
    assert next(molecule.deprotect(), None) is None
    assert list(molecule.deprotect(partial=True)) == []


# --- what the whole corpus does together ----------------------------------------------------------

def test_every_row_is_reachable_through_the_claim_walk():
    """Each row claims its own example, so no row is permanently shadowed by another.

    The direct acceptance test above still passes for a row that can never claim anything, so the two
    together are what prove a row usable.
    """
    unreachable = [rule.id for rule in RULES
                   if rule.name not in protective_groups(read_smiles(rule.protected))]
    assert not unreachable, f'{len(unreachable)} rows never claim their own example: {unreachable}'


def test_every_row_strips_its_example_end_to_end():
    """The same rows through `deprotect()` rather than their own template: the claim walk, the specificity
    sort and the multi-pass strip together give each row's documented answer.

    Each is asked for ITS OWN rule, since `cleaved` documents what this row reveals and an example may
    carry a second group an unrestricted strip would rightly take off too.  Claiming still runs over the
    whole table, so the shadowing this test exists to catch is caught.
    """
    wrong = []
    for rule in RULES:
        molecule = read_smiles(rule.protected)
        products = stripped(molecule, rule.name)
        if products is None or skeleton(read_smiles(rule.cleaved)) not in products:
            wrong.append(rule.id)
    assert not wrong, f'{len(wrong)} rows do not strip their own example: {wrong}'


def test_the_subset_count_is_two_to_the_k_sites():
    """The report counts sites, so `sum(report.values())` is the exponent -- three sites on two rules here,
    none symmetry-equivalent, so nothing is deduped away and the count is exact."""
    molecule = read_smiles('O=C(OC(C)(C)C)NCC(OC(C)=O)COC(C)=O')   # one Boc, two distinct acetates
    k = sum(protective_groups(molecule).values())
    assert k == 3
    assert sum(1 for _ in deprotect(molecule, partial=True)) == 2 ** k - 1
    assert sum(len(list(combinations(range(k), size))) for size in range(1, k + 1)) == 2 ** k - 1


def test_a_deprotection_comes_back_with_a_1_1_mapping():
    # The reactor imposes the mapping on every path it has, this one included: contiguous from 1 over
    # the atoms on both sides, and the group it cleaved is on one side only and stays 0.
    mol = read_smiles('c1ccccc1NC(=O)OC(C)(C)C')
    mol.canonicalize()
    rxn = next(mol.deprotect()).reaction
    reactant, product = rxn.reactants[0], rxn.products[0]

    numbers = sorted(product.map_number_of(n) for n in product.atom_numbers)
    assert numbers == list(range(1, len(numbers) + 1))
    forward = {reactant.map_number_of(n): n for n in reactant.atom_numbers
               if reactant.map_number_of(n)}
    assert len(forward) == len(numbers)                  # 1-1: no number used twice on the left
    for n in product.atom_numbers:                      # and paired atoms are the same element
        assert product.element_of(n) == reactant.element_of(forward[product.map_number_of(n)])
    view = rxn.modeling_view()
    assert view.collisions == {'reactants': (), 'products': ()}
    assert view.unmapped == {'reactants': len(reactant.atom_numbers) - len(numbers), 'products': 0}


def test_a_counter_ion_keeps_its_number_across_a_deprotection():
    # `_patch_within` re-unions the products so the working molecule stays one container, and a
    # REMAPPING union gives the counter-ion an atom number no reactant atom has -- the salt would come
    # back unmapped even though it is on both sides.
    mol = read_smiles('CC(C)(C)OC(=O)NCc1ccccc1.Cl')
    mol.canonicalize()
    rxn = next(mol.deprotect()).reaction
    reactant = rxn.reactants[0]
    product, chloride = next((p, n) for p in rxn.products for n in p.atom_numbers
                             if p.element_of(n) == 17)
    number = product.map_number_of(chloride)
    assert number
    source = next(n for n in reactant.atom_numbers if reactant.map_number_of(n) == number)
    assert reactant.element_of(source) == 17
