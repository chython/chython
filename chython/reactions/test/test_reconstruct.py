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
from collections import Counter
import pytest
from pytest import raises

from .._enumerate import _deprotect_toward
from ...core import ReactionContainer, read_smiles as smiles


def test_the_accessor_names_the_package_when_unregistered():
    from ...core._core import _reaction_reconstruct_fn, _set_reconstruct_fn
    try:
        kept = _reaction_reconstruct_fn()       # another test module may have registered it already
    except ImportError:
        kept = None
    _set_reconstruct_fn(None)
    try:
        with raises(ImportError, match='chython.reactions'):
            _reaction_reconstruct_fn()
    finally:
        _set_reconstruct_fn(kept)


def test_purification_is_the_identity_mapping():
    rxn = ReactionContainer([smiles('CCO'), smiles('O')], [smiles('CCO')])
    assert rxn.reconstruct_mapping() == ('purification',)
    product = rxn.products[0]
    numbers = {product.map_number_of(n) for n in product.atom_numbers}
    assert 0 not in numbers                                     # every product atom is numbered
    ethanol = rxn.reactants[0]
    assert numbers == {ethanol.map_number_of(n) for n in ethanol.atom_numbers}

    # PER ATOM, and not only as a set: the set assertion passes for any consistent permutation inside
    # ethanol.  `(element, degree, implicit_h)` separates all three of its atoms.
    source = {ethanol.map_number_of(n): n for n in ethanol.atom_numbers}
    for n in product.atom_numbers:
        m = source[product.map_number_of(n)]
        assert (product.atom(n).element, product.degree_of(n), product.atom(n).implicit_h) == \
               (ethanol.atom(m).element, ethanol.degree_of(m), ethanol.atom(m).implicit_h)


def test_a_multiproduct_record_is_refused():
    from ...core import Log                     # severities are `INFO`/`REFUSED` constants, not an enum
    rxn = ReactionContainer([smiles('CCOC(=O)CC'), smiles('O')], [smiles('CCO'), smiles('CCC(=O)O')])
    assert rxn.reconstruct_mapping() == ()
    log = rxn.log
    refused = log.refused()
    assert len(refused) == 1
    assert refused[0].rule == 'reconstruct:multiproduct'
    assert refused[0].stage == 'reconstruct'        # the pass names the stage; nothing above it does


def test_an_empty_side_is_refused():
    from ...core import Log
    rxn = ReactionContainer([smiles('CCO')], [])
    assert rxn.reconstruct_mapping() == ()
    log = rxn.log
    assert [r.rule for r in log.refused()] == ['reconstruct:empty']


def test_an_empty_reactant_side_is_refused():
    from ...core import Log
    rxn = ReactionContainer([], [smiles('CCO')])
    assert rxn.reconstruct_mapping() == ()
    log = rxn.log
    assert [r.rule for r in log.refused()] == ['reconstruct:empty']


def test_incoming_map_numbers_are_replaced_unconditionally():
    # 10/20/30 rather than 1/2/3 keeps the test non-vacuous: `_number_inputs` assigns 1..N, so the
    # incoming numbers cannot be mistaken for the derived ones even if `_clear` is removed.
    ethanol_in = smiles('CCO')
    for atom_n, mn in zip(list(ethanol_in.atom_numbers), [10, 20, 30]):
        ethanol_in.set_map_number(atom_n, mn)
    ethanol_out = smiles('CCO')
    for atom_n, mn in zip(list(ethanol_out.atom_numbers), [10, 20, 30]):
        ethanol_out.set_map_number(atom_n, mn)
    rxn = ReactionContainer([ethanol_in], [ethanol_out])
    assert rxn.reconstruct_mapping() == ('purification',)
    product = rxn.products[0]
    product_numbers = {product.map_number_of(n) for n in product.atom_numbers}
    assert 10 not in product_numbers
    assert 20 not in product_numbers
    assert 30 not in product_numbers
    assert 0 not in product_numbers   # every product atom has been numbered
    reactant = rxn.reactants[0]
    reactant_numbers = {reactant.map_number_of(n) for n in reactant.atom_numbers}
    assert product_numbers == reactant_numbers


def test_unexplained_record_emits_lost_and_clears_product():
    # Ethanol -> propane is not a reaction: an underivable mapping comes back with a LOST line rather
    # than silently disappearing, and the product's numbers are cleared and not replaced.
    from ...core import Log, LOST
    rxn = ReactionContainer([smiles('[CH3:1][CH2:2][OH:3]')], [smiles('[CH3:4][CH2:5][CH3:6]')])
    result = rxn.reconstruct_mapping()
    log = rxn.log
    assert result == ()
    lost = [r for r in log if r.severity == LOST]
    assert len(lost) == 1
    assert lost[0].rule == 'reconstruct:unexplained'
    product = rxn.products[0]
    product_numbers = {product.map_number_of(n) for n in product.atom_numbers}
    assert product_numbers == {0}  # all cleared, none replaced


def test_the_recorded_reaction_comes_back_canonical():
    # a Kekule aromatic input is rewritten in place: the answer is a normalized, mapped record
    rxn = ReactionContainer([smiles('C1=CC=CC=C1O')], [smiles('C1=CC=CC=C1O')])
    assert rxn.reconstruct_mapping() == ('purification',)
    assert all(rxn.reactants[0].atom(n).hybridization == 4
               for n in rxn.reactants[0].atom_numbers if rxn.reactants[0].atom(n).element == 6)


def test_translate_carries_the_inputs_numbering_onto_a_reactor_product():
    # The reactor numbers its own output; every rung of the ladder reads the INPUTS' numbering.
    # Asserted as set membership, not as a literal, so the test survives whichever scheme the reactor
    # imposes: the product's carbonyl carbon has to come back with a number the ACID input carries.
    from .._enumerate import _run
    from .._reconstruct import _number_inputs, _translate
    from .._tables import reaction_rules

    acid, amine = smiles('CCC(=O)O'), smiles('CCN')
    acid.canonicalize()
    amine.canonicalize()
    _number_inputs([acid, amine])
    outcome = next(o for o in _run([acid, amine], reaction_rules()) if o.name == 'amidation')
    product = _translate(outcome.reaction, [acid, amine])[0]

    acid_numbers = {acid.map_number_of(n) for n in acid.atom_numbers}
    amine_numbers = {amine.map_number_of(n) for n in amine.atom_numbers}
    carbonyl = next(n for n in product.atom_numbers
                    if product.element_of(n) == 6
                    and any(product.order_of(n, m) == 2 and product.element_of(m) == 8
                            for m in product.neighbors_of(n)))
    nitrogen = next(n for n in product.atom_numbers if product.element_of(n) == 7)
    assert product.map_number_of(carbonyl) in acid_numbers
    assert product.map_number_of(nitrogen) in amine_numbers


def test_a_direct_reaction_is_found_and_mapped():
    rxn = ReactionContainer([smiles('CCC(=O)O'), smiles('CCN')], [smiles('CCC(=O)NCC')])
    assert rxn.reconstruct_mapping() == ('react:amidation',)
    product = rxn.products[0]
    assert all(product.map_number_of(n) for n in product.atom_numbers
               if product.atom(n).element != 1)
    acid, amine = rxn.reactants
    acid_numbers = {acid.map_number_of(n) for n in acid.atom_numbers}
    amine_numbers = {amine.map_number_of(n) for n in amine.atom_numbers}
    assert acid_numbers.isdisjoint(amine_numbers)
    carried = {product.map_number_of(n) for n in product.atom_numbers}
    assert carried & acid_numbers and carried & amine_numbers
    # the carbonyl carbon carries the acid's number and not the amine's -- the literal claim the
    # assertion above only implies
    carbonyl_c = next(n for n in product.atom_numbers
                      if product.atom(n).element == 6
                      and any(product.bond(n, nb).order == 2 and product.atom(nb).element == 8
                              for nb in product.neighbors_of(n)))
    assert product.map_number_of(carbonyl_c) in acid_numbers
    assert product.map_number_of(carbonyl_c) not in amine_numbers


def test_a_spectator_input_does_not_block_the_match():
    # toluene is a solvent no row touches.  This is the subset enumeration: without it `_run` refuses
    # the whole list, every input having to be touched, and nothing is found.
    rxn = ReactionContainer([smiles('CCC(=O)O'), smiles('CCN'), smiles('Cc1ccccc1')],
                            [smiles('CCC(=O)NCC')])
    assert rxn.reconstruct_mapping() == ('react:amidation',)


def test_an_untouched_salt_component_survives_and_is_numbered():
    # the recorded product carries a chloride the template never sees; it must come back numbered
    # from the input it arrived in, not dropped and not left at zero
    rxn = ReactionContainer([smiles('CCC(=O)O'), smiles('CCN.Cl')], [smiles('CCC(=O)NCC.Cl')])
    assert rxn.reconstruct_mapping() == ('react:amidation',)
    amine_input = rxn.reactants[1]
    product = rxn.products[0]
    chloride = [n for n in product.atom_numbers if product.atom(n).element == 17]
    assert len(chloride) == 1
    input_cl = next(n for n in amine_input.atom_numbers if amine_input.atom(n).element == 17)
    assert product.map_number_of(chloride[0]) == amine_input.map_number_of(input_cl)


def test_the_winning_rule_id_is_logged():
    from ...core import Log
    rxn = ReactionContainer([smiles('CCC(=O)O'), smiles('CCN')], [smiles('CCC(=O)NCC')])
    assert rxn.reconstruct_mapping()
    log = rxn.log
    ids = [r.message for r in log if r.rule == 'reconstruct:rules']
    # pin the value and not the count, which is trivially 1 when only one row matches
    assert ids == ['explained by reactions:1']


def test_a_partial_product_emits_the_partial_log_line():
    from ...core import Log, INFO
    # the water in the recorded product came from nowhere, so `_number_product` cannot account for that
    # component and the orchestrator emits `reconstruct:partial`
    rxn = ReactionContainer([smiles('CCC(=O)O'), smiles('CCN')], [smiles('CCC(=O)NCC.O')])
    result = rxn.reconstruct_mapping()
    log = rxn.log
    assert result == ('react:amidation',)                     # rung fires despite the unaccounted component
    partial = [r for r in log if r.rule == 'reconstruct:partial']
    assert len(partial) == 1 and partial[0].severity == INFO
    product = rxn.products[0]
    components = list(product.split())
    water = next(c for c in components if len(list(c.atom_numbers)) == 1)
    amide = next(c for c in components if len(list(c.atom_numbers)) > 1)
    assert all(product.map_number_of(n) == 0 for n in water.atom_numbers)
    assert all(product.map_number_of(n) != 0 for n in amide.atom_numbers
               if product.atom(n).element != 1)


def test_a_standalone_deprotection_is_explained():
    # Boc removal: the recorded product is the recorded input, minus the protecting group
    rxn = ReactionContainer([smiles('CC(C)(C)OC(=O)NCc1ccccc1'), smiles('Cl')],
                            [smiles('NCc1ccccc1')])
    labels = rxn.reconstruct_mapping()
    assert labels and all(label.startswith('deprotect:') for label in labels)
    product = rxn.products[0]
    assert all(product.map_number_of(n) for n in product.atom_numbers)
    revealed = rxn.reactants[0]
    numbers = {revealed.map_number_of(n) for n in revealed.atom_numbers}
    assert {product.map_number_of(n) for n in product.atom_numbers} <= numbers

    # PER ATOM, and not only as a set: the set assertion passes for a swap of two ring CH atoms or of
    # the benzylic CH2 for a ring CH.  The revealed N is the one exception -- Boc removal changes its
    # degree and implicit_h -- so it is checked on element alone.
    source = {revealed.map_number_of(n): n for n in revealed.atom_numbers}
    for n in product.atom_numbers:
        m = source[product.map_number_of(n)]
        if product.atom(n).element == 7:      # the revealed atom; update if the substrate changes
            assert product.atom(n).element == revealed.atom(m).element
        else:
            assert (product.atom(n).element, product.degree_of(n), product.atom(n).implicit_h) == \
                   (revealed.atom(m).element, revealed.degree_of(m), revealed.atom(m).implicit_h)


def test_a_partial_deprotection_is_explained():
    # one of two Boc groups comes off: reachable only because the rung enumerates partial strips
    rxn = ReactionContainer([smiles('CC(C)(C)OC(=O)N(C(=O)OC(C)(C)C)Cc1ccccc1')],
                            [smiles('CC(C)(C)OC(=O)NCc1ccccc1')])
    labels = rxn.reconstruct_mapping()
    assert labels and all(label.startswith('deprotect:') for label in labels)


def test_partial_deprotection_strips_only_the_subsets_that_fit_the_product():
    # one Boc among twelve acetates: 2^13 site subsets, and one of them has the recorded formula
    protected = smiles('CC(C)(C)OC(=O)NC' + 'C(OC(C)=O)' * 12 + 'C')
    product = smiles('NC' + 'C(OC(C)=O)' * 12 + 'C')
    assert [o.names for o in _deprotect_toward(protected, [product])][1:] == [('boc',)]
    rxn = ReactionContainer([protected], [product])
    assert rxn.reconstruct_mapping() == ('deprotect:amine_boc',)


def test_deprotect_then_react_composes():
    # the amine arrives Boc-protected; strip it, then amidate.  Neither rung alone explains this.
    rxn = ReactionContainer([smiles('CCC(=O)O'), smiles('CC(C)(C)OC(=O)NCC')],
                            [smiles('CCC(=O)NCC')])
    assert rxn.reconstruct_mapping() == ('deprotect:amine_boc+react:amidation',)


def test_a_reactor_outcome_is_compared_in_its_canonical_form():
    # the patch leaves the benzotriazole NH where the input had it; the recorded side is canonical
    rxn = ReactionContainer([smiles('Clc1ccc2[nH]nnc2c1'), smiles('OB(O)c1ccccc1')],
                            [smiles('c1ccc(-c2ccc3[nH]nnc3c2)cc1')])
    assert rxn.reconstruct_mapping() == ('react:suzuki',)


def test_an_unmasked_azole_nitrogen_is_settled_before_comparison():
    # `deprotect()` leaves the revealed pyrrole-type N at H_UNKNOWN; `kekule()` is what closes it
    rxn = ReactionContainer([smiles('CC(C)(C)OC(=O)n1ccc2ccccc21')], [smiles('c1ccc2[nH]ccc2c1')])
    assert rxn.reconstruct_mapping() == ('deprotect:azole_boc',)


def test_a_stripped_azole_is_settled_before_the_corpus_fires():
    rxn = ReactionContainer([smiles('CC(C)(C)OC(=O)n1ccc2cc(Br)ccc21'), smiles('OB(O)c1ccccc1')],
                            [smiles('c1ccc(-c2ccc3[nH]ccc3c2)cc1')])
    assert rxn.reconstruct_mapping() == ('deprotect:azole_boc+react:suzuki',)


def test_a_direct_reaction_outranks_the_composed_rung():
    # NOT an ordering pin: CCN carries no protecting group, so `_deprotect_then_react` returns at once
    # whatever its position in `_PHASES`.  `test_phases_ladder_order_is_pinned` is the ordering pin.
    rxn = ReactionContainer([smiles('CCC(=O)O'), smiles('CCN')], [smiles('CCC(=O)NCC')])
    assert rxn.reconstruct_mapping() == ('react:amidation',)


def test_phases_ladder_order_is_pinned():
    # Protection right after purification: a label names conditions, and installing a named protecting
    # group is one family of them whatever supplies the group.  Behavioural counterpart:
    # `test_a_named_group_installed_is_a_protection_and_any_other_acylation_a_reaction`.
    from .._reconstruct import _PHASES
    assert [f.__name__ for f in _PHASES] == ['_purification', '_protect', '_react', '_deprotect',
                                             '_deprotect_then_react']


def test_two_surviving_candidates_yield_one_label():
    # Two inputs both deprotect to benzylamine, so `_deprotect` yields two candidates that reproduce the
    # recorded product.  Exactly one is applied: applying both writes two conflicting numberings.
    from ...core import Log
    rxn = ReactionContainer(
        [smiles('CC(C)(C)OC(=O)NCc1ccccc1'), smiles('O=C(OCc1ccccc1)NCc1ccccc1')],
        [smiles('NCc1ccccc1')],
    )
    labels = rxn.reconstruct_mapping()
    log = rxn.log
    assert len(labels) == 1 and labels[0].startswith('deprotect:')
    rules_logged = [r for r in log if r.rule == 'reconstruct:rules']
    assert len(rules_logged) == 1
    # pin the content and not only the count: the message must name a row and not a label
    assert rules_logged[0].message.startswith('explained by protective:')


def test_a_protection_is_explained():
    # Boc protection of benzylamine with Boc anhydride.  The rung deprotects the recorded product,
    # pairs the revealed fragment against the input amine, and numbers back.
    rxn = ReactionContainer([smiles('NCc1ccccc1'), smiles('CC(C)(C)OC(=O)OC(=O)OC(C)(C)C')],
                            [smiles('CC(C)(C)OC(=O)NCc1ccccc1')])
    labels = rxn.reconstruct_mapping()
    assert labels and all(label.startswith('protect:') for label in labels)
    product = rxn.products[0]
    amine = rxn.reactants[0]
    numbers = {amine.map_number_of(n) for n in amine.atom_numbers}
    carried = {product.map_number_of(n) for n in product.atom_numbers if product.map_number_of(n)}
    assert carried == numbers                       # the whole revealed fragment is numbered
    # the protecting group's own atoms are NEW and correctly carry no number
    assert any(product.map_number_of(n) == 0 for n in product.atom_numbers)

    # PER ATOM, and not only as a set.  `(element, degree, implicit_h)` separates the distinguishable
    # atoms; the ring CH positions share a triple and cannot be separated, being one symmetry orbit.
    source = {amine.map_number_of(n): n for n in amine.atom_numbers}
    for n in product.atom_numbers:
        mn = product.map_number_of(n)
        if mn == 0:        # protecting group atom, genuinely new -- no source to check against
            continue
        m = source[mn]
        if product.atom(n).element == 7:      # the reaction site; update if the substrate changes
            assert product.atom(n).element == amine.atom(m).element
        else:
            assert (product.atom(n).element, product.degree_of(n), product.atom(n).implicit_h) == \
                   (amine.atom(m).element, amine.degree_of(m), amine.atom(m).implicit_h)


@pytest.mark.parametrize('source', ['CC(=O)OC(C)=O', 'CC(=O)Cl', 'CC(=O)O'])
def test_a_named_group_installed_is_a_protection_whatever_supplies_it(source):
    """An acetate from Ac2O, AcCl or AcOH is one label, and the group's source is an agent: unmapped."""
    rxn = smiles(f'OCc1ccccc1.{source}>>CC(=O)OCc1ccccc1')
    assert rxn.reconstruct_mapping() == ('protect:hydroxyl_acyl',)
    acyl = next(m for m in rxn.reactants if 'O' in m and len(m) < 8)
    assert all(not acyl.map_number_of(n) for n in acyl.atom_numbers)


def test_a_named_group_installed_is_a_protection_and_any_other_acylation_a_reaction():
    benzoyl = ReactionContainer([smiles('OC(=O)c1ccccc1'), smiles('CCN')], [smiles('CCNC(=O)c1ccccc1')])
    assert benzoyl.reconstruct_mapping() == ('protect:amine_benzoate',)
    other = ReactionContainer([smiles('OC(=O)C1CCCCC1'), smiles('CCN')], [smiles('CCNC(=O)C1CCCCC1')])
    assert other.reconstruct_mapping() == ('react:amidation',)


def test_canonicalize_does_not_reconstruct():
    # NEVER IMPLICITLY COMPOSED: a caller asking for a canonical representation has not asked for a
    # mapping to be invented.  This would fail if `canonicalize()` called `reconstruct_mapping`, the
    # amidation row firing and writing nonzero numbers onto the product.
    rxn = ReactionContainer([smiles('CCC(=O)O'), smiles('CCN')], [smiles('CCC(=O)NCC')])
    rxn.canonicalize()
    product = rxn.products[0]
    assert not any(product.map_number_of(n) for n in product.atom_numbers)


def test_a_grossly_larger_product_is_refused():
    from ...core import Log
    # a 60-atom product from one 2-atom input: the corpus has nothing honest to say about this
    rxn = ReactionContainer([smiles('CO')], [smiles('C' * 60)])
    assert rxn.reconstruct_mapping() == ()
    log = rxn.log
    assert [r.rule for r in log.refused()] == ['reconstruct:unbalanced']
    # a bound on the SEARCH and not a rejection: the record still gets its ordinary unexplained line
    from ...core import LOST
    assert any(r.rule == 'reconstruct:unexplained' and r.severity == LOST for r in log)


def test_the_filter_is_off_below_the_floor():
    from ...core import Log
    # under `min_filter_size` the ratio is not consulted at all, however lopsided it looks
    rxn = ReactionContainer([smiles('CO')], [smiles('CCCCCCCCCC')])
    rxn.reconstruct_mapping()
    log = rxn.log
    assert not [r for r in log.refused() if r.rule == 'reconstruct:unbalanced']


def test_the_filter_can_be_disabled():
    from ...core import Log
    rxn = ReactionContainer([smiles('CO')], [smiles('C' * 60)])
    rxn.reconstruct_mapping(max_size_ratio=0.)
    log = rxn.log
    assert not [r for r in log.refused() if r.rule == 'reconstruct:unbalanced']


def test_a_protection_survives_the_filter():
    # Trityl protection of decanol.  THE THRESHOLDS ARE PASSED EXPLICITLY: at the defaults the filter
    # never engages and the test would pass with `_FILTER_EXEMPT` deleted.  Forced on, only the
    # exemption can let a `protect:` label out.
    rxn = ReactionContainer([smiles('OCCCCCCCCCC')],
                            [smiles('C(c1ccccc1)(c1ccccc1)(c1ccccc1)OCCCCCCCCCC')])
    assert len(rxn.products[0].atom_numbers) == 30       # 30 >= 10, and 30 >= 1.5 * 11
    labels = rxn.reconstruct_mapping(max_size_ratio=1.5, min_filter_size=10)
    assert labels and all(label.startswith('protect:') for label in labels)


def test_a_purification_survives_the_filter():
    # A purification's product IS one of its inputs, but the arithmetic does not know that: at a low
    # enough ratio it clears the bound.  Dropping `_purification` from `_FILTER_EXEMPT` turns this
    # answer into `()` with a `reconstruct:unexplained` line.
    from ...core import Log
    rxn = ReactionContainer([smiles('CCO'), smiles('O')], [smiles('CCO')])
    assert rxn.reconstruct_mapping(max_size_ratio=.5, min_filter_size=2) == ('purification',)
    log = rxn.log
    # confirm the filter genuinely engaged, or the test passes with `_purification` unreachable
    assert [r.rule for r in log.refused()] == ['reconstruct:unbalanced']


def test_a_reference_mapping_is_reproduced_exactly():
    from .._numbering import mapping_agrees
    for reference in _reference_records():
        reference.canonicalize()
        probe = reference.copy()
        assert probe.reconstruct_mapping(), 'nothing explained a record the corpus should explain'
        agreed, disagreed, missing = mapping_agrees(probe, reference)
        assert agreed, 'no product atom was traced back to an input'
        assert (disagreed, missing) == (0, 0)


def _reference_records():
    """Three public reactions, written with the mapping a chemist would draw.

    Hand-written rather than lifted from `mapping/golden.rdf`, a test not depending on a data file outside
    the installed package.  One per rung that can reproduce a whole product.

    THE COUPLING IS DELIBERATELY ASYMMETRIC -- 4-bromotoluene, not bromobenzene.  Biphenyl's two rings are
    one automorphism orbit of the PRODUCT while arriving from two different INPUTS, and `mapping_agrees`
    excuses a swap only within one input's orbits; for the symmetric spelling no mapping is the answer, so
    a reference naming one would state a convention rather than a fact.
    """
    return [
        ReactionContainer([smiles('[CH3:1][CH2:2][C:3](=[O:4])[OH:5]'), smiles('[CH3:6][CH2:7][NH2:8]')],
                          [smiles('[CH3:1][CH2:2][C:3](=[O:4])[NH:8][CH2:7][CH3:6]')]),
        ReactionContainer([smiles('[CH3:1][C:2]([CH3:3])([CH3:4])[O:5][C:6](=[O:7])'
                                  '[NH:8][CH2:9][c:10]1[cH:11][cH:12][cH:13][cH:14][cH:15]1')],
                          [smiles('[NH2:8][CH2:9][c:10]1[cH:11][cH:12][cH:13][cH:14][cH:15]1')]),
        ReactionContainer([smiles('[Br:1][c:2]1[cH:3][cH:4][c:5]([CH3:6])[cH:7][cH:8]1'),
                           smiles('[OH:9][B:10]([OH:11])[c:12]1[cH:13][cH:14][cH:15][cH:16][cH:17]1')],
                          [smiles('[CH3:6][c:5]1[cH:4][cH:3][c:2]([c:12]2[cH:13][cH:14][cH:15][cH:16]'
                                  '[cH:17]2)[cH:8][cH:7]1')]),
    ]


def test_the_reconstructed_mapping_is_1_1_and_starts_at_one():
    rxn = ReactionContainer([smiles('CCC(=O)O'), smiles('CCN')], [smiles('CCC(=O)NCC')])
    assert rxn.reconstruct_mapping() == ('react:amidation',)
    product = rxn.products[0]
    numbers = sorted(product.map_number_of(n) for n in product.atom_numbers)
    assert numbers == list(range(1, len(numbers) + 1))
    left = [m.map_number_of(n) for m in rxn.reactants for n in m.atom_numbers]
    assert sorted(x for x in left if x) == numbers        # 1-1, and the same set on both sides
    assert left.count(0) == 1                             # the acid's leaving OH, and only it
    assert rxn.modeling_view().collisions == {'reactants': (), 'products': ()}


def test_a_spectator_input_comes_back_unmapped():
    # A number on an atom the product never received is not part of a 1-1 mapping, and leaving it there
    # invites a reader to treat the toluene as a reagent that contributed atoms.
    rxn = ReactionContainer([smiles('CCC(=O)O'), smiles('CCN'), smiles('Cc1ccccc1')],
                            [smiles('CCC(=O)NCC')])
    assert rxn.reconstruct_mapping() == ('react:amidation',)
    toluene = rxn.reactants[2]
    assert all(toluene.map_number_of(n) == 0 for n in toluene.atom_numbers)


# --- what the walk is allowed to cost -------------------------------------------------------------

#: Solvents, a drying agent and a base that carry NO group the corpus names a slot for, so no row can
#: reach them: dichloromethane, THF, water, DMF, dioxane, hexane.  A record files its solvent among its
#: reagents, and the corpus median is 4 reagents, so this is the ordinary case and not a stress test.
_UNREACHABLE = ('ClCCl', 'C1CCOC1', 'O', 'CN(C)C=O', 'C1COCCO1', 'CCCCCC')


def _amidation(extra=()):
    return ReactionContainer([smiles('CCC(=O)O'), smiles('CCN'), *(smiles(s) for s in extra)],
                             [smiles('CCC(=O)NCC')])


def _pools(reaction):
    """The label, and how many times the walk applied a row to a pool of inputs to earn it."""
    from .. import _reconstruct

    offered = []
    kept = _reconstruct._outcomes

    def counting(rule, molecules):
        offered.append(len(molecules))
        return kept(rule, molecules)

    _reconstruct._outcomes = counting
    try:
        return reaction.reconstruct_mapping(), len(offered)
    finally:
        _reconstruct._outcomes = kept


def _mapping(reaction):
    """The written mapping, per molecule, as the thing the padding may not perturb."""
    return [sorted((n, m.map_number_of(n)) for n in m.atom_numbers) for m in reaction.molecules()]


def test_inert_inputs_do_not_change_the_answer():
    plain, padded = _amidation(), _amidation(_UNREACHABLE)
    assert plain.reconstruct_mapping() == padded.reconstruct_mapping() == ('react:amidation',)
    # PER ATOM, on the two molecules both records share: the acid, the amine and the amide
    written = _mapping(padded)
    assert _mapping(plain) == [*written[:2], written[-1]]


def test_inert_inputs_do_not_multiply_the_work():
    # A pool is one molecule per slot, drawn from the inputs carrying that slot's group, so the six inert
    # reagents are in no pool and cost nothing.
    plain_label, plain = _pools(_amidation())
    padded_label, padded = _pools(_amidation(_UNREACHABLE))
    assert plain_label == padded_label == ('react:amidation',)
    assert padded == plain, 'an input no row can reach must not be enumerated'


def test_a_plate_record_offers_each_row_one_pool_per_slot_choice():
    # One acid and a plate of ten amines, of which one reacted.  A pool is one molecule per slot, so the
    # acid-amine row is offered ten pools and the amine-amine urea row `C(10, 2)` pairs plus ten
    # symmetric ureas, each amine in excess.  Every subset of up to four inputs is `n**4/24` pools, which
    # on a 76-input record does not finish.
    from .. import _reconstruct

    offered = Counter()
    kept = _reconstruct._outcomes

    def counting(rule, molecules):
        offered[rule.id] += 1
        assert len(molecules) <= len(rule.groups)
        return kept(rule, molecules)

    amines = [smiles('C' * k + 'N') for k in range(2, 12)]
    rxn = ReactionContainer([smiles('CCC(=O)O'), *amines], [smiles('CCC(=O)NCC')])
    _reconstruct._outcomes = counting
    try:
        assert rxn.reconstruct_mapping() == ('react:amidation',)
    finally:
        _reconstruct._outcomes = kept
    assert offered['reactions:1'] == 10          # carboxylic_acid, primary_amine
    assert offered['reactions:182'] == 55        # primary_amine, primary_amine


def test_a_reagent_filed_once_fills_two_slots_as_an_excess():
    # The ester takes two methyls from one recorded MeMgBr: the first keeps the reagent's numbers, the
    # second equivalent is not in the record, so its carbon is left at 0 and the mapping stays 1-1.
    rxn = smiles('COC(=O)c1ccccc1.C[Mg]Br>>CC(C)(O)c1ccccc1')
    assert rxn.reconstruct_mapping() == ('react:grignard',)
    assert [r.rule for r in rxn.log if r.rule == 'reconstruct:excess'] == ['reconstruct:excess']
    product = rxn.products[0]
    numbers = [product.map_number_of(n) for n in product.atom_numbers]
    assert numbers.count(0) == 1
    assert len(set(numbers) - {0}) == len(numbers) - 1


def test_two_groups_on_one_input_are_one_molecule_not_two():
    # An amino acid fills an acid slot and an amine slot; that is never read as two molecules, so its
    # head-to-tail dimer is not an amidation of the monomer.
    rxn = smiles('NCC(=O)O>>NCC(=O)NCC(=O)O')
    assert rxn.reconstruct_mapping() == ()


def test_the_composed_rung_only_walks_what_the_stripping_reached():
    # `_deprotect_then_react` runs after `_react` on the same walk, so a pool of its stripped list that
    # holds nothing stripped is a pool `_react` already offered: here the acid alone, and the acid with
    # the isobutylene the Boc left behind.  Both rungs walk, so the count is per rung.
    from .._reconstruct import _applications

    rxn = ReactionContainer([smiles('CCC(=O)O'), smiles('CC(C)(C)OC(=O)NCC')], [smiles('CCC(=O)NCC')])
    stripped, unstripped = [], []
    for store, required in ((stripped, {1}), (unstripped, None)):
        pool = [smiles('CCC(=O)O'), smiles('CCN'), smiles('CC(C)=C')]
        store.extend(subset for subset, _ in _applications(rxn.products[0], pool, required=required))
    assert rxn.reconstruct_mapping() == ('deprotect:amine_boc+react:amidation',)
    assert len(stripped) < len(unstripped), 'a pool the stripping never reached must not be re-walked'


@pytest.mark.parametrize('reactants, product, label', [
    (('CC(C)(C)OC(=O)NCCC(=O)OC',), 'NCCC(=O)O', 'deprotect:amine_boc+ester_hydrolysis'),
    (('CC(C)(C)OC(=O)NCCC(=O)OC(C)(C)C',), 'NCCC(=O)O', 'deprotect:amine_boc+ester_hydrolysis'),
    (('COC(=O)CCC(=O)OC',), 'OC(=O)CCC(=O)O', 'deprotect:ester_hydrolysis+ester_hydrolysis'),
    (('COC(=O)c1ccc(Br)cc1', 'OB(O)c1ccccc1'), 'OC(=O)c1ccc(-c2ccccc2)cc1',
     'deprotect:ester_hydrolysis+react:suzuki'),
    (('COC(=O)c1ccccc1',), 'OC(=O)c1ccccc1', 'react:ester_hydrolysis'),
    (('CC(C)(C)OC(=O)c1ccccc1',), 'OC(=O)c1ccccc1', 'react:ester_hydrolysis'),
    (('COC(=O)c1ccc2cc[nH]c2n1', 'Ic1ccccc1'), 'OC(=O)c1ccc2ccn(-c3ccccc3)c2n1',
     'deprotect:ester_hydrolysis+react:ullmann_pyrrole'),
])
def test_an_ester_cleaved_to_its_acid_is_labelled_ester_hydrolysis(reactants, product, label):
    # the mechanism and not the protecting-group row: a methyl and a tert-butyl ester read alike
    rxn = ReactionContainer([smiles(x) for x in reactants], [smiles(product)])
    for molecule in rxn.molecules():
        molecule.canonicalize()
    assert rxn.reconstruct_mapping() == (label,)


def test_the_corpus_is_walked_once_for_both_strictnesses():
    # THE FLAT-PRODUCT RECORD, which only the loose walk explains: the strict walk enumerates every rung
    # and accepts none of it, so the loose walk is the one that earns the label.  What a rung enumerates
    # is a question about the inputs and the corpus, not about strictness, so it is asked once.
    from .. import _reconstruct

    rxn = ReactionContainer([smiles('C[C@H](O)CC'), smiles('CCC(=O)O')], [smiles('CC(CC)OC(=O)CC')])
    walks = []
    kept = _reconstruct._applications

    def counting(recorded, inputs, rules=None, required=None):
        walks.append(required)
        return kept(recorded, inputs, rules, required)

    _reconstruct._applications = counting
    try:
        assert rxn.reconstruct_mapping() == ('react:esterification',)
    finally:
        _reconstruct._applications = kept
    assert walks == [None], 'the corpus walk is paid once per call, not once per strictness'


def test_the_recorded_product_bounds_what_a_row_may_build():
    from .._reconstruct import _supplies
    from .._tables import reaction_rules

    rule = reaction_rules()['amidation'][0]
    created = max(len(template.created_atoms) for template in rule.templates)
    # the bound is over ELEMENTS, and one recorded component sufficing is enough
    assert _supplies(Counter({'C': 4, 'N': 1, 'O': 1}), [{'C': 4, 'N': 1, 'O': 1}], rule)
    assert _supplies(Counter({'C': 9}), [{'C': 40}, {'C': 4}], rule)
    # a record wanting an element no input carries is out of reach, up to what the row creates
    assert not _supplies(Counter({'C': 40}), [{'C': 4, 'Br': created + 1}], rule)
    # NECESSARY AND NOT SUFFICIENT: enough atoms of the right elements says nothing about the bonds
    assert _supplies(Counter({'C': 99}), [{'C': 4}], rule)


def test_an_unreachable_input_still_reports_as_a_spectator():
    # the pruning is of the SEARCH, so an input it skipped is still an input: it comes back at 0 like
    # any other spectator rather than disappearing from the record
    rxn = _amidation(_UNREACHABLE)
    assert rxn.reconstruct_mapping() == ('react:amidation',)
    assert len(rxn.reactants) + len(rxn.agents) == 2 + len(_UNREACHABLE)
    for molecule in (*rxn.reactants, *rxn.agents)[2:]:
        assert all(molecule.map_number_of(n) == 0 for n in molecule.atom_numbers)


# --- stereo strictness ----------------------------------------------------------------------------

def test_a_flat_product_is_explained_loosely():
    # THE SAD RECORD: a configured reactant and a product drawn flat.  Constitutionally the ester is
    # there, so the mapping is there, and the configuration the record failed to state is a log line.
    from ...core import INFO
    rxn = ReactionContainer([smiles('C[C@H](O)CC'), smiles('CCC(=O)O')], [smiles('CC(CC)OC(=O)CC')])
    assert rxn.reconstruct_mapping() == ('react:esterification',)
    mismatch = [r for r in rxn.log if r.rule == 'reconstruct:stereo-mismatch']
    assert len(mismatch) == 1 and mismatch[0].severity == INFO
    product = rxn.products[0]
    assert all(product.map_number_of(n) != 0 for n in product.atom_numbers
               if product.atom(n).element != 1)


def test_strict_refuses_what_loose_explains():
    from ...core import LOST
    rxn = ReactionContainer([smiles('C[C@H](O)CC'), smiles('CCC(=O)O')], [smiles('CC(CC)OC(=O)CC')])
    assert rxn.reconstruct_mapping(stereo='strict') == ()
    assert [r.rule for r in rxn.log if r.severity == LOST] == ['reconstruct:unexplained']


def test_the_strict_pass_runs_first_so_stereo_still_picks_the_row():
    # THE WHOLE REASON THE LOOSE PASS IS SECOND.  The corpus discriminates by configuration: one row
    # carries the centre through and another turns it over, and a record that states its product
    # configuration is entitled to the row that actually explains it.  Run blind, both rows fit.
    retention = ReactionContainer([smiles('C[C@H](O)CC'), smiles('CCC(=O)O')],
                                  [smiles('C[C@@H](CC)OC(=O)CC')])
    assert retention.reconstruct_mapping() == ('react:esterification',)
    inversion = ReactionContainer([smiles('C[C@H](O)CC'), smiles('CCC(=O)O')],
                                  [smiles('C[C@H](CC)OC(=O)CC')])
    assert inversion.reconstruct_mapping() == ('react:mitsunobu',)
    # and neither is a loose hit
    for rxn in (retention, inversion):
        assert not [r for r in rxn.log if r.rule == 'reconstruct:stereo-mismatch']


def test_an_inverted_purification_is_explained_loosely():
    # The lowest rung is stereo-strict too, so the enantiomer of the input is not the input; loosely it
    # is the same constitution and the record is a purification with a configuration disagreement.
    rxn = ReactionContainer([smiles('C[C@H](N)C(=O)O'), smiles('O')], [smiles('C[C@@H](N)C(=O)O')])
    assert rxn.reconstruct_mapping() == ('purification',)
    assert [r.rule for r in rxn.log if r.rule.startswith('reconstruct:stereo')] \
        == ['reconstruct:stereo-mismatch']


def test_the_mismatch_line_names_the_sites():
    rxn = ReactionContainer([smiles('C[C@H](N)C(=O)O'), smiles('O')], [smiles('C[C@@H](N)C(=O)O')])
    rxn.reconstruct_mapping()
    mismatch = next(r for r in rxn.log if r.rule == 'reconstruct:stereo-mismatch')
    product = rxn.products[0]
    centre = next(n for n in product.atom_numbers if product.parity_of(n))
    assert mismatch.atoms == (centre,)


def test_a_constitutional_mismatch_is_still_unexplained():
    # Loose relaxes the configuration and nothing else: the hydrogen counts, the kekule form and the
    # protonation are still what `isomorphism()` requires.
    rxn = ReactionContainer([smiles('[CH3:1][CH2:2][OH:3]')], [smiles('[CH3:4][CH2:5][CH3:6]')])
    assert rxn.reconstruct_mapping() == ()


def test_the_stereo_argument_is_validated():
    rxn = ReactionContainer([smiles('CCO')], [smiles('CCO')])
    with raises(ValueError, match="'loose' or 'strict'"):
        rxn.reconstruct_mapping(stereo='blind')


# --- healing the configuration --------------------------------------------------------------------

def _cip(molecule, number: int):
    """The CIP descriptor at the atom carrying map number `number`, or None where there is none.

    THE PARITY INTEGER IS NOT THE ASSERTION.  It is read against the atom's own refs, and a row that
    changes a substituent changes those, so a centre carried through retentively legitimately shows two
    different integers on the two sides of one record.
    """
    molecule.assign_cip()
    return next((d for n, d in molecule.atom_cips().items() if molecule.map_number_of(n) == number),
                None)


def _healed(record: str, **kwargs):
    """`reconstruct_mapping(**kwargs)` on `record`, then the label the STRICT walk gives the result.

    The ratchet for every repair: a healed record is one the strict walk accepts, which is a statement
    about the whole record rather than about one integer in one frame.
    """
    rxn = smiles(record)
    label = rxn.reconstruct_mapping(**kwargs)
    return rxn, label, smiles(format(rxn, 'm')).reconstruct_mapping(stereo='strict')


def test_heal_is_off_by_default():
    rxn, label, strict = _healed('CC(O)CC.CCC(=O)O>>C[C@H](CC)OC(=O)CC')
    assert label == ('react:esterification',)
    assert not [r for r in rxn.log if r.stage == 'heal']
    assert strict == ()                                   # the record is as inconsistent as it arrived


def test_the_product_configuration_fills_a_flat_input():
    # THE ORDINARY SAD RECORD, and the reason `'product'` is the direction a caller reaches for: the
    # product was copied from a catalogue with its configuration and the starting material drawn by hand
    # without one.  Fischer esterification never touches the carbinol bond, so the record's own centre is
    # the input's.
    rxn, label, strict = _healed('CC(O)CC.CCC(=O)O>>C[C@H](CC)OC(=O)CC', heal_stereo='product')
    assert label == ('react:esterification',) and strict == ('react:esterification',)
    alcohol = next(m for m in rxn.reactants if len(m) == 5)
    assert _cip(alcohol, 2) == _cip(rxn.products[0], 2) == 'R'
    assert [r.rule for r in rxn.log if r.stage == 'heal'] == ['heal:parity']


def test_the_other_enantiomer_fills_it_the_other_way():
    # The same record drawn S, to pin that the repair reads the record rather than a default.
    rxn, label, strict = _healed('CC(O)CC.CCC(=O)O>>C[C@@H](CC)OC(=O)CC', heal_stereo='product')
    assert strict == label == ('react:esterification',)
    alcohol = next(m for m in rxn.reactants if len(m) == 5)
    assert _cip(alcohol, 2) == _cip(rxn.products[0], 2) == 'S'


def test_the_product_configuration_turns_an_input_over():
    # The input states a configuration and the record contradicts it.  Every unit is two-state, so the
    # repair is the other state -- no inversion arithmetic and no second code path.
    rxn, label, strict = _healed('C[C@H](N)C(=O)O.O>>C[C@@H](N)C(=O)O', heal_stereo='product')
    assert label == ('purification',) and strict == ('purification',)
    alanine = next(m for m in rxn.reactants if len(m) > 1)
    assert _cip(alanine, 2) == _cip(rxn.products[0], 2)
    assert [r.rule for r in rxn.log if r.stage == 'heal'] == ['heal:parity']


def test_the_input_configuration_fills_a_flat_product():
    # `'reactant'`: the other direction, and the one that STATES a configuration where the record holds
    # none rather than turning one over.
    rxn, label, strict = _healed('C[C@H](O)CC.CCC(=O)O>>CC(CC)OC(=O)CC', heal_stereo='reactant')
    assert label == ('react:esterification',) and strict == ('react:esterification',)
    alcohol = next(m for m in rxn.reactants if len(m) == 5)
    assert _cip(rxn.products[0], 2) == _cip(alcohol, 2) == 'S'


def test_the_two_directions_repair_opposite_sides():
    # One record, both directions, and the answers differ: the geometry the record states in its product
    # against the geometry it states in its reactant.  Either way the record ends up consistent.
    record = r'C/C=C/CO.CCC(=O)O>>C/C=C\COC(=O)CC'
    from_product, label, strict = _healed(record, heal_stereo='product')
    assert label == ('react:esterification',) and strict == ('react:esterification',)
    from_reactant, label, strict = _healed(record, heal_stereo='reactant')
    assert label == ('react:esterification',) and strict == ('react:esterification',)
    assert format(from_product, 'm') != format(from_reactant, 'm')


def test_nothing_is_healed_where_the_record_states_nothing():
    # `'product'` onto a flat product asks for a configuration the record does not hold.  Silence and no
    # log line: there is no disagreement here, only an absence on the side that was named the truth.
    rxn, label, strict = _healed('C[C@H](O)CC.CCC(=O)O>>CC(CC)OC(=O)CC', heal_stereo='product')
    assert label == ('react:esterification',)
    assert not [r for r in rxn.log if r.stage == 'heal']
    assert strict == ()


def test_a_row_that_rebuilt_the_centre_refuses():
    # THE REFUSAL THAT MAKES THE REST SOUND.  An N-alkylation replaces the bond at the centre, so the
    # product's configuration says nothing about the halide's: a displacement turns the centre over and
    # the row does not state that it does.  Nothing is written and the line says what is missing.
    rxn, label, strict = _healed('CC(Cl)CC.CNC>>C[C@H](CC)N(C)C', heal_stereo='product')
    assert label == ('react:n_alkylation',)
    assert [r.rule for r in rxn.log if r.stage == 'heal'] == ['heal:course-unstated']
    halide = next(m for m in rxn.reactants if len(m) == 5)
    assert all(halide.parity_of(n) == 0 for n in halide.atom_numbers)
    assert strict == ()


def test_healing_is_frame_order_independent():
    # ONE RECORD, TWO ATOM ORDERS, one answer.  The parity integer differs between the two spellings and
    # the map numbers land on different atoms, so neither is the assertion: the healed alcohol is the
    # same molecule, `__eq__` being the canonical form and carrying the parity.
    one = _healed('CC(O)CC.CCC(=O)O>>C[C@H](CC)OC(=O)CC', heal_stereo='product')[0]
    two = _healed('CCC(O)C.OC(=O)CC>>CC[C@@H](C)OC(=O)CC', heal_stereo='product')[0]
    alcohols = [next(m for m in rxn.reactants if len(m) == 5) for rxn in (one, two)]
    assert alcohols[0] == alcohols[1]
    for alcohol in alcohols:
        alcohol.assign_cip()
        assert list(alcohol.atom_cips().values()) == ['R']


def test_every_heal_line_is_stage_heal():
    rxn = _healed('CC(O)CC.CCC(=O)O>>C[C@H](CC)OC(=O)CC', heal_stereo='product')[0]
    healed = [r for r in rxn.log if r.rule.startswith('heal:')]
    assert healed and all(r.stage == 'heal' for r in healed)
    assert all(r.stage == 'reconstruct' for r in rxn.log if not r.rule.startswith('heal:'))


def test_the_heal_line_names_the_atoms_of_the_molecule_it_wrote():
    rxn = _healed('CC(O)CC.CCC(=O)O>>C[C@H](CC)OC(=O)CC', heal_stereo='product')[0]
    line = next(r for r in rxn.log if r.rule == 'heal:parity')
    alcohol = next(m for m in rxn.reactants if len(m) == 5)
    assert line.atoms == tuple(n for n in alcohol.atom_numbers if alcohol.parity_of(n))


def test_heal_stereo_is_validated():
    rxn = ReactionContainer([smiles('CCO')], [smiles('CCO')])
    with raises(ValueError, match="False, 'product' or 'reactant'"):
        rxn.reconstruct_mapping(heal_stereo='both')


def test_a_collection_crosses_under_a_fresh_id():
    # A racemate is a statement about a configured unit, so it travels with the configuration.  The KIND
    # crosses and the NUMBER does not: two records numbering their `&1` differently must not have their
    # racemates merged by a heal.
    rxn, label, _ = _healed('C[C@H](O)CC.CCC(=O)O>>C[C@H](CC)OC(=O)CC |&1:1|', heal_stereo='reactant')
    assert label == ('react:mitsunobu',)
    product = rxn.products[0]
    assert [kind for kind, _ in product.stereo_groups()] \
        == [kind for kind, _ in rxn.reactants[0].stereo_groups()]
    assert [r.rule for r in rxn.log if r.stage == 'heal'] == ['heal:group']


def test_an_absolute_source_drops_a_collection():
    # The other way round, and the policy the argument states: the side the caller named wins outright.
    # A racemic input against a single-enantiomer product, read from the product, is absolute -- and the
    # line is there to be read, because this is the direction that discards a statement.
    rxn, label, _ = _healed('C[C@H](O)CC.CCC(=O)O>>C[C@H](CC)OC(=O)CC |&1:1|', heal_stereo='product')
    assert label == ('react:mitsunobu',)
    assert rxn.reactants[0].stereo_groups() == {}
    line = next(r for r in rxn.log if r.rule == 'heal:group')
    assert 'dropped' in line.message


def test_an_ester_hydrolysis_mapping_has_a_reaction_centre():
    # Both oxygens the row swaps are unmapped -- the alkoxy leaves, the OH arrives from an agent -- so the
    # centre is the acyl carbon's unmapped branch changing, and the reagents separate on it.
    rxn = ReactionContainer([smiles(x) for x in ('COC(=O)c1ccccc1', '[Li+]', '[OH-]', 'O')],
                            [smiles('OC(=O)c1ccccc1')])
    assert rxn.reconstruct_mapping(stereo='loose', heal_stereo='product') == ('react:ester_hydrolysis',)
    assert rxn.remove_reagents(keep_reagents=True) is True
    assert [m.smiles for m in rxn.reactants] == [smiles('COC(=O)c1ccccc1').smiles]
    assert sorted(m.smiles for m in rxn.agents) == sorted(smiles(x).smiles for x in ('[Li+]', '[OH-]', 'O'))


def test_a_counterion_carried_through_a_row_is_not_an_explanation():
    # Each row fires on a salt input and carries its metal through untouched; the recorded organic ion is
    # never built, so the metal matching itself explains nothing.
    for inputs, product in ((('COC(=O)c1ccc(F)cc1[N+](=O)[O-]', 'CC(C)(C)[O-].[K+]', 'O'),       # SNAr
                             '[O-]C(=O)c1ccc(F)cc1[N+](=O)[O-].[K+]'),
                            (('CCNCC', 'CCOS(=O)(=O)c1ccc(C)cc1.[K+].[Cl-]'),                  # nitration
                             'OS(=O)(=O)c1ccc(C)cc1.[K+]'),
                            (('CC(C)(C)OC(=O)CCC(=O)O', '[Ag+].[O-][N+](=O)[O-]', '[Na+].[OH-]'),  # protect
                             'CC(C)(C)OC(=O)CCC(=O)[O-].[Ag+]')):
        rxn = ReactionContainer([smiles(x) for x in inputs], [smiles(product)])
        assert rxn.reconstruct_mapping() == (), product


def test_a_protected_azole_pairs_with_its_input_tautomer():
    # The strip leaves the unmasked pyrazole nitrogen where the patch put it; settled, it is the input.
    for inputs, product, label in ((('Cc1cc[nH]n1', 'C1CCOC=C1', 'Cc1ccc(cc1)S(=O)(=O)O'),
                                    'Cc1ccn(n1)C1CCCCO1', 'protect:azole_thp'),
                                   (('Cc1cc[nH]n1', 'ClCOCC[Si](C)(C)C', '[Na+].[H-]'),
                                    'Cc1ccn(COCC[Si](C)(C)C)n1', 'protect:azole_sem')):
        rxn = ReactionContainer([smiles(x) for x in inputs], [smiles(product)])
        assert rxn.reconstruct_mapping() == (label,)
        assert rxn.remove_reagents(keep_reagents=True) is True


@pytest.mark.parametrize('record, label', [
    ('COC(=O)CC1=CCN(C(=O)OCc2ccccc2)CC1>>COC(=O)CC1CCNCC1', 'deprotect:amine_cbz+react:alkene_to_alkane'),
    ('OC(=O)CCCC(C(=O)OC)NC(=O)OC(C)(C)C>>O=C1CCCC(C(=O)OC)N1', 'deprotect:amine_boc+react:amidation'),
    ('COC(=O)CCC(=O)CCCCNC(=O)OC(C)(C)C>>COC(=O)CCC1=NCCCC1', 'deprotect:amine_boc+react:imine_condensation'),
    ('CC(C)(C)OC(=O)N1CCC(CC1)C=Cc1ccc(OCc2ccccc2)cc1>>CC(C)(C)OC(=O)N1CCC(CC1)CCc1ccc(O)cc1',
     'deprotect:phenol_benzyl+react:alkene_to_alkane'),
    ('CS(=O)(=O)OCCc1ccccc1.Oc1ccc(NC(=O)OC(C)(C)C)cc1>>Nc1ccc(OCCc2ccccc2)cc1',
     'deprotect:aniline_boc+react:williamson'),     # the mesylate is the leaving group, not a strip
])
def test_a_group_the_product_keeps_is_not_stripped(record, label):
    """The composed rung strips the groups the product lost, not every group it could: the ester or the Boc
    the product still carries stays on."""
    assert smiles(record).reconstruct_mapping() == (label,)


@pytest.mark.parametrize('record, label', [
    ('OC(=O)c1ccccc1.CO>>COC(=O)c1ccccc1', 'protect:carboxyl_methyl'),
    ('Oc1ccccc1.CI>>COc1ccccc1', 'protect:phenol_methyl'),
    ('OC(=O)c1ccccc1.BrCc1ccccc1>>O=C(OCc1ccccc1)c1ccccc1', 'protect:carboxyl_benzyl'),
    ('OCc1ccccc1.CC(=O)OC(C)=O>>CC(=O)OCc1ccccc1', 'protect:hydroxyl_acyl'),
])
def test_an_ester_installed_on_the_acid_is_named_by_the_acid(record, label):
    """One `hydroxyl_*` row reads an ether and an ester alike; the installs run under different conditions,
    so the ester made on the acid's oxygen is `carboxyl_*`."""
    assert smiles(record).reconstruct_mapping() == (label,)


TF2O, PHNTF2 = 'O=S(=O)(OS(=O)(=O)C(F)(F)F)C(F)(F)F', 'O=S(=O)(N(c1ccccc1)S(=O)(=O)C(F)(F)F)C(F)(F)F'


@pytest.mark.parametrize('record, label', [
    (f'Oc1ccccc1.{TF2O}>>O=S(=O)(Oc1ccccc1)C(F)(F)F', 'protect:phenol_triflyl'),
    (f'Oc1ccccc1.{PHNTF2}>>O=S(=O)(Oc1ccccc1)C(F)(F)F', 'protect:phenol_triflyl'),
    ('OCCc1ccccc1.CS(=O)(=O)Cl>>CS(=O)(=O)OCCc1ccccc1', 'protect:hydroxyl_mesyl'),
    ('OCCc1ccccc1.CS(=O)(=O)OS(C)(=O)=O>>CS(=O)(=O)OCCc1ccccc1', 'protect:hydroxyl_mesyl'),
    ('OCCc1ccccc1.Cc1ccc(cc1)S(=O)(=O)Cl>>Cc1ccc(cc1)S(=O)(=O)OCCc1ccccc1', 'protect:hydroxyl_tosyl'),
    (f'O=C1CCCCC1.{PHNTF2}>>O=S(=O)(OC1=CCCCC1)C(F)(F)F', 'protect:enol_triflyl'),
    ('O=S(=O)(Oc1ccccc1)C(F)(F)F.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1', 'react:suzuki'),
    ('CS(=O)(=O)OCCc1ccccc1.N1CCCCC1>>C(CN1CCCCC1)c1ccccc1', 'react:n_alkylation'),
])
def test_an_o_sulfonyl_pseudohalide_is_installed_as_a_group(record, label):
    """Tf, Ms and Ts on an OH are set by many reagents, so they are groups and the reagent an agent; the
    enol triflate of a ketone is the same group on the `enol` site, and using the pseudohalide is a reaction."""
    assert smiles(record).reconstruct_mapping() == (label,)


@pytest.mark.parametrize('record, label', [
    ('NCCc1ccccc1.CS(=O)(=O)Cl>>CS(=O)(=O)NCCc1ccccc1', 'protect:amine_mesyl'),
    (f'NCCc1ccccc1.{TF2O}>>O=S(=O)(NCCc1ccccc1)C(F)(F)F', 'protect:amine_triflyl'),
    ('NCCc1ccccc1.Cc1ccc(cc1)S(=O)(=O)Cl>>Cc1ccc(cc1)S(=O)(=O)NCCc1ccccc1', 'protect:amine_tosyl'),
    ('OCCc1ccccc1.O=S(=O)(Cl)c1ccccc1>>O=S(=O)(OCCc1ccccc1)c1ccccc1', 'protect:hydroxyl_phenylsulfonyl'),
    ('OCCc1ccccc1.O=S(=O)(Cl)c1ccccc1[N+](=O)[O-]>>O=S(=O)(OCCc1ccccc1)c1ccccc1[N+](=O)[O-]',
     'protect:hydroxyl_nosyl'),
    ('NCCc1ccccc1.O=S(=O)(Cl)c1ccc(F)cc1>>O=S(=O)(NCCc1ccccc1)c1ccc(F)cc1', 'react:sulfonamide_formation'),
])
def test_the_named_sulfonyl_groups_are_the_same_on_oxygen_and_nitrogen(record, label):
    """Ms, Tf, Ts, Ns and PhSO2 are groups on O and on N alike; a sulfonyl no row names is a reaction."""
    assert smiles(record).reconstruct_mapping() == (label,)


BOC2O = 'CC(C)(C)OC(=O)OC(=O)OC(C)(C)C'


@pytest.mark.parametrize('record, label', [
    (f'NCCc1ccccc1.{BOC2O}>>CC(C)(C)OC(=O)NCCc1ccccc1', 'protect:amine_boc'),
    (f'Nc1ccccc1.{BOC2O}>>CC(C)(C)OC(=O)Nc1ccccc1', 'protect:aniline_boc'),
    (f'O=C1CCCN1.{BOC2O}>>CC(C)(C)OC(=O)N1CCCC1=O', 'protect:amide_boc'),
    ('c1cn[nH]c1.C[Si](C)(C)CCOCCl>>C[Si](C)(C)CCOCn1cccn1', 'protect:azole_sem'),
    (f'Oc1ccccc1.{TF2O}>>O=S(=O)(Oc1ccccc1)C(F)(F)F', 'protect:phenol_triflyl'),
    (f'OCCc1ccccc1.{TF2O}>>O=S(=O)(OCCc1ccccc1)C(F)(F)F', 'protect:hydroxyl_triflyl'),
    ('OC(=O)c1ccccc1.CO>>COC(=O)c1ccccc1', 'protect:carboxyl_methyl'),
    ('COc1ccccc1>>Oc1ccccc1', 'deprotect:phenol_methyl'),
    ('CC(C)(C)OC(=O)Nc1ccccc1>>Nc1ccccc1', 'deprotect:aniline_boc'),
    ('CC(C)(C)OC(=O)NCCc1ccccc1>>NCCc1ccccc1', 'deprotect:amine_boc'),
    ('COC1=NCCCC1>>O=C1CCCCN1', 'deprotect:enol_methyl'),               # a lactim ether: the O sits on C=N
    ('CCOC=Cc1ccccc1>>O=CCc1ccccc1', 'deprotect:enol_ethyl'),           # an enol ether: the O sits on C=C
    ('COc1ccccn1>>O=c1cccc[nH]1', 'deprotect:phenol_methyl'),           # an alkoxyazine, the pyridone back
    ('CC(C)(C)[Si](C)(C)OC1=CCCCC1>>O=C1CCCCC1', 'deprotect:enol_tbs'),  # a ring enol, the ketone back
])
def test_a_group_is_named_for_the_atom_it_sits_on(record, label):
    """One group on an aniline, an amide N-H, an azole or a phenol is installed and removed under other
    conditions than on an alkylamine or an alcohol, so the label names the site: `aniline_boc`."""
    assert smiles(record).reconstruct_mapping() == (label,)
