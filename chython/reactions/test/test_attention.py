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
"""`ReactionContainer.attention_mapping()` end to end, and its agreement with a reference corpus.

MOST OF THIS FILE NEEDS THE MODEL and skips without it.  What does not is kept at the top: the
refusals -- an empty side and a hypervalent atom -- are decided before the weights are reached, and the
unregistered accessor is a question about the injection hook.

THE AGREEMENT FLOOR IS MEASURED, NOT ASSERTED FROM A DESIGN.  `test_the_golden_subset_agrees` runs 25
records of the public Golden benchmark and holds the numbers the first green run produced.  A floor
rather than an equality, so a better model or a better walk passes; a regression does not.

THE COMPARISON RUNS ON THE AROMATIC FORM, and that is not cosmetic.  `mapping_agrees` excuses a
disagreement when the two candidate atoms lie in one automorphism orbit, and a Kekulé ring has no
mirror automorphism -- alternating bond orders distinguish the two ortho carbons that the aromatic form
makes equivalent.  Measured on this subset, comparing the Kekulé structures scores 12 of 25 records
exact and comparing the aromatic ones scores 23, on identical mappings.  The 11 records between those
two numbers are ring-direction flips, which are the same mapping.
"""
from pathlib import Path

from pytest import mark, raises, skip

from ...core import MoleculeContainer, ReactionContainer, read_reaction_smiles, read_smiles as smiles
from .._numbering import mapping_agrees
from ..attention import attention_available


#: 25 records of the Golden benchmark, `index<TAB>reaction SMILES`.  Ships with the package: the
#: agreement claim above is only checkable where the corpus is.
GOLDEN = Path(__file__).resolve().parent / 'golden_subset.smi'

#: What the first green run produced, on the aromatic form.  Floors, not equalities.
EXACT_FLOOR = 23
AGREED_FLOOR = 551

#: Product atoms in the subset, and how many of them `mapping_agrees` can score.  THE TWO DIFFER BY ONE:
#: record 1258's reference gives a product nitrogen the number 19, which no atom on its reactant side
#: carries, so that atom has no `(input index, input atom id)` to compare and neither side scores it.
#: An incomplete reference is a property of the corpus, not something to repair here.
TOTAL_PRODUCT_ATOMS = 556
COMPARABLE_PRODUCT_ATOMS = 555


def _records():
    """`(index, reaction SMILES)` per corpus line, comments dropped."""
    out = []
    for line in GOLDEN.read_text(encoding='utf-8').splitlines():
        line = line.strip()
        if line and not line.startswith('#'):
            index, record = line.split('\t')
            out.append((int(index), record))
    return out


def _needs_model():
    if not attention_available():
        skip('needs `chython[mapping]`: onnxruntime and chython_rxnmap')


# --- decided before the model is reached ---------------------------------------------------------- #

@mark.parametrize('reactants, products', [(['CCO'], []), ([], ['CCO']), ([], [])])
def test_an_empty_side_is_refused_and_nothing_is_written(reactants, products):
    """No molecule on one side means no correspondence to find, and the refusal is logged as one."""
    rxn = ReactionContainer([smiles(s) for s in reactants], [smiles(s) for s in products])
    result = rxn.attention_mapping()
    assert result.skipped == 'empty'
    assert not result
    assert result.score == 0. and result.unplaced == ()

    refused = rxn.log.refused()
    assert [entry.rule for entry in refused] == ['attention:empty']
    assert refused[0].stage == 'attention_mapping'


def test_an_atom_past_fourteen_heavy_neighbours_is_refused():
    """Outside the domain the weights were trained on, so no map number is written at all.

    HEAVY DEGREE ALONE, and the 14 here is not the 14 the `neighbors` column clamps at: that clamp is on
    degree plus hydrogens and merely saturates a token the weights have seen, while this atom has no
    token at all.  A record with a 15-coordinate atom is accepted by the reader -- input is garbage by
    default -- so the mapper is where it has to be declined.
    """
    wide = MoleculeContainer()
    with wide.edit() as e:
        centre = e.add_atom('W')
        for _ in range(15):
            e.add_bond(centre, e.add_atom('C'), 1)
    assert wide.degree_of(centre) == 15

    rxn = ReactionContainer([wide], [wide.copy()])
    result = rxn.attention_mapping()
    assert result.skipped == 'hypervalent'
    assert not result
    assert all(not molecule.map_number_of(n) for molecule in rxn.molecules()
               for n in molecule.atom_numbers), 'a refusal wrote map numbers'

    refused = rxn.log.refused()
    assert [entry.rule for entry in refused] == ['attention:hypervalent']
    assert '15' not in refused[0].message or '14' in refused[0].message


def test_fourteen_heavy_neighbours_is_still_accepted():
    """The bound is exclusive, and a test that only proved the refusal would not pin which side."""
    _needs_model()
    wide = MoleculeContainer()
    with wide.edit() as e:
        centre = e.add_atom('W')
        for _ in range(14):
            e.add_bond(centre, e.add_atom('C'), 1)

    rxn = ReactionContainer([wide], [wide.copy()])
    assert rxn.attention_mapping().skipped is None


def test_the_accessor_names_the_package_when_unregistered():
    """The core owns the method name; `chython.reactions` registers the body at its own import."""
    from ...core._core import _reaction_attention_fn, _set_attention_fn
    try:
        kept = _reaction_attention_fn()          # another test module may have registered it already
    except ImportError:
        kept = None
    _set_attention_fn(None)
    try:
        with raises(ImportError, match='chython.reactions'):
            _reaction_attention_fn()
    finally:
        _set_attention_fn(kept)


# --- the model ------------------------------------------------------------------------------------ #

def test_an_amidation_maps_the_way_the_chemistry_reads():
    """One record asserted PER ATOM, because a count of placed atoms holds for a wrong mapping too."""
    _needs_model()
    rxn = read_reaction_smiles('CC(=O)O.CCN>>CC(=O)NCC.O')
    result = rxn.attention_mapping()
    assert result
    assert result.skipped is None and result.unplaced == ()

    # every atom, by the pair `(input index, input atom id)` it came from
    acid, amine = rxn.reactants
    amide, water = rxn.products
    source = {}
    for index, molecule in enumerate(rxn.reactants):
        for n in molecule.atom_numbers:
            source[molecule.map_number_of(n)] = (index, n)

    # the amine's nitrogen becomes the amide's nitrogen, and its two carbons stay its two carbons
    nitrogen = next(n for n in amine.atom_numbers if amine.atom(n).element == 7)
    amide_n = next(n for n in amide.atom_numbers if amide.atom(n).element == 7)
    assert source[amide.map_number_of(amide_n)] == (1, nitrogen)

    # the water is one of the acid's two oxygens -- the leaving group, and it comes from the acid
    water_o = water.atom_numbers[0]
    assert source[water.map_number_of(water_o)][0] == 0, 'the water did not come from the acid'
    assert acid.atom(source[water.map_number_of(water_o)][1]).element == 8


def test_the_mapping_written_is_a_partial_injection():
    """Two product atoms may not share a number, and every number used names a reactant atom.

    The property the greedy walk exists to hold, asserted over a record big enough for it to fail on.
    """
    _needs_model()
    rxn = read_reaction_smiles('CC(=O)Oc1ccccc1C(=O)O.O>>CC(=O)O.Oc1ccccc1C(=O)O')
    assert rxn.attention_mapping()

    reactant_numbers = {molecule.map_number_of(n) for molecule in rxn.reactants
                        for n in molecule.atom_numbers}
    assert 0 not in reactant_numbers, 'a reactant atom was left unnumbered'
    assert reactant_numbers == set(range(1, len(reactant_numbers) + 1)), 'reactants are not 1..N'

    written = [molecule.map_number_of(n) for molecule in rxn.products
               for n in molecule.atom_numbers if molecule.map_number_of(n)]
    assert len(written) == len(set(written)), 'two product atoms share one map number'
    assert set(written) <= reactant_numbers, 'a product atom names a number no reactant carries'


def test_the_elements_of_a_correspondence_always_match():
    """A carbon never maps to an oxygen, whatever the attention says -- the encoder's equality mask."""
    _needs_model()
    rxn = read_reaction_smiles('CC(=O)OCC.O>>CC(=O)O.CCO')
    assert rxn.attention_mapping()

    elements = {molecule.map_number_of(n): molecule.atom(n).element
                for molecule in rxn.reactants for n in molecule.atom_numbers}
    for molecule in rxn.products:
        for n in molecule.atom_numbers:
            number = molecule.map_number_of(n)
            if number:
                assert elements[number] == molecule.atom(n).element, \
                    f'atom {n} of a product maps to an atom of a different element'


def test_nothing_but_the_map_numbers_changes():
    """The structures, the ids and the bonds come back as they went in."""
    _needs_model()
    rxn = read_reaction_smiles('CC(=O)O.CCN>>CC(=O)NCC.O')
    before = [(molecule.atom_numbers, str(molecule)) for molecule in rxn.molecules()]
    assert rxn.attention_mapping()
    after = [(molecule.atom_numbers, str(molecule)) for molecule in rxn.molecules()]
    assert before == after


def test_a_product_atom_with_no_counterpart_keeps_zero_and_is_reported():
    """A bromine appearing only on the product side has nothing to map to.

    KEEPING 0 RATHER THAN TAKING A FRESH NUMBER is the point: a number above the reactant range would
    say "this atom is new", which is a claim the model never made -- it said nothing about this atom.
    The `unplaced` tuple is where the caller reads that, and the log records it as a loss.
    """
    _needs_model()
    rxn = ReactionContainer([smiles('CCO')], [smiles('CCBr')])
    result = rxn.attention_mapping()
    assert result

    product = rxn.products[0]
    bromine = next(n for n in product.atom_numbers if product.atom(n).element == 35)
    assert product.map_number_of(bromine) == 0
    assert result.unplaced == ((0, bromine),)
    assert [entry.rule for entry in rxn.log if entry.severity == 'lost'] == ['attention:unplaced']


def test_the_score_is_recorded_on_the_container_as_well_as_returned():
    """`molecule.log` is the one destination, so the number the caller reads is also written down."""
    _needs_model()
    rxn = read_reaction_smiles('CC(=O)O.CCN>>CC(=O)NCC.O')
    result = rxn.attention_mapping()
    scored = [entry for entry in rxn.log if entry.rule == 'attention:score']
    assert len(scored) == 1
    assert ('%.3f' % result.score) in scored[0].message


def test_mapping_twice_reports_no_second_change():
    """`changed` is measured against the numbers the record carried, so an idempotent run says so."""
    _needs_model()
    rxn = read_reaction_smiles('CC(=O)O.CCN>>CC(=O)NCC.O')
    assert rxn.attention_mapping().changed
    again = rxn.attention_mapping()
    assert not again.changed
    assert again.score > 0., 'the model still ran; only the write was a no-op'


def test_agents_are_numbered_last_and_never_modelled():
    """A catalyst gets a number above the reactant range, and its presence changes nothing else."""
    _needs_model()
    rxn = ReactionContainer([smiles('CCO')], [smiles('CC=O')])
    assert rxn.attention_mapping()
    without = [molecule.map_number_of(n) for molecule in rxn.molecules()
               for n in molecule.atom_numbers]

    withal = ReactionContainer([smiles('CCO')], [smiles('CC=O')], [smiles('[Pd]')])
    assert withal.attention_mapping()
    palladium = withal.agents[0]
    assert withal.agents[0].map_number_of(palladium.atom_numbers[0]) == 4, 'three reactant atoms first'
    assert [molecule.map_number_of(n) for molecule in (*withal.reactants, *withal.products)
            for n in molecule.atom_numbers] == without


def test_keep_reactant_mapping_leaves_the_reactant_numbers_alone():
    """For a record whose inputs are already mapped by something else."""
    _needs_model()
    rxn = read_reaction_smiles('[CH3:5][C:6](=[O:7])[OH:8].CCN>>CC(=O)NCC.O')
    acid = next(m for m in rxn.reactants if len(m) == 4)
    kept = {n: acid.map_number_of(n) for n in acid.atom_numbers}
    assert sorted(kept.values()) == [5, 6, 7, 8]

    assert rxn.attention_mapping(keep_reactant_mapping=True)
    assert {n: acid.map_number_of(n) for n in acid.atom_numbers} == kept

    # A REACTANT ATOM CARRYING NO NUMBER STILL GETS ONE, above the highest kept: leaving the hole and
    # numbering the product atom matched to it would write a correspondence to nothing.
    amine = next(m for m in rxn.reactants if len(m) != 4)
    assert sorted(amine.map_number_of(n) for n in amine.atom_numbers) == [9, 10, 11]

    written = {molecule.map_number_of(n) for molecule in rxn.products for n in molecule.atom_numbers}
    assert written <= set(range(5, 12)), 'a product atom names a number no reactant carries'


def test_the_multiplier_is_a_knob_and_not_a_constant():
    """Two values of `multiplier` reach different mappings on a record where the choice is close.

    Otherwise the argument is decoration.  WHICH mapping is better is not asserted -- that is what the
    Golden agreement below measures; this pins only that the walk actually reads the value.
    """
    _needs_model()
    record = 'CC(=O)Oc1ccccc1C(=O)O.O>>CC(=O)O.Oc1ccccc1C(=O)O'

    def numbers(multiplier):
        rxn = read_reaction_smiles(record)
        rxn.attention_mapping(multiplier=multiplier)
        return [molecule.map_number_of(n) for molecule in rxn.products for n in molecule.atom_numbers]

    assert numbers(1.75) != numbers(1.), 'the neighbourhood bonus changed no outcome at all'


def test_a_thread_count_changes_nothing_about_the_answer():
    """It is a runtime setting.  A second value costs a second loaded model and the same mapping."""
    _needs_model()
    def numbers(threads):
        rxn = read_reaction_smiles('CC(=O)O.CCN>>CC(=O)NCC.O')
        rxn.attention_mapping(threads=threads)
        return [molecule.map_number_of(n) for molecule in rxn.molecules()
                for n in molecule.atom_numbers]

    assert numbers(1) == numbers(2)


# --- the corpus ----------------------------------------------------------------------------------- #

def test_the_golden_subset_agrees():
    """25 public reference-mapped records, scored per record and per atom against their own mapping.

    Compared through `mapping_agrees` and never through container equality: `__eq__` excludes map
    numbers, so `probe == reference` holds for every record and would report a perfect score while
    measuring nothing.

    `thiele()` FIRST, ON BOTH SIDES.  It is what makes a ring's mirror symmetry an automorphism, which
    is what lets the orbit excuse apply to a ring-direction flip; see this module's docstring for the
    11 records it accounts for.  It is applied to the reference before the copy, so the mapper and the
    comparison see one structure.
    """
    _needs_model()
    exact = declined = agreed = disagreed = missing = 0
    imperfect = []
    for index, record in _records():
        reference = read_reaction_smiles(record)
        reference.thiele()
        probe = reference.copy()
        if not probe.attention_mapping():
            declined += 1
            continue
        a, d, m = mapping_agrees(probe, reference)
        agreed += a
        disagreed += d
        missing += m
        if d or m:
            imperfect.append('%d (%d agreed, %d disagreed, %d missing)' % (index, a, d, m))
        else:
            exact += 1

    assert not declined, 'the mapper declined a record of the benchmark'
    assert agreed + disagreed + missing == COMPARABLE_PRODUCT_ATOMS, 'the corpus changed size'
    assert exact >= EXACT_FLOOR, (
        'exact records fell to %d of 25 (floor %d).  Not exact:\n  %s'
        % (exact, EXACT_FLOOR, '\n  '.join(imperfect)))
    assert agreed >= AGREED_FLOOR, (
        'product atoms agreeing fell to %d of %d (floor %d)'
        % (agreed, COMPARABLE_PRODUCT_ATOMS, AGREED_FLOOR))


def test_the_corpus_file_is_the_size_the_floors_were_measured_on():
    """A floor is a number about a corpus, so the corpus is pinned too.

    Without this, dropping a record the mapper gets wrong raises the score and the floors still pass.
    """
    records = _records()
    assert len(records) == 25
    assert sum(len(m) for _, record in records
               for m in read_reaction_smiles(record).products) == TOTAL_PRODUCT_ATOMS
    assert [index for index, _ in records] == list(range(0, 25 * 74, 74)), \
        'the subset is every 74th record of the 1851, which is a slice and not a selection'
