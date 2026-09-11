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
"""The SMILES writer's constitution: the traversal, the atom and bond tokens, and the ring closures.

Stereo output is tested in test_smiles_write_stereo.py, which is where the creation-order sweep that
ruling F26 asked for lives; the sweep HERE covers the constitution only, and the two are separate
because a constitution sweep that passes proves nothing about parities.
"""
from itertools import permutations
from math import factorial
from random import Random, seed

from pytest import mark, raises

from chython.core import MoleculeContainer
from chython.core._core import (normalize_smiles_spec, read_smiles, smw_symbol_table,
                                smw_traversal, smv_valence_model, write_smiles)


# ------------------------------------------------------------------------------------------------
# FIXTURE PLUMBING.  Every molecule below is described as (atoms, bonds) with atoms indexed from 0,
# so that `build` and `build_in_order` can produce the SAME molecule from different creation orders
# -- which is the whole point of the sweep at the bottom of the file.
def build_in_order(atoms, bonds, order):
    """The molecule with its atoms created in `order` (a permutation of range(len(atoms)))."""
    m = MoleculeContainer()
    sids = {}
    for j in order:
        element, hydrogens = atoms[j]
        sids[j] = m.add_atom(element, implicit_h=hydrogens)
    for a, b, o in bonds:
        m.add_bond(sids[a], sids[b], o)
    return m


def build(atoms, bonds):
    return build_in_order(atoms, bonds, range(len(atoms)))


def one_atom(element, **kwargs):
    m = MoleculeContainer()
    m.add_atom(element, **kwargs)
    return m


ETHANOL = ([(6, 3), (6, 2), (8, 1)], [(0, 1, 1), (1, 2, 1)])
BENZENE = ([(6, 1)] * 6,
           [(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 0, 1)])
NAPHTHALENE = ([(6, 1)] * 8 + [(6, 0), (6, 0)],
               [(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 8, 1), (8, 4, 2), (4, 5, 1), (5, 6, 2),
                (6, 7, 1), (7, 9, 2), (9, 0, 1), (8, 9, 1)])
CUBANE = ([(6, 1)] * 8,
          [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 0, 1), (4, 5, 1), (5, 6, 1), (6, 7, 1), (7, 4, 1),
           (0, 4, 1), (1, 5, 1), (2, 6, 1), (3, 7, 1)])
# spiro[3.3]heptane: atom 3 is the spiro centre, in both rings and carrying no hydrogen.
SPIRO = ([(6, 2)] * 3 + [(6, 0)] + [(6, 2)] * 3,
         [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 0, 1), (3, 4, 1), (4, 5, 1), (5, 6, 1), (6, 3, 1)])
PYRIDINE = ([(7, 0)] + [(6, 1)] * 5,
            [(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 0, 1)])
# Two hubs each bonded to twelve bridges: eleven independent cycles, and eleven ring closures all
# opened at the same atom, which is the only way to reach `%10`.
_K2_12_ATOMS = [(6, 0)] * 14
_K2_12_BONDS = [(0, i, 1) for i in range(2, 14)] + [(1, i, 1) for i in range(2, 14)]
K2_12 = (_K2_12_ATOMS, _K2_12_BONDS)


# ------------------------------------------------------------------------------------------------
# THE TABLES.
def test_symbol_table_matches_the_arena_element_table():
    """SMW_SYMBOL is a second copy of `_elements.pxi`'s SYMBOLS, so it is checked against it.

    Through `add_atom`, which is the only Python door onto SYMBOL_TO_NUMBER and therefore onto
    SYMBOLS: if the writer's table disagreed by one entry the atomic number would come back wrong.
    """
    table = smw_symbol_table()
    assert len(table) == 119
    assert table[0] == 'R'  # element 0 is the fragment marker, not an unused slot
    for number in range(1, 119):
        m = MoleculeContainer()
        sid = m.add_atom(table[number])
        assert m.element_of(sid) == number, table[number]


def test_valence_model_is_the_documented_pair_of_sets():
    """The two valence sets, as measured against RDKit 2026.03.4 and OpenSMILES §3.1.5.

    Restated here rather than derived, because the point of the table is that the numbers were
    measured: a test that recomputed them from the same code would agree with any typo.
    """
    assert smv_valence_model() == {
        5: ((3,), (3,)),                    # B
        6: ((4,), (4,)),                    # C
        7: ((3,), (3, 5)),                  # N -- the models disagree above valence 3
        8: ((2,), (2,)),                    # O
        9: ((1,), (1,)),                    # F
        15: ((3, 5), (3, 5)),               # P
        16: ((2, 4, 6), (2, 4, 6)),         # S
        17: ((1,), (1,)),                   # Cl
        35: ((1,), (1,)),                   # Br
        53: ((1,), (1, 3, 5, 7))}           # I -- and here too


def test_the_two_valence_models_answer_differently():
    """`smv_default_h` and `_valence.pxi`'s `val_*` are NOT one table, and merging them breaks output.

    Mirrored by a test of the same name in test_valence.py, and it exists as a TEST rather than a
    comment on purpose: the person who is about to merge two tables reads a comment saying "these are
    not duplicates", decides it is stale, and deletes the comment along with the duplication.

    The witness is a neutral sulfur at valence 6, because it is unanswerable rather than merely true.
    Hexamethylsulfur is not a compound and `val_*` has no rule for it in ANY environment, while
    dimethyl sulfone -- same element, same charge, same valence -- is in V2's tables and gets zero
    hydrogens. So the chemical model's answer moves with the environment. This model's cannot even
    ask: a SMILES atom's syntax may not depend on what it is bonded to, only on its bond-order sum,
    so both sulfurs below are a bare `S` with no hydrogens. Different arity, therefore two tables.

    A gate on output would have written `[S]` or refused here, and the fidelity invariant says a
    library must be able to show a user the broken structure they actually handed it.
    """
    hexamethylsulfur = ([(16, 0)] + [(6, 3)] * 6, [(0, i, 1) for i in range(1, 7)])
    dimethylsulfone = ([(16, 0), (8, 0), (8, 0), (6, 3), (6, 3)],
                       [(0, 1, 2), (0, 2, 2), (0, 3, 1), (0, 4, 1)])
    assert write_smiles(build(*hexamethylsulfur)) == 'CS(C)(C)(C)(C)C', \
        'a valence-6 sulfur is spelled bare because the notation permits what the chemistry does not'
    assert write_smiles(build(*dimethylsulfone)) == 'CS(=O)(C)=O', \
        'and the sulfone gets the same treatment from the same input, which val_* cannot do'

    # The nitrogen half of the mirror, and it is one-sided in the opposite direction: `val_*` has no
    # rule for a neutral 5-valent nitrogen in any environment, while this model answers zero
    # hydrogens and drops the brackets, because 5 is above every narrow valence AND is a wide one.
    # Nitro written without charges is exactly that atom, and it is all over real corpora.
    uncharged_nitro = ([(7, 0), (8, 0), (8, 0), (6, 3)], [(0, 1, 2), (0, 2, 2), (0, 3, 1)])
    assert write_smiles(build(*uncharged_nitro)) == 'CN(=O)=O', \
        'a 5-valent neutral N is spelled bare here and has no chemical rule at all in val_*'
    # Where THIS model's two halves disagree is one bond lower, and that is a different fact: at sum
    # 4 the narrow set saturates to zero hydrogens and the wide set infers one, so the count is not
    # safe to omit and the atom brackets. Deleting either half of `smv_default_h` moves this line.
    tetramethylammonium_shaped = ([(7, 0)] + [(6, 3)] * 4, [(0, i, 1) for i in range(1, 5)])
    assert write_smiles(build(*tetramethylammonium_shaped)) == 'C[N](C)(C)C', \
        'at bond-order sum 4 the narrow and wide models disagree about N, so the count is not omitted'


def test_elements_outside_the_organic_subset_are_always_bracketed():
    assert write_smiles(one_atom(11)) == '[Na]'          # Na
    assert write_smiles(one_atom(26)) == '[Fe]'          # Fe
    assert write_smiles(one_atom(1)) == '[H]'            # even hydrogen: `H` alone is not an atom


# ------------------------------------------------------------------------------------------------
# THE BRACKET PREDICATE, one disjunct at a time.  Each pair below is (a molecule with exactly one
# reason to bracket, the same molecule with that reason removed): the second half is what makes the
# test fail if the corresponding clause is deleted, because deleting a clause turns the first half
# into the second half's answer and nothing else moves.
def test_bracket_because_the_element_is_not_in_the_subset():
    assert write_smiles(one_atom(11, implicit_h=0)) == '[Na]'
    assert write_smiles(one_atom(6, implicit_h=4)) == 'C'


def test_bracket_because_of_an_isotope():
    assert write_smiles(one_atom(6, implicit_h=4, isotope=13)) == '[13CH4]'
    assert write_smiles(one_atom(6, implicit_h=4)) == 'C'


def test_bracket_because_of_a_charge():
    assert write_smiles(one_atom(6, implicit_h=4, charge=1)) == '[CH4+]'
    assert write_smiles(one_atom(6, implicit_h=4, charge=2)) == '[CH4+2]'
    assert write_smiles(one_atom(6, implicit_h=4, charge=-1)) == '[CH4-]'
    assert write_smiles(one_atom(6, implicit_h=4, charge=-3)) == '[CH4-3]'
    # `!z` drops the charge, and with it the only reason to bracket -- so the clause is conditional
    # on the option and the test says which way.
    assert write_smiles(one_atom(6, implicit_h=4, charge=1), '!z') == 'C'


def test_bracket_because_of_a_radical():
    assert write_smiles(one_atom(6, implicit_h=4, radical=True)) == '[CH4] |^1:0|'
    assert write_smiles(one_atom(6, implicit_h=4)) == 'C'


def test_bracket_because_of_a_map_number():
    assert write_smiles(one_atom(6, implicit_h=4, map_number=5), 'm') == '[CH4:5]'
    # The map number is stored either way; without `m` it is not written and not a reason.
    assert write_smiles(one_atom(6, implicit_h=4, map_number=5)) == 'C'


def test_bracket_because_the_caller_asked_for_explicit_hydrogens():
    assert write_smiles(one_atom(6, implicit_h=4), 'h') == '[CH4]'
    assert write_smiles(one_atom(6, implicit_h=4)) == 'C'


def test_bracket_because_the_hydrogen_count_would_not_survive_omission():
    """The carbene case: a bare `C` carrying two hydrogens reads back as methane, so the count
    brackets."""
    assert write_smiles(one_atom(6, implicit_h=2)) == '[CH2]'
    assert write_smiles(one_atom(6, implicit_h=1)) == '[CH]'
    assert write_smiles(one_atom(6, implicit_h=0)) == '[C]'
    assert write_smiles(one_atom(6, implicit_h=4)) == 'C'


def test_bracket_because_the_two_valence_models_disagree():
    """A four-bonded N and a two-bonded I: `NH` to OpenSMILES, `N` to RDKit, and the mirror image.

    Neither may be written bare, and the reason is not any of the other seven clauses -- the count
    the writer holds equals the count ONE of the two models infers.
    """
    n, o, c1, c2 = (7, 1), (8, 0), (6, 3), (6, 3)
    amine_oxide = ([n, o, c1, c2], [(0, 1, 2), (0, 2, 1), (0, 3, 1)])
    assert write_smiles(build(*amine_oxide)) == 'C[NH](=O)C'
    diiodo = ([(53, 1), (6, 3), (6, 3)], [(0, 1, 1), (0, 2, 1)])
    assert write_smiles(build(*diiodo)) == 'C[IH]C'
    # A ONE-bonded iodine bearing no hydrogen agrees in both models -- both infer zero -- and is
    # written bare.  It is the same element and the same clause; only the bond-order sum moved.
    assert write_smiles(build([(53, 0), (6, 3)], [(0, 1, 1)])) == 'CI'


# ------------------------------------------------------------------------------------------------
# THE TRAVERSAL.
def test_every_atom_is_emitted_exactly_once():
    for name, (atoms, bonds) in [('ethanol', ETHANOL), ('benzene', BENZENE),
                                 ('naphthalene', NAPHTHALENE), ('cubane', CUBANE),
                                 ('spiro', SPIRO), ('K2,12', K2_12)]:
        m = build(atoms, bonds)
        order = smw_traversal(m)['order']
        assert len(order) == len(atoms), name
        assert set(order) == set(m.atom_numbers), name


def test_every_bond_is_classified_exactly_once():
    """Tree edges plus ring closures partition the bonds, and the closure count is the cycle rank."""
    for name, (atoms, bonds) in [('ethanol', ETHANOL), ('benzene', BENZENE),
                                 ('naphthalene', NAPHTHALENE), ('cubane', CUBANE),
                                 ('spiro', SPIRO), ('K2,12', K2_12)]:
        m = build(atoms, bonds)
        report = smw_traversal(m)
        seen = [frozenset(e) for e in report['tree']]
        seen += [frozenset(e[:2]) for e in report['closures']]
        assert len(seen) == len(set(seen)) == len(bonds), name
        # One component in each fixture, so the cycle rank is bonds - atoms + 1.
        assert len(report['closures']) == len(bonds) - len(atoms) + 1, name


def test_components_are_ordered_by_their_minimum_canonical_position():
    """Two components, and the one holding the globally-smallest canonical position goes first.

    Water and methane, built water-first: the string still starts with whichever component the
    canonical order puts first, so the assertion is on the ORDER of the two halves and not on which
    of them it is -- that is `canonical_order`'s business, not the writer's.
    """
    both = build([(8, 2), (6, 4)], [])
    text = write_smiles(both)
    assert text in ('O.C', 'C.O')
    order = smw_traversal(both)['order']
    positions = both.canonical_order()
    assert positions[order[0]] < positions[order[1]]
    # Built in the other creation order, the string is the same one.
    assert write_smiles(build([(6, 4), (8, 2)], [])) == text


def test_the_dot_separates_components_and_nothing_else():
    text = write_smiles(build([(8, 2), (6, 4), (17, 1)], []))
    assert text.count('.') == 2
    assert sorted(text.split('.')) == ['C', 'Cl', 'O']


# ------------------------------------------------------------------------------------------------
# BOND TOKENS.
def test_bond_orders_have_the_expected_spellings():
    assert write_smiles(build([(6, 3), (6, 3)], [(0, 1, 1)])) == 'CC'
    assert write_smiles(build([(6, 2), (6, 2)], [(0, 1, 2)])) == 'C=C'
    assert write_smiles(build([(6, 1), (6, 1)], [(0, 1, 3)])) == 'C#C'
    # Order 8 is chython 2's dialect any-bond; the arena can hold it, so the writer spells it.
    # An order-8 bond takes the atom out of the valence model entirely -- there is no count to
    # infer through a bond of unknown order -- so both atoms bracket.
    assert write_smiles(build([(6, 3), (6, 3)], [(0, 1, 8)])) == '[CH3]~[CH3]'


def test_no_bond_tokens_at_all_under_not_b():
    assert write_smiles(build([(6, 2), (6, 2)], [(0, 1, 2)]), '!b') == 'CC'


# ------------------------------------------------------------------------------------------------
# RING CLOSURES.
def test_a_ring_closure_number_is_reused_after_it_is_released():
    """Two separate rings joined by a single bond: the second ring gets number 1 again."""
    atoms = [(6, 2)] * 2 + [(6, 1)] + [(6, 2)] * 2 + [(6, 1)]
    bonds = [(0, 1, 1), (1, 2, 1), (2, 0, 1), (2, 5, 1), (3, 4, 1), (4, 5, 1), (5, 3, 1)]
    text = write_smiles(build(atoms, bonds))
    assert text.count('1') == 4                      # opened and closed twice, number 1 both times
    assert '2' not in text


def test_a_closure_number_is_not_released_before_the_atom_finishes_its_list():
    """The `C11` hazard: an atom that closes 1 and opens another closure must not write `11`.

    Spiro[3.3]heptane in STORED order, whose creation order puts the spiro centre fourth, so the
    first ring closes at exactly the atom the second ring opens at.  Releasing the number as soon as
    it closes would spell `C11`, which every reader takes as closure eleven -- one ring instead of
    two, and no error anywhere.
    """
    text = write_smiles(build(*SPIRO), 'i')
    assert text == 'C1CCC12CCC2'
    assert 'C11' not in text
    report = smw_traversal(build(*SPIRO), 'i')
    assert sorted(c[2] for c in report['closures']) == [1, 2]


def test_ten_simultaneous_closures_reach_the_percent_form():
    """Eleven closures opened at one atom, so numbers 10 and 11 are needed and must wear `%`."""
    text = write_smiles(build(*K2_12))
    assert '%10' in text and '%11' in text
    assert '%12' not in text
    report = smw_traversal(build(*K2_12))
    assert sorted(c[2] for c in report['closures']) == list(range(1, 12))


def test_running_out_of_closure_numbers_raises():
    """A hundred simultaneous closures: the writer refuses rather than writing an ambiguous string."""
    atoms = [(6, 0)] * 103
    bonds = [(0, i, 1) for i in range(2, 103)] + [(1, i, 1) for i in range(2, 103)]
    with raises(ValueError, match='closure numbers are exhausted'):
        write_smiles(build(atoms, bonds))


def test_a_ring_closure_carries_its_bond_order():
    text = write_smiles(build(*BENZENE))
    assert text == 'C1=CC=CC=C1'


# ------------------------------------------------------------------------------------------------
# THE STRINGS.  Checked against a hand-read expectation rather than against another writer, because
# the point of the epic is that this writer is the definition.
def test_known_strings():
    assert write_smiles(build(*ETHANOL)) == 'C(C)O'
    assert write_smiles(build(*BENZENE)) == 'C1=CC=CC=C1'
    assert write_smiles(build(*PYRIDINE)) == 'C1=CC=NC=C1'
    assert write_smiles(build(*NAPHTHALENE)) == 'C1=CC=2C(C=C1)=CC=CC=2'
    assert write_smiles(build(*CUBANE)) == 'C12C3C4C5C3C1C5C24'
    assert write_smiles(build(*SPIRO)) == 'C1CC2(CCC2)C1'


def test_an_empty_molecule_writes_an_empty_string():
    assert write_smiles(MoleculeContainer()) == ''
    assert smw_traversal(MoleculeContainer()) == {'order': (), 'tree': (), 'closures': (),
                                                  'directions': {}, 'lost': (), 'tokens': {},
                                                  'unknown_h': (), 'attachments': ()}


# ------------------------------------------------------------------------------------------------
# THE FORMAT SPEC.
def test_the_spec_keys_and_their_negations():
    m = build([(6, 3), (6, 1), (8, 0)], [(0, 1, 1), (1, 2, 2)])       # acetaldehyde
    assert write_smiles(m, '') == 'C(C)=O'
    assert write_smiles(m, '!b') == 'C(C)O'
    assert write_smiles(m, 'h') == '[CH]([CH3])=[O]'
    assert write_smiles(m, 'i') == 'CC=O'


def test_the_written_order_comes_back_beside_the_string():
    """`return_order=True` answers `(string, order)`, and the order describes THAT string.

    chython 2's `__format__(spec, _return_order=True)` and `smiles_atoms_order`, which exist for the
    callers that cannot use the string alone -- a reaction's CXSMILES tail indexes radicals by their
    position across every molecule in it.  Asserted against `smw_traversal`, which computes the same
    order by the same route, so the two cannot drift; and asserted to CHANGE with the spec, because
    the order is a function of the whole option set and a caller who fetched it under one spec must
    not reuse it under another.
    """
    m = build(*ETHANOL)
    plain = write_smiles(m)
    written, order = write_smiles(m, return_order=True)
    assert written == plain
    assert order == smw_traversal(m)['order']
    assert sorted(order) == sorted(m.atoms_order)       # a permutation of the stable ids
    stored, stored_order = write_smiles(m, 'i', return_order=True)
    assert stored == write_smiles(m, 'i')
    assert stored_order == smw_traversal(m, 'i')['order']
    assert stored_order != order, (order, stored_order)


def test_the_written_order_of_an_empty_molecule_is_empty():
    assert write_smiles(MoleculeContainer(), return_order=True) == ('', ())


def test_two_specs_normalize_equal_exactly_when_they_write_the_same_string():
    """THE CACHE KEY CONTRACT, and it is measured rather than declared.

    `normalize_smiles_spec` exists so that a cache can key on a spec without reimplementing the spec
    grammar -- and its whole value is the biconditional in the name of this test. So: take every
    permutation of a five-key spec, and every prefix length, normalize each, and require that two specs
    share a normal form if and only if they write the same string. The `only if` half is the one that
    would catch a normalizer that dropped a key it should have kept; the `if` half catches one that
    kept a key that makes no difference.

    The keys are swept in the spelling that is NOT the default -- `!s`, `!b` and the three positive
    ones -- because `s` and `b` are on to begin with and a spec of no-ops would test nothing.
    """
    m = spec_fixture()
    keys = ('!s', 'A', 'm', 'h', '!b')
    by_norm = {}
    for perm in permutations(keys):
        for k in range(len(perm) + 1):
            spec = ''.join(perm[:k])
            norm = normalize_smiles_spec(spec)
            assert normalize_smiles_spec(norm) == norm, (spec, norm)   # idempotent, so it IS a key
            by_norm.setdefault(norm, set()).add(write_smiles(m, spec))
    assert len(by_norm) == 2 ** len(keys), sorted(by_norm)        # every subset, one normal form each
    for norm, strings in by_norm.items():
        assert len(strings) == 1, (norm, sorted(strings))         # same normal form => same string
    # ... and no key is dropped on the way in: each one alone changes this fixture's output, so a
    # normalizer that forgot to render one would collapse two of the 32 groups and fail above.
    plain = write_smiles(m)
    for key in keys:
        assert write_smiles(m, key) != plain, key


def spec_fixture():
    """One molecule that every format key visibly changes: a stereocentre, an aromatic ring, a charge,
    a radical, a map number, and hydrogens in three different situations (a CH, an NH3+, an aromatic
    CH and a bare halogen)."""
    m = MoleculeContainer()
    sids = [m.add_atom(6, implicit_h=0, map_number=7),                # 0: the stereocentre
            m.add_atom(9), m.add_atom(17),                            # 1, 2: F, Cl
            m.add_atom(7, charge=1, implicit_h=3),                    # 3: NH3+
            m.add_atom(6, implicit_h=0), m.add_atom(6, implicit_h=0),  # 4, 5: ring
            m.add_atom(6, implicit_h=1), m.add_atom(6, implicit_h=0),  # 6, 7: ring
            m.add_atom(6, implicit_h=1), m.add_atom(6, implicit_h=1),  # 8, 9: ring
            m.add_atom(8, implicit_h=0, radical=True),                # 10: the phenoxyl radical
            m.add_atom(6, implicit_h=1), m.add_atom(8)]               # 11, 12: an aldehyde, for `!b`
    for a, b, o in ((0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1), (4, 5, 4), (5, 6, 4), (6, 7, 4),
                    (7, 8, 4), (8, 9, 4), (9, 4, 4), (7, 10, 1), (5, 11, 1), (11, 12, 2)):
        m.add_bond(sids[a], sids[b], o)
    m.set_parity(sids[0], 2)
    return m


def test_the_normalizer_refuses_what_the_writer_refuses():
    """One grammar, one refusal. A caller may normalize first and know the write will not fail on it."""
    for bad, message in (('Q', 'unknown format key'), ('!', 'ends with a bare'),
                         ('ir', 'both name the atom order'), ('ri', 'both name the atom order')):
        with raises(ValueError, match=message):
            normalize_smiles_spec(bad)
        with raises(ValueError, match=message):
            write_smiles(build(*ETHANOL), bad)


def test_an_alias_is_written_as_the_tail_s_label_field():
    # ONE ENTRY PER ATOM in the written order, trailing empties included, and the field first in the
    # block -- Marvin 25.1.3's own spelling, which is what keeps the string readable elsewhere.
    m = read_smiles('CCC')
    sids = [a.n for a in m.atoms()]
    m.set_aliases({sids[0]: b'Me'})
    assert write_smiles(m) == 'C(C)C |$;;Me$|'
    # no alias, no field; and `!x` suppresses the block the way it suppresses every other field
    assert write_smiles(read_smiles('CCC')) == 'C(C)C'
    assert write_smiles(m, '!x') == 'C(C)C'


def test_an_alias_survives_a_round_trip_through_the_string():
    # the atoms are re-ordered canonically, so what a round trip preserves is the PAIRING of a label
    # to its atom, which is what comparing the two molecules asserts
    for text in ('[Pol]CC[R3] |$;;;Resin$|', 'CCC |$Me;;OMe$|', 'C[C@H](N)O |$;lbl;;$|'):
        mol = read_smiles(text)
        back = read_smiles(write_smiles(mol))
        assert sorted(mol.aliases.values()) == sorted(back.aliases.values()), text
        assert mol == back, text
    # a marker's INDEX travels in the body, `[R3]`, and not as the tail's `_R3`
    assert write_smiles(read_smiles('[R3]C')) == '[R3]C'


def test_a_label_the_field_cannot_hold_is_written_as_a_character_reference():
    m = read_smiles('CC')
    sids = [a.n for a in m.atoms()]
    m.set_aliases({sids[0]: 'a;b|c$d&e f'.encode(), sids[1]: 'αβ'.encode()})
    text = write_smiles(m)
    assert text == 'CC |$&#945;&#946;;a&#59;b&#124;c&#36;d&#38;e&#32;f$|'
    assert sorted(read_smiles(text).aliases.values()) == sorted(m.aliases.values())


def test_unknown_and_refused_format_keys_raise():
    m = build(*ETHANOL)
    with raises(ValueError, match='unknown format key'):
        write_smiles(m, 'Q')
    with raises(ValueError, match="ends with a bare"):
        write_smiles(m, '!')


# ------------------------------------------------------------------------------------------------
# `r` -- A RANDOM ATOM ORDER.  Every fixture here is PARSED rather than built, because the assertion is
# a round trip and `spec_fixture` states no hydrogen count for its F, Cl and aldehyde O: those three are
# H_UNKNOWN, the writer spells them bare, and a reader answers 0 -- so it is unequal to its own string
# under any spec.
RANDOM_FIXTURE = 'C[C@H](N)/C=C/c1ccc([O])cc1[NH3+] |^1:8|'


def test_a_random_order_writes_the_same_molecule_many_ways():
    """`r` replaces the atom positions and nothing else, so every string it writes reads back equal.

    Both halves are the test: many DISTINCT strings (a stub that quietly kept the canonical order would
    give one) and every one of them the SAME molecule (an order that broke the traversal or the stereo
    signs would give a string that reads back as something else, or does not read at all).
    """
    m = read_smiles(RANDOM_FIXTURE)
    seen = {write_smiles(m, 'r') for _ in range(50)}
    assert len(seen) > 5, sorted(seen)
    for text in seen:
        assert read_smiles(text) == m, text


def test_a_random_order_carries_stereo_of_every_kind():
    """Tetrahedral, cis/trans and axial, over enough draws that each centre is written from several
    directions.  A sign that depended on the canonical order rather than on the written frame would
    survive the default spec and fail here."""
    for text in ('C[C@H](N)C(=O)O', 'C/C=C/C', 'CC=[C@]=CC', 'F[C@@H](Cl)[C@H](F)Br'):
        m = read_smiles(text)
        for _ in range(30):
            assert read_smiles(write_smiles(m, 'r')) == m, text


def test_a_random_order_is_seeded_from_the_random_module():
    """`random.seed()` reproduces a batch -- the property that makes `r` usable for augmentation, where
    a run has to be repeatable even though its strings are not predictable."""
    m = read_smiles(RANDOM_FIXTURE)
    seed(4)
    first = [write_smiles(m, 'r') for _ in range(8)]
    seed(4)
    assert [write_smiles(m, 'r') for _ in range(8)] == first
    assert len(set(first)) > 1                      # not one string repeated eight times


def test_the_two_atom_order_keys_cannot_be_combined():
    """`i` and `r` are two answers to the one question of where the order comes from, so a spec naming
    both raises either way round -- which is what keeps `normalize_smiles_spec` order-independent."""
    m = build(*ETHANOL)
    for spec in ('ir', 'ri'):
        with raises(ValueError, match='both name the atom order'):
            write_smiles(m, spec)
    assert normalize_smiles_spec('r') == 'r'
    assert normalize_smiles_spec('mr') == normalize_smiles_spec('rm') == 'rm'
    assert normalize_smiles_spec('i!ir') == 'r'     # resolved options, not the keys as written


# ------------------------------------------------------------------------------------------------
# CREATION-ORDER INVARIANCE (plan task 4).  The fixture the whole design exists for, in its
# constitution-only form.
#
# RULING F102 -- CAN THIS FIXTURE FAIL?  Yes, and it was made to.  Run against STORED-slot order
# (`write_smiles(m, 'i')`) instead of the canonical default, the same sweep produces
#
#   ethanol       4 distinct strings over the 6 creation orders
#   benzene       2 distinct strings over 720
#   pyridine     12
#   spiro         7
#   naphthalene  30
#   cubane       11
#
# against 1 apiece on the canonical path.  (Benzene's 2 is small because every atom is equivalent, so
# the only thing the creation order can move is which bond of the alternating pair the traversal
# enters -- the fixture is deliberately kept in the list as the weakest case that still fails.)
#
# so the assertion below is not vacuously true: it is the canonical path doing work, and
# test_stored_order_is_not_creation_order_invariant records that failure as a passing test so the
# evidence cannot rot.
SWEEP_LIMIT = 720


def _creation_orders(n):
    """Every permutation when there are few, a SEEDED SAMPLE when there are many.

    Sampled rather than truncated: `permutations` is lexicographic, so its first 720 entries of a
    10-atom molecule all share the same seven-atom prefix and vary only the tail -- a sample that
    tests almost nothing.  A fixed seed keeps the test reproducible.
    """
    if factorial(n) <= SWEEP_LIMIT:
        return list(permutations(range(n)))
    rng = Random(20260902)
    out = [tuple(range(n))]
    while len(out) < SWEEP_LIMIT:
        order = list(range(n))
        rng.shuffle(order)
        out.append(tuple(order))
    return out


def _distinct_over_creation_orders(atoms, bonds, spec=''):
    seen = set()
    orders = _creation_orders(len(atoms))
    for order in orders:
        seen.add(write_smiles(build_in_order(atoms, bonds, order), spec))
    return seen, len(orders)


def test_canonical_output_does_not_depend_on_the_creation_order():
    for name, (atoms, bonds) in [('ethanol', ETHANOL), ('benzene', BENZENE),
                                 ('pyridine', PYRIDINE), ('spiro', SPIRO),
                                 ('naphthalene', NAPHTHALENE), ('cubane', CUBANE)]:
        seen, count = _distinct_over_creation_orders(atoms, bonds)
        assert len(seen) == 1, (name, count, sorted(seen)[:4])
        assert count == min(factorial(len(atoms)), SWEEP_LIMIT), name


def test_stored_order_is_not_creation_order_invariant():
    """The could-have-failed evidence for the test above (ruling F102).

    Stored order is a function of the creation order BY DEFINITION, so the same sweep must produce
    more than one string.  If this test ever passes with one string, the sweep above has stopped
    measuring anything and both tests are wrong.
    """
    for name, (atoms, bonds) in [('ethanol', ETHANOL), ('benzene', BENZENE), ('spiro', SPIRO)]:
        seen, _ = _distinct_over_creation_orders(atoms, bonds, 'i')
        assert len(seen) > 1, name


def test_the_traversal_itself_is_creation_order_invariant_up_to_relabelling():
    """The shape of the traversal, not just the string: same tree, same closures, same numbers.

    Compared through canonical positions rather than stable ids -- the ids ARE the creation order,
    so comparing them would be comparing the input to itself.
    """
    atoms, bonds = NAPHTHALENE
    shapes = set()
    for order in _creation_orders(len(atoms)):
        m = build_in_order(atoms, bonds, order)
        pos = m.canonical_order()
        report = smw_traversal(m)
        # SORTED, because the probe lists edges in slot order and the slots ARE the creation
        # order -- an unsorted comparison would fail on the listing and not on the traversal.
        shapes.add((tuple(pos[s] for s in report['order']),
                    tuple(sorted((pos[a], pos[b]) for a, b in report['tree'])),
                    tuple(sorted((pos[a], pos[b], c) for a, b, c in report['closures']))))
    assert len(shapes) == 1


def test_the_arena_now_stores_an_aromatic_bond_so_the_hooks_are_live():
    """The premise the writer's three aromatic branches rest on, stated where they are tested.

    `add_bond(a, b, 4)` is accepted, so order 4 arrives from the ARENA rather than from perception, and
    it arrives on molecules the writer is asked to spell.  A build that refused the order instead would
    leave `smw_bond`, `smw_sticky_bond` and `smw_atom`'s aromatic branch inert and untestable.
    """
    m = MoleculeContainer()
    a = m.add_atom('C')
    b = m.add_atom('C')
    m.add_bond(a, b, 4)                       # accepted, not refused
    assert m.order_of(a, b) == 4 and m.aromatic_bond_count == 1
    write_smiles(m)                           # and the writer answers rather than crashing


def test_a_stored_aromatic_ring_is_written_aromatic():
    """`smw_bond`, `smw_sticky_bond` and `smw_atom` all read the stored order, so nothing is inferred.

    A writer deciding lowercase and `:` from an option instead spells aromatic benzene
    `[CH]1[CH][CH][CH][CH][CH]1` -- cyclohexane.  The rest of the aromatic surface is in
    test_smiles_write_aromatic.py.

    THE THREE HYDROGEN CASES ARE THE OTHER HALF, and no two of them may be confused for each other:
    one stated hydrogen is `c1ccccc1`, benzene; a STATED zero is `[c]1[c][c][c][c][c]1`, six bracketed
    aromatic carbons, because brackets are how SMILES states a count; and `add_atom('C')`, which states
    nothing and stores `H_UNKNOWN`, is the UNBRACKETED `c` -- SMILES for "the reader works it out".  A
    record that did not state a count must not be written as one that stated zero.
    """
    m = MoleculeContainer()
    ids = [m.add_atom('C', implicit_h=1) for _ in range(6)]
    for i in range(6):
        m.add_bond(ids[i], ids[(i + 1) % 6], 4)
    assert write_smiles(m) == 'c1ccccc1'

    zero = MoleculeContainer()
    ids = [zero.add_atom('C', implicit_h=0) for _ in range(6)]
    for i in range(6):
        zero.add_bond(ids[i], ids[(i + 1) % 6], 4)
    assert write_smiles(zero) == '[c]1[c][c][c][c][c]1', 'a STATED zero is bracketed'

    bare = MoleculeContainer()
    ids = [bare.add_atom('C') for _ in range(6)]
    for i in range(6):
        bare.add_bond(ids[i], ids[(i + 1) % 6], 4)
    assert write_smiles(bare) == 'c1ccccc1', 'an UNSTATED count is left to the reader, not called zero'


# ------------------------------------------------------------------------------------------------
# THE PUBLIC NAMES.  `chython.core`, not `chython.core._core`, is the import a caller writes.
def test_the_writers_entry_points_are_reachable_without_the_underscore_module():
    """Everything the rest of this file imports from `._core` is exported from the package.

    The tests reach into `chython.core._core` because that is where the symbols are DEFINED, and doing
    that everywhere hides the ordinary packaging mistake: a function that works perfectly and that no
    caller outside this repository can import.  `__all__` is checked as well as the attributes, because
    a name present but missing from `__all__` is invisible to `from chython.core import *` and to every
    documentation tool.

    `smw_symbol_table`, `smw_traversal` and `smw_stereo_seed_labels` are deliberately NOT here: they
    are probes this suite uses to look inside the writer, they have no caller-facing meaning, and
    exporting them would make three internal shapes part of the public surface.
    """
    import chython.core as core

    for name in ('write_smiles', 'normalize_smiles_spec', 'detached_smiles', 'DetachedSmiles'):
        assert hasattr(core, name), name
        assert name in core.__all__, name
    for probe in ('smw_symbol_table', 'smw_traversal', 'smw_stereo_seed_labels'):
        assert probe not in core.__all__, probe
    # and every name the package claims really resolves -- an `__all__` entry that does not is an
    # ImportError for anyone using the star form and nothing at all for anyone who is not
    for name in core.__all__:
        assert hasattr(core, name), name
    assert core.write_smiles is write_smiles
