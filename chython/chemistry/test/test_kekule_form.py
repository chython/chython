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
"""`standardize_kekule`: which Kekule form is stored decides what `thiele()` is able to see."""
from random import Random
from chython.chemistry._canonicalize import canonicalize
from chython.chemistry._implicit import check_valence
from chython.chemistry._kekule_form import standardize_kekule
from chython.core import INFO, MoleculeContainer, read_smiles


#: dibenzo[a,e]cyclooctatetraene, C16H12, in two of its five Kekule forms.  In the first the benzo rings
#: lend their fusion carbons' double bonds to the eight-ring between them, which `thiele()` reads as an
#: exocyclic double bond and refuses; in the second both benzo rings hold three of their own.
COT_LENT = 'C1=CC2=CC=C3C=CC=CC3=CC=C2C=C1'
COT_HELD = 'C1=CC=C2C=CC3=CC=CC=C3C=CC2=C1'

#: biphenylene, C12H8: the same shape with a four-ring between the two benzo rings.
BIPHENYLENE_LENT = 'C1=CC2=C3C=CC=CC3=C2C=C1'
BIPHENYLENE_HELD = 'C1=CC2=C(C=C1)C1=CC=CC=C21'

PAIRS = [(COT_LENT, COT_HELD, 'dibenzo[a,e]cyclooctatetraene'),
         (BIPHENYLENE_LENT, BIPHENYLENE_HELD, 'biphenylene')]

#: Every small ring already holds every double bond it has room for, so there is nothing to choose:
#: where the neighbouring ring is a candidate too, the fusion bond's double stays inside the set.  The
#: first two are naphthalene with its fusion bond written single and written double -- two forms, both
#: left alone: the phase stage sees `rings`, which is the SSSR, so the ten-atom perimeter alternation of
#: the first is not a ring it can offer to shift.
FULL = ['C1=CC2=CC=CC=C2C=C1', 'C1=CC2=C(C=C1)C=CC=C2',
        'C1=CC2=CC3=CC=CC=C3C=C2C=C1', 'C1=CC=C2C(=C1)C=CC1=CC=CC=C12',
        'C1=CC2=CC=C3C=CC4=CC=CC5=C4C3=C2C5=C1', 'C1=CC2=CC=CC=CC2=C1', 'C=CC1=CC=CC=C1',
        'C1=CC=CC=C1', 'C1=CC2=CC3=CC=C(N3)C=C4C=CC(=N4)C=C5C=CC(=N5)C=C1N2']

#: A ring atom whose double bond has nowhere else to go.  The oxygens and the exocyclic methylenes are
#: what `thiele()`'s exocyclic rule is for, and no Kekule form of these puts those bonds in the ring.
PINNED = ['O=C1C=CC(=O)C=C1', 'O=C1C=CC=CN1', 'C1=CC=CC=CC=C1', 'C=C1C=CC=C1', 'C=C1C=CC=CC1=C']


def state(molecule):
    return {n: (molecule.implicit_h_of(n), molecule.charge_of(n), molecule.radical_of(n))
            for n in molecule.atom_numbers}


def aromatic(molecule):
    molecule.thiele()
    return sum(1 for b in molecule.bonds() if b.order == 4)


def shuffled(source, seed):
    """The same molecule with its atoms added in another order: same compound, other stable ids."""
    atoms = [(a.n, a.atomic_symbol, source.charge_of(a.n), source.implicit_h_of(a.n))
             for a in source.atoms()]
    bonds = [(b.n, b.m, b.order) for b in source.bonds()]
    order = [n for n, *_ in atoms]
    Random(seed).shuffle(order)
    by = {n: rest for n, *rest in atoms}
    out = MoleculeContainer()
    new = {}
    with out.edit() as e:
        for n in order:
            element, charge, hydrogens = by[n]
            new[n] = e.add_atom(element, charge=charge, implicit_h=hydrogens)
        for u, v, bond in bonds:
            e.add_bond(new[u], new[v], bond)
    return out


def test_a_benzo_ring_gets_back_the_double_bonds_it_lent_its_neighbour():
    """The defect this pass exists for: the ring is a candidate and the ring it is fused to is not."""
    for lent, held, label in PAIRS:
        m = read_smiles(lent)
        assert standardize_kekule(m), f'{label}: {lent} was left in the form that hides its benzo rings'
        assert aromatic(m) == 12, f'{label}: {lent} did not aromatise after the form was chosen'
        assert not standardize_kekule(read_smiles(held)), f'{label}: {held} already holds them'


def test_both_kekule_forms_of_one_compound_canonicalize_alike():
    for lent, held, label in PAIRS:
        a, b = read_smiles(lent), read_smiles(held)
        canonicalize(a)
        canonicalize(b)
        assert a.canonical_bytes == b.canonical_bytes, f'{label}: {a} and {b} still differ'


def test_only_the_orders_move():
    """Another compound is not another form, so the hydrogens, charges and radicals are all still there."""
    for lent, _, label in PAIRS:
        m = read_smiles(lent)
        formula, before = m.brutto_formula, state(m)
        standardize_kekule(m)
        assert m.brutto_formula == formula, f'{label}: {formula} became {m.brutto_formula}'
        assert state(m) == before, f'{label}: a hydrogen, charge or radical moved with the orders'
        assert not check_valence(m), f'{label}: {check_valence(m)}'


def test_a_ring_holding_every_double_it_has_room_for_is_left_alone():
    for smiles in FULL:
        m = read_smiles(smiles)
        assert not standardize_kekule(m), f'{smiles} was rewritten with nothing to gain'


def test_a_double_bond_that_cannot_reach_the_ring_stays_where_it_is():
    """p-benzoquinone is the reason `thiele()` refuses an exocyclic double bond, and it still does."""
    for smiles in PINNED:
        m = read_smiles(smiles)
        assert not standardize_kekule(m), f'{smiles} was rewritten'
        assert not aromatic(m), f'{smiles} came out aromatic'


def test_an_aromatic_ring_is_not_this_passs_to_choose():
    m = read_smiles('c1ccccc1')
    assert not standardize_kekule(m)
    assert sum(1 for b in m.bonds() if b.order == 4) == 6


def test_the_answer_does_not_depend_on_the_atom_order():
    """A round takes the best trial by score and then by canonical bytes, so it cannot."""
    for lent, held, label in PAIRS:
        keys = set()
        for smiles in (lent, held):
            source = read_smiles(smiles)
            for seed in range(10):
                m = shuffled(source, seed)
                standardize_kekule(m)
                m.thiele()
                keys.add(m.canonical_bytes)
        assert len(keys) == 1, f'{label}: {len(keys)} forms over 20 drawings of one compound'


def test_it_is_idempotent():
    for lent, _, label in PAIRS:
        m = read_smiles(lent)
        assert standardize_kekule(m), label
        assert not standardize_kekule(m), f'{label}: the pass moved twice'


def test_a_configuration_on_a_shiftable_bond_does_not_pin_the_form():
    """The eight-ring's double bonds are an alternating cycle's to move, so a sign on one is the
    drawing's phase and not the compound's geometry: it does not stop the benzo rings being filled.
    """
    m = read_smiles(r'C1=C/C2=C/C=C3/C=CC=C/C/3=C/C=C\2/C=C1')
    assert [u['parity'] for u in m.stereo_units()] == [1, 1, 2, 2]
    assert not any(u['stereogenic'] for u in m.stereo_units())
    assert standardize_kekule(m)
    assert aromatic(m) == 12


def test_both_configured_drawings_of_one_compound_canonicalize_alike():
    """The lent and held spellings of dibenzo[a,e]cyclooctatetraene, each with its parities written."""
    a = read_smiles(r'C1=C/C2=C/C=C3/C=CC=C/C/3=C/C=C\2/C=C1')
    b = read_smiles(r'C1=CC=C2\C=C/C3=CC=CC=C3\C=C/C2=C1')
    canonicalize(a)
    canonicalize(b)
    assert a.canonical_bytes == b.canonical_bytes, f'{a} and {b} still differ'


def test_a_ring_too_big_to_aromatise_gets_one_of_its_two_alternations():
    """1,2-dimethylcyclooctatetraene, whose two bond-shift drawings nothing downstream collapses."""
    plain = ('CC1=CC=CC=CC=C1C', 'CC1=C(C)C=CC=CC=C1')
    one, two = read_smiles(plain[0]), read_smiles(plain[1])
    standardize_kekule(one)
    standardize_kekule(two)
    assert one.canonical_bytes == two.canonical_bytes, f'{one} and {two} chose apart'

    # and with the phase written as a configuration, which `canonicalize()` stops reading as one
    for a, b in [plain, (r'C/C1=C/C=C\C=C/C=C\1/C', r'C/C1=C(\C)/C=C\C=C/C=C\1')]:
        one, two = read_smiles(a), read_smiles(b)
        canonicalize(one)
        canonicalize(two)
        assert one.canonical_bytes == two.canonical_bytes, f'{one} and {two} still differ'


def test_a_ring_that_is_not_one_alternation_keeps_its_configuration():
    """(Z,Z)-1,5-cyclooctadiene: eight atoms and two double bonds no alternation joins, so the phase
    stage does not reach it and both parities are the compound's.
    """
    m = read_smiles(r'C1=C\CC/C=C\CC/1')
    assert [u['parity'] for u in m.stereo_units() if u['stereogenic']] == [2, 2]
    assert not standardize_kekule(m)
    canonicalize(m)
    assert [u['parity'] for u in m.stereo_units() if u['stereogenic']] == [2, 2]


def test_the_record_names_the_ring_and_is_information_rather_than_a_repair():
    m = read_smiles(COT_LENT)
    standardize_kekule(m)
    records = [r for r in m.log if r.stage == 'kekule-form']
    assert len(records) == 1
    assert records[0].rule == 'kekule-form:ring-filled'
    assert records[0].severity == INFO
    assert len(records[0].atoms) == 6


def test_the_phase_record_names_the_ring_it_shifted():
    m = read_smiles('CC1=CC=CC=CC=C1C')
    standardize_kekule(m)
    records = [r for r in m.log if r.stage == 'kekule-form']
    assert len(records) == 1
    assert records[0].rule == 'kekule-form:ring-phase'
    assert records[0].severity == INFO
    assert len(records[0].atoms) == 8
