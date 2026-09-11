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
"""The element and isotope tables: that the compiled arrays are the two TSVs, and that the TSVs
hold what they claim.

    elements.tsv + isotopes.tsv  ->  the arrays in _elements.pxi
    `-- hand-maintained data,        `-- generated; the sync test is what makes the TSVs
        and the authority                the authority rather than a copy of these numbers

The tables are checked-in data with no upstream, so "correct" is not a property anything can
re-derive -- what these tests hold is that the two files agree with each other and that the arrays
are their transposition.  Two tests at the foot use chython 2 as an independent witness that
nothing has been lost from its own tables -- through a subprocess, so this file does not import it;
nothing above depends on them and they skip when it is not provisioned.
"""
from math import fsum, isclose
from pytest import raises

from chython.core._core import (atomic_radius_table, element_period, element_symbols,
                                isotope_counts_table, isotope_data, isotope_offsets_table,
                                mdl_isotope_table, valence_electrons_table)
from .gen_element_tables import (ELEMENTS_HEADER, Element, Isotope, PXI, assert_invariants,
                                 compile_tables, read_elements, read_isotopes)


ELEMENTS = read_elements()
ISOTOPES = read_isotopes()


def test_the_compiled_tables_are_the_two_tsvs():
    """The gate the hand-written arrays never had.

    Without it, editing a TSV and forgetting to run `compile` changes nothing and says nothing --
    the .pxi keeps the old numbers and the file that claims to be the authority is not one.
    """
    assert compile_tables(ELEMENTS, ISOTOPES) in PXI.read_text(encoding='utf-8'), \
        'chython/core/_elements.pxi has drifted from the TSVs; run gen_element_tables.py compile'


def test_the_shipped_tables_satisfy_their_own_invariants():
    assert_invariants(ELEMENTS, ISOTOPES)


def test_the_symbol_column_of_isotopes_tsv_is_not_a_second_symbol_table():
    """It is there so a human can grep for an element, and `assert_invariants` checks it against
    elements.tsv rather than trusting it -- so a disagreement is a refusal to compile, not a second
    opinion about what element 26 is called."""
    broken = [Isotope(i.z, i.symbol, i.mass_number, i.mass, i.abundance) for i in ISOTOPES]
    broken[0].symbol = 'Xx'
    with raises(ValueError, match='one symbol table, not two'):
        assert_invariants(ELEMENTS, broken)


def test_offsets_are_the_prefix_sum_of_counts():
    """The relationship the two hand-written arrays left entirely unstated."""
    offsets = isotope_offsets_table()
    counts = isotope_counts_table()
    assert offsets[0] == 0 and counts[0] == 0        # index 0 unused: the index IS the number
    for z in range(1, 119):
        assert offsets[z] == sum(counts[1:z]), z
    assert offsets[118] + counts[118] == len(ISOTOPES)


def test_counts_are_the_run_lengths_of_the_flat_arrays():
    counts = isotope_counts_table()
    for z in range(1, 119):
        assert len(isotope_data(z)) == counts[z], z


def test_rows_are_sorted_by_mass_number_within_an_element():
    for z in range(1, 119):
        numbers = [a for a, _, _ in isotope_data(z)]
        assert numbers == sorted(numbers) == sorted(set(numbers)), z


def test_a_reordered_file_is_refused_because_the_sort_is_the_layout():
    with raises(ValueError, match='must be sorted'):
        assert_invariants(ELEMENTS, list(reversed(ISOTOPES)))


def test_weights_sum_to_one_or_to_zero_per_element():
    """There is no third case.  A partial set of weights would make `float(molecule)` answer with a
    number that is neither an average nor obviously wrong."""
    for z in range(1, 119):
        total = fsum(w for _, _, w in isotope_data(z))    # `fsum`, for the reason `assert_invariants` gives
        assert total == 0.0 or isclose(total, 1.0, abs_tol=1e-6), (z, total)


def test_a_partial_set_of_weights_is_refused():
    broken = [Isotope(i.z, i.symbol, i.mass_number, i.mass, i.abundance) for i in ISOTOPES]
    for i in broken:
        if i.z == 6 and i.mass_number == 12:
            i.abundance = 0.5                    # carbon-12 at half weight, carbon-13 unchanged
    with raises(ValueError, match='abundances sum to'):
        assert_invariants(ELEMENTS, broken)


# --- the defect the flat arrays hid ------------------------------------------------------------

def test_every_mdl_reference_mass_number_can_be_weighed():
    """`MDL_ISOTOPE[z]` is what an MDL file's mass-difference field is measured from, so a file is
    entitled to state exactly that mass number.  For 19 elements the isotope table had no row for
    it, and the atom then weighed nothing."""
    mdl = mdl_isotope_table()
    symbols = element_symbols()
    for z in range(1, 119):
        numbers = [a for a, _, _ in isotope_data(z)]
        assert mdl[z] in numbers, f'{symbols[z]}: MDL reference {mdl[z]} not in {numbers}'


def test_bromine_eighty_weighs_something():
    """The concrete case.  80 is what `MDL_ISOTOPE[35]` hands out -- 79.904 rounded -- and it used
    to be absent from bromine's rows, so `[80Br]` massed 0.0 while `[79Br]` massed 78.9."""
    masses = {a: m for a, m, _ in isotope_data(35)}
    assert 80 in masses
    assert 79 < masses[80] < 81


def test_a_missing_mdl_reference_row_is_refused():
    without_bromine_eighty = [i for i in ISOTOPES if not (i.z == 35 and i.mass_number == 80)]
    with raises(ValueError, match='no row for it'):
        assert_invariants(ELEMENTS, without_bromine_eighty)


def test_the_two_rows_nobody_has_a_mass_for_are_still_two():
    """Dubnium-270 and tennessine-297: named by MDL, weighed by nobody, so they compile to 0.0.
    Pinned by count so that a third one arriving is a failure and not a silent zero."""
    unweighable = [(i.symbol, i.mass_number) for i in ISOTOPES if i.mass is None]
    assert unweighable == [('Db', 270), ('Ts', 297)]
    assert [m for _, m, _ in isotope_data(105) if m == 0.0] == [0.0]


def test_a_massless_row_nothing_asks_for_is_refused():
    """A row with no mass and no weight does nothing but exist, and the one caller that needs a row
    to merely exist is the MDL reference lookup.  Any other massless row is a typo."""
    broken = [Isotope(i.z, i.symbol, i.mass_number, i.mass, i.abundance) for i in ISOTOPES]
    for i in broken:
        if i.z == 6 and i.mass_number == 14:     # carbon-14 is nobody's MDL reference
            i.mass = None
    with raises(ValueError, match='nothing needs the row'):
        assert_invariants(ELEMENTS, broken)


# --- the elements file -------------------------------------------------------------------------

def test_elements_tsv_covers_every_atomic_number_once():
    assert [e.z for e in ELEMENTS] == list(range(1, 119))
    with raises(ValueError, match='1..118, once each'):
        assert_invariants(ELEMENTS[:-1], ISOTOPES)


def test_a_repeated_symbol_is_refused():
    broken = [Element(e.z, e.symbol, e.mdl_isotope, e.valence_electrons, e.atomic_radius)
              for e in ELEMENTS]
    broken[1].symbol = 'H'
    with raises(ValueError, match='repeated symbol'):
        assert_invariants(broken, ISOTOPES)


def test_the_compiled_symbol_tuple_is_the_elements_file():
    symbols = element_symbols()
    assert symbols[0] == 'R'
    assert list(symbols[1:]) == [e.symbol for e in ELEMENTS]


def test_the_compiled_mdl_table_is_the_elements_file():
    mdl = mdl_isotope_table()
    assert mdl[0] == 0
    assert list(mdl[1:]) == [e.mdl_isotope for e in ELEMENTS]


# --- chython 2 as an independent witness -------------------------------------------------------
#
# Not an oracle: the TSVs are the authority, and these two tests neither define nor re-derive them.
# What they catch is a LOSS -- a nuclide or a mass that chython 2 carries and an edit here dropped
# or altered.  Nothing above depends on them.
#
# chython 2 is reached through `oracle`, which runs an INSTALLED copy in another interpreter, so
# this file imports no chython 2 and the witness survives V2 leaving the tree.  It is a genuine
# second opinion and not a copy of these numbers: V2's tables were compiled from the nuclear data
# independently of whoever wrote `isotopes.tsv`, which is exactly why a frozen snapshot of them
# would have been worth less -- it would be one more rendering of the same rows, checked against
# itself.  Unprovisioned, both tests skip and the self-consistency gates above still run.

V2_TABLES = """
from chython.periodictable.base import Element

masses = {}
mdl = {}
symbols = {}
radii = {}
for z in range(1, 119):
    e = Element.from_atomic_number(z)()
    symbols[str(z)] = e.atomic_symbol
    mdl[str(z)] = e.mdl_isotope
    radii[str(z)] = e.atomic_radius
    for a, mass in e.isotopes_masses.items():
        masses[f'{z}-{a}'] = (mass, e.isotopes_distribution.get(a, 0.0))
_emit({'masses': masses, 'mdl': mdl, 'symbols': symbols, 'radii': radii})
"""


def v2_tables():
    from .oracle import ask

    return ask(V2_TABLES)


# (z, mass_number) -> the weight chython 3 states instead of chython 2's, and why.  A MASS is never
# listed: a different mass is a different measurement and stays a failure.  A weight can differ for a
# reason that is not a measurement at all, and then the difference is declared here or it is a loss.
RESTATED_WEIGHTS = {
    (14, 28): (0.922296, 'silicon\'s three published abundances sum to 1.000001, which `compile` '
                         'refuses; the excess comes off the dominant nuclide -- see isotopes.tsv'),
}


def test_nothing_chython_two_carries_has_been_lost():
    """Exact float equality: both sides are IEEE doubles and JSON round-trips them exactly.

    A subset check, not equality, because the tables hold nuclides chython 2's do not -- every mass
    number `MDL_ISOTOPE` names has a row here, and for 19 elements chython 2 stopped short of its
    own reference value.  Extra rows are the point; changed and missing ones are the failure, unless
    the change is declared in `RESTATED_WEIGHTS` with its reason.
    """
    theirs = v2_tables()['masses']
    ours = {(i.z, i.mass_number): (i.mass, i.abundance) for i in ISOTOPES}

    for nuclide, (mass, abundance) in theirs.items():
        z, a = (int(p) for p in nuclide.split('-'))
        restated = RESTATED_WEIGHTS.get((z, a))
        expected = (mass, restated[0] if restated else abundance)
        assert ours.get((z, a)) == expected, nuclide
    # the comparison is evidence only if it looked at a real table: V2 carries about 350 nuclides,
    # and a bridge that handed back an empty dict would make this pass without comparing anything
    assert len(theirs) > 300, len(theirs)

    # the ratchet, in the other direction: an entry chython 2 turns out to agree with is deleted, not
    # left behind to describe a divergence that is not there
    agreed = sorted(k for k, (weight, _) in RESTATED_WEIGHTS.items()
                    if theirs.get('%d-%d' % k, (None, None))[1] == weight)
    assert not agreed, ('RESTATED_WEIGHTS declares a weight chython 2 already states: %s -- delete '
                        'the entry, the tables agree' % agreed)


def test_the_mdl_reference_table_agrees_with_chython_two():
    answer = v2_tables()
    for e in ELEMENTS:
        assert e.symbol == answer['symbols'][str(e.z)], e.z
        assert e.mdl_isotope == answer['mdl'][str(e.z)], e.symbol
    assert len(answer['symbols']) == 118


def test_the_compiled_isotope_rows_are_the_isotopes_file():
    flat = []
    for z in range(1, 119):
        for a, mass, weight in isotope_data(z):
            flat.append((z, a, mass, weight))
    assert flat == [(i.z, i.mass_number, 0.0 if i.mass is None else i.mass, i.abundance)
                    for i in ISOTOPES]


# --- the valence electron column ---------------------------------------------------------------


def test_the_valence_electron_column_is_the_group_number():
    """The convention D8 settles, spot-checked across all four blocks.

    Group number for groups 1-12 and group - 10 for groups 13-18, which is Kier and Hall's Zv for
    every main-group element: carbon 4, sulfur 6, chlorine 7.  Zinc is 12 and not 2 -- the filled
    d shell counts, which is the same convention the 18-electron rule uses.
    """
    by_symbol = {e.symbol: e.valence_electrons for e in ELEMENTS}
    assert by_symbol['H'] == 1
    assert by_symbol['He'] == 2
    assert by_symbol['C'] == 4
    assert by_symbol['N'] == 5
    assert by_symbol['O'] == 6
    assert by_symbol['F'] == 7
    assert by_symbol['Ne'] == 8
    assert by_symbol['Na'] == 1
    assert by_symbol['S'] == 6
    assert by_symbol['Cl'] == 7
    assert by_symbol['Sc'] == 3
    assert by_symbol['Fe'] == 8
    assert by_symbol['Zn'] == 12
    assert by_symbol['Ga'] == 3
    assert by_symbol['Br'] == 7
    assert by_symbol['Hf'] == 4
    assert by_symbol['Hg'] == 12
    assert by_symbol['Pb'] == 4
    assert by_symbol['Og'] == 8


def test_the_f_block_states_no_count_and_lanthanum_does():
    """`?`, the token isotopes.tsv already uses for a mass nobody has -- not a silent zero.

    The 4f and 5f electrons are neither reliably core nor reliably valence and every convention
    disagrees, so chython declines to pick one and `valence_electrons_count` refuses on such an
    atom.  Lanthanum and actinium are group 3 in every layout and are not part of the doubt.
    """
    by_symbol = {e.symbol: e.valence_electrons for e in ELEMENTS}
    assert by_symbol['La'] == 3
    assert by_symbol['Ac'] == 3
    assert by_symbol['Ce'] is None
    assert by_symbol['Gd'] is None
    assert by_symbol['Lu'] is None
    assert by_symbol['Th'] is None
    assert by_symbol['U'] is None
    assert by_symbol['Lr'] is None
    assert sum(1 for e in ELEMENTS if e.valence_electrons is None) == 28


def test_the_compiled_valence_electron_array_is_the_column():
    """The same gate the other five arrays have: the TSV is the authority or it is a copy."""
    table = valence_electrons_table()
    assert len(table) == 119
    assert table[0] == 0, 'index 0 is unused so that the index IS the atomic number'
    for e in ELEMENTS:
        expected = 0 if e.valence_electrons is None else e.valence_electrons
        assert table[e.z] == expected, f'{e.symbol}: {table[e.z]} != {expected}'


def test_an_unknown_count_compiles_to_the_reserved_zero():
    """0 is reserved: no element has zero valence electrons, so the sentinel costs no storage.

    It is a STORAGE spelling, exactly as H_UNKNOWN is, and the public surface never returns it --
    `valence_electrons_count` raises instead.
    """
    table = valence_electrons_table()
    assert table[58] == 0                      # cerium
    assert table[92] == 0                      # uranium
    assert all(table[z] for z in range(1, 58))
    assert all(table[z] for z in range(72, 90))
    assert all(table[z] for z in range(104, 119))


def test_the_period_boundaries_are_where_the_shells_close():
    """The period is arithmetic and NOT a column, so what is pinned is the six boundaries.

    2/3, 10/11, 18/19, 36/37, 54/55 and 86/87 are He/Li, Ne/Na, Ar/K, Kr/Rb, Xe/Cs and Rn/Fr -- the six
    places a period ends.  Asserting only that H is 1 and Og is 7 would pass on `1 + z // 17`, which is
    why the boundaries are the assertion and the ends are the sanity check.
    """
    for z in (2, 10, 18, 36, 54, 86):
        assert element_period(z) + 1 == element_period(z + 1), z
    assert element_period(1) == 1                               # hydrogen
    assert element_period(118) == 7                             # oganesson
    assert [element_period(z) for z in (6, 7, 8, 16, 17)] == [2, 2, 2, 3, 3]
    # non-decreasing, so no interior boundary can appear where no shell closes
    assert all(element_period(z) <= element_period(z + 1) for z in range(1, 118))


# --- the atomic radius column ------------------------------------------------------------------


def test_the_radius_column_is_the_calculated_radius_in_angstroms():
    """ONE RADIUS AND IT IS THE CALCULATED ONE: an SCF orbital measure, in angstroms.

    Not covalent and not van der Waals, and not the Hall-Kier alpha table in `_descriptors.pxi`, which
    is `r_cov / 0.77 - 1` per hybridization.  Spot-checked across the blocks; period 2 is the row that
    makes the unit unambiguous, every value there being under 2 A.
    """
    by_symbol = {e.symbol: e.atomic_radius for e in ELEMENTS}
    assert by_symbol['H'] == 0.53
    assert by_symbol['He'] == 0.31
    assert by_symbol['Li'] == 1.67
    assert by_symbol['C'] == 0.67
    assert by_symbol['N'] == 0.56
    assert by_symbol['O'] == 0.48
    assert by_symbol['F'] == 0.42
    assert by_symbol['Cl'] == 0.79
    assert by_symbol['Fe'] == 1.56
    assert by_symbol['I'] == 1.15
    assert by_symbol['Cs'] == 2.98


def test_every_element_states_a_radius():
    """No `?` in this column, unlike `valence_electrons`.

    The doubt the f block has about its valence electrons has no counterpart here -- a radius is a
    measure rather than a convention -- and the rows the published set does not reach state the group
    analogue's value instead, per the test below.
    """
    assert len([e for e in ELEMENTS if e.atomic_radius is not None]) == 118
    assert all(0.3 <= e.atomic_radius <= 3.0 for e in ELEMENTS)


def test_beyond_the_published_set_the_column_carries_the_group_analogue():
    """The calculated set ends at radon, and the last 32 rows take the value one period up.

    Fr takes Cs and Ra takes Ba; Rf through Og take Hf through Rn, 32 atomic numbers below each; and
    Ac through Lr take Lu's, the f block's own last row.  Pinned as a relationship so that a published
    actinide radius arriving is an edit to this test rather than a silent difference in a column of
    2.17s.
    """
    by_z = {e.z: e.atomic_radius for e in ELEMENTS}
    assert by_z[87] == by_z[55]                              # Fr takes Cs
    assert by_z[88] == by_z[56]                              # Ra takes Ba
    for z in range(104, 119):                                # Rf..Og take Hf..Rn
        assert by_z[z] == by_z[z - 32], z
    assert {by_z[z] for z in range(89, 104)} == {by_z[71]}    # Ac..Lr take Lu's


def test_the_column_may_not_say_it_does_not_know(tmp_path):
    """`?` is the token isotopes.tsv uses for a mass nobody has, and this column has no such row.

    Refused rather than compiled to a reserved zero: 0.0 is not a radius, and the two consumers of the
    column are geometric -- a sphere with no radius and a bond-perception threshold of zero are both
    wrong answers where a refusal names the row.
    """
    path = tmp_path / 'elements.tsv'
    path.write_text('\t'.join(ELEMENTS_HEADER) + '\n6\tC\t12\t4\t?\n')
    with raises(ValueError, match='states no atomic radius'):
        read_elements(path)


def test_a_radius_outside_the_plausible_range_is_refused():
    """The invariant that catches a misplaced decimal point, which is the one error this column has.

    167 for lithium is picometres in an angstrom column -- wider than any molecule -- and 1.67 is the
    value the rest of period 2 sits beside.
    """
    broken = [Element(e.z, e.symbol, e.mdl_isotope, e.valence_electrons, e.atomic_radius)
              for e in ELEMENTS]
    broken[2].atomic_radius = 167.0
    with raises(ValueError, match='atomic radius'):
        assert_invariants(broken, ISOTOPES)


def test_the_compiled_radius_array_is_the_column():
    """The same gate the other arrays have: the TSV is the authority or it is a copy."""
    table = atomic_radius_table()
    assert len(table) == 119
    assert table[0] == 0.0, 'index 0 is unused so that the index IS the atomic number'
    for e in ELEMENTS:
        assert table[e.z] == e.atomic_radius, e.symbol


def test_the_radius_column_agrees_with_chython_two_except_for_lithium():
    """The witness the mass and MDL columns have, over all 118 rows.

    ONE ROW DISAGREES DELIBERATELY.  chython 2.24 states 167 for lithium; V3 states 1.67, which is the
    value between helium's 0.31 and beryllium's 1.12 that the column's unit admits.  Named here because
    a witness that excused an unexplained disagreement would witness nothing.
    """
    theirs = v2_tables()['radii']
    assert len(theirs) == 118
    assert theirs['3'] == 167
    for e in ELEMENTS:
        if e.z == 3:
            assert e.atomic_radius == 1.67
        else:
            assert e.atomic_radius == theirs[str(e.z)], e.symbol
