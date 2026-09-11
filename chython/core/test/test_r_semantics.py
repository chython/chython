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
"""R matches nothing, reads as carbon for a neighbour's features, and is distinct for identity."""
from pytest import raises
from chython.core import MoleculeContainer, molecule_to_inchi, molecule_to_inchikey, read_smarts, read_smiles


def _toluene_with_r():
    """4-R-toluene: the R replaces the para hydrogen of toluene."""
    mol = MoleculeContainer()
    with mol.edit() as e:
        ring = [e.add_atom('C') for _ in range(6)]
        for i in range(6):
            e.add_bond(ring[i], ring[(i + 1) % 6], 4)
        methyl = e.add_atom('C')
        e.add_bond(ring[0], methyl, 1)
        r = e.add_atom('R')
        e.add_bond(ring[3], r, 1)
    mol.kekule()      # a hand-built molecule has no implicit H counts until kekule's fill runs
    mol.thiele()
    return mol, ring, methyl, r


def test_any_atom_query_does_not_reach_an_r():
    # [A] matches every non-R heavy atom; the R itself must be excluded.
    mol, ring, methyl, r = _toluene_with_r()
    hits = [m for m in read_smarts('[A]').get_mapping(mol)]
    matched = {next(iter(h.values())) for h in hits}
    assert r not in matched
    assert ring[0] in matched
    assert len(hits) == 7


def test_carbon_query_does_not_reach_an_r():
    # [C] is an element match; R (element 0) is not carbon and must not match.
    mol, ring, methyl, r = _toluene_with_r()
    hits = [m for m in read_smarts('[C]').get_mapping(mol)]
    assert not any(r in m.values() for m in hits)
    assert len(hits) == 7


def test_wildcard_star_does_not_reach_an_r():
    # [A,M] matches any organic or metal atom; R is a marker and must be excluded from both.
    mol, ring, methyl, r = _toluene_with_r()
    hits = [m for m in read_smarts('[A,M]').get_mapping(mol)]
    assert not any(r in m.values() for m in hits)
    assert len(hits) == 7


def test_a_query_through_the_r_bearing_atom_still_matches_its_ring():
    # [C;a:1] matches aromatic carbons including the one bearing the R substituent.
    mol, ring, methyl, r = _toluene_with_r()
    hits = list(read_smarts('[C;a:1]').get_mapping(mol))
    assert any(m[1] == ring[3] for m in hits)
    assert len(hits) == 6


def test_neighbour_of_an_r_counts_it_as_carbon():
    mol, ring, methyl, r = _toluene_with_r()
    carrier = mol.atom(ring[3])
    # Three heavy neighbours: two ring carbons and the R.  No implicit hydrogen left.
    assert carrier.neighbors == 3
    assert carrier.implicit_h == 0
    assert carrier.heteroatoms == 0


def test_the_same_shape_with_a_carbon_gives_the_same_neighbour_features():
    mol, ring, methyl, r = _toluene_with_r()
    para_xylene = read_smiles('Cc1ccc(C)cc1')
    r_atom = mol.atom(ring[3])
    found = False
    for atom in para_xylene.atoms():
        if atom.neighbors == 3 and atom.hybridization == 4:   # 4 is aromatic
            assert (atom.neighbors, atom.implicit_h, atom.heteroatoms) == \
                   (r_atom.neighbors, r_atom.implicit_h, r_atom.heteroatoms)
            found = True
    assert found, 'p-xylene has no substituted aromatic carbon'


def test_the_r_is_not_a_heteroatom_of_the_molecule():
    mol, ring, methyl, r = _toluene_with_r()
    assert mol.heteroatoms_count == 0


def test_an_r_next_to_nitrogen_leaves_the_nitrogen_one_hydrogen():
    # Aniline with the R on nitrogen: N keeps one H, as N-methylaniline's does.
    mol = MoleculeContainer()
    with mol.edit() as e:
        ring = [e.add_atom('C') for _ in range(6)]
        for i in range(6):
            e.add_bond(ring[i], ring[(i + 1) % 6], 4)
        n = e.add_atom('N')
        e.add_bond(ring[0], n, 1)
        r = e.add_atom('R')
        e.add_bond(n, r, 1)
    mol.kekule()
    mol.thiele()
    assert mol.atom(n).implicit_h == 1
    assert mol.atom(n).heteroatoms == 0


def test_an_r_neighbour_gives_the_same_hydrogen_count_as_a_carbon_one():
    """The hydrogen half of the rule, over the elements whose valence rows differ.

    `Ph-X-R` against `Ph-X-C`: the marker must not move X's count.  As and Sn are excluded because
    both sides report an unknown count there, which says nothing about the R.
    """
    def scaffold(element, substituent):
        mol = MoleculeContainer()
        with mol.edit() as e:
            ring = [e.add_atom('C') for _ in range(6)]
            for i in range(6):
                e.add_bond(ring[i], ring[(i + 1) % 6], 4)
            x = e.add_atom(element)
            e.add_bond(ring[0], x, 1)
            e.add_bond(x, e.add_atom(substituent), 1)
        mol.kekule()
        return mol.atom(x)

    for element, hydrogens in (('C', 2), ('N', 1), ('O', 0), ('S', 0), ('P', 1), ('B', 1),
                               ('Si', 2), ('Se', 0), ('Al', 1), ('Ge', 2), ('Te', 0)):
        with_r = scaffold(element, 'R')
        with_c = scaffold(element, 'C')
        assert with_r.implicit_h == with_c.implicit_h == hydrogens, element
        assert with_r.heteroatoms == with_c.heteroatoms == 0, element


def test_r_in_a_ring_is_not_a_heterocycle():
    # A six-membered saturated ring with one R member: the R is not a heteroatom and the ring is not
    # heterocyclic.  Both descriptors must agree.
    mol = MoleculeContainer()
    with mol.edit() as e:
        members = [e.add_atom('C') for _ in range(5)]
        members.append(e.add_atom('R'))
        for i in range(6):
            e.add_bond(members[i], members[(i + 1) % 6], 1)
    mol.kekule()      # a hand-built molecule has no implicit H counts until kekule's fill runs
    mol.thiele()
    assert mol.rings_count == 1
    assert mol.heteroatoms_count == 0
    assert mol.heterocycles_count == 0


def _ring_with(*symbols):
    """Benzene carrying one substituent per given symbol, in ring order from position 1."""
    mol = MoleculeContainer()
    with mol.edit() as e:
        ring = [e.add_atom('C') for _ in range(6)]
        for i in range(6):
            e.add_bond(ring[i], ring[(i + 1) % 6], 4)
        for position, symbol in enumerate(symbols):
            sub = e.add_atom(symbol)
            e.add_bond(ring[position], sub, 1)
    return mol


def test_an_r_is_not_a_carbon_in_the_canonical_form():
    assert _ring_with('R').canonical_bytes != _ring_with('C').canonical_bytes


def test_two_r_indices_are_distinct():
    assert _ring_with('R1').canonical_bytes != _ring_with('R2').canonical_bytes


def test_the_same_index_in_the_same_place_is_the_same_molecule():
    assert _ring_with('R1').canonical_bytes == _ring_with('R1').canonical_bytes


def _para_disubstituted(first, second):
    """Benzene with `first` at position 1 and `second` at position 4, returned with both atom ids."""
    mol = MoleculeContainer()
    with mol.edit() as e:
        ring = [e.add_atom('C') for _ in range(6)]
        for i in range(6):
            e.add_bond(ring[i], ring[(i + 1) % 6], 4)
        a = e.add_atom(first)
        e.add_bond(ring[0], a, 1)
        b = e.add_atom(second)
        e.add_bond(ring[3], b, 1)
    return mol, a, b


def test_two_identical_r_groups_stay_symmetric():
    # 1,4-di-R1-benzene keeps the automorphism a 1,4-disubstituted ring has: one index, one orbit.
    mol, a, b = _para_disubstituted('R1', 'R1')
    orbits = mol.automorphism_orbits()
    assert orbits[a] == orbits[b]


def test_two_different_r_indices_break_the_symmetry():
    # No automorphism swaps an R1 with an R2, so the two markers sit in different orbits.
    mol, a, b = _para_disubstituted('R1', 'R2')
    orbits = mol.automorphism_orbits()
    assert orbits[a] != orbits[b]


def test_an_r_carries_no_hydrogens_of_its_own():
    mol, ring, methyl, r = _toluene_with_r()
    assert mol.atom(r).implicit_h == 0
    assert mol.atom(r).total_h == 0
    # Derived and not unknown: nothing about an attachment point is underivable, so the formula is a
    # formula rather than a lower bound.
    assert mol.unknown_h_count == 0


def test_an_r_carries_no_mass():
    mol, ring, methyl, r = _toluene_with_r()
    toluene = read_smiles('Cc1ccccc1')
    # The R stands where toluene's para hydrogen stands, and contributes nothing of its own.
    assert round(float(toluene) - float(mol), 3) == 1.008


def test_element_counts_counts_the_marker():
    mol, ring, methyl, r = _toluene_with_r()
    # Keyed by atomic number, and an R is element 0.  The counts sum to the atom count.
    assert mol.element_counts == {6: 7, 0: 1}
    assert sum(mol.element_counts.values()) == len(list(mol.atoms()))


def test_brutto_names_the_marker():
    mol, ring, methyl, r = _toluene_with_r()
    assert mol.brutto == {'C': 7, 'H': 7, 'R': 1}
    # `brutto`'s order is the answer, and the marker sorts after every element.
    assert mol.brutto_formula == 'C7H7R'


def test_inchi_refuses_a_marker():
    mol, ring, methyl, r = _toluene_with_r()
    with raises(ValueError, match='R'):
        molecule_to_inchi(mol)
    with raises(ValueError, match='R'):
        molecule_to_inchikey(mol)


# `'R'` asks about ANY R and `'R7'` about one index -- the same two grains `atom == 'C'` has for the
# isotope.  `SYMBOL_TO_NUMBER['R'] == 0`, so a lookup that reads 0 as "unknown" answers False on both.

def test_the_bare_symbol_finds_any_indexed_marker():
    assert 'R' in read_smiles('[R1]c1ccccc1')


def test_the_indexed_symbol_finds_its_own_index():
    assert 'R1' in read_smiles('[R1]c1ccccc1')


def test_the_indexed_symbol_does_not_find_another_index():
    assert 'R2' not in read_smiles('[R1]c1ccccc1')


def test_the_bare_symbol_is_absent_from_a_molecule_with_no_marker():
    assert 'R' not in read_smiles('Cc1ccccc1')


def test_r0_is_not_a_spelling():
    # An unindexed R spells `R`, so `'R0'` is an unknown symbol with an answer rather than an error.
    assert 'R0' not in read_smiles('[R]C')


def test_an_index_past_the_domain_is_an_unknown_symbol():
    assert 'R100' not in read_smiles('[R]C')


def test_an_unknown_symbol_is_still_false():
    assert 'Xx' not in read_smiles('[R1]c1ccccc1')


def test_an_atom_equals_the_bare_symbol_and_its_own_index():
    atom = read_smiles('[R1]C').atom(1)
    assert atom == 'R'
    assert atom == 'R1'


def test_an_atom_does_not_equal_another_index():
    assert read_smiles('[R1]C').atom(1) != 'R2'


def test_an_unindexed_marker_equals_the_bare_symbol():
    assert read_smiles('[R]C').atom(1) == 'R'


def test_a_carbon_does_not_equal_the_marker():
    assert read_smiles('Cc1ccccc1').atom(1) != 'R'


def test_an_atom_does_not_equal_an_unknown_symbol():
    assert read_smiles('[R1]C').atom(1) != 'Xx'


def test_every_atom_equals_its_own_atomic_symbol():
    """The invariant both grains exist to keep: an R with no index spells `R`, an indexed one `R7`."""
    for smiles in ('[R]C', '[R1]C', '[R99]C', 'Cc1ccccc1'):
        for atom in read_smiles(smiles).atoms():
            assert atom == atom.atomic_symbol, (smiles, atom.atomic_symbol)
