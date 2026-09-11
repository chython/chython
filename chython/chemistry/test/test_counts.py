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
from pytest import mark, raises
from chython.chemistry import hydrogen_bond_acceptors_count, hydrogen_bond_donors_count, rotatable_bonds_count
from chython.chemistry._counts import hbond_atoms
from chython.core import read_smiles


ROTATABLE = [
    ('CCCC', 1),                              # butane: the central C-C only
    ('CC', 0),                                # ethane: both ends are terminal
    ('CC(C)Cc1ccc(cc1)C(C)C(=O)O', 4),        # ibuprofen
    ('CC(=O)NCCc1ccccc1', 3),                 # N-phenethylacetamide
    ('CN(C)C(=O)N(C)C', 0),                   # tetramethylurea: both amides excluded
    ('[O-][N+](=O)c1ccccc1', 1),              # nitrobenzene: the aryl-N bond.  Charged atoms are
                                              # counted, deliberately
    ('c1ccccc1', 0),                          # benzene: every bond is in a ring
    ('C1CCCCC1', 0),                          # cyclohexane: likewise
    ('CS(=O)(=O)N(C)C', 0),                   # a sulfonamide N-S, excluded
    ('CC#CC', 0),                             # 2-butyne: the sp carbons are D2 but the bond is #
    ('OCCO', 1),                              # ethylene glycol: only the central C-C
    ('c1ccccc1-c1ccccc1', 1),                 # biphenyl: the inter-ring bond.  The `-` is required,
                                              # see test_the_ar_ar_spelling_needs_the_repair_first
]


@mark.parametrize('smi,expected', ROTATABLE)
def test_rotatable_bonds_count_counts_bonds_not_mappings(smi, expected):
    assert rotatable_bonds_count(read_smiles(smi)) == expected


def test_the_container_property_agrees_with_the_function():
    m = read_smiles('CC(C)Cc1ccc(cc1)C(C)C(=O)O')
    assert m.rotatable_bonds_count == rotatable_bonds_count(m) == 4


def test_the_symmetric_pattern_maps_each_bond_twice_and_the_count_deduplicates():
    """The `found` set is load-bearing: row 1 is symmetric in `:1` and `:2` and the matcher emits both
    directions, so a pass counting the mapping stream would return exactly double.
    """
    from chython.chemistry._tables import rotatable_rules_by_role
    m = read_smiles('CC(C)Cc1ccc(cc1)C(C)C(=O)O')
    row, = rotatable_rules_by_role()['rotatable']
    assert sum(1 for _ in row.query.get_mapping(m)) == 8
    assert m.rotatable_bonds_count == 4


def test_the_ar_ar_spelling_needs_the_repair_pipeline_first():
    """An unspecified bond between two aromatic atoms is aromatic per OpenSMILES, so biphenyl written
    without the `-` has no acyclic single bond until the caller runs the repair pipeline.  Both
    spellings then converge to `c1cc(-c2ccccc2)ccc1`.
    """
    m = read_smiles('c1ccccc1c1ccccc1')
    assert m.rotatable_bonds_count == 0
    m.kekule()
    m.thiele()
    assert m.rotatable_bonds_count == 1


def test_the_count_is_renumbering_invariant():
    a = read_smiles('CC(=O)NCCc1ccccc1')
    b = read_smiles('c1ccccc1CCNC(C)=O')
    assert a.rotatable_bonds_count == b.rotatable_bonds_count == 3


HBOND = [
    # smiles,                                   donors, acceptors, name
    ('CCO', 1, 1),                              # ethanol
    ('CC(=O)O', 1, 2),                          # acetic acid: OH donates, both O accept
    ('CC(=O)OC', 0, 2),                         # methyl acetate: both O accept, nothing donates
    ('CC(=O)N', 1, 1),                          # acetamide: NH2 donates; the amide N does not accept
    ('CC(=O)NC', 1, 1),                         # N-methylacetamide: likewise
    ('c1ccccc1O', 1, 1),                        # phenol
    ('c1cc[nH]c1', 1, 0),                       # pyrrole: donates, does not accept
    ('c1ccncc1', 0, 1),                         # pyridine: accepts, does not donate
    ('CN', 1, 1),                               # methylamine
    ('C[N+](C)(C)C', 0, 0),                     # tetramethylammonium: neither
    ('[O-][N+](=O)c1ccccc1', 0, 0),             # nitrobenzene: a nitro oxygen accepts nothing
    ('C[N+](C)(C)[O-]', 0, 1),                  # trimethylamine N-oxide: the oxide oxygen accepts
    ('c1cc[n+]([O-])cc1', 0, 1),                # pyridine N-oxide: the oxide oxygen, not the ring N
    ('CS(=O)(=O)N', 1, 2),                      # methanesulfonamide: NH2 donates, the two S=O accept
    ('CC#N', 0, 1),                             # acetonitrile
    ('CCl', 0, 0),                              # a halogen is not an acceptor
    ('NC(=O)c1ccncc1', 1, 2),                   # isonicotinamide: C=O and the ring N accept
    ('Oc1ccccc1C(=O)O', 2, 3),                  # salicylic acid
    ('CC(=O)Nc1ccc(O)cc1', 2, 2),               # paracetamol
    ('NC(N)=O', 2, 1),                          # urea: two NH2 donate, only the O accepts
    ('O', 1, 1),                                # water
    ('CC(C)Cc1ccc(cc1)C(C)C(=O)O', 1, 2),       # ibuprofen
]


@mark.parametrize('smi,donors,acceptors', HBOND)
def test_hydrogen_bond_counts(smi, donors, acceptors):
    m = read_smiles(smi)
    assert m.hydrogen_bond_donors_count == donors, 'donors'
    assert m.hydrogen_bond_acceptors_count == acceptors, 'acceptors'


def test_the_amide_nitrogen_donates_and_does_not_accept():
    # not excluded by hybridization: an amide N is z1, exactly like an amine N
    m = read_smiles('CC(=O)NC')
    n = next(a.n for a in m.atoms() if a.element == 7)
    assert m.atom(n).hybridization == 1, 'the exclusion cannot be a z test; see step 3'
    assert n in hbond_atoms(m, 'donor')
    assert n not in hbond_atoms(m, 'acceptor')


def test_the_sulfonamide_nitrogen_donates_and_does_not_accept():
    # the other half of the same exclusion: the neighbour that disqualifies it is S, not a carbonyl
    m = read_smiles('CS(=O)(=O)N')
    n = next(a.n for a in m.atoms() if a.element == 7)
    assert n in hbond_atoms(m, 'donor')
    assert n not in hbond_atoms(m, 'acceptor')


def test_the_cations_donate_up_to_four_hydrogens_and_accept_none():
    # row 1's `H` list must run to H4: `[NH4+]` measures h=4, so a list stopping at H3 silently reads
    # the textbook donor cation as donating nothing.  The acceptor half is the charge rule -- an
    # unstated charge means neutral, so no acceptor row can reach a cation.
    for smi, donors in (('[NH4+]', 1), ('C[NH3+]', 1), ('CC[NH2+]C', 1)):
        m = read_smiles(smi)
        assert m.hydrogen_bond_donors_count == donors, smi
        assert m.hydrogen_bond_acceptors_count == 0, smi


def test_the_water_oxygen_is_both_and_has_no_heavy_neighbour():
    # D0: the row that types it may not demand a carbon neighbour, or water counts zero acceptors
    m = read_smiles('O')
    o = next(a.n for a in m.atoms() if a.element == 8)
    assert m.atom(o).degree == 0
    assert o in hbond_atoms(m, 'donor') and o in hbond_atoms(m, 'acceptor')


def test_a_dative_oxide_oxygen_accepts_and_a_nitro_oxygen_still_does_not():
    """The one oxide `standardize()` cannot spell neutrally, told apart from nitro by the cation's `x`.

    A charge-separated sulfoxide is the pipeline's to fix -- it neutralises to `CS(=O)C`, which row 2
    already types -- so it needs no row of its own to be answered about.  An amine or pyridine N-oxide
    has no neutral spelling at all: its oxygen reaches an acceptor row or nothing ever types it.  The
    separating primitive is `x` on the cationic centre, which counts one heteroatom neighbour for an
    oxide and two for a nitro group, so the table header's "a nitro O accepts nothing" survives
    unchanged rather than being re-argued.
    """
    for smi in ('C[N+](C)(C)[O-]', 'c1cc[n+]([O-])cc1', 'C[P+](C)(C)[O-]', 'C[S+]([O-])C'):
        m = read_smiles(smi)
        o = next(a.n for a in m.atoms() if a.element == 8)
        assert o in hbond_atoms(m, 'acceptor'), smi
        assert o not in hbond_atoms(m, 'donor'), smi


def test_a_cation_bearing_two_heteroatoms_is_left_out_of_the_oxide_row():
    """`x1` and not `x1,x2`, which is what keeps the nitro decision the header states.

    Azoxy rides along with it, and deliberately: it is the same shape -- a cationic nitrogen sharing its
    charge with a second heteroatom -- and the conservative reading is the one already chosen for nitro.
    """
    for smi in ('c1ccccc1[N+](=O)[O-]', 'C[N+](=NC)[O-]'):
        assert hbond_atoms(read_smiles(smi), 'acceptor') == frozenset(), smi


def test_hbond_atoms_returns_stable_ids_and_is_a_frozenset():
    m = read_smiles('CC(=O)O')
    ids = hbond_atoms(m, 'acceptor')
    assert isinstance(ids, frozenset)
    assert ids <= set(m.atom_numbers)


def test_an_unknown_role_is_refused_rather_than_answered_as_empty():
    with raises(ValueError, match='donor'):
        hbond_atoms(read_smiles('CCO'), 'halogen')


def test_counts_are_renumbering_invariant():
    assert (read_smiles('Oc1ccccc1C(=O)O').hydrogen_bond_donors_count
            == read_smiles('OC(=O)c1ccccc1O').hydrogen_bond_donors_count == 2)


def test_hbond_table_shape():
    from chython.chemistry._tables import HBOND_ROLES, hbond_rules, hbond_rules_by_role
    rows = hbond_rules()
    assert len(rows) == 13
    assert len(hbond_rules_by_role()['donor']) == 1
    assert len(hbond_rules_by_role()['acceptor']) == 12
    assert {r.role for r in rows} <= set(HBOND_ROLES)
    assert all(r.id.startswith('hbond:') for r in rows)
    assert all(r.description and r.description != '-' for r in rows)
    # every row must name its subject atom explicitly, or it types whichever atom the auto-numbering
    # reached first.  `map_numbers()` reports only the numbers actually written, which is what makes
    # this failable; `1 in r.numbers` cannot fail, `compile_smarts` numbering every atom from 1 up.
    assert all(1 in set(r.query.map_numbers().values()) for r in rows)


def test_the_subject_check_can_fail():
    # negative control for the check above, which was unfailable when written against `r.numbers`
    from chython.core import read_smarts
    assert 1 not in set(read_smarts('[N;*]-[C;*]').map_numbers().values())
    assert 1 in set(read_smarts('[N;*:1]-[C;*]').map_numbers().values())


def test_no_acceptor_row_admits_a_nitrogen_partner_on_the_carbonyl_row():
    # row 2's partner list is `[C,S,P;*]`; adding N makes nitrobenzene read one acceptor instead of 0.
    from chython.chemistry._tables import hbond_rules_by_role
    row = hbond_rules_by_role()['acceptor'][0]
    assert row.pattern == '[O;D1;z2;h0;!^:1]=[C,S,P;*]', row.pattern
