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
from pytest import approx, mark
from chython.chemistry import tpsa
from chython.chemistry._tables import TPSA_CLASSES
from chython.chemistry._tpsa import _POLAR, tpsa_contributions
from chython.core import read_smiles


# every expected value is the sum of tables/tpsa.tsv contributions, shown in the comment.
# Ertl, Rohde, Selzer, J. Med. Chem. 2000, 43, 3714.
CASES = [
    ('C', 0.0),                                     # methane: no polar atom
    ('CCCC', 0.0),                                  # butane
    ('c1ccccc1', 0.0),                              # benzene
    ('CCO', 20.23),                                 # ethanol: OH 20.23
    ('CC(=O)O', 37.30),                             # acetic acid: 17.07 + 20.23
    ('CC(=O)OC', 26.30),                            # methyl acetate: 17.07 + 9.23
    ('CN', 26.02),                                  # methylamine: NH2 26.02
    ('CNC', 12.03),                                 # dimethylamine: NH 12.03
    ('CN(C)C', 3.24),                               # trimethylamine
    ('C[N+](C)(C)C', 0.00),                         # tetramethylammonium: 0.00 by the table
    ('C[NH3+]', 27.64),                             # methylammonium
    ('CC#N', 23.79),                                # acetonitrile
    ('c1ccncc1', 12.89),                            # pyridine
    ('c1cc[nH]c1', 15.79),                          # pyrrole
    ('c1ccoc1', 13.14),                             # furan
    ('NC(N)=O', 69.11),                             # urea: 26.02 + 26.02 + 17.07
    ('CC(=O)Nc1ccc(O)cc1', 49.33),                  # paracetamol: 12.03 + 17.07 + 20.23
    ('CC(=O)Oc1ccccc1C(=O)O', 63.60),               # aspirin: 9.23 + 17.07 + 17.07 + 20.23
    ('CC(C)Cc1ccc(cc1)C(C)C(=O)O', 37.30),          # ibuprofen: 17.07 + 20.23
    ('Oc1ccccc1C(=O)O', 57.53),                     # salicylic acid: 20.23 + 17.07 + 20.23
    ('C1CO1', 12.53),                               # ethylene oxide: the epoxide row, not the ether
    ('C1CN1', 21.94),                               # aziridine: the ring NH row
    ('Nc1ccc(cc1)S(N)(=O)=O', 86.18),               # sulfanilamide: 26.02 + 26.02 + 17.07 + 17.07
    ('COP(=O)(OC)OC', 44.76),                       # trimethyl phosphate: 3 x 9.23 + 17.07
]


@mark.parametrize('smi,expected', CASES)
def test_tpsa_matches_the_table_sum(smi, expected):
    assert tpsa(read_smiles(smi)) == approx(expected, abs=1e-9)


def test_the_container_property_is_ertls_published_n_and_o_sum():
    m = read_smiles('CS(=O)C')  # DMSO: the sulfoxide S and its O
    assert m.tpsa == approx(17.07)                               # O only, as published
    assert tpsa(m, sulfur_phosphorus=True) == approx(36.28)      # + S 19.21


def test_the_sulfur_and_phosphorus_extension_sums_its_own_rows():
    """The three SP rows no N+O total can reach: the sulfone (36) and the phosphate (42).

    Each number below is a published N+O total plus one SP row's contribution, so a mistyped row names
    itself.  Row 35, the sulfoxide, is pinned by the DMSO assertion above and not repeated here.
    """
    # dimethyl sulfone: 2 x 17.07, + the sulfone S 8.38 -- half the sulfoxide's, which is Ertl's table
    # and not a transcription slip: the sulfone sulfur is the more buried of the two.
    assert tpsa(read_smiles('CS(=O)(=O)C'), sulfur_phosphorus=True) == approx(42.52)
    # sulfanilamide: 86.18 published, + the same sulfone S
    assert tpsa(read_smiles('Nc1ccc(cc1)S(N)(=O)=O'), sulfur_phosphorus=True) == approx(94.56)
    # trimethyl phosphate: 44.76 published, + the phosphate P 9.81
    assert tpsa(read_smiles('COP(=O)(OC)OC'), sulfur_phosphorus=True) == approx(54.57)


def test_the_charge_separated_and_pentavalent_nitro_differ_and_both_are_right():
    # IO does not mutate representation: a descriptor answers about the molecule as drawn.
    assert tpsa(read_smiles('[O-][N+](=O)c1ccccc1')) == approx(43.14)   # 23.06 + 17.07 + 3.01
    assert tpsa(read_smiles('CN(=O)=O')) == approx(45.82)               # 11.68 + 17.07 + 17.07


def test_an_atom_matching_no_row_contributes_zero_and_is_logged():
    # water has no heavy neighbour, so no Ertl pattern matches its oxygen
    m = read_smiles('O')
    assert tpsa(m) == 0.0
    assert len(m.log) == 1
    assert m.log[0].rule == 'tpsa:unmatched'


def test_the_default_does_not_log_the_sulfur_it_was_never_asked_about():
    # `sulfur_phosphorus=False` loads no SP row, so every S and P is unmatched by construction; a record
    # for one would report the caller's own choice back as a defect.  The element set is derived from
    # `wanted`.
    m = read_smiles('CSC')
    assert tpsa(m) == 0.0                                 # dimethyl sulfide: no N, no O
    assert m.log == []
    # asked for SP, the same S is typed by the thioether row -- so still nothing to report
    assert tpsa(m, sulfur_phosphorus=True) == approx(25.30)
    assert m.log == []
    # positive control, without which the assertions above pass on a path that never fires at all:
    # H2S is `D0 h2` and matches no sulfur row in the table, so asked for SP it is reported.
    m = read_smiles('S')
    assert tpsa(m, sulfur_phosphorus=True) == 0.0
    assert len(m.log) == 1 and m.log[0].rule == 'tpsa:unmatched'


def test_the_contributions_decompose_the_total():
    m = read_smiles('CC(=O)Oc1ccccc1C(=O)O')
    parts = tpsa_contributions(m)
    assert len(parts) == 4                            # four oxygens, no carbon in the dict
    # 63.60 spelled out, not `approx(m.tpsa)`: `tpsa()` is `sum(tpsa_contributions(...))`, so comparing
    # the two sides asserts X == X and scaling every contribution by two would still pass.  The number
    # is aspirin's published row in CASES, and duplicating it here is the point.
    assert sum(parts.values()) == approx(63.60)
    assert set(parts) <= set(m.atom_numbers)


def test_every_table_class_has_its_polar_elements():
    """`_POLAR` is the one thing in `_tpsa.py` a new `tpsa.tsv` class cannot bring with it.

    Which elements a class types is not derivable from the rows, and an absent entry would make the
    reporting branch raise; this turns that into a test failure the moment a class is added.
    """
    assert set(_POLAR) == set(TPSA_CLASSES)
    assert _POLAR['NO'] == frozenset((7, 8)) and _POLAR['SP'] == frozenset((15, 16))


def test_it_is_renumbering_invariant():
    assert (read_smiles('CC(=O)Nc1ccc(O)cc1').tpsa
            == read_smiles('Oc1ccc(NC(C)=O)cc1').tpsa == approx(49.33))
