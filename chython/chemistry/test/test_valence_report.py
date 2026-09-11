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
"""`check_valence`: a report whose two verdicts, `violation` and `unknown`, are not the same claim."""
from chython.chemistry._implicit import check_valence
from chython.core import H_UNKNOWN, read_smiles


def test_a_valid_molecule_reports_nothing():
    for string in ['CCO', 'C1=CC=CC=C1', 'C[N+](C)(C)C', 'O=S(=O)(O)O', '[Na+].[Cl-]']:
        assert check_valence(read_smiles(string)) == [], string


def test_a_violation_names_the_atom_and_only_that_atom():
    """Five bonds on a neutral nitrogen: the collection describes neutral N and no row accepts it."""
    mol = read_smiles('C[N](C)(C)C')
    assert check_valence(mol) == [(2, 'violation')]


def test_an_aromatic_molecule_is_CHECKED_and_not_shrugged_at():
    """An unkekulised ring is checkable, and kekulising must not CHANGE the verdict, so both are
    asserted.

    Order 4 needs no valence row of its own: `arom_classify_atom` decides whether the atom takes a
    ring double bond, and then the ordinary rows answer.  Kekulising only removes that step.
    """
    for string in ['c1ccccc1', 'c1ccncc1', 'c1cc[nH]c1', 'c1ccsc1', 'c1ccoc1', 'C[n+]1ccccc1',
                   'c1ccc2ccccc2c1', 'c1ccc2[nH]cnc2c1']:
        mol = read_smiles(string)
        assert check_valence(mol) == [], string
        mol.kekule()
        assert check_valence(mol) == [], f'{string} after kekule'


def test_unknown_survives_only_where_no_complete_QUESTION_can_be_put():
    """The two things `unknown` means.  Neither is a claim about the molecule, which is why it is not
    `'violation'`.
    """
    # (1) the collection describes nothing there: `[S+6]` is storable and undescribed (see
    # `_valence.pxi`'s header), so no row accepts it and none rejects it either.
    assert check_valence(read_smiles('C[S+6](C)C')) == [(2, 'unknown')]

    # (2) aromatic-ambiguous AND no count.  Pyrrole-versus-pyridine is the ring's choice, so there
    # are two candidate bond-order sums and picking one would invent the answer.  Erasing the count
    # is how a SMILES reaches the state a CTfile arrives in: MDL has no channel for an aromatic
    # nitrogen's hydrogens.
    mol = read_smiles('c1ccncc1')
    mol.set_hydrogens(4, H_UNKNOWN)
    assert check_valence(mol) == [(4, 'unknown')]

    # both halves are needed: erase a count whose aromatic class the classifier CAN settle and the
    # skeleton is still checkable, since benzene's carbon takes a ring double bond whatever its
    # hydrogens, so the sum is known and the collection has a row.
    mol = read_smiles('c1ccccc1')
    mol.set_hydrogens(1, H_UNKNOWN)
    assert check_valence(mol) == []


def test_an_unknown_count_is_not_read_as_a_count():
    """The sentinel is 15 and a hydrogen count is four bits, so it is a NUMBER unless someone looks.

    Passing the raw nibble asks the collection for a row with fifteen hydrogens and turns every
    unrecorded count into a violation; passing 0 claims the atom has no hydrogens.  Neither is a
    report about the molecule, so the honest verdict for `CCC` with an erased count is silence.
    """
    mol = read_smiles('CCC')
    mol.set_hydrogens(2, H_UNKNOWN)
    assert mol.implicit_h_of(2) is None
    assert check_valence(mol) == []


def test_the_report_never_edits_the_molecule():
    """It answers a question about what is stored, so storing something else would be a lie."""
    for string in ['C[N](C)(C)C', 'c1ccccc1', 'c1cc[nH]cc1']:
        mol = read_smiles(string)
        before = mol.canonical_bytes
        check_valence(mol)
        assert mol.canonical_bytes == before, string


def test_the_method_on_the_container_is_the_registered_report():
    """Registration is by INJECTION -- `chemistry` calls `_set_valence_fn`, `core` names no chemistry.

    `mol.kekule(); mol.check_valence()` is the triage sequence for a corpus, and it must not require
    knowing which package the verdicts live in.
    """
    mol = read_smiles('C[N](C)(C)C')
    assert mol.check_valence() == check_valence(mol) == [(2, 'violation')]
    assert read_smiles('CCO').check_valence() == []


def test_an_r_neighbour_does_not_make_a_sulfone_a_violation():
    # The collection enumerates hypervalent sulfur by neighbour element -- `-C -C =O =O` is dimethyl
    # sulfone -- and a marker reads as carbon for its neighbour, so the sulfur matches that row.
    assert check_valence(read_smiles('[R]S(C)(=O)=O')) == [(1, 'unknown')]


def test_the_marker_itself_is_still_unknown():
    # Element 0 is described by no row, and that is the honest answer for it.  Reading as carbon is
    # the NEIGHBOUR's rule, so it must not make the marker itself look like a described state.
    verdicts = dict(check_valence(read_smiles('[R]C')))
    assert verdicts == {1: 'unknown'}


def test_every_sulfur_oxidation_state_takes_the_marker():
    for probe in ('[R]S(C)=O', '[R]S(C)(=O)=O', '[R]S(=O)(=O)N', '[R]S(=O)(=O)Cl'):
        assert check_valence(read_smiles(probe)) == [(1, 'unknown')], probe


def test_a_carbon_in_the_markers_place_answers_the_same():
    # The rule is an equivalence, so state it as one: the only difference the marker may make to a
    # neighbour's verdict is its own `unknown` entry.
    for marked, plain in (('[R]S(C)(=O)=O', 'CS(C)(=O)=O'), ('[R]S(=O)(=O)Cl', 'CS(=O)(=O)Cl'),
                          ('[R][Si](C)(C)C', 'C[Si](C)(C)C'), ('[R]P(C)(C)=O', 'CP(C)(C)=O')):
        assert [v for i, v in check_valence(read_smiles(marked)) if v != 'unknown'] == \
               [v for i, v in check_valence(read_smiles(plain))], marked
