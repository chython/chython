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
"""The atropisomer axis rule: the summed ortho score against the ring pair's bar (ATROPISOMERS in
`_stereo.pxi`, the score and bar tables in `docs/stereo.rst`).

Every case is a candidate question -- is there a kind-3 unit at all -- and the score is spelled beside
it.  Stereogenicity is `mark_stereogenic`'s and is asserted only where a case is about it.
"""
import pytest

from chython.core import read_smiles

SU_ATROPISOMER = 3


def _axes(text):
    return [u for u in read_smiles(text).stereo_units() if u['kind'] == SU_ATROPISOMER]


@pytest.mark.parametrize('name,text,axis', [
    # 6C-6C, bar 5
    ('biphenyl', 'c1ccccc1-c1ccccc1', False),                                       # 0
    ("2,2'-dimethylbiphenyl", 'Cc1ccccc1-c1ccccc1C', False),                       # 2 + 2
    ("2-chloro-2'-iodobiphenyl", 'Ic1ccccc1-c1ccccc1Cl', True),                    # 3 + 2
    ("2-tert-butyl-2'-methylbiphenyl", 'CC(C)(C)c1ccccc1-c1ccccc1C', True),        # 3 + 2
    ("2-methoxy-6-methyl-2'-chlorobiphenyl", 'COc1cccc(C)c1-c1ccccc1Cl', True),    # 1 + 2 + 2
    ("2,6-dimethyl-2'-chlorobiphenyl", 'Cc1cccc(C)c1-c1ccccc1Cl', True),           # 2 + 2 + 2
    ('1,1\'-binaphthyl', 'c1ccc2c(c1)cccc2-c1cccc2ccccc12', True),                  # peri 3 + 3
    ('BINOL', 'Oc1ccc2ccccc2c1-c1c(O)ccc2ccccc12', True),                           # 3 + 2 + 3 + 2
    # 5N-6C, N-aryl azoles, bar 9
    ('1-phenylpyrrole', 'c1ccn(c1)-c1ccccc1', False),                               # 0
    ('1-(2-chloro-6-methylphenyl)-2,5-dimethylpyrrole', 'Cc1ccc(C)n1-c1c(C)cccc1Cl', False),        # 4 + 4
    ('1-(2-isopropyl-6-methylphenyl)-2,5-dimethylpyrrole', 'CC(C)c1cccc(C)c1-n1c(C)ccc1C', True),   # 4 + 5
    ('1-(2-methylphenyl)pyrazole', 'Cc1ccccc1-n1cccn1', False),                     # 0 + 2
    ('1-(2-tert-butyl-6-methylphenyl)-3,5-dimethylpyrazole',
     'Cc1cc(C)n(-c2c(C)cccc2C(C)(C)C)n1', False),                                   # 2 + 5: N2 is a lone pair
])
def test_the_ortho_score_against_the_ring_pair_bar(name, text, axis):
    assert bool(_axes(text)) is axis, name


@pytest.mark.parametrize('name,text', [
    ('fluorene', 'C1c2ccccc2-c2ccccc21'),
    ('9,10-dihydrophenanthrene', 'C1Cc2ccccc2-c2ccccc21'),
    ('4,5-dimethoxyfluorene', 'COc1cccc2c1-c1c(OC)cccc1C2'),
])
def test_a_bridge_of_four_atoms_or_fewer_makes_no_axis(name, text):
    """A ring of six or fewer through the pivot bond holds the two rings coplanar, however hindered."""
    assert _axes(text) == [], name


@pytest.mark.parametrize('name,text', [
    ('a dibenzocycloheptene, both bay positions methoxy', 'COc1cccc2c1-c1c(OC)cccc1CCC2'),
    ('a dibenzocyclooctadiene, both bay positions methoxy', 'COc1cccc2c1-c1c(OC)cccc1CCCC2'),
    ('a dibenzocyclooctadiene, both bay positions hydroxy', 'Oc1cccc2c1-c1c(O)cccc1CCCC2'),
    ('a twelve-ring strapped 2,2\'-dimethylbiphenyl', 'Cc1cccc2c1-c1c(C)cccc1CCCCCCCC2'),
])
def test_a_bridge_of_five_atoms_or_more_is_an_ortho_substituent(name, text):
    """The bond lies in a ring of seven or more and is judged as acyclic: the bridge's CH2 is a 3."""
    axes = _axes(text)
    assert len(axes) == 1 and axes[0]['stereogenic'], name


def test_a_macrocyclic_biaryl_bond_written_aromatic_is_still_an_axis():
    """OpenSMILES reads an unmarked bond between two aromatic atoms as aromatic, ring bond or not, so
    the bridged biaryl's pivot bond is order 4 here and the axis does not depend on that spelling."""
    text = 'COc1cccc2c1c1c(OC)cccc1CCCC2'
    assert len(_axes(text)) == 1


@pytest.mark.parametrize('aromatic,kekule', [
    ('Cc1ccccc1-c1ccccc1C', 'CC1=CC=CC=C1C1=CC=CC=C1C'),
    ('Oc1ccc2ccccc2c1-c1c(O)ccc2ccccc12', 'OC1=CC=C2C=CC=CC2=C1C1=C(O)C=CC2=CC=CC=C12'),
    ('Ic1ccccc1-c1ccccc1Cl', 'IC1=CC=CC=C1C1=CC=CC=C1Cl'),
])
def test_the_score_does_not_depend_on_the_kekule_form(aromatic, kekule):
    assert len(_axes(aromatic)) == len(_axes(kekule))


def test_a_candidate_with_a_symmetric_end_is_not_stereogenic():
    """2,6-dimethyl-2'-chlorobiphenyl clears the bar and is still no stereocentre: turning the
    dimethylphenyl ring over reproduces the molecule."""
    axis, = _axes('Cc1cccc(C)c1-c1ccccc1Cl')
    assert not axis['stereogenic']
