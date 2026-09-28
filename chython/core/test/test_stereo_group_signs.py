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
"""What a group's signs state, on named compounds (`docs/stereo.rst`, "Signs inside a group").

=====  ====================================================================================
Kind   Inverting every member's sign
=====  ====================================================================================
AND    the same racemate: one canonical phase, one string
OR     the other isolated isomer: signs kept as drawn
ABS    the enantiomer: signs kept as drawn, and an explicit label is not the same as none
=====  ====================================================================================

Only a mirror symmetry makes the inverted signs the same compound: a meso form under any kind, where the
label states nothing and the molecule is the unlabelled one.
"""
from pytest import mark

from chython.core._core import read_smiles as smiles


def inverted(text):
    """Every `@` sign of a SMILES inverted, the CXSMILES tail untouched."""
    return text.replace('@@', '#').replace('@', '@@').replace('#', '@')


def one_inverted(text):
    """The first `@` sign inverted and no other."""
    i = text.index('@')
    if text[i + 1] == '@':
        return text[:i] + text[i + 1:]
    return text[:i] + '@' + text[i:]


#: (name, SMILES, members, meso).  Each has two centres except butan-2-ol.
COMPOUNDS = [
    ('butan-2-ol', 'C[C@H](O)CC', '1', False),
    ('ephedrine', 'C[C@H](NC)[C@@H](O)c1ccccc1', '1,4', False),
    ('trans-cyclohexane-1,2-diol', 'O[C@@H]1CCCC[C@H]1O', '1,6', False),
    ('cis-cyclohexane-1,2-diol', 'O[C@@H]1CCCC[C@@H]1O', '1,6', True),
    ('(2R*,4R*)-pentane-2,4-diol', 'C[C@H](O)C[C@@H](O)C', '1,4', False),
    ('meso-pentane-2,4-diol', 'C[C@H](O)C[C@H](O)C', '1,4', True),
]
TWO_CENTRES = [c for c in COMPOUNDS if ',' in c[2]]


@mark.parametrize('name, text, members, meso', COMPOUNDS, ids=[c[0] for c in COMPOUNDS])
@mark.parametrize('kind', ['&1', 'o1', 'a'])
def test_inverting_every_member(kind, name, text, members, meso):
    drawn = smiles(f'{text} |{kind}:{members}|')
    other = smiles(f'{inverted(text)} |{kind}:{members}|')
    same = kind == '&1' or meso
    assert (drawn == other) is same, name
    assert (hash(drawn) == hash(other)) is same, name
    assert (str(drawn) == str(other)) is same, (name, str(drawn), str(other))


@mark.parametrize('name, text, members, meso', TWO_CENTRES, ids=[c[0] for c in TWO_CENTRES])
@mark.parametrize('kind', ['&1', 'o1', 'a'])
def test_inverting_one_member_is_the_other_diastereomer(kind, name, text, members, meso):
    """Ephedrine becomes pseudoephedrine, trans the cis diol, meso the racemic one: under every kind."""
    drawn = smiles(f'{text} |{kind}:{members}|')
    assert drawn != smiles(f'{one_inverted(text)} |{kind}:{members}|'), name


@mark.parametrize('name, text, members, meso', COMPOUNDS, ids=[c[0] for c in COMPOUNDS])
def test_the_three_kinds_and_no_label_are_four_identities(name, text, members, meso):
    """A meso form's label states nothing: all four forms are the unlabelled molecule."""
    forms = [smiles(text)] + [smiles(f'{text} |{kind}:{members}|') for kind in ('&1', 'o1', 'a')]
    for i, a in enumerate(forms):
        for b in forms[i + 1:]:
            assert (a == b) is meso, (name, str(a), str(b))
            assert (str(a) == str(b)) is meso, (name, str(a), str(b))


def test_an_and_pair_and_an_or_pair_in_one_molecule():
    """7,9-Diaminodecane-2,4-diol's four centres: the AND pair phases, the OR pair does not."""
    text = 'C[C@H](O)C[C@H](O)CC[C@H](N)C[C@H](N)C'
    tail = ' |&1:1,4,o1:8,11|'
    drawn = smiles(text + tail)
    and_inverted = text.replace('[C@H](O)', '[C@@H](O)') + tail
    or_inverted = text.replace('[C@H](N)', '[C@@H](N)') + tail
    assert smiles(and_inverted) == drawn
    assert str(smiles(and_inverted)) == str(drawn)
    assert smiles(or_inverted) != drawn


def test_one_group_per_kind_in_every_atom_order():
    """Racemic ephedrine drawn from either end and from either enantiomer writes one string."""
    forward = smiles('C[C@H](NC)[C@@H](O)c1ccccc1 |&1:1,4|')
    backward = smiles('O[C@@H](c1ccccc1)[C@H](C)NC |&1:1,8|')
    enantiomer = smiles('O[C@H](c1ccccc1)[C@@H](C)NC |&1:1,8|')
    assert forward == backward == enantiomer
    assert str(forward) == str(backward) == str(enantiomer)


@mark.parametrize('text, written', [
    ('CC(O)CC |o1:1|', 'C(C)C(O)C'),               # no sign
    ('CCCC |a:1|', 'C(C)CC'),                      # no centre
    ('C[C@H](C)C |&1:1|', 'C(C)(C)C'),             # a sign on a centre that is not stereogenic
])
def test_a_label_that_states_nothing_is_not_written(text, written):
    molecule = smiles(text)
    assert str(molecule) == written
    assert molecule == smiles(written)
