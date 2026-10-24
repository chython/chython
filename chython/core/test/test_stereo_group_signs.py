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
label states nothing and the molecule is the unlabelled one. Two labellings stating one mixture are one
compound and write one string after `canonicalize_stereo_groups()`; read, each is its stored labelling.
"""
import re
from itertools import product

from pytest import mark

from chython.core._core import read_smiles


def smiles(text):
    """`text` read and canonicalized: what these tests compare is the meaning, by the explicit rewrite."""
    m = read_smiles(text)
    m.canonicalize_stereo_groups()
    return m


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
    ('pseudoephedrine', 'C[C@H](NC)[C@@H](O)c1ccccc1', '1,4', False),
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
    """Pseudoephedrine becomes ephedrine, trans the cis diol, meso the racemic one: under every kind."""
    drawn = smiles(f'{text} |{kind}:{members}|')
    assert drawn != smiles(f'{one_inverted(text)} |{kind}:{members}|'), name


@mark.parametrize('name, text, members, meso', COMPOUNDS, ids=[c[0] for c in COMPOUNDS])
def test_the_three_kinds_and_no_label_are_four_identities(name, text, members, meso):
    """On a meso form the three kinds state one thing, the drawn diastereomer, as ABS does; unlabelled is
    unknown, so it is a fourth identity.  On a chiral form all four differ."""
    forms = [smiles(text)] + [smiles(f'{text} |{kind}:{members}|') for kind in ('&1', 'o1', 'a')]
    for i, a in enumerate(forms):
        for b in forms[i + 1:]:
            same = meso and i > 0
            assert (a == b) is same, (name, str(a), str(b))
            assert (str(a) == str(b)) is same, (name, str(a), str(b))


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
    """Racemic pseudoephedrine drawn from either end and from either enantiomer writes one string."""
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
    assert str(read_smiles(text)) != written, 'read, the stored label is written'
    molecule = smiles(text)
    assert str(molecule) == written
    assert molecule == smiles(written)


#: pentane-2,3,4-triol, centres at 1, 3 and 5
TRIOL = 'C[C%sH](O)[C%sH](O)[C%sH](O)C'
FLIP = {'@': '@@', '@@': '@'}


def triol_orders(signs, tail):
    """The triol with `signs` at C2, C3, C4 written forward and from the other end."""
    swap = {'1': '5', '3': '3', '5': '1'}
    backward = re.sub(r'(?<=[:,])\d+', lambda x: swap[x.group()], tail)
    return f'{TRIOL % signs} {tail}', f'{TRIOL % tuple(FLIP[x] for x in reversed(signs))} {backward}'


@mark.parametrize('signs', [('@', '@', '@'), ('@', '@@', '@')], ids=['meso-1', 'meso-2'])
def test_four_labellings_of_one_meso_triol_mixture_write_one_string(signs):
    """C2 and C4 as one AND pair give both meso forms: C3 in AND, in OR or in ABS, or the pair in ABS and C3
    in AND, all state that mixture and write its least labelling, the pair alone -- its other phase turns
    3r into 3s, so C3's sign states nothing."""
    tails = ['|&1:1,5,&2:3|', '|&1:1,5,o1:3|', '|a:3,&1:1,5|', '|a:1,5,&1:3|']
    forms = [smiles(t) for tail in tails for t in triol_orders(signs, tail)]
    for m in forms:
        assert m == forms[0] and hash(m) == hash(forms[0]) and str(m) == str(forms[0]), str(m)
    assert re.search(r' \|&1:\d+,\d+\|$', str(forms[0])), str(forms[0])     # the pair alone, one AND group


def test_an_or_pair_is_not_two_or_centres():
    """arabino-Pentane-2,3,4-triol, C2 absolute: C3 and C4 as one OR pair keep their relative
    configuration, as two OR groups they do not."""
    pair = [smiles(t) for t in triol_orders(('@@', '@', '@'), '|a:1,o1:3,5|')]
    apart = [smiles(t) for t in triol_orders(('@@', '@', '@'), '|a:1,o1:3,o2:5|')]
    assert pair[0] == pair[1] and str(pair[0]) == str(pair[1])
    assert apart[0] == apart[1] and str(apart[0]) == str(apart[1])
    assert pair[0] != apart[0] and str(pair[0]) != str(apart[0])


# --- exhaustive: every labelling of every stereoisomer ------------------------------------------------------

#: name: (template with one %s per centre, the centres' slots, meanings allowed more than one string)
EXHAUSTIVE = {
    'pentane-2,4-diol': ('C[C%sH](O)C[C%sH](O)C', (1, 4), 0),
    'ephedrine and pseudoephedrine': ('C[C%sH](NC)[C%sH](O)c1ccccc1', (1, 4), 0),
    'pentane-2,3,4-triol': ('C[C%sH](O)[C%sH](O)[C%sH](O)C', (1, 3, 5), 6),
    'hexane-2,3,4,5-tetrol': ('C[C%sH](O)[C%sH](O)[C%sH](O)[C%sH](O)C', (1, 3, 5, 7), 63),
}


def labellings(n):
    """Each centre in ABS, an AND group or an OR group; group ids numbered by first appearance per kind."""
    out = set()
    for raw in product(['a', '&1', '&2', '&3', '&4', 'o1', 'o2', 'o3', 'o4'], repeat=n):
        names, lab = {}, []
        for x in raw:
            if x != 'a' and x not in names:
                names[x] = x[0] + str(1 + sum(1 for y in names if y[0] == x[0]))
            lab.append(names.get(x, x))
        out.add(tuple(lab))
    return sorted(out)


def meaning(species, signs, lab):
    """(the drawn world, every world).  A world fixes each OR group's phase and is the set of species over
    the AND phases; the drawn world has every OR group as drawn, since an OR group keeps its sign."""
    ands = sorted({x for x in lab if x[0] == '&'})
    ors = sorted({x for x in lab if x[0] == 'o'})
    worlds, drawn = set(), None
    for op in product([0, 1], repeat=len(ors)):
        world = set()
        for ap in product([0, 1], repeat=len(ands)):
            flip = {g for g, p in zip(ors, op) if p} | {g for g, p in zip(ands, ap) if p}
            world.add(species(tuple(FLIP[x] if g in flip else x for x, g in zip(signs, lab))))
        worlds.add(frozenset(world))
        if not any(op):
            drawn = frozenset(world)
    return drawn, frozenset(worlds)


@mark.parametrize('name', list(EXHAUSTIVE))
def test_no_string_states_two_meanings(name):
    """One string never names two mixtures and reads back to itself; at most `split` mixtures keep more than
    one string, the labellings no relabelling move reaches."""
    template, slots, split = EXHAUSTIVE[name]
    cache = {}

    def species(signs):
        if signs not in cache:
            cache[signs] = str(smiles(template % signs))
        return cache[signs]

    drawings = {}
    for signs in product(['@', '@@'], repeat=len(slots)):
        drawings.setdefault(species(signs), signs)
    by_string, by_meaning = {}, {}
    for signs in drawings.values():
        for lab in labellings(len(slots)):
            groups = {}
            for slot, g in zip(slots, lab):
                groups.setdefault(g, []).append(str(slot))
            text = f"{template % signs} |{','.join(f'{g}:' + ','.join(v) for g, v in groups.items())}|"
            written = str(smiles(text))
            assert str(read_smiles(written)) == written, text
            key = meaning(species, signs, lab)
            by_string.setdefault(written, set()).add(key)
            by_meaning.setdefault(key, set()).add(written)
    assert [w for w, v in by_string.items() if len(v) > 1] == []
    assert sum(1 for v in by_meaning.values() if len(v) > 1) <= split
