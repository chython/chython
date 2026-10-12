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
"""One string per meaning of an enhanced-stereo drawing, against a brute-force reading of every labelling.

==========  =========================================================================================
Label       Reading
==========  =========================================================================================
cleared     unconfigured: either sign
unlabelled  the sign as drawn, ABS, OR or AND unknown: a component's unlabelled signs read together as each
ABS         the sign as drawn
AND g       the mixture over the group's two phases
OR g        one phase, unknown which; the drawn one is the drawn world
==========  =========================================================================================

A meaning is the set over the unlabelled readings of (drawn worlds, every world), with the components
holding a live label -- unlabelled is unknown, so on a meso pair `|a:1,4|` is not the unlabelled drawing; a
world is one choice of every OR group and cleared unit, as the set of pure compounds over the AND phases.
Before reading, an unlabelled or ABS sign whose inversion leaves the ABS reading unchanged is cleared, to a
fixpoint.
"""
from itertools import product
from random import Random
from pytest import mark
from time import perf_counter
from chython.core import inchi_library_loaded, read_smiles


class Reading:
    """A molecule with its units, and the meaning of any labelling of them."""

    def __init__(self, smiles, elements):
        self.base = read_smiles(smiles)
        self.elements = elements
        units = self.base.stereo_units()
        self.anchor = [e if isinstance(e, int) else
                       next(u['anchor'] for u in units if isinstance(u['owners'], tuple) and
                            sorted(u['owners']) == sorted(e)) for e in elements]
        component = {a: c for c, atoms in enumerate(self.base.connected_components) for a in atoms}
        self.component = [component[e if isinstance(e, int) else e[0]] for e in elements]
        self.k = k = len(elements)
        ids = {}
        self.pure = T = [ids.setdefault(str(self.build([1 + (mask >> i & 1) for i in range(k)], ['u'] * k)), len(ids))
                         for mask in range(1 << k)]
        # units are coupled when inverting one depends on the other, or inverting both restates the molecule
        adjacent = {i: set() for i in range(k)}
        for i in range(k):
            for j in range(k):
                bi, bj = 1 << i, 1 << j
                if i != j and any((T[m] == T[m ^ bi]) != (T[m ^ bj] == T[m ^ bi ^ bj]) or
                                  (T[m ^ bi ^ bj] == T[m] and T[m] != T[m ^ bi]) for m in range(1 << k)):
                    adjacent[i].add(j)
                    adjacent[j].add(i)
        self.clusters = []
        for i in range(k):
            members, stack = {i}, [i]
            while stack:
                for y in adjacent[stack.pop()] - members:
                    members.add(y)
                    stack.append(y)
            self.clusters.append(frozenset(members))

    def build(self, parities, labels):
        """'x' cleared, 'u' unlabelled, 'a' ABS, ('&', g) and ('o', g) groups."""
        m = self.base.copy()
        with m.edit() as e:
            for anchor, p, label in zip(self.anchor, parities, labels):
                e.set_parity(anchor, 0 if label == 'x' else p)
            for element, label in zip(self.elements, labels):
                if label == 'a':
                    e.set_stereo_group(element, 1, 0)
                elif isinstance(label, tuple):
                    e.set_stereo_group(element, 3 if label[0] == '&' else 2, label[1] + 1)
        return m

    def read(self, parities, labels):
        """(drawn worlds, every world), no unlabelled sign left."""
        mask = sum((p == 2) << i for i, p in enumerate(parities))
        groups = {}
        for i, label in enumerate(labels):
            if isinstance(label, tuple):
                groups[label] = groups.get(label, 0) | 1 << i
        ands = [m for g, m in groups.items() if g[0] == '&']
        ors = [m for g, m in groups.items() if g[0] == 'o']
        cleared = [i for i, label in enumerate(labels) if label == 'x']
        drawn, every = set(), set()
        for phases in product((0, 1), repeat=len(ors)):
            flip = 0
            for p, m in zip(phases, ors):
                flip ^= m * p
            for fill in product((0, 1), repeat=len(cleared)):
                m0 = mask ^ flip
                for i, b in zip(cleared, fill):
                    m0 = m0 & ~(1 << i) | b << i
                world = set()
                for more in product((0, 1), repeat=len(ands)):
                    m1 = m0
                    for p, m in zip(more, ands):
                        m1 ^= m * p
                    world.add(self.pure[m1])
                world = frozenset(world)
                every.add(world)
                if not any(phases):
                    drawn.add(world)
        return frozenset(drawn), frozenset(every)

    def as_abs(self, parities, labels):
        return self.read(parities, ['a' if label == 'u' else label for label in labels])

    def meaning(self, parities, labels):
        labels = list(labels)
        while True:
            ref = self.as_abs(parities, labels)
            inert = [i for i, label in enumerate(labels) if label in ('u', 'a') and
                     self.as_abs([3 - p if j == i else p for j, p in enumerate(parities)], labels) == ref]
            if not inert:
                break
            for i in inert:
                labels[i] = 'x'
        components = sorted({self.component[i] for i, label in enumerate(labels) if label == 'u'})
        out = set()
        for kinds in product('a&o', repeat=len(components)):
            read = list(labels)
            for c, kind in zip(components, kinds):
                for i, label in enumerate(labels):
                    if label == 'u' and self.component[i] == c:
                        read[i] = 'a' if kind == 'a' else (kind, 100 + c)
            out.add(self.read(parities, read))
        # unlabelled is unknown: a cluster whose labels, all cleared, move the reading keeps a label.  Symmetric
        # choices are one choice, so the set is named by the string of the unsigned molecule with every owner
        # atom of the stated clusters marked
        ref = self.as_abs(parities, labels)
        labelled = {i for i, label in enumerate(labels) if label not in ('x', 'u')}
        units = {j for c in {self.clusters[i] for i in labelled}
                 if self.as_abs(parities, ['x' if i in labelled and i in c else label
                                           for i, label in enumerate(labels)]) != ref for j in c}
        stated = ''
        if units:
            marked = self.base.copy()
            with marked.edit() as e:
                for i in units:
                    for atom in (self.elements[i] if isinstance(self.elements[i], tuple) else (self.elements[i],)):
                        e.set_isotope(atom, 13)
            stated = str(marked)
        return frozenset(out), stated


def labellings(k):
    """Per unit 'x', 'u', 'a' or a group, groups numbered by first use, each of either kind."""
    def partial(i, cur, blocks):
        if i == k:
            yield cur
            return
        for option in ['x', 'u', 'a', *range(blocks + 1)]:
            yield from partial(i + 1, cur + [option], blocks + (option == blocks))

    for shape in partial(0, [], 0):
        blocks = sorted({x for x in shape if isinstance(x, int)})
        for kinds in product('&o', repeat=len(blocks)):
            yield [x if not isinstance(x, int) else (kinds[x], x) for x in shape]


CASES = [
    ('pseudoephedrine', 'CC(NC)C(O)c1ccccc1', [2, 5]),
    ('butan-2-ol + pseudoephedrine', 'CC(O)CC.CC(NC)C(O)c1ccccc1', [2, 7, 10]),
    ('pentane-2,4-diol', 'CC(O)CC(O)C', [2, 5]),
    ('pentane-2,3,4-triol', 'CC(O)C(O)C(O)C', [2, 4, 6]),
    ('1,3-dimethylcyclobutane', 'CC1CC(C)C1', [2, 4]),
    ('1,4-dimethylcyclohexane', 'CC1CCC(C)CC1', [2, 5]),
    ('Fecht acid', 'OC(=O)C1CC2(C1)CC(C2)C(=O)O', [4, 6, 9]),
    ('2,6-dimethylspiro[3.3]heptane', 'CC1CC2(C1)CC(C)C2', [2, 4, 7]),
    ('2,7-dimethylspiro[4.4]nonane', 'CC1CCC2(C1)CCC(C)C2', [2, 5, 9]),
    ('spiro[4.4]nonane-1,6-dione', 'O=C1CCCC12CCCC2=O', [6]),
    ('4-methylcyclohexylideneacetic acid', 'CC1CCC(=CC(=O)O)CC1', [2, (5, 6)]),
    ("4,4'-dimethylbi(cyclohexylidene)", 'CC1CCC(=C2CCC(C)CC2)CC1', [2, (5, 6), 9]),
    ('muconic acid', 'OC(=O)C=CC=CC(=O)O', [(4, 5), (6, 7)]),
    ('penta-2,3-dien-1-ol + centre', 'CC(O)C=C=CC', [2, (4, 6)]),
    ('[6]helicene', 'c1ccc2c(c1)ccc1ccc3ccc4ccc5ccccc5c4c3c12', [(3, 22)]),
    ('2-(1-hydroxyethyl)[6]helicene', 'CC(O)c1ccc2ccc3ccc4ccc5ccc6ccccc6c5c4c3c2c1', [2, (23, 29)]),
    ('4-carboxy[2.2]paracyclophane', 'OC(=O)c1cc2ccc1CCc1ccc(cc1)CC2', [4]),
    ('4,12-dibromo[2.2]paracyclophane', 'Brc1cc2ccc1CCc1ccc(c(Br)c1)CC2', [14, 2]),
    ('(E)-cyclooctene', 'C1CCC/C=C/CC1', [(5, 6), 6]),
    ('PPFA', 'CN(C)C(C)c%30~%11c(P(c1ccccc1)c1ccccc1)~%12c~%13c~%14c~%15%30.'
             'c%31~%21c~%22c~%23c~%24c~%25%31.[Fe]%11%12%13%14%15%21%22%23%24%25', [4, 6, 7]),
]

INCHI_CASES = ('pseudoephedrine', 'pentane-2,3,4-triol', 'muconic acid')


@mark.parametrize('smiles,elements', [c[1:] for c in CASES], ids=[c[0] for c in CASES])
def test_every_meaning_writes_one_string_and_reads_back(smiles, elements):
    reading = Reading(smiles, elements)
    by_meaning, by_string = {}, {}
    for labels in labellings(reading.k):
        configured = [i for i, label in enumerate(labels) if label != 'x']
        for signs in product((1, 2), repeat=len(configured)):
            parities = [1] * reading.k
            for i, p in zip(configured, signs):
                parities[i] = p
            m = reading.build(parities, labels)
            text = str(m)
            meaning = reading.meaning(parities, labels)
            by_meaning.setdefault(meaning, set()).add(text)
            by_string.setdefault(text, set()).add(meaning)
            back = read_smiles(text)
            assert back == m and str(back) == text, (labels, parities, text)
    assert not [w for w in by_meaning.values() if len(w) > 1]
    assert not [s for s, v in by_string.items() if len(v) > 1]


FECHT = 'OC(=O)[C@@H]1C[C@]2(C[C@@H](C(=O)O)C2)C1'


@mark.parametrize('tail,written', [
    ('', FECHT),                                            # unlabelled: kept as drawn
    ('a:3,5,7', f'{FECHT} |a:5|'),                          # pure: ABS on the spiro atom, the ends unlabelled
    ('o1:3,5,7', f'{FECHT} |o1:5|'),                        # isolated, unknown which: OR on the spiro atom
    ('&1:3,5,7', 'OC(=O)C1C[C@]2(CC(C(=O)O)C2)C1 |&1:5|'),  # racemate: the spiro sign alone, the ends open
])
def test_the_fecht_acid_label_sits_on_the_spiro_atom(tail, written):
    m = read_smiles(f'{FECHT} |{tail}|' if tail else FECHT)
    assert str(m) == written
    assert read_smiles(written) == m


def test_a_spiro_sign_with_open_ends_states_nothing_but_in_and():
    # every fill of the open ends gives both enantiomers under an AND phase, and one isomer otherwise
    for tail in ('', ' |a:5|', ' |o1:5|'):
        m = read_smiles('OC(=O)C1C[C@]2(CC(C(=O)O)C2)C1' + tail)
        assert str(m) == 'OC(=O)C1CC2(CC(C(=O)O)C2)C1', tail
    assert read_smiles('OC(=O)C1C[C@]2(CC(C(=O)O)C2)C1 |&1:5|') == read_smiles(f'{FECHT} |&1:3,5,7|')


def test_the_bicyclohexylidene_racemate_is_the_double_bond_in_and():
    m = read_smiles('C/1(=C/2\\CC[C@@H](CC2)C)\\CC[C@@H](CC1)C |&1:0,5,12|')
    assert str(m) == 'C/1(=C\\2/CCC(CC2)C)\\CCC(CC1)C |&1:0|'


def test_an_and_group_over_one_axis_states_the_pure_compound():
    # inverting both units of 4-methylcyclohexylideneacetic acid restates it: the pair in AND is not the
    # racemate but the pure compound, as ABS states it, and a stated label is not the unlabelled drawing;
    # one unit in AND is the racemate
    plain = read_smiles('C(/C(O)=O)=C/1\\CC[C@H](C)CC1')
    pair = read_smiles('C(/C(O)=O)=C/1\\CC[C@H](C)CC1 |&1:0,7|')
    assert pair == read_smiles('C(/C(O)=O)=C/1\\CC[C@H](C)CC1 |a:0,7|') and pair != plain
    racemate = read_smiles('C(/C(O)=O)=C/1\\CC[C@H](C)CC1 |&1:7|')
    assert racemate != plain
    assert racemate == read_smiles('C(/C(O)=O)=C/1\\CC[C@H](C)CC1 |&1:0|')


@mark.parametrize('kind', ['&', 'o'])
def test_one_group_over_many_identical_components_is_exact_and_polynomial(kind):
    # the copy swaps preserve the group, so the canonical search prunes the n! labellings to one
    for n in (6, 12, 30):
        smiles = '.'.join(['C[C@H](O)CC'] * n)
        one = read_smiles(f"{smiles} |{kind}1:{','.join(str(5 * i + 1) for i in range(n))}|")
        each = read_smiles(f"{smiles} |{','.join(f'{kind}{i + 1}:{5 * i + 1}' for i in range(n))}|")
        start = perf_counter()
        text = str(one)
        back = read_smiles(text)
        assert back == one and str(back) == text and hash(back) == hash(one)
        assert one != each and str(one) != str(each)
        assert read_smiles(str(each)) == each
        assert perf_counter() - start < 2, n


def test_blocks_joined_by_no_group_are_searched_apart_and_exactly():
    # three copies of hexane-2,3,4,5-tetrol, four units each: twelve units in three blocks of four
    tetrol = 'CC(O)C(O)C(O)C(O)C'
    one = Reading(tetrol, [2, 4, 6, 8])
    spellings = {}
    for labels in labellings(4):
        configured = [i for i, label in enumerate(labels) if label != 'x']
        for signs in product((1, 2), repeat=len(configured)):
            parities = [1] * 4
            for i, p in zip(configured, signs):
                parities[i] = p
            spellings.setdefault(one.meaning(parities, labels), []).append((parities, labels))
    several = [m for m, s in spellings.items() if len(s) > 1]
    mixture = Reading('.'.join([tetrol] * 3), [2, 4, 6, 8, 12, 14, 16, 18, 22, 24, 26, 28])

    def placed(picks):
        parities, labels = [], []
        for c, (p, ls) in enumerate(picks):
            parities += p
            labels += [(x[0], x[1] + 4 * c) if isinstance(x, tuple) else x for x in ls]
        return parities, labels

    rng = Random(7)
    by_string = {}
    for _ in range(60):
        meanings = [rng.choice(several) for _ in range(3)]
        a = placed([rng.choice(spellings[m]) for m in meanings])
        b = placed([rng.choice(spellings[m]) for m in rng.sample(meanings, 3)])
        # the stated components are indices, which the permutation moves: the copies' own meanings are the
        # invariant, as a multiset
        meaning = mixture.meaning(*a)[0], tuple(sorted(map(repr, meanings)))
        assert mixture.meaning(*b)[0] == meaning[0]
        ma, mb = mixture.build(*a), mixture.build(*b)
        assert str(ma) == str(mb) and ma == mb and hash(ma) == hash(mb), (a, b)
        assert not ma.identity_approximate
        by_string.setdefault(str(ma), set()).add(meaning)
    assert not [s for s, v in by_string.items() if len(v) > 1]


def test_a_block_over_the_unit_limit_is_flagged():
    # nonane-2,3,4,5,6,7,8-heptol: one component, seven units joined by groups
    m = read_smiles('O[C@H](C)[C@H]([C@H]([C@H](O)[C@@H](O)[C@H]([C@H](C)O)O)O)O |&1:5,&2:3,&3:9,&4:1,&5:10,&6:4,7|')
    assert m.identity_approximate
    hash(m)
    assert [r.stage for r in m.log if r.rule == 'stereo:identity-approximate'] == ['stereo']
    assert not read_smiles('CC(O)C(O)C(O)C |&1:1,3,5|').identity_approximate


@mark.skipif(not inchi_library_loaded(), reason='libinchi not loaded')
@mark.parametrize('smiles,elements', [c[1:] for c in CASES if c[0] in INCHI_CASES],
                  ids=[c[0] for c in CASES if c[0] in INCHI_CASES])
def test_equal_molecules_write_one_inchi(smiles, elements):
    # an AND or OR member's sign is not written; a labelling with none writes its stored signs
    reading = Reading(smiles, elements)
    by_string = {}
    for labels in labellings(reading.k):
        configured = [i for i, label in enumerate(labels) if label != 'x']
        for signs in product((1, 2), repeat=len(configured)):
            parities = [1] * reading.k
            for i, p in zip(configured, signs):
                parities[i] = p
            m = reading.build(parities, labels)
            text = m.inchi
            assert m.inchi == text          # written twice, logged once
            by_string.setdefault(str(m), set()).add(text)
            stated = [r for r in m.log if r.rule == 'inchi:enhanced-stereo-not-stated']
            if not any(isinstance(label, tuple) for label in labels):
                plain = reading.build(parities, ['u' if label == 'a' else label for label in labels])
                assert text == plain.inchi and not stated, (labels, parities)
            assert len(stated) <= 1
    assert not [s for s, v in by_string.items() if len(v) > 1]


def test_past_the_atom_limit_the_labelling_is_not_searched():
    """Each combination costs a whole-molecule identity, so beside a 1000-carbon chain the ribo-triol's two
    labellings of one mixture are two identities, flagged and logged, and each reads back to itself."""
    left, right = 'C[C@H](O)[C@H](O)[C@H](O)C |&1:1,5,&2:3|', 'C[C@H](O)[C@H](O)[C@H](O)C |a:1,5,&1:3|'
    assert read_smiles(left) == read_smiles(right)
    big = [read_smiles(t.replace(' ', '.' + 'C' * 1000 + ' ')) for t in (left, right)]
    assert big[0] != big[1]
    for m in big:
        assert m.identity_approximate
        assert 'stereo:identity-approximate' in [r.rule for r in m.log]
        assert read_smiles(str(m)) == m
