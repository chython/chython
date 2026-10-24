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
"""Bridgeheads whose bridges are too short for in/out isomerism (`docs/stereo.rst`, "Bridgeheads").

The two bridgeheads of such a bicycle are both *out*: one sign implies the other, and a drawing that states the
other sign cannot be built, so neither of its signs is trusted.  Structures are PubChem isomeric SMILES (CID in
the row); a row without a CID is a signing of the named compound chosen out,out.
"""
import re
from pytest import importorskip, mark, skip
from chython.core import LOST, REPAIRED, inchi_library_loaded, read_smiles


def canon(text):
    """`text` read and canonicalized: a coupled set's label states the set after the explicit rewrite."""
    m = read_smiles(text)
    m.canonicalize_stereo_groups()
    return m


#: (name, PubChem CID, isomeric SMILES, bridgehead atom numbers in SMILES order).
COMPOUNDS = [
    # bridged: every bridge 1-3 atoms
    ('(+)-camphor', 159055, 'C[C@@]12CC[C@@H](C1(C)C)CC2=O', (2, 5)),
    ('(-)-borneol', 1201518, 'C[C@]12CC[C@H](C1(C)C)C[C@H]2O', (2, 5)),
    ('(-)-cocaine', 446220, 'CN1[C@H]2CC[C@@H]1[C@H]([C@H](C2)OC(=O)C3=CC=CC=C3)C(=O)OC', (3, 6)),
    ('(1R,2S,4S)-norbornane-2-carboxylic acid', 736412, 'C1C[C@@H]2C[C@H]1C[C@@H]2C(=O)O', (3, 5)),
    ('(1S,2R,4R)-norbornane-2-carboxylic acid', 11829530, 'C1C[C@H]2C[C@@H]1C[C@H]2C(=O)O', (3, 5)),
    ('(1R,2R,4S)-norbornane-2-carboxylic acid', 946560, 'C1C[C@@H]2C[C@H]1C[C@H]2C(=O)O', (3, 5)),
    ('(+)-fenchone', 1201521, 'C[C@]12CC[C@H](C1)C(C2=O)(C)C', (2, 5)),
    ('(+)-alpha-pinene', 82227, 'CC1=CC[C@@H]2C[C@H]1C2(C)C', (5, 7)),
    ('(-)-beta-pinene', 440967, 'CC1([C@H]2CCC(=C)[C@@H]1C2)C', (3, 8)),
    ('tropine', 449293, 'CN1[C@@H]2CC[C@H]1CC(C2)O', (3, 6)),
    ('(+)-camphene', 92221, 'CC1([C@@H]2CC[C@@H](C2)C1=C)C', (3, 6)),
    ('scopolamine', 3000322, 'CN1[C@@H]2CC(C[C@H]1[C@H]3[C@@H]2O3)OC(=O)[C@H](CO)C4=CC=CC=C4', (3, 7)),
    # PubChem signs C2 alone; the bridgeheads are signed out,out
    ('(2R)-bicyclo[2.2.2]octane-2-carboxylic acid', 6928769, 'C1C[C@H]2CC[C@@H]1C[C@H]2C(=O)O', (3, 6)),
    # fused: a zero-atom bridge, the two others at most 2 atoms
    ('tert-butyl (1R,3S)-1-fluoro-3-methyl-2-(2-methylpropyl)bicyclo[1.1.0]butane-2-carboxylate', 71070293,
     'CC(C)CC1([C@]2([C@@]1(C2)F)C)C(=O)OC(C)(C)C', (6, 7)),
    ('(1S,4S)-bicyclo[2.2.0]hexane-2-carboxylic acid', 145913449, 'C1C[C@H]2[C@@H]1CC2C(=O)O', (3, 4)),
    ('3,3-difluoro-1,4-bis(trimethylsilyloxy)bicyclo[2.2.0]hexan-2-one', 134919890,
     'C[Si](C)(C)O[C@]12CC[C@]1(C(C2=O)(F)F)O[Si](C)(C)C', (6, 9)),
    ('[5]-ladderane octanoic acid', 11659465,
     'C1C[C@@H]2[C@H]1[C@@H]3[C@H]2[C@@H]4[C@H]3[C@@H]5[C@H]4[C@H](C5)CCCCCCCC(=O)O', (3, 4)),
    ('tricyclo[4.2.0.0(2,5)]octane-3-carboxylic acid', 165987197, 'C1C[C@@H]2[C@H]1[C@@H]3[C@H]2C(C3)C(=O)O', (3, 4)),
    ('(1S,4R)-1-methylbicyclo[2.1.0]pentane', 101988110, 'C[C@@]12CC[C@@H]1C2', (2, 5)),
    ('(1S,4S)-1-(3,5-difluorophenyl)sulfonylbicyclo[2.1.0]pentane', 125434701,
     'C1C[C@]2([C@@H]1C2)S(=O)(=O)C3=CC(=CC(=C3)F)F', (3, 4)),
    ('3-methylbicyclo[2.1.0]pentane-1-carboxylic acid', 137935629, 'C[C@@H]1C[C@]2([C@@H]1C2)C(=O)O', (4, 5)),
    ('3-hydroxybicyclo[2.1.0]pentane-1-carboxylic acid', 138040499, 'C1[C@@H]2[C@]1(C[C@H]2O)C(=O)O', (2, 3)),
    ('eglumegad (LY354740), a cis-fused bicyclo[3.1.0]hexane', 213056,
     'C1C[C@]([C@H]2[C@@H]1[C@@H]2C(=O)O)(C(=O)O)N', (4, 5)),
]

#: achiral cages signed out,out: equal to the unsigned molecule.
ACHIRAL = [
    ('bicyclo[2.2.2]octane', 'C1C[C@H]2CC[C@@H]1CC2', 'C1CC2CCC1CC2'),
    ('adamantane', 'C1[C@H]2C[C@H]3C[C@@H]1C[C@@H](C2)C3', 'C1C2CC3CC1CC(C2)C3'),
    ('2-adamantanol', 'C1[C@H]2C[C@H]3C[C@@H]1C[C@@H](C2)C3O', 'C1C2CC3CC1CC(C2)C3O'),  # CID 64149
    ('memantine', 'C[C@]12C[C@@H]3C[C@](C1)(C[C@@](C3)(C2)N)C', 'CC12CC3CC(C1)(CC(C3)(C2)N)C'),  # CID 4054
    ('bicyclo[1.1.0]butane-1-carboxylic acid', 'C1[C@H]2[C@@]1(C2)C(=O)O', 'C1C2C1(C2)C(=O)O'),  # CID 20238048
]

_ATOM = re.compile(r'\[[^\]]*\]|Cl|Br|[BCNOPSFI]|[bcnops]|\*')


def _tokens(text):
    return list(_ATOM.finditer(text))


def _replace(text, n, fn):
    """The `n`-th atom token (1-based) replaced by `fn(token)`."""
    t = _tokens(text)[n - 1]
    return text[:t.start()] + fn(t.group()) + text[t.end():]


def drop_sign(text, n):
    return _replace(text, n, lambda t: t.replace('@', ''))


def invert_sign(text, n):
    return _replace(text, n, lambda t: t.replace('@@', '@') if '@@' in t else t.replace('@', '@@'))


def strip_signs(text, atoms=None):
    """`atoms` unsigned, every other sign kept; every sign dropped without `atoms`."""
    if atoms is None:
        return text.replace('@', '')
    for n in atoms:
        text = drop_sign(text, n)
    return text


#: PubChem's InChI of each COMPOUNDS row, by CID.
PUBCHEM_INCHI = {
    159055: 'InChI=1S/C10H16O/c1-9(2)7-4-5-10(9,3)8(11)6-7/h7H,4-6H2,1-3H3/t7-,10+/m1/s1',
    1201518: 'InChI=1S/C10H18O/c1-9(2)7-4-5-10(9,3)8(11)6-7/h7-8,11H,4-6H2,1-3H3/t7-,8+,10+/m0/s1',
    446220: 'InChI=1S/C17H21NO4/c1-18-12-8-9-13(18)15(17(20)21-2)14(10-12)22-16(19)11-6-4-3-5-7-11'
            '/h3-7,12-15H,8-10H2,1-2H3/t12-,13+,14-,15+/m0/s1',
    736412: 'InChI=1S/C8H12O2/c9-8(10)7-4-5-1-2-6(7)3-5/h5-7H,1-4H2,(H,9,10)/t5-,6+,7-/m0/s1',
    11829530: 'InChI=1S/C8H12O2/c9-8(10)7-4-5-1-2-6(7)3-5/h5-7H,1-4H2,(H,9,10)/t5-,6+,7-/m1/s1',
    946560: 'InChI=1S/C8H12O2/c9-8(10)7-4-5-1-2-6(7)3-5/h5-7H,1-4H2,(H,9,10)/t5-,6+,7+/m0/s1',
    1201521: 'InChI=1S/C10H16O/c1-9(2)7-4-5-10(3,6-7)8(9)11/h7H,4-6H2,1-3H3/t7-,10+/m1/s1',
    82227: 'InChI=1S/C10H16/c1-7-4-5-8-6-9(7)10(8,2)3/h4,8-9H,5-6H2,1-3H3/t8-,9-/m1/s1',
    440967: 'InChI=1S/C10H16/c1-7-4-5-8-6-9(7)10(8,2)3/h8-9H,1,4-6H2,2-3H3/t8-,9-/m0/s1',
    # PubChem writes /t6-,7+,8?: with C3 unsigned the meso pair's inversion restates the molecule, so it
    # states nothing and the InChI is CID 8424's
    449293: 'InChI=1S/C8H15NO/c1-9-6-2-3-7(9)5-8(10)4-6/h6-8,10H,2-5H2,1H3',
    92221: 'InChI=1S/C10H16/c1-7-8-4-5-9(6-8)10(7,2)3/h8-9H,1,4-6H2,2-3H3/t8-,9+/m0/s1',
    3000322: 'InChI=1S/C17H21NO4/c1-18-13-7-11(8-14(18)16-15(13)22-16)21-17(20)12(9-19)10-5-3-2-4-6-10'
             '/h2-6,11-16,19H,7-9H2,1H3/t11?,12-,13-,14+,15-,16+/m1/s1',
    6928769: 'InChI=1S/C9H14O2/c10-9(11)8-5-6-1-3-7(8)4-2-6/h6-8H,1-5H2,(H,10,11)/t6?,7?,8-/m1/s1',
    71070293: 'InChI=1S/C14H23FO2/c1-9(2)7-13(10(16)17-11(3,4)5)12(6)8-14(12,13)15/h9H,7-8H2,1-6H3/t12-,13?,14+/m0/s1',
    145913449: 'InChI=1S/C7H10O2/c8-7(9)6-3-4-1-2-5(4)6/h4-6H,1-3H2,(H,8,9)/t4-,5-,6?/m0/s1',
    134919890: 'InChI=1S/C12H22F2O3Si2/c1-18(2,3)16-10-7-8-11(10,17-19(4,5)6)12(13,14)9(10)15/h7-8H2,1-6H3/t10-,11-'
               '/m0/s1',
    11659465: 'InChI=1S/C20H30O2/c21-15(22)7-5-3-1-2-4-6-11-10-14-16(11)20-18-13-9-8-12(13)17(18)19(14)20'
              '/h11-14,16-20H,1-10H2,(H,21,22)/t11-,12-,13+,14-,16+,17+,18-,19-,20+/m0/s1',
    165987197: 'InChI=1S/C9H12O2/c10-9(11)7-3-6-4-1-2-5(4)8(6)7/h4-8H,1-3H2,(H,10,11)/t4-,5+,6+,7?,8-/m0/s1',
    101988110: 'InChI=1S/C6H10/c1-6-3-2-5(6)4-6/h5H,2-4H2,1H3/t5-,6+/m1/s1',
    213056: 'InChI=1S/C8H11NO4/c9-8(7(12)13)2-1-3-4(5(3)8)6(10)11/h3-5H,1-2,9H2,(H,10,11)(H,12,13)/t3-,4-,5-,8-/m0/s1',
    125434701: 'InChI=1S/C11H10F2O2S/c12-8-3-9(13)5-10(4-8)16(14,15)11-2-1-7(11)6-11/h3-5,7H,1-2,6H2/t7-,11-/m0/s1',
    137935629: 'InChI=1S/C7H10O2/c1-4-2-7(6(8)9)3-5(4)7/h4-5H,2-3H2,1H3,(H,8,9)/t4-,5-,7+/m1/s1',
    138040499: 'InChI=1S/C6H8O3/c7-4-2-6(5(8)9)1-3(4)6/h3-4,7H,1-2H2,(H,8,9)/t3-,4+,6-/m0/s1',
}


CONSTRAINED = [c for c in COMPOUNDS]


@mark.parametrize('name, cid, smiles, pair', CONSTRAINED, ids=[c[0] for c in CONSTRAINED])
def test_a_partner_sign_is_implied(name, cid, smiles, pair):
    full = read_smiles(smiles)
    for n in pair:
        half = read_smiles(drop_sign(smiles, n))
        assert half == full and str(half) == str(full)


@mark.parametrize('name, cid, smiles, pair', CONSTRAINED, ids=[c[0] for c in CONSTRAINED])
def test_an_impossible_pair_reads_as_unknown_and_is_logged(name, cid, smiles, pair):
    bad = read_smiles(invert_sign(smiles, pair[1]))
    unknown = read_smiles(strip_signs(smiles, pair))
    assert bad == unknown and str(bad) == str(unknown)
    assert hash(bad) == hash(unknown)
    assert [r.stage for r in bad.log if r.rule == 'stereo:bridgehead-unrealisable'] == ['stereo']
    assert not [r for r in read_smiles(smiles).log if r.rule == 'stereo:bridgehead-unrealisable']


@mark.parametrize('write', ['str', 'hash', 'inchi'])
def test_every_view_that_heals_an_impossible_pair_logs_it_once(write):
    """Camphor (CID 2537) with C4 inverted: str(), hash and `.inchi` each drop both signs and log it."""
    if write == 'inchi' and not inchi_library_loaded():
        skip('libinchi not loaded')
    bad = read_smiles('C[C@@]12CC[C@H](C1(C)C)CC2=O')
    for _ in range(2):
        str(bad) if write == 'str' else hash(bad) if write == 'hash' else bad.inchi
    assert [r.stage for r in bad.log if r.rule == 'stereo:bridgehead-unrealisable'] == ['stereo']
    good = read_smiles('C[C@@]12CC[C@@H](C1(C)C)CC2=O')
    str(good), hash(good)
    assert not [r for r in good.log if r.rule == 'stereo:bridgehead-unrealisable']


@mark.parametrize('name, signed, unsigned', ACHIRAL, ids=[c[0] for c in ACHIRAL])
def test_an_achiral_cage_signed_consistently_is_the_unsigned_molecule(name, signed, unsigned):
    assert read_smiles(signed) == read_smiles(unsigned)
    assert str(read_smiles(signed)) == str(read_smiles(unsigned))


def test_a_bridgehead_pair_carries_no_sign_of_its_own():
    """Bicyclo[2.2.2]octane-2-carboxylic acid has two stereoisomers: C2 alone names one."""
    only_c2 = read_smiles('C1CC2CCC1C[C@H]2C(=O)O')
    for text in ('C1C[C@H]2CC[C@@H]1C[C@H]2C(=O)O', 'C1C[C@@H]2CC[C@H]1C[C@H]2C(=O)O'):
        assert read_smiles(text) == only_c2


def test_in_out_bicycles_stay_free():
    """Bicyclo[4.4.4]tetradecane (CID 14450522) has in,out and out,out isomers: its bridgeheads are free."""
    out_out = read_smiles('C1CC[C@H]2CCCC[C@@H](C1)CCCC2')
    in_out = read_smiles('C1CC[C@H]2CCCC[C@H](C1)CCCC2')
    unsigned = read_smiles('C1CCC2CCCCC(C1)CCCC2')
    assert out_out != in_out and out_out != unsigned and in_out != unsigned
    assert not [r for r in in_out.log if r.rule == 'stereo:bridgehead-unrealisable']


def test_either_bridgehead_signed_writes_one_string():
    """(-)-cocaine, CID 446220, signed at C1 only, at C5 only and at both."""
    smiles = 'CN1[C@H]2CC[C@@H]1[C@H]([C@H](C2)OC(=O)C3=CC=CC=C3)C(=O)OC'
    full = read_smiles(smiles)
    texts = {str(full), str(read_smiles(str(full)))}
    texts.update(str(read_smiles(drop_sign(smiles, n))) for n in (3, 6))
    assert len(texts) == 1
    assert len([t for t in _tokens(texts.pop()) if '@' in t.group()]) == 4     # C1, C2, C3, C5


def _embeds(smiles):
    """Does RDKit's ETKDG place the molecule in 3D?  One seed, 100 iterations: an unrealisable pair
    fails within them, and every row of COMPOUNDS embeds."""
    chem = importorskip('rdkit.Chem')
    all_chem = importorskip('rdkit.Chem.AllChem')
    importorskip('rdkit').RDLogger.DisableLog('rdApp.*')
    params = all_chem.ETKDGv3()
    params.maxIterations = 100
    params.randomSeed = 7
    return all_chem.EmbedMolecule(chem.AddHs(chem.MolFromSmiles(smiles)), params) == 0


@mark.parametrize('name, cid, smiles, pair', CONSTRAINED, ids=[c[0] for c in CONSTRAINED])
def test_the_relation_agrees_with_a_3d_embedding(name, cid, smiles, pair):
    """The harness behind OUT_OUT_XOR, independent of the PubChem signs: the out,out drawing embeds and
    the pair with one sign inverted does not."""
    assert _embeds(smiles)
    assert not _embeds(invert_sign(smiles, pair[1]))


CAMPHOR = 'C[C@@]12CC[C@@H](C1(C)C)CC2=O'


@mark.parametrize('partial', ['|&1:1|', '|&1:4|', '|&1:1,o1:4|', '|o1:1,&1:4|', '|a:1,&1:4|', '|&1:1,a:4|'])
def test_a_coupled_set_takes_the_strongest_label_of_its_members(partial):
    """(+)-camphor, CID 159055: C1 and C4 are one choice, so one label states both, whichever member carries
    it, and AND outranks OR and ABS."""
    whole = canon(f'{CAMPHOR} |&1:1,4|')
    one = canon(f'{CAMPHOR} {partial}')
    assert one == whole and str(one) == str(whole) and hash(one) == hash(whole)


def test_an_or_group_on_one_bridgehead_is_the_or_group_on_both():
    one = canon(f'{CAMPHOR} |o1:1|')
    assert one == canon(f'{CAMPHOR} |o1:1,4|')
    assert one != canon(CAMPHOR) and one != canon(f'{CAMPHOR} |&1:1,4|')


def test_stereo_units_name_the_coupled_partner():
    units = {u['anchor']: u['coupled'] for u in read_smiles(CAMPHOR).stereo_units()}
    assert units[2] == (5,) and units[5] == (2,)
    assert not [a for a, c in units.items() if c and a not in (2, 5)]
    adamantane = read_smiles('C1[C@H]2C[C@H]3C[C@@H]1C[C@@H](C2)C3')
    assert adamantane.unit_of(2)['coupled'] == (4, 6, 8)
    assert read_smiles('C1CC[C@H]2CCCC[C@@H](C1)CCCC2').unit_of(4)['coupled'] == ()


def test_an_impossible_pair_beside_a_centre_it_decides_states_nothing():
    """Tropan-3-ol with C1 and C5 drawn alike: C3 is pseudo-asymmetric only between two bridgeheads."""
    bad = read_smiles('CN1[C@@H]2CC[C@@H]1C[C@@H](C2)O')
    assert bad == read_smiles('CN1C2CCC1CC(C2)O') and str(bad) == str(read_smiles('CN1C2CCC1CC(C2)O'))


@mark.skipif(not inchi_library_loaded(), reason='libinchi not loaded')
@mark.parametrize('name, cid, smiles, pair', CONSTRAINED, ids=[c[0] for c in CONSTRAINED])
def test_equal_bridgehead_drawings_write_one_inchi(name, cid, smiles, pair):
    full = read_smiles(smiles).inchi
    assert full == PUBCHEM_INCHI[cid]
    for n in pair:
        assert read_smiles(drop_sign(smiles, n)).inchi == full
    assert read_smiles(invert_sign(smiles, pair[1])).inchi == read_smiles(strip_signs(smiles, pair)).inchi


@mark.skipif(not inchi_library_loaded(), reason='libinchi not loaded')
@mark.parametrize('name, signed, unsigned', ACHIRAL, ids=[c[0] for c in ACHIRAL])
def test_an_achiral_cage_writes_the_unsigned_inchi(name, signed, unsigned):
    assert read_smiles(signed).inchi == read_smiles(unsigned).inchi
    assert '/t' not in read_smiles(signed).inchi


@mark.skipif(not inchi_library_loaded(), reason='libinchi not loaded')
def test_a_set_inverting_to_the_molecule_writes_unknown():
    """Saxagliptin, CID 11243969: its adamantyl bridgeheads write `?`, as PubChem's record does."""
    signed = read_smiles('C1[C@@H]2C[C@@H]2N([C@@H]1C#N)C(=O)[C@H]([C@]34C[C@@H]5C[C@@H](C3)C[C@](C5)(C4)O)N')
    assert signed.inchi == read_smiles(
        'C1[C@@H]2C[C@@H]2N([C@@H]1C#N)C(=O)[C@H](C34CC5CC(C3)CC(C5)(C4)O)N').inchi == (
        'InChI=1S/C18H25N3O2/c19-8-13-2-12-3-14(12)21(13)16(22)15(20)17-4-10-1-11(5-17)7-18(23,6-10)9-17'
        '/h10-15,23H,1-7,9,20H2/t10?,11?,12-,13+,14+,15-,17?,18?/m1/s1')


def test_a_refused_parity_beside_an_achiral_cage_is_judged():
    """2-Adamantylideneacetic acid: the exocyclic double bond is no cis/trans unit, and its bridgeheads are a
    set whose inversion restates the molecule, so they are no choice to enumerate.  A parity the wedge reader
    sets tentatively on the double bond states nothing."""
    m = read_smiles('OC(=O)C=C1C2CC3CC1CC(C2)C3')
    (unit,) = [u for u in m.stereo_units() if u['kind'] == 1]
    assert not unit['stereogenic']
    with m.edit():
        m.set_parity(unit['anchor'], 1)
    assert not m.stereo_is_live(unit['anchor'])
    assert str(m) == str(read_smiles('OC(=O)C=C1C2CC3CC1CC(C2)C3'))


@mark.parametrize('tail', ['|&1:7|', '|&1:2,5|', '|&1:2,&2:5|', '|o1:7|', '|&1:2,5,7|', '|a:7|'])
def test_a_tropane_labelling_reads_back_to_itself(tail):
    """Tropine, CID 449293: C3 is pseudo-asymmetric between the coupled bridgeheads.  A collection on the
    bridgeheads states the tropine and pseudotropine mixture with C3 open, so the written string reads back
    to the drawing it was written from."""
    m = read_smiles('CN1[C@@H]2CC[C@H]1C[C@@H](C2)O ' + tail)
    back = read_smiles(str(m))
    assert back == m and str(back) == str(m), str(m)
    for _ in range(20):
        spelling = read_smiles(format(m, 'r'))
        assert spelling == m and str(spelling) == str(m), format(m, 'r')


def test_the_tropine_and_pseudotropine_mixture_is_not_the_unsigned_compound():
    mixture = read_smiles('CN1[C@@H]2CC[C@H]1C[C@@H](C2)O |&1:7|')
    assert mixture != read_smiles('CN1C2CCC1CC(C2)O')
    assert mixture == read_smiles('CN1[C@@H]2CC[C@H]1C[C@H](C2)O |&1:7|')     # the drawn phase is no identity


def test_an_or_pair_of_tropane_bridgeheads_with_c3_open_states_nothing():
    """An OR set keeps its sign: with C3 open its worlds are tropine and pseudotropine alone, as unsigned."""
    m = canon('CN1[C@@H]2CC[C@H]1CC(C2)O |o1:2,5|')
    assert m == canon('CN1C2CCC1CC(C2)O')
    assert str(read_smiles(str(m))) == str(m)


def test_a_norbornene_imide_with_a_group_per_centre_reads_back_in_every_atom_order():
    """N-ethyl-5-norbornene-2,3-dicarboximide with each centre in its own AND group: the ring-fusion carbons
    racemic apart make the bridgehead signs state nothing, whichever atom order writes the string."""
    m = read_smiles('O=C1[C@H]2[C@@H]3C=C[C@@H](C3)[C@H]2C(=O)N1CC |&1:2,&2:3,&3:6,&4:8|')
    for _ in range(30):
        spelling = read_smiles(format(m, 'r'))
        assert spelling == m and str(spelling) == str(m), format(m, 'r')
    assert str(read_smiles(str(m))) == str(m)


#: camphor signed in,out at its bridgeheads, which no bicycle builds, and its realisable racemate (CID 159055)
CAMPHOR_IMPOSSIBLE = 'O=C1C[C@@H]2CC[C@]1(C)C2(C)C'
CAMPHOR_RACEMATE = 'C[C@@]12CC[C@@H](C1(C)C)CC2=O |&1:1,4|'


def test_validate_stereo_clears_an_unrealisable_pair_and_logs_it():
    m = read_smiles(CAMPHOR_IMPOSSIBLE)
    bridgeheads = sorted(u['anchor'] for u in m.stereo_units() if u['coupled'])
    assert m.validate_stereo() == bridgeheads
    assert not any(m.parity_of(x) for x in bridgeheads)
    assert [r.severity for r in m.log if r.rule == 'stereo:bridgehead-unrealisable'] == [LOST]
    assert m == read_smiles('O=C1CC2CCC1(C)C2(C)C')
    assert m.validate_stereo() == []


def test_validate_stereo_can_read_an_unrealisable_pair_as_the_racemate():
    m = read_smiles(CAMPHOR_IMPOSSIBLE)
    assert m.validate_stereo(recover_racemate=True) == []
    assert [r.severity for r in m.log if r.rule == 'stereo:bridgehead-racemate'] == [REPAIRED]
    assert m == read_smiles(CAMPHOR_RACEMATE)
    assert str(read_smiles(str(m))) == str(m)


def test_canonicalize_forwards_the_recovery():
    dropped, recovered = read_smiles(CAMPHOR_IMPOSSIBLE), read_smiles(CAMPHOR_IMPOSSIBLE)
    dropped.canonicalize()
    recovered.canonicalize(recover_racemate=True)
    racemate = read_smiles(CAMPHOR_RACEMATE)
    racemate.canonicalize()
    assert dropped != recovered and recovered == racemate


CAMPHOR_SPELLINGS = ['|&1:1|', '|&1:4|', '|&1:1,4|', '|&1:1,a:4|', '|a:1,&1:4|', '|&1:1,&2:4|', '|&2:1,&1:4|']


@mark.parametrize('tail', CAMPHOR_SPELLINGS)
def test_every_spelling_of_a_racemic_cage_writes_one_string_from_every_atom_order(tail):
    """(+)-camphor's racemate, CID 159055: one string, its label on both bridgeheads, whatever the spelling and
    the atom order."""
    whole = canon(f'{CAMPHOR} |&1:1,4|')
    m = canon(f'{CAMPHOR} {tail}')
    assert m == whole and hash(m) == hash(whole) and str(m) == str(whole)
    for _ in range(20):
        spelling = canon(format(m, 'r'))
        assert spelling == whole and str(spelling) == str(whole), format(m, 'r')
    back = read_smiles(str(m))
    assert back == whole and str(back) == str(m)
    assert [len(v) for v in back.stereo_groups().values()] == [2]


@mark.parametrize('kind', ['&', 'o', 'a'])
def test_a_cage_label_is_written_on_every_bridgehead(kind):
    for name, cid, text, pair in COMPOUNDS:
        if text.count('@') != 2:
            continue                                # the cage alone: no other centre to share a group
        a, b = (x - 1 for x in pair)
        label = 'a' if kind == 'a' else f'{kind}1'
        m = canon(f'{text} |{label}:{a}|')
        back = read_smiles(str(m))
        assert back == m and [len(v) for v in back.stereo_groups().values()] in ([2], []), (name, str(m))


#: (1R,2S,4S)-norbornane-2-carboxylic acid, CID 736412: C2 (SMILES atom 7) beside the C1-C4 cage (3, 5)
NORBORNANE_ACID = 'C1C[C@@H]2C[C@H]1C[C@@H]2C(=O)O'


@mark.parametrize('tail', ['|&1:2,4,6|', '|&1:2,6|', '|&1:4,6|', '|&1:2,6,&2:4|', '|&1:6,4,&2:2|'])
def test_two_groups_on_one_cage_are_one_group(tail):
    """The cage inverts with every group it carries, so a group on C2 and another on a bridgehead join: the
    racemate of the acid, its relative configuration stated."""
    whole = canon(f'{NORBORNANE_ACID} |&1:2,4,6|')
    m = canon(f'{NORBORNANE_ACID} {tail}')
    assert m == whole and str(m) == str(whole), str(m)
    for _ in range(10):
        assert str(canon(format(m, 'r'))) == str(whole)


def test_a_cage_apart_from_its_neighbour_centre_is_another_compound():
    """C2 racemic beside a racemic cage of its own group is the mixture of all four isomers' two
    diastereomeric racemates -- not the acid's racemate."""
    apart = canon(f'{NORBORNANE_ACID} |&1:6,&2:2,4|')
    assert apart != canon(f'{NORBORNANE_ACID} |&1:2,4,6|')
    assert apart == canon(f'{NORBORNANE_ACID} |&1:6,&2:2,&3:4|')       # the cage's two groups merge
    assert str(read_smiles(str(apart))) == str(apart)
    assert sorted(len(v) for v in read_smiles(str(apart)).stereo_groups().values()) == [1, 2]


def test_two_independent_centres_keep_both_labels():
    """trans-2-methylcyclohexan-1-ol, CID 23432: two free centres in one AND group state the racemate of one
    diastereomer, and one label alone states a mixture of two, so both are written."""
    pair = read_smiles('C[C@@H]1CCCC[C@H]1O |&1:1,6|')
    assert pair != read_smiles('C[C@@H]1CCCC[C@H]1O |&1:1|')
    assert [len(v) for v in read_smiles(str(pair)).stereo_groups().values()] == [2]


@mark.parametrize('tail', ['|&1:3,6|', '|&1:3,&2:6|'])
def test_validate_stereo_heals_an_unrealisable_set_stated_racemic(tail):
    """Camphor signed in,out with both bridgeheads in AND groups holding nothing else: the drawing states the
    racemate and no relation to anything else, so either realisable signing says the same -- healed without
    `recover_racemate`, logged REPAIRED."""
    m = read_smiles(f'{CAMPHOR_IMPOSSIBLE} {tail}')
    assert m.validate_stereo() == []
    assert m == read_smiles(CAMPHOR_RACEMATE)
    assert [r.severity for r in m.log if r.rule == 'stereo:bridgehead-racemate'] == ['repaired']


def test_an_unrealisable_set_in_an_and_group_beside_its_neighbour_is_not_guessed():
    """5-Norbornene-2,3-dicarboximide with the bridgeheads signed in,out and every centre in one AND group:
    each realisable signing puts the imide exo or endo to the ring fusion, a different compound, so the
    default clears the pair rather than choose."""
    m = read_smiles('O=C1NC(=O)[C@H]2[C@@H]1[C@H]1C=C[C@H]2C1 |&1:4,5,6,9|')
    assert m.validate_stereo()
    assert [r.rule for r in m.log if r.rule.startswith('stereo:bridgehead')] == ['stereo:bridgehead-unrealisable']
