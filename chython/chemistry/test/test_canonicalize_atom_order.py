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
"""`canonicalize()` returns one representation per compound, whatever the atom order and Kekule form.

Each drawing is a random atom order of the compound, written aromatic or re-kekulised in that order --
the second gives a different Kekule form per order wherever the compound has several.  The fullerenes
and porphyrinoids are where the stored ring basis, the kekuliser's search order and a tie between
mobile-hydrogen placements each followed the atom order; the benzenoids are controls.
"""
import random

from pytest import mark

from .. import canonicalize   # the import injects the method
from ...core import read_smiles as smiles


C60 = ('c12c3c4c5c1c1c6c7c2c2c8c3c3c9c4c4c%10c5c5c1c1c6c6c%11c7c2c2c7c8c3c3c8c9c4c4c9c%10c5c5c1c1c6c6c%11'
       'c2c2c7c3c3c8c4c4c9c5c1c1c6c2c3c41')
#: PCBM, [6,6]-phenyl-C61-butyric acid methyl ester, in one Kekule form.
PCBM = ('C1=2C3=C4C=5C=6C7=C8C=9C=5C5=C3C=3C=%10C%11%12C(C%13=CC=CC=C%13)(CCCC(OC)=O)C%115C=9C=5C=9C=%11'
        'C%13=C%14C=%15C=%16C%17=C%18C=%15C=%11C=5C%12=C%18C=%10C=5C(C1=3)=C1C(C%17=5)=C3C5=C%10C=%11C=%12'
        'C(C=6C4=C4C=2C1=C5C4=%12)=C1C=%11C(=C%13C(=C71)C=98)C%14=C%10C=%163')
#: phthalocyanine, CAS 574-93-6, and 5,10,15,20-tetraphenylporphyrin, CAS 917-23-7.
PHTHALOCYANINE = 'c1ccc2c(c1)c1nc2nc2[nH]c(nc3nc(nc4[nH]c(n1)c1ccccc41)c1ccccc31)c1ccccc21'
TPP = 'c1ccc(cc1)-c1c2ccc(n2)c(-c2ccccc2)c2ccc([nH]2)c(-c2ccccc2)c2ccc(n2)c(-c2ccccc2)c2ccc1[nH]2'
PORPHINE = 'c1cc2cc3ccc(cc4ccc(cc5ccc(cc1n2)[nH]5)n4)[nH]3'

CASES = {
    'C60': C60,
    'PCBM': PCBM,
    'phthalocyanine': PHTHALOCYANINE,
    'tetraphenylporphyrin': TPP,
    'porphine': PORPHINE,
    'naphthalene': 'c1ccc2ccccc2c1',
    'pyrene': 'c1cc2ccc3cccc4ccc(c1)c2c34',
    'triphenylene': 'c1ccc2c(c1)c1ccccc1c1ccccc21',
    'biphenylene': 'c1ccc2c(c1)-c1ccccc1-2',
}


def drawings(text, count):
    base = smiles(text)
    rng = random.getstate()
    random.seed(20260928)
    try:
        for i in range(count):
            order = smiles(format(base, 'r'))
            if i % 2:
                order.kekule()
                order = smiles(format(order, 'r'))
            yield order
    finally:
        random.setstate(rng)


@mark.parametrize('name', CASES)
def test_twenty_atom_orders_give_one_representation(name):
    forms = set()
    for molecule in drawings(CASES[name], 20):
        molecule.canonicalize()
        forms.add((str(molecule), molecule.canonical_bytes))
    assert len(forms) == 1, sorted(forms)


@mark.parametrize('name', ['C60', 'PCBM', 'phthalocyanine', 'tetraphenylporphyrin', 'porphine'])
def test_the_canonical_string_reads_back_to_itself(name):
    for molecule in drawings(CASES[name], 4):
        molecule.canonicalize()
        text = str(molecule)
        again = smiles(text)
        again.canonicalize()
        assert again.canonical_bytes == molecule.canonical_bytes
        assert str(again) == text


def test_every_face_of_c60_is_aromatic_from_every_order():
    for molecule in drawings(C60, 20):
        molecule.canonicalize()
        assert all(bond.order == 4 for bond in molecule.bonds())


@mark.parametrize('text', ['CC1=CC2=CC=C3C=CC=C4C=CC(=C1)[C@]2(C)[C@@]34C',
                           'CC=1C=C2C=CC3=CC=CC4=CC=C(C=1)[C@]2(C)[C@@]34C'])
def test_a_perimeter_thiele_does_not_take_keeps_its_bonds_and_its_configuration(text):
    """2,15,16-trimethyl-15,16-dihydropyrene: the 14-atom perimeter has two matchings and no aromatic form.

    Rewriting one matching onto the other on the same atoms keeps the stated parities while moving the
    bonds they are read against, which gives another stereoisomer.  So nothing is rewritten.
    """
    molecule = smiles(text)
    before = molecule.copy()
    molecule.canonicalize()
    assert molecule == before, 'canonicalize() changed the compound'
    assert not any(r.rule == 'kekule-form:canonical' for r in molecule.log)


def test_an_amidine_geometry_does_not_depend_on_the_atom_order():
    """A C=N is labile, so neither N-benzyl-N'-ethyl-N''-methylguanidine drawn E nor N,N-dimethylacetamidine
    N-methylimine keeps a geometry, and each lands on one string from every atom order."""
    for text in ('CN/C(/NCc1ccccc1)=N/CC', 'CN(C)/C(/C)=N/C'):
        m = smiles(text)
        out = set()
        for _ in range(20):
            r = smiles(format(m, 'r'))
            canonicalize(r)
            out.add(str(r))
        assert len(out) == 1, out
        assert '/' not in out.pop()


def test_folding_an_explicit_hydrogen_keeps_the_configuration_from_every_atom_order():
    """Bromochlorodeuteromethane drawn with an explicit [H]: folding it into the count re-expresses the sign
    against the implicit direction, so each enantiomer lands on one molecule from every atom order."""
    from .. import implicify_hydrogens
    for text, folded in (('[H][C@]([2H])(Cl)Br', '[C@H]([2H])(Cl)Br'), ('[H][C@@]([2H])(Cl)Br', '[C@@H]([2H])(Cl)Br')):
        m = smiles(text)
        for _ in range(20):
            r = smiles(format(m, 'r'))
            implicify_hydrogens(r)
            assert r == smiles(folded), (text, str(r))


def test_an_and_member_perception_alone_does_not_mark_survives_canonicalize():
    """3-Phenoxytropane with C3 in an AND group and its bridgeheads open is the exo and endo mixture: C3 is not
    stereogenic by itself, but its label states something, so the string `str()` writes canonicalizes as the
    drawing does."""
    m = smiles('CN1[C@H]2CC[C@@H]1C[C@@H](C2)Oc1ccccc1 |&1:7|')
    direct = m.copy()
    canonicalize(direct)
    written = smiles(str(m))
    canonicalize(written)
    assert str(direct) == str(written) and '&1' in str(direct)


@mark.parametrize('drawings', [('CN1CC(=O)N=C1N', 'CN1CC(=O)NC1=N'),          # creatinine, CID 588
                               ('NC1=NCCN1', 'N=C1NCCN1')])                # 2-aminoimidazoline, CID 66213
def test_a_cyclic_guanidine_or_amidine_is_one_molecule_from_every_tautomer_and_order(drawings):
    """A ring amidine with no C=C is placed like an acyclic one: an N-H beside a C=O keeps its hydrogen, and
    the canonical ranks decide the rest, so creatinine lands on its 2-imino form."""
    out = set()
    for text in drawings:
        m = smiles(text)
        for _ in range(10):
            r = smiles(format(m, 'r'))
            canonicalize(r)
            out.add(str(r))
    assert len(out) == 1, out
    if drawings[0].startswith('CN1CC(=O)'):
        lactam = smiles(out.pop())
        def carbonyl(x):
            return lactam.element_of(x) == 6 and any(lactam.order_of(x, y) == 2 and lactam.element_of(y) == 8
                                                     for y in lactam.neighbors_of(x))
        assert any(lactam.element_of(n) == 7 and lactam.implicit_h_of(n) and
                   any(carbonyl(x) for x in lactam.neighbors_of(n)) for n in lactam.atom_numbers)


def test_a_sulfonimidamide_is_one_molecule_from_either_tautomer():
    """N-methylbenzenesulfonimidamide: the S=N moves with the hydrogen as an amidine C=N does, the S=O stays,
    and a configured sulfur keeps its sign, since its four neighbours are the same atoms either way."""
    for pair in (('c1ccccc1S(=O)(=N)NC', 'c1ccccc1S(=O)(N)=NC'), ('c1ccccc1[S@](=O)(=N)NC', 'c1ccccc1[S@](=O)(N)=NC')):
        out = set()
        for text in pair:
            m = smiles(text)
            for _ in range(8):
                r = smiles(format(m, 'r'))
                canonicalize(r)
                out.add(str(r))
        assert len(out) == 1, out


@mark.parametrize('text', ['OC(=O)CN(CC(O)=O)CCN(CC(O)=O)CC(O)=O.[Na].[Na]',        # disodium EDTA, CAS 139-33-3
                           '[NH4+].[Na+].[Na+].OC(CC([O-])=O)(CC([O-])=O)C([O-])=O'])   # sodium ammonium citrate
def test_a_salt_charges_the_same_sites_from_every_atom_order(text):
    """Two of EDTA's four equivalent carboxyls take the sodium charges, and one of citrate's three carboxylates
    takes the ammonium proton -- which ones is the canonical order's choice and never the atom numbering's."""
    m = smiles(text)
    out = set()
    for _ in range(12):
        r = smiles(format(m, 'r'))
        canonicalize(r)
        out.add(str(r))
    assert len(out) == 1, out


#: (name, PubChem CID, tautomer drawings, the canonical string or None when only convergence is asserted)
AMIDINES = [
    ('creatinine', 588, ['CN1CC(=O)N=C1N', 'CN1CC(=O)NC1=N'], 'C1(NC(=N)N(C)C1)=O'),
    ('2-aminoimidazoline', 66213, ['NC1=NCCN1', 'N=C1NCCN1'], None),
    ('2-amino-1,4,5,6-tetrahydropyrimidine', 69625, ['NC1=NCCCN1', 'N=C1NCCCN1'], None),
    ('glycocyamidine', 10979, ['NC1=NC(=O)CN1', 'N=C1NC(=O)CN1', 'NC1=NCC(=O)N1'], 'N1C(=N)NCC1=O'),
    ('isocytosine', 66950, ['Nc1nccc(=O)[nH]1', 'NC1=NC(=O)C=CN1', 'N=C1NC=CC(=O)N1', 'Nc1nc(O)ccn1'],
     'O=C1NC(=NC=C1)N'),
    ('cytosine', 597, ['Nc1cc[nH]c(=O)n1', 'N=c1cc[nH]c(=O)[nH]1', 'NC1=NC(=O)NC=C1'], 'O=C1NC=CC(N)=N1'),
    ('2-aminopyridine', 10439, ['Nc1ccccn1', 'N=C1C=CC=CN1'], 'c1cccnc1N'),
    ('2-aminobenzimidazole', 13624, ['Nc1nc2ccccc2[nH]1', 'N=C1Nc2ccccc2N1'], 'c1ccc2c(c1)nc(N)[nH]2'),
    ('acetamidine', 67170, ['CC(=N)N', 'CC(N)=N'], None),
    ('N-methylacetamidine', 81059, ['CC(=N)NC', 'CC(=NC)N'], None),
    ('benzamidine', 2332, ['NC(=N)c1ccccc1', 'N=C(N)c1ccccc1'], None),
    ("N,N'-diphenylformamidine", 7535, ['C(=Nc1ccccc1)Nc1ccccc1', 'C(Nc1ccccc1)=Nc1ccccc1'], None),
    ('guanidine', 3520, ['NC(N)=N', 'N=C(N)N'], None),
    ('methylguanidine', 10111, ['CNC(N)=N', 'CN=C(N)N'], None),
    ('creatine', 586, ['CN(CC(=O)O)C(N)=N', 'CN(CC(O)=O)C(=N)N'], None),
    ('L-arginine', 6322, ['N[C@@H](CCCNC(N)=N)C(=O)O', 'N[C@@H](CCCN=C(N)N)C(=O)O'], None),
    ('metformin', 4091, ['CN(C)C(=N)NC(N)=N', 'CN(C)C(N)=NC(N)=N', 'CN(C)C(=N)N=C(N)N'], None),
    ('dicyandiamide', 10005, ['N#CN=C(N)N', 'N#CNC(N)=N'], None),
    ('amiloride', 16231, ['NC(=N)NC(=O)c1nc(Cl)c(N)nc1N', 'NC(N)=NC(=O)c1nc(Cl)c(N)nc1N'], None),
]


@mark.parametrize('name, cid, texts, expected', AMIDINES, ids=[row[0] for row in AMIDINES])
def test_an_amidine_or_guanidine_is_one_molecule_from_every_tautomer_and_drawing(name, cid, texts, expected):
    """Every tautomer, in random atom orders, aromatic and Kekule, gives one canonical string, which reads back
    to itself.  A conjugated ring takes the amino form -- an exocyclic imine hands its hydrogen back to the
    ring (`isomers:exocyclic-imine`) -- and a ring with no C=C keeps an N-H beside its C=O."""
    out = set()
    for text in texts:
        for m in drawings(text, 8):
            canonicalize(m)
            out.add(str(m))
    assert len(out) == 1, out
    [canonical] = out
    if expected is not None:
        assert canonical == expected
    back = smiles(canonical)
    canonicalize(back)
    assert str(back) == canonical


def test_an_exocyclic_imine_on_a_conjugated_ring_is_logged():
    m = smiles('N=C1C=CC=CN1')
    canonicalize(m)
    assert str(m) == 'c1cccnc1N'
    assert [r.atoms for r in m.log if r.rule == 'isomers:exocyclic-imine'] == [(1, 2, 7)]


def test_an_amidoxime_keeps_its_oxime():
    """A ring imine nitrogen bearing an oxygen is not handed over: 2-pyridone oxime stays an oxime."""
    m = smiles('ON=C1C=CC=CN1')
    canonicalize(m)
    assert not [r for r in m.log if r.rule == 'isomers:exocyclic-imine']
