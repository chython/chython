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
"""`standardize_isomers`: the placement is a choice, and two spellings must make the same one."""
from itertools import permutations, product
from chython.chemistry._canonicalize import canonicalize
from chython.chemistry._implicit import check_valence
from chython.chemistry._isomers import standardize_isomers
from chython.core import INFO, Log, read_smiles


def converge(a, b):
    """Both spellings through the pass; do they store the same molecule afterwards?"""
    ma, mb = read_smiles(a), read_smiles(b)
    standardize_isomers(ma)
    standardize_isomers(mb)
    return ma, mb, ma.canonical_bytes == mb.canonical_bytes


PAIRS = [('Cc1cnc[nH]1', 'Cc1c[nH]cn1', '4-methylimidazole'),
         ('c1cc2[nH]ncc2cn1', 'c1cc2n[nH]cc2cn1', 'pyrazolo[3,4-c]pyridine'),
         ('c1ccc2[nH]ncc2c1', 'c1ccc2n[nH]cc2c1', 'indazole'),
         ('Cc1n[nH]c2nc3[nH]nc(C)c3nc12', 'Cc1[nH]nc2nc3[nH]nc(C)c3nc12',
          'a bis-pyrazolo fused system with two mobile hydrogens'),
         ('Cc1cnc[n-]1', 'Cc1c[n-]cn1', '4-methylimidazolide -- a mobile CHARGE'),
         ('c1cc2[n-]ncc2cn1', 'c1cc2n[n-]cc2cn1', 'the same skeleton, anionic')]

PAIRS_SMILES = [(a, label) for a, _, label in PAIRS] + [(b, label) for _, b, label in PAIRS]


def test_two_annular_tautomers_of_one_compound_store_the_same_molecule():
    """The whole point.  Six pairs, each two valid drawings of one compound."""
    for a, b, label in PAIRS:
        ma, mb, same = converge(a, b)
        assert same, f'{label}: {a} and {b} still differ, {ma} vs {mb}'


def test_the_choice_does_not_depend_on_which_spelling_arrived():
    """The reference frame is a skeleton with every candidate stripped, so it cannot.

    Idempotence from either end -- not just equality of the two results -- is what proves the frame is
    spelling independent.
    """
    for a, b, label in PAIRS:
        ma, mb, _ = converge(a, b)
        assert not standardize_isomers(ma), f'{label}: not idempotent from {a}'
        assert not standardize_isomers(mb), f'{label}: not idempotent from {b}'


def test_a_molecule_with_nothing_to_place_is_untouched_and_says_so():
    # not an azaindole: `c1cc2[nH]ccc2cn1` has two sites, so the 1H and 7H forms are a real choice and
    # the pass does move its hydrogen.  Only molecules with nothing to decide belong here.
    for string in ['c1ccccc1', 'c1cc[nH]c1', 'c1ccncc1', 'c1ccoc1', 'CCO',
                   'c1ccc2[nH]ccc2c1']:                  # indole: one site, so no distribution
        mol = read_smiles(string)
        before = mol.canonical_bytes
        assert standardize_isomers(mol) is False, string
        assert mol.canonical_bytes == before, string


def test_only_a_placement_that_kekulises_is_ever_chosen():
    """Validity is proved by the complete backtracking kekuliser, never scored.

    1,2,4-triazole has three nitrogens and one hydrogen; not every placement has a Kekule form, and
    the pass must never return one that does not, however its ranks fall.
    """
    for string in ['c1nc[nH]n1', 'c1n[nH]cn1', 'c1[nH]ncn1']:
        mol = read_smiles(string)
        standardize_isomers(mol)
        assert mol.copy().kekule().unresolved == [], (string, str(mol))


def test_each_ring_system_is_decided_independently_of_the_others():
    """Limit independence: the reason groups are components rather than one pool.

    The eight-ring case must give each ring the same answer the one-ring case gives it, so no attempt
    budget may be shared across rings.
    """
    one = read_smiles('Cc1cc[nH]n1')
    standardize_isomers(one)
    many = read_smiles('.'.join(['Cc1cc[nH]n1'] * 8))
    standardize_isomers(many)
    # per component, not on `str(many)`: the SMILES writer picks its own start atom per component, so
    # the joined string differs even when every component stores the same molecule.
    components = [many.substructure(c) for c in many.connected_components]
    assert len(components) == 8
    for i, component in enumerate(components):
        assert component.canonical_bytes == one.canonical_bytes, f'ring {i} got a different answer'


def test_a_stereocentre_is_not_touched_and_needs_no_parity_restore():
    """Measured, and the measurement is why there is no parity restore in the pass.

    `set_hydrogens` and `set_charge` do not clear a parity -- unlike `delete_atom`, which is why
    `implicify_hydrogens` reads parities before its session and writes them back inside it.
    """
    for a, b in [('C[C@H](N)c1cnc[nH]1', 'C[C@H](N)c1c[nH]cn1'),
                 ('F[C@](Cl)(Br)c1cc2n[nH]cc2cn1', 'F[C@](Cl)(Br)c1cc2[nH]ncc2cn1')]:
        ma, mb = read_smiles(a), read_smiles(b)
        pa = {n: ma.parity_of(n) for n in ma.atom_numbers if ma.parity_of(n)}
        assert pa, a
        standardize_isomers(ma)
        standardize_isomers(mb)
        assert {n: ma.parity_of(n) for n in ma.atom_numbers if ma.parity_of(n)} == pa
        assert ma.canonical_bytes == mb.canonical_bytes, (a, b, str(ma), str(mb))


def test_an_epimer_does_not_become_its_own_mirror_image():
    """The other half of the stereo bar: unification must not reach across a stereocentre."""
    ma = read_smiles('C[C@H](N)c1cnc[nH]1')
    mb = read_smiles('C[C@@H](N)c1cnc[nH]1')
    standardize_isomers(ma)
    standardize_isomers(mb)
    assert ma.canonical_bytes != mb.canonical_bytes


def test_the_bool_means_changed_and_the_log_is_information_not_a_repair():
    """Both tautomers were valid molecules, so nothing here was wrong.  Severity says so."""
    mol = read_smiles('Cc1cnc[nH]1')
    log = mol.log
    assert standardize_isomers(mol) is True
    assert len(log) == 1
    assert log[0].severity is INFO
    assert 'placement' in log[0].rule


def test_nothing_is_logged_when_no_log_is_asked_for():
    mol = read_smiles('Cc1cnc[nH]1')
    assert standardize_isomers(mol) is True          # and it does not raise for want of a log


def test_the_kekule_form_is_decided_too_and_agrees_with_the_aromatic_one():
    """The aromatic form is not a precondition.

    `thiele()` refuses to aromatise a pyridone, so a lactam reaches this pass in its Kekule form however
    the pipeline is ordered.  The Kekule system is spelled aromatic on a working copy instead, which is
    also what makes the two answers the same answer -- `CC1=CN=CN1` is 4-methylimidazole drawn Kekule and
    must land where `Cc1cnc[nH]1` lands.
    """
    kekule, aromatic = read_smiles('CC1=CN=CN1'), read_smiles('Cc1cnc[nH]1')
    assert standardize_isomers(kekule) is True
    standardize_isomers(aromatic)
    # `thiele()` afterwards and not before: before, it would aromatise the Kekule one and there would be
    # nothing left to prove.  The pass answers where the hydrogen goes, not which Kekule form to draw.
    kekule.thiele()
    aromatic.thiele()
    assert kekule.canonical_bytes == aromatic.canonical_bytes, f'{kekule} vs {aromatic}'


#: Two Kekule drawings of one lactam, differing only in which ring nitrogen holds the hydrogen.  None of
#: these has an aromatic bond anywhere in it, so the aromatic `_sites` cannot see one site between them.
KEKULE_PAIRS = [
    ('CC1=CC=NC(=O)N1', 'CC1=NC(=O)NC=C1', '4-methylpyrimidin-2-one, N1-H and N3-H'),
    ('O=C1NC=CC(N)=N1', 'O=C1NC(N)=CC=N1', 'cytosine, N1-H and N3-H'),
    ('O=C1C=CN=CN1', 'O=C1C=CNC=N1', 'pyrimidin-4-one'),
    ('O=C1NC=NN1', 'O=C1NN=CN1', '1,2,4-triazol-3-one'),
    ('O=C1C=CNN1', 'O=C1NNC=C1', 'pyrazol-3-one'),
    ('O=C1NC=NC2=CC=CC=C12', 'O=C1N=CNC2=CC=CC=C12', 'quinazolin-4-one'),
    ('O=C1NC=NC2=C1NC=N2', 'O=C1N=CNC2=C1NC=N2', 'hypoxanthine, N1-H and N3-H'),
    ('NC1=NC2=C(N=CN2)C(=O)N1', 'NC1=NC2=C(NC=N2)C(=O)N1', 'guanine, the imidazole hydrogen moving'),
]


def test_two_kekule_tautomers_of_one_lactam_store_the_same_molecule():
    """The shape the pyridone preference put out of reach: a mobile hydrogen with no aromatic bond.

    Through `canonicalize()` and not the pass alone, deliberately.  A Kekule structure is not canonical
    by itself -- quinazolin-4-one drawn two ways differs in its benzo ring's double bonds as well as in
    its mobile hydrogen -- and choosing between two Kekule forms belongs to `thiele()` and `kekule()`.
    This pass decides only where the hydrogen goes, and the pipeline is where a caller reads a key.
    """
    for a, b, label in KEKULE_PAIRS:
        ma, mb = read_smiles(a), read_smiles(b)
        canonicalize(ma)
        canonicalize(mb)
        assert ma.canonical_bytes == mb.canonical_bytes, \
            f'{label}: {a} and {b} still differ, {ma} vs {mb}'


AROMATIC_PAIRS = [
    ('O=C1NC(=O)c2cnccc2N1', 'O=C1NC(=O)C2=CNC=CC2=N1', 'pyrido[4,3-d]pyrimidine-2,4-dione'),
    ('O=C1Nc2cc3ncccc3cc2N1', 'O=C1N=C2C=C3C(=CC=CN3)C=C2N1', 'imidazo[4,5-g]quinolin-2-one'),
]


def test_a_fused_lactam_keeps_its_pyridine_aromatic():
    """A hydrogen on a lactam nitrogen does not move onto a fused pyridine nitrogen, whichever form arrives.

    Both placements kekulise, so the canonical order alone would take whichever ranks lower; the imine
    leaves the pyridine without a sextet and is the form that loses.
    """
    for aromatic, quinoid, label in AROMATIC_PAIRS:
        ma, mb = read_smiles(aromatic), read_smiles(quinoid)
        canonicalize(ma)
        canonicalize(mb)
        assert ma.aromatic_rings_count == read_smiles(aromatic).aromatic_rings_count, f'{label}: {ma}'
        assert ma.canonical_bytes == mb.canonical_bytes, f'{label}: {ma} vs {mb}'


AZOLE_FIRST = [('c1cnc2[nH]ccc2c1', 'C1=CNC2=NC=CC2=C1', '7-azaindole'),
               ('OC(=O)c1ccc2cc[nH]c2n1', 'OC(=O)C1=CC=C2C=CN=C2N1', '7-azaindole-6-carboxylic acid'),
               ('NC(=O)c1ccc2cc[nH]c2n1', 'NC(=O)C1=CC=C2C=CN=C2N1', '7-azaindole-6-carboxamide')]


def test_an_azole_nh_outranks_an_azine_nh():
    """Both forms are aromatic, so the ring size decides: the hydrogen sits on the five-membered ring."""
    for azole, azine, label in AZOLE_FIRST:
        ma, mb = read_smiles(azole), read_smiles(azine)
        canonicalize(ma)
        canonicalize(mb)
        assert ma.canonical_bytes == mb.canonical_bytes, f'{label}: {ma} vs {mb}'
        nh = [n for n in ma.atom_numbers if ma.element_of(n) == 7 and ma.total_h_of(n)
              and len(tuple(ma.neighbors_of(n))) == 2]
        assert [5 in ma.ring_sizes_of(n) for n in nh] == [True], f'{label}: {ma}'


def test_a_kekule_placement_never_changes_the_formula():
    """The working copy is spelled aromatic and kekulised back, and `kekule()` repairs when it must.

    1,2,4-triazol-3-one is the case that taught this: protonating both nitrogens flanking its lone ring
    carbon leaves that carbon no partner for a double bond, and the kekuliser's relaxation hid the fact
    by taking the two hydrogens away -- a placement of a different compound.
    """
    for a, b, label in KEKULE_PAIRS:
        for string in (a, b):
            mol = read_smiles(string)
            before = mol.brutto_formula
            standardize_isomers(mol)
            assert mol.brutto_formula == before, f'{label}: {string} became {mol.brutto_formula}'


def test_the_kekule_path_is_idempotent_from_either_end():
    for a, b, label in KEKULE_PAIRS:
        ma, mb, _ = converge(a, b)
        assert not standardize_isomers(ma), f'{label}: not idempotent from {a}'
        assert not standardize_isomers(mb), f'{label}: not idempotent from {b}'


def test_a_ring_with_an_sp3_carbon_is_refused_whole():
    """The gate that keeps this module out of tautomerism it has no business doing.

    Spelling these aromatic would move a hydrogen off a carbon, which is a repair with a chemical opinion
    in it and belongs to `standardize()`'s table if it belongs anywhere.
    """
    for string, label in [('O=C1CC(=O)NC(=O)N1', 'barbituric acid keeps its CH2'),
                          ('O=C1CC=CC=C1', 'cyclohexa-2,4-dien-1-one does not become phenol'),
                          ('O=C1CN=C2C=CC=CC2=N1', 'quinoxalin-2-one keeps its CH2')]:
        mol = read_smiles(string)
        before = mol.canonical_bytes
        standardize_isomers(mol)
        assert mol.canonical_bytes == before, f'{label}: {string} became {mol}'


#: One amidine or guanidine written two ways.  No ring, so no kekuliser: which nitrogen may take the
#: double bond is decided by its valence and nothing else.
#:
#: The last three are the shape where two units share a nitrogen: the bridge is single-bonded to both
#: carbons, so each unit alone would accept it, and only one of them may.
AMIDINE_PAIRS = [
    ('CC(N)=NC', 'CC(=N)NC', 'N-methylacetamidine'),
    ('CCN=C(N)NC', 'CCNC(=NC)N', 'N-ethyl-N-methylguanidine'),
    ('COC(=N)NC', 'COC(N)=NC', 'O-methyl-N-methylisourea'),
    ('CNC(N)=NC(=N)NC', 'CNC(=N)NC(=N)NC', 'a biguanide, two amidines in one molecule'),
    ('CC(=NC)NCC', 'CC(NC)=NCC', 'N,N-disubstituted, so both nitrogens carry one hydrogen'),
    ('N=C(N)NC(=N)N', 'NC(=N)N=C(N)N', 'biguanide, the bridge single-bonded to both carbons'),
    ('CN(C)C(=N)NC(=N)N', 'CN(C)C(N)=NC(=N)N', 'metformin'),
    ('N=C(N)NC(=N)NCCc1ccccc1', 'NC(=N)N=C(N)NCCc1ccccc1', 'phenformin'),
]


def test_two_spellings_of_one_amidine_store_the_same_molecule():
    for a, b, label in AMIDINE_PAIRS:
        ma, mb, same = converge(a, b)
        assert same, f'{label}: {a} and {b} still differ, {ma} vs {mb}'


def test_the_amidine_path_conserves_the_formula_and_is_idempotent():
    for a, b, label in AMIDINE_PAIRS:
        for string in (a, b):
            mol = read_smiles(string)
            before = mol.brutto_formula
            standardize_isomers(mol)
            assert mol.brutto_formula == before, f'{label}: {string} became {mol.brutto_formula}'
            assert not standardize_isomers(mol), f'{label}: not idempotent from {string}'


def test_two_amidines_sharing_a_nitrogen_do_not_both_double_bond_it():
    """A bridging nitrogen may accept one C=N, so the two units are one placement problem.

    Decided as two, each unit picks the bridge -- it is the acceptor its own ranking prefers -- and the
    plan writes two double bonds onto one nitrogen.
    """
    for a, b, label in AMIDINE_PAIRS:
        for string in (a, b):
            mol = read_smiles(string)
            standardize_isomers(mol)
            assert not check_valence(mol), f'{label}: {string} became {mol}'


#: An amidine hung on a ring whose own hydrogen is mobile -- the two paths' gates do not overlap, but
#: their answers do, each mobile hydrogen sitting in the frame that ranks the other's placements.
COUPLED = [
    ('NC(=N)Nc1nc[nH]n1', 'NC(=N)Nc1n[nH]cn1', '3-guanidino-1H-1,2,4-triazole'),
    ('NC(=N)NC(=N)Nc1cc[nH]n1', 'NC(=N)NC(=N)Nc1n[nH]cc1', '1-(1H-pyrazol-3-yl)biguanide'),
]


def test_a_ring_and_an_amidine_on_it_are_decided_in_one_frame():
    """Two frames make the pass depend on its own output: idempotence is what catches that."""
    for a, b, label in COUPLED:
        for string in (a, b):
            mol = read_smiles(string)
            standardize_isomers(mol)
            assert not standardize_isomers(mol), f'{label}: not idempotent from {string}'
        ma, mb, same = converge(a, b)
        assert same, f'{label}: {a} and {b} still differ, {ma} vs {mb}'


def test_an_amide_is_not_an_amidine():
    """The double bond has to go to a nitrogen this group may move it to.

    An amide's is to oxygen, a nitro group's nitrogen carries its own, and a tertiary nitrogen has no
    room for one -- so none of these is a placement and all three are left alone.
    """
    for string, label in [('CC(=O)NC', 'an amide'),
                          ('CC(=O)N', 'a primary amide'),
                          ('CN(C)C(=N)N(C)C', 'both nitrogens tertiary but one'),
                          ('C[N+](=O)[O-]', 'nitromethane'),
                          ('CN=NC', 'an azo compound'),
                          ('N=C=NC', 'a carbodiimide, two double bonds on one carbon')]:
        mol = read_smiles(string)
        before = mol.canonical_bytes
        standardize_isomers(mol)
        assert mol.canonical_bytes == before, f'{label}: {string} became {mol}'


#: One ion, two drawings.  A charged nitrogen's hydrogen is as mobile as a neutral one's, and the charge
#: moves with it: which nitrogen of a guanidinium carries the `+` is a fact about the spelling and not
#: about the ion, so a placement has to choose the charge as well as the hydrogen.  Both paths, the ring
#: and the amidine, and both signs.
CHARGED_PAIRS = [
    ('NC(=[NH2+])N', '[NH3+]C(=N)N', 'guanidinium'),
    ('CNC(=[NH2+])N', 'CN=C([NH3+])N', 'N-methylguanidinium'),
    ('CNC(=[NH2+])N', 'C[NH2+]C(=N)N', 'N-methylguanidinium, the charge on the substituted nitrogen'),
    ('CC(=[NH2+])N', 'CC(=N)[NH3+]', 'acetamidinium'),
    ('CN(C)C(=[NH2+])NC(=N)N', 'CN(C)C(=N)NC(=[NH2+])N', 'metformin, protonated as it is at pH 7'),
    ('CN(C)C(=[NH2+])NC(=N)N', 'CN(C)C(=N)N=C(N)[NH3+]', 'metformin cation, the charge at the far end'),
    ('NC(=N)[NH-]', '[N-]=C(N)N', 'guanidide'),
    ('CC(=N)[NH-]', 'CC(N)=[N-]', 'acetamidate'),
    ('NC(=N)NC(=N)[NH-]', '[N-]=C(N)NC(=N)N', 'biguanidide, the anion two units away'),
    ('CN(C)C(=[NH2+])N', 'C[NH+](C)C(=N)N', 'N,N-dimethylguanidinium, the charge on the tertiary nitrogen'),
    ('CN(C)C(=[NH2+])NC(=N)N', 'C[NH+](C)C(=N)NC(=N)N', 'metformin cation, the charge on the tertiary N'),
    ('Cc1c[nH+]c[nH]1', 'Cc1c[nH]c[nH+]1', '4-methylimidazolium'),
    ('Cn1cc[nH+]c1', 'C[n+]1cc[nH]c1', '1-methylimidazolium: the charge is mobile where the H is not'),
    ('Cn1ccc[nH+]1', 'C[n+]1ccc[nH]1', '1-methylpyrazolium, the same on adjacent nitrogens'),
    ('c1[nH+]c[nH]n1', 'c1[nH]c[nH+]n1', '1,2,4-triazolium'),
    ('Cc1cc[nH+][nH]1', 'Cc1cc[nH][nH+]1', '3-methylpyrazolium'),
    ('c1cc2[nH+]ncc2cn1', 'c1cc2n[nH+]cc2cn1', 'pyrazolo[3,4-c]pyridinium'),
]

#: A charge this pass does not own, on an atom no placement may reach: too many neighbours for a site,
#: an element whose hydrogen is never mobile here, or a component with no site in it at all.
FIXED_CHARGES = [
    ('C[N+](C)(C)CCNC(=N)N', 'a quaternary ammonium beside a guanidine'),
    ('NC(=N)NCC(=O)[O-]', 'a carboxylate beside a guanidine -- oxygen to nitrogen is not this pass'),
    ('C[N+](=O)[O-]', 'nitromethane, a charge-separated group and not a placement'),
    ('CC(=O)[O-].NC(=N)N', 'a salt: the anion is in the component without the sites'),
]


def test_two_drawings_of_one_ion_store_the_same_molecule():
    """The point of the neutral pairs, for an ion: a guanidinium at pH 7 is one compound, not three."""
    for a, b, label in CHARGED_PAIRS:
        ma, mb, same = converge(a, b)
        assert same, f'{label}: {a} and {b} still differ, {ma} vs {mb}'


def test_a_charged_placement_is_idempotent_from_either_end():
    for a, b, label in CHARGED_PAIRS:
        for string in (a, b):
            mol = read_smiles(string)
            standardize_isomers(mol)
            assert not standardize_isomers(mol), f'{label}: not idempotent from {string}'


def test_a_charged_placement_conserves_the_formula_and_every_component_charge():
    """A charge may move between the sites of one system and nowhere else.

    Net charge is `neutralize()`'s business and this pass has none of it: the count of charges and their
    signs are read off the group and put back, so a component's charge is an invariant here.
    """
    for a, b, label in CHARGED_PAIRS + [('CC(=O)[O-].NC(=[NH2+])N', '', 'guanidinium acetate')]:
        for string in (a, b) if b else (a,):
            mol = read_smiles(string)
            formula = mol.brutto_formula
            charges = sorted(sum(mol.charge_of(n) for n in c) for c in mol.connected_components)
            standardize_isomers(mol)
            assert mol.brutto_formula == formula, f'{label}: {string} became {mol.brutto_formula}'
            assert sorted(sum(mol.charge_of(n) for n in c)
                          for c in mol.connected_components) == charges, \
                f'{label}: {string} moved a charge between components, {mol}'


def test_a_charged_placement_adds_no_valence_error():
    """Judged on the kekule form, because `valence_rules.tsv` has no row for an aromatic `[nH+]`:
    `check_valence` answers `unknown` for `c1cc[nH+]cc1` as read and `[]` for its kekule form, so the
    aromatic spelling cannot tell a placement a nitrogen can carry from one it cannot.

    The kinds and not the atoms, and no better than the input: a fused pyrazolium kekulises to a form the
    table already objects to whichever nitrogen holds its `+`.  A nitrogen asked for a hydrogen or a charge
    it has no room for is what this catches -- 1-methylpyrazolium handing its hydrogen to its methylated
    nitrogen, four bonds on a neutral nitrogen, which `kekule()` finds a form for and does not object to.
    """
    for a, b, label in CHARGED_PAIRS:
        for string in (a, b):
            mol, placed = read_smiles(string), read_smiles(string)
            mol.kekule()
            before = sorted(kind for _, kind in check_valence(mol))
            standardize_isomers(placed)
            placed.kekule()
            assert sorted(kind for _, kind in check_valence(placed)) == before, \
                f'{label}: {string} became {placed}, {check_valence(placed)}'


def test_a_charge_the_pass_does_not_own_stays_on_its_atom():
    """Mobile means mobile over the sites of one system, and a charge elsewhere is somebody else's."""
    for string, label in FIXED_CHARGES:
        mol = read_smiles(string)
        before = {n: mol.charge_of(n) for n in mol.atom_numbers if mol.charge_of(n)}
        standardize_isomers(mol)
        after = {n: mol.charge_of(n) for n in mol.atom_numbers if mol.charge_of(n)}
        assert after == before, f'{label}: {string} moved a charge, {before} -> {after}'


#: Species whose EVERY drawing has to collapse, enumerated rather than listed: a pair test only ever
#: covers the two spellings somebody thought of, and the ones that split were the ones nobody drew.  The
#: count is what `redrawings` finds today, asserted as a floor so a narrowed space cannot quietly make
#: this test pass by enumerating one drawing.
EXHAUSTIVE = [
    ('NC(=[NH2+])N', 2, 'guanidinium'),
    ('CC(=[NH2+])N', 2, 'acetamidinium'),
    ('CC(=N)[NH-]', 2, 'acetamidate'),
    ('NC(=N)[NH-]', 2, 'guanidide'),
    ('CC(N)=NC', 2, 'N-methylacetamidine'),
    ('NC(=N)NC(N)=N', 2, 'biguanide'),
    ('NC(=N)NC(N)=[NH2+]', 7, 'biguanidinium'),
    ('NC(=N)NC([NH-])=N', 6, 'biguanidide'),
    ('NC(=N)NC(=N)NC(=N)N', 5, 'triguanide, the flip-flop chain'),
    ('CN(C)C(=[NH2+])NC(N)=N', 16, 'metformin cation, a tertiary nitrogen among the sites'),
    ('C[N+](C)(C)CCNC(N)=N', 2, 'a quaternary ammonium beside a guanidine'),
    ('NC(=N)NCC(=O)[O-]', 2, 'a carboxylate beside a guanidine'),
    ('Cc1cnc[nH]1', 2, '4-methylimidazole'),
    ('Cc1c[nH+]c[nH]1', 2, '4-methylimidazolium'),
    ('Cc1cc[nH+][nH]1', 2, '3-methylpyrazolium'),
    ('Cn1cc[nH+]c1', 2, '1-methylimidazolium, a mobile charge on a substituted nitrogen'),
    ('Cn1ccc[nH+]1', 2, '1-methylpyrazolium, the same on adjacent nitrogens'),
    ('CCn1cc[n+](C)c1', 2, '1-ethyl-3-methylimidazolium: both nitrogens substituted, only the + moves'),
    ('c1[nH+]c[nH]n1', 3, '1,2,4-triazolium'),
    ('c1n[n-]cn1', 2, '1,2,4-triazolide'),
    ('c1cc2[nH]ncc2cn1', 3, 'pyrazolo[3,4-c]pyridine'),
    ('c1cc2[nH+]ncc2cn1', 9, 'pyrazolo[3,4-c]pyridinium'),
    ('Cc1n[nH]c2nc3[nH]nc(C)c3nc12', 7, 'a bis-pyrazolo fused system with two mobile hydrogens'),
]


def redrawings(string):
    """Every valid redrawing of one species: the same skeleton, the same charges over its nitrogens in any
    arrangement, the same nitrogen-bound hydrogen count in any distribution, and every double bond an
    acyclic amidine carbon could hold.  Aromatic as read, so a ring's bond orders are free too.

    Valid means the kekule form is no worse than the source's under `check_valence` -- which is the
    independent half, and the half that catches a hydrogen a nitrogen has no room for -- and that
    `kekule()` did not have to rewrite the counts to get there, since a form it repaired into existence is
    a form of another molecule.  Two exclusions, both deliberate: an arrangement changing the charge
    MULTISET invents a charge separation, which is `fix_resonance()`'s to remove, and one costing the
    source its aromatic ring is the flip out of an aromatic ring onto an sp3 nitrogen -- 4-methyl-
    imidazolium with both hydrogens on one nitrogen -- which no shape here decides.
    """
    base = read_smiles(string)
    reference = base.copy()
    reference.kekule()
    kinds = sorted(kind for _, kind in check_valence(reference))
    aromatic = base.aromatic_rings_count
    sites = tuple(n for n in base.atom_numbers if base.element_of(n) == 7 and not base.radical_of(n))
    charges = set(permutations([base.charge_of(n) for n in sites]))
    hydrogens = sum(base.implicit_h_of(n) or 0 for n in sites)
    amidines = []
    for c in base.atom_numbers:
        if base.element_of(c) != 6:
            continue
        near = tuple(n for n in base.neighbors_of(c) if n in sites and not base.bond_in_ring(c, n))
        if len(near) > 1:
            amidines.append((c, near))
    out = {}
    for qs, hs in product(charges, _spread(len(sites), hydrogens)):
        for choice in product(*([None, *near] for _, near in amidines)):
            taken = [n for n in choice if n is not None]
            if len(taken) != len(set(taken)):
                continue                  # one nitrogen cannot take two double bonds
            work = base.copy()
            with work.edit():
                for c, near in amidines:
                    for n in near:
                        work.set_order(c, n, 1)
                for (c, _), n in zip(amidines, choice):
                    if n is not None:
                        work.set_order(c, n, 2)
                for n, h, q in zip(sites, hs, qs):
                    work.set_charge(n, q)
                    work.set_hydrogens(n, h)
            asked = {n: (h, q) for n, h, q in zip(sites, hs, qs)}
            if work.kekule().unresolved:
                continue
            if any((work.implicit_h_of(n) or 0, work.charge_of(n)) != state for n, state in asked.items()):
                continue                  # kekulised only because it was repaired
            if sorted(kind for _, kind in check_valence(work)) != kinds:
                continue
            work.thiele()
            if work.aromatic_rings_count == aromatic:
                out[work.canonical_bytes] = work
    return list(out.values())


def _spread(sites, total):
    """Every way of distributing `total` hydrogens over `sites` nitrogens, at most three each."""
    if sites == 1:
        if total <= 3:
            yield (total,)
        return
    for h in range(min(3, total) + 1):
        for tail in _spread(sites - 1, total - h):
            yield (h, *tail)


def test_every_drawing_of_one_species_makes_the_same_placement():
    """The claim the pair tests only sample: over the whole space of valid drawings, one form comes out."""
    for string, count, label in EXHAUSTIVE:
        drawings = redrawings(string)
        assert len(drawings) >= count, f'{label}: {len(drawings)} drawings of {string}, expected {count}'
        placed = {}
        for mol in drawings:
            drawn = mol.smiles
            standardize_isomers(mol)
            placed.setdefault(mol.canonical_bytes, []).append(drawn)
        assert len(placed) == 1, f'{label}: {string} splits {len(placed)} ways, {list(placed.values())}'


def test_the_method_on_the_container_is_the_registered_pass():
    """Registration is by injection -- `chemistry` calls a setter, `core` never names `chemistry`."""
    mol = read_smiles('Cc1cnc[nH]1')
    assert mol.standardize_isomers() is True
    assert mol.standardize_isomers() is False


def test_the_method_forwards_the_log():
    mol = read_smiles('Cc1cnc[nH]1')
    log = mol.log
    mol.standardize_isomers()
    assert len(log) == 1


# Each shape is the smallest public molecule carrying the defect, and that is the only form one enters
# this file in.
#
# These are not tautomers, hence a regression net rather than `PAIRS`: they are placements with no
# Kekule form, made valid by `kekule()`'s two relaxations -- a neutral atom that must be a cation, and
# a site holding a hydrogen the ring cannot afford.
CORPUS_SHAPES = [
    ('c1ccn(C)cc1', 'N-methylpyridinium drawn without its charge'),
    ('c1ccn(CC(=O)N)cc1', 'an N-acyl-methyl pyridinium, the nicotinamide-conjugate shape'),
    ('NC(=O)c1cccn(C)c1', 'nicotinamide N-methylated -- the NAD(+) shape without the nucleotide'),
    ('c1ccn(O)cc1', 'N-hydroxypyridine, which is pyridine N-oxide charge separated'),
    ('c1ccn(N)cc1', '1-aminopyridinium drawn neutral'),
    ('Cn1cc[nH]c1', 'an N-methylimidazole carrying a second hydrogen it cannot afford'),
    ('Cc1cn(C)c[nH]1', 'the same, substituted'),
    ('c1cc[nH]cc1', 'a six-ring nitrogen with a hydrogen there is no room for'),
    ('Cn1cscc1', 'an N-methyl thiazolium drawn neutral'),
]


def test_every_corpus_shape_kekulises_after_the_pipeline():
    """The two relaxations, measured on the public analogue of every failing shape."""
    for string, label in CORPUS_SHAPES:
        mol = read_smiles(string)
        assert mol.kekule().unresolved == [], f'{label}: {string}'


def test_the_repairs_are_honest_about_valence():
    """A relaxation that invented a valid Kekule form by breaking a valence is not a repair.

    Checked on the Kekule form deliberately: an atom carrying an aromatic bond answers `'unknown'` for
    want of a table row, so an aromatic check cannot tell a violation from a gap.
    """
    for string, label in CORPUS_SHAPES:
        mol = read_smiles(string)
        mol.kekule()
        assert [(n, v) for n, v in check_valence(mol) if v == 'violation'] == [], f'{label}: {string}'


def test_every_corpus_shape_survives_a_canonical_smiles_round_trip():
    """Non-corruption: the canonical SMILES reads back as the same molecule."""
    for string, label in CORPUS_SHAPES:
        mol = read_smiles(string)
        canonicalize(mol)
        again = read_smiles(str(mol))
        assert again.canonical_bytes == mol.canonical_bytes, f'{label}: {mol}'


def test_the_whole_pipeline_is_idempotent_on_every_shape():
    for string, label in CORPUS_SHAPES:
        mol = read_smiles(string)
        canonicalize(mol)
        first = mol.canonical_bytes
        assert canonicalize(mol) is False, f'{label}: {string}'
        assert mol.canonical_bytes == first, f'{label}: {string}'


STEREO_CASES = [
    ('C[C@H](N)c1cnc[nH]1', 'a stereocentre beside a ring whose hydrogen moves'),
    ('F[C@](Cl)(Br)c1cc2n[nH]cc2cn1', 'a quaternary centre beside a fused mobile system'),
    ('C[C@H](N)c1cc2n[nH]cc2cn1.C[C@@H](O)c1c[nH]cn1', 'two components, both moving'),
    ('C/C=C/c1cc2n[nH]cc2cn1', 'a double-bond geometry beside a mobile system'),
]


def test_no_stereo_descriptor_is_lost_or_changed_by_the_pipeline():
    for string, label in STEREO_CASES:
        mol = read_smiles(string)
        before = {n: mol.parity_of(n) for n in mol.atom_numbers if mol.parity_of(n)}
        assert before, f'{label}: the fixture states no stereo, so it proves nothing'
        canonicalize(mol)
        after = {n: mol.parity_of(n) for n in mol.atom_numbers if mol.parity_of(n)}
        assert after == before, f'{label}: {before} became {after}'


def test_the_two_enantiomers_stay_two_compounds_through_the_pipeline():
    """Unification must reach across a tautomeric shift and never across a stereocentre."""
    for base in ['C[C@H](N)c1cnc[nH]1', 'F[C@](Cl)(Br)c1cc2n[nH]cc2cn1']:
        left = read_smiles(base)
        right = read_smiles(base.replace('[C@]', '[C@@]').replace('[C@H]', '[C@@H]'))
        canonicalize(left)
        canonicalize(right)
        assert left.canonical_bytes != right.canonical_bytes, base


def test_the_atom_count_and_the_formula_never_change():
    """`standardize_isomers` moves a hydrogen and must never add or remove one.

    Fixtures are the tautomer pairs and the stereo cases, deliberately not `CORPUS_SHAPES`, where
    `kekule()`'s surplus-hydrogen relaxation removes one on purpose and the formula must change.
    """
    for string in [s for s, _ in PAIRS_SMILES] + [s for s, _ in STEREO_CASES]:
        mol = read_smiles(string)
        heavy, formula = len(mol), mol.brutto
        canonicalize(mol)
        assert len(mol) == heavy, string
        assert mol.brutto == formula, string


# Do not add a `chython.standardize_isomers is standardize_isomers` test here:
# `test_dependency_direction.py` forbids anything under `chython/chemistry/` from importing the facade,
# tests included.  The facade re-export is checked from `chython/test/`.
