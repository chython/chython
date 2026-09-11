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
AMIDINE_PAIRS = [
    ('CC(N)=NC', 'CC(=N)NC', 'N-methylacetamidine'),
    ('CCN=C(N)NC', 'CCNC(=NC)N', 'N-ethyl-N-methylguanidine'),
    ('COC(=N)NC', 'COC(N)=NC', 'O-methyl-N-methylisourea'),
    ('CNC(N)=NC(=N)NC', 'CNC(=N)NC(=N)NC', 'a biguanide, two amidines in one molecule'),
    ('CC(=NC)NCC', 'CC(NC)=NCC', 'N,N-disubstituted, so both nitrogens carry one hydrogen'),
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
