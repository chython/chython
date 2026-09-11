"""
Benchmark: chython tautomer canonicalization vs RDKit TautomerEnumerator

Correctness: 62 test cases from chython's own test suite (test_isomers.py + test_groups.py)
Speed: throughput in mol/s over the full test set

Usage:
    python bench_tautomers.py
"""
import timeit
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')
from rdkit.Chem.MolStandardize import rdMolStandardize
from chython import smiles as chython_smiles, MoleculeContainer

# ---------------------------------------------------------------------------
# Test cases
# ---------------------------------------------------------------------------

# From test_isomers.py — canonical H/charge positioning (standardize_isomers)
# Comparison: standardize_isomers(inp) == prepare(expected)   [per test_isomers.py]
isomer_cases = [
    # fixed charge rules
    ('N1C=CN2[NH+]=CC=C12',     'N1C=CC2=[NH+]C=CN12',   'charge_fixed'),
    ('N(C)1C=CN2[N+](C)=CC=C12','N(C)1C=CC2=[N+](C)C=CN12','charge_fixed'),
    ('N1C=C[N+]2=C1C=CN2',     'N1C=CC2=[NH+]C=CN12',   'charge_fixed'),
    ('N1C=C2C=CN[N+]2=C1',     'N1C=CC2=C[NH+]=CN12',   'charge_fixed'),
    ('N1C=C2NC=C[N+]2=C1',     'N1C=CN2C=[NH+]C=C12',   'charge_fixed'),
    ('N1C2=[NH+]C=CC2=CC=C1',  'N1C=CC2=CC=C[NH+]=C12', 'charge_fixed'),
    ('N1C=C2C(=CC=[NH+]2)C=C1','N1C=CC2=CC=[NH+]C=C12', 'charge_fixed'),
    ('C1=CC=2C(=[NH+]1)C=CNC=2','N1C=CC2=C[NH+]=CC=C12','charge_fixed'),
    ('C=1C=[NH+]C=2C=1NC=CC=2','N1C=CC2=[NH+]C=CC=C12', 'charge_fixed'),
    ('C1=2C=[NH+]C=C1NC=CC=2', 'N1C=C2C=CC=[NH+]C2=C1','charge_fixed'),
    ('C1=2C=CNC=C1C=[NH+]C=2', 'N1C=C2C=C[NH+]=CC2=C1','charge_fixed'),
    # Morgan charge rules
    ('N1C=CC2=CC=[NH+]N12',    'N1C=CC2=CC=[NH+]N12',   'charge_morgan'),
    ('N1C=CC2=[N+]1NC=C2',     'N1C=CC2=CC=[NH+]N12',   'charge_morgan'),
    ('N1C=CN2C=C[NH+]=C12',    'N1C=CN2C=C[NH+]=C12',   'charge_morgan'),
    ('N1C=C[N+]2=C1NC=C2',     'N1C=CN2C=C[NH+]=C12',   'charge_morgan'),
    ('C=1N(C)C=[N+](CC)C=1',   'C=1N(C=[N+](C)C=1)CC', 'charge_morgan'),
    ('C=1N(C=[N+](C)C=1)CC',   'C=1N(C=[N+](C)C=1)CC', 'charge_morgan'),
    ('C=1N(C)C=[NH+]C=1',      '[N+]1(=CNC=C1)C',       'charge_morgan'),
    ('C=1N([N+](C)=CC=1)CC',   'C1=CC=[N+](CC)N1C',     'charge_morgan'),
    ('C1=CC=[N+](CC)N1C',      'C1=CC=[N+](CC)N1C',     'charge_morgan'),
    ('C1=CC=[N+](C)N1C',       'C1=CC=[N+](C)N1C',      'charge_morgan'),
    # ferrocene
    ('[CH-]1C=CC=C1.[Fe+2].[CH-]1C=CC=C1', '[CH-]1C=CC=C1.[Fe+2].[CH-]1C=CC=C1', 'ferrocene'),
    ('[CH-]1C=CC=C1.[Fe+2].C1=C[CH-]C=C1', '[CH-]1C=CC=C1.[Fe+2].[CH-]1C=CC=C1', 'ferrocene'),
    # fixed tautomer (triazole, tetrazole)
    ('N1C=NC=N1',  'N1C=NN=C1',      'triazole'),
    ('CC1=NC=NN1', 'N1=CNC(C)=N1',   'triazole'),
    ('CC1=NN=CN1', 'N1=CNC(C)=N1',   'triazole'),
    ('N1N=CN=N1',  'N1C=NN=N1',      'tetrazole'),
    ('CC1=NNN=N1', 'C=1(C)NN=NN=1',  'tetrazole'),
    ('CC1=NN=NN1', 'C=1(C)NN=NN=1',  'tetrazole'),
    # Morgan tautomer (pyrazole, imidazole, triazole)
    ('CC1=NNC=C1', 'N1C=CC(C)=N1',   'pyrazole'),
    ('CC1=CC=NN1', 'N1C=CC(C)=N1',   'pyrazole'),
    ('CC1=CN=CN1', 'C=1N=CNC=1C',    'imidazole'),
    ('CC1=CNC=N1', 'C=1N=CNC=1C',    'imidazole'),
    ('CC1=CN=NN1', 'N1N=NC(C)=C1',   'triazole_morgan'),
    ('CC1=CNN=N1', 'N1N=NC(C)=C1',   'triazole_morgan'),
    ('CC1=NNN=C1', 'N1N=NC(C)=C1',   'triazole_morgan'),
    # amidine/guanidine
    ('COC(=N)NC',  'COC(N)=NC',      'amidine'),
    ('CCN=C(N)NC', 'CCNC(=NC)N',     'amidine'),
    ('CCNC(=N)NC', 'CCNC(=NC)N',     'amidine'),
    ('CCNC(N)=NC', 'CCNC(=NC)N',     'amidine'),
    ('CNC(N)=NC(=N)NC', 'CNC(=N)NC(=N)NC', 'amidine'),
    ('CCN=CNC=NC', 'CCN=CN=CNC',     'amidine'),
]

# From test_groups.py — structural tautomers (canonicalize)
# Comparison: canonicalize(inp) == canonicalize(expected)
group_cases = [
    ('N=CO',   'NC=O',   'amide'),
    ('N=CS',   'NC=S',   'thioamide'),
    ('OC=C',   'O=CC',   'enol'),
    ('OC(C)=C','O=C(C)C','enol'),
    ('OC1=CC=NC=C1', 'O=C1C=CNC=C1', 'hydroxypyridine'),
    ('OC1=CC=CC=N1',    'O=C1NC=CC=C1',    'pyridone'),
    ('OC1=CC=NC=N1',    'O=C1NC=NC=C1',    'pyridone'),
    ('OC1=NC=CC=N1',    'O=C1N=CC=CN1',    'pyridone'),
    ('OC1=C(O)N=CC=N1', 'O=C1NC=CNC1=O',  'pyridone'),
    ('OC1=NC=CN=C1O',   'O=C1NC=CNC1=O',  'pyridone'),
    ('OC1=CC=NC(=O)N1', 'O=C1NC=CC(=O)N1','pyridone'),
    ('OC1=CC=NC(O)=N1', 'O=C1NC=CC(=O)N1','pyridone'),
    ('CN1C=CC(O)=N1', 'CN1NC(=O)C=C1', 'lactam_5'),
    ('CN1N=CC=C1O',   'CN1NC=CC1=O',   'lactam_5'),
    ('O=C1C=CCC=N1', 'O=C1NC=CC=C1',  'ring_keto'),
    ('O=C1CC=CC=N1', 'O=C1NC=CC=C1',  'ring_keto'),
    ('O=C1CC=NC=C1', 'O=C1C=CNC=C1',  'ring_keto'),
    ('C1C=NC=N1',    'N1C=CN=C1',     'imidazoline'),
    ('N=C1NC=CC(=O)N1', 'NC1=NC=CC(=O)N1', 'cytosine'),
    ('CN1N=CCC1=O', 'CN1NC=CC1=O',   'lactam_5'),
]

ALL_CASES = isomer_cases + group_cases

ISOMER_CATEGORIES = {'charge_fixed', 'charge_morgan', 'ferrocene', 'triazole', 'tetrazole',
                     'pyrazole', 'imidazole', 'triazole_morgan', 'amidine'}

# ---------------------------------------------------------------------------
# RDKit's own canonical tautomer test suite (canonTautomerData + testGithub3755)
# Source: Code/GraphMol/MolStandardize/testTautomer.cpp
# ---------------------------------------------------------------------------

rdkit_canon_cases = [
    # (input, rdkit_expected, category)
    # keto-enol
    ("C1(=CCCCC1)O",                    "O=C1CCCCC1",                     "keto_enol"),
    ("C1(CCCCC1)=O",                    "O=C1CCCCC1",                     "keto_enol"),
    ("C(=C)(O)C1=CC=CC=C1",             "CC(=O)c1ccccc1",                 "keto_enol"),
    ("CC(C)=O",                         "CC(C)=O",                        "keto_enol"),
    ("OC(C)=C(C)C",                     "CC(=O)C(C)C",                    "keto_enol"),
    ("c1(ccccc1)CC(=O)C",               "CC(=O)Cc1ccccc1",                "keto_enol"),
    ("C1(C=CCCC1)=O",                   "O=C1C=CCCC1",                    "keto_enol"),
    # imine-enamine
    ("C1(CCCCC1)=N",                    "N=C1CCCCC1",                     "imine_enamine"),
    ("C1(=CCCCC1)N",                    "N=C1CCCCC1",                     "imine_enamine"),
    ("C1(C=CC=CN1)=CC",                 "CCc1ccccn1",                     "imine_enamine"),
    ("C1(=NC=CC=C1)CC",                 "CCc1ccccn1",                     "imine_enamine"),
    # lactam-lactim
    ("O=c1cccc[nH]1",                   "O=c1cccc[nH]1",                  "lactam_lactim"),
    ("Oc1ccccn1",                       "O=c1cccc[nH]1",                  "lactam_lactim"),
    ("Oc1ncc[nH]1",                     "O=c1[nH]cc[nH]1",               "lactam_lactim"),
    ("OC(C)=NC",                        "CNC(C)=O",                       "lactam_lactim"),
    ("CNC(C)=O",                        "CNC(C)=O",                       "lactam_lactim"),
    ("Oc1ccncc1",                       "O=c1cc[nH]cc1",                  "lactam_lactim"),
    ("Oc1ncncc1",                       "O=c1cc[nH]cn1",                  "lactam_lactim"),
    ("Oc1c(cccc3)c3nc2ccncc12",         "O=c1c2ccccc2[nH]c2ccncc12",     "lactam_lactim"),
    ("C2(=C1C(=NC=N1)[NH]C(=N2)N)O",   "Nc1nc(=O)c2[nH]cnc2[nH]1",     "lactam_lactim"),
    ("C2(C1=C([NH]C=N1)[NH]C(=N2)N)=O","Nc1nc(=O)c2[nH]cnc2[nH]1",     "lactam_lactim"),
    ("O=c1nc2[nH]ccn2cc1",             "O=c1ccn2cc[nH]c2n1",            "lactam_lactim"),
    ("c1cc(=O)[nH]c2nccn12",           "O=c1ccn2cc[nH]c2n1",            "lactam_lactim"),
    ("c1cnc2ccc[nH]c12",               "c1cnc2cc[nH]c2c1",               "lactam_lactim"),
    ("C1=CC=C(O1)O",                   "Oc1ccco1",                       "lactam_lactim"),
    ("O=C1CC=CO1",                     "Oc1ccco1",                       "lactam_lactim"),
    ("Oc1nccc2cc[nH]c(=N)c12",         "Nc1nccc2cc[nH]c(=O)c12",        "lactam_lactim"),
    # amide/thioamide
    ("S=C(N)N",                         "NC(N)=S",                        "thioamide"),
    ("SC(N)=N",                         "NC(N)=S",                        "thioamide"),
    # N-heteroaromatic tautomers
    ("N=c1[nH]ccn(C)1",                "Cn1ccnc1N",                      "n_heteroarom"),
    ("CN=c1[nH]cncc1",                 "CNc1ccncn1",                     "n_heteroarom"),
    ("Cc1n[nH]c2ncnn12",               "Cc1n[nH]c2ncnn12",               "n_heteroarom"),
    ("Cc1nnc2nc[nH]n12",               "Cc1n[nH]c2ncnn12",               "n_heteroarom"),
    ("Oc1cccc2ccncc12",                "Oc1cccc2ccncc12",                 "n_heteroarom"),
    ("O=c1cccc2cc[nH]cc1-2",           "Oc1cccc2ccncc12",                "n_heteroarom"),
    ("Oc1n(C)ncc1",                    "Cn1[nH]ccc1=O",                  "n_heteroarom"),
    ("N=c1nc[nH]cc1",                  "Nc1ccncn1",                      "n_heteroarom"),
    ("N=c(c1)ccn2cc[nH]c12",          "Nc1ccn2ccnc2c1",                 "n_heteroarom"),
    ("CN=c1nc[nH]cc1",                 "CNc1ccncn1",                     "n_heteroarom"),
    ("c1ccc2[nH]c(-c3nc4ccccc4[nH]3)nc2c1", "c1ccc2[nH]c(-c3nc4ccccc4[nH]3)nc2c1", "n_heteroarom"),
    ("c1ccc2c(c1)NC(=C1N=c3ccccc3=N1)N2",   "c1ccc2[nH]c(-c3nc4ccccc4[nH]3)nc2c1", "n_heteroarom"),
    ("CNc1ccnc2ncnn21",                "CNc1ccnc2ncnn12",                "n_heteroarom"),
    ("CN=c1ccnc2nc[nH]n21",            "CNc1ccnc2ncnn12",                "n_heteroarom"),
    ("n1ccc2ccc[nH]c12",               "c1cnc2[nH]ccc2c1",               "n_heteroarom"),
    ("c1cnc2c[nH]ccc12",               "c1cc2cc[nH]c2cn1",               "n_heteroarom"),
    ("n1ccc2c[nH]ccc12",               "c1cc2[nH]ccc2cn1",               "n_heteroarom"),
    # special: p-quinol-like
    ("Nc1ccc(C=C2C=CC(=O)C=C2)cc1",   "Nc1ccc(C=C2C=CC(=O)C=C2)cc1",   "quinol"),
    ("N=C1C=CC(=Cc2ccc(O)cc2)C=C1",   "Nc1ccc(C=C2C=CC(=O)C=C2)cc1",   "quinol"),
    # nitroso-oxime
    ("CC(C)=NO",                       "CC(C)=NO",                       "nitroso_oxime"),
    ("CC(C)N=O",                       "CC(C)=NO",                       "nitroso_oxime"),
    ("O=Nc1ccc(O)cc1",                 "O=Nc1ccc(O)cc1",                 "nitroso_oxime"),
    ("O=C1C=CC(=NO)C=C1",              "O=Nc1ccc(O)cc1",                 "nitroso_oxime"),
    # cyanic acid / isocyanate
    ("C(#N)O",                         "N=C=O",                          "cyanic"),
    ("C(=N)=O",                        "N=C=O",                          "cyanic"),
    ("C#N",                            "C#N",                            "cyanic"),
    ("[C-]#[NH+]",                     "C#N",                            "cyanic"),
    # phosphorous acid
    ("[PH](=O)(O)(O)",                 "O=[PH](O)O",                     "phosphorous"),
    ("P(O)(O)O",                       "O=[PH](O)O",                     "phosphorous"),
    # nitro
    ("C([N+](=O)[O-])C",               "CC[N+](=O)[O-]",                 "nitro"),
    ("C(=[N+](O)[O-])C",               "CC[N+](=O)[O-]",                 "nitro"),
    # ketene
    ("CC=C=O",                         "CC=C=O",                         "ketene"),
    ("CC#CO",                          "CC=C=O",                         "ketene"),
    # amidine (github3755)
    ("NC(=N)C(N)CO",                   "N=C(N)C(N)CO",                   "amidine_rdkit"),
    ("NC(=N)NC(N)CO",                  "N=C(N)NC(N)CO",                  "amidine_rdkit"),
    # amino acid (github3755)
    ("OC(=O)C(N)CO",                   "NC(CO)C(=O)O",                   "amino_acid"),
    ("C([C@@H](C(=O)O)N)O",            "NC(CO)C(=O)O",                   "amino_acid"),
    ("OC(=O)C(N)CN",                   "NCC(N)C(=O)O",                   "amino_acid"),
    ("NC(=O)C(N)CO",                   "NC(=O)C(N)CO",                   "amino_acid"),
]

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def chython_prepare(smi):
    mol = chython_smiles(smi)
    mol.kekule()
    mol.thiele()
    return mol


def chython_canonicalize(smi):
    mol = chython_prepare(smi)
    mol.canonicalize()
    return mol


def chython_standardize_isomers(smi):
    mol = chython_prepare(smi)
    mol.standardize_isomers()
    return mol


_TE = rdMolStandardize.TautomerEnumerator()


def rdkit_to_chython(smi):
    """Apply RDKit tautomer canonicalization and return a MoleculeContainer (merged fragments)."""
    parts = smi.split('.')
    mols = []
    for p in parts:
        m = Chem.MolFromSmiles(p)
        if m is None:
            return None
        canonical = _TE.Canonicalize(m)
        Chem.SanitizeMol(canonical)
        mols.append(MoleculeContainer.from_rdkit(canonical))
    if len(mols) == 1:
        return mols[0]
    result = mols[0].copy()
    for m in mols[1:]:
        result = result | m
    return result


# ---------------------------------------------------------------------------
# Run comparison
# ---------------------------------------------------------------------------

results = []

for inp, expected, category in ALL_CASES:
    is_isomer = category in ISOMER_CATEGORIES

    # expected reference
    if is_isomer:
        expected_mol = chython_prepare(expected)  # kekule+thiele only (test_isomers.py style)
    else:
        expected_mol = chython_canonicalize(expected)

    # chython
    try:
        if is_isomer:
            chython_mol = chython_standardize_isomers(inp)
        else:
            chython_mol = chython_canonicalize(inp)
        chython_ok = (chython_mol == expected_mol)
        chython_out = str(chython_mol)
    except Exception as e:
        chython_ok = False
        chython_out = f'ERROR: {e}'

    # rdkit: convert with from_rdkit, then apply same chython post-processing.
    # For amidine: skip standardize_isomers (atom-numbering-dependent; compare the raw form).
    rdkit_raw = None
    try:
        rdkit_mol = rdkit_to_chython(inp)
        if rdkit_mol is None:
            rdkit_ok = False
            rdkit_out = 'PARSE_ERROR'
        else:
            rdkit_raw = str(rdkit_mol)
            if is_isomer and category != 'amidine':
                rdkit_mol.standardize_isomers()
            elif not is_isomer:
                rdkit_mol.canonicalize()
            rdkit_ok = (rdkit_mol == expected_mol)
            rdkit_out = rdkit_raw
    except Exception as e:
        rdkit_ok = False
        rdkit_out = rdkit_raw or f'ERROR: {e}'

    results.append({
        'input': inp,
        'expected': expected,
        'category': category,
        'chython_ok': chython_ok,
        'rdkit_ok': rdkit_ok,
        'chython_out': chython_out,
        'rdkit_out': rdkit_out,
    })

# ---------------------------------------------------------------------------
# Section 2: RDKit's own test suite — test chython against rdkit ground truth
# ---------------------------------------------------------------------------

results2 = []

for inp, rdkit_expected, category in rdkit_canon_cases:
    # RDKit ground truth: run rdkit, compare SMILES directly (rdkit->rdkit should always pass)
    try:
        rdkit_mol_raw = rdkit_to_chython(inp)
        if rdkit_mol_raw is None:
            rdkit_ok2 = False
            rdkit_out2 = 'PARSE_ERROR'
        else:
            # re-canonicalize expected through rdkit too, to get normalized form
            rdkit_expected_mol = rdkit_to_chython(rdkit_expected)
            if rdkit_expected_mol is None:
                rdkit_ok2 = False
                rdkit_out2 = str(rdkit_mol_raw)
            else:
                rdkit_ok2 = (rdkit_mol_raw == rdkit_expected_mol)
                rdkit_out2 = str(rdkit_mol_raw)
    except Exception as e:
        rdkit_ok2 = False
        rdkit_out2 = f'ERROR: {e}'

    # Chython: full canonicalize, compare against rdkit's expected (via chython)
    try:
        chython_mol2 = chython_canonicalize(inp)
        chython_out2 = str(chython_mol2)
        if rdkit_expected_mol is not None:
            rdkit_expected_canon = chython_canonicalize(rdkit_expected)
            chython_ok2 = (chython_mol2 == rdkit_expected_canon)
        else:
            chython_ok2 = False
    except Exception as e:
        chython_ok2 = False
        chython_out2 = f'ERROR: {e}'

    results2.append({
        'input': inp,
        'expected': rdkit_expected,
        'category': category,
        'chython_ok': chython_ok2,
        'rdkit_ok': rdkit_ok2,
        'chython_out': chython_out2,
        'rdkit_out': rdkit_out2,
    })

# ---------------------------------------------------------------------------
# Speed benchmark (combined set)
# ---------------------------------------------------------------------------

ALL_SPEED = list(ALL_CASES) + [(inp, exp, cat) for inp, exp, cat in rdkit_canon_cases]

def bench_chython():
    for inp, _, category in ALL_SPEED:
        try:
            if category in ISOMER_CATEGORIES:
                chython_standardize_isomers(inp)
            else:
                chython_canonicalize(inp)
        except Exception:
            pass


def bench_rdkit():
    for inp, _, _ in ALL_SPEED:
        try:
            for p in inp.split('.'):
                m = Chem.MolFromSmiles(p)
                if m:
                    _TE.Canonicalize(m)
        except Exception:
            pass


N = 30
n1 = len(ALL_CASES)
n2 = len(rdkit_canon_cases)
t_chython = timeit.timeit(bench_chython, number=N) / N
t_rdkit   = timeit.timeit(bench_rdkit,   number=N) / N
chython_speed = (n1 + n2) / t_chython
rdkit_speed   = (n1 + n2) / t_rdkit

# ---------------------------------------------------------------------------
# Report helpers
# ---------------------------------------------------------------------------

from collections import defaultdict


def print_section(title, results, note=''):
    cat_stats = defaultdict(lambda: {'chython': 0, 'rdkit': 0, 'total': 0})
    for r in results:
        cat_stats[r['category']]['total'] += 1
        if r['chython_ok']:
            cat_stats[r['category']]['chython'] += 1
        if r['rdkit_ok']:
            cat_stats[r['category']]['rdkit'] += 1

    total = len(results)
    chython_total = sum(1 for r in results if r['chython_ok'])
    rdkit_total   = sum(1 for r in results if r['rdkit_ok'])
    only_chython  = sum(1 for r in results if r['chython_ok'] and not r['rdkit_ok'])
    only_rdkit    = sum(1 for r in results if not r['chython_ok'] and r['rdkit_ok'])
    neither       = sum(1 for r in results if not r['chython_ok'] and not r['rdkit_ok'])

    print(f"\n{'=' * 72}")
    print(f"{title}")
    if note:
        print(f"({note})")
    print(f"{'=' * 72}")
    print(f"\nOVERALL  ({total} cases)")
    print(f"  Chython correct: {chython_total}/{total}  ({100*chython_total/total:.1f}%)")
    print(f"  RDKit correct:   {rdkit_total}/{total}  ({100*rdkit_total/total:.1f}%)")
    print(f"  Chython only:    {only_chython}")
    print(f"  RDKit only:      {only_rdkit}")
    print(f"  Neither:         {neither}")

    print(f"\nBY CATEGORY")
    print(f"  {'Category':<22} {'Chython':>9} {'RDKit':>9} {'n':>4}")
    print(f"  {'-'*22} {'-'*9} {'-'*9} {'-'*4}")
    for cat in sorted(cat_stats):
        s = cat_stats[cat]
        t = s['total']
        print(f"  {cat:<22} {s['chython']:>4}/{t:<4} {s['rdkit']:>4}/{t:<4} {t:>4}")

    chython_wins = [r for r in results if r['chython_ok'] and not r['rdkit_ok']]
    rdkit_wins   = [r for r in results if r['rdkit_ok'] and not r['chython_ok']]
    neither_list = [r for r in results if not r['chython_ok'] and not r['rdkit_ok']]

    if chython_wins:
        print(f"\nCHYTHON WINS ({len(chython_wins)})")
        for r in chython_wins:
            print(f"  [{r['category']}]  {r['input']}")
            print(f"    expected:  {r['expected']}")
            print(f"    RDKit out: {r['rdkit_out']}")

    if rdkit_wins:
        print(f"\nRDKIT WINS ({len(rdkit_wins)})")
        for r in rdkit_wins:
            print(f"  [{r['category']}]  {r['input']}")
            print(f"    expected:    {r['expected']}")
            print(f"    Chython out: {r['chython_out']}")

    if neither_list:
        print(f"\nNEITHER CORRECT ({len(neither_list)})")
        for r in neither_list:
            print(f"  [{r['category']}]  {r['input']}  -> expected: {r['expected']}")
            print(f"    Chython: {r['chython_out']}")
            print(f"    RDKit:   {r['rdkit_out']}")


# ---------------------------------------------------------------------------
# Print results
# ---------------------------------------------------------------------------

print(f"{'=' * 72}")
print(f"TAUTOMER CANONICALIZATION BENCHMARK")
print(f"RDKit {Chem.rdBase.rdkitVersion}  vs  chython")

print_section(
    "SECTION 1: CHYTHON TEST SUITE  (chython's ground truth)",
    results,
    "test_isomers.py + test_groups.py — chython's own expected canonical forms"
)

print_section(
    "SECTION 2: RDKIT TEST SUITE  (RDKit's ground truth)",
    results2,
    "canonTautomerData from testTautomer.cpp — RDKit's expected canonical forms"
)

print(f"\n{'=' * 72}")
print(f"SPEED  ({n1+n2} mol/pass, {N} passes each)")
print(f"  Chython: {chython_speed:>7.0f} mol/s")
print(f"  RDKit:   {rdkit_speed:>7.0f} mol/s")
if chython_speed > rdkit_speed:
    print(f"  -> Chython is {chython_speed/rdkit_speed:.1f}x faster")
else:
    print(f"  -> RDKit is {rdkit_speed/chython_speed:.1f}x faster")
print("=" * 72)
