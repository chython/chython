"""Deduplication harness: does each toolkit collapse a tautomer set to ONE key?

The metric is not "does the output match an expected string" -- that measures agreement with a
chosen spelling.  For deduplication the only question is whether every member of a set of
drawings of one compound lands on the same key, and whether two different compounds do not.

Each toolkit is measured in its own currency: chython by `canonicalize()` then `str()`
(a canonical SMILES), RDKit by `TautomerEnumerator().Canonicalize()` then `MolToSmiles`.
No interop conversion is involved, so neither is scored on the other's writer.

Usage:  python bench_tautomer_dedup.py
"""
import timeit
from collections import defaultdict

from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')
from rdkit.Chem.MolStandardize import rdMolStandardize

from chython import smiles as chython_smiles


# ---------------------------------------------------------------------------
# Corpus.  Every group is one compound written several ways.  `class` names the
# tautomerism involved; `src` says who asserts the group is one compound.
# ---------------------------------------------------------------------------

GROUPS = [
    # ---- annular (ring N-H / ring charge shift) -- what standardize_isomers is for
    ('annular', 'pyrazole-3-Me',        ['CC1=NNC=C1', 'CC1=CC=NN1', 'Cc1cc[nH]n1', 'Cc1ccn[nH]1']),
    ('annular', 'pyrazole-4-Me',        ['Cc1c[nH]nc1', 'Cc1cn[nH]c1']),
    ('annular', 'imidazole-4-Me',       ['CC1=CN=CN1', 'CC1=CNC=N1', 'Cc1cnc[nH]1', 'Cc1c[nH]cn1']),
    ('annular', '1,2,3-triazole-4-Me',  ['CC1=CN=NN1', 'CC1=CNN=N1', 'CC1=NNN=C1']),
    ('annular', '1,2,4-triazole-3-Me',  ['CC1=NC=NN1', 'CC1=NN=CN1']),
    ('annular', 'tetrazole-5-Me',       ['CC1=NNN=N1', 'CC1=NN=NN1']),
    ('annular', '1,2,4-triazole',       ['N1C=NC=N1', 'c1nc[nH]n1', 'c1[nH]cnn1']),
    ('annular', '1,2,3-triazole',       ['c1cn[nH]n1', 'c1c[nH]nn1']),
    ('annular', 'benzimidazole-5-Me',   ['Cc1ccc2[nH]cnc2c1', 'Cc1ccc2nc[nH]c2c1']),
    ('annular', 'indazole-5-Cl',        ['Clc1ccc2[nH]ncc2c1', 'Clc1ccc2n[nH]cc2c1']),
    ('annular', 'purine',               ['c1ncc2[nH]cnc2n1', 'c1ncc2nc[nH]c2n1']),
    ('annular', 'adenine',              ['Nc1ncnc2[nH]cnc12', 'Nc1ncnc2nc[nH]c12']),
    ('annular', 'pyrazolo[3,4-b]pyr',   ['Cc1n[nH]c2ncccc12', 'Cc1[nH]nc2ncccc12']),
    ('annular', '4-Me-imidazolium',     ['Cc1c[nH]c[nH+]1', 'Cc1c[nH+]c[nH]1']),
    ('annular', 'pyrazol-3-olate',      ['Cc1cc[n-]n1', 'Cc1ccn[n-]1']),
    ('annular', '4-nitroimidazole',     ['[O-][N+](=O)c1cnc[nH]1', '[O-][N+](=O)c1c[nH]cn1']),
    ('annular', 'bis-imidazole',        ['c1c[nH]cn1.Cc1c[nH]cn1', 'c1cnc[nH]1.Cc1cnc[nH]1']),
    ('annular', '8-fused-pyrazole-x8',  ['c1cc[nH]n1.' * 8, 'c1ccn[nH]1.' * 8]),

    # ---- lactam / lactim (2-pyridone family) -- standardize's SMARTS rules
    ('lactam', '2-pyridone',            ['Oc1ccccn1', 'O=c1cccc[nH]1', 'OC1=CC=CC=N1']),
    ('lactam', '4-pyridone',            ['Oc1ccncc1', 'O=c1cc[nH]cc1', 'OC1=CC=NC=C1']),
    ('lactam', '2-hydroxypyrimidine',   ['Oc1ncccn1', 'O=c1[nH]cccn1', 'OC1=NC=CC=N1']),
    ('lactam', 'uracil',                ['Oc1cc[nH]c(=O)n1', 'O=c1cc[nH]c(=O)[nH]1', 'OC1=CC=NC(=O)N1']),
    ('lactam', '2-hydroxyimidazole',    ['Oc1ncc[nH]1', 'O=c1[nH]cc[nH]1']),
    ('lactam', 'N-methylacetamide',     ['OC(C)=NC', 'CNC(C)=O']),
    ('lactam', 'formamide',             ['N=CO', 'NC=O']),
    ('lactam', 'thiourea',              ['S=C(N)N', 'SC(N)=N']),
    ('lactam', 'thioformamide',         ['N=CS', 'NC=S']),
    ('lactam', '4-quinolone',           ['Oc1ccnc2ccccc12', 'O=c1cc[nH]c2ccccc12']),
    ('lactam', 'guanine',               ['Nc1nc(O)c2[nH]cnc2n1', 'Nc1nc(=O)c2[nH]cnc2[nH]1']),
    ('lactam', 'cytosine',              ['N=C1NC=CC(=O)N1', 'NC1=NC=CC(=O)N1']),

    # ---- keto/enol
    ('keto_enol', 'cyclohexanone',      ['C1(=CCCCC1)O', 'O=C1CCCCC1']),
    ('keto_enol', 'acetophenone',       ['C(=C)(O)C1=CC=CC=C1', 'CC(=O)c1ccccc1']),
    ('keto_enol', 'acetaldehyde',       ['OC=C', 'O=CC']),
    ('keto_enol', 'acetone',            ['OC(C)=C', 'O=C(C)C']),
    ('keto_enol', 'MIBK-enol',          ['OC(C)=C(C)C', 'CC(=O)C(C)C']),
    ('keto_enol', 'cyclohex-2-enone',   ['C1(C=CCCC1)=O', 'OC1=CC=CCC1']),

    # ---- imine/enamine
    ('imine', 'cyclohexanimine',        ['C1(CCCCC1)=N', 'C1(=CCCCC1)N']),
    ('imine', '2-ethylpyridine',        ['C1(C=CC=CN1)=CC', 'C1(=NC=CC=C1)CC', 'CCc1ccccn1']),
    ('imine', '2-aminopyrimidine',      ['N=c1nc[nH]cc1', 'Nc1ccncn1']),
    ('imine', '2-methylaminopyrimidine', ['CN=c1[nH]cncc1', 'CNc1ccncn1']),

    # ---- amidine / guanidine
    ('amidine', 'O-Me-N-Me-isourea',    ['COC(=N)NC', 'COC(N)=NC']),
    ('amidine', 'N,N-diethylguanidine', ['CCN=C(N)NC', 'CCNC(=NC)N', 'CCNC(N)=NC']),
    ('amidine', 'biguanide-ish',        ['CNC(N)=NC(=N)NC', 'CNC(=N)NC(=N)NC']),

    # ---- classes chython does not claim (RDKit's enumerator does)
    ('nitroso_oxime', 'acetoxime',      ['CC(C)=NO', 'CC(C)N=O']),
    ('nitroso_oxime', 'p-nitrosophenol', ['O=Nc1ccc(O)cc1', 'O=C1C=CC(=NO)C=C1']),
    ('nitro', 'nitroethane',            ['C([N+](=O)[O-])C', 'C(=[N+](O)[O-])C']),
    ('cyanic', 'cyanic acid',           ['C(#N)O', 'C(=N)=O']),
    ('phosphorous', 'phosphorous acid', ['[PH](=O)(O)(O)', 'P(O)(O)O']),
    ('ketene', 'ketene',                ['CC=C=O', 'CC#CO']),
    ('ring_chain', 'glucose',           ['OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O',
                                         'OC[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)C=O']),
]

# Pairs that are DIFFERENT compounds and must land on different keys.  An over-merge is worse
# than an under-merge for a registry: it silently loses a compound.
MUST_NOT_MERGE = [
    ('constitution', '2- vs 3-pyridone',        'O=c1cccc[nH]1',      'Oc1cccnc1'),
    ('constitution', '4-Me- vs 5-Me-imid',      'Cc1cnc[nH]1',        'Cc1ncc[nH]1'),
    ('N-substituted', '1-Me vs 2-Me triazole', 'Cn1ccnn1',           'Cn1nccn1'),
    ('N-substituted', '1,3- vs 1,5-diMe-pyraz', 'Cc1ccn(C)n1',        'Cc1cc(C)nn1'),
    ('N-substituted', '1-Me-imid vs 4-Me-imid', 'Cn1ccnc1',           'Cc1cnc[nH]1'),
    ('scaffold', 'phenol vs cyclohexadienone',  'Oc1ccccc1',          'O=C1CC=CC=C1'),
    ('scaffold', 'aniline vs cyclohexadienimine', 'Nc1ccccc1',        'N=C1CC=CC=C1'),
    ('oxidation', 'pyridine vs pyridine-N-oxide', 'c1ccncc1',         '[O-][n+]1ccccc1'),
    ('tautomer-vs-isomer', 'acetamide vs Me-formamide', 'CC(N)=O',    'CNC=O'),
    ('regio', 'indazole 1H vs 2H N-Me',         'Cn1ncc2ccccc21',     'Cn1cc2ccccc2n1'),
]


# ---------------------------------------------------------------------------
# Keys
# ---------------------------------------------------------------------------

_TE = rdMolStandardize.TautomerEnumerator()


def chython_key(smi):
    """chython's dedup key: `canonical_bytes` after the documented pipeline.

    `canonical_bytes` and not `str()`: the docstring on `canonicalize()` names it, `__hash__` and
    `__eq__` are built on it, and a SMILES string is a writer's opinion about a canonical form
    rather than the form itself.
    """
    mol = chython_smiles(smi)
    mol.canonicalize()
    return mol.canonical_bytes


def rdkit_key(smi):
    """RDKit's dedup key: canonical tautomer per fragment, then canonical SMILES."""
    frags = []
    for p in smi.split('.'):
        m = Chem.MolFromSmiles(p)
        if m is None:
            raise ValueError(f'rdkit refused {p!r}')
        frags.append(Chem.MolToSmiles(_TE.Canonicalize(m)))
    return '.'.join(sorted(frags))


def _show(smi, key):
    """A human-readable stand-in for the byte key: the canonical SMILES of the same result."""
    if isinstance(key, str) and key.startswith('!'):
        return key
    try:
        m = chython_smiles(smi)
        m.canonicalize()
        return str(m)
    except Exception as e:
        return f'!{e}'


def keys_of(fn, members):
    out = []
    for smi in members:
        try:
            out.append(fn(smi))
        except Exception as e:
            out.append(f'!{type(e).__name__}: {e}')
    return out


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

rows = []
for cls, name, members in GROUPS:
    if len(members) < 2:
        continue
    ck, rk = keys_of(chython_key, members), keys_of(rdkit_key, members)
    rows.append({
        'class': cls, 'name': name, 'n': len(members), 'members': members,
        'chython_collapsed': len(set(ck)) == 1 and not any(isinstance(k, str) for k in ck),
        'rdkit_collapsed': len(set(rk)) == 1 and not any(k.startswith('!') for k in rk),
        'chython_keys': ck, 'rdkit_keys': rk,
    })

split_rows = []
for cls, name, a, b in MUST_NOT_MERGE:
    ca, cb = keys_of(chython_key, [a, b])
    ra, rb = keys_of(rdkit_key, [a, b])
    split_rows.append({
        'class': cls, 'name': name, 'a': a, 'b': b,
        'chython_kept': ca != cb, 'rdkit_kept': ra != rb,
        'chython_keys': (ca, cb), 'rdkit_keys': (ra, rb),
    })


# ---------------------------------------------------------------------------
# Speed, on the same corpus
# ---------------------------------------------------------------------------

FLAT = [s for _, _, m in GROUPS for s in m]


def run_chython():
    for s in FLAT:
        try:
            chython_key(s)
        except Exception:
            pass


def run_rdkit():
    for s in FLAT:
        try:
            rdkit_key(s)
        except Exception:
            pass


N = 20
t_c = timeit.timeit(run_chython, number=N) / N
t_r = timeit.timeit(run_rdkit, number=N) / N


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

print('=' * 78)
print(f'TAUTOMER DEDUPLICATION -- chython (this tree) vs RDKit {Chem.rdBase.rdkitVersion}')
print('=' * 78)

per_class = defaultdict(lambda: [0, 0, 0])
for r in rows:
    s = per_class[r['class']]
    s[2] += 1
    s[0] += r['chython_collapsed']
    s[1] += r['rdkit_collapsed']

print('\nCOLLAPSE: a set of drawings of one compound must give ONE key')
print(f"  {'class':<16} {'chython':>10} {'rdkit':>10}")
print(f"  {'-'*16} {'-'*10} {'-'*10}")
for cls in sorted(per_class):
    c, r, t = per_class[cls]
    print(f'  {cls:<16} {c:>5}/{t:<4} {r:>5}/{t:<4}')
tc = sum(r['chython_collapsed'] for r in rows)
tr = sum(r['rdkit_collapsed'] for r in rows)
print(f"  {'-'*16} {'-'*10} {'-'*10}")
print(f'  {"TOTAL":<16} {tc:>5}/{len(rows):<4} {tr:>5}/{len(rows):<4}')

print('\nSEPARATION: two different compounds must give TWO keys')
kc = sum(r['chython_kept'] for r in split_rows)
kr = sum(r['rdkit_kept'] for r in split_rows)
print(f'  chython kept apart: {kc}/{len(split_rows)}     rdkit kept apart: {kr}/{len(split_rows)}')
for r in split_rows:
    if not (r['chython_kept'] and r['rdkit_kept']):
        who = []
        if not r['chython_kept']:
            who.append('chython MERGED')
        if not r['rdkit_kept']:
            who.append('rdkit MERGED')
        print(f"    [{r['class']}] {r['name']}: {' / '.join(who)}")
        print(f"        {r['a']}  |  {r['b']}")
        print(f"        chython: {_show(r['a'], r['chython_keys'][0])}  |  "
              f"{_show(r['b'], r['chython_keys'][1])}")
        print(f"        rdkit:   {r['rdkit_keys'][0]}  |  {r['rdkit_keys'][1]}")

print('\nDETAIL: every group where the two toolkits differ, or where both fail')
for r in rows:
    if r['chython_collapsed'] and r['rdkit_collapsed']:
        continue
    tag = ('chython OK, rdkit split' if r['chython_collapsed'] else
           'rdkit OK, chython split' if r['rdkit_collapsed'] else 'both split')
    print(f"\n  [{r['class']}] {r['name']}  -- {tag}")
    for smi, ck, rk in zip(r['members'], r['chython_keys'], r['rdkit_keys']):
        print(f'      {smi}')
        print(f'        chython -> {_show(smi, ck)}')
        print(f'        rdkit   -> {rk}')

print('\n' + '=' * 78)
print(f'SPEED  ({len(FLAT)} structures per pass, {N} passes)')
print(f'  chython: {len(FLAT)/t_c:>8.0f} struct/s')
print(f'  rdkit:   {len(FLAT)/t_r:>8.0f} struct/s   -> chython {t_r/t_c:.1f}x')
print('=' * 78)
