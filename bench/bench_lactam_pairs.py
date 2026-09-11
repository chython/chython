"""Does chython collapse a pair of drawings that differ ONLY in which ring N holds the hydrogen?

That is the question the pyridone preference reframes.  `standardize_isomers` places a mobile ring
hydrogen, but only where both its bonds are aromatic.  chython canonicalises hydroxy-azines to the
oxo form and stores that ring NON-aromatic, so for the whole lactam family the placement pass never
sees a site -- the drawings must already agree, or the keys differ and a registry holds the same
compound twice.

Pairs are not hand-written.  RDKit's `TautomerEnumerator.Enumerate` preserves atom indices, so for
each input the enumerated forms can be filtered down to exactly those that differ from it at RING
NITROGEN ONLY -- no oxygen, no carbon, no exocyclic amine.  That is the mobile-ring-N-H class and
nothing else.  RDKit is used as a drawing GENERATOR here, never as an oracle for the direction; the
question asked of each toolkit is only whether it puts its own drawings on one key.

Each group is then labelled by what chython's canonical form actually looks like:

  AROMATIC   >= 2 sites `standardize_isomers` can act on   -> the pass owns it
  LACTAM     the mobile N sit on non-aromatic ring bonds   -> nothing in V3 places it
  MIXED      both shapes present in the group

Usage:  python bench_lactam_pairs.py [n_corpus]
"""
import sys
from collections import Counter, defaultdict
from pathlib import Path

from rdkit import Chem, RDConfig, RDLogger
RDLogger.DisableLog('rdApp.*')
from rdkit.Chem.MolStandardize import rdMolStandardize

from chython import smiles as chython_smiles

LIMIT = int(sys.argv[1]) if len(sys.argv) > 1 else 1500
SRC = Path(RDConfig.RDDataDir) / 'NCI' / 'first_5K.smi'
MOBILE = (7, 15, 33)
_TE = rdMolStandardize.TautomerEnumerator()

#: Public compounds of the family the preference concerns -- nucleobases and simple azinones,
#: every one of them in any pharmacopoeia or textbook.
CURATED = [
    ('2-pyridone', 'O=C1NC=CC=C1'),
    ('4-pyridone', 'O=C1C=CNC=C1'),
    ('uracil', 'O=C1NC(=O)C=CN1'),
    ('thymine', 'CC1=CNC(=O)NC1=O'),
    ('cytosine', 'NC1=NC(=O)NC=C1'),
    ('isocytosine', 'NC1=NC=CC(=O)N1'),
    ('guanine', 'NC1=NC2=C(N=CN2)C(=O)N1'),
    ('hypoxanthine', 'O=C1NC=NC2=C1NC=N2'),
    ('xanthine', 'O=C1NC(=O)C2=C(N1)NC=N2'),
    ('allopurinol', 'O=C1NC=NC2=C1C=NN2'),
    ('4-quinazolinone', 'O=C1NC=NC2=CC=CC=C12'),
    ('2-quinoxalinone', 'O=C1CN=C2C=CC=CC2=N1'),
    ('4-pyrimidinone', 'O=C1C=CN=CN1'),
    ('1,2,4-triazol-3-one', 'O=C1NN=CN1'),
    ('pyrazol-3-one', 'O=C1C=CNN1'),
    ('maleic hydrazide', 'O=C1C=CC(=O)NN1'),
    ('barbituric acid', 'O=C1CC(=O)NC(=O)N1'),
    ('purin-6-one', 'O=C1NC=NC2=C1NC=N2'),
    ('2-thiouracil', 'S=C1NC(=O)C=CN1'),
    ('cyanuric acid', 'O=C1NC(=O)NC(=O)N1'),
    ('phthalazin-1-one', 'O=C1NN=CC2=CC=CC=C12'),
    ('quinazoline-2,4-dione', 'O=C1NC(=O)C2=CC=CC=C2N1'),
    ('5-azacytosine', 'NC1=NC(=O)NN=C1'),
    ('imidazol-2-one', 'O=C1NC=CN1'),
    ('1,3,5-triazin-2-one', 'O=C1NC=NC=N1'),
]


def aromatic_sites(mol):
    """Exactly `_isomers._sites`: what `standardize_isomers` can act on."""
    out = []
    for n in mol.atoms_numbers:
        if mol.element_of(n) not in MOBILE or mol.radical_of(n):
            continue
        if mol.charge_of(n) not in (0, -1):
            continue
        h = mol.implicit_h_of(n)
        if h is None or h > 1:
            continue
        nb = tuple(mol.neighbors_of(n))
        if len(nb) != 2 or any(mol.order_of(n, m) != 4 for m in nb):
            continue
        out.append(n)
    return out


def lactam_mobile(mol):
    """Ring N of degree 2 on non-aromatic bonds, in a ring system bearing an exocyclic C=O/C=S/C=N,
    where the hydrogens are unevenly distributed -- so another placement exists."""
    rings = [set(r) for r in mol.sssr]
    if not rings:
        return []
    systems = []
    for r in rings:
        hit = [s for s in systems if s & r]
        if hit:
            merged = set(r)
            for s in hit:
                merged |= s
                systems.remove(s)
            systems.append(merged)
        else:
            systems.append(set(r))
    groups = []
    for system in systems:
        cand, carbonyl = [], False
        for n in system:
            for m in mol.neighbors_of(n):
                if m not in system and mol.order_of(n, m) == 2 \
                        and mol.element_of(m) in (7, 8, 16):
                    carbonyl = True
            if mol.element_of(n) not in MOBILE or mol.radical_of(n) or mol.charge_of(n):
                continue
            nb = tuple(mol.neighbors_of(n))
            if len(nb) != 2 or any(mol.order_of(n, m) == 4 for m in nb):
                continue
            cand.append((n, bool(mol.implicit_h_of(n))))
        if carbonyl and len(cand) >= 2:
            k = sum(h for _, h in cand)
            if 0 < k < len(cand):
                groups.append(cand)
    return groups


def ring_n_shifts(smi):
    """Every enumerated tautomer that differs from `smi` at RING NITROGEN ONLY.

    Atom indices are preserved by `Enumerate`, so the difference is read atom by atom.  A change on
    oxygen, carbon or an exocyclic nitrogen disqualifies the form: that is a different tautomerism
    and would confound the measurement.
    """
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None, []
    base = [a.GetTotalNumHs() for a in m.GetAtoms()]
    ok_idx = {a.GetIdx() for a in m.GetAtoms()
              if a.GetAtomicNum() == 7 and a.IsInRing() and a.GetDegree() == 2}
    out = []
    for t in _TE.Enumerate(m):
        if t.GetNumAtoms() != m.GetNumAtoms():
            continue
        diff = [a.GetIdx() for a in t.GetAtoms() if a.GetTotalNumHs() != base[a.GetIdx()]]
        if not diff or not set(diff) <= ok_idx:
            continue
        if sum(t.GetAtomWithIdx(i).GetTotalNumHs() - base[i] for i in diff):
            continue                     # net hydrogen must be conserved on the ring nitrogens
        out.append(Chem.MolToSmiles(t))
    return m, sorted(set(out))


def chython_state(smi):
    """(key, canonical smiles, n aromatic sites, n lactam groups) or an error string."""
    try:
        mol = chython_smiles(smi, log=[])
        mol.kekule()
        mol.canonicalize()
    except Exception as e:
        return f'!{type(e).__name__}: {e}'
    return (mol.canonical_bytes, str(mol), len(aromatic_sites(mol)), len(lactam_mobile(mol)))


def rdkit_key(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return '!refused'
    return Chem.MolToSmiles(_TE.Canonicalize(m))


def measure(label, smi):
    """One group: the input plus every ring-N-H shift of it."""
    m, shifts = ring_n_shifts(smi)
    if m is None or not shifts:
        return None
    members = [Chem.MolToSmiles(m)] + shifts
    states = [chython_state(s) for s in members]
    if any(isinstance(s, str) for s in states):
        return {'label': label, 'members': members, 'error': [s for s in states
                                                              if isinstance(s, str)]}
    ck = {s[0] for s in states}
    rk = {rdkit_key(s) for s in members}
    arom = any(s[2] >= 2 for s in states)
    lact = any(s[3] for s in states)
    shape = 'MIXED' if arom and lact else 'AROMATIC' if arom else 'LACTAM' if lact else 'NEITHER'
    return {'label': label, 'members': members, 'shape': shape,
            'chython_collapsed': len(ck) == 1, 'rdkit_collapsed': len(rk) == 1,
            'out': [s[1] for s in states], 'n': len(members)}


# ---------------------------------------------------------------------------
print('=' * 78)
print('MOBILE RING-N-H COLLAPSE -- drawings that differ only in which ring N carries the H')
print(f'chython (this tree) vs RDKit {Chem.rdBase.rdkitVersion}')
print('=' * 78)

for title, cases in (('CURATED public azinones / nucleobases', CURATED),
                     (f'NCI first_5K, first {LIMIT}', None)):
    if cases is None:
        lines = [l.split()[0] for l in SRC.read_text().splitlines() if l.strip()][:LIMIT]
        cases = [(f'nci:{i}', s) for i, s in enumerate(lines)]
    results = []
    for label, smi in cases:
        r = measure(label, smi)
        if r is not None:
            results.append(r)

    bad = [r for r in results if 'error' in r]
    results = [r for r in results if 'error' not in r]
    by_shape = defaultdict(lambda: [0, 0, 0])
    for r in results:
        s = by_shape[r['shape']]
        s[2] += 1
        s[0] += r['chython_collapsed']
        s[1] += r['rdkit_collapsed']

    print(f'\n--- {title}: {len(results)} groups with a ring-N-H shift'
          + (f' ({len(bad)} unreadable)' if bad else ''))
    print(f"    {'shape':<10} {'chython':>10} {'rdkit':>10}   {'drawings':>9}")
    for shape in sorted(by_shape):
        c, k, t = by_shape[shape]
        n = sum(r['n'] for r in results if r['shape'] == shape)
        print(f'    {shape:<10} {c:>5}/{t:<4} {k:>5}/{t:<4}   {n:>9}')
    c = sum(r['chython_collapsed'] for r in results)
    k = sum(r['rdkit_collapsed'] for r in results)
    print(f"    {'TOTAL':<10} {c:>5}/{len(results):<4} {k:>5}/{len(results):<4}"
          f"   {sum(r['n'] for r in results):>9}")

    show = [r for r in results if not r['chython_collapsed']]
    if cases and len(cases) <= 40:
        print(f'\n    every group chython does not collapse ({len(show)}):')
        for r in show[:40]:
            print(f"      [{r['shape']}] {r['label']}"
                  f"  rdkit={'collapsed' if r['rdkit_collapsed'] else 'split'}")
            for smi, out in zip(r['members'], r['out']):
                print(f'         {smi:<44} -> {out}')
    else:
        print(f'\n    chython does not collapse {len(show)}; first 8:')
        for r in show[:8]:
            print(f"      [{r['shape']}] {r['label']}"
                  f"  rdkit={'collapsed' if r['rdkit_collapsed'] else 'split'}")
            for smi, out in zip(r['members'], r['out']):
                print(f'         {smi:<44} -> {out}')
    if bad:
        print('\n    unreadable:')
        for r in bad[:5]:
            print(f"      {r['label']}: {r['error'][0][:90]}")
print('=' * 78)
