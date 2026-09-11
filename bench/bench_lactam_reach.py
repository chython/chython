"""How much of a real corpus does the pyridone preference put out of `standardize_isomers`' reach?

chython canonicalises a hydroxy-azine to the pyridone, and a pyridone ring is stored NON-aromatic
(`thiele()` is a no-op on `O=C1NC=CC=C1`).  `standardize_isomers` only considers atoms whose two
ring bonds are both order 4.  So every mobile hydrogen on a lactam ring is outside its reach by
construction, however many rules are added to it.

Measured over the public NCI 5K set that ships with RDKit:

  AROMATIC MOBILE   >= 2 sites `standardize_isomers` can act on today
  LACTAM MOBILE     a non-aromatic ring carrying an exocyclic C=O/C=S/C=N and >= 2 ring N of
                    degree 2, exactly one of which holds a hydrogen -- the shape whose hydrogen
                    has somewhere else it could equally sit, and which nothing in V3 places

Usage:  python bench_lactam_reach.py [n]
"""
import sys
from collections import Counter
from pathlib import Path

from rdkit import RDConfig
from chython import smiles as chython_smiles

LIMIT = int(sys.argv[1]) if len(sys.argv) > 1 else 5000
SRC = Path(RDConfig.RDDataDir) / 'NCI' / 'first_5K.smi'
MOBILE = (7, 15, 33)


def aromatic_sites(mol):
    """Exactly `_isomers._sites`: what the pass can act on today."""
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


def lactam_sites(mol):
    """Ring N of degree 2 in a non-aromatic ring that carries an exocyclic C=O / C=S / C=N.

    Grouped per ring system so that "two sites" means two sites the same hydrogen could occupy.
    Returns a list of groups, each a list of (atom, has_hydrogen).
    """
    rings = [set(r) for r in mol.sssr]
    if not rings:
        return []
    # ring systems: rings sharing an atom
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
        # non-aromatic part only: a fused system may be half aromatic (4-quinolone)
        cand = []
        carbonyl = False
        for n in system:
            for m in mol.neighbors_of(n):
                if m in system:
                    continue
                if mol.element_of(m) in (8, 16) and mol.order_of(n, m) == 2:
                    carbonyl = True
                elif mol.element_of(m) == 7 and mol.order_of(n, m) == 2:
                    carbonyl = True
            if mol.element_of(n) not in MOBILE or mol.radical_of(n) or mol.charge_of(n):
                continue
            nb = tuple(mol.neighbors_of(n))
            if len(nb) != 2:
                continue
            if any(mol.order_of(n, m) == 4 for m in nb):
                continue            # aromatic: `standardize_isomers` already owns it
            h = mol.implicit_h_of(n)
            cand.append((n, bool(h)))
        if carbonyl and len(cand) >= 2:
            k = sum(h for _, h in cand)
            if 0 < k < len(cand):   # a choice exists; all-H or no-H is not a placement
                groups.append(cand)
    return groups


lines = [l.split()[0] for l in SRC.read_text().splitlines() if l.strip()][:LIMIT]

stat = Counter()
examples = []
for smi in lines:
    try:
        mol = chython_smiles(smi, log=[])
        mol.kekule()
        mol.canonicalize()
    except Exception as e:
        stat[f'unreadable ({type(e).__name__})'] += 1
        continue
    stat['canonicalized'] += 1
    a = aromatic_sites(mol)
    groups = lactam_sites(mol)
    arom_multi = len([1 for g in [a] if len(g) >= 2])
    if arom_multi:
        stat['has >=2 aromatic mobile sites (reachable today)'] += 1
    if groups:
        stat['has a mobile lactam N-H (unreachable)'] += 1
        if len(examples) < 15:
            examples.append((smi, str(mol), [[n for n, _ in g] for g in groups]))
    if arom_multi and groups:
        stat['both'] += 1

print('=' * 78)
print(f'MOBILE-HYDROGEN REACH over {len(lines)} public NCI compounds')
print(f'source: {SRC}')
print('=' * 78)
total = stat['canonicalized']
for k, v in stat.most_common():
    pct = f'{100 * v / total:5.1f}%' if total and k != 'canonicalized' else ''
    print(f'  {k:<52} {v:>5}  {pct}')

print('\nexamples of the unreachable shape (group = ring N the hydrogen could sit on)')
for smi, out, groups in examples:
    print(f'  in : {smi}')
    print(f'  out: {out}      groups: {groups}')
print('=' * 78)
