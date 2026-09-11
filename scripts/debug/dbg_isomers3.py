"""Prototype of `standardize_isomers`: canonical PLACEMENT of mobile hydrogens.

PYTHONPATH=. python dbg_isomers3.py

The load-bearing idea: rank the atoms on a copy whose candidate hydrogens have all been STRIPPED.
`atoms_order` on the molecule as written depends on where the hydrogens sit, so it cannot be used to
choose where they should sit -- that is the circularity chython 2 never escaped.  On the stripped
skeleton the ranks are the same numbers whichever tautomer arrived, so the choice is spelling
independent by construction.

Validity is PROVED, not scored: a placement is admissible when `kekule()` on a copy leaves nothing
unresolved.  The complete backtracking kekuliser is the oracle.
"""
import sys
from itertools import combinations

from chython.core._core import read_smiles

MOBILE = frozenset((7, 15, 33))


def candidates(mol):
    """The ring heteroatoms whose hydrogen the ring, not the string, decides.

    Neutral, non-radical, two heavy neighbours, every bond aromatic: such an atom is a pyrrole-type
    donor with a hydrogen and a pyridine-type acceptor without one, and nothing in the graph says
    which.  An atom with three neighbours has no room for a double bond either way and is not mobile.
    """
    out = []
    for n in mol.atoms_numbers:
        if mol.element_of(n) not in MOBILE or mol.charge_of(n) or mol.radical_of(n):
            continue
        nbrs = tuple(mol.neighbors_of(n))
        if len(nbrs) != 2 or any(mol.order_of(n, m) != 4 for m in nbrs):
            continue
        out.append(n)
    return out


def stripped_ranks(mol, cands):
    """`atoms_order` of a copy with every candidate hydrogen removed -- the reference frame."""
    work = mol.copy()
    with work.edit():
        for n in cands:
            work.set_hydrogens(n, 0)
    return work.atoms_order


def kekulises(mol, cands, placement):
    work = mol.copy()
    with work.edit():
        for n in cands:
            work.set_hydrogens(n, 1 if n in placement else 0)
    return not work.kekule().unresolved


def choose(mol):
    """The canonical placement, or None when nothing has to move."""
    cands = candidates(mol)
    if len(cands) < 2:
        return None
    held = frozenset(n for n in cands if mol.implicit_h_of(n))
    if not held or len(held) == len(cands):
        return None                       # nothing to distribute: no choice to make
    ranks = stripped_ranks(mol, cands)
    best = None
    for placement in combinations(cands, len(held)):
        placement = frozenset(placement)
        if not kekulises(mol, cands, placement):
            continue
        key = sorted(ranks[n] for n in placement)
        if best is None or key < best[0]:
            best = (key, placement)
    if best is None or best[1] == held:
        return None
    return cands, held, best[1]


def apply(mol):
    answer = choose(mol)
    if answer is None:
        return False
    cands, held, placement = answer
    with mol.edit():
        for n in cands:
            mol.set_hydrogens(n, 1 if n in placement else 0)
    return True


PAIRS = [('Cc1cnc[nH]1', 'Cc1c[nH]cn1', '4-methylimidazole, the two N-H forms'),
         ('c1cc2[nH]ncc2cn1', 'c1cc2n[nH]cc2cn1', 'pyrazolo[3,4-c]pyridine'),
         ('Cc1n[nH]c2nc3[nH]nc(C)c3nc12', 'Cc1[nH]nc2nc3[nH]nc(C)c3nc12',
          'bis-pyrazolo fused, two mobile hydrogens'),
         ('c1ccc2[nH]ncc2c1', 'c1ccc2n[nH]cc2c1', 'indazole'),
         ('c1cnc2[nH]ccc2c1', 'c1cnc2[nH]ccc2c1', 'pyrrolopyridine, one H, one site'),
         ('c1c[nH]cn1', 'c1cnc[nH]1', 'imidazole itself -- symmetric, already one compound')]


def main():
    bad = 0
    for a, b, label in PAIRS:
        ma, mb = read_smiles(a), read_smiles(b)
        same_before = ma.canonical_bytes == mb.canonical_bytes
        apply(ma)
        apply(mb)
        same_after = ma.canonical_bytes == mb.canonical_bytes
        flag = 'OK ' if same_after else 'BAD'
        if not same_after:
            bad += 1
        print(f'{flag} {label}\n    {a} / {b}   equal before={same_before} after={same_after}')
        print(f'    -> {ma} / {mb}')

    # idempotence and non-corruption: a second pass must find nothing
    for a, _, label in PAIRS:
        m = read_smiles(a)
        apply(m)
        first = m.canonical_bytes
        changed = apply(m)
        if changed or m.canonical_bytes != first:
            bad += 1
            print(f'NOT IDEMPOTENT {label}: {a}')
    print(f'\n{bad} failures')


if __name__ == '__main__':
    sys.exit(main())
