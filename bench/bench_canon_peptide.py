"""Timing for the core's canonical labelling on peptides. Throwaway."""
from random import Random
from time import perf_counter

from chython.core import MoleculeContainer

# side chains as (element, bonds-to-previous-in-chain) walks rooted at CA
SIDE = {
    'G': [],
    'A': [('C', 0)],
    'V': [('C', 0), ('C', 1), ('C', 1)],
    'L': [('C', 0), ('C', 1), ('C', 2), ('C', 2)],
    'S': [('C', 0), ('O', 1)],
    'T': [('C', 0), ('O', 1), ('C', 1)],
    'F': [('C', 0), ('C', 1), ('C', 2), ('C', 3), ('C', 4), ('C', 5)],  # ring closed below
    'K': [('C', 0), ('C', 1), ('C', 2), ('C', 3), ('N', 4)],
    'D': [('C', 0), ('C', 1), ('O', 2), ('O', 2)],
}


def build(seq, cyclic=False):
    m = MoleculeContainer()
    with m.edit():
        prev_c = None
        first_n = None
        for res in seq:
            n = m.add_atom('N')
            ca = m.add_atom('C')
            c = m.add_atom('C')
            o = m.add_atom('O')
            m.add_bond(n, ca, 1)
            m.add_bond(ca, c, 1)
            m.add_bond(c, o, 2)
            if prev_c is None:
                first_n = n
            else:
                m.add_bond(prev_c, n, 1)
            prev_c = c
            walk = [ca]
            for element, parent in SIDE[res]:
                sid = m.add_atom(element)
                m.add_bond(walk[parent], sid, 1)
                walk.append(sid)
            if res == 'F':                      # close the phenyl ring
                m.add_bond(walk[2], walk[6], 1)
        if cyclic:
            m.add_bond(prev_c, first_n, 1)
        else:
            m.add_bond(prev_c, m.add_atom('O'), 1)   # C-term OH
    return m


def probe(label, m, repeats=5):
    n = len(m.stable_ids)
    classes = m.atoms_order_classes
    t = perf_counter()
    for _ in range(repeats):
        m.canonical_order()
    dt = (perf_counter() - t) / repeats
    t = perf_counter()
    orbits = len(set(m.automorphism_orbits().values()))
    dto = perf_counter() - t
    print(f'{label:<34} n={n:<5} refine_classes={classes:<5} '
          f'{"DISCRETE" if classes == n else "search":<9} '
          f'canonical={dt * 1000:8.3f} ms  orbits={orbits:<4} ({dto * 1000:.3f} ms)')


rng = Random(20260901)
alphabet = 'GAVLSTFKD'

print('--- realistic linear peptides, mixed sequence ---')
for k in (10, 20, 50, 100, 200):
    seq = ''.join(rng.choice(alphabet) for _ in range(k))
    probe(f'linear {k} residues (mixed)', build(seq))

print('\n--- homo-oligomers: every residue identical ---')
for k in (10, 50, 100):
    probe(f'linear poly-Gly {k}', build('G' * k))
for k in (10, 50, 100):
    probe(f'linear poly-Phe {k}', build('F' * k))

print('\n--- adversarial: cyclic homo-peptides (Cn symmetry) ---')
for k in (5, 10, 20, 40):
    probe(f'cyclo-(Gly){k}', build('G' * k, cyclic=True))
for k in (5, 10, 20):
    probe(f'cyclo-(Phe){k}', build('F' * k, cyclic=True))

print('\n--- cyclosporine-shaped: cyclic, mixed sequence ---')
for k in (11, 20):
    seq = ''.join(rng.choice(alphabet) for _ in range(k))
    probe(f'cyclic {k} residues (mixed)', build(seq, cyclic=True))
