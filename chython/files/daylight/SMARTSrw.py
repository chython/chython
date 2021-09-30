from collections import defaultdict
from ...containers import QueryContainer
from ...containers.bonds import QueryBond


class SMARTSRead:
    def _convert_query(self, data):
        atoms = data['atoms']
        bonds = defaultdict(dict)

        for n, m, value in data['query']:
            bonds[n][m] = bonds[m][n] = QueryBond(value)

        g = object.__new__(QueryContainer)
        g_atoms = {}
        g_bonds = {}
        plane = {}
        charges = {}
        radicals = {}
        atoms_stereo = {}
        allenes_stereo = {}
        cis_trans_stereo = {}
        neighbors = {}
        hybridizations = {}
        hydrogens = {}
        rings_sizes = {}
        heteroatoms = {}

        for n, atom in enumerate(atoms, 1):
            g_atoms[n] = ...
            g_bonds[n] = {}
            charges[n] = atom['charge']
            radicals[n] = ...
            plane[n] = (0., 0.)
        for n, m, b in data['bonds']:
            if m in bonds[n]:
                b = bonds[n][m]
            else:
                b = QueryBond(b)
            n += 1
            m += 1
            if n in g_bonds[m]:
                raise ValueError('atoms already bonded')
            g_bonds[n][m] = g_bonds[m][n] = b

        g.__setstate__({'atoms': g_atoms, 'bonds': g_bonds, 'plane': plane, 'charges': charges, 'radicals': radicals,
                        'neighbors': neighbors, 'hybridizations': hybridizations, 'hydrogens': hydrogens,
                        'rings_sizes': rings_sizes, 'heteroatoms': heteroatoms, 'atoms_stereo': atoms_stereo,
                        'allenes_stereo': allenes_stereo, 'cis_trans_stereo': cis_trans_stereo})
        return g


__all__ = ['SMARTSRead']
