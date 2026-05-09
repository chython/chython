# -*- coding: utf-8 -*-
#
#  Copyright 2022-2026 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2024 Philippe Gantzer <p.gantzer@icredd.hokudai.ac.jp>
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
from functools import cache
from itertools import repeat
from logging import getLogger, INFO
from numpy import array, argmax, clip, concatenate, int32, int64, isclose, ix_, mean, nan_to_num
from numpy import nonzero, ones, unravel_index, zeros
from scipy.sparse.csgraph import shortest_path

logger = getLogger('chython.attention')
logger.setLevel(INFO)

_max_distance = 10
_max_neighbors = 14


@cache
def _get_session():
    import onnxruntime as ort
    from os import cpu_count

    try:
        from chython_rxnmap import model_path
    except ImportError:
        raise ImportError('chython-rxnmap package is required for attention mapping. '
                          'Install it with: pip install chython-rxnmap')

    opts = ort.SessionOptions()
    opts.inter_op_num_threads = 1
    opts.intra_op_num_threads = min(cpu_count() or 4, 8)
    opts.graph_optimization_level = ort.GraphOptimizationLevel.ORT_ENABLE_ALL
    return ort.InferenceSession(model_path, opts, providers=['CPUExecutionProvider'])


def _run_model(reactants, products, p2r, r2p, equal_atoms):
    """Run ONNX model on reaction, return symmetrized attention filtered by atom type."""
    atoms, neighbors, distances, roles = _encode_reaction(reactants, products)
    am = _get_session().run(None, {
        'atoms': atoms[None].astype(int64),
        'neighbors': neighbors[None].astype(int64),
        'distances': distances[None].astype(int64),
        'roles': roles[None].astype(int64),
    })[0]
    return (am[p2r] + am[r2p].T) * equal_atoms


def _encode_reaction(reactants, products):
    """Encode reaction into model input tensors (atoms, neighbors, distances, roles)."""
    atoms_list = [array([0], dtype=int32)]  # rxn_cls
    neighbors_list = [array([0], dtype=int32)]
    roles_list = [1]  # rxn_cls role
    distances_list = []

    for mol in reactants:
        a, n, d = _encode_molecule(mol)
        atoms_list.append(a)
        neighbors_list.append(n)
        distances_list.append(d)
        roles_list.append(0)  # mol_cls (hidden)
        roles_list.extend(repeat(2, len(mol)))

    for mol in products:
        a, n, d = _encode_molecule(mol)
        atoms_list.append(a)
        neighbors_list.append(n)
        distances_list.append(d)
        roles_list.append(0)
        roles_list.extend(repeat(3, len(mol)))

    atoms = concatenate(atoms_list)
    neighbors = concatenate(neighbors_list)
    roles = array(roles_list, dtype=int32)

    total = len(roles)
    distances = zeros((total, total), dtype=int32)
    distances[0, 0] = 1  # rxn_cls self-loop
    i = 1
    for d in distances_list:
        j = i + d.shape[0]
        distances[i:j, i:j] = d
        i = j

    return atoms, neighbors, distances, roles


def _encode_molecule(mol):
    """Encode molecule into (atoms, neighbors, distances) tensors.

    Encoding:
        atoms: 0=padding, 1=mol_cls, 2+=atomic_number+2
        neighbors: 0=cls, 2+=neighbor_count+2
        distances: 0=padding, 1=unreachable/cross-component, 2+=shortest_path+2
    """
    n_atoms = len(mol)
    size = n_atoms + 1  # +1 for mol_cls

    atoms = zeros(size, dtype=int32)
    neighbors = zeros(size, dtype=int32)
    atoms[0] = 1  # mol_cls

    bonds = mol._bonds
    for i, (n, a) in enumerate(mol.atoms(), 1):
        atoms[i] = a.atomic_number + 2
        nb = len(bonds[n]) + (a.implicit_hydrogens or 0)
        if nb > _max_neighbors:
            nb = _max_neighbors
        neighbors[i] = nb + 2

    adj = mol.adjacency_matrix()
    dist = shortest_path(adj, method='FW', directed=False, unweighted=True)
    nan_to_num(dist, copy=False, posinf=1.0)  # unreachable -> 1 (cross-component attention)
    clip(dist, None, _max_distance, out=dist)
    dist = (dist + 2).astype(int32)

    distances = ones((size, size), dtype=int32)  # mol_cls distance=1 to all atoms
    distances[1:, 1:] = dist
    return atoms, neighbors, distances


def _prepare_masks(reactants, products, reagents):
    """Build atom maps, adjacency matrices, element-equality mask, and token-level index slices."""
    r_map = [n for m in reactants for n in m]
    p_map = [n for m in products for n in m]
    rg_map = [n for m in reagents for n in m]
    ra = len(r_map)
    pa = len(p_map)

    # token layout: [rxn_cls, mol_cls_1, atoms_1, ..., mol_cls_n, atoms_n, ...]
    ram = [False]  # rxn_cls
    r_atoms = []
    r_adj = zeros((ra, ra), dtype=bool)
    i = 0
    for m in reactants:
        ram.append(False)  # mol_cls
        ram.extend(repeat(True, len(m)))
        j = i + len(m)
        r_adj[i:j, i:j] = m.adjacency_matrix()
        i = j
        r_atoms.extend(a.atomic_number for _, a in m.atoms())

    pam = [False] * len(ram)
    p_atoms = []
    p_adj = zeros((pa, pa), dtype=bool)
    i = 0
    for m in products:
        pam.append(False)
        pam.extend(repeat(True, len(m)))
        j = i + len(m)
        p_adj[i:j, i:j] = m.adjacency_matrix()
        i = j
        p_atoms.extend(a.atomic_number for _, a in m.atoms())

    ram.extend(repeat(False, len(pam) - len(ram)))
    ram = array(ram, dtype=bool)
    pam = array(pam, dtype=bool)
    r_atoms = array(r_atoms, dtype=int)
    p_atoms = array(p_atoms, dtype=int)

    equal_atoms = p_atoms[:, None] == r_atoms
    return r_map, p_map, rg_map, equal_atoms, ix_(pam, ram), ix_(ram, pam), r_adj, p_adj


def _greedy_mapping(am, r_map, p_map, r_adj, p_adj, multiplier):
    """Greedy attention-based product-to-reactant atom assignment."""
    amc = am.copy()
    pa = len(p_map)
    mapping = {}
    scope = zeros(pa, dtype=bool)
    seen = ones(pa, dtype=bool)
    score = []

    for x in range(pa):
        if not x:
            i, j = unravel_index(argmax(am), am.shape)
        else:
            ams = am[scope]
            if ams.size:
                i, j = unravel_index(argmax(ams), ams.shape)
                i = nonzero(scope)[0][i]
            else:
                i, j = unravel_index(argmax(am), am.shape)

        if isclose(am[i, j], 0.):
            for n in set(p_map).difference(mapping):
                mapping[n] = 0
            break

        score.append(amc[i, j])
        mapping[p_map[i]] = r_map[j]
        am[ix_(p_adj[i], r_adj[j])] *= multiplier
        am[i] = am[:, j] = 0
        seen[i] = False
        scope[i] = False
        scope[p_adj[i] & seen] = True

    return mapping, float(mean(score)) if score else 0.


class Attention:
    __slots__ = ()

    def attention_mapping(self, *, return_score: bool = False, multiplier=1.75,
                          keep_reactants_numbering=False) -> bool | float:
        """Do atom-to-atom mapping. Return True if mapping changed."""
        if any(len(bs) > _max_neighbors for m in self.molecules() for bs in m._bonds.values()):
            logger.info('atom-to-atom mapping not supported for hypervalent compounds')
            return False

        fixed = self.reset_mapping()
        r_map, p_map, rg_map, equal_atoms, p2r, r2p, r_adj, p_adj = _prepare_masks(
            self.reactants, self.products, self.reagents
        )

        am = _run_model(self.reactants, self.products, p2r, r2p, equal_atoms)
        mapping, score = _greedy_mapping(am, r_map, p_map, r_adj, p_adj, multiplier)

        if any(n != m for n, m in mapping.items()):
            if keep_reactants_numbering or fixed:
                r_mapping = {n: n for n in r_map}
            else:
                r_mapping = {m: n for n, m in enumerate(r_map, 1)}
                for m in self.reactants:
                    m.remap(r_mapping)

            p_mapping = {}
            nm = max(r_mapping.values())
            for n, m in mapping.items():
                if m := r_mapping.get(m):
                    p_mapping[n] = m
                else:
                    nm += 1
                    p_mapping[n] = nm
            for m in self.products:
                m.remap(p_mapping)

            if not keep_reactants_numbering and not fixed:
                rg_mapping = {m: n for n, m in enumerate(rg_map, nm + 1)}
                for m in self.reagents:
                    m.remap(rg_mapping)

            self.flush_cache()
            fixed = True

        if self.fix_mapping():
            fixed = True
        return score if return_score else fixed


__all__ = ['Attention']
