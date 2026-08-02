# -*- coding: utf-8 -*-
#
#  Copyright 2019-2026 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2019, 2020 Dinar Batyrshin <batyrshin-dinar@mail.ru>
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
from importlib.resources import files
from typing import Literal
from ...exceptions import ImplementationError
from ...periodictable.base.vector import Vector

try:
    from py_mini_racer import MiniRacer

    ctx = MiniRacer()
    ctx.eval(files(__package__).joinpath('clean2d.js').read_text())
except (ImportError, RuntimeError):
    ctx = None


class Calculate2DMolecule:
    __slots__ = ()

    def clean2d(self, *, engine: Literal['rdkit', 'smilesdrawer', 'cdk', 'obabel', 'indigo'] = None):
        """
        Calculate 2d layout of graph.

        By default, https://pubs.acs.org/doi/10.1021/acs.jcim.7b00425 JS implementation is used.
        Can be changed globally with the `chython.clean2d_engine` parameter.

        :param engine: override globally set engine
        """
        if engine is None:
            from chython import clean2d_engine as engine

        plane = {}
        if engine == 'rdkit':
            from rdkit.Chem.AllChem import Compute2DCoords

            mol = self.to_rdkit(keep_mapping=False)
            Compute2DCoords(mol)
            # set coordinates from the first rdkit conformer. usually it's 2d layout
            for n, (x, y, _) in zip(self, mol.GetConformers()[0].GetPositions()):
                plane[n] = (x, y)
        elif engine == 'smilesdrawer':
            if ctx is None:
                raise ImportError('mini_racer is not installed or broken')
            # smiles-drawer normalizes the layout regardless of the tree root, so a single
            # deterministic layout pass is enough.
            tree, order = self.__clean2d_tree()
            try:
                xy = ctx.call('$.clean2d', tree)
            except Exception:
                raise ImplementationError

            shift_x, shift_y = xy[0]
            for n, (x, y) in zip(order, xy):
                plane[n] = (x - shift_x, shift_y - y)
        elif engine == 'cdk':
            from ..._java import get_cdk

            sdg = get_cdk().layout.StructureDiagramGenerator()
            sdg.setUseTemplates(False)
            sdg.setMolecule(self.to_cdk())
            sdg.generateCoordinates()
            mol = sdg.getMolecule()

            for i, n in enumerate(self.smiles_atoms_order):
                xy = mol.getAtom(i).getPoint2d()
                plane[n] = (xy.x, xy.y)
        elif engine == 'obabel':
            from openbabel import openbabel

            mol = self.to_openbabel()
            assert openbabel.OBOp.FindType('gen2D').Do(mol), 'OpenBabel failed to generate 2d layout'
            assert mol.NumAtoms() == len(self), 'OpenBabel modified molecule'

            # to_openbabel preserves atom order: OBMol index i (1-based) matches the i-th atom
            for i, n in enumerate(self, 1):
                xy = mol.GetAtom(i).GetVector()
                plane[n] = (xy.GetX(), xy.GetY())
        elif engine == 'indigo':
            mol = self.to_indigo()
            assert not mol.layout(), 'Indigo failed to generate 2d layout'

            # to_indigo preserves atom order: iterateAtoms() matches self order
            for n, a in zip(self, mol.iterateAtoms()):
                x, y, _ = a.xyz()
                plane[n] = (x, y)
        else: raise ValueError(f'Invalid clean2d engine: {engine}')

        atoms = self._atoms
        for n, xy in plane.items():
            atoms[n].xy = xy
        self.rescale2d()

        if self.connected_components_count > 1:
            shift_x = 0.
            for c in self.connected_components:
                shift_x = self._fix_plane_mean(shift_x, component=c) + .9
        self.__dict__.pop('__cached_method__repr_svg_', None)

    def rescale2d(self):
        """
        Rescale coordinates to average bond length 0.825.
        """
        bonds = []
        atoms = self._atoms
        for n, m, _ in self.bonds():
            bonds.append(float(atoms[n].xy - atoms[m].xy))
        if bonds:
            bond_reduce = sum(bonds) / len(bonds) / .825
            if bond_reduce > .5:  # check for singularity
                for a in atoms.values():
                    a.xy /= bond_reduce
                return True
        return False

    def _fix_plane_mean(self, shift_x: float, shift_y=0., component=None) -> float:
        atoms = self._atoms
        if component is None:
            component = atoms

        left_atom = atoms[min(component, key=lambda x: atoms[x].x)]
        right_atom = atoms[max(component, key=lambda x: atoms[x].x)]

        min_x = left_atom.x - shift_x
        if len(left_atom.atomic_symbol) == 2:
            min_x -= .2

        max_x = right_atom.x - min_x
        min_y = min(atoms[x].y for x in component)
        max_y = max(atoms[x].y for x in component)
        mean_y = (max_y + min_y) / 2 - shift_y
        delta = Vector(min_x, mean_y)
        for n in component:
            atoms[n].xy -= delta

        if -.18 <= right_atom.y <= .18:
            factor = right_atom.implicit_hydrogens
            if factor == 1:
                max_x += .15
            elif factor:
                max_x += .25
        return max_x

    def _fix_plane_min(self, shift_x: float, shift_y=0., component=None) -> float:
        atoms = self._atoms
        if component is None:
            component = atoms

        right_atom = atoms[max(component, key=lambda x: atoms[x].x)]
        min_x = min(atoms[x].x for x in component) - shift_x
        max_x = right_atom.x - min_x
        min_y = min(atoms[x].y for x in component) - shift_y
        delta = Vector(min_x, min_y)
        for n in component:
            atoms[n].xy -= delta

        if shift_y - .18 <= right_atom.y <= shift_y + .18:
            factor = right_atom.implicit_hydrogens
            if factor == 1:
                max_x += .15
            elif factor:
                max_x += .25
        return max_x

    def __clean2d_tree(self):
        """
        Build a smiles-drawer parse tree directly from the molecule graph. Any spanning
        tree is fine: the layout only needs valid connectivity (and is root-invariant), so
        this is a plain iterative DFS with ring-closure detection, not a canonical SMILES
        traversal. All neighbours are emitted as `branches`; `next` is used only to chain
        `.`-separated connected components. Returns `(tree, order)` where `order[i]` is the
        atom number of the i-th heavy atom created (matching smiles-drawer's atom index).

        See clean2d/README for the parse-tree node schema.
        """
        atoms = self._atoms
        bonds = self._bonds
        order = []
        nodes = {}  # atom number -> its node dict
        parent_of = {}  # atom number -> tree parent atom number (root: None)

        # each atom becomes one node. layout depends only on element and connectivity, so
        # `atom` is a bare element string (charge, isotope, hydrogens and even aromaticity
        # do not move atoms). bond order maps to a smiles-drawer token; orders without a
        # single/double/triple counterpart (aromatic 4, any 8) collapse to '-'.
        bond_symbol = {1: '-', 2: '=', 3: '#', 4: '-', 8: '-'}

        # spanning forest via iterative DFS. all neighbours become `branches`.
        components = []
        for root in self:
            if root in nodes:
                continue
            nodes[root] = rnode = {'atom': atoms[root].atomic_symbol, 'isBracket': False,
                                   'branches': [], 'branchCount': 0, 'ringbonds': [], 'ringbondCount': 0,
                                   'bond': '-', 'branchBond': '-', 'next': None, 'hasNext': False}
            order.append(root)
            components.append(rnode)
            parent_of[root] = None
            stack = [root]
            while stack:
                parent = stack[-1]
                for child in bonds[parent]:
                    if child not in nodes:
                        bs = bond_symbol[bonds[parent][child].order]
                        nodes[child] = cnode = {'atom': atoms[child].atomic_symbol, 'isBracket': False,
                                                'branches': [], 'branchCount': 0, 'ringbonds': [],
                                                'ringbondCount': 0, 'bond': bs, 'branchBond': bs,
                                                'next': None, 'hasNext': False}
                        order.append(child)
                        pnode = nodes[parent]
                        pnode['branches'].append(cnode)
                        pnode['branchCount'] += 1
                        parent_of[child] = parent
                        stack.append(child)
                        break
                else:
                    stack.pop()

        # ring closures: every non-tree edge gets a matching ringbond id on both ends.
        cycle = 0
        for n in order:
            for m in bonds[n]:
                if m <= n or parent_of.get(m) == n or parent_of.get(n) == m:
                    continue
                cycle += 1
                bs = bond_symbol[bonds[n][m].order]
                for k in (n, m):
                    nodes[k]['ringbonds'].append({'bond': bs, 'id': cycle})
                    nodes[k]['ringbondCount'] += 1

        # chain disconnected components through `next` with a '.' bond.
        for prev, comp in zip(components, components[1:]):
            comp['bond'] = '.'
            prev['next'] = comp
            prev['hasNext'] = True
        return components[0], order


__all__ = ['Calculate2DMolecule']
