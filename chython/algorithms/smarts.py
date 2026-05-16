# -*- coding: utf-8 -*-
#
#  Copyright 2025 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from chython._cache import cached_method
from typing import TYPE_CHECKING
from .smiles import Smiles, charge_str, order_str

if TYPE_CHECKING:
    from ..containers.query import QueryContainer

_aromatic_lower = {'B': 'b', 'C': 'c', 'N': 'n', 'O': 'o', 'P': 'p', 'S': 's', 'Se': 'se', 'Te': 'te'}
_hybridization_str = {1: '1', 2: '2', 3: '3', 4: '4'}


class Smarts(Smiles):
    """SMARTS writer mixin for QueryContainer.

    Provides ``str(qc)`` returning standard SMARTS and ``format(qc, spec)``
    for chython-specific or mapping-annotated variants.

    Format specs (combinable, order-independent):
        c — chython-specific extensions (z, x, M)
        m — include atom mapping numbers
    """
    __slots__ = ()

    def to_smarts(self: 'QueryContainer', chython_specific: bool = False,
                  mapping_renumber: dict = None) -> str:
        """Generate SMARTS string.

        :param chython_specific: Include z, x, M extensions.
        :param mapping_renumber: Optional dict {internal_atom_number: exported_mapping}.
            If provided, atom mapping is included using renumbered values.
        """
        kwargs = {}
        if chython_specific:
            kwargs['chython_specific'] = True
        if mapping_renumber is not None:
            kwargs['mapping'] = True
            kwargs['mapping_renumber'] = mapping_renumber

        if not kwargs:
            return str(self)

        smiles, order = self._smiles(self._smiles_order(), _return_order=True, **kwargs)
        if (cx := self._format_cxsmiles(order)) is not None:
            smiles.append(' ')
            smiles.append(cx)
        return ''.join(smiles)

    @cached_method
    def __str__(self: 'QueryContainer'):
        if self._smarts:
            return self._smarts
        smiles, order = self._smiles(self._smiles_order(), _return_order=True)
        if (cx := self._format_cxsmiles(order)) is not None:
            smiles.append(' ')
            smiles.append(cx)
        return ''.join(smiles)

    def __format__(self: 'QueryContainer', format_spec):
        if not format_spec:
            return str(self)
        kwargs = {}
        if 'c' in format_spec:
            kwargs['chython_specific'] = True
        if 'm' in format_spec:
            kwargs['mapping'] = True

        smiles, order = self._smiles(self._smiles_order(), _return_order=True, **kwargs)
        if (cx := self._format_cxsmiles(order)) is not None:
            smiles.append(' ')
            smiles.append(cx)
        return ''.join(smiles)

    def __repr__(self: 'QueryContainer'):
        return f'smarts({self})'

    def _smiles_order(self: 'QueryContainer', stereo=True):
        # QueryContainer has no Morgan ordering; use atom numbers directly.
        return lambda n: n

    def _format_atom(self: 'QueryContainer', n, adjacency, **kwargs):
        atom = self._atoms[n]
        chython_specific = kwargs.get('chython_specific', False)

        smi = ['[']

        # isotope
        isotope = getattr(atom, 'isotope', None)
        if isotope:
            smi.append(str(isotope))

        # element symbol — lowercase for aromatic
        symbol = atom.atomic_symbol
        hybridization = getattr(atom, 'hybridization', ())
        is_aromatic = hybridization == (4,)

        from ..periodictable.base.query import AnyElement, AnyMetal, ListElement
        if isinstance(atom, AnyMetal):
            if chython_specific:
                smi.append('M')
            else:
                smi.append('*')
        elif isinstance(atom, AnyElement):
            smi.append('*')
        elif isinstance(atom, ListElement):
            if is_aromatic:
                parts = symbol.split(',')
                smi.append(','.join(_aromatic_lower.get(p, p) for p in parts))
            else:
                smi.append(symbol)
        else:
            # QueryElement
            if is_aromatic and symbol in _aromatic_lower:
                smi.append(_aromatic_lower[symbol])
            else:
                smi.append(symbol)

        # stereo
        stereo = getattr(atom, 'stereo', None)
        if stereo is True:
            smi.append('@')
        elif stereo is False:
            smi.append('@@')

        # neighbors (D) — Daylight standard
        neighbors = atom.neighbors
        if neighbors:
            smi.append(';D')
            smi.append(',D'.join(str(x) for x in neighbors))

        # implicit hydrogens (h) — Daylight standard
        hydrogens = getattr(atom, 'implicit_hydrogens', ())
        if hydrogens:
            smi.append(';h')
            smi.append(',h'.join(str(x) for x in hydrogens))

        # ring sizes (r) — Daylight standard; (0,) means not in ring
        ring_sizes = getattr(atom, 'ring_sizes', ())
        if ring_sizes:
            if ring_sizes == (0,):
                smi.append(';!R')
            else:
                smi.append(';r')
                smi.append(',r'.join(str(x) for x in ring_sizes))

        # chython-specific extensions
        if chython_specific:
            # hybridization (z) — all values including aromatic
            if hybridization:
                smi.append(';z')
                smi.append(',z'.join(_hybridization_str[x] for x in hybridization))

            # heteroatoms (x)
            heteroatoms = getattr(atom, 'heteroatoms', ())
            if heteroatoms:
                smi.append(';x')
                smi.append(',x'.join(str(x) for x in heteroatoms))

            # masked (M)
            if getattr(atom, 'masked', False):
                smi.append(';M')

        # charge
        charge = getattr(atom, 'charge', 0)
        if isinstance(atom, AnyMetal):
            charge = 0
        if charge:
            smi.append(charge_str[charge])

        # mapping
        if kwargs.get('mapping', False):
            renumber = kwargs.get('mapping_renumber')
            m = renumber[n] if renumber else n
            smi.append(f':{m}')

        smi.append(']')
        return ''.join(smi)

    def _format_bond(self: 'QueryContainer', n, m, adjacency, **kwargs):
        bond = self._bonds[n][m]
        orders = bond.order
        if len(orders) == 1:
            return order_str[orders[0]]
        return ','.join(order_str[o] for o in orders)

    def _format_cxsmiles(self: 'QueryContainer', order):
        atoms = self._atoms
        rd = []
        for i, n in enumerate(order):
            a = atoms[n]
            if getattr(a, 'is_radical', False):
                rd.append(str(i))
        if rd:
            return '|^1:' + ','.join(rd) + '|'
        return None


def _query_smarts_body(q: 'QueryContainer'):
    """SMARTS body without trailing CX; returns (body, order) so callers
    can build a reaction-wide CX block instead of per-molecule."""
    smiles, order = q._smiles(q._smiles_order(), _return_order=True)
    return ''.join(smiles), order


__all__ = ['Smarts']
