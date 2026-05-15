# -*- coding: utf-8 -*-
#
#  Copyright 2017-2025 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from collections import defaultdict
from functools import reduce
from itertools import chain
from math import ceil
from operator import itemgetter, or_
from typing import Optional, Union
from collections.abc import Iterator, Sequence
from zlib import compress, decompress
from .bonds import DynamicBond
from .cgr import CGRContainer
from .cgr_query import QueryCGRContainer
from .molecule import MoleculeContainer
from .query import QueryContainer
from ..algorithms.smarts import _query_smarts_body
from ..periodictable import DynamicQueryElement, DynamicAnyElement
from ..periodictable.base.query import AnyElement
from ..algorithms.calculate2d import Calculate2DReaction
from ..algorithms.depict import DepictReaction
from ..algorithms.mapping import Mapping
from ..algorithms.standardize import StandardizeReaction


GraphContainer = Union[MoleculeContainer, QueryContainer, CGRContainer, QueryCGRContainer]


class ReactionContainer(StandardizeReaction, Mapping, Calculate2DReaction, DepictReaction):
    """
    Reaction storage. Contains reactants, products and reagents lists.

    Reaction storage hashable and comparable. based on reaction unique signature (SMILES/SMARTS).
    """
    __slots__ = ('_reactants', '_products', '_reagents', '_meta', '_name', '_arrow', '_signs', '_graph_cls', '__dict__')

    def __init__(self, reactants: Sequence[GraphContainer] = (), products: Sequence[GraphContainer] = (),
                 reagents: Sequence[GraphContainer] = (), meta: Optional[dict] = None, name: Optional[str] = None):
        """
        New reaction object creation

        :param reactants: list of MoleculeContainers (or other supported graph containers) in left side of reaction
        :param products: right side of reaction. see reactants
        :param reagents: middle side of reaction: solvents, catalysts, etc. see reactants
        :param meta: dictionary of metadata. like DTYPE-DATUM in RDF

        """
        reactants = tuple(reactants)
        products = tuple(products)
        reagents = tuple(reagents)
        if not reactants and not products and not reagents:
            raise ValueError('At least one graph object required')
        graphs = reactants + reagents + products
        graph_cls = type(graphs[0])
        allowed = (MoleculeContainer, QueryContainer, CGRContainer, QueryCGRContainer)
        if not issubclass(graph_cls, allowed):
            raise TypeError('MoleculeContainer or QueryContainer or CGRContainer or QueryCGRContainer expected')
        if not all(isinstance(x, graph_cls) for x in graphs):
            raise TypeError(f'{graph_cls.__name__} expected')

        self._reactants = reactants
        self._products = products
        self._reagents = reagents
        self._graph_cls = graph_cls
        if meta is None:
            self._meta = None
        else:
            self._meta = dict(meta)
        if name is None:
            self._name = None
        else:
            self._name = name
        self._arrow = None
        self._signs = None

    @property
    def reactants(self) -> tuple[GraphContainer, ...]:
        return self._reactants

    def _get_graph_cls(self):
        graph_cls = getattr(self, '_graph_cls', None)
        if graph_cls is None:
            graphs = self._reactants + self._reagents + self._products
            if not graphs:
                raise ValueError('At least one graph object required')
            graph_cls = type(graphs[0])
            allowed = (MoleculeContainer, QueryContainer, CGRContainer, QueryCGRContainer)
            if not issubclass(graph_cls, allowed):
                raise TypeError('MoleculeContainer or QueryContainer or CGRContainer or QueryCGRContainer expected')
            self._graph_cls = graph_cls
        return graph_cls

    @property
    def reagents(self) -> tuple[GraphContainer, ...]:
        return self._reagents

    @property
    def products(self) -> tuple[GraphContainer, ...]:
        return self._products

    def molecules(self) -> Iterator[GraphContainer]:
        """
        Iterator of all reaction graphs
        """
        return chain(self.reactants, self.reagents, self.products)

    @property
    def meta(self) -> dict:
        """
        Dictionary of metadata.
        Like DTYPE-DATUM in RDF
        """
        if self._meta is None:
            self._meta = {}  # lazy
        return self._meta

    @property
    def name(self) -> str:
        return self._name or ''

    @name.setter
    def name(self, name: str):
        if not isinstance(name, str):
            raise TypeError('name should be string up to 80 symbols')
        self._name = name

    def copy(self) -> 'ReactionContainer':
        """
        Get copy of object
        """
        copy = object.__new__(self.__class__)
        copy._reactants = tuple(x.copy() for x in self.reactants)
        copy._products = tuple(x.copy() for x in self.products)
        copy._reagents = tuple(x.copy() for x in self.reagents)
        copy._name = self._name
        copy._graph_cls = self._get_graph_cls()
        if self._meta is None:
            copy._meta = None
        else:
            copy._meta = self._meta.copy()
        copy._arrow = self._arrow
        copy._signs = self._signs
        return copy

    @classmethod
    def from_cgr(cls, cgr: 'CGRContainer') -> 'ReactionContainer':
        """
        Decompose CGR into reaction

        :param cgr: CGRContainer to decompose
        :return: ReactionContainer with separate reactant and product molecules
        """
        if not isinstance(cgr, CGRContainer):
            raise TypeError('CGR expected')
        r, p = cgr.decompose()
        reaction = object.__new__(cls)
        reaction._reactants = tuple(r.split())
        reaction._products = tuple(p.split())
        reaction._reagents = ()
        reaction._graph_cls = MoleculeContainer
        reaction._meta = None
        reaction._name = None
        reaction._arrow = None
        reaction._signs = None
        return reaction

    def compose(self, *, dynamic=True) -> Union[CGRContainer, QueryCGRContainer]:
        """
        Get CGR of reaction

        Reagents will be presented as unchanged molecules
        :return: CGRContainer for MoleculeContainer reactions, QueryCGRContainer for QueryContainer reactions
        """
        graph_cls = self._get_graph_cls()
        if issubclass(graph_cls, QueryContainer):
            return self._compose_query()
        if issubclass(graph_cls, (CGRContainer, QueryCGRContainer)):
            raise TypeError('CGR-based reactions are not composable')
        if not issubclass(graph_cls, MoleculeContainer):
            raise TypeError('Only MoleculeContainer or QueryContainer reactions are composable')

        rr = self.reagents + self.reactants
        if rr:
            r = reduce(or_, rr)
        else:
            r = MoleculeContainer()
        if self.products:
            p = reduce(or_, self.products)
        else:
            p = MoleculeContainer()
        return r.compose(p, dynamic)

    def _compose_query(self) -> QueryCGRContainer:
        """
        Compose QueryContainer reaction into QueryCGRContainer.

        Maps atoms by their atom mapping numbers between reactants and products.
        """
        rr = self.reagents + self.reactants
        if rr:
            r = reduce(or_, rr)
        else:
            r = QueryContainer('')
        if self.products:
            p = reduce(or_, self.products)
        else:
            p = QueryContainer('')

        h = QueryCGRContainer()
        ha = h._atoms
        hb = h._bonds

        common = r._atoms.keys() & p._atoms.keys()
        bonds = []
        adj = defaultdict(lambda: defaultdict(lambda: [None, None]))

        def _qb_order(qbond):
            orders = qbond.order
            if len(orders) == 1:
                return orders[0]
            return 8  # any bond

        def _make_dqe(atom):
            if isinstance(atom, AnyElement):
                return DynamicAnyElement(isotope=getattr(atom, 'isotope', None))
            return DynamicQueryElement.from_atom(atom)

        def _get_charge(atom):
            return getattr(atom, '_charge', getattr(atom, 'charge', 0))

        def _get_radical(atom):
            return getattr(atom, '_is_radical', getattr(atom, 'is_radical', False))

        # Cleavage atoms (only in reactant)
        for n in r._atoms.keys() - common:
            ra = r._atoms[n]
            ha[n] = _make_dqe(ra)
            ha[n]._attach_to_graph(h, n)
            hb[n] = {}
            h._plane[n] = r._plane.get(n, (0., 0.))
            charge = _get_charge(ra)
            radical = _get_radical(ra)
            h._charges[n] = charge
            h._radicals[n] = radical
            h._p_charges[n] = charge
            h._p_radicals[n] = radical
            h._neighbors[n] = ra.neighbors
            h._hybridizations[n] = ra.hybridization
            h._p_neighbors[n] = ra.neighbors
            h._p_hybridizations[n] = ra.hybridization

            for m, bond in r._bonds[n].items():
                if m not in ha:
                    order = _qb_order(bond)
                    if m in common:
                        bonds.append((n, m, DynamicBond(order, None)))
                    else:
                        bonds.append((n, m, DynamicBond(order, order)))

        # Coupling atoms (only in product)
        for n in p._atoms.keys() - common:
            pa = p._atoms[n]
            ha[n] = _make_dqe(pa)
            ha[n]._attach_to_graph(h, n)
            hb[n] = {}
            h._plane[n] = p._plane.get(n, (0., 0.))
            charge = _get_charge(pa)
            radical = _get_radical(pa)
            h._charges[n] = charge
            h._radicals[n] = radical
            h._p_charges[n] = charge
            h._p_radicals[n] = radical
            h._neighbors[n] = pa.neighbors
            h._hybridizations[n] = pa.hybridization
            h._p_neighbors[n] = pa.neighbors
            h._p_hybridizations[n] = pa.hybridization

            for m, bond in p._bonds[n].items():
                if m not in ha:
                    order = _qb_order(bond)
                    if m in common:
                        bonds.append((n, m, DynamicBond(None, order)))
                    else:
                        bonds.append((n, m, DynamicBond(order, order)))

        # Common atoms bond orders
        for n in common:
            an = adj[n]
            for m, bond in r._bonds[n].items():
                if m in common:
                    an[m][0] = _qb_order(bond)
            for m, bond in p._bonds[n].items():
                if m in common:
                    an[m][1] = _qb_order(bond)

        # Common atoms
        for n in common:
            ra = r._atoms[n]
            pa = p._atoms[n]
            ha[n] = _make_dqe(ra)
            ha[n]._attach_to_graph(h, n)
            hb[n] = {}
            h._plane[n] = r._plane.get(n, (0., 0.))
            h._charges[n] = _get_charge(ra)
            h._radicals[n] = _get_radical(ra)
            h._p_charges[n] = _get_charge(pa)
            h._p_radicals[n] = _get_radical(pa)
            h._neighbors[n] = ra.neighbors
            h._hybridizations[n] = ra.hybridization
            h._p_neighbors[n] = pa.neighbors
            h._p_hybridizations[n] = pa.hybridization

            for m, (o1, o2) in adj[n].items():
                if m not in ha:
                    bonds.append((n, m, DynamicBond(o1, o2)))

        for n, m, bond in bonds:
            hb[n][m] = hb[m][n] = bond

        return h

    def to_reactor(self, **kwargs) -> 'Reactor':
        """Convert query reaction to Reactor.

        Only works for reactions composed of QueryContainers.

        :param kwargs: Extra keyword arguments passed to Reactor.__init__
            (delete_atoms, one_shot, etc.).
        """
        if not issubclass(self._graph_cls, QueryContainer):
            raise TypeError('Only QueryContainer reactions can be converted to Reactor')
        from ..reactor import Reactor
        return Reactor(patterns=self.reactants, products=self.products, **kwargs)

    def flush_cache(self, keep_molecule_cache=False, **kwargs):
        self.__dict__.clear()
        if not keep_molecule_cache:
            for m in self.molecules():
                m.flush_cache(**kwargs)

    def pack(self, *, compressed=True, check=True) -> bytes:
        """
        Pack into compressed bytes.

        Note:
            * Same restrictions as in molecules pack.
            * reactants, reagents nad products should contain less than 256 molecules.

        Format specification:
        Big endian bytes order
        8 bit - header byte = 0x01 (current format specification)
        8 bit - reactants count
        8 bit - reagents count
        8 bit - products count
        x bit - concatenated molecules packs

        :param compressed: return zlib-compressed pack.
        :param check: check molecules for format restrictions.
        """
        if not issubclass(self._get_graph_cls(), MoleculeContainer):
            raise TypeError('packing supported only for MoleculeContainer reactions')

        data = b''.join((bytearray((1, len(self.reactants), len(self.reagents), len(self.products))),
                         *(m.pack(compressed=False, check=check) for m in self.molecules())))
        if compressed:
            return compress(data, 9)
        return data

    def pach(self, *, compressed=True, check=True) -> bytes:
        """
        Pack into compressed bytes.
        """
        return self.pack(compressed=compressed, check=check)

    @classmethod
    def pack_len(cls, data: bytes, /, *, compressed=True) -> tuple[list[int], list[int], list[int]]:
        """
        Returns reactants, reagents, products molecules atoms count in reaction pack.
        """
        if compressed:
            data = decompress(data)
        data = memoryview(data)
        if data[0] != 1:
            raise ValueError('invalid pack header')
        reactants, reagents, products = data[1], data[2], data[3]

        v = data[4]  # mol pack version
        shift = 5  # RH+RC+RC+PC+MH
        molecules = []
        for _ in range(reactants + reagents + products - 1):
            acs = int.from_bytes(data[shift: shift + 3], 'big')
            neighbors = 0
            ac = acs >> 12
            shift += 4  # AC+CC+AN
            for _ in range(ac):
                neighbors += data[shift] & 0x0f
                shift += 9
            neighbors //= 2
            if v == 2:
                shift += 3 * neighbors + ceil(neighbors * 3 / 8) + (acs & 0x0fff) * 4
            elif v == 0:
                shift += 3 * neighbors + ceil(neighbors / 5) * 2 + (acs & 0x0fff) * 4
            molecules.append(ac)
        if reactants or reagents or products:
            molecules.append(int.from_bytes(data[shift: shift + 3], 'big') >> 12)
        return molecules[:reactants], molecules[reactants: -products], molecules[-products:]

    @classmethod
    def unpack(cls, data: bytes, /, *, compressed=True) -> 'ReactionContainer':
        """
        Unpack from compressed bytes.

        :param compressed: decompress data before processing.
        """
        if compressed:
            data = decompress(data)
        data = memoryview(data)
        if data[0] != 1:
            raise ValueError('invalid pack header')

        reactants, reagents, products = data[1], data[2], data[3]
        molecules: list[MoleculeContainer] = []
        shift = 4
        for _ in range(reactants + reagents + products):
            m, pl = MoleculeContainer.unpack(data[shift:], compressed=False, _return_pack_length=True)
            molecules.append(m)
            shift += pl
        return cls(molecules[:reactants], molecules[-products:], molecules[reactants: -products])

    @classmethod
    def unpach(cls, data: bytes, /, *, compressed=True) -> 'ReactionContainer':
        """
        Unpack from compressed bytes.
        """
        return cls.unpack(data, compressed=compressed)

    def __bytes__(self):
        return self.pack()

    @cached_method
    def __invert__(self) -> Union[CGRContainer, QueryCGRContainer]:
        """
        Get CGR of reaction
        """
        return self.compose()

    def __eq__(self, other):
        return isinstance(other, ReactionContainer) and str(self) == str(other)

    @cached_method
    def __hash__(self):
        return hash(str(self))

    def __bool__(self):
        """
        Exists both reactants and products
        """
        return bool(self.reactants and self.products)

    @cached_method
    def __str__(self):
        return format(self)

    def __format__(self, format_spec):
        """
        :param format_spec:
            !c - Keep nested containers order.
            a - Generate asymmetric closures.
            !s - Disable stereo marks.
            A - Use aromatic bonds instead aromatic atoms.
            m - Set atom mapping.
            r - Generate random-ordered smiles.
            h - Show implicit hydrogens.
            !b - Disable bonds tokens.
            !x - Disable CXSMILES extension.
            !z - Disable charge representation.
            Query-only reactions return SMARTS and respect only !c (keep order).
        """
        graph_cls = self._get_graph_cls()
        if issubclass(graph_cls, QueryContainer):
            keep_order = '!c' in format_spec
            disable_cx = '!x' in format_spec
            sig = []
            radicals = []  # absolute-index `is_radical` flag for every atom
            for ml in (self.reactants, self.reagents, self.products):
                per_mol = []
                for m in ml:
                    body, order = _query_smarts_body(m)
                    radicals.extend(m._atoms[n].is_radical for n in order)
                    per_mol.append(body)
                if keep_order:
                    sig.append('.'.join(per_mol))
                else:
                    sig.append('.'.join(sorted(per_mol)))
            base = '>'.join(sig)
            if disable_cx:
                return base
            rad_idx = ','.join(str(i) for i, r in enumerate(radicals) if r)
            if rad_idx:
                return f"{base} |^1:{rad_idx}|"
            return base

        sig = []
        count = 0
        contract = []
        radicals = []

        for ml in (self.reactants, self.reagents, self.products):
            mso = [(m, *m.__format__(format_spec, _return_order=True)) for m in ml]
            if not format_spec or '!c' not in format_spec:
                mso.sort(key=itemgetter(1))

            ss = []
            for m, s, o in mso:
                components_count = getattr(m, 'connected_components_count', 1)
                if components_count > 1:
                    contract.append([str(x + count) for x in range(components_count)])
                    count += components_count
                else:
                    count += 1

                radicals.extend(m.atom(n).is_radical for n in o)
                ss.append(s)
            sig.append('.'.join(ss))

        if not format_spec or '!x' not in format_spec:
            cx = []
            if r := ','.join(str(n) for n, r in enumerate(radicals) if r):
                cx.append(f'^1:{r}')
            if contract:
                cx.append(f"f:{','.join('.'.join(x) for x in contract)}")
            if cx:
                return f"{'>'.join(sig)} |{','.join(cx)}|"
        return '>'.join(sig)

    def __len__(self):
        return len(self.reactants) + len(self.products) + len(self.reagents)


__all__ = ['ReactionContainer']
