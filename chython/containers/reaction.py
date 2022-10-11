# -*- coding: utf-8 -*-
#
#  Copyright 2017-2022 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from CachedMethods import cached_method
from functools import reduce
from hashlib import sha512
from itertools import chain
from math import ceil
from operator import itemgetter, or_
from typing import Dict, Iterable, Iterator, Optional, Tuple, List
from zlib import compress, decompress
from .cgr import CGRContainer
from .molecule import MoleculeContainer
from ..algorithms.calculate2d import Calculate2DReaction
from ..algorithms.depict import DepictReaction
from ..algorithms.mapping import Mapping
from ..algorithms.standardize import StandardizeReaction


class ReactionContainer(StandardizeReaction, Mapping, Calculate2DReaction, DepictReaction):
    """
    Reaction storage. Contains reactants, products and reagents lists.

    Reaction storage hashable and comparable. based on reaction unique signature (SMILES).
    """
    __slots__ = ('__reactants', '__products', '__reagents', '__meta', '__name', '_arrow', '_signs', '__dict__')
    __class_cache__ = {}

    def __init__(self, reactants: Iterable[MoleculeContainer] = (), products: Iterable[MoleculeContainer] = (),
                 reagents: Iterable[MoleculeContainer] = (), meta: Optional[Dict] = None, name: Optional[str] = None):
        """
        New reaction object creation

        :param reactants: list of MoleculeContainers in left side of reaction
        :param products: right side of reaction. see reactants
        :param reagents: middle side of reaction: solvents, catalysts, etc. see reactants
        :param meta: dictionary of metadata. like DTYPE-DATUM in RDF

        """
        reactants = tuple(reactants)
        products = tuple(products)
        reagents = tuple(reagents)
        if not reactants and not products and not reagents:
            raise ValueError('At least one graph object required')
        elif not all(isinstance(x, MoleculeContainer) for x in chain(reactants, products, reagents)):
            raise TypeError(f'MoleculeContainers expected')

        self.__reactants = reactants
        self.__products = products
        self.__reagents = reagents
        if meta is None:
            self.__meta = None
        else:
            self.__meta = dict(meta)
        if name is None:
            self.__name = None
        else:
            self.name = name
        self._arrow = None
        self._signs = None

    @property
    def reactants(self) -> Tuple[MoleculeContainer, ...]:
        return self.__reactants

    @property
    def reagents(self) -> Tuple[MoleculeContainer, ...]:
        return self.__reagents

    @property
    def products(self) -> Tuple[MoleculeContainer, ...]:
        return self.__products

    def molecules(self) -> Iterator[MoleculeContainer]:
        """
        Iterator of all reaction molecules
        """
        return chain(self.__reactants, self.__reagents, self.__products)

    @property
    def meta(self) -> Dict:
        """
        Dictionary of metadata.
        Like DTYPE-DATUM in RDF
        """
        if self.__meta is None:
            self.__meta = {}  # lazy
        return self.__meta

    @property
    def name(self) -> str:
        return self.__name or ''

    @name.setter
    def name(self, name: str):
        if not isinstance(name, str):
            raise TypeError('name should be string up to 80 symbols')
        self.__name = name

    def copy(self) -> 'ReactionContainer':
        """
        Get copy of object
        """
        copy = object.__new__(self.__class__)
        copy._ReactionContainer__reactants = tuple(x.copy() for x in self.__reactants)
        copy._ReactionContainer__products = tuple(x.copy() for x in self.__products)
        copy._ReactionContainer__reagents = tuple(x.copy() for x in self.__reagents)
        copy._ReactionContainer__name = self.__name
        if self.__meta is None:
            copy._ReactionContainer__meta = None
        else:
            copy._ReactionContainer__meta = self.__meta.copy()
        copy._arrow = self._arrow
        copy._signs = self._signs
        return copy

    @cached_method
    def compose(self) -> CGRContainer:
        """
        Get CGR of reaction

        Reagents will be presented as unchanged molecules
        :return: CGRContainer
        """
        rr = self.__reagents + self.__reactants
        if rr:
            r = reduce(or_, rr)
        else:
            r = MoleculeContainer()
        if self.__products:
            p = reduce(or_, self.__products)
        else:
            p = MoleculeContainer()
        return r ^ p

    def flush_cache(self):
        self.__dict__.clear()
        for m in self.molecules():
            m.flush_cache()

    def pack(self, *, compressed=True, check=True):
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
        data = b''.join((bytearray((1, len(self.__reactants), len(self.__reagents), len(self.__products))),
                         *(m.pack(compressed=False, check=check) for m in self.molecules())))
        if compressed:
            return compress(data, 9)
        return data

    @classmethod
    def pack_len(cls, data: bytes, /, *, compressed=True) -> Tuple[List[int], List[int], List[int]]:
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
        molecules = []
        shift = 4
        for _ in range(reactants + reagents + products):
            m, pl = MoleculeContainer.unpack(data[shift:], compressed=False, _return_pack_length=True)
            molecules.append(m)
            shift += pl
        return cls(molecules[:reactants], molecules[-products:], molecules[reactants: -products])

    def __invert__(self) -> CGRContainer:
        """
        Get CGR of reaction
        """
        return self.compose()

    def __eq__(self, other):
        return isinstance(other, ReactionContainer) and str(self) == str(other)

    @cached_method
    def __hash__(self):
        return hash(str(self))

    @cached_method
    def __bytes__(self):
        return sha512(str(self).encode()).digest()

    def __bool__(self):
        """
        Exists both reactants and products
        """
        return bool(self.__reactants and self.__products)

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
        """
        sig = []
        count = 0
        contract = []
        orders = []

        for ml in (self.__reactants, self.__reagents, self.__products):
            mso = [(m, *m.__format__(format_spec, _return_order=True)) for m in ml]
            if not format_spec or '!c' not in format_spec:
                mso.sort(key=itemgetter(1))

            ss = []
            for m, s, o in mso:
                if m.connected_components_count > 1:
                    contract.append([str(x + count) for x in range(m.connected_components_count)])
                    count += m.connected_components_count
                else:
                    count += 1

                orders.append((m, o))
                ss.append(s)
            sig.append('.'.join(ss))

        if not format_spec or '!x' not in format_spec:
            cx = []
            if r := ','.join(str(n) for n, (m, a) in enumerate((m, a) for m, o in orders for a in o) if m._radicals[a]):
                cx.append(f'^1:{r}')
            if contract:
                cx.append(f"f:{','.join('.'.join(x) for x in contract)}")
            if cx:
                return f"{'>'.join(sig)} |{','.join(cx)}|"
        return '>'.join(sig)

    @cached_method
    def __len__(self):
        return len(self.__reactants) + len(self.__products) + len(self.__reagents)

    def __getstate__(self):
        state = {'reactants': self.__reactants, 'products': self.__products, 'reagents': self.__reagents,
                 'meta': self.__meta, 'name': self.__name, 'arrow': self._arrow, 'signs': self._signs}
        from chython import pickle_cache

        if pickle_cache:
            state['cache'] = self.__dict__
        return state

    def __setstate__(self, state):
        self.__reactants = state['reactants']
        self.__products = state['products']
        self.__reagents = state['reagents']
        self.__meta = state['meta']
        self.__name = state['name']
        self._arrow = state['arrow']
        self._signs = state['signs']
        if 'cache' in state:
            self.__dict__.update(state['cache'])


__all__ = ['ReactionContainer']
