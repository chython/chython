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
from CachedMethods import cached_method
from functools import reduce
from itertools import chain
from math import ceil
from operator import itemgetter, or_
from typing import Dict, Iterator, Optional, Tuple, List, Sequence
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
    __slots__ = ('_reactants', '_products', '_reagents', '_meta', '_name', '_arrow', '_signs', '__dict__')

    def __init__(self, reactants: Sequence[MoleculeContainer] = (), products: Sequence[MoleculeContainer] = (),
                 reagents: Sequence[MoleculeContainer] = (), meta: Optional[Dict] = None, name: Optional[str] = None):
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
            raise TypeError('MoleculeContainers expected')

        self._reactants = reactants
        self._products = products
        self._reagents = reagents
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
    def reactants(self) -> Tuple[MoleculeContainer, ...]:
        return self._reactants

    @property
    def reagents(self) -> Tuple[MoleculeContainer, ...]:
        return self._reagents

    @property
    def products(self) -> Tuple[MoleculeContainer, ...]:
        return self._products

    def molecules(self) -> Iterator[MoleculeContainer]:
        """
        Iterator of all reaction molecules
        """
        return chain(self.reactants, self.reagents, self.products)

    @property
    def meta(self) -> Dict:
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
        if self._meta is None:
            copy._meta = None
        else:
            copy._meta = self._meta.copy()
        copy._arrow = self._arrow
        copy._signs = self._signs
        return copy

    def compose(self, *, dynamic=True) -> CGRContainer:
        """
        Get CGR of reaction

        Reagents will be presented as unchanged molecules
        :return: CGRContainer
        """
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
        molecules: List[MoleculeContainer] = []
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
        """
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
                if m.connected_components_count > 1:
                    contract.append([str(x + count) for x in range(m.connected_components_count)])
                    count += m.connected_components_count
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
