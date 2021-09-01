# -*- coding: utf-8 -*-
#
#  Copyright 2017-2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from operator import or_
from typing import Dict, Iterable, Iterator, Optional, Tuple
from .cgr import CGRContainer
from .molecule import MoleculeContainer
from ..algorithms.calculate2d import Calculate2DReaction
from ..algorithms.depict import DepictReaction
from ..algorithms.standardize import StandardizeReaction


class ReactionContainer(StandardizeReaction, Calculate2DReaction, DepictReaction):
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
        copy._ReactionContainer__meta = self.__meta.copy()
        copy._ReactionContainer__name = self.__name
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

    @classmethod
    def from_cgr(cls, cgr: CGRContainer) -> 'ReactionContainer':
        """
        Decompose CGR into reaction
        """
        if not isinstance(cgr, CGRContainer):
            raise TypeError('CGR expected')
        r, p = ~cgr
        return ReactionContainer(r, p)

    def flush_cache(self):
        self.__dict__.clear()
        for m in self.molecules():
            m.flush_cache()

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
        :param format_spec: see specification of nested containers.
            !c - Keep nested containers order
            !C - skip cxsmiles fragments contract
        """
        sig = []
        count = 0
        contract = []
        for ml in (self.__reactants, self.__reagents, self.__products):
            if not format_spec or '!c' not in format_spec:
                ml = sorted(ml, key=str)
            for m in ml:
                if m.connected_components_count > 1:
                    contract.append([str(x + count) for x in range(m.connected_components_count)])
                    count += m.connected_components_count
                else:
                    count += 1
            if format_spec:
                sig.append('.'.join(format(x, format_spec) for x in ml))
            else:
                sig.append('.'.join(str(x) for x in ml))
        if (not format_spec or '!C' not in format_spec) and contract:
            return f"{'>'.join(sig)} |f:{','.join('.'.join(x) for x in contract)}|"
        return '>'.join(sig)

    @cached_method
    def __len__(self):
        return len(self.__reactants) + len(self.__products) + len(self.__reagents)

    def __getstate__(self):
        state = {'reactants': self.__reactants, 'products': self.__products, 'reagents': self.__reagents,
                 'meta': self.__meta, 'name': self.__name, 'arrow': self._arrow, 'signs': self._signs}
        import chython
        if chython.pickle_cache:
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
