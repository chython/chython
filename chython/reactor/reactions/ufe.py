# -*- coding: utf-8 -*-
#
#  Copyright 2024 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from typing import Iterator
from ._amidation import template as amidation_template
from ._amine_isocyanate import template as amine_isocyanate_template
from ._buchwald_hartwig import template as buchwald_hartwig_template
from ._esterification import template as esterification_template
from ._macmillan import template as macmillan_template
from ._reductive_amination import template as reductive_amination_template
from ._sonogashira import template as sonogashira_template
from ._sulfonamidation import template as sulfonamidation_template
from ._suzuki_miyaura import template as suzuki_miyaura_template
from ..transformer import Transformer
from ... import MoleculeContainer, smarts
from ...periodictable import At


__all__ = ['PreparedUFE']
__all__.extend(k[:-9] for k, v in globals().items() if k.endswith('_template') and isinstance(v, dict) and v)
_cache = {}


class TransformerWrapper:
    def __init__(self, query, transformation, name):
        if isinstance(transformation, str):
            self.transformer = Transformer(smarts(query), smarts(transformation), copy_metadata=True,
                                           fix_aromatic_rings=False, fix_tautomers=False)
        else:
            self.query = smarts(query)
            self.mapping = transformation
            self.transformer = None
        self.name = name

    def __call__(self, molecule: MoleculeContainer) -> Iterator[MoleculeContainer]:
        if self.transformer is None:
            for mapping in self.query.get_mapping(molecule):
                n = mapping[self.mapping]
                copy = molecule.copy()
                copy._atoms[n].__class__ = At  # ad-hoc for masking leaving group
                copy._hydrogens[n] = 0
                copy.meta[self.name] = n
                yield copy
        else:
            for copy in self.transformer(molecule):
                copy.meta[self.name] = max(copy)
                yield copy


class PreparedUFE:
    def __init__(self, rules, name):
        self.name = name
        self.rules = rules
        self.transformations = []

        for n, rule in enumerate(rules['templates']):
            for g in 'AB':
                for s in rule[g]:
                    t = TransformerWrapper(s, rule['ufe'][g], f'{name}_{g}{n}')
                    self.transformations.append(t)

    def __call__(self, molecule: MoleculeContainer) -> Iterator[MoleculeContainer]:
        for transformer in self.transformations:
            yield from transformer(molecule)

    def __repr__(self):
        return f'{__name__}.{self.name}'

    def __str__(self):
        return f'UFE<{self.rules["name"]}>'


def __getattr__(name):
    try:
        return _cache[name]
    except KeyError:
        if name in __all__:
            _cache[name] = t = PreparedUFE(globals()[f'{name}_template'], name)
            return t
        raise AttributeError


def __dir__():
    return __all__
