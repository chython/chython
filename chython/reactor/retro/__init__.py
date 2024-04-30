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
from chython import Reactor, ReactionContainer, smarts
from typing import Iterator
from ._amidation import template as amidation_template
from ._aryl_amination import template as aryl_amination_template
from ._mitsunobu import template as mitsunobu_template
from ._sonogashira import template as sonogashira_template
from ._suzuki_miyaura import template as suzuki_miyaura_template


__all__ = ['PreparedReactor']
__all__.extend(k[:-9] for k, v in globals().items() if k.endswith('_template') and isinstance(v, dict) and v)
_cache = {}


class PreparedReactor:
    """
    Prepared retrosynthetic reactors with predefined sets of templates.
    """
    def __init__(self, rules, name):
        self.name = name
        self.rules = rules

        self.rxn = rxn = []
        for tmp in rules['templates']:
            p = smarts(tmp['product'])
            rs = [smarts(x) for x in tmp['reactants']]
            rxn.append(Reactor([p], rs, automorphism_filter=False))  # noqa

    def __repr__(self):
        return f'{__name__}.{self.name}'

    def __str__(self):
        return f'RetroReactor<{self.rules["name"]}>'

    def __call__(self, molecule) -> Iterator[ReactionContainer]:
        """
        :param molecule: Product molecule
        """
        seen = set()
        for rx in self.rxn:
            for r in rx(molecule):
                if str(r) in seen:
                    continue
                seen.add(str(r))
                yield r


def __getattr__(name):
    try:
        return _cache[name]
    except KeyError:
        if name in __all__:
            _cache[name] = t = PreparedReactor(globals()[f'{name}_template'], name)
            return t
        raise AttributeError


def __dir__():
    return __all__
