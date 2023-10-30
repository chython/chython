# -*- coding: utf-8 -*-
#
#  Copyright 2022, 2023 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from collections import deque
from itertools import product
from typing import Iterator, Optional, List
from .. import smarts, ReactionContainer
from .reactor import Reactor, fix_mapping_overlap

"""
Predefined reactors for common reactions.
"""

_esterification = {
    'name': 'Fischer esterification',
    'description': 'Esters formation from alcohols and acids',
    'templates': [
        # reactants sets fully mixable
        {
            'A': [
                '[O;D1;x0;z1:2]-[C;x2;z2:1]=[O;M]'
            ],
            'B': [
                '[O;D1;x0;z1:3]-[C;x1;z1;M]'
            ],
            'product': '[A:1]-[A:3]',
            # condition-specific untolerant groups
            'alerts': [
                '[S;D1;x0;z1][C;x1;z1]',  # thiol
                '[O,S;D1;z1][A;a]'  # [thia]phenol
            ]
        }
    ],
    'alerts': []  # global untolerant groups
}


#################
# Magic Factory #
#################

__all__ = ['PreparedReactor', 'prepare_reactor']
__all__.extend(k[1:] for k, v in globals().items() if k.startswith('_') and isinstance(v, dict) and v)


class PreparedReactor:
    """
    Prepared reactors with predefined sets of templates.
    """
    def __init__(self, rules, name):
        self.name = name
        self.rules = rules

        self.rxn_ms = []
        self.rxn_os = []
        self.alerts = []

        self.global_alerts = [smarts(x) for x in rules['alerts']]

        for c in rules['templates']:
            alerts = [smarts(x) for x in c['alerts']]
            p = smarts(c['product'])
            for rs in product(*([smarts(x) for x in c[x]] for x in 'ABCD' if x in c)):
                self.rxn_ms.append(Reactor(rs, [p], one_shot=False, automorphism_filter=False))  # noqa
                self.rxn_os.append(Reactor(rs, [p], one_shot=True, automorphism_filter=False))  # noqa
                self.alerts.append(alerts)

    def __repr__(self):
        return f'{__name__}.{self.name}'

    def __str__(self):
        return f'Reactor<{self.rules["name"]}>'

    def __call__(self, *molecules, one_shot=True, check_alerts: bool = True,
          excess: Optional[List[int]] = None) -> Iterator[ReactionContainer]:
        """
        %s

        :param molecules: Reactants molecules.
        :param one_shot: Generate only single stage products. Otherwise, all possible combinations, including products.
        :param check_alerts: Check structural alerts of reactants.
        :param excess: Molecules indices which can be involved in multistep synthesis. All by default.
        """
        if not molecules:
            raise ValueError('empty molecule list')
        if check_alerts and any(a < m for a, m in product(self.global_alerts, molecules)):
            return

        molecules = fix_mapping_overlap(molecules)
        seen = set()
        if one_shot:
            for rx, al in zip(self.rxn_os, self.alerts):
                if check_alerts and any(a < m for a, m in product(al, molecules)):
                    continue
                for r in rx(*molecules):
                    if str(r) in seen:
                        continue
                    seen.add(str(r))
                    yield r
            return

        excess = molecules if excess is None else [molecules[x] for x in excess]
        stack = deque([])
        for i, (rx, al) in enumerate(zip(self.rxn_ms, self.alerts)):
            if check_alerts and any(a < m for a, m in product(al, molecules)):
                continue
            x = self.rxn_ms.copy()
            del x[i]
            stack.appendleft((rx, molecules, x))

        while stack:
            rx, rct, nxt_rxn = stack.pop()
            for r in rx(*rct):
                if str(r) in seen:
                    continue
                seen.add(str(r))

                r = ReactionContainer([x.copy() for x in molecules], r.products)
                yield r

                x = excess.copy()
                for p in reversed(r.products):
                    x.insert(0, p.copy())
                x = fix_mapping_overlap(x)
                if excess is not molecules:
                    # expected that product can react with all excess molecules simultaneously.
                    # e.g. multicomponent reaction (Ugi)
                    for m, nrx in enumerate(nxt_rxn):
                        z = nxt_rxn.copy()
                        del z[m]
                        stack.append((nrx, x.copy(), z))
                else:  # drop one of the reactants
                    for n in range(len(r.products), len(x)):
                        y = x.copy()
                        del y[n]
                        for m, nrx in enumerate(nxt_rxn):
                            z = nxt_rxn.copy()
                            del z[m]
                            stack.append((nrx, y, z))


prepare_reactor = PreparedReactor  # backward compatibility
_cache = {}


def __getattr__(name):
    try:
        return _cache[name]
    except KeyError:
        if name in __all__:
            _cache[name] = t = PreparedReactor(globals()[f'_{name}'], name)
            return t
        raise AttributeError


def __dir__():
    return __all__
