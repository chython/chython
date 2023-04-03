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

_esterification =   {
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

__all__ = [k[1:] for k, v in globals().items() if k.startswith('_') and isinstance(v, dict) and v]
__all__.append('prepare_reactor')

_cache = {}


def prepare_reactor(rules, name):
    """
    Prepare reactor with predefined sets of templates.
    """
    rxn_ms = []
    rxn_os = []
    alerts = []
    g_alerts = [smarts(x) for x in rules['alerts']]
    for c in rules['templates']:
        c_alerts = [smarts(x) for x in c['alerts']] + g_alerts
        p = smarts(c['product'])
        for rs in product(*([smarts(x) for x in c[x]] for x in 'ABCD' if x in c)):
            rxn_ms.append(Reactor(rs, [p], one_shot=False, automorphism_filter=False))  # noqa
            rxn_os.append(Reactor(rs, [p], one_shot=True, automorphism_filter=False))  # noqa
            alerts.append(c_alerts)

    def w(*molecules, one_shot=True, check_alerts: bool = True,
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

        molecules = fix_mapping_overlap(molecules)
        seen = set()
        if one_shot:
            for rx, al in zip(rxn_os, alerts):
                if check_alerts and any(a < m for a, m in product(al, molecules)):
                    return
                for r in rx(*molecules):
                    if str(r) in seen:
                        continue
                    seen.add(str(r))
                    yield r
            return

        excess = molecules if excess is None else [molecules[x] for x in excess]
        stack = deque([])
        for i, (rx, al) in enumerate(zip(rxn_ms, alerts)):
            if check_alerts and any(a < m for a, m in product(al, molecules)):
                return
            x = rxn_ms.copy()
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

    w.__module__ = __name__
    w.__qualname__ = w.__name__ = name
    doc = w.__doc__ % rules['description']
    w.__doc__ = '\n'.join(x.lstrip() for x in doc.splitlines())
    return w


def __getattr__(name):
    try:
        return _cache[name]
    except KeyError:
        if name in __all__:
            _cache[name] = t = prepare_reactor(globals()[f'_{name}'], name)
            return t
        raise AttributeError


def __dir__():
    return __all__
