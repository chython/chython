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
                # C(=O)O
                '[O;D1;x0;z1:2]-[C;x2;z2:1]=[O;M]',
            ],
            'B': [
                # CO
                '[O;D1;x0;z1:3]-[C;x1;z1;M]'
            ],
            'product': '[A:1]-[A:3]',
            # condition-specific untolerant groups
            'alerts': [
                '[S;D1;x0;z1][C;x1;z1]',  # thiol
                '[O,S;D1;z1][A;a]'  # [thia]phenol
            ]
        },
    ],
    'alerts': []  # global untolerant groups
}

_amidation = {
    'name': 'Amidation Reaction',
    'description': 'Amides formation from acids and amines',
    'templates': [
        {
            'A': [
                # [H,R]COOH
                '[O;x0;z2;M]=[C;x2:1][O;D1:2]'
            ],
            'B': [
                # Ar-NH2
                '[N;D1;x0;z1:3][C;a;M]',
                # Alk-NH2
                '[N;D1;x0;z1:3][C;z1;x1;M]',
                # Ar-NH-Ar
                '[N;D2;x0;z1:3]([C;a;M])[C;a;M]',
                # Alk-NH-Ar
                '[N;D2;x0;z1:3]([C;a;M])[C;z1;x1;M]',
                # Alk2NH
                '[N;D2;x0;z1:3]([C;z1;x1;M])[C;z1;x1;M]'
            ],
            'product': '[A:1]-[A:3]',
            'alerts': []
        },
    ],
    'alerts': ['[O;D1;x0;z1][C;z1;x1]', '[O;D1;z1][C,N;a]']  # global untolerant groups
}

_amine_carbonyl_reductive_amination = {
    'name': 'Amine carbonyl reductive amination reaction',
    'description': 'Amines formation from carbonyls and amines',
    'templates': [
        {
            'A': [
                # O=CR2
                '[O;x0;z2:2]=[C;x1:1]'
            ],
            'B': [
                # Ar-NH2
                '[N;D1;x0;z1:3][C;a;M]',
                # Alk-NH2
                '[N;D1;x0;z1:3][C;z1;x1;M]',
                # Alk-NH-Ar
                '[N;D2;x0;z1:3]([C;a;M])[C;z1;x1;M]',
                # Alk2NH
                '[N;D2;x0;z1:3]([C;z1;x1;M])[C;z1;x1;M]'
            ],
            'product': '[A:1]-[A:3]',
            'alerts': []
        },
    ],
    'alerts': []
}

_suzuki = {
    'name': 'Suzuki-Miyaura reaction',
    'description': 'Suzuki-Miyaura C-C coupling reaction',
    'templates': [
        {
            'A': [
                # X-Ar
                '[Cl,Br,I;D1:1]-[C;a:2]'
            ],
            'B': [
                # Ar-B
                '[B;D3;x2;z1:4]([O;M])([O;M])-[C;a:3]',
                # C=C-B
                '[B;D3;x2;z1:4]([O;M])([O;M])-[C;x1;z2:3]=[C;x0;z2;M]',

            ],
            'product': '[A:2]-[A:3]',
            'alerts': []
        },
        {
            'A': [
                # X-C=C
                '[Cl,Br,I;D1:2]-[C;x1;z2:1]=[C;x0;z2;M]'
            ],
            'B': [
                # C=C-B
                '[B;D3;x2;z1:4]([O;M])([O;M])-[C;x1;z2:3]=[C;x0;z2;M]'
            ],
            'product': '[A:2]-[A:3]',
            'alerts': []
        },
    ],
    'alerts': []
}

_buchwald_hartwig = {
    'name': 'Buchwald-Hartwig reaction ',
    'description': 'Buchwald-Hartwig amination reaction, C-N coupling reaction',
    'templates': [
        {
            'A': [
                # Hal-Ar
                '[Cl,Br,I;D1:1]-[C;a:2]'
            ],
            'B': [
                # Ar-NH2
                '[N;D1;x0;z1:3][C;a;M]',
                # Alk-NH2
                '[N;D1;x0;z1:3][C;z1;x1;M]',
                # Alk-NH-Ar
                '[N;D2;x0;z1:3]([C;a;M])[C;z1;x1;M]',
                # Alk2NH
                '[N;D2;x0;z1:3]([C;z1;x1;M])[C;z1;x1;M]'
            ],
            'product': '[A:2]-[A:3]',
            'alerts': []
        },
    ],
    'alerts': []
}

_sulfonamidation = {
    'name': 'Sulfoamination reaction ',
    'description': 'Sulfoamination reaction, S-N coupling reaction',
    'templates': [
        {
            'A': [
                # RS(=O)(=O)X
                '[S;D4;x3;z3:1]([O,F,Cl,Br,I;D1:2])(=[O;M])(=[O;M])[C;M]'
            ],
            'B': [
                # Ar-NH2
                '[N;D1;x0;z1:3][C;a;M]',
                # Alk-NH2
                '[N;D1;x0;z1:3][C;z1;x1;M]',
                # Ar-NH-Ar
                '[N;D2;x0;z1:3]([C;a;M])[C;a;M]',
                # Alk-NH-Ar
                '[N;D2;x0;z1:3]([C;a;M])[C;z1;x1;M]',
                # Alk2NH
                '[N;D2;x0;z1:3]([C;z1;x1;M])[C;z1;x1;M]'
            ],
            'product': '[A:1]-[A:3]',
            'alerts': []
        },
    ],
    'alerts': []
}

_amine_isocyanate = {
    'name': 'Amine with isocyanate reaction ',
    'description': 'Amine with isocyanate reaction, C-N coupling reaction',
    'templates': [
        {
            'A': [
                # RN=C=O
                '[C;D2;x2;z3:1](=[O;M])=[N;D2;x0;z2:2]'
            ],
            'B': [
                # Ar-NH2
                '[N;D1;x0;z1:3][C;a;M]',
                # Alk-NH2
                '[N;D1;x0;z1:3][C;z1;x1;M]',
                # Ar-NH-Ar
                '[N;D2;x0;z1:3]([C;a;M])[C;a;M]',
                # Alk-NH-Ar
                '[N;D2;x0;z1:3]([C;a;M])[C;z1;x1;M]',
                # Alk2NH
                '[N;D2;x0;z1:3]([C;z1;x1;M])[C;z1;x1;M]'
            ],
            'product': '[A:1][A:2]-[A:3]',
            'alerts': []
        },
    ],
    'alerts': []
}

_sonogashira = {
    'name': 'Sonogashira reaction ',
    'description': 'Sonogashira reaction, C-C coupling reaction. It employs a palladium catalyst as well as copper '
                   'co-catalyst',
    'templates': [
        {
            'A': [
                # HC#C-R
                '[C;D1;x0;z3:1]#[C;D2;x0;M]'
            ],
            'B': [
                # Ar-Hal
                '[Cl,Br,I;D1:3]-[C;a:2]',
                # C=C-Hal
                '[Cl,Br,I;D1:3]-[C;x1;z2:2]=[C;x0;z2;M]',
                # R-C(=O)-Hal
                '[Cl,Br,I;D1:3]-[C;x2;z2:2]=[O;M]'
            ],
            'product': '[A:1]-[A:2]',
            'alerts': []
        },
    ],
    'alerts': []
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
