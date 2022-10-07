# -*- coding: utf-8 -*-
#
#  Copyright 2022 Ramil Nugmanov <nougmanoff@protonmail.com>
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

_amidation = (
    # alerts
    ['[C;z1;x1]-[O;D1]', '[C,N;a]-[O;D1]'],
    # [H,R]COOH + Ar-NH2
    (['[O;M:1]=[C;x2:2][O;D1:3]', '[N;D1:4][C;a;M:5]'], ['[A:2]-[A:4]']),  # 2 reactants, 1 product
    # Alk-NH2
    (['[O;M:1]=[C;x2:2][O;D1:3]', '[N;D1:4][C;z1;x1;M:5]'], ['[A:2]-[A:4]']),
    # Ar-NH-Ar
    (['[O;M:1]=[C;x2:2][O;D1:3]', '[N;D2:4]([C;a;M:5])[C;a;M:6]'], ['[A:2]-[A:4]']),
    # Alk-NH-Ar
    (['[O;M:1]=[C;x2:2][O;D1:3]', '[N;D2:4]([C;a;M:5])[C;z1;x1;M:6]'], ['[A:2]-[A:4]']),
    # Alk2NH
    (['[O;M:1]=[C;x2:2][O;D1:3]', '[N;D2:4]([C;z1;x1;M:5])[C;z1;x1;M:6]'], ['[A:2]-[A:4]']),
    # pyrrole
    (['[O;M:1]=[C;x2:2][O;D1:3]', '[N;a;h1:4]'], ['[A:2]-[A:4]'])
)

_amine_carbonyl_reductive_amination = (
    # alerts
    [],
    # R2C=O + Ar-NH2
    (['[C;x1;z2:1]=[O:2]', '[N;D1:3][C;a;M:4]'], ['[A:1]-[A:3]']),
    # Alk-NH2
    (['[C;x1;z2:1]=[O:2]', '[N;D1:3][C;z1;x1;M:4]'], ['[A:1]-[A:3]']),
    # Ar-NH-Ar
    # (['[C;x1;z2:1]=[O:2]', '[N;D2:3]([C;a;M:4])[C;a;M:5]'], ['[A:1]-[A:3]']),
    # Alk-NH-Ar
    (['[C;x1;z2:1]=[O:2]', '[N;D2:3]([C;a;M:4])[C;z1;x1;M:5]'], ['[A:1]-[A:3]']),
    # Alk2NH
    (['[C;x1;z2:1]=[O:2]', '[N;D2:3]([C;z1;x1;M:4])[C;z1;x1;M:5]'], ['[A:1]-[A:3]'])
)

_suzuki = (
    # alerts
    [],
    # Ar-X + Ar-B
    (['[C;a:1]-[Cl,Br,I;D1:2]', '[C;a:3]-[B;D3;z1:4]'], ['[A:1]-[A:3]']),
    # Ar-X + C=C-B
    (['[C;a:1]-[Cl,Br,I;D1:2]', '[C;x1;z2:3]-[B;D3;z1:4]'], ['[A:1]-[A:3]']),
    # C=C-X + C=C-B
    (['[C;x1;z2:1]-[Cl,Br,I;D1:2]', '[C;x1;z2:3]-[B;D3;z1:4]'], ['[A:1]-[A:3]']),
    # C=C-X + Ar-B
    (['[C;x1;z2:1]-[Cl,Br,I;D1:2]', '[C;a:3]-[B;D3;z1:4]'], ['[A:1]-[A:3]'])
)

_suzuki_amide = (
    # alerts
    [],
    (['[C;a:1]-[B;D3;z1:2]', '[C;a;M:3]-[C:4](=[O;M:5])-N([C;x1;z1])-C(=O)[C;x0;z1]'], ['[A:1]-[A:4]'])
)

_buchwald_hartwig = (
    # alerts
    [],
    # Ar-Hal + Ar-NH2
    (['[C;a:1]-[Cl,Br,I;D1:2]', '[N;D1:3][C;a;M:4]'], ['[A:1]-[A:3]']),
    # Alk-NH2
    (['[C;a:1]-[Cl,Br,I;D1:2]', '[N;D1:3][C;z1;x1;M:4]'], ['[A:1]-[A:3]']),
    # Ar-NH-Ar
    # (['[C;a:1]-[Cl,Br,I;D1:2]', '[N;D2:3]([C;a;M:4])[C;a;M:5]'], ['[A:1]-[A:3]']),
    # Alk-NH-Ar
    (['[C;a:1]-[Cl,Br,I;D1:2]', '[N;D2:3]([C;a;M:4])[C;z1;x1;M:5]'], ['[A:1]-[A:3]']),
    # Alk2NH
    (['[C;a:1]-[Cl,Br,I;D1:2]', '[N;D2:3]([C;z1;x1;M:4])[C;z1;x1;M:5]'], ['[A:1]-[A:3]'])
)

_sulfonamidation = (
    # alerts
    [],
    # RSOOX + Ar-NH2
    (['[O;M:1]=[S:2](=[O;M:3])([C;M:4])-[O,F,Cl,Br,I;D1:5]', '[N;D1:6][C;a;M:7]'], ['[A:2]-[A:6]']),
    # Alk-NH2
    (['[O;M:1]=[S:2](=[O;M:3])([C;M:4])-[O,F,Cl,Br,I;D1:5]', '[N;D1:6][C;z1;x1;M:7]'], ['[A:2]-[A:6]']),
    # Ar-NH-Ar
    (['[O;M:1]=[S:2](=[O;M:3])([C;M:4])-[O,F,Cl,Br,I;D1:5]', '[N;D2:6]([C;a;M:7])[C;a;M:8]'], ['[A:2]-[A:6]']),
    # Alk-NH-Ar
    (['[O;M:1]=[S:2](=[O;M:3])([C;M:4])-[O,F,Cl,Br,I;D1:5]', '[N;D2:6]([C;a;M:7])[C;z1;x1;M:8]'], ['[A:2]-[A:6]']),
    # Alk2NH
    (['[O;M:1]=[S:2](=[O;M:3])([C;M:4])-[O,F,Cl,Br,I;D1:5]', '[N;D2:6]([C;z1;x1;M:7])[C;z1;x1;M:8]'], ['[A:2]-[A:6]']),
    (['[O;M:1]=[S:2](=[O;M:3])([C;M:4])-[O,F,Cl,Br,I;D1:5]', '[N;a;h1:6]'], ['[A:2]-[A:6]'])
)

_amine_isocyanate = (
    # alerts
    [],
    # RN=C=O + Ar-NH2
    (['[C;M:1][N:2]=[C:3]=[O;M:4]', '[N;D1:5][C;a;M:6]'], ['[A:2][A:3]-[A:5]']),
    # Alk-NH2
    (['[C;M:1][N:2]=[C:3]=[O;M:4]', '[N;D1:5][C;z1;x1;M:6]'], ['[A:2][A:3]-[A:5]']),
    # Ar-NH-Ar
    (['[C;M:1][N:2]=[C:3]=[O;M:4]', '[N;D2:5]([C;a;M:6])[C;a;M:7]'], ['[A:2][A:3]-[A:5]']),
    # Alk-NH-Ar
    (['[C;M:1][N:2]=[C:3]=[O;M:4]', '[N;D2:5]([C;a;M:6])[C;z1;x1;M:7]'], ['[A:2][A:3]-[A:5]']),
    # Alk2NH
    (['[C;M:1][N:2]=[C:3]=[O;M:4]', '[N;D2:5]([C;z1;x1;M:6])[C;z1;x1;M:7]'], ['[A:2][A:3]-[A:5]'])
)

_sonogashira = (
    # alerts
    [],
    # Ar-Hal + HC#C-R
    (['[C;a:1]-[Cl,Br,I;D1:2]', '[C;D1:3]#[C;D2;M:4]'], ['[A:1]-[A:3]']),
    # C=C-Hal
    (['[C;x1;z2:1]-[Cl,Br,I;D1:2]', '[C;D1:3]#[C;D2;M:4]'], ['[A:1]-[A:3]']),
    # R-C(=O)-Hal
    (['[C;x2:1](=[O;M:2])-[Cl,Br,I;D1:3]', '[C;D1:4]#[C;D2;M:5]'], ['[A:1]-[A:4]'])
)


#################
# Magic Factory #
#################

__all__ = [k[1:] for k, v in globals().items() if k.startswith('_') and isinstance(v, tuple) and v]

_cache = {}


def _prepare_reactor(rules, name):
    rxn = [Reactor([smarts(r) for r in rs], [smarts(p) for p in ps], one_shot=False, automorphism_filter=False)  # noqa
           for rs, ps in rules[1:]]
    rxn_os = [Reactor([smarts(r) for r in rs], [smarts(p) for p in ps], automorphism_filter=False)  # noqa
              for rs, ps in rules[1:]]
    alerts = [smarts(x) for x in rules[0]]

    def w(*molecules, one_shot=False, rules_set: Optional[List[int]] = None,
          alerts_set: Optional[List[int]] = None, excess: Optional[List[int]] = None) -> Iterator[ReactionContainer]:
        """
        Generate reactions for given reactants.

        Rules:

        %s

        Alerts:

        %s

        :param molecules: Reactants molecules.
        :param one_shot: Generate only single stage products. Otherwise, all possible combinations, including products.
        :param rules_set: Select only specific rules.
        :param alerts_set: Check only specific structural alerts of reactants. Pass empty list to disable any checks.
        :param excess: Molecules indices which can be involved in multistep synthesis. All by default.
        """
        if not molecules:
            raise ValueError('empty molecule list')
        alerts_set = alerts if alerts_set is None else [alerts[x] for x in alerts_set]
        if any(a < m for a, m in product(alerts_set, molecules)):
            return

        molecules = fix_mapping_overlap(molecules)
        seen = set()
        if one_shot:
            picked_rxn = enumerate(rxn_os) if rules_set is None else [(x, rxn_os[x]) for x in rules_set]
            for i, rx in picked_rxn:
                for r in rx(*molecules):
                    if str(r) in seen:
                        continue
                    seen.add(str(r))
                    r.meta['MATCHED_RULES'] = [i]
                    yield r
            return

        excess = molecules if excess is None else [molecules[x] for x in excess]
        picked_rxn = list(enumerate(rxn)) if rules_set is None else [(x, rxn[x]) for x in rules_set]

        stack = deque([])
        for i, r in enumerate(picked_rxn):
            x = picked_rxn.copy()
            del x[i]
            stack.appendleft((r, [], molecules, x))

        while stack:
            (i, rx), chain, rct, nxt_rxn = stack.pop()
            chain = chain + [i]
            for r in rx(*rct):
                if str(r) in seen:
                    continue
                seen.add(str(r))

                r = ReactionContainer([x.copy() for x in molecules], r.products)
                r.meta['MATCHED_RULES'] = chain
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
                        stack.append((nrx, chain, x.copy(), z))
                else:  # drop one of the reactants
                    for n in range(len(r.products), len(x)):
                        y = x.copy()
                        del y[n]
                        for m, nrx in enumerate(nxt_rxn):
                            z = nxt_rxn.copy()
                            del z[m]
                            stack.append((nrx, chain, y, z))

    w.__module__ = __name__
    w.__qualname__ = w.__name__ = name
    doc = w.__doc__ % ('\n'.join(f'{n}: {".".join(r)}>>{".".join(p)}' for n, (r, p) in enumerate(rules[1:])),
                       '\n'.join(f'{n}: {r}' for n, r in enumerate(rules[0])) or '--')
    w.__doc__ = '\n'.join(x.lstrip() for x in doc.splitlines())
    return w


def __getattr__(name):
    try:
        return _cache[name]
    except KeyError:
        if name in __all__:
            _cache[name] = t = _prepare_reactor(globals()[f'_{name}'], name)
            return t
        raise AttributeError


def __dir__():
    return __all__
