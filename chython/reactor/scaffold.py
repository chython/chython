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
from .. import smarts, MoleculeContainer
from .transformer import Transformer

"""
Predefined transformers for common reactive groups cleavage.
"""

_alcohol = (
    ('[O;D1;z1;x0][C;z1;x1:1]', '[A:1]',  # rule
     'CCO', 'CC',  # match test
     'c1ccccc1O'),  # false-match test
)

_phenol = (
    ('[O;D1;z1;x0][C;a:1]', '[A:1]', 'c1ccccc1O', 'c1ccccc1', 'CCO'),
)

_chloro_aryl = (
    ('[Cl;D1;z1;x0][C;a:1]', '[A:1]', 'c1ccccc1Cl', 'c1ccccc1', 'c1ccccc1Br', 'CCCl'),
)

_bromo_aryl = (
    ('[Br;D1;z1;x0][C;a:1]', '[A:1]', 'c1ccccc1Br', 'c1ccccc1', 'c1ccccc1I', 'CCBr'),
)

_iodo_aryl = (
    ('[I;D1;z1;x0][C;a:1]', '[A:1]', 'c1ccccc1I', 'c1ccccc1', 'c1ccccc1Cl', 'CCI'),
)

_chloro_alkyl = (
    ('[Cl;D1;z1;x0][C;x1;z1:1]', '[A:1]', 'CCCl', 'CC', 'c1ccccc1Cl'),
)

_bromo_alkyl = (
    ('[Br;D1;z1;x0][C;x1;z1:1]', '[A:1]', 'CCBr', 'CC', 'c1ccccc1Br'),
)

_iodo_alkyl = (
    ('[I;D1;z1;x0][C;x1;z1:1]', '[A:1]', 'CCI', 'CC', 'c1ccccc1I'),
)

_carboxy = (
    ('[O;D1;z1;x0][C;D3;!R;x2;z2](=[O;D1])[C:1]', '[A:1]', 'CCC(=O)O', 'CC'),
)

_chloro_anhydride = (
    ('[Cl;D1;z1;x0][C;D3;!R;x2;z2](=[O;D1])[C:1]', '[A:1]', 'CCC(=O)Cl', 'CC'),
)

_amine_primary = (
    ('[N;D1;z1;x0][C;x1;z1:1]', '[A:1]', 'CCN', 'CC', 'CNC'),
    ('[N;D1;z1;x0][C;a:1]', '[A:1]', 'c1ccccc1N', 'c1ccccc1', 'c1ccccc1NC'),
)

#################
# Magic Factory #
#################

_groups = [k[1:] for k, v in globals().items() if k.startswith('_') and isinstance(v, tuple) and v]
__all__ = ['apply_all'] + _groups
_cache = {}


def _prepare_reactor(rules, name):
    rxn = [Transformer(smarts(r), smarts(p)) for r, p, *_ in rules]

    def w(molecule: MoleculeContainer, /) -> MoleculeContainer:
        """
        Remove reactive groups from the given molecule if applicable.
        """
        for r in rxn:
            while True:
                try:
                    molecule = next(r(molecule))
                except StopIteration:
                    break
        return molecule

    w.__module__ = __name__
    w.__qualname__ = w.__name__ = name
    return w


def apply_all(molecule: MoleculeContainer, /) -> MoleculeContainer:
    """
    Remove all found reactive groups from the given molecule.
    """
    for name in _groups:
        molecule = __getattr__(name)(molecule)
    return molecule


def __getattr__(name):
    try:
        return _cache[name]
    except KeyError:
        if name in _groups:
            _cache[name] = t = _prepare_reactor(globals()[f'_{name}'], name)
            return t
        raise AttributeError


def __dir__():
    return __all__
