# -*- coding: utf-8 -*-
#
#  Copyright 2021 Aleksandr Sizov <murkyrussian@gmail.com>
#  Copyright 2021, 2022 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from CachedMethods import cached_args_method
from functools import cached_property
from typing import Dict, Set, TYPE_CHECKING
from .linear import *
from .morgan import *
from .pharmacophore import *


if TYPE_CHECKING:
    from chython import MoleculeContainer


class Fingerprints(LinearFingerprint, MorganFingerprint, Pharmacophore):
    __slots__ = ()

    @cached_property
    def _screen_fingerprint(self) -> Dict[int, Set[int]]:
        """
        Fingerprint of available linear fragments with set of mapped atoms.
        Required for isomorphism tests filtering speedup.
        Parameters can be modified globally in `MoleculeContainer._fingerprint_config`.
        """
        from chython import fingerprint_config

        if fingerprint_config:
            return {hash(k): {x for x in v for x in x} for k, v in self._fragments(**fingerprint_config).items()}
        return {}

    @cached_args_method
    def _component_fingerprint(self: 'MoleculeContainer', component):
        """
        Fingerprint of specific component.
        """
        scope = set(self.connected_components[component])
        return {k: v & scope for k, v in self._screen_fingerprint.items() if not v.isdisjoint(scope)}


__all__ = ['Fingerprints']
