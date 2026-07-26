# -*- coding: utf-8 -*-
#
#  Copyright 2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from lazy_object_proxy import Proxy


# Role glossary: role name -> [(functional-group pattern name, cap template), ...].
# The cap template is a reactor product SMARTS. `[A:N]` preserves the source
# atom; the pattern's `:100` (and orphaned atoms) are deleted; the attachment
# atom `:1` is capped with a new `[At:20]`.
#
# NOTE: this content is a minimal illustrative starter set. It is curated by the
# maintainer from a chemical standpoint — extend `roles` here without touching
# the consuming code.
roles = {
    # aryl/alkyl/alkenyl/alkynyl halides (Cl/Br/I) — fluoride is a separate role
    'leaving_halogen': [
        ('aryl_chloride',    '[A:1]-[At:20]'),
        ('aryl_bromide',     '[A:1]-[At:20]'),
        ('aryl_iodide',      '[A:1]-[At:20]'),
        ('alkyl_chloride',   '[A:1]-[At:20]'),
        ('alkyl_bromide',    '[A:1]-[At:20]'),
        ('alkyl_iodide',     '[A:1]-[At:20]'),
        ('alkenyl_bromide',  '[A:1](=[A:2])-[At:20]'),
        ('alkynyl_bromide',  '[A:1](#[A:2])-[At:20]'),
    ],
    # fluoride kept separate from the other halides
    'fluoride': [
        ('aryl_fluoride',  '[A:1]-[At:20]'),
        ('alkyl_fluoride', '[A:1]-[At:20]'),
    ],
    # boron handles: B(OH)2 / Bpin / BF3(-) all cap the attachment carbon
    'boron': [
        ('aryl_boronic_acid',   '[A:1]-[At:20]'),
        ('aryl_boronic_ester',  '[A:1]-[At:20]'),
        ('alkyl_boronic_acid',  '[A:1]-[At:20]'),
        ('alkyl_boronic_ester', '[A:1]-[At:20]'),
    ],
    # acyl handle: keep the carbonyl O, cap the carbonyl C
    'acyl': [
        ('carboxylic_acid', '[A:1](=[A:2])-[At:20]'),
        ('acyl_chloride',   '[A:1](=[A:2])-[At:20]'),
    ],
    # amine handle: cap N, keep its carbon
    'amine': [
        ('primary_amine',   '[A:1](-[A:2])-[At:20]'),
        ('primary_aniline', '[A:1](-[A:2])-[At:20]'),
    ],
}


def _transformers():
    from ._functional import rules as functional_rules
    from ...reactor import Transformer
    from ... import smarts

    compiled = {}
    for role_name, entries in roles.items():
        built = []
        for fg_name, cap_template in entries:
            t = Transformer(functional_rules[fg_name], smarts(cap_template),
                            delete_atoms=True, automorphism_filter=True, canonicalize=True)
            built.append((fg_name, t))
        compiled[role_name] = built
    return compiled


transformers = Proxy(_transformers)


__all__ = ['roles', 'transformers']
