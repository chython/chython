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
"""The reaction corpus and the enumeration surface over it.

The corpus is TSV in `tables/`: `functional.tsv` names groups, `reactions.tsv` what happens to them,
`protective.tsv` the protecting groups, `roles.tsv` the handles a coupling cuts.  Importing this package
registers `react`, `functional_groups`, `functional_group_hits`, `protective_groups`,
`protective_group_hits`, `deprotect`, `sticky_fragments`, `sticky_linkers` and `@` onto
`MoleculeContainer` by injection -- the direction is `core <- reactions`, and
a special method must be compiled into the `cdef class` that owns the slot.

Mapping lives here too and is not part of the corpus: `reconstruct_mapping` derives a mapping by
rebuilding the recorded product from the recorded inputs through those templates, and `attention/` derives
one from a transformer's attention over the two sides.  The second imports neither the corpus nor
`onnxruntime` until it is called.

`functional_rules()`, `protective_rules()`, `reaction_rules()` and `roles()` are the corpus itself, on
the façade: every container method answers about one molecule, and these answer what the tables hold.
Each is keyed by the name its caller selects by -- `deprotect(protective=...)`, `react(reaction=...)`.
"""
from ._enumerate import (EnumeratedDeprotection, EnumeratedReaction, GroupHit, deprotect,
                         functional_group_hits, functional_groups, protective_group_hits,
                         protective_groups, react)
from .attention import attention_available, attention_mapping
from ._numbering import fast_mapping, mapping_agrees
from ._reconstruct import reconstruct_mapping
from ._stickers import StickyFragment, StickyLinker, sticky_fragments, sticky_linkers
from ._tables import (FunctionalGroup, ProtectiveGroup, ReactionRule, Role, SLOT_STRIDE,
                      compose_smirks, functional_rules, protective_rules, reaction_rules, read_table,
                      roles)
from ..core._core import _set_attention_fn, _set_reactions_fns, _set_reconstruct_fn


__all__ = ['EnumeratedDeprotection', 'EnumeratedReaction', 'FunctionalGroup', 'GroupHit',
           'ProtectiveGroup', 'ReactionRule', 'Role', 'StickyFragment', 'StickyLinker',
           # `attention_available` answers about the installation and not about a reaction, so it is a
           # function on the façade where every other name here is a container method or the corpus.
           'attention_available',
           'functional_rules', 'protective_rules', 'reaction_rules', 'roles']

_set_reactions_fns(deprotect=deprotect, functional_group_hits=functional_group_hits,
                   functional_groups=functional_groups, protective_group_hits=protective_group_hits,
                   protective_groups=protective_groups, react=react,
                   sticky_fragments=sticky_fragments, sticky_linkers=sticky_linkers)
_set_reconstruct_fn(reconstruct_mapping)
_set_attention_fn(attention_mapping)
