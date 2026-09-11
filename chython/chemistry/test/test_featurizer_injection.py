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
"""Chemistry-side injection test: importing chython.chemistry registers all ten names.

The core-side file (`chython/core/test/test_featurizer_injection.py`) proves the core owns the slots.
`PROPERTIES` and `METHODS` are duplicated there deliberately: a shared helper would have to live under
one layer and be imported from the other, which is the import this split exists to prevent.
"""
from importlib.util import find_spec

import pytest
from chython.core import read_smiles


PROPERTIES = ('rotatable_bonds_count', 'hydrogen_bond_donors_count',
              'hydrogen_bond_acceptors_count', 'tpsa', 'crippen_logp', 'crippen_mr', 'qed')
METHODS = ('maccs_keys', 'maccs_bit_set', 'pharmacophore_invariants')

# A name moves in here when its body lands.  The two lists above are the whole registered surface, so
# what is not in here is exactly what must still raise `NotImplementedError`.  Every one of the ten is
# in here now, which is what `test_the_registered_surface_is_fully_implemented` asserts.
IMPLEMENTED = frozenset({'rotatable_bonds_count', 'hydrogen_bond_donors_count',
                         'hydrogen_bond_acceptors_count', 'tpsa', 'pharmacophore_invariants',
                         'crippen_logp', 'crippen_mr', 'qed', 'maccs_keys', 'maccs_bit_set'})

#: The implemented names whose answer is a numpy array, and numpy is optional (`chython[ml]`).  THIS
#: TEST IS NOT SKIPPED WHEN NUMPY IS ABSENT -- it is the injection ratchet, and a minimal install is
#: exactly where a broken hook would go unnoticed.  What changes is only what the call is allowed to
#: raise; see the docstring below.  `maccs_bit_set` is in here although it answers a `frozenset`: it
#: builds the vector first, so it asks for numpy exactly like the array-valued names.
NEEDS_NUMPY = (frozenset({'pharmacophore_invariants', 'maccs_keys', 'maccs_bit_set'})
               if find_spec('numpy') is None else frozenset())


UNIMPLEMENTED = tuple(n for n in PROPERTIES + METHODS if n not in IMPLEMENTED)


def test_the_registered_surface_is_fully_implemented():
    """No registered name is a stub any more, and this is where that stops being an assumption.

    It replaces a per-name check on `NotImplementedError` messages that had nothing left to run over:
    a `parametrize` across an empty `UNIMPLEMENTED` asserts nothing, silently.  Should a name ever be
    registered ahead of its body again, this fails and names it -- and the message ratchet belongs back
    in the same commit.
    """
    assert UNIMPLEMENTED == (), UNIMPLEMENTED
    assert IMPLEMENTED == frozenset(PROPERTIES + METHODS)


def test_importing_chemistry_registers_every_name():
    """Importing chython.chemistry registers all ten names via `_set_featurizer_fns`.

    A name in `IMPLEMENTED` must answer without raising; one not yet in it must raise
    `NotImplementedError`, which proves the slot is wired and only the body is absent.  `ImportError`
    is the failure this distinction exists to catch: it means the injection hook broke.

    WITH NUMPY ABSENT, ONE IMPLEMENTED NAME ANSWERS WITH AN `ImportError` OF ITS OWN, and the two are
    told apart by what the message says rather than by skipping the check.  `pharmacophore_invariants`
    is a numpy array, so on a minimal install it raises the core's single refusal naming `chython[ml]`
    -- which still proves what this test is here for: the slot is wired, the call reached the body, and
    the body got as far as asking for its optional dependency.  A broken injection hook raises an
    `ImportError` that does NOT name the extra, so it still fails here, which is why the match matters
    and a bare `pytest.raises(ImportError)` would not do.
    """
    import chython.chemistry  # noqa: F401 -- the import is the registration
    m = read_smiles('CCO')
    for name in PROPERTIES:
        if name in NEEDS_NUMPY:
            with pytest.raises(ImportError, match=r'chython\[ml\]'):
                getattr(m, name)
        elif name in IMPLEMENTED:
            getattr(m, name)      # must not raise
        else:
            with pytest.raises(NotImplementedError):
                getattr(m, name)  # must reach the stub, not die on ImportError
    for name in METHODS:
        if name in NEEDS_NUMPY:
            with pytest.raises(ImportError, match=r'chython\[ml\]'):
                getattr(m, name)()
        elif name in IMPLEMENTED:
            getattr(m, name)()
        else:
            with pytest.raises(NotImplementedError):
                getattr(m, name)()
