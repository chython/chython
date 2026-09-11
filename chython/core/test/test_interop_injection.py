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
"""The interop injection: the names, the hook, the forwarding and the error path.

The core-side half, and NO TOOLKIT IS INSTALLED FOR ANY OF IT -- the bodies are stubs registered
through the hook, which is what lets `pytest chython/core/` prove the core owns `to_rdkit` and the
other five names on a machine with neither RDKit nor a JVM.  What the real converters answer is
`chython/interop/test/`'s subject, and a test here that imported one would be that file's test run
from the wrong side.

RESTORING THE REGISTRATION IS PART OF EVERY TEST THAT REPLACES IT.  `import chython.core` imports the
package `chython`, which imports `chython.interop`, so the real dispatchers ARE registered by the time
this module is collected; leaving a stub in place would break the interop suite in the same session.
"""
from contextlib import contextmanager

import pytest
from chython.core import MoleculeContainer, ReactionContainer, _core


NAMES = ('rdkit', 'indigo', 'openbabel', 'cdk', 'cdpkit')
METHODS = ('to_rdkit', 'to_indigo', 'to_openbabel', 'to_cdk', 'to_cdpkit')


@contextmanager
def registered(**fns):
    """Register `fns` as the whole interop surface, then put the real one back."""
    saved = {n: _core._reaction_interop_fn(n) for n in NAMES}
    saved['iupac'] = _core._reaction_interop_fn('iupac')
    _core._set_interop_fns(**fns)
    try:
        yield
    finally:
        _core._set_interop_fns(**saved)


@pytest.mark.parametrize('name', METHODS + ('iupac',))
def test_name_exists_on_the_sealed_container(name):
    assert hasattr(MoleculeContainer, name), name


def test_a_reaction_has_the_rdkit_method_and_only_that_one():
    """The only toolkit with a reaction form gets the only reaction method.

    The four absences are the claim: a reaction has no Indigo, OpenBabel, CDK or CDPKit shape, and a
    method that only ever raised would advertise one.
    """
    assert hasattr(ReactionContainer, 'to_rdkit')
    for name in ('to_indigo', 'to_openbabel', 'to_cdk', 'to_cdpkit', 'iupac'):
        assert not hasattr(ReactionContainer, name), name


def test_the_hook_is_exported():
    assert callable(_core._set_interop_fns)
    assert callable(_core._reaction_interop_fn)


def test_an_unregistered_converter_names_its_package_in_the_error():
    with registered():  # nothing at all registered
        for name, method in zip(NAMES, METHODS):
            with pytest.raises(ImportError, match='chython.interop'):
                getattr(MoleculeContainer(), method)()
            with pytest.raises(ImportError, match='chython.interop'):
                _core._reaction_interop_fn(name)
        with pytest.raises(ImportError, match='chython.interop'):
            MoleculeContainer().iupac
        with pytest.raises(ImportError, match='chython.interop'):
            ReactionContainer().to_rdkit()


@pytest.mark.parametrize('name,method', list(zip(NAMES, METHODS)))
def test_each_method_reaches_its_own_converter(name, method):
    """No transposition: the method calls the body registered under ITS toolkit's name.

    Every stub answers its own name, so a swapped pair fails here rather than handing a caller an
    Indigo object out of `to_rdkit()`.
    """
    mol = MoleculeContainer()
    with registered(**{n: (lambda answer: lambda x, **kw: (answer, x, kw))(n) for n in NAMES}):
        answer, got, kwargs = getattr(mol, method)()
    assert answer == name
    assert got is mol
    assert kwargs == {}


def test_keywords_are_forwarded_untouched():
    """The method spells no keyword of its own, so the converter's signature stays the one authority."""
    mol = MoleculeContainer()
    with registered(rdkit=lambda x, **kw: kw):
        assert mol.to_rdkit(keep_mapping=False, keep_numbers=True) == {'keep_mapping': False,
                                                                       'keep_numbers': True}
        with pytest.raises(TypeError):  # a wrong keyword fails at the body, naming it
            MoleculeContainer().to_rdkit(1)


def test_the_iupac_property_takes_no_arguments_and_is_not_cached():
    """A property, and read twice it asks twice: nothing on a mutable container caches a name."""
    calls = []
    mol = MoleculeContainer()
    with registered(iupac=lambda x: calls.append(x) or 'methane'):
        assert mol.iupac == 'methane'
        assert mol.iupac == 'methane'
    assert calls == [mol, mol]


def test_the_reaction_method_reaches_the_rdkit_converter():
    rxn = ReactionContainer()
    with registered(rdkit=lambda x, **kw: ('rdkit', x, kw)):
        assert rxn.to_rdkit(keep_mapping=False) == ('rdkit', rxn, {'keep_mapping': False})
