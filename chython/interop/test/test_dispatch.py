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
"""
The dispatch layer, tested without any toolkit installed.

Every callable picks its direction by asking whether the argument is a chython container.  The failure
mode is silent: a container that fails the predicate goes down the import path and blames the caller.
"""
from pytest import mark, raises

from chython.core import MoleculeContainer, QueryContainer, ReactionContainer, read_smarts, read_smiles
from chython.exceptions import DirectionNotImplemented
from chython.interop import cdk, cdpkit, indigo, iupac, is_container, openbabel, rdkit


CALLABLES = [rdkit, indigo, openbabel, cdk, cdpkit, iupac]

# What a caller gets from each half; implementing a direction costs one edit here.
#
# STUB   -- not built yet; raises `DirectionNotImplemented`.  Flip to LIVE when it is built.
# ABSENT -- deliberately never built; raises `DirectionNotImplemented` too, and does NOT flip.
# LIVE   -- built.  These tests stop probing it and its own test module owns it from here.
STUB, ABSENT, LIVE = 'stub', 'absent', 'live'

DIRECTIONS = [
    # callable     export  import
    (rdkit,        LIVE,   LIVE),
    (indigo,       LIVE,   LIVE),
    (openbabel,    LIVE,   LIVE),
    (cdk,          LIVE,   LIVE),
    (cdpkit,       LIVE,   ABSENT),   # export-only by decision, not by backlog
    (iupac,        LIVE,   LIVE),
]

# The toolkits that must not be imported as a side effect of choosing a direction, or of `import
# chython`.  Top-level distribution names, since that is what lands in `sys.modules`.
TOOLKIT_ROOTS = {'rdkit', 'indigo', 'openbabel', 'CDPL', 'jpype', 'openclatura'}

STUBBED_EXPORT = [fn for fn, e, _ in DIRECTIONS if e is not LIVE]
STUBBED_IMPORT = [fn for fn, _, i in DIRECTIONS if i is not LIVE]


def test_the_table_covers_every_callable():
    """
    The negative control for the table itself.

    Every test below iterates a filtered view of `DIRECTIONS`, so a callable missing from the table
    would be silently unprobed.
    """
    assert [fn for fn, _, _ in DIRECTIONS] == CALLABLES
    assert all(e in (STUB, ABSENT, LIVE) and i in (STUB, ABSENT, LIVE) for _, e, i in DIRECTIONS)


def test_every_container_class_passes_the_predicate():
    """
    Molecules, queries and reactions all reach the export half.

    Each is a separate entry in the predicate's `isinstance` tuple, so listing only the molecule would
    send queries and reactions down the import path instead.
    """
    assert is_container(MoleculeContainer())
    assert is_container(read_smiles('CCO'))
    assert is_container(QueryContainer())
    assert is_container(read_smarts('[C;D2]'))
    assert is_container(ReactionContainer([read_smiles('CCO')], [read_smiles('CC=O')]))


@mark.parametrize('value', ['CCO', 42, None, b'CCO', object(), ['CCO']])
def test_non_containers_are_not_containers(value):
    """
    Nothing that is not a chython container may pass the predicate.

    A string is in the list on purpose: `iupac` takes one on its import side.
    """
    assert not is_container(value)


@mark.parametrize('fn', STUBBED_EXPORT, ids=lambda f: f.__name__)
def test_export_direction_is_reached(fn):
    """
    A container reaches the export half of every callable whose export half is still a stub.

    The stub's message naming its direction is what makes this checkable; otherwise dispatch could be
    inverted with nothing red.
    """
    with raises(DirectionNotImplemented, match='export'):
        fn(read_smiles('CCO'))


@mark.parametrize('fn', STUBBED_IMPORT, ids=lambda f: f.__name__)
def test_import_direction_is_reached(fn):
    """
    A non-container reaches the import half of every callable whose import half is not built.

    `ABSENT` raises the same class as `STUB`, permanently rather than temporarily.
    """
    with raises(DirectionNotImplemented):
        fn(object())


def test_cdpkit_import_says_it_is_permanent():
    """CDPKit's refusal must read as a decision, not as an unfinished stub."""
    with raises(DirectionNotImplemented, match='export only'):
        cdpkit(object())


@mark.parametrize('fn', CALLABLES, ids=lambda f: f.__name__)
def test_the_argument_is_positional_only(fn):
    """
    No callable here may accept its subject by keyword.

    Any name would be the toolkit's or chython's, and wrong in one of the two directions.
    """
    with raises(TypeError):
        fn(x=read_smiles('CCO'))


def test_the_container_methods_are_these_very_callables():
    """`mol.to_rdkit()` is the export half of `interop.rdkit` and cannot become a second reading of it.

    Asserted by identity on what the core holds, so it is checkable with no toolkit installed: the
    hook was handed these six functions, so a method's answer is a dispatcher's answer by construction.
    A converter that grew a method-only keyword, or a method registered under the wrong name, fails
    here rather than at whichever toolkit is on the machine that day.
    """
    from chython.core import _core

    for fn in CALLABLES:
        assert _core._reaction_interop_fn(fn.__name__) is fn


def test_the_predicate_imports_no_toolkit():
    """
    Choosing a direction must not import a toolkit.

    Asked of `is_container` directly rather than inferred from a stub, so it stays measurable after the
    last stub is gone.
    """
    from sys import modules

    before = set(modules)
    for value in (read_smiles('CCO'), MoleculeContainer(), 'ethanol', 42, None, object(), b'CCO'):
        is_container(value)
    assert not {m.split('.')[0] for m in set(modules) - before} & TOOLKIT_ROOTS


def test_importing_chython_starts_no_jvm_and_loads_no_toolkit():
    """
    `import chython` alone starts no JVM and loads no toolkit.

    In a subprocess: by the time any test runs, the converter suites have already pulled a toolkit in.
    """
    from subprocess import run
    from sys import executable

    probe = ('import sys, chython\n'
             f'roots = {TOOLKIT_ROOTS!r}\n'
             'print(",".join(sorted({m.split(".")[0] for m in sys.modules} & roots)))\n')
    done = run([executable, '-c', probe], capture_output=True, text=True)
    assert not done.returncode, done.stderr
    assert not done.stdout.strip(), f'import chython pulled in: {done.stdout.strip()}'


@mark.parametrize('fn', STUBBED_EXPORT, ids=lambda f: f.__name__)
def test_no_toolkit_is_imported_on_the_way_to_a_stub(fn):
    """For a half still a stub, the whole path from the call to the raise ran without a toolkit."""
    from sys import modules

    before = set(modules)
    with raises(DirectionNotImplemented):
        fn(read_smiles('CCO'))
    # chython's own modules may be imported lazily here; a toolkit may not be.
    assert not {m.split('.')[0] for m in set(modules) - before} & TOOLKIT_ROOTS
