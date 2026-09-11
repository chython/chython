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
"""Fixtures for the XML tests.

Files under the repository's ``test/`` carry everything about a *document* -- entity policy, CML 1 and
CML 2 spellings, vendor exports.  Inline strings in the test that uses them carry everything about a
*rule*: one attribute, one child element, one malformation.
"""

from pathlib import Path

from pytest import fixture

from .._tree import available_engines


#: Every fixture this package's tests may open, named rather than globbed so a new file in `test/`
#: cannot silently change what a passing suite means.  Enforced by `data` below.
FIXTURES = ('cml_stereo.cml', 'cml_stereo.mol', 'cml_marvin.cml', 'cml_quirks.cml',
            'cml_damaged.cml', 'cml_entity_bomb.cml', 'cml_external_dtd.cml',
            'implicit.mrv', 'mrv_hydrogens.mrv', 'mrv_reaction.mrv',
            # One molecule, two official writers, two serialisations: 2-methyltetrahydropyran with one
            # `or1` centre from Marvin Sketch v25.1.3 (column form of `<atomArray>`) and from MarvinJS
            # (element form).  Exact bytes, unindented: the bytes are the evidence.
            'mrv_stereo_sketch.mrv', 'mrv_stereo_js.mrv')


def _root():
    """The repository root, found by walking up from this file to the directory holding ``test/``."""
    for parent in Path(__file__).resolve().parents:
        if (parent / 'test').is_dir() and (parent / 'chython').is_dir():
            return parent
    raise RuntimeError('cannot locate the repository root from ' + __file__)


@fixture(scope='session')
def root():
    return _root()


@fixture(scope='session')
def data(root):
    """``f(name) -> Path`` for a file under the repository's ``test/``, checked to exist.

    Asserted rather than skipped: a missing committed fixture is a broken checkout, and a skip would
    turn the entity-bomb test into a silent pass.  `name` must also be declared in :data:`FIXTURES`.
    """
    def path(name):
        assert name in FIXTURES, f'{name} is not in FIXTURES; declare it there before reading it'
        p = root / 'test' / name
        assert p.exists(), f'fixture {name} is missing from {root / "test"}'
        return p
    return path


@fixture(params=available_engines())
def engine(request):
    """Each XML backend in turn, so a test parameterized on this proves the policy and not a backend.

    ``defusedxml`` is defence in depth over a guarantee :mod:`.._tree` makes on its own, so the stdlib
    case must always be measured -- a preferred-backend-only test would pass for the wrong reason.
    """
    return request.param
