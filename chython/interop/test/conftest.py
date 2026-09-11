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
Shared fixtures and optional-toolkit skip guards for the toolkit converters.
"""
from importlib.util import find_spec
from pytest import fixture, mark, skip


def _requires(module: str):
    """
    Skip marker for an optional toolkit, testing importability without importing it.

    `find_spec` rather than `importorskip`: importing RDKit or starting a JVM during collection costs
    seconds on every run of the whole suite.
    """
    return mark.skipif(find_spec(module) is None, reason=f'{module} is not installed')


requires_rdkit = _requires('rdkit')
requires_indigo = _requires('indigo')
requires_openbabel = _requires('openbabel')
requires_cdpkit = _requires('CDPL')
requires_jpype = _requires('jpype')
requires_openclatura = _requires('openclatura')

# Inverse guard, for the tests that pin what a missing optional dependency does: where openclatura is
# installed, `to_iupac` walks past its `ImportError` and reports the next gap instead.
absent_openclatura = mark.skipif(find_spec('openclatura') is not None,
                                 reason='openclatura is installed; this test pins what its absence does')


@fixture
def cdk():
    """
    The CDK Java package, or a skip.

    A JVM that starts but cannot find the jar raises `ImportError` from `get_cdk`; that counts as a
    missing toolkit, not a failure.  `CDK_PATH` selects the jar.
    """
    if find_spec('jpype') is None:
        skip('jpype is not installed')
    from .._java import get_cdk

    try:
        return get_cdk()
    except ImportError as e:
        skip(str(e))


@fixture
def opsin():
    """
    The OPSIN `NameToStructure` instance, or a skip.  `OPSIN_PATH` selects the jar.
    """
    if find_spec('jpype') is None:
        skip('jpype is not installed')
    from .._java import get_opsin

    try:
        return get_opsin()
    except ImportError as e:
        skip(str(e))
