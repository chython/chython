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
"""Guards and a write cache for the conformance matrix.

The four oracle fixtures do one thing: skip the test when the toolkit is not installed on this machine,
so a checkout without Marvin still runs green.  `written` is the cache that keeps the matrix affordable
-- one fixture file is written once and read by five readers, and a writer that refused is remembered as
its exception so a repeated cell does not re-run a `molconvert` that already failed.

Nothing here builds a molecule from a string to stand in for a file: a reader is tested on bytes a
reference implementation wrote.
"""
from pytest import fixture

from . import oracles


@fixture(scope='session')
def rdkit_oracle():
    oracles.require('rdkit')


@fixture(scope='session')
def indigo_oracle():
    oracles.require('indigo')


@fixture(scope='session')
def cdk_oracle():
    oracles.require('cdk')


@fixture(scope='session')
def marvin_oracle():
    oracles.require('marvin')


@fixture(scope='session')
def oracle_versions():
    """`{toolkit: version or None}`, for the report header."""
    return oracles.versions()


@fixture(scope='session')
def written():
    """`{(toolkit, fmt, seed_key): text | Exception}` -- filled on demand by `oracles.write_once`."""
    return {}


@fixture(scope='session')
def parses(written):
    """`{(toolkit, fmt, hash(text)): Parsed}` -- filled on demand by `oracles.parse_once`.

    `oracles.warm` fills `written` first and then prebatches the Marvin readbacks over it: one
    `molconvert` per format instead of one per record, which is the difference between about a hundred
    JVM starts and about eight.
    """
    cache = {}
    oracles.warm(written, cache)
    return cache
