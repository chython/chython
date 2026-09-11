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
"""The conformance harness measured against itself.

A matrix that passes because it asserts nothing is worse than no matrix, so the four controls here are
about the harness and not about a format: every toolkit column exists, an absent toolkit skips instead of
failing, a present one with no version fails instead of skipping, and a file no reader can parse comes back
`REFUSED` with the bytes in the message.
"""
from pytest import mark, raises, skip

from . import oracles
from .test_conformance import ABSENT, Cell, REFUSED, classify_outbound


#: A V3000 CTAB cut off inside the atom block.  No reader can complete it, which is the point.
TRUNCATED_V3000 = '''
  chython

  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 7 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0 0 0 0
M  V30 2 C 0 -0.825 0 0
'''


def test_the_probe_set_is_the_same_in_all_five_toolkits():
    """A probe missing from one dict would drop that toolkit's whole column without a failure."""
    for toolkit, probes in oracles.PROBES.items():
        assert set(probes) == set(oracles.FEATURE_NAMES), \
            f'{toolkit} probes differ: {sorted(set(probes) ^ set(oracles.FEATURE_NAMES))}'
    assert len(oracles.FEATURE_NAMES) == len(set(oracles.FEATURE_NAMES))
    assert set(oracles.FEATURE_FORMATS) == set(oracles.FEATURE_NAMES)


@mark.parametrize('toolkit', oracles.TOOLKITS)
def test_an_absent_oracle_skips_and_does_not_fail(toolkit, monkeypatch):
    """The `core/test/oracle.py` posture: a machine without an oracle still gets a green suite."""
    monkeypatch.setitem(oracles._PROBES, toolkit, lambda: None)
    with raises(BaseException) as caught:
        oracles.require(toolkit)
    assert caught.typename == 'Skipped', f'an absent {toolkit} raised {caught.typename}'


def test_the_out_of_tree_probes_report_absent_when_the_path_is_not_there(monkeypatch):
    """The two oracles that are files and not imports, probed through their uncached bodies.

    `__wrapped__` rather than `cache_clear`: the cached answer belongs to the rest of the session, and
    a test that invalidates it makes every later cell re-probe.
    """
    monkeypatch.setattr(oracles, 'MOLCONVERT', oracles.Path('/nonexistent/molconvert'))
    monkeypatch.setattr(oracles, 'CDK_JAR', oracles.Path('/nonexistent/cdk.jar'))
    assert oracles.marvin_version.__wrapped__() is None
    assert oracles.cdk.__wrapped__() is None


def test_every_present_oracle_reports_a_version():
    """A provisioned oracle that states no version is a broken provision, and that fails."""
    found = {name: version for name, version in oracles.versions().items() if version is not None}
    if not found:
        skip('no oracle on this machine; the matrix skips entirely')
    for name, version in found.items():
        assert version.strip(), f'{name} is present and states an empty version'


def test_warming_the_caches_does_not_raise_when_the_marvin_binary_is_absent(monkeypatch):
    """`warm` runs from a session fixture, so an exception there is not one skip but every cell.

    MEASURED: the batch prewarm called `molconvert` unconditionally, and with `MARVIN_BIN` pointing at
    nothing `pytest chython/formats/test/` reported 927 errors on a machine that simply has no Marvin.
    """
    monkeypatch.setattr(oracles, 'MOLCONVERT', oracles.Path('/nonexistent/molconvert'))
    monkeypatch.setattr(oracles, 'marvin_version', lambda: None)
    # three seeds and not one: the write batch skips a batch of fewer than two, so a single seed would
    # never reach `molconvert` and the guard could be deleted with this test still green
    monkeypatch.setattr(oracles, 'SEEDS', oracles.SEEDS[:3])
    written, parses = {}, {}
    oracles.warm(written, parses)
    assert not [k for k in parses if k[0] == 'marvin'], 'a Marvin readback without Marvin'


@mark.parametrize('toolkit', ('chython',) + oracles.TOOLKITS)
def test_an_unmapped_reaction_states_no_mapping(toolkit):
    """The control behind the mapping cells.

    A probe that answered `((), True)` for an unmapped record would be truthy, the classifier would take
    it for a stated feature, and the cell would report `recovered` having compared one empty answer with
    another.  So an unmapped RXN must make every mapping probe falsy.
    """
    from ..ctfile import rxn
    from ...core import read_reaction_smiles

    if toolkit != 'chython':
        oracles.require(toolkit)
    text = rxn(read_reaction_smiles('CC(=O)O.OCC>[H+]>CC(=O)OCC.O'), version=3000)
    parsed = oracles.READERS[toolkit](text, 'rxn')
    assert parsed.error is None, f'{toolkit} refused an unmapped RXN: {parsed.error}'
    assert not oracles.PROBES[toolkit]['mapping'](parsed.obj), \
        f'{toolkit} reports a truthy mapping for a file that states none'


@mark.parametrize('toolkit', oracles.TOOLKITS)
def test_a_broken_fixture_is_classified_refused_and_not_recovered(toolkit, monkeypatch):
    """The control that keeps the classifier honest.

    A `classify` that swallowed every exception would report a green matrix having measured nothing, so a
    file cut off mid-CTAB must come back `REFUSED` -- or `ABSENT` for a reader that returns a partial
    container rather than an error -- and never `recovered`, with the bytes in the message either way.
    """
    oracles.require(toolkit)
    cell = Cell('implicit_h', 'v3000', toolkit, 'OUT', 'aspirin')
    written = {('chython', 'v3000', 'aspirin'): TRUNCATED_V3000}
    outcome, evidence = classify_outbound(cell, written, {})
    assert outcome in (REFUSED, ABSENT), f'a truncated CTAB classified {outcome}'
    assert 'M  V30 BEGIN CTAB' in evidence, 'the message must carry the bytes that produced it'
