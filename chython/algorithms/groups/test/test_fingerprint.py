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
from chython import smiles
from chython.algorithms.groups import fingerprint_schema, fingerprint_size
from chython.algorithms.groups._functional import rules as functional_rules
from chython.algorithms.groups._protective import rules as protective_rules


_PG_OFFSET = 35


def _bits(mol):
    """Set of bit indices raised in the molecule's fingerprint."""
    fp = mol.functional_fingerprint
    return {i for i in range(fingerprint_size() * 8) if fp[i >> 3] >> (i & 7) & 1}


def _mol(smi):
    m = smiles(smi)
    m.canonicalize()
    return m


def test_names_disjoint():
    # protective and functional names share one flat schema namespace, so a
    # collision would silently drop a bit.
    assert not (set(functional_rules) & set(protective_rules))


def test_schema_matches_vector():
    schema = fingerprint_schema()
    total = _PG_OFFSET + len(protective_rules) + len(functional_rules)

    assert len(schema) == total  # 35 descriptors + every rule, all named
    assert sorted(schema.values()) == list(range(total))  # unique, contiguous
    assert max(schema.values()) < fingerprint_size() * 8  # fits the vector


def test_schema_block_order():
    schema = fingerprint_schema()
    fg_offset = _PG_OFFSET + len(protective_rules)

    for offset, name in enumerate(protective_rules):
        assert schema[name] == _PG_OFFSET + offset
    for offset, name in enumerate(functional_rules):
        assert schema[name] == fg_offset + offset


def test_fingerprint_width():
    fp = _mol('CCCCO').functional_fingerprint
    assert isinstance(fp, bytes)
    assert len(fp) == fingerprint_size()


def test_protective_bit_set():
    schema = fingerprint_schema()
    mol = _mol('CC(C)(C)OC(=O)NCCc1ccccc1')  # Boc-protected phenethylamine
    fg_offset = _PG_OFFSET + len(protective_rules)

    assert mol.protective_groups == {'amine_boc': 1}
    pg_bits = {i for i in _bits(mol) if _PG_OFFSET <= i < fg_offset}
    assert pg_bits == {schema['amine_boc']}


def test_protective_bits_track_property():
    # bits must agree with the property, which filters overlapping sub-patterns
    schema = fingerprint_schema()
    mol = _mol('CC(C)(C)OC(=O)NCCO[Si](C)(C)C(C)(C)C')  # Boc amine + TBS ether
    fg_offset = _PG_OFFSET + len(protective_rules)

    pg_bits = {i for i in _bits(mol) if _PG_OFFSET <= i < fg_offset}
    assert pg_bits == {schema[n] for n in mol.protective_groups}


def test_no_protective_bits_when_unprotected():
    mol = _mol('CCCCO')
    fg_offset = _PG_OFFSET + len(protective_rules)

    assert not mol.protective_groups
    assert not [i for i in _bits(mol) if _PG_OFFSET <= i < fg_offset]


def test_functional_bits_at_new_offsets():
    schema = fingerprint_schema()
    mol = _mol('CCCCO')
    fg_offset = _PG_OFFSET + len(protective_rules)

    assert 'primary_alcohol' in mol.functional_groups
    fg_bits = {i for i in _bits(mol) if i >= fg_offset}
    assert fg_bits == {schema[n] for n in mol.functional_groups}
    assert schema['primary_alcohol'] in fg_bits


def test_descriptor_bits():
    # benzene: 1 ring, no hba/hbd, no sp3 carbon, MW 78, no rotatable bonds
    mol = _mol('c1ccccc1')
    assert {i for i in _bits(mol) if i < _PG_OFFSET} == {1, 6, 12, 18, 23, 29}
