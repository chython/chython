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
"""One `Log` per container and it is THE storage: what a reader recovered, what a pass repaired, what an
edit session lost.  `sgroup_log` and `cip_log` are views over it, not second storages.  No pass takes a
`log=` -- `mol.log` is where every record goes, whether or not anyone asked, and `KekuleResult.log` is a
copy of what its call put there rather than the only place it can be read."""
from chython.core import Log, ReactionContainer, read_smiles as smiles


def test_log_is_created_on_first_access():
    mol = smiles('CCO')
    assert isinstance(mol.log, Log) and len(mol.log) == 0


def test_both_containers_have_one():
    rxn = ReactionContainer([smiles('CCO')], [smiles('CC=O')])
    assert isinstance(rxn.log, Log)


def test_sgroup_loss_lands_in_the_one_log_and_in_the_view():
    """A `DAT` S-group on an atom that is then deleted.  One event, two ways to read it."""
    mol = smiles('CCO')
    n = mol.number_of(0)
    mol.set_sgroups([{'type': b'DAT', 'name': b'BATCH', 'atoms': (n,), 'data': [b'lot-42']}])
    with mol.edit() as e:
        e.delete_atom(n)
    assert mol.sgroup_log == ('1 sgroup record(s) lost a reference to a deleted atom',)
    assert [str(x) for x in mol.log.by_stage('edit:sgroup')] == list(mol.sgroup_log)
    assert mol.log.lost(), 'a lost reference is a LOST record, not an INFO one'


def test_kekule_records_on_the_molecule_as_well_as_on_its_result():
    """The result object is a convenience; the molecule is the storage.

    Pyridine written with the hydrogen an aromatic N cannot have, so `kekule()` repairs and says so;
    benzene kekulises silently and would have asserted nothing.  Nobody asked for a log here, which is
    the point -- `mol.log` is filled anyway, and `result.log` holds the same records.
    """
    mol = smiles('c1cc[nH]cc1')
    result = mol.kekule()
    assert result.log, 'the fixture has to produce a record for the rest of this to mean anything'
    assert all(x.rule.startswith('kekule:') for x in result.log)
    assert [str(x) for x in mol.log] == [str(x) for x in result.log]
    assert {x.stage for x in mol.log} == {'kekule'}, 'the pass names the stage it wrote in'
    assert mol.cip_log == () and mol.sgroup_log == ()


def test_thiele_records_on_the_molecule_as_well_as_on_its_result():
    mol = smiles('C1=CC=CC=CN1')   # 1H-azepine: a Kekule match at 8 pi, so thiele() declines and says so
    result = mol.thiele()
    assert result.log, 'the fixture has to produce a record for the rest of this to mean anything'
    assert all(x.rule.startswith('thiele:') for x in result.log)
    assert [str(x) for x in mol.log] == [str(x) for x in result.log]
    assert {x.stage for x in mol.log} == {'thiele'}


def test_a_second_kekulisation_repairs_nothing():
    """Every repair the kekuliser writes is triggered by an aromatic feature the first pass removed,
    and `thiele()` only re-aromatises rings that already had a Kekule form -- so `canonicalize()`'s
    fixed-point loop cannot grow a duplicate repair however many rounds it runs.  A comment would have
    asserted this; this asserts it.

    Each round does rewrite the representation, and each rewrite says so: `kekule:kekulized` and
    `thiele:aromatized` are per call and are what the round-tripping below is counted against."""
    routine = ('kekule:kekulized', 'thiele:aromatized')
    for smi in ('c1ccc-c1', 'c1cc[nH]cc1', 'c1cc[n+]([O-])cc1', 'Oc1[nH]cnc2nncc1-2', 'c1ccc2c(c1)cccc2'):
        mol = smiles(smi)
        mol.kekule()
        before = [str(x) for x in mol.log if x.rule not in routine]
        for _ in range(3):
            mol.thiele()
            assert not [x for x in mol.kekule().log if x.rule not in routine], smi
            assert [str(x) for x in mol.log if x.rule not in routine] == before, smi


def test_copy_leaves_the_log_behind():
    """`cip_log` is per handle -- `test_cip_storage.py:306` states why -- and it lives here now."""
    mol = smiles('CCO')
    mol.log.record('read something')
    assert len(mol.copy().log) == 0
