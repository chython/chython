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
"""Stereogenicity reads the constitution, never the drawing (`mark_stereogenic`).

| rule | refused | kept |
| --- | --- | --- |
| a Kekule phase is not constitution | fluorene C=N/C=C, benzo rings in opposite phases | 2-bromofluorenone oxime |
| an axis needs both ends unsymmetric | N-aryl imide, dihydropyridine axes and their partners | BINOL, nitro biaryl |
| an achiral cage's phase states nothing | norbornane, tropinone, adamantane bridgeheads | camphor, 2-norbornanol |
"""
import pytest

from chython.core import read_smiles


def _units(text):
    m = read_smiles(text)
    return sorted((u['kind'], u['owners']) for u in m.stereo_units() if u['stereogenic'])


@pytest.mark.parametrize('name,kekule,aromatic', [
    ('fluoren-9-one oxime, mirror phases', 'ON=C1C2=CC=CC=C2C2=CC=CC=C12', 'ON=C1c2ccccc2-c2ccccc12'),
    ('fluoren-9-one oxime, opposite phases', 'ON=C1C2=CC=CC=C2C2=C1C=CC=C2', 'ON=C1c2ccccc2-c2ccccc12'),
    ('9-benzylidenefluorene, opposite phases', 'C(=C1C2=CC=CC=C2C2=C1C=CC=C2)C1=CC=CC=C1',
     'C(=C1c2ccccc2-c2ccccc12)c1ccccc1'),
])
def test_the_kekule_phase_does_not_make_a_double_bond_stereogenic(name, kekule, aromatic):
    """A benzenoid ring's two Kekule phases are one constitution, so fluorene's halves stay equivalent."""
    assert _units(kekule) == _units(aromatic) == [], name
    drawn = read_smiles(kekule)
    drawn.thiele()
    assert drawn == read_smiles(aromatic), 'the two spellings are one compound'


@pytest.mark.parametrize('text', ['ON=C1C2=CC=C(Br)C=C2C2=C1C=CC=C2', 'ON=C1C2=C(C=C(Br)C=C2)C2=CC=CC=C12',
                                  'ON=C1c2ccccc2-c2ccc(Br)cc12'])
def test_a_substituted_fluorenone_oxime_keeps_its_e_z_in_every_spelling(text):
    assert [k for k, _ in _units(text)] == [1], text


@pytest.mark.parametrize('name,text', [
    ('N-(o-tolyl)-4-methylglutarimide', 'CC1CC(=O)N(c2ccccc2C)C(=O)C1'),
    ('N-(2-chlorophenyl)thiomorpholine-3,5-dione S-oxide', 'O=C1CS(=O)CC(=O)N1c1ccccc1Cl'),
    ('N-(2-nitrophenyl)-2,6-dimethyl-1,4-dihydropyridin-4-ylidene',
     'N#CC(C(=O)OC)=C1C=C(C)N(c2ccccc2[N+](=O)[O-])C(C)=C1'),
])
def test_an_axis_with_a_symmetric_end_is_no_unit_and_constrains_nothing(name, text):
    """The imide's two C=O, the dihydropyridine's two methyl-bearing halves: the end maps onto itself with
    its ortho pair exchanged, so there is no axis, and the ring's centre or exocyclic C=C stays equivalent."""
    assert _units(text) == [], name


@pytest.mark.parametrize('name,text,kinds', [
    ('BINOL', 'Oc1ccc2ccccc2c1-c1c(O)ccc2ccccc12', [3]),
    ("2,2'-dimethyl-6-nitrobiphenyl", 'Cc1cccc([N+](=O)[O-])c1-c1ccccc1C', [3]),
    ('4-methylcyclohexylideneacetic acid', 'CC1CCC(=CC(=O)O)CC1', [0, 1]),
    ('1-(o-tolyl)-4-tert-butylcyclohexane', 'Cc1ccccc1C1CCC(C(C)(C)C)CC1', [0, 0]),
])
def test_real_axes_and_ring_cis_trans_pairs_stay(name, text, kinds):
    assert sorted(k for k, _ in _units(text)) == kinds, name


@pytest.mark.parametrize('name,text', [
    ('norbornane', 'C1CC2CCC1C2'), ('bicyclo[2.2.2]octane', 'C1CC2CCC1CC2'),
    ('bicyclo[3.2.2]nonane', 'C1CC2CCC(C1)CC2'), ('tropinone', 'CN1C2CCC1CC(=O)C2'),
    ('adamantane', 'C1C2CC3CC1CC(C2)C3'), ('adamantan-2-ol', 'OC1C2CC3CC(C2)CC1C3'),
])
def test_the_bridgeheads_of_an_achiral_cage_are_not_stereogenic(name, text):
    """The cage realises one relative configuration, and an automorphism carries it onto its inverse: the
    set's phase, its one freedom, states nothing.  Tropinone's mirror exchanges the two bridgeheads."""
    assert _units(text) == [], name


@pytest.mark.parametrize('name,text', [('camphor', 'CC12CCC(C1(C)C)CC2=O'), ('2-norbornanol', 'OC1CC2CCC1C2')])
def test_the_bridgeheads_of_a_chiral_cage_stay(name, text):
    assert len([o for k, o in _units(text) if k == 0]) >= 2, name


def test_tropine_and_pseudotropine_stay_apart():
    """C3's configuration against the bridge is what tells them apart; the cage's mirror inverts C3 too."""
    tropine, pseudo = 'CN1[C@@H]2CC[C@H]1C[C@@H](O)C2', 'CN1[C@@H]2CC[C@H]1C[C@H](O)C2'
    assert read_smiles(tropine) != read_smiles(pseudo)
    assert len(_units('CN1C2CCC1CC(O)C2')) >= 2


_TWO_CAGES = '''
import sys
sys.path.insert(0, sys.argv[1])
from chython.core import read_smiles
m = read_smiles('N1C2CCC1CN(C2)CC1CC2CCC(C2)C1')
with m.edit() as e:
    for n, p in {2: 1, 5: 2, 10: 2, 12: 1, 15: 2}.items():
        e.set_parity(n, p)
    for k, n in enumerate((2, 5, 10, 12, 15), 1):
        e.set_stereo_group(n, 3, k)
c = m.copy()
c.canonicalize_stereo_groups()
print(str(c))
'''


def test_two_cages_in_singleton_and_groups_do_not_recurse():
    """3-(norbornan-2-ylmethyl)-3,8-diazabicyclo[3.2.1]octane, every signed centre its own AND group, built in
    this atom order: a combination's view asked the meaning about the same combination, which asked for its
    view again until the stack ran out.  Run in a child, so a regression fails here rather than killing pytest."""
    import subprocess
    import sys
    from pathlib import Path
    root = str(Path(__file__).resolve().parents[3])          # the tree holding this chython package
    run = subprocess.run([sys.executable, '-c', _TWO_CAGES, root], capture_output=True, text=True, timeout=120)
    assert run.returncode == 0, run.stderr[-500:]
    assert run.stdout.strip() == 'C1CC2CN(C[C@H]3CC4CCC(C4)C3)CC1N2 |&1:6|'
