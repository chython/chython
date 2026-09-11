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
"""`fix_resonance`: what it fixes, what it refuses, and that it terminates deterministically.

Every molecule below is a textbook one -- formamide, butadiene, ethylene, o-xylylene, malononitrile,
tetramethylammonium, methyl azide, tetramethylborate, a sulfilimine.
"""
import subprocess
import sys

import pytest

from .. import fix_resonance
from .._resonance import alternating_paths
from ...core import read_smiles


def run(smiles):
    """Apply the pass and hand back `(changed, canonical SMILES, log)`."""
    molecule = read_smiles(smiles)
    log = molecule.log
    changed = fix_resonance(molecule)
    return changed, format(molecule, ''), log


def rules(log):
    return {record.rule for record in log}


def test_dipole_is_neutralised():
    """N-methylformamide, drawn as its zwitterion.  The textbook case the pass exists for."""
    changed, result, log = run('[O-]C=[NH+]C')
    assert changed
    assert result == format(read_smiles('O=CNC'), '')
    assert rules(log) == {'resonance:donor_anion'}


def test_biradical_is_paired():
    """Two radicals four atoms apart become 1,3-butadiene."""
    changed, result, log = run('[CH2]C=C[CH2] |^1:0,3|')
    assert changed
    assert result == format(read_smiles('C=CC=C'), '')
    assert rules(log) == {'resonance:radical'}


def test_adjacent_biradical_is_paired():
    """The shortest possible case: two methyl radicals bonded to each other are ethylene."""
    changed, result, log = run('[CH2][CH2] |^1:0,1|')
    assert changed
    assert result == format(read_smiles('C=C'), '')


def test_amine_beside_a_cation_becomes_an_iminium():
    """An sp3 amine donates even though it is neutral."""
    changed, result, _ = run('N(C)(C)C=C[CH2+]')
    assert changed
    assert result == format(read_smiles('C=CC=[N+](C)C'), '')


def test_nitrile_takes_a_remote_carbanion_charge():
    """`N#X-[X-] >> [N-]=X=X`; a charge is better on nitrogen."""
    changed, result, _ = run('[CH2-]C#N')
    assert changed
    assert result == format(read_smiles('C=C=[N-]'), '')


def test_neutral_molecule_is_untouched():
    changed, result, log = run('CC(=O)NC')
    assert not changed
    assert result == format(read_smiles('CC(=O)NC'), '')
    assert not log


def test_nitro_group_is_untouched():
    """A nitro group's dipole is the correct drawing: the nitrogen cannot take a fourth bond."""
    changed, _, _ = run('C[N+](=O)[O-]')
    assert not changed


def test_aromatic_ring_never_gains_a_triple_bond():
    """o-xylylene as a biradical: no ring bond becomes a triple bond.

    Asserted on the bond orders, not on the return value: a pass that "changed nothing" while having
    written a triple bond into a benzene ring would satisfy `not changed` on some other numbering.
    """
    molecule = read_smiles('[CH2]c1ccccc1[CH2] |^1:0,7|')
    before = {(min(n, m), max(n, m)): molecule.order_of(n, m)
              for n in molecule.atom_numbers for m in molecule.neighbors_of(n)}
    log = molecule.log
    changed = fix_resonance(molecule)
    after = {(min(n, m), max(n, m)): molecule.order_of(n, m)
             for n in molecule.atom_numbers for m in molecule.neighbors_of(n)}

    assert 3 not in after.values(), 'a triple bond appeared; this is chython 2 defect 1'
    assert after == before
    assert not changed
    # the refusal is on the record, naming the pattern that offered the endpoint
    assert log
    assert rules(log) == {'resonance:radical'}
    assert all('refused' in record.message or 'aromatic' in record.message for record in log)


def test_radical_on_an_aromatic_atom_is_refused_by_name():
    """The same molecule with a radical drawn on a ring carbon: refused, and the log says why."""
    changed, _, log = run('[CH2]c1ccccc1[CH2] |^1:1,7|')
    assert not changed
    assert any('aromatic' in record.message for record in log)
    assert rules(log) == {'resonance:radical'}


def test_walk_refuses_to_cross_an_aromatic_bond():
    """The guard is in the walk, not only in the endpoint patterns.

    Every atom of benzene is handed to the walker as allowed and as a target -- what `_classify`
    refuses to do -- so a walk that did not stop at an order-4 bond would find paths.
    """
    benzene = read_smiles('c1ccccc1')
    allowed = frozenset(benzene.atom_numbers)
    targets = frozenset(n for n in benzene.atom_numbers if n != 1)
    assert not list(alternating_paths(benzene, 1, targets, allowed))


def test_walk_refuses_to_cross_a_dative_bond():
    """Order 8 is outside valence bookkeeping everywhere else, and it is impassable here too."""
    molecule = read_smiles('C[N]->[Zn]')
    assert molecule.order_of(2, 3) == 8, 'the probe stopped containing a dative bond'
    allowed = frozenset(molecule.atom_numbers)
    paths = list(alternating_paths(molecule, 1, frozenset({2, 3}), allowed))
    # the nitrogen is reachable, the zinc behind the dative bond is not
    assert paths
    assert all(path[-1][1] == 2 for path in paths)
    assert all(molecule.order_of(u, v) != 8 for path in paths for u, v, _ in path)


def test_anion_end_valence_is_checked():
    """The anion end is validated, not only the cation end.

    Neutralising this dipole would leave a neutral four-coordinate boron with a bond order sum of 4,
    which the valence collection has no rule for at all.
    """
    changed, result, log = run('C[B-](C)C=[NH+]C')
    assert not changed
    assert result == format(read_smiles('C[B-](C)C=[NH+]C'), '')
    assert log
    record, = log
    assert record.rule == 'resonance:donor_anion'
    assert 'no rule for at all' in record.message
    # the refused atom is the anion end
    assert 2 in record.atoms


def test_interior_valence_is_checked():
    """An interior atom is validated too.

    The sulfur keeps its charge, radical state and bond order sum -- only its environment moves, from
    two double bonds to one single and one triple -- so only a rule keyed on the environment refuses it.
    """
    changed, result, log = run('[CH2-]C=S=C=[NH2+]')
    assert not changed
    assert result == format(read_smiles('[CH2-]C=S=C=[NH2+]'), '')
    assert log
    assert any('environment' in record.message for record in log)
    # the refused atom is the interior sulfur, not either end
    atoms = {n for record in log for n in record.atoms}
    assert 3 in atoms


# the opt-outs: each is a row, and each refusal names it
@pytest.mark.parametrize('smiles,rule', [
    # a quaternary ammonium has no empty orbital; a carbanion three bonds away must not reach it
    ('[CH2-]C=C[N+](C)(C)C', 'resonance:veto_ammonium'),
    # methyl azide: the 1,3-dipole is the correct drawing
    ('CN=[N+]=[N-]', 'resonance:veto_azide'),
    # tetramethylborate: four substituents, no pair to give and no room for a fifth bond
    ('C=C[B-](C)(C)C', 'resonance:veto_borate'),
    ('C[S+]=NC', 'resonance:veto_sulfonium_ylide'),      # an S-methyl sulfilimine
    ('C[BH-](C)C', 'resonance:veto_hydroborate'),
    ('F[P-](F)(F)(F)(F)F', 'resonance:veto_hexacoordinate_p'),
    ('C[P+](C)(C)C', 'resonance:veto_phosphonium'),
    ('[CH2-]C=C[S+](C)C', 'resonance:veto_sulfonium'),
])
def test_opt_out_is_refused_and_named(smiles, rule):
    changed, result, log = run(smiles)
    assert not changed, f'{smiles} was rewritten; {rule} did not hold'
    assert result == format(read_smiles(smiles), '')
    assert rule in rules(log), f'{smiles} produced no record naming {rule}: {log}'
    record = next(r for r in log if r.rule == rule)
    assert 'vetoed by' in record.message
    assert len(record.atoms) == 1


def test_an_atom_no_row_wanted_is_not_logged():
    """A veto is logged only when an accept row also matched.  Silence is not a refusal."""
    _, _, log = run('CC(=O)NC')
    assert not log


def test_no_oscillation():
    """The malononitrile anion: two nitriles and one carbanion, which can loop.

    The first move is accepted -- a charge leaves carbon for nitrogen, strictly decreasing the
    potential -- and the return trip is refused with a record that says so.
    """
    changed, result, log = run('[CH2-](C#N)C#N')
    assert changed
    first = result
    refusals = [r for r in log if 'would not reduce' in r.message]
    assert refusals, f'nothing refused the return trip: {log}'

    # a fixed point: a second pass finds nothing
    molecule = read_smiles(result)
    assert not fix_resonance(molecule)
    assert format(molecule, '') == first


@pytest.mark.parametrize('smiles', [
    '[O-]C=[NH+]C', 'N(C)(C)C=C[CH2+]', '[CH2]C=C[CH2] |^1:0,3|', '[CH2][CH2] |^1:0,1|',
    '[CH2-](C#N)C#N', '[CH2-]C#N', '[CH2]c1ccccc1[CH2] |^1:0,7|', 'C[N+](=O)[O-]',
    'C[B-](C)C=[NH+]C', '[CH2-]C=S=C=[NH2+]', 'CN=[N+]=[N-]', 'C[N+](C)(C)C',
])
def test_idempotent(smiles):
    """The second call changes nothing."""
    molecule = read_smiles(smiles)
    fix_resonance(molecule)
    once = format(molecule, '')
    assert not fix_resonance(molecule)
    assert format(molecule, '') == once


# every case above, so determinism is checked over the whole surface rather than one molecule
CASES = ['[O-]C=[NH+]C', 'N(C)(C)C=C[CH2+]', '[CH2]C=C[CH2] |^1:0,3|', '[CH2][CH2] |^1:0,1|',
         '[CH2-](C#N)C#N', '[CH2-]C#N', '[CH2]c1ccccc1[CH2] |^1:0,7|', 'C[N+](=O)[O-]',
         'C[B-](C)C=[NH+]C', '[CH2-]C=S=C=[NH2+]', 'CN=[N+]=[N-]', 'C[N+](C)(C)C',
         '[CH2-]C=C[N+](C)(C)C', '[CH2-]C=C[S+](C)C', 'C=C[B-](C)(C)C',
         '[CH2]C=CC=C[CH2] |^1:0,5|', '[O-]C=C[NH+]=C']

PROBE = '''
import sys
from chython.core import read_smiles
from chython.chemistry import fix_resonance
for smiles in sys.argv[1:]:
    molecule = read_smiles(smiles)
    log = molecule.log
    changed = fix_resonance(molecule)
    print(changed, format(molecule, ''), [(r.rule, r.atoms) for r in log])
'''


def probe(seed):
    out = subprocess.run([sys.executable, '-c', PROBE, *CASES], capture_output=True, text=True,
                         env={'PYTHONHASHSEED': seed, 'PATH': '/usr/bin:/bin'}, check=True)
    return out.stdout


def test_deterministic_within_a_process():
    """Same molecule, same answer, log included -- twice in one interpreter."""
    for smiles in CASES:
        first = run(smiles)
        second = run(smiles)
        assert first[:2] == second[:2]
        assert [(r.rule, r.atoms, r.message) for r in first[2]] == \
               [(r.rule, r.atoms, r.message) for r in second[2]]


def test_deterministic_across_hash_seeds():
    """And across `PYTHONHASHSEED`, which catches a decision taken off set order: the answer must be a
    function of the atom numbering alone.
    """
    reference = probe('0')
    assert reference.strip()
    for seed in ('1', '2', '12345', 'random'):
        assert probe(seed) == reference, f'PYTHONHASHSEED={seed} gave a different answer'
