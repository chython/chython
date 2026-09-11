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
"""Asking for the aromatic spelling must not silently choose a tautomer.

`kekule()` repairs, because it is the boundary where input arrives; `thiele()` must not, because by
then the data is the library's own arena.  The hydrogen move chython 2 did inside `thiele` is owed to
`standardize_isomers` instead.  Both halves are pinned here.
"""
import pytest

from chython.chemistry.test._oracle import ask, needs_oracle
from chython.core import read_smiles
import chython.chemistry                                                         # noqa: F401


#: 4,7-dihydro-pyrrolo-pyridine drawn with the hydrogen on the wrong nitrogen.  chython 2's own
#: fixture (`algorithms/aromatics/test/test_thiele.py::test_tautomer_fix`), so the divergence is
#: anchored to their test rather than to ours.
SHIFTED = 'N1C=CC2=NC=CC2=C1'

#: `{atom: total hydrogens}` as written: atom 1 carries the NH, atom 5 is the bare aromatic N.
AS_WRITTEN = {1: 1, 2: 1, 3: 1, 4: 0, 5: 0, 6: 1, 7: 1, 8: 0, 9: 1}


def test_thiele_aromatises_the_condensed_ring():
    """The part `thiele` is for.  Without this the test below passes for the wrong reason."""
    molecule = read_smiles(SHIFTED)
    result = molecule.thiele()
    assert result.changed, 'nothing was aromatised, so the hydrogen check proves nothing'
    assert molecule.smiles == 'c12ccnc2cc[nH]c1', molecule.smiles


def test_thiele_does_not_move_a_hydrogen():
    """The ruling, as an assertion: every hydrogen count survives aromatisation unchanged."""
    molecule = read_smiles(SHIFTED)
    assert {n: molecule.total_h_of(n) for n in molecule.atom_numbers} == AS_WRITTEN
    molecule.thiele()
    assert {n: molecule.total_h_of(n) for n in molecule.atom_numbers} == AS_WRITTEN, (
        'thiele moved a hydrogen.  Aromatisation is a change of spelling; moving a hydrogen is a '
        'change of molecule, and belongs in standardize_isomers where a caller asks for it')


def test_thiele_reports_no_tautomer_decision_in_its_result():
    """`ThieleResult` has nowhere to log a hydrogen move.

    What it does report is the aromatisation itself, one line per ring system it accepted; nothing else,
    and `refused` comes back empty.
    """
    result = read_smiles(SHIFTED).thiele()
    assert [r.rule for r in result.log] == ['thiele:aromatized']
    assert result.refused == []


@needs_oracle
def test_chython_2_moved_the_hydrogen_and_that_is_the_defect():
    """Pins the divergence at its source, so nobody 'fixes' V3 back into chython 2's behaviour.

    Triaged: V3 is right.  Both `fix_tautomers` branches are asserted, since the flag being the only
    difference is what makes the hydrogen move a hidden second job.
    """
    script = """
from chython import smiles

for line in sys.stdin:
    source, _, flag = line.rstrip('\\n').partition('\\t')
    m = smiles(source)
    m.thiele(fix_tautomers=flag == 'on')
    sys.stdout.write('%s\\t%s\\t%s\\n' % (flag, format(m, 's'),
                     ','.join('%d:%d' % (n, a.total_hydrogens) for n, a in m.atoms())))
"""
    answers = dict(line.split('\t', 1) for line in ask(script, f'{SHIFTED}\ton\n{SHIFTED}\toff'))
    on_smiles, on_h = answers['on'].split('\t')
    off_smiles, off_h = answers['off'].split('\t')

    as_written = ','.join(f'{n}:{h}' for n, h in AS_WRITTEN.items())
    assert off_h == as_written, 'chython 2 leaves the hydrogens alone with the flag off'
    assert on_h != as_written, (
        'chython 2 no longer moves the hydrogen, so this divergence is stale -- re-triage it rather '
        'than deleting the test')
    assert on_smiles != off_smiles


@pytest.mark.xfail(strict=True, reason='owed: standardize_isomers has not been ported yet, so the '
                                       'hydrogen move chython 2 performed inside thiele is '
                                       'currently unavailable anywhere.  Delete this marker when '
                                       'the isomer pass lands -- a strict xfail is what makes that '
                                       'a required edit rather than an optional one')
def test_the_hydrogen_move_is_owed_to_standardize_isomers():
    """What the split costs until the isomer pass exists, written down as a failing test.

    The end state is the one chython 2 reached by combining the two jobs, and is what the ported rule
    must reproduce: the hydrogen on atom 5, not atom 1.
    """
    molecule = read_smiles(SHIFTED)
    molecule.thiele()
    molecule.standardize_isomers()
    assert molecule.total_h_of(1) == 0
    assert molecule.total_h_of(5) == 1


@pytest.mark.xfail(strict=True, reason='the core aromatises the six-membered ring and leaves the '
                                       'fused five-membered one in its Kekule form while reporting '
                                       'refused==[], so a declined system is invisible.  For the '
                                       'core epic: either apply the rule-based decision chython 2 '
                                       'makes for these, or report the refusal')
def test_a_declined_five_membered_ring_is_reported():
    """A partial aromatisation must say which system it gave up on.  Measured: it does not.

    chython 2 handles this shape -- a five-membered ring with three sp2 atoms fused to an aromatic
    ring -- unconditionally, with no hydrogen moved, so it is an aromatisation decision and belongs
    in `thiele` rather than in the isomer pass.
    """
    molecule = read_smiles('N1C=Cn2cccc12')
    molecule.kekule()
    result = molecule.thiele()
    assert result.refused or molecule.smiles == 'n12cc[nH]c1ccc2', (
        f'{molecule.smiles} is half aromatic and nothing was reported')
