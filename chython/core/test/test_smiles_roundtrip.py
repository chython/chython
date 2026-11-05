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
"""The writer and the reader against each other: `write` then `read` then `write` must not move.

A CLOSED loop, where the only two things in the circuit are ours, so a disagreement is provably a bug
in one of them and not a convention difference with a third-party parser.

WHAT THE LOOP CAN AND CANNOT PROVE.  `write_smiles(read_smiles(write_smiles(m, spec)), spec) ==
write_smiles(m, spec)` is a FIXED POINT claim about strings.  It catches a writer emitting something
the reader misreads and a reader building something the writer disagrees with.  It cannot catch a
misunderstanding both share, and it is not chemistry: two different strings can be the same
molecule.

AND THE INVARIANT CAN LIE IN ONE SHAPE, which is why `unknown_h` vetoes rather than merely explains:
the `A` dialect on a molecule whose hydrogen count nobody stated is a string fixed point AND a
changed molecule at the same time.  The last test in this file is that counterexample, worked in
full; the rule it establishes is that a fixed point is evidence and `unknown_h` is the veto.

`first_5K.smi` is 4999 public NCI records shipped inside RDKit, so nothing is fetched; a missing
corpus skips.  It carries no stereo, so the loop over it is a pure constitution loop -- the one that
would catch a ring-closure digit or a bracket rule that only the writer believes in.
"""
from pathlib import Path

from pytest import fixture, mark, skip

from chython.core import MoleculeContainer
from chython.core._core import read_smiles, smw_traversal, write_smiles


# `!s` is in the sweep although the corpus has no stereo: the spec still selects a different seed for
# the canonical order, so it is a different traversal reaching the same answer, and that is worth a
# column.  `r` is not here and never will be -- it raises.
SPECS = ['', 'A', 'i', '!s']


@fixture(scope='module')
def nci():
    """The NCI SMILES strings, read from the RDKit installation."""
    try:
        import rdkit
    except ImportError:
        skip('RDKit is not installed, and the NCI corpus ships inside it')
    path = Path(rdkit.__file__).resolve().parent / 'Data' / 'NCI' / 'first_5K.smi'
    if not path.is_file():
        skip('the RDKit NCI corpus is not installed (looked for %s)' % path)
    with path.open(encoding='utf-8') as f:
        return [line.split()[0] for line in f if line.split()]


def loop(mol, spec=''):
    """`(out, out_again)`. The reader is given OUR string, so a refusal is our bug, not bad input."""
    out = write_smiles(mol, spec)
    return out, write_smiles(read_smiles(out), spec)


# ------------------------------------------------------------------------------------------------
# The constitution loop.  The corpus carries no stereo, so nothing here can pass or fail on a
# configuration.
@mark.parametrize('spec', SPECS)
def test_every_NCI_string_reads_back_to_ITSELF(nci, spec):
    """4999 of 4999, on all four specs, with the CLOSED loop: our reader in, our writer out."""
    moved = []
    refused = []
    total = 0
    for text in nci:
        total += 1
        out = write_smiles(read_smiles(text), spec)
        try:
            again = write_smiles(read_smiles(out), spec)
        except Exception as e:
            refused.append((text, out, '%s: %s' % (type(e).__name__, e)))
            continue
        if out != again:
            moved.append((text, out, again))
    assert total > 4900, 'the corpus shrank: %d records' % total
    assert not refused, '%d strings our own reader refused, e.g. %r' % (len(refused), refused[:3])
    assert not moved, '%d of %d strings moved, e.g. %r' % (len(moved), total, moved[:3])


# ------------------------------------------------------------------------------------------------
# WHERE THE FIXED POINT IS NOT ENOUGH.  The one shape in which this whole file's invariant lies.
def test_a_string_fixed_point_does_NOT_prove_the_molecule_survived():
    """The `A` dialect closes the loop on a benzene whose hydrogen count nobody stated -- and the
    molecule changes anyway.  This is the counterexample that makes `unknown_h` the veto on a fixed
    point.

    The dialect brackets its aromatic atoms so the count is explicit, which is the right call when
    there IS a count.  When there is not, a bracket with no H term states zero -- a number the
    molecule never claimed -- and there is no honest alternative, because a bare upper-case atom on a
    `:` bond has no defined count either.  So the writer does the only thing left: it writes the
    lossy spelling and reports every atom in `unknown_h`.

    The trap is that the loop then LOOKS clean.  Read that string back and the atoms have a stated
    zero, write it in the same dialect and you get the same string, so string-fixed-point says pass.
    The default dialect is what gives the game away: `c1ccccc1` before, `[c]1[c][c][c][c][c]1` after,
    because a stated zero on an aromatic carbon has to be spelled and an unstated count does not.
    Hence the rule this file follows: a fixed point is evidence, `unknown_h` is the veto.
    """
    mol = MoleculeContainer()
    with mol.edit() as e:
        ids = [e.add_atom('C') for _ in range(6)]
        for i in range(6):
            e.add_bond(ids[i], ids[(i + 1) % 6], 4)

    assert smw_traversal(mol)['unknown_h'] == (2, 3, 4, 5, 6, 1), 'the arena stated a count'
    assert write_smiles(mol) == 'c1ccccc1'
    assert write_smiles(mol, 'A') == '[C]:1:[C]:[C]:[C]:[C]:[C]:1'

    out, again = loop(mol, 'A')
    assert out == again, 'the A dialect is a string fixed point here'
    back = read_smiles(out)
    assert smw_traversal(back)['unknown_h'] == (), 'the reader had to invent a count'
    assert write_smiles(back) == '[c]1[c][c][c][c][c]1', \
        'the DEFAULT dialect is what exposes the change the A loop hid'

    # the default dialect, on the same molecule, loses nothing: a bare atom degrades to the
    # valence-derived count, which is the count the caller most likely meant
    out, again = loop(mol)
    assert out == again == 'c1ccccc1'
    assert smw_traversal(read_smiles(out))['unknown_h'] == ()
