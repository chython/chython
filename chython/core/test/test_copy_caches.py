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
"""`copy()` may hand a cached answer to the copy, and may not hand it a STALE one.

Two caches travel with a copy because they are functions of the arena and the arena is immutable and
shared: `canonical_bytes` and the canonical SMILES `str()` returns.  Both are stored as a PAIR --
the value and the generation counter it belongs to -- and both are read back through the same guard,
`cache is not None and cache_gen == _gen`.

THE HAZARD IS IN THE HANDOFF, NOT IN THE GUARD.  A `copy()` that copies the value and stamps it with
`self._gen` re-validates a row the source itself would reject, because an edit bumps `_gen` and leaves
the cache generation behind:

    m.canonical_bytes       # fill it
    with m.edit(): ...      # invalidate it -- for `m`
    m.copy()                # ...and re-validate it, for the copy

That hands the copy the PRE-EDIT molecule's identity, sworn to be the post-edit arena's, so `==`,
`hash()` and `str()` answer for a molecule that does not exist -- and `split()` inherits the whole
thing, being `copy()` plus deletions.

`split()` is how it reaches a caller: in the SMIRKS patcher a product whose stereocentre the template
has just inverted reports the canonical bytes of the molecule the reaction started from, so the
patcher's structural dedupe key (design requirement N9) collapses two genuinely different outcomes
into one.  That is the third test here, and it is why a cache bug in `copy()` is worth four tests
instead of a one-line assertion.
"""
from chython.core import MoleculeContainer, read_smiles


def test_a_stale_identity_does_not_travel_with_the_copy():
    """The direct measurement: fill the cache, edit, copy, ask.

    The reference is a fresh parse of the edited molecule's own SMILES, so nothing here trusts a
    second cache to check the first one.

    THE COPY IS TAKEN BEFORE THE SOURCE IS ASKED AGAIN, and that ordering is the test.  Reading
    `mol.canonical_bytes` first refills the source's row validly and the copy then inherits a correct
    answer -- which is why the same assertions in the other order pass either way and measure nothing.
    """
    mol = read_smiles('CCO')
    stale = mol.canonical_bytes
    with mol.edit() as e:
        e.set_charge(3, -1)

    copy = mol.copy()
    # after the copy, so nothing the reference needs has refilled the source's row.  `str(mol)` and
    # not `'CC[O-]'`: a charge write does not recompute hydrogens, so the edited oxygen still holds
    # the one it had and the molecule is an ethanol anion with an H on it.
    fresh = read_smiles(str(mol)).canonical_bytes
    assert copy.canonical_bytes == fresh
    assert copy.canonical_bytes != stale
    assert mol.canonical_bytes == fresh, 'the source itself must reject its own stale row'
    assert copy == mol and hash(copy) == hash(mol)


def test_a_stale_smiles_does_not_travel_with_the_copy():
    """The same hazard in the other cache, and the visible one: `str()` answers the old molecule.

    Worth its own test because the two caches are separate fields filled by separate readers -- a guard
    on one and not the other leaves a wrong STRING behind, which is the more likely of the two to be
    believed.  The copy is taken before the source is asked again, for the reason the
    test above gives.
    """
    mol = read_smiles('CCO')
    stale = str(mol)
    with mol.edit() as e:
        e.set_charge(3, -1)

    copy = mol.copy()
    assert str(copy) != stale
    assert str(copy) == str(mol) != stale


def test_a_split_component_reports_its_own_identity():
    """`split()` is `copy()` plus deletions, so it inherits the re-validation.

    The parity flip is the change the stale row hid.  Nothing here reads a SMILES string as an
    identity -- the reference is a fresh parse of the configuration the edit produces.
    """
    mol = read_smiles('O[C@H]1CCC[C@H]1C')      # 2-methylcyclopentan-1-ol, both centres real
    before = mol.canonical_bytes
    work = mol.copy()
    with work.edit() as e:
        e.set_parity(2, 1 if work.parity_of(2) == 2 else 2)

    part = work.split()[0]
    assert part.canonical_bytes != before, 'the component must not report the input it came from'
    assert part.canonical_bytes == read_smiles(str(part)).canonical_bytes


def test_a_valid_cache_still_travels():
    """The negative control, because the cheap fix is to stop copying the caches at all.

    Copying them is the point: the arena is shared and immutable, so the copy's answer is the same
    answer and recomputing it costs a full canonicalisation.  A copy taken with no edit in between
    must therefore not recompute -- measured by identity of the returned objects, which is what a
    cache hit gives and a recompute does not.
    """
    mol = read_smiles('c1ccccc1C(=O)O')
    identity = mol.canonical_bytes
    text = str(mol)

    copy = mol.copy()
    assert copy.canonical_bytes is identity
    assert str(copy) is text


def test_a_copy_of_a_never_read_molecule_computes_its_own():
    """And the other end of it: an empty cache is not a stale one.

    A copy taken before anything read the source has nothing to inherit, which must be an ordinary
    cache miss rather than an inherited `None` treated as an answer.
    """
    mol = MoleculeContainer()
    with mol.edit() as e:
        e.add_atom(6)
        e.add_atom(8)
        e.add_bond(1, 2, 1)

    copy = mol.copy()
    assert copy.canonical_bytes == mol.canonical_bytes
    assert str(copy) == str(mol) == 'CO'
